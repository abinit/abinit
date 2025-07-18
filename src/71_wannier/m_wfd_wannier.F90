!!****m*ABINIT/m_wfd_wannier
!! NAME
!!  m_wfd_wannier
!!
!! FUNCTION
!!  The high level wfd_t inteface for building wannier functions
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2022 ABINIT group (hexu)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"


!===============================================================
! m_wfd_wannier
!> @description: Wannier function from wfd_t
!===============================================================

module m_wfd_wannier

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use defs_abitypes,     only : mpi_type
  use m_dtset,           only : dataset_type
  use m_hdr,             only : hdr_type
  use m_wfd,             only : wfd_t
  use m_crystal,         only : crystal_t
  use m_pawang,          only : pawang_type
  use m_pawrad,          only : pawrad_type
  use m_pawtab,          only : pawtab_type
  use m_pawcprj,         only : pawcprj_type
  use defs_datatypes,    only : pseudopotential_type
  use m_dtfil,           only : datafiles_type
  use m_abstract_wf,     only : abstract_wf, compute_iwav, write_cg_and_cprj, wann_ksetting_t, init_mywfc
  use m_mlwfovlp,        only : mlwfovlp
  use m_ebands,          only : ebands_t

  use defs_wannier90

  implicit none
  private
  public :: wfd_run_wannier
  integer,  parameter :: master=0

contains


!-----------------------------------------------------------------------------
!> @brief The high level wfd_t inteface for building wannier functions
!> @param[in] cryst: crystal_t type,  crystal structure
!> @param[inout] ebands: ebands_t type,   eigenvalues and eigenvectors
!> @param[inout] hdr: hdr_type type,   header
!> @param[in] mpi_enreg: mpi_type type,   MPI information
!> @param[in] ngfftc: integer,   FFT grid for coarse grid
!> @param[in] ngfftf: integer,   FFT grid for fine grid
!> @param[in] wfd: wfd_t type,   wavefunction data
!> @param[in] dtset: dataset_type type,   dataset
!> @param[in] dtfil: datafiles_type type,   datafiles
!> @param[in] pawang: pawang_type type,   PAW angular momentum
!> @param[in] pawrad: pawrad_type type,   PAW radial functions
!> @param[in] pawtab: pawtab_type type,   PAW tabulated functions
!> @param[in] psps: pseudopotential_type type,   pseudopotential
!> @param[in] kg: real(dp),   FFT grid, shouldn't be used
!> @param[in] cg: complex(dp),   wavefunction coefficients, shouldn't be used
!> @param[in] cprj: complex(dp),   wavefunction coefficients, shouldn't be used
!-----------------------------------------------------------------------------
  subroutine wfd_run_wannier(cryst, ebands, hdr, mpi_enreg, &
       & ngfftc, ngfftf,  wfd, dtset, dtfil,  &
       & pawang,  pawrad, pawtab, psps , kg, cg, cprj)

    type(crystal_t), intent(in) :: cryst
    type(ebands_t), intent(in) :: ebands
    type(hdr_type), intent(in) :: hdr
    integer, intent(in) :: ngfftc(18),ngfftf(18)
    type(dataset_type), intent(in) :: dtset
    type(datafiles_type),intent(in) :: dtfil
    type(mpi_type), intent(inout) :: mpi_enreg
    type(pseudopotential_type),intent(in) :: psps
    type(pawang_type),intent(in) :: pawang
    !type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
    type(pawrad_type),intent(in) :: pawrad(:)
    !type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
    type(pawtab_type),intent(in) :: pawtab(:)

    type(wfd_t), optional, intent(inout) :: wfd
    real(dp), optional, target, intent(in) :: cg(:, :)
    integer, optional, target, intent(in) :: kg(:, :)
    type(pawcprj_type), optional, target, intent(in) :: cprj(:, :)
    class(abstract_wf), pointer :: mywfc
    integer :: mgfftc
    integer :: nfft
    integer :: mcg, mcprj
    !real(dp), pointer:: ptr_cg(:, :)=>null()
    integer, pointer:: ptr_kg(:, :)=>null()
    !type(pawcprj_type), pointer:: ptr_cprj(:, :)=>null()
    integer :: mpw, nspinor, mband, mkmem, nsppol
    ! TODO: anything todo with nkpt: fullBZ
    ! TODO: ebands for fullBZ
    ! TODO: mcprj
    ! TODO: mgfftc: is it ngfft. NO
    ! TODO: ngfft: is it ngfftc, or ngfftf
    ! TODO: check mpw in wfd
    ! TODO: mkmem: note: the mkmem in mpi_enreg is not used.
    ! TODO: gather kg from wfd
    ! TODO: mcg
    ! TODO: mcprj
    integer :: exclude_bands(hdr%mband, hdr%nsppol)

    integer :: spaceComm, nprocs, rank, master

    !print *, "============================================================"
    !print *, "Starting WFD Wannier"
    !print *, "============================================================"

    spaceComm=MPI_enreg%comm_world
    nprocs = MPI_enreg%nproc
    rank= MPI_enreg%me
    master=0

    mgfftc=dtset%mgfft
    mpw=MAXVAL(hdr%npwarr)
    nspinor=hdr%nspinor
    mband=hdr%mband
    nsppol = hdr%nsppol
    !mcg=mpw*nspinor*mband* mkmem *nsppol
    !mcprj = nspinor*mband* mkmem*nsppol
    ! Is this correct for paw?
    nfft=product(ngfftc(1:3))
       if (present(cg)) then
         mkmem= dtset%mkmem
         call init_mywfc(mywfc=mywfc, ebands=ebands, cg=cg, cprj=cprj, &
           &  cryst=cryst,dtset=dtset, &
           & dtfil=dtfil, hdr=hdr, MPI_enreg=MPI_enreg, nprocs=nprocs, psps=psps, &
           & pawtab=pawtab, rank=rank, comm=spaceComm)
         ptr_kg => kg
          ! call mlwfovlp(mywfc=mywfc, crystal=cryst, ebands=ebands, hdr=hdr, atindx1=cryst%atindx1  &
          !      &,dtset=dtset,dtfil=dtfil, &
          !      & eigen=ebands%eig,gprimd=cryst%gprimd,kg=ptr_kg,&
          !      & mband=hdr%mband,mcg=mcg,mcprj=mcprj,mgfftc=mgfftc, &
          !      & mkmem=mkmem,mpi_enreg=mpi_enreg,mpw=mpw,natom=cryst%natom,&
          !      & nattyp=cryst%nattyp,nfft=nfft,ngfft=ngfftf,nkpt=hdr%nkpt,npwarr= hdr%npwarr , &
          !      &nsppol=dtset%nsppol,ntypat=cryst%ntypat,occ=ebands%occ,&
          !      &pawang=pawang,pawrad=pawrad,pawtab=pawtab,prtvol=dtset%prtvol,psps=psps, &
          !      &rprimd=cryst%rprimd,ucvol=cryst%ucvol, xred=cryst%xred)
       else
           call init_mywfc(mywfc=mywfc, ebands=ebands, wfd=wfd, &
             &  cryst=cryst,dtset=dtset, &
             & dtfil=dtfil, hdr=hdr, MPI_enreg=MPI_enreg, nprocs=nprocs, psps=psps, &
             & pawtab=pawtab, rank=rank, comm=spaceComm)

           ! TODO: distribute k
           mkmem = mywfc%kset%my_nkpt
           mpw=maxval(mywfc%hdr%npwarr)
           mcg=mpw*nspinor*mband* mkmem *nsppol
           mcprj = nspinor*mband* mkmem*nsppol
           ! block
           !   integer :: ik
           !   do ik =1, mkmem
           !     print *, ik, ":", mywfc%ebands%kptns(:, ik)
           !   end do
           ! end block

           if (present(kg)) then
             ptr_kg=> kg
           else
             ABI_MALLOC(ptr_kg, (3, mpw*mkmem))
             call mywfc%get_kgs(ptr_kg)
           end if
         end if

         call mlwfovlp(mywfc=mywfc, crystal=cryst, ebands=mywfc%ebands, hdr=mywfc%hdr, atindx1=cryst%atindx1  &
               &,dtset=mywfc%dtset,dtfil=dtfil, &
               & eigen=mywfc%ebands%eig,gprimd=cryst%gprimd,kg=ptr_kg,&
               & mband=mywfc%hdr%mband,mcg=mcg,mcprj=mcprj,mgfftc=mgfftc, &
               & mkmem=mkmem,mpi_enreg=mpi_enreg,mpw=mpw,natom=cryst%natom,&
               & nattyp=cryst%nattyp,nfft=nfft,ngfft=ngfftf,nkpt=mywfc%hdr%nkpt,npwarr= mywfc%hdr%npwarr , &
               &nsppol=dtset%nsppol,ntypat=cryst%ntypat,occ=mywfc%ebands%occ,&
               &pawang=pawang,pawrad=pawrad,pawtab=pawtab,prtvol=dtset%prtvol,psps=psps, &
               &rprimd=cryst%rprimd,ucvol=cryst%ucvol, xred=cryst%xred, exclude_bands=exclude_bands)

    if (.not. present(kg)) then
       ABI_FREE(ptr_kg)
    end if
    nullify(ptr_kg)
  end subroutine wfd_run_wannier
 !

end module m_wfd_wannier
