!!****m*ABINIT/m_wfd_wannie
!! NAME
!!  m_wfd_wannier
!!
!! FUNCTION
!!  The high level wfd_t inteface for building wannier functions
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

  use m_fstrings,      only : strcat, sjoin
  use m_numeric_tools, only : uniformrandom, simpson_int, c2r, l2int
  use m_build_info,      only :  abinit_version
  use defs_abitypes,     only : mpi_type
  use m_nctk,            only: nctk_try_fort_or_ncfile
  use m_fstrings,        only : ltoa, sjoin
  use m_fftcore,         only : print_ngfft
  use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
  use defs_wvltypes,  only : wvl_internal_type
  use m_dtset, only:dataset_type
  use m_hdr, only: hdr_type, fform_from_ext, hdr_init_lowlvl
  use m_wfd, only: wfd_t, wfd_init
  use m_crystal, only: crystal_t
  use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk, kpts_timrev_from_kptopt
  use m_ebands, only: ebands_from_hdr, ebands_print, ebands_expandk, ebands_free, ebands_ncwrite
  use m_fftcore,        only : get_kg
  use m_geometry,  only : wigner_seitz

  use m_pawang,          only : pawang_type
  use m_pawrad,          only : pawrad_type
  use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
  use m_pawcprj,         only : pawcprj_type
  use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, pawrhoij_inquire_dim
  use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
  use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init, pawfgrtab_print
  use defs_datatypes,    only : pseudopotential_type, ebands_t
  use m_dtfil,          only : datafiles_type
  use m_paw_overlap, only : smatrix_pawinit
  use m_evdw_wannier, only : evdw_wannier

  use m_abstract_wf,     only : abstract_wf, compute_iwav, write_cg_and_cprj, wann_ksetting_t, init_mywfc
  !use m_mlwfovlp,        only : mlwfovlp, mlwfovlp_pw, mlwfovlp_proj, mlwfovlp_projpaw, mlwfovlp_setup, mlwfovlp_seedname
  use m_mlwfovlp2,        only : mlwfovlp2

  use defs_wannier90
#ifdef HAVE_NETCDF
  use netcdf
#endif
  use m_nctk

  implicit none
  private
  integer,  parameter :: master=0
  public :: wfd_run_wannier
  !public :: wfd_mlwfovlp

contains


  subroutine wfd_run_wannier(cryst, ebands, hdr, mpi_enreg, &
       & ngfftc, ngfftf,  wfd, dtset, dtfil,  &
       & pawang,  pawrad, pawtab, psps , kg, cg, cprj)
    type(crystal_t), intent(in) :: cryst
    type(ebands_t), intent(inout) :: ebands
    type(hdr_type), intent(inout) :: hdr
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
    logical, allocatable :: bks_mask(:, :, :), keep_ur(:, :, :)
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

    integer :: spaceComm, nprocs, rank, master
    integer :: ik_ibz, spin, npw_k, ikg, iblk
    integer :: my_nbands,my_band_list(hdr%mband)

    print *, "============================================================"
    print *, "Starting WFD Wannier"
    print *, "============================================================"

    spaceComm=MPI_enreg%comm_cell
    nprocs=xmpi_comm_size(spaceComm)
    !rank=xmpi_comm_self
    ! TODO: is MPI_enreg  initialized?
    rank=MPI_enreg%me_kpt
    master=0

    mgfftc=dtset%mgfft
    mpw=MAXVAL(hdr%npwarr)
    nspinor=hdr%nspinor
    mband=hdr%mband
    nsppol = hdr%nsppol
    mcg=mpw*nspinor*mband* mkmem *nsppol
    mcprj = nspinor*mband* mkmem*nsppol
    ! Is this correct for paw?
    nfft=product(ngfftc(1:3))

    print *, "mpw:", mpw
    print *, "nspinor:", nspinor
    print *, "mband:", mband
    print *, "mkmem:", mkmem
    print *, "nsppol", nsppol
    print *, "mcg", mcg
    print *, "mcprj", mcprj
    print *, "nfft", nfft


    !print *, "kg=", kg

       if (present(cg)) then

         mkmem= dtset%mkmem
         print *, "mkmem:", mkmem
         call init_mywfc(mywfc=mywfc, ebands=ebands, cg=cg, cprj=cprj, &
           &  cryst=cryst,dtset=dtset, &
           & dtfil=dtfil, hdr=hdr, MPI_enreg=MPI_enreg, nprocs=nprocs, psps=psps, &
           & pawtab=pawtab, rank=rank, comm=spaceComm)
         ptr_kg => kg
          ! call mlwfovlp2(mywfc=mywfc, crystal=cryst, ebands=ebands, hdr=hdr, atindx1=cryst%atindx1  &
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
           mkmem = mywfc%hdr%nkpt
           mpw=maxval(mywfc%hdr%npwarr)
           print *, "mkmem:", mkmem
           ! TODO : why is mpw in hdr_bz larger than in hdr_ibz?
           print *, "mpw:", mpw
           mcg=mpw*nspinor*mband* mkmem *nsppol
           mcprj = nspinor*mband* mkmem*nsppol
           block
             integer :: ik
             do ik =1, mkmem
               print *, ik, ":", mywfc%ebands%kptns(:, ik)
             end do
           end block

           if (present(kg)) then
             ptr_kg=> kg
           else
             ABI_MALLOC(ptr_kg, (3, mpw*mkmem))
             call mywfc%get_kgs(ptr_kg)
           end if
         end if

         call mlwfovlp2(mywfc=mywfc, crystal=cryst, ebands=mywfc%ebands, hdr=mywfc%hdr, atindx1=cryst%atindx1  &
               &,dtset=mywfc%dtset,dtfil=dtfil, &
               & eigen=mywfc%ebands%eig,gprimd=cryst%gprimd,kg=ptr_kg,&
               & mband=mywfc%hdr%mband,mcg=mcg,mcprj=mcprj,mgfftc=mgfftc, &
               & mkmem=mkmem,mpi_enreg=mpi_enreg,mpw=mpw,natom=cryst%natom,&
               & nattyp=cryst%nattyp,nfft=nfft,ngfft=ngfftf,nkpt=mywfc%hdr%nkpt,npwarr= mywfc%hdr%npwarr , &
               &nsppol=dtset%nsppol,ntypat=cryst%ntypat,occ=mywfc%ebands%occ,&
               &pawang=pawang,pawrad=pawrad,pawtab=pawtab,prtvol=dtset%prtvol,psps=psps, &
               &rprimd=cryst%rprimd,ucvol=cryst%ucvol, xred=cryst%xred)

    if (.not. present(kg)) then
       ABI_FREE(ptr_kg)
    end if
    nullify(ptr_kg)
  end subroutine wfd_run_wannier
 !

end module m_wfd_wannier
