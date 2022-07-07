!!****m*ABINIT/m_wfd_wannier
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
  use defs_abitypes,     only : mpi_type
  use m_nctk,            only: nctk_try_fort_or_ncfile
  use m_io_tools,        only : open_file
  use m_fstrings,        only : ltoa, sjoin
  use m_fftcore,         only : print_ngfft
  use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq

  use m_dtset, only:dataset_type
  use m_hdr, only: hdr_type
  use m_wfd, only: wfd_t
  use m_crystal, only: crystal_t
  use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk, kpts_timrev_from_kptopt
  use m_ebands, only: ebands_from_hdr, ebands_print

  use m_pawang,          only : pawang_type
  use m_pawrad,          only : pawrad_type
  use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
  use m_pawcprj,         only : pawcprj_type
  use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, pawrhoij_inquire_dim
  use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
  use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init, pawfgrtab_print
  use defs_datatypes,    only : pseudopotential_type, ebands_t
  use m_dtfil,          only : datafiles_type

  use m_mlwfovlp,        only : mlwfovlp
  use m_wannier_io,      only : write_Amn, compute_and_write_unk, write_eigenvalues
  implicit none
  private
  integer,  parameter :: master=0
  public :: wfd_run_wannier
contains

  ! fix mpi_enreg for nproc=1
  subroutine mpi_enreg_init_seq(mpi_enreg)
    type(mpi_type), intent(inout) :: mpi_enreg
  end subroutine mpi_enreg_init_seq

  subroutine wfd_run_wannier(cryst, ebands, hdr, mpi_enreg, &
       & ngfftc, ngfftf,  wfd, dtset, dtfil,  &
       & pawang,  pawrad, pawtab, psps )
    type(crystal_t), intent(in) :: cryst
    type(ebands_t), intent(in) :: ebands
    type(hdr_type), intent(in) :: hdr
    integer, intent(in) :: ngfftc(18),ngfftf(18)
    type(wfd_t), intent(inout) :: wfd
    type(dataset_type), intent(in) :: dtset
    type(datafiles_type),intent(in) :: dtfil
    type(mpi_type), intent(in) :: mpi_enreg
    type(pseudopotential_type),intent(in) :: psps
    type(pawang_type),intent(in) :: pawang
    type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
    type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    integer :: nfft
    integer :: mcg, mcprj, mgfftc
    real(dp), allocatable :: cg(:, :)
    integer, allocatable :: kg(:, :)
    type(pawcprj_type), allocatable :: cprj(:, :)

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

    integer :: ik_ibz, spin, npw_k, ikg, iblk
    integer :: my_nbands,my_band_list(wfd%mband)

    print *, "============================================================"
    print *, "Starting WFD Wannier"
    print *, "============================================================"

    mpw=MAXVAL(wfd%npwarr)
    nspinor=hdr%nspinor
    mband=hdr%mband
    mkmem= hdr%nkpt ! FIXME: this is not right in parallel mode
    nsppol = hdr%nsppol
    mcg=mpw*nspinor*mband* mkmem *nsppol
    mcprj = nspinor*mband* mkmem*nsppol
    nfft=wfd%nfft

    print *, "mpw:", mpw
    print *, "nspinor:", nspinor
    print *, "mband:", mband
    print *, "mkmem:", mkmem
    print *, "nsppol", nsppol
    print *, "mcg", mcg
    print *, "mcprj", mcprj
    print *, "nfft", nfft

    !gather kg.
    !do ik_ibz = 1, mkmem
    !   print *, "npw ik_ibz=", ik_ibz, "npw:", wfd%npwarr(ik_ibz)
    !end do

       ABI_MALLOC(cg, (2, mcg))
       ABI_MALLOC(cprj, (cryst%natom, mcprj))
       ABI_MALLOC(kg, (3, mpw*mkmem))

       ikg=0
       do ik_ibz = 1, mkmem
          npw_k = wfd%npwarr(ik_ibz)
          !kg(:, (ik_ibz-1)*mpw+1:ik_ibz*mpw) = wfd%Kdata(ik_ibz)%kg_k(:,:)
          kg(:,1+ikg:npw_k+ikg)=wfd%Kdata(ik_ibz)%kg_k(:,:)
          ikg=ikg+ npw_k
       end do

   iblk=0
   do spin =1, nsppol
       !FIXME fill cg and cprj
      do ik_ibz=1, mkmem
         npw_k=wfd%npwarr(ik_ibz)
         call wfd%mybands(ik_ibz, spin, my_nbands, my_band_list(ik_ibz))
         !mcg=mpw*nspinor*mband* mkmem *nsppol
         print *, size(cg, dim=1)
         print *, size(my_band_list, dim=1)
         print *, wfd%mband
         print *, mpw*nspinor*wfd%mband
         call wfd%extract_cgblock(band_list=my_band_list, ik_ibz=ik_ibz, &
            & spin=spin, cgblock=cg(:,iblk+1: iblk+npw_k*nspinor*wfd%mband))
         iblk = iblk + mpw*nspinor*my_nbands
    end do
 end do

    call mlwfovlp(crystal=cryst, ebands=ebands, hdr=hdr, atindx1=cryst%atindx1  &
         &,cg=cg,cprj=cprj,dtset=dtset,dtfil=dtfil, &
         &eigen=ebands%eig,gprimd=cryst%gprimd,kg=kg,&
         & mband=wfd%mband,mcg=mcg,mcprj=mcprj,mgfftc=mgfftc, &
         &mkmem=mkmem,mpi_enreg=mpi_enreg,mpw=mpw,natom=cryst%natom,&
         & nattyp=cryst%nattyp,nfft=nfft,ngfft=ngfftf,nkpt=hdr%nkpt,npwarr= hdr%npwarr , &
         &nsppol=dtset%nsppol,ntypat=cryst%ntypat,occ=ebands%occ,&
         &pawang=pawang,pawrad=pawrad,pawtab=pawtab,prtvol=dtset%prtvol,psps=psps, &
         &rprimd=cryst%rprimd,ucvol=cryst%ucvol,xred=cryst%xred)
  end subroutine wfd_run_wannier


  subroutine mlwfovlp_wfd(cryst, ebands, hdr, wfd, dtset, dtfil, mpi_enreg, pawang, pawrad, pawtab, psps, ngfftc, ngfftf)
    type(crystal_t), intent(in) :: cryst
    type(ebands_t), intent(in) :: ebands
    type(hdr_type), intent(in) :: hdr
    integer, intent(in) :: ngfftc(18),ngfftf(18)
    type(wfd_t), intent(inout) :: wfd
    type(dataset_type), intent(in) :: dtset
    type(datafiles_type),intent(in) :: dtfil
    type(mpi_type), intent(in) :: mpi_enreg
    type(pseudopotential_type),intent(in) :: psps
    type(pawang_type),intent(in) :: pawang
    type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
    type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  end subroutine mlwfovlp_wfd



  end module m_wfd_wannier

