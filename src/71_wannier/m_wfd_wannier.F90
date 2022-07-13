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
  use m_io_tools,        only : open_file
  use m_fstrings,        only : ltoa, sjoin
  use m_fftcore,         only : print_ngfft
  use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
  use defs_wvltypes,  only : wvl_internal_type
  use m_dtset, only:dataset_type
  use m_hdr, only: hdr_type, fform_from_ext
  use m_wfd, only: wfd_t, wfd_init
  use m_crystal, only: crystal_t
  use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk, kpts_timrev_from_kptopt
  use m_ebands, only: ebands_from_hdr, ebands_print, ebands_expandk, ebands_free, ebands_ncwrite


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
  use m_wannier_io,      only : write_Amn, compute_and_write_unk, write_eigenvalues, write_mmn
  use m_io_tools,        only : delete_file, get_unit, open_file

  use defs_wannier90
#ifdef HAVE_NETCDF
  use netcdf
#endif
  use m_nctk

  implicit none
  private
  integer,  parameter :: master=0
  public :: wfd_run_wannier
  public :: wfd_mlwfovlp

contains

  subroutine wfd_mlwfovlp(cryst, ebands, hdr, mpi_enreg, &
       & ngfftc, ngfftf,  dtset, dtfil,  &
       & pawang,  pawrad, pawtab, psps, comm, my_rank,nprocs )
    type(crystal_t), target,  intent(in) :: cryst
    type(ebands_t), target, intent(in) :: ebands
    type(hdr_type), target, intent(in) :: hdr
    integer,  target,intent(in) :: ngfftc(18),ngfftf(18)
    type(dataset_type),  target,intent(in) :: dtset
    type(datafiles_type), target,intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: mpi_enreg
    type(pseudopotential_type), target,intent(in) :: psps
    type(pawang_type), target,intent(in) :: pawang
    type(pawrad_type), target,intent(in) :: pawrad(psps%ntypat*psps%usepaw)
    type(pawtab_type), target,intent(in) :: pawtab(psps%ntypat*psps%usepaw)
    type(wfd_t) :: wfd
    integer, intent(in) :: comm, my_rank, nprocs
    integer :: natom, nsppol, nspinor,  mband
    integer :: nkibz, nkbz
    !!  nkibz,nkbz = Number of points in IBZ and BZ, respectively.
    real(dp),pointer :: wtk(:),kibz(:,:), kbz(:,:)
    !!  wtk(nkibz) = weights of the k-points in the IBZ (normalized to one).
    !!  kibz(3,nkibz) = k-points in the IBZ.
    !!  kbz(3,nkbz) = k-points in the BZ.
    integer, allocatable :: bz2ibz(:,:)
    !!  [bz2ibz(6,nkbz)]=Mapping BZ --> IBZ


    type(ebands_t), target :: ebands_bz

    type(hdr_type) :: hdr_bz

    integer, allocatable :: indkk(:)

    integer,allocatable :: bstart_ks(:,:)
    ! bstart_ks(nkcalc, nsppol)
    ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
    ! Depends on spin because all degenerate states should be included when symmetries are used.

    integer,allocatable :: bstop_ks(:,:)
    ! bstop_ks(nkcalc, nsppol)

    integer,allocatable :: nbcalc_ks(:,:)
    ! nbcalc_ks(nkcalc, nsppol)
    ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
    ! Depends on spin because all degenerate states should be included when symmetries are used.

    integer,allocatable :: kcalc2ibz(:,:)
    !kcalc2ibz(nkcalc, 6))
    ! Mapping ikcalc --> IBZ as reported by listkk.

    integer :: my_nspin
    ! Number of spins treated by this MPI rank

    integer,allocatable :: my_spin(:)
    ! my_spins(my_nspins)
    ! Indirect table giving the spin indices treated by this rank.
    ! Used only the collinear case with nspinor == 1

    integer :: my_nkpt
    ! Number of spins treated by this MPI rank

    integer,allocatable :: my_kpt(:)
    ! my_spins(my_nspins)
    ! Indirect table giving the spin indices treated by this rank.
    ! Used only the collinear case with nspinor == 1


    integer :: my_nkcalc
    ! Number of k-points treated by this MPI rank

    integer,allocatable :: my_ikcalc(:)
    ! my_ikcalc(my_nkcalc)
    ! List of ikcalc indices treated by this pool if k-point parallelism is activated.

    integer,allocatable :: myq2ibz_k(:)
    ! myq2ibz_k(my_nqibz_k)
    ! Mapping my q-point index --> index in nqibz_k arrays (IBZ_k)
    ! Differs from nqibz_k only if imag with tetra because in this case we can introduce a cutoff.

    integer(i1b),allocatable :: itreat_qibz(:)
    ! itreat_qibz(nqibz)
    ! Table used to distribute potentials over q-points in the IBZ.
    ! The loop over qpts in the IBZ(k) is MPI distributed inside qpt_comm accordinging to this table.
    ! 0 if this IBZ point is not reated by this proc.
    ! 1 if this IBZ is treated.

    call run_all()

  contains
    subroutine run_all()
      call initialize()
      call finalize()
    end subroutine run_all

    subroutine initialize()
      natom = cryst%natom
      nsppol = ebands%nsppol
      nspinor = ebands%nspinor
      ! What if the wavefunction is not computed in the IBZ?
      mband = dtset%mband
      call prepare_ebands_and_kpt_bz()
      call prepare_mpi_enreg()
      call split_mpi_tasks()
      call read_wfd()
      ! quantities to be transformed into full BZ
      !ebands or eig
      ! hdr
    end subroutine initialize

    subroutine finalize()
      call ebands_free(ebands_bz)
      call hdr_bz%free()
      call wfd%free()
      ABI_FREE(my_spin)
      nullify(wtk)
      nullify(kibz)
      nullify(kbz)
      ABI_FREE(bz2ibz)
      ABI_FREE(indkk)
      ABI_FREE(my_kpt)
    end subroutine finalize

    subroutine prepare_kpts()
      ! NOTE: this subroutine is not used.
      ! TODO:  Remove this once we know that the prepare_ebands_and_kpts_bz() works.
      !integer, allocatable :: temp
      !call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, &
      !     & ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
      !     nkibz, kibz, wtk, nkbz, &
      !     & kbz, bz2ibz=bz2ibz)
      ! NOTE free bz2ibz because it will be allocated again in ebands_expandk
      ABI_FREE(bz2ibz)

      !TODO Below is copied from sigmaph.
      ! But it seems the same thing is done through listkk in ebands_expandk.
      ! Which will also be called. So it is safe to 
      !
      ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
      !ABI_MALLOC(temp, (6, new%nqbz))

      !qrank = krank_from_kptrlatt(new%nqibz, new%qibz, qptrlatt, compute_invrank=.False.)
      !call qrank%get_mapping(new%nqbz, new%qbz, dksqmax, cryst%gmet, temp, &
      !     cryst%nsym, cryst%symafm, cryst%symrec, 1, use_symrec=.True.)
      !call qrank%free()

      !if (dksqmax > tol12) then
      !   ABI_ERROR("Cannot map BZ to IBZ!")
      !end if

      !new%ind_qbz2ibz(1,:) = temp(1,:)
      !new%ind_qbz2ibz(2,:) = temp(2,:)
      !new%ind_qbz2ibz(3,:) = temp(6,:)
      !new%ind_qbz2ibz(4,:) = temp(3,:)
      !new%ind_qbz2ibz(5,:) = temp(4,:)
      !new%ind_qbz2ibz(6,:) = temp(5,:)
      !ABI_FREE(temp)
      print *, "nkbz:", nkbz
      print *, "nkibz:", nkibz
    end subroutine prepare_kpts

    subroutine split_mpi_tasks()
      ! Only parallel over kpoints. Here we split over
      ! the fullBZ kpoints, as the wannierization task is for
      ! each kpt there.
      integer :: i
      my_nspin=nsppol
      ABI_MALLOC(my_spin, (nsppol))
      do i=1, nsppol
         my_spin(i) = i
      end do

      ! split kpoints
      call xmpi_split_block(nkbz, comm, my_nkpt, my_kpt)

      print *, "my_nspin", my_nspin
      print *, "my_spin", my_spin
      print *, "my_nkpt", my_nkpt
      print *, "my_kpt", my_kpt
    end subroutine split_mpi_tasks


    subroutine prepare_mpi_enreg()
      ! TODO check if we need to
    end subroutine prepare_mpi_enreg


    subroutine read_wfd()
      logical:: bks_mask(mband, nkibz, nsppol)
      logical:: keep_ur(mband, nkibz, nsppol)
      integer :: spin, ikcalc, ik_ibz, band, ispin
      real(dp) :: ecut_eff
      character(len=fnlen) :: wfk0_path
      wfk0_path = dtfil%fnamewffk
      ecut_eff = dtset%ecut * dtset%dilatmx **2

      ! setup bks_mask
      bks_mask= .False.
      keep_ur=.False.
      do ispin=1, my_nspin
         spin = my_spin(ispin)
         do ik_ibz=1,nkibz
            do band=1, mband
               !eigk = ebands%eig(band, ik_ibz, spin)
               bks_mask(band, ik_ibz ,spin) = .True.
            end do
         end do
      end do
      call wfd_init(wfd,cryst,pawtab,psps,keep_ur,ebands%mband,ebands%nband,ebands%nkpt,dtset%nsppol,bks_mask,&
           dtset%nspden,dtset%nspinor,ecut_eff,dtset%ecutsm,dtset%dilatmx,hdr%istwfk,ebands%kptns,ngfftc,&
           dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
      call wfd%read_wfk(wfk0_path,IO_MODE_MPI)
    end subroutine read_wfd


    subroutine prepare_bks_mask()
      ! Each wfd contains:
      ! the k and the k+b points.
      ! the k
      integer :: spin, ikcalc, ik_ibz, band, ispin, ik_bz, iknn_bz, inn, ik_nn
      real(dp) :: eigk
      logical:: bks_mask(mband, nkibz, nsppol)
      logical:: keep_ur(mband, nkibz, nsppol)
      integer :: nntot
      bks_mask= .False.
      keep_ur=.False.

      do ispin=1, my_nspin
         spin = my_spin(ispin)
         do ik_bz=1,nkbz
            ik_ibz=indkk(ik_bz)
            ! TODO : find the neighbors
            do band=1, mband
               bks_mask(band, ik_ibz ,spin) = .True.
            end do
            do inn=1, nntot
               ik_ibz=indkk(ik_nn)
            end do
         end do
      end do
    end subroutine prepare_bks_mask


    subroutine prepare_ebands_and_kpt_bz()
      real(dp) :: ecut_eff, dksqmax
      character(len=200) :: msg
      ! NOTE: is this OK to assume so?
      ecut_eff = dtset%ecut * dtset%dilatmx **2
      call ebands_expandk(ebands, cryst, ecut_eff, .True., dksqmax, &
           & bz2ibz, ebands_bz )
      if (dksqmax > tol12) then
         write(msg, '(3a,es16.6,4a)' )&
              'At least one of the k points could not be generated from a symmetrical one.',ch10,&
              'dksqmax=',dksqmax,ch10,&
              'Action: check your WFK file and k-point input variables',ch10,&
              '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
         ABI_ERROR(msg)
      end if

      nkibz = ebands%nkpt
      nkbz= ebands_bz%nkpt

      print *, "nkbz:", nkbz
      print *, "nkibz", nkibz

      kibz => ebands%kptns
      kbz => ebands_bz%kptns
      wtk => ebands_bz%wtk
      ABI_MALLOC(indkk, (nkbz))
      indkk(:) = bz2ibz(:, 1)

    end subroutine prepare_ebands_and_kpt_bz

    !! Wrapper to the wfd_sym_ug_kg
    !! OUTPUT
    !!  istwf_kbz: Time-reversal flag associated to output wavefunctions
    !!  npw_kbz: Number of G-vectors in kk_bz G-sphere
    !!  kg_kbz: G-vectors in reduced coordinates.
    !!  cgs_kbz: Periodic part of wavefunctions at kk_bz
    subroutine get_cg_ks(bstart, kk_ibz, kk_bz, spin, nband, mpw , cgs_kbz, kg_kbz, work, work_ngfft, istwf_kbz, npw_kbz)
      integer, intent(in) :: bstart, nband,  spin, mpw
      real(dp), intent(in) :: kk_ibz(3), kk_bz(3)
      integer, intent(in) :: work_ngfft(18)
      real(dp) :: work(2, work_ngfft(4), work_ngfft(5), work_ngfft(6))
      integer,intent(out) :: istwf_kbz, npw_kbz
      integer,intent(out) :: kg_kbz(3, mpw)
      real(dp),intent(out) :: cgs_kbz(2, mpw*nspinor, nband)
      real(dp) :: ecut
      !call wfd_sym_ur(Wfd,Cryst,Kmesh,band,ik_bz,spin,ur_kbz,trans,with_umklp,ur_kibz)
      call wfd%sym_ug_kg( ecut, kk_bz, kk_ibz, bstart, nband, spin, mpw, indkk, cryst, &
           & work_ngfft, work, istwf_kbz, npw_kbz, kg_kbz, cgs_kbz)
    end subroutine get_cg_ks



    subroutine prepare_hdr_bz()
      ! see wfk_tofullbz in m_wfk
      type(wvl_internal_type) :: dummy_wvl
      integer:: kptopt3=3

      call hdr_init_lowlvl(hdr_bz,ebands_bz,psps,pawtab,dummy_wvl,abinit_version,&
           hdr%pertcase,hdr%natom,hdr%nsym,hdr%nspden,hdr%ecut,dtset%pawecutdg,hdr%ecutsm,dtset%dilatmx,&
           hdr%intxc,hdr%ixc,hdr%stmbias,hdr%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,hdr%so_psp,&
           hdr%qptn,cryst%rprimd,cryst%xred,hdr%symrel,hdr%tnons,hdr%symafm,hdr%typat,hdr%amu,hdr%icoulomb,&
           kptopt3,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD,dtset%ivalence,dtset%cellcharge(1),&
           dtset%kptrlatt_orig,dtset%kptrlatt,&
           dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk)
      ! End CP modified
      if (psps%usepaw == 1) call pawrhoij_copy(hdr%pawrhoij, hdr_bz%pawrhoij)
    end subroutine prepare_hdr_bz

  end subroutine wfd_mlwfovlp


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
    mkmem= dtset%mkmem
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
          ! TODO: Parallel
          ! TODO: IBZ->BZ
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


  subroutine mlwfovlp_wfd(cryst, ebands, hdr, wfd, dtset, dtfil, mpi_enreg, &
       & pawang, pawrad, pawtab, psps, ngfftc, ngfftf, cg, cprj, kg, occ)
    type(crystal_t), target, intent(in) :: cryst
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
    ! input, not set.i
    real(dp),intent(in) :: cg(:, :)
    type(pawcprj_type), intent(in) :: cprj(:, :)
    integer, intent(in) :: kg(:,:)!,nattyp(ntypat),ngfft(18),npwarr(nkpt)
    real(dp),intent(in) :: occ(hdr%mband*hdr%nkpt*hdr%nsppol)
    integer :: ngfft(18)
    ! message
    character(len=1000) :: message
    ! Labels
    logical :: gamma_only,leig,lmmn,lwannierrun,spinors !,have_disentangled
    integer :: lwanniersetup
    ! filenames
    character(len=fnlen) :: wfnname
    character(len=fnlen) :: seed_name(dtset%nsppol)
    character(len=fnlen) :: fname,filew90_win(dtset%nsppol), &
         & filew90_wout(dtset%nsppol),filew90_amn(dtset%nsppol), &
         & filew90_ramn(dtset%nsppol)
    character(len=fnlen) :: filew90_mmn(dtset%nsppol),filew90_eig(dtset%nsppol)

    ! mpi
    integer :: spaceComm, nprocs, rank, master


    ! 
    integer :: mband,mcg,mcprj,mgfftc,mkmem,mpw,natom, nspinor, nsppol
    integer :: nfft,nkpt, max_num_bands, num_bands(dtset%nsppol)
    integer :: ntypat,prtvol
    integer, pointer :: atindx1(:)

    ! w90 required
    real(dp) :: real_lattice(3,3), recip_lattice(3,3)
    integer :: nwan(dtset%nsppol), num_nnmax, nband_inc(dtset%nsppol), nntot, mwan
    integer :: g1temp(3),ngkpt(3)
    integer,allocatable :: g1(:,:,:),gbound(:,:),icg(:,:)
    integer,allocatable:: iwav(:,:,:),kg_k(:,:),ovikp(:,:)
    integer,allocatable :: proj_l(:,:),proj_m(:,:),proj_radial(:,:)
    integer,allocatable :: proj_s_loc(:)
    real(dp),allocatable :: cm1(:,:,:,:,:,:),cm2_paw(:,:,:),cwavef(:,:)
    real(dp),allocatable :: denpot(:,:,:)
    real(dp),allocatable :: eigenvalues_w(:,:,:),fofgout(:,:),fofr(:,:,:,:)
    real(dp),allocatable :: proj_site(:,:,:),proj_x(:,:,:),proj_z(:,:,:),proj_zona(:,:)

    complex(dpc),allocatable :: M_matrix(:,:,:,:,:)
    complex(dpc),pointer:: A_matrix(:,:,:,:)
    complex(dpc),allocatable :: U_matrix_opt(:,:,:,:), U_matrix(:,:,:,:)
    logical,allocatable :: lwindow(:,:,:)
    real(dp),allocatable :: wann_centres(:,:,:),wann_spreads(:,:),xcart(:,:)
    real(dp),allocatable :: proj_s_qaxis_loc(:,:)
    complex(dpc),allocatable :: A_paw(:,:,:,:)
    logical,allocatable :: band_in(:,:)
    character(len=3),allocatable :: atom_symbols(:)
    logical,allocatable::just_augmentation(:,:)
#if defined HAVE_WANNIER90
    real(dp) :: spreadw(3,dtset%nsppol)
    real(dp),allocatable :: csix(:,:,:,:)
    real(dpc),allocatable :: occ_arr(:,:,:),occ_wan(:,:,:)
    real(dp),allocatable :: tdocc_wan(:,:)

#ifdef HAVE_NETCDF
    integer :: ncid, ncerr, nrpts
    character(len=fnlen) :: abiwan_fname
    integer :: have_disentangled_spin(dtset%nsppol)
    integer,allocatable :: irvec(:,:),ndegen(:)
#endif

#if defined HAVE_WANNIER90
    real(dp) :: corrvdw
    complex(dpc) :: caux,caux2,caux3
#endif


#endif

    call set_parameters()
    call sanity_check()
    call wannier90_unit_convert()
    call setup_wannier90()
    !3) Write Eigenvalues (file seed_name.eig)
    if(leig) then
       call write_eigenvalues(filew90_eig, ebands%eig, band_in,  eigenvalues_w, &
            &  nsppol, nkpt, mband,  dtset, rank, master )
    end if !leig

    if(lmmn) then
       call compute_and_write_mmn()
    end if !lmmn
    ABI_FREE(ovikp)
    ABI_FREE(g1)
    ABI_FREE(iwav)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !5) Calculate initial projections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call compute_and_write_amn()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !6) write files for wannier function plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( dtset%w90prtunk>0) then
       call compute_and_write_unk(wfnname, psps%usepaw, dtset%w90prtunk, &
            & mpi_enreg, ngfft, nsppol, dtset%nspinor,  &
            & nkpt, mband,  mpw, mgfftc, mkmem,  nprocs, rank, hdr%npwarr, &
            & band_in,  dtset, kg, cg)
    end if !dtset%w90prtunk

    if(lwannierrun) then
       call run_wannier_and_use_wannier()
    end if


  contains
    subroutine set_parameters()
      !Some initialization and checks
      !
      lwanniersetup=1 ! 1 is mandatory ( 0 is for debug)
      !to use lwanniersetup=0, one would need
      !to define which bands to exclude.
      lwannierrun=.true.   ! .false. and .true. are possible
      lmmn=.true.          ! .false. and .true. are possible
      leig=.true.          ! .false. and .true. are possible
      !
      gamma_only=.false. !not yet implemented
      spinors=.false. !not yet implemented
      !
      !mpi initialization
      !
      spaceComm=MPI_enreg%comm_cell
      nprocs=xmpi_comm_size(spaceComm)
      rank=MPI_enreg%me_kpt
      master=0

      ngfft(:)=ngfftc(:)

      mband = hdr%mband
      mpw=MAXVAL(hdr%npwarr)
      nspinor=hdr%nspinor
      nsppol=hdr%nsppol
      mband=hdr%mband
      mkmem= dtset%mkmem
      nsppol = hdr%nsppol
      mcg=mpw*nspinor*mband* mkmem *nsppol
      mcprj = nspinor*mband* mkmem*nsppol
      nfft=wfd%nfft
      atindx1=> cryst%atindx1

      !Generate seed names for wannier90 files, and file names
      call mlwfovlp_seedname(dtfil%fnameabo_w90,filew90_win,filew90_wout,filew90_amn,&
           & filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)

    end subroutine set_parameters




    subroutine sanity_check()
      if(rank==master) then
         write(message, '(a,a,a,a)' ) ch10,&
              &   '   mlwfovlp:  you should give k-point in the full brillouin zone ',ch10,&
              &   '   with explicit k-points (or kptopt=3) and istwfk 1'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
      end if
      !
      if(MPI_enreg%paral_spinor==1) then
         message = ' Parallelization over spinorial components not yet available !'
         ABI_ERROR(message)
      end if

      if(nsppol==2) then
         write(message, '(3a)' ) ch10,&
              &   '   mlwfovlp:  Calculating matrices for both spin polarization  ',ch10
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
      end if

      if(psps%npsp/=psps%ntypat) then
         ABI_ERROR("prb npsp")
      end if

    end subroutine sanity_check

    subroutine wannier90_unit_convert()
      real_lattice(:,:)=Bohr_Ang*cryst%rprimd(:,:)
      recip_lattice(:,:)=two_pi*cryst%gprimd(:,:)/Bohr_Ang
    end subroutine wannier90_unit_convert

    subroutine setup_wannier90()
      integer :: isppol
      num_nnmax=12 !limit fixed for compact structure in wannier_setup.
      ABI_MALLOC(g1,(3,nkpt,num_nnmax))
      ABI_MALLOC(ovikp,(nkpt,num_nnmax))
      ABI_MALLOC(atom_symbols,(natom))
      ABI_MALLOC(xcart,(3,natom))
      ABI_MALLOC(band_in,(mband,nsppol))
      ABI_MALLOC(proj_site,(3,mband,nsppol))
      ABI_MALLOC(proj_l,(mband,nsppol))
      ABI_MALLOC(proj_m,(mband,nsppol))
      ABI_MALLOC(proj_radial,(mband,nsppol))
      ABI_MALLOC(proj_x,(3,mband,nsppol))
      ABI_MALLOC(proj_s_loc,(mband))
      ABI_MALLOC(proj_s_qaxis_loc,(3,mband))
      ABI_MALLOC(proj_z,(3,mband,nsppol))
      ABI_MALLOC(proj_zona,(mband,nsppol))
      nullify(A_matrix)

      call mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&
           &  g1,lwanniersetup,mband,natom,nband_inc,nkpt,&
           &  nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&
           &  proj_l,proj_m,proj_radial,proj_site,proj_s_loc, proj_s_qaxis_loc, proj_x,proj_z,proj_zona,&
           &  real_lattice,recip_lattice,cryst%rprimd,seed_name,spinors,xcart,cryst%xred)

      do isppol=1, nsppol
         write(message, '(6a)' ) ch10,&
              &   '   mlwfovlp :  mlwfovlp_setup done -',ch10,&
              &   '-  see ',trim(filew90_wout(isppol)),' for details.'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
      end do
      !
      !some allocations after wannier90 setup
      !
      max_num_bands=maxval(num_bands(:))
      mwan=maxval(nwan(:))
      ABI_MALLOC(eigenvalues_w,(max_num_bands,nkpt,nsppol))
      ABI_MALLOC(M_matrix,(max_num_bands,max_num_bands,nntot,nkpt,nsppol))
      ABI_MALLOC(A_matrix,(max_num_bands,mwan,nkpt,nsppol))
      ABI_MALLOC(iwav,(mband,nkpt,nsppol))
    end subroutine setup_wannier90


    subroutine compute_and_write_mmn()
      integer :: ierr, isppol, ikpt

      !  In case of parallelization write out cg for all k-points
      ! TODO: Remove this and use wfd
      if (nprocs > 1) then
         call write_cg_and_cprj(dtset, cg, cprj, dtfil, iwav, hdr%npwarr, mband, natom, &
              &nsppol, nkpt,  MPI_enreg, rank, psps, pawtab)
      end if !MPI nprocs>1
      !
      !  End of MPI preliminarities
      !  Calculate PW contribution of overlaps
      !
      ABI_MALLOC(cm1,(2,mband,mband,nntot,nkpt,nsppol))
      ! this loops over spin internally
      call mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,&
           &   mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&
           &   hdr%npwarr,dtset%nspinor,nsppol,ovikp,dtfil%fnametmp_cg)
      write(message, '(a,a)' ) ch10,&
           &   '   mlwfovlp : PW part of overlap computed   '
      call wrtout(std_out,  message,'COLL')
      !
      !  compute PAW Contribution and add it to PW contribution
      !
      if(psps%usepaw==1) then
         call compute_mmn_for_paw_part()
      end if ! usepaw
      !
      call xmpi_barrier(spaceComm)
      call xmpi_sum(cm1,spaceComm,ierr)


      call write_Mmn(filew90_mmn, band_in, cm1, ovikp, g1, M_matrix, &
           &  nkpt, nsppol, nntot, mband, num_bands,  message,&
           & iam_master=(rank==master))

      ABI_FREE(cm1)

      if (nprocs > 1) then
         if(prtvol>0) then
            write(message, '(3a)' ) ch10,&
                 &       '   mlwfovlp :  Removing temporary files with cg and cprj (PAW)',ch10
            call wrtout(ab_out,message,'COLL')
            call wrtout(std_out,  message,'COLL')
         end if
         !
         !    Just master  node will remove the files
         !
         if(rank==master) then
            do isppol=1,nsppol
               do ikpt=1,nkpt
                  write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cg),ikpt,isppol
                  call delete_file(wfnname,ierr)
                  if(psps%usepaw==1) then
                     write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cprj),ikpt,isppol
                     call delete_file(wfnname,ierr)
                  end if
               end do !ikpt
            end do !isppol
         end if
      end if !MPI nprocs>1
    end subroutine compute_and_write_mmn

    subroutine compute_mmn_for_paw_part()
      integer :: isppol,  intot, ikpt1, ikpt2
      write(message, '(a,a)' ) ch10,&
           &     '** smatrix_pawinit : PAW part of overlap  '
      call wrtout(std_out,  message,'COLL')
      ABI_MALLOC(cm2_paw,(2,mband,mband))
      do isppol=1,nsppol
         do ikpt1=1,nkpt
            !
            !        MPI:cycle over k-points not treated by this node
            !
            if (nprocs>1 ) then !sometimes we can have just one processor
               if ( ABS(MPI_enreg%proc_distrb(ikpt1,1,isppol)-rank)  /=0) CYCLE
            end if

            write(message, '(a,i6,a,2i6)' ) &
                 &         '   processor',rank,' computes PAW part for kpt and spin',ikpt1,isppol
            call wrtout(std_out,  message,'COLL')

            do intot=1,nntot
               ikpt2= ovikp(ikpt1,intot)
               g1temp(:)=g1(:,ikpt1,intot)
               call smatrix_pawinit(atindx1,cm2_paw,cprj,ikpt1,ikpt2,isppol,&
                    &           g1temp,cryst%gprimd,dtset%kpt,mband,mband,mkmem,mpi_enreg,&
                    &           natom,dtset%nband,nkpt,dtset%nspinor,nsppol,dtset%ntypat,pawang,pawrad,pawtab,cryst%rprimd,&
                    &           dtfil%fnametmp_cprj,dtset%typat,cryst%xred)
               !          cm1(:,:,:,intot,ikpt1,isppol)=four_pi*cm2_paw(:,:,:)
               !           write(6,*) "ikpt1=",ikpt1
               !           do iband=1,mband
               !             write(6,*) "iband=",iband
               !             write(6,*) "Wannier PW       overlap",cm1(:,iband,iband,intot,ikpt1,isppol)
               !             write(6,*) "Wannier PAW      overlap",four_pi*cm2_paw(:,iband,iband)
               !             write(6,*) "Wannier PW+PAW   overlap",cm1(:,iband,iband,intot,ikpt1,isppol)+four_pi*cm2_paw(:,iband,iband)
               !           enddo
               cm1(:,:,:,intot,ikpt1,isppol)=cm1(:,:,:,intot,ikpt1,isppol)+four_pi*cm2_paw(:,:,:)
            end do ! intot
         end do ! ikpt1
      end do ! isppol
      ABI_FREE(cm2_paw)
      write(message, '(a,a)' ) ch10,&
           &     '   mlwfovlp : PAW part of overlap computed '
      call wrtout(std_out,  message,'COLL')
    end subroutine compute_mmn_for_paw_part

    subroutine compute_and_write_Amn()
      if(dtset%w90iniprj/=0 )  then

         call compute_Amn()
         if(rank==master) then
            if(dtset%w90iniprj==1) then
               call write_Amn(A_matrix, filew90_ramn, nsppol, mband, nkpt, num_bands, nwan, band_in)
            else
               call write_Amn(A_matrix, filew90_amn, nsppol, mband, nkpt, num_bands, nwan, band_in)
            end if
         end if

      end if
    end subroutine compute_and_write_Amn

    subroutine compute_Amn()
      integer :: lproj, isppol, iwan, ierr

      !
      !  Set value for lproj (type of projections to be computed)
      !  In PAW, options 5 and 6 are not in use.
      !  5 means that there will be a contribution from inside the spheres and another from the PW part
      !  6 means that we take into account just the inside-spheres contribution
      !  2 means that PW part will be calculated

      !
      lproj=dtset%w90iniprj
      if(dtset%w90iniprj == 5 ) lproj=2 ! Necessary to calculate PW contribution
      !
      ABI_MALLOC(just_augmentation,(mwan,nsppol))
      just_augmentation(:,:)=.false.

      if( psps%usepaw==1 .and. (dtset%w90iniprj==2 .or. dtset%w90iniprj>4)) then
         if (dtset%w90iniprj==6) just_augmentation(:,:)=.true.
         if (dtset%w90iniprj==5) then
            do isppol=1,nsppol
               do iwan=1,nwan(isppol)
                  !
                  !          Trick to skip the planewave contribution for some Wannier functions
                  !          (Not in production).
                  !
                  if(proj_radial(iwan,isppol) > 4) then
                     just_augmentation(iwan,isppol)=.true.
                     proj_radial(iwan,isppol)=proj_radial(iwan,isppol)-3
                     write(message, '(2a,2i4)' ) &
                          &             '   ','Skiping planewave contribution for iwan, ispin=',iwan,isppol
                     call wrtout(std_out,  message,'COLL')
                  end if !proj_radial>4
               end do !iwan
            end do !isppol
         end if !w90iniprj == 5
      end if !paw
      !
      !  Call mlwfovlp_proj (plane waves part of projections)
      !
      if (dtset%w90iniprj/=6) then ! option 6 not yet in use
         call mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,cryst%gprimd,just_augmentation,kg,&
              &     lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,&
              &     cryst%nattyp,nkpt,hdr%npwarr,&
              &     dtset%nspinor,nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,&
              &     proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,cryst%ucvol)
         write(message, '(a,a,a,a)' ) ch10,&
              &     '   mlwfovlp:  mlwfovlp_proj done -',ch10,&
              &     '   Projectors computed.'
         call wrtout(std_out,  message,'COLL')
      end if !w90proj/=6
      !
      !  Calculate inside-sphere part of projections (PAW)
      !
      if (psps%usepaw ==1 .and. ( dtset%w90iniprj>4)) then
         ABI_MALLOC(A_paw,(max_num_bands,mwan,nkpt,nsppol))
         call mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&
              &     mwan,natom,dtset%nband,nkpt,&
              &     dtset%nspinor,nsppol,dtset%ntypat,nwan,pawrad,pawtab,&
              &     proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&
              &     cryst%rprimd,dtset%typat,cryst%xred)
         !
         write(message, '(a,a,a,a)' ) ch10,&
              &     '   mlwfovlp:  mlwfovlp_proj_paw done -',ch10,&
              &     '   Inside-spheres part of projectors computed.'
         call wrtout(std_out,  message,'COLL')
         !
         !    Add in-sphere contribution to A_matrix
         !
         !    w90iniprj==5. Plane waves + augmentation contributions
         !
         if(dtset%w90iniprj==5) A_matrix(:,:,:,:)=A_matrix(:,:,:,:)+A_paw(:,:,:,:)
         !
         !    w90iniprj==6. Just augmentation contribution
         !
         if(dtset%w90iniprj==6) A_matrix(:,:,:,:)=A_paw(:,:,:,:)
         !
         !    deallocations
         !
         ABI_FREE(A_paw)
      end if !usepaw==1

      ABI_FREE(just_augmentation)
      !
      call xmpi_barrier(spaceComm)
      call xmpi_sum(A_matrix,spaceComm,ierr)

      ABI_FREE(proj_site)
      ABI_FREE(proj_l)
      ABI_FREE(proj_m)
      ABI_FREE(proj_radial)
      ABI_FREE(proj_x)
      ABI_FREE(proj_z)
      ABI_FREE(proj_zona)
      ABI_FREE(proj_s_loc)
      ABI_FREE(proj_s_qaxis_loc)
    end subroutine compute_Amn

    subroutine run_wannier_and_use_wannier()

      if(lwanniersetup.ne.1) ABI_ERROR("lwanniersetup.ne.1")
      ABI_MALLOC(U_matrix,(mwan,mwan,nkpt,nsppol))
      ABI_MALLOC(U_matrix_opt,(max_num_bands,mwan,nkpt,nsppol))
      ABI_MALLOC(lwindow,(max_num_bands,nkpt,nsppol))
      ABI_MALLOC(wann_centres,(3,mwan,nsppol))
      ABI_MALLOC(wann_spreads,(mwan,nsppol))
!  Initialize
   U_matrix(:,:,:,:)=czero
   U_matrix_opt(:,:,:,:)=czero
   lwindow(:,:,:)=.false.
   wann_centres(:,:,:)=zero
   wann_spreads(:,:)=zero
!
!  write(std_out,*) seed_name
!  write(std_out,*) ngkpt
   ngkpt(1)=dtset%kptrlatt(1,1)
   ngkpt(2)=dtset%kptrlatt(2,2) ! ajouter test de verif que kptrlatt est bien diagonal
   ngkpt(3)=dtset%kptrlatt(3,3)
!  write(std_out,*) nkpt
!  write(std_out,*) rprimd*Bohr_Ang
!  write(std_out,*) two_pi*gprimd/Bohr_Ang
!  write(std_out,*) mband
!  write(std_out,*) "nwan",nwan
!  write(std_out,*) nntot
!  write(std_out,*) natom
!  write(std_out,*) atom_symbols
!  write(std_out,*) xcart
!  write(std_out,*) num_bands,num_bands,nntot,nkpt
!  write(std_out,*) wann_spreads
!  wann_spreads=2
!  do i=1, nkpt
!  do j=1, nntot
!  write(std_out,*) i,j
!  do k=1, num_bands
!  do l=1, num_bands
!  write(std_out,*) "m",M_matrix(l,k,j,i,1)
!  enddo
!  enddo
!  enddo
!  enddo
   !  FIXME: looks like there is no automatic test which goes through here: g95 bot did not catch
   !  the missing deallocations

   call run_wannier90()
   call write_to_abiwann_nc()
   call evaluate_vdw_energy_with_mlwf()

   ABI_FREE(wann_centres)
   ABI_FREE(wann_spreads)
   ABI_FREE(U_matrix)
   ABI_FREE(U_matrix_opt)
   ABI_FREE(lwindow)
end subroutine run_wannier_and_use_wannier


   subroutine run_wannier90()
#if defined HAVE_WANNIER90
     integer :: isppol, ierr
     do isppol=1,nsppol
        !    when nsppol>1, master runs isppol 1 and rank==1 runs isppol 2
        if(nprocs>1 .and. isppol==1.and.rank.ne.master) cycle
        if(nprocs>1 .and. isppol==2.and.rank.ne.1) cycle

        write(message, '(8a)' ) ch10,&
             &     '** mlwfovlp :   call wannier90 library subroutine wannier_run ',ch10,&
             &     '   Calculation is running         ',ch10,&
             &     '-  see ',trim(filew90_wout(isppol)),' for details.'
        call wrtout(std_out,  message,'COLL')
        !
        call wannier_run(trim(seed_name(isppol)),ngkpt,nkpt,&            !input
             &    real_lattice,recip_lattice,dtset%kpt,num_bands(isppol),& !input
             &    nwan(isppol),nntot,natom,atom_symbols,&                  !input
             &    xcart*Bohr_Ang,gamma_only,M_matrix(:,:,:,:,isppol),A_matrix(:,:,:,isppol),eigenvalues_w(:,:,isppol),& !input
             &    U_matrix(1:nwan(isppol),1:nwan(isppol),:,isppol),& !output
             &    U_matrix_opt(1:num_bands(isppol),1:nwan(isppol),:,isppol),& !output
             &    lwindow_loc=lwindow(1:num_bands(isppol),:,isppol),& !output
             &    wann_centres_loc=wann_centres(:,1:nwan(isppol),isppol),&     !output
             &    wann_spreads_loc=wann_spreads(1:nwan(isppol),isppol),spread_loc=spreadw(:,isppol))                            !output


        write(message, '(7a)' ) ch10,&
             &     '   mlwfovlp :  mlwfovlp_run completed -',ch10,&
             &     '-  see ',trim(filew90_wout(isppol)),' for details.',ch10
        call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,  message,'COLL')

     end do !isppol

     !  collect output of  wannier90 from different processors
     call xmpi_barrier(spaceComm)

     call xmpi_sum(U_matrix,spaceComm,ierr)
     call xmpi_sum(U_matrix_opt,spaceComm,ierr)
     call xmpi_lor(lwindow,spaceComm)
     call xmpi_sum(wann_centres,spaceComm,ierr)
     call xmpi_sum(wann_spreads,spaceComm,ierr)
#endif
   end subroutine run_wannier90


   subroutine write_to_abiwann_nc()
     integer :: isppol, ierr, ncid
     ! Output ABIWAN.nc file
#ifdef HAVE_WANNIER90
#ifdef HAVE_NETCDF
     if (dtset%kptopt == 0) then
        ABI_WARNING("Output of ABIWAN.nc requires kptopt /= 0. ABIWAN.nc file won't be produced!")
        ! Need kptrlatt in wigner_seitz and client code need to know the k-grid.
     end if
     if (rank == master .and. dtset%kptopt /= 0) then
        abiwan_fname = strcat(dtfil%filnam_ds(4), "_ABIWAN.nc")
        call wrtout(std_out, sjoin("Saving wannier90 ouput results in:", abiwan_fname))
        call wigner_seitz([zero, zero, zero], [2, 2, 2], dtset%kptrlatt, cryst%rmet, nrpts, irvec, ndegen)
        ! We know if disentanglement has been done by looking at the output values of lwindow
        ! Not elegant but it is the only way to avoid the parsing of the wannier input.
        ! In wannier_run lwindow is set to True if not disentanglement
        have_disentangled_spin = 0
        do isppol=1,nsppol
           !if nwan(isppol) < num_bands(isppol)
           if (.not. all(lwindow(:,:,isppol))) have_disentangled_spin(isppol) = 1
        end do

        NCF_CHECK(nctk_open_create(ncid, abiwan_fname, xmpi_comm_self))
        NCF_CHECK(hdr%ncwrite(ncid, fform_from_ext("ABIWAN"), nc_define=.True.))
        NCF_CHECK(cryst%ncwrite(ncid))
        NCF_CHECK(ebands_ncwrite(ebands, ncid))

        ncerr = nctk_def_dims(ncid, [ &
             nctkdim_t("mwan", mwan), &
             nctkdim_t("max_num_bands", max_num_bands), &
             nctkdim_t("nrpts", nrpts) &
             ], defmode=.True.)
        NCF_CHECK(ncerr)

        ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "nntot"])
        NCF_CHECK(ncerr)
        !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
        !NCF_CHECK(ncerr)

        ncerr = nctk_def_arrays(ncid, [ &
             nctkarr_t("nwan", "int", "number_of_spins"), &
             nctkarr_t("num_bands", "int", "number_of_spins"), &
             nctkarr_t("band_in_int", "int", "max_number_of_states, number_of_spins"), &
             nctkarr_t("lwindow_int", "int", "max_num_bands, number_of_kpoints, number_of_spins"), &
                                !nctkarr_t("exclude_bands", "int", "max_number_of_states, number_of_spins"), &
                                !nctkarr_t("eigenvalues_w", "int", "max_num_bands, number_of_kpoints, number_of_spins"), &
             nctkarr_t("spread", "dp", "three, number_of_spins"), &
                                !nctkarr_t("A_matrix", "dp", "two, max_num_bands, mwan, number_of_kpoints, number_of_spins"), &
             nctkarr_t("irvec", "int", "three, nrpts"), &
             nctkarr_t("ndegen", "int", "nrpts"), &
             nctkarr_t("have_disentangled_spin", "int", "number_of_spins"), &
             nctkarr_t("U_matrix", "dp", "two, mwan, mwan, number_of_kpoints, number_of_spins"), &
             nctkarr_t("U_matrix_opt", "dp", "two, max_num_bands, mwan, number_of_kpoints, number_of_spins"), &
             nctkarr_t("wann_centres", "dp", "three, mwan, number_of_spins"), &
             nctkarr_t("wann_spreads", "dp", "mwan, number_of_spins") &
             ])
        NCF_CHECK(ncerr)

        ! Write data.
        NCF_CHECK(nctk_set_datamode(ncid))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nntot"), nntot))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nwan"), nwan))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "num_bands"), num_bands))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "band_in_int"), l2int(band_in)))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "lwindow_int"), l2int(lwindow)))
        !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "exclude_bands"), exclude_bands))
        !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eigenvalues_w"), eigenvalues_w))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "spread"), spreadw))
        !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "A_matrix"), c2r(A_matrix)))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "irvec"), irvec))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ndegen"), ndegen))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "have_disentangled_spin"), have_disentangled_spin))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "U_matrix"), c2r(U_matrix)))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "U_matrix_opt"), c2r(U_matrix_opt)))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wann_centres"), wann_centres))
        NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wann_spreads"), wann_spreads))
        NCF_CHECK(nf90_close(ncid))

        ABI_FREE(irvec)
        ABI_FREE(ndegen)
     end if
#endif
#endif
   end subroutine write_to_abiwann_nc


   subroutine evaluate_vdw_energy_with_mlwf()

#ifdef HAVE_WANNIER90
     integer :: isppol, ikpt, iband, ii, jj,kk, iwan
!  CALL SILVESTRELLI'S APPROACH TO EVALUATE vdW INTERACTION ENERGY USING MLWF!!
!  ----------------------------------------------------------------------------------------------
   if (dtset%vdw_xc==10.or.dtset%vdw_xc==11.or.dtset%vdw_xc==12.or.dtset%vdw_xc==14.and.rank==master) then
!    vdw_xc==10,11,12,14 starts the
!    vdW interaction using MLWFs
     write(std_out,*) 'nwan(nsppol)=',ch10
     do ii=1,nsppol
       write(std_out,*) 'nsppol=',ii, 'nwan(nsppol)=',nwan(ii),ch10
     end do
     write(std_out,*) 'mwan=', mwan, ch10

     ABI_MALLOC(occ_arr,(mband,nkpt,isppol))
     ABI_MALLOC(occ_wan,(mwan,nkpt,nsppol))
     ABI_MALLOC(tdocc_wan,(mwan,nsppol))

     occ_arr(:,:,:)=zero
     occ_wan(:,:,:)=zero
     tdocc_wan(:,:)=zero
     jj = 0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iband=1,num_bands(isppol)
           jj = jj + 1
           occ_arr(iband,ikpt,isppol) = occ(jj)
         end do
       end do
     end do

     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iwan=1,nwan(isppol)
           caux=czero
           caux2=czero
           caux3=czero
           do iband=1,num_bands(isppol) !nband_inc(isppol) !nwan(isppol)
             do ii=1,nwan(isppol)
               caux=U_matrix(ii,iwan,ikpt,isppol)*U_matrix_opt(iband,ii,ikpt,isppol)
!              DEBUG
!              if(ISNAN(dble(caux))) then
!              write(std_out,*) 'NaN: caux(ikpt,iwan,iband,ii):',ikpt,iwan,iband,ii,ch10
!              end if
!              END DEBUG
               do kk=1,nwan(isppol)
                 caux2=conjg(U_matrix(kk,iwan,ikpt,isppol))*conjg(U_matrix_opt(iband,kk,ikpt,isppol))
                 caux3= caux3+caux*caux2*occ_arr(iband,ikpt,isppol) !take care here as exclude_bands case is not well
!                DEBUG
!                if(ISNAN(dble(caux2))) then
!                write(std_out,*) 'NaN: caux2(ikpt,iwan,iband,kk):',ikpt,iwan,iband,kk,ch10
!                end if
!                if(ISNAN(dble(caux3))) then
!                write(std_out,*) 'NaN: caux3(ikpt,iwan,iband,kk,jj):',ikpt,iwan,iband,kk,jj
!                end if
!                END DEBUG
               end do
             end do
           end do
           occ_wan(iwan,ikpt,isppol) = dble(caux3)
!          DEBUG
!          write(std_out,*) occ_wan(iwan,ikpt,isppol)
!          END DEBUG
!          end do
         end do
       end do
     end do

     write(std_out,*) ch10,'MLWFs Occupation Matrix diagonal terms:',ch10

     do jj=1,nsppol
        forall(iwan=1:nwan(jj)) tdocc_wan(iwan,jj) = &
             & sum(occ_wan(iwan,1:nkpt,jj)) / real(nkpt,dp)
       write(std_out,*) 'tdocc_wan(iwan),isppol:',ch10
       write(std_out,*) (tdocc_wan(iwan,jj),iwan=1,nwan(jj)),jj
     end do

     ABI_MALLOC(csix,(mwan,mwan,nsppol,nsppol))

     call evdw_wannier(csix,corrvdw,mwan,natom,nsppol, &
          & nwan,tdocc_wan,dtset%vdw_nfrag,&
          & dtset%vdw_supercell,dtset%vdw_typfrag,dtset%vdw_xc, &
          & cryst%rprimd,wann_centres,wann_spreads,cryst%xcart)

     ABI_FREE(csix)
     ABI_FREE(occ_arr)
     ABI_FREE(occ_wan)
     ABI_FREE(tdocc_wan)
   end if
#endif

   end subroutine evaluate_vdw_energy_with_mlwf

  end subroutine mlwfovlp_wfd
  !

end module m_wfd_wannier

