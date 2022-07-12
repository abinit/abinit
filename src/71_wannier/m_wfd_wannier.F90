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
  use m_build_info,      only :  abinit_version
  use defs_abitypes,     only : mpi_type
  use m_nctk,            only: nctk_try_fort_or_ncfile
  use m_io_tools,        only : open_file
  use m_fstrings,        only : ltoa, sjoin
  use m_fftcore,         only : print_ngfft
  use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
  use defs_wvltypes,  only : wvl_internal_type
  use m_dtset, only:dataset_type
  use m_hdr, only: hdr_type
  use m_wfd, only: wfd_t, wfd_init
  use m_crystal, only: crystal_t
  use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk, kpts_timrev_from_kptopt
  use m_ebands, only: ebands_from_hdr, ebands_print, ebands_expandk, ebands_free

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
         integer :: spin, ikcalc, ik_ibz, band, ispin, ik_bz
         real(dp) :: eigk
         logical:: bks_mask(mband, nkibz, nsppol)
         logical:: keep_ur(mband, nkibz, nsppol)
         bks_mask= .False.
         keep_ur=.False.
         do ispin=1, my_nspin
            spin = my_spin(ispin)
            do ik_bz=1,nkbz
               ik_ibz=indkk(ik_bz)
               do band=1, mband
                  bks_mask(band, ik_ibz ,spin) = .True.
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
         print *, "size(bz2ibz)", size(bz2ibz)
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

