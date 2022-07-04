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

  implicit none
  private
  integer,  parameter :: master=0

  public :: wfd_run_wannier

contains

  subroutine wfd_run_wannier(cryst, ebands, hdr, mpi_enreg, &
       & nfft, ngfftc, ngfftf,  wfd, dtset, dtfil,  &
       & pawang,  pawrad, pawtab, psps )
    type(crystal_t) :: cryst
    type(ebands_t) :: ebands
    type(hdr_type) :: hdr

    integer :: nfft
    integer :: ngfftc(18),ngfftf(18)
    type(wfd_t), intent(in) :: wfd

    type(dataset_type), intent(in) :: dtset
    type(datafiles_type),intent(in) :: dtfil
    type(mpi_type), intent(in) :: mpi_enreg
    type(pseudopotential_type),intent(inout) :: psps
    type(pawang_type),intent(inout) :: pawang
    type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
    type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

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
    ! TODO: mkmem
    ! TODO: gather kg from wfd
    ! TODO: mcg
    ! TODO: mcprj

    integer :: ikpt

    mcg=hdr%npwarr(ikpt)*hdr%nspinor*hdr%mband* mpi_enreg%mkmem(ikpt)*hdr%nsppol
    mcprj = hdr%nspinor*hdr%mband* mpi_enreg%mkmem(ikpt)*hdr%nsppol
    mkmem = mpi_enreg%mkmem(ikpt)

    ABI_MALLOC(cg, (2, mcg))
    ABI_MALLOC(cprj, (cryst%natom, mcprj))
    ABI_MALLOC(kg, (hdr%npwarr(ikpt), mpi_enreg%mkmem(ikpt)))

    !FIXME fill cg and cprj

    ! gather kg
    !kg(3,mpw*mkmem)


    call mlwfovlp(crystal=cryst, ebands=ebands, hdr=hdr, atindx1=cryst%atindx1  &
         &,cg=cg,cprj=cprj,dtset=dtset,dtfil=dtfil, &
         &eigen=ebands%eig,gprimd=cryst%gprimd,kg=kg,&
         & mband=wfd%mband,mcg=mcg,mcprj=mcprj,mgfftc=mgfftc, &
         &mkmem=mkmem,mpi_enreg=mpi_enreg,mpw=hdr%npwarr(ikpt),natom=cryst%natom,&
         & nattyp=cryst%nattyp,nfft=nfft,ngfft=ngfftf,nkpt=hdr%nkpt,npwarr= hdr%npwarr , &
         &nsppol=dtset%nsppol,ntypat=cryst%ntypat,occ=ebands%occ,&
         &pawang=pawang,pawrad=pawrad,pawtab=pawtab,prtvol=dtset%prtvol,psps=psps, &
         &rprimd=cryst%rprimd,ucvol=cryst%ucvol,xred=cryst%xred)

  end subroutine wfd_run_wannier



  !subroutine prepare_fft_grid(cryst, dtset)
  !  type(crystal_t) :: cryst
  !  type(dataset_type),intent(in) :: dtset
  !  type(pawfgr_type) :: pawfgr
  !  real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff,gsqcut_shp
  !  integer :: ngfftc(18),ngfftf(18)
  !  integer :: mgfftf,nfftf
  !  real(dp),parameter :: k0(3)=zero
  !  type(mpi_type) :: mpi_enreg

  !  ! input: dtset, k0, gmet
  !  call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
  !       gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=cryst%gmet,k0=k0)
  !  call print_ngfft(ngfftc, header='Coarse FFT mesh used for the wavefunctions')
  !  call print_ngfft(ngfftf, header='Dense FFT mesh used for densities and potentials')

  !  ! TODO hexu: this should be moved out of this subroutine
  !  ! Fake MPI_type for the sequential part.
  !  call initmpi_seq(mpi_enreg)
  !  call init_distribfft_seq(mpi_enreg%distribfft,'c',ngfftc(2),ngfftc(3),'all')
  !  call init_distribfft_seq(mpi_enreg%distribfft,'f',ngfftf(2),ngfftf(3),'all')
  !end subroutine prepare_fft_grid

  !! Is this a routine somewhere else?
  !subroutine prepare_paw(pawtab, cryst, psps, dtset)
  !  type(pseudopotential_type),intent(in) :: psps
  !  type(crystal_t), intent(inout) :: cryst
  !  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  !  type(dataset_type),intent(inout) :: dtset

  !  type(Pawrhoij_type),allocatable :: pawrhoij(:)
  !  integer :: cplex,cplex_dij,cplex_rhoij,ndij,nspden_rhoij,gnt_option

  ! call chkpawovlp(cryst%natom,cryst%ntypat,dtset%pawovlp,pawtab,cryst%rmet,cryst%typat,cryst%xred)

  ! cplex_dij=dtset%nspinor; cplex=1; ndij=1

  ! ABI_MALLOC(pawrhoij,(cryst%natom))
  ! ! TODO : Hexu: why spden and sporb is not present?
  ! call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
  !                           nspden=Dtset%nspden,spnorb=Dtset%pawspnorb,cpxocc=Dtset%pawcpxocc)
  ! call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,cryst%typat,pawtab=pawtab)

  ! ! Initialize values for several basic arrays
  ! gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

  ! ! Test if we have to call pawinit
  ! call paw_gencond(dtset,gnt_option,"test",call_pawinit)

  ! if (psp_gencond==1 .or. call_pawinit) then
  !   call timab(553,1,tsec)
  !   gsqcut_shp = two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
  !   call pawinit(dtset%effmass_free,gnt_option,gsqcut_shp,zero,dtset%pawlcutd,dtset%pawlmix,&
  !                psps%mpsang,dtset%pawnphi,cryst%nsym,dtset%pawntheta,pawang,Pawrad,&
  !                dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%ixc,dtset%usepotzero)
  !   call timab(553,2,tsec)

  !   ! Update internal values
  !   call paw_gencond(dtset,gnt_option,"save",call_pawinit)

  ! else
  !   if (pawtab(1)%has_kij  ==1) pawtab(1:cryst%ntypat)%has_kij  =2
  !   if (pawtab(1)%has_nabla==1) pawtab(1:cryst%ntypat)%has_nabla=2
  ! end if

  ! psps%n1xccc=MAXVAL(pawtab(1:cryst%ntypat)%usetcore)

  ! ! Initialize optional flags in pawtab to zero
  ! ! (Cannot be done in Pawinit since the routine is called only if some pars. are changed)
  ! pawtab(:)%has_nabla = 0
  ! pawtab(:)%usepawu   = 0
  ! pawtab(:)%useexexch = 0
  ! pawtab(:)%exchmix   =zero

  ! call setsym_ylm(cryst%gprimd,pawang%l_max-1,cryst%nsym,dtset%pawprtvol,cryst%rprimd,cryst%symrec,pawang%zarot)

  ! ! Initialize and compute data for DFT+U
  ! !paw_dmft%use_dmft=dtset%usedmft
  ! !call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%f4of2_sla,dtset%f6of2_sla,&
  ! !    .false.,dtset%jpawu,dtset%lexexch,dtset%lpawu,cryst%ntypat,pawang,dtset%pawprtvol,&
  ! !    Pawrad,pawtab,dtset%upawu,dtset%usedmft,dtset%useexexch,dtset%usepawu)
  ! !ABI_CHECK(paw_dmft%use_dmft==0,"DMFT not available")
  ! !call destroy_sc_dmft(paw_dmft)

  ! if (my_rank == master) call pawtab_print(pawtab, unit=std_out)

  ! ! Get Pawrhoij from the header of the WFK file.
  ! call pawrhoij_copy(wfk0_hdr%pawrhoij,pawrhoij)

  ! ! Variables/arrays related to the fine FFT grid.
  ! ABI_MALLOC(pawfgrtab,(cryst%natom))
  ! call pawtab_get_lsize(pawtab,l_size_atm,cryst%natom,cryst%typat)
  ! cplex=1
  ! call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat)
  ! ABI_FREE(l_size_atm)

  ! usexcnhat=maxval(pawtab(:)%usexcnhat)
  ! ! 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
  ! ! 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
  ! call wrtout(std_out,sjoin("using usexcnhat= ",itoa(usexcnhat)))
  ! !
  ! ! Identify parts of the rectangular grid where the density has to be calculated ===
  ! !optcut=0; optgr0=dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-dtset%pawstgylm
  ! !if (dtset%xclevel==2 .and. usexcnhat>0) optgr1=dtset%pawstgylm
  ! optcut=1; optgr0=1; optgr1=1; optgr2=1; optrad=1

  ! call nhatgrid(cryst%atindx1,cryst%gmet,cryst%natom,cryst%natom,cryst%nattyp,ngfftf,cryst%ntypat,&
  !              optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,cryst%rprimd,cryst%typat,cryst%ucvol,cryst%xred)

  ! call pawfgrtab_print(pawfgrtab,cryst%natom,unit=std_out,prtvol=dtset%pawprtvol)
  !end subroutine prepare_paw

  !subroutine wfd_prepare(wfk0_path, dtset, bks_mask, ngfft, comm)
  !  character(len=*), intent(inout) :: wfk0_path
  !  type(dataset_type),intent(in) :: dtset
  !  integer,intent(in) :: ngfft(18)
  !  logical,allocatable :: bks_mask(:,:,:)
  !  logical, allocatable :: keep_ur(:,:,:)
  !  integer,allocatable :: nband(:,:),wfd_istwfk(:)

  !  real(dp),pointer :: gs_eigen(:,:,:)
  !  type(hdr_type) :: wfk0_hdr
  !  type(crystal_t) :: cryst
  !  type(ebands_t) :: ebands

  !  integer :: nsppol, nspinor, nspden
  !  integer :: nkibz, mband, timrev
  !  real(dp) :: ecut


  !  integer :: comm, nprocs, my_rank
  !  integer :: ierr
  !  character(len=500) :: msg

  !  ! Copied from m_wfk_analysis
  !  comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

  !  if (my_rank == master) then
  !     ! Accept WFK file in Fortran or netcdf format.
  !     if (nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
  !        ABI_ERROR(sjoin("Cannot find GS WFK file:", ch10, msg))
  !     end if
  !  end if
  !  call xmpi_bcast(wfk0_path, master, comm, ierr)
  !  call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path))


  !  ! Costruct crystal and ebands from the GS WFK file.
  !  call wfk_read_eigenvalues(wfk0_path, gs_eigen, wfk0_hdr, comm) !,gs_occ)
  !  call wfk0_hdr%vs_dtset(dtset)



  !  cryst = wfk0_hdr%get_crystal()
  !  call cryst%print(header="Crystal structure from WFK file")

  !  ebands = ebands_from_hdr(wfk0_hdr, maxval(wfk0_hdr%nband), gs_eigen)
  !  nkibz = ebands%nkpt
  !  mband = ebands%mband
  !  timrev = kpts_timrev_from_kptopt(ebands%kptopt)

  !  nsppol = dtset%nsppol
  !  nspinor = ebands%nspinor
  !  nspden = dtset%nspden

  !  !call ebands_update_occ(ebands, spinmagntarget)
  !  call ebands_print(ebands,header="Ground state energies", prtvol=dtset%prtvol)
  !  ABI_FREE(gs_eigen)

  !  ! TODO : pawtab
  !  !call prepare_pawtab()


  !  ABI_MALLOC(nband, (nkibz, nsppol))
  !  ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
  !  ABI_MALLOC(keep_ur, (mband, nkibz ,nsppol))
  !  nband = mband; bks_mask = .False.; keep_ur = .False.


  !  ! build bsk_mask
  !  ! TODO: check Matteo ibte branch m_gstore gstore_t: parameters. fill_bks_mask. 
  !  ! TODO: adapt this to the need of Wannier  
  !  ! TODO: for k-point in this process, find their neighbors, and set the bsk_mask
  !  do mys=1,gams%my_nspins
  !     spin = gams%my_spins(mys)
  !     fs => fstab(spin)
  !     do ik_bz=1,fs%nkfs
  !        ! NOTE:  indkk_fs maps the kpt in bz to ibz.
  !        !listkk : legacy, slow but works. m_kpts.F90.  use_symrec=.False.
  !        ! mpw: m_gstore  get_kg, and get mpw. 
  !        !krank:
  !        ik_ibz = fs%indkk_fs(1, ik_bz)
  !        bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)
  !        bks_mask(bstart_k:bstart_k+nband_k-1, ik_ibz, spin) = .True.
  !     end do
  !  end do

  !  ! Impose istwfk = 1 for all k points. This is also done in respfn (see inkpts)
  !  ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
  !  ABI_MALLOC(wfd_istwfk, (nkibz))
  !  wfd_istwfk = 1

  !  call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
  !       nspden, nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
  !       dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

  !end subroutine wfd_prepare


  !subroutine wfd_wannier
  !  ! prepare

  !  ! Pieces of code to be used.

  !  ! To get the overlap.
  !  do isppol = 1, my_nsppol
  !     do ikpt =1, my_nkpt
  !        ! wfd%get_ug(ikpt, isppol)
  !        !call wfd_get_ug(Wfd, band, ik_ibz, spin, ug)
  !        wfd_ug_uk
  !        !call wfd_get_cprj(Wfd, band, ik_ibz, spin, Cryst, Cprj_out, sorted)
  !        ! paw_symcprj. (kmesh obj. ) ! see m_paw_sym
  !        !to init kmesh (skip. )
  !        do iknn =1, my_nntot(ikpt)
  !           ikpt2= ikpt_nn(iknn, ikpt)
  !           ! find the mapping to the ibz kpoint.
  !           ! use
  !           ! loop iband1
  !           ! loop iband2
  !           ! compute the overlaps
  !           ! sum up the overlaps
  !        end do
  !  end do
  !end subroutine wfd_wannier

  end module m_wfd_wannier

