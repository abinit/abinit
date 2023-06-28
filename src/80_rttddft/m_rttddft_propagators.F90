!!****m* ABINIT/m_rttddft_propagators
!! NAME
!!  m_rttddft_propagators
!!
!! FUNCTION
!!  Contains various propagators for the KS orbitals
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft_propagators

 use defs_basis
 use defs_abitypes,         only: MPI_type
 use defs_datatypes,        only: pseudopotential_type

 use m_bandfft_kpt,         only: bandfft_kpt, bandfft_kpt_type, &
                                & bandfft_kpt_set_ikpt,          &
                                & prep_bandfft_tabs
 use m_dtset,               only: dataset_type
 use m_energies,            only: energies_type, energies_init, energies_copy
 use m_gemm_nonlop,         only: make_gemm_nonlop
 use m_hamiltonian,         only: gs_hamiltonian_type, gspot_transgrid_and_pack
 use m_invovl,              only: make_invovl
 use m_kg,                  only: mkkin, mkkpg
 use m_mkffnl,              only: mkffnl
 use m_mpinfo,              only: proc_distrb_cycle
 use m_profiling_abi,       only: abimem_record
 use m_rttddft,             only: rttddft_init_hamiltonian
 use m_rttddft_exponential, only: rttddft_exp_taylor
 use m_rttddft_properties,  only: rttddft_calc_density, &
                                & rttddft_calc_occ,     &
                                & rttddft_calc_kin
 use m_rttddft_tdks,        only: tdks_type
 use m_specialmsg,          only: wrtout
 use m_xmpi,                only: xmpi_comm_rank, xmpi_sum, xmpi_max

 implicit none

 private
!!***

 public :: rttddft_propagator_er
 public :: rttddft_propagator_emr
!!***

contains
!!***

!!****f* m_rttddft/rttddft_propagator_er
!!
!! NAME
!!  rttddft_propagator_er
!!
!! FUNCTION
!!  Main subroutine to propagate the KS orbitals using
!!  the Exponential Rule (ER) propagator.
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  ham_k <type(gs_hamiltonian_type)> = Hamiltonian object
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!  calc_properties <logical> = logical governing the computation of
!!                              some properties (energies, occupations, eigenvalues..)
!!
!! NOTES
!!  Other propagators such as the Exponential Midpoint Rule (EMR)
!!  should usually be prefered over this one since the ER propagator
!!  alone violates time reversal symmetry. Using this propagator with
!!  the exponential approximated by Taylor expansion of order 1 leads
!!  to the famous Euler method which is fast and simple but unstable
!!  and thus insufficient for RT-TDDFT.
!!
!! SOURCE
subroutine rttddft_propagator_er(dtset, ham_k, istep, mpi_enreg, psps, tdks, calc_properties)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 logical,          optional, intent(in)    :: calc_properties
 type(dataset_type),         intent(inout) :: dtset
 type(gs_hamiltonian_type),  intent(inout) :: ham_k
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),    target, intent(inout) :: tdks

 !Local variables-------------------------------
 !scalars
 integer                        :: bdtot_index
 integer                        :: calc_forces
 integer                        :: dimffnl
 integer                        :: gemm_nonlop_ikpt_this_proc_being_treated
 integer                        :: iband
 integer                        :: ibg, icg
 integer                        :: ider, idir
 integer                        :: ierr, ilm
 integer                        :: ikpt, ikpt_loc, ikg
 integer                        :: isppol
 integer                        :: istwf_k
 integer                        :: me_distrb
 integer                        :: me_bandfft
 integer                        :: my_ikpt, my_nspinor
 integer                        :: nband_k, nband_k_mem
 integer                        :: npw_k, nkpg
 integer                        :: shift
 integer                        :: spaceComm_distrb
 integer                        :: n4, n5, n6
 logical                        :: with_vxctau
 logical                        :: lcalc_properties
 type(energies_type)            :: energies
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 !arrays
 integer,  allocatable          :: kg_k(:,:)
 real(dp), pointer              :: cg(:,:) => null()
 real(dp), pointer              :: cg0(:,:) => null()
 real(dp), allocatable          :: enl(:)
 real(dp), pointer              :: eig(:) => null()
 real(dp), allocatable          :: ffnl(:,:,:,:)
 real(dp), allocatable          :: kpg_k(:,:)
 real(dp)                       :: kpoint(3)
 real(dp), allocatable          :: kinpw(:)
 real(dp), pointer              :: occ(:) => null()
 real(dp), pointer              :: occ0(:) => null()
 real(dp), allocatable          :: ph3d(:,:,:)
 real(dp), allocatable          :: vlocal(:,:,:,:)
 real(dp), allocatable          :: vxctaulocal(:,:,:,:,:)
 real(dp), allocatable          :: ylm_k(:,:)
 logical                        :: lproperties(4)

! ***********************************************************************

 !Init MPI
 spaceComm_distrb=mpi_enreg%comm_cell
 if (mpi_enreg%paral_kgb==1) spaceComm_distrb=mpi_enreg%comm_kpt
 me_distrb=xmpi_comm_rank(spaceComm_distrb)

 !Do we calculate properties?
 !Governed by lproperties:
 !  lproperties(1) = compute energy contributions (kinetic)
 !  lproperties(2) = NL energy contribution in NC case
 !  lproperties(3) = eigenvalues
 !  lproperties(4) = occupations
 lproperties(:) = .false.
 lcalc_properties = .false.
 if (present(calc_properties)) then
   lcalc_properties = calc_properties
   if (lcalc_properties) then
      !compute energy contributions
      lproperties(1) = .true.
      !Init to zero different energies
      call energies_init(energies)
      energies%entropy=tdks%energies%entropy
      energies%e_corepsp=tdks%energies%e_corepsp
      energies%e_ewald=tdks%energies%e_ewald
      !including NL part in NC case?
      if (dtset%usepaw == 0) then
         lproperties(2) = .true.
      else
         ABI_MALLOC(enl,(0))
      end if
      !other properties
      ! eigenvalues
      if (dtset%prteig /= 0 .or. dtset%prtdos /= 0) then
         if (mod(istep-1,dtset%td_prtstr) == 0) then
            lproperties(3) = .true.
            tdks%eigen(:) = zero
         end if
      end if
      ! occupations
      lproperties(4) = .true.
      tdks%occ(:) = zero
   end if
 end if

 !Set "vtrial" and initialize the Hamiltonian
 call rttddft_init_hamiltonian(dtset,energies,ham_k,istep,mpi_enreg,psps,tdks)

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 n4=dtset%ngfft(4); n5=dtset%ngfft(5); n6=dtset%ngfft(6)
 ABI_MALLOC(vlocal,(n4,n5,n6,ham_k%nvloc))
 with_vxctau=(dtset%usekden/=0)
 if(with_vxctau) then
   ABI_MALLOC(vxctaulocal,(n4,n5,n6,ham_k%nvloc,4))
 end if
 if (dtset%ionmov/=0) then
   calc_forces=1
 else
   calc_forces=0
 end if

!FB: @MT Needed?
!has_vectornd = (with_vectornd .EQ. 1)
!if(has_vectornd) then
!  ABI_MALLOC(vectornd_pac,(n4,n5,n6,ham_k%nvloc,3))
!  vectornd_pac=zero
!end if

 icg=0; ibg=0
 bdtot_index=0

 !*** LOOP OVER SPINS
 do isppol=1,dtset%nsppol

   ikpt_loc=0
   ikg=0

   ! Set up local potential vlocal on the coarse FFT mesh from vtrial taking into account the spin.
   ! Also, continue to initialize the Hamiltonian.
   call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, tdks%nfftf, &
                               & dtset%nspden, ham_k%nvloc, 1, tdks%pawfgr, mpi_enreg, tdks%vtrial, vlocal)
   call ham_k%load_spin(isppol, vlocal=vlocal, with_nonlocal=.true.)

   if (with_vxctau) then
      call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, tdks%nfftf, &
                                  & dtset%nspden, ham_k%nvloc, 4, tdks%pawfgr, mpi_enreg, tdks%vxctau, vxctaulocal)
      call ham_k%load_spin(isppol, vxctaulocal=vxctaulocal)
   end if

!FB: @MT Needed?
!  ! if vectornd is present, set it up for addition to ham_k similarly to how it's done for
!  ! vtrial. Note that it must be done for the three Cartesian directions. Also, the following
!  ! code assumes explicitly and implicitly that nvloc = 1. This should eventually be generalized.
!  if(has_vectornd) then
!     do idir = 1, 3
!        ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
!        call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,idir))
!        call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
!        ABI_FREE(cgrvtrial)
!     end do
!     call ham_k%load_spin(isppol, vectornd=vectornd_pac)
!  end if

   !*** BIG FAT k POINT LOOP
   ikpt = 0
   do while (ikpt_loc < dtset%nkpt)

      ikpt_loc = ikpt_loc + 1
      ikpt = ikpt_loc
      my_ikpt = mpi_enreg%my_kpttab(ikpt)

      nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      if (mpi_enreg%paral_kgb==1) then
         nband_k_mem=mpi_enreg%bandpp
      else
         nband_k_mem=nband_k
      end if
      istwf_k=dtset%istwfk(ikpt)
      npw_k=tdks%npwarr(ikpt)

      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then
         bdtot_index=bdtot_index+nband_k
         cycle
      end if

      if (mpi_enreg%paral_kgb==1) my_bandfft_kpt => bandfft_kpt(my_ikpt)
      call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

      ABI_MALLOC(kg_k,(3,npw_k))
      ABI_MALLOC(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
      kg_k(:,1:npw_k)=tdks%kg(:,1+ikg:npw_k+ikg)
      if (psps%useylm==1) then
         do ilm=1,psps%mpsang*psps%mpsang
            ylm_k(1:npw_k,ilm)=tdks%ylm(1+ikg:npw_k+ikg,ilm)
         end do
      end if

      !** Set up the remaining k-dependent part of the Hamiltonian
      ! Kinetic energy - Compute (1/2) (2 Pi)**2 (k+G)**2:
      kpoint(:)=dtset%kptns(:,ikpt)
      ABI_MALLOC(kinpw,(npw_k))
      call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,tdks%gmet,kg_k,kinpw,kpoint,npw_k,0,0)
      ! Compute (k+G) vectors (only if useylm=1)
      nkpg=3*calc_forces*dtset%nloalg(3)
      ABI_MALLOC(kpg_k,(npw_k,nkpg))
      if ((mpi_enreg%paral_kgb/=1.or.istep<=tdks%first_step).and.nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
      ! Compute nonlocal form factors ffnl at all (k+G):
      ider=0;idir=0;dimffnl=1
      ABI_MALLOC(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
      if (mpi_enreg%paral_kgb/=1.or.istep<=tdks%first_step) then
        call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,tdks%gmet,tdks%gprimd, &
                  & ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,     &
                  & psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,       &
                  & psps%qgrid_ff,tdks%rmet,psps%usepaw,psps%useylm,ylm_k,tdks%ylmgr)
      end if

      !** Load k-dependent part in the Hamiltonian datastructure
      !**  - Compute 3D phase factors
      !**  - Prepare various tabs in case of band-FFT parallelism
      !**  - Load k-dependent quantities in the Hamiltonian
      ABI_MALLOC(ph3d,(2,npw_k,ham_k%matblk))
      call ham_k%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k, &
                      & ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=tdks%first_step),   &
                      & compute_gbound=(mpi_enreg%paral_kgb/=1))

      !** Load band-FFT tabs (transposed k-dependent arrays)
      if (mpi_enreg%paral_kgb==1) then
         if (istep<=tdks%first_step) call prep_bandfft_tabs(ham_k,ikpt,dtset%mkmem,mpi_enreg)
         call ham_k%load_k(npw_fft_k =my_bandfft_kpt%ndatarecv,    &
                         & gbound_k  =my_bandfft_kpt%gbound,       &
                         & kinpw_k   =my_bandfft_kpt%kinpw_gather, &
                         & kg_k      =my_bandfft_kpt%kg_k_gather,  &
                         & kpg_k     =my_bandfft_kpt%kpg_k_gather, &
                         & ffnl_k    =my_bandfft_kpt%ffnl_gather,  &
                         & ph3d_k    =my_bandfft_kpt%ph3d_gather)
      end if

      !** Build inverse of overlap matrix
      if(psps%usepaw == 1 .and. istep <= tdks%first_step) then
         call make_invovl(ham_k, dimffnl, ffnl, ph3d, mpi_enreg)
      end if

      ! Setup gemm_nonlop
      if (tdks%gemm_nonlop_use_gemm) then
         !set the global variable indicating to gemm_nonlop where to get its data from
         gemm_nonlop_ikpt_this_proc_being_treated = my_ikpt
         if (istep <= tdks%first_step) then
            !Init the arrays
            call make_gemm_nonlop(my_ikpt,ham_k%npw_fft_k,ham_k%lmnmax, &
                ham_k%ntypat, ham_k%indlmn, ham_k%nattyp, ham_k%istwf_k, &
                ham_k%ucvol,  ham_k%ffnl_k, &
                ham_k%ph3d_k, ham_k%kpt_k, ham_k%kg_k, ham_k%kpg_k)
         end if
      end if

      !** Compute the exp[(S^{-1})H]*cg using Taylor expansion to approximate the exponential
      cg => tdks%cg(:,icg+1:icg+nband_k*npw_k*my_nspinor)
      ! Compute properties "on-the-fly" if required
      if (lcalc_properties) then
         cg0 => tdks%cg0(:,icg+1:icg+nband_k*npw_k*my_nspinor)
         occ => tdks%occ(bdtot_index+1:bdtot_index+nband_k)
         occ0 => tdks%occ0(bdtot_index+1:bdtot_index+nband_k)
         if (dtset%paral_kgb /= 1) then
            shift = bdtot_index
         else
            me_bandfft = xmpi_comm_rank(mpi_enreg%comm_band)
            shift = bdtot_index+me_bandfft*mpi_enreg%bandpp
         end if
         ! kinetic energy
         if (lproperties(1)) then
            call rttddft_calc_kin(energies%e_kinetic,cg,dtset,ham_k,nband_k,npw_k,my_nspinor, &
                                & occ0,dtset%wtk(ikpt),mpi_enreg,my_bandfft_kpt)
         end if
         ! for NL PSP part
         if (lproperties(2)) then
            ABI_MALLOC(enl,(mpi_enreg%bandpp))
            enl = zero
         end if
         ! for eigenvalues
         if (lproperties(3)) then
            eig => tdks%eigen(1+shift:nband_k_mem+shift)
         end if
         ! occupations
         if (lproperties(4)) then
            !note that occupations are computed at istep-1 like energies and eigenvalues
            call rttddft_calc_occ(cg,cg0,dtset,ham_k,ikpt,ibg,isppol,mpi_enreg, &
                                & nband_k,npw_k,my_nspinor,occ,occ0,tdks)
         end if

         !Propagate cg and compute the requested properties
         call rttddft_exp_taylor(cg,dtset,ham_k,mpi_enreg,nband_k,npw_k,my_nspinor,enl=enl,eig=eig)

         !Finish computing NL PSP part for NC PSP
         if (lproperties(2)) then
            do iband = 1, mpi_enreg%bandpp
               energies%e_nlpsp_vfock=energies%e_nlpsp_vfock+dtset%wtk(ikpt)*tdks%occ0(shift+iband)*enl(iband)
            end do
            ABI_FREE(enl)
         end if
      else
         !Propagate cg only
         call rttddft_exp_taylor(cg,dtset,ham_k,mpi_enreg,nband_k,npw_k,my_nspinor)
      end if

      ABI_FREE(kg_k)
      ABI_FREE(ylm_k)
      ABI_FREE(kpg_k)
      ABI_FREE(kinpw)
      ABI_FREE(ffnl)
      ABI_FREE(ph3d)

      bdtot_index = bdtot_index+nband_k

      !** Also shift array memory if dtset%mkmem/=0
      if (dtset%mkmem/=0) then
         ibg=ibg+nband_k_mem
         icg=icg+npw_k*my_nspinor*nband_k
         ikg=ikg+npw_k
      end if

   end do !nkpt

 end do !nsppol

 ABI_FREE(vlocal)
 if(dtset%usekden/=0) then
   ABI_FREE(vxctaulocal)
 end if

 !Keep the computed energies in memory
 if (lcalc_properties) then
   call xmpi_sum(energies%e_kinetic,mpi_enreg%comm_kptband,ierr)
   call xmpi_sum(energies%e_nlpsp_vfock,mpi_enreg%comm_kptband,ierr)
   call energies_copy(energies,tdks%energies)
   if (lproperties(3)) call xmpi_sum(tdks%eigen,mpi_enreg%comm_kptband,ierr)
   if (lproperties(4)) call xmpi_sum(tdks%occ,mpi_enreg%comm_kpt,ierr)
 end if

end subroutine rttddft_propagator_er
!!***

!!****f* m_rttddft/rttddft_propagator_emr
!!
!! NAME
!!  rttddft_propagator_emr
!!
!! FUNCTION
!!  Main subroutine to propagate the KS orbitals using
!!  the Exponential Midpoint Rule (EMR) propagator
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ham_k <type(gs_hamiltonian_type)> = Hamiltonian object
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! NOTES
!!  This propagator is time reversible
!!  (if H(t+dt/2) and the exponential are computed exactly).
!!
!! SOURCE
subroutine rttddft_propagator_emr(dtset, ham_k, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(gs_hamiltonian_type),  intent(inout) :: ham_k
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks

 !Local variables-------------------------------
 integer :: me
 !scalars
 character(len=500) :: msg
 integer            :: ics
 integer            :: ierr
 logical            :: lconv
 !arrays
 real(dp)           :: cg(SIZE(tdks%cg(:,1)),SIZE(tdks%cg(1,:)))
 real(dp)           :: diff(SIZE(tdks%cg(:,1)),SIZE(tdks%cg(1,:)))
 real(dp)           :: max_diff(2)

! ***********************************************************************

 cg(:,:) = tdks%cg(:,:) !Psi(t)

 !** Predictor step
 ! predict psi(t+dt) using ER propagator
 call rttddft_propagator_er(dtset,ham_k,istep,mpi_enreg,psps,tdks,calc_properties=.true.)
 ! for convergence check
 diff = tdks%cg
 ! estimate psi(t+dt/2) = (psi(t)+psi(t+dt))/2
 tdks%cg(:,:) = 0.5_dp*(tdks%cg(:,:)+cg(:,:))
 ! calc associated density at t+dt/2
 call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)
 ! go back to time t ..
 tdks%cg(:,:) = cg(:,:)
 ! .. and evolve psi(t) using the EMR propagator with the estimated density at t+dt/2
 call rttddft_propagator_er(dtset,ham_k,istep,mpi_enreg,psps,tdks)

 ! check convergence
 diff = abs(diff-tdks%cg)
 me = xmpi_comm_rank(mpi_enreg%comm_world)
 call xmpi_max(maxval(diff(1,:)),max_diff(1),mpi_enreg%comm_world,ierr)
 call xmpi_max(maxval(diff(2,:)),max_diff(2),mpi_enreg%comm_world,ierr)
 lconv = (max_diff(1) < dtset%td_scthr .and. max_diff(2) < dtset%td_scthr)
 ics = 0
 if (mpi_enreg%me == 0) then
   write(msg,'(a,a,i3,a,3(es8.2,1x),l1,a)') ch10, 'SC Step', ics, ' - ', max_diff(1), max_diff(2), &
                                          & dtset%td_scthr, lconv, ch10
   if (do_write_log) call wrtout(std_out,msg)
 end if
 if (.not. lconv) then
   !** Corrector steps
   do ics = 1, dtset%td_scnmax
      ! for convergence check
      diff = tdks%cg
      ! estimate psi(t+dt/2) = (psi(t)+psi(t+dt))/2
      tdks%cg(:,:) = 0.5_dp*(tdks%cg(:,:)+cg(:,:))
      ! calc associated density at t+dt/2
      call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)
      ! Go back to time t ..
      tdks%cg(:,:) = cg(:,:)
      ! .. and evolve psi(t) using estimated density at t+dt/2
      call rttddft_propagator_er(dtset,ham_k,istep,mpi_enreg,psps,tdks)
      ! check convergence
      diff = abs(diff-tdks%cg)
      me = xmpi_comm_rank(mpi_enreg%comm_world)
      call xmpi_max(maxval(diff(1,:)),max_diff(1),mpi_enreg%comm_world,ierr)
      call xmpi_max(maxval(diff(2,:)),max_diff(2),mpi_enreg%comm_world,ierr)
      lconv = (max_diff(1) < dtset%td_scthr .and. max_diff(2) < dtset%td_scthr)
      if (mpi_enreg%me == 0) then
         write(msg,'(a,a,i3,a,3(es8.2,1x),l1,a)') ch10, 'SC Step', ics, ' - ', max_diff(1), max_diff(2), &
                                                & dtset%td_scthr, lconv, ch10
         if (do_write_log) call wrtout(std_out,msg)
      end if
      if (lconv) exit
   end do
 end if

 if (lconv) then
   write(msg,'(a,i4,a)') "Converged after ", ics, " self-consistent corrector steps"
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)
 else
   write(msg,'(a)') "Reached maximum number of corrector steps before convergence"
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)
 end if

end subroutine rttddft_propagator_emr
!!***

end module m_rttddft_propagators
!!***
