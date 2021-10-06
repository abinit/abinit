!!****m* ABINIT/m_rttddft_propagators
!! NAME
!!  m_rttddft_propagators
!!
!! FUNCTION
!!  Contains various propagators for the KS orbitals
!!
!! COPYRIGHT
!!  Copyright (C) 2021 ABINIT group (FB, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  m_rttddft_propagate
!!
!! CHILDREN
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
                                & bandfft_kpt_set_ikpt, &
                                & bandfft_kpt_get_ikpt, &
                                & prep_bandfft_tabs
 use m_cgtools,             only: dotprod_g
 use m_dtset,               only: dataset_type
 use m_energies,            only: energies_type, energies_init, energies_copy
 use m_gemm_nonlop,         only: make_gemm_nonlop
 use m_hamiltonian,         only: gs_hamiltonian_type, gspot_transgrid_and_pack
 use m_invovl,              only: make_invovl
 use m_kg,                  only: mkkin, mkkpg
 use m_mkffnl,              only: mkffnl
 use m_mpinfo,              only: proc_distrb_cycle
 use m_rttddft,             only: rttddft_init_hamiltonian, &
                                & rttddft_calc_density
 use m_rttddft_exponential, only: rttddft_exp_taylor
 use m_rttddft_types,       only: tdks_type 
 use m_spacepar,            only: meanvalue_g
 use m_specialmsg,          only: wrtout
 use m_xmpi,                only: xmpi_comm_rank

 implicit none

 private
!!***

 public :: rttddft_propagator_er
 public :: rttddft_propagator_emr
!!***

contains 

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
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)> = Hamiltonian object
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Other propagators such as the Exponential Midpoint Rule (EMR) 
!!  propagator should be prefered over this one since the ER propagator 
!!  violates time reversal symmetry. Using this propagator with the 
!!  exponential approximated by Taylor expansion of order 2 leads to 
!!  the famous Euler method which is fast and simple but extremely 
!!  unstable and thus usually insufficient for RT-TDDFT.
!!  
!!
!! PARENTS
!!  m_rttddft_propagate/rttddft_propagate_ele
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagator_er(dtset, gs_hamk, istep, mpi_enreg, psps, tdks, store_energies)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 logical,          optional, intent(in)    :: store_energies
 type(dataset_type),         intent(inout) :: dtset
 type(gs_hamiltonian_type),  intent(inout) :: gs_hamk
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks
 
 !Local variables-------------------------------
 !scalars
 integer                        :: bdtot_index
 integer                        :: blocksize
 integer                        :: calc_forces
 integer                        :: dimffnl
 integer                        :: gemm_nonlop_ikpt_this_proc_being_treated
 integer                        :: ibg, icg
 integer                        :: iblocksize, iblock, iband
 integer                        :: ider, idir, ilm
 integer                        :: ikpt, ikpt_loc, ikg
 integer                        :: ikpt_this_proc
 integer                        :: istwf_k
 integer                        :: me_distrb
 integer                        :: my_ikpt, my_nspinor
 integer                        :: nband_k, nband_cprj_k
 integer                        :: nblockbd
 integer                        :: npw_k, nkpg
 integer                        :: npw, nband
 integer                        :: spaceComm_distrb
 integer                        :: shift
 integer                        :: isppol
 integer                        :: n4, n5, n6
 logical                        :: with_vxctau
 real(dp)                       :: ar
 real(dp)                       :: dprod_i
 real(dp)                       :: enlx
 type(energies_type)            :: energies
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 !arrays
 integer,allocatable            :: kg_k(:,:)
 real(dp),allocatable           :: ffnl(:,:,:,:)
 real(dp),allocatable           :: gvnlxc(:,:)
 real(dp),allocatable           :: kpg_k(:,:)
 real(dp)                       :: kpoint(3)
 real(dp),allocatable           :: kinpw(:)
 real(dp),allocatable           :: occ_k(:)
 real(dp),allocatable           :: ph3d(:,:,:)
 real(dp),allocatable           :: vlocal(:,:,:,:)
 real(dp),allocatable           :: vxctaulocal(:,:,:,:,:)
 real(dp),allocatable           :: ylm_k(:,:)
 
! ***********************************************************************

 !Init MPI
 spaceComm_distrb=mpi_enreg%comm_cell
 me_distrb=xmpi_comm_rank(spaceComm_distrb)

 !Init to zero different energies
 call energies_init(energies)
 energies%e_corepsp=tdks%energies%e_corepsp
 energies%e_ewald=tdks%energies%e_ewald

 !Set "vtrial" and intialize the Hamiltonian
 call rttddft_init_hamiltonian(dtset,energies,gs_hamk,istep,mpi_enreg,psps,tdks)

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 n4=dtset%ngfft(4); n5=dtset%ngfft(5); n6=dtset%ngfft(6)
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamk%nvloc))
 with_vxctau=(dtset%usekden/=0)
 if(with_vxctau) then
   ABI_MALLOC(vxctaulocal,(n4,n5,n6,gs_hamk%nvloc,4))
 end if
 if (dtset%ionmov/=0) then
   calc_forces=1
 else
   calc_forces=0
 end if
 bdtot_index=0

!FB: @MT Needed?
!has_vectornd = (with_vectornd .EQ. 1)
!if(has_vectornd) then
!  ABI_MALLOC(vectornd_pac,(n4,n5,n6,gs_hamk%nvloc,3))
!  vectornd_pac=zero
!end if

 icg=0

 !*** LOOP OVER SPINS
 do isppol=1,dtset%nsppol

   ikpt_loc=0
   ikg=0

   ! Set up local potential vlocal on the coarse FFT mesh from vtrial taking into account the spin.
   ! Also, continue to initialize the Hamiltonian.
   call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, tdks%nfftf, &
                               & dtset%nspden, gs_hamk%nvloc, 1, tdks%pawfgr, mpi_enreg, tdks%vtrial, vlocal)
   call gs_hamk%load_spin(isppol, vlocal=vlocal, with_nonlocal=.true.)

   if (with_vxctau) then
      call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, tdks%nfftf, &
                                  & dtset%nspden, gs_hamk%nvloc, 4, tdks%pawfgr, mpi_enreg, tdks%vxctau, vxctaulocal)
      call gs_hamk%load_spin(isppol, vxctaulocal=vxctaulocal)
   end if

!FB: @MT Needed?
!  ! if vectornd is present, set it up for addition to gs_hamk similarly to how it's done for
!  ! vtrial. Note that it must be done for the three Cartesian directions. Also, the following
!  ! code assumes explicitly and implicitly that nvloc = 1. This should eventually be generalized.
!  if(has_vectornd) then
!     do idir = 1, 3
!        ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
!        call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,idir))
!        call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
!        ABI_FREE(cgrvtrial)
!     end do
!     call gs_hamk%load_spin(isppol, vectornd=vectornd_pac)
!  end if

   !*** BIG FAT k POINT LOOP
   ikpt = 0
   do while (ikpt_loc < dtset%nkpt)

      ikpt_loc = ikpt_loc + 1
      ikpt = ikpt_loc
      my_ikpt = mpi_enreg%my_kpttab(ikpt)

      nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      nband_cprj_k=nband_k/mpi_enreg%nproc_band
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

      kpoint(:)=dtset%kptns(:,ikpt)

      !** Set up the remaining k-dependent part of the Hamiltonian
      ! Kinetic energy - Compute (1/2) (2 Pi)**2 (k+G)**2:
      ABI_MALLOC(kinpw,(npw_k))
      call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,tdks%gmet,kg_k,kinpw,kpoint,npw_k,0,0)
      ! Compute (k+G) vectors (only if useylm=1)
      nkpg=3*calc_forces*dtset%nloalg(3)
      ABI_MALLOC(kpg_k,(npw_k,nkpg))
      if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
      ! Compute nonlocal form factors ffnl at all (k+G):
      ider=0;idir=0;dimffnl=1
      ABI_MALLOC(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
      if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
         call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,tdks%gmet,tdks%gprimd, &
                   & ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,     &
                   & psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,       &
                   & psps%qgrid_ff,tdks%rmet,psps%usepaw,psps%useylm,ylm_k,tdks%ylmgr)
      end if

      !** Load k-dependent part in the Hamiltonian datastructure
      !**  - Compute 3D phase factors
      !**  - Prepare various tabs in case of band-FFT parallelism
      !**  - Load k-dependent quantities in the Hamiltonian
      ABI_MALLOC(ph3d,(2,npw_k,gs_hamk%matblk))
      call gs_hamk%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,kinpw_k=kinpw,kg_k=kg_k,         &
                        & kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1), &
                        & compute_gbound=(mpi_enreg%paral_kgb/=1))

      !** Load band-FFT tabs (transposed k-dependent arrays)
      if (mpi_enreg%paral_kgb==1) then
         if (istep<=1) call prep_bandfft_tabs(gs_hamk,ikpt,dtset%mkmem,mpi_enreg)
         call gs_hamk%load_k(npw_fft_k=my_bandfft_kpt%ndatarecv,    &
                           & gbound_k =my_bandfft_kpt%gbound,       &
                           & kinpw_k  =my_bandfft_kpt%kinpw_gather, &
                           & kg_k     =my_bandfft_kpt%kg_k_gather,  &
                           & kpg_k    =my_bandfft_kpt%kpg_k_gather, &
                           & ffnl_k   =my_bandfft_kpt%ffnl_gather,  &
                           & ph3d_k   =my_bandfft_kpt%ph3d_gather)
      end if
   
      !** Build inverse of overlap matrix
      if(psps%usepaw == 1) then
         call make_invovl(gs_hamk, dimffnl, ffnl, ph3d, mpi_enreg)
      end if

      ! Setup gemm_nonlop
      if (tdks%gemm_nonlop_use_gemm) then
         !set the global variable indicating to gemm_nonlop where to get its data from
         gemm_nonlop_ikpt_this_proc_being_treated = my_ikpt
         if (istep <= 1) then
            !Init the arrays
            call make_gemm_nonlop(my_ikpt,gs_hamk%npw_fft_k,gs_hamk%lmnmax,gs_hamk%ntypat,        &
                               & gs_hamk%indlmn, gs_hamk%nattyp, gs_hamk%istwf_k, gs_hamk%ucvol, &
                               & gs_hamk%ffnl_k,gs_hamk%ph3d_k)
         end if
      end if
     
      if (dtset%paral_kgb == 1) then
         ikpt_this_proc = bandfft_kpt_get_ikpt()
         npw = bandfft_kpt(ikpt_this_proc)%ndatarecv
         nband = mpi_enreg%bandpp
      else
         npw = npw_k
         nband = nband_k
      end if

      ABI_MALLOC(gvnlxc, (2, npw*my_nspinor*nband))
      
      !** Compute the exp[(S^{-1})H]*cg using taylor expansion to approximate the exponential
      call rttddft_exp_taylor(tdks%cg(:,icg+1:),dtset%dtele,dtset,gs_hamk,gvnlxc,1,mpi_enreg,nband,npw,my_nspinor)

      !** Compute different contribution to the energy if needed
      if (PRESENT(store_energies)) then
         if (store_energies) then
            ABI_MALLOC(occ_k,(nband_k))
            occ_k(:)=tdks%occ(1+bdtot_index:nband_k+bdtot_index)
            nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
            blocksize=nband_k/nblockbd
            do iblock=1,nblockbd
               do iblocksize=1,blocksize
                  iband=(iblock-1)*blocksize+iblocksize
                  shift = npw_k*my_nspinor*(iband-1)
                  !Compute kinetic energy contribution 
                  call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
                                 & tdks%cg(:,1+shift+icg:shift+npw_k*my_nspinor+icg),&
                                 & tdks%cg(:,1+shift+icg:shift+npw_k*my_nspinor+icg),0)
                  energies%e_kinetic = energies%e_kinetic + dtset%wtk(ikpt)*occ_k(iband)*ar
                  !Compute non local psp energy contribution
                  if (gs_hamk%usepaw /= 1) then
                     call dotprod_g(enlx,dprod_i,gs_hamk%istwf_k,npw_k*my_nspinor,1,&
                                  & tdks%cg(:,1+shift+icg:shift+npw_k*my_nspinor+icg),&
                                  & gvnlxc(:,1+shift:shift+npw_k*my_nspinor),&
                                  & mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
                     energies%e_nlpsp_vfock=energies%e_nlpsp_vfock+dtset%wtk(ikpt)*occ_k(iband)*enlx                    
                  end if
               end do
            end do
            ABI_FREE(occ_k)
         end if
      end if

      ABI_FREE(gvnlxc)
      ABI_FREE(kg_k)
      ABI_FREE(ylm_k)
      ABI_FREE(kpg_k)
      ABI_FREE(kinpw)
      ABI_FREE(ffnl)
      ABI_FREE(ph3d)

      !** Also shift array memory if dtset%mkmem/=0
      if (dtset%mkmem/=0) then
         ibg=ibg+my_nspinor*nband_cprj_k
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
 if (PRESENT(store_energies)) then
   if (store_energies) call energies_copy(energies,tdks%energies)
 end if

 !FB: There should be some MPI reduce here I think.. ?!

 end subroutine rttddft_propagator_er

!!****f* m_rttddft/rttddft_propagator_emr
!!
!! NAME
!!  rttddft_propagator_emr
!!
!! FUNCTION
!!  Main subroutine to propagate the KS orbitals using 
!!  the Euler propagator
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)> = Hamiltonian object
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  This propagator is time reversible (if H(t+dt/2) and the exponential 
!!  would be computed exactly).
!!
!! PARENTS
!!  m_rttddft_propagate/rttddft_propagate_ele
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagator_emr(dtset, gs_hamk, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(gs_hamiltonian_type),  intent(inout) :: gs_hamk
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 integer              :: ics
 logical              :: lconv
 !arrays
 real(dp)             :: conv(2)
 real(dp)             :: cg(SIZE(tdks%cg(:,1)),SIZE(tdks%cg(1,:)))
 
! ***********************************************************************

 cg(:,:) = tdks%cg(:,:) !Psi(t)

 !** Predictor step
 ! predict psi(t+dt) using ER propagator
 call rttddft_propagator_er(dtset,gs_hamk,istep,mpi_enreg,psps,tdks,store_energies=.true.)
 ! for convergence check
 conv(1) = sum(abs(tdks%cg(1,:))); conv(2) = sum(abs(tdks%cg(2,:)))
 ! estimate psi(t+dt/2) = (psi(t)+psi(t+dt))/2
 tdks%cg(:,:) = 0.5_dp*(tdks%cg(:,:)+cg(:,:))
 ! calc associated density at t+dt/2
 call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)
 ! go back to time t ..
 tdks%cg(:,:) = cg(:,:)
 ! .. and evolve psi(t) using the EMR propagator with the estimated density at t+dt/2
 call rttddft_propagator_er(dtset,gs_hamk,istep,mpi_enreg,psps,tdks)
 
 ics = 0
 print*, "SC Step", ics, " - ", 100*(conv(1)-sum(abs(tdks%cg(1,:))))/conv(1), &
         & 100*(conv(2)-sum(abs(tdks%cg(2,:))))/conv(2), dtset%td_scthr, lconv
 ! Check convergence
 if (100*abs(conv(1)-sum(abs(tdks%cg(1,:))))/conv(1) < dtset%td_scthr .and. &
   & 100*abs(conv(2)-sum(abs(tdks%cg(2,:))))/conv(2) < dtset%td_scthr) then
   lconv = .true.
 else
   lconv = .false.
   !** Corrector steps
   do ics = 1, dtset%td_scnmax
      conv(1) = sum(abs(tdks%cg(1,:))); conv(2) = sum(abs(tdks%cg(2,:)))
      ! estimate psi(t+dt/2) = (psi(t)+psi(t+dt))/2
      tdks%cg(:,:) = 0.5_dp*(tdks%cg(:,:)+cg(:,:))
      ! calc associated density at t+dt/2
      call rttddft_calc_density(dtset,mpi_enreg,psps,tdks)
      ! start back from time t
      tdks%cg(:,:) = cg(:,:)
      ! evolve psi(t) using estimated density at t+dt/2
      call rttddft_propagator_er(dtset,gs_hamk,istep,mpi_enreg,psps,tdks)
      ! check convergence
      if (100*abs(conv(1)-sum(abs(tdks%cg(1,:))))/conv(1) < dtset%td_scthr .and. &
        & 100*abs(conv(2)-sum(abs(tdks%cg(2,:))))/conv(2) < dtset%td_scthr) then
         lconv = .true.
         print*, "SC Step", ics, " - ", 100*(conv(1)-sum(abs(tdks%cg(1,:))))/conv(1), &
               & 100*(conv(2)-sum(abs(tdks%cg(2,:))))/conv(2), dtset%td_scthr, lconv
         exit
     else
         print*, "SC Step", ics, " - ", 100*(conv(1)-sum(abs(tdks%cg(1,:))))/conv(1), &
               & 100*(conv(2)-sum(abs(tdks%cg(2,:))))/conv(2), dtset%td_scthr, lconv
        conv(1) = sum(abs(tdks%cg(1,:)))
        conv(2) = sum(abs(tdks%cg(2,:)))
     end if
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

end module m_rttddft_propagators
!!***
