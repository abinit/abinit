!!****m* ABINIT/m_rttddft_propagate
!! NAME
!!  m_rttddft_propagate
!!
!! FUNCTION
!!  Contains various subroutines to propagate the KS 
!!  orbitals and potentially also the nuclei in RTTDDFT
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
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft_propagate

 use defs_basis
 use defs_abitypes,   only: MPI_type
 use defs_datatypes,  only: pseudopotential_type
 
 use m_bandfft_kpt,   only: bandfft_kpt, bandfft_kpt_type, bandfft_kpt_set_ikpt, &
                            bandfft_kpt_get_ikpt, prep_bandfft_tabs
 use m_cgtools,       only: dotprod_g
 use m_dtset,         only: dataset_type
 use m_getghc,        only: getghc
 use m_gemm_nonlop,   only: make_gemm_nonlop
 use m_hamiltonian,   only: gs_hamiltonian_type, gspot_transgrid_and_pack
 use m_invovl,        only: make_invovl, apply_invovl
 use m_kg,            only: mkkin, mkkpg
 use m_mkffnl,        only: mkffnl
 use m_pawcprj,       only: pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_prep_kgb,      only: prep_getghc
 use m_rttddft,       only: rttddft_setup_ele_step
 use m_rttddft_types, only: tdks_type 
 use m_specialmsg,    only: wrtout
 use m_symtk,         only: symmetrize_xred
 use m_xmpi,          only: xmpi_sum

 implicit none

 private
!!***

 public :: rttddft_propagate_ele
 public :: rttddft_propagate_nuc

contains 

!!****f* m_rttddft/rttddft_propagate_ele
!!
!! NAME
!!  rttddft_propagate_ele
!!
!! FUNCTION
!!  Main subroutine to propagate time-dependent KS orbitals
!!
!! INPUTS
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagate_ele(tdks, dtset, istep, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),            intent(inout) :: tdks
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 character(len=500)             :: msg
 integer                        :: calc_forces
 integer                        :: dimffnl
 integer                        :: gemm_nonlop_ikpt_this_proc_being_treated
 integer                        :: ibg, icg
 integer                        :: ider, idir, ilm
 integer                        :: ikpt, ikpt_loc, ikg
 integer                        :: istwf_k
 integer                        :: my_ikpt, my_nspinor
 integer                        :: nband_k, nband_cprj_k
 integer                        :: npw_k, nkpg
 integer                        :: isppol
 integer                        :: n4, n5, n6
 logical                        :: with_vxctau
 type(gs_hamiltonian_type)      :: gs_hamk
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 !arrays
 integer,allocatable            :: kg_k(:,:)
 real(dp),allocatable           :: ffnl(:,:,:,:)
 real(dp),allocatable           :: kpg_k(:,:)
 real(dp)                       :: kpoint(3)
 real(dp),allocatable           :: kinpw(:)
 real(dp),allocatable           :: ph3d(:,:,:)
 real(dp),allocatable           :: vlocal(:,:,:,:)
 real(dp),allocatable           :: vxctaulocal(:,:,:,:,:)
 real(dp),allocatable           :: ylm_k(:,:)
 
! ***********************************************************************

 write(msg,'(2a,i5,a)') ch10,'--- Iteration',istep,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 ! Init/Update various quantities before performing propagating KS orbitals
 call rttddft_setup_ele_step(tdks,dtset,gs_hamk,istep,mpi_enreg,psps)

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

!FB: Needed?
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

!FB: Needed?
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
      npw_k=tdks%npwarr(ikpt)
      istwf_k=dtset%istwfk(ikpt)

      if (mpi_enreg%paral_kgb==1) my_bandfft_kpt => bandfft_kpt(my_ikpt)
      call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

       npw_k=tdks%npwarr(ikpt)
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
      if (dtset%usefock==1) then
         if(tdks%fock%fock_common%use_ACE/=0) then
            call gs_hamk%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,kinpw_k=kinpw,kg_k=kg_k, &
                              & kpg_k=kpg_k,ffnl_k=ffnl,fockACE_k=tdks%fock%fockACE(ikpt,isppol),ph3d_k=ph3d,  &
                              & compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1),compute_gbound=(mpi_enreg%paral_kgb/=1))
         end if
      else
         call gs_hamk%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,kinpw_k=kinpw,kg_k=kg_k,         &
                           & kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1), &
                           & compute_gbound=(mpi_enreg%paral_kgb/=1))
      end if

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
   
      !FB: Seems like this should be done at every step?
      !** Build inverse of overlap matrix
      if(psps%usepaw == 1 .and. istep <= 1) then
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
     
     !FB: Needed?
     !! Update the value of ikpt,isppol in fock_exchange and allocate the memory space to perform HF calculation.
     !if (usefock) call fock_updateikpt(fock%fock_common,ikpt,isppol)
     !if (psps%usepaw==1 .and. usefock) then
     !  if ((fock%fock_common%optfor).and.(usefock_ACE==0)) fock%fock_common%forces_ikpt=zero
     !end if
      
      !FB: Here we should now call the propagator to evolve the cg
      !1. Apply H (get_ghc) and possibly S^-1 (apply_invovl) as many times as needed
      call rttddft_update_cg_k(tdks%cg(:,icg+1:),dtset,gs_hamk,mpi_enreg,nband_k,npw_k,my_nspinor)

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

 !FB: the cprj and rhoij should also be updated once the cg have changed..
 !FB: probably to be done after propagating?

 end subroutine rttddft_propagate_ele

!!****f* m_rttddft/rttddft_update_cg_k
!!
!! NAME
!!  rttddft_update_cg_k
!!
!! FUNCTION
!!  Update the values of the cg coefficients at point k according to the propagator
!!
!! INPUTS
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_propagate_ele
!!
!! CHILDREN
!!
!! SOURCE
 subroutine rttddft_update_cg_k(cg,dtset,gs_hamk,mpi_enreg,nband_k,npw_k,nspinor)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)    :: nband_k
 integer,                   intent(in)    :: npw_k
 integer,                   intent(in)    :: nspinor
 type(dataset_type),        intent(in)    :: dtset
 type(gs_hamiltonian_type), intent(inout) :: gs_hamk
 type(MPI_type),            intent(inout) :: mpi_enreg
 !arrays
 real(dp),                  intent(inout) :: cg(:,:)
 
 !Local variables-------------------------------
 !scalars
 integer                         :: iband, icg, ierr
 integer                         :: ikpt_this_proc
 integer                         :: nband, npw
 integer                         :: shift
 integer                         :: sij_opt,cpopt
 integer                         :: tim_getghc = 5
 logical                         :: paw
 !arrays
 integer,            allocatable :: nline_bands(:)
 real(dp)                        :: cg_old(2)
 real(dp)                        :: dprod_r, dprod_i
 real(dp)                        :: eig(nband_k)
 real(dp), pointer               :: exp_operator(:,:)
 real(dp), target,   allocatable :: ghc(:,:)
 real(dp), target,   allocatable :: gsm1hc(:,:)
 real(dp),           allocatable :: gsc(:,:)
 real(dp),           allocatable :: gvnlxc(:,:)
 real(dp),           allocatable :: resids(:), residvec(:,:)
 type(pawcprj_type), allocatable :: cwaveprj(:,:)

 if (dtset%paral_kgb == 1) then
   ikpt_this_proc = bandfft_kpt_get_ikpt()
   npw = bandfft_kpt(ikpt_this_proc)%ndatarecv
   nband = mpi_enreg%bandpp
 else
    npw = npw_k
    nband = nband_k
 end if
 
 paw = gs_hamk%usepaw == 1
 if(paw) then
   ABI_MALLOC(cwaveprj, (gs_hamk%natom,nspinor*nband))
   call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
   !FB: Not sure about the values of sij_opt and cpopt
   sij_opt = 1 ! recompute S
   cpopt = 0   ! save cprojs
 else
    ABI_MALLOC(cwaveprj,(0,0))
   sij_opt = 0
   cpopt = -1
 end if

 !FB: Probably not needed to create so many arrrays
 ABI_MALLOC(ghc,    (2, npw*nspinor*nband))
 ABI_MALLOC(gsc,    (2, npw*nspinor*nband))
 ABI_MALLOC(gvnlxc, (2, npw*nspinor*nband))
 ABI_MALLOC(gsm1hc, (2, npw*nspinor*nband))

 !** Apply Hamiltonian to the cg (computes <G|H|C>)
 !In paral_kgb case, we are here in the representation where 
 !the procs have all the bands but only some G-points 
 if (dtset%paral_kgb /= 1) then
   call getghc(cpopt,cg,cwaveprj,ghc,gsc,gs_hamk,gvnlxc,1.0_dp,mpi_enreg,nband_k, &
             & dtset%prtvol,sij_opt,tim_getghc,0)
 else
    !FB: already_transposed put to false since we have not done any all_to_all before contrarily to chebfi
   call prep_getghc(cg,gs_hamk,gvnlxc,ghc,gsc,1.0_dp,nband_k,mpi_enreg,dtset%prtvol, &
                  & sij_opt,cpopt,cwaveprj,already_transposed=.false.)
 end if

 ! update eigenvalues and residuals
 ABI_MALLOC(resids, (nband))
 ABI_MALLOC(residvec, (2, npw*nspinor))
 ABI_MALLOC(nline_bands, (nband))
 do iband=1, nband
   shift = npw*nspinor*(iband-1)
   call dotprod_g(eig(iband),dprod_i,gs_hamk%istwf_k,npw*nspinor,1,ghc(:, shift+1:shift+npw*nspinor),&
                & cg(:, shift+1:shift+npw*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   if(paw) then
     call dotprod_g(dprod_r,dprod_i,gs_hamk%istwf_k,npw*nspinor,1,gsc(:, shift+1:shift+npw*nspinor),&
                  & cg(:, shift+1:shift+npw*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     eig(iband) = eig(iband)/dprod_r
   end if

   if(paw) then
     residvec = ghc(:, shift+1 : shift+npw*nspinor) - eig(iband)*gsc(:, shift+1 : shift+npw*nspinor)
   else
     residvec = ghc(:, shift+1 : shift+npw*nspinor) - eig(iband)*cg(:, shift+1 : shift+npw*nspinor)
   end if
   resids(iband) = SUM(residvec**2)
 end do
 call xmpi_sum(resids,mpi_enreg%comm_fft,ierr)

 print*, "TEST - eigenvalues", eig(:)
 print*, "TEST - resids", resids(:)

 ABI_FREE(resids)
 ABI_FREE(residvec)
 ABI_FREE(nline_bands)

 !** Also apply S^-1 in PAW case
 if(paw) then 
    call apply_invovl(gs_hamk, ghc, gsm1hc, cwaveprj, npw, nband, mpi_enreg, nspinor)
    exp_operator => gsm1hc
 else
    exp_operator => ghc
 end if

 write(97,*) cg(1,:)
 write(98,*) cg(2,:)
 !** Now update the cg 
 cg(1,:) = cg(1,:) + exp_operator(2,:)*dtset%dtele !real part
 cg(2,:) = cg(2,:) - exp_operator(1,:)*dtset%dtele !imaginary part
 ! Multiplication by a fixed phase - for tests
 !do icg = 1, size(cg(1,:))
 !   cg_old(1) = cg(1,icg)
 !   cg_old(2) = cg(2,icg)
 !   cg(1,icg) = cg_old(1)*cos(dtset%dtele/5.0_dp)-cg_old(2)*sin(dtset%dtele/5.0_dp) !real part
 !   cg(2,icg) = cg_old(1)*sin(dtset%dtele/5.0_dp)+cg_old(2)*cos(dtset%dtele/5.0_dp) !imaginary part
 !end do

 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(gvnlxc)
 ABI_FREE(gsm1hc)

 if(paw) call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 end subroutine rttddft_update_cg_k

!!****f* m_rttddft/rttddft_propagate_nuc
!!
!! NAME
!!  rttddft_propagate_nuc
!!
!! FUNCTION
!!  Main subroutine to propagate nuclei using Ehrenfest dynamics
!!
!! INPUTS
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_propagate_nuc(tdks, dtset, istep, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(in)    :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 !arrays
 
! ***********************************************************************

 write(msg,'(2a,i5,a)') ch10,'--- Iteration',istep,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 ! FB: Should we do this? 
 ! Eventually symmetrize atomic coordinates over space group elements:
 call symmetrize_xred(dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,tdks%xred,indsym=tdks%indsym)

 end subroutine rttddft_propagate_nuc

end module m_rttddft_propagate
!!***
