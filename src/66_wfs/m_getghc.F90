!!****m* ABINIT/m_getghc
!! NAME
!!  m_getghc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space;
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2025 ABINIT group (DCA, XG, GMR, LSI, MT, JB, JWZ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

module m_getghc

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_abi_linalg

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_fstrings,    only : sjoin, itoa
 use m_cgtools,     only : cg_copy_spin, cg_put_spin
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_getdim, pawcprj_copy
 use m_bandfft_kpt, only : bandfft_kpt, bandfft_kpt_get_ikpt
 use m_hamiltonian, only : gs_hamiltonian_type, KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME
 use m_fock,        only : fock_common_type, fock_get_getghc_call
 use m_fock_getghc, only : fock_getghc, fock_ACE_getghc
 use m_nonlop,      only : nonlop
 use m_gemm_nonlop_projectors, only : gemm_nonlop_use_gemm
 use m_fft,         only : fourwf
 use m_ompgpu_fourwf,      only : ompgpu_fourwf_work_mem
 use m_gemm_nonlop,        only : gemm_nonlop_ompgpu_work_mem

#if defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif
#ifdef HAVE_FFTW3_THREADS
 use m_fftw3,       only : fftw3_spawn_threads_here, fftw3_use_lib_threads
#endif
#if defined HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif
#if defined HAVE_YAKL
 use gator_mod
#endif
#ifdef HAVE_KOKKOS
 use m_manage_kokkos, only : assemble_energy_contribution_kokkos
#endif

 implicit none

 private
!!***

 public :: getghc      ! Compute <G|H|C> for input vector |C> expressed in reciprocal space
 public :: getgsc      ! Compute <G|S|C> for all input vectors |Cnk> at a given k-point
 public :: getghc_mGGA
 public :: multithreaded_getghc
 public :: getghc_nucdip ! compute <G|H_nucdip|C> for input vector |C> expressed in recip space
 public :: getghc_ompgpu_work_mem ! assess GPU memory requirements for running getghc with OpenMP GPU
!!***

contains
!!***

!!****f* ABINIT/getghc_ompgpu_work_mem
!! NAME
!! getghc_ompgpu_work_mem
!!
!! FUNCTION
!! Returns work memory requirement for getghc_ompgpu
!!
!! INPUTS
!! gs_ham <type(gs_hamiltonian_type)>=contains dimensions of FFT domain
!! ndat=size of batch for fourwf and nonlop processing
!!
!! OUTPUT
!! req_mem=amount in bytes of required memory for getghc_ompgpu

function getghc_ompgpu_work_mem(gs_ham, ndat) result(req_mem)

 type(gs_hamiltonian_type),intent(in),target :: gs_ham
 integer, intent(in) :: ndat
 integer(kind=c_size_t) :: req_mem, ghc_mem, nonlop_mem

 ! getghc use a GPU work buffer only when using fourwf
 ! Therefore, max GPU memory required by getghc is either:
 !   - the sum of getghc and fourwf work buffers memory requirements
 !   - the amount of memory required by gemm_nonlop_ompgpu work buffers
 ghc_mem = 0
 ghc_mem = int(2, c_size_t) * dp * gs_ham%n4 * gs_ham%n5 * gs_ham%n6 * ndat
 !ghc_mem = ghc_mem + ompgpu_fourwf_work_mem(gs_ham%ngfft, ndat)

 nonlop_mem = gemm_nonlop_ompgpu_work_mem(gs_ham%istwf_k, ndat, 0, gs_ham%npw_fft_k,&
 &               gs_ham%indlmn, gs_ham%nattyp, gs_ham%ntypat, gs_ham%lmnmax, 2, 11)

 req_mem = MAX(ghc_mem, nonlop_mem)

end function getghc_ompgpu_work_mem
!!***

!!****f* ABINIT/getghc
!! NAME
!! getghc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space;
!! Result is put in array ghc.
!! <G|Vnonlocal + VfockACE|C> is also returned in gvnlxc if either NLoc NCPP or FockACE.
!! If required, <G|S|C> is returned in gsc (S=overlap - PAW only)
!! Note that left and right k points can be different, i.e. ghc=<k^prime+G|H|C_k>.
!!
!! INPUTS
!! cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!       (same meaning as in nonlop.F90 routine)
!!       if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!       if cpopt= 0, <p_lmn|in> are computed here and saved
!!       if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!       if cpopt= 2  <p_lmn|in> are already in memory;
!!       if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!       if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!! cwavef(2,npw*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! cwavef_r(2,n4,n5,n6,nspinor) = wave function in real space
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFT to do in parallel
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+Fock+kinetic only (added to the existing Hamiltonian)
!!            3: local + kinetic only (added to the existing Hamiltonian)
!! ===== Optional inputs =====
!!   [kg_fft_k(3,:)]=optional, (k+G) vector coordinates to be used for the FFT transformation
!!                   instead of the one contained in gs_ham datastructure.
!!                   Typically used for real WF (in parallel) which are FFT-transformed 2 by 2.
!!   [kg_fft_kp(3,:)]=optional, (k^prime+G) vector coordinates to be used for the FFT transformation
!!   [select_k]=optional, option governing the choice of k points to be used.
!!             gs_ham datastructure contains quantities needed to apply Hamiltonian
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|H|k>       is applied [default]
!!             if select_k=2, <k|H|k^prime>       is applied
!!             if select_k=3, <k|H|k>             is applied
!!             if select_k=4, <k^prime|H|k^prime> is applied
!!
!! OUTPUT
!!  ghc(2,npw*my_nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                          or <G|H-lambda.S|C> (if sij_opt=-1)
!!  gvnlxc(2,npw*my_nspinor*ndat)=matrix elements <G|Vnonlocal+VFockACE|C> (if sij_opt>=0)
!!                                            or <G|Vnonlocal+VFockACE-lambda.S|C> (if sij_opt=-1)
!!      include Vnonlocal if NCPP and non-local Fock if associated(gs_ham%fockcommon)
!!  if (sij_opt=1)
!!    gsc(2,npw*my_nspinor*ndat)=matrix elements <G|S|C> (S=overlap).
!!
!! SIDE EFFECTS
!!  cwaveprj(natom,my_nspinor*(1+cpopt)*ndat)= wave function projected on nl projectors (PAW only)
!!
!! SOURCE

subroutine getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
                  prtvol,sij_opt,tim_getghc,type_calc,&
                  kg_fft_k,kg_fft_kp,select_k,cwavef_r,filter_dilatmx_loc) ! optional arguments

!Arguments ------------------------------------
!scalars
 logical,intent(in),optional :: filter_dilatmx_loc
 integer,intent(in) :: cpopt,ndat, prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 integer,intent(in),optional :: select_k
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
 real(dp),intent(out),target :: gsc(:,:)
 real(dp),intent(inout), target :: cwavef(:,:)
 real(dp),optional,intent(inout) :: cwavef_r(:,:,:,:,:)
 real(dp),intent(out), target :: ghc(:,:)
 real(dp),intent(out),target :: gvnlxc(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
 !MG: Passing these arrays assumed-shape has the drawback that client code is
 !forced to use vec(2, npw*ndat) instead of the more user-friendly vec(2,npw,ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=114, tim_fourwf=1
 integer :: choice,cplex,cpopt_here,i1,i2,i3,idat,idir,ierr,i0
 integer :: ig,igspinor,istwf_k_,ii,iispinor,ikpt_this_proc,ipw,ispinor,my_nspinor
 integer :: n4,n5,n6,ndat_,nnlout,npw_fft,npw_k1,npw_k2,nspinortot,option_fft
 integer :: paw_opt,select_k_,shift1,shift2,signs,tim_nonlop
 logical(kind=c_bool) :: k1_eq_k2
 logical :: double_rfft_trick,have_to_reequilibrate,has_fock,local_gvnlxc
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc,use_cwavef_r, filter_dilatmx_loc_
 real(dp) :: ghcim,ghcre,weight
#ifdef HAVE_OPENMP_OFFLOAD
 complex(dpc), parameter :: cminusone  = (-1._dp,0._dp)
#endif
 character(len=500) :: msg
!arrays
 integer,  pointer                :: gbound_k1(:,:), gbound_k2(:,:)
 integer,  pointer                :: kg_k1(:,:), kg_k2(:,:)
 integer,  ABI_CONTIGUOUS pointer :: indices_pw_fft(:), kg_k_fft(:,:)
 integer,  ABI_CONTIGUOUS pointer :: recvcount_fft(:), recvdisp_fft(:)
 integer,  ABI_CONTIGUOUS pointer :: sendcount_fft(:), senddisp_fft(:)
 integer,  allocatable:: dimcprj(:)
 real(dp)                         :: enlout(ndat), lambda_ndat(ndat), tsec(2)
 real(dp), target                 :: nonlop_dum(1,1)
 real(dp), allocatable            :: buff_wf(:,:), cwavef1(:,:), cwavef2(:,:), cwavef_fft(:,:), cwavef_fft_tr(:,:)
 real(dp), allocatable            :: cwavef_spin(:,:), gvnlxc_spin(:,:)
 real(dp), allocatable            :: ghc1(:,:), ghc2(:,:), ghc3(:,:),  ghc4(:,:)
 real(dp), allocatable            :: ghc_mGGA(:,:), ghc_vectornd(:,:)

#if defined HAVE_GPU && defined HAVE_YAKL
 real(c_double), ABI_CONTIGUOUS pointer :: gvnlc(:,:)
 real(c_double), ABI_CONTIGUOUS pointer :: gvnlxc_(:,:)
#else
 real(dp), allocatable            :: gvnlc(:,:)
 real(dp), pointer                :: gvnlxc_(:,:)
#endif
 real(dp), allocatable            :: vlocal_tmp(:,:,:), work(:,:,:,:)
 real(dp), pointer                :: kinpw_k1(:), kinpw_k2(:), kpt_k1(:), kpt_k2(:)
 real(dp), pointer                :: gsc_ptr(:,:)
 type(fock_common_type),pointer :: fock
 type(pawcprj_type),pointer :: cwaveprj_fock(:,:),cwaveprj_idat(:,:),cwaveprj_nonlop(:,:)
 logical :: transfer_ghc,transfer_gsc,transfer_cwavef,transfer_gvnlxc
 real(c_double), parameter        :: hugevalue = huge(zero)*1.d-11
! *********************************************************************

 DBG_ENTER("COLL")

 !Keep track of total time spent in getghc:
 call timab(350+tim_getghc,1,tsec)
 ABI_NVTX_START_RANGE(NVTX_GETGHC)

 !Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' getghc : enter, debugging '
   call wrtout(std_out,msg,'PERS')
 end if

 !Select k-dependent objects according to select_k input parameter
 select_k_=KPRIME_H_K;if (present(select_k)) select_k_=select_k
 select case (select_k_)
 case (KPRIME_H_K)
   ! <k^prime|H|k>
   npw_k1    =  gs_ham%npw_fft_k ; npw_k2    =  gs_ham%npw_fft_kp
   kpt_k1    => gs_ham%kpt_k     ; kpt_k2    => gs_ham%kpt_kp
   kg_k1     => gs_ham%kg_k      ; kg_k2     => gs_ham%kg_kp
   gbound_k1 => gs_ham%gbound_k  ; gbound_k2 => gs_ham%gbound_kp
   kinpw_k1  => gs_ham%kinpw_k   ; kinpw_k2  => gs_ham%kinpw_kp
 case (K_H_KPRIME)
   ! <k|H|k^prime>
   npw_k1    =  gs_ham%npw_fft_kp; npw_k2    =  gs_ham%npw_fft_k
   kpt_k1    => gs_ham%kpt_kp    ; kpt_k2    => gs_ham%kpt_k
   kg_k1     => gs_ham%kg_kp     ; kg_k2     => gs_ham%kg_k
   gbound_k1 => gs_ham%gbound_kp ; gbound_k2 => gs_ham%gbound_k
   kinpw_k1  => gs_ham%kinpw_kp  ; kinpw_k2  => gs_ham%kinpw_k
 case (K_H_K)
   ! <k|H|k>
   npw_k1    =  gs_ham%npw_fft_k ; npw_k2    =  gs_ham%npw_fft_k
   kpt_k1    => gs_ham%kpt_k     ; kpt_k2    => gs_ham%kpt_k
   kg_k1     => gs_ham%kg_k      ; kg_k2     => gs_ham%kg_k
   gbound_k1 => gs_ham%gbound_k  ; gbound_k2 => gs_ham%gbound_k
   kinpw_k1  => gs_ham%kinpw_k   ; kinpw_k2  => gs_ham%kinpw_k
 case (KPRIME_H_KPRIME)
   ! <k^prime|H|k^prime>
   npw_k1    =  gs_ham%npw_fft_kp; npw_k2    =  gs_ham%npw_fft_kp
   kpt_k1    => gs_ham%kpt_kp    ; kpt_k2    => gs_ham%kpt_kp
   kg_k1     => gs_ham%kg_kp     ; kg_k2     => gs_ham%kg_kp
   gbound_k1 => gs_ham%gbound_kp ; gbound_k2 => gs_ham%gbound_kp
   kinpw_k1  => gs_ham%kinpw_kp  ; kinpw_k2  => gs_ham%kinpw_kp
 case default
   ABI_ERROR(sjoin("Invalid select_k: ", itoa(select_k_)))
 end select

 k1_eq_k2=(all(abs(kpt_k1(:)-kpt_k2(:))<tol8))

 ! Check sizes
 my_nspinor=max(1,gs_ham%nspinor/mpi_enreg%nproc_spinor)
 ABI_CHECK_IGEQ(size(cwavef), 2*npw_k1*my_nspinor*ndat, 'wrong size for cwavef!')
 ABI_CHECK_IGEQ(size(ghc), 2*npw_k2*my_nspinor*ndat, 'wrong size for ghc!')
 if (sij_opt==1) then
   ABI_CHECK_IGEQ(size(gsc), 2*npw_k2*my_nspinor*ndat, 'wrong size for gsc!')
 end if
 if (gs_ham%usepaw==1.and.cpopt>=0) then
   ABI_CHECK_IGEQ(size(cwaveprj), gs_ham%natom*my_nspinor*ndat, 'wrong size for cwaveprj!')
 end if
 if (any(type_calc == [0, 2, 3])) then
   local_gvnlxc = size(gvnlxc)<=1
   if (local_gvnlxc) then
     if(gs_ham%gpu_option==ABI_GPU_KOKKOS) then
#if defined HAVE_GPU && defined HAVE_YAKL
       ABI_MALLOC_MANAGED(gvnlxc_, (/2,npw_k2*my_nspinor*ndat/))
#endif
     else
       ABI_MALLOC(gvnlxc_,(2,npw_k2*my_nspinor*ndat))
     end if
   else
     gvnlxc_ => gvnlxc
   end if
   ABI_CHECK_IGEQ(size(gvnlxc_), 2*npw_k2*my_nspinor*ndat, 'wrong size for gvnlxc!')
 end if

 use_cwavef_r=present(cwavef_r)
 n4=gs_ham%n4 ; n5=gs_ham%n5 ; n6=gs_ham%n6
 nspinortot=gs_ham%nspinor

 if (use_cwavef_r) then
   ABI_CHECK_IEQ(size(cwavef_r,1), 2, 'wrong size for cwavef_r (dimension 1)')
   ABI_CHECK_IEQ(size(cwavef_r,2), n4, 'wrong size for cwavef_r (dimension 2)')
   ABI_CHECK_IEQ(size(cwavef_r,3), n5, 'wrong size for cwavef_r (dimension 3)')
   ABI_CHECK_IEQ(size(cwavef_r,4), n6, 'wrong size for cwavef_r (dimension 4)')
   ABI_CHECK_IEQ(size(cwavef_r,5), nspinortot, 'wrong size for cwavef_r (dimension 5)')
 end if

 !Eventually overwrite plane waves data for FFT
 if (present(kg_fft_k)) then
   kg_k1 => kg_fft_k ; kg_k2 => kg_fft_k
   npw_k1=size(kg_k1,2) ; npw_k2=size(kg_k2,2)
 end if
 if (present(kg_fft_kp)) then
   kg_k2 => kg_fft_kp ; npw_k2=size(kg_k2,2)
 end if

 !paral_kgb constraint
 if (mpi_enreg%paral_kgb==1.and.(.not.k1_eq_k2)) then
   ABI_BUG('paral_kgb=1 not allowed for k/=k_^prime!')
 end if

 !Do we add Fock exchange term ?
 has_fock = associated(gs_ham%fockcommon)
 if (has_fock) fock => gs_ham%fockcommon

 !Parallelization over spinors management
 if (mpi_enreg%paral_spinor==0) then
   shift1=npw_k1;shift2=npw_k2
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(nspinortot==2)
 else
   shift1=0;shift2=0
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

 filter_dilatmx_loc_ = .true.; if ( present(filter_dilatmx_loc) ) filter_dilatmx_loc_ = filter_dilatmx_loc
 transfer_ghc = .false.; transfer_gsc = .false.; transfer_gvnlxc = .false.; transfer_cwavef = .false.
 if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
 transfer_ghc =  .not. xomp_target_is_present(c_loc(ghc))
 transfer_gsc =  .not. xomp_target_is_present(c_loc(gsc))
 transfer_gvnlxc =  .not. xomp_target_is_present(c_loc(gvnlxc_))
 transfer_cwavef =  .not. xomp_target_is_present(c_loc(cwavef))

 !$OMP TARGET ENTER DATA MAP(alloc:ghc)     IF(transfer_ghc)
 !$OMP TARGET ENTER DATA MAP(alloc:gsc)     IF(transfer_gsc)
 !$OMP TARGET ENTER DATA MAP(alloc:gvnlxc_) IF(transfer_gvnlxc)
 !$OMP TARGET ENTER DATA MAP(to:cwavef)     IF(transfer_cwavef)
 if (type_calc == 2) then
   !$OMP TARGET UPDATE TO(ghc)     IF(transfer_ghc)
   !$OMP TARGET UPDATE TO(gsc)     IF(transfer_gsc .and. gs_ham%usepaw==1)
 end if
#endif
 end if

 !============================================================
 ! Application of the local potential
 !============================================================

 ABI_NVTX_START_RANGE(NVTX_GETGHC_LOCPOT)

 if (any(type_calc == [0, 1, 3])) then

   ! Need a Vlocal
   ABI_CHECK(associated(gs_ham%vlocal), "We need vlocal in gs_ham!")

   ! fourwf can only process with one value of istwf_k
   if (gs_ham%use_gbt == 0) then
     ABI_CHECK(k1_eq_k2, 'vlocal (fourwf) cannot be computed with k/=k^prime!')
   end if

   ! Eventually adjust load balancing for FFT (by changing FFT distrib)
   ! Not that have_to_reequilibrate can be true only if mpi_enreg%nproc_fft>1
   have_to_reequilibrate=.false.
   if (mpi_enreg%paral_kgb==1) then
     ikpt_this_proc=bandfft_kpt_get_ikpt()
     have_to_reequilibrate=bandfft_kpt(ikpt_this_proc)%have_to_reequilibrate
   end if
   ndat_             = ndat
   istwf_k_          = gs_ham%istwf_k
   double_rfft_trick = istwf_k_==2.and.ndat>1.and.mpi_enreg%paral_kgb==1.and.gs_ham%gpu_option==ABI_GPU_DISABLED
   ! LB-08-2024 : double_rfft_trick works only for paral_kgb=1, but I don't know why...
   ! Note that the trick can be activated only if nspinortot=1 (if =2 then istwf_k=1), so gs_ham%nvloc=1 too
   if (double_rfft_trick) then
     ABI_CHECK(mpi_enreg%nproc_fft == 1, 'double_rfft_trick inside getghc is not implemented for npfft>1')
     istwf_k_ = 1
     ndat_ = ndat / 2
     if (modulo(ndat,2)/=0) ndat_=ndat_+1
     npw_fft=2*npw_k1
     i0=1
     if (mpi_enreg%me_g0_fft==1) then ! Do not include G=(0,0,0) twice
       npw_fft=npw_fft-1
       i0=2
     end if
     ABI_MALLOC(kg_k_fft,(3,npw_fft))
     ABI_MALLOC(cwavef_fft,(2,npw_fft*ndat_))
     kg_k_fft(:,1:npw_k1) = kg_k1(:,1:npw_k1)
     kg_k_fft(:,npw_k1+1:npw_fft) = -kg_k1(:,i0:npw_k1)
     call cwavef_double_rfft_trick_pack(cwavef,cwavef_fft,mpi_enreg%me_g0_fft,ndat,npw_k1)
   end if

   if (have_to_reequilibrate.and.double_rfft_trick) then
     ABI_ERROR("In getghc: have_to_reequilibrate cannot be activated with double_rfft_trick")
   end if

   if (have_to_reequilibrate) then
     ! Note: for this case we have ndat_=ndat
     npw_fft =  bandfft_kpt(ikpt_this_proc)%npw_fft
     sendcount_fft  => bandfft_kpt(ikpt_this_proc)%sendcount_fft(:)
     recvcount_fft  => bandfft_kpt(ikpt_this_proc)%recvcount_fft(:)
     senddisp_fft   => bandfft_kpt(ikpt_this_proc)%senddisp_fft(:)
     recvdisp_fft   => bandfft_kpt(ikpt_this_proc)%recvdisp_fft(:)
     indices_pw_fft => bandfft_kpt(ikpt_this_proc)%indices_pw_fft(:)
     kg_k_fft       => bandfft_kpt(ikpt_this_proc)%kg_k_fft(:,:)
     ABI_MALLOC(buff_wf,(2,npw_k1*ndat) )
     ABI_MALLOC(cwavef_fft,(2,npw_fft*ndat) )
     if(ndat>1) then
       ABI_MALLOC(cwavef_fft_tr, (2,npw_fft*ndat))
     end if
     do idat=1, ndat
       do ipw = 1 ,npw_k1
         buff_wf(1:2, idat + ndat*(indices_pw_fft(ipw)-1) ) = cwavef(1:2,ipw + npw_k1*(idat-1))
       end do
     end do
     if(ndat > 1) then
       call xmpi_alltoallv(buff_wf,2*ndat*sendcount_fft,2*ndat*senddisp_fft,  &
        cwavef_fft_tr,2*ndat*recvcount_fft, 2*ndat*recvdisp_fft, mpi_enreg%comm_fft,ierr)
       ! We need to transpose data
       do idat=1,ndat
         do ipw = 1 ,npw_fft
           cwavef_fft(1:2,  ipw + npw_fft*(idat-1)) = cwavef_fft_tr(1:2,  idat + ndat*(ipw-1))
         end do
       end do
     else
       call xmpi_alltoallv(buff_wf,2*sendcount_fft,2*senddisp_fft,  &
        cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, mpi_enreg%comm_fft,ierr)
     end if
   end if

   ! Apply the local potential to the wavefunction
   ! Start from wavefunction in reciprocal space cwavef
   ! End with function ghc in reciprocal space also.
   ABI_MALLOC(work,(2,gs_ham%n4,gs_ham%n5,gs_ham%n6*ndat))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:work) IF(gs_ham%gpu_option==ABI_GPU_OPENMP)
#endif
   weight=one
   if (.not.use_cwavef_r) then
     option_fft=2
     if (nspinortot==2) then
       ! Note: for this case we have ndat_=ndat
       ABI_MALLOC(cwavef1,(2,npw_k1*ndat))
       ABI_MALLOC(cwavef2,(2,npw_k1*ndat))
       if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET ENTER DATA MAP(alloc:cwavef1,cwavef2)
         !$OMP TARGET TEAMS DISTRIBUTE MAP(to:cwavef1,cwavef2,cwavef)
         do idat=1,ndat
           !$OMP PARALLEL DO COLLAPSE(2)
           do ipw=1,npw_k1
             do ig=1,2
               cwavef1(ig,ipw+(idat-1)*npw_k1)=cwavef(ig,ipw+(idat-1)*my_nspinor*npw_k1)
               cwavef2(ig,ipw+(idat-1)*npw_k1)=cwavef(ig,ipw+(idat-1)*my_nspinor*npw_k1+shift1)
             end do
           end do
         end do
#endif
       else
         do idat=1,ndat
           do ipw=1,npw_k1
             cwavef1(1:2,ipw+(idat-1)*npw_k1)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k1)
             cwavef2(1:2,ipw+(idat-1)*npw_k1)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k1+shift1)
           end do
         end do
       end if
     end if
   else
     option_fft=3
     if (nspinortot==2) then
       ABI_MALLOC(cwavef1,(0,0))
       ABI_MALLOC(cwavef2,(0,0))
     end if
   end if

   if (gs_ham%nvloc==1) then
     !  Treat scalar local potentials

     if (nspinortot==1) then

       if (use_cwavef_r) then
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               work(1,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(1,i1,i2,i3,1)
               work(2,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(2,i1,i2,i3,1)
             end do
           end do
         end do
       end if
       if (have_to_reequilibrate.or.double_rfft_trick) then
         call fourwf(1,gs_ham%vlocal,cwavef_fft,cwavef_fft,work,gbound_k1,gbound_k2,&
          istwf_k_,kg_k_fft,kg_k_fft,gs_ham%mgfft,mpi_enreg,ndat_,gs_ham%ngfft,&
          npw_fft,npw_fft,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
          weight,weight,gpu_option=gs_ham%gpu_option)
       else
         call fourwf(1,gs_ham%vlocal,cwavef,ghc,work,gbound_k1,gbound_k2,&
          istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat_,gs_ham%ngfft,&
          npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
          weight,weight,gpu_option=gs_ham%gpu_option)
       end if

     else
       ! nspinortot==2

       ! Note: for this case we have ndat_=ndat
       if (nspinor1TreatedByThisProc) then
         if (use_cwavef_r) then
           do i3=1,n6
             do i2=1,n5
               do i1=1,n4
                 work(1,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(1,i1,i2,i3,1)
                 work(2,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(2,i1,i2,i3,1)
               end do
             end do
           end do
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET UPDATE TO(work) IF(gs_ham%gpu_option==ABI_GPU_OPENMP)
#endif
         end if
         ABI_MALLOC(ghc1,(2,npw_k2*ndat))
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET ENTER DATA MAP(alloc:ghc1) IF(gs_ham%gpu_option==ABI_GPU_OPENMP)
#endif
         call fourwf(1,gs_ham%vlocal,cwavef1,ghc1,work,gbound_k1,gbound_k2,&
          istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
          npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
          weight,weight,gpu_option=gs_ham%gpu_option)
         if(gs_ham%gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc1)
           do idat=1,ndat
             !$OMP PARALLEL DO COLLAPSE(2)
             do ipw =1, npw_k2
               do ig =1,2
                 ghc(ig,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(ig,ipw+(idat-1)*npw_k2)
               end do
             end do
           end do
#endif
         else
           do idat=1,ndat
             do ipw =1, npw_k2
               ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(1:2,ipw+(idat-1)*npw_k2)
             end do
           end do
         end if
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET EXIT DATA MAP(delete:ghc1) IF(gs_ham%gpu_option==ABI_GPU_OPENMP)
#endif
         ABI_FREE(ghc1)
       end if ! spin 1 treated by this proc

       if (nspinor2TreatedByThisProc) then
         if (use_cwavef_r) then
           do i3=1,n6
             do i2=1,n5
               do i1=1,n4
                 work(1,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(1,i1,i2,i3,2)
                 work(2,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(2,i1,i2,i3,2)
               end do
             end do
           end do
         end if
         ABI_MALLOC(ghc2,(2,npw_k2*ndat))
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET ENTER DATA MAP(alloc:ghc2) IF(gs_ham%gpu_option==ABI_GPU_OPENMP)
#endif
         call fourwf(1,gs_ham%vlocal,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
           istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
           npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
           gpu_option=gs_ham%gpu_option)

         if(gs_ham%gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc2)
           do idat=1,ndat
             !$OMP PARALLEL DO COLLAPSE(2)
             do ipw=1,npw_k2
               do ig =1,2
                 ghc(ig,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc2(ig,ipw+(idat-1)*npw_k2)
               end do
             end do
           end do
#endif
         else
           do idat=1,ndat
             do ipw=1,npw_k2
               ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc2(1:2,ipw+(idat-1)*npw_k2)
             end do
           end do
         end if
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET EXIT DATA MAP(delete:ghc2) IF(gs_ham%gpu_option==ABI_GPU_OPENMP)
#endif
         ABI_FREE(ghc2)
       end if ! spin 2 treated by this proc

     end if ! nspinortot

   else if (gs_ham%nvloc==4) then
     ! Treat non-collinear local potentials
     ! Note: for this case we have ndat_=ndat
     ABI_MALLOC(ghc1,(2,npw_k2*ndat))
     ABI_MALLOC(ghc2,(2,npw_k2*ndat))
     ABI_MALLOC(ghc3,(2,npw_k2*ndat))
     ABI_MALLOC(ghc4,(2,npw_k2*ndat))
     if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET ENTER DATA MAP(alloc:ghc1,ghc2,ghc3,ghc4)
       call gpu_set_to_zero(ghc1, int(2,c_size_t)*npw_k2*ndat)
       call gpu_set_to_zero(ghc2, int(2,c_size_t)*npw_k2*ndat)
       call gpu_set_to_zero(ghc3, int(2,c_size_t)*npw_k2*ndat)
       call gpu_set_to_zero(ghc4, int(2,c_size_t)*npw_k2*ndat)
#endif
     else
       ghc1(:,:)=zero; ghc2(:,:)=zero; ghc3(:,:)=zero ;  ghc4(:,:)=zero
     end if
     if (use_cwavef_r) then
       ABI_MALLOC(vlocal_tmp,(0,0,0))
     else
       ABI_MALLOC(vlocal_tmp,(gs_ham%n4,gs_ham%n5,gs_ham%n6))
     end if
     ! ghc1=v11*phi1
     if (nspinor1TreatedByThisProc) then
       if (use_cwavef_r) then
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               work(1,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(1,i1,i2,i3,1)
               work(2,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(2,i1,i2,i3,1)
             end do
           end do
         end do
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET UPDATE TO(work) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
       else
         ! LB,07/22:
         ! Weird segmentation fault encountered here if called with multithreaded_getghc for big systems.
         ! Using an explicit loop instead of fortran syntax seems to solve the problem, I don't understand why...
         !vlocal_tmp(:,:,:)=gs_ham%vlocal(:,:,:,1)
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               vlocal_tmp(i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,1)
             end do
           end do
         end do
       end if
       call fourwf(1,vlocal_tmp,cwavef1,ghc1,work,gbound_k1,gbound_k2,&
         istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
         gpu_option=gs_ham%gpu_option)
     end if
     ! ghc2=v22*phi2
     if (nspinor2TreatedByThisProc) then
       if (use_cwavef_r) then
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               work(1,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,2)*cwavef_r(1,i1,i2,i3,2)
               work(2,i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,2)*cwavef_r(2,i1,i2,i3,2)
             end do
           end do
         end do
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET UPDATE TO(work) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
       else
         ! LB,07/22:
         ! Weird segmentation fault encountered here if called with multithreaded_getghc for big systems.
         ! Using an explicit loop instead of fortran syntax seems to solve the problem, I don't understand why...
         !vlocal_tmp(:,:,:)=gs_ham%vlocal(:,:,:,2)
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               vlocal_tmp(i1,i2,i3) = gs_ham%vlocal(i1,i2,i3,2)
             end do
           end do
         end do
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET UPDATE TO(work) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
       end if
       call fourwf(1,vlocal_tmp,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
         istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
         gpu_option=gs_ham%gpu_option)
     end if
     ABI_FREE(vlocal_tmp)
     cplex=2
     if (use_cwavef_r) then
       ABI_MALLOC(vlocal_tmp,(0,0,0))
     else
       ABI_MALLOC(vlocal_tmp,(cplex*gs_ham%n4,gs_ham%n5,gs_ham%n6))
     end if
     ! ghc3=(re(v12)-im(v12))*phi1
     if (nspinor1TreatedByThisProc) then
       if (use_cwavef_r) then
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               work(1,i1,i2,i3)= gs_ham%vlocal(i1,i2,i3,3)*cwavef_r(1,i1,i2,i3,1)+gs_ham%vlocal(i1,i2,i3,4)*cwavef_r(2,i1,i2,i3,1)
               work(2,i1,i2,i3)=-gs_ham%vlocal(i1,i2,i3,4)*cwavef_r(1,i1,i2,i3,1)+gs_ham%vlocal(i1,i2,i3,3)*cwavef_r(2,i1,i2,i3,1)
             end do
           end do
         end do
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET UPDATE TO(work) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
       else
         do i3=1,gs_ham%n6
           do i2=1,gs_ham%n5
             do i1=1,gs_ham%n4
               vlocal_tmp(2*i1-1,i2,i3)= gs_ham%vlocal(i1,i2,i3,3)
               vlocal_tmp(2*i1  ,i2,i3)=-gs_ham%vlocal(i1,i2,i3,4)
             end do
           end do
         end do
       end if
       call fourwf(cplex,vlocal_tmp,cwavef1,ghc3,work,gbound_k1,gbound_k2,&
         istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
         gpu_option=gs_ham%gpu_option)
     end if
     ! ghc4=(re(v12)+im(v12))*phi2
     if (nspinor2TreatedByThisProc) then
       if (use_cwavef_r) then
         do i3=1,n6
           do i2=1,n5
             do i1=1,n4
               work(1,i1,i2,i3)= gs_ham%vlocal(i1,i2,i3,3)*cwavef_r(1,i1,i2,i3,2)-gs_ham%vlocal(i1,i2,i3,4)*cwavef_r(2,i1,i2,i3,2)
               work(2,i1,i2,i3)= gs_ham%vlocal(i1,i2,i3,4)*cwavef_r(1,i1,i2,i3,2)+gs_ham%vlocal(i1,i2,i3,3)*cwavef_r(2,i1,i2,i3,2)
             end do
           end do
         end do
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET UPDATE TO(work) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
       else
         do i3=1,gs_ham%n6
           do i2=1,gs_ham%n5
             do i1=1,gs_ham%n4
               vlocal_tmp(2*i1-1,i2,i3)= gs_ham%vlocal(i1,i2,i3,3)
               vlocal_tmp(2*i1  ,i2,i3)= gs_ham%vlocal(i1,i2,i3,4)
             end do
           end do
         end do
       end if
       call fourwf(cplex,vlocal_tmp,cwavef2,ghc4,work,gbound_k1,gbound_k2,&
         istwf_k_,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
         gpu_option=gs_ham%gpu_option)
     end if
     ABI_FREE(vlocal_tmp)
     ! Build ghc from pieces
     ! (v11,v22,Re(v12)+iIm(v12);Re(v12)-iIm(v12))(psi1;psi2): matrix product
     if (mpi_enreg%paral_spinor==0) then
       if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc1,ghc2,ghc3,ghc4)
         do idat=1,ndat
           !$OMP PARALLEL DO
           do ipw=1,npw_k2
             ghc(1,ipw+(idat-1)*my_nspinor*npw_k2)       =ghc1(1,ipw+(idat-1)*npw_k2)+ghc4(1,ipw+(idat-1)*npw_k2)
             ghc(2,ipw+(idat-1)*my_nspinor*npw_k2)       =ghc1(2,ipw+(idat-1)*npw_k2)+ghc4(2,ipw+(idat-1)*npw_k2)
             ghc(1,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1,ipw+(idat-1)*npw_k2)+ghc2(1,ipw+(idat-1)*npw_k2)
             ghc(2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(2,ipw+(idat-1)*npw_k2)+ghc2(2,ipw+(idat-1)*npw_k2)
           end do
         end do
#endif
       else
         do idat=1,ndat
           do ipw=1,npw_k2
             ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2)       =ghc1(1:2,ipw+(idat-1)*npw_k2)+ghc4(1:2,ipw+(idat-1)*npw_k2)
             ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1:2,ipw+(idat-1)*npw_k2)+ghc2(1:2,ipw+(idat-1)*npw_k2)
           end do
         end do
       end if
     else
       call xmpi_sum(ghc4,mpi_enreg%comm_spinor,ierr)
       call xmpi_sum(ghc3,mpi_enreg%comm_spinor,ierr)
       if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         if (nspinor1TreatedByThisProc) then
           !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc1,ghc4)
           do idat=1,ndat
             !$OMP PARALLEL DO
             do ipw=1,npw_k2
               ghc(1,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(1,ipw+(idat-1)*npw_k2)+ghc4(1,ipw+(idat-1)*npw_k2)
               ghc(2,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(2,ipw+(idat-1)*npw_k2)+ghc4(2,ipw+(idat-1)*npw_k2)
             end do
           end do
         else if (nspinor2TreatedByThisProc) then
           !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc2,ghc3)
           do idat=1,ndat
             !$OMP PARALLEL DO
             do ipw=1,npw_k2
               ghc(1,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1,ipw+(idat-1)*npw_k2)+ghc2(1,ipw+(idat-1)*npw_k2)
               ghc(2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(2,ipw+(idat-1)*npw_k2)+ghc2(2,ipw+(idat-1)*npw_k2)
             end do
           end do
         end if
#endif
       else
         if (nspinor1TreatedByThisProc) then
           do idat=1,ndat
             do ipw=1,npw_k2
               ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(1:2,ipw+(idat-1)*npw_k2)+ghc4(1:2,ipw+(idat-1)*npw_k2)
             end do
           end do
         else if (nspinor2TreatedByThisProc) then
           do idat=1,ndat
             do ipw=1,npw_k2
               ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1:2,ipw+(idat-1)*npw_k2)+ghc2(1:2,ipw+(idat-1)*npw_k2)
             end do
           end do
         end if
       end if
     end if
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET EXIT DATA MAP(delete:ghc1,ghc2,ghc3,ghc4) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
     ABI_FREE(ghc1)
     ABI_FREE(ghc2)
     ABI_FREE(ghc3)
     ABI_FREE(ghc4)
   end if ! nvloc

   if (nspinortot==2)  then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET EXIT DATA MAP(delete:cwavef1,cwavef2) IF (gs_ham%gpu_option == ABI_GPU_OPENMP .and. .not.use_cwavef_r)
#endif
     ABI_FREE(cwavef1)
     ABI_FREE(cwavef2)
   end if

#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:work) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
   ABI_FREE(work)

   if (double_rfft_trick) then
     call cwavef_double_rfft_trick_unpack(ghc,cwavef_fft,mpi_enreg%me_g0_fft,ndat,npw_k1)
     ABI_FREE(kg_k_fft)
     ABI_FREE(cwavef_fft)
   end if

   ! Retrieve eventually original FFT distrib
   if(have_to_reequilibrate) then
     ! Note: for this case we have ndat_=ndat
     if(ndat > 1 ) then
       do idat=1,ndat
         do ipw = 1 ,npw_fft
           cwavef_fft_tr(1:2,  idat + ndat*(ipw-1)) = cwavef_fft(1:2,  ipw + npw_fft*(idat-1))
         end do
       end do
       call xmpi_alltoallv(cwavef_fft_tr,2*ndat*recvcount_fft, 2*ndat*recvdisp_fft, &
         buff_wf,2*ndat*sendcount_fft,2*ndat*senddisp_fft, mpi_enreg%comm_fft,ierr)
     else
       call xmpi_alltoallv(cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, &
         buff_wf,2*sendcount_fft,2*senddisp_fft, mpi_enreg%comm_fft,ierr)
     end if
     do idat=1,ndat
       do ipw = 1 ,npw_k2
         ghc(1:2,ipw + npw_k2*(idat-1)) = buff_wf(1:2, idat + ndat*(indices_pw_fft(ipw)-1))
       end do
     end do
     ABI_FREE(buff_wf)
     ABI_FREE(cwavef_fft)
     if(ndat > 1) then
       ABI_FREE(cwavef_fft_tr)
     end if
   end if

   ! Add metaGGA contribution
   if (associated(gs_ham%vxctaulocal)) then
     ABI_CHECK(k1_eq_k2, 'metaGGA not allowed for k/=k_^prime!')
     ABI_CHECK_IEQ(size(gs_ham%vxctaulocal), gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*4, 'wrong sizes for vxctaulocal!')

     ABI_MALLOC(ghc_mGGA,(2,npw_k2*my_nspinor*ndat))
     call getghc_mGGA(cwavef,ghc_mGGA,gbound_k1,gs_ham%gprimd,istwf_k_,kg_k1,kpt_k1,&
       gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
       gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vxctaulocal,gs_ham%gpu_option)
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE FROM(ghc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_mGGA(1:2,1:npw_k2*my_nspinor*ndat)
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE TO(ghc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
     ABI_FREE(ghc_mGGA)
   end if

   !  Add nuclear dipole moment contribution
   if (associated(gs_ham%vectornd)) then
     ABI_CHECK(k1_eq_k2, 'nuclear dipole vector potential not allowed for k/=k_^prime!')
     if (size(gs_ham%vectornd)/=gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*3) then
       ABI_BUG('wrong sizes for vectornd in getghc!')
     end if
     ABI_MALLOC(ghc_vectornd,(2,npw_k2*my_nspinor*ndat))
     ghc_vectornd=zero
     call getghc_nucdip(cwavef,ghc_vectornd,gbound_k1,istwf_k_,kg_k1,kpt_k1,&
       gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
       gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vectornd,gs_ham%vlocal,gs_ham%zora,gs_ham%gpu_option)

#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE FROM(ghc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_vectornd(1:2,1:npw_k2*my_nspinor*ndat)
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE TO(ghc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif

     ABI_FREE(ghc_vectornd)
   end if

   ! If only local part is applied, still we have to filter the result because of dilatmx.
   ! Otherwise it is done when adding kinetic term.
   if (type_calc==1.and.filter_dilatmx_loc_) then

     if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,kinpw_k2)
       do idat=1,ndat
         !$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2)
         do ispinor=1,my_nspinor
           do ig=1,npw_k2
             igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
             if(kinpw_k2(ig)>huge(zero)*1.d-11) ghc(:,igspinor)=zero
           end do ! ig
         end do ! ispinor
       end do
#endif
     else
       do idat=1,ndat
         !$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2) IF(gemm_nonlop_use_gemm)
         do ispinor=1,my_nspinor
           do ig=1,npw_k2
             igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
             if(kinpw_k2(ig)>huge(zero)*1.d-11) ghc(:,igspinor)=zero
           end do ! ig
         end do ! ispinor
       end do
     end if
   end if

 end if ! type_calc

 ABI_NVTX_END_RANGE()

 if (any(type_calc == [0, 2, 3])) then
   !============================================================
   ! Application of the non-local potential and the Fock potential
   !============================================================

   ABI_NVTX_START_RANGE(NVTX_GETGHC_NLOCPOT)

   if (type_calc==0 .or. type_calc==2) then
     signs=2 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=1
     cpopt_here=-1;if (gs_ham%usepaw==1) cpopt_here=cpopt
     if (has_fock) then
       if (gs_ham%usepaw==1) then
         cpopt_here=max(cpopt,0)
         if (cpopt<2) then
           ABI_MALLOC(cwaveprj_fock,(gs_ham%natom,my_nspinor*ndat))
           ABI_MALLOC(dimcprj,(gs_ham%natom))
           call pawcprj_getdim(dimcprj,gs_ham%natom,gs_ham%nattyp,gs_ham%ntypat,gs_ham%typat,fock%pawtab,'O')
           call pawcprj_alloc(cwaveprj_fock,0,dimcprj)
           ABI_FREE(dimcprj)
         else
           cwaveprj_fock=>cwaveprj
         end if
         cwaveprj_nonlop=>cwaveprj_fock
       else
         cwaveprj_nonlop=>cwaveprj
         cwaveprj_fock=>cwaveprj
       end if
     else
       cwaveprj_nonlop=>cwaveprj
     end if
     paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3
     lambda_ndat = lambda

     if (gs_ham%use_gbt == 0) then
       if (gs_ham%usepaw==0) gsc_ptr => nonlop_dum
       if (gs_ham%usepaw==1) gsc_ptr => gsc

       call nonlop(choice,cpopt_here,cwaveprj_nonlop,enlout,gs_ham,idir,lambda_ndat,mpi_enreg,ndat,&
          nnlout,paw_opt,signs,gsc_ptr,tim_nonlop,cwavef,gvnlxc_,select_k=select_k_)
     else
       ! GBT case
       gs_ham%nspinor = 1

       ! Split cwavef and gvnlxc
       ABI_MALLOC(cwavef_spin, (2, npw_k1*ndat))
       ABI_MALLOC(gvnlxc_spin, (2, npw_k1*ndat))

       if (gs_ham%usepaw==0) gsc_ptr => nonlop_dum
       if (gs_ham%usepaw==1) gsc_ptr => gsc

       ! Apply Vnl{k-q/2} to u^up
       call cg_copy_spin(1, npw_k1, nspinortot, ndat, cwavef, cwavef_spin)
       call nonlop(choice, cpopt_here, cwaveprj, enlout, gs_ham, idir, lambda_ndat, mpi_enreg, ndat, &
                   nnlout, paw_opt, signs, gsc_ptr, tim_nonlop, cwavef_spin, gvnlxc_spin, select_k=K_H_K)
       ! Insert results in the right position.
       call cg_put_spin(1, npw_k1, nspinortot, ndat, gvnlxc_spin, gvnlxc_)

       ! Apply H_{k+q/2} to u^down
       call cg_copy_spin(2, npw_k1, nspinortot, ndat, cwavef, cwavef_spin)
       call nonlop(choice, cpopt_here, cwaveprj, enlout, gs_ham, idir, lambda_ndat, mpi_enreg, ndat, &
                   nnlout, paw_opt, signs, gsc_ptr, tim_nonlop, cwavef_spin, gvnlxc_spin, select_k=KPRIME_H_KPRIME)
       ! Insert results in the right position.
       call cg_put_spin(2, npw_k1, nspinortot, ndat, gvnlxc_spin, gvnlxc_)

       gs_ham%nspinor = 2
       ABI_FREE(cwavef_spin)
       ABI_FREE(gvnlxc_spin)
     end if ! use_gbt

     if (gs_ham%usepaw==1 .and. has_fock)then
       if (fock_get_getghc_call(fock)==1) then
         if(gs_ham%gpu_option==ABI_GPU_KOKKOS) then
#if defined HAVE_GPU && defined HAVE_YAKL
           ABI_MALLOC_MANAGED(gvnlc, (/2,npw_k2*my_nspinor*ndat/))
#endif
         else
           ABI_MALLOC(gvnlc, (2,npw_k2*my_nspinor*ndat))
#if defined HAVE_GPU && defined HAVE_YAKL
           !$OMP TARGET ENTER DATA MAP(to:gvnlc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
         end if
         if(gs_ham%gpu_option==ABI_GPU_OPENMP) then
           call gpu_copy(gvnlc, gvnlxc_, int(2,c_size_t)*npw_k2*my_nspinor*ndat)
         else
           gvnlc=gvnlxc_
         end if
       endif
     endif

     ! Calculation of the Fock exact exchange contribution from the Fock or ACE operator
     if (has_fock) then
       if (fock_get_getghc_call(fock)==1) then
         if (gs_ham%usepaw==0) cwaveprj_idat => cwaveprj
         if (fock%use_ACE==0) then
#if defined HAVE_GPU && defined HAVE_YAKL
           !$OMP TARGET UPDATE FROM(cwavef,gvnlxc_) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
           call timab(360,1,tsec)
           call fock_getghc(cwavef,cwaveprj,gvnlxc_,gs_ham,mpi_enreg,ndat)
           call timab(360,2,tsec)
#if defined HAVE_GPU && defined HAVE_YAKL
           !$OMP TARGET UPDATE TO(cwavef,gvnlxc_) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
         else
           call fock_ACE_getghc(cwavef,gvnlxc_,gs_ham,mpi_enreg,ndat)
         end if
       end if
     end if

   else if (type_calc == 3) then
     ! for kinetic and local only, nonlocal and vfock should be zero
     if(gs_ham%gpu_option==ABI_GPU_OPENMP) then
       call gpu_set_to_zero(gvnlxc_, int(2,c_size_t)*npw_k2*my_nspinor*ndat)
     else
       gvnlxc_(:,:) = zero
     end if
   end if ! if(type_calc...

   ABI_NVTX_END_RANGE()

   !============================================================
   ! Assemble kinetic, local, nonlocal and Fock contributions
   !============================================================

   ABI_NVTX_START_RANGE(NVTX_GETGHC_KIN)

#ifdef FC_NVHPC
   !FIXME This Kokkos kernel seems to cause issues under NVHPC so it is disabled
   if (.false.) then
#else
   if (gs_ham%gpu_option == ABI_GPU_KOKKOS) then
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS)
     call assemble_energy_contribution_kokkos(c_loc(ghc), &
       & c_loc(gsc), c_loc(kinpw_k2), c_loc(cwavef), c_loc(gvnlxc_), &
       & ndat, my_nspinor, npw_k2, sij_opt, k1_eq_k2, hugevalue)
     ! sync device so that data can be reused safely on host
     ! will probably be moved elsewhere once all the scf loop runs on device
     call gpu_device_synchronize()
#endif

   else

#ifdef FC_NVHPC
#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS)
     !Related to FIXME above
     if (gs_ham%gpu_option == ABI_GPU_KOKKOS) call gpu_device_synchronize()
#endif
#endif

     ! Assemble modified kinetic, local and nonlocal contributions
     ! to <G|H|C(n,k)>. Take also into account build-in debugging.
     if(prtvol/=-level)then

       if (gs_ham%gpu_option == ABI_GPU_OPENMP) then
         ! OpenMP GPU
#ifdef HAVE_OPENMP_OFFLOAD
         if (k1_eq_k2) then
           !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(to:ghc,kinpw_k2,gvnlxc_,gsc,cwavef) MAP(tofrom:kinpw_k2)
           do idat=1,ndat
             do ispinor=1,my_nspinor
               !$OMP PARALLEL DO PRIVATE(igspinor)
               do ig=1,npw_k2
                 igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
                 if(kinpw_k2(ig)<huge(zero)*1.d-11)then
                   ghc(:,igspinor) = ghc(:,igspinor) + kinpw_k2(ig)*cwavef(:,igspinor) + gvnlxc_(:,igspinor)
                 else
                   ghc(:,igspinor)=zero
                   if (sij_opt==1) gsc(:,igspinor)=zero
                 end if
               end do ! ig
             end do ! ispinor
           end do ! idat
         else
           !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(to:ghc,gvnlxc_,gsc) MAP(tofrom:kinpw_k2)
           do idat=1,ndat
             do ispinor=1,my_nspinor
               !$OMP PARALLEL DO PRIVATE(igspinor)
               do ig=1,npw_k2
                 igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
                 if(kinpw_k2(ig)<huge(zero)*1.d-11)then
                   ghc(:,igspinor)= ghc(:,igspinor) + gvnlxc_(:,igspinor)
                 else
                   ghc(:,igspinor)=zero
                   if (sij_opt==1) gsc(:,igspinor)=zero
                 end if
               end do ! ig
             end do ! ispinor
           end do ! idat
         end if
#endif

       else

         !CPU (+ Kokkos eventually)
         if (k1_eq_k2) then
           !$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2) IF(gemm_nonlop_use_gemm)
           do idat=1,ndat
             do ispinor=1,my_nspinor
               do ig=1,npw_k2
                 igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
                 if(kinpw_k2(ig)<huge(zero)*1.d-11)then
                   ghc(:,igspinor) = ghc(:,igspinor) + kinpw_k2(ig)*cwavef(:,igspinor) + gvnlxc_(:,igspinor)
                 else
                   ghc(:,igspinor)=zero
                   if (sij_opt==1) gsc(:,igspinor)=zero
                 end if
               end do ! ig
             end do ! ispinor
           end do ! idat
         else
           if (gs_ham%use_gbt == 0) then
             !$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2) IF(gemm_nonlop_use_gemm)
             do idat=1,ndat
               do ispinor=1,my_nspinor
                 do ig=1,npw_k2
                   igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
                   if(kinpw_k2(ig)<huge(zero)*1.d-11)then
                     ghc(:,igspinor)= ghc(:,igspinor) + gvnlxc_(:,igspinor)
                   else
                     ghc(:,igspinor)=zero
                     if (sij_opt==1) gsc(:,igspinor)=zero
                   end if
                 end do ! ig
               end do ! ispinor
             end do ! idat
           else
             !$OMP PARALLEL DO PRIVATE(igspinor, iispinor) COLLAPSE(2) IF(gemm_nonlop_use_gemm)
             do idat=1,ndat
               do ispinor=1,my_nspinor
                 iispinor=ispinor;if (mpi_enreg%paral_spinor==1) iispinor=mpi_enreg%me_spinor+1
                 if (iispinor == 1) then
                   do ig=1,npw_k2
                     igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
                     if(gs_ham%kinpw_k(ig)<huge(zero)*1.d-11)then
                       ghc(:,igspinor) = ghc(:,igspinor) + gs_ham%kinpw_k(ig)*cwavef(:,igspinor) + gvnlxc_(:,igspinor)
                     else
                       ghc(:,igspinor)=zero
                       if (sij_opt==1) gsc(:,igspinor)=zero
                     end if
                   end do ! ig
                else
                   do ig=1,npw_k2
                     igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
                     if(gs_ham%kinpw_kp(ig)<huge(zero)*1.d-11)then
                       ghc(:,igspinor) = ghc(:,igspinor) + gs_ham%kinpw_kp(ig)*cwavef(:,igspinor) + gvnlxc_(:,igspinor)
                     else
                       ghc(:,igspinor)=zero
                       if (sij_opt==1) gsc(:,igspinor)=zero
                     end if
                   end do ! ig
                 end if
               end do ! ispinor
             end do ! idat
           end if
         end if

       end if ! gs_ham%gpu_option

     else
       ! Here, debugging section
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET UPDATE FROM(ghc,gsc,cwavef,gvnlxc_) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
       call wrtout(std_out,' getghc : components of ghc ','PERS')
       write(msg,'(a)')&
         'icp ig ispinor igspinor re/im     ghc        kinpw         cwavef      glocc        gvnlxc  gsc'
       call wrtout(std_out,msg,'PERS')
       do idat=1,ndat
         do ispinor=1,my_nspinor
           do ig=1,npw_k2
             igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
             if(kinpw_k2(ig)<huge(zero)*1.d-11)then
               if (k1_eq_k2) then
                 ghcre=kinpw_k2(ig)*cwavef(1,igspinor)+ghc(1,igspinor)+gvnlxc_(1,igspinor)
                 ghcim=kinpw_k2(ig)*cwavef(2,igspinor)+ghc(2,igspinor)+gvnlxc_(2,igspinor)
               else
                 ghcre=ghc(1,igspinor)+gvnlxc_(1,igspinor)
                 ghcim=ghc(2,igspinor)+gvnlxc_(2,igspinor)
               end if
             else
               ghcre=zero
               ghcim=zero
               if (sij_opt==1) gsc(:,igspinor)=zero
             end if
             iispinor=ispinor;if (mpi_enreg%paral_spinor==1) iispinor=mpi_enreg%me_spinor+1
             if (sij_opt == 1) then
               write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
                 kinpw_k2(ig),cwavef(1,igspinor),ghc(1,igspinor),gvnlxc_(1,igspinor), gsc(1,igspinor)
               call wrtout(std_out,msg,'PERS')
               write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
                 kinpw_k2(ig),cwavef(2,igspinor),ghc(2,igspinor),gvnlxc_(2,igspinor), gsc(2,igspinor)
               call wrtout(std_out,msg,'PERS')
             else
               write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
                 kinpw_k2(ig),cwavef(1,igspinor),ghc(1,igspinor),gvnlxc_(1,igspinor)
               call wrtout(std_out,msg,'PERS')
               write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
                 kinpw_k2(ig),cwavef(2,igspinor),ghc(2,igspinor),gvnlxc_(2,igspinor)
               call wrtout(std_out,msg,'PERS')
             end if
             ghc(:,igspinor) = [ghcre, ghcim]
           end do ! ig
         end do ! ispinor
       end do ! idat
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET UPDATE TO(ghc,gsc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
     end if
   end if ! gs_ham%gpu_option

   ABI_NVTX_END_RANGE()

!  Special case of PAW + Fock : only return Fock operator contribution in gvnlxc_
   if (gs_ham%usepaw==1 .and. has_fock) then
     if(gs_ham%gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET DATA USE_DEVICE_ADDR(gvnlc,gvnlxc_)
       call abi_gpu_xaxpy(1, 2*npw_k2*my_nspinor*ndat, cminusone, c_loc(gvnlc), 1, c_loc(gvnlxc_), 1)
       !$OMP END TARGET DATA
       !$OMP TARGET EXIT DATA MAP(delete:gvnlc) IF(gs_ham%gpu_option == ABI_GPU_OPENMP)
#endif
     else
       gvnlxc_=gvnlxc_-gvnlc
     end if
     if(gs_ham%gpu_option==ABI_GPU_KOKKOS) then
#if defined HAVE_GPU && defined HAVE_YAKL
       ABI_FREE_MANAGED(gvnlc)
#endif
     else
       ABI_FREE(gvnlc)
     end if
   endif

#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE FROM(ghc)     IF(transfer_ghc)
   !$OMP TARGET UPDATE FROM(gsc)     IF(transfer_gsc)
   !$OMP TARGET UPDATE FROM(cwavef)  IF(transfer_cwavef)
   !$OMP TARGET UPDATE FROM(gvnlxc_) IF(transfer_gvnlxc .and. .not. local_gvnlxc)
   !$OMP TARGET EXIT DATA MAP(delete:gvnlxc_) IF(transfer_gvnlxc)
#endif
   if (local_gvnlxc) then
     if(gs_ham%gpu_option==ABI_GPU_KOKKOS) then
#if defined HAVE_GPU && defined HAVE_YAKL
       ABI_FREE_MANAGED(gvnlxc_)
#endif
     else
       ABI_FREE(gvnlxc_)
     end if
   end if

   ! Structured debugging: if prtvol=-level, stop here.
   if (prtvol == -level) then
     ABI_ERROR(sjoin(' getghc: exit prtvol=-',itoa(level),', debugging mode => stop '))
   end if

   if (type_calc==0.or.type_calc==2) then
     if (has_fock.and.gs_ham%usepaw==1.and.cpopt<2) then
       call pawcprj_free(cwaveprj_fock)
       ABI_FREE(cwaveprj_fock)
     end if
   end if

 end if ! type_calc

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET EXIT DATA MAP(delete:ghc)     IF(transfer_ghc)
 !$OMP TARGET EXIT DATA MAP(delete:gsc)     IF(transfer_gsc)
 !$OMP TARGET EXIT DATA MAP(delete:cwavef)  IF(transfer_cwavef)
#endif
 call timab(350+tim_getghc,2,tsec)

 DBG_EXIT("COLL")
 ABI_NVTX_END_RANGE()

end subroutine getghc
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/getghc_nucdip
!!
!! NAME
!! getghc_nucdip
!!
!! FUNCTION
!! Compute magnetic nuclear dipole moment contribution to <G|H|C>
!! for input vector |C> expressed in reciprocal space.
!!
!! INPUTS
!! cwavef(2,npw_k*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! gbound_k(2*mgfft+4)=sphere boundary info
!! gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!! istwf_k=input parameter that describes the storage of wfs
!! kg_k(3,npw_k)=G vec coordinates wrt recip lattice transl.
!! kpt(3)=current k point
!! mgfft=maximum single fft dimension
!! mpi_enreg=information about MPI parallelization
!! my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!! ndat=number of FFTs to perform in parall
!! ngfft(18)=contain all needed information about 3D FFT
!! npw_k=number of planewaves in basis for given k point.
!! nvloc=number of spin components of vxctaulocal
!! n4,n5,n6=for dimensioning of vxctaulocal
!! gpu_option= GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!! vectornd(n4,n5,n6,nvloc,3)= local potential corresponding to the vector potential of the array
!!  of nuclear magnetic dipoles, in real space, on the augmented fft grid.
!!
!! OUTPUT
!!  ghc_vectornd(2,npw_k*my_nspinor*ndat)=A.p contribution to <G|H|C> for array of nuclear dipoles
!!
!! SIDE EFFECTS
!!
!! NOTES
!! this code is a copied, simplified version of getghc_mGGA (see below) and should eventually be
!! integrated into that code, to simplify maintenance
!!
!! SOURCE

subroutine getghc_nucdip(cwavef,ghc_vectornd,gbound_k,istwf_k,kg_k,kpt,mgfft,mpi_enreg,&
&                      ndat,ngfft,npw_k,nvloc,n4,n5,n6,my_nspinor,vectornd,vlocal,zora,gpu_option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,my_nspinor,ndat,npw_k,nvloc,n4,n5,n6,gpu_option,zora
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: gbound_k(2*mgfft+4),kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(inout) :: cwavef(2,npw_k*my_nspinor*ndat)
 real(dp),intent(inout) :: ghc_vectornd(2,npw_k*my_nspinor*ndat)
 real(dp),intent(inout) :: vectornd(n4,n5,n6,nvloc,3),vlocal(n4,n5,n6,nvloc)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf=1
 integer :: idat,idir,ipw,iv1,iv2,nspinortot,shift
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc,usezora
 real(dp), parameter :: HalfFineStruct2=half/InvFineStruct**2
 real(dp) :: weight=one
 !arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: gcwavef(:,:,:),gcwavef1(:,:,:),gcwavef2(:,:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),kgkpk(:,:),vectornd_dir(:,:,:,:)
 real(dp),allocatable :: work(:,:,:,:),zk(:,:,:)
! *********************************************************************

 ghc_vectornd(:,:)=zero
 if (nvloc/=1) return

 nspinortot=min(2,(1+mpi_enreg%paral_spinor)*my_nspinor)
 if (mpi_enreg%paral_spinor==0) then
   shift=npw_k
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(nspinortot==2)
 else
   shift=0
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

 usezora=((zora.EQ.1).OR.(zora.EQ.3))
 if(usezora) then
   ABI_MALLOC(zk,(n4,n5,n6))
   zk(1:n4,1:1:n5,1:n6)=1.0/(1.0-HalfFineStruct2*vlocal(1:n4,1:n5,1:n6,nvloc))
 end if

 ABI_MALLOC(work,(2,n4,n5,n6*ndat))

 if (nspinortot==1) then

    ABI_MALLOC(ghc1,(2,npw_k*ndat))

    !  Do it in 2 STEPs:
    !  STEP1: Compute grad of cwavef
    ABI_MALLOC(gcwavef,(2,npw_k*ndat,3))

    gcwavef = zero

    ! compute k + G. Note these are in reduced coords
    ABI_MALLOC(kgkpk,(npw_k,3))
    do ipw = 1, npw_k
       kgkpk(ipw,:) = kpt(:) + kg_k(:,ipw)
    end do

    ! make 2\pi(k+G)c(G)|G> by element-wise multiplication
    do idir = 1, 3
       do idat = 1, ndat
          iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
          gcwavef(1,iv1:iv2,idir) = cwavef(1,iv1:iv2)*kgkpk(1:npw_k,idir)
          gcwavef(2,iv1:iv2,idir) = cwavef(2,iv1:iv2)*kgkpk(1:npw_k,idir)
       end do
    end do
    ABI_FREE(kgkpk)
    gcwavef = gcwavef*two_pi

    !  STEP2: Compute sum of (grad components of vectornd)*(grad components of cwavef)
    ABI_MALLOC(vectornd_dir,(n4,n5,n6,nvloc))
    do idir=1,3
      if (usezora) then
        vectornd_dir(1:n4,1:n5,1:n6,nvloc)=zk(1:n4,1:n5,1:n6)*vectornd(1:n4,1:n5,1:n6,nvloc,idir)
      else
        vectornd_dir(1:n4,1:n5,1:n6,nvloc)=vectornd(1:n4,1:n5,1:n6,nvloc,idir)
      end if
      call fourwf(1,vectornd_dir,gcwavef(:,:,idir),ghc1,work,gbound_k,gbound_k,&
           istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
           &     tim_fourwf,weight,weight,gpu_option=gpu_option)
!      call fourwf(1,vectornd(:,:,:,:,idir),gcwavef(:,:,idir),ghc1,work,gbound_k,gbound_k,&
!           istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
!           &     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
      ghc_vectornd=ghc_vectornd+ghc1
    end do ! idir
    ABI_FREE(vectornd_dir)
    ABI_FREE(gcwavef)
    ABI_FREE(ghc1)

 else ! nspinortot==2

    ABI_MALLOC(cwavef1,(2,npw_k*ndat))
    ABI_MALLOC(cwavef2,(2,npw_k*ndat))
    do idat=1,ndat
       iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
       cwavef1(1:2,iv1:iv2) = cwavef(1:2,1+(idat-1)*my_nspinor*npw_k:npw_k+(idat-1)*my_nspinor*npw_k)
       cwavef2(1:2,iv1:iv2) = &
         & cwavef(1:2,1+(idat-1)*my_nspinor*npw_k+shift:npw_k+(idat-1)*my_nspinor*npw_k+shift)
    end do

    ! compute k + G. Note these are in reduced coords
    ABI_MALLOC(kgkpk,(npw_k,3))
    do ipw = 1, npw_k
       kgkpk(ipw,:) = kpt(:) + kg_k(:,ipw)
    end do

    if (nspinor1TreatedByThisProc) then

       ABI_MALLOC(ghc1,(2,npw_k*ndat))

       !  Do it in 2 STEPs:
       !  STEP1: Compute grad of cwavef
       ABI_MALLOC(gcwavef1,(2,npw_k*ndat,3))

       gcwavef1 = zero
       ! make 2\pi(k+G)c(G)|G> by element-wise multiplication
       do idir = 1, 3
         do idat = 1, ndat
           iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
           gcwavef1(1,iv1:iv2,idir) = cwavef1(1,iv1:iv2)*kgkpk(1:npw_k,idir)
           gcwavef1(2,iv1:iv2,idir) = cwavef1(2,iv1:iv2)*kgkpk(1:npw_k,idir)
         end do
       end do
       gcwavef1 = gcwavef1*two_pi

       !  STEP2: Compute sum of (grad components of vectornd)*(grad components of cwavef)
       ABI_MALLOC(vectornd_dir,(n4,n5,n6,nvloc))
       do idir=1,3
         if (usezora) then
           vectornd_dir(1:n4,1:n5,1:n6,nvloc)=zk(1:n4,1:n5,1:n6)*vectornd(1:n4,1:n5,1:n6,nvloc,idir)
         else
           vectornd_dir(1:n4,1:n5,1:n6,nvloc)=vectornd(1:n4,1:n5,1:n6,nvloc,idir)
         end if
         call fourwf(1,vectornd_dir,gcwavef1(:,:,idir),ghc1,work,gbound_k,gbound_k,&
           & istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
           & tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
         do idat=1,ndat
           iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
           ghc_vectornd(1:2,iv1:iv2)=ghc_vectornd(1:2,iv1:iv2)+&
             &  ghc1(1:2,iv1:iv2)
         end do
       end do ! idir
       ABI_FREE(vectornd_dir)
       ABI_FREE(gcwavef1)
       ABI_FREE(ghc1)

    end if ! end spinor 1

    if (nspinor2TreatedByThisProc) then

       ABI_MALLOC(ghc2,(2,npw_k*ndat))

       !  Do it in 2 STEPs:
       !  STEP1: Compute grad of cwavef
       ABI_MALLOC(gcwavef2,(2,npw_k*ndat,3))
       gcwavef2 = zero
       ! make 2\pi(k+G)c(G)|G> by element-wise multiplication
       do idir = 1, 3
         do idat = 1, ndat
           iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
           gcwavef2(1,iv1:iv2,idir) = cwavef2(1,iv1:iv2)*kgkpk(1:npw_k,idir)
           gcwavef2(2,iv1:iv2,idir) = cwavef2(2,iv1:iv2)*kgkpk(1:npw_k,idir)
          end do
       end do
       gcwavef2 = gcwavef2*two_pi

       !  STEP2: Compute sum of (grad components of vectornd)*(grad components of cwavef)
       ABI_MALLOC(vectornd_dir,(n4,n5,n6,nvloc))
       do idir=1,3
         if (usezora) then
           vectornd_dir(1:n4,1:n5,1:n6,nvloc)=zk(1:n4,1:n5,1:n6)*vectornd(1:n4,1:n5,1:n6,nvloc,idir)
         else
           vectornd_dir(1:n4,1:n5,1:n6,nvloc)=vectornd(1:n4,1:n5,1:n6,nvloc,idir)
         end if
         call fourwf(1,vectornd_dir,gcwavef2(:,:,idir),ghc2,work,gbound_k,gbound_k,&
           & istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
           & tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
         do idat=1,ndat
           iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
           ghc_vectornd(1:2,iv1+shift:iv2+shift)=ghc_vectornd(1:2,iv1+shift:iv2+shift)+&
             & ghc2(1:2,iv1:iv2)
         end do
       end do ! idir
       ABI_FREE(vectornd_dir)
       ABI_FREE(gcwavef2)
       ABI_FREE(ghc2)

    end if ! end spinor 2

    ABI_FREE(cwavef1)
    ABI_FREE(cwavef2)
    ABI_FREE(kgkpk)

 end if ! nspinortot

 ABI_FREE(work)
 if(usezora) then
   ABI_FREE(zk)
 end if

end subroutine getghc_nucdip
!!***

!!****f* ABINIT/getghc_mGGA
!!
!! NAME
!! getghc_mGGA
!!
!! FUNCTION
!! Compute metaGGA contribution to <G|H|C> for input vector |C> expressed in reciprocal space.
!!
!! INPUTS
!! cwavef(2,npw_k*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! gbound_k(2*mgfft+4)=sphere boundary info
!! gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!! istwf_k=input parameter that describes the storage of wfs
!! kg_k(3,npw_k)=G vec coordinates wrt recip lattice transl.
!! kpt(3)=current k point
!! mgfft=maximum single fft dimension
!! mpi_enreg=information about MPI parallelization
!! my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!! ndat=number of FFTs to perform in parall
!! ngfft(18)=contain all needed information about 3D FFT
!! npw_k=number of planewaves in basis for given k point.
!! nvloc=number of spin components of vxctaulocal
!! n4,n5,n6=for dimensionning of vxctaulocal
!! gpu_option= GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!! vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!  kinetic energy density, in real space, on the augmented fft grid.
!!  This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!
!! OUTPUT
!!  ghc_mGGA(2,npw_k*my_nspinor*ndat)=metaGGA contribution to <G|H|C>
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine getghc_mGGA(cwavef,ghc_mGGA,gbound_k,gprimd,istwf_k,kg_k,kpt,mgfft,mpi_enreg,&
&                      ndat,ngfft,npw_k,nvloc,n4,n5,n6,my_nspinor,vxctaulocal,gpu_option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,my_nspinor,ndat,npw_k,nvloc,n4,n5,n6,gpu_option
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: gbound_k(2*mgfft+4),kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),kpt(3)
 real(dp),intent(inout) :: cwavef(2,npw_k*my_nspinor*ndat)
 real(dp),intent(inout) :: ghc_mGGA(2,npw_k*my_nspinor*ndat)
 real(dp),intent(inout) :: vxctaulocal(n4,n5,n6,nvloc,4)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf=1
 integer :: idat,idir,ipw,nspinortot,shift
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: weight=one
!arrays
 real(dp) :: kg_k_cart_vec(3)
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: gcwavef(:,:,:),gcwavef1(:,:,:),gcwavef2(:,:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:)
 real(dp),allocatable :: lcwavef(:,:),lcwavef1(:,:),lcwavef2(:,:)
 real(dp),allocatable :: work(:,:,:,:)
! *********************************************************************

 ghc_mGGA(:,:)=zero

 if (nvloc/=1) return

 nspinortot=min(2,(1+mpi_enreg%paral_spinor)*my_nspinor)
 if (mpi_enreg%paral_spinor==0) then
   shift=npw_k
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(nspinortot==2)
 else
   shift=0
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

 ABI_MALLOC(work,(2,n4,n5,n6*ndat))

 if (nspinortot==1) then

   ABI_MALLOC(ghc1,(2,npw_k*ndat))

!  Do it in 3 STEPs:
!  STEP1: Compute grad of cwavef and Laplacian of cwavef
   ABI_MALLOC(gcwavef,(2,npw_k*ndat,3))
   ABI_MALLOC(lcwavef,(2,npw_k*ndat))
!!$OMP PARALLEL DO
   gcwavef = zero; lcwavef = zero
   do idat=1,ndat
     do ipw=1,npw_k
       ! convert k + G from reduced coords to Cartesian
       kg_k_cart_vec = two_pi*MATMUL(gprimd,kpt(1:3)+kg_k(1:3,ipw))
       ! form \grad\psi = i(k + G) \psi in Cartesian frame
       gcwavef(1,ipw+(idat-1)*npw_k,1:3)=  cwavef(2,ipw+(idat-1)*npw_k)*kg_k_cart_vec(1:3)
       gcwavef(2,ipw+(idat-1)*npw_k,1:3)= -cwavef(1,ipw+(idat-1)*npw_k)*kg_k_cart_vec(1:3)
       ! form \nabla^2\psi = -|k + G|^2 \psi
       lcwavef(1:2,ipw+(idat-1)*npw_k)=&
         &lcwavef(1:2,ipw+(idat-1)*npw_k)-cwavef(1:2,ipw+(idat-1)*npw_k)*DOT_PRODUCT(kg_k_cart_vec,kg_k_cart_vec)
     end do
   end do
!  STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
   call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef,ghc1,work,gbound_k,gbound_k,&
&   istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
&   tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
   do idat=1,ndat
     do ipw=1,npw_k
        ghc_mGGA(:,ipw+(idat-1)*npw_k)=ghc_mGGA(:,ipw+(idat-1)*npw_k)-half*ghc1(:,ipw+(idat-1)*npw_k)
     end do
   end do
   ABI_FREE(lcwavef)
!  STEP3: Compute sum of (grad components of vxctaulocal)*(grad components of cwavef)
!  note: since grad cwavef is in Cart frame, evidently grad vxc is also
   do idir=1,3
     call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef(:,:,idir),ghc1,work,gbound_k,gbound_k,&
     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
&     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
     do idat=1,ndat
       do ipw=1,npw_k
          ghc_mGGA(:,ipw+(idat-1)*npw_k)=ghc_mGGA(:,ipw+(idat-1)*npw_k)-half*ghc1(:,ipw+(idat-1)*npw_k)
       end do
     end do
   end do ! idir
   ABI_FREE(gcwavef)
   ABI_FREE(ghc1)

 else ! nspinortot==2

   ABI_MALLOC(cwavef1,(2,npw_k*ndat))
   ABI_MALLOC(cwavef2,(2,npw_k*ndat))
   do idat=1,ndat
     do ipw=1,npw_k
       cwavef1(1:2,ipw+(idat-1)*npw_k)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k)
       cwavef2(1:2,ipw+(idat-1)*npw_k)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k+shift)
     end do
   end do
!  call cg_zcopy(npw*ndat,cwavef(1,1),cwavef1)
!  call cg_zcopy(npw*ndat,cwavef(1,1+shift),cwavef2)


   if (nspinor1TreatedByThisProc) then

     ABI_MALLOC(ghc1,(2,npw_k*ndat))

!    Do it in 3 STEPs:
!    STEP1: Compute grad of cwavef and Laplacian of cwavef
     ABI_MALLOC(gcwavef1,(2,npw_k*ndat,3))
     ABI_MALLOC(lcwavef1,(2,npw_k*ndat))
     gcwavef1 = zero; lcwavef1 = zero
!!$OMP PARALLEL DO
      do idat=1,ndat
        do ipw=1,npw_k
          ! convert k + G from reduced coords to Cartesian
          kg_k_cart_vec = two_pi*MATMUL(gprimd,kpt(1:3)+kg_k(1:3,ipw))
          ! form \grad\psi = i(k + G) \psi in Cartesian frame
          gcwavef1(1,ipw+(idat-1)*npw_k,1:3)=  cwavef1(2,ipw+(idat-1)*npw_k)*kg_k_cart_vec(1:3)
          gcwavef1(2,ipw+(idat-1)*npw_k,1:3)= -cwavef1(1,ipw+(idat-1)*npw_k)*kg_k_cart_vec(1:3)
          ! form \nabla^2\psi = -|k + G|^2 \psi
          lcwavef1(1:2,ipw+(idat-1)*npw_k)=&
            &lcwavef1(1:2,ipw+(idat-1)*npw_k)-cwavef1(1:2,ipw+(idat-1)*npw_k)*DOT_PRODUCT(kg_k_cart_vec,kg_k_cart_vec)
        end do
      end do
!    STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
     call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef1,ghc1,work,gbound_k,gbound_k,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
&     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
     do idat=1,ndat
       do ipw=1,npw_k
         ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)=ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)-half*ghc1(:,ipw+(idat-1)*npw_k)
       end do
     end do
     ABI_FREE(lcwavef1)
!    STEP3: Compute (grad components of vxctaulocal)*(grad components of cwavef)
     do idir=1,3
       call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef1(:,:,idir),ghc1,work,gbound_k,gbound_k,&
       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
&      tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
       do idat=1,ndat
         do ipw=1,npw_k
           ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k) = ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)-half*ghc1(:,ipw+(idat-1)*npw_k)
         end do
       end do
     end do ! idir
     ABI_FREE(gcwavef1)
     ABI_FREE(ghc1)

   end if ! spin 1 treated by this proc

   if (nspinor2TreatedByThisProc) then

     ABI_MALLOC(ghc2,(2,npw_k*ndat))

!    Do it in 3 STEPs:
!    STEP1: Compute grad of cwavef and Laplacian of cwavef
     ABI_MALLOC(gcwavef2,(2,npw_k*ndat,3))
     ABI_MALLOC(lcwavef2,(2,npw_k*ndat))
!!$OMP PARALLEL DO
     gcwavef2 = zero; lcwavef2 = zero
     do idat=1,ndat
       do ipw=1,npw_k
         ! convert k + G from reduced coords to Cartesian
         kg_k_cart_vec = two_pi*MATMUL(gprimd,kpt(1:3)+kg_k(1:3,ipw))
         ! form \grad\psi = i(k + G) \psi in Cartesian frame
         gcwavef2(1,ipw+(idat-1)*npw_k,1:3)=  cwavef2(2,ipw+(idat-1)*npw_k)*kg_k_cart_vec(1:3)
         gcwavef2(2,ipw+(idat-1)*npw_k,1:3)= -cwavef2(1,ipw+(idat-1)*npw_k)*kg_k_cart_vec(1:3)
         ! form \nabla^2\psi = -|k + G|^2 \psi
         lcwavef2(1:2,ipw+(idat-1)*npw_k)=&
           &lcwavef2(1:2,ipw+(idat-1)*npw_k)-cwavef2(1:2,ipw+(idat-1)*npw_k)*DOT_PRODUCT(kg_k_cart_vec,kg_k_cart_vec)
       end do
     end do
!    STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
     call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef2,ghc2,work,gbound_k,gbound_k,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
&     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
     do idat=1,ndat
        do ipw=1,npw_k
           ! original code
           ! ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)=ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)-half*ghc2(:,ipw+(idat-1)*npw_k)
           ! but this stores the spinor2 result in the spinor1 location. Should be stored with shift
           ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k+shift)=ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k+shift)&
                & -half*ghc2(:,ipw+(idat-1)*npw_k)
       end do
     end do
     ABI_FREE(lcwavef2)
!    STEP3: Compute sum of (grad components of vxctaulocal)*(grad components of cwavef)
     do idir=1,3
       call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef2(:,:,idir),ghc2,work,gbound_k,gbound_k,&
       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
&      tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
       do idat=1,ndat
         do ipw=1,npw_k
           ! original code
           ! ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)=ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k)-half*ghc2(:,ipw+(idat-1)*npw_k)
           ! but this stores the spinor2 result in the spinor1 location. Should be stored with shift
            ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k+shift)=ghc_mGGA(:,ipw+(idat-1)*my_nspinor*npw_k+shift)&
                 & -half*ghc2(:,ipw+(idat-1)*npw_k)
         end do
       end do
     end do ! idir

     ABI_FREE(gcwavef2)
     ABI_FREE(ghc2)

   end if ! spin 2 treated by this proc

   ABI_FREE(cwavef1)
   ABI_FREE(cwavef2)

 end if ! nspinortot

 ABI_FREE(work)

end subroutine getghc_mGGA
!!***

!!****f* ABINIT/cwavef_double_rfft_trick_pack
!!
!! NAME
!! cwavef_double_rfft_trick_pack
!!
!! FUNCTION
!!
!! We have C(G)=C(-G)^* and D(G)=D(-G)^* and only G components are in memory (istwfk=2)
!! The Fourier transform is:
!! C(r) = sum_G e^(iGr) C(G) = sum_(G_z>=0,G/=0) 2 Re[e^(iGr) C(G)] + C(0)
!! so C(r) is a real function (same for D)
!! Here we construct:
!! E( G) = C(G)   + i D(G)
!! E(-G) = C(G)^* + i D(G)^* (G/=0)
!! so:
!! E(r) = C(r) + i D(r)
!! In short, one can do only one FFT on E (with istwfk=1) and obtains the FFT of C and D (istwfk=2)
!!
!! INPUTS
!! cwavef(2,npw_k*my_nspinor*(2*ndat))=planewave coefficients of wavefunctioni (istwfk=2)
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFTs to perform in parall
!! npw_k=number of planewaves in basis for given k point (istwfk=2).
!!
!! OUTPUT
!! cwavef_fft(2,npw_fft*my_nspinor*ndat)=planewave coefficients of wavefunction (istwfk=1)
!!
!! SOURCE
!!
subroutine cwavef_double_rfft_trick_pack(cwavef,cwavef_fft,me_g0,ndat,npw_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,npw_k,me_g0
!arrays
 real(dp),intent(in) :: cwavef(:,:)
 real(dp),intent(out) :: cwavef_fft(:,:)
!Local variables-------------------------------
!scalars
 integer :: idat,i0,ib1,ib2,npw_fft,ndat_
 logical :: ndat_is_odd

 ndat_ = ndat / 2
 ndat_is_odd=.False.
 if (modulo(ndat,2)/=0) then
   ndat_=ndat_+1
   ndat_is_odd=.True.
 end if

 npw_fft=2*npw_k
 i0=1
 if (me_g0==1) then! Do not include G=(0,0,0) twice
   npw_fft=npw_fft-1
   i0=2
 end if

 if (size(cwavef,1)/=2) then
   ABI_BUG('wrong size for cwavef (dim 1)')
 end if
 if (ndat_is_odd) then
   if (size(cwavef,2)/=npw_k*(2*ndat_-1)) then
     ABI_BUG('wrong size for cwavef (dim 2) (odd)')
   end if
 else
   if (size(cwavef,2)/=npw_k*2*ndat_) then
     ABI_BUG('wrong size for cwavef (dim 2)')
   end if
 end if
 if (size(cwavef_fft,1)/=2.or.size(cwavef_fft,2)/=npw_fft*ndat_) then
   ABI_BUG('wrong size for cwavef_fft')
 end if

 do idat=1,ndat_
   ib1=(idat-1)*npw_fft ! band shift for cwavef_fft
   ib2=(idat-1)*2*npw_k ! band shift for cwavef
   if (.not.ndat_is_odd.or.idat<ndat_) then
     ! E(G) = C(G) + i D(G)
     cwavef_fft(1,1+ib1:npw_k+ib1) = cwavef(1,1+ib2      :  npw_k+ib2) &
                                    -cwavef(2,1+npw_k+ib2:2*npw_k+ib2)
     cwavef_fft(2,1+ib1:npw_k+ib1) = cwavef(2,1+ib2      :  npw_k+ib2) &
                                    +cwavef(1,1+npw_k+ib2:2*npw_k+ib2)
     ! E(-G) = C(G)^* + i D(G)^* (G/=0)
     cwavef_fft(1,1+npw_k+ib1:npw_fft+ib1) = cwavef(1,i0      +ib2:  npw_k+ib2) &
                                            +cwavef(2,i0+npw_k+ib2:2*npw_k+ib2)
     cwavef_fft(2,1+npw_k+ib1:npw_fft+ib1) =-cwavef(2,i0      +ib2:  npw_k+ib2) &
                                            +cwavef(1,i0+npw_k+ib2:2*npw_k+ib2)
   else ! idat=ndat_ and ndat_is_odd : the vector D does not exist
     ! E(G) = C(G)
     cwavef_fft(1,1+ib1:npw_k+ib1) = cwavef(1,1+ib2:npw_k+ib2)
     cwavef_fft(2,1+ib1:npw_k+ib1) = cwavef(2,1+ib2:npw_k+ib2)
     ! E(-G) = C(G)^* (G/=0)
     cwavef_fft(1,1+npw_k+ib1:npw_fft+ib1) = cwavef(1,i0+ib2:npw_k+ib2)
     cwavef_fft(2,1+npw_k+ib1:npw_fft+ib1) =-cwavef(2,i0+ib2:npw_k+ib2)
   end if
 end do

end subroutine cwavef_double_rfft_trick_pack
!!***

!!****f* ABINIT/cwavef_double_rfft_trick_unpack
!!
!! NAME
!! cwavef_double_rfft_trick_unpack
!!
!! FUNCTION
!!
!! From the "cwavef_double_rfft_trick_pack" routine we have:
!! E( G) = C(G)   + i D(G)
!! E(-G) = C(G)^* + i D(G)^* (G/=0)
!! Here we compute:
!! C(G) = ( E(G) + E(-G)^* ) / 2   (G/=0)
!! D(G) = ( iE(-G)^* - iE(G)) / 2  (G/=0)
!! and:
!! C(0) = Re(E(0))
!! D(0) = Im(E(0))
!!
!! INPUTS
!! cwavef_fft(2,npw_fft*my_nspinor*ndat_)=planewave coefficients of wavefunction (istwfk=1)
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFTs to perform in parall
!! npw_k=number of planewaves in basis for given k point (istwfk=2).
!!
!! OUTPUT
!! cwavef(2,npw_fft*my_nspinor*(2*ndat))=planewave coefficients of wavefunction (istwfk=2)
!!
!! SOURCE
!!
subroutine cwavef_double_rfft_trick_unpack(cwavef,cwavef_fft,me_g0,ndat,npw_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,npw_k,me_g0
!arrays
 real(dp),intent(out) :: cwavef(:,:)
 real(dp),intent(in) :: cwavef_fft(:,:)
!Local variables-------------------------------
!scalars
 integer :: idat,i0,ib1,ib2,npw_fft,ndat_
 logical :: ndat_is_odd

 ndat_ = ndat / 2
 ndat_is_odd=.False.
 if (modulo(ndat,2)/=0) then
   ndat_=ndat_+1
   ndat_is_odd=.True.
 end if

 npw_fft=2*npw_k
 i0=1
 if (me_g0==1) then! Do not include G=(0,0,0) twice
   npw_fft=npw_fft-1
   i0=2
 end if

 if (size(cwavef,1)/=2) then
   ABI_BUG('wrong size for cwavef (dim 1)')
 end if
 if (ndat_is_odd) then
   if (size(cwavef,2)/=npw_k*(2*ndat_-1)) then
     ABI_BUG('wrong size for cwavef (dim 2) (odd)')
   end if
 else
   if (size(cwavef,2)/=npw_k*2*ndat_) then
     ABI_BUG('wrong size for cwavef (dim 2)')
   end if
 end if
 if (size(cwavef_fft,1)/=2.or.size(cwavef_fft,2)/=npw_fft*ndat_) then
   ABI_BUG('wrong size for cwavef_fft')
 end if

 do idat=1,ndat_
   ib1=(idat-1)*npw_fft  ! band shift for cwavef_fft
   ib2=(idat-1)*2*npw_k ! band shift for cwavef
   ! C(G) = ( E(G) + E(-G)^* ) / 2 (factor 1/2 will be applied later)
   cwavef(1,i0+ib2:npw_k+ib2) = cwavef_fft(1,i0      +ib1:npw_k  +ib1) & !+Re(E( G))
                               +cwavef_fft(1,1 +npw_k+ib1:npw_fft+ib1)   !+Re(E(-G))
   cwavef(2,i0+ib2:npw_k+ib2) = cwavef_fft(2,i0      +ib1:npw_k  +ib1) & !+Im(E( G))
                               -cwavef_fft(2,1 +npw_k+ib1:npw_fft+ib1)   !-Im(E(-G))
   if (.not.ndat_is_odd.or.idat<ndat_) then
     ! D(G) = ( iE(-G)^* - iE(G) ) / 2 (factor 1/2 will be applied later)
     cwavef(1,i0+npw_k+ib2:2*npw_k+ib2) = cwavef_fft(2,i0      +ib1:npw_k  +ib1) & !+Im(E( G))
                                         +cwavef_fft(2,1 +npw_k+ib1:npw_fft+ib1)   !+Im(E(-G))
     cwavef(2,i0+npw_k+ib2:2*npw_k+ib2) =-cwavef_fft(1,i0      +ib1:npw_k  +ib1) & !-Re(E( G))
                                         +cwavef_fft(1,1 +npw_k+ib1:npw_fft+ib1)   !+Re(E(-G))
   end if
   if (me_g0==1) then
     ! Compute C(G=0) and D(G=0) and multiply by 2 as we apply 1/2 to the whole array shortly afterwards
     ! C(G=0) = Re(E(G=0))
     cwavef(1,1+ib2) = two*cwavef_fft(1,1+ib1)
     cwavef(2,1+ib2) = zero
     if (.not.ndat_is_odd.or.idat<ndat_) then
       ! D(G=0) = Im(E(G=0))
       cwavef(1,1+npw_k+ib2) = two*cwavef_fft(2,1+ib1)
       cwavef(2,1+npw_k+ib2) = zero
     end if
   end if

 end do

 cwavef(:,:) = half*cwavef(:,:)

end subroutine cwavef_double_rfft_trick_unpack
!!***

!!****f* ABINIT/getgsc
!! NAME
!! getgsc
!!
!! FUNCTION
!! Compute <G|S|C> for all input vectors |Cnk> at a given k-point,
!!              OR for one input vector |Cnk>.
!! |Cnk> are expressed in reciprocal space.
!! S is the overlap operator between |Cnk> (used for PAW).
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  cprj(natom,mcprj)= wave functions projected with non-local projectors: cprj=<p_i|Cnk>
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj (beginning of current k-point)
!!  icg=shift to be applied on the location of data in the array cg (beginning of current k-point)
!!  igsc=shift to be applied on the location of data in the array gsc (beginning of current k-point)
!!  ikpt,isppol=indexes of current (spin.kpoint)
!!  mcg=second dimension of the cg array
!!  mcprj=second dimension of the cprj array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=information about MPI parallelization
!!  ndat=number of bands to compute in parallel
!!  natom=number of atoms in unit cell.
!!  nband= if positive: number of bands at this k point for that spin polarization
!!         if negative: abs(nband) is the index of the only band to be computed
!!  npw_k=number of planewaves in basis for given k point.
!!  nspinor=number of spinorial components of the wavefunctions
!! [select_k]=optional, option governing the choice of k points to be used.
!!             gs_ham datastructure contains quantities needed to apply overlap operator
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|S|k>       is applied [default]
!!             if select_k=2, <k|S|k^prime>       is applied
!!             if select_k=3, <k|S|k>             is applied
!!             if select_k=4, <k^prime|S|k^prime> is applied
!!
!! OUTPUT
!!  gsc(2,mgsc)= <g|S|Cnk> or <g|S^(1)|Cnk> (S=overlap)
!!
!! SOURCE

subroutine getgsc(cg,cprj,gs_ham,gsc,ibg,icg,igsc,ikpt,isppol,&
&                 mcg,mcprj,mgsc,mpi_enreg,ndat,natom,nband,npw_k,nspinor,select_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ibg,icg,igsc,ikpt,isppol,mcg,mcprj
 integer,intent(in) :: mgsc,natom,nband,npw_k,nspinor,ndat
!TODO : may be needed to distribute cprj over band procs
! integer,intent(in) :: mband_mem
 integer,intent(in),optional :: select_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 real(dp),intent(in),  target :: cg(2,mcg)
 real(dp),intent(out), target :: gsc(2,mgsc)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimenl1,dimenl2,iband,iband1,iband2,index_cg,index_cprj
 integer :: index_gsc,me,my_nspinor,my_ndat,paw_opt,select_k_,signs,tim_nonlop,useylm
 !character(len=500) :: msg
!arrays
 real(dp) :: enlout_dum(ndat),tsec(2)
 real(dp), ABI_CONTIGUOUS pointer :: cwavef(:,:),scwavef(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
! *********************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 my_ndat=ndat; if(nband<1) my_ndat=1
 if(gs_ham%usepaw==0) then
   ABI_BUG('Only compatible with PAW (usepaw=1) !')
 end if
 if(nband<0.and.(mcg<npw_k*my_nspinor.or.mgsc<npw_k*my_nspinor.or.mcprj<my_nspinor)) then
   ABI_BUG('Invalid value for mcg, mgsc or mcprj !')
 end if

!Keep track of total time spent in getgsc:
 call timab(565,1,tsec)

 if(gs_ham%gpu_option==ABI_GPU_DISABLED) then
   gsc = zero
 else  if(gs_ham%gpu_option==ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
   call gpu_set_to_zero(gsc,int(2,c_size_t)*mgsc)
#endif
 end if

!Prepare some data
 if (gs_ham%usecprj==1) then
   ABI_MALLOC(cwaveprj,(natom,my_nspinor*my_ndat))
   call pawcprj_alloc(cwaveprj,0,gs_ham%dimcprj)
 else
   ABI_MALLOC(cwaveprj,(0,0))
 end if
 dimenl1=gs_ham%dimekb1;dimenl2=natom;tim_nonlop=0
 choice=1;signs=2;cpopt=-1+3*gs_ham%usecprj;paw_opt=3;useylm=1
 select_k_=1;if (present(select_k)) select_k_=select_k
 me=mpi_enreg%me_kpt

!Loop over bands
 index_cprj=ibg;index_cg=icg;index_gsc=igsc
 if (nband>0) then
   iband1=1;iband2=nband
!only do 1 band in case nband < 0 (the |nband|th one)
 else if (nband<0) then
   iband1=abs(nband);iband2=iband1
   index_cprj=index_cprj+(iband1-1)*my_nspinor
   index_cg  =index_cg  +(iband1-1)*npw_k*my_nspinor
   index_gsc =index_gsc +(iband1-1)*npw_k*my_nspinor
 end if

 do iband=iband1,iband2,my_ndat

   if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me.and.nband>0) then
! No longer needed 28/03/2020 to parallelize memory
!     gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=zero
!     index_gsc=index_gsc+npw_k*my_nspinor
     !index_cprj=index_cprj+my_nspinor
     !index_cg=index_cg+npw_k*my_nspinor

     cycle
   end if

!  Retrieve WF at (n,k)
   cwavef(1:2,1:npw_k*my_nspinor*my_ndat)  => cg(:,1+index_cg:npw_k*my_nspinor*my_ndat+index_cg)
   scwavef(1:2,1:npw_k*my_nspinor*my_ndat) => gsc(:,1+index_gsc:npw_k*my_nspinor*my_ndat+index_gsc)
   if (gs_ham%usecprj==1) then
     call pawcprj_copy(cprj(:,1+index_cprj:my_nspinor*my_ndat+index_cprj),cwaveprj)
   end if

!  Compute <g|S|Cnk>
   call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_ham,0,(/zero/),mpi_enreg,my_ndat,1,paw_opt,&
&   signs,scwavef,tim_nonlop,cwavef,cwavef,select_k=select_k_)


!  End of loop over bands
   index_cprj=index_cprj+my_nspinor*my_ndat
   index_cg=index_cg+npw_k*my_nspinor*my_ndat
   index_gsc=index_gsc+npw_k*my_nspinor*my_ndat
 end do

!Memory deallocation
 if (gs_ham%usecprj==1) then
   call pawcprj_free(cwaveprj)
 end if
 ABI_FREE(cwaveprj)

 call timab(565,2,tsec)

 DBG_EXIT("COLL")

end subroutine getgsc
!!***

!!****f* ABINIT/multithreaded_getghc
!!
!! NAME
!! multithreaded_getghc
!!
!! FUNCTION
!!
!! INPUTS
!! cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!       (same meaning as in nonlop.F90 routine)
!!       if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!       if cpopt= 0, <p_lmn|in> are computed here and saved
!!       if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!       if cpopt= 2  <p_lmn|in> are already in memory;
!!       if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!       if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!! cwavef(2,npw*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFT to do in parallel
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!!            0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+kinetic only (added to the existing Hamiltonian)
!!            3: local + kinetic only (added to the existing Hamiltonian)
!! ===== Optional inputs =====
!!   [kg_fft_k(3,:)]=optional, (k+G) vector coordinates to be used for the FFT transformation
!!                   instead of the one contained in gs_ham datastructure.
!!                   Typically used for real WF (in parallel) which are FFT-transformed 2 by 2.
!!   [kg_fft_kp(3,:)]=optional, (k^prime+G) vector coordinates to be used for the FFT transformation
!!   [select_k]=optional, option governing the choice of k points to be used.
!!             gs_ham datastructure contains quantities needed to apply Hamiltonian
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|H|k>       is applied [default]
!!             if select_k=2, <k|H|k^prime>       is applied
!!             if select_k=3, <k|H|k>             is applied
!!             if select_k=4, <k^prime|H|k^prime> is applied
!!
!! OUTPUT
!!  ghc(2,npw*my_nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                          or <G|H-lambda.S|C> (if sij_opt=-1)
!!  gvnlxc(2,npw*my_nspinor*ndat)=matrix elements <G|Vnonlocal|C> (if sij_opt>=0)
!!                                            or <G|Vnonlocal-lambda.S|C> (if sij_opt=-1)
!!  if (sij_opt=1)
!!    gsc(2,npw*my_nspinor*ndat)=matrix elements <G|S|C> (S=overlap).
!!
!! SIDE EFFECTS
!!  cwaveprj(natom,my_nspinor*(1+cpopt)*ndat)= wave function projected on nl projectors (PAW only)
!!
!! SOURCE

subroutine multithreaded_getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
&                 prtvol,sij_opt,tim_getghc,type_calc,&
&                 kg_fft_k,kg_fft_kp,select_k,filter_dilatmx_loc) ! optional arguments

#ifdef HAVE_OPENMP
   use omp_lib
#endif

!Arguments ------------------------------------
!scalars
 logical,intent(in),optional :: filter_dilatmx_loc
 integer,intent(in) :: cpopt,ndat, prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 integer,intent(in),optional :: select_k
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
 real(dp),intent(out),target :: gsc(:,:)
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(out) :: ghc(:,:),gvnlxc(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 logical :: filter_dilatmx_loc_
 integer :: firstelt, firstprj, lastelt, lastprj,usegvnlxc,usegsc
 integer :: nthreads,fftalga
 integer :: ithread
 integer :: chunk
 integer :: residuchunk
 integer :: firstband
 integer :: lastband
 integer :: spacedim, spacedim_prj
 logical :: fftw3_use_lib_threads_sav
 integer :: select_k_default
 ! *************************************************************************

 select_k_default = 1; if (present(select_k)) select_k_default = select_k
 filter_dilatmx_loc_ = .true.; if ( present(filter_dilatmx_loc) ) filter_dilatmx_loc_ = filter_dilatmx_loc

 spacedim     = size(cwavef  ,dim=2)/ndat
 spacedim_prj = size(cwaveprj,dim=2)/ndat

 nthreads = xomp_get_num_threads(open_parallel=.True.)
 fftalga = gs_ham%ngfft(7)/100
 if (fftalga==FFT_SG.and.nthreads>1.and.ndat>1) then
   ABI_ERROR("fftalg=1XX is not thread-safe, so it cannot be used in multi-threaded hamiltonian with nthreads>1 and ndat>1.")
 end if

 ! Disabling multithreading for GPU variants (getghc_ompgpu is not thread-safe for now)
 !$omp parallel default (none) &
 !$omp& private(ithread,nthreads,chunk,firstband,lastband,residuchunk,firstelt,lastelt), &
 !$omp& private(firstprj,lastprj,usegvnlxc,usegsc,fftw3_use_lib_threads_sav), &
 !$omp& shared(cwavef,ghc,gsc, gvnlxc,spacedim,spacedim_prj,ndat,kg_fft_k,kg_fft_kp,gs_ham,cwaveprj,mpi_enreg), &
 !$omp& shared(gemm_nonlop_use_gemm), &
 !$omp& firstprivate(cpopt,lambda,prtvol,sij_opt,tim_getghc,type_calc,select_k_default,filter_dilatmx_loc_) &
 !$omp& IF(gs_ham%gpu_option==ABI_GPU_DISABLED .and. .not. gemm_nonlop_use_gemm)
 ithread = 0
 nthreads = 1
 fftw3_use_lib_threads_sav = .false.
 if(gs_ham%gpu_option==ABI_GPU_DISABLED .and. .not. gemm_nonlop_use_gemm) then
#ifdef HAVE_OPENMP
   ithread = omp_get_thread_num()
   nthreads = omp_get_num_threads()
!Ensure that libs are used without threads (mkl, openblas, fftw3, ...)
#ifdef HAVE_LINALG_MKL_THREADS
   call mkl_set_num_threads(1)
#endif
!LB-23/07/24: OpenBLAS detects parallel sections automatically. To comment this line improves performances for some cases.
!#ifdef HAVE_LINALG_OPENBLAS_THREADS
!   call openblas_set_num_threads(1)
!#endif
#ifdef HAVE_LINALG_NVPL_THREADS
   call nvpl_blas_set_num_threads(1)
#endif
#ifdef HAVE_FFTW3_THREADS
   fftw3_use_lib_threads_sav=(.not.fftw3_spawn_threads_here(nthreads,nthreads))
   call fftw3_use_lib_threads(.false.)
#endif
#endif
 end if
 chunk = ndat/nthreads ! Divide by 2 to construct chunk of even number of bands
 residuchunk = ndat - nthreads*chunk
 if ( ithread < nthreads-residuchunk ) then
   firstband = ithread*chunk+1
   lastband = (ithread+1)*chunk
 else
   firstband = (nthreads-residuchunk)*chunk + ( ithread -(nthreads-residuchunk) )*(chunk+1) +1
   lastband = firstband+chunk
 end if
 usegvnlxc=1
 if (size(gvnlxc)<=1) usegvnlxc=0
 usegsc=0
 if (gs_ham%usepaw==1.and.sij_opt==1) usegsc=1

 if ( lastband /= 0 ) then
   firstelt = (firstband-1)*spacedim+1
   firstprj = (firstband-1)*spacedim_prj+1
   lastelt = lastband*spacedim
   lastprj = lastband*spacedim_prj
      ! Don't know how to manage optional arguments .... :(
   if ( present(kg_fft_k) ) then
     if (present(kg_fft_kp)) then
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*usegsc),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_k=select_k_default,kg_fft_k=kg_fft_k,kg_fft_kp=kg_fft_kp,filter_dilatmx_loc=filter_dilatmx_loc_)
     else
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*usegsc),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_k=select_k_default,kg_fft_k=kg_fft_k,filter_dilatmx_loc=filter_dilatmx_loc_)
     end if
   else
     if (present(kg_fft_kp)) then
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*usegsc),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_k=select_k_default,kg_fft_kp=kg_fft_kp,filter_dilatmx_loc=filter_dilatmx_loc_)
     else
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*usegsc),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_K=select_k_default,filter_dilatmx_loc=filter_dilatmx_loc_)
     end if
   end if
 end if
 if(gs_ham%gpu_option==ABI_GPU_DISABLED .and. .not. gemm_nonlop_use_gemm) then
#ifdef HAVE_OPENMP
  !Restore libs behavior (mkl, openblas, fftw3, ...)
#ifdef HAVE_LINALG_MKL_THREADS
   call mkl_set_num_threads(nthreads)
#endif
!LB-23/07/24: OpenBLAS detects parallel sections automatically. To comment this line improves performances for some cases.
!#ifdef HAVE_LINALG_OPENBLAS_THREADS
!   call openblas_set_num_threads(nthreads)
!#endif
#ifdef HAVE_LINALG_NVPL_THREADS
   call nvpl_blas_set_num_threads(nthreads)
#endif
#ifdef HAVE_FFTW3_THREADS
   call fftw3_use_lib_threads(fftw3_use_lib_threads_sav)
#endif
#endif
 end if
!$omp end parallel

end subroutine multithreaded_getghc
!!***

end module m_getghc
!!***
