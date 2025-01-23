!!****m* ABINIT/m_getghc_ompgpu
!! NAME
!!  m_getghc_ompgpu
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space - OpenMP GPU version;
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

module m_getghc_ompgpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_abi_linalg
 use, intrinsic :: iso_c_binding

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_getdim, pawcprj_copy
 use m_hamiltonian, only : gs_hamiltonian_type, KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME
 use m_fock,        only : fock_common_type, fock_get_getghc_call
 use m_fock_getghc, only : fock_getghc, fock_ACE_getghc
 use m_nonlop,      only : nonlop
 use m_fft,         only : fourwf

 !FIXME Keep those in these modules or moves them together ?
 use m_ompgpu_fourwf,      only : ompgpu_fourwf_work_mem
 use m_gemm_nonlop_ompgpu, only : gemm_nonlop_ompgpu_work_mem

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#endif

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

 implicit none

 private
!!***

 public :: getghc_ompgpu      ! Compute <G|H|C> for input vector |C> expressed in reciprocal space

 ! These routines allows to manage work buffers outside of getghc_ompgpu
 public :: alloc_getghc_ompgpu_buffers
 public :: free_getghc_ompgpu_buffers

 ! This routine is here to assess memory requirements
 public :: getghc_ompgpu_work_mem
!!***
#ifdef HAVE_OPENMP_OFFLOAD
 integer, save :: mod__n4=0, mod__n5=0, mod__n6=0, mod__nspinor=0, mod__ndat=0, mod__npw=0
 integer, save :: buf_initialized=0
 real(dp),allocatable :: work(:,:,:,:)
#endif

contains
!!***


!!****f* ABINIT/alloc_getghc_ompgpu_buffers
!! NAME
!! alloc_getghc_ompgpu_buffers
!!
!! FUNCTION
!! Allocate getghc_ompgpu work buffer
!!
subroutine alloc_getghc_ompgpu_buffers(npw, nspinor, ndat, n4, n5, n6)
  integer, intent(in) :: npw, nspinor, ndat, n4, n5, n6

#ifdef HAVE_OPENMP_OFFLOAD
  if (buf_initialized == 1) then
    call free_getghc_ompgpu_buffers()
  end if

 ABI_MALLOC(work,(2, n4, n5, n6*ndat))
 mod__n4 = n4
 mod__n5 = n5
 mod__n6 = n6
 mod__nspinor = nspinor
 mod__ndat = ndat
 mod__npw = npw
 !FIXME Smater buffer management ?
#ifdef HAVE_GPU_HIP
  !Work buffer allocated once to save time in HIP (alloc costful)
 !$OMP TARGET ENTER DATA MAP(alloc:work)
#endif

 buf_initialized = 1

#else
 ABI_UNUSED((/npw,nspinor,ndat,n4,n5,n6/))
#endif
end subroutine alloc_getghc_ompgpu_buffers
!!***


!!****f* ABINIT/free_getghc_ompgpu_buffers
!! NAME
!! free_getghc_ompgpu_buffers
!!
!! FUNCTION
!! Free getghc_ompgpu work buffer
!!
subroutine free_getghc_ompgpu_buffers

#ifdef HAVE_OPENMP_OFFLOAD
 if(allocated(work)) then
   !FIXME Smater buffer management ?
#ifdef HAVE_GPU_HIP
   !$OMP TARGET EXIT DATA MAP(delete:work)
#endif
   ABI_FREE(work)
 end if

 buf_initialized = 0

#endif
end subroutine free_getghc_ompgpu_buffers
!!***


!!****f* ABINIT/getghc_ompgpu_work_mem
!! NAME
!! getghc_ompgpu_work_mem
!!
!! FUNCTION
!! Returns work memory requirement for getghc_ompgpu
!!
!! INPUTS
!!
!! gs_ham <type(gs_hamiltonian_type)>=contains dimensions of FFT domain
!! ndat=size of batch for fourwf and nonlop processing
!!
!! OUTPUT
!!
!! req_mem=amount in bytes of required memory for getghc_ompgpu
function getghc_ompgpu_work_mem(gs_ham, ndat) result(req_mem)

 implicit none

 type(gs_hamiltonian_type),intent(in),target :: gs_ham
 integer, intent(in) :: ndat
 integer(kind=c_size_t) :: req_mem, ghc_mem, nonlop_mem

 ! getghc use a GPU work buffer only when using fourwf
 ! Therefore, max GPU memory required by getghc is either:
 !   - the sum of getghc and fourwf work buffers memory requirements
 !   - the amount of memory required by gemm_nonlop_ompgpu work buffers
 ghc_mem = 0
 ghc_mem = int(2, c_size_t) * dp * gs_ham%n4 * gs_ham%n5 * gs_ham%n6 * ndat
 ghc_mem = ghc_mem + ompgpu_fourwf_work_mem(gs_ham%ngfft, ndat)

 nonlop_mem = gemm_nonlop_ompgpu_work_mem(gs_ham%istwf_k, ndat, gs_ham%npw_fft_k,&
 &               gs_ham%indlmn, gs_ham%nattyp, gs_ham%ntypat, gs_ham%lmnmax)

 req_mem = MAX(ghc_mem, nonlop_mem)

end function getghc_ompgpu_work_mem
!!***

!!****f* ABINIT/getghc_ompgpu
!! NAME
!! getghc_ompgpu
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
!!   [kg_fft_k(3,:)]=optional, (k+G) vector coordinates to be used for the FFT tranformation
!!                   instead of the one contained in gs_ham datastructure.
!!                   Typically used for real WF (in parallel) which are FFT-transformed 2 by 2.
!!   [kg_fft_kp(3,:)]=optional, (k^prime+G) vector coordinates to be used for the FFT tranformation
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
!! PARENTS
!!      m_cgwf,m_cgwf_cprj,m_chebfi,m_dfpt_cgwf,m_dft_energy,m_getghc
!!      m_gwls_hamiltonian,m_ksdiago,m_lobpcgwf_old,m_orbmag,m_rf2,m_rmm_diis
!!
!! CHILDREN
!!      getghc,mkl_set_num_threads,omp_set_nested
!!
!! SOURCE

subroutine getghc_ompgpu(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
                  prtvol,sij_opt,tim_getghc,type_calc,&
                  kg_fft_k,kg_fft_kp,select_k,cwavef_r) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,ndat, prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 integer,intent(in),optional :: select_k
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
 real(dp),intent(out),target :: gsc(:,:)
 real(dp),intent(inout),target :: cwavef(:,:)
 real(dp),optional,intent(inout) :: cwavef_r(:,:,:,:,:)
 real(dp),intent(out),target :: ghc(:,:)
 real(dp),intent(out),target :: gvnlxc(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
 !MG: Passing these arrays assumed-shape has the drawback that client code is
 !forced to use vec(2, npw*ndat) instead of the more user-friendly vec(2,npw,ndat)

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! - AMD GPUs    : FC_LLVM + HIP
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

!Local variables-------------------------------
!scalars
 integer,parameter :: level=114,re=1,im=2,tim_fourwf=1
 integer :: choice,cplex,cpopt_here,i1,i2,i3,idat,idir,ierr
 integer :: ig,igspinor,ii,iispinor,ikpt_this_proc,ipw,ispinor,my_nspinor
 integer :: n4,n5,n6,nnlout,npw_fft,npw_k1,npw_k2,nspinortot,option_fft
 integer :: paw_opt,select_k_,shift1,shift2,signs,tim_nonlop
 logical :: k1_eq_k2,has_fock,local_gvnlxc
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc,use_cwavef_r
 real(dp) :: ghcim,ghcre,weight
 character(len=500) :: msg
!arrays
 integer, pointer :: gbound_k1(:,:),gbound_k2(:,:),kg_k1(:,:),kg_k2(:,:)
 integer, ABI_CONTIGUOUS pointer :: indices_pw_fft(:),kg_k_fft(:,:)
 integer, ABI_CONTIGUOUS pointer :: recvcount_fft(:),recvdisp_fft(:)
 integer, ABI_CONTIGUOUS pointer ::  sendcount_fft(:),senddisp_fft(:)
 integer, allocatable:: dimcprj(:)
 real(dp) :: enlout(ndat),lambda_ndat(ndat),tsec(2)
 real(dp),target :: nonlop_dum(1,1)
 real(dp),allocatable :: buff_wf(:,:),cwavef1(:,:),cwavef2(:,:),cwavef_fft(:,:),cwavef_fft_tr(:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),ghc3(:,:),ghc4(:,:),ghc_mGGA(:,:),ghc_vectornd(:,:)
 real(dp),allocatable :: gvnlc(:,:),vlocal_tmp(:,:,:)
 real(dp),pointer :: gvnlxc_(:,:),kinpw_k1(:),kinpw_k2(:),kpt_k1(:),kpt_k2(:)
 real(dp),pointer :: gsc_ptr(:,:)
 type(fock_common_type),pointer :: fock
 type(pawcprj_type),pointer :: cwaveprj_fock(:,:),cwaveprj_idat(:,:),cwaveprj_nonlop(:,:)
 logical :: transfer_ghc,transfer_gsc,transfer_cwavef,transfer_gvnlxc,fourwf_on_cpu
 integer :: tmp2i(2)

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
 select_k_=1;if (present(select_k)) select_k_=select_k
 if (select_k_==KPRIME_H_K) then
!  <k^prime|H|k>
   npw_k1    =  gs_ham%npw_fft_k ; npw_k2    =  gs_ham%npw_fft_kp
   kpt_k1    => gs_ham%kpt_k     ; kpt_k2    => gs_ham%kpt_kp
   kg_k1     => gs_ham%kg_k      ; kg_k2     => gs_ham%kg_kp
   gbound_k1 => gs_ham%gbound_k  ; gbound_k2 => gs_ham%gbound_kp
   kinpw_k1  => gs_ham%kinpw_k   ; kinpw_k2  => gs_ham%kinpw_kp
 else if (select_k_==K_H_KPRIME) then
!  <k|H|k^prime>
   npw_k1    =  gs_ham%npw_fft_kp; npw_k2    =  gs_ham%npw_fft_k
   kpt_k1    => gs_ham%kpt_kp    ; kpt_k2    => gs_ham%kpt_k
   kg_k1     => gs_ham%kg_kp     ; kg_k2     => gs_ham%kg_k
   gbound_k1 => gs_ham%gbound_kp ; gbound_k2 => gs_ham%gbound_k
   kinpw_k1  => gs_ham%kinpw_kp  ; kinpw_k2  => gs_ham%kinpw_k
 else if (select_k_==K_H_K) then
!  <k|H|k>
   npw_k1    =  gs_ham%npw_fft_k ; npw_k2    =  gs_ham%npw_fft_k
   kpt_k1    => gs_ham%kpt_k     ; kpt_k2    => gs_ham%kpt_k
   kg_k1     => gs_ham%kg_k      ; kg_k2     => gs_ham%kg_k
   gbound_k1 => gs_ham%gbound_k  ; gbound_k2 => gs_ham%gbound_k
   kinpw_k1  => gs_ham%kinpw_k   ; kinpw_k2  => gs_ham%kinpw_k
 else if (select_k_==KPRIME_H_KPRIME) then
!  <k^prime|H|k^prime>
   npw_k1    =  gs_ham%npw_fft_kp; npw_k2    =  gs_ham%npw_fft_kp
   kpt_k1    => gs_ham%kpt_kp    ; kpt_k2    => gs_ham%kpt_kp
   kg_k1     => gs_ham%kg_kp     ; kg_k2     => gs_ham%kg_kp
   gbound_k1 => gs_ham%gbound_kp ; gbound_k2 => gs_ham%gbound_kp
   kinpw_k1  => gs_ham%kinpw_kp  ; kinpw_k2  => gs_ham%kinpw_kp
 end if
 k1_eq_k2=(all(abs(kpt_k1(:)-kpt_k2(:))<tol8))

!Check sizes
 my_nspinor=max(1,gs_ham%nspinor/mpi_enreg%nproc_spinor)
 if (size(cwavef)<2*npw_k1*my_nspinor*ndat) then
   ABI_BUG('wrong size for cwavef!')
 end if
 if (size(ghc)<2*npw_k2*my_nspinor*ndat) then
   ABI_BUG('wrong size for ghc!')
 end if
 if (sij_opt==1) then
   if (size(gsc)<2*npw_k2*my_nspinor*ndat) then
     ABI_BUG('wrong size for gsc!')
   end if
 end if
 if (gs_ham%usepaw==1.and.cpopt>=0) then
   if (size(cwaveprj)<gs_ham%natom*my_nspinor*ndat) then
     ABI_BUG('wrong size for cwaveprj!')
   end if
 end if
 if (any(type_calc == [0, 2, 3])) then
   local_gvnlxc = size(gvnlxc)<=1
   if (local_gvnlxc) then
     ABI_MALLOC(gvnlxc_,(2,npw_k2*my_nspinor*ndat))
   else
     gvnlxc_ => gvnlxc
   end if
   if (size(gvnlxc_)<2*npw_k2*my_nspinor*ndat) then
     ABI_BUG('wrong size for gvnlxc!')
   end if
 end if
 fourwf_on_cpu=.false.; if((mpi_enreg%me_g0==1 .or. mpi_enreg%me_g0_fft==1) .and. gs_ham%istwf_k==2) fourwf_on_cpu=.true.
 use_cwavef_r=present(cwavef_r)
 n4=gs_ham%n4 ; n5=gs_ham%n5 ; n6=gs_ham%n6
 nspinortot=gs_ham%nspinor
 if (use_cwavef_r) then
   if (size(cwavef_r,1)/=2) then
     ABI_BUG('wrong size for cwavef_r (dimension 1)')
   end if
   if (size(cwavef_r,2)/=n4) then
     ABI_BUG('wrong size for cwavef_r (dimension 2)')
   end if
   if (size(cwavef_r,3)/=n5) then
     ABI_BUG('wrong size for cwavef_r (dimension 3)')
   end if
   if (size(cwavef_r,4)/=n6) then
     ABI_BUG('wrong size for cwavef_r (dimension 4)')
   end if
   if (size(cwavef_r,5)/=nspinortot) then
     ABI_BUG('wrong size for cwavef_r (dimension 5)')
   end if
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

has_fock=.false.
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

 paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3
 if (gs_ham%usepaw==0) gsc_ptr => nonlop_dum
 if (gs_ham%usepaw==1) gsc_ptr => gsc


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

!============================================================
! Application of the local potential
!============================================================
 ABI_NVTX_START_RANGE(NVTX_GETGHC_LOCPOT)

 if (any(type_calc == [0, 1, 3])) then

   if(fourwf_on_cpu) then
     gs_ham%gpu_option=ABI_GPU_DISABLED
     !$OMP TARGET UPDATE FROM(cwavef)
   end if
!  Need a Vlocal
   if (.not.associated(gs_ham%vlocal)) then
     ABI_BUG("We need vlocal in gs_ham!")
   end if

!  fourwf can only process with one value of istwf_k
   if (.not.k1_eq_k2) then
     ABI_BUG('vlocal (fourwf) cannot be computed with k/=k^prime!')
   end if

!  Apply the local potential to the wavefunction
!  Start from wavefunction in reciprocal space cwavef
   if(buf_initialized==0 .or. mod__n4/=gs_ham%n4 .or. mod__n5/=gs_ham%n5 &
  &   .or. mod__n6/=gs_ham%n6 .or. mod__nspinor/=my_nspinor .or. mod__ndat/=ndat .or. npw_k1/=mod__npw) then
      call alloc_getghc_ompgpu_buffers(npw_k1, my_nspinor, ndat, gs_ham%n4, gs_ham%n5, gs_ham%n6)
   end if
!  End with function ghc in reciprocal space also.
#ifndef HAVE_GPU_HIP
   !Work buffer allocated at each call to save memory (but too costful on HIP)
   !$OMP TARGET ENTER DATA MAP(alloc:work)
#endif
   weight=one
   if (.not.use_cwavef_r) then
     option_fft=2
     if (nspinortot==2) then
       ABI_MALLOC(cwavef1,(2,npw_k1*ndat))
       ABI_MALLOC(cwavef2,(2,npw_k1*ndat))
       !$OMP TARGET ENTER DATA MAP(alloc:cwavef1,cwavef2)
       !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(idat) MAP(to:cwavef1,cwavef2,cwavef)
       do idat=1,ndat
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ipw,ig)
         do ipw=1,npw_k1
           do ig=1,im
             cwavef1(ig,ipw+(idat-1)*npw_k1)=cwavef(ig,ipw+(idat-1)*my_nspinor*npw_k1)
             cwavef2(ig,ipw+(idat-1)*npw_k1)=cwavef(ig,ipw+(idat-1)*my_nspinor*npw_k1+shift1)
           end do
         end do
       end do
!      call cg_zcopy(npw_k1*ndat,cwavef(1,1),cwavef1)
!      call cg_zcopy(npw_k1*ndat,cwavef(1,1+shift1),cwavef2)
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
       !$OMP TASKWAIT
       call fourwf(1,gs_ham%vlocal,cwavef,ghc,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
&       weight,weight,gpu_option=gs_ham%gpu_option)

     else
       ! nspinortot==2

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
           !$OMP TARGET UPDATE TO(work)
         end if
         ABI_MALLOC(ghc1,(2,npw_k2*ndat))
         !$OMP TARGET ENTER DATA MAP(alloc:ghc1)
         call fourwf(1,gs_ham%vlocal,cwavef1,ghc1,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
&         weight,weight,gpu_option=gs_ham%gpu_option)
         !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(idat) MAP(to:ghc,ghc1)
         do idat=1,ndat
           !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ipw,ig)
           do ipw =1, npw_k2
             do ig =1, im
               ghc(ig,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(ig,ipw+(idat-1)*npw_k2)
             end do
           end do
         end do
         !$OMP TARGET EXIT DATA MAP(delete:ghc1)
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
           !$OMP TARGET UPDATE TO(work)
         end if
         ABI_MALLOC(ghc2,(2,npw_k2*ndat))
         !$OMP TARGET ENTER DATA MAP(alloc:ghc2)
         call fourwf(1,gs_ham%vlocal,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&         gpu_option=gs_ham%gpu_option)
         !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(idat) MAP(to:ghc,ghc2)
         do idat=1,ndat
           !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ipw,ig)
           do ipw=1,npw_k2
             do ig =1, im
               ghc(ig,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc2(ig,ipw+(idat-1)*npw_k2)
             end do
           end do
         end do
         !$OMP TARGET EXIT DATA MAP(delete:ghc2)
         ABI_FREE(ghc2)
       end if ! spin 2 treated by this proc

     end if ! nspinortot

   else if (gs_ham%nvloc==4) then
!    Treat non-collinear local potentials

     ABI_MALLOC(ghc1,(2,npw_k2*ndat))
     ABI_MALLOC(ghc2,(2,npw_k2*ndat))
     ABI_MALLOC(ghc3,(2,npw_k2*ndat))
     ABI_MALLOC(ghc4,(2,npw_k2*ndat))
     !$OMP TARGET ENTER DATA MAP(alloc:ghc1,ghc2,ghc3,ghc4)
     call gpu_set_to_zero(ghc1, int(2,c_size_t)*npw_k2*ndat)
     call gpu_set_to_zero(ghc2, int(2,c_size_t)*npw_k2*ndat)
     call gpu_set_to_zero(ghc3, int(2,c_size_t)*npw_k2*ndat)
     call gpu_set_to_zero(ghc4, int(2,c_size_t)*npw_k2*ndat)

     if (use_cwavef_r) then
       ABI_MALLOC(vlocal_tmp,(0,0,0))
     else
       ABI_MALLOC(vlocal_tmp,(gs_ham%n4,gs_ham%n5,gs_ham%n6))
     end if
!    ghc1=v11*phi1
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
         !$OMP TARGET UPDATE TO(work)
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
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&       gpu_option=gs_ham%gpu_option)
     end if
!    ghc2=v22*phi2
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
         !$OMP TARGET UPDATE TO(work)
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
         !$OMP TARGET UPDATE TO(work)
       end if
       call fourwf(1,vlocal_tmp,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&       gpu_option=gs_ham%gpu_option)
     end if
     ABI_FREE(vlocal_tmp)
     cplex=2
     if (use_cwavef_r) then
       ABI_MALLOC(vlocal_tmp,(0,0,0))
     else
       ABI_MALLOC(vlocal_tmp,(cplex*gs_ham%n4,gs_ham%n5,gs_ham%n6))
     end if
!    ghc3=(re(v12)-im(v12))*phi1
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
         !$OMP TARGET UPDATE TO(work)
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
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&       gpu_option=gs_ham%gpu_option)
     end if
!    ghc4=(re(v12)+im(v12))*phi2
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
         !$OMP TARGET UPDATE TO(work)
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
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&       gpu_option=gs_ham%gpu_option)
     end if
     ABI_FREE(vlocal_tmp)
!    Build ghc from pieces
!    (v11,v22,Re(v12)+iIm(v12);Re(v12)-iIm(v12))(psi1;psi2): matrix product
     if (mpi_enreg%paral_spinor==0) then
       !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc1,ghc2,ghc3,ghc4) PRIVATE(idat)
       do idat=1,ndat
         !$OMP PARALLEL DO PRIVATE(ipw)
         do ipw=1,npw_k2
           ghc(1,ipw+(idat-1)*my_nspinor*npw_k2)       =ghc1(1,ipw+(idat-1)*npw_k2)+ghc4(1,ipw+(idat-1)*npw_k2)
           ghc(2,ipw+(idat-1)*my_nspinor*npw_k2)       =ghc1(2,ipw+(idat-1)*npw_k2)+ghc4(2,ipw+(idat-1)*npw_k2)
           ghc(1,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1,ipw+(idat-1)*npw_k2)+ghc2(1,ipw+(idat-1)*npw_k2)
           ghc(2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(2,ipw+(idat-1)*npw_k2)+ghc2(2,ipw+(idat-1)*npw_k2)
         end do
       end do
     else
       call xmpi_sum(ghc4,mpi_enreg%comm_spinor,ierr,use_omp_map=.true.)
       call xmpi_sum(ghc3,mpi_enreg%comm_spinor,ierr,use_omp_map=.true.)
       if (nspinor1TreatedByThisProc) then
         !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc1,ghc4) PRIVATE(idat)
         do idat=1,ndat
           !$OMP PARALLEL DO PRIVATE(ipw)
           do ipw=1,npw_k2
             ghc(1,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(1,ipw+(idat-1)*npw_k2)+ghc4(1,ipw+(idat-1)*npw_k2)
             ghc(2,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(2,ipw+(idat-1)*npw_k2)+ghc4(2,ipw+(idat-1)*npw_k2)
           end do
         end do
       else if (nspinor2TreatedByThisProc) then
         !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ghc,ghc2,ghc3) PRIVATE(idat)
         do idat=1,ndat
           !$OMP PARALLEL DO PRIVATE(ipw)
           do ipw=1,npw_k2
             ghc(1,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1,ipw+(idat-1)*npw_k2)+ghc2(1,ipw+(idat-1)*npw_k2)
             ghc(2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(2,ipw+(idat-1)*npw_k2)+ghc2(2,ipw+(idat-1)*npw_k2)
           end do
         end do
       end if
     end if
     !$OMP TARGET EXIT DATA MAP(delete:ghc1,ghc2,ghc3,ghc4)
     ABI_FREE(ghc1)
     ABI_FREE(ghc2)
     ABI_FREE(ghc3)
     ABI_FREE(ghc4)
   end if ! nvloc

   if (nspinortot==2)  then
     !$OMP TARGET EXIT DATA MAP(delete:cwavef1,cwavef2) IF (.not.use_cwavef_r)
     ABI_FREE(cwavef1)
     ABI_FREE(cwavef2)
   end if

!  Add metaGGA contribution
   if (associated(gs_ham%vxctaulocal)) then
     if (.not.k1_eq_k2) then
       ABI_BUG('metaGGA not allowed for k/=k_^prime!')
     end if
     if (size(gs_ham%vxctaulocal)/=gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*4) then
       ABI_BUG('wrong sizes for vxctaulocal!')
     end if
     ABI_MALLOC(ghc_mGGA,(2,npw_k2*my_nspinor*ndat))
     call getghc_mGGA_ompgpu(cwavef,ghc_mGGA,gbound_k1,gs_ham%gprimd,gs_ham%istwf_k,kg_k1,kpt_k1,&
     &    gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
     &    gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vxctaulocal,gs_ham%gpu_option)
     !$OMP TARGET UPDATE FROM(ghc)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_mGGA(1:2,1:npw_k2*my_nspinor*ndat)
     !$OMP TARGET UPDATE TO(ghc)
     ABI_FREE(ghc_mGGA)
   end if

   !  Add nuclear dipole moment contribution
   if (associated(gs_ham%vectornd)) then
     if (.not.k1_eq_k2) then
       ABI_BUG('nuclear dipole vector potential not allowed for k/=k_^prime!')
     end if
     if (size(gs_ham%vectornd)/=gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*3) then
       ABI_BUG('wrong sizes for vectornd in getghc!')
     end if
     ABI_MALLOC(ghc_vectornd,(2,npw_k2*my_nspinor*ndat))

     call getghc_nucdip_ompgpu(cwavef,ghc_vectornd,gbound_k1,gs_ham%istwf_k,kg_k1,kpt_k1,&
     &    gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
     &    gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vectornd,gs_ham%gpu_option)

     !$OMP TARGET UPDATE FROM(ghc)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_vectornd(1:2,1:npw_k2*my_nspinor*ndat)
     !$OMP TARGET UPDATE TO(ghc)

     ABI_FREE(ghc_vectornd)
   end if

   if (type_calc==1 .and. .not. fourwf_on_cpu) then
     !$OMP TARGET UPDATE FROM(ghc) nowait IF(transfer_ghc)
   end if

#ifndef HAVE_GPU_HIP
   !$OMP TARGET EXIT DATA MAP(delete:work)
#endif
   if(fourwf_on_cpu) then
     gs_ham%gpu_option=ABI_GPU_OPENMP
     !$OMP TARGET UPDATE TO(ghc)
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
           call pawcprj_getdim(dimcprj,gs_ham%natom,gs_ham%nattyp,gs_ham%ntypat,&
&           gs_ham%typat,fock%pawtab,'O')
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

     call nonlop(choice,cpopt_here,cwaveprj_nonlop,enlout,gs_ham,idir,lambda_ndat,mpi_enreg,ndat,&
&     nnlout,paw_opt,signs,gsc_ptr,tim_nonlop,cwavef,gvnlxc_,select_k=select_k_)

     if (gs_ham%usepaw==1 .and. has_fock)then
       if (fock_get_getghc_call(fock)==1) then
         ABI_MALLOC(gvnlc, (2,npw_k2*my_nspinor*ndat))
         !$OMP TARGET UPDATE FROM(gvnlxc_)
         gvnlc=gvnlxc_
         !$OMP TARGET ENTER DATA MAP(to:gvnlc)
       endif
     endif

!    Calculation of the Fock exact exchange contribution from the Fock or ACE operator
     if (has_fock) then
       if (fock_get_getghc_call(fock)==1) then
         if (gs_ham%usepaw==0) cwaveprj_idat => cwaveprj
         if (fock%use_ACE==0) then
           !$OMP TARGET UPDATE FROM(cwavef,gvnlxc_)
           call timab(360,1,tsec)
           call fock_getghc(cwavef,cwaveprj_idat,gvnlxc_,gs_ham,mpi_enreg,ndat)
           call timab(360,2,tsec)
           !$OMP TARGET UPDATE TO(cwavef,gvnlxc_)
         else
           call fock_ACE_getghc(cwavef,gvnlxc_,gs_ham,mpi_enreg,ndat,gpu_option=ABI_GPU_OPENMP)
         end if
       end if
     end if

   else if (type_calc == 3) then
     ! for kinetic and local only, nonlocal and vfock should be zero
     call gpu_set_to_zero(gvnlxc_, int(2,c_size_t)*npw_k2*my_nspinor*ndat)
   end if ! if(type_calc...

   ABI_NVTX_END_RANGE()

!============================================================
! Assemble kinetic, local, nonlocal and Fock contributions
!============================================================

   ABI_NVTX_START_RANGE(NVTX_GETGHC_KIN)
!  Assemble modified kinetic, local and nonlocal contributions
   !  to <G|H|C(n,k)>. Take also into account build-in debugging.
   if(prtvol/=-level)then
     if (k1_eq_k2) then
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& PRIVATE(idat,ispinor) &
       !$OMP& MAP(to:ghc,kinpw_k2,gvnlxc_,gsc,cwavef) MAP(tofrom:kinpw_k2)
       do idat=1,ndat
         do ispinor=1,my_nspinor
           !$OMP PARALLEL DO PRIVATE(ig,igspinor)
           do ig=1,npw_k2
             igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
             if(kinpw_k2(ig)<huge(zero)*1.d-11)then
               ghc(re,igspinor)= ghc(re,igspinor) + kinpw_k2(ig)*cwavef(re,igspinor) + gvnlxc_(re,igspinor)
               ghc(im,igspinor)= ghc(im,igspinor) + kinpw_k2(ig)*cwavef(im,igspinor) + gvnlxc_(im,igspinor)
             else
               ghc(re,igspinor)=zero
               ghc(im,igspinor)=zero
               if (sij_opt==1) then
                 gsc(re,igspinor)=zero
                 gsc(im,igspinor)=zero
               end if
             end if
           end do ! ig
         end do ! ispinor
       end do ! idat
     else
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& PRIVATE(idat,ispinor) &
       !$OMP& MAP(to:ghc,gvnlxc_,gsc) MAP(tofrom:kinpw_k2)
       do idat=1,ndat
         do ispinor=1,my_nspinor
           !$OMP PARALLEL DO PRIVATE(ig,igspinor)
           do ig=1,npw_k2
             igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
             if(kinpw_k2(ig)<huge(zero)*1.d-11)then
               ghc(re,igspinor)= ghc(re,igspinor) + gvnlxc_(re,igspinor)
               ghc(im,igspinor)= ghc(im,igspinor) + gvnlxc_(im,igspinor)
             else
               ghc(re,igspinor)=zero
               ghc(im,igspinor)=zero
               if (sij_opt==1) then
                 gsc(re,igspinor)=zero
                 gsc(im,igspinor)=zero
               end if
             end if
           end do ! ig
         end do ! ispinor
       end do ! idat
     end if
   else
     !$OMP TARGET UPDATE FROM(ghc)
     !$OMP TARGET UPDATE FROM(gsc_ptr)
!    Here, debugging section
     call wrtout(std_out,' getghc : components of ghc ','PERS')
     write(msg,'(a)')&
&     'icp ig ispinor igspinor re/im     ghc        kinpw         cwavef      glocc        gvnlxc  gsc'
     call wrtout(std_out,msg,'PERS')
     do idat=1,ndat
       do ispinor=1,my_nspinor
         do ig=1,npw_k2
           igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
           if(kinpw_k2(ig)<huge(zero)*1.d-11)then
             if (k1_eq_k2) then
               ghcre=kinpw_k2(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlxc_(re,igspinor)
               ghcim=kinpw_k2(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlxc_(im,igspinor)
             else
               ghcre=ghc(re,igspinor)+gvnlxc_(re,igspinor)
               ghcim=ghc(im,igspinor)+gvnlxc_(im,igspinor)
             end if
           else
             ghcre=zero
             ghcim=zero
             if (sij_opt==1) then
               gsc(re,igspinor)=zero
               gsc(im,igspinor)=zero
             end if
           end if
           iispinor=ispinor;if (mpi_enreg%paral_spinor==1) iispinor=mpi_enreg%me_spinor+1
           if (sij_opt == 1) then
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw_k2(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlxc_(re,igspinor), gsc(re,igspinor)
             call wrtout(std_out,msg,'PERS')
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw_k2(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlxc_(im,igspinor), gsc(im,igspinor)
             call wrtout(std_out,msg,'PERS')
           else
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw_k2(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlxc_(re,igspinor)
             call wrtout(std_out,msg,'PERS')
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw_k2(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlxc_(im,igspinor)
             call wrtout(std_out,msg,'PERS')
           end if
           ghc(re,igspinor)=ghcre
           ghc(im,igspinor)=ghcim
         end do ! ig
       end do ! ispinor
     end do ! idat
   end if
   ABI_NVTX_END_RANGE()

!  Special case of PAW + Fock : only return Fock operator contribution in gvnlxc_
   if (gs_ham%usepaw==1 .and. has_fock) then
     !$OMP TARGET UPDATE FROM(gvnlxc_,gvnlc)
     gvnlxc_=gvnlxc_-gvnlc
     !$OMP TARGET UPDATE TO(gvnlxc_) if(.not. local_gvnlxc)
     !$OMP TARGET EXIT DATA MAP(delete:gvnlc)
     ABI_FREE(gvnlc)
   endif

   !$OMP TARGET UPDATE FROM(ghc)     nowait IF(transfer_ghc)
   !$OMP TARGET UPDATE FROM(gsc)     nowait IF(transfer_gsc)
   !$OMP TARGET UPDATE FROM(cwavef)  nowait IF(transfer_cwavef)
   if (.not. local_gvnlxc) then
     !$OMP TARGET UPDATE FROM(gvnlxc_) nowait IF(transfer_gvnlxc)
   end if
   !$OMP TASKWAIT

!  Structured debugging : if prtvol=-level, stop here.
   if(prtvol==-level)then
     write(msg,'(a,i0,a)')' getghc : exit prtvol=-',level,', debugging mode => stop '
     ABI_ERROR(msg)
   end if

 end if ! type_calc

 !$OMP TARGET EXIT DATA MAP(delete:ghc)     IF(transfer_ghc)
 !$OMP TARGET EXIT DATA MAP(delete:gsc)     IF(transfer_gsc)
 !$OMP TARGET EXIT DATA MAP(delete:cwavef)  IF(transfer_cwavef)
 !$OMP TARGET EXIT DATA MAP(delete:gvnlxc_) IF(transfer_gvnlxc)
 if (local_gvnlxc .and. any(type_calc == [0, 2, 3])) then
   ABI_FREE(gvnlxc_)
 end if

 call timab(350+tim_getghc,2,tsec)

 DBG_EXIT("COLL")

 ABI_NVTX_END_RANGE()
#else

 ABI_UNUSED((/cpopt,ndat,prtvol,sij_opt,tim_getghc/))
 ABI_UNUSED((/kg_fft_k,kg_fft_kp,type_calc,select_k/))
 ABI_UNUSED((/gsc,cwavef,cwavef_r,ghc,gvnlxc,lambda/))
 ABI_UNUSED_A(mpi_enreg)
 ABI_UNUSED_A(cwaveprj)
 ABI_UNUSED_A(gs_ham)
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")

#endif
end subroutine getghc_ompgpu
!!***

!!****f* ABINIT/getghc_nucdip_ompgpu
!!
!! NAME
!! getghc_nucdip_ompgpu
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
!! n4,n5,n6=for dimensionning of vxctaulocal
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
!! this code is a copied, simplied version of getghc_mGGA (see below) and should eventually be
!! integrated into that code, to simplify maintenance
!!
!! SOURCE

subroutine getghc_nucdip_ompgpu(cwavef,ghc_vectornd,gbound_k,istwf_k,kg_k,kpt,mgfft,mpi_enreg,&
&                      ndat,ngfft,npw_k,nvloc,n4,n5,n6,my_nspinor,vectornd,gpu_option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,my_nspinor,ndat,npw_k,nvloc,n4,n5,n6,gpu_option
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: gbound_k(2*mgfft+4),kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(inout) :: cwavef(2,npw_k*my_nspinor*ndat)
 real(dp),intent(inout) :: ghc_vectornd(2,npw_k*my_nspinor*ndat)
 real(dp),intent(inout) :: vectornd(n4,n5,n6,nvloc,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf=1
 integer :: icmplx,idat,idir,ipw,iv1,iv2,nspinortot,shift
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: scale_conversion,weight=one
 !arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: gcwavef(:,:,:),gcwavef1(:,:,:),gcwavef2(:,:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),kgkpk(:,:)
 real(dp),allocatable :: dx(:),dy(:)

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

 ! scale conversion from SI to atomic units,
 ! here \alpha^2 where \alpha is the fine structure constant
 scale_conversion = FineStructureConstant2

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
    ABI_MALLOC(dx,(npw_k))
    ABI_MALLOC(dy,(npw_k))
    do idir=1,3
      call fourwf(1,vectornd(:,:,:,:,idir),gcwavef(:,:,idir),ghc1,work,gbound_k,gbound_k,&
           istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
           &     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
       ! DAXPY is a BLAS routine for y -> A*x + y, here x = ghc1, A = scale_conversion, and y = ghc_vectornd
       ! should be faster than explicit loop over ipw as npw_k gets large
      do idat=1,ndat
        iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
        do icmplx=1,2
          dx=ghc1(icmplx,iv1:iv2)
          dy=ghc_vectornd(icmplx,iv1:iv2)
          call DAXPY(npw_k,scale_conversion,dx,1,dy,1)
          ghc_vectornd(icmplx,iv1:iv2)=dy
        end do
      end do
    end do ! idir
    ABI_FREE(dx)
    ABI_FREE(dy)
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
       ABI_MALLOC(dx,(npw_k))
       ABI_MALLOC(dy,(npw_k))
       do idir=1,3
          call fourwf(1,vectornd(:,:,:,:,idir),gcwavef1(:,:,idir),ghc1,work,gbound_k,gbound_k,&
               istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
               &     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
          ! DAXPY is a BLAS routine for y -> A*x + y, here x = ghc1, A = scale_conversion, and y = ghc_vectornd
          ! should be faster than explicit loop over ipw as npw_k gets large
          do idat=1,ndat
            iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
            do icmplx=1,2
              dx=ghc1(icmplx,iv1:iv2)
              dy=ghc_vectornd(icmplx,iv1:iv2)
              call DAXPY(npw_k,scale_conversion,dx,1,dy,1)
              ghc_vectornd(icmplx,iv1:iv2)=dy
            end do
          end do
       end do ! idir
       ABI_FREE(dx)
       ABI_FREE(dy)
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
       ABI_MALLOC(dx,(npw_k))
       ABI_MALLOC(dy,(npw_k))
       do idir=1,3
          call fourwf(1,vectornd(:,:,:,:,idir),gcwavef2(:,:,idir),ghc2,work,gbound_k,gbound_k,&
               istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
               &     tim_fourwf,weight,weight,gpu_option=gpu_option)
!!$OMP PARALLEL DO
          ! DAXPY is a BLAS routine for y -> A*x + y, here x = ghc1, A = scale_conversion, and y = ghc_vectornd
          ! should be faster than explicit loop over ipw as npw_k gets large
          do idat=1,ndat
            iv1=1+(idat-1)*npw_k; iv2=-1+iv1+npw_k
            do icmplx=1,2
              dx=ghc2(icmplx,iv1:iv2)
              dy=ghc_vectornd(icmplx,iv1+shift:iv2+shift)
              call DAXPY(npw_k,scale_conversion,dx,1,dy,1)
              ghc_vectornd(icmplx,iv1+shift:iv2+shift)=dy
            end do
          end do
       end do ! idir
       ABI_FREE(dx)
       ABI_FREE(dy)
       ABI_FREE(gcwavef2)
       ABI_FREE(ghc2)

    end if ! end spinor 2

    ABI_FREE(cwavef1)
    ABI_FREE(cwavef2)
    ABI_FREE(kgkpk)

 end if ! nspinortot

end subroutine getghc_nucdip_ompgpu
!!***

!!****f* ABINIT/getghc_mGGA_ompgpu
!!
!! NAME
!! getghc_mGGA_ompgpu
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

subroutine getghc_mGGA_ompgpu(cwavef,ghc_mGGA,gbound_k,gprimd,istwf_k,kg_k,kpt,mgfft,mpi_enreg,&
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
 real(dp) :: gp2pi1,gp2pi2,gp2pi3,kpt_cart,kg_k_cart,weight=one
!arrays
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
   do idat=1,ndat
     do ipw=1,npw_k
       gcwavef(:,ipw+(idat-1)*npw_k,1:3)=zero
       lcwavef(:,ipw+(idat-1)*npw_k)  =zero
     end do
   end do
   do idir=1,3
     gp2pi1=gprimd(idir,1)*two_pi
     gp2pi2=gprimd(idir,2)*two_pi
     gp2pi3=gprimd(idir,3)*two_pi
     kpt_cart=gp2pi1*kpt(1)+gp2pi2*kpt(2)+gp2pi3*kpt(3)
!    Multiplication by 2pi i (G+k)_idir for gradient
!    Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
     do idat=1,ndat
       do ipw=1,npw_k
         kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
         gcwavef(1,ipw+(idat-1)*npw_k,idir)= cwavef(2,ipw+(idat-1)*npw_k)*kg_k_cart
         gcwavef(2,ipw+(idat-1)*npw_k,idir)=-cwavef(1,ipw+(idat-1)*npw_k)*kg_k_cart
         lcwavef(1,ipw+(idat-1)*npw_k)=lcwavef(1,ipw+(idat-1)*npw_k)-cwavef(1,ipw+(idat-1)*npw_k)*kg_k_cart**2
         lcwavef(2,ipw+(idat-1)*npw_k)=lcwavef(2,ipw+(idat-1)*npw_k)-cwavef(2,ipw+(idat-1)*npw_k)*kg_k_cart**2
       end do
     end do
   end do ! idir
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
!!$OMP PARALLEL DO
     do idat=1,ndat
       do ipw=1,npw_k
         gcwavef1(:,ipw+(idat-1)*npw_k,1:3)=zero
         lcwavef1(:,ipw+(idat-1)*npw_k)=zero
       end do
     end do
     do idir=1,3
       gp2pi1=gprimd(idir,1)*two_pi
       gp2pi2=gprimd(idir,2)*two_pi
       gp2pi3=gprimd(idir,3)*two_pi
       kpt_cart=gp2pi1*kpt(1)+gp2pi2*kpt(2)+gp2pi3*kpt(3)
!      Multiplication by 2pi i (G+k)_idir for gradient
!      Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
       do idat=1,ndat
         do ipw=1,npw_k
           kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
           gcwavef1(1,ipw+(idat-1)*npw_k,idir)= cwavef1(2,ipw+(idat-1)*npw_k)*kg_k_cart
           gcwavef1(2,ipw+(idat-1)*npw_k,idir)=-cwavef1(1,ipw+(idat-1)*npw_k)*kg_k_cart
           lcwavef1(1,ipw+(idat-1)*npw_k)=lcwavef1(1,ipw+(idat-1)*npw_k)-cwavef1(1,ipw+(idat-1)*npw_k)*kg_k_cart**2
           lcwavef1(2,ipw+(idat-1)*npw_k)=lcwavef1(2,ipw+(idat-1)*npw_k)-cwavef1(2,ipw+(idat-1)*npw_k)*kg_k_cart**2
         end do
       end do
     end do ! idir
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
     do idat=1,ndat
       do ipw=1,npw_k
         gcwavef2(:,ipw+(idat-1)*npw_k,1:3)=zero
         lcwavef2(:,ipw+(idat-1)*npw_k)  =zero
       end do
     end do
     do idir=1,3
       gp2pi1=gprimd(idir,1)*two_pi
       gp2pi2=gprimd(idir,2)*two_pi
       gp2pi3=gprimd(idir,3)*two_pi
       kpt_cart=gp2pi1*kpt(1)+gp2pi2*kpt(2)+gp2pi3*kpt(3)
!      Multiplication by 2pi i (G+k)_idir for gradient
!      Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
       do idat=1,ndat
         do ipw=1,npw_k
           kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
           gcwavef2(1,ipw+(idat-1)*npw_k,idir)= cwavef2(2,ipw+(idat-1)*npw_k)*kg_k_cart
           gcwavef2(2,ipw+(idat-1)*npw_k,idir)=-cwavef2(1,ipw+(idat-1)*npw_k)*kg_k_cart
           lcwavef2(1,ipw+(idat-1)*npw_k)=lcwavef2(1,ipw+(idat-1)*npw_k)-cwavef2(1,ipw+(idat-1)*npw_k)*kg_k_cart**2
           lcwavef2(2,ipw+(idat-1)*npw_k)=lcwavef2(2,ipw+(idat-1)*npw_k)-cwavef2(2,ipw+(idat-1)*npw_k)*kg_k_cart**2
         end do
       end do
     end do ! idir
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

end subroutine getghc_mGGA_ompgpu
!!***

end module m_getghc_ompgpu
