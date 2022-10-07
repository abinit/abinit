!!****m* ABINIT/m_getghc
!! NAME
!!  m_getghc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space;
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2022 ABINIT group (DCA, XG, GMR, LSI, MT, JB, JWZ)
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

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_getdim, pawcprj_copy
 use m_bandfft_kpt, only : bandfft_kpt, bandfft_kpt_get_ikpt
 use m_hamiltonian, only : gs_hamiltonian_type, KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME
 use m_fock,        only : fock_common_type, fock_get_getghc_call
 use m_fock_getghc, only : fock_getghc, fock_ACE_getghc
 use m_nonlop,      only : nonlop
 use m_fft,         only : fourwf

#if defined HAVE_YAKL
 use gator_mod
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_GPU_NVTX_V3)
 use m_nvtx_data
#endif

 implicit none

 private
!!***

 public :: getghc      ! Compute <G|H|C> for input vector |C> expressed in reciprocal space
 public :: getgsc      ! Compute <G|S|C> for all input vectors |Cnk> at a given k-point
 public :: getghc_mGGA
 public :: multithreaded_getghc
 public :: getghc_nucdip ! compute <G|H_nucdip|C> for input vector |C> expressed in recip space
!!***

contains
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
!! SOURCE

subroutine getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
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
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),optional,intent(inout) :: cwavef_r(:,:,:,:,:)
 real(dp),intent(out) :: ghc(:,:)
 real(dp),intent(out),target :: gvnlxc(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
 !MG: Passing these arrays assumed-shape has the drawback that client code is
 !forced to use vec(2, npw*ndat) instead of the more user-friendly vec(2,npw,ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=114,re=1,im=2,tim_fourwf=1
 integer :: choice,cplex,cpopt_here,i1,i2,i3,idat,idir,ierr
 integer :: ig,igspinor,ii,iispinor,ikpt_this_proc,ipw,ispinor,my_nspinor
 integer :: n4,n5,n6,nnlout,npw_fft,npw_k1,npw_k2,nspinortot,option_fft
 integer :: paw_opt,select_k_,shift1,shift2,signs,tim_nonlop
 logical :: k1_eq_k2,have_to_reequilibrate,has_fock,local_gvnlxc
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc,use_cwavef_r
 real(dp) :: ghcim,ghcre,weight
 character(len=500) :: msg
!arrays
 integer,  pointer                :: gbound_k1(:,:)
 integer,  pointer                :: gbound_k2(:,:)
 integer,  pointer                :: kg_k1(:,:)
 integer,  pointer                :: kg_k2(:,:)
 integer,  ABI_CONTIGUOUS pointer :: indices_pw_fft(:), kg_k_fft(:,:)
 integer,  ABI_CONTIGUOUS pointer :: recvcount_fft(:), recvdisp_fft(:)
 integer,  ABI_CONTIGUOUS pointer :: sendcount_fft(:), senddisp_fft(:)
 integer,  allocatable:: dimcprj(:)
 real(dp)                         :: enlout(ndat), lambda_ndat(ndat), tsec(2)
 real(dp), target                 :: nonlop_dum(1,1)
 real(dp), allocatable            :: buff_wf(:,:)
 real(dp), allocatable            :: cwavef1(:,:)
 real(dp), allocatable            :: cwavef2(:,:)
 real(dp), allocatable            :: cwavef_fft(:,:)
 real(dp), allocatable            :: cwavef_fft_tr(:,:)

 real(dp), allocatable            :: ghc1(:,:)
 real(dp), allocatable            :: ghc2(:,:)
 real(dp), allocatable            :: ghc3(:,:)
 real(dp), allocatable            :: ghc4(:,:)
 real(dp), allocatable            :: ghc_mGGA(:,:)
 real(dp), allocatable            :: ghc_vectornd(:,:)

#if defined HAVE_GPU && defined HAVE_YAKL
 real(c_double), ABI_CONTIGUOUS pointer :: gvnlc(:,:)
#else
 real(dp), allocatable            :: gvnlc(:,:)
#endif

 real(dp), allocatable            :: vlocal_tmp(:,:,:)
 real(dp), allocatable            :: work(:,:,:,:)

#if defined HAVE_GPU && defined HAVE_YAKL
 real(c_double), ABI_CONTIGUOUS pointer :: gvnlxc_(:,:)
#else
 real(dp), pointer                :: gvnlxc_(:,:)
#endif

 real(dp), pointer                :: kinpw_k1(:)
 real(dp), pointer                :: kinpw_k2(:)
 real(dp), pointer                :: kpt_k1(:)
 real(dp), pointer                :: kpt_k2(:)
 real(dp), pointer                :: gsc_ptr(:,:)
 type(fock_common_type),pointer :: fock
 type(pawcprj_type),pointer :: cwaveprj_fock(:,:),cwaveprj_idat(:,:),cwaveprj_nonlop(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in getghc:
 call timab(350+tim_getghc,1,tsec)

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
   local_gvnlxc = size(gvnlxc)==0
   if (local_gvnlxc) then
#if defined HAVE_GPU && defined HAVE_YAKL
     ABI_MALLOC_MANAGED(gvnlxc_, (/2,npw_k2*my_nspinor*ndat/))
#else
     ABI_MALLOC(gvnlxc_,(2,npw_k2*my_nspinor*ndat))
#endif
   else
     gvnlxc_ => gvnlxc
   end if
   if (size(gvnlxc_)<2*npw_k2*my_nspinor*ndat) then
     ABI_BUG('wrong size for gvnlxc!')
   end if
 end if
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

!============================================================
! Application of the local potential
!============================================================
 ABI_NVTX_START_RANGE(NVTX_GETGHC_LOCPOT)

 if (any(type_calc == [0, 1, 3])) then

!  Need a Vlocal
   if (.not.associated(gs_ham%vlocal)) then
     ABI_BUG("We need vlocal in gs_ham!")
   end if

!  fourwf can only process with one value of istwf_k
   if (.not.k1_eq_k2) then
     ABI_BUG('vlocal (fourwf) cannot be computed with k/=k^prime!')
   end if

!  Eventually adjust load balancing for FFT (by changing FFT distrib)
   have_to_reequilibrate=.false.
   if (mpi_enreg%paral_kgb==1) then
     ikpt_this_proc=bandfft_kpt_get_ikpt()
     have_to_reequilibrate=bandfft_kpt(ikpt_this_proc)%have_to_reequilibrate
   end if
   if (have_to_reequilibrate) then
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
&       cwavef_fft_tr,2*ndat*recvcount_fft, 2*ndat*recvdisp_fft, mpi_enreg%comm_fft,ierr)
!      We need to transpose data
       do idat=1,ndat
         do ipw = 1 ,npw_fft
           cwavef_fft(1:2,  ipw + npw_fft*(idat-1)) = cwavef_fft_tr(1:2,  idat + ndat*(ipw-1))
         end do
       end do
     else
       call xmpi_alltoallv(buff_wf,2*sendcount_fft,2*senddisp_fft,  &
&       cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, mpi_enreg%comm_fft,ierr)
     end if
   end if

!  Apply the local potential to the wavefunction
!  Start from wavefunction in reciprocal space cwavef
!  End with function ghc in reciprocal space also.
   ABI_MALLOC(work,(2,gs_ham%n4,gs_ham%n5,gs_ham%n6*ndat))
   weight=one
   if (.not.use_cwavef_r) then
     option_fft=2
     if (nspinortot==2) then
       ABI_MALLOC(cwavef1,(2,npw_k1*ndat))
       ABI_MALLOC(cwavef2,(2,npw_k1*ndat))
       do idat=1,ndat
         do ipw=1,npw_k1
           cwavef1(1:2,ipw+(idat-1)*npw_k1)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k1)
           cwavef2(1:2,ipw+(idat-1)*npw_k1)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k1+shift1)
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
       if (have_to_reequilibrate) then
         call fourwf(1,gs_ham%vlocal,cwavef_fft,cwavef_fft,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k_fft,kg_k_fft,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_fft,npw_fft,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_ham%use_gpu_cuda)
       else
         call fourwf(1,gs_ham%vlocal,cwavef,ghc,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_ham%use_gpu_cuda)
       end if

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
         end if
         ABI_MALLOC(ghc1,(2,npw_k2*ndat))
         call fourwf(1,gs_ham%vlocal,cwavef1,ghc1,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw =1, npw_k2
             ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(1:2,ipw+(idat-1)*npw_k2)
           end do
         end do
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
         call fourwf(1,gs_ham%vlocal,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw=1,npw_k2
             ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc2(1:2,ipw+(idat-1)*npw_k2)
           end do
         end do
         ABI_FREE(ghc2)
       end if ! spin 2 treated by this proc

     end if ! nspinortot

   else if (gs_ham%nvloc==4) then
!    Treat non-collinear local potentials

     ABI_MALLOC(ghc1,(2,npw_k2*ndat))
     ABI_MALLOC(ghc2,(2,npw_k2*ndat))
     ABI_MALLOC(ghc3,(2,npw_k2*ndat))
     ABI_MALLOC(ghc4,(2,npw_k2*ndat))
     ghc1(:,:)=zero; ghc2(:,:)=zero; ghc3(:,:)=zero ;  ghc4(:,:)=zero
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
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
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
       end if
       call fourwf(1,vlocal_tmp,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,option_fft,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
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
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
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
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
     ABI_FREE(vlocal_tmp)
!    Build ghc from pieces
!    (v11,v22,Re(v12)+iIm(v12);Re(v12)-iIm(v12))(psi1;psi2): matrix product
     if (mpi_enreg%paral_spinor==0) then
       do idat=1,ndat
         do ipw=1,npw_k2
           ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2)       =ghc1(1:2,ipw+(idat-1)*npw_k2)+ghc4(1:2,ipw+(idat-1)*npw_k2)
           ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc3(1:2,ipw+(idat-1)*npw_k2)+ghc2(1:2,ipw+(idat-1)*npw_k2)
         end do
       end do
     else
       call xmpi_sum(ghc4,mpi_enreg%comm_spinor,ierr)
       call xmpi_sum(ghc3,mpi_enreg%comm_spinor,ierr)
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
     ABI_FREE(ghc1)
     ABI_FREE(ghc2)
     ABI_FREE(ghc3)
     ABI_FREE(ghc4)
   end if ! nvloc

   if (nspinortot==2)  then
     ABI_FREE(cwavef1)
     ABI_FREE(cwavef2)
   end if

   ABI_FREE(work)

!  Retrieve eventually original FFT distrib
   if(have_to_reequilibrate) then
     if(ndat > 1 ) then
       do idat=1,ndat
         do ipw = 1 ,npw_fft
           cwavef_fft_tr(1:2,  idat + ndat*(ipw-1)) = cwavef_fft(1:2,  ipw + npw_fft*(idat-1))
         end do
       end do
       call xmpi_alltoallv(cwavef_fft_tr,2*ndat*recvcount_fft, 2*ndat*recvdisp_fft, &
&       buff_wf,2*ndat*sendcount_fft,2*ndat*senddisp_fft, mpi_enreg%comm_fft,ierr)
     else
       call xmpi_alltoallv(cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, &
&       buff_wf,2*sendcount_fft,2*senddisp_fft, mpi_enreg%comm_fft,ierr)
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

!  Add metaGGA contribution
   if (associated(gs_ham%vxctaulocal)) then
     if (.not.k1_eq_k2) then
       ABI_BUG('metaGGA not allowed for k/=k_^prime!')
     end if
     if (size(gs_ham%vxctaulocal)/=gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*4) then
       ABI_BUG('wrong sizes for vxctaulocal!')
     end if
     ABI_MALLOC(ghc_mGGA,(2,npw_k2*my_nspinor*ndat))
     call getghc_mGGA(cwavef,ghc_mGGA,gbound_k1,gs_ham%gprimd,gs_ham%istwf_k,kg_k1,kpt_k1,&
&     gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
&     gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vxctaulocal,gs_ham%use_gpu_cuda)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_mGGA(1:2,1:npw_k2*my_nspinor*ndat)
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
     ghc_vectornd=zero
     call getghc_nucdip(cwavef,ghc_vectornd,gbound_k1,gs_ham%istwf_k,kg_k1,kpt_k1,&
&     gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
&     gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vectornd,gs_ham%use_gpu_cuda)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_vectornd(1:2,1:npw_k2*my_nspinor*ndat)
     ABI_FREE(ghc_vectornd)
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

     if (gs_ham%usepaw==0) gsc_ptr => nonlop_dum
     if (gs_ham%usepaw==1) gsc_ptr => gsc
     call nonlop(choice,cpopt_here,cwaveprj_nonlop,enlout,gs_ham,idir,lambda_ndat,mpi_enreg,ndat,&
&     nnlout,paw_opt,signs,gsc_ptr,tim_nonlop,cwavef,gvnlxc_,select_k=select_k_)

     if (gs_ham%usepaw==1 .and. has_fock)then
       if (fock_get_getghc_call(fock)==1) then
#if defined HAVE_GPU && defined HAVE_YAKL
         ABI_MALLOC_MANAGED(gvnlc, (/2,npw_k2*my_nspinor*ndat/))
#else
         ABI_MALLOC(gvnlc, (2,npw_k2*my_nspinor*ndat))
#endif
         gvnlc=gvnlxc_
       endif
     endif

!    Calculation of the Fock exact exchange contribution from the Fock or ACE operator
     if (has_fock) then
       if (fock_get_getghc_call(fock)==1) then
         if (gs_ham%usepaw==0) cwaveprj_idat => cwaveprj
         if (fock%use_ACE==0) then
           call timab(360,1,tsec)
           do idat=1,ndat
             if (gs_ham%usepaw==1) cwaveprj_idat => cwaveprj_fock(:,(idat-1)*my_nspinor+1:idat*my_nspinor)
             call fock_getghc(cwavef(:,1+(idat-1)*npw_k1*my_nspinor:idat*npw_k1*my_nspinor),cwaveprj_idat,&
&             gvnlxc_(:,1+(idat-1)*npw_k2*my_nspinor:idat*npw_k2*my_nspinor),gs_ham,mpi_enreg)
           end do ! idat
           call timab(360,2,tsec)
         else
           do idat=1,ndat
             call fock_ACE_getghc(cwavef(:,1+(idat-1)*npw_k1*my_nspinor:idat*npw_k1*my_nspinor),&
&             gvnlxc_(:,1+(idat-1)*npw_k2*my_nspinor:idat*npw_k2*my_nspinor),gs_ham,mpi_enreg)
           end do ! idat
         end if
       end if
     end if

   else if (type_calc == 3) then
     ! for kinetic and local only, nonlocal and vfock should be zero
     gvnlxc_(:,:) = zero
   end if ! if(type_calc...
   ABI_NVTX_END_RANGE()

!============================================================
! Assemble kinetic, local, nonlocal and Fock contributions
!============================================================
   ABI_NVTX_START_RANGE(NVTX_GETGHC_KIN)

!  Assemble modified kinetic, local and nonlocal contributions
   !  to <G|H|C(n,k)>. Take also into account build-in debugging.
   if(prtvol/=-level)then
     do idat=1,ndat
       if (k1_eq_k2) then
!      !!$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2)
         do ispinor=1,my_nspinor
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

       else
!      !!$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2)
         do ispinor=1,my_nspinor
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
       end if
     end do ! idat
   else
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
     gvnlxc_=gvnlxc_-gvnlc
#if defined HAVE_GPU && defined HAVE_YAKL
     ABI_FREE_MANAGED(gvnlc)
#else
     ABI_FREE(gvnlc)
#endif
   endif

   if (local_gvnlxc) then
#if defined HAVE_GPU && defined HAVE_YAKL
     ABI_FREE_MANAGED(gvnlxc_)
#else
     ABI_FREE(gvnlxc_)
#endif
   end if

!  Structured debugging : if prtvol=-level, stop here.
   if(prtvol==-level)then
     write(msg,'(a,i0,a)')' getghc : exit prtvol=-',level,', debugging mode => stop '
     ABI_ERROR(msg)
   end if

   if (type_calc==0.or.type_calc==2) then
     if (has_fock.and.gs_ham%usepaw==1.and.cpopt<2) then
       call pawcprj_free(cwaveprj_fock)
       ABI_FREE(cwaveprj_fock)
     end if
   end if

 end if ! type_calc

 call timab(350+tim_getghc,2,tsec)

 DBG_EXIT("COLL")

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
!! n4,n5,n6=for dimensionning of vxctaulocal
!! use_gpu_cuda=1 if Cuda (GPU) is on
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

subroutine getghc_nucdip(cwavef,ghc_vectornd,gbound_k,istwf_k,kg_k,kpt,mgfft,mpi_enreg,&
&                      ndat,ngfft,npw_k,nvloc,n4,n5,n6,my_nspinor,vectornd,use_gpu_cuda)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,my_nspinor,ndat,npw_k,nvloc,n4,n5,n6,use_gpu_cuda
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
 integer :: idat,idir,ipw,nspinortot,shift
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: scale_conversion,weight=one
 !arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: gcwavef(:,:,:),gcwavef1(:,:,:),gcwavef2(:,:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),kgkpk(:,:)
 real(dp),allocatable :: work(:,:,:,:)

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

 ABI_MALLOC(work,(2,n4,n5,n6*ndat))

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
          gcwavef(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k,idir) = &
               & cwavef(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k)*kgkpk(1:npw_k,idir)
          gcwavef(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k,idir) = &
               & cwavef(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k)*kgkpk(1:npw_k,idir)
       end do
    end do
    ABI_FREE(kgkpk)
    gcwavef = gcwavef*two_pi

    !  STEP2: Compute sum of (grad components of vectornd)*(grad components of cwavef)
    do idir=1,3
      call fourwf(1,vectornd(:,:,:,:,idir),gcwavef(:,:,idir),ghc1,work,gbound_k,gbound_k,&
           istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
           &     tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
!!$OMP PARALLEL DO
       ! DAXPY is a BLAS routine for y -> A*x + y, here x = ghc1, A = scale_conversion, and y = ghc_vectornd
       ! should be faster than explicit loop over ipw as npw_k gets large
      do idat=1,ndat
        call DAXPY(npw_k,scale_conversion,ghc1(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1,&
             & ghc_vectornd(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1)
        call DAXPY(npw_k,scale_conversion,ghc1(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1,&
             & ghc_vectornd(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1)
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
             gcwavef1(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k,idir) = &
                  & cwavef1(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k)*kgkpk(1:npw_k,idir)
             gcwavef1(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k,idir) = &
                  & cwavef1(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k)*kgkpk(1:npw_k,idir)
          end do
       end do
       gcwavef1 = gcwavef1*two_pi

       !  STEP2: Compute sum of (grad components of vectornd)*(grad components of cwavef)
       do idir=1,3
          call fourwf(1,vectornd(:,:,:,:,idir),gcwavef1(:,:,idir),ghc1,work,gbound_k,gbound_k,&
               istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
               &     tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
!!$OMP PARALLEL DO
          ! DAXPY is a BLAS routine for y -> A*x + y, here x = ghc1, A = scale_conversion, and y = ghc_vectornd
          ! should be faster than explicit loop over ipw as npw_k gets large
          do idat=1,ndat
             call DAXPY(npw_k,scale_conversion,ghc1(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1,&
                  & ghc_vectornd(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1)
             call DAXPY(npw_k,scale_conversion,ghc1(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1,&
                  & ghc_vectornd(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1)
          end do
       end do ! idir
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
             gcwavef2(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k,idir) = &
                  & cwavef2(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k)*kgkpk(1:npw_k,idir)
             gcwavef2(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k,idir) = &
                  & cwavef2(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k)*kgkpk(1:npw_k,idir)
          end do
       end do
       gcwavef2 = gcwavef2*two_pi

       !  STEP2: Compute sum of (grad components of vectornd)*(grad components of cwavef)
       do idir=1,3
          call fourwf(1,vectornd(:,:,:,:,idir),gcwavef2(:,:,idir),ghc2,work,gbound_k,gbound_k,&
               istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,ngfft,npw_k,npw_k,n4,n5,n6,2,&
               &     tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
!!$OMP PARALLEL DO
          ! DAXPY is a BLAS routine for y -> A*x + y, here x = ghc1, A = scale_conversion, and y = ghc_vectornd
          ! should be faster than explicit loop over ipw as npw_k gets large
          do idat=1,ndat
             call DAXPY(npw_k,scale_conversion,ghc2(1,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1,&
                  & ghc_vectornd(1,1+(idat-1)*npw_k+shift:npw_k+(idat-1)*npw_k+shift),1)
             call DAXPY(npw_k,scale_conversion,ghc2(2,1+(idat-1)*npw_k:npw_k+(idat-1)*npw_k),1,&
                  & ghc_vectornd(2,1+(idat-1)*npw_k+shift:npw_k+(idat-1)*npw_k+shift),1)
          end do
       end do ! idir
       ABI_FREE(gcwavef2)
       ABI_FREE(ghc2)

    end if ! end spinor 2

    ABI_FREE(cwavef1)
    ABI_FREE(cwavef2)
    ABI_FREE(kgkpk)

 end if ! nspinortot

 ABI_FREE(work)

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
!! use_gpu_cuda=1 if Cuda (GPU) is on
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
&                      ndat,ngfft,npw_k,nvloc,n4,n5,n6,my_nspinor,vxctaulocal,use_gpu_cuda)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,my_nspinor,ndat,npw_k,nvloc,n4,n5,n6,use_gpu_cuda
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
&   tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
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
&     tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
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
&     tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
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
&      tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
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
&     tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
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
&      tim_fourwf,weight,weight,use_gpu_cuda=use_gpu_cuda)
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
&                 mcg,mcprj,mgsc,mpi_enreg,natom,nband,npw_k,nspinor,select_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ibg,icg,igsc,ikpt,isppol,mcg,mcprj
 integer,intent(in) :: mgsc,natom,nband,npw_k,nspinor
!TODO : may be needed to distribute cprj over band procs
! integer,intent(in) :: mband_mem
 integer,intent(in),optional :: select_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: gsc(2,mgsc)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimenl1,dimenl2,iband,iband1,iband2,index_cg,index_cprj
 integer :: index_gsc,me,my_nspinor,paw_opt,select_k_,signs,tim_nonlop,useylm
 !character(len=500) :: msg
!arrays
 real(dp) :: enlout_dum(1),tsec(2)
 real(dp),allocatable :: cwavef(:,:),scwavef(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 if(gs_ham%usepaw==0) then
   ABI_BUG('Only compatible with PAW (usepaw=1) !')
 end if
 if(nband<0.and.(mcg<npw_k*my_nspinor.or.mgsc<npw_k*my_nspinor.or.mcprj<my_nspinor)) then
   ABI_BUG('Invalid value for mcg, mgsc or mcprj !')
 end if

!Keep track of total time spent in getgsc:
 call timab(565,1,tsec)

 gsc = zero

!Prepare some data
 ABI_MALLOC(cwavef,(2,npw_k*my_nspinor))
 ABI_MALLOC(scwavef,(2,npw_k*my_nspinor))
 if (gs_ham%usecprj==1) then
   ABI_MALLOC(cwaveprj,(natom,my_nspinor))
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

 do iband=iband1,iband2

   if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me.and.nband>0) then
! No longer needed 28/03/2020 to parallelize memory
!     gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=zero
!     index_gsc=index_gsc+npw_k*my_nspinor
     !index_cprj=index_cprj+my_nspinor
     !index_cg=index_cg+npw_k*my_nspinor

     cycle
   end if

!  Retrieve WF at (n,k)
   cwavef(:,1:npw_k*my_nspinor)=cg(:,1+index_cg:npw_k*my_nspinor+index_cg)
   if (gs_ham%usecprj==1) then
     call pawcprj_copy(cprj(:,1+index_cprj:my_nspinor+index_cprj),cwaveprj)
   end if

!  Compute <g|S|Cnk>
   call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_ham,0,(/zero/),mpi_enreg,1,1,paw_opt,&
&   signs,scwavef,tim_nonlop,cwavef,cwavef,select_k=select_k_)

   gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=scwavef(:,1:npw_k*my_nspinor)

!  End of loop over bands
   index_cprj=index_cprj+my_nspinor
   index_cg=index_cg+npw_k*my_nspinor
   index_gsc=index_gsc+npw_k*my_nspinor
 end do

!Memory deallocation
 ABI_FREE(cwavef)
 ABI_FREE(scwavef)
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
&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments

#ifdef HAVE_OPENMP
   use omp_lib
#endif

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
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(out) :: ghc(:,:),gvnlxc(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: firstelt, firstprj, lastelt, lastprj,usegvnlxc
 integer :: nthreads
 integer :: ithread
 integer :: chunk
 integer :: residuchunk
 integer :: firstband
 integer :: lastband
 integer :: spacedim, spacedim_prj
#ifdef HAVE_OPENMP
 logical :: is_nested
#ifdef HAVE_FFTW3_THREADS
 logical ::  fftw3_use_lib_threads_sav
#endif
#endif

 integer :: select_k_default

 ! *************************************************************************

 select_k_default = 1; if ( present(select_k) ) select_k_default = select_k

 spacedim     = size(cwavef  ,dim=2)/ndat
 spacedim_prj = size(cwaveprj,dim=2)/ndat


#ifdef HAVE_GPU
 usegvnlxc=1
 if (size(gvnlxc)==0) usegvnlxc=0

 if ( present(kg_fft_k) ) then
   if (present(kg_fft_kp)) then
     call getghc(cpopt,cwavef,cwaveprj,&
       &      ghc,gsc(:,1:ndat*spacedim*gs_ham%usepaw),&
       &      gs_ham,gvnlxc(:,1:ndat*spacedim*usegvnlxc),lambda, mpi_enreg,ndat,&
       &      prtvol,sij_opt,tim_getghc,type_calc,&
       &      select_k=select_k_default,kg_fft_k=kg_fft_k,kg_fft_kp=kg_fft_kp)
   else
     call getghc(cpopt,cwavef,cwaveprj,&
       &      ghc,gsc(:,1:ndat*spacedim*gs_ham%usepaw),&
       &      gs_ham,gvnlxc(:,1:ndat*spacedim*usegvnlxc),lambda, mpi_enreg,ndat,&
       &      prtvol,sij_opt,tim_getghc,type_calc,&
       &      select_k=select_k_default,kg_fft_k=kg_fft_k)
   end if
 else
   if (present(kg_fft_kp)) then
     call getghc(cpopt,cwavef,cwaveprj,&
       &      ghc,gsc(:,1:ndat*spacedim*gs_ham%usepaw),&
       &      gs_ham,gvnlxc(:,1:ndat*spacedim*usegvnlxc),lambda, mpi_enreg,ndat,&
       &      prtvol,sij_opt,tim_getghc,type_calc,&
       &      select_k=select_k_default,kg_fft_kp=kg_fft_kp)
   else
     call getghc(cpopt,cwavef,cwaveprj,&
       &      ghc,gsc(:,1:ndat*spacedim*gs_ham%usepaw),&
       &      gs_ham,gvnlxc(:,1:ndat*spacedim*usegvnlxc),lambda, mpi_enreg,ndat,&
       &      prtvol,sij_opt,tim_getghc,type_calc,&
       &      select_K=select_k_default)
   end if
 end if

 ! end GPU version
#else

    !$omp parallel default (none) &
    !$omp& private(ithread,nthreads,chunk,firstband,lastband,residuchunk,firstelt,lastelt,firstprj,lastprj,is_nested,usegvnlxc), &
    !$omp& shared(cwavef,ghc,gsc, gvnlxc,spacedim,spacedim_prj,ndat,kg_fft_k,kg_fft_kp,gs_ham,cwaveprj,mpi_enreg), &
    !$omp& firstprivate(cpopt,lambda,prtvol,sij_opt,tim_getghc,type_calc,select_k_default)
#ifdef HAVE_OPENMP
 ithread = omp_get_thread_num()
 nthreads = omp_get_num_threads()
! is_nested = omp_get_nested()
 is_nested = .false.
! call omp_set_nested(.false.)
!Ensure that libs are used without threads (mkl, openblas, fftw3, ...)
#ifdef HAVE_LINALG_MKL_THREADS
 call mkl_set_num_threads(1)
#endif
#ifdef HAVE_LINALG_OPENBLAS_THREADS
 call openblas_set_num_threads(1)
#endif
#ifdef HAVE_FFTW3_THREADS
 fftw3_use_lib_threads_sav=(.not.fftw3_spawn_threads_here(nthreads,nthreads))
 call fftw3_use_lib_threads(.false.)
#endif
#else
 ithread = 0
 nthreads = 1
#endif
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
 if (size(gvnlxc)==0) usegvnlxc=0

 if ( lastband /= 0 ) then
   firstelt = (firstband-1)*spacedim+1
   firstprj = (firstband-1)*spacedim_prj+1
   lastelt = lastband*spacedim
   lastprj = lastband*spacedim_prj
      ! Don't know how to manage optional arguments .... :(
   if ( present(kg_fft_k) ) then
     if (present(kg_fft_kp)) then
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_k=select_k_default,kg_fft_k=kg_fft_k,kg_fft_kp=kg_fft_kp)
     else
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_k=select_k_default,kg_fft_k=kg_fft_k)
     end if
   else
     if (present(kg_fft_kp)) then
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_k=select_k_default,kg_fft_kp=kg_fft_kp)
     else
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj(:,firstprj:lastprj),&
&      ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
&      gs_ham,gvnlxc(:,firstelt:lastelt*usegvnlxc),lambda, mpi_enreg,lastband-firstband+1,&
&      prtvol,sij_opt,tim_getghc,type_calc,&
&      select_K=select_k_default)
     end if
   end if
 end if
#ifdef HAVE_OPENMP
! call omp_set_nested(is_nested)
!Restire libs behavior (mkl, openblas, fftw3, ...)
#ifdef HAVE_LINALG_MKL_THREADS
 call mkl_set_num_threads(nthreads)
#endif
#ifdef HAVE_LINALG_OPENBLAS_THREADS
 call openblas_set_num_threads(nthreads)
#endif
#ifdef HAVE_FFTW3_THREADS
 call fftw3_use_lib_threads(fftw3_use_lib_threads_sav)
#endif
#endif
    !$omp end parallel

! end HAVE_GPU
#endif

end subroutine multithreaded_getghc
!!***

end module m_getghc
!!***
