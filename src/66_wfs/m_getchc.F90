!!****m* ABINIT/m_getchc
!! NAME
!!  m_getchc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space;
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LSI, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_getchc

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_getdim, pawcprj_copy
 use m_bandfft_kpt, only : bandfft_kpt, bandfft_kpt_get_ikpt
 use m_hamiltonian, only : gs_hamiltonian_type, KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME
 use m_nonlop,      only : nonlop
 use m_fock,        only : fock_common_type, fock_get_getghc_call
 use m_fock_getghc, only : fock_getghc, fock_ACE_getghc
 use m_fft,         only : fourwf
 use m_cgtools,     only : dotprod_g
 use m_getghc,      only : getghc,getghc_mGGA,getghc_nucdip
 !LTEST
 use testing
 !LTEST

 implicit none

 private
!!***

 public :: getchc     ! Compute <CP|H|C> for input vector |C> expressed in reciprocal space
 public :: getcsc     ! Compute <CP|S|C> for all input vectors |Cnk> at a given k-point
! public :: multithreaded_getchc
!!***

contains
!!***

!!****f* ABINIT/getchc
!!
!! NAME
!! getchc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space;
!! Result is put in array ghc.
!! <G|Vnonlocal + VfockACE|C> is also returned in gvnlxc if either NLoc NCPP or FockACE.
!! if required, <G|S|C> is returned in gsc (S=overlap - PAW only)
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
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFT to do in parallel
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getchc=timing code of the calling subroutine(can be set to 0 if not attributed)
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
!!      cgwf,chebfi,dfpt_cgwf,gwls_hamiltonian,ks_ddiago,lobpcgwf,m_io_kss
!!      m_rf2,mkresi,multithreaded_getchc
!!
!! CHILDREN
!!      fourwf
!!
!! SOURCE

subroutine getchc(chc_re,chc_im,cpopt,cwavef,cwavef_left,cwaveprj,cwaveprj_left,gs_ham,lambda,mpi_enreg,ndat,&
&                 npw,nspinor,prtvol,sij_opt,tim_getchc,type_calc,&
&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,ndat, prtvol
 integer,intent(in) :: npw,nspinor,sij_opt,tim_getchc,type_calc
 integer,intent(in),optional :: select_k
 real(dp),intent(in) :: lambda
 real(dp),intent(out) :: chc_re,chc_im
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
 real(dp),intent(inout) :: cwavef(:,:),cwavef_left(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:),cwaveprj_left(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=114,re=1,im=2,tim_fourwf=33
 integer :: choice,cplex,cpopt_here,i1,i2,i3,idat,idir,ierr
 integer :: ig,igspinor,ii,iispinor,ikpt_this_proc,ipw,ispinor,my_nspinor
 integer :: nnlout,npw_fft,npw_k1,npw_k2,nspinortot
 integer :: paw_opt,select_k_,shift1,shift2,signs,tim_nonlop
 logical :: k1_eq_k2,have_to_reequilibrate,has_fock
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: ghcim,ghcre,weight
 character(len=500) :: msg
!arrays
 integer, pointer :: gbound_k1(:,:),gbound_k2(:,:),kg_k1(:,:),kg_k2(:,:)
 integer, ABI_CONTIGUOUS pointer :: indices_pw_fft(:),kg_k_fft(:,:)
 integer, ABI_CONTIGUOUS pointer :: recvcount_fft(:),recvdisp_fft(:)
 integer, ABI_CONTIGUOUS pointer ::  sendcount_fft(:),senddisp_fft(:)
 integer, allocatable:: dimcprj(:)
 real(dp) :: enlout(ndat),enlout_im(ndat),lambda_ndat(ndat),tsec(2),dotr,doti
 real(dp),target :: nonlop_dum(1,1)
 real(dp),allocatable :: buff_wf(:,:),cwavef1(:,:),cwavef2(:,:),cwavef_fft(:,:),cwavef_fft_tr(:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),ghc3(:,:),ghc4(:,:),ghc_mGGA(:,:),ghc_vectornd(:,:)
 real(dp),allocatable :: gvnlc(:,:),vlocal_tmp(:,:,:),work(:,:,:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:)
 real(dp), pointer :: kinpw_k1(:),kinpw_k2(:),kpt_k1(:),kpt_k2(:)
 real(dp), pointer :: gsc_ptr(:,:)
 type(fock_common_type),pointer :: fock
 type(pawcprj_type),pointer :: cwaveprj_fock(:,:),cwaveprj_idat(:,:),cwaveprj_nonlop(:,:)
 !LTEST
 real(dp) :: chc_re_tmp,chc_im_tmp
 !LTEST

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in getchc:
 call timab(1370,1,tsec)

!!Structured debugging if prtvol==-level
! if(prtvol==-level)then
!   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' getchc : enter, debugging '
!   call wrtout(std_out,msg,'PERS')
! end if

 ABI_ALLOCATE(ghc,(2,npw*nspinor))
 ABI_ALLOCATE(gvnlxc,(2,npw*nspinor))
 ABI_ALLOCATE(gsc,(2,npw*nspinor))

!LTEST
! call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
!&                 prtvol,sij_opt,tim_getchc,type_calc,&
!&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments
!
! call dotprod_g(chc_re_tmp,chc_im_tmp,gs_ham%istwf_k,npw*nspinor,2,cwavef_left,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!LTEST

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
   msg='wrong size for cwavef!'
   MSG_BUG(msg)
 end if
! if (size(ghc)<2*npw_k2*my_nspinor*ndat) then
!   msg='wrong size for ghc!'
!   MSG_BUG(msg)
! end if
! if (size(gvnlxc)<2*npw_k2*my_nspinor*ndat) then
!   msg='wrong size for gvnlxc!'
!   MSG_BUG(msg)
! end if
! if (sij_opt==1) then
!   if (size(gsc)<2*npw_k2*my_nspinor*ndat) then
!     msg='wrong size for gsc!'
!     MSG_BUG(msg)
!   end if
! end if
 if (gs_ham%usepaw==1.and.cpopt>=0) then
   if (size(cwaveprj)<gs_ham%natom*my_nspinor*ndat) then
     msg='wrong size for cwaveprj!'
     MSG_BUG(msg)
   end if
 end if
 if (gs_ham%usepaw==1) then
   if (size(cwaveprj_left)<gs_ham%natom*my_nspinor*ndat) then
     msg='wrong size for cwaveprj_left!'
     MSG_BUG(msg)
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
   msg='paral_kgb=1 not allowed for k/=k_^prime!'
   MSG_BUG(msg)
 end if

!Do we add Fock exchange term ?
 has_fock=(associated(gs_ham%fockcommon))
 if (has_fock) then
   MSG_BUG('Fock not implemented yet')
 end if
! if (has_fock) fock => gs_ham%fockcommon

!Parallelization over spinors management
 nspinortot=gs_ham%nspinor
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

 if ((type_calc==0).or.(type_calc==1).or.(type_calc==3)) then

!  Need a Vlocal
   if (.not.associated(gs_ham%vlocal)) then
     MSG_BUG("We need vlocal in gs_ham!")
   end if

!  fourwf can only process with one value of istwf_k
   if (.not.k1_eq_k2) then
     MSG_BUG('vlocal (fourwf) cannot be computed with k/=k^prime!')
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
     ABI_ALLOCATE(buff_wf,(2,npw_k1*ndat) )
     ABI_ALLOCATE(cwavef_fft,(2,npw_fft*ndat) )
     if(ndat>1) then
       ABI_ALLOCATE(cwavef_fft_tr, (2,npw_fft*ndat))
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
   ABI_ALLOCATE(work,(2,gs_ham%n4,gs_ham%n5,gs_ham%n6*ndat))
   weight=one

   if (nspinortot==2) then
     ABI_ALLOCATE(cwavef1,(2,npw_k1*ndat))
     ABI_ALLOCATE(cwavef2,(2,npw_k1*ndat))
     do idat=1,ndat
       do ipw=1,npw_k1
         cwavef1(1:2,ipw+(idat-1)*npw_k1)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k1)
         cwavef2(1:2,ipw+(idat-1)*npw_k1)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k1+shift1)
       end do
     end do
!    call cg_zcopy(npw_k1*ndat,cwavef(1,1),cwavef1)
!    call cg_zcopy(npw_k1*ndat,cwavef(1,1+shift1),cwavef2)
   end if

!  Treat scalar local potentials
   if (gs_ham%nvloc==1) then

     if (nspinortot==1) then

       if (have_to_reequilibrate) then
         call fourwf(1,gs_ham%vlocal,cwavef_fft,cwavef_fft,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k_fft,kg_k_fft,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_fft,npw_fft,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_ham%use_gpu_cuda)
       else
         call fourwf(1,gs_ham%vlocal,cwavef,ghc,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_ham%use_gpu_cuda)
       end if

     else ! nspinortot==2

       if (nspinor1TreatedByThisProc) then
         ABI_ALLOCATE(ghc1,(2,npw_k2*ndat))
         call fourwf(1,gs_ham%vlocal,cwavef1,ghc1,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw =1, npw_k2
             ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2)=ghc1(1:2,ipw+(idat-1)*npw_k2)
           end do
         end do
         ABI_DEALLOCATE(ghc1)
       end if ! spin 1 treated by this proc

       if (nspinor2TreatedByThisProc) then
         ABI_ALLOCATE(ghc2,(2,npw_k2*ndat))
         call fourwf(1,gs_ham%vlocal,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
&         gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw=1,npw_k2
             ghc(1:2,ipw+(idat-1)*my_nspinor*npw_k2+shift2)=ghc2(1:2,ipw+(idat-1)*npw_k2)
           end do
         end do
         ABI_DEALLOCATE(ghc2)
       end if ! spin 2 treated by this proc

     end if ! npsinortot

!    Treat non-collinear local potentials
   else if (gs_ham%nvloc==4) then

     ABI_ALLOCATE(ghc1,(2,npw_k2*ndat))
     ABI_ALLOCATE(ghc2,(2,npw_k2*ndat))
     ABI_ALLOCATE(ghc3,(2,npw_k2*ndat))
     ABI_ALLOCATE(ghc4,(2,npw_k2*ndat))
     ghc1(:,:)=zero; ghc2(:,:)=zero; ghc3(:,:)=zero ;  ghc4(:,:)=zero
     ABI_ALLOCATE(vlocal_tmp,(gs_ham%n4,gs_ham%n5,gs_ham%n6))
!    ghc1=v11*phi1
     vlocal_tmp(:,:,:)=gs_ham%vlocal(:,:,:,1)
     if (nspinor1TreatedByThisProc) then
       call fourwf(1,vlocal_tmp,cwavef1,ghc1,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
!    ghc2=v22*phi2
     vlocal_tmp(:,:,:)=gs_ham%vlocal(:,:,:,2)
     if (nspinor2TreatedByThisProc) then
       call fourwf(1,vlocal_tmp,cwavef2,ghc2,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
     cplex=2
     ABI_ALLOCATE(vlocal_tmp,(cplex*gs_ham%n4,gs_ham%n5,gs_ham%n6))
!    ghc3=(re(v12)-im(v12))*phi1
     do i3=1,gs_ham%n6
       do i2=1,gs_ham%n5
         do i1=1,gs_ham%n4
           vlocal_tmp(2*i1-1,i2,i3)= gs_ham%vlocal(i1,i2,i3,3)
           vlocal_tmp(2*i1  ,i2,i3)=-gs_ham%vlocal(i1,i2,i3,4)
         end do
       end do
     end do
     if (nspinor1TreatedByThisProc) then
       call fourwf(cplex,vlocal_tmp,cwavef1,ghc3,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
!    ghc4=(re(v12)+im(v12))*phi2
     if (nspinor2TreatedByThisProc) then
       do i3=1,gs_ham%n6
         do i2=1,gs_ham%n5
           do i1=1,gs_ham%n4
             vlocal_tmp(2*i1,i2,i3)=-vlocal_tmp(2*i1,i2,i3)
           end do
         end do
       end do
       call fourwf(cplex,vlocal_tmp,cwavef2,ghc4,work,gbound_k1,gbound_k2,&
&       gs_ham%istwf_k,kg_k1,kg_k2,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw_k1,npw_k2,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
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
     ABI_DEALLOCATE(ghc1)
     ABI_DEALLOCATE(ghc2)
     ABI_DEALLOCATE(ghc3)
     ABI_DEALLOCATE(ghc4)
   end if ! nvloc

   if (nspinortot==2)  then
     ABI_DEALLOCATE(cwavef1)
     ABI_DEALLOCATE(cwavef2)
   end if
   ABI_DEALLOCATE(work)

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
     ABI_DEALLOCATE(buff_wf)
     ABI_DEALLOCATE(cwavef_fft)
     if(ndat > 1) then
       ABI_DEALLOCATE(cwavef_fft_tr)
     end if
   end if

!  Add metaGGA contribution
   if (associated(gs_ham%vxctaulocal)) then
     if (.not.k1_eq_k2) then
       MSG_BUG('metaGGA not allowed for k/=k_^prime!')
     end if
     if (size(gs_ham%vxctaulocal)/=gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*4) then
       MSG_BUG('wrong sizes for vxctaulocal!')
     end if
     ABI_ALLOCATE(ghc_mGGA,(2,npw_k2*my_nspinor*ndat))
     call getghc_mGGA(cwavef,ghc_mGGA,gbound_k1,gs_ham%gprimd,gs_ham%istwf_k,kg_k1,kpt_k1,&
&     gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
&     gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vxctaulocal,gs_ham%use_gpu_cuda)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_mGGA(1:2,1:npw_k2*my_nspinor*ndat)
     ABI_DEALLOCATE(ghc_mGGA)
   end if

   !  Add nuclear dipole moment contribution
   if (associated(gs_ham%vectornd)) then
     if (.not.k1_eq_k2) then
       MSG_BUG('nuclear dipole vector potential not allowed for k/=k_^prime!')
     end if
     if (size(gs_ham%vectornd)/=gs_ham%n4*gs_ham%n5*gs_ham%n6*gs_ham%nvloc*3) then
       MSG_BUG('wrong sizes for vectornd in getchc!')
     end if
     ABI_ALLOCATE(ghc_vectornd,(2,npw_k2*my_nspinor*ndat))
     call getghc_nucdip(cwavef,ghc_vectornd,gbound_k1,gs_ham%istwf_k,kg_k1,kpt_k1,&
&     gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
&     gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vectornd,gs_ham%use_gpu_cuda)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_vectornd(1:2,1:npw_k2*my_nspinor*ndat)
     ABI_DEALLOCATE(ghc_vectornd)
   end if

   call dotprod_g(chc_re,chc_im,gs_ham%istwf_k,npw*nspinor,2,cwavef_left,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

 else
   chc_re = zero
   chc_im = zero
 end if ! type_calc

 if ((type_calc==0).or.(type_calc==2).or.(type_calc==3)) then

!============================================================
! Application of the non-local potential and the Fock potential
!============================================================

   if ((type_calc==0).or.(type_calc==2)) then
     signs=1 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=15
     cpopt_here=-1;if (gs_ham%usepaw==1) cpopt_here=cpopt
!     if (has_fock) then
!       if (gs_ham%usepaw==1) then
!         cpopt_here=max(cpopt,0)
!         if (cpopt<2) then
!           ABI_DATATYPE_ALLOCATE(cwaveprj_fock,(gs_ham%natom,my_nspinor*ndat))
!           ABI_ALLOCATE(dimcprj,(gs_ham%natom))
!           call pawcprj_getdim(dimcprj,gs_ham%natom,gs_ham%nattyp,gs_ham%ntypat,&
!&           gs_ham%typat,fock%pawtab,'O')
!           call pawcprj_alloc(cwaveprj_fock,0,dimcprj)
!           ABI_DEALLOCATE(dimcprj)
!         else
!           cwaveprj_fock=>cwaveprj
!         end if
!         cwaveprj_nonlop=>cwaveprj_fock
!       else
!         cwaveprj_nonlop=>cwaveprj
!         cwaveprj_fock=>cwaveprj
!       end if
!     else
!     end if
     paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3
     lambda_ndat = lambda

     !LTEST
!     signs=2 
!     call nonlop(choice,cpopt_here,cwaveprj,enlout,gs_ham,idir,lambda_ndat,mpi_enreg,ndat,&
!&     nnlout,paw_opt,signs,gsc,tim_nonlop,cwavef,gvnlxc,select_k=select_k_,cprjin_left=cwaveprj_left)
!
!     call dotprod_g(dotr,doti,gs_ham%istwf_k,npw*nspinor,2,cwavef_left,gvnlxc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     signs=1
!     !LTEST

     call nonlop(choice,cpopt_here,cwaveprj,enlout,gs_ham,idir,lambda_ndat,mpi_enreg,ndat,&
&     nnlout,paw_opt,signs,gsc,tim_nonlop,cwavef,gvnlxc,select_k=select_k_,&
&     cprjin_left=cwaveprj_left,enlout_im=enlout_im)

!     !LTEST
!     call writeout(999,'dotr    ',dotr)
!     call writeout(999,'dotr bis',enlout(1))
!     call writeout(999,'doti    ',doti)
!     call writeout(999,'doti bis',enlout_im(1))
!     !LTEST
     do idat=1,ndat
       chc_re = chc_re + enlout(idat)
       chc_im = chc_im + enlout_im(idat)
     end do

!     if (gs_ham%usepaw==1 .and. has_fock)then
!       if (fock_get_getghc_call(fock)==1) then
!         ABI_ALLOCATE(gvnlc,(2,npw_k2*my_nspinor*ndat))
!         gvnlc=gvnlxc
!       endif
!     endif

!    Calculation of the Fock exact exchange contribution from the Fock or ACE operator
!     if (has_fock) then
!       if (fock_get_getghc_call(fock)==1) then
!         if (gs_ham%usepaw==0) cwaveprj_idat => cwaveprj
!         do idat=1,ndat
!           if (fock%use_ACE==0) then
!             if (gs_ham%usepaw==1) cwaveprj_idat => cwaveprj_fock(:,(idat-1)*my_nspinor+1:idat*my_nspinor)
!             call fock_getghc(cwavef(:,1+(idat-1)*npw_k1*my_nspinor:idat*npw_k1*my_nspinor),cwaveprj_idat,&
!&             gvnlxc(:,1+(idat-1)*npw_k2*my_nspinor:idat*npw_k2*my_nspinor),gs_ham,mpi_enreg)
!           else
!             call fock_ACE_getghc(cwavef(:,1+(idat-1)*npw_k1*my_nspinor:idat*npw_k1*my_nspinor),&
!&             gvnlxc(:,1+(idat-1)*npw_k2*my_nspinor:idat*npw_k2*my_nspinor),gs_ham,mpi_enreg)
!           end if
!         end do ! idat
!       end if
!     end if

   end if ! if(type_calc...

!============================================================
! Assemble kinetic, local, nonlocal and Fock contributions
!============================================================

!  Add modified kinetic contributions
   !  to <CP|H|C(n,k)>.
   do idat=1,ndat
!    !!$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2)
     do ispinor=1,my_nspinor
       do ig=1,npw_k2
         igspinor=ig+npw_k2*(ispinor-1)+npw_k2*my_nspinor*(idat-1)
         if(kinpw_k2(ig)<huge(zero)*1.d-11)then
           chc_re = chc_re +  kinpw_k2(ig)*cwavef(re,igspinor)*cwavef_left(re,igspinor)
           chc_re = chc_re +  kinpw_k2(ig)*cwavef(im,igspinor)*cwavef_left(im,igspinor)
           chc_im = chc_im +  kinpw_k2(ig)*cwavef(im,igspinor)*cwavef_left(re,igspinor)
           chc_im = chc_im -  kinpw_k2(ig)*cwavef(re,igspinor)*cwavef_left(im,igspinor)
         end if
       end do ! ig
     end do ! ispinor
   end do ! idat

!  Special case of PAW + Fock : only return Fock operator contribution in gvnlxc
!   if (gs_ham%usepaw==1 .and. has_fock)then
!     gvnlxc=gvnlxc-gvnlc
!     ABI_DEALLOCATE(gvnlc)
!   endif
!
!   if ((type_calc==0).or.(type_calc==2)) then
!     if (has_fock.and.gs_ham%usepaw==1.and.cpopt<2) then
!       call pawcprj_free(cwaveprj_fock)
!       ABI_DATATYPE_DEALLOCATE(cwaveprj_fock)
!     end if
!   end if

 end if ! type_calc

 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlxc)
 ABI_DEALLOCATE(gsc)

!LTEST
! if (abs(chc_re-chc_re_tmp)>tol8) then
!   call writeout(999,'chc_re (dif)',chc_re-chc_re_tmp)
!   call writeout(999,'chc_im (dif)',chc_im-chc_im_tmp)
!   MSG_ERROR('ERROR in getchc : dif re is too large')
! end if
! if (abs(chc_im-chc_im_tmp)>tol8) then
!   call writeout(999,'chc_re (dif)',chc_re-chc_re_tmp)
!   call writeout(999,'chc_im (dif)',chc_im-chc_im_tmp)
!   MSG_ERROR('ERROR in getchc : dif im is too large')
! end if
! call writeout(999,'chc_re_tmp',chc_re_tmp)
! call writeout(999,'chc_im_tmp',chc_im_tmp)
!LTEST
!
 call timab(1370,2,tsec)

 DBG_EXIT("COLL")

end subroutine getchc
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/getcsc
!! NAME
!! getcsc
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
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,timab,xmpi_sum
!!
!! SOURCE

subroutine getcsc(csc,cpopt,cwavef,cwavef_left,cprj,cprj_left,gs_ham,mpi_enreg,ndat,&
&                 npw,nspinor,prtvol,tim_getcsc,&
&                 select_k) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,ndat,prtvol
 integer,intent(in) :: npw,nspinor,tim_getcsc
 integer,intent(in),optional :: select_k
 real(dp),intent(out) :: csc(2*ndat)
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout) :: gs_ham
!arrays
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(inout),target :: cwavef_left(:,:)
 type(pawcprj_type),intent(inout) :: cprj(:,:)
 type(pawcprj_type),intent(inout),target :: cprj_left(:,:)

!Local variables-------------------------------
!scalars
 integer :: band_shift,choice,dimenl1,dimenl2,iband,idat,idir,ierr,index_cg,index_cprj
 integer :: index_gsc,me,my_nspinor,paw_opt,select_k_,signs,tim_nonlop,useylm,nnlout
 !character(len=500) :: msg
!arrays
 real(dp) :: enlout(1),enlout_im(1),tsec(2)
 real(dp),allocatable :: gsc(:,:),gvnlxc(:,:)
 !LTEST
 real(dp) :: csc_re_tmp,csc_im_tmp
 real(dp), pointer :: cwavef_left_oneband(:,:)
 type(pawcprj_type),pointer :: cprj_left_oneband(:,:)
 !LTEST
! *********************************************************************

 DBG_ENTER("COLL")

 call timab(1360+tim_getcsc,1,tsec)

!!Compatibility tests
! my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
! if(gs_ham%usepaw==0) then
!   MSG_BUG('Only compatible with PAW (usepaw=1) !')
! end if
! if(nband<0.and.(mcg<npw_k*my_nspinor.or.mgsc<npw_k*my_nspinor.or.mcprj<my_nspinor)) then
!   MSG_BUG('Invalid value for mcg, mgsc or mcprj !')
! end if

 !LTEST
! ABI_ALLOCATE(gsc,(2,nspinor*npw))
! ABI_ALLOCATE(gvnlxc,(2,nspinor*npw))
!
! signs=2
! call nonlop(choice,cpopt,cwaveprj,enlout,gs_ham,idir,(/zero/),mpi_enreg,ndat,&
!& nnlout,paw_opt,signs,gsc,tim_nonlop,cwavef,gvnlxc,select_k=select_k_)
!
! call dotprod_g(csc_re_tmp,csc_im_tmp,gs_ham%istwf_k,npw*nspinor,2,cwavef_left,gsc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
! ABI_DEALLOCATE(gsc)
! ABI_DEALLOCATE(gvnlxc)
 !LTEST

 do idat=1,ndat
   band_shift = (idat-1)*npw*nspinor
   cwavef_left_oneband => cwavef_left(:,1+band_shift:npw*nspinor+band_shift)
   call timab(1361,1,tsec)
   call dotprod_g(csc(2*idat-1),csc(2*idat),gs_ham%istwf_k,npw*nspinor,2,cwavef_left_oneband,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   call timab(1361,2,tsec)
 end do

 if (gs_ham%usepaw==1) then
   select_k_=1;if (present(select_k)) select_k_=select_k
   choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=16 ; paw_opt=3
   ABI_ALLOCATE(gsc,(0,0))
   ABI_ALLOCATE(gvnlxc,(0,0))
   signs=1
   do idat=1,ndat
     cprj_left_oneband => cprj_left(:,idat:idat)
     call nonlop(choice,cpopt,cprj,enlout,gs_ham,idir,(/zero/),mpi_enreg,1,&
&     nnlout,paw_opt,signs,gsc,tim_nonlop,cwavef,gvnlxc,select_k=select_k_,&
&     cprjin_left=cprj_left_oneband,enlout_im=enlout_im)
     csc(2*idat-1) = csc(2*idat-1) + enlout(1)
     csc(2*idat  ) = csc(2*idat  ) + enlout_im(1)
   end do
   ABI_DEALLOCATE(gsc)
   ABI_DEALLOCATE(gvnlxc)
 end if
! !LTEST
! call writeout(999,'dotr    ',dotr)
! call writeout(999,'dotr bis',enlout(1))
! call writeout(999,'doti    ',doti)
! call writeout(999,'doti bis',enlout_im(1))
! !LTEST

!LTEST
! if (abs(csc_re-csc_re_tmp)>tol8) then
!   call writeout(999,'csc_re (dif)',csc_re-csc_re_tmp)
!   call writeout(999,'csc_im (dif)',csc_im-csc_im_tmp)
!   MSG_ERROR('ERROR in getcsc : dif re is too large')
! end if
! if (abs(csc_im-csc_im_tmp)>tol8) then
!   call writeout(999,'csc_re (dif)',csc_re-csc_re_tmp)
!   call writeout(999,'csc_im (dif)',csc_im-csc_im_tmp)
!   MSG_ERROR('ERROR in getcsc : dif im is too large')
! end if
 !LTEST
! call writeout(999,'chc_re_tmp',chc_re_tmp)
!!Prepare some data
! ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
! ABI_ALLOCATE(scwavef,(2,npw_k*my_nspinor))
! if (gs_ham%usecprj==1) then
!   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor))
!   call pawcprj_alloc(cwaveprj,0,gs_ham%dimcprj)
! else
!   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
! end if
! dimenl1=gs_ham%dimekb1;dimenl2=natom;tim_nonlop=0
! choice=1;signs=2;cpopt=-1+3*gs_ham%usecprj;paw_opt=3;useylm=1
! select_k_=1;if (present(select_k)) select_k_=select_k
! me=mpi_enreg%me_kpt
!
!!Loop over bands
! index_cprj=ibg;index_cg=icg;index_gsc=igsc
! if (nband>0) then
!   iband1=1;iband2=nband
! else if (nband<0) then
!   iband1=abs(nband);iband2=iband1
!   index_cprj=index_cprj+(iband1-1)*my_nspinor
!   index_cg  =index_cg  +(iband1-1)*npw_k*my_nspinor
!   index_gsc =index_gsc +(iband1-1)*npw_k*my_nspinor
! end if
!
! do iband=iband1,iband2
!
!   if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me.and.nband>0) then
!     gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=zero
!     index_cprj=index_cprj+my_nspinor
!     index_cg=index_cg+npw_k*my_nspinor
!     index_gsc=index_gsc+npw_k*my_nspinor
!     cycle
!   end if
!
!!  Retrieve WF at (n,k)
!   cwavef(:,1:npw_k*my_nspinor)=cg(:,1+index_cg:npw_k*my_nspinor+index_cg)
!   if (gs_ham%usecprj==1) then
!     call pawcprj_copy(cprj(:,1+index_cprj:my_nspinor+index_cprj),cwaveprj)
!   end if
!
!!  Compute <g|S|Cnk>
!   call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_ham,0,(/zero/),mpi_enreg,1,1,paw_opt,&
!&   signs,scwavef,tim_nonlop,cwavef,cwavef,select_k=select_k_)
!
!   gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=scwavef(:,1:npw_k*my_nspinor)
!
!!  End of loop over bands
!   index_cprj=index_cprj+my_nspinor
!   index_cg=index_cg+npw_k*my_nspinor
!   index_gsc=index_gsc+npw_k*my_nspinor
! end do
!
!!Reduction in case of parallelization
! if ((xmpi_paral==1)) then
!   call timab(48,1,tsec)
!   call xmpi_sum(gsc,mpi_enreg%comm_band,ierr)
!   call timab(48,2,tsec)
! end if
!
!!Memory deallocation
! ABI_DEALLOCATE(cwavef)
! ABI_DEALLOCATE(scwavef)
! if (gs_ham%usecprj==1) then
!   call pawcprj_free(cwaveprj)
! end if
! ABI_DATATYPE_DEALLOCATE(cwaveprj)
!
 call timab(1360+tim_getcsc,2,tsec)

 DBG_EXIT("COLL")

end subroutine getcsc
!!***

!!****f* ABINIT/multithreaded_getchc
!!
!! NAME
!! multithreaded_getchc
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2016-2020 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!! mpi_enreg=informations about MPI parallelization
!! ndat=number of FFT to do in parallel
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getchc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+kinetic only (added to the exixting Hamiltonian)
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
!! PARENTS
!!      m_lobpcgwf,prep_getchc
!!
!! CHILDREN
!!      getchc,mkl_set_num_threads,omp_set_nested
!!
!! SOURCE

!subroutine multithreaded_getchc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
!&                 prtvol,sij_opt,tim_getchc,type_calc,&
!&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments
!
!#ifdef HAVE_OPENMP
!   use omp_lib
!#endif
!
!!Arguments ------------------------------------
!!scalars
! integer,intent(in) :: cpopt,ndat, prtvol
! integer,intent(in) :: sij_opt,tim_getchc,type_calc
! integer,intent(in),optional :: select_k
! real(dp),intent(in) :: lambda
! type(MPI_type),intent(in) :: mpi_enreg
! type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!!arrays
! integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
! real(dp),intent(out),target :: gsc(:,:)
! real(dp),intent(inout) :: cwavef(:,:)
! real(dp),intent(out) :: ghc(:,:),gvnlxc(:,:)
! type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)
!
!!Local variables-------------------------------
!!scalars
! integer :: firstelt, lastelt
! integer :: nthreads
! integer :: ithread
! integer :: chunk
! integer :: residuchunk
! integer :: firstband
! integer :: lastband
! integer :: spacedim
!#ifdef HAVE_OPENMP
! logical :: is_nested
!#endif
!
! integer :: select_k_default
!
! ! *************************************************************************
!
! select_k_default = 1; if ( present(select_k) ) select_k_default = select_k
!
! spacedim = size(cwavef,dim=2)/ndat
!
!    !$omp parallel default (none) private(ithread,nthreads,chunk,firstband,lastband,residuchunk,firstelt,lastelt, is_nested), &
!    !$omp& shared(cwavef,ghc,gsc, gvnlxc,spacedim,ndat,kg_fft_k,kg_fft_kp,gs_ham,cwaveprj,mpi_enreg), &
!    !$omp& firstprivate(cpopt,lambda,prtvol,sij_opt,tim_getchc,type_calc,select_k_default)
!#ifdef HAVE_OPENMP
! ithread = omp_get_thread_num()
! nthreads = omp_get_num_threads()
! is_nested = omp_get_nested()
! call omp_set_nested(.false.)
!#ifdef HAVE_LINALG_MKL_THREADS
! call mkl_set_num_threads(1)
!#endif
!#else
! ithread = 0
! nthreads = 1
!#endif
! chunk = ndat/nthreads ! Divide by 2 to construct chunk of even number of bands
! residuchunk = ndat - nthreads*chunk
! if ( ithread < nthreads-residuchunk ) then
!   firstband = ithread*chunk+1
!   lastband = (ithread+1)*chunk
! else
!   firstband = (nthreads-residuchunk)*chunk + ( ithread -(nthreads-residuchunk) )*(chunk+1) +1
!   lastband = firstband+chunk
! end if
!
! if ( lastband /= 0 ) then
!   firstelt = (firstband-1)*spacedim+1
!   lastelt = lastband*spacedim
!      ! Don't know how to manage optional arguments .... :(
!   if ( present(kg_fft_k) ) then
!     if (present(kg_fft_kp)) then
!       call getchc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
!       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getchc,type_calc,&
!       select_k=select_k_default,kg_fft_k=kg_fft_k,kg_fft_kp=kg_fft_kp)
!     else
!       call getchc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
!       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getchc,type_calc,&
!       select_k=select_k_default,kg_fft_k=kg_fft_k)
!     end if
!   else
!     if (present(kg_fft_kp)) then
!       call getchc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
!       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getchc,type_calc,&
!       select_k=select_k_default,kg_fft_kp=kg_fft_kp)
!     else
!       call getchc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
!       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getchc,type_calc,&
!       select_K=select_k_default)
!     end if
!   end if
! end if
!#ifdef HAVE_OPENMP
! call omp_set_nested(is_nested)
!#ifdef HAVE_LINALG_MKL_THREADS
! call mkl_set_num_threads(nthreads)
!#endif
!#endif
!    !$omp end parallel
!
!end subroutine multithreaded_getchc
!!***

end module m_getchc
!!***
