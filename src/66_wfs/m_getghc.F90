!!****m* ABINIT/m_getghc
!! NAME
!!  m_getghc
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

 implicit none

 private
!!***

 public :: getghc     ! Compute <G|H|C> for input vector |C> expressed in reciprocal space
 public :: getgsc     ! Compute <G|S|C> for all input vectors |Cnk> at a given k-point
 public :: multithreaded_getghc
!!***

contains
!!***

!!****f* ABINIT/getghc
!!
!! NAME
!! getghc
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
!!      m_cgwf,m_chebfi,m_dfpt_cgwf,m_dft_energy,m_getghc,m_gwls_hamiltonian
!!      m_io_kss,m_ksdiago,m_lobpcgwf_old,m_orbmag,m_rf2
!!
!! CHILDREN
!!      getghc,mkl_set_num_threads,omp_set_nested
!!
!! SOURCE

subroutine getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlxc,lambda,mpi_enreg,ndat,&
&                 prtvol,sij_opt,tim_getghc,type_calc,&
&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments

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
 integer,parameter :: level=114,re=1,im=2,tim_fourwf=1
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
 real(dp) :: enlout(ndat),lambda_ndat(ndat),tsec(2)
 real(dp),target :: nonlop_dum(1,1)
 real(dp),allocatable :: buff_wf(:,:),cwavef1(:,:),cwavef2(:,:),cwavef_fft(:,:),cwavef_fft_tr(:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),ghc3(:,:),ghc4(:,:),ghc_mGGA(:,:),ghc_vectornd(:,:)
 real(dp),allocatable :: gvnlc(:,:),vlocal_tmp(:,:,:),work(:,:,:,:)
 real(dp), pointer :: kinpw_k1(:),kinpw_k2(:),kpt_k1(:),kpt_k2(:)
 real(dp), pointer :: gsc_ptr(:,:)
 type(fock_common_type),pointer :: fock
 type(pawcprj_type),pointer :: cwaveprj_fock(:,:),cwaveprj_idat(:,:),cwaveprj_nonlop(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in getghc:
 call timab(200+tim_getghc,1,tsec)

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
   msg='wrong size for cwavef!'
   MSG_BUG(msg)
 end if
 if (size(ghc)<2*npw_k2*my_nspinor*ndat) then
   msg='wrong size for ghc!'
   MSG_BUG(msg)
 end if
 if (size(gvnlxc)<2*npw_k2*my_nspinor*ndat) then
   msg='wrong size for gvnlxc!'
   MSG_BUG(msg)
 end if
 if (sij_opt==1) then
   if (size(gsc)<2*npw_k2*my_nspinor*ndat) then
     msg='wrong size for gsc!'
     MSG_BUG(msg)
   end if
 end if
 if (gs_ham%usepaw==1.and.cpopt>=0) then
   if (size(cwaveprj)<gs_ham%natom*my_nspinor*ndat) then
     msg='wrong size for cwaveprj!'
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
 if (has_fock) fock => gs_ham%fockcommon

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
       MSG_BUG('wrong sizes for vectornd in getghc!')
     end if
     ABI_ALLOCATE(ghc_vectornd,(2,npw_k2*my_nspinor*ndat))
     call getghc_nucdip(cwavef,ghc_vectornd,gbound_k1,gs_ham%istwf_k,kg_k1,kpt_k1,&
&     gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,npw_k1,gs_ham%nvloc,&
&     gs_ham%n4,gs_ham%n5,gs_ham%n6,my_nspinor,gs_ham%vectornd,gs_ham%use_gpu_cuda)
     ghc(1:2,1:npw_k2*my_nspinor*ndat)=ghc(1:2,1:npw_k2*my_nspinor*ndat)+ghc_vectornd(1:2,1:npw_k2*my_nspinor*ndat)
     ABI_DEALLOCATE(ghc_vectornd)
   end if

 end if ! type_calc


 if ((type_calc==0).or.(type_calc==2).or.(type_calc==3)) then

!============================================================
! Application of the non-local potential and the Fock potential
!============================================================

   if ((type_calc==0).or.(type_calc==2)) then
     signs=2 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=1
     cpopt_here=-1;if (gs_ham%usepaw==1) cpopt_here=cpopt
     if (has_fock) then
       if (gs_ham%usepaw==1) then
         cpopt_here=max(cpopt,0)
         if (cpopt<2) then
           ABI_DATATYPE_ALLOCATE(cwaveprj_fock,(gs_ham%natom,my_nspinor*ndat))
           ABI_ALLOCATE(dimcprj,(gs_ham%natom))
           call pawcprj_getdim(dimcprj,gs_ham%natom,gs_ham%nattyp,gs_ham%ntypat,&
&           gs_ham%typat,fock%pawtab,'O')
           call pawcprj_alloc(cwaveprj_fock,0,dimcprj)
           ABI_DEALLOCATE(dimcprj)
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
&     nnlout,paw_opt,signs,gsc_ptr,tim_nonlop,cwavef,gvnlxc,select_k=select_k_)

     if (gs_ham%usepaw==1 .and. has_fock)then
       if (fock_get_getghc_call(fock)==1) then
         ABI_ALLOCATE(gvnlc,(2,npw_k2*my_nspinor*ndat))
         gvnlc=gvnlxc
       endif
     endif

!    Calculation of the Fock exact exchange contribution from the Fock or ACE operator
     if (has_fock) then
       if (fock_get_getghc_call(fock)==1) then
         if (gs_ham%usepaw==0) cwaveprj_idat => cwaveprj
         do idat=1,ndat
           if (fock%use_ACE==0) then
             if (gs_ham%usepaw==1) cwaveprj_idat => cwaveprj_fock(:,(idat-1)*my_nspinor+1:idat*my_nspinor)
             call fock_getghc(cwavef(:,1+(idat-1)*npw_k1*my_nspinor:idat*npw_k1*my_nspinor),cwaveprj_idat,&
&             gvnlxc(:,1+(idat-1)*npw_k2*my_nspinor:idat*npw_k2*my_nspinor),gs_ham,mpi_enreg)
           else
             call fock_ACE_getghc(cwavef(:,1+(idat-1)*npw_k1*my_nspinor:idat*npw_k1*my_nspinor),&
&             gvnlxc(:,1+(idat-1)*npw_k2*my_nspinor:idat*npw_k2*my_nspinor),gs_ham,mpi_enreg)
           end if
         end do ! idat
       end if
     end if

   else if (type_calc == 3) then ! for kinetic and local only, nonlocal and vfock should be zero

     gvnlxc(:,:) = zero

   end if ! if(type_calc...

!============================================================
! Assemble kinetic, local, nonlocal and Fock contributions
!============================================================

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
               ghc(re,igspinor)= ghc(re,igspinor) + kinpw_k2(ig)*cwavef(re,igspinor) + gvnlxc(re,igspinor)
               ghc(im,igspinor)= ghc(im,igspinor) + kinpw_k2(ig)*cwavef(im,igspinor) + gvnlxc(im,igspinor)
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
               ghc(re,igspinor)= ghc(re,igspinor) + gvnlxc(re,igspinor)
               ghc(im,igspinor)= ghc(im,igspinor) + gvnlxc(im,igspinor)
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
               ghcre=kinpw_k2(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlxc(re,igspinor)
               ghcim=kinpw_k2(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlxc(im,igspinor)
             else
               ghcre=ghc(re,igspinor)+gvnlxc(re,igspinor)
               ghcim=ghc(im,igspinor)+gvnlxc(im,igspinor)
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
&             kinpw_k2(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlxc(re,igspinor), gsc(re,igspinor)
             call wrtout(std_out,msg,'PERS')
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw_k2(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlxc(im,igspinor), gsc(im,igspinor)
             call wrtout(std_out,msg,'PERS')
           else
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw_k2(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlxc(re,igspinor)
             call wrtout(std_out,msg,'PERS')
             write(msg,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw_k2(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlxc(im,igspinor)
             call wrtout(std_out,msg,'PERS')
           end if
           ghc(re,igspinor)=ghcre
           ghc(im,igspinor)=ghcim
         end do ! ig
       end do ! ispinor
     end do ! idat
   end if

!  Special case of PAW + Fock : only return Fock operator contribution in gvnlxc
   if (gs_ham%usepaw==1 .and. has_fock)then
     gvnlxc=gvnlxc-gvnlc
     ABI_DEALLOCATE(gvnlc)
   endif

!  Structured debugging : if prtvol=-level, stop here.
   if(prtvol==-level)then
     write(msg,'(a,i0,a)')' getghc : exit prtvol=-',level,', debugging mode => stop '
     MSG_ERROR(msg)
   end if

   if ((type_calc==0).or.(type_calc==2)) then
     if (has_fock.and.gs_ham%usepaw==1.and.cpopt<2) then
       call pawcprj_free(cwaveprj_fock)
       ABI_DATATYPE_DEALLOCATE(cwaveprj_fock)
     end if
   end if

 end if ! type_calc

 call timab(200+tim_getghc,2,tsec)

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
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LSI, MT, JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cwavef(2,npw_k*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! gbound_k(2*mgfft+4)=sphere boundary info
!! gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!! istwf_k=input parameter that describes the storage of wfs
!! kg_k(3,npw_k)=G vec coordinates wrt recip lattice transl.
!! kpt(3)=current k point
!! mgfft=maximum single fft dimension
!! mpi_enreg=informations about MPI parallelization
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
!! PARENTS
!!      m_getghc
!!
!! CHILDREN
!!      getghc,mkl_set_num_threads,omp_set_nested
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

 ABI_ALLOCATE(work,(2,n4,n5,n6*ndat))

 ! scale conversion from SI to atomic units,
 ! here \alpha^2 where \alpha is the fine structure constant
 scale_conversion = FineStructureConstant2
 ! JWZ debug
 ! scale_conversion = zero

 if (nspinortot==1) then

    ABI_ALLOCATE(ghc1,(2,npw_k*ndat))

    !  Do it in 2 STEPs:
    !  STEP1: Compute grad of cwavef
    ABI_ALLOCATE(gcwavef,(2,npw_k*ndat,3))

    gcwavef = zero

    ! compute k + G. Note these are in reduced coords
    ABI_ALLOCATE(kgkpk,(npw_k,3))
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
    ABI_DEALLOCATE(kgkpk)
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
    ABI_DEALLOCATE(gcwavef)
    ABI_DEALLOCATE(ghc1)

    ! JWZ debug blank this term for now
    ! ghc_vectornd = zero

 else ! nspinortot==2

    ABI_ALLOCATE(cwavef1,(2,npw_k*ndat))
    ABI_ALLOCATE(cwavef2,(2,npw_k*ndat))
    do idat=1,ndat
       do ipw=1,npw_k
          cwavef1(1:2,ipw+(idat-1)*npw_k)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k)
          cwavef2(1:2,ipw+(idat-1)*npw_k)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k+shift)
       end do
    end do

    ! compute k + G. Note these are in reduced coords
    ABI_ALLOCATE(kgkpk,(npw_k,3))
    do ipw = 1, npw_k
       kgkpk(ipw,:) = kpt(:) + kg_k(:,ipw)
    end do

    if (nspinor1TreatedByThisProc) then

       ABI_ALLOCATE(ghc1,(2,npw_k*ndat))

       !  Do it in 2 STEPs:
       !  STEP1: Compute grad of cwavef
       ABI_ALLOCATE(gcwavef1,(2,npw_k*ndat,3))

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
       ABI_DEALLOCATE(gcwavef1)
       ABI_DEALLOCATE(ghc1)

    end if ! end spinor 1

    if (nspinor2TreatedByThisProc) then

       ABI_ALLOCATE(ghc2,(2,npw_k*ndat))

       !  Do it in 2 STEPs:
       !  STEP1: Compute grad of cwavef
       ABI_ALLOCATE(gcwavef2,(2,npw_k*ndat,3))

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
       ABI_DEALLOCATE(gcwavef2)
       ABI_DEALLOCATE(ghc2)

    end if ! end spinor 2

    ABI_DEALLOCATE(cwavef1)
    ABI_DEALLOCATE(cwavef2)
    ABI_DEALLOCATE(kgkpk)

 end if ! nspinortot

 ABI_DEALLOCATE(work)

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
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LSI, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cwavef(2,npw_k*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! gbound_k(2*mgfft+4)=sphere boundary info
!! istwf_k=input parameter that describes the storage of wfs
!! kg_k(3,npw_k)=G vec coordinates wrt recip lattice transl.
!! kpt(3)=current k point
!! mgfft=maximum single fft dimension
!! mpi_enreg=informations about MPI parallelization
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
!! PARENTS
!!      m_getghc
!!
!! CHILDREN
!!      getghc,mkl_set_num_threads,omp_set_nested
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

 ABI_ALLOCATE(work,(2,n4,n5,n6*ndat))

 if (nspinortot==1) then

   ABI_ALLOCATE(ghc1,(2,npw_k*ndat))

!  Do it in 3 STEPs:
!  STEP1: Compute grad of cwavef and Laplacian of cwavef
   ABI_ALLOCATE(gcwavef,(2,npw_k*ndat,3))
   ABI_ALLOCATE(lcwavef,(2,npw_k*ndat))
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
   ABI_DEALLOCATE(lcwavef)
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
   ABI_DEALLOCATE(gcwavef)
   ABI_DEALLOCATE(ghc1)

 else ! nspinortot==2

   ABI_ALLOCATE(cwavef1,(2,npw_k*ndat))
   ABI_ALLOCATE(cwavef2,(2,npw_k*ndat))
   do idat=1,ndat
     do ipw=1,npw_k
       cwavef1(1:2,ipw+(idat-1)*npw_k)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k)
       cwavef2(1:2,ipw+(idat-1)*npw_k)=cwavef(1:2,ipw+(idat-1)*my_nspinor*npw_k+shift)
     end do
   end do
!  call cg_zcopy(npw*ndat,cwavef(1,1),cwavef1)
!  call cg_zcopy(npw*ndat,cwavef(1,1+shift),cwavef2)


   if (nspinor1TreatedByThisProc) then

     ABI_ALLOCATE(ghc1,(2,npw_k*ndat))

!    Do it in 3 STEPs:
!    STEP1: Compute grad of cwavef and Laplacian of cwavef
     ABI_ALLOCATE(gcwavef1,(2,npw_k*ndat,3))
     ABI_ALLOCATE(lcwavef1,(2,npw_k*ndat))
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
     ABI_DEALLOCATE(lcwavef1)
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
     ABI_DEALLOCATE(gcwavef1)
     ABI_DEALLOCATE(ghc1)

   end if ! spin 1 treated by this proc

   if (nspinor2TreatedByThisProc) then

     ABI_ALLOCATE(ghc2,(2,npw_k*ndat))

!    Do it in 3 STEPs:
!    STEP1: Compute grad of cwavef and Laplacian of cwavef
     ABI_ALLOCATE(gcwavef2,(2,npw_k*ndat,3))
     ABI_ALLOCATE(lcwavef2,(2,npw_k*ndat))
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
     ABI_DEALLOCATE(lcwavef2)
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

     ABI_DEALLOCATE(gcwavef2)
     ABI_DEALLOCATE(ghc2)

   end if ! spin 2 treated by this proc

   ABI_DEALLOCATE(cwavef1)
   ABI_DEALLOCATE(cwavef2)

 end if ! npsinortot

 ABI_DEALLOCATE(work)

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
!! PARENTS
!!      m_dfpt_vtowfk,m_dfptnl_pert
!!
!! CHILDREN
!!      getghc,mkl_set_num_threads,omp_set_nested
!!
!! SOURCE

subroutine getgsc(cg,cprj,gs_ham,gsc,ibg,icg,igsc,ikpt,isppol,&
&                 mcg,mcprj,mgsc,mpi_enreg,natom,nband,npw_k,nspinor,select_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ibg,icg,igsc,ikpt,isppol,mcg,mcprj
 integer,intent(in) :: mgsc,natom,nband,npw_k,nspinor
 integer,intent(in),optional :: select_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: gsc(2,mgsc)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimenl1,dimenl2,iband,iband1,iband2,ierr,index_cg,index_cprj
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
   MSG_BUG('Only compatible with PAW (usepaw=1) !')
 end if
 if(nband<0.and.(mcg<npw_k*my_nspinor.or.mgsc<npw_k*my_nspinor.or.mcprj<my_nspinor)) then
   MSG_BUG('Invalid value for mcg, mgsc or mcprj !')
 end if

!Keep track of total time spent in getgsc:
 call timab(565,1,tsec)

 gsc = zero

!Prepare some data
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(scwavef,(2,npw_k*my_nspinor))
 if (gs_ham%usecprj==1) then
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor))
   call pawcprj_alloc(cwaveprj,0,gs_ham%dimcprj)
 else
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
 end if
 dimenl1=gs_ham%dimekb1;dimenl2=natom;tim_nonlop=0
 choice=1;signs=2;cpopt=-1+3*gs_ham%usecprj;paw_opt=3;useylm=1
 select_k_=1;if (present(select_k)) select_k_=select_k
 me=mpi_enreg%me_kpt

!Loop over bands
 index_cprj=ibg;index_cg=icg;index_gsc=igsc
 if (nband>0) then
   iband1=1;iband2=nband
 else if (nband<0) then
   iband1=abs(nband);iband2=iband1
   index_cprj=index_cprj+(iband1-1)*my_nspinor
   index_cg  =index_cg  +(iband1-1)*npw_k*my_nspinor
   index_gsc =index_gsc +(iband1-1)*npw_k*my_nspinor
 end if

 do iband=iband1,iband2

   if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me.and.nband>0) then
     gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=zero
     index_cprj=index_cprj+my_nspinor
     index_cg=index_cg+npw_k*my_nspinor
     index_gsc=index_gsc+npw_k*my_nspinor
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

!Reduction in case of parallelization
 if ((xmpi_paral==1)) then
   call timab(48,1,tsec)
   call xmpi_sum(gsc,mpi_enreg%comm_band,ierr)
   call timab(48,2,tsec)
 end if

!Memory deallocation
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(scwavef)
 if (gs_ham%usecprj==1) then
   call pawcprj_free(cwaveprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

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
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
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
!!      m_lobpcgwf,m_prep_kgb
!!
!! CHILDREN
!!      getghc,mkl_set_num_threads,omp_set_nested
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
 integer :: firstelt, lastelt
 integer :: nthreads
 integer :: ithread
 integer :: chunk
 integer :: residuchunk
 integer :: firstband
 integer :: lastband
 integer :: spacedim
#ifdef HAVE_OPENMP
 logical :: is_nested
#endif

 integer :: select_k_default

 ! *************************************************************************

 select_k_default = 1; if ( present(select_k) ) select_k_default = select_k

 spacedim = size(cwavef,dim=2)/ndat

    !$omp parallel default (none) private(ithread,nthreads,chunk,firstband,lastband,residuchunk,firstelt,lastelt, is_nested), &
    !$omp& shared(cwavef,ghc,gsc, gvnlxc,spacedim,ndat,kg_fft_k,kg_fft_kp,gs_ham,cwaveprj,mpi_enreg), &
    !$omp& firstprivate(cpopt,lambda,prtvol,sij_opt,tim_getghc,type_calc,select_k_default)
#ifdef HAVE_OPENMP
 ithread = omp_get_thread_num()
 nthreads = omp_get_num_threads()
 is_nested = omp_get_nested()
 call omp_set_nested(.false.)
#ifdef HAVE_LINALG_MKL_THREADS
 call mkl_set_num_threads(1)
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

 if ( lastband /= 0 ) then
   firstelt = (firstband-1)*spacedim+1
   lastelt = lastband*spacedim
      ! Don't know how to manage optional arguments .... :(
   if ( present(kg_fft_k) ) then
     if (present(kg_fft_kp)) then
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
       select_k=select_k_default,kg_fft_k=kg_fft_k,kg_fft_kp=kg_fft_kp)
     else
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
       select_k=select_k_default,kg_fft_k=kg_fft_k)
     end if
   else
     if (present(kg_fft_kp)) then
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
       select_k=select_k_default,kg_fft_kp=kg_fft_kp)
     else
       call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
       gs_ham,gvnlxc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
       select_K=select_k_default)
     end if
   end if
 end if
#ifdef HAVE_OPENMP
 call omp_set_nested(is_nested)
#ifdef HAVE_LINALG_MKL_THREADS
 call mkl_set_num_threads(nthreads)
#endif
#endif
    !$omp end parallel

end subroutine multithreaded_getghc
!!***

! !!****f* ABINIT/getghcnd
! !!
! !! NAME
! !! getghcnd
! !!
! !! FUNCTION
! !! Compute <G|H_ND|C> for input vector |C> expressed in reciprocal space
! !! Result is put in array ghcnc. H_ND is the Hamiltonian due to magnetic dipoles
! !! on the nuclear sites.
! !!
! !! INPUTS
! !! cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
! !! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
! !! my_nspinor=number of spinorial components of the wavefunctions (on current proc)
! !! ndat=number of FFT to do in //
! !!
! !! OUTPUT
! !! ghcnd(2,npw*my_nspinor*ndat)=matrix elements <G|H_ND|C>
! !!
! !! NOTES
! !!  This routine applies the Hamiltonian due to an array of magnetic dipoles located
! !!  at the atomic nuclei to the input wavefunction. Strategy below is to take advantage of
! !!  Hermiticity to store H_ND in triangular form and then use a BLAS call to zhpmv to apply to
! !!  input vector in one shot.
! !! Application of <k^prime|H|k> or <k|H|k^prime> not implemented!
! !!
! !! PARENTS
! !!      getghc
! !!
! !! CHILDREN
! !!      zhpmv
! !!
! !! SOURCE

! subroutine getghcnd(cwavef,ghcnd,gs_ham,my_nspinor,ndat)

! !Arguments ------------------------------------
! !scalars
!  integer,intent(in) :: my_nspinor,ndat
!  type(gs_hamiltonian_type),intent(in) :: gs_ham
! !arrays
!  real(dp),intent(in) :: cwavef(2,gs_ham%npw_k*my_nspinor*ndat)
!  real(dp),intent(out) :: ghcnd(2,gs_ham%npw_k*my_nspinor*ndat)

! !Local variables-------------------------------
! !scalars
!  integer :: cwavedim
!  character(len=500) :: message
!  !arrays
!  complex(dpc),allocatable :: inwave(:),hggc(:)

! ! *********************************************************************

!  if (gs_ham%matblk /= gs_ham%natom) then
!    write(message,'(a,i4,a,i4)')' gs_ham%matblk = ',gs_ham%matblk,' but natom = ',gs_ham%natom
!    MSG_ERROR(message)
!  end if
!  if (ndat /= 1) then
!    write(message,'(a,i4,a)')' ndat = ',ndat,' but getghcnd requires ndat = 1'
!    MSG_ERROR(message)
!  end if
!  if (my_nspinor /= 1) then
!    write(message,'(a,i4,a)')' nspinor = ',my_nspinor,' but getghcnd requires nspinor = 1'
!    MSG_ERROR(message)
!  end if
!  if (any(abs(gs_ham%kpt_k(:)-gs_ham%kpt_kp(:))>tol8)) then
!    message=' not allowed for kpt(left)/=kpt(right)!'
!    MSG_BUG(message)
!  end if

!  if (any(abs(gs_ham%nucdipmom_k)>tol8)) then
!    cwavedim = gs_ham%npw_k*my_nspinor*ndat
!    ABI_ALLOCATE(hggc,(cwavedim))
!    ABI_ALLOCATE(inwave,(cwavedim))

!    inwave(1:gs_ham%npw_k) = cmplx(cwavef(1,1:gs_ham%npw_k),cwavef(2,1:gs_ham%npw_k),kind=dpc)

!     ! apply hamiltonian hgg to input wavefunction inwave, result in hggc
!     ! ZHPMV is a level-2 BLAS routine, does Matrix x Vector multiplication for double complex
!     ! objects, with the matrix as Hermitian in packed storage
!    call ZHPMV('L',cwavedim,cone,gs_ham%nucdipmom_k,inwave,1,czero,hggc,1)

!    ghcnd(1,1:gs_ham%npw_k) = real(hggc)
!    ghcnd(2,1:gs_ham%npw_k) = aimag(hggc)

!    ABI_DEALLOCATE(hggc)
!    ABI_DEALLOCATE(inwave)

!  else

!    ghcnd(:,:) = zero

!  end if

! end subroutine getghcnd
! !!***

end module m_getghc
!!***
