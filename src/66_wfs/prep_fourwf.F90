!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_fourwf
!! NAME
!! prep_fourwf
!!
!! FUNCTION
!! this routine prepares the data to the call of fourwf.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FBottin,MT,GZ,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cwavef(2,npw*ndat)=planewave coefficients of wavefunction (one spinorial component?).
!!  dtfil <type(datafiles_type)>=variables related to files
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=informations about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndat=number of FFT to do in //
!!  ngfft(18)= contain all needed information about 3D FFT
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  n4,n5,n6 used for dimensionning of vlocal
!!  option_fourwf=option for fourwf (see fourwf.F90)
!!  prtvol=control print volume and debugging output
!!  ucvol=unit cell volume
!!  [bandfft_kpt_tab]= (optional) if present, contains tabs used to implement
!!                     the "band-fft" parallelism
!!                      if not present, the bandfft_kpt global variable is used
!!  [use_gpu_cuda]= (optional) 1 if Cuda (GPU) is on
!!
!! OUTPUT
!!  gwavef=(2,npw*ndat)=matrix elements <G|H|C>.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mkrho,posdoppler,vtowfk
!!
!! CHILDREN
!!      fourwf,gpu_fourwf,prep_index_wavef_bandpp,prep_wavef_sym_do
!!      prep_wavef_sym_undo,timab,xmpi_alltoallv,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,iblock,istwf_k,mgfft,&
&          mpi_enreg,nband_k,ndat,ngfft,npw_k,n4,n5,n6,occ_k,option_fourwf,ucvol,wtk,&
&          bandfft_kpt_tab,use_gpu_cuda) ! Optional arguments

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_time,        only : timab
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_bandfft_kpt, only : bandfft_kpt_type,bandfft_kpt,bandfft_kpt_get_ikpt

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_fourwf'
 use interfaces_53_ffts
 use interfaces_66_wfs, except_this_one => prep_fourwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,iblock,istwf_k,mgfft,n4,n5,n6,nband_k,ndat,npw_k
 integer,intent(in) :: option_fourwf
 integer,intent(in),optional :: use_gpu_cuda
 real(dp),intent(in) :: ucvol,wtk
 type(bandfft_kpt_type),optional,target,intent(in) :: bandfft_kpt_tab
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: occ_k(nband_k)
 real(dp),intent(out) :: rhoaug(n4,n5,n6)
 real(dp),intent(in) :: cwavef(2,npw_k*blocksize)
 real(dp),target,intent(inout) :: wfraug(2,n4,n5,n6*ndat)

!Local variables-------------------------------
!scalars
 integer :: bandpp,bandpp_sym,ier,iibandpp,ikpt_this_proc,ind_occ,ind_occ1,ind_occ2,ipw
 integer :: istwf_k_,jjbandpp,me_fft,nd3,nproc_band,nproc_fft,npw_fft
 integer :: old_me_g0=0,spaceComm=0,tim_fourwf,use_gpu_cuda_
 integer,pointer :: idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 logical :: flag_inv_sym,have_to_reequilibrate
 real(dp) :: weight,weight1,weight2
 type(bandfft_kpt_type),pointer :: bandfft_kpt_ptr
!arrays
 integer,ABI_CONTIGUOUS pointer :: indices_pw_fft(:),kg_k_fft(:,:),kg_k_gather(:,:),kg_k_gather_sym(:,:)
 integer,ABI_CONTIGUOUS pointer :: rdispls(:),rdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: recvcounts(:),recvcount_fft(:),recvcounts_sym(:),recvcounts_sym_tot(:)
 integer,ABI_CONTIGUOUS pointer :: recvdisp_fft(:),sdispls(:),sdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: sendcounts(:),sendcount_fft(:),sendcounts_sym(:),sendcounts_sym_all(:)
 integer,ABI_CONTIGUOUS pointer :: senddisp_fft(:),tab_proc(:)
 integer,allocatable :: rdisplsloc(:)
 integer,allocatable :: recvcountsloc(:),sdisplsloc(:)
 integer,allocatable :: sendcountsloc(:)
 integer,allocatable :: index_wavef_band(:),index_wavef_send(:)
 integer,pointer :: gbound_(:,:)
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: buff_wf(:,:),cwavef_alltoall1(:,:),cwavef_alltoall2(:,:)
 real(dp),allocatable :: cwavef_fft(:,:), cwavef_fft_tr(:,:)
 real(dp),allocatable :: weight_t(:),weight1_t(:),weight2_t(:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:),wfraug_ptr(:,:,:,:)

! *************************************************************************

 ABI_CHECK((option_fourwf/=3),'Option=3 (FFT r->g) not implemented')
 ABI_CHECK((mpi_enreg%bandpp==ndat),'BUG: bandpp/=ndat')

 spaceComm=mpi_enreg%comm_band
 nproc_band = mpi_enreg%nproc_band
 nproc_fft  = mpi_enreg%nproc_fft
 bandpp     = mpi_enreg%bandpp
 me_fft     = mpi_enreg%me_fft

 use_gpu_cuda_=0;if (present(use_gpu_cuda)) use_gpu_cuda_=use_gpu_cuda

 if (present(bandfft_kpt_tab)) then
   bandfft_kpt_ptr => bandfft_kpt_tab
 else
   ikpt_this_proc=bandfft_kpt_get_ikpt()
   bandfft_kpt_ptr => bandfft_kpt(ikpt_this_proc)
 end if

 istwf_k_=istwf_k
 flag_inv_sym = (istwf_k_==2 .and. any(ngfft(7) == [401,402,312]))
 if (option_fourwf==0) flag_inv_sym=((flag_inv_sym).and.(use_gpu_cuda_==0))

 if (flag_inv_sym) then
   istwf_k_       = 1
   if (modulo(bandpp,2)==0) then
     bandpp_sym   = bandpp/2
   else
     bandpp_sym   = bandpp
   end if
 end if

!====================================================================================
 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 recvcounts   =>bandfft_kpt_ptr%recvcounts(:)
 sendcounts   =>bandfft_kpt_ptr%sendcounts(:)
 rdispls      =>bandfft_kpt_ptr%rdispls   (:)
 sdispls      =>bandfft_kpt_ptr%sdispls   (:)
 ndatarecv    =>bandfft_kpt_ptr%ndatarecv

 kg_k_gather  =>bandfft_kpt_ptr%kg_k_gather(:,:)
 gbound_      =>bandfft_kpt_ptr%gbound(:,:)

 if (flag_inv_sym ) then
   idatarecv0           =>bandfft_kpt_ptr%idatarecv0
   ndatarecv_tot        =>bandfft_kpt_ptr%ndatarecv_tot
   ndatasend_sym        =>bandfft_kpt_ptr%ndatasend_sym
   kg_k_gather_sym      =>bandfft_kpt_ptr%kg_k_gather_sym(:,:)
   rdispls_sym          =>bandfft_kpt_ptr%rdispls_sym(:)
   recvcounts_sym       =>bandfft_kpt_ptr%recvcounts_sym(:)
   recvcounts_sym_tot   =>bandfft_kpt_ptr%recvcounts_sym_tot(:)
   sdispls_sym          =>bandfft_kpt_ptr%sdispls_sym(:)
   sendcounts_sym       =>bandfft_kpt_ptr%sendcounts_sym(:)
   sendcounts_sym_all   =>bandfft_kpt_ptr%sendcounts_sym_all(:)
   tab_proc             =>bandfft_kpt_ptr%tab_proc(:)
 end if

 ABI_ALLOCATE(cwavef_alltoall2,(2,ndatarecv*bandpp))
 if ( ((.not.flag_inv_sym) .and. (bandpp>1) ) .or. flag_inv_sym )then
   ABI_ALLOCATE(cwavef_alltoall1,(2,ndatarecv*bandpp))
 end if

 recvcountsloc(:)=recvcounts(:)*2*bandpp
 rdisplsloc(:)=rdispls(:)*2*bandpp
 sendcountsloc(:)=sendcounts(:)*2
 sdisplsloc(:)=sdispls(:)*2

 call timab(547,1,tsec)
 call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall2,&
& recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(547,2,tsec)

!If me_fft==0, I have the G=0 vector, but keep for the record the old value
 if (me_fft==0) then
   old_me_g0=mpi_enreg%me_g0
   mpi_enreg%me_g0=1
 end if

 tim_fourwf=16

!Eventually adjust load balancing for FFT (by changing FFT distrib)
 have_to_reequilibrate = bandfft_kpt_ptr%have_to_reequilibrate
 if(have_to_reequilibrate) then
   npw_fft =  bandfft_kpt_ptr%npw_fft
   sendcount_fft  => bandfft_kpt_ptr%sendcount_fft(:)
   recvcount_fft  => bandfft_kpt_ptr%recvcount_fft(:)
   senddisp_fft   => bandfft_kpt_ptr%senddisp_fft(:)
   recvdisp_fft   => bandfft_kpt_ptr%recvdisp_fft(:)
   indices_pw_fft => bandfft_kpt_ptr%indices_pw_fft(:)
   kg_k_fft       => bandfft_kpt_ptr%kg_k_fft(:,:)
   ABI_ALLOCATE( buff_wf, (2,ndatarecv*bandpp) ) ! for sorting cgwavef
   ABI_ALLOCATE( cwavef_fft, (2,npw_fft*bandpp) )
   if(bandpp>1) then
     ABI_ALLOCATE( cwavef_fft_tr, (2,npw_fft*bandpp))
   end if
 end if

 if (option_fourwf==0) wfraug(:,:,:,:)=zero

!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then

!  Compute the index of the band
   ind_occ = (iblock-1)*blocksize + mpi_enreg%me_band + 1

   if(abs(occ_k(ind_occ))>=tol8.or.option_fourwf==0) then

!    Compute the weight of the band
     weight=occ_k(ind_occ)*wtk/ucvol

     if(have_to_reequilibrate) then
!      filling of sorted send buffers before exchange
       do ipw = 1 ,ndatarecv
         buff_wf(1:2, indices_pw_fft(ipw) ) = cwavef_alltoall2(1:2,ipw)
       end do
       call xmpi_alltoallv(buff_wf,2*sendcount_fft,2*senddisp_fft,  &
&       cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, mpi_enreg%comm_fft,ier)
       call fourwf(1,rhoaug,cwavef_fft,dummy,wfraug,gbound_,gbound_,&
&       istwf_k_,kg_k_fft,kg_k_fft,mgfft,mpi_enreg,1,&
&       ngfft,npw_fft,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=use_gpu_cuda_)
     else
       call fourwf(1,rhoaug,cwavef_alltoall2,dummy,wfraug,gbound_,gbound_,&
&       istwf_k_,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&       ngfft,ndatarecv,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=use_gpu_cuda_)
     end if
     if (option_fourwf==0.and.nproc_fft>1) then
       if (me_fft>0) then
         nd3=(ngfft(3)-1)/nproc_fft+1
         wfraug(:,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=wfraug(:,:,:,1:nd3)
         wfraug(:,:,:,1:nd3)=zero
       end if
       call xmpi_sum(wfraug,mpi_enreg%comm_fft,ier)
     end if
   end if

!====================================================================
 else if ((.not.(flag_inv_sym)) .and. (bandpp>1) ) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the wave functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall1(:,:) = cwavef_alltoall2(:,index_wavef_band)

   if(have_to_reequilibrate) then
!    filling of sorted send buffers before exchange
     do iibandpp=1,bandpp
       do ipw = 1 ,ndatarecv
         buff_wf(1:2, iibandpp + bandpp*(indices_pw_fft(ipw)-1)) = &
&         cwavef_alltoall1(1:2,ipw + ndatarecv*(iibandpp-1))
       end do
     end do
     call xmpi_alltoallv(buff_wf,2*bandpp*sendcount_fft,2*bandpp*senddisp_fft,  &
&     cwavef_fft_tr,2*bandpp*recvcount_fft, 2*bandpp*recvdisp_fft, mpi_enreg%comm_fft,ier)
     do iibandpp=1,bandpp
       do ipw = 1 ,npw_fft
         cwavef_fft(1:2,  ipw + npw_fft*(iibandpp -1)) = cwavef_fft_tr(1:2,  iibandpp + bandpp*(ipw-1))
       end do
     end do
   end if

!  -------------------
!  Fourier calculation
!  -------------------
!  Cuda version
   if(use_gpu_cuda_==1) then
     ABI_ALLOCATE(weight_t,(bandpp))
     do iibandpp=1,bandpp
!      Compute the index of the band
       ind_occ = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + iibandpp
!      Compute the weight of the band
       weight_t(iibandpp)=occ_k(ind_occ)*wtk/ucvol
       if(abs(occ_k(ind_occ)) < tol8) weight_t(iibandpp) = zero
     end do
!    Accumulate time because it is not done in gpu_fourwf
     call timab(240+tim_fourwf,1,tsec)
#if defined HAVE_GPU_CUDA
     call gpu_fourwf(1,rhoaug,&
&     cwavef_alltoall1,&
&     dummy,wfraug,gbound_,gbound_,&
&     istwf_k_,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,bandpp,&
&     ngfft,ndatarecv,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,&
&     tim_fourwf,weight_t,weight_t)
#endif
     call timab(240+tim_fourwf,2,tsec)
     ABI_DEALLOCATE(weight_t)

!  Standard version
   else
     do iibandpp=1,bandpp
!      Compute the index of the band
       ind_occ = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + iibandpp
!      Compute the weight of the band
       weight=occ_k(ind_occ)*wtk/ucvol
       if (option_fourwf==0) then
         wfraug_ptr => wfraug(:,:,:,(iibandpp-1)*n6+1:iibandpp*n6)
       else
         wfraug_ptr => wfraug
       end if
       if (abs(occ_k(ind_occ)) >=tol8.or.option_fourwf==0) then
         if(have_to_reequilibrate) then
           call fourwf(1,rhoaug, &
&           cwavef_fft(:,(npw_fft*(iibandpp-1))+1:(npw_fft*iibandpp)), &
&           dummy,wfraug_ptr,gbound_,gbound_,&
&           istwf_k_,kg_k_fft,kg_k_fft,mgfft,mpi_enreg,1,&
&           ngfft,npw_fft,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=use_gpu_cuda_)
         else
           call fourwf(1,rhoaug,&
&           cwavef_alltoall1(:,(ndatarecv*(iibandpp-1))+1:(ndatarecv*iibandpp)),&
&           dummy,wfraug_ptr,gbound_,gbound_,&
&           istwf_k_,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&           ngfft,ndatarecv,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,&
&           tim_fourwf,weight,weight)
         end if
         if (option_fourwf==0.and.nproc_fft>1) then
           if (me_fft>0) then
             nd3=(ngfft(3)-1)/nproc_fft+1
             wfraug_ptr(:,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=wfraug_ptr(:,:,:,1:nd3)
             wfraug_ptr(:,:,:,1:nd3)=zero
           end if
           call xmpi_sum(wfraug_ptr,mpi_enreg%comm_fft,ier)
         end if
       end if
     end do
   end if ! (use_gpu_cuda==1)

!  -----------------------------------------------------
!  Sorting waves functions below the processors
!  -----------------------------------------------------
!  cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)   ! NOT NEEDED
   ABI_DEALLOCATE(index_wavef_band)

!====================================================================
 else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting the wave functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall1(:,:) = cwavef_alltoall2(:,index_wavef_band)

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,bandpp,1,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,tab_proc,&
&   cwavef_alltoall1,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   ewavef_alltoall_sym,&
&   index_wavef_send)

!  ------------------------------------------------------------
!  Fourier calculation
!  ------------------------------------------------------------
!  Cuda version
   if (use_gpu_cuda_==1) then
     ABI_ALLOCATE(weight1_t,(bandpp_sym))
     ABI_ALLOCATE(weight2_t,(bandpp_sym))
     do iibandpp=1,bandpp_sym
       if (bandpp/=1) then
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp-1)
         ind_occ2 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp  )
       else
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + 1
         ind_occ2 = ind_occ1
       end if
       weight1_t(iibandpp) = occ_k(ind_occ1)*wtk/ucvol
       weight2_t(iibandpp) = occ_k(ind_occ2)*wtk/ucvol
     end do
     call timab(240+tim_fourwf,1,tsec)
#if defined HAVE_GPU_CUDA
     call gpu_fourwf(1,rhoaug,&
&     ewavef_alltoall_sym,&
&     dummy,wfraug,gbound_,gbound_,&
&     istwf_k_,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,bandpp_sym,&
&     ngfft,ndatarecv_tot,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,&
&     tim_fourwf,weight1_t,weight2_t)
#endif
     call timab(240+tim_fourwf,2,tsec)
     ABI_DEALLOCATE(weight1_t)
     ABI_DEALLOCATE(weight2_t)

!  Standard version
   else
     if (option_fourwf==0.and.bandpp>1) then
       ABI_ALLOCATE(wfraug_ptr,(2,n4,n5,n6))
     else
       wfraug_ptr => wfraug
     end if
     do iibandpp=1,bandpp_sym
       if (bandpp/=1) then
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp-1)
         ind_occ2 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp  )
       else
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + 1
         ind_occ2 = ind_occ1
       end if
       weight1 = occ_k(ind_occ1)*wtk/ucvol
       weight2 = occ_k(ind_occ2)*wtk/ucvol
       call fourwf(1,rhoaug,&
&       ewavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&       dummy,wfraug_ptr,gbound_,gbound_,&
&       istwf_k_,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,1,&
&       ngfft,ndatarecv_tot,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,&
&       tim_fourwf,weight1,weight2)
       if (option_fourwf==0) then
         if (modulo(bandpp,2)==0) then
           jjbandpp=2*iibandpp-1
           wfraug(1,:,:,(jjbandpp-1)*n6+1:jjbandpp*n6)=wfraug_ptr(1,:,:,1:n6)
           wfraug(1,:,:,(jjbandpp)*n6+1:(jjbandpp+1)*n6)=wfraug_ptr(2,:,:,1:n6)
         else if (bandpp>1) then
           wfraug(1,:,:,(iibandpp-1)*n6+1:iibandpp*n6)=wfraug_ptr(1,:,:,1:n6)
         end if
         if (nproc_fft>1) then
           if (me_fft>0) then
             nd3=(ngfft(3)-1)/nproc_fft+1
             wfraug(1,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=wfraug(1,:,:,1:nd3)
             wfraug(1,:,:,1:nd3)=zero
           end if
           call xmpi_sum(wfraug,mpi_enreg%comm_fft,ier)
         end if
       end if
     end do
     if (option_fourwf==0.and.bandpp>1) then
       ABI_DEALLOCATE(wfraug_ptr)
     end if
   end if ! (use_gpu_cuda==1)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,bandpp,1,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,idatarecv0,&
&   cwavef_alltoall1,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   ewavef_alltoall_sym,&
&   index_wavef_send)

   ABI_DEALLOCATE(ewavef_alltoall_sym)
   ABI_DEALLOCATE(index_wavef_send)

!  -------------------------------------------------------
!  Sorting waves functions below the processors
!  -------------------------------------------------------
!  cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:) ! NOT NEEDED

   ABI_DEALLOCATE(index_wavef_band)

 end if

!====================================================================
 if (me_fft==0) mpi_enreg%me_g0=old_me_g0
 if(have_to_reequilibrate) then
   ABI_DEALLOCATE(buff_wf)
   ABI_DEALLOCATE(cwavef_fft)
   if(bandpp > 1) then
     ABI_DEALLOCATE(cwavef_fft_tr)
   end if
 end if
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall2)
 if ( ((.not.flag_inv_sym) .and. (bandpp>1) ) .or. flag_inv_sym ) then
   ABI_DEALLOCATE(cwavef_alltoall1)
 end if

end subroutine prep_fourwf
!!***


!!****f* ABINIT/prep_wavef_sym_do
!! NAME
!! prep_wavef_sym_do
!!
!! FUNCTION
!! this routine associates waves functions by two as following
!!      E(G)  = C(G) + D(G)
!!      E(-G) = C*(G) + iD*(G)
!! the values are distributed on the processors in function of
!! the value of mpi_enreg%distribfft%tab_fftwf2_distrib( (-kg_k_gather(2,i) )
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg          = informations about mpi parallelization
!!  bandpp             = number of couple of waves functions
!!  nspinor            = number of spin
!!  ndatarecv          = number of values received by the processor and sended
!!                       by the other processors band
!!  ndatarecv_tot      = total number of received values
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited
!!                       planewave coordinates
!!  tab_proc           = positions of opposited planewave coordinates in the list of the
!!                       processors fft
!!  cwavef_alltoall    = planewave coefficients of wavefunction
!!                      ( initial of the processor + sended by other processors band)
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor from each processor fft
!!  rdispls_sym        = postions of the received values by the processor from each processor fft
!!
!! OUTPUT
!!  ewavef_alltoall_sym = planewave coefficients of wavefunction
!!                        initial of the processor +
!!                        sended by other processors band +
!!                        sended by other processors fft  +
!!                        and compisited if bandpp >1
!!  index_wavef_send    = index to send the values in blocks to the other processor fft
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prep_fourwf,prep_getghc
!!
!! CHILDREN
!!      xmpi_alltoallv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,tab_proc,&
&     cwavef_alltoall,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     ewavef_alltoall_sym,&
&     index_wavef_send)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_wavef_sym_do'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,intent(in) :: nspinor
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,allocatable,intent(out) :: index_wavef_send(:)
 integer,intent(in) :: rdispls_sym(:),recvcounts_sym(:)
 integer,intent(in) :: sdispls_sym(:),sendcounts_sym(:)
 integer,intent(in) :: tab_proc(:)
 real(dp),intent(inout) :: cwavef_alltoall(2,ndatarecv*nspinor*bandpp)
 real(dp),pointer :: ewavef_alltoall_sym(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp_sym,ibandpp,idatarecv,ideb_loc,idebc,idebd,idebe
 integer :: ier,ifin_loc,ifinc,ifind,ifine,iproc,jbandpp,jsendloc
 integer :: kbandpp,newspacecomm,nproc_fft
 logical :: flag_compose
!arrays
 integer,allocatable :: rdispls_sym_loc(:),recvcounts_sym_loc(:)
 integer,allocatable :: sdispls_sym_loc(:),sendcounts_sym_loc(:)
 real(dp),allocatable :: ewavef_alltoall_loc(:,:),ewavef_alltoall_send(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_wavef_sym_do : enter '
!ENDDEBUG

!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft

 newspacecomm = mpi_enreg%comm_fft

 if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
   flag_compose = .TRUE.
 else
   bandpp_sym   = bandpp
   flag_compose = .FALSE.
 end if

!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(ewavef_alltoall_sym     ,(2,ndatarecv_tot*bandpp_sym))
 ABI_ALLOCATE(ewavef_alltoall_loc     ,(2,ndatarecv    *bandpp_sym))
 ABI_ALLOCATE(ewavef_alltoall_send    ,(2,ndatasend_sym*bandpp_sym))
 ABI_ALLOCATE(index_wavef_send        ,(  ndatasend_sym*bandpp_sym))

 ABI_ALLOCATE(sendcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(sdispls_sym_loc       ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(rdispls_sym_loc       ,(nproc_fft))


!Initialisation
!--------------
 ewavef_alltoall_sym(:,:) =0.
 ewavef_alltoall_loc(:,:) =0.

 sendcounts_sym_loc(:) =0
 sdispls_sym_loc(:)    =0
 recvcounts_sym_loc(:) =0
 rdispls_sym_loc(:)    =0

 index_wavef_send(:)   =0


!-------------------------------------------------
!We are bandpp blocks which we want to :
!associate by two      (band_sym==bandpp/2)
!or not associate by two  (band_sym==bandpp)
!
!So We'll have got bandpp_sym blocks
!So we loop on the bandpp_sym blocks
!--------------------------------------------------

 do kbandpp=1,bandpp_sym

!  position of the two blocks
!  --------------------------
   ibandpp = (kbandpp-1) * 2
   jbandpp =  ibandpp    + 1

   idebe = (kbandpp-1) * ndatarecv_tot + 1
   ifine = idebe       + ndatarecv     - 1

   idebc = ibandpp * ndatarecv     + 1
   ifinc = idebc   + ndatarecv     - 1

   idebd = jbandpp * ndatarecv     + 1
   ifind = idebd   + ndatarecv     - 1

   ideb_loc = (kbandpp-1) * ndatarecv  + 1
   ifin_loc = ideb_loc    + ndatarecv  - 1


   if (flag_compose) then


!    calcul ewavef(G)
!    ----------------
     ewavef_alltoall_sym(1,idebe:ifine) =    &
&     cwavef_alltoall(1,idebc:ifinc) &
&     - cwavef_alltoall(2,idebd:ifind)

     ewavef_alltoall_sym(2,idebe:ifine) =    &
&     cwavef_alltoall(2,idebc:ifinc) &
&     + cwavef_alltoall(1,idebd:ifind)

!    calcul ewavef_loc(-G)
!    ---------------------
     ewavef_alltoall_loc(1,ideb_loc:ifin_loc) =  &
&     cwavef_alltoall(1,idebc:ifinc) &
&     + cwavef_alltoall(2,idebd:ifind)

     ewavef_alltoall_loc(2,ideb_loc:ifin_loc) =  &
&     - cwavef_alltoall(2,idebc:ifinc) &
&     + cwavef_alltoall(1,idebd:ifind)
   else

!    calcul ewavef(G)
!    ----------------
     ewavef_alltoall_sym(1,idebe:ifine)   = cwavef_alltoall(1,idebc:ifinc)
     ewavef_alltoall_sym(2,idebe:ifine)   = cwavef_alltoall(2,idebc:ifinc)

!    calcul ewavef_loc(-G)
!    ---------------------
     ewavef_alltoall_loc(1,ideb_loc:ifin_loc) =   cwavef_alltoall(1,idebc:ifinc)
     ewavef_alltoall_loc(2,ideb_loc:ifin_loc) = - cwavef_alltoall(2,idebc:ifinc)

   end if

 end do



!------------------------------------------------------------------------
!Creation of datas blocks for each processor fft from ewavef_alltoall_loc
!to send datas by blocks with a alltoall...
!------------------------------------------------------------------------

!Position of the blocks
!----------------------
 jsendloc=0
 do ibandpp=1,bandpp_sym
   do iproc=1,nproc_fft
     do idatarecv=1,ndatarecv
       if (tab_proc(idatarecv)==(iproc-1)) then
         jsendloc=jsendloc+1
         index_wavef_send(jsendloc)  = idatarecv + ndatarecv * (ibandpp-1)
       end if
     end do
   end do
 end do

!Classment
!----------
 ewavef_alltoall_send(:,:)=ewavef_alltoall_loc(:,index_wavef_send)


!-------------------------------------------------
!Calcul of the number of received and sended datas
!-------------------------------------------------
 sendcounts_sym_loc = sendcounts_sym*2
 recvcounts_sym_loc = recvcounts_sym*2

!------------------------------------------
!Exchange of the datas ewavef_allto_all_loc
!------------------------------------------
 do ibandpp=1,bandpp_sym

!  ------------------------------------------------
!  Deplacment of the sended datas because of bandpp
!  ------------------------------------------------
   sdispls_sym_loc(:) = sdispls_sym(:) + ndatasend_sym * (ibandpp-1)
   sdispls_sym_loc    = sdispls_sym_loc   *2

!  --------------------------------------------------
!  Deplacment of the received datas because of bandpp
!  --------------------------------------------------
   rdispls_sym_loc(:) = rdispls_sym(:) + ndatarecv_tot * (ibandpp-1)
   rdispls_sym_loc    = rdispls_sym_loc   *2


   call xmpi_alltoallv(&
&   ewavef_alltoall_send(:,:) ,sendcounts_sym_loc,sdispls_sym_loc,&
&   ewavef_alltoall_sym(:,:)  ,recvcounts_sym_loc,rdispls_sym_loc,&
&   newspacecomm,ier)

 end do

!-----------------------
!Desallocation
!-----------------------

 ABI_DEALLOCATE(sendcounts_sym_loc)
 ABI_DEALLOCATE(recvcounts_sym_loc)
 ABI_DEALLOCATE(sdispls_sym_loc)
 ABI_DEALLOCATE(rdispls_sym_loc)

 ABI_DEALLOCATE(ewavef_alltoall_loc)
 ABI_DEALLOCATE(ewavef_alltoall_send)

end subroutine prep_wavef_sym_do
!!***


!!****f* ABINIT/prep_wavef_sym_undo
!! NAME
!! prep_wavef_sym_undo
!!
!! FUNCTION
!! this routine dissociates each wave function in two waves functions as following
!!      C(G) =   ( E*(-G) + E(G))/2
!!      D(G) = i*( E*(-G) - E(G))/2
!! the values are redistributed on the processors in function of
!! the value of mpi_enreg%distribfft%tab_fftwf2_distrib( (-kg_k_gather(2,i) )
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg          = informations about mpi parallelization
!!  bandpp             = number of groups of couple of waves functions
!!  nspinor            = number of spin
!!  ndatarecv          = number of values received by the processor and sended
!!                       by the other processors band
!!  ndatarecv_tot      = total number of received values
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited
!!                       planewave coordinates
!!  idatarecv0         = position of the planewave coordinates (0,0,0)
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor to each processor fft
!!!  rdispls_sym        = postions of the received values by the processor to each processor fft
!!
!!  gwavef_alltoall_sym = planewave coefficients of wavefunction
!!                        initial of the processor +
!!                        sended by other processors band +
!!                        sended by other processors fft  +
!!                        and composited if bandpp >1
!!  index_wavef_send    = index to send the values by block to the other processor fft
!!
!! OUTPUT
!!  gwavef_alltoall     = planewave coefficients of wavefunction
!!                        ( for of the processor + to send to other processors band)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prep_fourwf,prep_getghc
!!
!! CHILDREN
!!      xmpi_alltoallv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,idatarecv0,&
&     gwavef_alltoall,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     gwavef_alltoall_sym,&
&     index_wavef_send)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_wavef_sym_undo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,intent(in) :: nspinor
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: index_wavef_send(:),rdispls_sym(:),recvcounts_sym(:)
 integer,intent(in) :: sdispls_sym(:),sendcounts_sym(:)
 real(dp),intent(inout) :: gwavef_alltoall(2,ndatarecv*nspinor*bandpp)
 real(dp),intent(inout) :: gwavef_alltoall_sym(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp_sym,ibandpp,ideb_loc,idebc,idebd
 integer :: idebe,ier,ifin_loc,ifinc,ifind,ifine
 integer :: jbandpp,kbandpp,newspacecomm,nproc_fft
 logical :: flag_compose
!arrays
 integer,allocatable :: rdispls_sym_loc(:),recvcounts_sym_loc(:)
 integer,allocatable :: sdispls_sym_loc(:),sendcounts_sym_loc(:)
 real(dp),allocatable :: gwavef_alltoall_loc(:,:),gwavef_alltoall_rcv(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_wavef_sym_undo : enter '
!ENDDEBUG


!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft

 newspacecomm = mpi_enreg%comm_fft

 if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
   flag_compose = .TRUE.
 else
   bandpp_sym   = bandpp
   flag_compose = .FALSE.
 end if

!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(gwavef_alltoall_loc     ,(2,ndatarecv     *bandpp_sym))
 ABI_ALLOCATE(gwavef_alltoall_rcv     ,(2,ndatasend_sym *bandpp_sym))

 ABI_ALLOCATE(sendcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(sdispls_sym_loc       ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(rdispls_sym_loc       ,(nproc_fft))


!---------------------------------------------
!Initialisation
!---------------------------------------------
 gwavef_alltoall_loc(:,:) =0.

 sendcounts_sym_loc(:) =0
 sdispls_sym_loc(:)    =0
 recvcounts_sym_loc(:) =0
 rdispls_sym_loc(:)    =0


!-------------------------------------------------
!Calcul of number of the sended and received datas
!-------------------------------------------------
 sendcounts_sym_loc = sendcounts_sym*2
 recvcounts_sym_loc = recvcounts_sym*2

!----------------------------------------------------
!Sending of the values
!----------------------------------------------------
 do ibandpp = 1,bandpp_sym

!  -------------------------------------------------
!  Deplacment of the sended values because of bandpp
!  -------------------------------------------------
   sdispls_sym_loc(:) = sdispls_sym(:) + ndatasend_sym * (ibandpp-1)
   sdispls_sym_loc    = sdispls_sym_loc   *2

!  ---------------------------------------------------
!  Deplacment of the received values because of bandpp
!  ---------------------------------------------------
   rdispls_sym_loc(:) = rdispls_sym(:) + ndatarecv_tot * (ibandpp-1)
   rdispls_sym_loc    = rdispls_sym_loc   *2


   call xmpi_alltoallv(&
&   gwavef_alltoall_sym(:,:) ,recvcounts_sym_loc,rdispls_sym_loc,&
&   gwavef_alltoall_rcv(:,:) ,sendcounts_sym_loc,sdispls_sym_loc,&
&   newspacecomm,ier)

 end do


!----------------------
!Dispatching the blocks
!----------------------
 gwavef_alltoall_loc(:,index_wavef_send(:)) = gwavef_alltoall_rcv(:,:)

!----------------------
!Case  -kg = [0 0 0]
!----------------------
 if (idatarecv0/=-1) then
   do kbandpp=1,bandpp_sym
     gwavef_alltoall_loc(:,(kbandpp-1)*ndatarecv     + idatarecv0)= &
     gwavef_alltoall_sym(:,(kbandpp-1)*ndatarecv_tot + idatarecv0)
   end do
 end if

!---------------------------------------------------
!Build of hwavef_alltoall
!
!We have got :
!bandpp_sym blocks to dissociate
!or   bandpp_sym blokcs to not dissociate
!--------------------------------------------------
 do kbandpp=1,bandpp_sym

!  position of the 2 blocks
!  ----------------------------------
   ibandpp = (kbandpp-1) * 2
   jbandpp =  ibandpp    + 1

   idebe = (kbandpp-1) * ndatarecv_tot + 1
   ifine = idebe       + ndatarecv     - 1

   idebc = ibandpp * ndatarecv     + 1
   ifinc = idebc   + ndatarecv     - 1

   idebd = jbandpp * ndatarecv     + 1
   ifind = idebd   + ndatarecv     - 1

   ideb_loc = (kbandpp-1) * ndatarecv  + 1
   ifin_loc = ideb_loc    + ndatarecv  - 1


   if (flag_compose) then

!    calcul cwavef(G)
!    ----------------
     gwavef_alltoall(1,idebc:ifinc) =   gwavef_alltoall_sym(1,idebe:ifine)  &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebc:ifinc) =   gwavef_alltoall_sym(2,idebe:ifine)  &
&     - gwavef_alltoall_loc(2,ideb_loc:ifin_loc)

!    calcul dwavef(G)
!    ------------------
     gwavef_alltoall(1,idebd:ifind) =   gwavef_alltoall_sym(2,idebe:ifine) &
&     + gwavef_alltoall_loc(2,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebd:ifind) = - gwavef_alltoall_sym(1,idebe:ifine) &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
   else

!    calcul cwavef(G)
!    ----------------
     gwavef_alltoall(1,idebc:ifinc) =   gwavef_alltoall_sym(1,idebe:ifine)  &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebc:ifinc) =   gwavef_alltoall_sym(2,idebe:ifine)  &
&     - gwavef_alltoall_loc(2,ideb_loc:ifin_loc)
   end if

 end do

!We divise by two
 gwavef_alltoall(:,:)    = gwavef_alltoall(:,:)/2

!-----------------------
!Desallocation
!-----------------------

 ABI_DEALLOCATE(sendcounts_sym_loc)
 ABI_DEALLOCATE(recvcounts_sym_loc)
 ABI_DEALLOCATE(sdispls_sym_loc)
 ABI_DEALLOCATE(rdispls_sym_loc)

 ABI_DEALLOCATE(gwavef_alltoall_loc)
 ABI_DEALLOCATE(gwavef_alltoall_rcv)

end subroutine prep_wavef_sym_undo
!!***

!!****f* ABINIT/prep_index_wavef_bandpp
!! NAME
!! prep_index_wavef_bandpp
!!
!! FUNCTION
!! this routine sorts the waves functions by bandpp and by processors
!! after the alltoall
!!
!! COPYRIGHT
!!
!! INPUTS
!!  nproc_band = number of processors below the band
!!  bandpp     = number of groups of couple of waves functions
!!  nspinor    = number of spin
!!  ndatarecv  = total number of values received by the processor and sended
!!               by the other processors band
!!  recvcounts = number of values sended by each processor band and received
!!               by the processor
!!  rdispls    = positions of the values received by the processor and
!!               sended by each processor band
!!
!! OUTPUT
!!  index_wavef_band = position of the sorted values
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      chebfi,prep_fourwf,prep_getghc,prep_nonlop
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&
                             nspinor,ndatarecv,&
                             recvcounts,rdispls,&
                             index_wavef_band)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_index_wavef_bandpp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,allocatable,intent(out) :: index_wavef_band(:)

!Local variables-------------------------------
!scalars
 integer :: delta,idebc,idebe,ifinc,ifine,iindex,iproc,kbandpp,nb

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_index_wavef_banpp : enter '
!write(std_out,*) 'ndatarecv = ', ndatarecv
!write(std_out,*) 'rdispls(:) = ', rdispls(:)
!write(std_out,*) 'recvcounts(:) = ', recvcounts(:)
!ENDDEBUG


!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(index_wavef_band ,(bandpp*nspinor*ndatarecv))
 index_wavef_band(:)   =0

!---------------------------------------------
!Calcul : loops on bandpp and processors band
!---------------------------------------------
 nb = sum(recvcounts(1:nproc_band))
 do kbandpp=1,bandpp

   do iproc=1,nproc_band

     idebe = (rdispls(iproc) + 1)  + (kbandpp-1) * ndatarecv*nspinor
     ifine = idebe + recvcounts(iproc) -1

     if (iproc==1) then
       idebc =   (kbandpp-1)* recvcounts(iproc)*nspinor + 1
     else
       idebc = (bandpp)  * sum(recvcounts(1:iproc-1))*nspinor &
       + (kbandpp-1)* recvcounts(iproc)*nspinor &
       + 1
     end if
     ifinc = idebc + recvcounts(iproc) -1
     index_wavef_band(idebe:ifine) = (/( iindex,iindex=idebc,ifinc)/)
     delta=ifine-idebe
     if (nspinor==2) then
       index_wavef_band(idebe+nb :idebe+nb +delta)=(/( iindex,iindex=ifinc+1,ifinc+1+delta)/)
     end if
   end do
 end do

end subroutine prep_index_wavef_bandpp
!!***


!!****f* ABINIT/prep_sort_wavef_spin
!!
!! NAME
!! prep_sort_wavef_spin
!!
!! FUNCTION
!! Compute index used to sort a spinorial wave-function by spin
!! Sort to have all nspinor=1 fisrt, then all nspinor=2
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (MD)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nproc_band=size of "band" communicator
!!  nspinor=number of spinorial components of the wavefunction
!!  ndatarecv=total number of values on all processors
!!  recvcounts(nproc_band)= number of received values by the processor
!!  rdispls(nproc_band)= offsets of the received values by the processor
!!
!! OUTPUT
!!  index_wavef(:)=array containing the sorted indexes (pointer, allocated in this routine)
!!
!! PARENTS
!!      prep_getghc,prep_nonlop
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_sort_wavef_spin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,allocatable,intent(out) :: index_wavef(:)

!Local variables-------------------------------
!scalars
 integer :: isft,isft1,iproc,iindex
!arrays
 integer,allocatable :: recvcountsloc(:),rdisplsloc(:)

! *********************************************************************

 ABI_ALLOCATE(index_wavef,(ndatarecv*nspinor))

 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc,(nproc_band))
 recvcountsloc(:)=recvcounts(:)*2*nspinor
 rdisplsloc(:)=rdispls(:)*2*nspinor

!---------------------------------------------
!Loops on bandpp and processors band
!---------------------------------------------
 isft=0
 do iproc=1,nproc_band

!  ===== Spin up
   if (iproc==1) then
     isft= 0
   else
     isft= sum(recvcounts(1: (iproc-1)))
   end if
   isft1 = 0.5*rdisplsloc(iproc)

   index_wavef(1+isft:isft+recvcounts(iproc))= &
&   (/(iindex,iindex=isft1+1,isft1+recvcounts(iproc))/)

!  =====Spin down
   if (iproc==1)then
     isft=sum(recvcounts(1:nproc_band))
   else
     isft=sum(recvcounts(1:nproc_band)) &
&     +sum(recvcounts(1:iproc-1))
   end if
   isft1 = 0.5 * rdisplsloc(iproc) + recvcounts(iproc)

   index_wavef(1+isft:isft+recvcounts(iproc))= &
&   (/(iindex,iindex=isft1+1,isft1+ recvcounts(iproc))/)

 end do

 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)

end subroutine prep_sort_wavef_spin
!!***
