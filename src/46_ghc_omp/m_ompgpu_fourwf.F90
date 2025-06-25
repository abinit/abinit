!!****m* ABINIT/m_ompgpu_fourwf
!! NAME
!!  m_ompgpu_fourwf
!!
!! FUNCTION
!!  Fake module to dupe the build system and allow it to include cuda files
!!   in the chain of dependencies.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2025 ABINIT group (MT)
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

module m_ompgpu_fourwf

 use defs_basis
 use m_abicore
 use m_abi_linalg
 use m_errors
 use m_xomp

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#endif

#ifdef HAVE_GPU
 use m_gpu_toolbox
#endif

 implicit none

 private

#ifdef HAVE_OPENMP_OFFLOAD

! STATIC vars to avoid too much cuda overhead
integer :: fourdp_initialized=0
integer :: fourwf_initialized=0

#define FOURDP_ID 1
#define FOURWF_ID 0

! FFT plans
integer :: fft_size_fourdp=-1
integer :: ndat_fourdp=-1
integer :: fft_size_fourwf=-1
integer :: ndat_fourwf=-1

! GPU ffts buffers (fourwf)
real(dp),allocatable,target :: work_gpu(:,:,:,:)

#endif

 public :: ompgpu_fourdp
 !public :: alloc_ompgpu_fourdp
 !public :: free_ompgpu_fourdp
 public :: ompgpu_fourwf
 public :: alloc_ompgpu_fourwf
 public :: free_ompgpu_fourwf

 ! This routine is here to assess memory requirements
 public :: ompgpu_fourwf_work_mem
contains

function ompgpu_fourwf_work_mem(ngfft, ndat) result(req_mem)

 integer, intent(in) :: ngfft(:), ndat
 integer(kind=c_size_t) :: req_mem

 req_mem = int(2, c_size_t) * dp * ngfft(1) * ngfft(2) * ngfft(3) * ndat

end function ompgpu_fourwf_work_mem

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! - AMD GPUs    : FC_LLVM + HIP
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

subroutine ompgpu_fourdp(cplex,ngfft,ldx,ldy,ldz,ndat,isign,fofg,fofr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ngfft(18),ldx,ldy,ldz,ndat,isign
!arrays
 real(dp),intent(inout) :: fofg(2*ldx*ldy*ldz*ndat)
 real(dp),intent(inout) :: fofr(cplex*ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer      :: n1,n2,n3,nfft_tot
 complex(dpc) :: norm
 logical      :: transfer_fofr, transfer_fofg
 character(len=500) :: msg

! *************************************************************************

 n1=ngfft(1);
 n2=ngfft(2);
 n3=ngfft(3);
 nfft_tot=n1*n2*n3;
 norm=dcmplx(one/(n1*n2*n3), 0.0_dp)

 !*********** CHECK some compatibilities **************
 if( (n1/=ldx) .or. (n2/=ldy) .or. (n3/=ldz)) then
   write(msg,"(a,a,i5,i5,i5,a,i5,i5,i5,a)") "FFT SIZE ERROR: \n when gpu mode is on the fft grid must not be augmented",&
    "(n1,n2,n3) = (", n1,n2,n3,") whereas (ldx,ldy,ldz) = (", ldx,ldy,ldz, ")"
   ABI_ERROR(msg)
 end if

 ! ***********  GPU ALLOCATIONS  ***********************

 if (fourdp_initialized == 0) then
   call alloc_ompgpu_fourdp(ngfft,ndat)
 endif !end of initialisation

 ! If fft size has changed, we realloc our buffers
 if((nfft_tot/=fft_size_fourdp) .or. (ndat/=ndat_fourdp)) then
   call free_ompgpu_fourdp
   call alloc_ompgpu_fourdp(ngfft,ndat)
 end if !end if "fft size changed"

 transfer_fofg=   .not. xomp_target_is_present(c_loc(fofg))
 transfer_fofr=   .not. xomp_target_is_present(c_loc(fofr))

 !$OMP TARGET ENTER DATA MAP(alloc:fofg)    IF(transfer_fofg)
 !$OMP TARGET ENTER DATA MAP(alloc:fofr)    IF(transfer_fofr)
 !$OMP TARGET UPDATE TO(fofg) IF(transfer_fofg .and. isign==FFT_INVERSE)
 !$OMP TARGET UPDATE TO(fofr) IF(transfer_fofr .and. isign==FFT_FORWARD)

 select case (cplex)
 case (2)
   ! Complex to Complex.
   select case (isign)
   case (FFT_INVERSE) ! +1
     !$OMP TARGET DATA USE_DEVICE_ADDR(fofg,fofr)
     call gpu_fft_exec_z2z(FOURDP_ID, c_loc(fofg), c_loc(fofr), FFT_INVERSE)
     !$OMP END TARGET DATA
     call gpu_fft_stream_synchronize(FOURDP_ID)
   case (FFT_FORWARD) ! -1
     !$OMP TARGET DATA USE_DEVICE_ADDR(fofg,fofr)
     call gpu_fft_exec_z2z(FOURDP_ID, c_loc(fofr), c_loc(fofg), FFT_FORWARD)
     !$OMP END TARGET DATA
     call gpu_fft_stream_synchronize(FOURDP_ID)
     ! Normalize here
     call abi_gpu_xscal(2,ldx*ldy*ldz*ndat,norm,fofg,1)
   case default
     ABI_BUG("Wrong isign")
   end select
 case (1)
   ! Real case.
   ABI_BUG("Real case for OpenMP GPU fourdp not handled yet !")
   !select case (isign)
   !case (ABI_FFTW_BACKWARD) ! +1
   !  ! +1; G --> R
   !  call fftw3_c2r_op(nx,ny,nz,ldx,ldy,ldz,ndat,fofg,fofr)
   !case (ABI_FFTW_FORWARD)  ! -1
   !  ! -1; R --> G
   !  call fftw3_r2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,fofr,fofg)
   !case default
   !  ABI_BUG("Wrong isign")
   !end select
 end select

 !$OMP TARGET UPDATE FROM(fofg) IF(transfer_fofg .and. isign==FFT_FORWARD)
 !$OMP TARGET UPDATE FROM(fofr) IF(transfer_fofr .and. isign==FFT_INVERSE)
 !$OMP TARGET EXIT DATA MAP(delete:fofg)    IF(transfer_fofg)
 !$OMP TARGET EXIT DATA MAP(delete:fofr)    IF(transfer_fofr)

end subroutine ompgpu_fourdp

subroutine ompgpu_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,me_g0,ndat,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i,&
&  use_ndo,fofginb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,ldx,ldy,ldz,me_g0,ndat,npwin,npwout,option,mgfft
 real(dp),intent(in) :: weight_i(ndat),weight_r(ndat)
 integer,intent(in),optional :: use_ndo
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),target,intent(inout) :: denpot(cplex*ldx,ldy,ldz)
 real(dp),target,intent(in)    :: fofgin(2,npwin*ndat)
 real(dp),target,intent(inout) :: fofr(2,ldx,ldy,ldz*ndat)
 real(dp),target,intent(out)   :: fofgout(2,npwout*ndat)
 real(dp),target,intent(inout),optional :: fofginb(2,npwin*ndat)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
 real(dp), allocatable :: fofrb(:,:,:,:)
 logical :: l_use_ndo

 real(dp) :: xnorm,one

 integer :: n1,n2,n3,nfft_tot,npwmin
 integer :: cfft_size
 integer :: shift_inv1,shift_inv2,shift_inv3
 integer :: i1,i2,i3,ipw,idat;
 integer :: i1inv,i2inv,i3inv
 logical :: transfer_fofgin, transfer_fofginb, transfer_fofgout, transfer_denpot, transfer_fofr
 integer(C_SIZE_T) :: byte_count

#if defined HAVE_GPU_HIP && defined FC_LLVM
 real(dp), ABI_CONTIGUOUS pointer :: fofr_amdref(:,:,:,:)
#endif

 npwmin=1; if(me_g0==1 .and. istwf_k==2) npwmin=2
 l_use_ndo=.false.; if(present(use_ndo)) l_use_ndo=use_ndo==1

 n1=ngfft(1);
 n2=ngfft(2);
 n3=ngfft(3);
 nfft_tot=n1*n2*n3;

 !*********** CHECK some compatibilities **************
 if( (n1/=ldx) .or. (n2/=ldy) .or. (n3/=ldz)) then
   write(msg,"(a,a,i5,i5,i5,a,i5,i5,i5,a)") "FFT SIZE ERROR: \n when gpu mode is on the fft grid must not be augmented",&
    "(n1,n2,n3) = (", n1,n2,n3,") whereas (ldx,ldy,ldz) = (", ldx,ldy,ldz, ")"
   ABI_ERROR(msg)
 end if

 !*************** CUDA INITIALISATION STAGE ****
 if (fourwf_initialized == 0) then
   call alloc_ompgpu_fourwf(ngfft,ndat)
 endif !end of initialisation


 ! ***********  GPU ALLOCATIONS  ***********************

 ! If fft size has changed, we realloc our buffers
 if((nfft_tot/=fft_size_fourwf) .or. (ndat/=ndat_fourwf)) then
   call free_ompgpu_fourwf
   call alloc_ompgpu_fourwf(ngfft,ndat)
 end if !end if "fft size changed"

 transfer_fofgin= .not. xomp_target_is_present(c_loc(fofgin))  .and. (option/=3)
 transfer_fofgout=.not. xomp_target_is_present(c_loc(fofgout)) .and. (option==2 .or. option==3)
 transfer_denpot =.not. xomp_target_is_present(c_loc(denpot))  .and. (option==2 .or. option==1)
 transfer_fofr=   .not. xomp_target_is_present(c_loc(fofr))
 if(l_use_ndo) then
   transfer_fofginb= .not. xomp_target_is_present(c_loc(fofginb))  .and. (option/=3)
 end if

 !$OMP TARGET ENTER DATA MAP(to:fofgin)     IF(transfer_fofgin)
 !$OMP TARGET ENTER DATA MAP(alloc:fofgout) IF(transfer_fofgout)
 !$OMP TARGET ENTER DATA MAP(alloc:denpot)  IF(transfer_denpot)
 !$OMP TARGET ENTER DATA MAP(alloc:fofr)    IF(transfer_fofr)
 if(l_use_ndo) then
   !$OMP TARGET ENTER DATA MAP(to:fofginb)     IF(transfer_fofginb)
   ABI_MALLOC(fofrb, (2,n1,n2,n3*ndat))
   !$OMP TARGET ENTER DATA MAP(alloc:fofrb)
 end if

 !$OMP TARGET ENTER DATA MAP(alloc:work_gpu)

 if(option==3) then
   !$OMP TARGET UPDATE TO(fofr)
 endif

 if(option==1 .or. option==2) then
   ! We launch async transfert of denpot
   !FIXME This async transfer might be better handled through CUDA/HIP after all...
   ! Issues randomly occurs when using Cray compiler, seems fine with NVHPC.
#ifdef FC_CRAY
   !$OMP TARGET UPDATE TO(denpot) IF(transfer_denpot)
#else
   !$OMP TARGET UPDATE TO(denpot) NOWAIT IF(transfer_denpot)
#endif
   if(option == 1) then
#ifdef FC_CRAY
     !$OMP TARGET ENTER DATA MAP(to:weight_r,weight_i)
#else
     !$OMP TARGET ENTER DATA MAP(to:weight_r,weight_i) NOWAIT
#endif
   endif
 endif

 if(option/=3) then

#if defined HAVE_GPU_HIP && defined FC_LLVM
   !FIXME Work-around for AOMP v15.0.3 (AMD Flang fork)
   ! For some reason, fofr won't be processed normally when passed as argument
   ! of FFT routine within TARGET DATA directives
   fofr_amdref => fofr
#endif

   ! GPU_SPHERE_IN

   cfft_size = 2*n1*n2*n3*ndat

   call gpu_set_to_zero(work_gpu,int(2,c_size_t)*n1*n2*n3*ndat)

   ! During GPU calculation we do some pre-calculation on symetries
   if((istwf_k==2) .or. (istwf_k==4) .or. (istwf_k==6) .or. (istwf_k==8)) then
     shift_inv1 = n1;
   else
     shift_inv1 = n1-1;
   end if

   if((istwf_k>=2) .and. (istwf_k<=5))  then
     shift_inv2 = n2;
   else
     shift_inv2 = n2-1;
   end if

   if((istwf_k==2) .or. (istwf_k==3) .or. (istwf_k==6) .or. (istwf_k==7)) then
     shift_inv3 = n3;
   else
     shift_inv3 = n3-1;
   end if


   !!$OMP TARGET TEAMS LOOP &
   !$OMP TARGET TEAMS DISTRIBUTE &
   !$OMP& PRIVATE(idat) MAP(to:work_gpu,kg_kin,fofgin)
   do idat = 1, ndat
     !$OMP PARALLEL DO PRIVATE(ipw,i1,i2,i3)
     do ipw = 1, npwin
       i1=kg_kin(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
       i2=kg_kin(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
       i3=kg_kin(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
       ! We write cfft(i1,i2,i3)
       ! (double2): cfft[i1 + n1*(i2 + n2*(i3 + n3*idat))] = cg[ipw + npw*idat]
       work_gpu(1, i1, i2, i3+n3*(idat-1)) = fofgin(1, (ipw + npwin*(idat-1)))
       work_gpu(2, i1, i2, i3+n3*(idat-1)) = fofgin(2, (ipw + npwin*(idat-1)))
     end do
   end do

   if(istwf_k > 1) then
     !!$OMP TARGET TEAMS LOOP &
     if (istwf_k==2 .and. me_g0==1) then
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
       !$OMP& PRIVATE(idat) MAP(to:work_gpu,fofgin)
       do idat = 1, ndat
         work_gpu(1, 1, 1, 1+n3*(idat-1)) = fofgin(1, (1 + npwin*(idat-1)))
         work_gpu(2, 1, 1, 1+n3*(idat-1)) = zero
       end do
     end if
     !$OMP TARGET TEAMS DISTRIBUTE &
     !$OMP& PRIVATE(idat) MAP(to:work_gpu,kg_kin,fofgin)
     do idat = 1, ndat
       !$OMP PARALLEL DO PRIVATE(ipw,i1,i2,i3,i1inv,i2inv,i3inv)
       do ipw = npwmin, npwin
         i1=kg_kin(1,ipw); if(i1<0)i1=i1+n1;
         i2=kg_kin(2,ipw); if(i2<0)i2=i2+n2;
         i3=kg_kin(3,ipw); if(i3<0)i3=i3+n3;
#if defined HAVE_GPU_CUDA
         i1inv = modulo(shift_inv1 - i1, n1) + 1
         i2inv = modulo(shift_inv2 - i2, n2) + 1
         i3inv = modulo(shift_inv3 - i3, n3) + 1
#elif defined HAVE_GPU_HIP
         i1inv = (shift_inv1-i1) - ( ((shift_inv1-i1)/n1) * n1 ) + 1
         i2inv = (shift_inv2-i2) - ( ((shift_inv2-i2)/n2) * n2 ) + 1
         i3inv = (shift_inv3-i3) - ( ((shift_inv3-i3)/n3) * n3 ) + 1
#endif
         work_gpu(1, i1inv, i2inv, i3inv+n3*(idat-1)) =  fofgin(1, (ipw + npwin*(idat-1)))
         work_gpu(2, i1inv, i2inv, i3inv+n3*(idat-1)) = -fofgin(2, (ipw + npwin*(idat-1)))
       end do
     end do
   end if

   ! call backward fourrier transform on gpu work_gpu => fofr_gpu
#if defined HAVE_GPU_HIP && defined FC_LLVM
   !$OMP TARGET DATA USE_DEVICE_ADDR(work_gpu,fofr_amdref)
   call gpu_fft_exec_z2z(FOURWF_ID, c_loc(work_gpu), c_loc(fofr_amdref), FFT_INVERSE)
   !$OMP END TARGET DATA
#else
   !$OMP TARGET DATA USE_DEVICE_ADDR(work_gpu,fofr)
   call gpu_fft_exec_z2z(FOURWF_ID, c_loc(work_gpu), c_loc(fofr), FFT_INVERSE)
   !$OMP END TARGET DATA
#endif
   call gpu_fft_stream_synchronize(FOURWF_ID)

   ! In non-diagonal specific use-case, perform same operation with fofginb:
   ! - box-to-sphere : fofginb => work_gpu
   ! - backward FFT  : work_gpu => fofrb
   if(l_use_ndo) then

     call gpu_set_to_zero(work_gpu,int(2,c_size_t)*n1*n2*n3*ndat)

     !$OMP TARGET TEAMS DISTRIBUTE &
     !$OMP& PRIVATE(idat) MAP(to:work_gpu,kg_kin,fofginb)
     do idat = 1, ndat
       !$OMP PARALLEL DO PRIVATE(ipw,i1,i2,i3)
       do ipw = 1, npwin
         i1=kg_kin(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
         i2=kg_kin(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
         i3=kg_kin(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
         ! We write cfft(i1,i2,i3)
         ! (double2): cfft[i1 + n1*(i2 + n2*(i3 + n3*idat))] = cg[ipw + npw*idat]
         work_gpu(1, i1, i2, i3+n3*(idat-1)) = fofginb(1, (ipw + npwin*(idat-1)))
         work_gpu(2, i1, i2, i3+n3*(idat-1)) = fofginb(2, (ipw + npwin*(idat-1)))
       end do
     end do

     !$OMP TARGET DATA USE_DEVICE_ADDR(work_gpu,fofrb)
     call gpu_fft_exec_z2z(FOURWF_ID, c_loc(work_gpu), c_loc(fofrb), FFT_INVERSE)
     !$OMP END TARGET DATA
     call gpu_fft_stream_synchronize(FOURWF_ID)
   end if
 end if

 if(option==0) then
   ! We copy back fofr
   !$OMP TARGET UPDATE from(fofr)    IF(transfer_fofr)
 end if

 if(option==1) then
   ! We finish denpot and weight transferts
   !!$OMP TASKWAIT depend(in:denpot,weight_r,weight_i)
   !$OMP TASKWAIT

   ! Perform density accumulation

   ! Non-diagonal case
   if(l_use_ndo) then

     !!$OMP TARGET TEAMS LOOP &
     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
     !$OMP& PRIVATE(idat,i1,i2,i3) MAP(to:fofr,fofrb,denpot,weight_r,weight_i)
     do i3=1, n3
       do i2=1, n2
         do i1=1, n1
           do idat = 1, ndat
             denpot(i1,i2,i3) = denpot(i1,i2,i3) + weight_r(idat) *&
             &    (fofr(1, i1, i2, i3+(idat-1)*n3)*fofrb(1, i1, i2, i3+(idat-1)*n3)&
             &    +fofr(2, i1, i2, i3+(idat-1)*n3)*fofrb(2, i1, i2, i3+(idat-1)*n3))
             denpot(i1,i2,i3) = denpot(i1,i2,i3) + weight_i(idat) *&
             &    (fofrb(2, i1, i2, i3+(idat-1)*n3)*fofr(1, i1, i2, i3+(idat-1)*n3)&
             &    -fofrb(1, i1, i2, i3+(idat-1)*n3)*fofr(2, i1, i2, i3+(idat-1)*n3))
           end do
         end do
       end do
     end do

   ! General case
   else

     !!$OMP TARGET TEAMS LOOP &
     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
     !$OMP& PRIVATE(idat,i1,i2,i3) MAP(to:fofr,denpot,weight_r,weight_i)
     do i3=1, n3
       do i2=1, n2
         do i1=1, n1
           do idat = 1, ndat
             denpot(i1,i2,i3) = denpot(i1,i2,i3) + fofr(1, i1, i2, i3+(idat-1)*n3)*fofr(1, i1, i2, i3+(idat-1)*n3)*weight_r(idat)
             denpot(i1,i2,i3) = denpot(i1,i2,i3) + fofr(2, i1, i2, i3+(idat-1)*n3)*fofr(2, i1, i2, i3+(idat-1)*n3)*weight_i(idat)
           end do
         end do
       end do
     end do
   end if

   !$OMP TARGET UPDATE from(denpot) IF(transfer_denpot)

 end if

 if(option==2) then
   ! We finished denpot transfert
   !!$OMP TASKWAIT depend(in:denpot)
   !$OMP TASKWAIT

   ! call gpu routine to  Apply local potential
   if(cplex==1) then
     !!$OMP TARGET TEAMS LOOP &
     !$OMP TARGET TEAMS DISTRIBUTE &
     !$OMP& PRIVATE(idat) MAP(to:denpot,fofr)
     do idat = 1, ndat
       !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i1,i2,i3)
       do i3=1, n3
         do i2=1, n2
           do i1=1, n1
             fofr(1, i1, i2, i3+(idat-1)*n3) = fofr(1, i1, i2, i3+(idat-1)*n3) * denpot(i1,i2,i3)
             fofr(2, i1, i2, i3+(idat-1)*n3) = fofr(2, i1, i2, i3+(idat-1)*n3) * denpot(i1,i2,i3)
           end do
         end do
       end do
     end do
   else ! cplex==2
     call gpu_copy(work_gpu, fofr, int(2,c_size_t)*n1*n2*n3*ndat)

     !!$OMP TARGET TEAMS LOOP &
     !$OMP TARGET TEAMS DISTRIBUTE &
     !$OMP& PRIVATE(idat) MAP(to:denpot,fofr,work_gpu)
     do idat = 1, ndat
       !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i1,i2,i3)
       do i3=1, n3
         do i2=1, n2
           do i1=1, n1
             fofr(1, i1, i2, i3+(idat-1)*n3) = &
             &    work_gpu(1, i1, i2, i3+(idat-1)*n3) * denpot(2*i1-1,i2,i3)&
             &    - work_gpu(2, i1, i2, i3+(idat-1)*n3) * denpot(2*i1,i2,i3)
             fofr(2, i1, i2, i3+(idat-1)*n3) = &
             &    work_gpu(2, i1, i2, i3+(idat-1)*n3) * denpot(2*i1-1,i2,i3)&
             &    + work_gpu(1, i1, i2, i3+(idat-1)*n3) * denpot(2*i1,i2,i3)
           end do
         end do
       end do
     end do
   end if
 end if

 if(option==2 .or. option==3) then

   ! call forward fourier transform on gpu: fofr_gpu ==> work_gpu
#if defined HAVE_GPU_HIP && defined FC_LLVM
   !$OMP TARGET DATA USE_DEVICE_ADDR(work_gpu,fofr_amdref)
   call gpu_fft_exec_z2z(FOURWF_ID, c_loc(fofr_amdref), c_loc(work_gpu), FFT_FORWARD)
   !$OMP END TARGET DATA
#else
   !$OMP TARGET DATA USE_DEVICE_ADDR(work_gpu,fofr)
   call gpu_fft_exec_z2z(FOURWF_ID, c_loc(fofr), c_loc(work_gpu), FFT_FORWARD)
   !$OMP END TARGET DATA
#endif
   call gpu_fft_stream_synchronize(FOURWF_ID)

   one=1
   xnorm=one/dble(n1*n2*n3)
   if (istwf_k==2 .and. me_g0==1) then
     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
     !$OMP& PRIVATE(idat) MAP(to:work_gpu,fofgout)
     do idat = 1, ndat
       fofgout (1, 1 + npwout*(idat-1)) = work_gpu(1, 1, 1, 1+n3*(idat-1)) * xnorm
       fofgout (2, 1 + npwout*(idat-1)) = zero
     end do
   end if
   !!$OMP TARGET TEAMS LOOP &
   !$OMP TARGET TEAMS DISTRIBUTE &
   !$OMP& PRIVATE(idat) MAP(to:work_gpu,kg_kout,fofgout)
   do idat = 1, ndat
     !$OMP PARALLEL DO PRIVATE(ipw,i1,i2,i3)
     do ipw = npwmin, npwout
       i1=kg_kout(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
       i2=kg_kout(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
       i3=kg_kout(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

       ! We write cg(ig)
       fofgout (1, ipw + npwout*(idat-1)) = work_gpu(1, i1, i2, i3+n3*(idat-1)) * xnorm
       fofgout (2, ipw + npwout*(idat-1)) = work_gpu(2, i1, i2, i3+n3*(idat-1)) * xnorm
     end do
   end do

   !$OMP TARGET UPDATE FROM(fofgout)   IF(transfer_fofgout)
 end if

 if(option==1 .or. option==2) then
   !$OMP TARGET EXIT DATA MAP(delete:denpot)  IF(transfer_denpot)
   if(option == 1) then
     !$OMP TARGET EXIT DATA MAP(delete:weight_r,weight_i)
   endif
 endif

 !$OMP TARGET EXIT DATA MAP(delete:work_gpu)

 !$OMP TARGET EXIT DATA MAP(delete:fofgin)  IF(transfer_fofgin)
 !$OMP TARGET EXIT DATA MAP(delete:fofgout) IF(transfer_fofgout)
 !$OMP TARGET EXIT DATA MAP(delete:fofr)    IF(transfer_fofr)

 if(l_use_ndo) then
   !$OMP TARGET EXIT DATA MAP(delete:fofginb)     IF(transfer_fofginb)
   !$OMP TARGET EXIT DATA MAP(delete:fofrb)
   ABI_FREE(fofrb)
 end if

end subroutine ompgpu_fourwf
!!***

!!****************************************************************************************************************
!!****************************************************************************************************************
!!****************************************************************************************************************

! Memory allocation routine
subroutine alloc_ompgpu_fourdp(ngfft, ndat)
 integer, intent(in) :: ngfft(18), ndat

 integer :: n1,n2,n3, ldx, ldy, ldz
 integer, target :: t_fft(3)

 fourdp_initialized = 1

 ! Size modification if too little memory allocation
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 ldx = ngfft(4)
 ldy = ngfft(5)
 ldz = ngfft(6)
 fft_size_fourdp=n1*n2*n3
 ndat_fourdp = ndat

 ! Initialisation des plans FFT
 t_fft(1) = n3;
 t_fft(2) = n2;
 t_fft(3) = n1;

 ! Creation du plan
 call gpu_fft_plan_many(FOURDP_ID, 3, c_loc(t_fft), c_null_ptr, 1, ldx*ldy*ldz, c_null_ptr, 1, ldx*ldy*ldz, FFT_Z2Z, ndat);

end subroutine alloc_ompgpu_fourdp

subroutine free_ompgpu_fourdp()

 if(fourdp_initialized==0) return

 ! On detruit l'ancien plan
 call gpu_fft_plan_destroy(FOURDP_ID)

 fourdp_initialized = 0;

end subroutine free_ompgpu_fourdp

subroutine alloc_ompgpu_fourwf(ngfft, ndat)
 integer, intent(in) :: ngfft(18), ndat

 integer :: n1,n2,n3, ldx, ldy, ldz
 integer, target :: t_fft(3),embed(3)

 fourwf_initialized = 1

 ! Size modification if too little memory allocation
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 ldx = ngfft(4)
 ldy = ngfft(5)
 ldz = ngfft(6)
 fft_size_fourwf=n1*n2*n3
 ndat_fourwf = ndat

 ! Initialisation des plans FFT
 t_fft(1) = n3;
 t_fft(2) = n2;
 t_fft(3) = n1;

 embed(1) = ldz
 embed(2) = ldy
 embed(3) = ldx

 ! Creation du plan
 call gpu_fft_plan_many(FOURWF_ID, 3, c_loc(t_fft), c_loc(embed), 1, ldx*ldy*ldz, c_loc(embed), 1, ldx*ldy*ldz, FFT_Z2Z, ndat);

 ABI_MALLOC(work_gpu, (2,n1,n2,n3*ndat))

end subroutine alloc_ompgpu_fourwf

subroutine free_ompgpu_fourwf()

 if(fourwf_initialized==0) return

 ! On detruit l'ancien plan
 call gpu_fft_plan_destroy(FOURWF_ID)

 ABI_FREE(work_gpu)

 fourwf_initialized = 0

end subroutine free_ompgpu_fourwf

#else
! interface for unsupported compilers
subroutine ompgpu_fourdp(cplex,ngfft,ldx,ldy,ldz,ndat,isign,fofg,fofr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ngfft(18),ldx,ldy,ldz,ndat,isign
!arrays
 real(dp),intent(inout) :: fofg(2*ldx*ldy*ldz*ndat)
 real(dp),intent(inout) :: fofr(cplex*ldx*ldy*ldz*ndat)

! *************************************************************************

 ABI_UNUSED((/cplex,ldx,ldy,ldz,ndat,isign/))
 ABI_UNUSED((/ngfft/))
 ABI_UNUSED((/fofg,fofr/))

end subroutine ompgpu_fourdp

subroutine ompgpu_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,ldx,ldy,ldz,ndat,npwin,npwout,option,mgfft
 real(dp),intent(in) :: weight_i(:),weight_r(:)
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),target,intent(inout) :: denpot(cplex*ldx,ldy,ldz)
 real(dp),intent(in) :: fofgin(2,npwin*ndat)
 real(dp),intent(inout),target :: fofr(2,ldx,ldy,ldz*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

 ! Marking unused variables to please code quality checks
 ABI_UNUSED((/cplex,istwf_k,ldx,ldy,ldz,ndat,npwin,npwout,option,mgfft/))
 ABI_UNUSED((/weight_i,weight_r/))
 ABI_UNUSED((/gboundin,gboundout/))
 ABI_UNUSED((/kg_kin,kg_kout,ngfft/))
 ABI_UNUSED((/denpot,fofgin,fofgout,fofr/))
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
end subroutine ompgpu_fourwf

subroutine alloc_ompgpu_fourwf(ngfft, ndat)
!Arguments ------------------------------------
 integer, intent(in) :: ngfft(6), ndat

 ABI_UNUSED((/ngfft,ndat/))
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
end subroutine alloc_ompgpu_fourwf

subroutine free_ompgpu_fourwf()
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
end subroutine free_ompgpu_fourwf

#endif

end module m_ompgpu_fourwf
!!***
