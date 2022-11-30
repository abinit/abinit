!!****m* ABINIT/m_ompgpu_fourwf
!! NAME
!!  m_ompgpu_fourwf
!!
!! FUNCTION
!!  Fake module to dupe the build system and allow it to include cuda files
!!   in the chain of dependencies.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2022 ABINIT group (MT)
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
integer :: fourwf_initialized=0

! FFT plan
integer :: fft_size=-1
integer :: ndat_loc=-1
integer :: npw=-1

! GPU ffts buffers
real(dp),allocatable,target :: work_gpu(:,:,:,:)

#endif

 public :: ompgpu_fourwf
 public :: alloc_ompgpu_fourwf
 public :: free_ompgpu_fourwf
contains

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

subroutine ompgpu_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,ldx,ldy,ldz,ndat,npwin,npwout,option,mgfft
 real(dp),intent(in) :: weight_i(ndat),weight_r(ndat)
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),target,intent(inout) :: denpot(cplex*ldx,ldy,ldz)
 real(dp),intent(in) :: fofgin(2,npwin*ndat)
 real(dp),intent(inout),target :: fofr(2,ldx,ldy,ldz*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: me_g0=1,ndat1=1
 integer :: nx,ny,nz,fftalg,fftalga,fftalgc,fftcache,dat,ptg,ptr,ptgin,ptgout,nthreads
 character(len=500) :: msg

 real(dp) :: xnorm,one

 !Cuda return code
 integer fft_state
 integer rc

 !Local integer
 integer n1,n2,n3,nfft_tot,i,offset

 integer :: cfft_size
 integer :: shift_inv1,shift_inv2,shift_inv3
 integer :: i1,i2,i3,ipw,idat;
 integer :: i1inv,i2inv,i3inv
 logical :: transfer_fofgin, transfer_fofgout, transfer_fofr
 integer(C_SIZE_T) :: byte_count

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

 if(cplex==2)  then
   ABI_ERROR("gpu_fourwf: cplex == 2 is not treated in GPU mode")
 end if

 !*************** CUDA INITIALISATION STAGE ****
 if (fourwf_initialized == 0) then
   call alloc_ompgpu_fourwf(ngfft,ndat)
 endif !end of initialisation


 ! ***********  GPU ALLOCATIONS  ***********************

 ! If fft size has changed, we realloc our buffers
 if((nfft_tot/=fft_size) .or. (ndat/=ndat_loc)) then
   call free_ompgpu_fourwf
   call alloc_ompgpu_fourwf(ngfft,ndat)
 end if !end if "fft size changed"

 transfer_fofgin=.not. xomp_target_is_present(c_loc(fofgin))
 transfer_fofgout=.not. xomp_target_is_present(c_loc(fofgout))
 transfer_fofr=.not. xomp_target_is_present(c_loc(fofr))

 !$OMP TARGET ENTER DATA MAP(to:fofgin)     IF(transfer_fofgin)
 !$OMP TARGET ENTER DATA MAP(alloc:fofgout) IF(transfer_fofgout)
 !$OMP TARGET ENTER DATA MAP(alloc:fofr)    IF(transfer_fofr)

 !$OMP TARGET ENTER DATA MAP(alloc:work_gpu)

 if(option==3) then
   !$OMP TARGET UPDATE TO(fofr)
 endif

 if(option/=3) then
   ! We launch async transfert of denpot
   !$OMP TARGET ENTER DATA MAP(alloc:denpot)
   !$OMP TARGET UPDATE TO(denpot) NOWAIT DEPEND(OUT:denpot)
   if(option == 1) then
     !$OMP TARGET ENTER DATA MAP(to:weight_r,weight_i) NOWAIT DEPEND(OUT:weight_r,weight_i)
   endif

   ! GPU_SPHERE_IN

   cfft_size = 2*n1*n2*n3*ndat

   byte_count=sizeof(work_gpu)
   !$OMP TARGET DATA USE_DEVICE_PTR(work_gpu)
   call gpu_memset_omp(c_loc(work_gpu), zero, byte_count)
   !$OMP END TARGET DATA

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
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
   !$OMP& PRIVATE(idat,ipw,i1,i2,i3) MAP(to:work_gpu,kg_kin,fofgin)
   do idat = 1, ndat
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
     !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
     !$OMP& PRIVATE(idat,ipw,i1,i2,i3,i1inv,i2inv,i3inv) MAP(to:work_gpu,kg_kin,fofgin)
     do idat = 1, ndat
       do ipw = 1, npwin
         i1=kg_kin(1,ipw); if(i1<0)i1=i1+n1;
         i2=kg_kin(2,ipw); if(i2<0)i2=i2+n2;
         i3=kg_kin(3,ipw); if(i3<0)i3=i3+n3;
#ifdef HAVE_GPU_CUDA
         i1inv = modulo(shift_inv1 - i1, n1) + 1
         i2inv = modulo(shift_inv2 - i2, n2) + 1
         i3inv = modulo(shift_inv3 - i3, n3) + 1
#endif
         work_gpu(1, i1inv, i2inv, i3inv+n3*(idat-1)) =  fofgin(1, (ipw + npwin*(idat-1)))
         work_gpu(2, i1inv, i2inv, i3inv+n3*(idat-1)) = -fofgin(2, (ipw + npwin*(idat-1)))
       end do
     end do
   end if

   ! call backward fourrier transform on gpu work_gpu => fofr_gpu
   !$OMP TARGET DATA USE_DEVICE_PTR(work_gpu,fofr)
   call gpu_fft_exec_z2z(c_loc(work_gpu), c_loc(fofr), FFT_INVERSE)
   !$OMP END TARGET DATA
   call gpu_fft_stream_synchronize()
 end if

 if(option==0) then
   ! We copy back fofr
   !$OMP TARGET UPDATE from(fofr)
 end if

 if(option==1) then
   ! We finish denpot and weight transferts
   !!$OMP TASKWAIT depend(in:denpot,weight_r,weight_i)
   !$OMP TASKWAIT

   ! call density accumulation routine on gpu
   !!$OMP TARGET TEAMS LOOP &
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
   !$OMP& PRIVATE(idat,i1,i2,i3) MAP(to:fofr,denpot,weight_r,weight_i) &
   !$OMP& DEPEND(in:denpot,weight_r,weight_i)
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
   !$OMP TARGET UPDATE from(denpot)

 end if

 if(option==2) then
   ! We finished denpot transfert
   !!$OMP TASKWAIT depend(in:denpot)

   ! call gpu routine to  Apply local potential
   !!$OMP TARGET TEAMS LOOP &
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
   !$OMP& PRIVATE(idat,i1,i2,i3) MAP(to:denpot,fofr) DEPEND(in:denpot)
   do idat = 1, ndat
     do i3=1, n3
       do i2=1, n2
         do i1=1, n1
           fofr(1, i1, i2, i3+(idat-1)*n3) = fofr(1, i1, i2, i3+(idat-1)*n3) * denpot(i1,i2,i3)
           fofr(2, i1, i2, i3+(idat-1)*n3) = fofr(2, i1, i2, i3+(idat-1)*n3) * denpot(i1,i2,i3)
         end do
       end do
     end do
   end do
 end if

 if(option==2 .or. option==3) then

   ! call forward fourier transform on gpu: fofr_gpu ==> work_gpu
   !$OMP TARGET DATA USE_DEVICE_PTR(work_gpu,fofr)
   call gpu_fft_exec_z2z(c_loc(fofr), c_loc(work_gpu), FFT_FORWARD)
   !$OMP END TARGET DATA
   call gpu_fft_stream_synchronize()

   one=1
   xnorm=one/dble(n1*n2*n3)
   !!$OMP TARGET TEAMS LOOP &
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
   !$OMP& PRIVATE(idat,ipw,i1,i2,i3) MAP(to:work_gpu,kg_kout,fofgout)
   do idat = 1, ndat
     do ipw = 1, npwout
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

 if(option/=3) then
   ! We launch async transfert of denpot
   !!$OMP TARGET UPDATE FROM(denpot)
   !$OMP TARGET EXIT DATA MAP(release:denpot)
   if(option == 1) then
     !$OMP TARGET EXIT DATA MAP(release:weight_r,weight_i)
   endif
 endif

 !$OMP TARGET EXIT DATA MAP(release:work_gpu)

 !$OMP TARGET EXIT DATA MAP(release:fofgin)  IF(transfer_fofgin)
 !$OMP TARGET EXIT DATA MAP(release:fofgout) IF(transfer_fofgout)
 !$OMP TARGET EXIT DATA MAP(release:fofr)    IF(transfer_fofr)

end subroutine ompgpu_fourwf
!!***

!!****************************************************************************************************************
!!****************************************************************************************************************
!!****************************************************************************************************************

! Memory allocation routine
subroutine alloc_ompgpu_fourwf(ngfft, ndat)
 implicit none
 integer, intent(in) :: ngfft(6), ndat

 integer :: n1,n2,n3, ldx, ldy, ldz
 integer, target :: t_fft(3)
 character(len=500) :: msg
 logical :: pinnedFlag

 !Cuda return code
 integer fft_state

 fourwf_initialized = 1

 ! Size modification if too little memory allocation
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 ldx = ngfft(4)
 ldy = ngfft(5)
 ldz = ngfft(6)
 fft_size=n1*n2*n3
 ndat_loc = ndat

 ! Initialisation des plans FFT
 t_fft(1) = n3;
 t_fft(2) = n2;
 t_fft(3) = n1;

 ! Creation du plan
 call gpu_fft_plan_many(3, c_loc(t_fft), c_null_ptr, 1, ldx*ldy*ldz, c_null_ptr, 1, ldx*ldy*ldz, FFT_Z2Z, ndat);

 ABI_MALLOC(work_gpu, (2,n1,n2,n3*ndat))
 !FIXME Smarter buffer management ?
 !!$OMP TARGET ENTER DATA MAP(alloc:work_gpu)

end subroutine alloc_ompgpu_fourwf

subroutine free_ompgpu_fourwf()
 character(len=500) :: msg

 !Cuda return code
 integer fft_state

 ! On detruit l'ancien plan
 call gpu_fft_plan_destroy()

 !FIXME Smarter buffer management ?
 !!$OMP TARGET EXIT DATA MAP(release:work_gpu)
 ABI_FREE(work_gpu)

 fourwf_initialized = 0;
end subroutine free_ompgpu_fourwf

#else
! interface for unsupported compilers
subroutine ompgpu_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
 implicit none

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

subroutine alloc_ompgpu_fourwf(ngfft, ndat, npwin, npwout)
 implicit none
!Arguments ------------------------------------
 integer, intent(in) :: ngfft(6), ndat, npwin, npwout

 ABI_UNUSED((/ngfft,ndat,npwin,npwout/))
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
end subroutine alloc_ompgpu_fourwf

subroutine free_ompgpu_fourwf()
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
end subroutine free_ompgpu_fourwf

#endif

end module m_ompgpu_fourwf
!!***
