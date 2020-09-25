!!****m* ABINIT/m_hidecudarec
!! NAME
!! m_hidecudarec
!!
!! FUNCTION
!!  Call the C-cu program to make recursion on GPU
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2020 ABINIT group (MMancini)
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
#if defined HAVE_GPU_CUDA
#include "cuda_common.h"
#endif

#include "abi_common.h"


module m_hidecudarec

 use defs_basis
 use defs_rectypes
 use m_abicore
 use m_fft,        only : fourdp

#if defined HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif

 implicit none

 private

#if defined HAVE_GPU_CUDA
 private ::  prt_mem_usage          ! Print memory usage
#endif

#if defined HAVE_GPU_CUDA
 public ::  InitRecGPU_0        ! Initialize recGPU_type
 public ::  InitRecGPU          ! InitRecGPU
 public ::  cudarec             ! Make the recursion on GPU
#endif
 public :: CleanRecGPU            ! deallocate all pointers.


CONTAINS !===========================================================
!!***


!!****f* m_initcuda/prt_mem_usage
!! NAME
!! prt_mem_usage
!!
!! FUNCTION
!! Print information about allocation on GPU device during recursion
!!
!! INPUTS
!! nptrec=number of vectors allocated on device
!! nfft=size of the grid (and so of a vector)
!! PARENTS
!!      m_hidecudarec
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE
#if defined HAVE_GPU_CUDA
subroutine prt_mem_usage(nptrec,nfft)

  implicit none
!Arguments ------------------------------------
  integer,intent(in) :: nptrec,nfft
!Local ---------------------------
  integer :: ii
  integer(kind=i4b) :: largeur,clargeur
  real(dp):: totmem,rpart
  character(500) :: msg
! *********************************************************************

  largeur  = cudap*nfft
  clargeur = cudap*nfft*2
  !for CUDA version <3.0 :
  !   rpart = 3.d0*real(largeur,dp)/1048576.d0*real(nptrec,dp)
  !   totmem = rpart+real(2*clargeur+largeur+(2*cudap+i2b)*nptrec,dp)/1048576.d0
  !for CUDA version 3.0 :
    rpart = 6.d0*real(largeur,dp)/1048576.d0*real(nptrec,dp)
    totmem = rpart+real(clargeur+largeur+(2*cudap+i2b)*nptrec,dp)/1048576.d0

  write(msg,'(a,80a)')' ',('_',ii=1,80)
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a18,a44,a18,a)')'_________________',&
&  '  Allocated Memory on Device for Recursion ','___________________' ,ch10
  call wrtout(std_out,msg,'COLL')

  write (msg,'(2(a32,i10,a),2(a32,i10,a6,a),2(a32,f10.2,a7,a))')&
    & '   Number of Points            ',nfft  ,ch10, &
    & '   Number of Vectors           ',nptrec,ch10, &
    & '   Size Real Vectors           ',largeur ,'bytes',ch10, &
    & '   Size Complex Vectors        ',clargeur,'bytes',ch10, &
    & '   Size Matrix of Vectors      ',real(largeur*nptrec,dp)/1048576.d0,'Mbytes',ch10, &
    & '   Allocated Memory on GPU     ',totmem,'Mbytes',ch10
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,80a)')' ',('_',ii=1,80)
  call wrtout(std_out,msg,'COLL')
end subroutine prt_mem_usage

#endif
!!***

!!****f* m_hidecudarec/InitRecGPU_0
!! NAME
!! InitRecGPU_0
!!
!! FUNCTION
!!  recGPU_type is initialized
!!
!! INPUTS
!! mpi_ab=MPI information
!!
!! OUTPUT
!! recgpu=initialisation of GPU variables for recursion
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE
#if defined HAVE_GPU_CUDA

subroutine InitRecGPU_0(recgpu,mpi_ab)

 implicit none

!Arguments ------------------------------------
 type(MPI_type),intent(in) :: mpi_ab
 type(recGPU_type),intent(inout) :: recgpu
!Local variables-------------------------------
! *************************************************************************
 recgpu%nptrec = 0
 nullify(recgpu%map)
 ABI_ALLOCATE(recgpu%map,(0:mpi_ab%nproc-1))
 recgpu%map = -1       !--Initial guess no gpu

end subroutine InitRecGPU_0
!!***
#endif

!!****f* m_hidecudarec/InitRecGPU
!! NAME
!! InitRecGPU
!!
!! FUNCTION
!!  If there are devices available then the recGPU_type is initialized
!!
!! INPUTS
!!  rset<recusion_type>= contains information of recusion
!!  gpuinfo<devGPU_type>=contains information of GPU
!!  calc_type=if 0 takes the possible max for nptrec (to test the
!!  completly full graphic card). 1 after test to calculate the min
!!  possible value for nptrec

!!
!! OUTPUT
!!  recgpuinfo<recGPU_type>=contains information of recursion with GPU
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE
#if defined HAVE_GPU_CUDA

subroutine InitRecGPU(rset,nfft,gratio,gpudevice,calc_type)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nfft,gpudevice
 integer,intent(in) :: gratio
 integer,intent(in) :: calc_type
 type(recursion_type),intent(inout) :: rset
!Local variables-------------------------------
 integer :: pos_size,resto,nfftc
! real(dp) ::
 type(devGPU_type) :: gpuinfo
 character(len=500) :: msg
! *************************************************************************
 nfftc = nfft/(gratio**3)
 pos_size = 1
 rset%gpudevice = gpudevice


 call InitGPU(gpuinfo,gpudevice)
 !-- look if it is possible to set devices CUDA compatible
 call set_dev(gpudevice)
 if(gpudevice>-1)then
   !--Take the approximate use of memory to compute the number of points on any GPU
   if(rset%tronc)then
     !for CUDA version <3.0 :
     !      pos_size = (.90d0*real(gpuinfo%maxmemdev(0))/real(cudap)&
     !        &           -real(nfft+4*rset%nfftrec))/real((4*rset%nfftrec+15+2))
     ! for CUDA version 3.0 with batched FFT:
    pos_size = (.50d0*real(gpuinfo%maxmemdev(0))/real(cudap)&
      &           -real(nfft+2*rset%nfftrec))/real((6*rset%nfftrec+15+2))


     else
       !for CUDA version <3.0 :
       !       pos_size = (.90d0*real(gpuinfo%maxmemdev(0))/real(cudap)&
       !         &           -real(5*rset%nfftrec))/real((3*rset%nfftrec+15)+2)
       ! for CUDA version 3.0 with batched FFT:
      pos_size = (.5d0*real(gpuinfo%maxmemdev(0))/real(cudap)-real(3&
        &*rset%nfftrec))/real((5*rset%nfftrec)+15+2)

   endif
   !--The nbr of points has to be bigger than 1 and smaller than
   !  rset%par%npt (which is the number of points given to any proc to compute
   !  it is smaller than nfftrec)
   pos_size = min(pos_size,nfftc)

   !--if production and not timing test
   if(calc_type==1) pos_size = min(pos_size,rset%GPU%par%npt)

   if(pos_size<1) then
     write(msg,'(a)')' ERROR NO SUFFICENT MEMORY ON DEVICE'
     call wrtout(std_out,msg,'PERS')
   end if


   !--For GPU calculation it is better to have a number of point
   !  proportional to the half-warp size (16)

   if(pos_size >16 )then
     resto = mod(pos_size,16)
     if(resto /=0) then
       pos_size = pos_size-resto
       if(pos_size<nfftc) pos_size = pos_size+16
     endif
   endif

   rset%GPU%nptrec = pos_size
   if(rset%mpi%me==0) then
     call prt_mem_usage(pos_size,rset%nfftrec)
   end if
 endif
 call CleanGPU(gpuinfo)

end subroutine InitRecGPU
!!***
#endif

!!****f* m_hidecudarec/cudarec
!! NAME
!! cudarec
!!
!! FUNCTION
!! Make recursion on a GPU device
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vtorhorec
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE
#if defined HAVE_GPU_CUDA

subroutine cudarec(rset,exppot,an,bn2,beta,trotter,tolrec,gratio,ngfft,max_rec)

 implicit none

!Arguments ------------------------------------
 integer,intent(in)     :: trotter,gratio
 real(dp),intent(in)    :: beta,tolrec
 integer,intent(inout)  :: max_rec
 type(recursion_type),intent(inout) :: rset
 integer,intent(in)         :: ngfft(1:3)
 real(dp), intent(in)       :: exppot(0:product(ngfft)-1)
 real(cudap), intent(inout) :: an(0:rset%GPU%par%npt-1,0:rset%min_nrec)
 real(cudap), intent(inout) :: bn2(0:rset%GPU%par%npt-1,0:rset%min_nrec)
!Local variables-------------------------------
 ! character(len=500) :: msg
 !integer  ::  maxpt,ipt,ii,jj,kk
 real(dp) :: T_p(0:rset%nfftrec-1)
 ! **integer***********************************************************************


!DEBUG
! write (std_out,*) ' m_hidecudarec/cudarec : enter'
!ENDDEBUG

 call fourdp(1,rset%ZT_p,T_p,1,rset%mpi,rset%nfftrec,1,rset%ngfftrec,0)
 T_p = (one/rset%nfftrec)*T_p

 if(.not.(rset%tronc)) then
   call cuda_rec_cal(trotter,&
     &               gratio,&
     &               rset%GPU%par%npt,&
     &               rset%min_nrec,&
     &               rset%GPU%nptrec,&
     &               max_rec,&
     &               real(beta,cudap),&
     &               real(rset%efermi,cudap),&
     &               real(tolrec,cudap),&
     &               real(rset%inf%ucvol,cudap),&
     &               rset%GPU%par%pt0,&
     &               rset%GPU%par%pt1,&
     &               rset%ngfftrec(1:3),&
     &               real(T_p,cudap),&
     &               real(exppot,cudap),&
     &               an,bn2)

 else


   call cuda_rec_cal_cut(trotter,&
     &                   gratio,&
     &                   rset%GPU%par%npt,&
     &                   rset%min_nrec,&
     &                   rset%GPU%nptrec,&
     &                   max_rec,&
     &                   real(beta,cudap),&
     &                   real(rset%efermi,cudap),&
     &                   real(tolrec,cudap),&
     &                   real(rset%inf%ucvol,cudap),&
     &                   rset%GPU%par%pt0,&
     &                   rset%GPU%par%pt1,&
     &                   ngfft,&
     &                   rset%ngfftrec(1:3),&
     &                   real(T_p,cudap),&
     &                   real(exppot,cudap),&
     &                   an,bn2)

 endif

!DEBUG
!write (std_out,*) ' m_hidecudarec/cudarec : exit'
!ENDDEBUG

end subroutine cudarec
!!***
#endif

!!****f* m_hidecudarec/CleanRecGPU
!! NAME
!! CleanRecGPU
!!
!! FUNCTION
!!  If there are devices availeble than the recGPU_type is initialized
!!
!! INPUTS
!!  load=marks allocation of some arrays
!!  recgpu<type(devGPU_type)>=contains information of GPU
!!
!! OUTPUT
!! nptrec(ndevice)=number of points for recursion on GPU
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

subroutine CleanRecGPU(recgpu,load)

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: load
 type(recGPU_type),intent(inout) :: recgpu
! *************************************************************************

 recgpu%nptrec = 0

 if(associated(recgpu%map))  then
   ABI_DEALLOCATE(recgpu%map)
 end if
 if(load==1)then
   if(allocated(recgpu%par%displs)) then
     ABI_DEALLOCATE(recgpu%par%displs)
   end if
   if(allocated(recgpu%par%vcount)) then
     ABI_DEALLOCATE(recgpu%par%vcount)
   end if
 endif
 call unset_dev()

end subroutine CleanRecGPU
!!***

end module m_hidecudarec
!!***
