!!****m* ABINIT/m_initcuda
!! NAME
!! m_initcuda
!!
!! FUNCTION
!!  Module containing all variables concerning GPU device
!!  and the functions needed to extract them
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2025 ABINIT group (MMancini, MT, FDahm)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  Is an experimental development
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#if defined HAVE_GPU_CUDA
#include "cuda_common.h"
#endif

#include "abi_common.h"

module m_initcuda

 use defs_basis
 use m_abicore
 use m_xomp
 use m_xmpi, only: xmpi_world,xmpi_comm_rank,xmpi_comm_size,xmpi_abort,xmpi_sum

#ifdef HAVE_KOKKOS
 use m_kokkos_utils
#endif

#ifdef HAVE_YAKL
 use gator_mod
#endif

 implicit none

#if defined HAVE_GPU_CUDA
 integer,parameter,public :: cudap=kind(CUDA_KIND)
#endif

!Structures
!!***

!!****t* m_initcuda/devGPU_type
!! NAME
!! devGPU_type
!!
!! FUNCTION
!! This structured datatype used to contains GPU properties
!!
!!
!! SOURCE
 type,public :: devGPU_type
  integer :: ndevice  !--number of available devices
  real(dp),allocatable :: maxmemdev(:)  !--max global memory on any device
 end type devGPU_type
!!***

 private

 private ::            &
   prt_device_info !, &    ! To print information about GPU

 public ::             &
   InitGPU,            & ! Initialise GPU
   Get_Mem_Dev,        & ! To obtain the max memory available on GPU device
   Get_ndevice,        & ! Number of devices of Capability>1.2
   CleanGPU,           & ! Clean devGPU_type variables
   setdevice_cuda,     & ! Set device, print info, ...
   unsetdevice_cuda      ! Unset device


CONTAINS !===========================================================
!!***


!!****f* m_initcuda/prt_device_info
!! NAME
!! prt_device_info
!!
!! FUNCTION
!! Print information about GPU device
!!
!! SOURCE

 subroutine prt_device_info(device)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: device
!Local variables ------------------------------
!scalars
 integer :: gflops,constmem,sharemem
 integer :: ii,regist,lenname,ncores,nprocs,nprocs_per_gpu
 real(sp) :: globalmem,clockRate
 integer  :: uuid(4)
 character(20)  :: name
 character(20)  :: formatdev
 character(60)  :: gflops_stg
 character(500) :: msg
!arrays
 integer :: vers(0:1)
! *********************************************************************
#if defined HAVE_GPU
 write(msg,'(a,80a)')' ',('_',ii=1,80)
 call wrtout(std_out,msg,'PERS')
 write(msg,'(a25,a25,a31,a)')  '________________________',&
&  ' Graphic Card Properties ','_______________________________' ,ch10
 call wrtout(std_out,msg,'PERS')

 call get_dev_info(device,name,lenname,vers,globalmem,clockRate,gflops,constmem,sharemem,regist,nprocs,ncores,uuid)
 call gpu_get_nprocs_per_gpu(nprocs_per_gpu);
 if (gflops<0) then
   gflops_stg="undefined (add new def. in version_2_cores function)"
 else
   write(gflops_stg,'(i7,a)') gflops,' GFP'
 end if

 write(formatdev,'(a12,i4,a)') '(a23,i4,a3,a',lenname,')'
 write (msg,formatdev)&
       & '  Device             ',device,' : ',name(1:lenname)
 call wrtout(std_out,msg,'PERS')
 write (msg,'(a,2(i1,a),a,i9,a,a,a,f7.1,a,a,a,i9,a,i9,4a,2(a,i9,2a),a,i9,a,a,4(z8,a1),a)')&
       & ' Revision number:                   ',vers(0),'.',vers(1),ch10, &
       & ' Total amount of global memory: ',nint(globalmem),' Mbytes',ch10, &
       & ' Clock rate:                    ',clockRate,' GHz',ch10, &
       & ' Number of processors/cores:    ',nprocs,'/',ncores,ch10, &
       & ' Max FP64 GFLOPS:               ',trim(gflops_stg),ch10, &
       & ' Total  constant memory:        ',constmem,' bytes',ch10, &
       & ' Shared memory per block:       ',sharemem,' bytes',ch10, &
       & ' Number of registers per block: ',regist,ch10,&
       & ' UUID:                              ',uuid(1),'-',uuid(2),'-',uuid(3),'-',uuid(4),' ',ch10
 call wrtout(std_out,msg,'PERS')
 if(nprocs_per_gpu > 1) then
   write (msg,'(a,i9,a,a)')&
       & ' Global memory per task:        ',nint(globalmem)/nprocs_per_gpu,' Mbytes',ch10
   call wrtout(std_out,msg,'PERS')
 end if
 if(device == -1)then
   write(msg,'(a)')' no cuda-GPU devices found'
   call wrtout(std_out,msg,'PERS')
 end if
 write(msg,'(a,80a)')' ',('_',ii=1,80)
 call wrtout(std_out,msg,'PERS')
#endif
 end subroutine prt_device_info
!!***


!!****f* m_initcuda/InitGPU
!! NAME
!! InitGPU
!!
!! FUNCTION
!! Print information about GPU device
!!
!! SOURCE

 subroutine InitGPU(gpuinfo,device)

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: device
 type(devGPU_type),intent(inout) :: gpuinfo
!Local variables ------------------------------
!scalars
 real(sp) :: locmax
! *********************************************************************
 gpuinfo%ndevice = 0
#if defined HAVE_GPU
!--Initialization
 if(device>-1)then
   !--Get the number of device for this proc
   gpuinfo%ndevice = 1
   ABI_MALLOC(gpuinfo%maxmemdev,(0:1))
   call get_GPU_max_mem(device,locmax)
   gpuinfo%maxmemdev(0:1) = locmax
   call  prt_device_info(device)
 endif
#endif
 end subroutine InitGPU
!!***


!****f* m_initcuda/Get_ndevice
!! NAME
!! Get_ndevice
!!
!! FUNCTION
!! Give the number of device with capability>=1.2
!!
!! SOURCE

 subroutine Get_ndevice(ndevice)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ndevice
! *********************************************************************
#if defined HAVE_GPU
!--Get the number of device for this proc
 call c_get_ndevice(ndevice)
#endif
 end subroutine Get_ndevice
!!***



!!****f* m_initcuda/Get_Mem_Dev
!! NAME
!! Get_Mem_Dev
!!
!! FUNCTION
!! Get the max memory availeble on device
!!
!! INPUTS
!! device  device number
!!
!! OUTPUT
!! max_mem_dev
!!
!! SOURCE

subroutine Get_Mem_Dev(device,max_mem_dev)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: device
 real(sp),intent(out) :: max_mem_dev
!Local variables ------------------------------
! *********************************************************************
#if defined HAVE_GPU
 call get_GPU_max_mem(device,max_mem_dev)
#endif
end subroutine Get_Mem_Dev
!!***


!!****f* m_initcuda/CleanGPU
!! NAME
!! CleanGPU
!!
!! FUNCTION
!! Print information about GPU device
!!
!! SOURCE

 subroutine CleanGPU(gpuinfo)

!Arguments ------------------------------------
!scalars
 type(devGPU_type),intent(inout) :: gpuinfo
! *********************************************************************
#if defined HAVE_GPU
 if (allocated(gpuinfo%maxmemdev))  then
   ABI_FREE(gpuinfo%maxmemdev)
 end if
#endif

 end subroutine CleanGPU
!!***


!!****f* m_initcuda/setdevice_cuda
!! NAME
!! setdevice_cuda
!!
!! FUNCTION
!! Detect and activate a GPU device from current CPU core
!!
!! INPUTS
!!  gpu_devices(12)= list of GPU devices to choose on one node (in case of multiple devices);
!!                   if set to 20*-1, will choose the devices by order of performances.
!!
!! SIDE EFFECTS
!!  gpu_option= which GPU implementation is used (None, CUDA, OpenMP, Kokkos)
!!
!! SOURCE

 subroutine setdevice_cuda(gpu_devices_node,gpu_option)

#ifdef FC_NAG
 use f90_unix_proc
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: gpu_option
!arrays
 integer, intent(in) :: gpu_devices_node(12)
!Local variables ------------------------------
!scalars
 integer :: device,ii,jj,me,nb_devices,nproc,ierr,nprocs_per_gpu
 logical :: testopen
 character(len=500) :: msg
 type(devGPU_type) :: gpuinfo
!arrays
 integer,allocatable :: fastest_devices(:),uuid_array(:,:)
! *********************************************************************

 if (gpu_option==ABI_GPU_DISABLED) return

 nproc=xmpi_comm_size(xmpi_world)
 me=xmpi_comm_rank(xmpi_world)

#if defined HAVE_GPU
 device=-1
 call c_get_ndevice(nb_devices)
 write(msg,'(a,i2,a)') ch10,nb_devices,' GPU device(s) have been detected on the current node'
 call wrtout(std_out,msg,'PERS')

 !nb_devices=min(nb_devices,20)
 if(nb_devices>0) then
   if(nb_devices==1) then
     device=0
   else if(all(gpu_devices_node(1:nb_devices)==-1)) then
     nprocs_per_gpu=max(1,nproc/nb_devices)
     device=me/nprocs_per_gpu
   else
     jj=nb_devices
     do ii=jj,2,-1
       if(gpu_devices_node(ii)==-1) nb_devices=ii-1
     end do
     device=gpu_devices_node(1+mod(me,nb_devices))
   end if

   ! Initialize Kokkos and YAKL if requested
   if(gpu_option==ABI_GPU_KOKKOS) then
#ifdef HAVE_KOKKOS
     ! initialize kokkos
     if (xmpi_comm_rank(xmpi_world) == 0) then
       write(std_out,*)'initializinging kokkos in MPI process ', xmpi_comm_rank(xmpi_world)
     end if
     call kokkos_initialize()

     ! only master MPI process print kokkos config
     if (xmpi_comm_rank(xmpi_world) == 0) then
       call abinit_kokkos_print_config()
     endif
#endif

#ifdef HAVE_YAKL
     call gator_init()
#endif
   end if

   call set_dev(device)
   call check_context(nb_devices,msg)
   ABI_MALLOC(uuid_array,(4,nproc))
   uuid_array(:,:) = 0; nprocs_per_gpu = 0
   call get_gpu_uuid(device,uuid_array(:,me+1))
   call xmpi_sum(uuid_array,xmpi_world,ierr)
   do ii=1,nproc
     if(all(uuid_array(:,ii) == uuid_array(:,me+1))) nprocs_per_gpu = nprocs_per_gpu + 1
   end do
   call gpu_set_nprocs_per_gpu(nprocs_per_gpu)
   if(gpu_option==ABI_GPU_OPENMP) then
     call xomp_set_default_device(device)
   end if
   if(nb_devices==1) then !allocation succeed
     write(msg, '(4a,i1,2a)' ) ch10,&
&     ' setdevice_cuda : COMMENT -',ch10,&
&     '  GPU ',device,' has been properly initialized, continuing...',ch10
     call wrtout(std_out,msg,'PERS')
     !Advertise when there are more than one GPU per task
     !When running on NVIDIA GPU, warn if MPS isn't used
     if(nprocs_per_gpu > 1) then
       write(msg, '(4a,i1,2a)' ) ch10,&
&       ' setdevice_cuda : COMMENT -',ch10,&
&       '  There are ',nprocs_per_gpu,' MPI tasks assigned per GPU.',ch10
       call wrtout(std_out,msg,'PERS')
#ifdef HAVE_GPU_CUDA
       call gpu_get_mps_status(device, ierr);
       if(ierr==0) then
         write(msg, '(10a,i3,3a)' ) ch10,&
&         ' setdevice_cuda : WARNING -',ch10,&
&         "  NVIDIA MPS doesn't seem to be enabled while multiple MPI tasks were assigned per GPU.", ch10,&
&         "  This won't break the code but will harm performance.",ch10,&
&         "  In most environments, MPS is disabled by default and must be requested or enabled by user.",ch10,&
&         "  A GPU percent-share should be given to MPS for further tuning, which would be ideally ", 100/nprocs_per_gpu,&
&         "%", ch10, "  in current MPI tasks per GPU configuration."
         call wrtout(std_out,msg,'PERS')
       else if(ierr==-1) then ! unknow status (when CUDA < 12.3)
         write(msg, '(4a)' ) ch10,&
&         ' setdevice_cuda : WARNING -',ch10,&
&         "  Multiple MPI tasks were assigned per GPU. Make sure NVIDIA MPS has been enabled for best performance !"
         call wrtout(std_out,msg,'PERS')
       end if
#endif
     end if
   else !gpu allocation failed we print error message returned and exit
     device=-1
     call wrtout(std_out,msg,'COLL')
     call xmpi_abort()
     inquire(std_out,OPENED=testopen)
     if (testopen) close(std_out)
#if defined FC_NAG
     call exit(-1)
#elif defined HAVE_FC_EXIT
     call exit(1)
#else
      stop 1
#endif
   end if
   call InitGPU(gpuinfo,device)
   call CleanGPU(gpuinfo)
 else
   gpu_option=ABI_GPU_DISABLED
 end if
#endif
 end subroutine setdevice_cuda
!!***


!!****f* m_initcuda/unsetdevice_cuda
!! NAME
!! unsetdevice_cuda
!!
!! FUNCTION
!! Deactivate a GPU device from current CPU core
!!
!! INPUTS
!!  gpu_option= which GPU implementation is used (None, CUDA, OpenMP, Kokkos)
!!
!! SOURCE

 subroutine unsetdevice_cuda(gpu_option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: gpu_option
!Local variables ------------------------------
!scalars
 character(len=500) :: msg
! *********************************************************************

 if (gpu_option==ABI_GPU_DISABLED) return

#if defined HAVE_GPU

 ! Closing YAKL and Kokkos if opened
 if (gpu_option==ABI_GPU_KOKKOS) then
#ifdef HAVE_YAKL
   call gator_finalize()
   write(std_out,*)'yakl gator finalized'
#endif
#ifdef HAVE_KOKKOS
   ! finalize kokkos
   call kokkos_finalize()
   write(std_out,*)'kokkos finalized'
#endif
 !kokkos_finalize already reset GPU context
 !if (gpu_option/=ABI_GPU_KOKKOS) call unset_dev()
 end if

 if (gpu_option==ABI_GPU_LEGACY) then
   call unset_dev()
 end if

#endif
 end subroutine unsetdevice_cuda
!!***

!!****f* m_initcuda/get_fastest_devices
!! NAME
!! get_fastest_devices
!!
!! FUNCTION
!! In case of multiple devices, sort them by performances
!! and output the resulting list of devices.
!!
!! SOURCE

 subroutine get_fastest_devices(devices,nb_devices)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb_devices
 integer,intent(out) :: devices(:)
!Local variables ------------------------------
!scalars
 integer :: ii,nproc
 character(len=500) :: msg
#if defined HAVE_GPU
 integer :: constmem,gflops,jj,lenname,nprocs,ncores,regist,sharemem
 real(sp) :: clockRate,globalmem
 character(len=20) :: name
#endif
!arrays
#if defined HAVE_GPU
 integer :: vers(0:1)
 integer :: uuid(4)
 integer,allocatable :: isort(:)
 real(dp),allocatable :: flops(:),mem(:)
#endif

! *********************************************************************

 nproc=xmpi_comm_size(xmpi_world)
 if (size(devices)/=nproc) stop 'wrong size for devices array!'

!Default
 do ii=0,nproc-1
   devices(ii+1) = MOD(ii,nb_devices)
 end do
 if (nb_devices==1) return

 write(msg,'(a,i2,a)') ch10,nb_devices,' GPU device(s) have been detected on the current node:'
 call wrtout(std_out,msg,'PERS')

#if defined HAVE_GPU
!Check device(s) properties
 ABI_MALLOC(flops,(nb_devices))
 ABI_MALLOC(mem,  (nb_devices))
 do ii=0,nb_devices-1
   call set_dev(ii)
   call get_dev_info(ii,name,lenname,vers,globalmem,clockRate,gflops,constmem,&
&                    sharemem,regist,nprocs,ncores,uuid)
   flops(ii+1)=dble(gflops) ; mem(ii+1)=dble(globalmem)
   call unset_dev()
   write(msg,'(a,i2,3a,i1,a,i1,a,i6,a,f7.1,a,i7,a,i4,a,i4,a)') &
&   '  Device ',ii,': ',trim(name(1:lenname)),', v',vers(0),'.',vers(1),', Mem=',nint(globalmem),&
&   ' Mbytes, Clock=',clockrate,' GHz, ',gflops,' GFLOPS, ',nprocs,' processors, ',ncores,' cores'
   call wrtout(std_out,msg,'PERS')
 end do

!Sort devices (first by flops, then by memory)
 ABI_MALLOC(isort,(nb_devices))
 isort(:)=(/(ii,ii=1,nb_devices)/)
 call my_sort(flops,mem,isort)

!Distribute cards among procs
 do ii=0,nproc-1
   jj=MOD(ii,nb_devices)
   devices(ii+1) = isort(jj+1)-1
 end do

 ABI_FREE(isort)
 ABI_FREE(flops)
 ABI_FREE(mem)
#endif

contains
!!***

!!****f* m_initcuda/my_sort
!! NAME
!! my_sort
!!
!! FUNCTION
!!  Small sorting routine: change iperm array
!!  according to list1 values then list2 values
!!
!! SOURCE

 subroutine my_sort(list1,list2,iperm)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iperm(:)
 real(dp),intent(in) :: list1(:),list2(:)
!Local variables ------------------------------
!scalars
 integer :: ll,mm,nn,pp
 real(dp) :: xx
!arrays
 real(dp),allocatable :: llist(:)

! *********************************************************************

 nn=size(iperm)
 ABI_MALLOC(llist,(nn))
 llist(:)=list1(:)
 do ll=1,nn-1
   do mm=ll+1,nn
     if (llist(mm)>llist(ll)) then
       xx=llist(ll);llist(ll)=llist(mm);llist(mm)=xx
       pp=iperm(ll);iperm(ll)=iperm(mm);iperm(mm)=pp
     end if
   end do
 end do
 do ll=1,nn-1
   do mm=ll+1,nn
     if (abs(llist(mm)-llist(ll))<tol8) then
       if (list2(iperm(mm))>list2(iperm(ll))) then
         xx=llist(ll);llist(ll)=llist(mm);llist(mm)=xx
         pp=iperm(ll);iperm(ll)=iperm(mm);iperm(mm)=pp
       end if
     end if
   end do
 end do
 ABI_FREE(llist)

 end subroutine my_sort
!!***

 end subroutine get_fastest_devices
!!***

end module m_initcuda
!!***
