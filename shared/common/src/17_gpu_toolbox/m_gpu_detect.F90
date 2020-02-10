!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gpu_detect
!! NAME
!! m_gpu_detect
!!
!! FUNCTION
!! Detects the GPU associated to any cpu and associates a GPU, if
!! possible, to any proc
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (MMancini)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_gpu_detect

 use m_abicore

 use defs_basis
 use m_xmpi

#if defined HAVE_GPU_CUDA
 use m_initcuda, only     : Get_ndevice
#endif


 implicit none

 private

 public ::            &
  find_set_gpu,       &  !Calc. the number of point,GPU,for any proc
  get_topo               !Put the topology of machine in an integer
CONTAINS  !===========================================================
!!***


!!****f* m_gpu_detect/find_set_gpu
!! NAME
!! find_set_gpu
!!
!! FUNCTION
!! Calculate the number of point,GPU,for any proc
!!
!! INPUTS
!!  nproc= number of processor
!!  commcart= mpi communicator
!!
!! OUTPUT
!! ngpu=total number of gpu distributed on all node <=nproc
!! gpu_map(0,nproc-1)=contains for any proc the associated device
!! number. -1 if no gpu is associated
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!
!! SOURCE

 subroutine find_set_gpu(nproc,commcart,gpu_map,ngpu)

  implicit none


!Arguments ------------------------------------
  integer,intent(in) :: nproc,commcart
  integer,intent(out) :: ngpu
  integer,intent(out) :: gpu_map(0:nproc-1)
!Local ---------------------------
  integer :: ierr,ndev,avail_gpu
  integer :: me,icpu,cpu_map_me
  character(20) :: name_ch
  character(20) :: nodes(0:nproc-1)
  character(500) :: msg
! *********************************************************************

  ngpu = 0
  gpu_map = -1
  ndev = 0
  me = 0

#if defined HAVE_GPU_CUDA
  me = xmpi_comm_rank(commcart)

  !--Get the number of compatible device on this CPU
  call Get_ndevice(ndev)

  if(nproc == 1) then
    if(ndev /= 0) gpu_map(0) = 0
    ngpu = count(gpu_map>-1)
    return
  end if

  !--Get the name of the node
  call xmpi_name(name_ch,ierr)

  !--Array containing the number of gpu seen by any cpu
  call  xmpi_allgather(ndev,gpu_map,commcart,ierr)
  !   write(std_out,*)' me,nedevice ',gpu_map

  !--Array containing the name of the cpu
  call  xmpi_allgather(name_ch,nodes,commcart,ierr)

  !--Printing Nodes name
  write(msg,'(3a)')&
    & ' -Node names---------------',ch10,&
    & '   me                name  '
  call wrtout(std_out,msg,'COLL')
  do icpu=0,nproc-1
    write(msg,'(i5,a22)') icpu,trim(nodes(icpu))
    call wrtout(std_out,msg,'COLL')
  end do

  !--research of the cpu on the same node of this cpu
  !   write(std_out,*)'ndev ',ndev
  icpu = 0
  avail_gpu = ndev
  cpu_map_me = -1
  do while(avail_gpu /= 0 .and. icpu <= me )
    if( trim(nodes(icpu)) == trim(name_ch)) then
      !--yes on the same node
      if(me == icpu) cpu_map_me = ndev-avail_gpu
      avail_gpu = avail_gpu -1
    endif
    icpu = icpu +1
  end do

  !--All cpu know the cpu with associated gpu (and which gpu on the node)
  !--Now gpu_map contains the number of the device which is associated
  !  with any cpu (-1 if not)
  call  xmpi_allgather(cpu_map_me,gpu_map,commcart,ierr)

  !--Count the total number of gpu
  ngpu = count(gpu_map>-1)
  !write(std_out,*)'total gpu',ngpu

#endif

end subroutine find_set_gpu
!!***

!!****f* m_gpu_detect/get_topo
!! NAME
!! get_topo
!!
!! FUNCTION
!! Put the topology of machine in an integer
!!
!! INPUTS
!!  nproc= number of processor
!!  ngpu = mpi communicator
!!
!! OUTPUT
!!  topo= 0: 1 cpu;
!!        1: n cpu;
!!        2: 1 cpu 1 gpu;
!!        3: n cpu n gpu
!!        4: n cpu > m gpu;
!!        5: n cpu < m gpu
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!
!! SOURCE

 subroutine get_topo(nproc,ngpu,topo)

  implicit none

!Arguments ------------------------------------
  integer,intent(in)  :: nproc,ngpu
  integer,intent(out) :: topo
!Local ---------------------------
  integer :: ierr,ndev,avail_gpu
! *********************************************************************
  topo = 0
  if(nproc>1) topo = 1 !ncpu>1
  if(ngpu==0) return   !no gpu

  if(ngpu==nproc)then
    topo = 2             !1cpu,1gpu
    if (nproc>1)topo = 3 !ncpu,ngpu
    return
  else
    topo = 4               !ncpu>ngpu
    if(nproc<ngpu)topo = 5 !ncpu<ngpu
  endif

end subroutine get_topo
!!***


end module m_gpu_detect
!!***
