!!****m* ABINIT/m_gpu_detect
!! NAME
!! m_gpu_detect
!!
!! FUNCTION
!! Provide profiling helper routine to annotate execution ranges on both CPU and GPU
!! The code below is just a wrapper around the Nvidia NVTX (version 3) library.
!! It is borrowed from https://developer.nvidia.com/blog/customize-cuda-fortran-profiling-nvtx/
!! This module should (TBC) only be activated when GPU execution is enabled.
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2021 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
module m_nvtx

  use iso_c_binding
  implicit none

  integer,private,parameter :: nbcol=14
  integer,private :: col(nbcol) = [ &
       & Z'0000ff00', & ! GREEN
       & Z'000000ff', & ! BLUE
       & Z'00ffff00', & ! YELLOW
       & Z'00ff00ff', & ! PURPLE
       & Z'0000ffff', & ! CYAN
       & Z'00ff0000', & ! READ
       & Z'00ff8000', & ! ORANGE
       & Z'000080ff', & ! LIGHT BLUE
       & Z'00ff80ff', & ! PINK
       & Z'0080ff80', & ! LIGHT GREEN
       & Z'00b832ff', &
       & Z'00f9fa7d', & ! LIGHT YELLOW
       & Z'00f96c56', &
       & Z'0094b5dc' ]
  character,private,target :: tempName(256)

  type, bind(C):: nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1
     type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes

  interface nvtxRangePush
     ! push range with custom label and standard color
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR) :: name(256)
     end subroutine nvtxRangePushA

     ! push range with custom label and custom color
     subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import:: nvtxEventAttributes
       type(nvtxEventAttributes):: event
     end subroutine nvtxRangePushEx
  end interface nvtxRangePush

  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop

contains

  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
    type(nvtxEventAttributes):: event
    character(kind=c_char,len=256) :: trimmed_name
    integer:: i

    trimmed_name=trim(name)//c_null_char

    ! move scalar trimmed_name into character array tempName
    do i=1,LEN(trim(name)) + 1
       tempName(i) = trimmed_name(i:i)
    enddo


    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,nbcol)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
  end subroutine nvtxStartRange

  subroutine nvtxEndRange
    call nvtxRangePop
  end subroutine nvtxEndRange

end module m_nvtx
