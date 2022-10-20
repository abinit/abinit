!!****m* ABINIT/m_xgTools
!! NAME
!!  m_xgTools
!!
!! FUNCTION
!! This is a module to manage and help developer with 2D arrays for low level routines.
!! Particularly, it manages memory for allocations and deallocations (see xg_routines),
!! It handles MPI, complex and real values (*8 kind only) automatically.
!! It is also possible to build sub-block of an array and work on it very easily (see xgBlock_routines)
!! Several routines are also available for performing blas/lapack which again
!! manage the type and MPI (and openmp if needed)
!! Almost all routines are timed by abinit timers
!! This is a starting point and has to be improved/developed
!! An example of how to use those types and routines can be found in
!! 30_diago/m_lobpcg2.F90. This is a full rewrite of LOBPCG algorithm which uses
!! only these types to perfom calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2022 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xg

  use, intrinsic :: iso_c_binding, only: c_loc, c_double, c_double_complex, c_int32_t, c_size_t

  use m_errors
  use m_abicore
  use defs_basis
  use m_time, only : timab
  use m_xmpi, only : xmpi_sum

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
  use m_gpu_toolbox
  use m_xg_kokkos
#endif

#if defined HAVE_YAKL
 use gator_mod
#endif

  implicit none

  private

  integer, parameter, public :: SPACE_R = 1
  integer, parameter, public :: SPACE_C = 2
  integer, parameter, public :: SPACE_CR = 3

  integer, parameter :: tim_potrf = 1670
  integer, parameter :: tim_trsm  = 1671
  integer, parameter :: tim_gemm  = 1672
  integer, parameter :: tim_set   = 1673
  integer, parameter :: tim_get   = 1674
  integer, parameter :: tim_heev  = 1675
  integer, parameter :: tim_heevd = 1676
  integer, parameter :: tim_hpev  = 1677
  integer, parameter :: tim_hpevd = 1678
  integer, parameter :: tim_hegv  = 1679
  integer, parameter :: tim_hegvx = 1680
  integer, parameter :: tim_hegvd = 1681
  integer, parameter :: tim_hpgv  = 1682
  integer, parameter :: tim_hpgvx = 1683
  integer, parameter :: tim_hpgvd = 1684
  integer, parameter :: tim_copy  = 1685
  integer, parameter :: tim_cshift= 1686
  integer, parameter :: tim_pack  = 1687

  integer, save, private :: lrwork = 0
  integer, save, private :: lcwork = 0
  integer, save, private :: liwork = 0
#if defined HAVE_GPU && defined HAVE_YAKL
  integer(kind=c_int32_t),        ABI_CONTIGUOUS pointer, save, private :: iwork(:) => null()
  real(kind=c_double),            ABI_CONTIGUOUS pointer, save, private :: rwork(:) => null()
  complex(kind=c_double_complex), ABI_CONTIGUOUS pointer, save, private :: cwork(:) => null()
#else
  integer, allocatable,  save, private :: iwork(:)
  real(kind=c_double),            allocatable, save, private :: rwork(:)
  complex(kind=c_double_complex), allocatable, save, private :: cwork(:)
#endif

  type, public :: xgBlock_t
    integer, private :: space
    integer, private :: rows
    integer, private :: LDim
    integer, private :: cols
    character, public :: trans
    character, public :: normal
    integer, private :: spacedim_comm
    real(kind=c_double)           , pointer, private :: vecR(:,:) => null()
    complex(kind=c_double_complex), pointer, private :: vecC(:,:) => null()
  end type xgBlock_t

  type, public :: xg_t
    integer, private :: space
    integer, private :: rows
    integer, private :: cols
    character, public :: trans
    character, public :: normal
    integer, private :: spacedim_comm
#if defined HAVE_GPU && defined HAVE_YAKL
    real(kind=c_double)            , ABI_CONTIGUOUS pointer, private :: vecR(:,:) => null()
    complex(kind=c_double_complex) , ABI_CONTIGUOUS pointer, private :: vecC(:,:) => null()
#else
    real(kind=c_double)            , allocatable, private :: vecR(:,:)
    complex(kind=c_double_complex) , allocatable, private :: vecC(:,:)
#endif
    type(xgBlock_t), public :: self
  end type xg_t

  interface xgBlock_gemm
    module procedure xgBlock_gemmR
    module procedure xgBlock_gemmC
  end interface xgBlock_gemm

  interface xgBlock_saxpy
    module procedure xgBlock_saxpyR
    module procedure xgBlock_saxpyC
  end interface xgBlock_saxpy

  interface xgBlock_colwiseMul
    module procedure xgBlock_colwiseMulR
    module procedure xgBlock_colwiseMulC
  end interface xgBlock_colwiseMul

  interface xgBlock_trsm
    module procedure xgBlock_trsmR
    module procedure xgBlock_trsmC
  end interface xgBlock_trsm

  interface xgBlock_scale
    module procedure xgBlock_scaleR
    module procedure xgBlock_scaleC
  end interface xgBlock_scale

  interface checkResize
    module procedure checkResizeI
    module procedure checkResizeR
    module procedure checkResizeC
  end interface checkResize

  public :: space
  public :: cols
  public :: rows
  public :: comm
  public :: xgBlock_setComm
  private :: getClocR
  private :: getClocC
  private :: checkResize

  public :: xg_init
  public :: xg_set
  public :: xg_get
  public :: xg_setBlock
  public :: xg_free

  public :: xg_associated

  public :: xgBlock_setBlock
  public :: xgBlock_set
  public :: xgBlock_map
  public :: xgBlock_reverseMap
  public :: xgBlock_prefetch_async
  public :: xgBlock_get
  public :: xgBlock_copy
  public :: xgBlock_pack
  public :: xgBlock_getSize

  public :: xgBlock_potrf
  public :: xgBlock_trsm

  public :: xgBlock_heev
  public :: xgBlock_heevd

  public :: xgBlock_hpev
  public :: xgBlock_hpevd

  public :: xgBlock_hegv
  public :: xgBlock_hegvx
  public :: xgBlock_hegvd

  public :: xgBlock_hpgv
  public :: xgBlock_hpgvx
  public :: xgBlock_hpgvd

  public :: xgBlock_gemm
  public :: xgBlock_add
  public :: xgBlock_cshift
  public :: xgBlock_colwiseNorm2
  public :: xgBlock_colwiseDotProduct
  public :: xgBlock_colwiseDivision
  public :: xgBlock_colwiseCymax
  public :: xgBlock_saxpy
  public :: xgBlock_colwiseMul
  public :: xgBlock_scale

  public :: xgBlock_zero
  public :: xgBlock_one
  public :: xgBlock_diagonal
  public :: xgBlock_diagonalOnly

  public :: xgBlock_average
  public :: xgBlock_deviation

  public :: xgBlock_reshape
  public :: xgBlock_print
  public :: xg_finalize

contains
  !!***

  !!****f* m_xg/checkResizeI
  !!
  !! NAME
  !! checkResizeI

  subroutine checkResizeI(array,current_dim,asked_dim)

#if defined HAVE_GPU && defined HAVE_YAKL

    integer(kind=c_int32_t), ABI_CONTIGUOUS pointer, intent(inout) :: array(:)
    integer, intent(inout)  :: current_dim
    integer, intent(in   )  :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( associated(array) ) then
        ABI_FREE_MANAGED(array)
      end if
      ABI_MALLOC_MANAGED(array,(/asked_dim/))
    end if

#else

    integer, allocatable, intent(inout) :: array(:)
    integer, intent(inout)  :: current_dim
    integer, intent(in   )  :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
    end if

#endif

  end subroutine checkResizeI
  !!***

  !!****f* m_xg/checkResizeR
  !!
  !! NAME
  !! checkResizeR

  subroutine checkResizeR(array,current_dim,asked_dim)

#if defined HAVE_GPU && defined HAVE_YAKL

    real(kind=c_double), ABI_CONTIGUOUS pointer, intent(inout) :: array(:)
    integer, intent(inout) :: current_dim
    integer, intent(in   ) :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( associated(array) ) then
        ABI_FREE_MANAGED(array)
      end if
      ABI_MALLOC_MANAGED(array,(/asked_dim/))
    end if

#else

    double precision, allocatable, intent(inout) :: array(:)
    integer, intent(inout) :: current_dim
    integer, intent(in   )  :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
    end if

#endif

  end subroutine checkResizeR
  !!***

  !!****f* m_xg/checkResizeC
  !!
  !! NAME
  !! checkResizeC

  subroutine checkResizeC(array,current_dim,asked_dim)

#if defined HAVE_GPU && defined HAVE_YAKL

    complex(kind=c_double_complex), ABI_CONTIGUOUS pointer, intent(inout) :: array(:)
    integer, intent(inout)  :: current_dim
    integer, intent(in   )  :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( associated(array) ) then
        ABI_FREE_MANAGED(array)
      end if
      ABI_MALLOC_MANAGED(array,(/asked_dim/))
    end if

#else

    complex(kind=8), allocatable, intent(inout) :: array(:)
    integer, intent(inout) :: current_dim
    integer, intent(in   )  :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
    end if

#endif

  end subroutine checkResizeC
  !!***

  !!****f* m_xg/getClocR
  !!
  !! NAME
  !! getClocR

  function getClocR(rows,cols,array) result(cptr)
    use, intrinsic :: iso_c_binding
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    double precision, target, intent(inout) :: array(rows,cols)
    type(c_ptr) :: cptr
    cptr = c_loc(array)
  end function getClocR
  !!***

  !!****f* m_xg/getClocC
  !!
  !! NAME
  !! getClocC

  function getClocC(rows,cols,array) result(cptr)
    use, intrinsic :: iso_c_binding
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    complex(kind=8), target, intent(inout) :: array(rows,cols)
    type(c_ptr) :: cptr
    cptr = c_loc(array)
  end function getClocC
  !!***

  !!****f* m_xg/xg_init
  !!
  !! NAME
  !! xg_init

  subroutine xg_init(xg, space, rows, cols, comm)

    type(xg_t), intent(inout) :: xg
    integer   , intent(in   ) :: space
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    integer   , optional, intent(in) :: comm

    integer(kind=c_int32_t), parameter :: izero = 0
    integer(kind=c_size_t)             :: size_bytes

    if ( rows < 1 ) then
      ABI_ERROR("rows < 1 ")
    endif
    if ( cols < 1 ) then
      ABI_ERROR("cols < 1 ")
    end if

    ! MG: Initialize arrays with zero to avoid SIGFPE in xmpi_sum
#if defined HAVE_GPU && defined HAVE_YAKL
    select case (space)
    case (SPACE_R,SPACE_CR)
      if ( associated(xg%vecR) ) then
        ABI_FREE_MANAGED(xg%vecR)
      end if
      ABI_MALLOC_MANAGED_BOUNDS(xg%vecR,(/rows,cols/), (/1,1/))
      !xg%vecR(:,:) = zero
      size_bytes = rows*cols*dp
      call gpu_memset(c_loc(xg%vecR), izero, size_bytes)
      xg%trans = 't'
    case (SPACE_C)
      if ( associated(xg%vecC) ) then
        ABI_FREE_MANAGED(xg%vecC)
      end if
      ABI_MALLOC_MANAGED_BOUNDS(xg%vecC,(/rows,cols/), (/1,1/))
      !xg%vecC(:,:) = zero
      size_bytes = rows*cols*dpc
      call gpu_memset(c_loc(xg%vecC), izero, size_bytes)
      xg%trans = 'c'
    case default
      ABI_ERROR("Invalid space")
    end select
#else
    select case (space)
    case (SPACE_R,SPACE_CR)
      if ( allocated(xg%vecR) ) then
        ABI_FREE(xg%vecR)
      end if
      ABI_MALLOC(xg%vecR,(1:rows,1:cols))
      xg%vecR(:,:) = zero
      xg%trans = 't'
    case (SPACE_C)
      if ( allocated(xg%vecC) ) then
        ABI_FREE(xg%vecC)
      end if
      ABI_MALLOC(xg%vecC,(1:rows,1:cols))
      xg%vecC(:,:) = zero
      xg%trans = 'c'
    case default
      ABI_ERROR("Invalid space")
    end select
#endif

    xg%space = space
    xg%normal = 'n'
    xg%cols = cols
    xg%rows = rows
    xg%spacedim_comm = -1
    if ( present(comm) ) xg%spacedim_comm = comm

    call xg_setBlock(xg,xg%self,1,rows,cols)

  end subroutine xg_init
  !!***

  !!****f* m_xg/xg_set
  !!
  !! NAME
  !! xg_set

  subroutine xg_set(xg,array,shift_col,rows)

    type(xg_t), intent(inout) :: xg
    double precision, intent(in) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_set,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xg%cols ) then
      ABI_WARNING("Ignore some columns, input array to large")
    endif

    select case (xg%space)
    case (SPACE_R)
      do col = 1, min(cols,xg%cols-shift_col)
        xg%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
      end do
    case (SPACE_CR)
      if ( xg%rows /= 2*rows ) then
        ABI_ERROR("Bad number of rows")
      end if

      do col = 1, min(cols,xg%cols-shift_col)
        xg%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
        xg%vecR(rows+1:2*rows,shift_col+col) = array(2,(col-1)*rows+1:col*rows)
      end do
    case (SPACE_C)
      do col = 1, min(cols,xg%cols-shift_col)
        xg%vecC(1:rows,shift_col+col) = dcmplx(array(1,(col-1)*rows+1:col*rows), &
          array(2,(col-1)*rows+1:col*rows))
      end do
    end select

    call timab(tim_set,2,tsec)

  end subroutine xg_set
  !!***

  !!****f* m_xg/xgBlock_set
  !!
  !! NAME
  !! xgBlock_set

  subroutine xgBlock_set(xgBlock,array,shift_col,rows)

    type(xgBlock_t), intent(inout) :: xgBlock
    double precision, intent(in) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_set,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xgBlock%cols ) then
      ABI_WARNING("Block Ignore some columns, input array to large")
    endif

    select case (xgBlock%space)
    case (SPACE_R)
      do col = 1, min(cols,xgBlock%cols-shift_col)
        xgBlock%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
      end do
    case (SPACE_CR)
      if ( xgBlock%rows /= 2*rows ) then
        ABI_ERROR("Bad number of rows")
      end if

      do col = 1, min(cols,xgBlock%cols-shift_col)
        xgBlock%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
        xgBlock%vecR(rows+1:2*rows,shift_col+col) = array(2,(col-1)*rows+1:col*rows)
      end do
    case (SPACE_C)
      do col = 1, min(cols,xgBlock%cols-shift_col)
        xgBlock%vecC(1:rows,shift_col+col) = dcmplx(array(1,(col-1)*rows+1:col*rows), &
          array(2,(col-1)*rows+1:col*rows))
      end do
    end select

    call timab(tim_set,2,tsec)

  end subroutine xgBlock_set
  !!***

  !!****f* m_xg/xgBlock_map
  !!
  !! NAME
  !! xgBlock_map

  subroutine xgBlock_map(xgBlock,array,space,rows,cols,comm)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(inout) :: xgBlock
    double precision, intent(inout) :: array(:,:)
    integer   , intent(in   ) :: space
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    integer   , optional, intent(in) :: comm
    integer :: fullsize
    type(c_ptr) :: cptr

    fullsize = size(array)
    select case (space)
    case ( SPACE_R,SPACE_CR )
      if ( fullsize < cols*rows .or. mod(fullsize,rows) /= 0) then
        ABI_ERROR("Bad size for real array")
      end if
      cptr = getClocR(size(array,dim=1),size(array,dim=2),array)
      call c_f_pointer(cptr,xgBlock%vecR,(/ rows, cols /))
      xgBlock%trans = 't'
    case ( SPACE_C )
      if ( fullsize/2 < cols*rows .or. mod(fullsize/2,rows) /= 0) then
        ABI_ERROR("Bad size for complex array")
      end if
      cptr = getClocR(size(array,dim=1),size(array,dim=2),array)
      call c_f_pointer(cptr,xgBlock%vecC,(/ rows, cols /))
      xgBlock%trans = 'c'
    end select

    xgBlock%space = space
    xgBlock%rows = rows
    xgBlock%LDim = rows
    xgBlock%cols = cols
    xgBlock%normal = 'n'
    if ( present(comm) ) xgBlock%spacedim_comm = comm

  end subroutine xgBlock_map
  !!***

  !!****f* m_xg/xgBlock_reverseMap
  !!
  !! NAME
  !! xgBlock_reverseMap

  subroutine xgBlock_reverseMap(xgBlock,array,rows,cols)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(inout) :: xgBlock
    double precision, pointer, intent(inout) :: array(:,:)
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    type(c_ptr) :: cptr

    select case (xgBlock%space)
    case ( SPACE_R,SPACE_CR )
      if ( xgBlock%cols*xgBlock%Ldim < cols*rows ) then
        write(std_out,*) xgBlock%cols,xgBlock%Ldim,cols,rows
        write(std_out,*) xgBlock%cols*xgBlock%Ldim,cols*rows
        ABI_ERROR("Bad reverseMapping")
      end if
      cptr = getClocR(xgBlock%Ldim,xgBlock%cols,xgBlock%vecR(:,:))
      call c_f_pointer(cptr,array,(/ rows, cols /))
    case ( SPACE_C )
      if ( xgBlock%cols*xgBlock%Ldim < cols*rows ) then
        ABI_ERROR("Bad complex reverseMapping")
      end if
      cptr = getClocC(xgBlock%Ldim,xgBlock%cols,xgBlock%vecC(:,:))
      call c_f_pointer(cptr,array,(/ 2*rows, cols /))
    end select

  end subroutine xgBlock_reverseMap
  !!***

  !!****f* m_xg/xgBlock_prefetch_async
  !!
  !! Help the compiler to upload / dowload data to/from GPU memory
  !!
  !! if deviceId is CPU_DEVICE_ID (= -1, defined in 17_gpu_toolbox/m_gpu_toolbox.F90)
  !! then data is prefetch on host, else it is prefetch to GPU memory
  !!
  !! NAME
  !! xgBlock_prefetch_async

  subroutine xgBlock_prefetch_async(xgBlock, deviceId)
    use iso_c_binding
    type(xgBlock_t) ,             intent(inout) :: xgBlock
    integer(C_INT32_T), optional, intent(in)    :: deviceId

    real(dp), pointer                           :: array(:,:)
    integer                                     :: blockdim
    integer                                     :: spacedim
    integer                                     :: ldim
    integer(C_SIZE_T)                           :: byte_count

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)

    ! get the array pointer underneath the xgBlock
    call xgBlock_getSize(xgBlock,spacedim,blockdim,ldim)
    call xgBlock_reverseMap(xgBlock,array,1,spacedim*blockdim)

    select case (xgBlock%space)
    case ( SPACE_R,SPACE_CR )
      byte_count = ldim*blockdim*dp
    case ( SPACE_C )
      byte_count = ldim*blockdim*dpc
    end select

    ! now we can call the memory prefetch
    if (present(deviceId)) then
      call gpu_data_prefetch_async(c_loc(array), byte_count , deviceId)
    else
      call gpu_data_prefetch_async(c_loc(array), byte_count)
    end if

#endif

  end subroutine xgBlock_prefetch_async
  !!***

  !!****f* m_xg/xg_get
  !!
  !! NAME
  !! xg_get

  subroutine xg_get(xg,array,shift_col,rows)

    type(xg_t), intent(inout) :: xg
    double precision, intent(out) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_get,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xg%cols ) then
      ABI_WARNING("Ignore some columns, input array to large")
    endif

    select case (xg%space)
    case (SPACE_R)
      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xg%vecR(1:rows,shift_col+col)
      end do
    case (SPACE_CR)
      if ( xg%rows /= 2*rows ) then
        ABI_ERROR("Bad number of rows")
      end if

      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xg%vecR(1:rows,shift_col+col)
        array(2,(col-1)*rows+1:col*rows) = xg%vecR(rows+1:2*rows,shift_col+col)
      end do
    case (SPACE_C)
      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = dble(xg%vecC(1:rows,shift_col+col))
        array(2,(col-1)*rows+1:col*rows) = aimag(xg%vecC(1:rows,shift_col+col))
      end do
    end select

    call timab(tim_get,2,tsec)

  end subroutine xg_get
  !!***

  !!****f* m_xg/xgBlock_get
  !!
  !! NAME
  !! xgBlock_get

  subroutine xgBlock_get(xgBlock,array,shift_col,rows)

    type(xgBlock_t), intent(in   ) :: xgBlock
    double precision, intent(out) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_get,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xgBlock%cols ) then
      ABI_ERROR("Ignore some columns, input array to large")
    endif

    select case (xgBlock%space)
    case (SPACE_R)
      do col = 1, min(cols,xgBlock%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xgBlock%vecR(1:rows,shift_col+col)
      end do
    case (SPACE_CR)
      if ( xgBlock%rows /= 2*rows ) then
        ABI_ERROR("Bad number of rows")
      end if

      do col = 1, min(cols,xgBlock%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xgBlock%vecR(1:rows,shift_col+col)
        array(2,(col-1)*rows+1:col*rows) = xgBlock%vecR(rows+1:2*rows,shift_col+col)
      end do
    case (SPACE_C)
      do col = 1, min(cols,xgBlock%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = dble(xgBlock%vecC(1:rows,shift_col+col))
        array(2,(col-1)*rows+1:col*rows) = aimag(xgBlock%vecC(1:rows,shift_col+col))
      end do
    end select

    call timab(tim_get,2,tsec)

  end subroutine xgBlock_get
  !!***

  !!****f* m_xg/xg_setBlock
  !!
  !! NAME
  !! xg_setBlock

  subroutine xg_setBlock(xg, Xgblock, fcol, rows, cols)
    use, intrinsic :: iso_c_binding
    type(xg_t), intent(inout) :: xg
    type(xgBlock_t), intent(inout) :: xgBlock
    integer, intent(in) :: fcol
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    type(c_ptr) :: cptr

    if ( (fcol+cols-1 ) > xg%cols ) then
      ABI_ERROR("Too many columns")
    endif
    if ( rows > xg%rows ) then
      ABI_ERROR("Too many rows")
    end if

    xgBlock%space = xg%space
    xgBlock%rows = rows
    xgBlock%LDim = xg%rows
    xgBlock%cols = cols
    xgBlock%trans = xg%trans
    xgBlock%normal = xg%normal
    xgBlock%spacedim_comm= xg%spacedim_comm

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xg%rows,xg%cols,xg%vecR(:,fcol:fcol+cols-1))
      call c_f_pointer(cptr,xgBlock%vecR,(/ xgBlock%LDim,cols /))
    case(SPACE_C)
      cptr = getClocC(xg%rows,xg%cols,xg%vecC(:,fcol:fcol+cols-1))
      call c_f_pointer(cptr,xgBlock%vecC,(/ xgBlock%LDim,cols /))
    end select

  end subroutine xg_setBlock
  !!***

  !!****f* m_xg/xgBlock_setBlock
  !!
  !! NAME
  !! xgBlock_setBlock

  subroutine xgBlock_setBlock(xgBlockA,xgBlockB, fcol, rows, cols)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB
    integer, intent(in) :: fcol
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    type(c_ptr) :: cptr

    if ( (fcol+cols-1 ) > xgblockA%cols ) then
      ABI_ERROR("Too many columns")
    endif
    if ( rows > xgblockA%rows ) then
      ABI_ERROR("Too many rows")
    end if

    xgBlockB%space = xgBlockA%space
    xgBlockB%rows = rows
    xgBlockB%LDim = xgBlockA%LDim
    xgBlockB%cols = cols
    xgBlockB%trans = xgBlockA%trans
    xgBlockB%normal = xgBlockA%normal
    xgBlockB%spacedim_comm= xgBlockA%spacedim_comm

    select case(xgBlockA%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xgBlockA%LDim,xgBlockA%cols,xgBlockA%vecR(:,fcol:fcol+cols-1))
      call c_f_pointer(cptr,xgBlockB%vecR,(/ xgBlockB%LDim,cols /))
    case(SPACE_C)
      cptr = getClocC(xgBlockA%LDim,xgBlockA%cols,xgBlockA%vecC(:,fcol:fcol+cols-1))
      call c_f_pointer(cptr,xgBlockB%vecC,(/ xgBlockB%LDim,cols /))
    end select

  end subroutine xgBlock_setBlock
  !!***

  !!****f* m_xg/xg_free
  !!
  !! NAME
  !! xg_free

  subroutine xg_free(xg)

    type(xg_t), intent(inout) :: xg

#if defined HAVE_GPU && defined HAVE_YAKL

    if ( associated(xg%vecR) ) then
      ABI_FREE_MANAGED(xg%vecR)
    end if

    if ( associated(xg%vecC) ) then
      ABI_FREE_MANAGED(xg%vecC)
    end if

#else

    if ( allocated(xg%vecR) ) then
      ABI_FREE(xg%vecR)
    end if

    if ( allocated(xg%vecC) ) then
      ABI_FREE(xg%vecC)
    end if

#endif

  end subroutine xg_free
  !!***

  !!****f* m_xg/space
  !!
  !! NAME
  !! space

  function space(xgBlock)

    type(xgBlock_t), intent(in) :: xgBlock
    integer :: space
    space = xgBlock%space
  end function space
  !!***

  !!****f* m_xg/comm
  !!
  !! NAME
  !! comm

  function comm(xgBlock)
    type(xgBlock_t), intent(in) :: xgBlock
    integer :: comm
    comm = xgBlock%spacedim_comm
  end function comm
  !!***

  !!****f* m_xg/setComm
  !!
  !! NAME
  !! setComm

  subroutine xgBlock_setComm(xgBlock,comm)
    type(xgBlock_t), intent(inout) :: xgBlock
    integer :: comm
    xgBlock%spacedim_comm = comm
  end subroutine xgBlock_setComm
  !!***

  !!****f* m_xg/cols
  !!
  !! NAME
  !! cols

  function cols(xgBlock)

    type(xgBlock_t), intent(in) :: xgBlock
    integer :: cols
    cols = xgBlock%cols
  end function cols
  !!***

  !!****f* m_xg/rows
  !!
  !! NAME
  !! rows

  function rows(xgBlock)
    type(xgBlock_t), intent(in) :: xgBlock
    integer :: rows
    rows = xgBlock%rows
    if ( rows /= xgBlock%ldim ) then
      ABI_WARNING("rows/ldim ! Be very careful at what you are doing")
    end if
  end function rows
  !!***

  !!****f* m_xg/xgBlock_copy
  !!
  !! NAME
  !! xgBlock_copy

  subroutine xgBlock_copy(xgBlockA, xgBlockB, inc1, inc2, use_gpu_cuda)

    type(xgBlock_t),   intent(inout) :: xgBlockA
    type(xgBlock_t),   intent(inout) :: xgBlockB
    integer, optional, intent(in   ) :: inc1
    integer, optional, intent(in   ) :: inc2
    integer, optional, intent(in   ) :: use_gpu_cuda

    integer :: incx
    integer :: incy
    integer :: l_use_gpu_cuda = 0

    integer :: size1
    integer :: size2
    integer :: size
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)
    incx = 1; if ( present(inc1) ) incx = inc1
    incy = 1; if ( present(inc2) ) incy = inc2
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if
    !if ( xgBlockA%LDim*xgBlockA%cols/incx /= xgBlockB%LDim*xgBlockB%cols/incy ) then
    !  ABI_ERROR("Number of element different")
    !end if

    size1 = xgBlockA%LDim*xgBlockA%cols/incx ; if ( size1 * incx < xgBlockA%LDim*xgBlockA%cols ) size1 = size1+1
    size2 = xgBlockB%LDim*xgBlockB%cols/incy ; if ( size2 * incy < xgBlockB%LDim*xgBlockB%cols ) size2 = size2+1
    size = min(size1,size2)

    if (l_use_gpu_cuda==1) then

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call gpu_xcopy(1, size, c_loc(xgBlockA%vecR), incx, c_loc(xgBlockB%vecR), incy)
      case(SPACE_C)
        call gpu_xcopy(2, size, c_loc(xgBlockA%vecC), incx, c_loc(xgBlockB%vecC), incy)
      end select

    else

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call dcopy(size,xgBlockA%vecR,incx,xgBlockB%vecR,incy)
      case(SPACE_C)
        call zcopy(size,xgBlockA%vecC,incx,xgBlockB%vecC,incy)
      end select

    end if

    call timab(tim_copy,2,tsec)

  end subroutine xgBlock_copy
  !!***

  !!****f* m_xg/xgBlock_pack
  !!
  !! NAME
  !! xgBlock_pack

  subroutine xgBlock_pack(xgBlockA,xgBlockB,uplo)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB
    character, intent(in) :: uplo
    integer :: j
    integer :: i
    integer :: col
    type(c_ptr) :: cptr
    double precision, pointer :: subR(:)
    complex(kind=8), pointer :: subC(:)
    double precision :: tsec(2)

    call timab(tim_pack,1,tsec)
    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Both blocks must be the same space")
    end if

    if ( xgBlockA%Ldim /= xgBlockA%rows ) then
      ABI_ERROR("Cannot pack when ldim /= rows")
    end if

    if ( xgBlockA%Ldim /= xgBlockA%cols ) then
      ABI_ERROR("Cannot pack when cols /= rows")
    end if

    if ( xgBlockA%rows*(xgBlockA%rows+1)/2 > xgBlockB%Ldim*xgBlockB%cols ) then
      ABI_ERROR("Not enought memory in destination")
    end if

    ! make a fake pointer to pack in a 1-D array instead of 2
    ! This allows to directly use the lapack conventions of transformation
    select case(xgBlockA%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xgBlockB%Ldim,xgBlockB%cols,xgBlockB%vecR(:,:))
      call c_f_pointer(cptr,subR,(/ (xgBlockB%cols*(xgBlockB%cols+1))/2 /))
    case (SPACE_C)
      cptr = getClocC(xgBlockB%Ldim,xgBlockB%cols,xgBlockB%vecC(:,:))
      call c_f_pointer(cptr,subC,(/ (xgBlockB%cols*(xgBlockB%cols+1))/2 /))
    end select

    select case(uplo)
    case ('u','U')
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        do j = 1, xgBlockA%cols
          col = (j*(j-1))/2
          do i = 1, j
            subR(i+col) = xgBlockA%vecR(i,j)
          end do
        end do
      case (SPACE_C)
        do j = 1, xgBlockA%cols
          col = (j*(j-1))/2
          do i = 1, j
            subC(i+col) = xgBlockA%vecC(i,j)
          end do
        end do
      end select

    case ('l','L')
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        do j = 1, xgBlockA%cols
          col = ((2*xgBlockA%cols-j)*(j-1))/2
          do i = j, xgBlockA%cols
            subR(i+col) = xgBlockA%vecR(i,j)
          end do
        end do
      case (SPACE_C)
        do j = 1, xgBlockA%cols
          col = ((2*xgBlockA%cols-j)*(j-1))/2
          do i = j, xgBlockA%cols
            subC(i+col) = xgBlockA%vecC(i,j)
          end do
        end do
      end select
    case default
      ABI_ERROR("Error for packing matrix")
    end select
    call timab(tim_pack,2,tsec)

  end subroutine xgBlock_pack
  !!***

  !!****f* m_xg/xgBlock_gemmR
  !!
  !! NAME
  !! xgBlock_gemmR

  subroutine xgBlock_gemmR(transa, transb, alpha, xgBlockA, xgBlockB, beta, xgBlockW, use_gpu_cuda)

    character,        intent(in   )           :: transa
    character,        intent(in   )           :: transb
    double precision, intent(in   )           :: alpha
    type(xgBlock_t),  intent(in   )           :: xgBlockA
    type(xgBlock_t),  intent(in   )           :: xgBlockB
    double precision, intent(in   )           :: beta
    type(xgBlock_t),  intent(inout)           :: xgBlockW
    integer        ,  intent(in   ), optional :: use_gpu_cuda

    complex(kind=8)  :: calpha
    complex(kind=8)  :: cbeta
    integer          :: K
    double precision :: tsec(2)
    integer          :: l_use_gpu_cuda = 0

    call timab(tim_gemm,1,tsec)

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockB%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if

    if ( transa == 'n' ) then
      K = xgBlockA%cols
    else
      K = xgBlockA%rows
    end if

    calpha = dcmplx(alpha,0.d0)
    cbeta  = dcmplx(beta, 0.d0)

    select case(xgBlockA%space)

    case (SPACE_R,SPACE_CR)
      if (l_use_gpu_cuda==1) then
        call gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          calpha, &
          c_loc(xgBlockA%vecR), xgBlockA%LDim, &
          c_loc(xgBlockB%vecR), xgBlockB%LDim, &
          cbeta, &
          c_loc(xgBlockW%vecR), xgBlockW%LDim)
      else
        call dgemm(transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          alpha, &
          xgBlockA%vecR, xgBlockA%LDim, &
          xgBlockB%vecR, xgBlockB%LDim, &
          beta, &
          xgBlockW%vecR, xgBlockW%LDim)

      end if

      if ( transa == xgBlockA%trans .and. (beta) < 1d-10) then
        if (l_use_gpu_cuda==1) then
          ! CPU waits for GPU to finish before doing MPI communications
#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
          call gpu_device_synchronize()
#endif
        end if
        call xmpi_sum(xgBlockW%vecR,xgBlockW%spacedim_comm,K)
      end if

    case(SPACE_C)

      if (l_use_gpu_cuda==1) then
        call gpu_xgemm(2, transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          calpha, &
          c_loc(xgBlockA%vecC), xgBlockA%LDim, &
          c_loc(xgBlockB%vecC), xgBlockB%LDim, &
          cbeta, &
          c_loc(xgBlockW%vecC), xgBlockW%LDim)
      else
        call zgemm(transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          calpha, &
          xgBlockA%vecC, xgBlockA%LDim, &
          xgBlockB%vecC, xgBlockB%LDim, &
          cbeta, &
          xgBlockW%vecC, xgBlockW%LDim)
      end if

      if ( xgBlockW%spacedim_comm/= -1 .and. transa == xgBlockW%trans .and. abs(beta) < 1d-10 ) then
        if (l_use_gpu_cuda==1) then
          ! CPU waits for GPU to finish before doing MPI communications
#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
          call gpu_device_synchronize()
#endif
        end if
        call xmpi_sum(xgBlockW%vecC,xgBlockW%spacedim_comm,K)
      end if

    end select

    call timab(tim_gemm,2,tsec)

  end subroutine xgBlock_gemmR
  !!***

  !!****f* m_xg/xgBlock_gemmC
  !!
  !! NAME
  !! xgBlock_gemmC

  subroutine xgBlock_gemmC(transa, transb, alpha, xgBlockA, xgBlockB, beta, xgBlockW, use_gpu_cuda)

    character,       intent(in   )           :: transa
    character,       intent(in   )           :: transb
    complex(kind=8), intent(in   )           :: alpha
    type(xgBlock_t), intent(in   )           :: xgBlockA
    type(xgBlock_t), intent(in   )           :: xgBlockB
    complex(kind=8), intent(in   )           :: beta
    type(xgBlock_t), intent(inout)           :: xgBlockW
    integer        , intent(in   ), optional :: use_gpu_cuda

    integer          :: K
    double precision :: tsec(2)
    integer          :: l_use_gpu_cuda = 0

    call timab(tim_gemm,1,tsec)

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockB%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockA%space /= SPACE_C ) then
      ABI_ERROR("Not correct space")
    end if

    if ( transa == 'n' ) then
      K = xgBlockA%cols
    else
      K = xgBlockA%rows
    end if

    if (l_use_gpu_cuda==1) then
      call gpu_xgemm(2, transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
        alpha, &
        c_loc(xgBlockA%vecC), xgBlockA%LDim, &
        c_loc(xgBlockB%vecC), xgBlockB%LDim, &
        beta, &
        c_loc(xgBlockW%vecC), xgBlockW%LDim)
    else
      call zgemm(transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
        alpha, &
        xgBlockA%vecC, xgBlockA%LDim, &
        xgBlockB%vecC, xgBlockB%LDim, &
        beta, &
        xgBlockW%vecC, xgBlockW%LDim)
    end if
    if ( xgBlockW%spacedim_comm/= -1 .and. transa == xgBlockA%trans .and. abs(beta) < 1.d-10 ) then
      if (l_use_gpu_cuda==1) then
        ! CPU waits for GPU to finish before doing MPI communications
#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
        call gpu_device_synchronize()
#endif
      end if
      call xmpi_sum(xgBlockW%vecC,xgBlockW%spacedim_comm,K)
    end if

    call timab(tim_gemm,2,tsec)

  end subroutine xgBlock_gemmC
  !!***

  !!****f* m_xg/xgBlock_potrf
  !!
  !! NAME
  !! xgBlock_potrf

  subroutine xgBlock_potrf(xgBlock,uplo,info)

    type(xgBlock_t), intent(inout) :: xgBlock
    character      , intent(in   ) :: uplo
    integer        , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_potrf,1,tsec)

    if ( xgBlock%rows /= xgBlock%cols ) then
      ABI_ERROR("Matrix should be a square matrixx")
    endif

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      call dpotrf(uplo,xgBlock%rows,xgBlock%vecR,xgBlock%LDim,info)
    case (SPACE_C)
      call zpotrf(uplo,xgBlock%rows,xgBlock%vecC,xgBlock%LDim,info)
    end select

    call timab(tim_potrf,2,tsec)

  end subroutine xgBlock_potrf
  !!***

  !!****f* m_xg/xgBlock_heev
  !!
  !! NAME
  !! xgBlock_heev


  !===================================================
  != Hermitian Full Matrix diago
  !===================================================
  subroutine xgBlock_heev(jobz,uplo,xgBlockA,xgBlockW,info)

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_heev,1,tsec)

    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    select case(xgBlockA%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,8*xgBlockA%rows)

      call dsyev(jobz,uplo,xgBlockA%cols, &
        xgBlockA%vecR,xgBlockA%LDim, &
        xgBlockW%vecR, &
        rwork, lrwork,info)


    case (SPACE_C)
      call checkResize(rwork,lrwork,3*xgBlockA%cols-2)
      call checkResize(cwork,lcwork,lrwork)

      call zheev(jobz,uplo,xgBlockA%cols, &
        xgBlockA%vecC,xgBlockA%LDim, &
        xgBlockW%vecR, &
        cwork, lrwork, rwork, info)

    end select

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    call timab(tim_heev,2,tsec)

  end subroutine xgBlock_heev
  !!***

  !!****f* m_xg/xgBlock_heevd
  !!
  !! NAME
  !! xgBlock_heevd

  subroutine xgBlock_heevd(jobz,uplo,xgBlockA,xgBlockW, info)

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_heevd,1,tsec)

    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlockA%rows+3)

    select case(xgBlockA%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+6*xgBlockA%rows+1)

      call dsyevd(jobz,uplo,xgBlockA%cols, &
        xgBlockA%vecR,xgBlockA%LDim, &
        xgBlockW%vecR, rwork, lrwork, &
        iwork, liwork,info)

    case (SPACE_C)
      call checkResize(cwork,lcwork,xgBlockA%rows*xgBlockA%rows+2*xgBlockA%rows)
      call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+5*xgBlockA%rows+1)

      call zheevd(jobz,uplo,xgBlockA%cols, &
        xgBlockA%vecC,xgBlockA%LDim, &
        xgBlockW%vecR, &
        cwork, lcwork, rwork, lrwork, iwork, liwork, info)

      if ( int(cwork(1)) > lcwork ) then
        !write(std_out,*) "Allocate work from", int(lcwork), "to", int(cwork(1))
        call checkResize(cwork,lcwork,int(cwork(1)))
      end if

    end select

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    if ( iwork(1) > liwork ) then
      !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
      call checkResize(iwork,liwork,int(iwork(1)))
    end if

    call timab(tim_heevd,2,tsec)

  end subroutine xgBlock_heevd
  !!***

  !!****f* m_xg/xgBlock_hpev
  !!
  !! NAME
  !! xgBlock_hpev

  !===================================================
  != Hermitian Packed Matrix diago
  !===================================================
  subroutine xgBlock_hpev(jobz,uplo,xgBlockAP,xgBlockW,xgBlockZ,info)

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockAP
    type(xgBlock_t) , intent(inout) :: xgBlockW
    type(xgBlock_t) , intent(inout) :: xgBlockZ
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpev,1,tsec)

    if ( xgBlockAP%space /= xgBlockZ%space ) then
      ABI_ERROR("Not same space")
    end if

    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    select case(xgBlockAP%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,3*xgBlockZ%cols)

      call dspev(jobz,uplo,xgBlockZ%cols, &
        xgBlockAP%vecR, xgBlockW%vecR, xgBlockZ%vecR, xgBlockZ%Ldim, &
        rwork, info)


    case (SPACE_C)
      call checkResize(cwork,lcwork,2*xgBlockZ%cols-1)
      call checkResize(rwork,lrwork,3*xgBlockZ%cols-2)

      call zhpev(jobz,uplo,xgBlockZ%cols, &
        xgBlockAP%vecC, xgBlockW%vecR, xgBlockZ%vecC, xgBlockZ%Ldim, &
        cwork, rwork, info)

      if ( int(cwork(1)) > lcwork ) then
        !write(std_out,*) "Allocate cwork from", lcwork, "to", int(cwork(1))
        call checkResize(cwork,lcwork,int(cwork(1)))
      end if

    end select

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    call timab(tim_hpev,2,tsec)

  end subroutine xgBlock_hpev
  !!***

  !!****f* m_xg/xgBlock_hpevd
  !!
  !! NAME
  !! xgBlock_hpevd

  subroutine xgBlock_hpevd(jobz,uplo,xgBlockAP,xgBlockW,xgBlockZ,info)

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockAP
    type(xgBlock_t) , intent(inout) :: xgBlockW
    type(xgBlock_t) , intent(inout) :: xgBlockZ
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpevd,1,tsec)

    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    if ( xgBlockAP%space /= xgBlockZ%space ) then
      ABI_ERROR("Block 1 and 3 must have the same space")
    end if

    call checkResize(iwork,liwork,5*xgBlockZ%rows+3)
    select case(xgBlockAP%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,xgBlockZ%rows*xgBlockZ%rows+6*xgBlockZ%rows+1)

      call dspevd(jobz,uplo,xgBlockZ%cols, &
        xgBlockAP%vecR, xgBlockW%vecR, xgBlockZ%vecR, xgBlockZ%Ldim, &
        rwork, lrwork, iwork, liwork,info)


    case (SPACE_C)
      call checkResize(cwork,lcwork,2*xgBlockZ%rows)
      call checkResize(rwork,lrwork,2*xgBlockZ%rows*xgBlockZ%rows+5*xgBlockZ%rows+1)

      call zhpevd(jobz,uplo,xgBlockZ%cols, &
        xgBlockAP%vecC, xgBlockW%vecR, xgBlockZ%vecC, xgBlockZ%Ldim, &
        cwork, lcwork, rwork, lrwork, iwork, liwork, info)

      if ( int(cwork(1)) > lcwork ) then
        !write(std_out,*) "Allocate work from", lcwork, "to", int(cwork(1))
        call checkResize(cwork,lcwork,int(cwork(1)))
      end if

    end select

    if ( int(rwork(1)) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    if ( iwork(1) > liwork ) then
      !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
      call checkResize(iwork,liwork,int(iwork(1)))
    end if

    call timab(tim_hpevd,2,tsec)

  end subroutine xgBlock_hpevd
  !!***

  !!****f* m_xg/xgBlock_hgev
  !!
  !! NAME
  !! xgBlock_hgev

  !===================================================
  != Hermitian Full Generalized Matrix diago
  !===================================================

  subroutine xgBlock_hegv(itype, jobz, uplo, xgBlockA, xgBlockB, xgBlockW, info)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockB
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hegv,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    select case(xgBlockA%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+6*xgBlockA%rows+1)

      call dsygv(itype, jobz, uplo, xgBlockA%rows, xgBlockA%vecR, xgBlockA%ldim, &
        xgBlockB%vecR, xgBlockB%ldim, xgBlockW%vecR, rwork, lrwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,2*xgBlockA%rows-1)
      call checkResize(rwork,lrwork,3*xgBlockA%rows-2)

      call zhegv(itype, jobz, uplo, xgBlockA%rows, xgBlockA%vecC, xgBlockA%ldim,&
        xgBlockB%vecC, xgBlockB%ldim, xgBlockW%vecR, cwork, lcwork, &
        rwork, info)

      if ( int(cwork(1)) > lcwork ) then
        !write(std_out,*) "Allocate work from", lcwork, "to", int(cwork(1))
        call checkResize(cwork,lcwork,int(cwork(1)))
      end if

    end select

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate rwork from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    call timab(tim_hegv,2,tsec)

  end subroutine xgBlock_hegv
  !!***

  !!****f* m_xg/xgBlock_hegvx
  !!
  !! NAME
  !! xgBlock_hegvx

  subroutine xgBlock_hegvx(itype,jobz,range,uplo,xgBlockA,xgBlockB,vl,vu,il,iu,abstol,xgBlockW,xgBlockZ,info)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: range
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockB
    double precision, intent(in   ) :: vl
    double precision, intent(in   ) :: vu
    integer         , intent(in   ) :: il
    integer         , intent(in   ) :: iu
    double precision, intent(in   ) :: abstol
    type(xgBlock_t) , intent(inout) :: xgBlockW
    type(xgBlock_t) , intent(inout) :: xgBlockZ
    integer         , intent(  out) :: info
    integer :: neigen
    integer, allocatable :: ifail(:)
    double precision :: tsec(2)

    call timab(tim_hegvx,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockA%space /= xgBlockZ%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlockA%rows)

    ABI_MALLOC(ifail,(xgBlockA%rows))
    ifail = 0

    select case(xgBlockA%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,8*xgBlockA%rows)

      call dsygvx(itype,jobz,range,uplo,xgBlockA%rows, &
        xgBlockA%vecR,xgBlockA%LDim,xgBlockB%vecR,xgBlockB%LDim, &
        vl,vu,il,iu,abstol,&
        neigen,xgBlockW%vecR, xgBlockZ%vecR, xgBlockZ%LDim, &
        rwork, lrwork,iwork,ifail,info)

    case (SPACE_C)
      call checkResize(rwork,lrwork,7*xgBlockA%rows)
      call checkResize(cwork,lcwork,lrwork)

      call zhegvx(itype,jobz,range,uplo,xgBlockA%rows, &
        xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim, &
        vl,vu,il,iu,abstol,&
        neigen,xgBlockW%vecR, xgBlockZ%vecC, xgBlockZ%LDim, &
        cwork, lcwork, rwork, iwork,ifail,info)

    end select
    ABI_FREE(ifail)

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if


    call timab(tim_hegvx,2,tsec)

  end subroutine xgBlock_hegvx
  !!***

  !!****f* m_xg/xgBlock_hegvd
  !!
  !! NAME
  !! xgBlock_hegvd

  subroutine xgBlock_hegvd(itype, jobz, uplo, xgBlockA, xgBlockB, xgBlockW, info, use_gpu_cuda)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockB
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer         , intent(  out) :: info
    integer         , intent(in   ), optional :: use_gpu_cuda

    double precision :: tsec(2)
    integer          :: l_use_gpu_cuda = 0
    integer          :: bufferSize = -1

    call timab(tim_hegvd,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlockA%space)

      case (SPACE_R,SPACE_CR)

        ! probe needed bufferSize
        call gpu_xsygvd_buffersize(1, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             c_loc(xgBlockA%vecR), xgBlockA%ldim, &
          &             c_loc(xgBlockB%vecR), xgBlockB%ldim, &
          &             c_loc(xgBlockW%vecR), &
          &             bufferSize)


        ! resize rwork if needed
        call checkResize(rwork,lrwork,bufferSize)

        ! and compute
        call gpu_xsygvd(1, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             c_loc(xgBlockA%vecR), xgBlockA%ldim, &
          &             c_loc(xgBlockB%vecR), xgBlockB%ldim, &
          &             c_loc(xgBlockW%vecR), &
          &             c_loc(rwork), bufferSize, info)

      case (SPACE_C)

        ! probe needed bufferSize
        call gpu_xsygvd_buffersize(2, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             c_loc(xgBlockA%vecC), xgBlockA%ldim, &
          &             c_loc(xgBlockB%vecC), xgBlockB%ldim, &
          &             c_loc(xgBlockW%vecR), &
          &             bufferSize)

        ! resize cwork if needed
        call checkResize(cwork,lcwork,bufferSize)

        ! and compute
        call gpu_xsygvd(2, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             c_loc(xgBlockA%vecC), xgBlockA%ldim, &
          &             c_loc(xgBlockB%vecC), xgBlockB%ldim, &
          &             c_loc(xgBlockW%vecR), &
          &             c_loc(cwork), bufferSize, info)

      end select
#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
      call wrtout(std_out,"We shouldn't be here : abinit was not compiled with GPU/CUDA support (Kokkos+YAKL).")
      call abi_abort('COLL')
#endif

    else

      call checkResize(iwork,liwork,5*xgBlockA%rows+3)

      select case(xgBlockA%space)

      case (SPACE_R,SPACE_CR)

        call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+6*xgBlockA%rows+1)

        call dsygvd(itype, jobz, uplo, xgBlockA%rows, xgBlockA%vecR, xgBlockA%ldim, &
          xgBlockB%vecR, xgBlockB%ldim, xgBlockW%vecR, rwork, lrwork, iwork, liwork, info)

      case (SPACE_C)

        call checkResize(cwork,lcwork,xgBlockA%rows*xgBlockA%rows+2*xgBlockA%rows)
        call checkResize(rwork,lrwork,2*(xgBlockA%rows*xgBlockA%rows)+5*xgBlockA%rows+1)

        call zhegvd(itype, jobz, uplo, xgBlockA%rows, xgBlockA%vecC, xgBlockA%ldim,&
          xgBlockB%vecC, xgBlockB%ldim, xgBlockW%vecR, cwork, lcwork, &
          rwork, lrwork, iwork, liwork, info)

        if ( int(cwork(1)) > lcwork ) then
          !write(std_out,*) "Allocate work from", lcwork, "to", int(cwork(1))
          call checkResize(cwork,lcwork,int(cwork(1)))
        end if

      end select

      if ( rwork(1) > lrwork ) then
        !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
        call checkResize(rwork,lrwork,int(rwork(1)))
      end if


      if ( iwork(1) > liwork ) then
        !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
        call checkResize(iwork,liwork,int(iwork(1)))
      end if

    end if

    call timab(tim_hegvd,2,tsec)

  end subroutine xgBlock_hegvd

  !!***

  !!****f* m_xg/xgBlock_hpgv
  !!
  !! NAME
  !! xgBlock_hpgv

  !===================================================
  != Hermitian Full Generalized Matrix diago
  !===================================================

  subroutine xgBlock_hpgv(itype, jobz, uplo, xgBlockAP, xgBlockBP, xgBlockW, xgBlockZ,info)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockAP
    type(xgBlock_t) , intent(inout) :: xgBlockBP
    type(xgBlock_t) , intent(inout) :: xgBlockW
    type(xgBlock_t) , intent(inout) :: xgBlockZ
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpgv,1,tsec)

    if ( xgBlockAP%space /= xgBlockBP%space ) then
      ABI_ERROR("Not same space")
    end if

    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    select case(xgBlockAP%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,3*xgBlockZ%rows)

      call dspgv(itype, jobz, uplo, xgBlockZ%rows, xgBlockAP%vecR, xgBlockBP%vecR, &
        xgBlockW%vecR, xgBlockZ%vecR, xgBlockZ%Ldim, rwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,2*xgBlockZ%rows-1)
      call checkResize(rwork,lrwork,3*xgBlockZ%rows-2)

      call zhpgv(itype, jobz, uplo, xgBlockAP%rows, xgBlockAP%vecC, xgBlockBP%vecC, &
        xgBlockW%vecR, xgBlockZ%vecC, xgBlockZ%Ldim, cwork, rwork, info)

      if ( int(cwork(1)) > lcwork ) then
        !write(std_out,*) "Allocate work from", lcwork, "to", int(cwork(1))
        call checkResize(cwork,lcwork,int(cwork(1)))
      end if

    end select

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate rwork from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    call timab(tim_hpgv,2,tsec)

  end subroutine xgBlock_hpgv
  !!***

  !!****f* m_xg/xgBlock_hpgvx
  !!
  !! NAME
  !! xgBlock_hpgvx

  subroutine xgBlock_hpgvx(itype,jobz,range,uplo,xgBlockAP,xgBlockBP,vl,vu,il,iu,abstol,xgBlockW,xgBlockZ,info)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: range
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockAP
    type(xgBlock_t) , intent(inout) :: xgBlockBP
    double precision, intent(in   ) :: vl
    double precision, intent(in   ) :: vu
    integer         , intent(in   ) :: il
    integer         , intent(in   ) :: iu
    double precision, intent(in   ) :: abstol
    type(xgBlock_t) , intent(inout) :: xgBlockW
    type(xgBlock_t) , intent(inout) :: xgBlockZ
    integer         , intent(  out) :: info
    integer :: neigen
    integer, allocatable :: ifail(:)
    double precision :: tsec(2)

    call timab(tim_hpgvx,1,tsec)

    if ( xgBlockAP%space /= xgBlockBP%space .or. xgBlockAP%space /= xgBlockZ%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlockZ%rows)

    ABI_MALLOC(ifail,(xgBlockZ%rows))
    ifail = 0

    select case(xgBlockAP%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,8*xgBlockZ%rows)

      call dspgvx(itype,jobz,range,uplo,xgBlockZ%rows, &
        xgBlockAP%vecR,xgBlockBP%vecR, &
        vl,vu,il,iu,abstol,&
        neigen,xgBlockW%vecR, xgBlockZ%vecR, xgBlockZ%LDim, &
        rwork, iwork, ifail, info)

    case (SPACE_C)
      call checkResize(rwork,lrwork,7*xgBlockAP%rows)
      call checkResize(cwork,lcwork,2*xgBlockZ%rows)

      call zhpgvx(itype,jobz,range,uplo,xgBlockZ%rows, &
        xgBlockAP%vecC,xgBlockBP%vecC, &
        vl,vu,il,iu,abstol,&
        neigen,xgBlockW%vecR, xgBlockZ%vecC, xgBlockZ%LDim, &
        cwork, rwork, iwork, ifail, info)

    end select
    ABI_FREE(ifail)

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if


    call timab(tim_hpgvx,2,tsec)

  end subroutine xgBlock_hpgvx
  !!***

  !!****f* m_xg/xgBlock_hpgvd
  !!
  !! NAME
  !! xgBlock_hpgvd

  subroutine xgBlock_hpgvd(itype, jobz, uplo, xgBlockAP, xgBlockBP, xgBlockW, xgBlockZ, info)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockAP
    type(xgBlock_t) , intent(inout) :: xgBlockBP
    type(xgBlock_t) , intent(inout) :: xgBlockW
    type(xgBlock_t) , intent(inout) :: xgBlockZ
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpgvd,1,tsec)

    if ( xgBlockAP%space /= xgBlockBP%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlockZ%rows+3)

    select case(xgBlockAP%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlockZ%rows*xgBlockZ%rows+6*xgBlockZ%rows+1)

      call dspgvd(itype, jobz, uplo, xgBlockZ%rows, xgBlockAP%vecR, xgBlockBP%vecR, &
        xgBlockW%vecR, xgBlockZ%vecR, xgBlockZ%Ldim, &
        rwork, lrwork, iwork, liwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,2*xgBlockZ%rows)
      call checkResize(rwork,lrwork,2*(xgBlockZ%rows*xgBlockZ%rows)+5*xgBlockZ%rows+1)

      call zhpgvd(itype, jobz, uplo, xgBlockZ%rows, xgBlockAP%vecC, xgBlockBP%vecC, &
        xgBlockW%vecR, xgBlockZ%vecC, xgBlockZ%Ldim, &
        cwork, lcwork, rwork, lrwork, iwork, liwork, info)

      if ( int(cwork(1)) > lcwork ) then
        !write(std_out,*) "Allocate work from", lcwork, "to", int(cwork(1))
        call checkResize(cwork,lcwork,int(cwork(1)))
      end if

    end select

    if ( rwork(1) > lrwork ) then
      !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
      call checkResize(rwork,lrwork,int(rwork(1)))
    end if

    if ( iwork(1) > liwork ) then
      !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
      call checkResize(iwork,liwork,int(iwork(1)))
    end if

    call timab(tim_hpgvd,2,tsec)

  end subroutine xgBlock_hpgvd
  !!***

  !!****f* m_xg/xgBlock_trsmR
  !!
  !! NAME
  !! xgBlock_trsmR

  subroutine xgBlock_trsmR(side,uplo,transa,diag,alpha, xgBlockA,xgBlockB)

    character       , intent(in   ) :: side
    character       , intent(in   ) :: uplo
    character       , intent(in   ) :: transa
    character       , intent(in   ) :: diag
    double precision, intent(in   ) :: alpha
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockB
    complex(kind=8) :: calpha
    double precision :: tsec(2)

    call timab(tim_trsm,1,tsec)
    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if

    select case(xgBlockA%space)
    case (SPACE_R,SPACE_CR)
      call dtrsm(side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
        alpha,xgBlockA%vecR,xgBlockA%LDim,xgBlockB%vecR,xgBlockB%LDim)
    case (SPACE_C)
      calpha = dcmplx(alpha,0.d0)
      call ztrsm(side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
        calpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)
    end select

    call timab(tim_trsm,2,tsec)

  end subroutine xgBlock_trsmR
  !!***

  !!****f* m_xg/xgBlock_trsmC
  !!
  !! NAME
  !! xgBlock_trsmC

  subroutine xgBlock_trsmC(side,uplo,transa,diag,alpha, xgBlockA,xgBlockB)

    character      , intent(in   ) :: side
    character      , intent(in   ) :: uplo
    character      , intent(in   ) :: transa
    character      , intent(in   ) :: diag
    complex(kind=8), intent(in   ) :: alpha
    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB
    double precision :: tsec(2)

    call timab(tim_trsm,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockA%space /= SPACE_C) then
      ABI_ERROR("Not same space")
    end if

    call ztrsm(side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
      alpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)

    call timab(tim_trsm,2,tsec)

  end subroutine xgBlock_trsmC
  !!***

  !!****f* m_xg/xgBlock_colwiseCymax
  !!
  !! NAME
  !! xgBlock_colwiseCymax

  subroutine xgBlock_colwiseCymax(xgBlockA, da, xgBlockB, xgBlockW, use_gpu_cuda)

    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(in   ) :: da
    type(xgBlock_t), intent(in   ) :: xgBlockB
    type(xgBlock_t), intent(in   ) :: xgBlockW
    integer        , intent(in   ), optional :: use_gpu_cuda

    integer :: iblock
    integer          :: l_use_gpu_cuda = 0


    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockA%space /= xgBlockW%space ) then
      ABI_ERROR("Must be same space for caxmy")
    end if
    if ( xgBlockA%LDim /= xgBlockB%LDim .or. xgBlockA%LDim /= xgBlockW%LDim) then
      ABI_ERROR("Must have same LDim for caxmy")
    end if
    if ( xgBlockA%cols /= xgBlockB%cols .or. xgBlockA%cols /= xgBlockW%cols ) then
      ABI_ERROR("Must have same cols for caxmy")
    end if
    if ( da%rows /= xgBlockA%cols ) then
      ABI_ERROR("Must have same cols for caxmy")
    end if

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call compute_colwiseCymax_scalar(c_loc(xgBlockA%vecR), c_loc(da%vecR), c_loc(xgBlockB%vecR), &
          &                              c_loc(xgBlockW%vecR), xgBlockA%rows, xgBlockA%cols, xgBlockA%ldim)
      case (SPACE_C)
        call compute_colwiseCymax_cplx  (c_loc(xgBlockA%vecC), c_loc(da%vecR), c_loc(xgBlockB%vecC), &
          &                              c_loc(xgBlockW%vecC), xgBlockA%rows, xgBlockA%cols, xgBlockA%ldim)
      end select

#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
      call wrtout(std_out,"We shouldn't be here : abinit was not compiled with GPU/CUDA support (Kokkos+YAKL).")
      call abi_abort('COLL')
#endif

    else

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do shared(da,xgBlockB,xgBlockW,xgBlockA), &
        !$omp& schedule(static)
        do iblock = 1, xgBlockA%cols
          xgBlockA%vecR(:,iblock) = - da%vecR(iblock,1) * xgBlockB%vecR(:,iblock) + xgBlockW%vecR(:,iblock)
        end do
        !$omp end parallel do
      case (SPACE_C)
        !$omp parallel do shared(da,xgBlockB,xgBlockW,xgBlockA), &
        !$omp& schedule(static)
        do iblock = 1, xgBlockA%cols
          xgBlockA%vecC(:,iblock) = - da%vecR(iblock,1) * xgBlockB%vecC(:,iblock) + xgBlockW%vecC(:,iblock)
        end do
        !$omp end parallel do
      end select

    end if

  end subroutine xgBlock_colwiseCymax
  !!***

  !!****f* m_xg/xgBlock_colwiseMulR
  !!
  !! NAME
  !! xgBlock_colwiseMulR

  subroutine xgBlock_colwiseMulR(xgBlock, vec, shift, use_gpu_cuda)

    type(xgBlock_t) , intent(inout)           :: xgBlock
    double precision, intent(in   ), target   :: vec(:)
    integer,          intent(in   )           :: shift
    integer,          intent(in   ), optional :: use_gpu_cuda

    integer :: rows
    integer :: iblock
    integer :: l_use_gpu_cuda = 0

    rows = size(vec,dim=1)

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        call compute_colwiseMul_scalar_scalar(c_loc(xgBlock%vecR), c_loc(vec), &
          &                                   shift, xgBlock%rows, xgBlock%cols, &
          &                                   xgBlock%ldim, rows)
      case (SPACE_C)
        call compute_colwiseMul_cplx_scalar(c_loc(xgBlock%vecC), c_loc(vec), &
          &                                 shift, xgBlock%rows, xgBlock%cols, &
          &                                 xgBlock%ldim, rows)
      end select

#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
      call wrtout(std_out,"We shouldn't be here : abinit was not compiled with GPU/CUDA support (Kokkos+YAKL).")
      call abi_abort('COLL')
#endif

    else

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do shared(xgBlock,vec), &
        !$omp& schedule(static)
        do iblock = 1, xgBlock%cols
          xgBlock%vecR(shift+1:min(xgBlock%rows,shift+rows),iblock) = &
            xgBlock%vecR(shift+1:min(xgBlock%rows,shift+rows),iblock) * vec(1:min(xgBlock%rows-shift,rows))
        end do
      case (SPACE_C)
        !$omp parallel do shared(xgBlock,vec), &
        !$omp& schedule(static)
        do iblock = 1, xgBlock%cols
          xgBlock%vecC(shift+1:min(xgBlock%rows,shift+rows),iblock) = &
            xgBlock%vecC(shift+1:min(xgBlock%rows,shift+rows),iblock) * vec(1:min(xgBlock%rows-shift,rows))
        end do
      end select

    end if

  end subroutine xgBlock_colwiseMulR
  !!***

  !!****f* m_xg/xgBlock_colwiseMulC
  !!
  !! NAME
  !! xgBlock_colwiseMulC

  subroutine xgBlock_colwiseMulC(xgBlock, vec, shift, use_gpu_cuda)

    type(xgBlock_t), intent(inout)           :: xgBlock
    complex(kind=8), intent(in   ), target   :: vec(:)
    integer,         intent(in   )           :: shift
    integer,         intent(in   ), optional :: use_gpu_cuda

    integer :: rows
    integer :: iblock
    integer :: l_use_gpu_cuda = 0

    rows = size(vec,dim=1)

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        ABI_ERROR("Error colwiseMulC")
      case (SPACE_C)
        call compute_colwiseMul_cplx_cplx(c_loc(xgBlock%vecC), c_loc(vec), &
          &                               shift, xgBlock%rows, xgBlock%cols, &
          &                               xgBlock%ldim, rows)
      end select

#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
      call wrtout(std_out,"We shouldn't be here : abinit was not compiled with GPU/CUDA support (Kokkos+YAKL).")
      call abi_abort('COLL')
#endif

    else

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        ABI_ERROR("Error colwiseMulC")
      case (SPACE_C)
        !$omp parallel do shared(xgBlock,vec), &
        !$omp& schedule(static)
        do iblock = 1, xgBlock%cols
          xgBlock%vecC(shift+1:min(xgBlock%rows,shift+rows),iblock) = &
            xgBlock%vecC(shift+1:min(xgBlock%rows,shift+rows),iblock) * vec(1:min(xgBlock%rows-shift,rows))
        end do
      end select

    end if

  end subroutine xgBlock_colwiseMulC
  !!***

  !!****f* m_xg/xgBlock_saxpyR
  !!
  !! NAME
  !! xgBlock_saxpyR

  subroutine xgBlock_saxpyR(xgBlock1, da, xgBlock2, use_gpu_cuda)

    type(xgBlock_t),  intent(inout) :: xgBlock1
    double precision, intent(in   ) :: da
    type(xgBlock_t),  intent(in   ) :: xgBlock2
    integer        ,  intent(in   ), optional :: use_gpu_cuda

    integer      :: l_use_gpu_cuda = 0
    complex(dpc) :: da_cplx

    da_cplx = dcmplx(da,0.0_dp)

    if ( xgBlock1%space /= xgBlock2%space ) then
      ABI_ERROR("Must be same space for saxpy")
    end if
    if ( xgBlock1%LDim /= xgBlock2%LDim ) then
      ABI_ERROR("Must have same LDim for saxpy")
    end if
    if ( xgBlock1%cols /= xgBlock2%cols ) then
      ABI_ERROR("Must have same cols for saxpy")
    end if

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        call gpu_xaxpy(1, xgBlock1%cols*xgBlock1%LDim, da_cplx, c_loc(xgBlock2%vecR),1,c_loc(xgBlock1%vecR),1)
        !call daxpy(xgBlock1%cols*xgBlock1%LDim,da,xgBlock2%vecR,1,xgBlock1%vecR,1)
      case (SPACE_C)
        call gpu_xaxpy(2, xgBlock1%cols*xgBlock1%LDim, da_cplx, c_loc(xgBlock2%vecC),1,c_loc(xgBlock1%vecC),1)
        !call zaxpy(xgBlock1%cols*xgBlock1%LDim,dcmplx(da,0.d0),xgBlock2%vecC,1,xgBlock1%vecC,1)
      end select

    else

      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        call daxpy(xgBlock1%cols*xgBlock1%LDim,da,xgBlock2%vecR,1,xgBlock1%vecR,1)
      case (SPACE_C)
        call zaxpy(xgBlock1%cols*xgBlock1%LDim,dcmplx(da,0.d0),xgBlock2%vecC,1,xgBlock1%vecC,1)
      end select

    end if

  end subroutine xgBlock_saxpyR
  !!***

  !!****f* m_xg/xgBlock_saxpyC
  !!
  !! NAME
  !! xgBlock_saxpyC

  subroutine xgBlock_saxpyC(xgBlock1, da, xgBlock2, use_gpu_cuda)

    type(xgBlock_t), intent(inout) :: xgBlock1
    double complex,  intent(in   ) :: da
    type(xgBlock_t), intent(in   ) :: xgBlock2
    integer        , intent(in   ), optional :: use_gpu_cuda

    integer      :: l_use_gpu_cuda = 0

    if ( xgBlock1%space /= xgBlock2%space ) then
      ABI_ERROR("Must be same space for Saxpy")
    end if
    if ( xgBlock1%LDim /= xgBlock2%LDim ) then
      ABI_ERROR("Must have same LDim for Saxpy")
    end if
    if ( xgBlock1%cols /= xgBlock2%cols ) then
      ABI_ERROR("Must have same cols for Saxpy")
    end if
    if ( xgBlock1%space /= SPACE_C ) then
      ABI_ERROR("Not correct space")
    end if

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then
      call gpu_xaxpy(2, xgBlock1%cols*xgBlock1%LDim, da, c_loc(xgBlock2%vecC), 1, c_loc(xgBlock1%vecC), 1)
    else
      call zaxpy(xgBlock1%cols*xgBlock1%LDim, da, xgBlock2%vecC, 1, xgBlock1%vecC, 1)
    end if

  end subroutine xgBlock_saxpyC
  !!***

  !!****f* m_xg/xgBlock_add
  !!
  !! NAME
  !! xgBlock_add

  subroutine xgBlock_add(xgBlockA, xgBlockB)

    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB
    integer :: col
    integer :: row

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Must be same space for add")
    end if
    if ( xgBlockA%rows /= xgBlockB%rows ) then
      ABI_ERROR("Must have same LDim for add")
    end if
    if ( xgBlockA%cols /= xgBlockB%cols ) then
      ABI_ERROR("Must have same cols for add")
    end if

    select case(xgBlockA%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do schedule(static)
      do col = 1, xgBlockB%cols
        do row = 1, xgBlockB%rows
          xgBlockA%vecR(row,col) = xgBlockA%vecR(row,col) + xgBlockB%vecR(row,col)
        end do
      end do
      !call daxpy(xgBlockA%cols*xgBlockA%LDim,1.d0,xgBlockB%vecR,1,xgBlockA%vecR1)
    case (SPACE_C)
      !$omp parallel do schedule(static)
      do col = 1, xgBlockB%cols
        do row = 1, xgBlockB%rows
          xgBlockA%vecC(row,col) = xgBlockA%vecC(row,col) + xgBlockB%vecC(row,col)
        end do
      end do
      !call zaxpy(xgBlockA%cols*xgBlockA%LDim,1.d0,xgBlockB%vecR,1,xgBlockA%vecR1)
    end select

  end subroutine xgBlock_add
  !!***

  !!****f* m_xg/xgBlock_cshift
  !!
  !! NAME
  !! xgBlock_cshift

  subroutine xgBlock_cshift(xgBlock,nshift,shiftdim)

    type(xgBlock_t), intent(inout) :: xgBlock
    integer        , intent(in   ) :: nshift
    integer        , intent(in   ) :: shiftdim
    double precision :: tsec(2)

    call timab(tim_cshift,1,tsec)
    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      xgBlock%vecR(:,:) = cshift(xgBlock%vecR(:,:),nshift,dim=shiftdim) ! Bottom 2*blockdim lines are now at the top
    case (SPACE_C)
      xgBlock%vecC(:,:) = cshift(xgBlock%vecC(:,:),nshift,dim=shiftdim) ! Bottom 2*blockdim lines are now at the top
    end select
    call timab(tim_cshift,2,tsec)

  end subroutine xgBlock_cshift
  !!***

  !!****f* m_xg/xgBlock_colwiseNorm2
  !!
  !! NAME
  !! xgBlock_colwiseNorm2

  subroutine xgBlock_colwiseNorm2(xgBlock, dot, max_val, max_elt, min_val, min_elt, use_gpu_cuda)

    type(xgBlock_t) , intent(in   ) :: xgBlock
    type(xgBlock_t) , intent(inout) :: dot
    double precision, intent(  out), optional :: max_val
    integer         , intent(  out), optional :: max_elt
    double precision, intent(  out), optional :: min_val
    integer         , intent(  out), optional :: min_elt
    integer         , intent(in   ), optional :: use_gpu_cuda

    integer :: icol, ierr
    double precision,external :: ddot
    integer      :: l_use_gpu_cuda = 0


    if ( dot%space /= SPACE_R ) then
      ABI_ERROR("error space")
    end if

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlock%space)
      case(SPACE_R,SPACE_CR)
        call computeBatchedDotProduct_scalar(c_loc(xgBlock%vecR), c_loc(xgBlock%vecR), &
          & c_loc(dot%vecR), xgBlock%rows, xgBlock%cols, xgBlock%ldim)

      case(SPACE_C)
        call computeBatchedDotProduct_cplx_scalar(c_loc(xgBlock%vecC), c_loc(xgBlock%vecC), &
          & c_loc(dot%vecR), xgBlock%rows, xgBlock%cols, xgBlock%ldim)

      end select
      call xmpi_sum(dot%vecR,xgBlock%spacedim_comm,ierr)

      ! do reductions
      if ( present(max_val) ) then
        call computeMax_scalar(c_loc(dot%vecR(1,1)), xgBlock%cols, max_val)
      end if
      if ( present(min_val) ) then
        call computeMin_scalar(c_loc(dot%vecR(1,1)), xgBlock%cols, min_val)
      end if
      if ( present(max_elt) ) then
        call computeMaxloc_scalar(c_loc(dot%vecR(1,1)), xgBlock%cols, max_elt)
      end if
      if ( present(min_elt) ) then
        call computeMinloc_scalar(c_loc(dot%vecR(1,1)), xgBlock%cols, min_elt)
      end if

#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
#endif

    else

      select case(xgBlock%space)
      case(SPACE_R,SPACE_CR)
        !$omp parallel do shared(dot,xgBlock), &
        !$omp& schedule(static)
        do icol = 1, xgBlock%cols
          dot%vecR(icol,1) = ddot(xgBlock%rows,xgBlock%vecR(:,icol),1,xgBlock%vecR(:,icol),1)
        end do
        !$omp end parallel do
      case(SPACE_C)
        !$omp parallel do shared(dot,xgBlock), &
        !$omp& schedule(static)
        do icol = 1, xgBlock%cols
          ! Instead of calling a complex function to get only the real part of the
          ! result
          !dot%vecR(icol,1) = dble(zdotc(xgBlock%rows,xgBlock%vecC(:,icol),1,xgBlock%vecC(:,icol),1))
          ! Directely call a real function which gives what we want.
          dot%vecR(icol,1) = ddot(2*xgBlock%rows,xgBlock%vecC(:,icol),1,xgBlock%vecC(:,icol),1)
        end do
        !$omp end parallel do
      end select
      call xmpi_sum(dot%vecR,xgBlock%spacedim_comm,ierr)

      if ( present(max_val) ) then
        max_val = maxval(dot%vecR(1:xgBlock%cols,1))
      end if
      if ( present(min_val) ) then
        min_val = minval(dot%vecR(1:xgBlock%cols,1))
      end if
      if ( present(max_elt) ) then
        max_elt = maxloc(dot%vecR(1:xgBlock%cols,1),dim=1)
      end if
      if ( present(min_elt) ) then
        min_elt = minloc(dot%vecR(1:xgBlock%cols,1),dim=1)
      end if

    end if ! if l_use_gpu_cuda==1

  end subroutine xgBlock_colwiseNorm2
  !!***

  !!****f* m_xg/xgBlock_colwiseDotProduct
  !!
  !! NAME
  !! xgBlock_colwiseDotProduct

  subroutine xgBlock_colwiseDotProduct(xgBlockA,xgBlockB,dot,max_val,max_elt,min_val,min_elt,use_gpu_cuda)

    type(xgBlock_t)  , intent(in   ) :: xgBlockA
    type(xgBlock_t)  , intent(in   ) :: xgBlockB
    type(xgBlock_t)  , intent(inout) :: dot
    double precision , intent(  out), optional :: max_val
    integer          , intent(  out), optional :: max_elt
    double precision , intent(  out), optional :: min_val
    integer          , intent(  out), optional :: min_elt
    integer          , intent(in   ), optional :: use_gpu_cuda
    integer :: icol
    double precision,external :: ddot
    double complex,external :: zdotc !conjugated dot product
    integer :: l_use_gpu_cuda = 0

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        call computeBatchedDotProduct_scalar(c_loc(xgBlockA%vecR), c_loc(xgBlockB%vecR), &
          & c_loc(dot%vecR), xgBlockA%rows, xgBlockA%cols, xgBlockA%ldim)

        ! do reductions
        if ( present(max_val) ) then
          call computeMax_scalar(c_loc(dot%vecR(1,1)), xgBlockA%cols, max_val)
        end if
        if ( present(min_val) ) then
          call computeMin_scalar(c_loc(dot%vecR(1,1)), xgBlockA%cols, min_val)
        end if
        if ( present(max_elt) ) then
          call computeMaxloc_scalar(c_loc(dot%vecR(1,1)), xgBlockA%cols, max_elt)
        end if
        if ( present(min_elt) ) then
          call computeMinloc_scalar(c_loc(dot%vecR(1,1)), xgBlockA%cols, min_elt)
        end if

      case(SPACE_C)
        call computeBatchedDotProduct_cplx(c_loc(xgBlockA%vecC), c_loc(xgBlockB%vecC), &
          & c_loc(dot%vecC), xgBlockA%rows, xgBlockA%cols, xgBlockA%ldim)

        ! do reductions
        if ( present(max_val) ) then
          call computeMax_complex(c_loc(dot%vecC(1,1)), xgBlockA%cols, max_val)
        end if
        if ( present(min_val) ) then
          call computeMin_complex(c_loc(dot%vecC(1,1)), xgBlockA%cols, min_val)
        end if
        if ( present(max_elt) ) then
          call computeMaxloc_complex(c_loc(dot%vecC(1,1)), xgBlockA%cols, max_elt)
        end if
        if ( present(min_elt) ) then
          call computeMinloc_scalar(c_loc(dot%vecC(1,1)), xgBlockA%cols, min_elt)
        end if

      end select

#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
#endif

    else

      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        !$omp parallel do shared(dot,xgBlockA,xgBlockB) &
        !$omp& schedule(static)
        do icol = 1, xgBlockA%cols
          dot%vecR(icol,1) = ddot(xgBlockA%rows,xgBlockA%vecR(:,icol),1,xgBlockB%vecR(:,icol),1)
        end do
        !$omp end parallel do

        if ( present(max_val) ) then
          max_val = maxval(dot%vecR(1:xgBlockA%cols,1))
        end if
        if ( present(min_val) ) then
          min_val = minval(dot%vecR(1:xgBlockA%cols,1))
        end if
        if ( present(max_elt) ) then
          max_elt = maxloc(dot%vecR(1:xgBlockA%cols,1),dim=1)
        end if
        if ( present(min_elt) ) then
          min_elt = minloc(dot%vecR(1:xgBlockA%cols,1),dim=1)
        end if

      case(SPACE_C)
        !$omp parallel do shared(dot,xgBlockA,xgBlockB), &
        !$omp& schedule(static)
        do icol = 1, xgBlockA%cols
          dot%vecC(icol,1) = zdotc(xgBlockA%rows,xgBlockA%vecC(:,icol),1,xgBlockB%vecC(:,icol),1)
        end do
        !$omp end parallel do

        if ( present(max_val) ) then
          max_val = maxval(dble(dot%vecC(1:xgBlockA%cols,1)))
        end if
        if ( present(min_val) ) then
          min_val = minval(dble(dot%vecC(1:xgBlockA%cols,1)))
        end if
        if ( present(max_elt) ) then
          max_elt = maxloc(dble(dot%vecC(1:xgBlockA%cols,1)),dim=1)
        end if
        if ( present(min_elt) ) then
          min_elt = minloc(dble(dot%vecC(1:xgBlockA%cols,1)),dim=1)
        end if

      end select

    end if ! use_gpu_cuda

  end subroutine xgBlock_colwiseDotProduct
  !!***

  !!****f* m_xg/xgBlock_colwiseDivision
  !!
  !! NAME
  !! xgBlock_colwiseDivision

  subroutine xgBlock_colwiseDivision(xgBlockA, xgBlockB, divResult, &
    & max_val, max_elt, min_val, min_elt, use_gpu_cuda)

    type(xgBlock_t) ,      intent(in   )           :: xgBlockA
    type(xgBlock_t) ,      intent(in   )           :: xgBlockB
    type(xgBlock_t) ,      intent(inout)           :: divResult
    double precision,      intent(inout), optional :: max_val
    integer, dimension(2), intent(inout), optional, target :: max_elt
    double precision,      intent(inout), optional :: min_val
    integer, dimension(2), intent(inout), optional, target :: min_elt
    integer,               intent(in   ), optional :: use_gpu_cuda

    integer :: irow
    integer :: l_use_gpu_cuda = 0

    ! TODO: evaluate if total_size should be a 64 bit integer, i.e.
    ! does spacedim * neigenpairs be larger than 2^31 = 2. 10^9
    integer(kind=c_int32_t)  :: total_size

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      total_size = xgBlockA%rows * xgBlockA%cols

      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        call computeColwiseDivision_scalar(c_loc(xgBlockA%vecR(1,1)), &
          &                                c_loc(xgBlockB%vecR(1,1)), &
          &                                total_size,                &
          &                                c_loc(divResult%vecR(1,1)))

        ! do reductions
        if ( present(max_val) ) then
          call computeMax_scalar(c_loc(divResult%vecR(1,1)), total_size, max_val)
        end if
        if ( present(min_val) ) then
          call computeMin_scalar(c_loc(divResult%vecR(1,1)), total_size, min_val)
        end if
        if ( present(max_elt) ) then
          call computeMaxloc_scalar_2d(c_loc(divResult%vecR(1,1)), xgBlockA%rows, xgBlockA%cols, c_loc(max_elt))
        end if
        if ( present(min_elt) ) then
          call computeMinloc_scalar_2d(c_loc(divResult%vecR(1,1)), xgBlockA%rows, xgBlockA%cols, c_loc(min_elt))
        end if

      case(SPACE_C)
        call computeColwiseDivision_complex(c_loc(xgBlockA%vecC(1,1)), &
          &                                 c_loc(xgBlockB%vecC(1,1)), &
          &                                 total_size,                &
          &                                 c_loc(divResult%vecC(1,1)))

        if ( present(max_val) ) then
          call computeMax_complex(c_loc(divResult%vecC(1,1)), total_size, max_val)
        end if
        if ( present(min_val) ) then
          call computeMin_complex(c_loc(divResult%vecC(1,1)), total_size, min_val)
        end if
        if ( present(max_elt) ) then
          call computeMaxloc_complex_2d(c_loc(divResult%vecC(1,1)), xgBlockA%rows, xgBlockA%cols, c_loc(max_elt))
        end if
        if ( present(min_elt) ) then
          call computeMinloc_complex_2d(c_loc(divResult%vecC(1,1)), xgBlockA%rows, xgBlockA%cols, c_loc(min_elt))
        end if

      end select

#else
      ! we shouldn't be here, it means use_gpu_cuda was wrongly set to 1 in
      ! input parameter file
#endif

    else

      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        !$omp parallel do shared(divResult,xgBlockA,xgBlockB), &
        !$omp& schedule(static)
        do irow = 1, xgBlockA%rows
          divResult%vecR(irow,:) = xgBlockA%vecR(irow,:)/xgBlockB%vecR(irow,:)
        end do
        !$omp end parallel do

        if ( present(max_val) ) then
          max_val = maxval(dble(divResult%vecR))
        end if
        if ( present(min_val) ) then
          min_val = minval(dble(divResult%vecR))
        end if
        if ( present(max_elt) ) then
          max_elt = maxloc(dble(divResult%vecR(1:xgBlockA%rows,1:xgBlockA%cols)))
        end if
        if ( present(min_elt) ) then
          min_elt = minloc(dble(divResult%vecR(1:xgBlockA%rows,1:xgBlockA%cols)))
        end if

      case(SPACE_C)

        !$omp parallel do shared(divResult,xgBlockA,xgBlockB), &
        !$omp& schedule(static)
        do irow = 1, xgBlockA%rows
          divResult%vecC(irow,:) = xgBlockA%vecC(irow,:)/xgBlockB%vecC(irow,:)
        end do
        !$omp end parallel do

        if ( present(max_val) ) then
          max_val = maxval(dble(divResult%vecC))
        end if
        if ( present(min_val) ) then
          min_val = minval(dble(divResult%vecC))
        end if
        if ( present(max_elt) ) then
          max_elt = maxloc(dble(divResult%vecC(1:xgBlockA%rows,1:xgBlockA%cols)))
        end if
        if ( present(min_elt) ) then
          min_elt = minloc(dble(divResult%vecC(1:xgBlockA%rows,1:xgBlockA%cols)))
        end if
      end select

    end if ! use_gpu_cuda

  end subroutine xgBlock_colwiseDivision
  !!***

  !!****f* m_xg/xgBlock_scaleR
  !!
  !! NAME
  !! xgBlock_scaleR

  subroutine xgBlock_scaleR(xgBlock, val, inc, use_gpu_cuda)

    type(xgBlock_t) , intent(inout)           :: xgBlock
    double precision, intent(in   )           :: val
    integer         , intent(in   )           :: inc
    integer         , intent(in   ), optional :: use_gpu_cuda

    integer      :: i
    integer      :: l_use_gpu_cuda = 0
    complex(dpc) :: valc

    valc = dcmplx(val,0.0_dp)

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          call gpu_xscal(1, xgBlock%ldim*xgBlock%cols/inc, valc, c_loc(xgBlock%vecR), inc)
        case (SPACE_C)
          call gpu_xscal(2, xgBlock%ldim*xgBlock%cols/inc, valc, c_loc(xgBlock%vecC), inc)
        end select
      else

        ! TODO (PK) : evaluate if it is really necessary to deal with this case
        ABI_ERROR("Scaling a xgBlock when xgBlock%ldim != xgBlock%rows is not implemented for GPU. FIX ME if needed.")

      end if

    else

      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          call dscal(xgBlock%ldim*xgBlock%cols/inc,val,xgBlock%vecR,inc)
        case (SPACE_C)
          call zdscal(xgBlock%ldim*xgBlock%cols/inc,val,xgBlock%vecC,inc)
        end select
      else
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          !$omp parallel do
          do i=1,xgBlock%cols
            call dscal(xgBlock%rows/inc,val,xgBlock%vecR(:,i),inc)
          end do
        case (SPACE_C)
          !$omp parallel do
          do i=1,xgBlock%cols
            call zdscal(xgBlock%rows/inc,val,xgBlock%vecC(:,i),inc)
          end do
        end select
      end if

    end if

  end subroutine xgBlock_scaleR
  !!***

  !!****f* m_xg/xgBlock_scaleC
  !!
  !! NAME
  !! xgBlock_scaleC

  subroutine xgBlock_scaleC(xgBlock, val, inc, use_gpu_cuda)

    type(xgBlock_t), intent(inout)           :: xgBlock
    complex(kind=8), intent(in   )           :: val
    integer        , intent(in   )           :: inc
    integer        , intent(in   ), optional :: use_gpu_cuda

    integer :: i
    integer :: l_use_gpu_cuda = 0

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          ABI_ERROR("Scaling real vector with a complex not possible")
        case (SPACE_C)
          call gpu_xscal(2, xgBlock%ldim*xgBlock%cols/inc, val, c_loc(xgBlock%vecC), inc)
        end select
      else

        ! TODO (PK) : evaluate if it is really necessary to deal with this case
        ABI_ERROR("Scaling a xgBlock when xgBlock%ldim != xgBlock%rows is not implemented for GPU. FIX ME if needed.")

      end if

    else

      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          ABI_ERROR("Scaling real vector with a complex not possible")
        case (SPACE_C)
          call zscal(xgBlock%ldim*xgBlock%cols/inc,val,xgBlock%vecC,inc)
        end select
      else
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          ABI_ERROR("Scaling real vector with a complex not possible")
        case (SPACE_C)
          !$omp parallel do
          do i=1,xgBlock%cols
            call zscal(xgBlock%rows/inc,val,xgBlock%vecC(:,i),inc)
          end do
        end select
      end if

    end if

  end subroutine xgBlock_scaleC
  !!***

  !!****f* m_xg/xgBlock_getSize
  !!
  !! NAME
  !! xgBlock_getSize

  subroutine xgBlock_getSize(xgBlock, rows, cols, ldim)

    type(xgBlock_t)  , intent(in   ) :: xgBlock
    integer          , intent(  out) :: rows
    integer          , intent(  out) :: cols
    integer, optional, intent(  out) :: ldim

    rows = xgBlock%rows
    cols = xgBlock%cols

    if (present(ldim)) then
      ldim = xgBlock%ldim
    end if

  end subroutine xgBlock_getSize
  !!***

  !!****f* m_xg/xgBlock_reshape
  !!
  !! NAME
  !! xgBlock_reshape

  subroutine xgBlock_reshape(xgBlock,newShape)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t), intent(inout) :: xgBlock
    integer        , intent(in   ) :: newShape(2)
    type(c_ptr) :: cptr

    if ( xgBLock%rows*xgBlock%cols /= newShape(1)*newShape(2) ) then
      write(std_out,*) "xgBLock%rows", xgBLock%rows
      write(std_out,*) "xgBlock%cols", xgBlock%cols
      write(std_out,*) "newShape(1)", newShape(1)
      write(std_out,*) "newShape(2)", newShape(2)
      write(std_out,*) "xgBLock%rows*xgBlock%cols", xgBLock%rows*xgBlock%cols
      write(std_out,*) "newShape(1)*newShape(2)", newShape(1)*newShape(2)
      ABI_ERROR("Bad shape")
    end if

    xgBlock%LDim = newShape(1)+( (xgBlock%LDim-xgBLock%rows)* xgBlock%cols)/newShape(2)
    xgBlock%rows = newShape(1)
    xgBlock%cols = newShape(2)
    select case(xgBLock%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xgBlock%LDim,xgBlock%cols,xgBlock%vecR)
      call c_f_pointer(cptr,xgBlock%vecR,newshape)
    case (SPACE_C)
      cptr = getClocC(xgBlock%LDim,xgBlock%cols,xgBlock%vecC)
      call c_f_pointer(cptr,xgBlock%vecC,newshape)
    end select
  end subroutine xgBlock_reshape
  !!***

  !!****f* m_xg/xgBlock_zero
  !!
  !! NAME
  !! xgBlock_zero

  subroutine xgBlock_zero(xgBlock, use_gpu_cuda)

    type(xgBlock_t), intent(inout) :: xgBlock
    integer        , intent(in   ), optional :: use_gpu_cuda

    integer :: i
    integer :: l_use_gpu_cuda = 0
    integer(C_SIZE_T) :: byte_count

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    if (present(use_gpu_cuda)) then
      l_use_gpu_cuda = use_gpu_cuda
    end if

    if (l_use_gpu_cuda==1) then

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        byte_count = xgBlock%ldim * xgBlock%cols * dp
        call gpu_memset(c_loc(xgBlock%vecR), 0, byte_count)
      case (SPACE_C)
        byte_count = xgBlock%ldim * xgBlock%cols * dpc
        call gpu_memset(c_loc(xgBlock%vecC), 0, byte_count)
      end select

    else

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do
        do i = 1, xgBlock%cols
          xgBlock%vecR(:,i) = 0.d0
        end do
      case (SPACE_C)
        !$omp parallel do
        do i = 1, xgBlock%cols
          xgBlock%vecC(:,i) = dcmplx(0.d0)
        end do
      end select
    end if

  end subroutine xgBlock_zero
  !!***

  !!****f* m_xg/xgBlock_one
  !!
  !! NAME
  !! xgBlock_one

  subroutine xgBlock_one(xgBlock)

    type(xgBlock_t), intent(inout) :: xgBlock
    integer :: i

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do
      do i = 1, min(xgBlock%rows,xgBlock%cols)
        xgBlock%vecR(i,i) = 1.d0
      end do
    case (SPACE_C)
      !$omp parallel do
      do i = 1, min(xgBlock%rows,xgBlock%cols)
        xgBlock%vecC(i,i) = dcmplx(1.d0)
      end do
    end select

  end subroutine xgBlock_one
  !!***

  !!****f* m_xg/xgBlock_diagonal
  !!
  !! NAME
  !! xgBlock_diagonal

  subroutine xgBlock_diagonal(xgBlock,diag)

    type(xgBlock_t), intent(inout) :: xgBlock
    type(xgBlock_t), intent(in   ) :: diag
    integer :: i

    if ( diag%cols /= 1 .or. diag%rows/= min(xgBlock%rows,xgBlock%cols) ) then
      ABI_ERROR("Bad diagonal")
    end if

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      select case(diag%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecR(i,i) = diag%vecR(i,1)
        end do
      case (SPACE_C)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecR(i,i) = dble(diag%vecC(i,1))
        end do
      end select
    case (SPACE_C)
      select case(diag%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecC(i,i) = dcmplx(diag%vecR(i,1))
        end do
      case (SPACE_C)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecC(i,i) = diag%vecR(i,1)
        end do
      end select
    end select

  end subroutine xgBlock_diagonal
  !!***

  !!****f* m_xg/xgBlock_diagonalOnly
  !!
  !! NAME
  !! xgBlock_diagonalOnly

  subroutine xgBlock_diagonalOnly(xgBlock)

    type(xgBlock_t) , intent(inout) :: xgBlock
    type(xg_t) :: diag
    integer :: i

    if ( xgBlock%rows /= xgBlock%cols) then
      ABI_ERROR("Bad xgBlock shape")
    end if

    call xg_init(diag,space(xgBlock),xgBlock%rows,1,xgBlock%spacedim_comm)
    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do
      do i = 1, xgBlock%cols
        diag%vecR(i,1) = xgBlock%vecR(i,i)
      end do
    case (SPACE_C)
      !$omp parallel do
      do i = 1, xgBlock%cols
        diag%vecC(i,1) = xgBlock%vecC(i,i)
      end do
    end select
    call xgBlock_zero(xgBlock)
    call xgBlock_diagonal(xgBlock,diag%self)
    call xg_free(diag)
  end subroutine xgBlock_diagonalOnly
  !!***

  !!****f* m_xg/xgBlock_average
  !!
  !! NAME
  !! xgBlock_average

  subroutine xgBlock_average(xgBlock,average)

    type(XgBlock_t) , intent(in)  :: xgBlock
    double precision, intent(out) :: average
    complex(kind=8) :: averageC
    integer :: i

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      average = 0.d0
      do i = 1, xgBlock%cols
        average = average + sum(xgBlock%vecR(1:xgBlock%rows,i))
      end do
      average = average / dble(xgBlock%cols*xgBlock%rows)
    case (SPACE_C)
      averageC = dcmplx(0.d0,0.d0)
      do i = 1, xgBlock%cols
        averageC = averageC + sum(xgBlock%vecC(1:xgBlock%rows,i))
      end do
      averageC = averageC / dble(xgBlock%cols*xgBlock%rows)
      average = dble(averageC)
    end select

  end subroutine xgBlock_average
  !!***

  !!****f* m_xg/xgBlock_deviation
  !!
  !! NAME
  !! xgBlock_deviation

  subroutine xgBlock_deviation(xgBlock,deviation)

    type(XgBlock_t) , intent(in)  :: xgBlock
    double precision, intent(out) :: deviation
    complex(kind=8) :: deviationC
    double precision :: average
    integer :: i

    call xgBlock_average(xgBlock,average)
    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      deviation = 0.d0
      do i = 1, xgBlock%cols
        deviation = deviation + sum((xgBlock%vecR(1:xgBlock%rows,i)-average)*(xgBlock%vecR(1:xgBlock%rows,i)-average))
      end do
      deviation = sqrt( deviation / dble(xgBlock%cols*xgBlock%rows) )
    case (SPACE_C)
      deviationC = dcmplx(0.d0,0.d0)
      do i = 1, xgBlock%cols
        deviationC = deviationC + sum((xgBlock%vecC(1:xgBlock%rows,i)-average)*(xgBlock%vecC(1:xgBlock%rows,i)-average))
      end do
      deviationC = deviationC / dble(xgBlock%cols*xgBlock%rows)
      deviation = abs(deviationC)
    end select
  end subroutine xgBlock_deviation
  !!***

  !!****f* m_xg/xgBlock_print
  !!
  !! NAME
  !! xgBlock_print

  subroutine xgBlock_print(xgBlock,outunit)

    type(xgBlock_t), intent(in) :: xgBlock
    integer, intent(in) :: outunit
    integer :: i, j
    character(len=4) :: ccols
    character(len=50) :: fstring


    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      write(ccols,'(i4)') xgBlock%cols
      fstring = '(1x,'//trim(adjustl(ccols))//'ES22.14)'
      do i = 1, xgBlock%rows
        write(outunit,fstring) (/ (xgBlock%vecR(i,j), j = 1, xgBlock%cols) /)
      end do
    case (SPACE_C)
      write(ccols,'(i4)') xgBlock%cols
      fstring = '(1x,2(1x,'//trim(adjustl(ccols))//'ES22.14))'
      do i = 1, xgBlock%rows
        write(outunit,fstring) (/ (xgBlock%vecC(i,j), j = 1, xgBlock%cols) /)
      end do
    end select
  end subroutine xgBlock_print
  !!***

  !!****f* m_xg/xg_finalize
  !!
  !! NAME
  !! xg_finalize

  subroutine xg_finalize()

#if defined HAVE_GPU && defined HAVE_YAKL

    if ( associated(iwork) ) then
      ABI_FREE_MANAGED(iwork)
    end if
    if ( associated(rwork) ) then
      ABI_FREE_MANAGED(rwork)
    end if
    if ( associated(cwork) ) then
      ABI_FREE_MANAGED(cwork)
    end if

#else

    if ( allocated(iwork) ) then
      ABI_FREE(iwork)
    end if
    if ( allocated(rwork) ) then
      ABI_FREE(rwork)
    end if
    if ( allocated(cwork) ) then
      ABI_FREE(cwork)
    end if

#endif

    liwork = 0
    lrwork = 0
    lcwork = 0
  end subroutine xg_finalize
  !!***

  !!****f* m_xg/xg_associated
  !!
  !! NAME
  !! xg_associated

  function xg_associated(xgB) result (tf)

    type(xgBlock_t), intent(inout) :: xgB
    logical :: tf

    if ( associated(xgB%vecR) ) then
      tf = .TRUE.
    end if

    if ( associated(xgB%vecC) ) then
      tf = .FALSE.
    end if

  end function xg_associated
  !!***

end module m_xg
!!***
