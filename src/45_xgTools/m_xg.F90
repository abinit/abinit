!{\src2tex{textfont=tt}}
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
!!  Copyright (C) 2016-2018 ABINIT group (J. Bieder)
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

#include "abi_common.h"

module m_xg

  use m_errors
  use m_abicore

  use defs_basis, only : std_err, std_out
  use m_time,     only: timab
  use m_xmpi,     only : xmpi_sum


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
  integer, allocatable,  save, private :: iwork(:)
  double precision, allocatable, save, private :: rwork(:)
  complex(kind= 8), allocatable, save, private :: cwork(:)

  type, public :: xgBlock_t
    integer, private :: space
    integer, private :: rows
    integer, private :: LDim
    integer, private :: cols
    character, public :: trans
    character, public :: normal
    integer, private :: spacedim_comm
    double precision, pointer, private :: vecR(:,:) => null()
    complex(kind=8) , pointer, private :: vecC(:,:) => null()
  end type xgBlock_t

  type, public :: xg_t
    integer, private :: space
    integer, private :: rows
    integer, private :: cols
    character, public :: trans
    character, public :: normal
    integer, private :: spacedim_comm
    double precision, allocatable, private :: vecR(:,:)
    complex(kind=8) , allocatable, private :: vecC(:,:)
    type(xgBlock_t), public :: self
  end type xg_t

  interface xgBlock_gemm
    module procedure xgBlock_gemmR
    module procedure xgBlock_gemmC
  end interface

  interface xgBlock_colwiseMul
    module procedure xgBlock_colwiseMulR
    module procedure xgBlock_colwiseMulC
  end interface

  interface xgBlock_trsm
    module procedure xgBlock_trsmR
    module procedure xgBlock_trsmC
  end interface

  interface xgBlock_scale
    module procedure xgBlock_scaleR
    module procedure xgBlock_scaleC
  end interface

  interface checkResize
    module procedure checkResizeI
    module procedure checkResizeR
    module procedure checkResizeC
  end interface

  public :: space
  public :: cols
  private :: getClocR
  private :: getClocC
  private :: checkResize

  public :: xg_init
  public :: xg_set
  public :: xg_get
  public :: xg_setBlock
  public :: xg_free

  public :: xgBlock_setBlock
  public :: xgBlock_set
  public :: xgBlock_map
  public :: xgBlock_reverseMap
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
  public :: xgBlock_colwiseCaxmy
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'checkResizeI'
!End of the abilint section

    integer, allocatable, intent(inout) :: array(:)
    integer, intent(inout) :: current_dim
    integer, intent(in   )  :: asked_dim

    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
    end if
  end subroutine checkResizeI
!!***

!!****f* m_xg/checkResizeR
!!
!! NAME
!! checkResizeR

  subroutine checkResizeR(array,current_dim,asked_dim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'checkResizeR'
!End of the abilint section

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
  end subroutine checkResizeR
!!***

!!****f* m_xg/checkResizeC
!!
!! NAME
!! checkResizeC

  subroutine checkResizeC(array,current_dim,asked_dim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'checkResizeC'
!End of the abilint section

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
  end subroutine checkResizeC
!!***

!!****f* m_xg/getClocR
!!
!! NAME
!! getClocR

  function getClocR(rows,cols,array) result(cptr)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getClocR'
!End of the abilint section

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
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getClocC'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xg_init'
!End of the abilint section

    type(xg_t), intent(inout) :: xg
    integer   , intent(in   ) :: space
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    integer   , optional, intent(in) :: comm

    if ( rows < 1 ) then
      MSG_ERROR("rows < 1 ")
    endif
    if ( cols < 1 ) then
      MSG_ERROR("cols < 1 ")
    end if

    select case (space)
    case (SPACE_R,SPACE_CR)
      if ( allocated(xg%vecR) ) then
        ABI_FREE(xg%vecR)
      end if
      ABI_MALLOC(xg%vecR,(1:rows,1:cols))
      !xg%vecR(:,:) = 0.d0
      xg%trans = 't'
    case (SPACE_C)
      if ( allocated(xg%vecC) ) then
        ABI_FREE(xg%vecC)
      end if
      ABI_MALLOC(xg%vecC,(1:rows,1:cols))
      !xg%vecC(:,:) = dcmplx(0.d0)
      xg%trans = 'c'
    case default
      MSG_ERROR("Invalid space")
    end select

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xg_set'
!End of the abilint section

    type(xg_t), intent(inout) :: xg
    double precision, intent(in) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_set,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      MSG_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xg%cols ) then
      MSG_WARNING("Ignore some columns, input array to large")
    endif

    select case (xg%space)
    case (SPACE_R)
      do col = 1, min(cols,xg%cols-shift_col)
        xg%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
      end do
    case (SPACE_CR)
      if ( xg%rows /= 2*rows ) then
        MSG_ERROR("Bad number of rows")
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_set'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    double precision, intent(in) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_set,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      MSG_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xgBlock%cols ) then
      MSG_WARNING("Block Ignore some columns, input array to large")
    endif

    select case (xgBlock%space)
    case (SPACE_R)
      do col = 1, min(cols,xgBlock%cols-shift_col)
        xgBlock%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
      end do
    case (SPACE_CR)
      if ( xgBlock%rows /= 2*rows ) then
        MSG_ERROR("Bad number of rows")
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
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_map'
!End of the abilint section

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
        MSG_ERROR("Bad size for real array")
      end if
      cptr = getClocR(size(array,dim=1),size(array,dim=2),array)
      call c_f_pointer(cptr,xgBlock%vecR,(/ rows, cols /))
    xgBlock%trans = 't'
    case ( SPACE_C )
      if ( fullsize/2 < cols*rows .or. mod(fullsize/2,rows) /= 0) then
        MSG_ERROR("Bad size for complex array")
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
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_reverseMap'
!End of the abilint section

    type(xgBlock_t) , intent(inout) :: xgBlock
    double precision, pointer, intent(inout) :: array(:,:)
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    type(c_ptr) :: cptr

    select case (xgBlock%space)
    case ( SPACE_R,SPACE_CR )
      if ( xgBlock%cols*xgBlock%Ldim < cols*rows ) then
          MSG_ERROR("Bad reverseMapping")
      end if
      cptr = getClocR(xgBlock%Ldim,xgBlock%cols,xgBlock%vecR(:,:))
      call c_f_pointer(cptr,array,(/ rows, cols /))
    case ( SPACE_C )
      if ( xgBlock%cols*xgBlock%Ldim < cols*rows ) then
          MSG_ERROR("Bad complex reverseMapping")
      end if
      cptr = getClocC(xgBlock%Ldim,xgBlock%cols,xgBlock%vecC(:,:))
      call c_f_pointer(cptr,array,(/ 2*rows, cols /))
    end select

  end subroutine xgBlock_reverseMap
!!***

!!****f* m_xg/xg_get
!!
!! NAME
!! xg_get

  subroutine xg_get(xg,array,shift_col,rows)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xg_get'
!End of the abilint section

    type(xg_t), intent(inout) :: xg
    double precision, intent(out) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_get,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      MSG_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xg%cols ) then
      MSG_WARNING("Ignore some columns, input array to large")
    endif

    select case (xg%space)
    case (SPACE_R)
      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xg%vecR(1:rows,shift_col+col)
      end do
    case (SPACE_CR)
      if ( xg%rows /= 2*rows ) then
        MSG_ERROR("Bad number of rows")
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_get'
!End of the abilint section

    type(xgBlock_t), intent(in   ) :: xgBlock
    double precision, intent(out) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
    double precision :: tsec(2)

    call timab(tim_get,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      MSG_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xgBlock%cols ) then
      MSG_ERROR("Ignore some columns, input array to large")
    endif

    select case (xgBlock%space)
    case (SPACE_R)
      do col = 1, min(cols,xgBlock%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xgBlock%vecR(1:rows,shift_col+col)
      end do
    case (SPACE_CR)
      if ( xgBlock%rows /= 2*rows ) then
        MSG_ERROR("Bad number of rows")
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

  subroutine xg_setBlock(xg,Xgblock, fcol, rows, cols)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xg_setBlock'
!End of the abilint section

    type(xg_t), intent(inout) :: xg
    type(xgBlock_t), intent(inout) :: xgBlock
    integer, intent(in) :: fcol
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    type(c_ptr) :: cptr

    if ( (fcol+cols-1 ) > xg%cols ) then
      MSG_ERROR("Too many columns")
    endif
    if ( rows > xg%rows ) then
      MSG_ERROR("Too many rows")
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

  subroutine xgBlock_setBlock(xgBlock1,xgBlock2, fcol, rows, cols)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_setBlock'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock1
    type(xgBlock_t), intent(inout) :: xgBlock2
    integer, intent(in) :: fcol
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    type(c_ptr) :: cptr

    if ( (fcol+cols-1 ) > xgblock1%cols ) then
      MSG_ERROR("Too many columns")
    endif
    if ( rows > xgblock1%rows ) then
      MSG_ERROR("Too many rows")
    end if

    xgBlock2%space = xgBlock1%space
    xgBlock2%rows = rows
    xgBlock2%LDim = xgBlock1%LDim
    xgBlock2%cols = cols
    xgBlock2%trans = xgBlock1%trans
    xgBlock2%normal = xgBlock1%normal
    xgBlock2%spacedim_comm= xgBlock1%spacedim_comm

    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xgBlock1%LDim,xgBlock1%cols,xgBlock1%vecR(:,fcol:fcol+cols-1))
      call c_f_pointer(cptr,xgBlock2%vecR,(/ xgBlock2%LDim,cols /))
    case(SPACE_C)
      cptr = getClocC(xgBlock1%LDim,xgBlock1%cols,xgBlock1%vecC(:,fcol:fcol+cols-1))
      call c_f_pointer(cptr,xgBlock2%vecC,(/ xgBlock2%LDim,cols /))
    end select

  end subroutine xgBlock_setBlock
!!***

!!****f* m_xg/xg_free
!!
!! NAME
!! xg_free

  subroutine xg_free(xg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xg_free'
!End of the abilint section

    type(xg_t), intent(inout) :: xg

    if ( allocated(xg%vecR) ) then
      ABI_FREE(xg%vecR)
    end if

    if ( allocated(xg%vecC) ) then
      ABI_FREE(xg%vecC)
    end if
  end subroutine xg_free
!!***

!!****f* m_xg/space
!!
!! NAME
!! space

  function space(xgBlock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'space'
!End of the abilint section

    type(xgBlock_t), intent(in) :: xgBlock
    integer :: space
    space = xgBlock%space
  end function space
!!***

!!****f* m_xg/cols
!!
!! NAME
!! cols

  function cols(xgBlock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cols'
!End of the abilint section

    type(xgBlock_t), intent(in) :: xgBlock
    integer :: cols
    cols = xgBlock%cols
  end function cols
!!***

!!****f* m_xg/xgBlock_copy
!!
!! NAME
!! xgBlock_copy

  subroutine xgBlock_copy(xgBlock1, xgBlock2, inc1, inc2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_copy'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock1
    type(xgBlock_t), intent(inout) :: xgBlock2
    integer, optional :: inc1
    integer, optional :: inc2
    integer :: incx
    integer :: incy
    integer :: size1
    integer :: size2
    integer :: size
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)
    incx = 1; if ( present(inc1) ) incx = inc1
    incy = 1; if ( present(inc2) ) incy = inc2
    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if
    !if ( xgBlock1%LDim*xgBlock1%cols/incx /= xgBlock2%LDim*xgBlock2%cols/incy ) then
    !  MSG_ERROR("Number of element different")
    !end if

    size1 = xgBlock1%LDim*xgBlock1%cols/incx ; if ( size1 * incx < xgBlock1%LDim*xgBlock1%cols ) size1 = size1+1
    size2 = xgBlock2%LDim*xgBlock2%cols/incy ; if ( size2 * incy < xgBlock2%LDim*xgBlock2%cols ) size2 = size2+1
    size = min(size1,size2)

    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      call dcopy(size,xgBlock1%vecR,incx,xgBlock2%vecR,incy)
    case(SPACE_C)
      call zcopy(size,xgBlock1%vecC,incx,xgBlock2%vecC,incy)
    end select
    call timab(tim_copy,2,tsec)

  end subroutine xgBlock_copy
!!***

!!****f* m_xg/xgBlock_pack
!!
!! NAME
!! xgBlock_pack

  subroutine xgBlock_pack(xgBlock1,xgBlock2,uplo)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_pack'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock1
    type(xgBlock_t), intent(inout) :: xgBlock2
    character, intent(in) :: uplo
    integer :: j
    integer :: i
    integer :: col
    type(c_ptr) :: cptr
    double precision, pointer :: subR(:)
    complex(kind=8), pointer :: subC(:)
    double precision :: tsec(2)

    call timab(tim_pack,1,tsec)
    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Both blocks must be the same space")
    end if

    if ( xgBlock1%Ldim /= xgBlock1%rows ) then
      MSG_ERROR("Cannot pack when ldim /= rows")
    end if

    if ( xgBlock1%Ldim /= xgBlock1%cols ) then
      MSG_ERROR("Cannot pack when cols /= rows")
    end if

    if ( xgBlock1%rows*(xgBlock1%rows+1)/2 > xgBlock2%Ldim*xgBlock2%cols ) then
      MSG_ERROR("Not enought memory in destination")
    end if

    ! make a fake pointer to pack in a 1-D array instead of 2
    ! This allows to directly use the lapack conventions of transformation
    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xgBlock2%Ldim,xgBlock2%cols,xgBlock2%vecR(:,:))
      call c_f_pointer(cptr,subR,(/ (xgBlock2%cols*(xgBlock2%cols+1))/2 /))
    case (SPACE_C)
      cptr = getClocC(xgBlock2%Ldim,xgBlock2%cols,xgBlock2%vecC(:,:))
      call c_f_pointer(cptr,subC,(/ (xgBlock2%cols*(xgBlock2%cols+1))/2 /))
    end select

    select case(uplo)
    case ('u','U')
      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        do j = 1, xgBlock1%cols
          col = (j*(j-1))/2
          do i = 1, j
            subR(i+col) = xgBlock1%vecR(i,j)
          end do
        end do
      case (SPACE_C)
        do j = 1, xgBlock1%cols
          col = (j*(j-1))/2
          do i = 1, j
            subC(i+col) = xgBlock1%vecC(i,j)
          end do
        end do
      end select

    case ('l','L')
      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        do j = 1, xgBlock1%cols
          col = ((2*xgBlock1%cols-j)*(j-1))/2
          do i = j, xgBlock1%cols
            subR(i+col) = xgBlock1%vecR(i,j)
          end do
        end do
      case (SPACE_C)
        do j = 1, xgBlock1%cols
          col = ((2*xgBlock1%cols-j)*(j-1))/2
          do i = j, xgBlock1%cols
            subC(i+col) = xgBlock1%vecC(i,j)
          end do
        end do
      end select
    case default
      MSG_ERROR("Error for packing matrix")
    end select
    call timab(tim_pack,2,tsec)

  end subroutine xgBlock_pack
!!***

!!****f* m_xg/xgBlock_gemmR
!!
!! NAME
!! xgBlock_gemmR

  subroutine xgBlock_gemmR(transa, transb, alpha, xgBlock1, xgBlock2, beta, xgBlock3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_gemmR'
!End of the abilint section

    character, intent(in) :: transa
    character, intent(in) :: transb
    double precision, intent(in) :: alpha
    type(xgBlock_t), intent(in) :: xgBlock1
    type(xgBlock_t), intent(in) :: xgBlock2
    double precision, intent(in) :: beta
    type(xgBlock_t), intent(inout) :: xgBlock3
    complex(kind=8) :: calpha
    complex(kind=8) :: cbeta
    integer :: K
    double precision :: tsec(2)

    call timab(tim_gemm,1,tsec)
    if ( xgBlock1%space /= xgBlock2%space .or. xgBlock2%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if

    if ( transa == 'n' ) then
      K = xgBlock1%cols
    else
      K = xgBlock1%rows
    end if

    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      call dgemm(transa,transb,xgBlock3%rows, xgBlock3%cols, K, &
        alpha,xgBlock1%vecR, xgBlock1%LDim, &
        xgBlock2%vecR, xgBlock2%LDim, beta,xgBlock3%vecR,xgBlock3%LDim)
      if ( transa == xgBlock1%trans .and. (beta) < 1d-10) then
        call xmpi_sum(xgBlock3%vecR,xgBlock3%spacedim_comm,K)
      end if
    case(SPACE_C)
      calpha = dcmplx(alpha,0.d0)
      cbeta = dcmplx(beta,0.d0)
      call zgemm(transa,transb,xgBlock3%rows, xgBlock3%cols, K, &
        calpha,xgBlock1%vecC, xgBlock1%LDim, &
        xgBlock2%vecC, xgBlock2%LDim, cbeta,xgBlock3%vecC,xgBlock3%LDim)
      if ( xgBlock3%spacedim_comm/= -1 .and. transa == xgBlock3%trans .and. abs(beta) < 1d-10 ) then
        call xmpi_sum(xgBlock3%vecC,xgBlock3%spacedim_comm,K)
      end if
    end select

    call timab(tim_gemm,2,tsec)

  end subroutine xgBlock_gemmR
!!***

!!****f* m_xg/xgBlock_gemmC
!!
!! NAME
!! xgBlock_gemmC

  subroutine xgBlock_gemmC(transa, transb, alpha, xgBlock1, xgBlock2, beta, xgBlock3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_gemmC'
!End of the abilint section

    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(kind=8), intent(in) :: alpha
    type(xgBlock_t), intent(in) :: xgBlock1
    type(xgBlock_t), intent(in) :: xgBlock2
    complex(kind=8), intent(in) :: beta
    type(xgBlock_t), intent(inout) :: xgBlock3
    integer :: K
    double precision :: tsec(2)

    call timab(tim_gemm,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space .or. xgBlock2%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if
    if ( xgBlock1%space /= SPACE_C ) then
      MSG_ERROR("Not correct space")
    end if

    if ( transa == 'n' ) then
      K = xgBlock1%cols
    else
      K = xgBlock1%rows
    end if

    call zgemm(transa,transb,xgBlock3%rows, xgBlock3%cols, K, &
      alpha,xgBlock1%vecC, xgBlock1%LDim, &
      xgBlock2%vecC, xgBlock2%LDim, beta,xgBlock3%vecC,xgBlock3%LDim)
    if ( xgBlock3%spacedim_comm/= -1 .and. transa == xgBlock1%trans .and. abs(beta) < 1.d-10 ) then
      call xmpi_sum(xgBlock3%vecC,xgBlock3%spacedim_comm,K)
    end if

    call timab(tim_gemm,2,tsec)

  end subroutine xgBlock_gemmC
!!***

!!****f* m_xg/xgBlock_potrf
!!
!! NAME
!! xgBlock_potrf

  subroutine xgBlock_potrf(xgBlock,uplo,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_potrf'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    character      , intent(in   ) :: uplo
    integer        , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_potrf,1,tsec)

    if ( xgBlock%rows /= xgBlock%cols ) then
      MSG_ERROR("Matrix should be a square matrixx")
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
  subroutine xgBlock_heev(jobz,uplo,xgBlock1,xgBlock3,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_heev'
!End of the abilint section

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock3
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_heev,1,tsec)

    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,8*xgBlock1%rows)

      call dsyev(jobz,uplo,xgBlock1%cols, &
        xgBlock1%vecR,xgBlock1%LDim, &
      xgBlock3%vecR, &
      rwork, lrwork,info)


    case (SPACE_C)
      call checkResize(rwork,lrwork,3*xgBlock1%cols-2)
      call checkResize(cwork,lcwork,lrwork)

      call zheev(jobz,uplo,xgBlock1%cols, &
        xgBlock1%vecC,xgBlock1%LDim, &
      xgBlock3%vecR, &
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

  subroutine xgBlock_heevd(jobz,uplo,xgBlock1,xgBlock3, info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_heevd'
!End of the abilint section

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock3
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_heevd,1,tsec)

    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlock1%rows+3)

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlock1%rows*xgBlock1%rows+6*xgBlock1%rows+1)

      call dsyevd(jobz,uplo,xgBlock1%cols, &
        xgBlock1%vecR,xgBlock1%LDim, &
      xgBlock3%vecR, rwork, lrwork, &
      iwork, liwork,info)

    case (SPACE_C)
      call checkResize(cwork,lcwork,xgBlock1%rows*xgBlock1%rows+2*xgBlock1%rows)
      call checkResize(rwork,lrwork,2*xgBlock1%rows*xgBlock1%rows+5*xgBlock1%rows+1)

      call zheevd(jobz,uplo,xgBlock1%cols, &
        xgBlock1%vecC,xgBlock1%LDim, &
      xgBlock3%vecR, &
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
  subroutine xgBlock_hpev(jobz,uplo,xgBlock1,xgBlock3,xgBlock4,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hpev'
!End of the abilint section

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock3
    type(xgBlock_t) , intent(inout) :: xgBlock4
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpev,1,tsec)

    if ( xgBlock1%space /= xgBlock4%space ) then
      MSG_ERROR("Not same space")
    end if

    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,3*xgBlock4%cols)

      call dspev(jobz,uplo,xgBlock4%cols, &
        xgBlock1%vecR, xgBlock3%vecR, xgBlock4%vecR, xgBlock4%Ldim, &
      rwork, info)


    case (SPACE_C)
      call checkResize(cwork,lcwork,2*xgBlock4%cols-1)
      call checkResize(rwork,lrwork,3*xgBlock4%cols-2)

      call zhpev(jobz,uplo,xgBlock4%cols, &
        xgBlock1%vecC, xgBlock3%vecR, xgBlock4%vecC, xgBlock4%Ldim, &
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

  subroutine xgBlock_hpevd(jobz,uplo,xgBlock1,xgBlock3,xgBlock4,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hpevd'
!End of the abilint section

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock3
    type(xgBlock_t) , intent(inout) :: xgBlock4
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpevd,1,tsec)

    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    if ( xgBlock1%space /= xgBlock4%space ) then
      MSG_ERROR("Block 1 and 3 must have the same space")
    end if

    call checkResize(iwork,liwork,5*xgBlock4%rows+3)
    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,xgBlock4%rows*xgBlock4%rows+6*xgBlock4%rows+1)

      call dspevd(jobz,uplo,xgBlock4%cols, &
        xgBlock1%vecR, xgBlock3%vecR, xgBlock4%vecR, xgBlock4%Ldim, &
        rwork, lrwork, iwork, liwork,info)


    case (SPACE_C)
      call checkResize(cwork,lcwork,2*xgBlock4%rows)
      call checkResize(rwork,lrwork,2*xgBlock4%rows*xgBlock4%rows+5*xgBlock4%rows+1)

      call zhpevd(jobz,uplo,xgBlock4%cols, &
        xgBlock1%vecC, xgBlock3%vecR, xgBlock4%vecC, xgBlock4%Ldim, &
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

  subroutine xgBlock_hegv(itype, jobz, uplo, xgBlock1, xgBlock2, xgBlock3, info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hegv'
!End of the abilint section

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    type(xgBlock_t) , intent(inout) :: xgBlock3
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hegv,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if
    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlock1%rows*xgBlock1%rows+6*xgBlock1%rows+1)

      call dsygv(itype, jobz, uplo, xgBlock1%rows, xgBlock1%vecR, xgBlock1%ldim, &
        xgBlock2%vecR, xgBlock2%ldim, xgBlock3%vecR, rwork, lrwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,2*xgBlock1%rows-1)
      call checkResize(rwork,lrwork,3*xgBlock1%rows-2)

      call zhegv(itype, jobz, uplo, xgBlock1%rows, xgBlock1%vecC, xgBlock1%ldim,&
        xgBlock2%vecC, xgBlock2%ldim, xgBlock3%vecR, cwork, lcwork, &
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

  subroutine xgBlock_hegvx(itype,jobz,range,uplo,xgBlock1,xgBlock2,vl,vu,il,iu,abstol,xgBlock3,xgBlock4,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hegvx'
!End of the abilint section

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: range
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    double precision, intent(in   ) :: vl
    double precision, intent(in   ) :: vu
    integer         , intent(in   ) :: il
    integer         , intent(in   ) :: iu
    double precision, intent(in   ) :: abstol
    type(xgBlock_t) , intent(inout) :: xgBlock3
    type(xgBlock_t) , intent(inout) :: xgBlock4
    integer         , intent(  out) :: info
    integer :: neigen
    integer, allocatable :: ifail(:)
    double precision :: tsec(2)

    call timab(tim_hegvx,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space .or. xgBlock1%space /= xgBlock4%space ) then
      MSG_ERROR("Not same space")
    end if
    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlock1%rows)

    ABI_MALLOC(ifail,(xgBlock1%rows))
    ifail = 0

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,8*xgBlock1%rows)

      call dsygvx(itype,jobz,range,uplo,xgBlock1%rows, &
        xgBlock1%vecR,xgBlock1%LDim,xgBlock2%vecR,xgBlock2%LDim, &
      vl,vu,il,iu,abstol,&
      neigen,xgBlock3%vecR, xgBlock4%vecR, xgBlock4%LDim, &
      rwork, lrwork,iwork,ifail,info)

    case (SPACE_C)
      call checkResize(rwork,lrwork,7*xgBlock1%rows)
      call checkResize(cwork,lcwork,lrwork)

      call zhegvx(itype,jobz,range,uplo,xgBlock1%rows, &
        xgBlock1%vecC,xgBlock1%LDim,xgBlock2%vecC,xgBlock2%LDim, &
      vl,vu,il,iu,abstol,&
      neigen,xgBlock3%vecR, xgBlock4%vecC, xgBlock4%LDim, &
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

  subroutine xgBlock_hegvd(itype, jobz, uplo, xgBlock1, xgBlock2, xgBlock3, info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hegvd'
!End of the abilint section

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    type(xgBlock_t) , intent(inout) :: xgBlock3
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hegvd,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if
    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlock1%rows+3)

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlock1%rows*xgBlock1%rows+6*xgBlock1%rows+1)

      call dsygvd(itype, jobz, uplo, xgBlock1%rows, xgBlock1%vecR, xgBlock1%ldim, &
        xgBlock2%vecR, xgBlock2%ldim, xgBlock3%vecR, rwork, lrwork, iwork, liwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,xgBlock1%rows*xgBlock1%rows+2*xgBlock1%rows)
      call checkResize(rwork,lrwork,2*(xgBlock1%rows*xgBlock1%rows)+5*xgBlock1%rows+1)

      call zhegvd(itype, jobz, uplo, xgBlock1%rows, xgBlock1%vecC, xgBlock1%ldim,&
        xgBlock2%vecC, xgBlock2%ldim, xgBlock3%vecR, cwork, lcwork, &
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

  subroutine xgBlock_hpgv(itype, jobz, uplo, xgBlock1, xgBlock2, xgBlock3, xgBlock4,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hpgv'
!End of the abilint section

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    type(xgBlock_t) , intent(inout) :: xgBlock3
    type(xgBlock_t) , intent(inout) :: xgBlock4
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpgv,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if

    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,3*xgBlock4%rows)

      call dspgv(itype, jobz, uplo, xgBlock4%rows, xgBlock1%vecR, xgBlock2%vecR, &
        xgBlock3%vecR, xgBlock4%vecR, xgBlock4%Ldim, rwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,2*xgBlock4%rows-1)
      call checkResize(rwork,lrwork,3*xgBlock4%rows-2)

      call zhpgv(itype, jobz, uplo, xgBlock1%rows, xgBlock1%vecC, xgBlock2%vecC, &
        xgBlock3%vecR, xgBlock4%vecC, xgBlock4%Ldim, cwork, rwork, info)

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

  subroutine xgBlock_hpgvx(itype,jobz,range,uplo,xgBlock1,xgBlock2,vl,vu,il,iu,abstol,xgBlock3,xgBlock4,info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hpgvx'
!End of the abilint section

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: range
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    double precision, intent(in   ) :: vl
    double precision, intent(in   ) :: vu
    integer         , intent(in   ) :: il
    integer         , intent(in   ) :: iu
    double precision, intent(in   ) :: abstol
    type(xgBlock_t) , intent(inout) :: xgBlock3
    type(xgBlock_t) , intent(inout) :: xgBlock4
    integer         , intent(  out) :: info
    integer :: neigen
    integer, allocatable :: ifail(:)
    double precision :: tsec(2)

    call timab(tim_hpgvx,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space .or. xgBlock1%space /= xgBlock4%space ) then
      MSG_ERROR("Not same space")
    end if
    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlock4%rows)

    ABI_MALLOC(ifail,(xgBlock4%rows))
    ifail = 0

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,8*xgBlock4%rows)

      call dspgvx(itype,jobz,range,uplo,xgBlock4%rows, &
        xgBlock1%vecR,xgBlock2%vecR, &
      vl,vu,il,iu,abstol,&
      neigen,xgBlock3%vecR, xgBlock4%vecR, xgBlock4%LDim, &
      rwork, iwork, ifail, info)

    case (SPACE_C)
      call checkResize(rwork,lrwork,7*xgBlock1%rows)
      call checkResize(cwork,lcwork,2*xgBlock4%rows)

      call zhpgvx(itype,jobz,range,uplo,xgBlock4%rows, &
        xgBlock1%vecC,xgBlock2%vecC, &
      vl,vu,il,iu,abstol,&
      neigen,xgBlock3%vecR, xgBlock4%vecC, xgBlock4%LDim, &
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

  subroutine xgBlock_hpgvd(itype, jobz, uplo, xgBlock1, xgBlock2, xgBlock3, xgBlock4, info)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_hpgvd'
!End of the abilint section

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    type(xgBlock_t) , intent(inout) :: xgBlock3
    type(xgBlock_t) , intent(inout) :: xgBlock4
    integer         , intent(  out) :: info
    double precision :: tsec(2)

    call timab(tim_hpgvd,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if
    if ( xgBlock3%space /= SPACE_R ) then
      MSG_ERROR("Block3 must be real")
    end if

    call checkResize(iwork,liwork,5*xgBlock4%rows+3)

    select case(xgBlock1%space)

    case (SPACE_R,SPACE_CR)
      call checkResize(rwork,lrwork,2*xgBlock4%rows*xgBlock4%rows+6*xgBlock4%rows+1)

      call dspgvd(itype, jobz, uplo, xgBlock4%rows, xgBlock1%vecR, xgBlock2%vecR, &
        xgBlock3%vecR, xgBlock4%vecR, xgBlock4%Ldim, &
        rwork, lrwork, iwork, liwork, info)

    case (SPACE_C)

      call checkResize(cwork,lcwork,2*xgBlock4%rows)
      call checkResize(rwork,lrwork,2*(xgBlock4%rows*xgBlock4%rows)+5*xgBlock4%rows+1)

      call zhpgvd(itype, jobz, uplo, xgBlock4%rows, xgBlock1%vecC, xgBlock2%vecC, &
        xgBlock3%vecR, xgBlock4%vecC, xgBlock4%Ldim, &
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

  subroutine xgBlock_trsmR(side,uplo,transa,diag,alpha, xgBlock1,xgBlock2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_trsmR'
!End of the abilint section

    character       , intent(in   ) :: side
    character       , intent(in   ) :: uplo
    character       , intent(in   ) :: transa
    character       , intent(in   ) :: diag
    double precision, intent(in   ) :: alpha
    type(xgBlock_t) , intent(inout) :: xgBlock1
    type(xgBlock_t) , intent(inout) :: xgBlock2
    complex(kind=8) :: calpha
    double precision :: tsec(2)

    call timab(tim_trsm,1,tsec)
    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Not same space")
    end if

    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      call dtrsm(side,uplo,transa,diag,xgBlock2%rows,xgBlock2%cols, &
        alpha,xgBlock1%vecR,xgBlock1%LDim,xgBlock2%vecR,xgBlock2%LDim)
    case (SPACE_C)
      calpha = dcmplx(alpha,0.d0)
      call ztrsm(side,uplo,transa,diag,xgBlock2%rows,xgBlock2%cols, &
        calpha,xgBlock1%vecC,xgBlock1%LDim,xgBlock2%vecC,xgBlock2%LDim)
    end select

    call timab(tim_trsm,2,tsec)

  end subroutine xgBlock_trsmR
!!***

!!****f* m_xg/xgBlock_trsmC
!!
!! NAME
!! xgBlock_trsmC

  subroutine xgBlock_trsmC(side,uplo,transa,diag,alpha, xgBlock1,xgBlock2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_trsmC'
!End of the abilint section

    character      , intent(in   ) :: side
    character      , intent(in   ) :: uplo
    character      , intent(in   ) :: transa
    character      , intent(in   ) :: diag
    complex(kind=8), intent(in   ) :: alpha
    type(xgBlock_t), intent(inout) :: xgBlock1
    type(xgBlock_t), intent(inout) :: xgBlock2
    double precision :: tsec(2)

    call timab(tim_trsm,1,tsec)

    if ( xgBlock1%space /= xgBlock2%space .or. xgBlock1%space /= SPACE_C) then
      MSG_ERROR("Not same space")
    end if

    call ztrsm(side,uplo,transa,diag,xgBlock2%rows,xgBlock2%cols, &
      alpha,xgBlock1%vecC,xgBlock1%LDim,xgBlock2%vecC,xgBlock2%LDim)

    call timab(tim_trsm,2,tsec)

  end subroutine xgBlock_trsmC
!!***

!!****f* m_xg/xgBlock_colwiseCaxmy
!!
!! NAME
!! xgBlock_colwiseCaxmy

  subroutine xgBlock_colwiseCaxmy(xgBlock1, da, xgBlock2,xgBlock3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_colwiseCaxmy'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock1
    type(xgBlock_t), intent(in   ) :: da
    type(xgBlock_t), intent(in   ) :: xgBlock2
    type(xgBlock_t), intent(in   ) :: xgBlock3
    integer :: iblock

    if ( xgBlock1%space /= xgBlock2%space .or. xgBlock1%space /= xgBlock3%space ) then
      MSG_ERROR("Must be same space for caxmy")
    end if
    if ( xgBlock1%LDim /= xgBlock2%LDim .or. xgBlock1%LDim /= xgBlock3%LDim) then
      MSG_ERROR("Must have same LDim for caxmy")
    end if
    if ( xgBlock1%cols /= xgBlock2%cols .or. xgBlock1%cols /= xgBlock3%cols ) then
      MSG_ERROR("Must have same cols for caxmy")
    end if
    if ( da%rows /= xgBlock1%cols ) then
      MSG_ERROR("Must have same cols for caxmy")
    end if

    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do shared(da,xgBlock2,xgBlock3,xgBlock1), &
      !$omp& schedule(static)
      do iblock = 1, xgBlock1%cols
        xgBlock1%vecR(:,iblock) = - da%vecR(iblock,1) * xgBlock2%vecR(:,iblock) + xgBlock3%vecR(:,iblock)
      end do
      !$omp end parallel do
    case (SPACE_C)
      !$omp parallel do shared(da,xgBlock2,xgBlock3,xgBlock1), &
      !$omp& schedule(static)
      do iblock = 1, xgBlock1%cols
        xgBlock1%vecC(:,iblock) = - da%vecR(iblock,1) * xgBlock2%vecC(:,iblock) + xgBlock3%vecC(:,iblock)
      end do
      !$omp end parallel do
    end select

  end subroutine xgBlock_colwiseCaxmy
!!***

!!****f* m_xg/xgBlock_colwiseMulR
!!
!! NAME
!! xgBlock_colwiseMulR

  subroutine xgBlock_colwiseMulR(xgBlock, vec, shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_colwiseMulR'
!End of the abilint section

    type(xgBlock_t) , intent(inout) :: xgBlock
    double precision, intent(in   ) :: vec(:)
    integer, intent(in   )          :: shift
    integer :: rows
    integer :: iblock

    rows = size(vec,dim=1)

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

  end subroutine xgBlock_colwiseMulR
!!***

!!****f* m_xg/xgBlock_colwiseMulC
!!
!! NAME
!! xgBlock_colwiseMulC

  subroutine xgBlock_colwiseMulC(xgBlock, vec, shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_colwiseMulC'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    complex(kind=8), intent(in   ) :: vec(:)
    integer, intent(in   )         :: shift
    integer :: rows
    integer :: iblock

    rows = size(vec,dim=1)

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      MSG_ERROR("Error colwiseMulC")
    case (SPACE_C)
      !$omp parallel do shared(xgBlock,vec), &
      !$omp& schedule(static)
      do iblock = 1, xgBlock%cols
        xgBlock%vecC(shift+1:min(xgBlock%rows,shift+rows),iblock) = &
        xgBlock%vecC(shift+1:min(xgBlock%rows,shift+rows),iblock) * vec(1:min(xgBlock%rows-shift,rows))
      end do
    end select

  end subroutine xgBlock_colwiseMulC
!!***

!!****f* m_xg/xgBlock_add
!!
!! NAME
!! xgBlock_add

  subroutine xgBlock_add(xgBlock1, xgBlock2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_add'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock1
    type(xgBlock_t), intent(inout) :: xgBlock2
    integer :: col
    integer :: row

    if ( xgBlock1%space /= xgBlock2%space ) then
      MSG_ERROR("Must be same space for add")
    end if
    if ( xgBlock1%rows /= xgBlock2%rows ) then
      MSG_ERROR("Must have same LDim for add")
    end if
    if ( xgBlock1%cols /= xgBlock2%cols ) then
      MSG_ERROR("Must have same cols for add")
    end if

    select case(xgBlock1%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do schedule(static)
      do col = 1, xgBlock2%cols
        do row = 1, xgBlock2%rows
          xgBlock1%vecR(row,col) = xgBlock1%vecR(row,col) + xgBlock2%vecR(row,col)
        end do
      end do
      !call daxpy(xgBlock1%cols*xgBlock1%LDim,1.d0,xgBlock2%vecR,1,xgBlock1%vecR1)
    case (SPACE_C)
      !$omp parallel do schedule(static)
      do col = 1, xgBlock2%cols
        do row = 1, xgBlock2%rows
          xgBlock1%vecC(row,col) = xgBlock1%vecC(row,col) + xgBlock2%vecC(row,col)
        end do
      end do
      !call zaxpy(xgBlock1%cols*xgBlock1%LDim,1.d0,xgBlock2%vecR,1,xgBlock1%vecR1)
    end select

  end subroutine xgBlock_add
!!***

!!****f* m_xg/xgBlock_cshift
!!
!! NAME
!! xgBlock_cshift

  subroutine xgBlock_cshift(xgBlock,nshift,shiftdim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_cshift'
!End of the abilint section

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

  subroutine xgBlock_colwiseNorm2(xgBlock,dot,max_val,max_elt,min_val,min_elt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_colwiseNorm2'
!End of the abilint section

    type(xgBlock_t) , intent(in   ) :: xgBlock
    type(xgBlock_t) , intent(inout) :: dot
    double precision, intent(  out), optional :: max_val
    integer         , intent(  out), optional :: max_elt
    double precision, intent(  out), optional :: min_val
    integer         , intent(  out), optional :: min_elt
    integer :: icol
    double precision,external :: ddot


    if ( dot%space /= SPACE_R ) then
      MSG_ERROR("error space")
    end if

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
        ! resuld
        !dot%vecR(icol,1) = dble(zdotc(xgBlock%rows,xgBlock%vecC(:,icol),1,xgBlock%vecC(:,icol),1))
        ! Directely call a real function which gives what we want.
        dot%vecR(icol,1) = ddot(2*xgBlock%rows,xgBlock%vecC(:,icol),1,xgBlock%vecC(:,icol),1)
      end do
      !$omp end parallel do
    end select
    call xmpi_sum(dot%vecR,xgBlock%spacedim_comm,icol)

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
  end subroutine xgBlock_colwiseNorm2
!!***

!!****f* m_xg/xgBlock_scaleR
!!
!! NAME
!! xgBlock_scaleR

  subroutine xgBlock_scaleR(xgBlock, val, inc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_scaleR'
!End of the abilint section

    type(xgBlock_t) , intent(inout) :: xgBlock
    double precision, intent(in   ) :: val
    integer         , intent(in   ) :: inc
    integer :: i

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

  end subroutine xgBlock_scaleR
!!***

!!****f* m_xg/xgBlock_scaleC
!!
!! NAME
!! xgBlock_scaleC

  subroutine xgBlock_scaleC(xgBlock, val, inc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_scaleC'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    complex(kind=8), intent(in   ) :: val
    integer         , intent(in   ) :: inc
    integer :: i

    if ( xgBlock%ldim .eq. xgBlock%rows ) then
      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        MSG_ERROR("Scaling real vector with a complex not possible")
      case (SPACE_C)
        call zscal(xgBlock%ldim*xgBlock%cols/inc,val,xgBlock%vecC,inc)
      end select
    else
      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        MSG_ERROR("Scaling real vector with a complex not possible")
      case (SPACE_C)
        !$omp parallel do
        do i=1,xgBlock%cols
          call zscal(xgBlock%rows/inc,val,xgBlock%vecC(:,i),inc)
        end do
      end select
    end if

  end subroutine xgBlock_scaleC
!!***

!!****f* m_xg/xgBlock_getSize
!!
!! NAME
!! xgBlock_getSize

  subroutine xgBlock_getSize(xgBlock, rows, cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_getSize'
!End of the abilint section

    type(xgBlock_t), intent(in   ) :: xgBlock
    integer        , intent(  out) :: rows
    integer        , intent(  out) :: cols

    rows = xgBlock%rows
    cols = xgBlock%cols
  end subroutine xgBlock_getSize
!!***

!!****f* m_xg/xgBlock_reshape
!!
!! NAME
!! xgBlock_reshape

  subroutine xgBlock_reshape(xgBlock,newShape)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_reshape'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    integer        , intent(in   ) :: newShape(2)
    type(c_ptr) :: cptr

    if ( xgBLock%rows*xgBlock%cols /= newShape(1)*newShape(2) ) then
      MSG_ERROR("Bad shape")
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

  subroutine xgBlock_zero(xgBlock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_zero'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    integer :: i

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do
      do i = 1, xgBlock%cols
        xgBlock%vecR(:,i) = 0.d0
      end do
    case (SPACE_C)
      !$omp parallel do
      do i = 1, xgBlock%cols
        xgBlock%vecC = dcmplx(0.d0)
      end do
    end select

  end subroutine xgBlock_zero
!!***

!!****f* m_xg/xgBlock_one
!!
!! NAME
!! xgBlock_one

  subroutine xgBlock_one(xgBlock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_one'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_diagonal'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: xgBlock
    type(xgBlock_t), intent(in   ) :: diag
    integer :: i

    if ( diag%cols /= 1 .or. diag%rows/= min(xgBlock%rows,xgBlock%cols) ) then
      MSG_ERROR("Bad diagonal")
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_diagonalOnly'
!End of the abilint section

    type(xgBlock_t) , intent(inout) :: xgBlock
    type(xg_t) :: diag
    integer :: i

    if ( xgBlock%rows /= xgBlock%cols) then
      MSG_ERROR("Bad xgBlock shape")
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_average'
!End of the abilint section

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
      averageC = cmplx(0.d0,0.d0)
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_deviation'
!End of the abilint section

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
      deviationC = cmplx(0.d0,0.d0)
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgBlock_print'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xg_finalize'
!End of the abilint section

    if ( allocated(iwork) ) then
      ABI_FREE(iwork)
    end if
    if ( allocated(rwork) ) then
      ABI_FREE(rwork)
    end if
    if ( allocated(cwork) ) then
      ABI_FREE(cwork)
    end if
    liwork = 0
    lrwork = 0
    lcwork = 0
  end subroutine xg_finalize
!!***

end module m_xg
!!***
