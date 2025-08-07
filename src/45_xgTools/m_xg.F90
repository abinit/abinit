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
!!  Copyright (C) 2016-2025 ABINIT group (J. Bieder, MS, L. Baguet)
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

  use, intrinsic :: iso_c_binding, only: c_loc, c_double, c_double_complex, c_int32_t, c_size_t, c_ptr

  use m_errors
  use m_abicore
  use defs_basis
  use m_time, only : timab
  USE_MPI
  use m_xmpi
  use m_xomp
  use m_abi_linalg

#if defined(HAVE_GPU)
  use m_gpu_toolbox
#endif

#if defined(HAVE_KOKKOS)
  use m_xg_kokkos
#endif

#if defined HAVE_YAKL
 use gator_mod
#endif

  implicit none

  private

  integer, parameter, public :: SPACE_R =  1
  integer, parameter, public :: SPACE_C =  2
  integer, parameter, public :: SPACE_CR = 3

  integer, parameter, public :: SMALL2BIG =  1
  integer, parameter, public :: BIG2SMALL = -1

  integer, parameter, public :: COLS2ROWS =  1
  integer, parameter, public :: ROWS2COLS = -1

  integer, parameter :: tim_gemm_blas  = 1670
  integer, parameter :: tim_trsm       = 1671
  integer, parameter :: tim_potrf      = 1672
  integer, parameter :: tim_zero       = 1673
  integer, parameter :: tim_zero_im_g0 = 1674
!  integer, parameter :: tim_set        = 1673
!  integer, parameter :: tim_get        = 1674
  integer, parameter :: tim_heev       = 1675
  integer, parameter :: tim_heevd      = 1676
  integer, parameter :: tim_hpev       = 1677
  integer, parameter :: tim_hpevd      = 1678
  integer, parameter :: tim_hegv       = 1679
  integer, parameter :: tim_hegvx      = 1680
  integer, parameter :: tim_hegvd      = 1681
  integer, parameter :: tim_hpgv       = 1682
  integer, parameter :: tim_hpgvx      = 1683
  integer, parameter :: tim_hpgvd      = 1684
  integer, parameter :: tim_copy       = 1685
  integer, parameter :: tim_cshift     = 1686
  integer, parameter :: tim_pack       = 1687
  integer, parameter :: tim_gemm_mpi   = 1688
  integer, parameter :: tim_apply_diag = 1689
  integer, parameter :: tim_invertri   = 1696

  integer, parameter :: tim_scale       = 2000
  integer, parameter :: tim_colw_dot    = 2001
  integer, parameter :: tim_colw_mul    = 2002
  integer, parameter :: tim_colw_cymax  = 2003
  integer, parameter :: tim_colw_div    = 2004
  integer, parameter :: tim_colw_norm2  = 2005
  integer, parameter :: tim_saxpy       = 2006
  integer, parameter :: tim_minmax      = 2007
  integer, parameter :: tim_partialcopy = 2008
  integer, parameter :: tim_gemmcyclic  = 2009
  integer, parameter :: tim_yxmax       = 2010
  integer, parameter :: tim_ymax        = 2011
  integer, parameter :: tim_add         = 2012
  integer, parameter :: tim_add_diag    = 2013
  integer, parameter :: tim_invert      = 2014
  integer, parameter :: tim_invert_sy   = 2015

  integer, save, private :: lrwork = 0
  integer, save, private :: lcwork = 0
  integer, save, private :: liwork = 0
  integer,                        allocatable, save, private :: iwork(:)
  real(kind=c_double),            allocatable, save, private :: rwork(:)
  complex(kind=c_double_complex), allocatable, save, private :: cwork(:)

  type, public :: xgBlock_t
    integer, private :: space
    integer, private :: rows
    integer, private :: LDim
    integer, private :: cols
    integer, private :: me_g0
    character, public :: trans
    character, public :: normal
    integer, private :: spacedim_comm
    integer, private :: gpu_option
    real(kind=c_double)            , ABI_CONTIGUOUS pointer, private :: vecR(:,:) => null()
    complex(kind=c_double_complex) , ABI_CONTIGUOUS pointer, private :: vecC(:,:) => null()
  end type xgBlock_t

  type, public :: xg_t
    integer, private :: space
    integer, private :: rows
    integer, private :: cols
    integer, private :: me_g0
    character, public :: trans
    character, public :: normal
    integer, private :: spacedim_comm
    !FIXME Settle this
    real(kind=c_double)            , ABI_CONTIGUOUS pointer, private :: vecR(:,:) => null()
    complex(kind=c_double_complex) , ABI_CONTIGUOUS pointer, private :: vecC(:,:) => null()
    integer, private :: gpu_option
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

  interface xgBlock_reverseMap
    module procedure xgBlock_reverseMapR
    module procedure xgBlock_reverseMapC
  end interface xgBlock_reverseMap

  interface xgBlock_reverseMap_1d
    module procedure xgBlock_reverseMap_1dR
    module procedure xgBlock_reverseMap_1dC
  end interface xgBlock_reverseMap_1d

  interface checkResize
    module procedure checkResizeI
    module procedure checkResizeR
    module procedure checkResizeC
  end interface checkResize

  public :: space
  public :: cols
  public :: rows
  public :: comm
  public :: gpu_option
  public :: me_g0
  public :: xgBlock_setComm
  private :: getClocR
  private :: getClocC
  private :: checkResize

  public :: xg_init
  public :: xg_set ! LB-06/03/24: Be careful, this routine is not used (so not tested)
  public :: xg_get ! LB-06/03/24: Be careful, this routine is not used (so not tested)
  public :: xg_setBlock
  public :: xg_free

  public :: xg_associated

  public :: xgBlock_setBlock
  public :: xgBlock_set ! LB-06/03/24: Be careful, this routine is not used (so not tested)
  public :: xgBlock_map
  public :: xgBlock_map_1d
  public :: xgBlock_reverseMap
  public :: xgBlock_reverseMap_1d
  public :: xgBlock_prefetch_async
  public :: xgBlock_get ! LB-06/03/24: Be careful, this routine is not used (so not tested)
  public :: xgBlock_copy
  public :: xgBlock_partialcopy
  public :: xgBlock_pack
  public :: xgBlock_getSize

  public :: xgBlock_check
  public :: xgBlock_check_gpu_option

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
  public :: xgBlock_trmmR
  public :: xgBlock_add
  public :: xgBlock_cshift
  public :: xgBlock_colwiseNorm2
  public :: xgBlock_colwiseDotProduct
  public :: xgBlock_colwiseDivision
  public :: xgBlock_ymax
  public :: xgBlock_yxmax
  public :: xgBlock_colwiseCymax
  public :: xgBlock_saxpy
  public :: xgBlock_colwiseMul
  public :: xgBlock_scale
  public :: xgBlock_transpose
  public :: xgBlock_r2c
  public :: xgBlock_c2r

  public :: xgBlock_apply_diag
  public :: xgBlock_add_diag

  public :: xgBlock_mpi_sum
  public :: xgBlock_mpi_send
  public :: xgBlock_mpi_isend
  public :: xgBlock_mpi_recv
  public :: xgBlock_gemm_mpi_cyclic_permutation

  public :: xgBlock_invert
  public :: xgBlock_invert_sy
  public :: xgBlock_invert_tri
  public :: xgBlock_yxpa

  public :: xgBlock_zero
  public :: xgBlock_zerotri
  public :: xgBlock_zero_im_g0
  public :: xgBlock_one
  public :: xgBlock_diagonal
  public :: xgBlock_diagonalOnly

  public :: xgBlock_minmax
  public :: xgBlock_average
  public :: xgBlock_deviation

  public :: xgBlock_reshape
  public :: xgBlock_reshape_spinor
  public :: xgBlock_free_reshape
  public :: xgBlock_print
  public :: xgBlock_getId
  public :: xgBlock_get_im_g0
  public :: xgBlock_copy_from_gpu
  public :: xgBlock_copy_to_gpu
  public :: xg_finalize

contains
  !!***

  !!****f* m_xg/checkResizeI
  !!
  !! NAME
  !! checkResizeI

  subroutine checkResizeI(array,current_dim,asked_dim)

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

  end subroutine checkResizeI
  !!***

  !!****f* m_xg/checkResizeR
  !!
  !! NAME
  !! checkResizeR

  subroutine checkResizeR(array,current_dim,asked_dim)

    double precision, allocatable, intent(inout) :: array(:)
    integer, intent(inout) :: current_dim
    integer, intent(in   ) :: asked_dim

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

    complex(kind=8), allocatable, intent(inout) :: array(:)
    integer, intent(inout)  :: current_dim
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
    use, intrinsic :: iso_c_binding
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    double precision, target, intent(in) :: array(rows,cols)
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
    complex(kind=8), target, intent(in) :: array(rows,cols)
    type(c_ptr) :: cptr
    cptr = c_loc(array)
  end function getClocC
  !!***

  !!****f* m_xg/xg_init
  !!
  !! NAME
  !! xg_init

  subroutine xg_init(xg, space, rows, cols, comm, me_g0, gpu_option)

    type(xg_t), target, intent(inout) :: xg
    integer   , intent(in   ) :: space
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    integer   , optional, intent(in) :: comm, me_g0, gpu_option
    integer                   :: l_gpu_option,fact
#if defined HAVE_GPU
    integer(kind=c_int32_t), parameter :: izero = 0
#endif
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), pointer :: xg__vecC(:,:)
    real(dp), pointer :: xg__vecR(:,:)
#endif

    if ( rows < 1 ) then
      ABI_ERROR("rows < 1 ")
    endif
    if ( cols < 1 ) then
      ABI_ERROR("cols < 1 ")
    end if

    ! if optional parameter is present, use it
    ! else use default value, i.e. don't use GPU
    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    fact = 1 ; if (space==SPACE_CR) fact = 2

    if (l_gpu_option==ABI_GPU_KOKKOS) then
#if defined HAVE_GPU && defined HAVE_YAKL
      select case (space)
      case (SPACE_R,SPACE_CR)
        if ( associated(xg%vecR) ) then
          ABI_FREE_MANAGED(xg%vecR)
        end if
        ABI_MALLOC_MANAGED_BOUNDS(xg%vecR,(/fact*rows,cols/), (/1,1/))
        xg%trans = 't'
      case (SPACE_C)
        if ( associated(xg%vecC) ) then
          ABI_FREE_MANAGED(xg%vecC)
        end if
        ABI_MALLOC_MANAGED_BOUNDS(xg%vecC,(/rows,cols/), (/1,1/))
        xg%trans = 'c'
      case default
        ABI_ERROR("Invalid space")
      end select
#endif

    else if (l_gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      select case (space)

      case (SPACE_R,SPACE_CR)
        if ( associated(xg%vecR) ) then
          !$OMP TARGET EXIT DATA MAP(delete:xg%vecR)
          ABI_FREE(xg%vecR)
        end if
        ABI_MALLOC(xg%vecR,(1:fact*rows,1:cols))
        xg%trans = 't'
#if defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        !$OMP TARGET ENTER DATA MAP(alloc:xg%vecR)
#else
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        xg__vecR => xg%vecR
        !$OMP TARGET ENTER DATA MAP(alloc:xg__vecR)
#endif

      case (SPACE_C)
        if ( associated(xg%vecC) ) then
          !$OMP TARGET EXIT DATA MAP(delete:xg%vecC)
          ABI_FREE(xg%vecC)
        end if
        ABI_MALLOC(xg%vecC,(1:rows,1:cols))
        xg%trans = 'c'
#if defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        !$OMP TARGET ENTER DATA MAP(alloc:xg%vecC)
#else
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        xg__vecC => xg%vecC
        !$OMP TARGET ENTER DATA MAP(alloc:xg__vecC)
#endif

      case default
        ABI_ERROR("Invalid space")
      end select
#endif

    else if ( l_gpu_option==ABI_GPU_DISABLED .or. l_gpu_option==ABI_GPU_LEGACY ) then

      select case (space)
      case (SPACE_R,SPACE_CR)
        if ( associated(xg%vecR) ) then
          ABI_FREE(xg%vecR)
        end if
        ABI_MALLOC(xg%vecR,(1:fact*rows,1:cols))
        xg%trans = 't'
      case (SPACE_C)
        if ( associated(xg%vecC) ) then
          ABI_FREE(xg%vecC)
        end if
        ABI_MALLOC(xg%vecC,(1:rows,1:cols))
        xg%trans = 'c'
      case default
        ABI_ERROR("Invalid space")
      end select

    else
        ABI_ERROR("Invalid gpu_option")
    end if

    xg%space = space
    xg%normal = 'n'
    xg%cols = cols
    xg%rows = rows
    xg%spacedim_comm = xmpi_comm_null
    xg%gpu_option = l_gpu_option
    xg%me_g0 = -1

    if ( present(comm) ) xg%spacedim_comm = comm
    if ( present(me_g0) ) then
      if (me_g0/=-1.and.me_g0/=0.and.me_g0/=1) then
        ABI_ERROR('Bad value of me_g0 in xg_init')
      else
        xg%me_g0 = me_g0
      end if
    end if

    call xg_setBlock(xg,xg%self,rows,cols)
    call xgBlock_zero(xg%self)

  end subroutine xg_init
  !!***

  !LB-06/03/24: Be careful, xg_set is not used anywhere in the code, so not tested...
  !!****f* m_xg/xg_set
  !!
  !! NAME
  !! xg_set

  subroutine xg_set(xg,array,shift_col,rows)

    type(xg_t), target, intent(inout) :: xg
    double precision, intent(in) :: array(:,:)
    integer, intent(in) :: shift_col
    integer, intent(in) :: rows
    integer :: cols
    integer :: col
!    double precision :: tsec(2)

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    complex(dpc), pointer :: xg__vecC(:,:)
    real(dp), pointer :: xg__vecR(:,:)
#endif

!    call timab(tim_set,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xg%cols ) then
      ABI_WARNING("Ignore some columns, input array to large")
    endif

    if(xg%gpu_option == ABI_GPU_OPENMP) then
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      select case (xg%space)
      case (SPACE_R)
        xg__vecR => xg%vecR
        do col = 1, min(cols,xg%cols-shift_col)
          xg%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
        end do
        !$OMP TARGET UPDATE TO(xg__vecR)
      case (SPACE_CR)
        xg__vecR => xg%vecR
        if ( xg%rows /= 2*rows ) then
          ABI_ERROR("Bad number of rows")
        end if

        do col = 1, min(cols,xg%cols-shift_col)
          xg%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
          xg%vecR(rows+1:2*rows,shift_col+col) = array(2,(col-1)*rows+1:col*rows)
        end do
        !$OMP TARGET UPDATE TO(xg__vecR)
      case (SPACE_C)
        xg__vecC => xg%vecC
        do col = 1, min(cols,xg%cols-shift_col)
          xg%vecC(1:rows,shift_col+col) = dcmplx(array(1,(col-1)*rows+1:col*rows), &
            array(2,(col-1)*rows+1:col*rows))
        end do
        !$OMP TARGET UPDATE TO(xg__vecC)
      end select
#endif
    else
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
    end if

!    call timab(tim_set,2,tsec)

  end subroutine xg_set
  !!***

  !LB-06/03/24: Be careful, xgBlock_set is not used anywhere in the code, so not tested...
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
!    double precision :: tsec(2)

!    call timab(tim_set,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xgBlock%cols ) then
      ABI_WARNING("Block Ignore some columns, input array to large")
    endif

    if(xgBlock%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      select case (xgBlock%space)
      case (SPACE_R)
        call xgBlock_copy_from_gpu(xgBlock)
        do col = 1, min(cols,xgBlock%cols-shift_col)
          xgBlock%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
        end do
        call xgBlock_copy_to_gpu(xgBlock)
      case (SPACE_CR)
        if ( xgBlock%rows /= 2*rows ) then
          ABI_ERROR("Bad number of rows")
        end if

        call xgBlock_copy_from_gpu(xgBlock)
        do col = 1, min(cols,xgBlock%cols-shift_col)
          xgBlock%vecR(1:rows,shift_col+col) = array(1,(col-1)*rows+1:col*rows)
          xgBlock%vecR(rows+1:2*rows,shift_col+col) = array(2,(col-1)*rows+1:col*rows)
        end do
        call xgBlock_copy_to_gpu(xgBlock)
      case (SPACE_C)
        call xgBlock_copy_from_gpu(xgBlock)
        do col = 1, min(cols,xgBlock%cols-shift_col)
          xgBlock%vecC(1:rows,shift_col+col) = dcmplx(array(1,(col-1)*rows+1:col*rows), &
            array(2,(col-1)*rows+1:col*rows))
        end do
        call xgBlock_copy_to_gpu(xgBlock)
      end select
#endif
    else
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
    end if

!    call timab(tim_set,2,tsec)

  end subroutine xgBlock_set
  !!***

  !!****f* m_xg/xgBlock_map
  !!
  !! NAME
  !! xgBlock_map

  subroutine xgBlock_map(xgBlock,array,space,rows,cols,comm,me_g0,gpu_option)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(inout) :: xgBlock
    double precision, target, intent(in) :: array(:,:)
    integer   , intent(in   ) :: space
    integer   , intent(in   ) :: rows
    integer   , intent(in   ) :: cols
    integer   , optional, intent(in) :: comm
    integer   , optional, intent(in) :: me_g0
    integer   , optional, intent(in) :: gpu_option
    integer :: fullsize,fact
    type(c_ptr) :: cptr

    fullsize = size(array)
    fact = 1 ; if (space==SPACE_CR) fact = 2
    select case (space)
    case ( SPACE_R,SPACE_CR )
      if ( fullsize < fact*cols*rows .or. mod(fullsize,fact*rows) /= 0) then
        ABI_ERROR("Bad size for real array")
      end if
      cptr = getClocR(size(array,dim=1),size(array,dim=2),array)
      call c_f_pointer(cptr,xgBlock%vecR,(/ fact*rows, cols /))
      xgBlock%trans = 't'
    case ( SPACE_C )
      if ( fullsize/2 < cols*rows .or. mod(fullsize/2,rows) /= 0) then
        ABI_ERROR("Bad size for complex array")
      end if
      cptr = getClocR(size(array,dim=1),size(array,dim=2),array)
      call c_f_pointer(cptr,xgBlock%vecC,(/ rows, cols /))

      xgBlock%trans = 'c'
    case default
      ABI_ERROR('Bad space in xgBlock_map')
    end select

    xgBlock%space = space
    xgBlock%rows = rows
    xgBlock%LDim = rows
    xgBlock%cols = cols
    xgBlock%normal = 'n'
    xgBlock%spacedim_comm = xmpi_comm_null
    xgBlock%me_g0 = -1
    xgBlock%gpu_option = ABI_GPU_DISABLED

    if ( present(comm) ) xgBlock%spacedim_comm = comm

    if ( present(me_g0) ) then
      if (me_g0/=-1.and.me_g0/=0.and.me_g0/=1) then
        ABI_ERROR('Bad value of me_g0 in xg_init')
      else
        xgBlock%me_g0 = me_g0
      end if
    end if

    if ( present(gpu_option) ) xgBlock%gpu_option = gpu_option
#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
    if ( xgBlock%gpu_option == ABI_GPU_OPENMP ) then
      ABI_CHECK(xomp_target_is_present(c_loc(array)), "Mapped array isn't... mapped with OpenMP")
    end if
#endif
    if ( xgBlock%gpu_option /= ABI_GPU_DISABLED .and. xgBlock%gpu_option /= ABI_GPU_LEGACY .and. &
         xgBlock%gpu_option /= ABI_GPU_OPENMP   .and. xgBlock%gpu_option /= ABI_GPU_KOKKOS ) then
       ABI_ERROR('Bad GPU option in xgBlock_map')
    end if

  end subroutine xgBlock_map
  !!***

  !!****f* m_xg/xgBlock_map_1d
  !!
  !! NAME
  !! xgBlock_map_1d

    subroutine xgBlock_map_1d(xgBlock,array,space,rows,comm,me_g0,gpu_option)
      use iso_c_binding
      type(xgBlock_t) , intent(inout) :: xgBlock
      double precision, intent(in), target :: array(:)
      integer   , intent(in   ) :: space
      integer   , intent(in   ) :: rows
      integer   , optional, intent(in) :: comm
      integer   , optional, intent(in) :: me_g0
      integer   , optional, intent(in) :: gpu_option

      integer :: comm_,me_g0_,gpu_option_
      type(c_ptr) :: cptr
      real(dp), pointer :: array_ptr(:,:) => NULL()

      ! Trick the with C to change rank of arrays (:) to (:,:)
      cptr = c_loc(array)
      call c_f_pointer(cptr,array_ptr,(/ rows,1 /))

      comm_ =xmpi_comm_null
      if (present(comm)) then
        comm_=comm
      end if

      me_g0_=-1
      if (present(me_g0)) then
        me_g0_=me_g0
      end if

      gpu_option_=ABI_GPU_DISABLED
      if (present(gpu_option)) then
        gpu_option_=gpu_option
      end if

      call xgBlock_map(xgBlock,array_ptr,space,rows,1,comm=comm_,me_g0=me_g0_,gpu_option=gpu_option_)

    end subroutine xgBlock_map_1d
  !!***

  !!****f* m_xg/xgBlock_reverseMapR
  !!
  !! NAME
  !! xgBlock_reverseMapR

  subroutine xgBlock_reverseMapR(xgBlock,array,rows,cols)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(in) :: xgBlock
    double precision, pointer, intent(inout) :: array(:,:)
    integer,optional,intent(in) :: rows
    integer,optional,intent(in) :: cols
    type(c_ptr) :: cptr

    integer :: rows_,cols_,fact

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    rows_ = fact*xgBlock%ldim
    cols_ = xgBlock%cols
    if (present(rows)) rows_=fact*rows
    if (present(cols)) cols_=cols

    select case (xgBlock%space)
    case ( SPACE_R,SPACE_CR )
      if ( xgBlock%cols*fact*xgBlock%Ldim < cols_*rows_ ) then
        write(std_out,*) xgBlock%cols,xgBlock%Ldim,cols,rows
        write(std_out,*) xgBlock%cols*xgBlock%Ldim,cols*rows
        ABI_ERROR("Bad reverseMapping")
      end if
      cptr = getClocR(fact*xgBlock%Ldim,xgBlock%cols,xgBlock%vecR(:,:))
      call c_f_pointer(cptr,array,(/ rows_, cols_ /))
    case ( SPACE_C )
      if ( xgBlock%cols*xgBlock%Ldim < cols_*rows_ ) then
        ABI_ERROR("Bad complex reverseMapping")
      end if
      cptr = getClocC(xgBlock%Ldim,xgBlock%cols,xgBlock%vecC(:,:))
      call c_f_pointer(cptr,array,(/ 2*rows_, cols_ /))
    end select

  end subroutine xgBlock_reverseMapR
  !!***

  !!****f* m_xg/xgBlock_reverseMap_1dR
  !!
  !! NAME
  !! xgBlock_reverseMap_1dR

  subroutine xgBlock_reverseMap_1dR(xgBlock,array,array_dim)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(in) :: xgBlock
    double precision, pointer, intent(inout) :: array(:)
    integer,optional,intent(in) :: array_dim
    type(c_ptr) :: cptr

    integer :: dim_,fact

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    dim_ = fact*xgBlock%ldim*xgBlock%cols
    if (present(array_dim)) dim_=fact*array_dim

    select case (xgBlock%space)
    case ( SPACE_R,SPACE_CR )
      if ( xgBlock%cols*fact*xgBlock%Ldim < dim_ ) then
        write(std_out,*) xgBlock%cols,xgBlock%Ldim,dim_
        write(std_out,*) xgBlock%cols*xgBlock%Ldim,dim_
        ABI_ERROR("Bad reverseMapping")
      end if
      cptr = getClocR(fact*xgBlock%Ldim,xgBlock%cols,xgBlock%vecR(:,:))
      call c_f_pointer(cptr,array,(/ dim_ /))
    case ( SPACE_C )
      if ( xgBlock%cols*xgBlock%Ldim < dim_ ) then
        ABI_ERROR("Bad complex reverseMapping")
      end if
      cptr = getClocC(xgBlock%Ldim,xgBlock%cols,xgBlock%vecC(:,:))
      call c_f_pointer(cptr,array,(/ 2*dim_ /))
    end select

  end subroutine xgBlock_reverseMap_1dR
  !!***

  !!****f* m_xg/xgBlock_reverseMapC
  !!
  !! NAME
  !! xgBlock_reverseMapC

  subroutine xgBlock_reverseMapC(xgBlock,array,rows,cols)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(in) :: xgBlock
    complex(dpc), pointer, intent(inout) :: array(:,:)
    integer,optional,intent(in) :: rows
    integer,optional,intent(in) :: cols
    type(c_ptr) :: cptr

    integer :: rows_,cols_

    if (xgBlock%space/=SPACE_C) then
      ABI_ERROR('space(xgBlock) should be SPACE_C')
    end if
    rows_ = xgBlock%ldim
    cols_ = xgBlock%cols
    if (present(rows)) rows_=rows
    if (present(cols)) cols_=cols

    if ( xgBlock%cols*xgBlock%Ldim < cols_*rows_ ) then
      ABI_ERROR("Bad complex reverseMapping")
    end if
    cptr = getClocC(xgBlock%Ldim,xgBlock%cols,xgBlock%vecC(:,:))
    call c_f_pointer(cptr,array,(/ rows_, cols_ /))

  end subroutine xgBlock_reverseMapC
  !!***

  !!****f* m_xg/xgBlock_reverseMap_1dC
  !!
  !! NAME
  !! xgBlock_reverseMap_1dC

  subroutine xgBlock_reverseMap_1dC(xgBlock,array,array_dim)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t) , intent(in) :: xgBlock
    complex(dpc), pointer, intent(inout) :: array(:)
    integer,optional,intent(in) :: array_dim
    type(c_ptr) :: cptr

    integer :: dim_

    if (xgBlock%space/=SPACE_C) then
      ABI_ERROR('space(xgBlock) should be SPACE_C')
    end if
    dim_ = xgBlock%ldim*xgBlock%cols
    if (present(array_dim)) dim_=array_dim

    if ( xgBlock%cols*xgBlock%Ldim < dim_ ) then
      ABI_ERROR("Bad complex reverseMapping")
    end if
    cptr = getClocC(xgBlock%Ldim,xgBlock%cols,xgBlock%vecC(:,:))
    call c_f_pointer(cptr,array,(/ dim_ /))

  end subroutine xgBlock_reverseMap_1dC
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

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)

    real(dp), pointer                           :: array(:,:)
    integer                                     :: blockdim
    integer                                     :: spacedim
    integer                                     :: ldim
    integer(C_SIZE_T)                           :: byte_count

    ! get the array pointer underneath the xgBlock
    call xgBlock_getSize(xgBlock,spacedim,blockdim,ldim)
    call xgBlock_reverseMap(xgBlock,array,1,spacedim*blockdim)

    select case (xgBlock%space)
    case ( SPACE_R )
      byte_count = ldim*blockdim*dp
    case ( SPACE_CR )
      byte_count = 2*ldim*blockdim*dp
    case ( SPACE_C )
      byte_count = ldim*blockdim*2*dpc ! Note the factor 2, needed here!
    end select

    ! now we can call the memory prefetch
    if (present(deviceId)) then
      call gpu_data_prefetch_async(c_loc(array), byte_count , deviceId)
    else
      call gpu_data_prefetch_async(c_loc(array), byte_count)
    end if

#else
    ABI_UNUSED(deviceId)
    ABI_UNUSED_A(xgBlock)
#endif

  end subroutine xgBlock_prefetch_async
  !!***

  !LB-06/03/24: Be careful, xg_get is not used anywhere in the code, so not tested...
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
!    double precision :: tsec(2)

!    call timab(tim_get,1,tsec)

    if ( size(array,dim=1) /= 2 ) then
      ABI_ERROR("First dim must be 2")
    end if

    cols = size(array,dim=2)/rows
    if ( shift_col+cols > xg%cols ) then
      ABI_WARNING("Ignore some columns, input array to large")
    endif

    select case (xg%space)
    case (SPACE_R)
      !!$OMP TARGET UPDATE FROM(xg%vecR)
      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xg%vecR(1:rows,shift_col+col)
      end do
    case (SPACE_CR)
      if ( xg%rows /= 2*rows ) then
        ABI_ERROR("Bad number of rows")
      end if

      !!$OMP TARGET UPDATE FROM(xg%vecR)
      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = xg%vecR(1:rows,shift_col+col)
        array(2,(col-1)*rows+1:col*rows) = xg%vecR(rows+1:2*rows,shift_col+col)
      end do
    case (SPACE_C)
      !!$OMP TARGET UPDATE FROM(xg%vecC)
      do col = 1, min(cols,xg%cols-shift_col)
        array(1,(col-1)*rows+1:col*rows) = dble(xg%vecC(1:rows,shift_col+col))
        array(2,(col-1)*rows+1:col*rows) = aimag(xg%vecC(1:rows,shift_col+col))
      end do
    end select

!    call timab(tim_get,2,tsec)

  end subroutine xg_get
  !!***

  !LB-06/03/24: Be careful, xgBlock_get is not used anywhere in the code, so not tested...
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
!    double precision :: tsec(2)

!    call timab(tim_get,1,tsec)

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

!    call timab(tim_get,2,tsec)

  end subroutine xgBlock_get
  !!***

  !!****f* m_xg/xg_setBlock
  !!
  !! NAME
  !! xg_setBlock

  subroutine xg_setBlock(xg, Xgblock, rows, cols, fcol)
    use, intrinsic :: iso_c_binding
    type(xg_t), intent(inout) :: xg
    type(xgBlock_t), intent(inout) :: xgBlock
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    integer, optional, intent(in) :: fcol
    type(c_ptr) :: cptr
    integer :: fcol_,fact

    fcol_=1
    if (present(fcol)) fcol_=fcol

    if ( (fcol_+cols-1 ) > xg%cols ) then
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
    xgBlock%me_g0 = xg%me_g0
    xgBlock%gpu_option = xg%gpu_option

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2
    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xg%rows,xg%cols,xg%vecR(:,fcol_:fcol_+cols-1))
      call c_f_pointer(cptr,xgBlock%vecR,(/ fact*xgBlock%LDim,cols /))
    case(SPACE_C)
      cptr = getClocC(xg%rows,xg%cols,xg%vecC(:,fcol_:fcol_+cols-1))
      call c_f_pointer(cptr,xgBlock%vecC,(/ xgBlock%LDim,cols /))
    end select

  end subroutine xg_setBlock
  !!***

  !!****f* m_xg/xgBlock_setBlock
  !!
  !! NAME
  !! xgBlock_setBlock

  subroutine xgBlock_setBlock(xgBlockA,xgBlockB, rows, cols, fcol)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t), intent(in   ) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    integer,optional,intent(in) :: fcol
    type(c_ptr) :: cptr
    integer :: fcol_,fact

    fcol_=1
    if (present(fcol)) fcol_=fcol

    if ( (fcol_+cols-1 ) > xgblockA%cols ) then
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
    xgBlockB%me_g0 = xgBlockA%me_g0
    xgBlockB%gpu_option= xgBlockA%gpu_option

    fact = 1 ; if (xgBlockA%space==SPACE_CR) fact = 2
    select case(xgBlockA%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(xgBlockA%LDim,xgBlockA%cols,xgBlockA%vecR(:,fcol_:fcol_+cols-1))
      call c_f_pointer(cptr,xgBlockB%vecR,(/ fact*xgBlockB%LDim,cols /))
    case(SPACE_C)
      cptr = getClocC(xgBlockA%LDim,xgBlockA%cols,xgBlockA%vecC(:,fcol_:fcol_+cols-1))
      call c_f_pointer(cptr,xgBlockB%vecC,(/ xgBlockB%LDim,cols /))
    end select

  end subroutine xgBlock_setBlock
  !!***

  !!****f* m_xg/xg_free
  !!
  !! NAME
  !! xg_free

  subroutine xg_free(xg)

    type(xg_t),target, intent(inout) :: xg

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    complex(dpc), pointer :: xg__vecC(:,:)
    real(dp), pointer :: xg__vecR(:,:)
#endif

    if(xg%gpu_option==ABI_GPU_KOKKOS) then
#if defined HAVE_GPU && defined HAVE_YAKL
      if ( associated(xg%vecR) ) then
        ABI_FREE_MANAGED(xg%vecR)
      end if
      if ( associated(xg%vecC) ) then
        ABI_FREE_MANAGED(xg%vecC)
      end if
#endif

    else
      if(xg%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
        if ( associated(xg%vecR) ) then
          xg__vecR => xg%vecR
          !$OMP TARGET EXIT DATA MAP(delete:xg__vecR)
        end if
        if ( associated(xg%vecC) ) then
          xg__vecC => xg%vecC
          !$OMP TARGET EXIT DATA MAP(delete:xg__vecC)
        end if
#endif
      end if
      if ( associated(xg%vecR) ) then
        ABI_FREE(xg%vecR)
      end if
      if ( associated(xg%vecC) ) then
        ABI_FREE(xg%vecC)
      end if

    end if

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

  !!****f* m_xg/me_g0
  !!
  !! NAME
  !! me_g0

  function me_g0(xgBlock)
    type(xgBlock_t), intent(in) :: xgBlock
    integer :: me_g0
    me_g0 = xgBlock%me_g0
  end function me_g0
  !!***

  !!****f* m_xg/gpu_option
  !!
  !! NAME
  !! gpu_option

  function gpu_option(xgBlock)
    type(xgBlock_t), intent(in) :: xgBlock
    integer :: gpu_option
    gpu_option = xgBlock%gpu_option
  end function gpu_option
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
    !LB-31/05/24 : this warning is too verbose when activated... (even with correct coding!)
    !if ( rows /= xgBlock%ldim ) then
    !  ABI_WARNING("rows/ldim ! Be very careful at what you are doing")
    !end if
  end function rows
  !!***

  !!****f* m_xg/xgBlock_copy
  !!
  !! NAME
  !! xgBlock_copy

  subroutine xgBlock_copy(xgBlockA, xgBlockB, inc1, inc2)

    type(xgBlock_t),   intent(in   ) :: xgBlockA
    type(xgBlock_t),   intent(inout) :: xgBlockB
    integer, optional, intent(in   ) :: inc1
    integer, optional, intent(in   ) :: inc2

    integer :: incx
    integer :: incy

    integer :: size1
    integer :: size2
    integer :: size
    double precision :: tsec(2)
    integer :: l_gpu_option

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:)
#endif

    call timab(tim_copy,1,tsec)

    if (xgBlockA%gpu_option==xgBlockB%gpu_option) then
      l_gpu_option = xgBlockA%gpu_option
    else if ((xgBlockA%gpu_option==ABI_GPU_DISABLED.and.xgBlockB%gpu_option==ABI_GPU_LEGACY) &
       .or.  (xgBlockA%gpu_option==ABI_GPU_LEGACY  .and.xgBlockB%gpu_option==ABI_GPU_DISABLED)) then
      l_gpu_option = ABI_GPU_DISABLED
    else if (xgBlockA%gpu_option==ABI_GPU_DISABLED.and.xgBlockB%gpu_option==ABI_GPU_OPENMP) then
      l_gpu_option = ABI_GPU_DISABLED
      call xgBlock_copy_from_gpu(xgBlockB)
    else if (xgBlockB%gpu_option==ABI_GPU_DISABLED.and.xgBlockA%gpu_option==ABI_GPU_OPENMP) then
      l_gpu_option = ABI_GPU_DISABLED
      call xgBlock_copy_from_gpu(xgBlockA)
    else
      ABI_ERROR('xgA%gpu_option/=xgB%gpu_option is possible only with ABI_GPU_OPENMP, ABI_GPU_LEGACY or ABI_GPU_DISABLED')
    end if

    incx = 1; if ( present(inc1) ) incx = inc1
    incy = 1; if ( present(inc2) ) incy = inc2

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if
    !if ( xgBlockA%LDim*xgBlockA%cols/incx /= xgBlockB%LDim*xgBlockB%cols/incy ) then
    !  ABI_ERROR("Number of element different")
    !end if

    size1 = xgBlockA%LDim*xgBlockA%cols/incx ; if ( size1 * incx < xgBlockA%LDim*xgBlockA%cols ) size1 = size1+1
    size2 = xgBlockB%LDim*xgBlockB%cols/incy ; if ( size2 * incy < xgBlockB%LDim*xgBlockB%cols ) size2 = size2+1
    size = min(size1,size2)
    if (xgBlockA%space==SPACE_CR) then
      size = 2*size
    end if

    xgBlockB%me_g0 = xgBlockA%me_g0

    if (l_gpu_option==ABI_GPU_KOKKOS .or. l_gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call abi_gpu_xcopy(1, size, xgBlockA%vecR, incx, xgBlockB%vecR, incy)
      case(SPACE_C)
        call abi_gpu_xcopy(2, size, xgBlockA%vecC, incx, xgBlockB%vecC, incy)
      end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockB__vecR)
        call abi_gpu_xcopy(1, size, c_loc(xgBlockA__vecR), incx, c_loc(xgBlockB__vecR), incy)
        !$OMP END TARGET DATA
      case(SPACE_C)
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecC,xgBlockB__vecC)
        call abi_gpu_xcopy(2, size, c_loc(xgBlockA__vecC), incx, c_loc(xgBlockB__vecC), incy)
        !$OMP END TARGET DATA
      end select
#endif
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

!!****f* m_xg/xgBlock_partialcopy
!!
!! NAME
!! xgBlock_partialcopy

  subroutine xgBlock_partialcopy(xgBlock_in, xgBlock_out, shift_row, shift_col, option)

    type(xgBlock_t), intent(in   ) :: xgBlock_in
    type(xgBlock_t), intent(inout) :: xgBlock_out
    integer,intent(in) :: shift_row,shift_col
    integer,intent(in) :: option

    integer :: icol,ncols_small,ncols_big
    integer :: nrows_small,nrows_big
    integer :: shift_col_big,end_row,max_col
    double precision :: tsec(2)

    call timab(tim_partialcopy,1,tsec)

    if (xgBlock_in%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlock_in,xgBlock_out)

    if (option==SMALL2BIG) then
      ncols_big   = xgBlock_out%cols
      nrows_big   = xgBlock_out%rows
      ncols_small = xgBlock_in%cols
      nrows_small = xgBlock_in%rows
    else if (option==BIG2SMALL) then
      ncols_big   = xgBlock_in%cols
      nrows_big   = xgBlock_in%rows
      ncols_small = xgBlock_out%cols
      nrows_small = xgBlock_out%rows
    else
      ABI_ERROR('Bad option')
    end if
    if (xgBlock_in%space/=xgBlock_out%space) then
      ABI_ERROR('xgBlockA%space/=xgBlockB%space')
    end if
    if (ncols_small>ncols_big) then
      ABI_ERROR('ncols_small>ncols_big')
    end if
    if (nrows_small>nrows_big) then
      ABI_ERROR('nrows_small>nrows_big')
    end if
    if (shift_row<0) then
      ABI_ERROR('start_row<0')
    end if
    if (shift_col<0) then
      ABI_ERROR('start_col<0')
    end if

    end_row = shift_row + nrows_small
    if (end_row>nrows_big) then
      ABI_ERROR('end_row>nrows_big')
    end if

    max_col = shift_col + ncols_small
    if (max_col>ncols_big) then
      ABI_ERROR('max_col>nrows_big')
    end if

    if (option==SMALL2BIG) then ! copy small matrix into the big one
      select case(xgBlock_in%space)
      case (SPACE_R)
        do icol=1,ncols_small
          shift_col_big = shift_col + (icol-1)
          xgBlock_out%vecR(1+shift_row:end_row,1+shift_col_big) = xgBlock_in%vecR(1:nrows_small,icol)
        end do
      case (SPACE_CR)
        do icol=1,ncols_small
          shift_col_big = shift_col + (icol-1)
          xgBlock_out%vecR(1+2*shift_row:2*end_row,1+shift_col_big) = xgBlock_in%vecR(1:2*nrows_small,icol)
        end do
      case(SPACE_C)
        do icol=1,ncols_small
          shift_col_big = shift_col + (icol-1)
          xgBlock_out%vecC(1+shift_row:end_row,1+shift_col_big) = xgBlock_in%vecC(1:nrows_small,icol)
        end do
      end select
    else if (option==BIG2SMALL) then ! copy a part of the big matrix in the small one
      select case(xgBlock_in%space)
      case (SPACE_R)
        do icol=1,ncols_small
          shift_col_big = shift_col + (icol-1)
          xgBlock_out%vecR(1:nrows_small,icol) = xgBlock_in%vecR(1+shift_row:end_row,1+shift_col_big)
        end do
      case (SPACE_CR)
        if (shift_row>0) xgBlock_out%me_g0 = 0
        do icol=1,ncols_small
          shift_col_big = shift_col + (icol-1)
          xgBlock_out%vecR(1:2*nrows_small,icol) = xgBlock_in%vecR(1+2*shift_row:2*end_row,1+shift_col_big)
        end do
      case(SPACE_C)
        do icol=1,ncols_small
          shift_col_big = shift_col + (icol-1)
          xgBlock_out%vecC(1:nrows_small,icol) = xgBlock_in%vecC(1+shift_row:end_row,1+shift_col_big)
        end do
      end select
    else
      ABI_ERROR('Bad option')
    end if

    call timab(tim_partialcopy,2,tsec)

  end subroutine xgBlock_partialcopy
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

  subroutine xgBlock_gemmR(transa, transb, alpha, xgBlockA, xgBlockB, beta, xgBlockW, comm, timing)

    character,        intent(in   )           :: transa
    character,        intent(in   )           :: transb
    double precision, intent(in   )           :: alpha
    type(xgBlock_t),  intent(in   )           :: xgBlockA
    type(xgBlock_t),  intent(in   )           :: xgBlockB
    double precision, intent(in   )           :: beta
    type(xgBlock_t),  intent(inout)           :: xgBlockW
    integer,optional, intent(in)              :: comm
    logical,optional, intent(in)              :: timing

    real(dp)       :: alpha_
    complex(dpc)   :: calpha
    complex(dpc)   :: cbeta
    character(kind=1) :: transa_,transb_
    integer           :: K
    double precision  :: tsec(2)
    logical :: timing_

#if defined HAVE_OPENMP_OFFLOAD
#if !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:),xgBlockW__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:),xgBlockW__vecR(:,:)
#endif
#if !defined HAVE_MPI2_INPLACE
    complex(kind=c_double_complex), ABI_CONTIGUOUS pointer :: vecC_buf(:,:)
    real(kind=c_double), ABI_CONTIGUOUS pointer :: vecR_buf(:,:)
#endif
#endif

    timing_ = .true.
    if (present(timing)) then
      timing_ = timing
    end if
    if (timing_) call timab(tim_gemm_blas,1,tsec)

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)

    if (transa /= 'n' .and. transa /= 't') then
      ABI_ERROR("transa should be 'n' or 't'")
    end if
    if (transb /= 'n' .and. transb /= 't') then
      ABI_ERROR("transb should be 'n' or 't'")
    end if

    if ( transa == 'n' ) then
      K = xgBlockA%cols
      if ( xgBlockA%rows /= xgBlockW%rows ) then
        ABI_ERROR("rows(A)/=rows(W)")
      end if
    else
      K = xgBlockA%rows
      if ( xgBlockA%cols /= xgBlockW%rows ) then
        ABI_ERROR("cols(A)/=rows(W)")
      end if
    end if
    if ( transb == 'n' ) then
      if ( xgBlockB%cols /= xgBlockW%cols ) then
        ABI_ERROR("cols(B)/=cols(W)")
      end if
    else
      if ( xgBlockB%rows /= xgBlockW%cols ) then
        ABI_ERROR("rows(B)/=cols(W)")
      end if
    end if

    calpha = dcmplx(alpha,0.d0)
    cbeta  = dcmplx(beta, 0.d0)

    if ( xgBlockA%space == xgBlockB%space ) then

      if ( xgBlockA%space /= SPACE_CR .and. xgBlockW%space /= xgBlockA%space ) then
        ABI_ERROR("Not same space for A and W")
      end if
      if ( xgBlockA%space == SPACE_CR ) then
        if ( transa/= 't' .or. transb /='n' ) then
          ABI_ERROR("if space(A)==SPACE_CR, transa shoulbe 't' and transb shoulb 'n'")
        end if
        if ( xgBlockW%space /= SPACE_R) then
          ABI_ERROR("space(W) should be SPACE_R")
        end if
        if (xgBlockA%me_g0<0) then
          ABI_ERROR("xgBlockA me_g0 is not initialized")
        end if
        if (xgBlockB%me_g0<0) then
          ABI_ERROR("xgBlockB me_g0 is not initialized")
        end if
        if (xgBlockA%me_g0/=xgBlockB%me_g0) then
          ABI_ERROR("xgBlockA and xgBlockB should have same me_g0")
        end if
      end if
      ! CALL GEMM
      select case(xgBlockA%space)
      case (SPACE_R)
        if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          calpha, &
          xgBlockA%vecR, xgBlockA%LDim, &
          xgBlockB%vecR, xgBlockB%LDim, &
          cbeta, &
          xgBlockW%vecR, xgBlockW%LDim)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        xgBlockW__vecR => xgBlockW%vecR
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockB__vecR,xgBlockW__vecR)
        call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          calpha, &
          c_loc(xgBlockA__vecR), xgBlockA%LDim, &
          c_loc(xgBlockB__vecR), xgBlockB%LDim, &
          cbeta, &
          c_loc(xgBlockW__vecR), xgBlockW%LDim)
        !$OMP END TARGET DATA
#endif
      else
        call dgemm(transa, transb, xgBlockW%rows, xgBlockW%cols, K, &
          alpha, &
          xgBlockA%vecR, xgBlockA%LDim, &
          xgBlockB%vecR, xgBlockB%LDim, &
          beta, &
          xgBlockW%vecR, xgBlockW%LDim)
      end if

      case (SPACE_CR)
        calpha = dcmplx(2*alpha,0.d0)
        alpha_ = 2.0d0 * alpha
        if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
          call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, 2*K, &
            calpha, &
            xgBlockA%vecR, 2*xgBlockA%LDim, &
            xgBlockB%vecR, 2*xgBlockB%LDim, &
            cbeta, &
            xgBlockW%vecR, xgBlockW%LDim)
          if (xgBlockA%me_g0 == 1) then
            calpha = dcmplx(-2*alpha,0.d0)
            call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, 2, &
              calpha, &
              xgBlockA%vecR, 2*xgBlockA%LDim, &
              xgBlockB%vecR, 2*xgBlockB%LDim, &
              cone, &
              xgBlockW%vecR, xgBlockW%LDim)
            calpha = dcmplx(alpha,0.d0)
            call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, 1, &
              calpha, &
              xgBlockA%vecR, 2*xgBlockA%LDim, &
              xgBlockB%vecR, 2*xgBlockB%LDim, &
              cone, &
              xgBlockW%vecR, xgBlockW%LDim)
          end if
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
          xgBlockA__vecR => xgBlockA%vecR
          xgBlockB__vecR => xgBlockB%vecR
          xgBlockW__vecR => xgBlockW%vecR
          !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockB__vecR,xgBlockW__vecR)
          call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, 2*K, &
            calpha, &
            c_loc(xgBlockA__vecR), 2*xgBlockA%LDim, &
            c_loc(xgBlockB__vecR), 2*xgBlockB%LDim, &
            cbeta, &
            c_loc(xgBlockW__vecR), xgBlockW%LDim)
          if (xgBlockA%me_g0 == 1) then
            calpha = dcmplx(-2*alpha,0.d0)
            call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, 2, &
              calpha, &
              c_loc(xgBlockA__vecR), 2*xgBlockA%LDim, &
              c_loc(xgBlockB__vecR), 2*xgBlockB%LDim, &
              cone, &
              c_loc(xgBlockW__vecR), xgBlockW%LDim)
            calpha = dcmplx(alpha,0.d0)
            call abi_gpu_xgemm(1, transa, transb, xgBlockW%rows, xgBlockW%cols, 1, &
              calpha, &
              c_loc(xgBlockA__vecR), 2*xgBlockA%LDim, &
              c_loc(xgBlockB__vecR), 2*xgBlockB%LDim, &
              cone, &
              c_loc(xgBlockW__vecR), xgBlockW%LDim)
          end if
          !$OMP END TARGET DATA
#endif
        else
          call dgemm(transa, transb, xgBlockW%rows, xgBlockW%cols, 2*K, &
            alpha_, &
            xgBlockA%vecR, 2*xgBlockA%LDim, &
            xgBlockB%vecR, 2*xgBlockB%LDim, &
            beta, &
            xgBlockW%vecR, xgBlockW%LDim)
          if (xgBlockA%me_g0 == 1) then
            alpha_ = - 2*alpha
            call dgemm(transa,transb,xgBlockW%rows, xgBlockW%cols, 2, &
              alpha_, &
              xgBlockA%vecR, 2*xgBlockA%LDim, &
              xgBlockB%vecR, 2*xgBlockB%LDim, &
              one, &
              xgBlockW%vecR,xgBlockW%LDim)
            alpha_ = alpha
            call dgemm(transa,transb,xgBlockW%rows, xgBlockW%cols, 1, &
              alpha_, &
              xgBlockA%vecR, 2*xgBlockA%LDim, &
              xgBlockB%vecR, 2*xgBlockB%LDim, &
              one, &
              xgBlockW%vecR,xgBlockW%LDim)
          end if
        end if

      case(SPACE_C)

        transa_=transa
        if (transa=='t') transa_ = 'c'
        transb_=transb
        if (transb=='t') transb_ = 'c'

        if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
          call abi_gpu_xgemm(2, transa_, transb_, xgBlockW%rows, xgBlockW%cols, K, &
            calpha, &
            xgBlockA%vecC, xgBlockA%LDim, &
            xgBlockB%vecC, xgBlockB%LDim, &
            cbeta, &
            xgBlockW%vecC, xgBlockW%LDim)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
          xgBlockA__vecC => xgBlockA%vecC
          xgBlockB__vecC => xgBlockB%vecC
          xgBlockW__vecC => xgBlockW%vecC
          !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecC,xgBlockB__vecC,xgBlockW__vecC)
          call abi_gpu_xgemm(2, transa_, transb_, xgBlockW%rows, xgBlockW%cols, K, &
            calpha, &
            c_loc(xgBlockA__vecC), xgBlockA%LDim, &
            c_loc(xgBlockB__vecC), xgBlockB%LDim, &
            cbeta, &
            c_loc(xgBlockW__vecC), xgBlockW%LDim)
          !$OMP END TARGET DATA
#endif
        else
          call zgemm(transa_, transb_, xgBlockW%rows, xgBlockW%cols, K, &
            calpha, &
            xgBlockA%vecC, xgBlockA%LDim, &
            xgBlockB%vecC, xgBlockB%LDim, &
            cbeta, &
            xgBlockW%vecC, xgBlockW%LDim)
        end if

      end select

    else ! not same space for A and B

      if (xgBlockA%space==SPACE_CR.and.xgBlockB%space==SPACE_R) then
        if (transa/='n'.or.transb/='n') then
          ABI_ERROR('Not implemented')
        end if
        if (xgBlockW%space/=SPACE_CR) then
          ABI_ERROR('space(W) should be SPACE_CR')
        end if
        if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
          call abi_gpu_xgemm(1, transa, transb, 2*xgBlockW%rows, xgBlockW%cols, K, &
            calpha, &
            xgBlockA%vecR, 2*xgBlockA%LDim, &
            xgBlockB%vecR, xgBlockB%LDim, &
            cbeta, &
            xgBlockW%vecR, 2*xgBlockW%LDim)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
          xgBlockA__vecR => xgBlockA%vecR
          xgBlockB__vecR => xgBlockB%vecR
          xgBlockW__vecR => xgBlockW%vecR
          !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockB__vecR,xgBlockW__vecR)
          call abi_gpu_xgemm(1, transa, transb, 2*xgBlockW%rows, xgBlockW%cols, K, &
            calpha, &
            c_loc(xgBlockA__vecR), 2*xgBlockA%LDim, &
            c_loc(xgBlockB__vecR), xgBlockB%LDim, &
            cbeta, &
            c_loc(xgBlockW__vecR), 2*xgBlockW%LDim)
          !$OMP END TARGET DATA
#endif
        else ! CPU
          call dgemm(transa,transb,2*xgBlockW%rows, xgBlockW%cols,K, &
            alpha,xgBlockA%vecR, 2*xgBlockA%LDim, &
            xgBlockB%vecR, xgBlockB%LDim, beta,xgBlockW%vecR,2*xgBlockW%LDim)
        end if
      else
        ABI_ERROR('Not implemented')
      end if

    end if

    if (timing_) call timab(tim_gemm_blas,2,tsec)
    ! END CALL GEMM

    ! MPI SUM
    if ( present(comm) ) then
      if (timing_) call timab(tim_gemm_mpi,1,tsec)
      call xgBlock_mpi_sum(xgBlockW,comm=comm)
      if (timing_) call timab(tim_gemm_mpi,2,tsec)
    end if

  end subroutine xgBlock_gemmR
  !!***

  !!****f* m_xg/xgBlock_gemmC
  !!
  !! NAME
  !! xgBlock_gemmC

  subroutine xgBlock_gemmC(transa, transb, alpha, xgBlockA, xgBlockB, beta, xgBlockW, comm, timing)

    character,       intent(in   ) :: transa
    character,       intent(in   ) :: transb
    complex(kind=8), intent(in   ) :: alpha
    type(xgBlock_t), intent(in   ) :: xgBlockA
    type(xgBlock_t), intent(in   ) :: xgBlockB
    complex(kind=8), intent(in   ) :: beta
    type(xgBlock_t), intent(inout) :: xgBlockW
    integer,optional,intent(in)    :: comm
    logical,optional,intent(in)    :: timing

    integer          :: K
    double precision :: tsec(2)
    character(kind=1) :: transa_,transb_
    logical :: timing_

    timing_ = .true.
    if (present(timing)) then
      timing_ = timing
    end if
    if (timing_) call timab(tim_gemm_blas,1,tsec)

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockB%space /= xgBlockW%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockA%space /= SPACE_C ) then
      ABI_ERROR("Not correct space")
    end if

    if (transa /= 'n' .and. transa /= 't') then
      ABI_ERROR("transa should be 'n' or 't'")
    end if
    if (transb /= 'n' .and. transb /= 't') then
      ABI_ERROR("transb should be 'n' or 't'")
    end if

    if ( transa == 'n' ) then
      K = xgBlockA%cols
    else
      K = xgBlockA%rows
    end if

    transa_=transa
    if (transa=='t') transa_ = 'c'
    transb_=transb
    if (transb=='t') transb_ = 'c'

    ! CALL GEMM
    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
      call abi_gpu_xgemm(2, transa_, transb_, xgBlockW%rows, xgBlockW%cols, K, &
        alpha, &
        xgBlockA%vecC, xgBlockA%LDim, &
        xgBlockB%vecC, xgBlockB%LDim, &
        beta, &
        xgBlockW%vecC, xgBlockW%LDim)
    else
      call zgemm(transa_, transb_, xgBlockW%rows, xgBlockW%cols, K, &
        alpha, &
        xgBlockA%vecC, xgBlockA%LDim, &
        xgBlockB%vecC, xgBlockB%LDim, &
        beta, &
        xgBlockW%vecC, xgBlockW%LDim)
    end if
    ! END CALL GEMM
    if (timing_) call timab(tim_gemm_blas,2,tsec)

    ! MPI SUM
    if ( present(comm) ) then
      if (timing_) call timab(tim_gemm_mpi,1,tsec)
      call xgBlock_mpi_sum(xgBlockW,comm=comm)
      if (timing_) call timab(tim_gemm_mpi,2,tsec)
    end if

  end subroutine xgBlock_gemmC
  !!***

!!****f* m_xg/xgBlock_trmmR
!!
!! NAME
!! xgBlock_trmmR

  subroutine xgBlock_trmmR(side, uplo, transa, diag, alpha, xgBlockA, xgBlockB)

    character, intent(in) :: transa,side,uplo,diag
    double precision, intent(in) :: alpha
    type(xgBlock_t), intent(in) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB

    complex(kind=8) :: calpha

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockB%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if

    select case(xgBlockA%space)
    case (SPACE_R)
      call dtrmm(side,uplo,transa,diag,transa,xgBlockB%rows,xgBlockB%cols,&
        alpha,xgBlockA%vecR, xgBlockA%LDim, &
        xgBlockB%vecR, xgBlockB%LDim)
    case (SPACE_CR)
      ABI_ERROR("Not implemented")
      !call dtrmm(side,uplo,transa,diag,transa,xgBlockB%rows,xgBlockB%cols,&
      !  alpha,xgBlockA%vecR, xgBlockA%LDim, &
      !  xgBlockB%vecR, xgBlockB%LDim)
    case(SPACE_C)
      calpha = dcmplx(alpha,0.d0)
      call ztrmm(side,uplo,transa,diag,transa,xgBlockB%rows,xgBlockB%cols,&
        calpha,xgBlockA%vecC, xgBlockA%LDim, &
        xgBlockB%vecC, xgBlockB%LDim)
    end select

  end subroutine xgBlock_trmmR
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

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)
#endif

    call timab(tim_potrf,1,tsec)

    if ( xgBlock%rows /= xgBlock%cols ) then
      ABI_ERROR("Matrix should be a square matrixx")
    endif

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS .or. xgBlock%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      select case(xgBlock%space)
      case (SPACE_R)
        call abi_gpu_xpotrf(1,uplo,xgBlock%rows,xgBlock%vecR,xgBlock%LDim,info)
      case (SPACE_C)
        call abi_gpu_xpotrf(2,uplo,xgBlock%rows,xgBlock%vecC,xgBlock%LDim,info)
      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
      end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      select case(xgBlock%space)
      case (SPACE_R)
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock__vecR)
        call abi_gpu_xpotrf(1,uplo,xgBlock%rows,c_loc(xgBlock__vecR),xgBlock%LDim,info)
        !$OMP END TARGET DATA
      case (SPACE_C)
        xgBlock__vecC => xgBlock%vecC
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock__vecC)
        call abi_gpu_xpotrf(2,uplo,xgBlock%rows,c_loc(xgBlock__vecC),xgBlock%LDim,info)
        !$OMP END TARGET DATA
      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
      end select
#endif
      if(xgBlock%gpu_option==ABI_GPU_KOKKOS) call gpu_device_synchronize()

    else
      select case(xgBlock%space)
      case (SPACE_R)
        call dpotrf(uplo,xgBlock%rows,xgBlock%vecR,xgBlock%LDim,info)
      case (SPACE_C)
        call zpotrf(uplo,xgBlock%rows,xgBlock%vecC,xgBlock%LDim,info)
      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
      end select
    end if

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

  subroutine xgBlock_heevd(jobz, uplo, xgBlockA, xgBlockW, info)

    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer         , intent(  out) :: info
    double precision :: tsec(2)

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockW__vecR(:,:)
#endif

    call timab(tim_heevd,1,tsec)

    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)

    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      select case(xgBlockA%space)

      case (SPACE_R)
        call abi_gpu_xheevd(1,jobz,uplo,xgBlockA%cols, &
            xgBlockA%vecR,xgBlockA%LDim, &
            xgBlockW%vecR,info)

      case (SPACE_C)
        call abi_gpu_xheevd(2,jobz,uplo,xgBlockA%cols, &
            xgBlockA%vecC,xgBlockA%LDim, &
            xgBlockW%vecR,info)

      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')

      end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      select case(xgBlockA%space)

      case (SPACE_R)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockW__vecR => xgBlockW%vecR
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockW__vecR)
        call abi_gpu_xheevd(1,jobz,uplo,xgBlockA%cols, &
            c_loc(xgBlockA__vecR),xgBlockA%LDim, &
            c_loc(xgBlockW__vecR),info)
        !$OMP END TARGET DATA
        !call checkResize(iwork,liwork,5*xgBlockA%rows+3)
        !call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+6*xgBlockA%rows+1)
        !!$OMP TARGET UPDATE FROM(xgBlockA__vecR)
        !call dsyevd(jobz,uplo,xgBlockA%cols, &
        !  xgBlockA%vecR,xgBlockA%LDim, &
        !  xgBlockW%vecR, rwork, lrwork, &
        !  iwork, liwork,info)
        !!$OMP TARGET UPDATE TO(xgBlockA__vecR)
        !!$OMP TARGET UPDATE TO(xgBlockW__vecR)
        !if ( rwork(1) > lrwork ) then
        !  !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
        !  call checkResize(rwork,lrwork,int(rwork(1)))
        !end if
        !if ( iwork(1) > liwork ) then
        !  !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
        !  call checkResize(iwork,liwork,int(iwork(1)))
        !end if

      case (SPACE_C)
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockW__vecR => xgBlockW%vecR
        !!$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecC,xgBlockW__vecR)
        !call abi_gpu_xheevd(2,jobz,uplo,xgBlockA%cols, &
        !    c_loc(xgBlockA__vecC),xgBlockA%LDim, &
        !    c_loc(xgBlockW__vecR),info)
        !!$OMP END TARGET DATA
        call checkResize(iwork,liwork,5*xgBlockA%rows+3)
        call checkResize(cwork,lcwork,xgBlockA%rows*xgBlockA%rows+2*xgBlockA%rows)
        call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+5*xgBlockA%rows+1)
        !$OMP TARGET UPDATE FROM(xgBlockA__vecC)
        call zheevd(jobz,uplo,xgBlockA%cols, &
          xgBlockA%vecC,xgBlockA%LDim, &
        xgBlockW%vecR, &
        cwork, lcwork, rwork, lrwork, iwork, liwork, info)
        !$OMP TARGET UPDATE TO(xgBlockA__vecC)
        !$OMP TARGET UPDATE TO(xgBlockW__vecR)
        if ( int(cwork(1)) > lcwork ) then
          !write(std_out,*) "Allocate work from", int(lcwork), "to", int(cwork(1))
          call checkResize(cwork,lcwork,int(cwork(1)))
        end if
        if ( rwork(1) > lrwork ) then
          !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
          call checkResize(rwork,lrwork,int(rwork(1)))
        end if
        if ( iwork(1) > liwork ) then
          !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
          call checkResize(iwork,liwork,int(iwork(1)))
        end if

      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')

      end select
#endif

      if(xgBlockA%gpu_option==ABI_GPU_KOKKOS) call gpu_device_synchronize()

    else

      call checkResize(iwork,liwork,5*xgBlockA%rows+3)

      select case(xgBlockA%space)

      case (SPACE_R)
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

      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

  subroutine xgBlock_hegvd(itype, jobz, uplo, xgBlockA, xgBlockB, xgBlockW, info)

    integer         , intent(in   ) :: itype
    character       , intent(in   ) :: jobz
    character       , intent(in   ) :: uplo
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockB
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer         , intent(  out) :: info

    double precision :: tsec(2)

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:),xgBlockW__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:),xgBlockW__vecR(:,:)
#endif

    call timab(tim_hegvd,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Not same space")
    end if
    if ( xgBlockW%space /= SPACE_R ) then
      ABI_ERROR("Block3 must be real")
    end if

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then

      select case(xgBlockA%space)

      case (SPACE_R)

#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        call abi_gpu_xhegvd(1, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             xgBlockA%vecR, xgBlockA%ldim, &
          &             xgBlockB%vecR, xgBlockB%ldim, &
          &             xgBlockW%vecR, &
          &             info)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        xgBlockW__vecR => xgBlockW%vecR
        !!$OMP TARGET UPDATE FROM(xgBlockA__vecR)
        !!$OMP TARGET UPDATE FROM(xgBlockB__vecR)
        !call checkResize(iwork,liwork,5*xgBlockA%rows+3)
        !call checkResize(rwork,lrwork,2*xgBlockA%rows*xgBlockA%rows+6*xgBlockA%rows+1)
        !call dsygvd(itype, jobz, uplo, xgBlockA%rows, xgBlockA%vecR, xgBlockA%ldim, &
        !  xgBlockB%vecR, xgBlockB%ldim, xgBlockW%vecR, rwork, lrwork, iwork, liwork, info)
        !!$OMP TARGET UPDATE TO(xgBlockA__vecR)
        !!$OMP TARGET UPDATE TO(xgBlockB__vecR)
        !!$OMP TARGET UPDATE TO(xgBlockW__vecR)
        !if ( rwork(1) > lrwork ) then
        !  !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
        !  call checkResize(rwork,lrwork,int(rwork(1)))
        !end if
        !if ( iwork(1) > liwork ) then
        !  !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
        !  call checkResize(iwork,liwork,int(iwork(1)))
        !end if
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockB__vecR,xgBlockW__vecR)
        call abi_gpu_xhegvd(1, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             c_loc(xgBlockA__vecR), xgBlockA%ldim, &
          &             c_loc(xgBlockB__vecR), xgBlockB%ldim, &
          &             c_loc(xgBlockW__vecR), &
          &             info)
        !$OMP END TARGET DATA
#endif

      case (SPACE_C)

        !call xgBlock_prefetch_async(xgBlockA, 0)
        !call xgBlock_prefetch_async(xgBlockB, 0)
        !call xgBlock_prefetch_async(xgBlockW, 0)
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        call abi_gpu_xhegvd(2, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             xgBlockA%vecC, xgBlockA%ldim, &
          &             xgBlockB%vecC, xgBlockB%ldim, &
          &             xgBlockW%vecR, &
          &             info)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        xgBlockW__vecR => xgBlockW%vecR
        !call checkResize(iwork,liwork,5*xgBlockA%rows+3)
        !call checkResize(cwork,lcwork,xgBlockA%rows*xgBlockA%rows+2*xgBlockA%rows)
        !call checkResize(rwork,lrwork,2*(xgBlockA%rows*xgBlockA%rows)+5*xgBlockA%rows+1)

        !!$OMP TARGET UPDATE FROM(xgBlockA__vecC)
        !!$OMP TARGET UPDATE FROM(xgBlockB__vecC)
        !call zhegvd(itype, jobz, uplo, xgBlockA%rows, xgBlockA%vecC, xgBlockA%ldim,&
        !  xgBlockB%vecC, xgBlockB%ldim, xgBlockW%vecR, cwork, lcwork, &
        !  rwork, lrwork, iwork, liwork, info)
        !!$OMP TARGET UPDATE TO(xgBlockA__vecC)
        !!$OMP TARGET UPDATE TO(xgBlockB__vecC)
        !!$OMP TARGET UPDATE TO(xgBlockW__vecR)

        !if ( int(cwork(1)) > lcwork ) then
        !  !write(std_out,*) "Allocate work from", lcwork, "to", int(cwork(1))
        !  call checkResize(cwork,lcwork,int(cwork(1)))
        !end if
        !if ( rwork(1) > lrwork ) then
        !  !write(std_out,*) "Allocate work from", lrwork, "to", int(rwork(1))
        !  call checkResize(rwork,lrwork,int(rwork(1)))
        !end if
        !if ( iwork(1) > liwork ) then
        !  !write(std_out,*) "Allocate work from", liwork, "to", int(iwork(1))
        !  call checkResize(iwork,liwork,int(iwork(1)))
        !end if
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecC,xgBlockB__vecC,xgBlockW__vecR)
        call abi_gpu_xhegvd(2, itype, jobz, uplo, &
          &             xgBlockA%rows, &
          &             c_loc(xgBlockA__vecC), xgBlockA%ldim, &
          &             c_loc(xgBlockB__vecC), xgBlockB%ldim, &
          &             c_loc(xgBlockW__vecR), &
          &             info)
        !$OMP END TARGET DATA
#endif

      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')

      end select

      if(xgBlockA%gpu_option==ABI_GPU_KOKKOS) call gpu_device_synchronize()

    else

      call checkResize(iwork,liwork,5*xgBlockA%rows+3)

      select case(xgBlockA%space)

      case (SPACE_R)

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

      case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

    case (SPACE_R)
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

    case (SPACE_CR)
      ABI_ERROR('Not implemented for SPACE_CR')

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

  subroutine xgBlock_trsmR(side,uplo,transa,diag,alpha,xgBlockA,xgBlockB)

    character       , intent(in   ) :: side
    character       , intent(in   ) :: uplo
    character       , intent(in   ) :: transa
    character       , intent(in   ) :: diag
    double precision, intent(in   ) :: alpha
    type(xgBlock_t) , intent(inout) :: xgBlockA
    type(xgBlock_t) , intent(inout) :: xgBlockB
    complex(kind=8) :: calpha
    double precision :: tsec(2)
    integer :: fact

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:)
#endif

    call timab(tim_trsm,1,tsec)
    if ( xgBlockB%space/=SPACE_CR ) then
      if ( xgBlockA%space /= xgBlockB%space ) then
        ABI_ERROR("Not same space")
      end if
    else
      if ( xgBlockA%space /= SPACE_R ) then
        ABI_ERROR("If space(B)=SPACE_CR, space(A) should be space(R)")
      end if
    end if

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)

    calpha = dcmplx(alpha,0.d0)

    fact = 1 ; if (xgBlockB%space==SPACE_CR) fact = 2

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call abi_gpu_xtrsm(1,side,uplo,transa,diag,fact*xgBlockB%rows,xgBlockB%cols, &
          calpha,xgBlockA%vecR,xgBlockA%LDim,xgBlockB%vecR,fact*xgBlockB%LDim)
      case (SPACE_C)
        call abi_gpu_xtrsm(2,side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
          calpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)
      end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecR,xgBlockB__vecR)
        call abi_gpu_xtrsm(1,side,uplo,transa,diag,fact*xgBlockB%rows,xgBlockB%cols, &
          calpha,c_loc(xgBlockA__vecR),xgBlockA%LDim,c_loc(xgBlockB__vecR),fact*xgBlockB%LDim)
        !$OMP END TARGET DATA
      case (SPACE_C)
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecC,xgBlockB__vecC)
        call abi_gpu_xtrsm(2,side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
          calpha,c_loc(xgBlockA__vecC),xgBlockA%LDim,c_loc(xgBlockB__vecC),xgBlockB%LDim)
        !$OMP END TARGET DATA
      end select
#endif

    else
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call dtrsm(side,uplo,transa,diag,fact*xgBlockB%rows,xgBlockB%cols, &
          alpha,xgBlockA%vecR,xgBlockA%LDim,xgBlockB%vecR,fact*xgBlockB%LDim)
      case (SPACE_C)
        calpha = dcmplx(alpha,0.d0)
        call ztrsm(side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
          calpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)
      end select

    end if

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

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:)
#endif

    call timab(tim_trsm,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space .or. xgBlockA%space /= SPACE_C) then
      ABI_ERROR("Space should be SPACE_C for xgBlockA and xgBlockB")
    end if

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)

    if (xgBlockA%gpu_option == ABI_GPU_KOKKOS .or. xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      call abi_gpu_xtrsm(2,side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
        alpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlockA__vecC,xgBlockB__vecC)
      xgBlockA__vecC => xgBlockA%vecC
      xgBlockB__vecC => xgBlockB%vecC
      call abi_gpu_xtrsm(2,side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
        alpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)
      !$OMP END TARGET DATA
#endif

    else
      call ztrsm(side,uplo,transa,diag,xgBlockB%rows,xgBlockB%cols, &
        alpha,xgBlockA%vecC,xgBlockA%LDim,xgBlockB%vecC,xgBlockB%LDim)
    end if

    call timab(tim_trsm,2,tsec)

  end subroutine xgBlock_trsmC
  !!***

  !!****f* m_xg/xgBlock_ymax
  !!
  !! NAME
  !! xgBlock_ymax
  subroutine xgBlock_ymax(xgBlockA, da, shift, nblocks)

    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(in   ) :: da
    integer, intent(in) :: shift,nblocks

    integer :: iblock,ncols,irow,nrows,fact
    double precision :: tsec(2)

    call timab(tim_ymax,1,tsec)

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,da)

    nrows = xgBlockA%rows
    ncols = xgBlockA%cols

    if ( da%rows /= nblocks*ncols ) then
      ABI_ERROR("rows(da)/=nblocks*ncols")
    end if
    if ( shift<0 ) then
      ABI_ERROR("shift<0")
    end if
    if ( shift+ncols > da%rows ) then
      ABI_ERROR("shift+xgBlockA%cols > da%rows")
    end if

    fact = 1 ; if (xgBlockA%space==SPACE_CR) fact = 2

    if (space(da)==SPACE_R) then
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do collapse(2) shared(da,xgBlockA) private(irow,iblock)
        do iblock = 1, ncols
          do irow = 1, fact*nrows
            xgBlockA%vecR(irow,iblock) = - da%vecR(iblock+shift,1) &
             & * xgBlockA%vecR(irow,iblock)
          end do
        end do
        !$omp end parallel do
      case (SPACE_C)
        !$omp parallel do collapse(2) shared(da,xgBlockA) private(irow,iblock)
        do iblock = 1, ncols
          do irow = 1, nrows
            xgBlockA%vecC(irow,iblock) = - da%vecR(iblock+shift,1) &
             & * xgBlockA%vecC(irow,iblock)
          end do
        end do
        !$omp end parallel do
      end select
    else if (space(da)==SPACE_C) then
      if (xgBlockA%space/=SPACE_C) then
        ABI_ERROR('If space(da)=SPACE_C, space(xgBlockA) has to be SPACE_C')
      end if
      !$omp parallel do collapse(2) shared(da,xgBlockA) private(irow,iblock)
      do iblock = 1, ncols
        do irow = 1, nrows
          xgBlockA%vecC(irow,iblock) = - da%vecC(iblock+shift,1) &
           & * xgBlockA%vecC(irow,iblock)
        end do
      end do
      !$omp end parallel do
    else
      ABI_ERROR('Only SPACE_R or SPACE_C (for da) are implemented.')
    end if

    call timab(tim_ymax,2,tsec)

  end subroutine xgBlock_ymax
  !!***

  !!****f* m_xg/xgBlock_colwiseCymax
  !!
  !! NAME
  !! xgBlock_colwiseCymax

  subroutine xgBlock_colwiseCymax(xgBlockA, da, xgBlockB, xgBlockW)

    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(in   ) :: da
    type(xgBlock_t), intent(in   ) :: xgBlockB
    type(xgBlock_t), intent(in   ) :: xgBlockW

    integer :: iblock,fact,rows,cols,jblock
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:),xgBlockW__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:),xgBlockW__vecR(:,:),da__vecR(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_colw_cymax,1,tsec)

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

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)
    call xgBlock_check_gpu_option(xgBlockA,da)

    fact = 1 ; if (xgBlockA%space==SPACE_CR) fact = 2

    rows = fact*xgBlockA%rows; cols = xgBlockA%cols

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS) then

#if defined HAVE_GPU && defined HAVE_KOKKOS

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        call compute_colwiseCymax_scalar(c_loc(xgBlockA%vecR), c_loc(da%vecR), c_loc(xgBlockB%vecR), &
          &                              c_loc(xgBlockW%vecR), fact*xgBlockA%rows, xgBlockA%cols, fact*xgBlockA%ldim)
      case (SPACE_C)
        call compute_colwiseCymax_cplx  (c_loc(xgBlockA%vecC), c_loc(da%vecR), c_loc(xgBlockB%vecC), &
          &                              c_loc(xgBlockW%vecC), xgBlockA%rows, xgBlockA%cols, xgBlockA%ldim)
      end select

#endif

    else if (xgBlockA%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        xgBlockW__vecR => xgBlockW%vecR
        da__vecR => da%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !$OMP& MAP(to:xgBlockA__vecR,xgBlockB__vecR,xgBlockW__vecR,da__vecR)
        do iblock = 1, cols
          do jblock = 1, rows
            xgBlockA__vecR(jblock,iblock) = - da__vecR(iblock,1) * xgBlockB__vecR(jblock,iblock) &
                + xgBlockW__vecR(jblock,iblock)
          end do
        end do
      case (SPACE_C)
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        xgBlockW__vecC => xgBlockW%vecC
        da__vecR => da%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !$OMP& MAP(to:xgBlockA__vecC,xgBlockB__vecC,xgBlockW__vecC,da__vecR)
        do iblock = 1, cols
          do jblock = 1, rows
            xgBlockA__vecC(jblock,iblock) = - da__vecR(iblock,1) * xgBlockB__vecC(jblock,iblock) &
                + xgBlockW__vecC(jblock,iblock)
          end do
        end do
      end select

#endif

    else

      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do collapse(2) shared(da,xgBlockB,xgBlockW,xgBlockA) private(iblock,jblock)
        do iblock = 1, cols
          do jblock = 1, rows
            xgBlockA%vecR(jblock,iblock) = - da%vecR(iblock,1) * xgBlockB%vecR(jblock,iblock) + xgBlockW%vecR(jblock,iblock)
          end do
        end do
        !$omp end parallel do
      case (SPACE_C)
        !$omp parallel do collapse(2) shared(da,xgBlockB,xgBlockW,xgBlockA) private(iblock,jblock)
        do iblock = 1, cols
          do jblock = 1, rows
            xgBlockA%vecC(jblock,iblock) = - da%vecR(iblock,1) * xgBlockB%vecC(jblock,iblock) + xgBlockW%vecC(jblock,iblock)
          end do
        end do
        !$omp end parallel do
      end select

    end if

    call timab(tim_colw_cymax,2,tsec)

  end subroutine xgBlock_colwiseCymax
  !!***

!!****f* m_xg/xgBlock_yxmax
!!
!! NAME
!! xgBlock_yxmax

  subroutine xgBlock_yxmax(xgBlockA, da, xgBlockB)

    type(xgBlock_t), intent(inout) :: xgBlockA
    type(xgBlock_t), intent(in   ) :: da
    type(xgBlock_t), intent(in   ) :: xgBlockB

    integer :: iblock,irow,cols,rows,fact
    double precision :: tsec(2)

    call timab(tim_yxmax,1,tsec)

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,da)
    call xgBlock_check_gpu_option(xgBlockA,da)

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Must be same space for ymax")
    end if
    if ( xgBlockA%rows /= xgBlockB%rows ) then
      ABI_ERROR("Must have same rows for ymax")
    end if
    if ( xgBlockA%cols /= xgBlockB%cols ) then
      ABI_ERROR("Must have same cols for ymax")
    end if
    if ( da%rows /= xgBlockA%cols ) then
      ABI_ERROR("Must have same cols for ymax")
    end if

    cols = xgBlockA%cols
    rows = xgBlockA%rows
    fact = 1 ; if (xgBlockA%space==SPACE_CR) fact = 2

    select case(xgBlockA%space)
    case (SPACE_R,SPACE_CR)
      !$omp parallel do collapse(2) shared(da,xgBlockB,xgBlockA) private(iblock,irow)
      do iblock = 1, cols
        do irow = 1, fact*rows
          xgBlockA%vecR(irow,iblock) = xgBlockA%vecR(irow,iblock) - da%vecR(iblock,1) * xgBlockB%vecR(irow,iblock)
        end do
      end do
      !$omp end parallel do
    case (SPACE_C)
      !$omp parallel do collapse(2) shared(da,xgBlockB,xgBlockA) private(iblock,irow)
      do iblock = 1, cols
        do irow = 1, rows
          xgBlockA%vecC(irow,iblock) = xgBlockA%vecC(irow,iblock) - da%vecR(iblock,1) * xgBlockB%vecC(irow,iblock)
        end do
      end do
      !$omp end parallel do
    end select

    call timab(tim_yxmax,2,tsec)

  end subroutine xgBlock_yxmax
!!***

  !!****f* m_xg/xgBlock_apply_diag
  !!
  !! NAME
  !! xgBlock_apply_diag

  subroutine xgBlock_apply_diag(X, diag, nspinor, Y)

    type(xgBlock_t) , intent(inout) :: X
    type(xgBlock_t) , intent(in)    :: diag
    integer,          intent(in)    :: nspinor
    type(xgBlock_t) , optional, intent(inout) :: Y

    type(xgBlock_t) :: X_spinor, Y_spinor
    real(dp)    , pointer :: array(:)
    complex(dpc), pointer :: arrayc(:)
    double precision :: tsec(2)

    call timab(tim_apply_diag,1,tsec)

    if (X%rows/=nspinor*diag%rows) then
      ABI_ERROR('xgBlock%rows/=nspinor*xgBlock_diag%rows')
    end if
    if (diag%cols/=1) then
      ABI_ERROR('diag should have one column')
    end if
    if (diag%space/=SPACE_R.and.diag%space/=SPACE_C) then
      ABI_ERROR('space(diag) should be SPACE_C or SPACE_R')
    end if
    if (X%space==SPACE_R) then
      if (diag%space/=SPACE_R) then
        ABI_ERROR('If space(X)==SPACE_R, space(diag) should be SPACE_R')
      end if
    end if

    if (present(Y)) then
      call xgBlock_check(Y,X)
      call xgBlock_copy(X,Y)
      call xgBlock_reshape_spinor(Y,Y_spinor,nspinor,ROWS2COLS)
    else
      call xgBlock_reshape_spinor(X,X_spinor,nspinor,ROWS2COLS)
      Y_spinor = X_spinor
    end if

    if (space(diag)==SPACE_R) then
      call xgBlock_reverseMap_1dR(diag,array,array_dim=diag%rows)
      call xgBlock_colwiseMulR(Y_spinor,array)
    else if (space(diag)==SPACE_C) then
      call xgBlock_reverseMap_1dC(diag,arrayc,array_dim=diag%rows)
      call xgBlock_colwiseMulC(Y_spinor,arrayc)
    end if

    call timab(tim_apply_diag,2,tsec)

  end subroutine xgBlock_apply_diag
  !!***

!****f* m_xg/xgBlock_add_diag
!
! NAME
! xgBlock_add_diag

  subroutine xgBlock_add_diag(X, diag, nspinor, Y)

    type(xgBlock_t) , intent(in)    :: X
    type(xgBlock_t) , intent(in)    :: diag
    integer,          intent(in)    :: nspinor
    type(xgBlock_t) , intent(inout) :: Y

    integer :: iblock,irow,rows,cols
    type(xgBlock_t) :: X_spinor, Y_spinor
    double precision :: tsec(2)

    call timab(tim_add_diag,1,tsec)

    if (X%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(X,diag)
    call xgBlock_check_gpu_option(X,Y)

    if (X%rows/=nspinor*diag%rows) then
      ABI_ERROR('xgBlock%rows/=nspinor*xgBlock_diag%rows')
    end if
    if (diag%cols/=1) then
      ABI_ERROR('xgBlock_diag should have one column')
    end if
    if (diag%space==SPACE_CR) then
      ABI_ERROR('space(diag) should be SPACE_C or SPACE_R')
    end if
    if (X%space==SPACE_R) then
      if (diag%space/=SPACE_R) then
        ABI_ERROR('If space(X)==SPACE_R, space(diag) should be SPACE_R')
      end if
    end if

    call xgBlock_check(Y,X)

    call xgBlock_reshape_spinor(X,X_spinor,nspinor,ROWS2COLS)
    call xgBlock_reshape_spinor(Y,Y_spinor,nspinor,ROWS2COLS)

    rows = X_spinor%rows
    cols = X_spinor%cols

    select case(X%space)
    case (SPACE_R)
      !$omp parallel do collapse(2) shared(X_spinor,Y_spinor,diag) private(iblock,irow)
      do iblock = 1, cols
        do irow=1,rows
          Y_spinor%vecR(irow,iblock) = Y_spinor%vecR(irow,iblock) &
            & + X_spinor%vecR(irow,iblock) * diag%vecR(irow,1)
        end do
      end do
    case (SPACE_CR)
      if (diag%space==SPACE_R) then
        !$omp parallel do collapse(2) shared(X_spinor,Y_spinor,diag) private(iblock,irow)
        do iblock = 1, cols
          do irow=1,rows
            Y_spinor%vecR(2*irow-1,iblock) = Y_spinor%vecR(2*irow-1,iblock) &
              & + X_spinor%vecR(2*irow-1,iblock) * diag%vecR(irow,1)
            Y_spinor%vecR(2*irow  ,iblock) = Y_spinor%vecR(2*irow  ,iblock) &
              & + X_spinor%vecR(2*irow  ,iblock) * diag%vecR(irow,1)
          end do
        end do
      else
        ABI_ERROR('Not implemented')
      end if
    case (SPACE_C)
      if (diag%space==SPACE_C) then
        !$omp parallel do collapse(2) shared(X_spinor,Y_spinor,diag) private(iblock,irow)
        do iblock = 1, cols
          do irow=1,rows
            Y_spinor%vecC(irow,iblock) = Y_spinor%vecC(irow,iblock) &
              & + X_spinor%vecC(irow,iblock) * diag%vecC(irow,1)
          end do
        end do
      else if (diag%space==SPACE_R) then
        !$omp parallel do collapse(2) shared(X_spinor,Y_spinor,diag) private(iblock,irow)
        do iblock = 1, cols
          do irow=1,rows
            Y_spinor%vecC(irow,iblock) = Y_spinor%vecC(irow,iblock) &
              & + X_spinor%vecC(irow,iblock) * diag%vecR(irow,1)
          end do
        end do
      else
        ABI_ERROR('Not implemented')
      end if
    end select

    call timab(tim_add_diag,2,tsec)

  end subroutine xgBlock_add_diag
!!***

  !!****f* m_xg/xgBlock_mpi_sum
  !!
  !! NAME
  !! xgBlock_mpi_sum

  subroutine xgBlock_mpi_sum(xgBlock,comm)

    type(xgBlock_t) , intent(inout) :: xgBlock
    integer,intent(in),optional :: comm

    integer :: ierr,comm_

    if (.not.present(comm)) then
      comm_ = xgBlock%spacedim_comm
    else
      comm_ = comm
    end if

    if ( xmpi_comm_size(comm_) > 1) then
      if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then
        ! CPU waits for GPU to finish before doing MPI communications
        call gpu_device_synchronize()
      end if

      select case(xgBlock%space)

        case (SPACE_R,SPACE_CR)
          call xmpi_sum(xgBlock%vecR,comm_,ierr,use_omp_map=(xgBlock%gpu_option==ABI_GPU_OPENMP))
        case (SPACE_C)
          call xmpi_sum(xgBlock%vecC,comm_,ierr,use_omp_map=(xgBlock%gpu_option==ABI_GPU_OPENMP))
      end select
    end if ! xmpi_comm_size>1

  end subroutine xgBlock_mpi_sum
  !!***

  !!****f* m_xg/xgBlock_mpi_send
  !!
  !! NAME
  !! xgBlock_mpi_send

  subroutine xgBlock_mpi_send(xgBlock,dest,tag,comm)

    type(xgBlock_t) , intent(in) :: xgBlock
    integer,intent(in) :: dest,tag
    integer,intent(in),optional :: comm

    integer :: ierr,comm_
    real(dp), pointer :: vec(:,:)

    if (xgBlock%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if

    if (.not.present(comm)) then
      comm_ = xgBlock%spacedim_comm
    else
      comm_ = comm
    end if

    call xgBlock_reverseMap(xgBlock,vec)
    call xmpi_send(vec,dest,tag,comm_,ierr)

  end subroutine xgBlock_mpi_send
  !!***

  !!****f* m_xg/xgBlock_mpi_isend
  !!
  !! NAME
  !! xgBlock_mpi_isend

  subroutine xgBlock_mpi_isend(xgBlock,dest,tag,request,comm)

    type(xgBlock_t) , intent(in) :: xgBlock
    integer,intent(in) :: dest,tag
    integer,intent(inout) :: request
    integer,intent(in),optional :: comm

    integer :: ierr,comm_
    real(dp), pointer :: vec(:,:)

    if (xgBlock%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if

    if (.not.present(comm)) then
      comm_ = xgBlock%spacedim_comm
    else
      comm_ = comm
    end if

    call xgBlock_reverseMap(xgBlock,vec)
    call xmpi_isend(vec,dest,tag,comm_,request,ierr)

  end subroutine xgBlock_mpi_isend
  !!***

  !!****f* m_xg/xgBlock_mpi_recv
  !!
  !! NAME
  !! xgBlock_mpi_recv

  subroutine xgBlock_mpi_recv(xgBlock,source,tag,comm)

    type(xgBlock_t) , intent(inout) :: xgBlock
    integer,intent(in) :: source,tag
    integer,intent(in),optional :: comm

    integer :: ierr,comm_
    real(dp), pointer :: vec(:,:)

    if (xgBlock%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if

    if (.not.present(comm)) then
      comm_ = xgBlock%spacedim_comm
    else
      comm_ = comm
    end if

    call xgBlock_reverseMap(xgBlock,vec)
    call xmpi_recv(vec,source,tag,comm_,ierr)

  end subroutine xgBlock_mpi_recv
  !!***

  !!****f* m_xg/xgBlock_gemm_mpi_cyclic_permutation
  !!
  !! NAME
  !! xgBlock_gemm_mpi_cyclic_permutation

  subroutine xgBlock_gemm_mpi_cyclic_permutation(xgBlockA,xgBlockB,xgBlockW,me_comm,blocksize,comm)

    type(xgBlock_t) , intent(in)    :: xgBlockA,xgBlockB
    type(xgBlock_t) , intent(inout) :: xgBlockW
    integer,intent(in) :: me_comm
    integer,intent(in),optional :: blocksize,comm

    logical :: multiblock
    integer :: ierr,blocksize_,comm_,source,dest,tag,request
    integer :: iblock_left,iblock_right,iblock_mpi,nblocks_mpi,nblocks_left,nblocks_right
    integer :: shift_col,shift_row,shift_col_mpi,shift_row_mpi
    double precision :: tsec(2)
    type(xg_t) :: xg_mpi_work
    type(xg_t) :: subB,subB_mpi

    call timab(tim_gemmcyclic,1,tsec)

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)

    if (.not.present(comm)) then
      comm_ = xgBlockA%spacedim_comm
    else
      comm_ = comm
    end if

    nblocks_mpi = xmpi_comm_size(comm_)

    if (xgBlockA%rows/=xgBlockW%rows) then
      ABI_ERROR('rows(xgBlockA)/=rows(xgBlockW)')
    end if
    if (xgBlockB%rows/=nblocks_mpi*xgBlockA%cols) then
      ABI_ERROR('rows(xgBlockB)/=nblocks_mpi*cols(xgBlockA)')
    end if
    if (xgBlockB%cols/=nblocks_mpi*xgBlockW%cols) then
      ABI_ERROR('cols(xgBlockB)/=nblocks_mpi*cols(xgBlockW)')
    end if

    blocksize_ = xgBlockA%cols
    if (present(blocksize)) then
      if (mod(xgBlockA%cols,blocksize)/=0) then
        ABI_ERROR('invalid blocksize')
      end if
      if (mod(xgBlockB%cols/nblocks_mpi,blocksize)/=0) then
        ABI_ERROR('invalid blocksize')
      end if
      blocksize_ = blocksize
    end if

    if (nblocks_mpi==1) then

      ! If only one mpi process, use timing from gemm routine
      call timab(tim_gemmcyclic,2,tsec)
      call xgBlock_gemm('n','n',1.0d0,xgBlockA,xgBlockB,1.d0,xgBlockW)
      call timab(tim_gemmcyclic,1,tsec)

    else

      nblocks_left  = xgBlockA%cols / blocksize_
      nblocks_right = xgBlockB%cols / (blocksize_*nblocks_mpi)
      multiblock = .false.
      if (nblocks_left>1.or.nblocks_right>1) then
        multiblock = .true.
      end if

      call xg_init(xg_mpi_work,xgBlockA%space,xgBlockA%rows,xgBlockA%cols,xmpi_comm_null)
      call xg_init(subB_mpi,xgBlockB%space,xgBlockB%rows/nblocks_mpi,xgBlockB%cols/nblocks_mpi,xmpi_comm_null)
      if (multiblock) then
        call xg_init(subB,xgBlockB%space,blocksize_,blocksize_,xmpi_comm_null)
      end if

      do iblock_mpi=1,nblocks_mpi

        shift_row_mpi = mod((iblock_mpi-1)+me_comm,nblocks_mpi) * blocksize_
        shift_col_mpi = me_comm * blocksize_
        if (.not.multiblock) then
          call xgBlock_partialcopy(xgBlockB,subB_mpi%self,shift_row_mpi,shift_col_mpi,BIG2SMALL)
        else
          do iblock_right=1,nblocks_right
            do iblock_left=1,nblocks_left
              shift_row = shift_row_mpi + (iblock_left-1) * blocksize_*nblocks_mpi
              shift_col = shift_col_mpi + (iblock_right-1) * blocksize_*nblocks_mpi
              call xgBlock_partialcopy(xgBlockB,subB%self,shift_row,shift_col,BIG2SMALL)
              shift_row = (iblock_left-1) * blocksize_
              shift_col = (iblock_right-1) * blocksize_
              call xgBlock_partialcopy(subB%self,subB_mpi%self,shift_row,shift_col,SMALL2BIG)
            end do
          end do
        end if

        if (iblock_mpi==1) then
          call xgBlock_gemm('n','n',1.0d0,xgBlockA,subB_mpi%self,1.d0,xgBlockW,timing=.false.)
        else
          tag = iblock_mpi
          dest = mod(me_comm-(iblock_mpi-1),nblocks_mpi)
          if (dest<0) dest=dest+nblocks_mpi
          call xgBlock_mpi_isend(xgBlockA,dest,tag,request,comm=comm_)
          source = mod(me_comm+(iblock_mpi-1),nblocks_mpi)
          call xgBlock_mpi_recv(xg_mpi_work%self,source,tag,comm=comm_)
          call xgBlock_gemm('n','n',1.0d0,xg_mpi_work%self,subB_mpi%self,1.d0,xgBlockW,timing=.false.)
        end if

        if (iblock_mpi>1) call xmpi_wait(request,ierr)

      end do

      call xg_free(xg_mpi_work)
      call xg_free(subB_mpi)
      if (multiblock) then
        call xg_free(subB)
      end if

    end if

    call timab(tim_gemmcyclic,2,tsec)

  end subroutine xgBlock_gemm_mpi_cyclic_permutation
  !!***

  !!****f* m_xg/xgBlock_colwiseMulR
  !!
  !! NAME
  !! xgBlock_colwiseMulR

  subroutine xgBlock_colwiseMulR(xgBlock, vec)

    type(xgBlock_t) , intent(inout)           :: xgBlock
    double precision, intent(in   ), target   :: vec(:)

    integer :: rows
    integer :: iblock,irow

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    integer :: cols
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_colw_mul,1,tsec)

    rows = size(vec,dim=1)

    if (xgBlock%rows/=rows) then
      ABI_ERROR('dim(vec)/=xgBlock%rows')
    end if

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlock%space)
      case (SPACE_R)
        call compute_colwiseMul_scalar_scalar(c_loc(xgBlock%vecR), c_loc(vec), &
          &                                   0, xgBlock%rows, xgBlock%cols, &
          &                                   xgBlock%ldim, rows)
      case (SPACE_CR)
         ABI_ERROR('Not implemented')
      case (SPACE_C)
        call compute_colwiseMul_cplx_scalar(c_loc(xgBlock%vecC), c_loc(vec), &
          &                                 0, xgBlock%rows, xgBlock%cols, &
          &                                 xgBlock%ldim, rows)
      end select

#endif

    else if (xgBlock%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD

      !$OMP TARGET ENTER DATA MAP(to:vec)
      cols=xgBlock%cols
      select case(xgBlock%space)
      case (SPACE_R)
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlock__vecR) MAP(to:vec) PRIVATE(iblock,irow)
        do iblock = 1, cols
          do irow = 1, rows
            xgBlock__vecR(irow,iblock) = xgBlock__vecR(irow,iblock) * vec(irow)
          end do
        end do
      case (SPACE_CR)
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlock__vecR) MAP(to:vec) PRIVATE(iblock,irow)
        do iblock = 1, cols
          do irow = 1, rows
            xgBlock__vecR(2*irow-1,iblock) = xgBlock__vecR(2*irow-1,iblock) * vec(irow)
            xgBlock__vecR(2*irow  ,iblock) = xgBlock__vecR(2*irow  ,iblock) * vec(irow)
          end do
        end do
      case (SPACE_C)
        xgBlock__vecC => xgBlock%vecC
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlock__vecC) MAP(to:vec) PRIVATE(iblock,irow)
        do iblock = 1, cols
          do irow = 1, rows
            xgBlock__vecC(irow,iblock) = xgBlock__vecC(irow,iblock) * vec(irow)
          end do
        end do
      end select
      !$OMP TARGET EXIT DATA MAP(delete:vec)
#endif

    else

      select case(xgBlock%space)
      case (SPACE_R)
        !$omp parallel do collapse(2) shared(xgBlock,vec) private(iblock,irow)
        do iblock = 1, xgBlock%cols
          do irow = 1, rows
            xgBlock%vecR(irow,iblock) = xgBlock%vecR(irow,iblock) * vec(irow)
          end do
        end do
      case (SPACE_CR)
        !$omp parallel do collapse(2) shared(xgBlock,vec) private(iblock,irow)
        do iblock = 1, xgBlock%cols
          do irow = 1, rows
            xgBlock%vecR(2*irow-1,iblock) = xgBlock%vecR(2*irow-1,iblock) * vec(irow)
            xgBlock%vecR(2*irow  ,iblock) = xgBlock%vecR(2*irow  ,iblock) * vec(irow)
          end do
        end do
      case (SPACE_C)
        !$omp parallel do collapse(2) shared(xgBlock,vec) private(iblock,irow)
        do iblock = 1, xgBlock%cols
          do irow = 1, rows
            xgBlock%vecC(irow,iblock) = xgBlock%vecC(irow,iblock) * vec(irow)
          end do
        end do
      end select

    end if

    call timab(tim_colw_mul,2,tsec)

  end subroutine xgBlock_colwiseMulR
  !!***

  !!****f* m_xg/xgBlock_colwiseMulC
  !!
  !! NAME
  !! xgBlock_colwiseMulC

  subroutine xgBlock_colwiseMulC(xgBlock, vec)

    type(xgBlock_t), intent(inout)           :: xgBlock
    complex(kind=8), intent(in   ), target   :: vec(:)

    integer :: rows
    integer :: iblock,irow

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    integer :: cols
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_colw_mul,1,tsec)

    ABI_UNUSED((/irow/)) ! Use in OpenMP GPU
    rows = size(vec,dim=1)

    if (xgBlock%rows/=rows) then
      ABI_ERROR('dim(vec)/=xgBlock%rows')
    end if

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        ABI_ERROR("Error colwiseMulC")
      case (SPACE_C)
        call compute_colwiseMul_cplx_cplx(c_loc(xgBlock%vecC), c_loc(vec), &
          &                               0, xgBlock%rows, xgBlock%cols, &
          &                               xgBlock%ldim, rows)
      end select

#endif

    else if (xgBlock%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        ABI_ERROR("Error colwiseMulC")
      case (SPACE_C)
        cols = xgBlock%cols
        xgBlock__vecC => xgBlock%vecC
        !$OMP TARGET ENTER DATA MAP(to:vec)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlock__vecC) MAP(to:vec) PRIVATE(iblock,irow)
        do iblock = 1, cols
          do irow = 1, rows
            xgBlock__vecC(irow,iblock) = xgBlock__vecC(irow,iblock) * vec(irow)
          end do
        end do
      end select
      !$OMP TARGET EXIT DATA MAP(delete:vec)

#endif

    else

      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        ABI_ERROR("Error colwiseMulC")
      case (SPACE_C)
        !$omp parallel do collapse(2) shared(xgBlock,vec) private(iblock,irow)
        do iblock = 1, xgBlock%cols
          do irow=1, rows
            xgBlock%vecC(irow,iblock) = xgBlock%vecC(irow,iblock) * vec(irow)
          end do
        end do
      end select

    end if

    call timab(tim_colw_mul,2,tsec)

  end subroutine xgBlock_colwiseMulC
  !!***

  !!****f* m_xg/xgBlock_saxpyR
  !!
  !! NAME
  !! xgBlock_saxpyR

  subroutine xgBlock_saxpyR(xgBlock1, da, xgBlock2)

    type(xgBlock_t),  intent(inout) :: xgBlock1
    double precision, intent(in   ) :: da
    type(xgBlock_t),  intent(in   ) :: xgBlock2

    integer :: fact
    complex(dpc) :: da_cplx
#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock1__vecC(:,:),xgBlock2__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock1__vecR(:,:),xgBlock2__vecR(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_saxpy,1,tsec)

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

    call xgBlock_check_gpu_option(xgBlock1,xgBlock2)

    fact = 1 ; if (xgBlock1%space==SPACE_CR) fact = 2

    if (xgBlock1%gpu_option==ABI_GPU_KOKKOS .or. xgBlock1%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        call abi_gpu_xaxpy(1, xgBlock1%cols*fact*xgBlock1%LDim, da_cplx, xgBlock2%vecR,1,xgBlock1%vecR,1)
      case (SPACE_C)
        call abi_gpu_xaxpy(2, xgBlock1%cols*xgBlock1%LDim, da_cplx, xgBlock2%vecC,1,xgBlock1%vecC,1)
      end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        xgBlock1__vecR => xgBlock1%vecR
        xgBlock2__vecR => xgBlock2%vecR
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock1__vecR,xgBlock2__vecR)
        call abi_gpu_xaxpy(1, xgBlock1%cols*fact*xgBlock1%LDim, da_cplx, c_loc(xgBlock2__vecR),1,c_loc(xgBlock1__vecR),1)
        !$OMP END TARGET DATA
      case (SPACE_C)
        xgBlock1__vecC => xgBlock1%vecC
        xgBlock2__vecC => xgBlock2%vecC
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock1__vecC,xgBlock2__vecC)
        call abi_gpu_xaxpy(2, xgBlock1%cols*xgBlock1%LDim, da_cplx, c_loc(xgBlock2__vecC),1,c_loc(xgBlock1__vecC),1)
        !$OMP END TARGET DATA
      end select
#endif

    else
      select case(xgBlock1%space)
      case (SPACE_R,SPACE_CR)
        call daxpy(xgBlock1%cols*fact*xgBlock1%LDim,da,xgBlock2%vecR,1,xgBlock1%vecR,1)
      case (SPACE_C)
        call zaxpy(xgBlock1%cols*xgBlock1%LDim,dcmplx(da,0.d0),xgBlock2%vecC,1,xgBlock1%vecC,1)
      end select

    end if

    call timab(tim_saxpy,2,tsec)

  end subroutine xgBlock_saxpyR
  !!***

  !!****f* m_xg/xgBlock_saxpyC
  !!
  !! NAME
  !! xgBlock_saxpyC

  subroutine xgBlock_saxpyC(xgBlock1, da, xgBlock2)

    type(xgBlock_t), intent(inout) :: xgBlock1
    double complex,  intent(in   ) :: da
    type(xgBlock_t), intent(in   ) :: xgBlock2

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock1__vecC(:,:),xgBlock2__vecC(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_saxpy,1,tsec)

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

    call xgBlock_check_gpu_option(xgBlock1,xgBlock2)

    if (xgBlock1%gpu_option==ABI_GPU_KOKKOS .or. xgBlock2%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      call abi_gpu_xaxpy(2, xgBlock1%cols*xgBlock1%LDim, da, xgBlock2%vecC, 1, xgBlock1%vecC, 1)
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      xgBlock1__vecC => xgBlock1%vecC
      xgBlock2__vecC => xgBlock2%vecC
      !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock1__vecC,xgBlock2__vecC)
      call abi_gpu_xaxpy(2, xgBlock1%cols*xgBlock1%LDim, da, c_loc(xgBlock2__vecC),1,c_loc(xgBlock1__vecC),1)
      !$OMP END TARGET DATA
#endif

    else
      call zaxpy(xgBlock1%cols*xgBlock1%LDim, da, xgBlock2%vecC, 1, xgBlock1%vecC, 1)
    end if

    call timab(tim_saxpy,2,tsec)

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
    integer :: fact
    integer :: rows,cols

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:), xgBlockB__vecR(:,:)
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:), xgBlockB__vecC(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_add,1,tsec)

    if ( xgBlockA%space /= xgBlockB%space ) then
      ABI_ERROR("Must be same space for add")
    end if
    if ( xgBlockA%rows /= xgBlockB%rows ) then
      ABI_ERROR("Must have same LDim for add")
    end if
    if ( xgBlockA%cols /= xgBlockB%cols ) then
      ABI_ERROR("Must have same cols for add")
    end if

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)

    fact = 1 ; if (xgBlockA%space==SPACE_CR) fact = 2
    rows=fact*xgBlockB%rows
    cols=xgBlockB%cols

    if (xgBlockA%gpu_option==ABI_GPU_OPENMP) then
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlockA__vecR,xgBlockB__vecR)
        do col = 1, cols
          do row = 1, rows
            xgBlockA__vecR(row,col) = xgBlockA__vecR(row,col) + xgBlockB__vecR(row,col)
          end do
        end do
        !call daxpy(xgBlockA%cols*xgBlockA%LDim,1.d0,xgBlockB%vecR,1,xgBlockA%vecR1)
      case (SPACE_C)
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlockA__vecC,xgBlockB__vecC)
        do col = 1, cols
          do row = 1, rows
            xgBlockA__vecC(row,col) = xgBlockA__vecC(row,col) + xgBlockB__vecC(row,col)
          end do
        end do
        !call zaxpy(xgBlockA%cols*xgBlockA%LDim,1.d0,xgBlockB%vecR,1,xgBlockA%vecR1)
      end select
#endif
    else
      select case(xgBlockA%space)
      case (SPACE_R,SPACE_CR)
        !$omp parallel do collapse(2) shared(xgBlockA,xgBlockB) private(col,row)
        do col = 1, cols
          do row = 1, rows
            xgBlockA%vecR(row,col) = xgBlockA%vecR(row,col) + xgBlockB%vecR(row,col)
          end do
        end do
        !call daxpy(xgBlockA%cols*xgBlockA%LDim,1.d0,xgBlockB%vecR,1,xgBlockA%vecR1)
      case (SPACE_C)
        !$omp parallel do collapse(2) shared(xgBlockA,xgBlockB) private(col,row)
        do col = 1, cols
          do row = 1, rows
            xgBlockA%vecC(row,col) = xgBlockA%vecC(row,col) + xgBlockB%vecC(row,col)
          end do
        end do
        !call zaxpy(xgBlockA%cols*xgBlockA%LDim,1.d0,xgBlockB%vecR,1,xgBlockA%vecR1)
      end select
    end if

    call timab(tim_add,2,tsec)

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
    case (SPACE_R)
      xgBlock%vecR(:,:) = cshift(xgBlock%vecR(:,:),nshift,dim=shiftdim) ! Bottom 2*blockdim lines are now at the top
    case (SPACE_C)
      xgBlock%vecC(:,:) = cshift(xgBlock%vecC(:,:),nshift,dim=shiftdim) ! Bottom 2*blockdim lines are now at the top
    case (SPACE_CR)
      ABI_ERROR('Not implemented')
    end select
    call timab(tim_cshift,2,tsec)

  end subroutine xgBlock_cshift
  !!***

  !!****f* m_xg/xgBlock_colwiseNorm2
  !!
  !! NAME
  !! xgBlock_colwiseNorm2

  subroutine xgBlock_colwiseNorm2(xgBlock, dot, max_val, max_elt, min_val, min_elt, comm_loc)

    type(xgBlock_t) , intent(in   ) :: xgBlock
    type(xgBlock_t) , intent(inout) :: dot
    double precision, intent(  out), optional :: max_val
    integer         , intent(  out), optional :: max_elt
    double precision, intent(  out), optional :: min_val
    integer         , intent(  out), optional :: min_elt
    integer         , intent(in   ), optional :: comm_loc

    integer :: icol, ierr, fact, comm_
    double precision,external :: ddot
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    integer :: cols,rows
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:),dot__vecR(:,:)
#endif

#if (defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD) || defined FC_CRAY
    integer :: ii
    double precision :: tmp
#endif
    double precision :: tsec(2)

    call timab(tim_colw_norm2,1,tsec)

    if ( dot%space /= SPACE_R ) then
      ABI_ERROR("space(dot) should be SPACE_R")
    end if
    if ( dot%cols /= 1 ) then
      ABI_ERROR("cols(dot) should be 1")
    end if
    if ( dot%rows /= xgBlock%cols ) then
      ABI_ERROR("rows(dot) should be cols(xgBlock)")
    end if
    comm_=comm(xgBlock)
    if (present(comm_loc)) then
      comm_ = comm_loc
    end if

    if (xgBlock%space==SPACE_CR.and.xgBlock%me_g0<0) then
      ABI_ERROR("xgBlock me_g0 is not initialized")
    end if
    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlock%space)
      case(SPACE_R,SPACE_CR)
        if (xgBlock%space==SPACE_CR) then
          ABI_ERROR('Not implemented for GPU with KOKKOS')
        end if
        call computeBatchedDotProduct_scalar(c_loc(xgBlock%vecR), c_loc(xgBlock%vecR), &
          & c_loc(dot%vecR), fact*xgBlock%rows, xgBlock%cols, fact*xgBlock%ldim)

      case(SPACE_C)
        call computeBatchedDotProduct_cplx_scalar(c_loc(xgBlock%vecC), c_loc(xgBlock%vecC), &
          & c_loc(dot%vecR), xgBlock%rows, xgBlock%cols, xgBlock%ldim)

      end select
      call xmpi_sum(dot%vecR,comm_,ierr)

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
      ! we shouldn't be here, it means gpu_option was wrongly set to 1 in
      ! input parameter file
#endif

    else if (xgBlock%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD

      cols=xgBlock%cols
      rows=fact*xgBlock%rows
      dot__vecR => dot%vecR
      select case(xgBlock%space)
      case(SPACE_R,SPACE_CR)
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecR,xgBlock__vecR) PRIVATE(icol,tmp)
        do icol = 1, cols
          tmp=0
          !$OMP PARALLEL DO REDUCTION(+:tmp) PRIVATE(ii)
          do ii = 1,rows
            tmp = tmp + fact*xgBlock__vecR(ii,icol)*xgBlock__vecR(ii,icol)
          end do
          dot__vecR(icol,1)=tmp
        end do
        if (xgBlock%me_g0==1) then
          !$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecR,xgBlock__vecR) PRIVATE(icol)
          do icol = 1, cols
            dot__vecR(icol,1) = dot__vecR(icol,1) - xgBlock__vecR(1,icol)*xgBlock__vecR(1,icol) &
&                                                 - xgBlock__vecR(2,icol)*xgBlock__vecR(2,icol)
          end do
        end if
      case(SPACE_C)
        xgBlock__vecC => xgBlock%vecC
        !$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecR,xgBlock__vecC) PRIVATE(icol,tmp)
        do icol = 1, cols
          tmp=0
          !$OMP PARALLEL DO REDUCTION(+:tmp) PRIVATE(ii)
          do ii = 1,rows
            tmp = tmp + dconjg(xgBlock__vecC(ii,icol))*xgBlock__vecC(ii,icol)
          end do
          dot__vecR(icol,1)=tmp
        end do
      end select
      !FIXME This should happen inplace ideally
      !$OMP TARGET UPDATE FROM(dot__vecR)
      call xmpi_sum(dot%vecR,comm_,icol)
      !$OMP TARGET UPDATE TO(dot__vecR)

      ! do reductions
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

#endif

    else

      select case(xgBlock%space)
      case(SPACE_R,SPACE_CR)
        !$omp parallel do shared(dot,xgBlock)
        do icol = 1, xgBlock%cols
          dot%vecR(icol,1) = fact*ddot(fact*xgBlock%rows,xgBlock%vecR(:,icol),1,xgBlock%vecR(:,icol),1)
        end do
        !$omp end parallel do
        if (xgBlock%me_g0==1) then
          !$omp parallel do shared(dot,xgBlock)
          do icol = 1, xgBlock%cols
            dot%vecR(icol,1) = dot%vecR(icol,1) - ddot(2,xgBlock%vecR(:,icol),1,xgBlock%vecR(:,icol),1)
          end do
          !$omp end parallel do
        end if
      case(SPACE_C)
#if defined(FC_CRAY)
!FIXME zdotc call goes wrong with NVHPC (NVHPC 22.11, MKL 22.3) or CRAY
        !$omp parallel do private(ii,tmp)
        do icol = 1, xgBlock%cols
          tmp=0
          do ii = 1, xgBlock%rows
            tmp = tmp + dconjg(xgBlock%vecC(ii,icol))*xgBlock%vecC(ii,icol)
          end do
          dot%vecR(icol,1)=tmp
        end do
        !$omp end parallel do
#else
        !$omp parallel do shared(dot,xgBlock)
        do icol = 1, xgBlock%cols
          ! Instead of calling a complex function to get only the real part of the
          ! result
          !dot%vecR(icol,1) = dble(zdotc(xgBlock%rows,xgBlock%vecC(:,icol),1,xgBlock%vecC(:,icol),1))
          ! Directely call a real function which gives what we want.
          dot%vecR(icol,1) = ddot(2*xgBlock%rows,xgBlock%vecC(:,icol),1,xgBlock%vecC(:,icol),1)
        end do
        !$omp end parallel do
#endif
      end select
      call xmpi_sum(dot%vecR,comm_,ierr)

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

    end if ! if gpu_option==ABI_GPU_KOKKOS

    call timab(tim_colw_norm2,2,tsec)

  end subroutine xgBlock_colwiseNorm2
  !!***

  !!****f* m_xg/xgBlock_colwiseDotProduct
  !!
  !! NAME
  !! xgBlock_colwiseDotProduct

  subroutine xgBlock_colwiseDotProduct(xgBlockA,xgBlockB,dot,max_val,max_elt,min_val,min_elt,comm_loc)

    type(xgBlock_t)  , intent(in   ) :: xgBlockA
    type(xgBlock_t)  , intent(in   ) :: xgBlockB
    type(xgBlock_t)  , intent(inout) :: dot
    integer, intent(in), optional :: comm_loc
    double precision , intent(  out), optional :: max_val
    integer          , intent(  out), optional :: max_elt
    double precision , intent(  out), optional :: min_val
    integer          , intent(  out), optional :: min_elt
    integer :: icol,fact,comm_
    double precision,external :: ddot
    double complex,external :: zdotc !conjugated dot product

#if (defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD) || defined(FC_NVHPC) || defined(FC_CRAY)
    integer :: rows,cols,ii,me_g0
    double precision :: tmp
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:),dot__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:),dot__vecR(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_colw_dot,1,tsec)

    call xgBlock_check(xgBlockA,xgBlockB)

    comm_=comm(xgBlockA)
    if (present(comm_loc)) then
      comm_ = comm_loc
    else
      if (comm(xgBlockA)/=comm(xgBlockB)) then
        ABI_ERROR('xgBlockA and xgBlockB should have the same comm')
      end if
    end if
    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,dot)

    if (xgBlockA%space/=SPACE_CR) then
      if (dot%space/=xgBlockA%space) then
        ABI_ERROR('xgBlockA and dot should have the same space')
      end if
    else ! space(A) = space(B) = SPACE_CR
      if (xgBlockA%me_g0<0) then
        ABI_ERROR("xgBlockA me_g0 is not initialized")
      end if
      if (xgBlockB%me_g0<0) then
        ABI_ERROR("xgBlockB me_g0 is not initialized")
      end if
      if (xgBlockA%me_g0/=xgBlockB%me_g0) then
        ABI_ERROR('xgBlockA and xgBlockB should have the same me_g0')
      end if
      if (dot%space/=SPACE_R) then
        ABI_ERROR('if space(A)=SPACE_CR, space(dot) should be SPACE_R')
      end if
    end if
    if ( dot%cols /= 1 ) then
      ABI_ERROR("cols(dot) should be 1")
    end if
    if ( dot%rows /= xgBlockA%cols ) then
      ABI_ERROR("rows(dot) should be cols(xgBlockA)")
    end if

    fact = 1 ; if (xgBlockA%space==SPACE_CR) fact = 2

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        if (xgBlockA%space==SPACE_CR) then
          ABI_ERROR('Not implemented for GPU with KOKKOS')
        end if
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
      ! we shouldn't be here, it means gpu_option was wrongly set to 1 in
      ! input parameter file
#endif

    else if (xgBlockA%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      rows = xgBlockA%rows; cols = xgBlockA%cols
      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        dot__vecR => dot%vecR
        me_g0 = xgBlockA%me_g0
#if defined FC_NVHPC
        !$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecR,xgBlockA__vecR,xgBlockB__vecR) PRIVATE(icol,tmp)
        do icol = 1, cols
          tmp=0
          !$OMP PARALLEL DO REDUCTION(+:tmp) PRIVATE(ii)
          do ii = 1, fact*rows
            tmp = tmp + fact*xgBlockA__vecR(ii,icol)*xgBlockB__vecR(ii,icol)
          end do
          dot__vecR(icol,1)=tmp
        end do
        if (me_g0==1) then
          !$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecR,xgBlockA__vecR,xgBlockB__vecR) PRIVATE(icol)
          do icol = 1, cols
            dot__vecR(icol,1) = dot__vecR(icol,1) - xgBlockA__vecR(1,icol)*xgBlockB__vecR(1,icol)
          end do
        end if
        !$OMP TARGET UPDATE FROM(dot__vecR)
#else
!FIXME For several compilers, this section doesnt work properly
        !$OMP TARGET UPDATE FROM(dot__vecR,xgBlockA__vecR,xgBlockB__vecR)
        !!$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecR,xgBlockA__vecR,xgBlockB__vecR) PRIVATE(icol,tmp)
        do icol = 1, cols
          tmp=0
          !!$OMP PARALLEL DO REDUCTION(+:tmp) PRIVATE(ii)
          do ii = 1, fact*rows
            tmp = tmp + fact*xgBlockA__vecR(ii,icol)*xgBlockB__vecR(ii,icol)
          end do
          dot__vecR(icol,1)=tmp
        end do
        if (me_g0==1) then
          do icol = 1, cols
            dot__vecR(icol,1) = dot__vecR(icol,1) - xgBlockA__vecR(1,icol)*xgBlockB__vecR(1,icol)
          end do
        end if
        !$OMP TARGET UPDATE TO(dot__vecR)
#endif

        !TODO Port this to GPU (reductions)
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
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        dot__vecC => dot%vecC
#if defined FC_NVHPC
        !$OMP TARGET TEAMS DISTRIBUTE MAP(to:dot__vecC,xgBlockA__vecC,xgBlockB__vecC) PRIVATE(icol,tmp)
        do icol = 1, cols
          tmp=0
          !$OMP PARALLEL DO REDUCTION(+:tmp) PRIVATE(ii)
          do ii = 1, rows
            tmp = tmp + dconjg(xgBlockA__vecC(ii,icol))*xgBlockB__vecC(ii,icol)
          end do
          dot__vecC(icol,1)=tmp
        end do
        !$OMP TARGET UPDATE FROM(dot__vecC)
#else
!FIXME For several compilers, this section doesnt work properly
        !$OMP TARGET UPDATE FROM(xgBlockA__vecC,xgBlockB__vecC)
        do icol = 1, cols
          tmp=0
          !$OMP PARALLEL DO REDUCTION(+:tmp) PRIVATE(ii)
          do ii = 1, rows
            tmp = tmp + dconjg(xgBlockA__vecC(ii,icol))*xgBlockB__vecC(ii,icol)
          end do
          dot__vecC(icol,1)=tmp
        end do
        !$OMP TARGET UPDATE TO(dot__vecC)
#endif

        !TODO Port this to GPU (reductions)
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

#endif

    else

      select case(xgBlockA%space)
      case(SPACE_R,SPACE_CR)
        !$omp parallel do shared(dot,xgBlockA,xgBlockB)
        do icol = 1, xgBlockA%cols
          dot%vecR(icol,1) = fact*ddot(fact*xgBlockA%rows,xgBlockA%vecR(:,icol),1,xgBlockB%vecR(:,icol),1)
        end do
        !$omp end parallel do
        if (xgBlockA%me_g0==1) then
          !$omp parallel do shared(dot,xgBlockA,xgBlockB) &
          !$omp& schedule(static)
          do icol = 1, xgBlockA%cols
            dot%vecR(icol,1) = dot%vecR(icol,1) - ddot(2,xgBlockA%vecR(:,icol),1,xgBlockB%vecR(:,icol),1)
          end do
          !$omp end parallel do
        end if

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
#if defined(FC_NVHPC) || defined(FC_CRAY)
!FIXME zdotc call goes wrong with NVHPC (NVHPC 22.11, MKL 22.3) or CRAY
        !$omp parallel do private(ii,tmp) shared(dot,xgBlockA,xgBlockB)
        do icol = 1, xgBlockA%cols
          tmp=0
          do ii = 1, xgBlockA%rows
            tmp = tmp + dconjg(xgBlockA%vecC(ii,icol))*xgBlockB%vecC(ii,icol)
          end do
          dot%vecC(icol,1)=tmp
        end do
        !$omp end parallel do
#else
        !$omp parallel do shared(dot,xgBlockA,xgBlockB)
        do icol = 1, xgBlockA%cols
          dot%vecC(icol,1) = zdotc(xgBlockA%rows,xgBlockA%vecC(:,icol),1,xgBlockB%vecC(:,icol),1)
        end do
        !$omp end parallel do
#endif

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

    end if ! gpu_option

    call xgBlock_mpi_sum(dot,comm=comm_)

    call timab(tim_colw_dot,2,tsec)

  end subroutine xgBlock_colwiseDotProduct
  !!***

  !!****f* m_xg/xgBlock_colwiseDivision
  !!
  !! NAME
  !! xgBlock_colwiseDivision

  subroutine xgBlock_colwiseDivision(xgBlockA, xgBlockB, divResult, &
    & max_val, max_elt, min_val, min_elt)

    type(xgBlock_t) ,      intent(in   )           :: xgBlockA
    type(xgBlock_t) ,      intent(in   )           :: xgBlockB
    type(xgBlock_t) ,      intent(inout)           :: divResult
    double precision,      intent(inout), optional :: max_val
    integer, dimension(2), intent(inout), optional, target :: max_elt
    double precision,      intent(inout), optional :: min_val
    integer, dimension(2), intent(inout), optional, target :: min_elt

    integer :: irow,icol,rows,cols

#if defined HAVE_GPU
    ! TODO: evaluate if total_size should be a 64 bit integer, i.e.
    ! does spacedim * neigenpairs be larger than 2^31 = 2. 10^9
    integer(kind=c_int32_t)  :: total_size
#if defined HAVE_OPENMP_OFFLOAD
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlockA__vecC(:,:),xgBlockB__vecC(:,:),divResult__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlockA__vecR(:,:),xgBlockB__vecR(:,:),divResult__vecR(:,:)
#endif
#endif
    double precision :: tsec(2)

    call timab(tim_colw_div,1,tsec)

    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)
    call xgBlock_check_gpu_option(xgBlockA,divResult)

    rows = xgBlockA%rows; cols = xgBlockA%cols

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)

      total_size = xgBlockA%rows * xgBlockA%cols

      select case(xgBlockA%space)
      case(SPACE_R)
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
      case(SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
      end select

#else
      ! we shouldn't be here, it means gpu_option was wrongly set to 1 in
      ! input parameter file
#endif

    else if (xgBlockA%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD

      total_size = xgBlockA%rows * xgBlockA%cols
      select case(xgBlockA%space)
      case(SPACE_R)
        xgBlockA__vecR => xgBlockA%vecR
        xgBlockB__vecR => xgBlockB%vecR
        divResult__vecR => divResult%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlockA__vecR,xgBlockB__vecR,divResult__vecR)
        do irow = 1, rows
          do icol = 1, cols
            divResult__vecR(irow,icol) = xgBlockA__vecR(irow,icol)/xgBlockB__vecR(irow,icol)
          end do
        end do
        !FIXME Port this on GPU to avoid copy below ?
        !$OMP TARGET UPDATE FROM(divResult__vecR)
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
        xgBlockA__vecC => xgBlockA%vecC
        xgBlockB__vecC => xgBlockB%vecC
        divResult__vecC => divResult%vecC
#if !defined FC_LLVM
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlockA__vecC,xgBlockB__vecC,divResult__vecC)
        do irow = 1, rows
          do icol = 1, cols
            divResult__vecC(irow,icol) = xgBlockA__vecC(irow,icol)/xgBlockB__vecC(irow,icol)
          end do
        end do
        !FIXME Port this on GPU to avoid copy below ?
        !$OMP TARGET UPDATE FROM(divResult__vecC)
#else
        !FIXME LLVM AOMP 16 doesn't support complex division inside OpenMP !?
        !$OMP TARGET UPDATE FROM(xgBlockA__vecC,xgBlockB__vecC)
        do irow = 1, rows
          do icol = 1, cols
            divResult__vecC(irow,icol) = xgBlockA__vecC(irow,icol)/xgBlockB__vecC(irow,icol)
          end do
        end do
        !$OMP TARGET UPDATE TO(divResult__vecC)
#endif
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
      case(SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
      end select

#endif

    else

      select case(xgBlockA%space)
      case(SPACE_R)
        !$omp parallel do collapse(2) shared(divResult,xgBlockA,xgBlockB) private(icol,irow)
        do icol = 1, cols
          do irow = 1, rows
            divResult%vecR(irow,icol) = xgBlockA%vecR(irow,icol)/xgBlockB%vecR(irow,icol)
          end do
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

        !$omp parallel do collapse(2) shared(divResult,xgBlockA,xgBlockB) private(icol,irow)
        do icol = 1, cols
          do irow = 1, rows
            divResult%vecC(irow,icol) = xgBlockA%vecC(irow,icol)/xgBlockB%vecC(irow,icol)
          end do
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
      case(SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
      end select

    end if ! gpu_option

    call timab(tim_colw_div,2,tsec)

  end subroutine xgBlock_colwiseDivision
  !!***

  !!****f* m_xg/xgBlock_scaleR
  !!
  !! NAME
  !! xgBlock_scaleR

  subroutine xgBlock_scaleR(xgBlock, val, inc)

    type(xgBlock_t) , intent(inout)           :: xgBlock
    double precision, intent(in   )           :: val
    integer         , intent(in   )           :: inc

    integer      :: i,fact
    complex(dpc) :: valc

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_scale,1,tsec)

    valc = dcmplx(val,0.0_dp)

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS .or. xgBlock%gpu_option==ABI_GPU_OPENMP) then

      if ( xgBlock%ldim .eq. xgBlock%rows ) then
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          call abi_gpu_xscal(1, fact*xgBlock%ldim*xgBlock%cols/inc, valc, xgBlock%vecR, inc)
        case (SPACE_C)
          call abi_gpu_xscal(2, xgBlock%ldim*xgBlock%cols/inc, valc, xgBlock%vecC, inc)
        end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          xgBlock__vecR => xgBlock%vecR
          !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock__vecR)
          call abi_gpu_xscal(1, fact*xgBlock%ldim*xgBlock%cols/inc, valc, c_loc(xgBlock__vecR), inc)
          !$OMP END TARGET DATA
        case (SPACE_C)
          xgBlock__vecC => xgBlock%vecC
          !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock__vecC)
          call abi_gpu_xscal(2, xgBlock%ldim*xgBlock%cols/inc, valc, c_loc(xgBlock__vecC), inc)
          !$OMP END TARGET DATA
        end select
#endif

      else
        !FIXME Do loop that calls scal on each column sequentially, might be improved
#if defined HAVE_KOKKOS || defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          do i=1,xgBlock%cols
            call abi_gpu_xscal(1, fact*xgBlock%rows/inc, valc, xgBlock%vecR(:,i), inc)
          end do
        case (SPACE_C)
          do i=1,xgBlock%cols
            call abi_gpu_xscal(2, xgBlock%rows/inc, valc, xgBlock%vecC(:,i), inc)
          end do
        end select
#elif defined HAVE_OPENMP_OFFLOAD
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          xgBlock__vecR => xgBlock%vecR
          do i=1,xgBlock%cols
            !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock__vecR)
            call abi_gpu_xscal(1, fact*xgBlock%rows/inc, valc, c_loc(xgBlock__vecR(1,i)), inc)
            !$OMP END TARGET DATA
          end do
        case (SPACE_C)
          xgBlock__vecC => xgBlock%vecC
          do i=1,xgBlock%cols
            !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock__vecC)
            call abi_gpu_xscal(2, xgBlock%rows/inc, valc, c_loc(xgBlock__vecC(1,i)), inc)
            !$OMP END TARGET DATA
          end do
        end select
#endif
      end if

    else

      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          call dscal(fact*xgBlock%ldim*xgBlock%cols/inc,val,xgBlock%vecR,inc)
        case (SPACE_C)
          call zdscal(xgBlock%ldim*xgBlock%cols/inc,val,xgBlock%vecC,inc)
        end select
      else
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          !$omp parallel do
          do i=1,xgBlock%cols
            call dscal(fact*xgBlock%rows/inc,val,xgBlock%vecR(:,i),inc)
          end do
        case (SPACE_C)
          !$omp parallel do
          do i=1,xgBlock%cols
            call zdscal(xgBlock%rows/inc,val,xgBlock%vecC(:,i),inc)
          end do
        end select
      end if

    end if

    call timab(tim_scale,2,tsec)

  end subroutine xgBlock_scaleR
  !!***

  !!****f* m_xg/xgBlock_scaleC
  !!
  !! NAME
  !! xgBlock_scaleC

  subroutine xgBlock_scaleC(xgBlock, val, inc)

    type(xgBlock_t), intent(inout)           :: xgBlock
    complex(kind=8), intent(in   )           :: val
    integer        , intent(in   )           :: inc

    integer :: i
    double precision :: tsec(2)

    call timab(tim_scale,1,tsec)

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          ABI_ERROR("Scaling real vector with a complex not possible")
        case (SPACE_C)
          call abi_gpu_xscal(2, xgBlock%ldim*xgBlock%cols/inc, val, xgBlock%vecC, inc)
        end select
      else

        ! TODO (PK) : evaluate if it is really necessary to deal with this case
        ABI_ERROR("Scaling a xgBlock when xgBlock%ldim != xgBlock%rows is not implemented for GPU. FIX ME if needed.")

      end if
#endif

    else if (xgBlock%gpu_option==ABI_GPU_OPENMP) then

#if defined(HAVE_GPU) && defined(HAVE_OPENMP_OFFLOAD)
      if ( xgBlock%ldim .eq. xgBlock%rows ) then
        select case(xgBlock%space)
        case (SPACE_R,SPACE_CR)
          ABI_ERROR("Scaling real vector with a complex not possible")
        case (SPACE_C)
          call abi_gpu_xscal(2, xgBlock%ldim*xgBlock%cols/inc, val, xgBlock%vecC, inc)
        end select
      else

        ! TODO (PK) : evaluate if it is really necessary to deal with this case
        ABI_BUG("Scaling a xgBlock when xgBlock%ldim != xgBlock%rows is not implemented for GPU. FIX ME if needed.")

      end if
#endif
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

    call timab(tim_scale,2,tsec)

  end subroutine xgBlock_scaleC
  !!***

  !!****f* m_xg/xgBlock_transpose
  !!
  !! NAME
  !! xgBlock_transpose

  subroutine xgBlock_transpose(xgBlockI,xgBlockO)

    type(xgBlock_t) , intent(inout) :: xgBlockI,xgBlockO
    integer :: nrows,ncols

    if (xgBlockI%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockI,xgBlockO)

    nrows = xgBlockI%rows
    ncols = xgBlockI%cols
    if (nrows/=xgBlockO%cols) then
      ABI_ERROR('nrowsI/=ncolsO')
    end if
    if (ncols/=xgBlockO%rows) then
      ABI_ERROR('ncolsI/=nrowsO')
    end if

    select case(xgBlockI%space)
    case (SPACE_R)
      xgBlockO%vecR = TRANSPOSE(xgBlockI%vecR)
    case (SPACE_CR)
      ABI_ERROR("Not implemented")
    case (SPACE_C)
      xgBlockO%vecC = TRANSPOSE(CONJG(xgBlockI%vecC))
    end select

  end subroutine xgBlock_transpose
  !!***

  !!****f* m_xg/xgBlock_r2c
  !!
  !! NAME
  !! xgBlock_r2c

  subroutine xgBlock_r2c(xgBlockR,xgBlockC,nspinor)

    integer         , intent(in) :: nspinor
    type(xgBlock_t) , intent(in) :: xgBlockR
    type(xgBlock_t) , intent(inout) :: xgBlockC
    integer :: nrows,ncols,col
    double precision,allocatable :: zeros(:)

    if (xgBlockR%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockR,xgBlockC)

    if (space(xgBlockR)/=SPACE_R) then
      ABI_ERROR('space(xgBlockR)/=SPACE_R')
    end if
    if (space(xgBlockC)/=SPACE_C .and. space(xgBLockC)/=SPACE_CR) then
      ABI_ERROR('space(xgBlockC)/=SPACE_C')
    end if
    if (space(xgBlockC)==SPACE_CR .and. nspinor/=1) then
      ABI_ERROR('This should not happen')
    end if

    nrows = xgBlockR%rows
    ncols = xgBlockR%cols
    if (nspinor*nrows/=xgBlockC%rows) then
      ABI_ERROR('nspinor*nrowsR/=nrowsC')
    end if
    if (nspinor*ncols/=xgBlockC%cols) then
      ABI_ERROR('nspinor*ncolsR/=ncolsC')
    end if

    if (space(xgBlockC)==SPACE_C) then
      ABI_MALLOC(zeros,(nrows))
      zeros=zero
      do col = 1,ncols
        xgBlockC%vecC(1:nrows,col) = dcmplx(xgBlockR%vecR(1:nrows,col),zeros)
        if (nspinor==2) then
          xgBlockC%vecC(1+nrows:2*nrows,col      ) = dcmplx(zeros,zeros)
          xgBlockC%vecC(1      :  nrows,col+ncols) = dcmplx(zeros,zeros)
          xgBlockC%vecC(1+nrows:2*nrows,col+ncols) = dcmplx(xgBlockR%vecR(1:nrows,col),zeros)
        end if
      end do
      ABI_FREE(zeros)
    else ! space(C)==SPACE_CR
      do col = 1,ncols
        xgBlockC%vecR(1:2*nrows-1:2,col) = xgBlockR%vecR(1:nrows,col)
        xgBlockC%vecR(2:2*nrows  :2,col) = zero
      end do
    end if

  end subroutine xgBlock_r2c
  !!***

  !!****f* m_xg/xgBlock_c2r
  !!
  !! NAME
  !! xgBlock_c2r

  subroutine xgBlock_c2r(xgBlockC,xgBlockR)

    type(xgBlock_t) , intent(in) :: xgBlockC
    type(xgBlock_t) , intent(inout) :: xgBlockR
    integer :: nrows,ncols,col

    if (xgBlockR%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockR,xgBlockC)

    if (space(xgBlockR)/=SPACE_R) then
      ABI_ERROR('space(xgBlockR)/=SPACE_R')
    end if
    if (space(xgBlockC)/=SPACE_C .and. space(xgBLockC)/=SPACE_CR) then
      ABI_ERROR('space(xgBlockC)/=SPACE_C')
    end if

    nrows = xgBlockR%rows
    ncols = xgBlockR%cols
    if (nrows/=xgBlockC%rows) then
      ABI_ERROR('nrowsR/=nrowsC')
    end if
    if (ncols/=xgBlockC%cols) then
      ABI_ERROR('ncolsR/=ncolsC')
    end if

    if (space(xgBlockC)==SPACE_C) then
      do col = 1,ncols
        xgBlockR%vecR(1:nrows,col) = dble(xgBlockC%vecC(1:nrows,col))
      end do
    else ! space(C)==SPACE_CR
      do col = 1,ncols
        xgBlockR%vecR(1:nrows,col) = xgBlockC%vecR(1:2*nrows-1:2,col)
      end do
    end if

  end subroutine xgBlock_c2r
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

  !!****f* m_xg/xgBlock_check
  !!
  !! NAME
  !! xgBlock_check

  subroutine xgBlock_check(X, Y, fact_col)

    type(xgBlock_t) , intent(in) :: X
    type(xgBlock_t) , intent(in) :: Y
    integer,optional, intent(in) :: fact_col

    integer :: fact_col_

    fact_col_ = 1
    if (present(fact_col)) then
      fact_col_ = fact_col
    end if
    if (X%space/=Y%space) then
      ABI_ERROR('X%space/=Y%space')
    end if
    if (X%rows/=Y%rows) then
      ABI_ERROR('X%rows/=Y%rows')
    end if
    if (fact_col_*X%cols/=Y%cols) then
      ABI_ERROR('X%cols/=Y%cols')
    end if

  end subroutine xgBlock_check
  !!***

  !!****f* m_xg/xgBlock_check_gpu_option
  !!
  !! NAME
  !! xgBlock_check_gpu_option

  subroutine xgBlock_check_gpu_option(X, Y)

    type(xgBlock_t) , intent(in) :: X
    type(xgBlock_t) , intent(in) :: Y

    if (X%gpu_option/=Y%gpu_option) then
      ABI_ERROR('X%gpu_option /= Y%gpu_option')
    end if

  end subroutine xgBlock_check_gpu_option
  !!***

  !!****f* m_xg/xgBlock_copy_to_gpu
  !!
  !! NAME
  !! xgBlock_copy_to_gpu

  subroutine xgBlock_copy_to_gpu(xgBlock)
    type(xgBlock_t), target, intent(in   ) :: xgBlock
#if defined(HAVE_GPU) && defined(HAVE_OPENMP_OFFLOAD)
    integer(c_size_t) :: size
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      xgBlock__vecR => xgBlock%vecR
      !$OMP TARGET UPDATE TO(xgBlock__vecR)
    case (SPACE_C)
      xgBlock__vecC => xgBlock%vecC
      !$OMP TARGET UPDATE TO(xgBlock__vecC)
    end select
#else
    ABI_UNUSED_A(xgBlock)
#endif

  end subroutine xgBlock_copy_to_gpu
  !!***

  !!****f* m_xg/xgBlock_copy_from_gpu
  !!
  !! NAME
  !! xgBlock_copy_from_gpu

  subroutine xgBlock_copy_from_gpu(xgBlock)
    type(xgBlock_t), target, intent(in   ) :: xgBlock
#if defined(HAVE_GPU) && defined(HAVE_OPENMP_OFFLOAD)
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      xgBlock__vecR => xgBlock%vecR
      !$OMP TARGET UPDATE FROM(xgBlock__vecR)
    case (SPACE_C)
      xgBlock__vecC => xgBlock%vecC
      !$OMP TARGET UPDATE FROM(xgBlock__vecC)
    end select
#else
    ABI_UNUSED_A(xgBlock)
#endif

  end subroutine xgBlock_copy_from_gpu
  !!***

  !!****f* m_xg/xgBlock_reshape
  !!
  !! NAME
  !! xgBlock_reshape

  subroutine xgBlock_reshape(xgBlock,newrows,newcols)
    use, intrinsic :: iso_c_binding
    type(xgBlock_t), intent(inout) :: xgBlock
    integer        , intent(in   ) :: newrows
    integer        , intent(in   ) :: newcols
    integer :: fact,newshape(2)
    type(c_ptr) :: cptr

    if ( xgBlock%rows*xgBlock%cols /= newrows*newcols ) then
      write(std_out,*) "xgBlock%rows", xgBlock%rows
      write(std_out,*) "xgBlock%cols", xgBlock%cols
      write(std_out,*) "newrows", newrows
      write(std_out,*) "newcols", newcols
      write(std_out,*) "xgBlock%rows*xgBlock%cols", xgBlock%rows*xgBlock%cols
      write(std_out,*) "newrows*newcols", newrows*newcols
      ABI_ERROR("Bad shape")
    end if

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    xgBlock%LDim = newrows+( (xgBlock%LDim-xgBlock%rows)* xgBlock%cols)/newcols
    xgBlock%rows = newrows
    xgBlock%cols = newcols
    newshape(1) = newrows
    newshape(2) = newcols
    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      cptr = getClocR(fact*xgBlock%LDim,xgBlock%cols,xgBlock%vecR)
      call c_f_pointer(cptr,xgBlock%vecR,newshape)
    case (SPACE_C)
      cptr = getClocC(xgBlock%LDim,xgBlock%cols,xgBlock%vecC)
      call c_f_pointer(cptr,xgBlock%vecC,newshape)
    end select

  end subroutine xgBlock_reshape
  !!***

  !!****f* m_xg/xgBlock_reshape_spinor
  !!
  !! NAME
  !! xgBlock_reshape_spinor

  subroutine xgBlock_reshape_spinor(xgBlock,xgBlock_spinor,nspinor,option)
    use iso_c_binding
    integer, intent(in   ) :: nspinor,option
    type(xgBlock_t), intent(in   ) :: xgBlock
    type(xgBlock_t), intent(inout) :: xgBlock_spinor

    integer :: nrows,ncols

    if (nspinor/=1.and.nspinor/=2) then
      ABI_ERROR('It should not happen : nspinor must be 1 or 2')
    end if
    if (xgBlock%space==SPACE_CR.and.nspinor==2) then
      ABI_ERROR('It should not happen : space_CR cannot be used with nspinor=2')
    end if

    nrows = rows(xgBlock)
    ncols = cols(xgBlock)

    if (option==COLS2ROWS) then
      if (modulo(ncols,nspinor)/=0) then
        ABI_ERROR('nspinor should divide the number of cols')
      end if
      call xgBlock_setBlock(xgBlock,xgBlock_spinor,nrows,ncols)
      if (nspinor>1) call xgBlock_reshape(xgBlock_spinor,nrows*nspinor,ncols/nspinor)
    else if (option==ROWS2COLS) then
      if (modulo(nrows,nspinor)/=0) then
        ABI_ERROR('nspinor should divide the number of rows')
      end if
      call xgBlock_setBlock(xgBlock,xgBlock_spinor,nrows,ncols)
      if (nspinor>1) call xgBlock_reshape(xgBlock_spinor,nrows/nspinor,ncols*nspinor)
    else
      ABI_ERROR('bad option value')
    end if

  end subroutine xgBlock_reshape_spinor
  !!***

  !!****f* m_xg/xgBlock_free_reshape
  !!
  !! NAME
  !! xgBlock_free_reshape
  subroutine xgBlock_free_reshape(xgBlock,newrows,newcols,newldim,new_me_g0)
    use iso_c_binding
    type(xgBlock_t) , intent(inout) :: xgBlock
    integer         , intent(in   ) :: newrows,newcols
    integer,optional, intent(in   ) :: newldim
    integer,optional, intent(in   ) :: new_me_g0

    integer :: newshape(2)
    type(c_ptr) :: cptr

    if (newrows<1.or.newcols<1) then
      ABI_ERROR("Bad new shape")
    end if

    if ( xgBLock%Ldim*xgBlock%cols < newrows*newcols ) then
      write(std_out,*) "xgBLock%rows", xgBLock%rows
      write(std_out,*) "xgBlock%cols", xgBlock%cols
      write(std_out,*) "newrows", newrows
      write(std_out,*) "newcols", newcols
      write(std_out,*) "xgBLock%rows*xgBlock%cols", xgBLock%rows*xgBlock%cols
      write(std_out,*) "newrows*newcols", newrows*newcols
      ABI_ERROR("Bad shape (ldim*cols<newrows*newcols")
    end if

    xgBlock%LDim = newrows
    if (present(newldim)) then
      if (newldim<newrows) then
        ABI_ERROR("newldim<newrows")
      end if
      if ( xgBLock%Ldim*xgBlock%cols < newldim*newcols ) then
        ABI_ERROR("Bad shape (ldim*cols<newldim*newcols")
      end if
      xgBlock%LDim = newldim
    end if
    xgBlock%rows = newrows
    xgBlock%cols = newcols
    newshape(1)  = newrows
    newshape(2)  = newcols
    select case(xgBLock%space)
    case (SPACE_R)
      cptr = getClocR(xgBlock%LDim,xgBlock%cols,xgBlock%vecR)
      call c_f_pointer(cptr,xgBlock%vecR,newshape)
    case (SPACE_CR)
      cptr = getClocR(2*xgBlock%LDim,xgBlock%cols,xgBlock%vecR)
      newshape(1) = 2*newshape(1)
      call c_f_pointer(cptr,xgBlock%vecR,newshape)
      if (present(new_me_g0)) then
        xgBlock%me_g0=new_me_g0
      end if
    case (SPACE_C)
      cptr = getClocC(xgBlock%LDim,xgBlock%cols,xgBlock%vecC)
      call c_f_pointer(cptr,xgBlock%vecC,newshape)
    end select
  end subroutine xgBlock_free_reshape
  !!***

  !!****f* m_xg/xgBlock_zero
  !!
  !! NAME
  !! xgBlock_zero

  subroutine xgBlock_zero(xgBlock)

    type(xgBlock_t), intent(inout) :: xgBlock

    integer :: i,fact
#if defined HAVE_GPU
    integer(C_SIZE_T) :: byte_count
#endif

#if defined HAVE_OPENMP_OFFLOAD && !defined HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
    complex(dpc), ABI_CONTIGUOUS pointer :: xgBlock__vecC(:,:)
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)
    integer :: rows,cols,iblock,jblock
#endif
    double precision :: tsec(2)

    call timab(tim_zero,1,tsec)

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then

#if defined HAVE_GPU && defined HAVE_KOKKOS
      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        byte_count = fact * xgBlock%ldim * xgBlock%cols * dp
        call gpu_memset(c_loc(xgBlock%vecR), 0, byte_count)
      case (SPACE_C)
        byte_count = xgBlock%ldim * xgBlock%cols * 2 * dpc ! Note the factor 2, needed here!
        call gpu_memset(c_loc(xgBlock%vecC), 0, byte_count)
      end select
#endif

    else if (xgBlock%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_OPENMP_OFFLOAD
#ifdef HAVE_OPENMP_OFFLOAD_DATASTRUCTURE
      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        byte_count = int(fact, c_size_t) * xgBlock%ldim * xgBlock%cols * dp
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock%vecR)
        call gpu_memset(c_loc(xgBlock%vecR), 0, byte_count)
        !$OMP END TARGET DATA
      case (SPACE_C)
        byte_count = int(xgBlock%ldim, c_size_t) * xgBlock%cols * 2 * dpc ! Note the factor 2, needed here!
        !$OMP TARGET DATA USE_DEVICE_ADDR(xgBlock%vecC)
        call gpu_memset(c_loc(xgBlock%vecC), 0, byte_count)
        !$OMP END TARGET DATA
      end select
#else
!FIXME For several compilers, OMP doesn't work correctly with structured types, so use pointers
      rows = xgBlock%rows; cols = xgBlock%cols
      select case(xgBlock%space)
      case (SPACE_R,SPACE_CR)
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlock__vecR)
        do iblock = 1, cols
          do jblock = 1, fact * rows
            xgBlock__vecR(jblock,iblock) = zero
          end do
        end do
      case (SPACE_C)
        xgBlock__vecC => xgBlock%vecC
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(to:xgBlock__vecC)
        do iblock = 1, cols
          do jblock = 1, fact * rows
            xgBlock__vecC(jblock,iblock) = dcmplx(0,0)
          end do
        end do
      end select
#endif
#endif

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

    call timab(tim_zero,2,tsec)

  end subroutine xgBlock_zero
  !!***

  !!****f* m_xg/xgBlock_zerotri
  !!
  !! NAME
  !! xgBlock_zerotri

  subroutine xgBlock_zerotri(xgBlockA,uplo)
    use iso_c_binding
    type(xgBlock_t), intent(inout) :: xgBlockA
    character, intent(in) :: uplo
    integer :: j
    integer :: i
    integer :: nn
    integer :: col

    nn = xgBlockA%cols
    if (xgBlockA%rows/=nn) then
      ABI_ERROR('rows should be equal to cols!')
    end if

    if (xgBlockA%gpu_option==ABI_GPU_KOKKOS) then
      ABI_ERROR('Not implemented for GPU Kokkos')
    end if
    if (xgBlockA%gpu_option==ABI_GPU_OPENMP) then
      call xgBlock_copy_from_gpu(xgBlockA) !FIXME Avoid that transfer
    end if

    select case(uplo)
    case ('u','U')
      select case(xgBlockA%space)
      case (SPACE_R)
        do j = 1, nn
          col = (j*(j-1))/2
          do i = j+1,nn
            xgBlockA%vecR(i,j) = zero
          end do
        end do
      case (SPACE_CR)
        ABI_ERROR("Not implemented")
      case (SPACE_C)
        do j = 1, nn
          col = (j*(j-1))/2
          do i = j+1,nn
            xgBlockA%vecC(i,j) = czero
          end do
        end do
      end select

    case ('l','L')
      select case(xgBlockA%space)
      case (SPACE_R)
        do j = 1, nn
          col = ((2*xgBlockA%cols-j)*(j-1))/2
          do i = 1, j-1
            xgBlockA%vecR(i,j) = zero
          end do
        end do
      case (SPACE_CR)
        ABI_ERROR("Not implemented")
      case (SPACE_C)
        do j = 1, nn
          col = ((2*xgBlockA%cols-j)*(j-1))/2
          do i = 1, j-1
            xgBlockA%vecC(i,j) = czero
          end do
        end do
      end select
    case default
      ABI_ERROR("Error for zerotri")
    end select

    if (xgBlockA%gpu_option==ABI_GPU_OPENMP) then
      call xgBlock_copy_to_gpu(xgBlockA) !FIXME Avoid that transfer
    end if

  end subroutine xgBlock_zerotri
  !!***

  !!****f* m_xg/xgBlock_zero_im_g0
  !!
  !! NAME
  !! xgBlock_zero_im_g0

  subroutine xgBlock_zero_im_g0(xgBlock)

    type(xgBlock_t), intent(inout) :: xgBlock

    integer :: ii,cols

#if defined HAVE_OPENMP_OFFLOAD
    real(dp), ABI_CONTIGUOUS pointer :: xgBlock__vecR(:,:)
#endif
    double precision :: tsec(2)

    call timab(tim_zero_im_g0,1,tsec)
    if (xgBlock%space==SPACE_CR) then

      cols = xgBlock%cols

      if (xgBlock%me_g0<0) then
        ABI_ERROR("xgBlock me_g0 is not initialized")
      else if (xgBlock%me_g0==1) then

        if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then

          ABI_ERROR('Not implemented')

        else if (xgBlock%gpu_option==ABI_GPU_OPENMP) then

#if defined HAVE_OPENMP_OFFLOAD
          xgBlock__vecR => xgBlock%vecR
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(to:xgBlock__vecR)
          do ii = 1, cols
            xgBlock__vecR(2,ii) = zero
          end do
#endif

        else

          !$omp parallel do
          do ii = 1, cols
            xgBlock%vecR(2,ii) = zero
          end do

        end if ! gpu_option

      end if ! me_g0>=0
    end if ! SPACE_CR
    call timab(tim_zero_im_g0,2,tsec)

  end subroutine xgBlock_zero_im_g0
  !!***

  !!****f* m_xg/xgBlock_invert
  !!
  !! NAME
  !! xgBlock_invert

  subroutine xgBlock_invert(xgBlockA,xgBlockW,xg_input)

    type(xgBlock_t), intent(inout) :: xgBlockA,xgblockW
    type(xgBlock_t), optional,intent(in) :: xg_input

    integer :: nn,nrows,info,ldim
    integer,allocatable :: ipiv(:)
    double precision :: tsec(2)

    call timab(tim_invert,1,tsec)

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)
    if (present(xg_input)) then
      call xgBlock_check_gpu_option(xgBlockA,xg_input)
    end if

    call xgBlock_getsize(xgBlockA,nrows,nn,ldim=ldim)
    if (nrows/=nn) then
      ABI_ERROR('nrows/=ncols')
    end if
    if (cols(xgBlockW)/=nn) then
      ABI_ERROR('cols(xgBlockW)/=ncols')
    end if
    if (rows(xgBlockW)/=nn) then
      ABI_ERROR('rows(xgBlockW)/=nrows')
    end if

    if (present(xg_input)) then
      if (cols(xg_input)/=nn) then
        ABI_ERROR('cols(xgBlockB)/=ncols')
      end if
      if (rows(xg_input)/=nn) then
        ABI_ERROR('rows(xgBlockB)/=nrows')
      end if
      call xgBlock_copy(xg_input,xgBlockA)
    end if

    ABI_MALLOC(ipiv,(nn))

    select case(xgBlockA%space)
    case (SPACE_R)
      call DGETRF(nn,nn,xgBlockA%vecR,ldim,ipiv,info)
      if (info==0) then
        call DGETRI(nn,xgBlockA%vecR,ldim,ipiv,xgBlockW%vecR,nn,info)
      else
        ABI_ERROR('info/=0 : something bad happened in xgetrf')
      end if
    case (SPACE_CR)
      ABI_ERROR("Not implemented")
      !call DGETRF(nn,nn,xgBlockA%vecR,ldim,ipiv,info)
      !if (info==0) then
      !  call DGETRI(nn,xgBlockA%vecR,ldim,ipiv,xgBlockW%vecR,nn,info)
      !else
      !  ABI_ERROR('info/=0 : something bad happened in xgetrf')
      !end if
    case (SPACE_C)
      call ZGETRF(nn,nn,xgBlockA%vecC,ldim,ipiv,info)
      if (info==0) then
        call ZGETRI(nn,xgBlockA%vecC,ldim,ipiv,xgBlockW%vecC,nn,info)
      else
        ABI_ERROR('info/=0 : something bad happened in xgetrf')
      end if
    end select

    ABI_FREE(ipiv)

    if (info/=0) then
      ABI_ERROR('info/=0 : something bad happened in xgetri')
    end if

    call timab(tim_invert,2,tsec)

  end subroutine xgBlock_invert
  !!***

  !!****f* m_xg/xgBlock_invert_sy
  !!
  !! NAME
  !! xgBlock_invert_sy

  subroutine xgBlock_invert_sy(xgBlockA,xgBlockW,xg_input)

    type(xgBlock_t), intent(inout)       :: xgBlockA,xgblockW
    type(xgBlock_t), optional,intent(in) :: xg_input

    integer :: nn,nrows,ldim,info
    integer :: ii,jj
    integer,allocatable :: ipiv(:)
    double precision :: tsec(2)

    call timab(tim_invert_sy,1,tsec)

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,xgBlockW)
    if (present(xg_input)) then
      call xgBlock_check_gpu_option(xgBlockA,xg_input)
    end if

    call xgBlock_getsize(xgBlockA,nrows,nn,ldim=ldim)
    if (nrows/=nn) then
      ABI_ERROR('nrows/=ncols')
    end if
    if (cols(xgBlockW)/=nn) then
      ABI_ERROR('cols(xgBlockW)/=ncols')
    end if
    if (rows(xgBlockW)/=nn) then
      ABI_ERROR('rows(xgBlockW)/=nrows')
    end if

    if (present(xg_input)) then
      if (cols(xg_input)/=nn) then
        ABI_ERROR('cols(xg_input)/=ncols')
      end if
      if (rows(xg_input)/=nn) then
        ABI_ERROR('rows(xg_input)/=nrows')
      end if
      call xgBlock_copy(xg_input,xgBlockA)
    end if

    ABI_MALLOC(ipiv,(nn))

    select case(xgBlockA%space)
    case (SPACE_R)
      call DSYTRF('U',nn,xgBlockA%vecR,ldim,ipiv,xgBlockW%vecR,nn,info)
      if (info==0) then
        call DSYTRI('U',nn,xgBlockA%vecR,ldim,ipiv,xgBlockW%vecR,info)
        if (info/=0) then
          ABI_ERROR('info/=0 : something bad happened in dsytri')
        end if
      else
        ABI_ERROR('info/=0 : something bad happened in dsytrf')
      end if
      ! complete the matrix
      do ii=1, nn
        do jj=1, ii-1
           xgBlockA%vecR(ii,jj) = xgBlockA%vecR(jj,ii)
        end do
      end do
    case (SPACE_CR)
      ABI_ERROR("Not implemented")
    case (SPACE_C)
      call ZHETRF('U',nn,xgBlockA%vecC,ldim,ipiv,xgBlockW%vecC,nn,info)
      if (info==0) then
        call ZHETRI('U',nn,xgBlockA%vecC,ldim,ipiv,xgBlockW%vecC,info)
        if (info/=0) then
          ABI_ERROR('info/=0 : something bad happened in zhetri')
        end if
      else
        ABI_ERROR('info/=0 : something bad happened in zhetrf')
      end if
      ! complete the matrix
      do ii=1, nn
        do jj=1, ii-1
           xgBlockA%vecC(ii,jj) = CONJG(xgBlockA%vecC(jj,ii))
        end do
      end do
    end select

    ABI_FREE(ipiv)

    call timab(tim_invert_sy,2,tsec)

  end subroutine xgBlock_invert_sy
  !!***

  !!****f* m_xg/xgBlock_invert_tri
  !!
  !! NAME
  !! xgBlock_invert_tri

  subroutine xgBlock_invert_tri(uplo,diag,xgBlock)

    type(xgBlock_t), intent(inout) :: xgBlock
    character, intent(in) :: uplo,diag

    integer :: info
    double precision :: tsec(2)

    call timab(tim_invertri,1,tsec)

    if (xgBlock%rows/=xgBlock%cols) then
      ABI_ERROR('nrows/=ncols')
    end if

    if (xgBlock%gpu_option==ABI_GPU_KOKKOS) then
      ABI_ERROR('Not implemented for GPU Kokkos')
    end if
    if (xgBlock%gpu_option==ABI_GPU_OPENMP) then
      call xgBlock_copy_from_gpu(xgBlock) !FIXME Avoid that transfer
    end if

    select case(xgBlock%space)
    case (SPACE_R)
      call DTRTRI(uplo,diag,xgBlock%rows,xgBlock%vecR,xgBlock%LDim,info)
    case (SPACE_CR)
      ABI_ERROR("Not implemented")
    case (SPACE_C)
      call ZTRTRI(uplo,diag,xgBlock%rows,xgBlock%vecC,xgBlock%LDim,info)
    end select

    if (xgBlock%gpu_option==ABI_GPU_OPENMP) then
      call xgBlock_copy_to_gpu(xgBlock) !FIXME Avoid that transfer
    end if

    if (info/=0) then
      ABI_ERROR('info/=0 : something bad happened in xtrtri')
    end if

    call timab(tim_invertri,2,tsec)

  end subroutine xgBlock_invert_tri
  !!***

  !!****f* m_xg/xgBlock_yxpa
  !!
  !! NAME
  !! xgBlock_yxpa

  subroutine xgBlock_yxpa(xgBlockA,xgBlockB,aa)

    double precision,intent(in) :: aa
    type(xgBlock_t), intent(in) :: xgBlockA
    type(xgBlock_t), intent(inout) :: xgBlockB

    integer :: nrows

    if (xgBlockA%gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(xgBlockA,xgBlockB)

    nrows = rows(xgBlockA)

    select case(xgBlockA%space)
    case (SPACE_R)
      xgBlockB%vecR(1:nrows,:) = xgBlockA%vecR(1:nrows,:) + aa
    case (SPACE_CR)
      xgBlockB%vecR(1:2*nrows,:) = xgBlockA%vecR(1:2*nrows,:) + aa
    case (SPACE_C)
      xgBlockB%vecC(1:nrows,:) = xgBlockA%vecC(1:nrows,:) + aa*(1.0d0,1.0d0)
    end select

  end subroutine xgBlock_yxpa
  !!***

  !!****f* m_xg/xgBlock_one
  !!
  !! NAME
  !! xgBlock_one

  subroutine xgBlock_one(xgBlock)

    type(xgBlock_t), intent(inout) :: xgBlock
    integer :: i

    select case(xgBlock%space)
    case (SPACE_R)
      !$omp parallel do
      do i = 1, min(xgBlock%rows,xgBlock%cols)
        xgBlock%vecR(i,i) = 1.d0
      end do
    case (SPACE_CR)
      !$omp parallel do
      do i = 1, min(2*xgBlock%rows,xgBlock%cols)
        xgBlock%vecR(2*i-1,i) = 1.d0
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
    case (SPACE_R)
      select case(diag%space)
      case (SPACE_R)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecR(i,i) = diag%vecR(i,1)
        end do
      case (SPACE_CR)
          ABI_ERROR('Not implemented for SPACE_CR')
      case (SPACE_C)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecR(i,i) = dble(diag%vecC(i,1))
        end do
      end select
    case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
    case (SPACE_C)
      select case(diag%space)
      case (SPACE_R)
        !$omp parallel do
        do i = 1, min(xgBlock%rows,xgBlock%cols)
          xgBlock%vecC(i,i) = dcmplx(diag%vecR(i,1))
        end do
      case (SPACE_CR)
          ABI_ERROR('Not implemented for SPACE_CR')
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
    case (SPACE_R)
      !$omp parallel do
      do i = 1, xgBlock%cols
        diag%vecR(i,1) = xgBlock%vecR(i,i)
      end do
    case (SPACE_CR)
        ABI_ERROR('Not implemented for SPACE_CR')
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

  !!****f* m_xg/xgBlock_minmax
  !!
  !! NAME
  !! xgBlock_minmax

  subroutine xgBlock_minmax(xgBlock,minimum,maximum,row_bound)

    type(XgBlock_t) , intent(in)  :: xgBlock
    double precision, intent(out) :: minimum,maximum
    integer,optional, intent(in)  :: row_bound

    integer :: row_bound_,fact
    double precision :: tsec(2)

    call timab(tim_minmax,1,tsec)

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    row_bound_ = fact*xgBlock%rows
    if (present(row_bound)) then
      if (row_bound<1.or.row_bound>fact*xgBlock%rows) then
        ABI_ERROR('Bad row_bound')
      else
        row_bound_ = fact*row_bound
      end if
    end if

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      minimum = minval(xgBlock%vecR(:row_bound_,:))
      maximum = maxval(xgBlock%vecR(:row_bound_,:))
    case (SPACE_C)
      minimum = minval(abs(xgBlock%vecC(:row_bound_,:)))
      maximum = maxval(abs(xgBlock%vecC(:row_bound_,:)))
    end select

    call timab(tim_minmax,2,tsec)

  end subroutine xgBlock_minmax
  !!***

  !!****f* m_xg/xgBlock_average
  !!
  !! NAME
  !! xgBlock_average

  subroutine xgBlock_average(xgBlock,average)

    type(XgBlock_t) , intent(in)  :: xgBlock
    double precision, intent(out) :: average
    complex(kind=8) :: averageC
    integer :: i,fact

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      average = 0.d0
      do i = 1, xgBlock%cols
        average = average + sum(xgBlock%vecR(1:fact*xgBlock%rows,i))
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
    integer :: i,fact

    fact = 1 ; if (xgBlock%space==SPACE_CR) fact = 2

    call xgBlock_average(xgBlock,average)
    select case(xgBlock%space)
    case (SPACE_R,SPACE_CR)
      deviation = 0.d0
      do i = 1, xgBlock%cols
        deviation = deviation + sum((xgBlock%vecR(1:fact*xgBlock%rows,i)-average)*(xgBlock%vecR(1:fact*xgBlock%rows,i)-average))
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
    real(dp), allocatable :: vecR_tmp(:)

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    complex(dpc), pointer :: xgBlock__vecC(:,:)
    real(dp), pointer :: xgBlock__vecR(:,:)
#endif

    select case(xgBlock%space)
    case (SPACE_R)
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      if (xgBlock%gpu_option==ABI_GPU_OPENMP) then
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET UPDATE FROM(xgBlock__vecR)
      end if
#endif
      write(ccols,'(i4)') xgBlock%cols
      !fstring = '(1x,'//trim(adjustl(ccols))//'ES22.14)'
      fstring = '(1x,'//trim(adjustl(ccols))//'f24.14)'
      do i = 1, xgBlock%rows
        write(outunit,fstring) (/ (xgBlock%vecR(i,j), j = 1, xgBlock%cols) /)
      end do
    case (SPACE_CR)
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      if (xgBlock%gpu_option==ABI_GPU_OPENMP) then
        xgBlock__vecR => xgBlock%vecR
        !$OMP TARGET UPDATE FROM(xgBlock__vecR)
      end if
#endif
      ABI_MALLOC(vecR_tmp,(2*xgBlock%cols))
      write(ccols,'(i4)') 2*xgBlock%cols
      !fstring = '(1x,2(1x,'//trim(adjustl(ccols))//'ES22.14))'
      fstring = '(1x,2(1x,'//trim(adjustl(ccols))//'f24.14))'
      do i = 1, xgBlock%rows
        do j = 1, xgBlock%cols
           vecR_tmp(2*j-1) = xgBlock%vecR(2*i-1,j)
           vecR_tmp(2*j  ) = xgBlock%vecR(2*i  ,j)
        end do
        write(outunit,fstring) (/ (vecR_tmp(j), j = 1, 2*xgBlock%cols) /)
      end do
      ABI_FREE(vecR_tmp)
    case (SPACE_C)
#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
      if (xgBlock%gpu_option==ABI_GPU_OPENMP) then
        xgBlock__vecC => xgBlock%vecC
        !$OMP TARGET UPDATE FROM(xgBlock__vecC)
      end if
#endif
      write(ccols,'(i4)') xgBlock%cols
      !fstring = '(1x,2(1x,'//trim(adjustl(ccols))//'ES22.14))'
      fstring = '(1x,2(1x,'//trim(adjustl(ccols))//'f24.14))'
      do i = 1, xgBlock%rows
        write(outunit,fstring) (/ (xgBlock%vecC(i,j), j = 1, xgBlock%cols) /)
      end do
    end select
  end subroutine xgBlock_print
  !!***

  !!****f* m_xg/xgBlock_getid
  !!
  !! NAME
  !! xgBlock_getid

  function xgBlock_getid(xgBlock,comm) result (id)

    type(xgBlock_t), intent(in) :: xgBlock
    integer, intent(in),optional :: comm

    real(dp) :: id
    integer :: ierr,comm_

    if (xgBlock%gpu_option/=ABI_GPU_DISABLED) then
      call xgBlock_copy_from_gpu(xgBlock)
    end if
    select case(xgBlock%space)
      case (SPACE_R)
        id = sum(abs(xgBlock%vecR(:,:)))
      case (SPACE_CR)
        if (xgBlock%me_g0<0) then
          ABI_ERROR("xgBlock me_g0 is not initialized")
        end if
        id = 2*sum(abs(xgBlock%vecR(:,:)))
        if (xgBlock%me_g0==1) then
          id = id - sum(abs(xgBlock%vecR(1,:)))
        end if
      case (SPACE_C)
        id = sum(abs(dble(xgBlock%vecC(:,:))))+sum(abs(dimag(xgBlock%vecC(:,:))))
    end select
    comm_=xgBlock%spacedim_comm
    if (present(comm)) then
      comm_=comm
    end if
    if (xmpi_comm_size(comm_)>1) call xmpi_sum(id,comm_,ierr)

  end function xgBlock_getid
  !!***

  !!****f* m_xg/xgBlock_get_im_g0
  !!
  !! NAME
  !! xgBlock_get_im_g0

  function xgBlock_get_im_g0(xgBlock,comm) result (im_g0)

    type(xgBlock_t), intent(in) :: xgBlock
    integer, intent(in),optional :: comm

    real(dp) :: im_g0
    integer :: ierr,comm_

    if (xgBlock%gpu_option/=ABI_GPU_DISABLED) then
      call xgBlock_copy_from_gpu(xgBlock)
    end if
    select case(xgBlock%space)
      case (SPACE_R)
        im_g0 = zero
      case (SPACE_CR)
        if (xgBlock%me_g0<0) then
          ABI_ERROR("xgBlock me_g0 is not initialized")
        end if
        im_g0 = zero
        if (xgBlock%me_g0==1) then
          im_g0 = im_g0 + sum(abs(xgBlock%vecR(2,:)))
        end if
      case (SPACE_C)
        im_g0 = zero
    end select
    comm_=xgBlock%spacedim_comm
    if (present(comm)) then
      comm_=comm
    end if
    if (xmpi_comm_size(comm_)>1) call xmpi_sum(im_g0,comm_,ierr)

  end function xgBlock_get_im_g0
  !!***

  !!****f* m_xg/xg_finalize
  !!
  !! NAME
  !! xg_finalize

  subroutine xg_finalize()

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
