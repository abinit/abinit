!!****m* ABINIT/m_numeric_tools
!! NAME
!!  m_numeric_tools
!!
!! FUNCTION
!!  This module contains basic tools for numeric computations.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, GMR, MJV, XG, MVeithen, NH, FJ, MT, DCS, FrD, Olevano, Reining, Sottile, AL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_numeric_tools

 use defs_basis
 use m_errors
 use m_abicore
 use m_linalg_interfaces

 use m_fstrings,   only : itoa, sjoin

 implicit none

 private

 public :: arth                  ! Return an arithmetic progression
 public :: linspace              ! Similar to the above but with start, stop and num of division
 public :: geop                  ! Return a geometric progression
 public :: reverse               ! Reverse a 1D array *IN PLACE*
 public :: set2unit              ! Set the matrix to be a unit matrix (if it is square)
 public :: get_trace             ! Calculate the trace of a square matrix
 public :: get_diag              ! Return the diagonal of a matrix as a vector
 public :: isdiagmat             ! True if matrix is diagonal
 public :: l2int                 ! convert logical data to int array
 public :: r2c,c2r               ! Transfer complex data stored in a real array to a complex array and vice versa
 public :: iseven                ! True if int is even
 public :: isinteger             ! True if all elements of rr differ from an integer by less than tol
 public :: is_zero               ! True if all elements of rr differ from zero by less than tol
 public :: inrange               ! True if (int/float) is inside an interval.
 public :: bisect                ! Given a monotonic array A and x find j such that A(j)>x>A(j+1) using bisection
 public :: imax_loc              ! Index of maxloc on an array returned as scalar instead of array-valued quantity
 public :: imin_loc              ! Index of minloc on an array returned as scalar instead of array-valued quantity
 public :: lfind                 ! Find the index of the first occurrence of .True. in a logical array.
 public :: list2blocks           ! Given a list of integers, find the number of contiguos groups of values.
 public :: mask2blocks           ! Find groups of .TRUE. elements in a logical mask.
 public :: linfit                ! Perform a linear fit, y=ax+b, of data
 public :: llsfit_svd            ! Linear least squares fit with SVD of an user-defined set of functions
 public :: polyn_interp          ! Polynomial interpolation with Nevilles"s algorithms, error estimate is reported
 public :: quadrature            ! Driver routine for performing quadratures in finite domains using different algorithms
 public :: cspint                ! Estimates the integral of a tabulated function.
 public :: ctrap                 ! Corrected trapezoidal integral on uniform grid of spacing hh.
 public :: coeffs_gausslegint    ! Compute the coefficients (supports and weights) for Gauss-Legendre integration.
 public :: simpson_cplx          ! Integrate a complex function via extended Simpson's rule.
 public :: hermitianize          ! Force a square matrix to be hermitian
 public :: mkherm                ! Make the complex array(2,ndim,ndim) hermitian, by adding half of it to its hermitian conjugate.
 public :: hermit                ! Rdefine diagonal elements of packed matrix to impose Hermiticity.
 public :: symmetrize            ! Force a square matrix to be symmetric
 public :: pack_matrix           ! Packs a matrix into hermitian format
 public :: print_arr             ! Print a vector/array
 public :: pade, dpade           ! Functions for Pade approximation (complex case)
 public :: newrap_step           ! Apply single step Newton-Raphson method to find root of a complex function
 public :: OPERATOR(.x.)         ! Cross product of two 3D vectors
 public :: l2norm                ! Return the length (ordinary L2 norm) of a vector
 public :: remove_copies         ! Find the subset of inequivalent items in a list.
 public :: denominator           ! Return the denominator of a rational number.
 public :: mincm                 ! Return the minimum common multiple of two integers.
 public :: continued_fract       ! Routine to calculate the continued fraction (see description).
 public :: cmplx_sphcart         ! Convert an array of cplx numbers from spherical to Cartesian coordinates or vice versa.
 public :: pfactorize            ! Factorize a number in terms of an user-specified set of prime factors.
 public :: isordered             ! Check the ordering of a sequence.
 public :: wrap2_zero_one        ! Transforms a real number in a reduced number in the interval [0,1[ where 1 is not included (tol12)
 public :: wrap2_pmhalf          ! Transforms a real number in areduced number in the interval ]-1/2,1/2] where -1/2 is not included (tol12)
 public :: interpol3d            ! Linear interpolation in 3D
 public :: interpol3d_indices    ! Computes the indices in a cube which are neighbors to the point to be interpolated in interpol3d
 public :: interpolate_denpot    ! Liner interpolation of scalar field e.g. density of potential
 public :: simpson_int           ! Simpson integral of a tabulated function. Returns arrays with integrated values
 public :: simpson               ! Simpson integral of a tabulated function. Returns scalar with the integral on the full mesh.
 public :: rhophi                ! Compute the phase and the module of a complex number.
 public :: smooth                ! Smooth data.
 public :: nderiv                ! Compute first or second derivative of input function y(x) on a regular grid.
 public :: central_finite_diff   ! Coefficients of the central differences, for several orders of accuracy.
 public :: uniformrandom         ! Returns a uniform random deviate between 0.0 and 1.0.
 public :: findmin               ! Compute the minimum of a function whose value and derivative are known at two points.
 public :: kramerskronig         ! check or apply the Kramers Kronig relation
 public :: invcb                 ! Compute a set of inverse cubic roots as fast as possible.
 public :: safe_div              ! Performs 'save division' that is to prevent overflow, underflow, NaN or infinity errors
 public :: bool2index            ! Allocate and return array with the indices in the input boolean array that evaluates to .True.

 !MG FIXME: deprecated: just to avoid updating refs while refactoring.
 public :: dotproduct

 interface arth
   module procedure arth_int
   module procedure arth_rdp
 end interface arth

 interface reverse
   module procedure reverse_int
   module procedure reverse_rdp
 end interface reverse

 interface set2unit
   module procedure unit_matrix_int
   module procedure unit_matrix_rdp
   module procedure unit_matrix_cdp
 end interface set2unit

 interface get_trace
   module procedure get_trace_int
   module procedure get_trace_rdp
   module procedure get_trace_cdp
 end interface get_trace

 !interface cart_prod33
 !  module procedure cart_prod33_int
 !  module procedure cart_prod33_rdp
 !  module procedure cart_prod33_cdp
 !end interface cart_prod33

 interface get_diag
   module procedure get_diag_int
   module procedure get_diag_rdp
   module procedure get_diag_cdp
 end interface get_diag

 interface isdiagmat
   module procedure isdiagmat_int
   module procedure isdiagmat_rdp
   !module procedure isdiagmat_cdp
 end interface isdiagmat

 interface inrange
   module procedure inrange_int
   module procedure inrange_dp
 end interface inrange

 interface l2int
   module procedure l2int_1D
   module procedure l2int_2D
   module procedure l2int_3D
 end interface l2int

 interface r2c
   module procedure rdp2cdp_1D
   module procedure rdp2cdp_2D
   module procedure rdp2cdp_3D
   module procedure rdp2cdp_4D
   module procedure rdp2cdp_5D
   module procedure rdp2cdp_6D
 end interface r2c

 interface c2r
   module procedure cdp2rdp_1D
   module procedure cdp2rdp_2D
   module procedure cdp2rdp_3D
   module procedure cdp2rdp_4D
   module procedure cdp2rdp_5D
 end interface c2r

 interface isinteger
   module procedure is_integer_0d
   module procedure is_integer_1d
 end interface isinteger

 interface is_zero
   module procedure is_zero_rdp_0d
   module procedure is_zero_rdp_1d
 end interface is_zero

 interface bisect
   module procedure bisect_rdp
   module procedure bisect_int
 end interface bisect

 interface imax_loc
   module procedure imax_loc_int
   module procedure imax_loc_rdp
 end interface imax_loc

 interface imin_loc
   module procedure imin_loc_int
   module procedure imin_loc_rdp
 end interface imin_loc

 interface linfit
   module procedure linfit_rdp
   module procedure linfit_spc
   module procedure linfit_dpc
 end interface linfit

 interface hermitianize
   module procedure hermitianize_spc
   module procedure hermitianize_dpc
 end interface hermitianize

 interface symmetrize
   module procedure symmetrize_spc
   module procedure symmetrize_dpc
 end interface symmetrize

 interface print_arr  !TODO add prtm
   module procedure print_arr1d_spc
   module procedure print_arr1d_dpc
   module procedure print_arr2d_spc
   module procedure print_arr2d_dpc
 end interface print_arr

 interface operator (.x.)
   module procedure cross_product_int
   module procedure cross_product_rdp
 end interface

 interface l2norm
   module procedure l2norm_rdp
 end interface l2norm

 interface isordered
   module procedure isordered_rdp
 end interface isordered
!!***

!----------------------------------------------------------------------

!!****t* m_numeric_tools/stats_t
!! NAME
!! stats_t
!!
!! FUNCTION
!!  Statistical parameters of a data distribution.
!!
!! SOURCE

 type,public :: stats_t
   real(dp) :: mean
   real(dp) :: stdev
   real(dp) :: min
   real(dp) :: max
 end type stats_t

 public :: stats_eval  ! Calculate the statistical parameters of a data distribution.
!!***

!----------------------------------------------------------------------

!!****t* m_numeric_tools/vdiff_t
!! NAME
!! vdiff_t
!!
!! FUNCTION
!!  Estimate the "distance" between two functions tabulated on a homogeneous grid.
!!  Use `vidff` function to construct the object.
!!
!! SOURCE

 type,public :: vdiff_t

   real(dp) :: int_adiff = zero     ! \int |f1-f2| dr
   real(dp) :: mean_adiff = zero    ! Mean {|f1-f2|}
   real(dp) :: stdev_adiff = zero   ! Standard deviation of {|f1-f2|}
   real(dp) :: min_adiff = zero     ! Min {|f1-f2|}
   real(dp) :: max_adiff = zero     ! Max {|f1-f2|}
   real(dp) :: l1_rerr = zero       ! (\int |f1-f2| dr) / (\int |f2| dr)

 end type vdiff_t

 public :: vdiff_eval         ! Estimate the "distance" between two functions tabulated on a homogeneous grid.
 public :: vdiff_print        ! Print vdiff_t to formatted file.
!!***


CONTAINS  !===========================================================
!!***

!!****f* m_numeric_tools/arth_int
!! NAME
!!  arth_int
!!
!! FUNCTION
!!  Returns an array of length nn containing an arithmetic progression whose
!!  starting value is start and whose step is step.
!!
!! INPUTS
!!  start=initial point
!!  step=the increment
!!  nn=the number of points
!!
!! OUTPUT
!!  arth(nn)=the progression
!!
!! SOURCE

pure function arth_int(start, step, nn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 integer,intent(in) :: start,step
 integer :: arth_int(nn)

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 select case (nn)

 case (1:)
   arth_int(1)=start
   do ii=2,nn
    arth_int(ii)=arth_int(ii-1)+step
   end do

 case (0)
   return
 end select

end function arth_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/arth_rdp
!! NAME
!!  arth_rdp
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function arth_rdp(start, step, nn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp),intent(in) :: start,step
 real(dp) :: arth_rdp(nn)

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 select case (nn)
 case (1:)
   arth_rdp(1)=start
   do ii=2,nn
    arth_rdp(ii)=arth_rdp(ii-1)+step
   end do

 case (0)
   return
 end select

end function arth_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/linspace
!! NAME
!!  linspace
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function linspace(start, stop, nn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp),intent(in) :: start, stop
 real(dp) :: linspace(nn)

!Local variables-------------------------------
 real(dp) :: length
 integer :: ii
! *********************************************************************

 select case (nn)
 case (1:)
   length = stop-start
   do ii=1,nn
     linspace(ii)=start+length*(ii-1)/(nn-1)
   end do

 case (0)
   return
 end select

end function linspace
!!***


!----------------------------------------------------------------------

!!****f* m_numeric_tools/geop
!! NAME
!!  geop
!!
!! FUNCTION
!!  Returns an array of length nn containing a geometric progression whose
!!  starting value is start and whose factor is factor!
!!
!! INPUTS
!!  start=initial point
!!  factor=the factor of the geometric progression
!!  nn=the number of points
!!
!! OUTPUT
!!  geop(nn)=the progression
!!
!! SOURCE


pure function geop(start,factor,nn) result(res)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: start,factor
 integer,intent(in) :: nn
 real(dp) :: res(nn)

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 if (nn>0) res(1)=start
 do ii=2,nn
   res(ii)=res(ii-1)*factor
 end do

end function geop
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/reverse_int
!! NAME
!!  reverse_int
!!
!! FUNCTION
!!   Reverse a 1D array *IN PLACE*. Target: INT arrays
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine reverse_int(arr)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: arr(:)
!arrays
 integer :: ii,nn,swap
! *************************************************************************

 nn = SIZE(arr)
 if (nn <= 1) return

 do ii=1,nn/2
   swap = arr(ii)
   arr(ii) = arr(nn-ii+1)
   arr(nn-ii+1) = swap
 end do

end subroutine reverse_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/reverse_rdp
!! NAME
!!  reverse_rdp
!!
!! FUNCTION
!!   Reverse a 1D array *IN PLACE*. Target: DP arrays
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine reverse_rdp(arr)

!Arguments ------------------------------------
!scalars
 real(dp),intent(inout) :: arr(:)
!arrays
 integer :: ii,nn
 real(dp) :: swap
! *************************************************************************

 nn = SIZE(arr)
 if (nn <= 1) return

 do ii=1,nn/2
   swap = arr(ii)
   arr(ii) = arr(nn-ii+1)
   arr(nn-ii+1) = swap
 end do

end subroutine reverse_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/unit_matrix_int
!! NAME
!!  unit_matrix_int
!!
!! FUNCTION
!!  Set the matrix matrix to be a unit matrix (if it is square).
!!
!! SIDE EFFECTS
!!  matrix(:,:)=set to unit on exit
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine unit_matrix_int(matrix)

!Arguments ------------------------------------
 integer,intent(inout) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
! *********************************************************************

 nn=MIN(SIZE(matrix,DIM=1),SIZE(matrix,DIM=2))
 matrix(:,:)=0
 do ii=1,nn
  matrix(ii,ii)=1
 end do

end subroutine unit_matrix_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/unit_matrix_rdp
!! NAME
!!  unit_matrix_rdp
!!
!! FUNCTION
!!  Set the matrix matrix to be a unit matrix (if it is square).
!!
!! SIDE EFFECTS
!!  matrix(:,:)=set to unit on exit
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine unit_matrix_rdp(matrix)

!Arguments ------------------------------------
 real(dp),intent(inout) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
! *********************************************************************

 nn=MIN(SIZE(matrix,DIM=1),SIZE(matrix,DIM=2))
 matrix(:,:)=zero
 do ii=1,nn
   matrix(ii,ii)=one
 end do

end subroutine unit_matrix_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/unit_matrix_cdp
!! NAME
!!  unit_matrix_cdp
!!
!! FUNCTION
!!  Set the matrix matrix to be a unit matrix (if it is square).
!!
!! SIDE EFFECTS
!!  matrix(:,:)=set to unit on exit
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine unit_matrix_cdp(matrix)

!Arguments ------------------------------------
 complex(dpc),intent(inout) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
! *********************************************************************

 nn=MIN(SIZE(matrix,DIM=1),SIZE(matrix,DIM=2))
 matrix=czero
 do ii=1,nn
   matrix(ii,ii)=cone
 end do

end subroutine unit_matrix_cdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/get_trace_int
!! NAME
!!  get_trace_int
!!
!! FUNCTION
!!  Calculate the trace of a square matrix
!!
!! INPUTS
!!  matrix(:,:)
!!
!! OUTPUT
!!  trace=the trace
!!
!! SOURCE

pure function get_trace_int(matrix) result(trace)

!Arguments ------------------------------------
 integer :: trace
 integer,intent(in) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii
! *********************************************************************

 trace=0
 do ii=1,size(matrix,dim=1)
   trace=trace+matrix(ii,ii)
 end do

end function get_trace_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/get_trace_rdp
!! NAME
!!  get_trace_int
!!
!! FUNCTION
!!  Calculate the trace of a square matrix (real(dp) version)
!!
!! INPUTS
!!  matrix(:,:)
!!
!! OUTPUT
!!  trace=the trace
!!
!! SOURCE

pure function get_trace_rdp(matrix) result(trace)

!Arguments ------------------------------------
 real(dp) :: trace
 real(dp),intent(in) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii
! *********************************************************************

 trace=zero
 do ii=1,size(matrix,dim=1)
    trace=trace+matrix(ii,ii)
 end do

end function get_trace_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/get_trace_cdp
!! NAME
!!  get_trace_cdp
!!
!! FUNCTION
!!  Calculate the trace of a square matrix (complex(dpc) version)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function get_trace_cdp(matrix) result(trace)

!Arguments ------------------------------------
 complex(dpc) :: trace
 complex(dpc),intent(in) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii
! *********************************************************************

 trace=czero
 do ii=1,size(matrix,dim=1)
    trace=trace+matrix(ii,ii)
 end do

end function get_trace_cdp
!!***

!!****f* m_numeric_tools/get_diag_int
!! NAME
!!  get_diag_int
!!
!! FUNCTION
!!  Return the diagonal of a square matrix as a vector
!!
!! INPUTS
!!  matrix(:,:)
!!
!! OUTPUT
!!  diag(:)=the diagonal
!!
!! SOURCE

function get_diag_int(mat) result(diag)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mat(:,:)
 integer :: diag(SIZE(mat,1))

!Local variables-------------------------------
 integer :: ii
! *************************************************************************

 ii = assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',__FILE__,__LINE__)

 do ii=1,SIZE(mat,1)
   diag(ii)=mat(ii,ii)
 end do

end function get_diag_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/get_diag_rdp
!! NAME
!!  get_diag_rdp
!!
!! FUNCTION
!!  Return the diagonal of a square matrix as a vector
!!
!! INPUTS
!!  matrix(:,:)
!!
!! OUTPUT
!!  diag(:)=the diagonal
!!
!! SOURCE

function get_diag_rdp(mat) result(diag)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: mat(:,:)
 real(dp) :: diag(SIZE(mat,1))

!Local variables-------------------------------
 integer :: ii
! *************************************************************************

 ABI_CHECK(SIZE(mat,1) == SIZE(mat,2), 'Matrix not square')

 do ii=1,SIZE(mat,1)
   diag(ii) = mat(ii,ii)
 end do

end function get_diag_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/get_diag_cdp
!! NAME
!!  get_diag_cdp
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function get_diag_cdp(cmat) result(cdiag)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: cmat(:,:)
 complex(dpc) :: cdiag(SIZE(cmat,1))

!Local variables-------------------------------
 integer :: ii
! *************************************************************************

 ABI_CHECK(SIZE(cmat,1) == SIZE(cmat,2), 'Matrix not square')

 do ii=1,SIZE(cmat,1)
   cdiag(ii)=cmat(ii,ii)
 end do

end function get_diag_cdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/isdiagmat_int
!! NAME
!!  isdiagmat_int
!!
!! FUNCTION
!!  True if matrix mat is diagonal
!!
!! SOURCE

pure logical function isdiagmat_int(mat) result(ans)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mat(:,:)

!Local variables-------------------------------
 integer :: ii,jj
! *************************************************************************

 ans = .True.
 do jj=1,size(mat,dim=2)
   do ii=1,size(mat,dim=1)
     if (ii == jj) cycle
     if (mat(ii,jj) /= 0) then
       ans = .False.; return
     end if
   end do
 end do

end function isdiagmat_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/isdiagmat_rdp
!! NAME
!!  isdiagmat_rdp
!!
!! FUNCTION
!!  True if matrix mat is diagonal withing the given absolute tolerance (default: tol12)
!!
!! SOURCE

pure logical function isdiagmat_rdp(mat, atol) result(ans)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: mat(:,:)
 real(dp),optional,intent(in) :: atol

!Local variables-------------------------------
 integer :: ii,jj
 real(dp) :: my_atol
! *************************************************************************

 my_atol = tol12; if (present(atol)) my_atol = atol

 ans = .True.
 do jj=1,size(mat,dim=2)
   do ii=1,size(mat,dim=1)
     if (ii == jj) cycle
     if (abs(mat(ii,jj)) > my_atol) then
       ans = .False.; return
     end if
   end do
 end do

end function isdiagmat_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/l2int_1D
!! NAME
!!  l2int_1D
!!
!! FUNCTION
!!  Convert a logical array into an int array (True --> 1, False --> 0)
!!
!! INPUTS
!!  larr(:)=the input logical array
!!
!! SOURCE

pure function l2int_1D(larr) result(int_arr)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: larr(:)
 integer :: int_arr(size(larr))

! *********************************************************************

 where (larr)
   int_arr = 1
 elsewhere
   int_arr = 0
 end where

end function l2int_1D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/l2int_2D
!! NAME
!!  l2int_2D
!!
!! FUNCTION
!!  Convert a logical array into an int array (True --> 1, False --> 0)
!!
!! INPUTS
!!  larr(:)=the input logical array
!!
!! SOURCE

pure function l2int_2D(larr) result(int_arr)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: larr(:,:)
 integer :: int_arr(size(larr,1), size(larr,2))

! *********************************************************************

 where (larr)
   int_arr = 1
 elsewhere
   int_arr = 0
 end where

end function l2int_2D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/l2int_3D
!! NAME
!!  l2int_3D
!!
!! FUNCTION
!!  Convert a logical array into an int array (True --> 1, False --> 0)
!!
!! INPUTS
!!  larr(:)=the input logical array
!!
!! SOURCE

pure function l2int_3D(larr) result(int_arr)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: larr(:,:,:)
 integer :: int_arr(size(larr,1), size(larr,2), size(larr,3))

! *********************************************************************

 where (larr)
   int_arr = 1
 elsewhere
   int_arr = 0
 end where

end function l2int_3D
!!***

!----------------------------------------------------------------------

!!***!!****f* m_numeric_tools/rdp2cdp_1D
!! NAME
!!  rdp2cdp_1D
!!
!! FUNCTION
!!  Create a complex array starting from a real array containing real and imaginary part
!!
!! INPUTS
!!  rr(:)=the real array
!!
!! OUTPUT
!!  cc(:)=the complex array
!!
!! SOURCE

pure function rdp2cdp_1D(rr) result(cc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rr(:,:)
 complex(dpc) :: cc(SIZE(rr,2))

! *********************************************************************

 cc(:)=CMPLX(rr(1,:),rr(2,:),kind=dpc)

end function rdp2cdp_1D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/rdp2cdp_2D
!! NAME
!!  rdp2cdp_2D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function rdp2cdp_2D(rr) result(cc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rr(:,:,:)
 complex(dpc) :: cc(SIZE(rr,2),SIZE(rr,3))

! *********************************************************************

 cc(:,:)=CMPLX(rr(1,:,:),rr(2,:,:), kind=dpc)

end function rdp2cdp_2D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/rdp2cdp_3D
!! NAME
!!  rdp2cdp_3D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function rdp2cdp_3D(rr) result(cc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rr(:,:,:,:)
 complex(dpc) :: cc(SIZE(rr,2),SIZE(rr,3),SIZE(rr,4))

! *********************************************************************

 cc(:,:,:)=CMPLX(rr(1,:,:,:),rr(2,:,:,:), kind=dpc)

end function rdp2cdp_3D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/rdp2cdp_4D
!! NAME
!!  rdp2cdp_4D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function rdp2cdp_4D(rr) result(cc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rr(:,:,:,:,:)
 complex(dpc) :: cc(SIZE(rr,2),SIZE(rr,3),SIZE(rr,4),SIZE(rr,5))

! *********************************************************************

 cc(:,:,:,:)=CMPLX(rr(1,:,:,:,:),rr(2,:,:,:,:), kind=dpc)

end function rdp2cdp_4D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/rdp2cdp_5D
!! NAME
!!  rdp2cdp_5D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function rdp2cdp_5D(rr) result(cc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rr(:,:,:,:,:,:)
 complex(dpc) :: cc(SIZE(rr,2),SIZE(rr,3),SIZE(rr,4),SIZE(rr,5),SIZE(rr,6))

! *********************************************************************

 cc(:,:,:,:,:)=CMPLX(rr(1,:,:,:,:,:),rr(2,:,:,:,:,:), kind=dpc)

end function rdp2cdp_5D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/rdp2cdp_6D
!! NAME
!!  rdp2cdp_6D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function rdp2cdp_6D(rr) result(cc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rr(:,:,:,:,:,:,:)
 complex(dpc) :: cc(SIZE(rr,2),SIZE(rr,3),SIZE(rr,4),SIZE(rr,5),SIZE(rr,6),SIZE(rr,7))

! *********************************************************************

 cc(:,:,:,:,:,:)=CMPLX(rr(1,:,:,:,:,:,:),rr(2,:,:,:,:,:,:), kind=dpc)

end function rdp2cdp_6D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cdp2rdp_1D
!! NAME
!!  cdp2rdp_1D
!!
!! FUNCTION
!!  Create a real array containing real and imaginary part starting from a complex array
!!
!! INPUTS
!!  cc(:)=the input complex array
!!
!! OUTPUT
!!  rr(2,:)=the real array
!!
!! SOURCE

pure function cdp2rdp_1D(cc) result(rr)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: cc(:)
 real(dp) :: rr(2,SIZE(cc))

! *********************************************************************

 rr(1,:)=REAL (cc(:))
 rr(2,:)=AIMAG(cc(:))

end function cdp2rdp_1D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cdp2rdp_2D
!! NAME
!!  cdp2rdp_2D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function cdp2rdp_2D(cc) result(rr)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: cc(:,:)
 real(dp) :: rr(2,SIZE(cc,1),SIZE(cc,2))
! *********************************************************************

 rr(1,:,:)=REAL (cc(:,:))
 rr(2,:,:)=AIMAG(cc(:,:))

end function cdp2rdp_2D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cdp2rdp_3D
!! NAME
!!  cdp2rdp_3D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function cdp2rdp_3D(cc) result(rr)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: cc(:,:,:)
 real(dp) :: rr(2,SIZE(cc,1),SIZE(cc,2),SIZE(cc,3))

! *********************************************************************

 rr(1,:,:,:)=REAL (cc(:,:,:))
 rr(2,:,:,:)=AIMAG(cc(:,:,:))

end function cdp2rdp_3D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cdp2rdp_4D
!! NAME
!!  cdp2rdp_4D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function cdp2rdp_4D(cc) result(rr)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: cc(:,:,:,:)
 real(dp) :: rr(2,SIZE(cc,1),SIZE(cc,2),SIZE(cc,3),SIZE(cc,4))
! *********************************************************************

 rr(1,:,:,:,:)=REAL (cc(:,:,:,:))
 rr(2,:,:,:,:)=AIMAG(cc(:,:,:,:))

end function cdp2rdp_4D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cdp2rdp_5D
!! NAME
!!  cdp2rdp_5D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function cdp2rdp_5D(cc) result(rr)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: cc(:,:,:,:,:)
 real(dp) :: rr(2,SIZE(cc,1),SIZE(cc,2),SIZE(cc,3),SIZE(cc,4),SIZE(cc,5))

! *********************************************************************

 rr(1,:,:,:,:,:)=REAL (cc(:,:,:,:,:))
 rr(2,:,:,:,:,:)=AIMAG(cc(:,:,:,:,:))

end function cdp2rdp_5D
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/iseven
!! NAME
!!  iseven
!!
!! FUNCTION
!!  Return .TRUE. if the given integer is even
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function iseven(nn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 logical :: iseven
! *********************************************************************

 iseven=.FALSE.; if ((nn/2)*2==nn) iseven=.TRUE.

end function iseven
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/is_integer_0D
!! NAME
!!  is_integer_0D
!!
!! FUNCTION
!!  Return .TRUE. if all elements differ from an integer by less that tol
!!
!! INPUTS
!!  rr=the set of real values to be checked
!!  tol=tolerance on the difference between real and integer
!!
!! SOURCE

pure function is_integer_0d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr

! *************************************************************************

 ans=(ABS(rr-NINT(rr))<tol)

end function is_integer_0d
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/is_integer_1D
!! NAME
!!  is_integer_1D
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function is_integer_1d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr(:)

! *************************************************************************

 ans=ALL((ABS(rr-NINT(rr))<tol))

end function is_integer_1d
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/is_zero_rdp_0D
!! NAME
!!  is_zero_rdp_0D
!!
!! FUNCTION
!!  Return .TRUE. if all elements differ from zero by less that tol
!!
!! INPUTS
!!  rr=the set of real values to be checked
!!  tol=tolerance
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function is_zero_rdp_0d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr
! *************************************************************************

 ans=(ABS(rr)<tol)

end function is_zero_rdp_0d
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/is_zero_rdp_1d
!! NAME
!!  is_zero_rdp_1d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function is_zero_rdp_1d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr(:)
! *************************************************************************

 ans=ALL(ABS(rr(:))<tol)

end function is_zero_rdp_1d
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/inrange_int
!! NAME
!!  inrange_int
!!
!! FUNCTION
!!  True if int `xval` is inside the interval [win(1), win(2)]
!!
!! SOURCE

pure logical function inrange_int(xval, win)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: xval,win(2)
! *************************************************************************

 inrange_int = (xval >= win(1) .and. xval <= win(2))

end function inrange_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/inrange_dp
!! NAME
!!  inrange_dp
!!
!! FUNCTION
!!  True if float `xval` is inside the interval [win(1), win(2)]
!!
!! SOURCE

pure logical function inrange_dp(xval, win)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: xval, win(2)
! *************************************************************************

 inrange_dp = (xval >= win(1) .and. xval <= win(2))

end function inrange_dp
!!***

!!****f* m_numeric_tools/bisect_rdp
!! NAME
!!  bisect_rdp
!!
!! FUNCTION
!!  Given an array AA(1:N), and a value x, returns the index j such that AA(j)<=x<= AA(j + 1).
!!  AA must be monotonic, either increasing or decreasing. j=0 or
!!  j=N is returned to indicate that x is out of range.
!!
!! SOURCE

pure function bisect_rdp(AA,xx) result(loc)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: AA(:)
 real(dp),intent(in) :: xx
 integer :: loc

!Local variables-------------------------------
 integer :: nn,jl,jm,ju
 logical :: ascnd
! *********************************************************************

 nn=SIZE(AA); ascnd=(AA(nn)>=AA(1))
 !
 ! Initialize lower and upper limits
 jl=0; ju=nn+1
 do
   if (ju-jl<=1) EXIT
   jm=(ju+jl)/2  ! Compute a midpoint,
   if (ascnd.EQV.(xx>=AA(jm))) then
     jl=jm ! Replace lower limit
   else
     ju=jm ! Replace upper limit
   end if
 end do
 !
 ! Set the output, being careful with the endpoints
 if (xx==AA(1)) then
   loc=1
 else if (xx==AA(nn)) then
   loc=nn-1
 else
   loc=jl
 end if

end function bisect_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/bisect_int
!! NAME
!!  bisect_int
!!
!! FUNCTION
!!  Given an array AA(1:N), and a value x, returns the index j such that AA(j)<=x<= AA(j + 1).
!!  AA must be monotonic, either increasing or decreasing. j=0 or
!!  j=N is returned to indicate that x is out of range.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
pure function bisect_int(AA,xx) result(loc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: AA(:)
 integer,intent(in) :: xx
 integer :: loc

!Local variables-------------------------------
 integer :: nn,jl,jm,ju
 logical :: ascnd
! *********************************************************************

 nn=SIZE(AA) ; ascnd=(AA(nn)>=AA(1))

 ! Initialize lower and upper limits
 jl=0 ; ju=nn+1
 do
  if (ju-jl<=1) EXIT
  jm=(ju+jl)/2  ! Compute a midpoint
  if (ascnd.EQV.(xx>=AA(jm))) then
   jl=jm ! Replace lower limit
  else
   ju=jm ! Replace upper limit
  end if
 end do
 !
 ! Set the output, being careful with the endpoints
 if (xx==AA(1)) then
  loc=1
 else if (xx==AA(nn)) then
  loc=nn-1
 else
  loc=jl
 end if

end function bisect_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/imax_loc_int
!! NAME
!!  imax_loc_int
!!
!! FUNCTION
!!  Index of maxloc on an array returned as scalar instead of array-valued
!!
!! SOURCE

pure function imax_loc_int(iarr,mask)

!Arguments ------------------------------------
!scalars
 integer :: imax_loc_int
!arrays
 integer,intent(in) :: iarr(:)
 logical,optional,intent(in) :: mask(:)

!Local variables-------------------------------
 integer :: imax(1)
! *************************************************************************

 if (PRESENT(mask)) then
  imax=MAXLOC(iarr,MASK=mask)
 else
  imax=MAXLOC(iarr)
 end if
 imax_loc_int=imax(1)

end function imax_loc_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/imax_loc_rdp
!! NAME
!!  imax_loc_rdp
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
pure function imax_loc_rdp(arr,mask)

!Arguments ------------------------------------
!scalars
 integer :: imax_loc_rdp
!arrays
 real(dp),intent(in) :: arr(:)
 logical,optional,intent(in) :: mask(:)

!Local variables-------------------------------
 integer :: imax(1)
! *************************************************************************

 if (PRESENT(mask)) then
  imax=MAXLOC(arr,MASK=mask)
 else
  imax=MAXLOC(arr)
 end if
 imax_loc_rdp=imax(1)

end function imax_loc_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/imin_loc_int
!! NAME
!!  imin_loc_int
!!
!! FUNCTION
!!  Index of minloc on an array returned as scalar instead of array-valued
!!
!! SOURCE

pure function imin_loc_int(arr,mask)

!Arguments ------------------------------------
!scalars
 integer :: imin_loc_int
!arrays
 integer,intent(in) :: arr(:)
 logical,optional,intent(in) :: mask(:)

!Local variables-------------------------------
 integer :: imin(1)
! *************************************************************************

 if (PRESENT(mask)) then
  imin=MINLOC(arr,MASK=mask)
 else
  imin=MINLOC(arr)
 end if
 imin_loc_int=imin(1)

end function imin_loc_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/imin_loc_rdp
!! NAME
!!  imin_loc_rdp
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function imin_loc_rdp(arr,mask)

!Arguments ------------------------------------
!scalars
 integer :: imin_loc_rdp
!arrays
 real(dp),intent(in) :: arr(:)
 logical,optional,intent(in) :: mask(:)

!Local variables-------------------------------
 integer :: imin(1)
! *************************************************************************

 if (PRESENT(mask)) then
  imin=MINLOC(arr,MASK=mask)
 else
  imin=MINLOC(arr)
 end if

 imin_loc_rdp=imin(1)

end function imin_loc_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/lfind
!! NAME
!!  lfind
!!
!! FUNCTION
!!  Find the index of the first occurrence of .True. in a logical array.
!!  Return -1 if not found. If back is True, the search starts from the
!!  last element of the array (default: False).
!!
!! INPUTS
!!  mask(:)=Input logical mask
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function lfind(mask, back)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: mask(:)
 logical,optional,intent(in) :: back
!arrays

!Local variables-------------------------------
!scalars
 integer :: ii,nitems
 logical :: do_back

!************************************************************************

 do_back = .False.; if (present(back)) do_back = back
 lfind = -1; nitems = size(mask); if (nitems == 0) return

 if (do_back) then
   ! Backward search
   do ii=nitems,1,-1
     if (mask(ii)) then
       lfind = ii; return
     end if
   end do
 else
   ! Forward search.
   do ii=1,nitems
     if (mask(ii)) then
       lfind = ii; return
     end if
   end do
 end if

end function lfind
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/list2blocks
!! NAME
!!  list2blocks
!!
!! FUNCTION
!!  Given a list of integers, find the number of contiguos groups of values.
!!  and returns the set of indices that can be used to loop over these groups
!!  Example list = [1,2,3,5,6] --> blocks = [[1,3], [4,5]]
!!
!! INPUTS
!!  list(:)=List of integers
!!
!! OUTPUTS
!! nblocks=Number of blocks
!! blocks(2,nblocks)=
!!    allocatable array in input
!!    in output:
!!       blocks(1,i) gives the start of the i-th block
!!       blocks(2,i) gives the end of the i-th block
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine list2blocks(list,nblocks,blocks)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: nblocks
 integer,intent(in) :: list(:)
!arrays
 integer,intent(out),allocatable :: blocks(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nitems
!arrays
 integer :: work(2,size(list))

!************************************************************************

 nitems = size(list)

 ! Handle nitems == 1 case
 if (nitems == 1) then
   ABI_MALLOC(blocks, (2,1))
   blocks = 1
   return
 end if

 nblocks = 1; work(1,1) = 1

 do ii=2,nitems
   if (list(ii) /= (list(ii-1) + 1)) then
     work(2,nblocks) = ii - 1
     nblocks = nblocks + 1
     work(1,nblocks) = ii
   end if
 end do

 work(2,nblocks) = nitems

 ABI_MALLOC(blocks, (2,nblocks))
 blocks = work(:,1:nblocks)

end subroutine list2blocks
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/mask2blocks
!! NAME
!!  mask2blocks
!!
!! FUNCTION
!!  Give a logical mask, find the number of contiguos groups of .TRUE. values.
!!  and return the set of indices that can be used to loop over these groups
!!
!! INPUTS
!!  mask(:)=Input logical mask
!!
!! OUTPUTS
!! nblocks=Number of blocks
!!
!! SIDE EFFECTS
!! blocks(:,:)= Null pointer in input. blocks(2,nblocks) in output where
!!   blocks(1,i) gives the start of the i-th block
!!   blocks(2,i) gives the end of the i-th block
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine mask2blocks(mask,nblocks,blocks)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: nblocks
 logical,intent(in) :: mask(:)
!arrays
 integer,allocatable :: blocks(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nitems,start
 logical :: inblock
!arrays
 integer :: work(2,SIZE(mask))

!************************************************************************

 ! Find first element.
 nitems = size(mask); start = 0
 do ii=1,nitems
   if (mask(ii)) then
     start = ii
     exit
   end if
 end do

 ! Handle no true element or just one.
 if (start == 0) then
   nblocks = 0
   ABI_MALLOC(blocks, (0,0))
   return
 end if
 if (start /= 0 .and. nitems == 1) then
   nblocks = 1
   ABI_MALLOC(blocks, (2,1))
   blocks(:,1) = [1,1]
 end if

 nblocks = 1; work(1,1) = start; inblock = .True.

 do ii=start+1,nitems
   if (.not.mask(ii)) then
     if (inblock) then
       inblock = .False.
       work(2,nblocks) = ii - 1
     end if
   else
     if (.not. inblock) then
       inblock = .True.
       nblocks = nblocks + 1
       work(1,nblocks) = ii
     end if
   end if
 end do

 if (mask(nitems) .and. inblock) work(2,nblocks) = nitems

 ABI_MALLOC(blocks, (2,nblocks))
 blocks = work(:,1:nblocks)

end subroutine mask2blocks
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/linfit_rdp
!! NAME
!!  linfit_rdp
!!
!! FUNCTION
!!  Perform a linear fit, y=ax+b, of data
!!
!! INPUTS
!!  xx(nn)=xx coordinates
!!  yy(nn)=yy coordinates
!!
!! OUTPUT
!!  aa=coefficient of linear term of fit
!!  bb=coefficient of constant term of fit
!!  function linfit=root mean square of differences between data and fit
!!
!! SOURCE

function linfit_rdp(nn,xx,yy,aa,bb) result(res)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: res
 real(dp),intent(out) :: aa,bb
!arrays
 real(dp),intent(in) :: xx(nn),yy(nn)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: msrt,sx2,sx,sxy,sy,tx,ty
! *************************************************************************

 sx=zero ; sy=zero ; sxy=zero ; sx2=zero
 do ii=1,nn
  tx=xx(ii)
  ty=yy(ii)
  sx=sx+tx
  sy=sy+ty
  sxy=sxy+tx*ty
  sx2=sx2+tx*tx
 end do

 aa=(nn*sxy-sx*sy)/(nn*sx2-sx*sx)
 bb=sy/nn-sx*aa/nn

 msrt=zero
 do ii=1,nn
  tx=xx(ii)
  ty=yy(ii)
  msrt=msrt+(ty-aa*tx-bb)**2
 end do
 msrt=SQRT(msrt/nn) ; res=msrt

end function linfit_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/linfit_spc
!! NAME
!!  linfit_spc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function linfit_spc(nn,xx,zz,aa,bb) result(res)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: res
 real(dp),intent(in) :: xx(nn)
 complex(spc),intent(in) :: zz(nn)
 complex(spc),intent(out) :: aa,bb
!arrays

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: sx,sx2,msrt
 complex(dpc) :: sz,sxz
! *************************************************************************

 sx=zero ; sx2=zero ; msrt=zero
 sz=czero ; sxz=czero
 do ii=1,nn
  sx=sx+xx(ii)
  sz=sz+zz(ii)
  sxz=sxz+xx(ii)*zz(ii)
  sx2=sx2+xx(ii)*xx(ii)
 end do

 aa=CMPLX((nn*sxz-sx*sz)/(nn*sx2-sx*sx), kind=spc)
 bb=CMPLX(sz/nn-sx*aa/nn, kind=spc)

 do ii=1,nn
  msrt=msrt+ABS(zz(ii)-aa*xx(ii)-bb)**2
 end do
 msrt=SQRT(msrt) ; res=msrt

end function linfit_spc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/linfit_dpc
!! NAME
!!  linfit_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function linfit_dpc(nn,xx,zz,aa,bb) result(res)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: res
 real(dp),intent(in) :: xx(nn)
 complex(dpc),intent(in) :: zz(nn)
 complex(dpc),intent(out) :: aa,bb
!arrays

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: sx,sx2,msrt
 complex(dpc) :: sz,sxz
! *************************************************************************

 sx=zero  ; sx2=zero ; msrt=zero
 sz=czero ; sxz=czero
 do ii=1,nn
  sx=sx+xx(ii)
  sz=sz+zz(ii)
  sxz=sxz+xx(ii)*zz(ii)
  sx2=sx2+xx(ii)*xx(ii)
 end do

 aa=(nn*sxz-sx*sz)/(nn*sx2-sx*sx)
 bb=sz/nn-sx*aa/nn

 do ii=1,nn
  msrt=msrt+ABS(zz(ii)-aa*xx(ii)-bb)**2
 end do
 msrt=SQRT(msrt) ; res=msrt

end function linfit_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/llsfit_svd
!! NAME
!!  llsfit_svd
!!
!! FUNCTION
!!  Given a set of N data points (x,y) with individual standard deviations sigma_i,
!!  use chi-square minimization to determine the M coefficients, par, of a function that
!!  depends linearly on nfuncs functions, i.e f(x) = \sum_i^{nfuncs} par_i * func_i(x).
!!  Solve the fitting equations using singular value decomposition of the design matrix as in Eq 14.3.17
!!  of Numerical Recipes. The program returns values for the M fit parameters par, and chi-square.
!!  The user supplies a subroutine funcs(x,nfuncs) that returns the M basis functions evaluated at xx.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine llsfit_svd(xx,yy,sigma,nfuncs,funcs,chisq,par,var,cov,info)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfuncs
 integer,intent(out) :: info
 real(dp),intent(out) :: chisq
!arrays
 real(dp),intent(in) :: xx(:),yy(:),sigma(:)
 real(dp),intent(out) :: par(:),var(:),cov(:,:)

 interface
  function funcs(xx,nf)
  use defs_basis
  implicit none
  real(dp),intent(in) :: xx
  integer,intent(in) :: nf
  real(dp) :: funcs(nf)
  end function funcs
 end interface

!Local variables-------------------------------
 integer,parameter :: PAD_=50
 integer :: ii,npts,lwork
 real(dp),parameter :: TOL_=1.0e-5_dp
!arrays
 real(dp),dimension(SIZE(xx)) :: bb,sigm1
 real(dp),dimension(SIZE(xx),nfuncs) :: dmat,dmat_save
 real(dp) :: tmp(nfuncs)
 real(dp),allocatable :: work(:),Vt(:,:),U(:,:),S(:)
! *************************************************************************

 npts = assert_eq(SIZE(xx),SIZE(yy),SIZE(sigma),'Wrong size in xx,yy,sigma', __FILE__, __LINE__)
 call assert((npts>=nfuncs),'No. of functions must greater than no. of points', __FILE__, __LINE__)
 ii = assert_eq(nfuncs,SIZE(cov,1),SIZE(cov,2),SIZE(var),'Wrong size in covariance', __FILE__, __LINE__)

 !
 ! === Calculate design matrix and b vector ===
 ! * dmat_ij=f_j(x_i)/sigma_i, b_i=y_i/sigma_i
 sigm1(:)=one/sigma(:) ; bb(:)=yy(:)*sigm1(:)
 do ii=1,npts
  dmat_save(ii,:)=funcs(xx(ii),nfuncs)
 end do
 dmat=dmat_save*SPREAD(sigm1,DIM=2,ncopies=nfuncs)
 dmat_save(:,:)=dmat(:,:)
 !
 ! === Singular value decomposition ===
 lwork=MAX(3*MIN(npts,nfuncs)+MAX(npts,nfuncs),5*MIN(npts,nfuncs)-4)+PAD_
 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(U,(npts,npts))
 ABI_MALLOC(S,(nfuncs))
 ABI_MALLOC(Vt,(nfuncs,nfuncs))

 call DGESVD('A','A',npts,nfuncs,dmat,npts,S,U,npts,Vt,nfuncs,work,lwork,info)
 ABI_FREE(work)
 GOTO 10
 !
 ! === Set to zero small singular values according to TOL_ and find coefficients ===
 WHERE (S>TOL_*MAXVAL(S))
  tmp=MATMUL(bb,U)/S
 ELSEWHERE
  S  =zero
  tmp=zero
 END WHERE
 par(:)=MATMUL(tmp,Vt)
 !
 ! === Evaluate chi-square ===
 chisq=l2norm(MATMUL(dmat_save,par)-bb)**2
 !
 ! === Calculate covariance and variance ===
 ! C_jk = V_ji V_ki / S_i^2
 WHERE (S/=zero) S=one/(S*S)

 ! check this but should be correct
 cov(:,:)=Vt*SPREAD(S,DIM=2,ncopies=nfuncs)
 cov(:,:)=MATMUL(TRANSPOSE(Vt),cov)
 var(:)=SQRT(get_diag(cov))

10 continue
 ABI_FREE(U)
 ABI_FREE(S)
 ABI_FREE(Vt)

end subroutine llsfit_svd
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/polyn_interp
!! NAME
!!  polyn_interp
!!
!! FUNCTION
!!  Given arrays xa and ya of length N, and given a value x, return a value y, and an error estimate dy.
!!  If P(x) is the polynomial of degree N-1 such that P(xai)=yai, i=1,...,N, then the returned value y=P(x).
!!
!! INPUTS
!!  xa(:)=abscissas in ascending order
!!  ya(:)=ordinates
!!  x=the point where the set of data has to be interpolated
!!
!! OUTPUT
!!  y=the interpolated value
!!  dy=error estimate
!!
!! NOTES
!!  Based on the polint routine reported in Numerical Recipies
!!
!! PARENTS
!!      m_numeric_tools
!!
!! CHILDREN
!!
!! SOURCE

subroutine polyn_interp(xa,ya,x,y,dy)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: xa(:),ya(:)
 real(dp),intent(in) :: x
 real(dp),intent(out) :: y,dy
!Local variables-------------------------------
!scalars
 integer :: m,n,ns
!arrays
 real(dp),dimension(SIZE(xa)) :: c,d,den,ho
! *************************************************************************

 n = assert_eq(SIZE(xa),SIZE(ya),'Different size in xa and ya',__FILE__,__LINE__)

 ! === Initialize the tables of c and d ===
 c(:)=ya(:) ; d(:)=ya(:) ; ho(:)=xa(:)-x
 ! === Find closest table entry and initial approximation to y ===
 ns=imin_loc(ABS(x-xa)) ; y=ya(ns)
 ns=ns-1
 !
 ! === For each column of the tableau loop over current c and d and up-date them ===
 do m=1,n-1
  den(1:n-m)=ho(1:n-m)-ho(1+m:n)
  if (ANY(den(1:n-m)==zero)) then
   MSG_ERROR('Two input xa are identical')
  end if

  den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
  d(1:n-m)=ho(1+m:n)*den(1:n-m) ! Update c and d
  c(1:n-m)=ho(1:n-m)*den(1:n-m)

  if (2*ns<n-m) then  ! Now decide which correction, c or d, we want to add to the
   dy=c(ns+1)         ! accumulating value of y, The last dy added is the error indication.
  else
   dy=d(ns)
   ns=ns-1
  end if

  y=y+dy
 end do

end subroutine polyn_interp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/trapezoidal_
!! NAME
!!  trapezoidal_ (PRIVATE)
!!
!! FUNCTION
!!  Compute the n-th stage of refinement of an extended trapezoidal rule
!!  adding 2^(n-2) additional interior point in the finite range of integration
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  nn=integer defining the refinement of the mesh, each call adds 2^(n-2) additional interior points
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  quad=the integral at the n-th stage.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!  When called with nn=1, the routine returns the crudest estimate of the integral
!!  Subsequent calls with nn=2,3,... (in that sequential order) will improve the accuracy
!!  by adding 2^(n-2) additional interior points. Note that quad should not be modified between sequential calls.
!!  Subroutine is defined as recursive to allow multi-dimensional integrations
!!
!! SOURCE

recursive subroutine trapezoidal_(func,nn,xmin,xmax,quad)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 !real(dp),external :: func
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(inout) :: quad

 interface
   function func(x)
     use defs_basis
     real(dp),intent(in) :: x
     real(dp) :: func
   end function func
 end interface

 !interface
 ! function func(x)
 !  use defs_basis
 !  real(dp),intent(in) :: x(:)
 !  real(dp) :: func(SIZE(x))
 ! end function func
 !end interface

!Local variables-------------------------------
!scalars
 integer :: npt,ix
 real(dp) :: space,new,yy
 character(len=500) :: msg
!arrays
 !real(dp),allocatable :: xx(:)
!************************************************************************

 select case (nn)

 case (1)
   ! === Initial crude estimate (xmax-xmin)(f1+f2)/2 ===
   !quad=half*(xmax-xmin)*SUM(func((/xmin,xmax/)))
   quad=half*(xmax-xmin)*(func(xmin)+func(xmax))

 case (2:)
   ! === Add npt interior points of spacing space ===
   npt=2**(nn-2) ; space=(xmax-xmin)/npt
   ! === The new sum is combined with the old integral to give a refined integral ===
   !new=SUM(func(arth(xmin+half*space,space,npt))) !PARALLEL version
   !allocate(xx(npt))
   !xx(:)=arth(xmin+half*space,space,npt)
   !xx(1)=xmin+half*space
   !do ii=2,nn
   ! xx(ii)=xx(ii-1)+space
   !end do
   new=zero
   yy=xmin+half*space
   do ix=1,npt
    !new=new+func(xx(ix))
    new=new+func(yy)
    yy=yy+space
   end do
   !deallocate(xx)
   quad=half*(quad+space*new)
   !write(std_out,*) 'trapezoidal',quad

 case (:0)
   write(msg,'(a,i3)')'Wrong value for nn ',nn
   MSG_BUG(msg)
 end select

end subroutine trapezoidal_
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/midpoint_
!! NAME
!!  midpoint_ (PRIVATE)
!!
!! FUNCTION
!!  This routine computes the n-th stage of refinement of an extended midpoint rule.
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  nn=integer defining the refinement of the mesh, each call adds (2/3)*3n-1 additional
!!   interior points between xmin ans xmax
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  quad=the integral at the n-th stage.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!  When called with nn=1, the routine returns as quad the crudest estimate of the integral
!!  Subsequent calls with nn=2,3,... (in that sequential order) will improve the accuracy of quad by adding
!!  (2/3)*3n-1 additional interior points. quad should not be modified between sequential calls.
!!  Subroutine is defined as recursive to allow multi-dimensional integrations
!!
!! SOURCE

 recursive subroutine midpoint_(func,nn,xmin,xmax,quad)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 !real(dp),external :: func
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(inout) :: quad

 interface
   function func(x)
     use defs_basis
     real(dp),intent(in) :: x
     real(dp) :: func
   end function func
 end interface

 !interface
 ! function func(x)
 !  use defs_basis
 !  real(dp),intent(in) :: x(:)
 !  real(dp) :: func(SIZE(x))
 ! end function func
 !end interface

!Local variables-------------------------------
!scalars
 integer  :: npt,ix
 real(dp) :: space
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: xx(:)

!************************************************************************

 select case (nn)

 case (1)
   ! === Initial crude estimate done at the middle of the interval
   !quad=(xmax-xmin)*SUM(func((/half*(xmin+xmax)/))) !PARALLEL version
   quad=(xmax-xmin)*func(half*(xmin+xmax))

 case (2:)
   ! === Add npt interior points, they alternate in spacing between space and 2*space ===
   ABI_MALLOC(xx,(2*3**(nn-2)))
   npt=3**(nn-2) ; space=(xmax-xmin)/(three*npt)
   xx(1:2*npt-1:2)=arth(xmin+half*space,three*space,npt)
   xx(2:2*npt:2)=xx(1:2*npt-1:2)+two*space
   ! === The new sum is combined with the old integral to give a refined integral ===
   !quad=quad/three+space*SUM(func(xx))  !PARALLEL version
   quad=quad/three
   do ix=1,SIZE(xx)
     quad=quad+space*func(xx(ix))
   end do
   ABI_FREE(xx)

 case (:0)
   write(msg,'(a,i3)')' wrong value for nn ',nn
   MSG_BUG('Wrong value for nn')
 end select

end subroutine midpoint_
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/quadrature
!! NAME
!!  quadrature
!!
!! FUNCTION
!!  Driver routine to perform quadratures in finite domains using different techniques.
!!  The routine improves the resolution of the grid until a given accuracy is reached
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  npts=Initial number of points, only for Gauss-Legendre. At each step this number is doubled
!!  accuracy=fractional accuracy required
!!  ntrial=Max number of attempts
!!  qopt=integer flag defining the algorithm for the quadrature:
!!    1 for Trapezoidal rule, closed, O(1/N^2)
!!    2 for Simpson based on trapezoidal,closed, O(1/N^4)
!!    3 for Midpoint rule, open, O(1/N^2)
!!    4 for midpoint rule with cancellation of leading error, open, O(1/N^4)
!!    5 for Romberg integration (closed form) and extrapolation for h-->0 (order 10 is hard-coded)
!!    6 for Romberg integration with midpoint rule and extrapolation for h-->0 (order 10 is hard-coded)
!!    7 for Gauss-Legendre
!!
!! OUTPUT
!!  quad=the integral
!!  ierr=0 if quadrature converged.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

recursive subroutine quadrature(func,xmin,xmax,qopt,quad,ierr,ntrial,accuracy,npts)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: qopt
 integer,intent(out) :: ierr
 integer,optional,intent(in) :: ntrial,npts
 real(dp),intent(in) :: xmin,xmax
 real(dp),optional,intent(in) :: accuracy
 real(dp),intent(out) :: quad

 interface
   function func(x)
     use defs_basis
     real(dp),intent(in) :: x
     real(dp) :: func
   end function func
 end interface

 !interface
 ! function func(x)
 !  use defs_basis
 !  real(dp),intent(in) :: x(:)
 !  real(dp) :: func(SIZE(x))
 ! end function func
 !end interface

!Local variables-------------------------------
!scalars
 integer :: K,KM,NT,NX,NX0,it,ix
 real(dp) :: EPS,old_st,st,old_quad,dqromb
 real(dp) :: TOL
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: h(:),s(:)
 real(dp),allocatable :: wx(:),xx(:)
! *************************************************************************

 ierr = 0
 TOL  =tol12
 EPS  =tol6  ; if (PRESENT(accuracy)) EPS=accuracy
 NT   =20    ; if (PRESENT(ntrial  )) NT=ntrial
 quad =zero

 select case (qopt)

 case (1)
   ! === Trapezoidal, closed form, O(1/N^2)
   do it=1,NT
     call trapezoidal_(func,it,xmin,xmax,quad)
     if (it>5) then ! Avoid spurious early convergence
       if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
     end if
     old_quad=quad
   end do

 case (2)
  ! === Extended Simpson rule based on trapezoidal O(1/N^4) ===
   do it=1,NT
     call trapezoidal_(func,it,xmin,xmax,st)
     if (it==1) then
       quad=st
     else
       quad=(four*st-old_st)/three
     end if
     if (it>5) then ! Avoid spurious early convergence
       if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
     end if
     old_quad=quad
     old_st=st
   end do

 case (3)
  ! === Midpoint rule, open form, O(1/N^2) ===
  do it=1,NT
    call midpoint_(func,it,xmin,xmax,quad)
    if (it>4) then ! Avoid spurious early convergence
      if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
    end if
    old_quad=quad
  end do

 case (4)
   ! === Midpoint rule with cancellation of leading 1/N^2 term, open form, O(1/N^4) ===
   do it=1,NT
     call midpoint_(func,it,xmin,xmax,st)
     if (it==1) then
       quad=st
     else
       quad=(nine*st-old_st)/eight
     end if
     if (it>4) then ! Avoid spurious early convergence
       if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
     end if
     old_quad=quad
     old_st=st
   end do

 case (5)
   ! === Romberg Integration, closed form ===
   K=5 ; KM=K-1 ! Order 10
   ABI_MALLOC(h,(NT+1))
   ABI_MALLOC(s,(NT+1))
   h=zero
   s=zero
   h(1)=one
   do it=1,NT
     call trapezoidal_(func,it,xmin,xmax,s(it))
     !write(std_out,*) ' romberg-trap at ',ncall,it,s(it)
     if (it>=K) then
       call polyn_interp(h(it-KM:it),s(it-KM:it),zero,quad,dqromb)
       if (ABS(dqromb)<EPS*ABS(quad)) then
         ABI_FREE(h)
         ABI_FREE(s)
         RETURN
       end if
     end if
     s(it+1)=s(it)
     h(it+1)=quarter*h(it) ! Quarter makes the extrapolation a polynomial in h^2,
   end do                 ! This is required to use the Euler-Maclaurin formula
   ABI_FREE(h)
   ABI_FREE(s)

 case (6)
   ! === Romberg Integration, closed form ===
   K=5 ; KM=K-1 ! Order 10
   ABI_MALLOC(h,(NT+1))
   ABI_MALLOC(s,(NT+1))
   h=zero
   s=zero
   h(1)=one
   do it=1,NT
     call midpoint_(func,it,xmin,xmax,s(it))
     if (it>=K) then
       call polyn_interp(h(it-KM:it),s(it-KM:it),zero,quad,dqromb)
       !write(std_out,*) quad,dqromb
       if (ABS(dqromb)<EPS*ABS(quad)) then
         ABI_FREE(h)
         ABI_FREE(s)
         RETURN
       end if
     end if
     s(it+1)=s(it)
     h(it+1)=ninth*h(it) ! factor is due to step tripling in midpoint and even error series
   end do
   ABI_FREE(h)
   ABI_FREE(s)

 case (7)
   ! === Gauss-Legendre ===
   NX0=5 ; if (PRESENT(npts)) NX0=npts
   NX=NX0
   do it=1,NT
     ABI_MALLOC(wx,(NX))
     ABI_MALLOC(xx,(NX))
     call coeffs_gausslegint(xmin,xmax,xx,wx,NX)
     quad=zero
     do ix=1,NX
       quad=quad+wx(ix)*func(xx(ix))
     end do
     ABI_FREE(wx)
     ABI_FREE(xx)
     if (it>1) then
       !write(std_out,*) quad
       if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
     end if
     old_quad=quad
     NX=NX+NX0
     !NX=2*NX
   end do

 case default
   write(msg,'(a,i3)')'Wrong value for qopt',qopt
   MSG_BUG(msg)
 end select

 write(msg,'(a,i0,2(a,es14.6))')&
&  "Results are not converged within the given accuracy. ntrial= ",NT,"; EPS= ",EPS,"; TOL= ",TOL
 MSG_WARNING(msg)
 ierr = -1

end subroutine quadrature
!!***

!!****f* m_numeric_tools/ctrap
!! NAME
!! ctrap
!!
!! FUNCTION
!! Do corrected trapezoidal integral on uniform grid of spacing hh.
!!
!! INPUTS
!!  imax=highest index of grid=grid point number of upper limit
!!  ff(imax)=integrand values
!!  hh=spacing between x points
!!
!! OUTPUT
!!  ans=resulting integral by corrected trapezoid
!!
!! NOTES
!!
!! PARENTS
!!      m_phonons,m_psp6,m_psptk,m_unittests,m_upf2abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine ctrap(imax,ff,hh,ans)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imax
 real(dp),intent(in) :: hh
 real(dp),intent(out) :: ans
!arrays
 real(dp),intent(in) :: ff(imax)

!Local variables-------------------------------
!scalars
 integer :: ir,ir2
 real(dp) :: endpt,sum

! *************************************************************************

 if (imax>=10)then

!  endpt=end point correction terms (low and high ends)
   endpt  = (23.75d0*(ff(1)+ff(imax  )) &
&   + 95.10d0*(ff(2)+ff(imax-1)) &
&   + 55.20d0*(ff(3)+ff(imax-2)) &
&   + 79.30d0*(ff(4)+ff(imax-3)) &
&   + 70.65d0*(ff(5)+ff(imax-4)))/ 72.d0
   ir2 = imax - 5
   sum=0.00d0
   if (ir2 > 5) then
     do ir=6,ir2
       sum = sum + ff(ir)
     end do
   end if
   ans = (sum + endpt ) * hh

 else if (imax>=8)then
   endpt  = (17.0d0*(ff(1)+ff(imax  )) &
&   + 59.0d0*(ff(2)+ff(imax-1)) &
&   + 43.0d0*(ff(3)+ff(imax-2)) &
&   + 49.0d0*(ff(4)+ff(imax-3)) )/ 48.d0
   sum=0.0d0
   if(imax==9)sum=ff(5)
   ans = (sum + endpt ) * hh

 else if (imax==7)then
   ans = (17.0d0*(ff(1)+ff(imax  )) &
&   + 59.0d0*(ff(2)+ff(imax-1)) &
&   + 43.0d0*(ff(3)+ff(imax-2)) &
&   + 50.0d0* ff(4)                )/ 48.d0  *hh

 else if (imax==6)then
   ans = (17.0d0*(ff(1)+ff(imax  )) &
&   + 59.0d0*(ff(2)+ff(imax-1)) &
&   + 44.0d0*(ff(3)+ff(imax-2)) )/ 48.d0  *hh

 else if (imax==5)then
   ans = (     (ff(1)+ff(5)) &
&   + four*(ff(2)+ff(4)) &
&   + two * ff(3)         )/ three  *hh

 else if (imax==4)then
   ans = (three*(ff(1)+ff(4)) &
&   + nine *(ff(2)+ff(3))  )/ eight  *hh

 else if (imax==3)then
   ans = (     (ff(1)+ff(3)) &
&   + four* ff(2)         )/ three *hh

 else if (imax==2)then
   ans = (ff(1)+ff(2))/ two  *hh

 else if (imax==1)then
   ans = ff(1)*hh

 end if

end subroutine ctrap
!!***

!!****f* m_numeric_tools/cspint
!! NAME
!!  cspint
!!
!! FUNCTION
!!  Estimates the integral of a tabulated function.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!!    The routine is given the value of a function F(X) at a set of
!!    nodes XTAB, and estimates
!!
!!      Integral ( A <= X <= B ) F(X) DX
!!
!!    by computing the cubic natural spline S(X) that interpolates
!!    F(X) at the nodes, and then computing
!!
!!      Integral ( A <= X <= B ) S(X) DX
!!
!!    exactly.
!!
!!    Other output from the program includes the definite integral
!!    from X(1) to X(I) of S(X), and the coefficients necessary for
!!    the user to evaluate the spline S(X) at any point.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Parameters:
!!
!!    Input, real (dp) FTAB(NTAB), contains the tabulated values of
!!    the function, FTAB(I) = F(XTAB(I)).
!!
!!    Input, real (dp) XTAB(NTAB), contains the points at which the
!!    function was evaluated.  The XTAB's must be distinct and
!!    in ascending order.
!!
!!    Input, integer NTAB, the number of entries in FTAB and
!!    XTAB.  NTAB must be at least 3.
!!
!!    Input, real (dp) A, lower limit of integration.
!!
!!    Input, real (dp) B, upper limit of integration.
!!
!!    Output, real (dp) Y(3,NTAB), will contain the coefficients
!!    of the interpolating natural spline over each subinterval.
!!
!!    For XTAB(I) <= X <= XTAB(I+1),
!!
!!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!!                   + Y(2,I)*(X-XTAB(I))**2
!!                   + Y(3,I)*(X-XTAB(I))**3
!!
!!    Output, real (dp) E(NTAB), E(I) = the definite integral from
!!    XTAB(1) to XTAB(I) of S(X).
!!
!!    Workspace, real (dp) WORK(NTAB).
!!
!!    Output, real (dp) RESULT, the estimated value of the integral.
!!
!!
!! PARENTS
!!      m_xc_vdw,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: ntab
  real(dp), intent(in) :: a
  real(dp), intent(in) :: b
  real(dp), intent(inout) :: e(ntab)
  real(dp), intent(in) :: ftab(ntab)
  real(dp), intent(inout) :: work(ntab)
  real(dp), intent(in) :: xtab(ntab)
  real(dp), intent(inout) :: y(3,ntab)
  real(dp), intent(out) :: result

!Local variables ------------------------------
!scalars
  integer :: i
  integer :: j
  real(dp) :: r
  real(dp) :: s
  real(dp) :: term
  real(dp) :: u

!************************************************************************

  if ( ntab < 3 ) then
    write(std_out,'(a)' ) ' '
    write(std_out,'(a)' ) 'CSPINT - Fatal error!'
    write(std_out,'(a,i6)' ) '  NTAB must be at least 3, but input NTAB = ',ntab
    MSG_ERROR("Aborting now")
  end if

  do i = 1, ntab-1

    if ( xtab(i+1) <= xtab(i) ) then
      write(std_out,'(a)' ) ' '
      write(std_out,'(a)' ) 'CSPINT - Fatal error!'
      write(std_out,'(a)' ) '  Nodes not in strict increasing order.'
      write(std_out,'(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
      write(std_out,'(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
      write(std_out,'(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
      MSG_ERROR("Aborting now")
    end if

  end do

  s = zero
  do i = 1, ntab-1
    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
    y(2,i) = r - s
    s = r
  end do

  result = zero
  s = zero
  r = zero
  y(2,1) = zero
  y(2,ntab) = zero

  do i = 2, ntab-1
    y(2,i) = y(2,i) + r * y(2,i-1)
    work(i) = two * ( xtab(i-1) - xtab(i+1) ) - r * s
    s = xtab(i+1) - xtab(i)
    r = s / work(i)
  end do

  do j = 2, ntab-1
    i = ntab+1-j
    y(2,i) = ( ( xtab(i+1) - xtab(i) ) * y(2,i+1) - y(2,i) ) / work(i)
  end do

  do i = 1, ntab-1
    s = xtab(i+1) - xtab(i)
    r = y(2,i+1) - y(2,i)
    y(3,i) = r / s
    y(2,i) = three * y(2,i)
    y(1,i) = ( ftab(i+1) - ftab(i) ) / s - ( y(2,i) + r ) * s
  end do

  e(1) = 0.0D+00
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    term = ((( y(3,i) * quarter * s + y(2,i) * third ) * s &
      + y(1,i) * half ) * s + ftab(i) ) * s
    e(i+1) = e(i) + term
  end do
!
!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!
  r = a
  u = one

  do j = 1, 2
!
!  The endpoint is less than or equal to XTAB(1).
!
    if ( r <= xtab(1) ) then
      result = result-u*((r-xtab(1))*y(1,1)*half +ftab(1))*(r-xtab(1))
!
!  The endpoint is greater than or equal to XTAB(NTAB).
!
    else if ( xtab(ntab) <= r ) then

      result = result -u * ( e(ntab) + ( r - xtab(ntab) ) &
        * ( ftab(ntab) + half * ( ftab(ntab-1) &
        + ( xtab(ntab) - xtab(ntab-1) ) * y(1,ntab-1) ) &
        * ( r - xtab(ntab) )))
!
!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!
    else

      do i = 1, ntab-1

        if ( r <= xtab(i+1) ) then
          r = r-xtab(i)
          result = result-u*(e(i)+(((y(3,i)*quarter*r+y(2,i)*third)*r &
            +y(1,i)*half )*r+ftab(i))*r)
          go to 120
        end if

      end do

    end if

  120   continue

    u = -one
    r = b

  end do

end subroutine cspint
!!***

!!****f* m_numeric_tools/coeffs_gausslegint
!! NAME
!!  coeffs_gausslegint
!!
!! FUNCTION
!! Compute the coefficients (supports and weights) for Gauss-Legendre integration.
!! Inspired by a routine due to G. Rybicki.
!!
!! INPUTS
!! xmin=lower bound of integration
!! xmax=upper bound of integration
!! n=order of integration
!!
!! OUTPUT
!! x(n)=array of support points
!! weights(n)=array of integration weights
!!
!! PARENTS
!!      m_bader,m_numeric_tools,m_screening_driver,m_sigc
!!
!! CHILDREN
!!
!! SOURCE

subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(out) :: x(n),weights(n)

!Local variables ------------------------------
!scalars
 integer :: i,j
 real(dp),parameter :: tol=1.d-13
 real(dp),parameter :: pi=4.d0*atan(1.d0)
 real(dp) :: z,z1,xmean,p1,p2,p3,pp,xl

!************************************************************************

 xl=(xmax-xmin)*0.5d0
 xmean=(xmax+xmin)*0.5d0

 do i=1,(n+1)/2
  z=cos(pi*(i-0.25d0)/(n+0.5d0))

  do
    p1=1.d0
    p2=0.d0

    do j=1,n
     p3=p2
     p2=p1
     p1=((2.d0*j - 1.d0)*z*p2 - (j-1.d0)*p3)/j
    end do

    pp=n*(p2-z*p1)/(1.0d0-z**2)
    z1=z
    z=z1-p1/pp

    if(abs(z-z1) < tol) exit
  end do

  x(i)=xmean-xl*z
  x(n+1-i)=xmean+xl*z
  weights(i)=2.d0*xl/((1.d0-z**2)*pp**2)
  weights(n+1-i)=weights(i)
 end do

end subroutine coeffs_gausslegint
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/simpson_cplx
!! NAME
!!  simpson_cplx
!!
!! FUNCTION
!!  Integrate a complex function using extended Simpson's rule.
!!
!! INPUTS
!!  npts=Number of points.
!!  step=Step of the mesh.
!!  ff(npts)=Values of the integrand.
!!
!! OUTPUT
!!  simpson_cplx=Integral of ff.
!!
!! NOTES
!!  If npts is odd, the integration is done with the extended Simpson's rule (error = O^(step^4))
!!  If npts is even, the last 4 four points are integrated separately via Simpson's 3/8 rule. Error = O(step^5)
!!  while the first npts-3 points are integrared with the extended Simpson's rule.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function simpson_cplx(npts,step,ff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
 real(dp),intent(in)  :: step
 complex(dpc),intent(in) :: ff(npts)
 complex(dpc) :: simpson_cplx

!Local variables ------------------------------
!scalars
 integer :: ii,my_n
 complex(dpc) :: sum_even, sum_odd

!************************************************************************

 my_n=npts; if ((npts/2)*2 == npts) my_n=npts-3

 if (my_n<2) then
   MSG_ERROR("Too few points")
 end if

 sum_odd=czero
 do ii=2,my_n-1,2
   sum_odd = sum_odd + ff(ii)
 end do

 sum_even=zero
 do ii=3,my_n-2,2
   sum_even = sum_even + ff(ii)
 end do

 ! Eq 25.4.6 Abramowitz. Error is O(step^4)
 simpson_cplx = step/three * (ff(1) + four*sum_odd + two*sum_even + ff(my_n))

 if (my_n/=npts) then ! Simpson's 3/8 rule. Eq 25.4.13 Abramowitz. Error is O(step^5)
  simpson_cplx = simpson_cplx + three*step/eight * (ff(npts-3) + 3*ff(npts-2) + 3*ff(npts-1) + ff(npts))
 end if

end function simpson_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/hermitianize_spc
!! NAME
!!  hermitianize_spc
!!
!! FUNCTION
!!  Force a square matrix to be hermitian
!!
!! INPUTS
!!  uplo=String describing which part of the matrix has been calculated.
!!    Only the first character is tested (no case sensitive). Possible values are:
!!    "All"= Full matrix is supplied in input
!!    "Upper"=Upper triangle is in input. Lower triangle is reconstructed by symmetry.
!!    "Lower"=Lower triangle is in input. Upper triangle is reconstructed by symmetry.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mat(:,:)=complex input matrix, hermitianized at output
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine hermitianize_spc(mat,uplo)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: uplo
!arrays
 complex(spc),intent(inout) :: mat(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii,jj
!arrays
 complex(spc),allocatable :: tmp(:)
! *************************************************************************

 nn = assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',__FILE__,__LINE__)

 select case (uplo(1:1))

 case ("A","a") ! Full matrix has been calculated.
   ABI_MALLOC(tmp,(nn))
   do ii=1,nn
     do jj=ii,nn
       ! reference half constant is dp not sp
       tmp(jj)=real(half)*(mat(ii,jj)+CONJG(mat(jj,ii)))
     end do
     mat(ii,ii:nn)=tmp(ii:nn)
     mat(ii:nn,ii)=CONJG(tmp(ii:nn))
   end do
   ABI_FREE(tmp)

 case ("U","u") ! Only the upper triangle is used.
   do jj=1,nn
     do ii=1,jj
       if (ii/=jj) then
         mat(jj,ii) = CONJG(mat(ii,jj))
       else
         mat(ii,ii) = CMPLX(REAL(mat(ii,ii)),0.0_sp)
       end if
     end do
   end do

 case ("L","l") ! Only the lower triangle is used.
  do jj=1,nn
    do ii=1,jj
      if (ii/=jj) then
        mat(ii,jj) = CONJG(mat(jj,ii))
      else
        mat(ii,ii) = CMPLX(REAL(mat(ii,ii)),0.0_sp)
      end if
    end do
  end do

 case default
   MSG_ERROR("Wrong uplo"//TRIM(uplo))
 end select

end subroutine hermitianize_spc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/hermitianize_dpc
!! NAME
!!  hermitianize_dpc
!!
!! FUNCTION
!!  Force a square matrix to be hermitian
!!
!! INPUTS
!!  uplo=String describing which part of the matrix has been calculated.
!!    Only the first character is tested (no case sensitive). Possible values are:
!!    "All"= Full matrix is supplied in input
!!    "Upper"=Upper triangle is in input. Lower triangle is reconstructed by symmetry.
!!    "Lower"=Lower triangle is in input. Upper triangle is reconstructed by symmetry.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mat(:,:)=complex input matrix, hermitianized in output
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine hermitianize_dpc(mat,uplo)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: uplo
!arrays
 complex(dpc),intent(inout) :: mat(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii,jj
!arrays
 complex(dpc),allocatable :: tmp(:)
! *************************************************************************

 nn = assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',__FILE__,__LINE__)

 select case (uplo(1:1))

 case ("A","a") ! Full matrix has been calculated.
   ABI_MALLOC(tmp,(nn))
   do ii=1,nn
     do jj=ii,nn
       tmp(jj)=half*(mat(ii,jj)+DCONJG(mat(jj,ii)))
     end do
     mat(ii,ii:nn)=tmp(ii:nn)
     mat(ii:nn,ii)=DCONJG(tmp(ii:nn))
   end do
   ABI_FREE(tmp)

 case ("U","u") ! Only the upper triangle is used.
   do jj=1,nn
     do ii=1,jj
       if (ii/=jj) then
         mat(jj,ii) = DCONJG(mat(ii,jj))
       else
         mat(ii,ii) = CMPLX(DBLE(mat(ii,ii)),zero, kind=dpc)
       end if
     end do
   end do

 case ("L","l") ! Only the lower triangle is used.
  do jj=1,nn
    do ii=1,jj
      if (ii/=jj) then
        mat(ii,jj) = DCONJG(mat(jj,ii))
      else
        mat(ii,ii) = CMPLX(REAL(mat(ii,ii)),zero, kind=dpc)
      end if
    end do
  end do

 case default
   MSG_ERROR("Wrong uplo"//TRIM(uplo))
 end select

end subroutine hermitianize_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/mkherm
!! NAME
!! mkherm
!!
!! FUNCTION
!! Make the complex array(ndim,ndim) hermitian,
!! by adding half of it to its hermitian conjugate.
!!
!! INPUTS
!! ndim=dimension of the matrix
!! array= complex matrix
!!
!! SIDE EFFECTS
!! array= hermitian matrix made by adding half of array to its hermitian conjugate
!!
!! PARENTS
!!      anaddb,dfpt_phfrq,m_ddb,m_phgamma
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine mkherm(array,ndim)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ndim
!arrays
 real(dp),intent(inout) :: array(2,ndim,ndim)

!Local variables -------------------------
!scalars
 integer :: i1,i2

! *********************************************************************

 do i1=1,ndim
   do i2=1,i1
     array(1,i1,i2)=(array(1,i1,i2)+array(1,i2,i1))*half
     array(2,i1,i2)=(array(2,i1,i2)-array(2,i2,i1))*half
     array(1,i2,i1)=array(1,i1,i2)
     array(2,i2,i1)=-array(2,i1,i2)
   end do
 end do

end subroutine mkherm
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/hermit
!! NAME
!! hermit
!!
!! FUNCTION
!! Take a matrix in hermitian storage mode (lower triangle stored)
!! and redefine diagonal elements to impose Hermiticity
!! (diagonal terms have to be real).
!! If abs(Im(H(i,i)))>4096*machine precision, print error warning.
!! (Typical 64 bit machine precision is 2^-52 or 2.22e-16)
!!
!! INPUTS
!!  chmin(n*n+n)=complex hermitian matrix with numerical noise possibly
!!   rendering Im(diagonal elements) approximately 1e-15 or so
!!  ndim=size of complex hermitian matrix
!!
!! OUTPUT
!!  chmout(n*n+n)=redefined matrix with strictly real diagonal elements.
!!   May be same storage location as chmin.
!!  ierr=0 if no problem, 1 if the imaginary part of some element
!!   too large (at present, stop in this case).
!!
!! TODO
!!  Name is misleading, perhaps hermit_force_diago?
!!  Interface allows aliasing
!!
!! PARENTS
!!      m_cgtools,m_extraprho
!!
!! CHILDREN
!!
!! SOURCE

subroutine hermit(chmin,chmout,ierr,ndim)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim
 integer,intent(out) :: ierr
!arrays
 real(dp),intent(inout) :: chmin(ndim*ndim+ndim)
 real(dp),intent(inout) :: chmout(ndim*ndim+ndim)

!Local variables-------------------------------
!scalars
 integer,save :: mmesgs=20,nmesgs=0
 integer :: idim,merrors,nerrors
 real(dp),parameter :: eps=epsilon(0.0d0)
 real(dp) :: ch_im,ch_re,moduls,tol
 character(len=500) :: message

! *************************************************************************

 tol=4096.0d0*eps

 ierr=0
 merrors=0

!Copy matrix into possibly new location
 chmout(:)=chmin(:)

!Loop over diagonal elements of matrix (off-diag not altered)
 do idim=1,ndim

   ch_im=chmout(idim*idim+idim  )
   ch_re=chmout(idim*idim+idim-1)

!  check for large absolute Im part and print warning when
!  larger than (some factor)*(machine precision)
   nerrors=0
   if( abs(ch_im) > tol .and. abs(ch_im) > tol8*abs(ch_re)) nerrors=2
   if( abs(ch_im) > tol .or. abs(ch_im) > tol8*abs(ch_re)) nerrors=1

   if( (abs(ch_im) > tol .and. nmesgs<mmesgs) .or. nerrors==2)then
     write(message, '(3a,i0,a,es20.12,a,es20.12,a)' )&
&     'Input Hermitian matrix has nonzero relative Im part on diagonal:',ch10,&
&     'for component',idim,' Im part is',ch_im,', Re part is',ch_re,'.'
     call wrtout(std_out,message,'PERS')
     nmesgs=nmesgs+1
   end if

   if( ( abs(ch_im) > tol8*abs(ch_re) .and. nmesgs<mmesgs) .or. nerrors==2)then
     write(message, '(3a,i0,a,es20.12,a,es20.12,a)' )&
&     'Input Hermitian matrix has nonzero relative Im part on diagonal:',ch10,&
&     'for component',idim,' Im part is',ch_im,', Re part is',ch_re,'.'
     call wrtout(std_out,message,'PERS')
     nmesgs=nmesgs+1
   end if

!  compute modulus $= (\Re^2+\Im^2)^{1/2}$
   moduls=sqrt(ch_re**2+ch_im**2)

!  set Re part to modulus with sign of original Re part
   chmout(idim*idim+idim-1)=sign(moduls,ch_re)

!  set Im part to 0
   chmout(idim*idim+idim)=zero

   merrors=max(merrors,nerrors)

 end do

 if(merrors==2)then
   ierr=1
   write(message, '(3a)' )&
&   'Imaginary part(s) of diagonal Hermitian matrix element(s) is too large.',ch10,&
&   'See previous messages.'
   MSG_BUG(message)
 end if

end subroutine hermit
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/symmetrize_spc
!! NAME
!!  symmetrize_spc
!!
!! FUNCTION
!!  Force a square matrix to be symmetric.
!!
!! INPUTS
!!  uplo=String describing which part of the matrix has been calculated.
!!    Only the first character is tested (no case sensitive). Possible values are:
!!    "All"= Full matrix is supplied in input
!!    "Upper"=Upper triangle is in input. Lower triangle is reconstructed by symmetry.
!!    "Lower"=Lower triangle is in input. Upper triangle is reconstructed by symmetry.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mat(:,:)=complex input matrix, symmetrized at output
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine symmetrize_spc(mat,uplo)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: uplo
!arrays
 complex(spc),intent(inout) :: mat(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii,jj
!arrays
 complex(spc),allocatable :: tmp(:)
! *************************************************************************

 nn = assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',__FILE__,__LINE__)

 select case (uplo(1:1))

 case ("A","a") ! Full matrix has been calculated.
   ABI_MALLOC(tmp,(nn))
   do ii=1,nn
     do jj=ii,nn
       tmp(jj)=REAL(half)*(mat(ii,jj)+mat(jj,ii))
     end do
     mat(ii,ii:nn)=tmp(ii:nn)
     mat(ii:nn,ii)=tmp(ii:nn)
   end do
   ABI_FREE(tmp)

 case ("U","u") ! Only the upper triangle is used.
   do jj=1,nn
     do ii=1,jj-1
       mat(jj,ii) = mat(ii,jj)
     end do
   end do

 case ("L","l") ! Only the lower triangle is used.
  do jj=1,nn
    do ii=1,jj-1
      mat(ii,jj) = mat(jj,ii)
    end do
  end do

 case default
   MSG_ERROR("Wrong uplo"//TRIM(uplo))
 end select

end subroutine symmetrize_spc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/symmetrize_dpc
!! NAME
!!  symmetrize_dpc
!!
!! FUNCTION
!!  Force a square matrix to be symmetric.
!!
!! INPUTS
!!  uplo=String describing which part of the matrix has been calculated.
!!    Only the first character is tested (no case sensitive). Possible values are:
!!    "All"= Full matrix is supplied in input
!!    "Upper"=Upper triangle is in input. Lower triangle is reconstructed by symmetry.
!!    "Lower"=Lower triangle is in input. Upper triangle is reconstructed by symmetry.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mat(:,:)=complex input matrix, symmetrized in output
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine symmetrize_dpc(mat,uplo)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: uplo
!arrays
 complex(dpc),intent(inout) :: mat(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii,jj
!arrays
 complex(dpc),allocatable :: tmp(:)
! *************************************************************************

 nn = assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',__FILE__,__LINE__)

 select case (uplo(1:1))
 case ("A","a") ! Full matrix has been calculated.
   ABI_MALLOC(tmp,(nn))
   do ii=1,nn
     do jj=ii,nn
       tmp(jj)=half*(mat(ii,jj)+mat(jj,ii))
     end do
     mat(ii,ii:nn)=tmp(ii:nn)
     mat(ii:nn,ii)=tmp(ii:nn)
   end do
   ABI_FREE(tmp)

 case ("U","u") ! Only the upper triangle is used.
   do jj=1,nn
     do ii=1,jj-1
       mat(jj,ii) = mat(ii,jj)
     end do
   end do

 case ("L","l") ! Only the lower triangle is used.
  do jj=1,nn
    do ii=1,jj-1
      mat(ii,jj) = mat(jj,ii)
    end do
  end do

 case default
   MSG_ERROR("Wrong uplo"//TRIM(uplo))
 end select

end subroutine symmetrize_dpc
!!***

!!****f* m_numeric_tools/pack_matrix
!! NAME
!! pack_matrix
!!
!! FUNCTION
!! Packs a matrix into hermitian format
!!
!! INPUTS
!! N: size of matrix
!! cplx: 2 if matrix is complex, 1 for real matrix.
!! mat_in(cplx, N*N)= matrix to be packed
!!
!! OUTPUT
!! mat_out(cplx*N*N+1/2)= packed matrix (upper triangle)
!!
!! PARENTS
!!      m_rayleigh_ritz
!!
!! CHILDREN
!!
!! SOURCE

subroutine pack_matrix(mat_in, mat_out, N, cplx)

 integer, intent(in) :: N, cplx
 real(dp), intent(in) :: mat_in(cplx, N*N)
 real(dp), intent(out) :: mat_out(cplx*N*(N+1)/2)

!local variables
 integer :: isubh, i, j

 ! *************************************************************************

 isubh = 1
 do j=1,N
   do i=1,j
     mat_out(isubh)    = mat_in(1, (j-1)*N+i)
     ! bad for vectorization, but it's not performance critical, so ...
     if(cplx == 2) then
       mat_out(isubh+1)  = mat_in(2, (j-1)*N+i)
     end if
     isubh=isubh+cplx
   end do
 end do

end subroutine pack_matrix
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/print_arr1d_spc
!! NAME
!!  print_arr1d_spc
!!
!! FUNCTION
!! Print an array using a nice (?) format
!!
!! INPUTS
!!  arr(:)=vector/matrix to be printed
!!  mode_paral(optional)=parallel mode, DEFAULT is "COLL"
!!   "COLL" if all procs are calling the routine with the same message to be written only once
!!   "PERS" if the procs are calling the routine with different mesgs each to be written,
!!          or if one proc is calling the routine
!!  unit(optional)=the unit number of the file, DEFAULT=std_out
!!  max_r,max_c(optional)=Max number of rows and columns to be printed
!!   (DEFAULT is 9, output format assumes to be less that 99, but there might be
!!    problems with wrtout if message size exceeds 500 thus max number of elements should be ~60)
!!
!! OUTPUT
!!  (only printing)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_arr1d_spc(arr,max_r,unit,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,max_r
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 complex(spc),intent(in) :: arr(:)

!Local variables-------------------------------
!scalars
 integer :: unt,ii,nr,mr
 character(len=4) :: mode
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1
! *************************************************************************

 unt=std_out ; if (PRESENT(unit      )) unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 mr=15       ; if (PRESENT(max_r     )) mr=max_r

 if (mode/='COLL'.and.mode/='PERS') then
  write(msg,'(2a)')' Wrong value of mode_paral ',mode
  MSG_BUG(msg)
 end if
 !
 ! === Print out matrix ===
 nr=SIZE(arr,DIM=1) ; if (mr>nr) mr=nr

 write(fmth,*)'(6x,',mr,'(i2,6x))'
 write(fmt1,*)'(3x,',mr,'f8.3)'

 write(msg,fmth)(ii,ii=1,mr)
 call wrtout(unt,msg,mode) !header
 write(msg,fmt1)REAL (arr(1:mr))
 call wrtout(unt,msg,mode) !real part
 write(msg,fmt1)AIMAG(arr(1:mr))
 call wrtout(unt,msg,mode) !imag part

end subroutine print_arr1d_spc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/print_arr1d_dpc
!! NAME
!!  print_arr1d_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_arr1d_dpc(arr,max_r,unit,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,max_r
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 complex(dpc),intent(in) :: arr(:)

!Local variables-------------------------------
!scalars
 integer :: unt,ii,nr,mr
 character(len=4) :: mode
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1
! *************************************************************************

 unt=std_out ; if (PRESENT(unit      )) unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 mr=15       ; if (PRESENT(max_r     )) mr=max_r

 if (mode/='COLL'.and.mode/='PERS') then
  write(msg,'(2a)')' Wrong value of mode_paral ',mode
  MSG_BUG(msg)
 end if
 !
 ! === Print out matrix ===
 nr=SIZE(arr,DIM=1) ; if (mr>nr) mr=nr

 write(fmth,*)'(6x,',mr,'(i2,6x))'
 write(fmt1,*)'(3x,',mr,'f8.3)'

 write(msg,fmth)(ii,ii=1,mr)
 call wrtout(unt,msg,mode) !header
 write(msg,fmt1)REAL (arr(1:mr))
 call wrtout(unt,msg,mode) !real part
 write(msg,fmt1)AIMAG(arr(1:mr))
 call wrtout(unt,msg,mode) !imag part

end subroutine print_arr1d_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/print_arr2d_spc
!! NAME
!!  print_arr2d_spc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_arr2d_spc(arr,max_r,max_c,unit,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,max_r,max_c
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 complex(spc),intent(in) :: arr(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,ii,jj,nc,nr,mc,mr
 character(len=4) :: mode
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1,fmt2
! *************************************************************************

 unt =std_out; if (PRESENT(unit      )) unt =unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 mc  =9      ; if (PRESENT(max_c     )) mc  =max_c
 mr  =9      ; if (PRESENT(max_r     )) mr  =max_r

 if (mode/='COLL'.and.mode/='PERS') then
   write(msg,'(2a)')'Wrong value of mode_paral ',mode
   MSG_BUG(msg)
 end if
 !
 ! === Print out matrix ===
 nr=SIZE(arr,DIM=1); if (mr>nr) mr=nr
 nc=SIZE(arr,DIM=2); if (mc>nc) mc=nc

 write(fmth,*)'(6x,',mc,'(i2,6x))'
 write(fmt1,*)'(3x,i2,',mc,'f8.3)'
 write(fmt2,*)'(5x   ,',mc,'f8.3,a)'

 write(msg,fmth)(jj,jj=1,mc)
 call wrtout(unt,msg,mode) !header
 do ii=1,mr
   write(msg,fmt1)ii,REAL(arr(ii,1:mc))
   call wrtout(unt,msg,mode) !real part
   write(msg,fmt2)  AIMAG(arr(ii,1:mc)),ch10
   call wrtout(unt,msg,mode) !imag part
 end do

end subroutine print_arr2d_spc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/print_arr2d_dpc
!! NAME
!!  print_arr2d_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_arr2d_dpc(arr,max_r,max_c,unit,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,max_r,max_c
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 complex(dpc),intent(in) :: arr(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,ii,jj,nc,nr,mc,mr
 character(len=4) :: mode
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1,fmt2
! *************************************************************************

 unt =std_out; if (PRESENT(unit      )) unt =unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 mc  =9      ; if (PRESENT(max_c     )) mc  =max_c
 mr  =9      ; if (PRESENT(max_r     )) mr  =max_r

 if (mode/='COLL'.and.mode/='PERS') then
   write(msg,'(2a)')'Wrong value of mode_paral ',mode
   MSG_BUG(msg)
 end if
 !
 ! === Print out matrix ===
 nr=SIZE(arr,DIM=1); if (mr>nr) mr=nr
 nc=SIZE(arr,DIM=2); if (mc>nc) mc=nc

 write(fmth,*)'(6x,',mc,'(i2,6x))'
 write(fmt1,*)'(3x,i2,',mc,'f8.3)'
 write(fmt2,*)'(5x   ,',mc,'f8.3,a)'

 write(msg,fmth)(jj,jj=1,mc)
 call wrtout(unt,msg,mode) !header
 do ii=1,mr
   write(msg,fmt1)ii,REAL(arr(ii,1:mc))
   call wrtout(unt,msg,mode) !real part
   write(msg,fmt2)  AIMAG(arr(ii,1:mc)),ch10
   call wrtout(unt,msg,mode) !imag part
 end do

end subroutine print_arr2d_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/pade
!! NAME
!!  pade
!!
!! FUNCTION
!!  Calculate the pade approximant in zz of the function f calculated at the n points z
!!
!! SOURCE

function pade(n,z,f,zz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: zz
 complex(dpc) :: pade
!arrays
 complex(dpc),intent(in) :: z(n),f(n)

!Local variables-------------------------------
!scalars
 complex(dpc) :: a(n)
 complex(dpc) :: Az(0:n), Bz(0:n)
 integer :: i
! *************************************************************************

 call calculate_pade_a(a,n,z,f)

 Az(0)=czero
 Az(1)=a(1)
 Bz(0)=cone
 Bz(1)=cone

 do i=1,n-1
   Az(i+1)=Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
   Bz(i+1)=Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
 end do
 !write(std_out,*) 'Bz(n)',Bz(n)
 if (REAL(Bz(n))==zero.and.AIMAG(Bz(n))==zero) write(std_out,*) ' Bz(n) ',Bz(n)
 pade=Az(n)/Bz(n)
 !write(std_out,*) 'pade_approx ', pade_approx

end function pade
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/dpade
!! NAME
!!  dpade
!!
!! FUNCTION
!!  Calculate the derivative of the pade approximant in zz of the function f calculated at the n points z
!!
!! SOURCE

function dpade(n,z,f,zz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: zz
 complex(dpc) :: dpade
!arrays
 complex(dpc),intent(in) :: z(n),f(n)

!Local variables-------------------------------
!scalars
 integer :: i
!arrays
 complex(dpc) :: a(n)
 complex(dpc) :: Az(0:n), Bz(0:n)
 complex(dpc) :: dAz(0:n), dBz(0:n)
! *************************************************************************

 call calculate_pade_a(a,n,z,f)

 Az(0)=czero
 Az(1)=a(1)
 Bz(0)=cone
 Bz(1)=cone
 dAz(0)=czero
 dAz(1)=czero
 dBz(0)=czero
 dBz(1)=czero

 do i=1,n-1
   Az(i+1)=Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
   Bz(i+1)=Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
   dAz(i+1)=dAz(i)+a(i+1)*Az(i-1)+(zz-z(i))*a(i+1)*dAz(i-1)
   dBz(i+1)=dBz(i)+a(i+1)*Bz(i-1)+(zz-z(i))*a(i+1)*dBz(i-1)
 end do
 !write(std_out,*) 'Bz(n)', Bz(n)
 if (REAL(Bz(n))==zero.and.AIMAG(Bz(n))==zero) write(std_out,*) 'Bz(n)',Bz(n)
 !pade_approx = Az(n) / Bz(n)
 dpade=dAz(n)/Bz(n) -Az(n)*dBz(n)/(Bz(n)*Bz(n))
 !write(std_out,*) 'pade_approx ', pade_approx

end function dpade
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/calculate_pade_a
!! NAME
!!  calculate_pade_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_numeric_tools
!!
!! CHILDREN
!!
!! SOURCE

subroutine calculate_pade_a(a,n,z,f)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: z(n),f(n)
 complex(dpc),intent(out) :: a(n)

!Local variables-------------------------------
!scalars
 integer :: i,j
!arrays
 complex(dpc) :: g(n,n)
! *************************************************************************

 g(1,1:n)=f(1:n)

 do i=2,n
   do j=i,n
     if (REAL(g(i-1,j))==zero.and.AIMAG(g(i-1,j))==zero) write(std_out,*) 'g_i(z_j)',i,j,g(i,j)
     g(i,j)=(g(i-1,i-1)-g(i-1,j)) / ((z(j)-z(i-1))*g(i-1,j))
     !write(std_out,*) 'g_i(z_j)',i,j,g(i,j)
   end do
 end do
 do i=1,n
   a(i)=g(i,i)
 end do
 !write(std_out,*) 'a ',a(:)

end subroutine calculate_pade_a
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/newrap_step
!! NAME
!!  newrap_step
!!
!! FUNCTION
!!  Apply single step newton-raphson method to find the root of a complex function
!!   z_k+1=z_k-f(z_k)/(df/dz(z_k))
!!
!! SOURCE

function newrap_step(z,f,df)

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: z,f,df
 complex(dpc) :: newrap_step

!Local variables-------------------------------
!scalars
 real(dp) :: dfm2
! *************************************************************************

 dfm2=ABS(df)*ABS(df)

 newrap_step= z - (f*CONJG(df))/dfm2
 !& z-one/(ABS(df)*ABS(df)) * CMPLX( REAL(f)*REAL(df)+AIMAG(f)*AIMAG(df), -REAL(f)*AIMAG(df)+AIMAG(f)*EAL(df) )

end function newrap_step
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cross_product_int
!! NAME
!!  cross_product_int
!!
!! FUNCTION
!!  Return the cross product of two vectors with integer components.
!!
pure function cross_product_int(vec1,vec2) result(res)

!Arguments ------------------------------------
 integer,intent(in) :: vec1(3),vec2(3)
 integer :: res(3)
! *************************************************************************

 res(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
 res(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
 res(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

end function cross_product_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cross_product_rdp
!! NAME
!!  cross_product_rdp
!!
!! FUNCTION
!!  Return the cross product of two vectors with real double precision components.
!!
pure function cross_product_rdp(vec1,vec2) result(res)

!Arguments ------------------------------------
 real(dp),intent(in) :: vec1(3),vec2(3)
 real(dp) :: res(3)
! *************************************************************************

 res(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
 res(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
 res(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

end function cross_product_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/l2norm_rdp
!! NAME
!!  l2norm_rdp
!!
!! FUNCTION
!!  Return the length (ordinary L2 norm) of a vector.
!!

pure function l2norm_rdp(vec) result(res)

!Arguments ------------------------------------
 real(dp),intent(in) :: vec(:)
 real(dp) :: res
! *************************************************************************

 res=SQRT(DOT_PRODUCT(vec,vec))

end function l2norm_rdp
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/remove_copies
!! NAME
!!  remove_copies
!!
!! FUNCTION
!!  Given an initial set of elements, set_in, return the subset of inequivalent items
!!  packed in the firt n_out positions of set_in. Use the logical function is_equal
!!  to define whether two items are equivalent.
!!
!! INPUTS
!! n_in=Initial number of elements.
!! is_equal=logical function used to discern if two items are equal.
!!
!! OUTPUT
!! n_out=Number of inequivalent items found.
!!
!! SIDE EFFECTS
!! set_in(3,n_in)=
!!  In input the initial set of n_in elements
!!  In output set_in(3,1:n_out) contains the inequivalent elements found.
!!
!! NOTES
!! The routines only deals with arrays of 3D-vectors although generalizing the
!! algorithm to nD-space is straightforward.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine remove_copies(n_in,set_in,n_out,is_equal)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n_in
 integer,intent(out) :: n_out

!arrays
 real(dp),target,intent(inout) :: set_in(3,n_in)

 interface
   function is_equal(k1,k2)
     use defs_basis
     real(dp),intent(in) :: k1(3),k2(3)
     logical :: is_equal
   end function is_equal
 end interface

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 logical :: isnew
!arrays
 type rdp1d_pt
  integer :: idx
  real(dp),pointer :: rpt(:)
 end type rdp1d_pt
 type(rdp1d_pt),allocatable :: Ap(:)

! *************************************************************************

 ABI_DATATYPE_ALLOCATE(Ap,(n_in))
 Ap(1)%idx = 1
 Ap(1)%rpt => set_in(:,1)

 n_out=1
 do ii=2,n_in

   isnew=.TRUE.
   do jj=1,n_out
     if (is_equal(set_in(:,ii),Ap(jj)%rpt(:))) then
       isnew=.FALSE.
       exit
     end if
   end do

   if (isnew) then
     n_out=n_out+1
     Ap(n_out)%rpt => set_in(:,ii)
     Ap(n_out)%idx = ii
   end if
 end do

 ! The n_out inequivalent items are packed first.
 if (n_out/=n_in) then
   do ii=1,n_out
     jj=Ap(ii)%idx
     set_in(:,ii) = set_in(:,jj)
     !write(std_out,*) Ap(ii)%idx,Ap(ii)%rpt(:)
   end do
 end if

 ABI_DATATYPE_DEALLOCATE(Ap)

end subroutine remove_copies
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/denominator
!! NAME
!!   denominator
!!
!! FUNCTION
!!  Return the denominator of the rational number dd, sign is not considered.
!!
!! INPUTS
!!  dd=The rational number
!!  tolerance=Absolute tolerance
!!
!! OUTPUT
!!  ierr=If /=0  the input number is not rational within the given tolerance.
!!
!! SOURCE

integer function denominator(dd,ierr,tolerance)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ierr
 real(dp),intent(in) :: dd
 real(dp),optional,intent(in) :: tolerance

!Local variables ------------------------------
!scalars
 integer,parameter :: largest_integer = HUGE(1)
 integer :: ii
 real(dp) :: my_tol

!************************************************************************

 ii=1
 my_tol=0.0001 ; if (PRESENT(tolerance)) my_tol=ABS(tolerance)
 do
   if (ABS(dd*ii-NINT(dd*ii))<my_tol) then
     denominator=ii
     ierr=0
     RETURN
   end if
   ! Handle the case in which dd is not rational within my_tol.
   if (ii==largest_integer) then
     denominator=ii
     ierr=-1
     RETURN
   end if
   ii=ii+1
 end do

end function denominator
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/mincm
!! NAME
!!  mincm
!!
!! FUNCTION
!!  Return the minimum common multiple of ii and jj.
!!
!! SOURCE

integer function mincm(ii,jj)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,jj

!************************************************************************

 if (ii==0.or.jj==0) then
   MSG_BUG('ii==0 or jj==0')
 end if

 mincm=MAX(ii,jj)
 do
   if ( ((mincm/ii)*ii)==mincm .and. ((mincm/jj)*jj)==mincm ) RETURN
   mincm=mincm+1
 end do

end function mincm
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/continued_fract
!! NAME
!!  continued_fract
!!
!! FUNCTION
!!  This routine calculates the continued fraction:
!!
!!                        1
!! f(z) =  _______________________________
!!           z - a1 -        b1^2
!!                   _____________________
!!                     z - a2 -    b2^2
!!                             ___________
!!                                z -a3 -    ........
!!
!! INPUTS
!!  nlev=Number of "levels" in the continued fraction.
!!  term_type=Type of the terminator.
!!    0 --> No terminator.
!!   -1 --> Assume constant coefficients for a_i and b_i for i>nlev with a_inf = a(nlev) and b_inf = b(nleb)
!!    1 --> Same as above but a_inf and b_inf are obtained by averaging over the nlev values.
!!  aa(nlev)=Set of a_i coefficients.
!!  bb(nlev)=Set of b_i coefficients.
!!  nz=Number of points on the z-mesh.
!!  zpts(nz)=z-mesh.
!!
!! OUTPUT
!!  spectrum(nz)=Contains f(z) on the input mesh.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine continued_fract(nlev,term_type,aa,bb,nz,zpts,spectrum)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nlev,term_type,nz
!arrays
 real(dp),intent(in) ::  bb(nlev)
 complex(dpc),intent(in) :: aa(nlev)
 complex(dpc),intent(in) :: zpts(nz)
 complex(dpc),intent(out) :: spectrum(nz)

!Local variables ------------------------------
!scalars
 integer :: it
 real(dp) ::  bb_inf,bg,bu,swap
 complex(dpc) :: aa_inf
 character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: div(:),den(:)

!************************************************************************

 ABI_MALLOC(div,(nz))
 ABI_MALLOC(den,(nz))

 select case (term_type)
 case (0) ! No terminator.
   div=czero
 case (-1,1)
   if (term_type==-1) then
     bb_inf=bb(nlev)
     aa_inf=aa(nlev)
   else
     bb_inf=SUM(bb)/nlev
     aa_inf=SUM(aa)/nlev
   end if
   ! Be careful with the sign of the SQRT.
   div(:) = half*(bb(nlev)/(bb_inf))**2 * ( zpts-aa_inf - SQRT((zpts-aa_inf)**2 - four*bb_inf**2) )
 case (2)
   MSG_ERROR("To be tested")
   div = zero
   if (nlev>4) then
     bg=zero; bu=zero
     do it=1,nlev,2
       if (it+2<nlev) bg = bg + bb(it+2)
       bu = bu + bb(it)
     end do
     bg = bg/(nlev/2+MOD(nlev,2))
     bu = bg/((nlev+1)/2)
     !if (iseven(nlev)) then
     if (.not.iseven(nlev)) then
       swap = bg
       bg = bu
       bu = bg
     end if
     !write(std_out,*)nlev,bg,bu
     !Here be careful with the sign of SQRT
     do it=1,nz
       div(it) = half/zpts(it) * (bb(nlev)/bu)**2 * &
&        ( (zpts(it)**2 +bu**2 -bg**2) - SQRT( (zpts(it)**2+bu**2-bg**2)**2 -four*(zpts(it)*bu)**2) )
     end do
   end if

 case default
   write(msg,'(a,i0)')" Wrong value for term_type : ",term_type
   MSG_ERROR(msg)
 end select

 do it=nlev,2,-1
   den(:) = zpts(:) - aa(it) - div(:)
   div(:) = (bb(it-1)**2 )/ den(:)
 end do

 den = zpts(:) - aa(1) - div(:)
 div = one/den(:)

 spectrum = div
 ABI_FREE(div)
 ABI_FREE(den)

end subroutine continued_fract
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/cmplx_sphcart
!! NAME
!! cmplx_sphcart
!!
!! FUNCTION
!! Convert an array of complex values stored in spherical coordinates
!! to Cartesian coordinates with real and imaginary part or vice versa.
!!
!! INPUTS
!!  from=Option specifying the format used to store the complex values. See below.
!!  [units]=Option to specify if angles are given in  "Radians" (default) or "Degrees".
!!
!! SIDE EFFECTS
!!  carr(:,:):
!!    input:  array with complex values in Cartesian form if from="C" or spherical form if from="S"
!!    output: array with values converted to the new representation.
!!
!! PARENTS
!!      m_ptgroups
!!
!! CHILDREN
!!
!! SOURCE

subroutine cmplx_sphcart(carr, from, units)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: from
 character(len=*),optional,intent(in) :: units
!arrays
 complex(dpc),intent(inout) :: carr(:,:)

!Local variables-------------------------------
!scalars
 integer :: jj,ii
 real(dp) :: rho,theta,fact
 character(len=500) :: msg

! *************************************************************************

 select case (from(1:1))

 case ("S","s") ! Spherical --> Cartesian

   fact = one
   if (PRESENT(units)) then
     if (units(1:1) == "D" .or. units(1:1) == "d") fact = two_pi/360_dp
   end if

   do jj=1,SIZE(carr,DIM=2)
     do ii=1,SIZE(carr,DIM=1)
        rho  = DBLE(carr(ii,jj))
        theta= AIMAG(carr(ii,jj)) * fact
        carr(ii,jj) = CMPLX(rho*DCOS(theta), rho*DSIN(theta), kind=dpc)
     end do
   end do

 case ("C","c") ! Cartesian --> Spherical \theta = 2 arctan(y/(rho+x))

   fact = one
   if (PRESENT(units)) then
     if (units(1:1) == "D" .or. units(1:1) == "d") fact = 360_dp/two_pi
   end if

   do jj=1,SIZE(carr,DIM=2)
     do ii=1,SIZE(carr,DIM=1)
        rho  = SQRT(ABS(carr(ii,jj)))
        if (rho > tol16) then
          theta= two * ATAN( AIMAG(carr(ii,jj)) / (DBLE(carr(ii,jj)) + rho) )
        else
          theta= zero
        end if
        carr(ii,jj) = CMPLX(rho, theta*fact, kind=dpc)
     end do
   end do

 case default
   msg = " Wrong value for from: "//TRIM(from)
   MSG_BUG(msg)
 end select

end subroutine cmplx_sphcart
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/pfactorize
!! NAME
!!  pfactorize
!!
!! FUNCTION
!!  Factorize a number in terms of an user-specified set of prime factors
!!  nn = alpha * Prod_i p^i   1)
!!
!! INPUTS
!!  nn=The number to be factorized.
!!  nfactors=The number of factors
!!  pfactors(nfactors)=The list of prime number e.g. (/ 2, 3, 5, 7, 11 /)
!!
!! OUTPUT
!!  powers(nfactors+1)=
!!   The first nfactors entries are the powers i in Eq.1. powers(nfactors+1) is alpha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pfactorize(nn,nfactors,pfactors,powers)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn,nfactors
 integer,intent(in) :: pfactors(nfactors)
 integer,intent(out) :: powers (nfactors+1)

!Local variables ------------------------------
!scalars
 integer :: tnn,ifc,fact,ipow,maxpwr

! *************************************************************************

 powers=0; tnn=nn

 fact_loop: do ifc=1,nfactors
   fact = pfactors (ifc)
   maxpwr = NINT ( LOG(DBLE(tnn))/LOG(DBLE(fact) ) ) + 1
   do ipow=1,maxpwr
     if (tnn==1) EXIT fact_loop
     if ( MOD(tnn,fact)==0 ) then
       tnn=tnn/fact
       powers(ifc)=powers(ifc) + 1
     end if
   end do
 end do fact_loop

 if ( nn /= tnn * PRODUCT( pfactors**powers(1:nfactors)) ) then
   MSG_BUG('nn/=tnn!')
 end if

 powers(nfactors+1) = tnn

end subroutine pfactorize
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/isordered
!! NAME
!!  isordered
!!
!! FUNCTION
!!  Return .TRUE. if values in array arr are ordered.
!!  Consider that two double precision numbers within tolerance tol are equal.
!!
!! INPUTS
!!  nn=Size of arr.
!!  arr(nn)=The array with real values to be tested.
!!  direction= ">" for ascending numerical order.
!!             ">" for decreasing numerical order.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function isordered_rdp(nn,arr,direction,tol) result(isord)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp),intent(in) :: tol
 logical :: isord
 character(len=*),intent(in) :: direction
!arrays
 real(dp),intent(in) :: arr(nn)

!Local variables ------------------------------
!scalars
 integer :: ii
 real(dp) :: prev
 character(len=500) :: msg

! *************************************************************************

 prev = arr(1); isord =.TRUE.

 SELECT CASE (direction(1:1))
 CASE(">")
 ii=2;
 do while (ii<=nn .and. isord)
   if (ABS(arr(ii)-prev) > tol) isord = (arr(ii) >= prev)
   prev = arr(ii)
   ii = ii +1
 end do

 CASE("<")
 ii=2;
 do while (ii<=nn .and. isord)
   if (ABS(arr(ii)-prev) > tol) isord = (arr(ii) <= prev)
   prev = arr(ii)
   ii = ii +1
 end do

 CASE DEFAULT
   msg = "Wrong direction: "//TRIM(direction)
   MSG_ERROR(msg)
 END SELECT

end function isordered_rdp
!!***

!----------------------------------------------------------------------


!!****f* m_numeric_tools/stats_eval
!! NAME
!!  stats_eval
!!
!! FUNCTION
!!  Helper function used to calculate the statistical parameters of a data set.
!!
!! INPUT
!!  arr(:)=Array with the values.
!!
!! OUTPUT
!!  stats<stats_t>=Data type storing the parameters of the data set.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!
!! SOURCE

pure function stats_eval(arr) result(stats)

!Arguments ------------------------------------
!scalars
 type(stats_t) :: stats
!arrays
 real(dp),intent(in) :: arr(:)

!Local variables ------------------------------
!scalars
 integer :: ii,nn
 real(dp) :: xx,x2_sum

! *************************************************************************

 !@stats_t
 stats%min   = +HUGE(one)
 stats%max   = -HUGE(one)
 stats%mean  = zero

 nn = SIZE(arr)

 do ii=1,nn
   xx = arr(ii)
   stats%max  = MAX(stats%max, xx)
   stats%min  = MIN(stats%min, xx)
   stats%mean = stats%mean + xx
 end do

 stats%mean  = stats%mean/nn

 ! Two-pass algorithm for the variance (more stable than the single-pass one).
 x2_sum = zero
 do ii=1,nn
   xx     = arr(ii)
   x2_sum = x2_sum + (xx - stats%mean)*(xx - stats%mean)
 end do

 if (nn>1) then
   stats%stdev  = x2_sum/(nn-1)
   stats%stdev = SQRT(ABS(stats%stdev))
 else
   stats%stdev = zero
 end if

end function stats_eval
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/wrap2_zero_one
!! NAME
!! wrap2_zero_one
!!
!! FUNCTION
!! Transforms a real number (num) in its corresponding reduced number
!! (red) in the interval [0,1[ where 1 is not included (tol12)
!! num=red+shift
!!
!! INPUTS
!!  num=real number
!!
!! OUTPUT
!! red=reduced number of num in the interval [0,1[ where 1 is not included
!! shift=num-red
!!
!! PARENTS
!!      exc_plot,k_neighbors,lin_interpq_gam,m_nesting,m_paw_pwaves_lmn
!!      pawmkaewf
!!
!! CHILDREN
!!
!! SOURCE

elemental subroutine wrap2_zero_one(num, red, shift)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: num
 real(dp),intent(out) :: red,shift

! *************************************************************************

 if (num>zero) then
   red=mod((num+tol12),one)-tol12
 else
   red=-mod(-(num-one+tol12),one)+one-tol12
 end if
 if(abs(red)<tol12)red=0.0_dp
 shift=num-red

end subroutine wrap2_zero_one
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/wrap2_pmhalf
!! NAME
!! wrap2_pmhalf
!!
!! FUNCTION
!! Transforms a real number (num) in its corresponding reduced number
!! (red) in the interval ]-1/2,1/2] where -1/2 is not included (tol12)
!! num=red+shift
!!
!! INPUTS
!! num=real number
!!
!! OUTPUT
!! red=reduced number of num in the interval ]-1/2,1/2] where -1/2 is not included
!! shift=num-red
!!
!! PARENTS
!!      canat9,dist2,elphon,ep_setupqpt,get_full_kgrid,interpolate_gkk
!!      m_bz_mesh,m_cprj_bspline,m_gamma,m_io_gkk,m_optic_tools,m_shirley
!!      mkfskgrid,mkfsqgrid,mkph_linwid,mkphbs,order_fs_kpts,printvtk,read_gkk
!!      smpbz,littlegroup_q
!!
!! CHILDREN
!!
!! SOURCE

elemental subroutine wrap2_pmhalf(num,red,shift)

!Arguments -------------------------------
!scalars
 real(dp),intent(in) :: num
 real(dp),intent(out) :: red,shift

! *********************************************************************

 if (num>zero) then
   red=mod((num+half-tol12),one)-half+tol12
 else
   red=-mod(-(num-half-tol12),one)+half+tol12
 end if
 if(abs(red)<tol12)red=0.0d0
 shift=num-red

end subroutine wrap2_pmhalf
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/interpol3d
!! NAME
!! interpol3d
!!
!! FUNCTION
!! Computes the value at any point r by linear interpolation
!! inside the eight vertices of the surrounding cube
!! r is presumed to be normalized, in a unit cube for the full grid
!!
!! INPUTS
!! r(3)=point coordinate
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! grid(nr1,nr2,nr3)=grid matrix
!!
!! OUTPUT
!! res=Interpolated value
!!
!! PARENTS
!!      integrate_gamma_alt,lin_interpq_gam,lineint,m_nesting,m_qparticles
!!      planeint,pointint,volumeint
!!
!! CHILDREN
!!
!! SOURCE

pure function interpol3d(r, nr1, nr2, nr3, grid) result(res)

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1, nr2, nr3
 real(dp) :: res
!arrays
 real(dp),intent(in) :: grid(nr1,nr2,nr3),r(3)

!Local variables--------------------------------------------------------
!scalars
 integer :: ir1,ir2,ir3,pr1,pr2,pr3
 real(dp) :: x1,x2,x3

! *************************************************************************

 call interpol3d_indices (r,nr1,nr2,nr3,ir1,ir2,ir3, pr1,pr2,pr3)

!weight
 x1=one+r(1)*nr1-real(ir1)
 x2=one+r(2)*nr2-real(ir2)
 x3=one+r(3)*nr3-real(ir3)

!calculation of the density value
 res=zero
 res=res + grid(ir1, ir2, ir3) * (one-x1)*(one-x2)*(one-x3)
 res=res + grid(pr1, ir2, ir3) * x1*(one-x2)*(one-x3)
 res=res + grid(ir1, pr2, ir3) * (one-x1)*x2*(one-x3)
 res=res + grid(ir1, ir2, pr3) * (one-x1)*(one-x2)*x3
 res=res + grid(pr1, pr2, ir3) * x1*x2*(one-x3)
 res=res + grid(ir1, pr2, pr3) * (one-x1)*x2*x3
 res=res + grid(pr1, ir2, pr3) * x1*(one-x2)*x3
 res=res + grid(pr1, pr2, pr3) * x1*x2*x3

end function interpol3d
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/interpol3d_indices
!! NAME
!! interpol3d_indices
!!
!! FUNCTION
!! Computes the indices in a cube which are neighbors to the point to be
!!  interpolated in interpol3d
!!
!! INPUTS
!! r(3)=point coordinate
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!!
!! OUTPUT
!! ir1,ir2,ir3 = bottom left neighbor
!! pr1,pr2,pr3 = top right neighbor
!!
!! PARENTS
!!      interpol3d,k_neighbors
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine interpol3d_indices (r,nr1,nr2,nr3,ir1,ir2,ir3,pr1,pr2,pr3)

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1,nr2,nr3
 integer,intent(out) :: ir1,ir2,ir3
 integer,intent(out) :: pr1,pr2,pr3
!arrays
 real(dp),intent(in) :: r(3)

!Local variables-------------------------------
 real(dp) :: d1,d2,d3

! *************************************************************************

!grid density
 d1=one/nr1
 d2=one/nr2
 d3=one/nr3

!lower left
 ir1=int(r(1)/d1)+1
 ir2=int(r(2)/d2)+1
 ir3=int(r(3)/d3)+1

!upper right
 pr1=mod(ir1+1,nr1)
 pr2=mod(ir2+1,nr2)
 pr3=mod(ir3+1,nr3)

 if(ir1==0) ir1=nr1
 if(ir2==0) ir2=nr2
 if(ir3==0) ir3=nr3

 if(ir1>nr1) ir1=ir1-nr1
 if(ir2>nr2) ir2=ir2-nr2
 if(ir3>nr3) ir3=ir3-nr3

 if(pr1==0) pr1=nr1
 if(pr2==0) pr2=nr2
 if(pr3==0) pr3=nr3

end subroutine interpol3d_indices
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/interpolate_denpot
!! NAME
!! interpolate_denpot
!!
!! FUNCTION
!!  Linear interpolation of density/potential given on the real space FFT mesh.
!!  Assumes array on full mesh i.e. no MPI-FFT.
!!
!! INPUTS
!!  cplex=1 for real, 2 for complex data.
!!  in_ngfft(3)=Mesh divisions of input array
!!  nspden=Number of density components.
!!  in_rhor(cplex * in_nfftot * nspden)=Input array
!!  out_ngfft(3)=Mesh divisions of output array
!!
!! OUTPUT
!!  outrhor(cplex * out_nfftot * nspden)=Output array with interpolated data.
!!
!! PARENTS
!!      m_ioarr
!!
!! CHILDREN
!!
!! SOURCE

subroutine interpolate_denpot(cplex, in_ngfft, nspden, in_rhor, out_ngfft, out_rhor)

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden
!arrays
 integer,intent(in) :: in_ngfft(3), out_ngfft(3)
 real(dp),intent(in) :: in_rhor(cplex, product(in_ngfft), nspden)
 real(dp),intent(out) :: out_rhor(cplex, product(out_ngfft), nspden)

!Local variables--------------------------------------------------------
!scalars
 integer :: ispden, ir1, ir2, ir3, ifft
 real(dp) :: rr(3)
 real(dp),allocatable :: re(:,:),im(:,:)

! *************************************************************************

 if (cplex == 2) then
   ! copy slices for efficiency reasons (the best would be to have stride option in interpol3d)
   ABI_MALLOC(re, (product(in_ngfft), nspden))
   ABI_MALLOC(im, (product(in_ngfft), nspden))
   re = in_rhor(1, :, :)
   im = in_rhor(2, :, :)
 end if

 ! Linear interpolation.
 do ispden=1,nspden
   do ir3=0,out_ngfft(3)-1
     rr(3) = DBLE(ir3)/out_ngfft(3)
     do ir2=0,out_ngfft(2)-1
       rr(2) = DBLE(ir2)/out_ngfft(2)
       do ir1=0,out_ngfft(1)-1
         rr(1) = DBLE(ir1)/out_ngfft(1)
         ifft = 1 + ir1 + ir2*out_ngfft(1) + ir3*out_ngfft(1)*out_ngfft(2)
         if (cplex == 1) then
           out_rhor(1, ifft, ispden) = interpol3d(rr, in_ngfft(1), in_ngfft(2), in_ngfft(3), in_rhor(1, :, ispden))
         else
           out_rhor(1, ifft, ispden) = interpol3d(rr, in_ngfft(1), in_ngfft(2), in_ngfft(3), re(:, ispden))
           out_rhor(2, ifft, ispden) = interpol3d(rr, in_ngfft(1), in_ngfft(2), in_ngfft(3), im(:, ispden))
         end if
       end do
     end do
   end do
 end do

 if (cplex == 2) then
   ABI_FREE(re)
   ABI_FREE(im)
 end if

end subroutine interpolate_denpot
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/simpson_int
!! NAME
!! simpson_int
!!
!! FUNCTION
!!   Simpson integral of input function
!!
!! INPUTS
!!  npts=max number of points on grid for integral
!!  step = space between integral arguments
!!  values(npts)=integrand function.
!!
!! OUTPUT
!!  int_values(npts)=integral of values.
!!
!! PARENTS
!!      m_a2ftr,m_ebands,m_eliashberg_1d,m_elphon,m_evdw_wannier,m_exc_spectra
!!      m_integrals,m_mlwfovlp,m_numeric_tools,m_outscfcv,m_phgamma,m_phonons
!!      m_rta,m_tdep_psij,m_xc_vdw
!!
!! CHILDREN
!!
!! SOURCE

subroutine simpson_int(npts, step, values, int_values)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
 real(dp),intent(in) :: step
!arrays
 real(dp),intent(in) :: values(npts)
 real(dp),intent(out) :: int_values(npts)

!Local variables -------------------------
!scalars
 integer :: ii
 real(dp),parameter :: coef1 = 0.375_dp                          !9.0_dp  / 24.0_dp
 real(dp),parameter :: coef2 = 1.166666666666666666666666667_dp  !28.0_dp / 24.0_dp
 real(dp),parameter :: coef3 = 0.958333333333333333333333333_dp  !23.0_dp / 24.0_dp
 character(len=500) :: msg

! *********************************************************************

 if (npts < 6) then
   write(msg,"(a,i0)")"Number of points in integrand function must be >=6 while it is: ",npts
   MSG_ERROR(msg)
 end if

!-----------------------------------------------------------------
!Simpson integral of input function
!-----------------------------------------------------------------

!first point is 0: don t store it
!do integration equivalent to Simpson O(1/N^4) from NumRec in C p 134  NumRec in Fortran p 128
 int_values(1) =               coef1*values(1)
 int_values(2) = int_values(1) + coef2*values(2)
 int_values(3) = int_values(2) + coef3*values(3)

 do ii=4,npts-3
   int_values(ii) = int_values(ii-1) + values(ii)
 end do

 int_values(npts-2) = int_values(npts-3) + coef3*values(npts-2)
 int_values(npts-1) = int_values(npts-2) + coef2*values(npts-1)
 int_values(npts  ) = int_values(npts-1) + coef1*values(npts  )

 int_values(:) = int_values(:) * step

end subroutine simpson_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/simpson
!! NAME
!! simpson
!!
!! FUNCTION
!!   Simpson integral of input function
!!
!! INPUTS
!!  step = space between integral arguments
!!  values(npts)=integrand function.
!!
!! OUTPUT
!!  integral of values on the full mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function simpson(step, values) result(res)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: step
 real(dp) :: res
!arrays
 real(dp),intent(in) :: values(:)

!Local variables -------------------------
!scalars
 real(dp) :: int_values(size(values))

! *********************************************************************

 call simpson_int(size(values),step,values,int_values)
 res = int_values(size(values))

end function simpson
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/rhophi
!! NAME
!! rhophi
!!
!! FUNCTION
!! Compute the phase and the module of a complex number.
!! The phase angle is fold into the interval [-pi,pi]
!!
!! INPUTS
!!  cx(2) = complex number
!!
!! OUTPUT
!!  phi = phase of cx fold into [-pi,pi]
!!  rho = modul of cx
!!
!! PARENTS
!!      berryphase_new,etheta,linemin
!!
!! SOURCE

pure subroutine rhophi(cx,phi,rho)

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: phi,rho
!arrays
 real(dp),intent(in) :: cx(2)


! ***********************************************************************

 rho = sqrt(cx(1)*cx(1) + cx(2)*cx(2))

 if (abs(cx(1)) > tol8) then

   phi = atan(cx(2)/cx(1))

!  phi is an element of [-pi,pi]
   if (cx(1) < zero) then
     if (phi < zero) then
       phi = phi + pi
     else
       phi = phi - pi
     end if
   end if

 else

   if (cx(2) > tol8) then
     phi = pi*half
   else if (cx(2) < tol8) then
     phi = -0.5_dp*pi
   else
     phi = 0
   end if

 end if

end subroutine rhophi
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/vdiff_eval
!! NAME
!! vdiff_eval
!!
!! FUNCTION
!!  Estimate the "distance" between two functions tabulated on a homogeneous grid.
!!  See vdiff_t
!!
!! INPUTS
!!  cplex=1 if f1 and f2 are real, 2 for complex.
!!  nr=Number of points in the mesh.
!!  f1(cplex,nr), f2(cplex,nr)=Vectors with values
!!  [vd_max]= Compute max value of the different entries.
!!
!! OUTPUT
!!  vdiff_t object
!!
!! PARENTS
!!
!! SOURCE

type(vdiff_t) function vdiff_eval(cplex, nr, f1, f2, volume, vd_max) result(vd)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nr
 real(dp),intent(in) :: volume
 type(vdiff_t),optional,intent(inout) :: vd_max
!arrays
 real(dp),intent(in) :: f1(cplex,nr),f2(cplex,nr)

!Local variables-------------------------------
!scalars
 integer :: ir
 real(dp) :: num,den,dr
 type(stats_t) :: stats
!arrays
 real(dp) :: abs_diff(nr)
! *********************************************************************

 dr = volume / nr

 if (cplex == 1) then
   abs_diff = abs(f1(1,:) - f2(1,:))
   num = sum(abs_diff)
   den = sum(abs(f2(1,:)))

 else if (cplex == 2) then
   do ir=1,nr
     abs_diff(ir) = sqrt((f1(1,ir) - f2(1,ir))**2 + (f1(2,ir) - f2(2,ir))**2)
   end do
   num = sum(abs_diff)
   den = zero
   do ir=1,nr
     den = den + sqrt(f2(1,ir)**2 + f2(2,ir)**2)
   end do
 end if

 vd%int_adiff = num * dr
 call safe_div(num,den,zero,vd%l1_rerr)

 stats = stats_eval(abs_diff)
 vd%mean_adiff = stats%mean
 vd%stdev_adiff = stats%stdev
 vd%min_adiff = stats%min
 vd%max_adiff = stats%max

 if (present(vd_max)) then
   vd_max%int_adiff =   max(vd_max%int_adiff, vd%int_adiff)
   vd_max%mean_adiff =  max(vd_max%mean_adiff, vd%mean_adiff)
   vd_max%stdev_adiff = max(vd_max%stdev_adiff, vd%stdev_adiff)
   vd_max%min_adiff =   max(vd_max%min_adiff, vd%min_adiff)
   vd_max%max_adiff =   max(vd_max%max_adiff, vd%max_adiff)
   vd_max%l1_rerr =     max(vd_max%l1_rerr, vd%L1_rerr)
 end if

end function vdiff_eval
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/vdiff_print
!! NAME
!! vdiff_print
!!
!! FUNCTION
!!  Print vdiff_t to unit
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine vdiff_print(vd, unit)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit
 type(vdiff_t),intent(in) :: vd

!Local variables-------------------------------
!scalars
 integer :: unt
! *********************************************************************

 unt = std_out; if (present(unit)) unt = unit
 write(unt,"(a,es10.3,a)")"  L1_rerr: ", vd%l1_rerr, ","
 write(unt,"(a,es10.3,a)")"  'Integral |f1-f2|dr': ", vd%int_adiff, ","
 write(unt,"(a,es10.3,a)")"  'min {|f1-f2|}': ", vd%min_adiff, ","
 write(unt,"(a,es10.3,a)")"  'Max {|f1-f2|}': ", vd%max_adiff, ","
 write(unt,"(a,es10.3,a)")"  'mean {|f1-f2|}': ", vd%mean_adiff, ","
 write(unt,"(a,es10.3,a)")"  'stdev {|f1-f2|}': ", vd%stdev_adiff, ","

end subroutine vdiff_print
!!***

!!****f* m_numeric_tools/smooth
!! NAME
!!  smooth data.
!!
!! FUNCTION
!!  smooth
!!
!! INPUTS
!!  mesh=Number of points.
!!  it=Number of iterations.
!!
!! SIDE EFFECTS
!!  a(mesh)=Input values, smoothed in output
!!
!! PARENTS
!!      m_psp6,m_upf2abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine smooth(a,mesh,it)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: it,mesh
 real(dp), intent(inout) :: a(mesh)
!Local variables-------------------------------
 integer :: i,k
 real(dp) :: asm(mesh)
! *********************************************************************

 do k=1,it
   asm(1)=1.0d0/3.0d0*(a(1)+a(2)+a(3))
   asm(2)=0.25d0*(a(1)+a(2)+a(3)+a(4))
   asm(3)=0.2d0*(a(1)+a(2)+a(3)+a(4)+a(5))
   asm(4)=0.2d0*(a(2)+a(3)+a(4)+a(5)+a(6))
   asm(5)=0.2d0*(a(3)+a(4)+a(5)+a(6)+a(7))
   asm(mesh-4)=0.2d0*(a(mesh-2)+a(mesh-3)+a(mesh-4)+&
&                   a(mesh-5)+a(mesh-6))
   asm(mesh-3)=0.2d0*(a(mesh-1)+a(mesh-2)+a(mesh-3)+&
&                   a(mesh-4)+a(mesh-5))
   asm(mesh-2)=0.2d0*(a(mesh)+a(mesh-1)+a(mesh-2)+&
&                   a(mesh-3)+a(mesh-4))
   asm(mesh-1)=0.25d0*(a(mesh)+a(mesh-1)+a(mesh-2)+a(mesh-3))
   asm(mesh)=1.0d0/3.0d0*(a(mesh)+a(mesh-1)+a(mesh-2))

   do i=6,mesh-5
     asm(i)=0.1d0*a(i)+0.1d0*(a(i+1)+a(i-1))+&
&             0.1d0*(a(i+2)+a(i-2))+&
&             0.1d0*(a(i+3)+a(i-3))+&
&             0.1d0*(a(i+4)+a(i-4))+&
&             0.05d0*(a(i+5)+a(i-5))
   end do

   do i=1,mesh
     a(i)=asm(i)
   end do
 end do

end subroutine smooth
!!***

!!****f* m_numeric_tools/nderiv
!! NAME
!! nderiv
!!
!! FUNCTION
!! Given an input function y(x) on a regular grid,
!! compute its first or second derivative.
!!
!! INPUTS
!!  hh= radial step
!!  ndim= radial mesh size
!!  yy(ndim)= input function
!!  norder= order of derivation (1 or 2)
!!
!! OUTPUT
!!  zz(ndim)= first or second derivative of y
!!
!! PARENTS
!!      m_upf2abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine nderiv(hh,yy,zz,ndim,norder)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndim,norder
 real(dp),intent(in) :: hh
!arrays
 real(dp),intent(in) :: yy(ndim)
 real(dp),intent(out) :: zz(ndim)

!Local variables ---------------------------------------
!scalars
 integer :: ier,ii
 real(dp) :: aa,bb,cc,h1,y1

! *********************************************************************

!Initialization (common to 1st and 2nd derivative)
 h1=one/(12.d0*hh)
 y1=yy(ndim-4)

!FIRST DERIVATIVE
!================
 if (norder==1) then

!  Prepare differentiation loop
   bb=h1*(-25.d0*yy(1)+48.d0*yy(2)-36.d0*yy(3)+16.d0*yy(4)-3.d0*yy(5))
   cc=h1*(-3.d0*yy(1)-10.d0*yy(2)+18.d0*yy(3)-6.d0*yy(4)+yy(5))
!  Start differentiation loop
   do ii=5,ndim
     aa=bb;bb=cc
     cc=h1*(yy(ii-4)-yy(ii)+8.d0*(yy(ii-1)-yy(ii-3)))
     zz(ii-4)=aa
   end do
!  Normal exit
   ier=0
   aa=h1*(-y1+6.d0*yy(ndim-3)-18.d0*yy(ndim-2)+10.d0*yy(ndim-1)+3.d0*yy(ndim))
   zz(ndim)=h1*(3.d0*y1-16.d0*yy(ndim-3)+36.d0*yy(ndim-2) -48.d0*yy(ndim-1)+25.d0*yy(ndim))
   zz(ndim-1)=aa
   zz(ndim-2)=cc
   zz(ndim-3)=bb

!  SECOND DERIVATIVE
!  =================
 else
   h1=h1/hh
!  Prepare differentiation loop
   bb=h1*(35.d0*yy(1)-104.d0*yy(2)+114.d0*yy(3)-56.d0*yy(4)+11.d0*yy(5))
   cc=h1*(11.d0*yy(1)-20.d0*yy(2)+6.d0*yy(3)+4.d0*yy(4)-yy(5))
!  Start differentiation loop
   do ii=5,ndim
     aa=bb;bb=cc
     cc=h1*(-yy(ii-4)-yy(ii)+16.d0*(yy(ii-1)+yy(ii-3))-30.d0*yy(ii-2))
     zz(ii-4)=aa
   end do
!  Normal exit
   ier=0
   aa=h1*(-y1+4.d0*yy(ndim-3)+6.d0*yy(ndim-2)-20.d0*yy(ndim-1)+11.d0*yy(ndim))
   zz(ndim)=h1*(11.d0*y1-56.d0*yy(ndim-3)+114.d0*yy(ndim-2) -104.d0*yy(ndim-1)+35.d0*yy(ndim))
   zz(ndim-1)=aa
   zz(ndim-2)=cc
   zz(ndim-3)=bb

 end if !norder

end subroutine nderiv
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/central_finite_diff
!! NAME
!! central_finite_diff
!!
!! FUNCTION
!! Coefficients of the central differences, for several orders of accuracy.
!! See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
!!
!! INPUTS
!!  order=Derivative order.
!!  ipos=Index of the point must be in [1,npts]
!!  npts=Number of points used in finite difference, origin at npts/2 + 1
!!
!! OUTPUT
!!  coeffient for central finite difference
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

real(dp) function central_finite_diff(order, ipos, npts) result(fact)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipos,order,npts

!Local variables ---------------------------------------
!scalars
 real(dp),parameter :: empty=huge(one)
! 1st derivative.
 real(dp),parameter :: d1(9,4) = reshape([ &
  [-1/2._dp, 0._dp, 1/2._dp, empty, empty, empty, empty, empty, empty], &
  [ 1/12._dp, -2/3._dp, 0._dp, 2/3._dp, -1/12._dp, empty, empty, empty, empty], &
  [-1/60._dp, 3/20._dp, -3/4._dp, 0._dp, 3/4._dp, -3/20._dp, 1/60._dp, empty, empty], &
  [ 1/280._dp, -4/105._dp, 1/5._dp, -4/5._dp, 0._dp, 4/5._dp, -1/5._dp, 4/105._dp, -1/280._dp]], [9,4])
! 2nd derivative.
 real(dp),parameter :: d2(9,4) = reshape([ &
   [ 1._dp, -2._dp, 1._dp, empty, empty, empty, empty, empty, empty], &
   [-1/12._dp, 4/3._dp, -5/2._dp, 4/3._dp, -1/12._dp, empty, empty, empty, empty], &
   [ 1/90._dp, -3/20._dp, 3/2._dp, -49/18._dp, 3/2._dp, -3/20._dp, 1/90._dp, empty, empty], &
   [-1/560._dp, 8/315._dp, -1/5._dp, 8/5._dp, -205/72._dp, 8/5._dp, -1/5._dp, 8/315._dp, -1/560._dp]], [9,4])
! 3th derivative.
 real(dp),parameter :: d3(9,3) = reshape([ &
   [-1/2._dp, 1._dp, 0._dp, -1._dp, 1/2._dp, empty, empty, empty, empty], &
   [ 1/8._dp, -1._dp, 13/8._dp, 0._dp, -13/8._dp, 1._dp, -1/8._dp, empty, empty], &
   [ -7/240._dp, 3/10._dp, -169/120._dp, 61/30._dp, 0._dp, -61/30._dp, 169/120._dp, -3/10._dp, 7/240._dp]], &
   [9,3])
! 4th derivative.
 real(dp),parameter :: d4(9,3) = reshape([ &
   [ 1._dp, -4._dp, 6._dp, -4._dp, 1._dp, empty, empty, empty, empty], &
   [ -1/6._dp, 2._dp, -13/2._dp, 28/3._dp, -13/2._dp, 2._dp, -1/6._dp, empty, empty], &
   [ 7/240._dp, -2/5._dp, 169/60._dp, -122/15._dp, 91/8._dp, -122/15._dp, 169/60._dp, -2/5._dp, 7/240._dp]], [9,3])
! 5th derivative.
 real(dp),parameter :: d5(7) = [ -1/2._dp, 2._dp, -5/2._dp, 0._dp, 5/2._dp, -2._dp, 1/2._dp]
! 6th derivative.
 real(dp),parameter :: d6(7) = [ 1._dp, -6._dp, 15._dp, -20._dp, 15._dp, -6._dp, 1._dp]
! *********************************************************************

 select case (order)
 case (1)
   if (ipos < 1 .or. ipos > 9 .or. npts < 1 .or. npts > 9) goto 10
   fact = d1(ipos, npts/2)
 case (2)
   if (ipos < 1 .or. ipos > 9 .or. npts < 1 .or. npts > 9) goto 10
   fact = d2(ipos, npts/2)
 case (3)
   if (ipos < 1 .or. ipos > 9 .or. npts < 1 .or. npts > 9) goto 10
   fact = d3(ipos, npts/2)
 case (4)
   if (ipos < 1 .or. ipos > 9 .or. npts < 1 .or. npts > 9) goto 10
   fact = d4(ipos, npts/2 - 1)
 case (5)
   if (ipos < 1 .or. ipos > 7 .or. npts /= 7) goto 10
   fact = d5(ipos)
 case (6)
   if (ipos < 1 .or. ipos > 7 .or. npts /= 7) goto 10
   fact = d6(ipos)
 case default
   MSG_ERROR(sjoin("No entry for ipos:",itoa(ipos),"order", itoa(order), "npts", itoa(npts)))
 end select

 if (fact == empty) then
   MSG_ERROR(sjoin("Invalid ipos:",itoa(ipos),"for order", itoa(order), "npts", itoa(npts)))
 end if
 return

10 MSG_ERROR(sjoin("No entry for ipos:",itoa(ipos),"order", itoa(order), "npts", itoa(npts)))

end function central_finite_diff
!!***

!!****f* m_numeric_tools/uniformrandom
!! NAME
!!  uniformrandom
!!
!! FUNCTION
!! Returns a uniform random deviate between 0.0 and 1.0.
!! Set seed to any value < 0 to initialize or reinitialize sequence.
!! Parameters are chosen from integer overflow=2**23 (conservative).
!! For some documentation, see Numerical Recipes, 1986, p196.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function uniformrandom(seed)

!Arguments ------------------------------------
!scalars
 real(dp) :: uniformrandom
 integer,intent(inout) :: seed

!Local variables ---------------------------------------
 integer, parameter :: im1=11979,ia1= 430,ic1=2531
 integer, parameter :: im2= 6655,ia2= 936,ic2=1399
 integer, parameter :: im3= 6075,ia3=1366,ic3=1283
 integer, save :: init=0
 integer, save :: ii1,ii2,ii3
 integer :: kk
 real(dp) :: im1inv,im2inv
 real(dp), save :: table(97)
 character(len=500) :: message

! *********************************************************************

 im1inv=1.0d0/im1 ; im2inv=1.0d0/im2

!Initialize on first call or when seed<0:
 if (seed<0.or.init==0) then
   seed=-abs(seed)

!  First generator
   ii1=mod(ic1-seed,im1)
   ii1=mod(ia1*ii1+ic1,im1)
!  Second generator
   ii2=mod(ii1,im2)
   ii1=mod(ia1*ii1+ic1,im1)
!  Third generator
   ii3=mod(ii1,im3)

!  Fill table
   do kk=1,97
     ii1=mod(ia1*ii1+ic1,im1)
     ii2=mod(ia2*ii2+ic2,im2)
     table(kk)=(dble(ii1)+dble(ii2)*im2inv)*im1inv
   enddo

   init=1 ; seed=1
 end if

!Third generator gives index
 ii3=mod(ia3*ii3+ic3,im3)
 kk=1+(97*ii3)/im3
 if (kk<1.or.kk>97) then
   write(message,'(a,2i0,a)' ) ' trouble in uniformrandom; ii3,kk=',ii3,kk,' =>stop'
   MSG_ERROR(message)
 end if
 uniformrandom=table(kk)

!Replace old value, based on generators 1 and 2
 ii1=mod(ia1*ii1+ic1,im1)
 ii2=mod(ia2*ii2+ic2,im2)
 table(kk)=(dble(ii1)+dble(ii2)*im2inv)*im1inv

end function uniformrandom
!!***

!!****f* m_numeric_tools/findmin
!!
!! NAME
!! findmin
!!
!! FUNCTION
!! Compute the minimum of a function whose value and derivative are known at two points.
!! Also deduce different quantities at this predicted point, and at the two other points
!! It uses a quartic interpolation, with the supplementary
!! condition that the second derivative vanishes at one and
!! only one point. See Schlegel, J. Comp. Chem. 3, 214 (1982) [[cite:Schlegel1982]].
!! For this option, lambda_1 must be 1 (new point),
!! and lambda_2 must be 0 (old point).
!! Also, if the derivative at the new point is more negative
!! than the derivative at the old point, the predicted
!! point cannot correspond to a minimum, but will be lambda=2.5_dp,
!! if the energy of the second point is lower than the energy
!! of the first point.
!!
!! INPUTS
!! etotal_1=first value of the function
!! etotal_2=second value of the function
!! dedv_1=first value of the derivative
!! dedv_2=second value of the derivative
!! lambda_1=first value of the argument
!! lambda_2=second value of the argument
!!
!! OUTPUT
!! dedv_predict=predicted value of the derivative (usually zero,
!!  except if choice=4, if it happens that a minimum cannot be located,
!!  and a trial step is taken)
!! d2edv2_predict=predicted value of the second derivative (not if choice=4)
!! d2edv2_1=first value of the second derivative (not if choice=4)
!! d2edv2_2=second value of the second derivative (not if choice=4)
!! etotal_predict=predicted value of the function
!! lambda_predict=predicted value of the argument
!! status= 0 if everything went normally ;
!!         1 if negative second derivative
!!         2 if some other problem
!!
!! PARENTS
!!      m_bfgs
!!
!! CHILDREN
!!
!! SOURCE

subroutine findmin(dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: status
 real(dp),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(dp),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(dp),intent(out) :: etotal_predict,lambda_predict

!Local variables-------------------------------
!scalars
 real(dp) :: aa,bb,bbp,cc,ccp,d_lambda,dd
 real(dp) :: discr,ee,eep,lambda_shift,sum1,sum2,sum3,uu
 real(dp) :: uu3,vv,vv3
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' findmin : enter'
!write(std_out,*)' choice,lambda_1,lambda_2=',choice,lambda_1,lambda_2
!ENDDEBUG

 status=0
 d_lambda=lambda_1-lambda_2

!DEBUG
!do choice=3,1,-1
!ENDDEBUG

 if(abs(lambda_1-1.0_dp)>tol12 .or. abs(lambda_2)>tol12) then
   message = '  For choice=4, lambda_1 must be 1 and lambda_2 must be 0.'
   MSG_BUG(message)
 end if

!Evaluate quartic interpolation
!etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3 + ee * lambda**4
!Impose positive second derivative everywhere, with
!one point where it vanishes :  3*dd**2=8*cc*ee
 aa=etotal_2
 bb=dedv_2
 sum1=etotal_1-aa-bb
 sum2=dedv_1-bb
 sum3=sum2-2.0_dp*sum1

!Build the discriminant of the associated 2nd degree equation
 discr=sum2**2-3.0_dp*sum3**2
 if(discr<0.0_dp .or. sum2<0.0_dp)then

! jmb init
   d2edv2_2=0.0
   d2edv2_1=0.0
   d2edv2_predict=0.0

!  Even if there is a problem, try to keep going ...
   message = 'The 2nd degree equation has no positive root (choice=4).'
   MSG_WARNING(message)
   status=2
   if(etotal_1<etotal_2)then
     write(message, '(a,a,a)' )&
&     'Will continue, since the new total energy is lower',ch10,&
&     'than the old. Take a larger step in the same direction.'
     MSG_COMMENT(message)
     lambda_predict=2.5_dp
   else
     write(message, '(a,a,a,a,a)' )&
&     'There is a problem, since the new total energy is larger',ch10,&
&     'than the old (choice=4).',ch10,&
&     'I take a point between the old and new, close to the old .'
     MSG_COMMENT(message)
     lambda_predict=0.25_dp
   end if
!  Mimick a zero-gradient lambda, in order to avoid spurious
!  action of the inverse hessian (the next line would be a realistic estimation)
   dedv_predict=0.0_dp
!  dedv_predict=dedv_2+lambda_predict*(dedv_1-dedv_2)
!  Uses the energies, and the gradient at lambda_2
   etotal_predict=etotal_2+dedv_2*lambda_predict&
&   +(etotal_1-etotal_2-dedv_2)*lambda_predict**2

 else

!  Here, there is an acceptable solution to the 2nd degree equation
   discr=sqrt(discr)
!  The root that gives the smallest ee corresponds to  -discr
!  This is the one to be used: one aims at modelling the
!  behaviour of the function as much as possible with the
!  lowest orders of the polynomial, not the quartic term.
   ee=(sum2-discr)*0.5_dp
   dd=sum3-2.0_dp*ee
   cc=sum1-dd-ee

!  DEBUG
!  write(std_out,*)'aa,bb,cc,dd,ee',aa,bb,cc,dd,ee
!  ENDDEBUG

!  Now, must find the unique root of
!  0 = bb + 2*cc * lambda + 3*dd * lambda^2 + 4*ee * lambda^3
!  This root is unique because it was imposed that the second derivative
!  of the quartic polynomial is everywhere positive.
!  First, remove the quadratic term, by a shift of lambda
!  lambdap=lambda-lambda_shift
!  0 = bbp + ccp * lambdap + eep * lambdap^3
   eep=4.0_dp*ee
   lambda_shift=-dd/(4.0_dp*ee)
   ccp=2.0_dp*cc-12.0_dp*ee*lambda_shift**2
   bbp=bb+ccp*lambda_shift+eep*lambda_shift**3

!  DEBUG
!  write(std_out,*)'bbp,ccp,eep,lambda_shift',bbp,ccp,eep,lambda_shift
!  ENDDEBUG

!  The solution of a cubic polynomial equation is as follows :
   discr=(bbp/eep)**2+(4.0_dp/27.0_dp)*(ccp/eep)**3
!  In the present case, discr will always be positive
   discr=sqrt(discr)
   uu3=0.5_dp*(-bbp/eep+discr) ; uu=sign((abs(uu3))**(1.0_dp/3.0_dp),uu3)
   vv3=0.5_dp*(-bbp/eep-discr) ; vv=sign((abs(vv3))**(1.0_dp/3.0_dp),vv3)
   lambda_predict=uu+vv

!  Restore the shift
   lambda_predict=lambda_predict+lambda_shift
   etotal_predict=aa+bb*lambda_predict+cc*lambda_predict**2+&
&   dd*lambda_predict**3+ee*lambda_predict**4
   dedv_predict=bb+2.0_dp*cc*lambda_predict+3.0_dp*dd*lambda_predict**2+&
&   4.0_dp*ee*lambda_predict**3
   d2edv2_1=2*cc+6*dd*lambda_1+12*ee*lambda_1**2
   d2edv2_2=2*cc+6*dd*lambda_2+12*ee*lambda_2**2
   d2edv2_predict=2*cc+6*dd*lambda_predict+12*ee*lambda_predict**2

 end if

 write(message, '(a,i3)' )'   line minimization, algorithm ',4
 call wrtout(std_out,message,'COLL')
 write(message, '(a,a)' )'                        lambda      etotal ','           dedv        d2edv2    '
 call wrtout(std_out,message,'COLL')
 write(message, '(a,es12.4,es18.10,2es12.4)' )'   old point         :',lambda_2,etotal_2,dedv_2,d2edv2_2
 call wrtout(std_out,message,'COLL')
 write(message, '(a,es12.4,es18.10,2es12.4)' )'   new point         :',lambda_1,etotal_1,dedv_1,d2edv2_1
 call wrtout(std_out,message,'COLL')
 write(message, '(a,es12.4,es18.10,2es12.4)' )'   predicted point   :',lambda_predict,etotal_predict,dedv_predict,d2edv2_predict
 call wrtout(std_out,message,'COLL')
 write(message, '(a)' ) ' '
 call wrtout(std_out,message,'COLL')

end subroutine findmin
!!***

!!****f* m_numeric_tools/kramerskronig
!! NAME
!! kramerskronig
!!
!! FUNCTION
!! check or apply the Kramers Kronig relation:
!!  Re \epsilon(\omega) = 1 + \frac{2}{\pi}
!!  \int_0^\infty d\omega' frac{\omega'}{\omega'^2 - \omega^2} Im \epsilon(\omega')
!!
!! INPUTS
!!  nomega=number of real frequencies
!!  omega(nomega)= real frequencies
!!  eps(nomega)= function on the frequency grid (both real and imaginary part)
!!   real part can be used to check whether the K-K relation is satisfied or not
!!  method=method used to perform the integration
!!   0= naive integration
!!   1=simpson rule
!!  only_check= if /=0 the real part of eps is checked against the imaginary part,
!!                a final report in written but the array eps is not modified
!!              if ==0 the real part of eps is overwritten using the
!!              results obtained using the Kramers-Kronig relation
!!
!! OUTPUT
!!
!! NOTES
!! Inspired to check_kramerskronig of the DP code
!!
!! PARENTS
!!      m_paw_optics
!!
!! CHILDREN
!!
!! SOURCE

subroutine kramerskronig(nomega,omega,eps,method,only_check)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega,only_check
!arrays
 real(dp),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: eps(nomega)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ii,ip
 real(dp) :: acc,domega,eav,kkdif,kkrms,ww,wwp
 character(len=500) :: msg
!arrays
 real(dp) :: e1kk(nomega),intkk(nomega),kk(nomega)

! *************************************************************************

!Check whether the frequency grid is linear or not
 domega = (omega(nomega) - omega(1)) / (nomega-1)
 do ii=2,nomega
   if (ABS(domega-(omega(ii)-omega(ii-1))) > 0.001) then
     if (only_check/=1) then
       MSG_WARNING("Check cannot be performed since the frequency step is not constant")
       RETURN
     else
       MSG_ERROR('Cannot perform integration since frequency step is not constant')
     end if
   end if
 end do

!Check whether omega(1) is small or not
 if (omega(1) > 0.1/Ha_eV) then
   if (only_check/=1) then
     MSG_WARNING('Check cannot be performed since first frequency on the grid > 0.1 eV')
     RETURN
   else
     MSG_ERROR('Cannot perform integration since first frequency on the grid > 0.1 eV')
   end if
 end if

!If eps(nomega) is not 0 warn
 if (AIMAG(eps(nomega)) > 0.1 .and. enough<50) then
   enough=enough+1
   write(msg,'(a,f8.4,3a,f8.2,2a)')&
&   'Im epsilon for omega = ',omega(nomega)*Ha_eV,' eV',ch10,&
&   'is not yet zero, epsilon_2 = ',AIMAG(eps(nomega)),ch10,&
&   'Kramers Kronig could give wrong results'
   MSG_WARNING(msg)
   if (enough==50) then
     write(msg,'(3a)')' sufficient number of WARNINGS-',ch10,' stop writing '
     call wrtout(std_out,msg,'COLL')
   end if
 end if


!Perform Kramers-Kronig using naive integration
 select case (method)
 case (0)

   do ii=1,nomega
     ww = omega(ii)
     acc = 0.0_dp
     do ip=1,nomega
       if (ip == ii) CYCLE
       wwp = omega(ip)
       acc = acc + wwp/(wwp**2-ww**2) *AIMAG(eps(ip))
     end do
     e1kk(ii) = one + two/pi*domega* acc
   end do

!    Perform Kramers-Kronig using Simpson integration
!    Simpson O(1/N^4), from NumRec in C p 134  NumRec in Fortran p 128
 case (1)

   kk=zero

   do ii=1,nomega
     ww=omega(ii)
     do ip=1,nomega
       if (ip == ii) CYCLE
       wwp = omega(ip)
       kk(ip) = wwp/(wwp**2-ww**2) *AIMAG(eps(ip))
     end do
     call simpson_int(nomega,domega,kk,intkk)
     e1kk(ii) = one + two/pi * intkk(nomega)
   end do

 case default
   write(msg,'(a,i0)')' Wrong value for method ',method
   MSG_BUG(msg)
 end select

!at this point real part is in e1kk, need to put it into eps
 do ii=1,nomega
   eps(ii)=CMPLX(e1kk(ii),AIMAG(eps(ii)), kind=dpc)
 end do

!Verify Kramers-Kronig
 eav   = zero
 kkdif = zero
 kkrms = zero

 do ii=1,nomega
   kkdif = kkdif + ABS(REAL(eps(ii)) - e1kk(ii))
   kkrms = kkrms + (REAL(eps(ii)) - e1kk(ii))*(REAL(eps(ii)) - e1kk(ii))
   eav = eav + ABS(REAL(eps(ii)))
 end do

 eav = eav/nomega
 kkdif = (kkdif/nomega) / eav
 kkrms = (kkrms/nomega) / (eav*eav)

 kk = ABS(REAL(eps(1)) - e1kk(1)) / REAL(eps(1))

!Write data
 write(msg,'(a,f7.2,a)')' Kramers-Kronig transform is verified within ',MAXVAL(kk)*100,"%"
 call wrtout(std_out,msg,'COLL')

end subroutine kramerskronig
!!***

!!****f* ABINIT/dotproduct
!! NAME
!! dotproduct
!!
!! FUNCTION
!! scalar product of two vectors
!!
!! INPUTS
!! v1 and v2: two real(dp) vectors
!!
!! OUTPUT
!! scalar product of the two vectors
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! vector size is not checked
!!
!! NOTES
!! I've benchmarked this to be speedier than the intrinsic dot_product even on
!! big vectors. The point is that less check is performed.
!!
!! MG: FIXME: Well, optized blas1 is for sure better than what you wrote!
!! Now I dont' have time to update ref files
!!
!! PARENTS
!! cgpr,brent
!!
!! CHILDREN
!!
!!
!! SOURCE

function dotproduct(nv1,nv2,v1,v2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nv1,nv2
 real(dp) :: dotproduct
!arrays
 real(dp),intent(in) :: v1(nv1,nv2),v2(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: i,j

! *************************************************************************
 dotproduct=zero
 do j=1,nv2
  do i=1,nv1
   dotproduct=dotproduct+v1(i,j)*v2(i,j)
  end do
 end do

end function dotproduct
!!***

!!****f* m_numeric_tools/invcb
!! NAME
!! invcb
!!
!! FUNCTION
!! Compute a set of inverse cubic roots as fast as possible :
!! rspts(:)=rhoarr(:)$^\frac{-1}{3}$
!!
!! INPUTS
!!  npts=number of real space points on which density is provided
!!  rhoarr(npts)=input data
!!
!! OUTPUT
!!  rspts(npts)=inverse cubic root of rhoarr
!!
!! PARENTS
!!      m_drivexc,m_gammapositron,m_xchcth,m_xclda,m_xcpbe,m_xcpositron
!!
!! CHILDREN
!!
!! SOURCE

subroutine invcb(rhoarr,rspts,npts)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
!arrays
 real(dp),intent(in) :: rhoarr(npts)
 real(dp),intent(out) :: rspts(npts)

!Local variables-------------------------------
!scalars
 integer :: ii,ipts
 real(dp),parameter :: c2_27=2.0e0_dp/27.0e0_dp,c5_9=5.0e0_dp/9.0e0_dp
 real(dp),parameter :: c8_9=8.0e0_dp/9.0e0_dp,m1thrd=-third
 real(dp) :: del,prod,rho,rhom1,rhomtrd
 logical :: test
!character(len=500) :: message

! *************************************************************************

!Loop over points : here, brute force algorithm
!do ipts=1,npts
!rspts(ipts)=sign( (abs(rhoarr(ipts)))**m1thrd,rhoarr(ipts))
!end do
!

 rhomtrd=sign( (abs(rhoarr(1)))**m1thrd, rhoarr(1) )
 rhom1=one/rhoarr(1)
 rspts(1)=rhomtrd
 do ipts=2,npts
   rho=rhoarr(ipts)
   prod=rho*rhom1
!  If the previous point is too far ...
   if(prod < 0.01_dp .or. prod > 10._dp )then
     rhomtrd=sign( (abs(rho))**m1thrd , rho )
     rhom1=one/rho
   else
     del=prod-one
     do ii=1,5
!      Choose one of the two next lines, the last one is more accurate
!      rhomtrd=((one+third*del)/(one+two_thirds*del))*rhomtrd
       rhomtrd=((one+c5_9*del)/(one+del*(c8_9+c2_27*del)))*rhomtrd
       rhom1=rhomtrd*rhomtrd*rhomtrd
       del=rho*rhom1-one
!      write(std_out,*)rhomtrd,del
       test = del*del < 1.0e-24_dp
       if(test) exit
     end do
     if( .not. test) then
       rhomtrd=sign( (abs(rho))**m1thrd , rho )
     end if
   end if
   rspts(ipts)=rhomtrd
 end do

end subroutine invcb
!!***

!!****f* ABINIT/safe_div
!! NAME
!! safe_div
!!
!! FUNCTION
!!  Subroutine safe_div performs "safe division", that is to prevent overflow,
!!  underflow, NaN, or infinity errors.  An alternate value is returned if the
!!  division cannot be performed. (bmy, 2/26/08)
!!
!!  For more information, see the discussion on:
!!  http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/8b367f44c419fa1d/
!!
!!  Taken by HM from:
!!  http://wiki.seas.harvard.edu/geos-chem/index.php/Floating_point_math_issues#Safe_floating-point_division
!!
!! INPUTS
!!  n    : Numerator for the division
!!  d    : Divisor for the division
!!  altv : Alternate value to be returned if the division can't be done
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental subroutine safe_div(n, d, altv, q)

!Arguments ----------------------------------------------
!scalars
 real(dp),intent(in) :: n, d, altv
 real(dp),intent(out) :: q

! *********************************************************************

 if ( exponent(n) - exponent(d) >= maxexponent(n) .or. d==zero ) then
    q = altv
 else
    q = n / d
 endif

end subroutine safe_div
!!***

!!****f* ABINIT/bool2index
!! NAME
!! bool2index
!!
!! FUNCTION
!!  Allocate and return array with the indices in the input boolean array `bool_list` that evaluates to .True.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine bool2index(bool_list, out_index)

!Arguments ----------------------------------------------
!scalars
 logical,intent(in) :: bool_list(:)
 integer,allocatable,intent(inout) :: out_index(:)

!Local variables-------------------------------
 integer :: ii, cnt
! *********************************************************************

 cnt = count(bool_list)
 ABI_REMALLOC(out_index, (cnt))
 cnt = 0
 do ii=1,size(bool_list)
   if (bool_list(ii)) then
     cnt = cnt + 1
     out_index(cnt) = ii
   end if
 end do

end subroutine bool2index
!!***

END MODULE m_numeric_tools
!!***
