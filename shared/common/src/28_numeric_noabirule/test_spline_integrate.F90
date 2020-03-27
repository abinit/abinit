#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program test_spline_integrate

use defs_basis
use m_splines
implicit none

real(dp), parameter :: int_tol = 5.0e-3_dp
real(dp), parameter :: f1(6) = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
real(dp), parameter :: f2(6) = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0]
real(dp), parameter :: f3(6) = [0.0, 1.0, 8.0, 27.0, 64.0, 125.0]

logical :: f1_ok, f2_ok, f3_ok, test_ok
integer :: npts
real(dp) :: dx, cal_int, ref_int, int_err

#if defined DEBUG_MODE
write(std_out,'(a,a)') "Unitary test: spline_integrate", ch10
write(std_out,'(1x,4(1x,a12),1x,a5)') "Function", "Analytical", "Calculated", &
& "Error (rel)", "Check"
write(std_out,'(1x,4(1x,a12),1x,a5)') "------------", "------------", &
& "------------", "------------", "-----"
#endif

! Integrate f1
ref_int = 12.5_dp
npts = size(f1)
dx = f1(2) - f1(1)
call spline_integrate(cal_int, npts, dx, f1)
int_err = abs((cal_int - ref_int) / ref_int)
f1_ok = ( int_err < int_tol )
#if defined DEBUG_MODE
write(std_out,'(2x,a12,1x,e12.5,1x,e12.5,1x,e12.5,1x,5l)') "f(x) = x  ", &
& ref_int, cal_int, int_err, f1_ok
#endif

! Integrate f2
ref_int = 125.0_dp / 3.0_dp
npts = size(f2)
dx = f2(2) - f2(1)
call spline_integrate(cal_int, npts, dx, f2)
int_err = abs((cal_int - ref_int) / ref_int)
f2_ok = ( int_err < int_tol )
#if defined DEBUG_MODE
write(std_out,'(2x,a12,1x,e12.5,1x,e12.5,1x,e12.5,1x,5l)') "f(x) = x^2", &
& ref_int, cal_int, int_err, f2_ok
#endif

! Integrate f3
ref_int = 625.0_dp / 4.0_dp
npts = size(f3)
dx = f3(2) - f3(1)
call spline_integrate(cal_int, npts, dx, f3)
int_err = abs((cal_int - ref_int) / ref_int)
f3_ok = ( int_err < int_tol )
#if defined DEBUG_MODE
write(std_out,'(2x,a12,1x,e12.5,1x,e12.5,1x,e12.5,1x,5l)') "f(x) = x^3", &
& ref_int, cal_int, int_err, f3_ok
#endif

! Report test result
test_ok = ( f1_ok .and. f2_ok .and. f3_ok )
#if defined DEBUG_MODE
write(std_out,*)
if ( test_ok ) then
  write(std_out,'(a,a)') "TEST OK", ch10
else
  write(std_out,'(a,a)') "TEST FAILED", ch10
end if
#else
if ( .not. test_ok ) then
  write(std_out,'(a,a)') "TEST FAILED", ch10
end if
#endif

end program test_spline_integrate
