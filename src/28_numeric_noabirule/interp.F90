#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function interp(n,z,f,z0,zz)

 use defs_basis
 use m_numeric_tools, only : linfit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'interp'
!End of the abilint section

implicit none
integer :: n
complex(gwpc) :: z(n), f(n), z0, zz
complex(gwpc) :: interp
complex(gwpc) :: a, b
real(dp) :: x(n),smrt
if(.false.)write(std_out,*)z0
x(:) = real(z(:))
smrt=linfit(n,x,f,a,b)
if( smrt>0.1/Ha_eV) then
  write(std_out,*) '**WARNING: the values are not linear'
  write(std_out,*) smrt,a,b,f
endif
interp = a * zz + b
return
end function interp


function dinterp(n,z,f,z0,zz)

 use defs_basis
 use m_numeric_tools, only : linfit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dinterp'
!End of the abilint section

implicit none
integer n
complex(gwpc) :: z(n), f(n), z0, zz
complex(gwpc) :: dinterp
complex(gwpc) a, b
real(dp) :: x(n), smrt
if(.false.)write(std_out,*)z0,zz
x(:) = real(z(:))
smrt=linfit(n,x,f,a,b)
if( smrt>0.1/Ha_eV) then
  write(std_out,*) '**WARNING: the values are not linear'
  write(std_out,*) smrt,a,b,f
endif
dinterp = a
return
end function dinterp


function taylor_interp(n,z,f,z0,zz)

use defs_basis
! calculate in zz the Taylor polinomial around z0 that interpolate the function f
! at the n points z

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'taylor_interp'
 use interfaces_28_numeric_noabirule, except_this_one => taylor_interp
!End of the abilint section

implicit none
integer n
complex(gwpc) :: z(n), f(n), z0, zz
complex(gwpc) :: taylor_interp
integer i
complex(gwpc) :: c(n), zz_minus_z0_to_i

call calculate_taylor_c(n,z,f,z0,c)

! calculate the Taylor polinomial in zz
zz_minus_z0_to_i = 1
taylor_interp = 0
do i = 1, n
  taylor_interp = taylor_interp + c(i) * zz_minus_z0_to_i
  zz_minus_z0_to_i = zz_minus_z0_to_i * (zz - z0)
enddo
return
end function taylor_interp


function dtaylor_interp(n,z,f,z0,zz)
use defs_basis
! calculate in zz the Taylor polinomial derivative

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtaylor_interp'
 use interfaces_28_numeric_noabirule, except_this_one => dtaylor_interp
!End of the abilint section

implicit none
integer n
complex(gwpc) :: z(n), f(n), z0, zz
complex(gwpc) :: dtaylor_interp
integer i
complex(gwpc) :: c(n), zz_minus_z0_to_i

call calculate_taylor_c(n,z,f,z0,c)

! calculate the Taylor polinomial in zz
zz_minus_z0_to_i = 1
dtaylor_interp = 0
do i = 2, n
  dtaylor_interp = dtaylor_interp + (i-1) * c(i) * zz_minus_z0_to_i
  zz_minus_z0_to_i = zz_minus_z0_to_i * (zz - z0)
enddo
return
end function dtaylor_interp


subroutine calculate_taylor_c(n,z,f,z0,c)

 use m_profiling_abi
 use defs_basis
 use m_errors

! calculate the Taylor coefficients polinomial expansion around z0
! for the function f at the n points z

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calculate_taylor_c'
!End of the abilint section

implicit none
integer n
complex(gwpc) :: z(n), f(n), z0, c(n)
integer, parameter :: nrhs=1
integer i,j,info
complex(gwpc) :: a(n,n), bb(n,nrhs)

! calculate the coefficient matrix
a(:,1) = 1
do i = 1, n
  do j = 2, n
    a(i,j) = a(i,j-1) * (z(i) - z0)
  enddo
enddo

! solve the linear system to find the Taylor expansion coefficients
bb(:,1) = f(:)
MSG_ERROR('introduce cgesv')
!call cgesv(n,nrhs,a,n,ipiv,bb,n,info)
if(info/=0) then
  write(std_out,*) info
  MSG_ERROR('cgesv failed')
endif
c(:) = bb(:,1)
return
end subroutine calculate_taylor_c
