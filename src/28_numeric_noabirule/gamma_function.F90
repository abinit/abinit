!{\src2tex{textfont=tt}}
!!****f* ABINIT/GAMMA_FUNCTION
!! NAME
!!  GAMMA_FUNCTION
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2016 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      pspnl_hgh_rec,pspnl_operat_rec,vso_realspace_local
!!
!! CHILDREN
!!      gsl_f90_sf_gamma
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine GAMMA_FUNCTION(X,GA)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'GAMMA_FUNCTION'
!End of the abilint section

  implicit none

#ifdef HAVE_MATH_GSL
! in case we have gsl, no need to use explicit function, just wrap the
!  call to the GSL C function in 01_gsl_ext/

  ! arguments

  real(dp),intent(in) :: x
  real(dp),intent(out) :: ga

  call gsl_f90_sf_gamma(x,ga)

#else
!       ====================================================
!       Purpose: This program computes the gamma function
!                Gamma(x) using subroutine GAMMA
!       Examples:
!                   x          Gamma(x)
!                ----------------------------
!                  1/3       2.678938534708
!                  0.5       1.772453850906
!                 -0.5      -3.544907701811
!                 -1.5       2.363271801207
!                  5.0      24.000000000000
!       ====================================================
!
!  This routine was downloaded from UIUC:
!  http://jin.ece.uiuc.edu/routines/routines.html
!
!  The programs appear to accompany a book "Computation of Special
!  Functions" (1996) John Wiley and Sons, but are distributed online
!  by the authors. Exact copyright should be checked.
!
!  Authors / copyright:
!     Shanjie Zhang and Jianming Jin
!     Proposed contact is:  j-jin1@uiuc.edu
!
!  20 October 2008:
!     Incorporated into ABINIT by M. Verstraete
!  
!
!
!       ==================================================
!       Purpose: Compute the gamma function Gamma(x)
!       Input :  x  --- Argument of Gamma(x)
!                       ( x is not equal to 0,-1,-2, etc )
!       Output:  GA --- Gamma(x)
!       ==================================================
!

  ! arguments

  real(dp),intent(in) :: x
  real(dp),intent(out) :: ga

  ! local variables
  integer :: k,m
  real(dp) :: m1,z,r,gr
  real(dp) :: G(26)

  ! source code:

  ! initialization of reference data
  G=(/1.0D0,0.5772156649015329D0, &
     &  -0.6558780715202538D0, -0.420026350340952D-1, &
     &   0.1665386113822915D0,-.421977345555443D-1, &
     &  -.96219715278770D-2, .72189432466630D-2, &
     &  -.11651675918591D-2, -.2152416741149D-3, &
     &   .1280502823882D-3, -.201348547807D-4, &
     &  -.12504934821D-5, .11330272320D-5, &
     &  -.2056338417D-6, .61160950D-8, &
     &   .50020075D-8, -.11812746D-8, &
     &   .1043427D-9, .77823D-11, &
     &  -.36968D-11, .51D-12, &
     &  -.206D-13, -.54D-14, .14D-14, .1D-15/)


  ! for the integer case, do explicit factorial
  if (X==int(X)) then
    if (X > 0.0D0) then
      GA=1.0D0
      M1=X-1
      do K=2,int(M1)
        GA=GA*K
      end do
    else
      GA=1.0D+300
    end if
  ! for the integer case, do explicit factorial
  else
    if (abs(X) > 1.0D0) then
      Z=abs(X)
      M=int(Z)
      R=1.0D0
      do K=1,M
        R=R*(Z-K)
      end do
      Z=Z-M
    else
      Z=X
    end if
    GR=G(26)
    do K=25,1,-1
      GR=GR*Z+G(K)
    end do
    GA=1.0D0/(GR*Z)
    if (abs(X) > 1.0D0) then
      GA=GA*R
      if (X < 0.0D0) GA=-PI/(X*GA*SIN(PI*X))
    end if
  end if
  return

#endif
!  end preproc for presence of GSL

end subroutine GAMMA_FUNCTION
!!***
