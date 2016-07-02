!{\src2tex{textfont=tt}}
!!****f* ABINIT/phim
!! NAME
!! phim
!!
!! FUNCTION
!! Computes Phi_m[theta]=Sqrt[2] cos[m theta],      if m>0
!!                       Sqrt[2] sin[Abs(m) theta], if m<0
!!                       1                        , if m=0
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  costeta= cos(theta)  (theta= input angle)
!!  mm = index m
!!  sinteta= sin(theta)  (theta= input angle)
!!
!! OUTPUT
!!  phim= Phi_m(theta) (see above)
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!     setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function phim(costheta,sintheta,mm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phim'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mm
 real(dp) :: phim
 real(dp),intent(in) :: costheta,sintheta

!Local variables-------------------------------

! *********************************************************************

 if (mm==0)  phim=one
 if (mm==1)  phim=sqrt2*costheta
 if (mm==-1) phim=sqrt2*sintheta
 if (mm==2)  phim=sqrt2*(costheta*costheta-sintheta*sintheta)
 if (mm==-2) phim=sqrt2*two*sintheta*costheta
 if (mm==3)  phim=sqrt2*&
& (costheta*(costheta*costheta-sintheta*sintheta)&
& -sintheta*two*sintheta*costheta)
 if (mm==-3) phim=sqrt2*&
& (sintheta*(costheta*costheta-sintheta*sintheta)&
& +costheta*two*sintheta*costheta)

 end function phim

!!***
