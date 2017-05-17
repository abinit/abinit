!{\src2tex{textfont=tt}}
!!****f* ABINIT/sincos
!! NAME
!! sincos
!!
!! FUNCTION
!! Update the sine and cosine values, needed inside the
!! pseudopotential routines psp1lo and psp1nl.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iq  = number of current wavevector q
!!  irmax = number of values  of r on the radial grid to be computed
!!  mmax = dimension of pspwk and rad
!!  pspwk(:,1,1) and pspwk(:,2,1) : sine and cosine of 2$\pi$ dq * rad
!!  pspwk(:,1,2) and pspwk(:,2,2) : sine and cosine of 2$\pi$ previous q * rad
!!  rad(mmax) radial grid
!!  tpiq = 2 $\pi$ * current wavevector q
!!
!! OUTPUT
!!  pspwk(*,1,2) and pspwk(*,2,2) : sine and cosine of 2$\pi$ current q * rad
!!
!! NOTES
!! The speed was a special concern, so iterative computation
!! based on addition formula is possible. Interestingly,
!! this algorithm places strong constraints on accuracy,
!! so this routine is machine-dependent.
!!
!! PARENTS
!!      psp1lo,psp1nl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sincos(iq,irmax,mmax,pspwk,rad,tpiq)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sincos'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,irmax,mmax
 real(dp),intent(in) :: tpiq
!arrays
 real(dp),intent(in) :: rad(mmax)
 real(dp),intent(inout) :: pspwk(mmax,2,2)

!Local variables-------------------------------
!scalars
 integer :: ir,nstep
 real(dp) :: prevcos,prevsin
 logical :: testmipspro
#if defined HAVE_LINALG_MLIB
 real(dp) :: halfpi
#endif


! *************************************************************************

#if defined HAVE_LINALG_MLIB
 halfpi=asin(1.0d0)
#endif

 if(iq==2)then

!  Here set up the sin and cos at iq=2
   do ir=2,irmax
     pspwk(ir,1,2)=pspwk(ir,1,1)
     pspwk(ir,2,2)=pspwk(ir,2,1)
   end do

 else
!  
!  The sensitivity of the algorithm to changes of nstep
!  has been tested : for all the machines except SGI - R10000 ,
!  either using only the hard way, or
!  using up to nstep=40 causes changes at the level
!  of 1.0d-16 in the total energy. Larger values of
!  nstep might be possible, but the associated residual
!  is already very small ! The accelerated computation of
!  sine and cosine is essential for a good speed on IBM, but,
!  fortunately, on the SGI - R10000 the normal computation is fast enough.

   testmipspro=.false.
#ifdef FC_MIPSPRO
   testmipspro=.true.
#endif
   nstep=40
   if(iq-(iq/nstep)*nstep == 0 .or. testmipspro)then

!    Every nstep steps, uses the hard way
     do ir=2,irmax
#if defined HAVE_LINALG_MLIB
!      There is a bug in the hp library !! Sine is slightly inaccurate !
       pspwk(ir,1,2)=cos(tpiq*rad(ir)-halfpi)
#else
       pspwk(ir,1,2)=sin(tpiq*rad(ir))
#endif
       pspwk(ir,2,2)=cos(tpiq*rad(ir))
     end do

   else

!    Here the fastest way, iteratively
     do ir=2,irmax
       prevsin=pspwk(ir,1,2)
       prevcos=pspwk(ir,2,2)
       pspwk(ir,1,2)=prevsin*pspwk(ir,2,1)+prevcos*pspwk(ir,1,1)
       pspwk(ir,2,2)=prevcos*pspwk(ir,2,1)-prevsin*pspwk(ir,1,1)
     end do

   end if 

 end if ! iq==2

end subroutine sincos
!!***
