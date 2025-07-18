!!****m* ABINIT/m_qplusg
!! NAME
!!  m_qplusg
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2025 ABINIT group (MG, FB, BG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_qplusg

 use defs_basis

 implicit none

 private
!!***

 public :: cmod_qpg
!!***

CONTAINS
!!***

!!****f* m_qplusg/cmod_qpg
!! NAME
!! cmod_qpg
!!
!! FUNCTION
!! Set up table of lengths |q+G| for the Coulomb potential.
!!
!! INPUTS
!! gvec(3,npwvec)=Reduced coordinates of the G vectors.
!! gprimd(3,3)=Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!! iq=Index specifying the q point.
!! npwvec=Number of planewaves
!! nq=Number of q points.
!! q=Coordinates of q points.
!!
!! OUTPUT
!! qplusg(npwvec)=Norm of q+G vector
!!
!! SOURCE

subroutine cmod_qpg(nq, iq, q, npwvec, gvec, gprimd, qplusg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,npwvec,nq
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 real(dp),intent(in) :: gprimd(3,3),q(3,nq)
 real(dp),intent(out) :: qplusg(npwvec)

!Local variables ------------------------------
!scalars
 integer :: ig
!arrays
 real(dp) :: gmet(3,3),gpq(3)

!************************************************************************

 ! Compute reciprocal space metrics
 gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

 if (ALL(ABS(q(:,iq)) < tol3)) then
   ! Treat q as if it were zero except when G=0
   qplusg(1)=two_pi*SQRT(DOT_PRODUCT(q(:,iq),MATMUL(gmet,q(:,iq))))
   do ig=2,npwvec
     gpq(:)=gvec(:,ig)
     qplusg(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
   end do
 else
   do ig=1,npwvec
     gpq(:)=gvec(:,ig)+q(:,iq)
     qplusg(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
   end do
 end if

end subroutine cmod_qpg
!!***

end module m_qplusg
!!***
