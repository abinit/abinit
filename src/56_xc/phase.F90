!{\src2tex{textfont=tt}}
!!****f* ABINIT/phase
!! NAME
!! phase
!!
!! FUNCTION
!! Compute ph(ig)=$\exp(\pi\ i \ n/ngfft)$ for n=0,...,ngfft/2,-ngfft/2+1,...,-1
!! while ig runs from 1 to ngfft.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ngfft=number of points
!! OUTPUT
!!  ph(2*ngfft)=phase array (complex)
!!
!! NOTES
!! XG 990504 : changed the formulation, in order to preserve
!! the invariance between n and -n, that was broken for n=ngfft/2 if ngfft even.
!! Simply suppresses the corresponding sine.
!!
!! PARENTS
!!      xcden,xcpot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine phase(ngfft,ph)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phase'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft
!arrays
 real(dp),intent(out) :: ph(2*ngfft)

!Local variables-------------------------------
!scalars
 integer :: id,ig,nn
 real(dp) :: arg,fac

! *************************************************************************

 id=ngfft/2+2
 fac=pi/dble(ngfft)
 do ig=1,ngfft
   nn=ig-1-(ig/id)*ngfft
   arg=fac*dble(nn)
   ph(2*ig-1)=cos(arg)
   ph(2*ig)  =sin(arg)

 end do
!XG 990504 Here zero the corresponding sine
 if((ngfft/2)*2==ngfft) ph(2*(id-1))=zero

end subroutine phase
!!***
