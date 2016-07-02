!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_dos_1band_m
!! NAME
!! get_dos_1band_m
!!
!! FUNCTION
!! calculate DOS from tetrahedron method for 1 band and 1 sppol
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dos_fractions_m=fractional DOS at each irred kpoint
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! ndosfraction_m=number of different fractional DOSs
!! tweight=sum of tetrahedron weights for each irred kpoint
!! dtweightde=energy derivative of tweight
!!
!! OUTPUT
!!  partial_dos_m(nene,ndosfraction_m)=partial DOS, for each different channel
!!  integ_dos_m_m(nene,ndosfraction_m)=integrated DOS, for each different channel
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine get_dos_1band_m (dos_fractions_m,enemin,enemax,&
&            integ_dos_m,nene,nkpt,ndosfraction_m,&
&            partial_dos_m,tweight,dtweightde)


 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dos_1band_m'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndosfraction_m,nene,nkpt
 real(dp),intent(in) :: enemax,enemin
!arrays
 real(dp),intent(in) :: dos_fractions_m(nkpt,ndosfraction_m)
 real(dp),intent(in) :: dtweightde(nkpt,nene),tweight(nkpt,nene)
 real(dp),intent(out) :: integ_dos_m(nene,ndosfraction_m)
 real(dp),intent(out) :: partial_dos_m(nene,ndosfraction_m)

!Local variables-------------------------------
!scalars
 integer :: ieps,ifract,ikpt
 real(dp) :: deltaene,eps
!arrays

! *********************************************************************

 partial_dos_m(:,:) = zero
 integ_dos_m(:,:) = zero

 deltaene = (enemax-enemin) / (nene-1)
!
!Calculate parameters of DOS at each point eps in [epsmin,epsmax]
!
 eps=enemin
 do ieps=1,nene

   do ifract=1,ndosfraction_m
     do ikpt=1,nkpt
       partial_dos_m(ieps,ifract) = partial_dos_m(ieps,ifract)+&
&       dtweightde(ikpt,ieps)*dos_fractions_m(ikpt,ifract)
       integ_dos_m(ieps,ifract) = integ_dos_m(ieps,ifract)+&
&       tweight(ikpt,ieps)*dos_fractions_m(ikpt,ifract)
     end do
   end do
   eps = eps + deltaene
 end do

end subroutine get_dos_1band_m
!!***
