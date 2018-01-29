!{\src2tex{textfont=tt}}
!!****f* ABINIT/gensymshub4
!! NAME
!! gensymshub4
!!
!! FUNCTION
!! Assigns the Bravais magnetic translations for shubnikov type IV
!! symmetry groups starting from the Fedorov space group
!! and the translation, generator of the anti-ferromagnetic
!! operations (as input). It will double nsym and tnons. In the end both
!! symrel and tnons are completely determined.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! genafm(3) = translation, generator of the anti-ferromagnetic symmatry operations
!! msym = maximum number of symmetry operations
!!
!! OUTPUT
!! symafm(msym)= (anti)ferromagnetic part of symmetry operations
!!
!! SIDE EFFECTS
!! nsym=number of symmetry operations, without magnetic operations at input,
!!  and with magnetic operations at output
!! symrel(3,3,msym)=symmetry operations in real space in terms
!!  of primitive translations, without magnetic operations at input,
!!  and with magnetic operations at output
!! tnons(3,msym)=nonsymmorphic translations for symmetry operations
!!  without magnetic operations at input,
!!  and with magnetic operations at output
!!
!! NOTES
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gensymshub4(genafm,msym,nsym,symafm,symrel,tnons)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gensymshub4'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym
 integer,intent(inout) :: nsym
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(in) :: genafm(3)
 real(dp),intent(inout) :: tnons(3,msym)

!Local variables ------------------------------
!scalars
 integer :: ii
 character(len=500) :: message

! *************************************************************************

 if(msym<2*nsym)then
   write(message, '(3a)' )&
&   'The number of symmetries in the Shubnikov type IV space group',ch10,&
&   'is larger than the maximal allowed number of symmetries.'
   MSG_ERROR(message)
 end if

 do ii=1,nsym
   tnons(:,nsym+ii)=tnons(:,ii)+genafm(:)
   symrel(:,:,nsym+ii)=symrel(:,:,ii)
   symafm(ii)=1
   symafm(nsym+ii)=-1
 end do
 nsym=nsym*2

end subroutine gensymshub4

!!***
