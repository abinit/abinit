!{\src2tex{textfont=tt}}
!!****f* ABINIT/symredcart
!! NAME
!! symredcart
!!
!! FUNCTION
!! Convert a symmetry operation from reduced coordinates (integers)
!! to cartesian coordinates (reals)
!! Can operate in real or reciprocal space
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! symred(3,3)=symmetry matrice in reduced coordinates (integers) (real or reciprocal space)
!! aprim(3,3)=real or reciprocal space dimensional primitive translations (see below)
!! bprim(3,3)=real or reciprocal space dimensional primitive translations (see below)
!!
!! OUTPUT
!! symcart(3,3)=symmetry matrice in cartesian coordinates (reals)
!!
!! NOTES
!! When aprim=rprimd and bprim=gprimd, the routine operates in real space (on a real space symmetry)
!! When aprim=gprimd and bprim=rprimd, the routine operates in reciprocal space (on a real space symmetry)
!!
!! PARENTS
!!      m_matlu,m_phonons,symrhg
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symredcart(aprim,bprim,symcart,symred)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symredcart'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: symred(3,3)
 real(dp),intent(in) :: aprim(3,3),bprim(3,3)
 real(dp),intent(out) :: symcart(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: symtmp
!arrays
 real(dp) :: work(3,3)

! *************************************************************************

 work=zero
 do kk=1,3
   do jj=1,3
     symtmp=dble(symred(jj,kk))
     do ii=1,3
       work(ii,jj)=work(ii,jj)+bprim(ii,kk)*symtmp
     end do
   end do
 end do

 symcart=zero
 do kk=1,3
   do jj=1,3
     symtmp=work(jj,kk)
     do ii=1,3
       symcart(ii,jj)=symcart(ii,jj)+aprim(ii,kk)*symtmp
     end do
   end do
 end do

end subroutine symredcart
!!***
