!{\src2tex{textfont=tt}}
!!****f* ABINIT/symrelrot
!! NAME
!! symrelrot
!!
!! FUNCTION
!! Transform the symmetry matrices symrel expressed in the coordinate system rprimd,
!! to symmetry matrices symrel expressed in the new coordinate system rprimd_new
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! rprimd_new(3,3)=new dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! symrel(3,3,nsym)=symmetry operations in real space in terms
!! of primitive translations rprimd at input and rprimd_new at output
!!
!! PARENTS
!!      ingeo,m_esymm,symbrav,symlatt,symspgr
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrelrot'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(inout) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),rprimd_new(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj
 real(dp) :: val
 character(len=500) :: message
!arrays
 real(dp) :: coord(3,3),coordinvt(3,3),matr1(3,3),matr2(3,3),rprimd_invt(3,3)

!**************************************************************************

!Compute the coordinates of rprimd_new in the system defined by rprimd(:,:)
 call matr3inv(rprimd,rprimd_invt)
 do ii=1,3
   coord(:,ii)=rprimd_new(1,ii)*rprimd_invt(1,:)+ &
&   rprimd_new(2,ii)*rprimd_invt(2,:)+ &
&   rprimd_new(3,ii)*rprimd_invt(3,:)
 end do

!Transform symmetry matrices in the system defined by rprimd_new
 call matr3inv(coord,coordinvt)
 do isym=1,nsym
   do ii=1,3
     matr1(:,ii)=symrel(:,1,isym)*coord(1,ii)+&
&     symrel(:,2,isym)*coord(2,ii)+&
&     symrel(:,3,isym)*coord(3,ii)
   end do
   do ii=1,3
     matr2(:,ii)=coordinvt(1,:)*matr1(1,ii)+&
&     coordinvt(2,:)*matr1(2,ii)+&
&     coordinvt(3,:)*matr1(3,ii)
   end do

!  Check that the new symmetry matrices are made of integers, and store them
   do ii=1,3
     do jj=1,3
       val=matr2(ii,jj)
!      Need to allow for twice tolsym, in case of centered Bravais lattices (but do it for all lattices ...)
       if(abs(val-nint(val))>two*tolsym)then
         write(message,'(2a,a,i3,a,a,3es14.6,a,a,3es14.6,a,a,3es14.6)')&
&         'One of the components of symrel is non-integer,',ch10,&
&         '  for isym=',isym,ch10,&
&         '  symrel=',matr2(:,1),ch10,&
&         '         ',matr2(:,2),ch10,&
&         '         ',matr2(:,3)
         MSG_ERROR_CLASS(message, "TolSymError")
       end if
       symrel(ii,jj,isym)=nint(val)
     end do
   end do
!  End loop isym
 end do

end subroutine symrelrot
!!***
