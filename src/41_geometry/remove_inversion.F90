!{\src2tex{textfont=tt}}
!!****f* ABINIT/remove_inversion
!! NAME
!! remove_inversion
!!
!! FUNCTION
!!  Remove the inversion symmetry from a symmetry set as well
!!  all the improper rotations (if present)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsym=initial number of symmetries
!!  symrel(3,3,nsym)=Initial set of symmetry operarations in real space
!!  tnons(3,nsym)=Initial fractional translations
!!
!! OUTPUT
!!  nsym_out=Number of symmetries in the set without improper rotation
!!  symrel_out(:,:) [pointer] = output symmetries without improper rotations
!!  tnons_out(:) [pointer] = fractional translations associated to symrel_out
!!  pinv=-1 if the inversion has been removed, 1 otherwise
!!
!! NOTES
!!  Note the use of pointers, memory is allocated inside the procedure and passed back
!!  to the caller. Thus memory deallocation is relegated to the caller. To be on the safe side
!!  the pointers should be nullified before entering.
!!
!! PARENTS
!!      m_crystal,m_io_kss,outkss
!!
!! CHILDREN
!!      set2unit,symdet,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine remove_inversion(nsym,symrel,tnons,nsym_out,symrel_out,tnons_out,pinv)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_numeric_tools, only : isinteger, set2unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'remove_inversion'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => remove_inversion
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: nsym_out,pinv
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,pointer :: symrel_out(:,:,:)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),pointer :: tnons_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: is,is2,is_discarded,is_inv,is_retained,nsym2
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: determinant(nsym),inversion(3,3),symrel2(3,3,nsym)
 real(dp) :: dtnons(3),tnons2(3,nsym)

! *********************************************************************

 msg = 'Removing inversion related symmetrie from initial set'
 MSG_WARNING(msg)
!
!=== Find the occurence of the inversion symmetry ===
 call set2unit(inversion) ; inversion=-inversion

 is_inv=0 ; found=.FALSE.
 do while (is_inv<nsym .and. .not.found)
   is_inv=is_inv+1 ; found=ALL(symrel(:,:,is_inv)==inversion)
 end do
 if (found) then
   write(msg,'(a,i3)')' The inversion is symmetry operation no. ',is_inv
 else
   write(msg,'(a)')' The inversion was not found in the symmetries list.'
 end if
 call wrtout(std_out,msg,'COLL')
!
!=== Find the symmetries that are related through the inversion symmetry ===
 call symdet(determinant,nsym,symrel)
 nsym2=0
 do is=1,nsym-1
   do is2=is+1,nsym

     dtnons(:)=tnons(:,is2)-tnons(:,is)-tnons(:,is_inv)
     found=ALL(symrel(:,:,is)==-symrel(:,:,is2)).and.isinteger(dtnons,tol8)

     if (found) then
       nsym2=nsym2+1
!      * Retain symmetries with positive determinant
       if (ALL(tnons(:,is2)<tol8).and.ALL(tnons(:,is)<tol8)) then
         is_retained=is2 ; is_discarded=is
         if (determinant(is)==1) then
           is_retained=is  ; is_discarded=is2
         end if
       else if (ALL(tnons(:,is2)<tol8)) then
         is_retained=is2 ; is_discarded=is
       else
         is_retained=is ;  is_discarded=is2
       end if

       symrel2(:,:,nsym2)=symrel(:,:,is_retained)
       tnons2   (:,nsym2)=tnons   (:,is_retained)
       write(msg,'(a,i3,a,i3,3a,i3,a)')&
&       ' Symmetry operations no. ',is,' and no. ',is2,&
&       ' are related through the inversion.',ch10,&
&       ' Symmetry operation no. ',is_discarded,' will be suppressed.'
       call wrtout(std_out,msg,'COLL')
     end if ! found

   end do !is2
 end do !is

 if (nsym2/=(nsym/2).or.nsym==1) then
   write(msg,'(a)')' Program uses the original set of symmetries '
   call wrtout(std_out,msg,'COLL')
   nsym_out=nsym
   ABI_ALLOCATE(symrel_out,(3,3,nsym))
   ABI_ALLOCATE(tnons_out,(3,nsym))
   symrel_out(:,:,:)=symrel(:,:,1:nsym)
   tnons_out(:,:)=tnons(:,1:nsym)
   pinv=1
 else
   write(msg,'(a)')' Inversion related operations have been suppressed from symmetries list.'
   call wrtout(std_out,msg,'COLL')
   nsym_out=nsym2
   ABI_ALLOCATE(symrel_out,(3,3,nsym2))
   ABI_ALLOCATE(tnons_out,(3,nsym2))
   symrel_out(:,:,:)=symrel2(:,:,1:nsym2)
   tnons_out(:,:)=tnons(:,1:nsym2)
   pinv=-1
 end if

end subroutine remove_inversion
!!***
