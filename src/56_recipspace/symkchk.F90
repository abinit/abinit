!{\src2tex{textfont=tt}}
!!****f* ABINIT/symkchk
!! NAME
!! symkchk
!!
!! FUNCTION
!! Checks that the set of k points chosen for a response function
!! calculation has the full space group symmetry, modulo time reversal if appropriate.
!! Returns ierr/=0 with error message if not satisfied
!! Currently used only when strain perturbation is treated. Based on symkpt.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (DRH, XG, LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! kptns(3,nkpt)= k vectors in reciprocal space
!! nkpt = number of k-points whose weights are wtk
!! nsym=number of space group symmetries
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! timrev: if 1, the time reversal operation has to be taken into account
!! if 0, no time reversal symmetry.
!!
!! OUTPUT
!!  msg=Error message if ierr /= 0
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


integer function symkchk(kptns,nkpt,nsym,symrec,timrev,errmsg) result(ierr)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_fstrings,    only : itoa, sjoin

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symkchk'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym,timrev
 character(len=*),intent(out) :: errmsg
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: kptns(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: identi,ii,ikpt,ikpt2,imatch,isym,jj,tident
 real(dp) :: difk,reduce
 character(len=500) :: message
!arrays
 real(dp) :: ksym(3)

! *********************************************************************
 ierr = 0

 if(timrev/=1 .and. timrev/=0)then
   write(errmsg, '(3a,i0,a)' )&
&   'timrev should be 0 or 1, while',ch10,&
&   'it is equal to ',timrev,'.'
   ierr = 1; return
 end if

 if(nsym/=1)then
!  Find the identity symmetry operation
   do isym=1,nsym
     tident=1
     do jj=1,3
       if(symrec(jj,jj,isym)/=1)tident=0
       do ii=1,3
         if( ii/=jj .and.&
&         symrec(ii,jj,isym)/=0)tident=0
       end do
     end do
     if(tident==1)then
       identi=isym
       call wrtout(std_out,sjoin(' symkchk: found identity with number:', itoa(identi)))
       exit
     end if
   end do
   if(tident==0)then
     errmsg = 'Did not found the identity operation.'
     ierr = 1; return
   end if
 end if

!Here begins the serious business
!The length sorting, etc. of symkpt have been dropped because the
!computational cost is estimated to be negligible.

 if(nsym>1 .or. timrev==1)then

!  Outer loop over kpts
   do ikpt=1,nkpt-1

!    Loop on the symmetries
!    For each k-point and each symmetry transformation, a matching
!    k-pointpt must be found, modulo time reversal if appropriate
     do isym=1,nsym

!      Get the symmetric of the vector
       do ii=1,3
         ksym(ii)= kptns(1,ikpt)*symrec(ii,1,isym)&
&         +kptns(2,ikpt)*symrec(ii,2,isym)&
&         +kptns(3,ikpt)*symrec(ii,3,isym)
       end do

!      Second loop k-points
       do ikpt2=1,nkpt

!        Test for match of symmetric and any vector (including original)
         imatch=1
         do ii=1,3
           difk= ksym(ii)-kptns(ii,ikpt2)
           reduce=difk-anint(difk)
           if(abs(reduce)>tol8)imatch=0
         end do
         if(imatch==1)exit

!        Test for match with time reversal
         if(timrev==1)then
           imatch=1
           do ii=1,3
             difk= ksym(ii)+kptns(ii,ikpt2)
             reduce=difk-anint(difk)
             if(abs(reduce)>tol8)imatch=0
           end do
           if(imatch==1)exit
         end if

       end do ! End secondary loop over k-points
       if (imatch/=1) then
         write(errmsg, '(a,a,a,i4,a,i4,a,a,a,a)' )&
&         'k-point set must have full space-group symmetry',ch10,&
&         'there is no match for kpt',ikpt,' transformed by symmetry',isym,ch10,&
&         'Action: change kptopt to 2 or 3 and/or change or use shiftk',ch10,&
&         'shiftk = 0 0 0 is always a safe choice.'
         ierr = 2; return
       end if

     end do ! End loop on isym
   end do ! End primary loop over k-points

   write(message,'(a)')' symkchk : k-point set has full space-group symmetry.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

end function symkchk
!!***
