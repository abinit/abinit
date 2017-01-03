!{\src2tex{textfont=tt}}
!!****f* ABINIT/littlegroup_q
!! NAME
!! littlegroup_q
!!
!! FUNCTION
!! Determines the symmetry operations by which reciprocal vector q is preserved,
!! modulo a primitive reciprocal lattice vector, and the time-reversal symmetry.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=number of space group symmetries
!! qpt(3)= vector in reciprocal space
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! [prtvol]=integer flag defining the verbosity of output. =0 if no output is provided.
!! prtgkk= integer flag. If 1 provide output of electron-phonon "gkk" matrix elements, for further
!!     treatment by mrggkk utility or anaddb utility. If 0 no output is provided. 
!!
!! OUTPUT
!! symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!     fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!     second index is one without time-reversal symmetry, two with time-reversal symmetry
!! timrev=1 if the time-reversal symmetry preserves the wavevector, modulo a reciprocal lattice vector (in principle, see below).
!!
!! NOTES
!! The condition is :
!!    $q =  O  S(q) - G$
!! with O being either the identity or the time reversal symmetry (= inversion in reciprocal space)
!! and G being a primitive vector of the reciprocal lattice.
!! If the time-reversal (alone) also preserves q, modulo a lattice vector, then timrev is set to 1, otherwise 0.
!!
!! TODO
!! timrev is put to 1 only for Gamma.
!! Better handling should be provided in further version.
!!
!! PARENTS
!!      get_npert_rbz,m_bz_mesh,m_ddb,m_dvdb,m_dynmat,m_esymm,m_gkk,m_phgamma
!!      m_sigmaph,memory_eval,read_gkk,respfn
!!
!! CHILDREN
!!      wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine littlegroup_q(nsym,qpt,symq,symrec,symafm,timrev,prtvol,use_sym)

 use defs_basis

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'littlegroup_q'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(in),optional :: prtvol,use_sym
 integer,intent(out) :: timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(in) :: symafm(nsym)
 integer,intent(out) :: symq(4,2,nsym)
 real(dp),intent(in) :: qpt(3)

!Local variables -------------------------
!scalars
 integer :: ii,isign,isym,itirev,my_prtvol
 real(dp),parameter :: tol=2.d-8
 real(dp) :: reduce
 character(len=500) :: message
!arrays
 real(dp) :: difq(3),qsym(3),shift(3)

! *********************************************************************

 my_prtvol=0 ; if (PRESENT(prtvol)) my_prtvol=prtvol

! Initialise the array symq
 symq(:,:,:)=0

 isym = symafm(1) ! just to fool abirules and use symafm for the moment

 do isym=1,nsym
!   if (symafm(isym) /= 1) cycle ! skip afm symops
! TODO: check how much of the afm syms are coded in the rf part of the code. cf
! test v3 / 12
   do itirev=1,2
     isign=3-2*itirev  ! isign is 1 without time-reversal, -1 with time-reversal

!    Get the symmetric of the vector
     do ii=1,3
       qsym(ii)=qpt(1)*isign*symrec(ii,1,isym)&
&       +qpt(2)*isign*symrec(ii,2,isym)&
&       +qpt(3)*isign*symrec(ii,3,isym)
     end do

!    Get the difference between the symmetric and the original vector

     symq(4,itirev,isym)=1
     do ii=1,3
       difq(ii)=qsym(ii)-qpt(ii)
!      Project modulo 1 in the interval ]-1/2,1/2] such that difq = reduce + shift
       call wrap2_pmhalf(difq(ii),reduce,shift(ii))
       if(abs(reduce)>tol)symq(4,itirev,isym)=0
     end do

!    SP: When prtgkk is asked (GKK matrix element will be output), one has to
!    disable symmetries. There is otherwise a jauge problem with the unperturbed
!    and the perturbed wavefunctions. This leads to a +- 5% increase in computational 
!    cost but provide the correct GKKs (i.e. the same as without the use of
!    symmerties.)  

     if (PRESENT(use_sym)) then
       if (use_sym == 0) then   
         symq(4,itirev,isym)=0
         symq(4,itirev,1)=1 
       end if
     end if
     
!    If the operation succeded, change shift from real(dp) to integer, then exit loop
     if(symq(4,itirev,isym)/=0)then
       if (my_prtvol>0) then
         if(itirev==1)write(message,'(a,i4,a)')' littlegroup_q : found symmetry',isym,' preserves q '
         if(itirev==2)write(message,'(a,i4,a)')' littlegroup_q : found symmetry ',isym,' + TimeReversal preserves q '
         call wrtout(std_out,message,'COLL')
       end if
!      Uses the mathematical function NINT = nearest integer
       do ii=1,3
         symq(ii,itirev,isym)=nint(shift(ii))
       end do
     end if

   end do !itirev
 end do !isym

!Test time-reversal symmetry
 timrev=1
 do ii=1,3
!  Unfortunately, this version does not work yet ...
!  call wrap2_pmhalf(2*qpt(ii),reduce,shift(ii))
!  if(abs(reduce)>tol)timrev=0
!  So, this is left ...
   if(abs(qpt(ii))>tol)timrev=0
 end do

 if(timrev==1.and.my_prtvol>0)then
   write(message, '(a,a,a)' )&
&   ' littlegroup_q : able to use time-reversal symmetry. ',ch10,&
&   '  (except for gamma, not yet able to use time-reversal symmetry)'
   call wrtout(std_out,message,'COLL')
 end if

end subroutine littlegroup_q
!!***
