!{\src2tex{textfont=tt}}
!!****f* ABINIT/littlegroup_pert
!!
!! NAME
!! littlegroup_pert
!!
!! FUNCTION
!! If syuse==0 and rfmeth==2, determines the set of symmetries that leaves a perturbation invariant.
!! (Actually, all symmetries that leaves a q-wavevector invariant should be used to reduce the number
!! of k-points for all perturbations. Unfortunately, one has to take into account the sign reversal of the
!! perturbation under the symmetry operations, which makes GS routines not usable for the respfn code.
!! The intermediate choice was to select only those that keep also the perturbation invariant.
!! Note that the wavevector of the perturbation must also be invariant, a translation vector in real space is NOT allowed ).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!! idir=direction of the perturbation
!! indsym(4,nsym,natom)=indirect indexing of atom labels--see subroutine symatm for definition (if nsym>1)
!! iout=if non-zero, output on unit iout
!! ipert=characteristics of the perturbation
!! natom= number of atoms
!! nsym=number of space group symmetries
!! rfmeth =
!!   1 if non-stationary block
!!   2 if stationary block
!!   3 if third order derivatives
!! symq(4,2,nsym)= Table computed by littlegroup_q.
!!   three first numbers define the G vector;
!!   fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!   second index is one without time-reversal symmetry, two with time-reversal symmetry
!! symafm(nsym)=(anti)ferromagnetic part of the symmetry operations
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! syuse= flag to use the symmetries or not. If 0 usei it, if 1 do not use it.
!! tnons(3,nsym)=nonsymmorphic translations of space group in terms
!!  of real space primitive translations (may be 0)
!! [unit]=By default the routine writes to std_out and this is very annoying if we are inside a big loop.
!!   Use unit=dev_null or a negative integer to disable writing.
!!
!! OUTPUT
!! nsym1 =number of space group symmetries that leaves the perturbation invariant
!! symaf1(nsym1)=(anti)ferromagnetic part of the corresponding symmetry operations
!! symrl1(3,3,nsym1)=corresponding 3x3 matrices of the group symmetries (real space)
!! tnons1(3,nsym1)=corresponding nonsymmorphic translations of space group in terms
!!   of real space primitive translations (may be 0)!!
!!
!! PARENTS
!!      dfpt_looppert,get_npert_rbz,m_dvdb,read_gkk
!!
!! CHILDREN
!!      stresssym,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine littlegroup_pert(gprimd,idir,indsym,iout,ipert,natom,nsym,nsym1, &
&    rfmeth,symafm,symaf1,symq,symrec,symrel,symrl1,syuse,tnons,tnons1, &
&    unit) ! Optional

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'littlegroup_pert'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => littlegroup_pert
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: idir,iout,ipert,natom,nsym,rfmeth,syuse
 integer,intent(in),optional :: unit
 integer,intent(out) :: nsym1
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym),symq(4,2,nsym)
 integer,intent(in) :: symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(out) :: symaf1(nsym),symrl1(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),tnons(3,nsym)
 real(dp),intent(out) :: tnons1(3,nsym)

!Local variables -------------------------
!scalars
 integer :: idir1,ii,istr,isym,jj,nsym_test,tok,ount
 character(len=500) :: msg
!arrays
 integer :: sym_test(3,3,2)
 real(dp) :: str_test(6)

! *********************************************************************

 ount = std_out; if (present(unit)) ount = unit

 nsym1=0
 if((ipert==natom+3 .or. ipert==natom+4) .and. syuse==0 .and. rfmeth==2) then
!  Strain perturbation section
!  Use ground state routine which symmetrizes cartesian stress as a quick
!  and dirty test for the invariance of the strain (ipert,idir) under
!  each candidate symmetry
!  I am presently assuming that translations are acceptable because I dont
!  see why not.

   istr=3*(ipert-natom-3)+idir
   nsym_test=2
!  Store identity as first element for test
   sym_test(:,:,1)=0
   sym_test(1,1,1)=1; sym_test(2,2,1)=1; sym_test(3,3,1)=1
   do isym=1,nsym
     sym_test(:,:,2)=symrec(:,:,isym)
     str_test(:)=0.0_dp
     str_test(istr)=1.0_dp
     call stresssym(gprimd,nsym_test,str_test,sym_test)
     if(abs(str_test(istr)-1.0_dp)<tol8)then
!      The test has been successful !
       nsym1=nsym1+1
       symaf1(nsym1)=symafm(isym)
       do ii=1,3
         tnons1(ii,nsym1)=tnons(ii,isym)
         do jj=1,3
           symrl1(ii,jj,nsym1)=symrel(ii,jj,isym)
         end do
       end do
     end if
   end do

 else if(ipert>natom .or. syuse/=0 .or. rfmeth/=2)then

!  Not yet coded for d/dk or electric field perturbations
   nsym1=1
   do ii=1,3
     tnons1(ii,1)=0._dp
     symaf1(1)=1
     do jj=1,3
       symrl1(ii,jj,1)=0
       if(ii==jj)symrl1(ii,jj,1)=1
     end do
   end do

 else

   do isym=1,nsym
!    Check that the symmetry operation preserves the wavevector
!    (a translation is NOT allowed)
     if(symq(4,1,isym)==1 .and.&
&     symq(1,1,isym)==0 .and.&
&     symq(2,1,isym)==0 .and.&
&     symq(3,1,isym)==0          )then
!      Check that the symmetry operation preserves the atom
       if(ipert==indsym(4,isym,ipert))then
!        Check if the direction is preserved
         tok=1
         do idir1=1,3
           if((idir1==idir.and.symrec(idir,idir1,isym)/=1) .or.&
&           (idir1/=idir.and.symrec(idir,idir1,isym)/=0))then
             tok=0
           end if
         end do
         if(tok==1)then
!          All the tests have been successful !
           nsym1=nsym1+1
           symaf1(nsym1)=symafm(isym)
           do ii=1,3
             tnons1(ii,nsym1)=tnons(ii,isym)
             do jj=1,3
               symrl1(ii,jj,nsym1)=symrel(ii,jj,isym)
             end do
           end do
         end if

       end if
     end if
   end do
 end if

 if (nsym1<1) then
   write(msg,'(a,i0,a)')&
&   ' The number of selected symmetries should be > 0, while it is nsym=',nsym1,'.'
   MSG_BUG(msg)
 end if

 if (nsym1 /= 1) then
   if (iout /= ount .and. iout > 0) then
     write(msg,'(a,i5,a)')' Found ',nsym1,' symmetries that leave the perturbation invariant.'
     call wrtout(iout,msg,'COLL')
   end if
   write(msg,'(a,i5,a)')' littlegroup_pert : found ',nsym1,' symmetries that leave the perturbation invariant :'
   call wrtout(ount,msg,'COLL')
 else
   if (iout /= ount .and. iout > 0) then
     write(msg,'(a,a)')' The set of symmetries contains',' only one element for this perturbation.'
     call wrtout(iout,msg,'COLL')
   end if
   write(msg,'(a)')' littlegroup_pert : only one element in the set of symmetries for this perturbation :'
   call wrtout(ount,msg,'COLL')
 end if

 if (ount > 0) then
   do isym=1,nsym1
     write(msg, '(9i4)' )((symrl1(ii,jj,isym),ii=1,3),jj=1,3)
     call wrtout(ount,msg,'COLL')
   end do
 end if

end subroutine littlegroup_pert
!!***
