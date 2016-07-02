!{\src2tex{textfont=tt}}
!!****f* ABINIT/irreducible_set_pert
!!
!! NAME
!! irreducible_set_pert
!!
!! FUNCTION
!! Determines a set of perturbations that form a basis
!! in that, using symmetry, they can be used to generate
!! all other perturbations that are asked to be calculated (target).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array described above: for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of iper
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  rfdir(3)=direction for the perturbations
!!  rfpert(mpert)=information on the perturbations
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!  symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!   fourth number is 0 if the q-vector is not preserved,
!!              is 1 otherwise
!!   second index is one without time-reversal symmetry,
!!                two with time-reversal symmetry
!!
!! OUTPUT
!!   pertsy(3,mpert)= the target perturbation is described by the two last indices (idir, and ipert),
!!                    the value is 0, 1 or -1, see notes.
!!
!! NOTES
!! Output will be in the pertsy array,
!!   0 for non-target perturbations
!!   1 for basis perturbations
!!  -1 for perturbations that can be found from basis perturbations
!!
!! PARENTS
!!      get_npert_rbz,m_dvdb,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine irreducible_set_pert(indsym,mpert,natom,nsym,pertsy,rfdir,rfpert,symq,symrec,symrel)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'irreducible_set_pert'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),rfdir(3),rfpert(mpert)
 integer,intent(in) :: symq(4,2,nsym),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(out) :: pertsy(3,mpert)

!Local variables -------------------------
!scalars
 integer :: found,idir1,idisy1,ii,ipert1,ipesy1,isign,isym,itirev,jj
!arrays
 integer :: sym1(3,3)

! *********************************************************************

!Zero pertsy
 pertsy(:,:)=0

 do ipert1=1,mpert
   do idir1=1,3
     if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then

!      DEBUG
!      write(std_out,*)' for candidate idir =',idir1,&
!      &            ' ipert = ',ipert1
!      ENDDEBUG

!      Loop on all symmetries, including time-reversal
       do isym=1,nsym
         do itirev=1,2
           isign=3-2*itirev

           if(symq(4,itirev,isym)/=0)then

             found=1

!            Here select the symmetric of ipert1
             if(ipert1<=natom)then
               ipesy1=indsym(4,isym,ipert1)
               do ii=1,3
                 do jj=1,3
                   sym1(ii,jj)=symrec(ii,jj,isym)
                 end do
               end do
             else if(ipert1==(natom+2) .or. ipert1==(natom+6))then
               ipesy1=ipert1
               do ii=1,3
                 do jj=1,3
                   sym1(ii,jj)=symrel(ii,jj,isym)
                 end do
               end do
             else
               found=0
             end if

!            Now that a symmetric perturbation has been obtained,
!            including the expression of the symmetry matrix, see
!            if the symmetric perturbations are available
             if( found==1 ) then

               do idisy1=1,3
                 if(sym1(idir1,idisy1)/=0)then
                   if(pertsy(idisy1,ipesy1)==0)then
                     found=0
                     exit
                   end if
                 end if
               end do
             end if

!            Now, if still found, then it is a symmetric
!            of some linear combination of existing perturbations
             if(found==1)then

!              DEBUG
!              write(std_out,*)' all found !  isym, isign= ',isym,isign
!              write(std_out,1010)((sym1(ii,jj),ii=1,3),jj=1,3)
!              write(std_out,1010)((sym2(ii,jj),ii=1,3),jj=1,3)
!              write(std_out,*)sumr,sumi
!              1010    format(9i4)
!              ENDDEBUG

               pertsy(idir1,ipert1)=-1
!              Exit loop on symmetry operations
               exit

             end if

!            End loop on all symmetries + time-reversal
           end if
         end do
       end do

!      Now that all symmetries have been examined,
!      if still not symmetric of a linear combination
!      of basis perturbations, then it is a basis
!      perturbation
       if(pertsy(idir1,ipert1)/=-1)then
         pertsy(idir1,ipert1)=1
       end if

!      DEBUG
!      write(std_out,'(a,3i5)' ) ' irreducible_set_pert :',idir1,ipert1,pertsy(idir1,ipert1)
!      ENDDEBUG

!      End big loop on all elements
     end if
   end do
 end do

!DEBUG
!write(std_out,*)' irreducible_set_pert : exit '
!ENDDEBUG

end subroutine irreducible_set_pert
!!***
