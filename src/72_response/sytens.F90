!{\src2tex{textfont=tt}}
!!****f* ABINIT/sytens
!!
!! NAME
!! sytens
!!
!! FUNCTION
!! Determines the set of irreductible elements of the non-linear
!! optical susceptibility and Raman tensors
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (MVeithen)
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
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!    At the input :
!!       1   ->   element has to be computed explicitely
!!    At the output :
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!      -2   ->   element is zero by symmetry
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sytens(indsym,mpert,natom,nsym,rfpert,symrec,symrel)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sytens'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(inout) :: rfpert(3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: flag,found,i1dir,i1dir_,i1pert,i1pert_,i2dir,i2dir_,i2pert,i2pert_
 integer :: i3dir,i3dir_,i3pert,i3pert_,idisy1,idisy2,idisy3,ipesy1,ipesy2
 integer :: ipesy3,isym,pert1,pert2,pert3
!arrays
 integer :: sym1(3,3),sym2(3,3),sym3(3,3)
 integer,allocatable :: pertsy(:,:,:,:,:,:)

!***********************************************************************

!DEBUG
!write(std_out,*)'sytens : enter'
!write(std_out,*)'indsym = '
!write(std_out,*)indsym
!stop
!ENDDEBUG

 ABI_ALLOCATE(pertsy,(3,mpert,3,mpert,3,mpert))
 pertsy(:,:,:,:,:,:) = 0

!Loop over perturbations

 do i1pert_ = 1, mpert
   do i2pert_ = 1, mpert
     do i3pert_ = 1, mpert

       do i1dir_ = 1, 3
         do i2dir_ = 1, 3
           do i3dir_ = 1, 3

             i1pert = (mpert - i1pert_ + 1)
             if (i1pert <= natom) i1pert = natom + 1 - i1pert
             i2pert = (mpert - i2pert_ + 1)
             if (i2pert <= natom) i2pert = natom + 1 - i2pert
             i3pert = (mpert - i3pert_ + 1)
             if (i3pert <= natom) i3pert = natom + 1 - i3pert

             if (i1pert <= natom) then
               i1dir = i1dir_ ; i2dir = i2dir_ ; i3dir = i3dir_
             else if (i2pert <= natom) then
               i1dir = i2dir_ ; i2dir = i1dir_ ; i3dir = i3dir_
             else if (i3pert <= natom) then
               i1dir = i3dir_ ; i2dir = i2dir_ ; i3dir = i1dir_
             else
               i1dir = i1dir_ ; i2dir = i2dir_ ; i3dir = i3dir_
             end if

             if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) /= 0) then

!              Loop over all symmetries

               flag = 0
               do isym = 1, nsym

                 found = 1

!                Select the symmetric element of i1pert,i2pert,i3pert

                 if (i1pert <= natom) then
                   ipesy1 = indsym(4,isym,i1pert)
                   sym1(:,:) = symrec(:,:,isym)
                 else if (i1pert == natom + 2) then
                   ipesy1 = i1pert
                   sym1(:,:) = symrel(:,:,isym)
                 else
                   found = 0
                 end if

                 if (i2pert <= natom) then
                   ipesy2 = indsym(4,isym,i2pert)
                   sym2(:,:) = symrec(:,:,isym)
                 else if (i2pert == natom + 2) then
                   ipesy2 = i2pert
                   sym2(:,:) = symrel(:,:,isym)
                 else
                   found = 0
                 end if

                 if (i3pert <= natom) then
                   ipesy3 = indsym(4,isym,i3pert)
                   sym3(:,:) = symrec(:,:,isym)
                 else if (i3pert == natom + 2) then
                   ipesy3 = i3pert
                   sym3(:,:) = symrel(:,:,isym)
                 else
                   found = 0
                 end if

!                See if the symmetric element is available and check if some
!                of the elements may be zeor. In the latter case, they do not need
!                to be computed.


                 if ((flag /= -1).and.&
&                 (ipesy1==i1pert).and.(ipesy2==i2pert).and.(ipesy3==i3pert)) then
                   flag = sym1(i1dir,i1dir)*sym2(i2dir,i2dir)*sym3(i3dir,i3dir)
                 end if


                 do idisy1 = 1, 3
                   do idisy2 = 1, 3
                     do idisy3 = 1, 3

                       if ((sym1(i1dir,idisy1) /= 0).and.(sym2(i2dir,idisy2) /= 0).and.&
&                       (sym3(i3dir,idisy3) /= 0)) then
                         if (pertsy(idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3) == 0) then
                           found = 0
!                          exit      ! exit loop over symmetries
                         end if
                       end if


                       if ((flag == -1).and.&
&                       ((idisy1/=i1dir).or.(idisy2/=i2dir).or.(idisy3/=i3dir))) then
                         if ((sym1(i1dir,idisy1)/=0).and.(sym2(i2dir,idisy2)/=0).and.&
&                         (sym3(i3dir,idisy3)/=0)) then
                           flag = 0
                         end if
                       end if



                     end do
                   end do
                 end do


                 if (found == 1) then
                   pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = -1
                 end if


!                In case a symmetry operation only changes the sign of an
!                element, this element has to be equal to zero

                 if (flag == -1) then
                   pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = -2
                   exit
                 end if

               end do    ! close loop on symmetries



!              If the elemetn i1pert,i2pert,i3pert is not symmetric
!              to a basis element, it is a basis element

               if (pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) > -1) then
                 pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
               end if

             end if ! rfpert /= 0

           end do        ! close loop over perturbations
         end do
       end do
     end do
   end do
 end do


!Now, take into account the permutation of (i1pert,i1dir)
!and (i3pert,i3dir)

!LTEST
! if (nsym>1) then
!!LTEST
! do i1pert = 1, mpert
!   do i2pert = 1, mpert
!     do i3pert = 1, mpert

!       do i1dir = 1, 3
!         do i2dir = 1, 3
!           do i3dir = 1, 3

!             if ((i1pert /= i3pert).or.(i1dir /= i3dir)) then

!               if ((pertsy(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) == 1).and.&
!&               (pertsy(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) == 1)) then
!                 pertsy(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) = -1
!               end if

!             end if

!           end do
!         end do
!       end do

!     end do
!   end do
! end do
!!LTEST
! end if
!LTEST

 rfpert(:,:,:,:,:,:) = pertsy(:,:,:,:,:,:)

 ABI_DEALLOCATE(pertsy)


end subroutine sytens
!!***
