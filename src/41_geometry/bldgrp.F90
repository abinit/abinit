!{\src2tex{textfont=tt}}
!!****f* ABINIT/bldgrp
!! NAME bldgrp
!! bldgrp
!!
!! FUNCTION
!! Yields all the symmetry operations starting from the generators.
!! Applies all the generators onto themselves, and obtains all the other operations.
!! Iterates until it reaches nsym.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! msym = default number of symmetry operations
!! nsym = number of symmetry operations
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym) = 3D matrix containg symmetry operations
!! tnons(3,msym) = 2D matrix containing translations of the symmery operations
!!
!! OUTPUT
!!
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym) = 3D matrix containg symmetry operations
!! tnons(3,msym) = 2D matrix containing translations of the symmery operations
!!
!! SIDE EFFECTS
!! nogen = number of generators, number of operations to be applied onto
!!  themselves
!!
!! NOTES
!!
!!
!! PARENTS
!!      symsgcube,symsghexa,symsgortho,symsgtetra
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine bldgrp(msym,nogen,nsym,symafm,symrel,tnons)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bldgrp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym
 integer,intent(inout) :: nogen
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(inout) :: tnons(3,msym)

!Local variables ------------------------------
!matrintoper(3,3) & matrinttransl(3) are intermediate arrays of the new
!      symmetry operations obtained, in order to check their uniqueness.
!flagop,flagtr = flags used during the checking of the similarity between
!      the obtained operation and the already existent ones
!ii,ijk,ijkl,jjj,kk = counters in the cycles
!scalars
 integer :: flagop,flagtr,ii,ijk,ijkl,jjj,kk,matrintsymafm,nogen_new
 real(dp) :: nastyzero
 character(len=500) :: message
!arrays
 integer :: bcksymafm(2*msym),bcksymrel(3,3,2*msym),matrintoper(3,3)
 real(dp) :: bcktnons(3,2*msym),matrinttransl(3)

! *************************************************************************

 nastyzero=0.1

!DEBUG
! write(std_out,*)' bldgrp : enter, builds the space group symmetry '
! write(std_out,*)' bldgrp : number of generators : ',nogen
! write(std_out,*)' bldgrp : nsym,msym=',nsym,msym
!ENDDEBUG

 if (nogen<1) then
   write(message, '(a,i4,a,a,a,a,a)' )&
&   'The number of generators nogen is ',nogen,&
&   'and it should be greater than one',ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

!Transfer the generators to bcksymrel
 do ii=1,nogen
   bcksymrel(:,:,ii)=symrel(:,:,ii)
   bcktnons(:,ii)=tnons(:,ii)
   bcksymafm(ii)=symafm(ii)
 end do

!DEBUG
 write(std_out,*)' Describe the different generators (index,symrel,tnons,symafm)'
 do ii=1,nogen
   write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
 end do
!ENDDEBUG

!Simply iterate until the group is complete
 do ijkl=1,nsym

!  DEBUG
!  write(std_out,*)' bldgrp : in loop, ijkl,nogen=',ijkl,nogen
!  ENDDEBUG

   nogen_new=nogen

   do jjj=2,nogen
     do kk=2,nogen

!      Computing block of the new symmetry operation according to:
!      !   $ { R1 | v1 }{ R2 | v2 } = { R1.R2 | v1+R1.v2 } $
       matrintoper(:,:) = matmul(bcksymrel(:,:,jjj),bcksymrel(:,:,kk))
       matrinttransl(:) = bcktnons(:,jjj)+matmul(bcksymrel(:,:,jjj),bcktnons(:,kk))
       matrintsymafm    = bcksymafm(jjj)*bcksymafm(kk)

!      Rescaling translation between 0 and 1
       do ii=1,3
         if (matrinttransl(ii)>=0.9) then
           do while (matrinttransl(ii)>=0.9)
             matrinttransl(ii)=matrinttransl(ii)-1.0
           end do
         end if
         if (matrinttransl(ii)<0.0) then
           do while (matrinttransl(ii)<0.0)
             matrinttransl(ii)=matrinttransl(ii)+1.0
           end do
         end if
         if ( abs(matrinttransl(ii))<nastyzero) matrinttransl(ii)=0.0
         if ( abs(matrinttransl(ii)-1.0)<nastyzero) matrinttransl(ii)=0.0
       end do

!      Cheking block to validate the new symmetry operation
       do ijk=1,nogen_new

         flagop=0 ; flagtr=0

!        Check for rotation similarity
         if(sum((matrintoper-bcksymrel(:,:,ijk))**2)==0)flagop=1

!        Check for translation similarity
         if(maxval((matrinttransl-bcktnons(:,ijk))**2)<nastyzero**2)flagtr=1

         if(flagop+flagtr==2)exit

       end do

!      Add the new determined symmetry if it is unique
       if (flagtr+flagop<2) then
         nogen_new=nogen_new+1
!        DEBUG
!         write(std_out,*)' added one more symmetry : nogen_new=',nogen_new
!        ENDDEBUG
         bcksymrel(:,:,nogen_new)=matrintoper(:,:)
         bcktnons(:,nogen_new)=matrinttransl(:)
         bcksymafm(nogen_new)=matrintsymafm
       end if

     end do
   end do

   nogen=nogen_new

   if(nogen==nsym)exit

 end do

!Transfer of the calculated symmetry to the routine output
 if (nogen==nsym) then
   symrel(:,:,1:nsym)=bcksymrel(:,:,1:nsym)
   tnons(:,1:nsym)=bcktnons(:,1:nsym)
   symafm(1:nsym)=bcksymafm(1:nsym)
 else
!  Problem with the generation of the symmetry operations
   write(message, '(a,i7,a,a,i7)' )&
&   'The symmetries obtained are  ',nogen,ch10,&
&   'and they should be ',nsym
   MSG_BUG(message)
 end if

!DEBUG
!write(std_out,*)' bldgrp : exit with  ',nogen,' operation symmetries'
!ENDDEBUG

end subroutine bldgrp
!!***
