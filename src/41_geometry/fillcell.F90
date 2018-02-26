!{\src2tex{textfont=tt}}
!!****f* ABINIT/fillcell
!! NAME
!! fillcell
!!
!! FUNCTION
!! Computes the atomic position of all the atoms in the unit cell starting
!! with the symmetry operations and the atoms from the asymetric unit cell.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natrd = number of atoms in the assymetric unit cell
!!  natom = total number of atoms (to be checked)
!!  nsym = number of symmetry operations
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  tolsym=tolerance on symmetries
!!  typat(1:natrd)=type integer for each atom in cell
!!  xred(3,1:natrd)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At input, for the assymetric unit cell
!!  nucdipmom(3,1:natrd)=nuclear magnetic dipole moments of the atoms
!!  spinat(3,1:natrd)=spin-magnetization of the atoms
!!  typat(1:natrd)=type integer for each atom in cell
!!  xred(3,1:natrd)=reduced dimensionless atomic coordinates
!!
!!  At output, for the complete unit cell
!!  nucdipmom(3,1:natom)=nuclear magnetic dipole moments of the atoms
!!  spinat(3,1:natom)=spin-magnetization of the atoms
!!  typat(1:natom)=type integer for each atom in cell
!!  xred(3,1:natom)=reduced dimensionless atomic coordinates
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


subroutine fillcell(natom,natrd,nsym,nucdipmom,spinat,symafm,symrel,tnons,tolsym,typat,xred)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fillcell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,natrd,nsym
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,intent(inout) :: typat(natom)
 real(dp),intent(in) :: tolsym
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(inout) :: nucdipmom(3,natom),spinat(3,natom),xred(3,natom)

!Local variables ------------------------------
!scalars
 integer :: curat,flagch,flageq,ii,iij,jj,kk
 character(len=500) :: message
!arrays
 integer :: bcktypat(nsym*natrd)
 real(dp) :: bckat(3),bcknucdipmom(3,nsym*natrd)
 real(dp) :: bckspinat(3,nsym*natrd),bckxred(3,nsym*natrd)

! *************************************************************************

!DEBUG
!write(std_out,*)' fillcell : enter with nsym, natrd= ',nsym,natrd
!write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!do ii=1,nsym
!write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!end do
!write(std_out,*)' Describe the input atoms (index,typat,xred,spinat)'
!do jj=1,natrd
!write(std_out,'(i3,2x,i3,6es12.2)')jj,typat(jj),xred(:,jj),spinat(:,jj)
!end do
!ENDDEBUG

 curat=0

!Cycle over all the symmetry operations
 do ii=1,nsym

!  Cycle over all the atoms in the assymetric unit cell
   do jj=1,natrd

!    Symmetry operation application
     bckat(:)=matmul(symrel(:,:,ii),xred(:,jj))+tnons(:,ii)

!    Normalization of the coordinates in [0,1)
     do iij=1,3
       do while (bckat(iij)<-tolsym)
         bckat(iij)=bckat(iij)+1.0d0
       end do
       do while (bckat(iij)>=1.0d0-tolsym)
         bckat(iij)=bckat(iij)-1.0d0
       end do
     end do

!    Check for duplicate atoms
     flagch=0
     do kk=1,curat
       flageq=0
       if ( abs(bckxred(1,kk)-bckat(1))<tolsym  .and. &
&       abs(bckxred(2,kk)-bckat(2))<tolsym  .and. &
&       abs(bckxred(3,kk)-bckat(3))<tolsym       ) exit
       flagch=flagch+1
     end do

     if (flagch==curat) then
!      Add the obtained atom to the bckxred list
       curat=curat+1
       bckxred(:,curat)=bckat
       bcktypat(curat)=typat(jj)
       bcknucdipmom(:,curat)=nucdipmom(:,jj)
       bckspinat(:,curat)=spinat(:,jj)*symafm(ii)
     end if

   end do

 end do

!DEBUG
!write(std_out,*)' fillcell : Proposed coordinates ='
!do ii=1,curat
!write(std_out,'(i4,3es16.6)' )ii,bckxred(:,ii)
!end do
!ENDDEBUG

 if (curat>natom) then
   write(message, '(a,i3,a,a,i7,a,a,a,a)' )&
&   'The number of atoms obtained from symmetries, ',curat,ch10,&
&   'is greater than the input number of atoms, natom=',natom,ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify natom or the symmetry data in the input file.'
   MSG_ERROR(message)
 end if

 if (curat<natom) then
   write(message, '(a,i3,a,a,i7,a,a,a,a)' )&
&   'fillcell : The number of atoms obtained from symmetries, ',curat,ch10,&
&   'is lower than the input number of atoms, natom=',natom,ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify natom or the symmetry data in the input file.'
   MSG_ERROR(message)
 end if

!Assignment of symmetry to xred
 xred(:,1:natom)=bckxred(:,1:natom)
 typat(1:natom)=bcktypat(1:natom)
 nucdipmom(1:3,1:natom)=bcknucdipmom(1:3,1:natom)
 spinat(1:3,1:natom)=bckspinat(1:3,1:natom)

!DEBUG
!write(std_out,*)' fillcell : exit with natom=',natom
!write(std_out,*)' Describe the output atoms (index,typat,xred,spinat)'
!do jj=1,natom
!write(std_out,'(i3,2x,i3,6es12.2)')jj,typat(jj),xred(:,jj),spinat(:,jj)
!end do
!ENDDEBUG

end subroutine fillcell
!!***
