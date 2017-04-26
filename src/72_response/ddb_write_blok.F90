!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddb_write_blok
!!
!! NAME
!! ddb_write_blok
!!
!! FUNCTION
!! This routine writes blocks of data in the DDBs.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! choice= (2 => write), (3 => write minimal info )
!! mpert =maximum number of ipert
!! msize=maximum size of the arrays flags and values
!! nunit=unit number for the data block file
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! ddb = ddb block datastructure
!! ddb%typ=type of the block:
!!   0 => total energy
!!   1 => second-order energy derivatives, non-stationary block
!!   2 => second-order energy derivatives, stationary block
!!   3 => third-order energy derivatives
!!   4 => first-order energy derivatives: forces, stresses and polarization
!!   5 => second-order eigenvalue derivatives
!! ddb%flg(msize)=flag for every matrix element (0=> the element is
!!  not in the data block), (1=> the element is in the data blok)
!! ddb%qpt(9)=wavevector of the perturbation(s). The elements from
!!  1 to 3 are used if we are dealing with the 2nd derivative of
!!  total energy (only one wavevector), while all elements are
!!  used in case of a third order derivative of total energy
!!  (three wavevector could be present)
!! ddb%nrm(3)=normalization factors for the three allowed wavevectors.
!! ddb%val(2,msize)=real(dp), complex, value of the
!!  matrix elements that are present in the data block
!! blkval2(2,msize,mband,nkpt) = value of the matrix elements
!!  that are present in a block of EIGR2D/EIGI2D
!!
!! NOTES
!! only executed by one processor.
!!
!! PARENTS
!!      dfptnl_doutput,gstate,mblktyp1,mblktyp5
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ddb_write_blok(ddb,iblok,choice,mband,mpert,msize,nkpt,nunit,&
&     blkval2,kpt) !optional

 use defs_basis
 use m_profiling_abi
 use m_ddb

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_write_blok'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,mband,mpert,msize,nkpt,nunit
 integer,intent(in) :: iblok
 type(ddb_type),intent(in) :: ddb
!arrays
 real(dp),intent(in),optional :: kpt(3,nkpt)
 real(dp),intent(in),optional :: blkval2(2,msize,mband,nkpt)

!Local variables -------------------------
!scalars
 integer :: iband,idir1,idir2,idir3,ii,ikpt,ipert1,ipert2,ipert3
 integer :: nelmts

! *********************************************************************
 

!Count the number of elements
 nelmts=0
 do ii=1,msize
   if(ddb%flg(ii,iblok)==1)nelmts=nelmts+1
 end do

!Write the block type and number of elements
 write(nunit,*)' '
 if (ddb%typ(iblok) == 0) then
   write(nunit, '(a,i8)' )&
&   ' Total energy                 - # elements :',nelmts
 else if (ddb%typ(iblok)==1) then
   write(nunit, '(a,i8)' )&
&   ' 2nd derivatives (non-stat.)  - # elements :',nelmts
 else if(ddb%typ(iblok)==2) then
   write(nunit, '(a,i8)' )&
&   ' 2nd derivatives (stationary) - # elements :',nelmts
 else if(ddb%typ(iblok)==3) then
   write(nunit, '(a,i8)' )&
&   ' 3rd derivatives              - # elements :',nelmts
 else if (ddb%typ(iblok) == 4) then
   write(nunit, '(a,i8)' )&
&   ' 1st derivatives              - # elements :',nelmts
 else if (ddb%typ(iblok) == 5) then
   write(nunit, '(a,i8)' )&
&   ' 2nd eigenvalue derivatives   - # elements :',nelmts
 end if

!Write the 2nd derivative block
 if(ddb%typ(iblok)==1.or.ddb%typ(iblok)==2)then

!  Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

!  Write the matrix elements
   if(choice==2)then
     ii=0
     do ipert2=1,mpert
       do idir2=1,3
         do ipert1=1,mpert
           do idir1=1,3
             ii=ii+1
             if(ddb%flg(ii,iblok)==1)then
               write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,&
&               ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
             end if
           end do
         end do
       end do
     end do
   end if

!  Write the 3rd derivative block
 else if(ddb%typ(iblok)==3)then

!  Write the phonon wavevectors
   write(nunit, '(a,3es16.8,f6.1)' )&
&   ' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
   write(nunit, '(a,3es16.8,f6.1)' )&
&   '    ',(ddb%qpt(ii,iblok),ii=4,6),ddb%nrm(2,iblok)
   write(nunit, '(a,3es16.8,f6.1)' )&
&   '    ',(ddb%qpt(ii,iblok),ii=7,9),ddb%nrm(3,iblok)

!  Write the matrix elements
   if(choice==2)then
     ii=0
     do ipert3=1,mpert
       do idir3=1,3
         do ipert2=1,mpert
           do idir2=1,3
             do ipert1=1,mpert
               do idir1=1,3
                 ii=ii+1
                 if(ddb%flg(ii,iblok)==1)then
                   write(nunit, '(6i4,2d22.14)' )&
&                   idir1,ipert1,idir2,ipert2,idir3,ipert3,&
&                   ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
                 end if
               end do
             end do
           end do
         end do
       end do
     end do
   end if

!  Write total energy
 else if (ddb%typ(iblok) == 0) then
   if (choice == 2) then
     write(nunit,'(2d22.14)')ddb%val(1,1,iblok),ddb%val(2,1,iblok)
   end if

!  Write the 1st derivative blok
 else if (ddb%typ(iblok) == 4) then
   if (choice == 2) then
     ii = 0
     do ipert1 = 1, mpert
       do idir1 = 1, 3
         ii = ii + 1
         if (ddb%flg(ii,iblok) == 1) then
           write(nunit,'(2i4,2d22.14)')idir1,ipert1,&
&           ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
         end if
       end do
     end do
   end if

 else if (ddb%typ(iblok)==5) then
!  Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
!  Write the matrix elements
   if(choice==2)then
     if(present(blkval2).and.present(kpt))then
       do ikpt=1,nkpt
         write(nunit,'(a,3es16.8)')' K-point:',(kpt(ii,ikpt),ii=1,3)
         do iband=1,mband
           write(nunit,'(a,i3)')' Band:',iband
           ii=0
           do ipert2=1,mpert
             do idir2=1,3
               do ipert1=1,mpert
                 do idir1=1,3
                   ii=ii+1
                   if(ddb%flg(ii,iblok)==1)then
                     write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,&
&                     blkval2(1,ii,iband,ikpt),blkval2(2,ii,iband,ikpt)
                   end if
                 end do !idir1
               end do  !ipert1
             end do   !idir2
           end do    !ipert2
         end do     !iband
       end do      !ikpt
     end if !blkval2
   end if !choice
 end if !ddb%typ(iblok)

end subroutine ddb_write_blok
!!***
