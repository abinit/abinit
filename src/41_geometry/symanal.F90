!{\src2tex{textfont=tt}}
!!****f* ABINIT/symanal
!!
!! NAME
!! symanal
!!
!! FUNCTION
!! Find the space group, Bravais lattice, including Shubnikov characteristics
!! from the list of symmetries (including magnetic characteristics), and lattice parameters
!! Warning : the recognition of the space group might not yet work for the
!! Shubnikov group of type IV
!! 
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! chkprim= if 1 then stop if the cell is not primitive
!! msym=default maximal number of symmetries
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,1:msym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!! tolsym=tolerance for the symmetry operations
!!
!! OUTPUT
!! bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!! genafm(3)=magnetic translation generator (in case of Shubnikov group type IV)
!! ptgroupma = magnetic point group number
!! spgroup=symmetry space group
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      ingeo,m_ab7_symmetry,m_use_ga
!!
!! CHILDREN
!!      chkprimit,getptgroupma,symbrav,symlatt,symptgroup,symspgr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symanal'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => symanal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkprim,msym,nsym
 integer,intent(out) :: ptgroupma,spgroup
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(out) :: bravais(11)
 integer,intent(in) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: tnons(3,msym)
 real(dp),intent(out) :: genafm(3)

!Local variables-------------------------------
!scalars
 integer :: iholohedry_nomagn,isym,isym_nomagn,multi
 integer :: nptsym,nsym_nomagn,shubnikov
 character(len=5) :: ptgroup,ptgroupha
 character(len=500) :: message
!arrays
 integer :: identity(3,3)
 integer,allocatable :: ptsymrel(:,:,:),symrel_nomagn(:,:,:)
 real(dp),allocatable :: tnons_nomagn(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' symanal : enter '
!call flush(6)
!stop
!ENDDEBUG

!This routine finds the Bravais characteristics, without actually
!looking at the symmetry operations.
 ABI_ALLOCATE(ptsymrel,(3,3,msym))
 call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
 ABI_DEALLOCATE(ptsymrel)

!Check whether the cell is primitive or not.
 call chkprimit(chkprim,multi,nsym,symafm,symrel)

 spgroup=0 ; ptgroupma=0 ; genafm(:)=zero

 if(multi>1)then !  Modify bravais if the cell is not primitive ; no determination of the space group
   bravais(1)=-bravais(1)
 else

!  The cell is primitive, so that the space group can be
!  determined. Need to distinguish Fedorov and Shubnikov groups.
!  Do not distinguish Shubnikov types I and II.
!  Also identify genafm, in case of Shubnikov type IV
   identity(:,:)=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
   shubnikov=1
   do isym=1,nsym
     if(symafm(isym)==-1)then
       shubnikov=3
       if(sum(abs(symrel(:,:,isym)-identity(:,:)))==0)then
         shubnikov=4
         genafm(:)=tnons(:,isym)
!        DEBUG
!        write(std_out,*)' isym=',isym
!        write(std_out,*)' symrel(:,:,isym)',symrel(:,:,isym)
!        write(std_out,*)' tnons(:,isym)',tnons(:,isym)
!        write(std_out,*)' symafm(isym)',symafm(isym)
!        ENDDEBUG
         exit
       end if
     end if
   end do

   if(shubnikov/=1)then
     if(shubnikov==3)write(message, '(a)' )' Shubnikov space group type III'
     if(shubnikov==4)write(message, '(a)' )' Shubnikov space group type IV'
     call wrtout(std_out,message,'COLL')
   end if

   if(shubnikov==1 .or. shubnikov==3)then
!    Find the correct Bravais characteristics and point group
!    Should also be used for Shubnikov groups of type IV ...
     call symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym)
!    Find the space group
     call symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)
   end if

   if(shubnikov/=1)then

!    Determine nonmagnetic symmetry operations
     nsym_nomagn=nsym/2
     ABI_ALLOCATE(symrel_nomagn,(3,3,nsym_nomagn))
     ABI_ALLOCATE(tnons_nomagn,(3,nsym_nomagn))
     isym_nomagn=0
     do isym=1,nsym
       if(symafm(isym)==1)then
         isym_nomagn=isym_nomagn+1
         symrel_nomagn(:,:,isym_nomagn)=symrel(:,:,isym)
         tnons_nomagn(:,isym_nomagn)=tnons(:,isym)
       end if
     end do

     if(shubnikov==3)then

!      DEBUG
!      write(std_out,*)' symanal : will enter symbrav with halved symmetry set'
!      write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!      do isym=1,nsym_nomagn
!      write(std_out,'(i3,2x,9i3,3es12.2,i3)')isym,symrel_nomagn(:,:,isym),tnons_nomagn(:,isym)
!      end do
!      ENDDEBUG

!      Find the point group of the halved symmetry set
       call symptgroup(iholohedry_nomagn,nsym_nomagn,ptgroupha,symrel_nomagn)

!      Deduce the magnetic point group (ptgroupma) from ptgroup and ptgroupha
       call getptgroupma(ptgroup,ptgroupha,ptgroupma)

     else if(shubnikov==4)then

!      Find the Fedorov space group of the halved symmetry set
       call symspgr(bravais,nsym_nomagn,spgroup,symrel_nomagn,tnons_nomagn,tolsym)

!      The magnetic translation generator genafm has already been determined

!      DEBUG
!      write(std_out,*)' genafm =',genafm
!      write(std_out,*)' spgroup=',spgroup
!      ENDDEBUG

     end if

     ABI_DEALLOCATE(symrel_nomagn)
     ABI_DEALLOCATE(tnons_nomagn)

   end if ! Shubnikov groups

 end if

end subroutine symanal
!!***
