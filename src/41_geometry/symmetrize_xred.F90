!{\src2tex{textfont=tt}}
!!****f* ABINIT/symmetrize_xred
!! NAME
!! symmetrize_xred
!!
!! FUNCTION
!! Symmetrize atomic coordinates using input symmetry matrices symrel
!! which are expressed in terms of the basis of real space primitive
!! translations (array elements are integers).
!! Input array indsym(4,isym,iatom) gives label of atom into which iatom
!! is rotated by INVERSE of symmetry element isym and also gives primitive
!! translation to get back to unit cell.
!! This version uses improvement in algorithm suggested by Andrew
!! Horsfield (see symatm.f).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! indsym(4,nsym,natom)=indirect indexing array giving label of atom
!!   into which iatom is rotated by symmetry element isym
!! natom=number of atoms
!! nsym=number of symmetries in group
!! symrel(3,3,nsym)=symmetry matrices in terms of real space
!!   primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for symmetries
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! xred(3,natom)=
!!  (input) atomic coordinates in terms of real space translations
!!  (output) symmetrized atomic coordinates in terms
!!    of real space translations
!!
!! PARENTS
!!      ingeo,mover,nonlinear,respfn,scfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symmetrize_xred(indsym,natom,nsym,symrel,tnons,xred)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symmetrize_xred'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer  :: iatom,ib,isym
 integer  :: ii,jj
 real(dp) :: fc1,fc2,fc3
 real(dp) :: diff
 logical  :: dissimilar
!arrays
 real(dp) :: tsum(3),tt(3)
 real(dp),allocatable :: xredsym(:,:)
 real(dp) :: transl(3) ! translation vector

! *************************************************************************
!
!Check whether group contains more than identity;
!if not then simply return
 if (nsym>1) then

!  loop over atoms
   ABI_ALLOCATE(xredsym,(3,natom))
   do iatom=1,natom
     tsum(1)=0.0d0
     tsum(2)=0.0d0
     tsum(3)=0.0d0
!    
!    loop over symmetries
     do isym=1,nsym
!      atom ib is atom into which iatom is rotated by inverse of
!      symmetry isym (inverse of symrel(mu,nu,isym))
       ib=indsym(4,isym,iatom)
!      Find the reduced coordinates after translation=t(indsym)+transl
       fc1=xred(1,ib)+dble(indsym(1,isym,iatom))
       fc2=xred(2,ib)+dble(indsym(2,isym,iatom))
       fc3=xred(3,ib)+dble(indsym(3,isym,iatom))
!      Compute [S * (x(indsym)+transl) ] + tnonsymmorphic
       tt(:)=dble(symrel(:,1,isym))*fc1+&
&       dble(symrel(:,2,isym))*fc2+&
&       dble(symrel(:,3,isym))*fc3+ tnons(:,isym)

!      Average over nominally equivalent atomic positions
       tsum(:)=tsum(:)+tt(:)
     end do
!    
!    Set symmetrized result to sum over number of terms
     xredsym(:,iatom)=tsum(:)/dble(nsym)

!    End loop over iatom
   end do

   transl(:)=xredsym(:,1)-nint(xredsym(:,1))

!  Compute the smallest translation to an integer
   do jj=2,natom
     do ii=1,3
       diff=xredsym(ii,jj)-nint(xredsym(ii,jj))
       if (diff<transl(ii)) transl(ii)=diff
     end do
   end do

!  Test if the translation on each direction is small
!  Tolerance 1E-13
   do ii=1,3
     if (abs(transl(ii))>1e-13) transl(ii)=0.0
   end do

!  Execute translation
   do jj=1,natom
     do ii=1,3
       xredsym(ii,jj)=xredsym(ii,jj)-transl(ii)
     end do
   end do

!  Test if xredsym is too similar to xred
!  Tolerance 1E-15
   dissimilar=.FALSE.
   do jj=1,natom
     do ii=1,3
       if (abs(xredsym(ii,jj)-xred(ii,jj))>1E-15) dissimilar=.TRUE.
     end do
   end do

   if (dissimilar) xred(:,:)=xredsym(:,:)
   ABI_DEALLOCATE(xredsym)

!  End condition of nsym/=1
 end if

end subroutine symmetrize_xred
!!***
