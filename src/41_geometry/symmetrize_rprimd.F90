!{\src2tex{textfont=tt}}
!!****f* ABINIT/symmetrize_rprimd
!! NAME
!! symmetrize_rprimd
!!
!! FUNCTION
!! Supposing the input rprimd does not preserve the length and angles
!! following the symmetries, will generates a new set rprimd,
!! on the basis of the expected characteristics of the conventional cell,
!! as specified in bravais(:)
!!
!! COPYRIGHT
!! Copyright (C) 2015-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! nsym=actual number of symmetries
!! symrel(3,3,1:nsym)=symmetry operations in real space in terms of primitive translations
!! tolsym=tolerance for the symmetry operations (only for checking purposes, the new set rprimd will
!!     be coherent with the symmetry operations at a much accurate level).
!!
!! SIDE EFFECTS
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      chkorthsy,holocell,matr3inv,metric
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symmetrize_rprimd(bravais,nsym,rprimd,symrel,tolsym)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_geometry,     only : metric

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symmetrize_rprimd'
 use interfaces_32_util
 use interfaces_41_geometry, except_this_one => symmetrize_rprimd
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: bravais(11),symrel(3,3,nsym)
 real(dp),intent(inout) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: foundc,iexit,ii,jj
 real(dp):: ucvol,reldiff
 character(len=500) :: msg
!arrays
 real(dp):: aa(3,3),ait(3,3),cell_base(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd_new(3,3)

! *************************************************************************

!Build the conventional cell basis vectors in cartesian coordinates
 aa(:,1)=bravais(3:5)
 aa(:,2)=bravais(6:8)
 aa(:,3)=bravais(9:11)
!Inverse transpose
 call matr3inv(aa,ait)
 do ii=1,3
   cell_base(:,ii)=ait(ii,1)*rprimd(:,1)+ait(ii,2)*rprimd(:,2)+ait(ii,3)*rprimd(:,3)
 end do

!Enforce the proper holohedry on the conventional cell vectors.
 call holocell(cell_base,1,foundc,bravais(1),tolsym)

!Reconstruct the dimensional primitive vectors
 do ii=1,3
   rprimd_new(:,ii)=aa(1,ii)*cell_base(:,1)+aa(2,ii)*cell_base(:,2)+aa(3,ii)*cell_base(:,3)
 end do

!Check whether the modification make sense
 do ii=1,3
   do jj=1,3
     reldiff=(rprimd_new(ii,jj)-rprimd(ii,jj))/sqrt(sum(rprimd(:,jj)**2))
!    Allow for twice tolsym
     if(abs(reldiff)>two*tolsym)then
       write(msg,'(a,6(2a,3es14.6))')&
&       'Failed rectification of lattice vectors to comply with Bravais lattice identification, modifs are too large',ch10,&
&       '  rprimd    =',rprimd(:,1),ch10,&
&       '             ',rprimd(:,2),ch10,&
&       '             ',rprimd(:,3),ch10,&
&       '  rprimd_new=',rprimd_new(:,1),ch10,&
&       '             ',rprimd_new(:,2),ch10,&
&       '             ',rprimd_new(:,3)
       MSG_ERROR_CLASS(msg, "TolSymError")
     end if
   end do
 end do

 rprimd(:,:)=rprimd_new(:,:)

!Check whether the symmetry operations are consistent with the lattice vectors
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 iexit=1

 call chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)

end subroutine symmetrize_rprimd
!!***

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
