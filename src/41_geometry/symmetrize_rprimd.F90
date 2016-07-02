!{\src2tex{textfont=tt}}
!!****f* ABINIT/symmetrize_rprimd
!! NAME 
!! symmetrize_rprimd
!!
!! FUNCTION
!! Supposing the input rprimd does not preserve the length and angles following the symmetries,
!! will generates a new set rprimd, on the basis of the expected characteristics of the conventional cell,
!! as specified in bravais(:) 
!!
!! COPYRIGHT
!! Copyright (C) 2015-2016 ABINIT group (XG)
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
 real(dp) :: reldiff,tolsym
 character(len=500):: msg
!arrays
 integer,intent(in) :: bravais(11),symrel(3,3,nsym)
 real(dp),intent(inout) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: foundc,iexit,ii,jj
 real(dp):: ucvol
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
