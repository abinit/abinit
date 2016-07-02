!{\src2tex{textfont=tt}}
!!****f* ABINIT/matpointsym
!! NAME
!! matpointsym
!!
!! FUNCTION
!! For given order of point group, symmetrizes a 3x3 input matrix using the
!! point symmetry of the input atom
!!
!! COPYRIGHT
!! Copyright (C) 2007-2016 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! iatom=index of atom to symmetrize around
!! natom=number of atoms in cell
!! nsym=order of group
!! rprimd(3,3)= real space primitive vectors
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! tnons(3,nsym) = nonsymmorphic translations
!! xred(3,natom)=locations of atoms in reduced coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! mat3(3,3) = matrix to be symmetrized, in cartesian frame
!!
!! PARENTS
!!      make_efg_el,make_efg_ion,make_efg_onsite
!!
!! CHILDREN
!!      dgemm,dgemv,mati3inv,matrginv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine matpointsym(iatom,mat3,natom,nsym,rprimd,symrel,tnons,xred)

 use defs_basis
 use m_profiling_abi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matpointsym'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom,natom,nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym),xred(3,natom)
 real(dp),intent(inout) :: mat3(3,3)

!Local variables-------------------------------
!scalars
 integer :: cell_index,cell_indexp,ii,isym,nsym_point
 real(dp) :: xreddiff
!arrays
 integer :: symrel_it(3,3)
 real(dp) :: mat3_tri(3,3),mat3_tri_sym(3,3),rprimd_inv(3,3),tmp_mat(3,3)
 real(dp) :: xredp(3)

!**************************************************************************

!DEBUG
!write(std_out,*)' matpointsym : enter '
!enddo
!ENDDEBUG

!copy rprimd input and construct inverse
 rprimd_inv = rprimd
 call matrginv(rprimd_inv,3,3)

!transform input mat3 to triclinic frame with rprimd^{-1} * mat3 * rprimd
 call dgemm('N','N',3,3,3,one,rprimd_inv,3,mat3,3,zero,tmp_mat,3)
 call dgemm('N','N',3,3,3,one,tmp_mat,3,rprimd,3,zero,mat3_tri,3)

!loop over symmetry elements to obtain symmetrized input matrix
 mat3_tri_sym = zero
 nsym_point = 0
 do isym = 1, nsym

! skip any nonsymmorphic symmetry elements, want to consider point elements only
   if(dot_product(tnons(:,isym),tnons(:,isym))>tol8) cycle

! for current symmetry element, find transformed reduced coordinates of target atom
! via xredp = symrel * xred
   call dgemv('N',3,3,one,dble(symrel(:,:,isym)),3,xred(:,iatom),1,zero,xredp,1)


! shift xredp into the same unit cell as xred, for comparison
! label cells as 0..1:0 1..2:1 2..3:2 and -1..0:-1 -2..-1:-2 and so forth
   do ii = 1, 3

     cell_index = int(xred(ii,iatom))
     if(xred(ii,iatom) < zero) cell_index = cell_index - 1
     cell_indexp = int(xredp(ii))
     if(xredp(ii) < zero) cell_indexp = cell_indexp - 1
     
     do while (cell_indexp < cell_index)
       xredp(ii) = xredp(ii)+one
       cell_indexp = cell_indexp + 1
     end do
     do while (cell_indexp > cell_index)
       xredp(ii) = xredp(ii)-one
       cell_indexp = cell_indexp - 1
     end do

   end do

! now compare xredp to xred
   xreddiff = dot_product(xredp-xred(:,iatom),xredp-xred(:,iatom))

   if (xreddiff < tol8) then

!  accumulate symrel^{-1}*mat3_tri*symrel into mat3_tri_sym iff xredp = xred + L,
!  where is a lattice vector, so symrel leaves the target atom invariant

!  mati3inv gives the inverse transpose of symrel
     call mati3inv(symrel(:,:,isym),symrel_it)
     call dgemm('N','N',3,3,3,one,mat3_tri,3,dble(symrel(:,:,isym)),3,zero,tmp_mat,3)
     call dgemm('T','N',3,3,3,one,dble(symrel_it),3,tmp_mat,3,one,mat3_tri_sym,3)
     nsym_point = nsym_point + 1
   end if

 end do

!normalize by number of point symmetry operations
 mat3_tri_sym = mat3_tri_sym/dble(nsym_point)

!transform mat3_tri_sym to cartesian frame with rprimd * mat3_tri_sym * rprimd^{-1}

 call dgemm('N','N',3,3,3,one,mat3_tri_sym,3,rprimd_inv,3,zero,tmp_mat,3)
 call dgemm('N','N',3,3,3,one,rprimd,3,tmp_mat,3,zero,mat3,3)

!DEBUG
!write(std_out,*)' matpointsym : exit '
!ENDDEBUG

end subroutine matpointsym
!!***
