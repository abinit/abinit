!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcart2deloc
!! NAME
!! xcart2deloc
!!
!! FUNCTION
!!  Calculate values of delocalized coordinates as a function of
!!  cartesian ones. First primitive internals, then B matrix,
!!  then F, then U then delocalized internals.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! deloc <type(delocint)>=Important variables for
!!   |                           pred_delocint
!!   |
!!   | nang     = Number of angles
!!   | nbond    = Number of bonds
!!   | ncart    = Number of cartesian directions
!!   |             (used for constraints)
!!   | ndihed   = Number of dihedrals
!!   | nrshift  = Dimension of rshift
!!   | ninternal= Number of internal coordinates
!!   |            ninternal=nbond+nang+ndihed+ncart
!!   |
!!   | angs(2,3,nang)  = Indexes to characterize angles
!!   | bonds(2,2,nbond)= For a bond between iatom and jatom
!!   |                   bonds(1,1,nbond) = iatom
!!   |                   bonds(2,1,nbond) = icenter
!!   |                   bonds(1,2,nbond) = jatom
!!   |                   bonds(2,2,nbond) = irshift
!!   | carts(2,ncart)  = Index of total primitive internal,
!!   |                   and atom (carts(2,:))
!!   | dihedrals(2,4,ndihed)= Indexes to characterize dihedrals
!!   |
!!   | rshift(3,nrshift)= Shift in xred that must be done to find
!!   |                    all neighbors of a given atom within a
!!   |                    given number of neighboring shells
!! natom = Number of atoms
!! rprimd(3,3) = Dimensional real space primitive translations
!!               (bohr)
!! xcart(3,natom) = Cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! bt_inv_matrix(3*(natom-1),3*natom) = Inverse of B^{T} matrix
!! deloc_int(3*(natom-1)) = Delocalized internal coordinates
!! prim_int(ninternal) = Primitive internal coordinates
!!
!! SIDE EFFECTS
!! u_matrix(ninternal,3*(natom-1)) = Eigenvectors of BB^T matrix
!!
!! NOTES
!!
!! PARENTS
!!      deloc2xcart,pred_delocint,xfh_recover_deloc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xcart2deloc(deloc,natom,rprimd,xcart,&
& bt_inv_matrix,u_matrix,deloc_int,prim_int)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcart2deloc'
 use interfaces_45_geomoptim, except_this_one => xcart2deloc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(inout) :: deloc
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
 real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
 real(dp),intent(out) :: deloc_int(3*(natom-1))
 real(dp),intent(out) :: prim_int(deloc%ninternal)

!Local variables-------------------------------
!scalars
integer :: ii
logical :: DEBUG=.FALSE.
!arrays
 real(dp) :: b_matrix(deloc%ninternal,3*natom)

! ******************************************************************

 call calc_prim_int(deloc,natom,rprimd,xcart,prim_int)
 if (DEBUG)then
   write(std_out,*) 'Primitive Internals'
   do ii=1,deloc%ninternal
     write(std_out,*) prim_int(ii)
   end do
 end if

 call calc_b_matrix(deloc,natom,rprimd,xcart,b_matrix)
 if (DEBUG)then
   write(std_out,*) 'B Matrix'
   do ii=1,deloc%ninternal
     write(std_out,*) b_matrix(:,ii)
   end do
 end if

 call calc_btinv_matrix(b_matrix,natom,deloc%ninternal,&
& bt_inv_matrix,u_matrix)
 if (DEBUG)then
   write(std_out,*) 'BT Inverse Matrix'
   do ii=1,3*natom
     write(std_out,*) bt_inv_matrix(:,ii)
   end do
 end if

!calculate value of delocalized internals

 call dgemv('T',deloc%ninternal,3*(natom-1),one,&
& u_matrix,deloc%ninternal,prim_int,1,zero,deloc_int,1)

end subroutine xcart2deloc
!!***


!!****f* ABINIT/calc_btinv_matrix
!! NAME
!! calc_btinv_matrix
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      xcart2deloc
!!
!! NOTES
!!   bt_inv_matrix is inverse transpose of the delocalized
!!    coordinate B matrix. b_matrix is the primitive internal B matrix
!!
!! CHILDREN
!!
!! SOURCE

 subroutine calc_btinv_matrix(b_matrix,natom,ninternal,bt_inv_matrix,u_matrix)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_btinv_matrix'
 use interfaces_45_geomoptim, except_this_one => calc_btinv_matrix
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ninternal,natom
 real(dp),intent(in) :: b_matrix(ninternal,3*natom)
 real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))

!Local variables ------------------------------------
!scalars
 integer :: ii,info,lwork
!arrays
 real(dp) :: f_eigs(3*natom),f_matrix(3*natom,3*natom)
 real(dp) :: s_matrix(3*natom,3*natom)
 real(dp) :: s_red(3*natom,3*(natom-1))
 real(dp) :: u_matrix_old(ninternal,3*(natom-1))
 real(dp),allocatable :: work(:)

!******************************************************************

!f matrix = B^{T} B
 call dgemm('T','N',3*natom,3*natom,ninternal,one,&
& b_matrix,ninternal,b_matrix,ninternal,zero,f_matrix,3*natom)

 lwork = max(1,3*3*natom-1)
 ABI_ALLOCATE(work,(lwork))
 s_matrix(:,:) = f_matrix(:,:)

 call dsyev('V','L',3*natom,s_matrix,3*natom,&
& f_eigs,work,lwork,info)

 ABI_DEALLOCATE(work)

 if (abs(f_eigs(1)) + abs(f_eigs(2)) + abs(f_eigs(3)) > tol10 ) then
   write(std_out,*) 'Error: 3 lowest eigenvalues are not zero'
   write(std_out,*) '  internal coordinates do NOT span the full degrees of freedom !'
   write(std_out,'(6E16.6)') f_eigs
   MSG_ERROR("Aborting now")
 end if
 if ( abs(f_eigs(4)) < tol10 ) then
   write(std_out,*) 'Error: fourth eigenvalue is zero'
   write(std_out,*) '  internal coordinates do NOT span the full degrees of freedom !'
   write(std_out,'(6E16.6)') f_eigs
   MSG_ERROR("Aborting now")
 end if

!calculate U matrix from U = B * S_red * lambda^{-1/2}
 do ii=1,3*(natom-1)
   s_red(:,ii) = s_matrix(:,ii+3)/sqrt(f_eigs(ii+3))
 end do

 u_matrix_old(:,:) = u_matrix(:,:)

 call dgemm('N','N',ninternal,3*(natom-1),3*natom,one,&
& b_matrix,ninternal,s_red,3*natom,zero,u_matrix,ninternal)


!align eigenvectors, to preserve a form of continuity in convergences
!!!! eigenvalues are no longer in increasing order!!! but only s_red is reordered
!so that btinv is correct.
 call align_u_matrices(natom,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)

!calculate B_deloc^{-1} matrix for transformation of forces to deloc coord.
!(B^{T}_deloc)^{-1} = (B_deloc B^{T}_deloc)^{-1} B_deloc = lambda^{-3/2} S^{T} F
!= ( S lambda^{3/2} )^{T} F

!! DEFINITION
!! real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)

!even better: B_deloc^{-1} = lambda^{-1/2} S^{T}
 do ii=1,3*(natom-1)
!  s_red(:,ii) = s_matrix(:,ii+3)*sqrt(f_eigs(ii+3))
   bt_inv_matrix(ii,:) = s_matrix(:,ii+3)/sqrt(f_eigs(ii+3))
 end do

end subroutine calc_btinv_matrix
!!***

!!****f* ABINIT/align_u_matrices
!! NAME
!! align_u_matrices
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      xcart2deloc
!!
!! CHILDREN
!!
!! SOURCE

 subroutine align_u_matrices(natom,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'align_u_matrices'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ninternal,natom
!arrays
 real(dp),intent(in) :: u_matrix_old(ninternal,3*(natom-1))
 real(dp),intent(inout) :: f_eigs(3*natom)
 real(dp),intent(inout) :: s_matrix(3*natom,3*natom)
 real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))

!Local variables ------------------------------
!scalars
 integer :: ii,iint1,imax
 real(dp) :: ss
!arrays
 integer :: eigv_flag(3*(natom-1)),eigv_ind(3*(natom-1))
 real(dp) :: tmps(3*natom,3*natom)
 real(dp) :: tmpu(ninternal,3*(natom-1))
 real(dp) :: tmpf(3*natom)

!******************************************************************

 eigv_flag(:) = 0
 eigv_ind(:) = 0

!just permit a change in sign
 do iint1=1,3*(natom-1)
   ss = zero
   do ii=1,ninternal
     ss = ss + u_matrix_old(ii,iint1)*u_matrix(ii,iint1)
   end do
   if (ss < -tol12) then
     imax = -iint1
   else
     imax = iint1
   end if
   eigv_ind(iint1) = imax
   eigv_flag(abs(imax)) = 1
 end do

 tmpu(:,:) = u_matrix
 tmps(:,:) = s_matrix
 tmpf(:) = f_eigs
!exchange eigenvectors...
 do iint1=1,3*(natom-1)
   ss = one
   if (eigv_ind(iint1) < 0) ss = -one

   imax = abs(eigv_ind(iint1))

   tmpu(:,imax) = ss*u_matrix(:,iint1)

   tmps(:,imax+3) = ss*s_matrix(:,iint1+3)

   tmpf(imax+3) = f_eigs(iint1+3)
 end do

 u_matrix(:,:) = tmpu(:,:)
 s_matrix(:,:) = tmps(:,:)
 f_eigs(:) = tmpf(:)

end subroutine align_u_matrices
!!***
