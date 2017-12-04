!{\src2tex{textfont=tt}}
!!****f* ABINIT/mknucdipmom_k
!! NAME
!! mknucdipmom_k
!!
!! FUNCTION
!! compute Hamiltonian in reciprocal space due to array of nuclear
!! dipole moments, at a given k point
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  nucdipmom_k(2,npw*(npw+1)/2) = nuclear dipole moment Hamiltonian matrix, in
!!                                 lower diagonal Hermitian packed storage
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mknucdipmom_k(gmet,kg,kpt,natom,nucdipmom,nucdipmom_k,npw,rprimd,ucvol,xred)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mknucdipmom_k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: natom,npw
 real(dp),intent(in) :: ucvol

 !arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),kpt(3),nucdipmom(3,natom),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: nucdipmom_k(2,npw*(npw+1)/2)
 
!Local variables-------------------------------
!scalars
 integer :: atom_nd_tot,col,iatom,ndp_index,row
 real(dp) :: crossfac,dg2,permeability,permfac,phasefac
 !arrays
 integer :: atom_nd(natom)
 real(dp) :: cprod(3),cprod_cart(3),dgp_red(3), gpk_red(3)

! *************************************************************************
 !

! magnetic permeability mu_0/four_pi in atomic units
! this constant is also used in m_pawdij.F90/pawdijnd, if you change it here,
! change it there also for consistency
 permeability=5.325135453D-5
 ! will need 4*pi*i*(\mu_0/four\pi)
 permfac = four_pi*permeability

 ! make list of atoms with non-zero nuclear magnetic dipoles
 atom_nd_tot = 0
 do iatom = 1, natom
    if(any(abs(nucdipmom(:,iatom))>tol8)) then
       atom_nd_tot = atom_nd_tot + 1
       atom_nd(atom_nd_tot) = iatom
    end if
 end do

 ndp_index = 0
 do col=1,npw ! enumerate plane waves G
    ! form k + G at this k point for current plane wave (this is the ket |k+G> )
    ! in reduced coordinates
   gpk_red(1)=dble(kg(1,col))+kpt(1)
   gpk_red(2)=dble(kg(2,col))+kpt(2)
   gpk_red(3)=dble(kg(3,col))+kpt(3)

   do row=col,npw ! enumerate lower diagonal from 1 to G
      ! index of the current matrix element, in lower triangular packed storage
      ! "packed sequentially, column by column"
      ndp_index = ndp_index + 1
      nucdipmom_k(:,ndp_index) = zero
      
      ! form G-G' = \Delta G at this k pt (this is the bra <k+G'| )
      ! in reduced coordinates
      dgp_red(1)=dble(kg(1,col)-kg(1,row))
      dgp_red(2)=dble(kg(2,col)-kg(2,row))
      dgp_red(3)=dble(kg(3,col)-kg(3,row))

      ! compute |\Delta G|^2
      ! must use gmet metric because G's are in reduced coords in reciprocal space
      dg2 = DOT_PRODUCT(dgp_red,MATMUL(gmet,dgp_red))
      ! if \Delta G = 0, Hamiltonian term is zero and move on to next one
      if (abs(dg2)<tol8) then
         nucdipmom_k(1:2,ndp_index)=zero
         cycle
      end if

      ! compute cross product \Delta G \times (k + G)
      ! notice that \Delta G and (k + G) are in reduced coords in reciprocal space
      cprod(1) = dgp_red(2)*gpk_red(3) - dgp_red(3)*gpk_red(2)
      cprod(2) = dgp_red(3)*gpk_red(1) - dgp_red(1)*gpk_red(3)
      cprod(3) = dgp_red(1)*gpk_red(2) - dgp_red(2)*gpk_red(1)

      ! proper cross product must account for reduced coords as follows:
      ! gprimd*dgp \times gprimd*gpk = (det gprimd)*(gprimd^{-1,T})*(dgp \times gpk)
      ! = rprimd * (dgp \times gpk)/ucvol
      ! final vector also includes the division by |\Delta G|^2
      cprod_cart = MATMUL(rprimd,cprod)/(ucvol*dg2)

      ! loop over the atoms with non-zero nuclear dipoles
      ! phase factors exp(i*\Delta G*I) where I is ion position,
      ! might be retrievable from ph1d, need to check
      do iatom = 1, atom_nd_tot
         phasefac = two_pi*DOT_PRODUCT(dgp_red,xred(:,atom_nd(iatom)))
         crossfac = DOT_PRODUCT(nucdipmom(:,iatom),cprod_cart)
         nucdipmom_k(1,ndp_index) = nucdipmom_k(1,ndp_index) - permfac*crossfac*sin(phasefac)
         nucdipmom_k(2,ndp_index) = nucdipmom_k(2,ndp_index) + permfac*crossfac*cos(phasefac)
      end do ! end loop over atoms with nonzero dipoles

   end do ! end loop over G' = G to npw

 end do ! end loop over G = 1 to npw

end subroutine mknucdipmom_k
!!***
