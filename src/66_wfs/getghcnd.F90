!{\src2tex{textfont=tt}}
!!****f* ABINIT/getghcnd
!!
!! NAME
!! getghcnd
!!
!! FUNCTION
!! Compute <G|H_ND|C> for input vector |C> expressed in reciprocal space
!! Result is put in array ghcnc. H_ND is the Hamiltonian due to magnetic dipoles
!! on the nuclear sites.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!! my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!! ndat=number of FFT to do in //
!!
!! OUTPUT
!! ghcnd(2,npw*my_nspinor*ndat)=matrix elements <G|H_ND|C>
!!
!! SIDE EFFECTS
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!!
!! NOTES
!! Application of <k^prime|H|k> or <k|H|k^prime> not implemented!
!!
!! PARENTS
!!      getghc
!!
!! CHILDREN
!!      matr3inv,zhpmv
!!
!! NOTES
!!  This routine applies the Hamiltonian due to an array of magnetic dipoles located
!!  at the atomic nuclei to the input wavefunction.  The original term is $H_ND = -q_e/m_e A.p$,
!!  where $A=(\mu_0/4\pi) m \times (r-I)/|r-I|^3$ is the vector potential due to a magnetic dipole
!!  of strength $m$ (vector input in cartesian space) located at atomic site $I$. This
!!  term is expressed in reciprocal space to obtain
!!  $ <k+G'| H_ND(k) |k+G> = (\mu_0/4\pi)*4\pi i/ucvol exp(i\Delta G.I) m.\Delta G \times (k+G)/|\Delta G|^2$,
!!  where $\Delta G$ is G - G', q_e = -1, and m_e = +1. Note the "non-locality" in k space:
!!  output at G' has contributions from all input G. Strategy below is to take advantage of
!!  Hermiticity to store H_ND in triangular form and then use a BLAS call to zhpmv to apply to
!!  input vector in one shot.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getghcnd(cwavef,ghcnd,gs_ham,my_nspinor,ndat)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getghcnd'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_nspinor,ndat
 type(gs_hamiltonian_type),intent(in),target :: gs_ham
!arrays
 real(dp),intent(in) :: cwavef(2,gs_ham%npw_k*my_nspinor*ndat)
 real(dp),intent(out) :: ghcnd(2,gs_ham%npw_k*my_nspinor*ndat)

!Local variables-------------------------------
!scalars
 integer :: atom_nd_tot,cwavedim,hggindx,iatom,iatom_nd_tot,ig,igp
 real(dp) :: ldgvec,permeability
 complex(dpc) :: cterm,imfac,phgvec,phgvecp
 character(len=500) :: message
!arrays
 integer,allocatable :: atom_nd(:)
 real(dp) :: crossp(3),dgvec(3),kplusg(3),mured(3),rprimd(3,3)
 complex(dpc),allocatable :: hgg(:),hggc(:),inwave(:)

! *********************************************************************

 if (gs_ham%matblk /= gs_ham%natom) then
   write(message,'(a,i4,a,i4)')' gs_ham%matblk = ',gs_ham%matblk,' but natom = ',gs_ham%natom
   MSG_ERROR(message)
 end if
 if (ndat /= 1) then
   write(message,'(a,i4,a)')' ndat = ',ndat,' but getghcnd requires ndat = 1'
   MSG_ERROR(message)
 end if
 if (my_nspinor /= 1) then
   write(message,'(a,i4,a)')' nspinor = ',my_nspinor,' but getghcnd requires nspinor = 1'
   MSG_ERROR(message)
 end if
 if (any(abs(gs_ham%kpt_k(:)-gs_ham%kpt_kp(:))>tol8)) then
   message=' not allowed for kpt(left)/=kpt(right)!'
   MSG_BUG(message)
 end if

! will need rprimd
 call matr3inv(gs_ham%gprimd,rprimd)

! magnetic permeability mu_0/four_pi in atomic units
! this constant is also used in m_pawdij.F90/pawdijnd, if you change it here,
! change it there also for consistency
 permeability=5.325135453D-5

! size of incoming wavefunction
 cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ABI_ALLOCATE(hgg,(cwavedim*(cwavedim+1)/2))
 ABI_ALLOCATE(hggc,(cwavedim))
 ABI_ALLOCATE(inwave,(cwavedim))

! construct list of atoms that have non-zero magnetic dipole moments
 ABI_ALLOCATE(atom_nd,(gs_ham%natom))
 atom_nd_tot = 0
 do iatom = 1, gs_ham%natom
   if (any(abs(gs_ham%nucdipmom(:,iatom))>tol8)) then
     atom_nd_tot = atom_nd_tot + 1
     atom_nd(atom_nd_tot) = iatom
   end if
 end do

 imfac = cmplx(zero,permeability*four_pi/(gs_ham%ucvol*gs_ham%ucvol),kind=dpc)

 hgg = czero
! loop over only atoms that have non-zero dipoles
 do iatom_nd_tot = 1, atom_nd_tot
   iatom = atom_nd(iatom_nd_tot)

!  nucdipmom is input w.r.t cartesian components. However, the G and k data used
!  below are given to us in reduced units with respect to the primitive translations
!  in reciprocal space. The matrix B (otherwise known as gprimd in abinit) converts
!  back to Cartesian coordinates as k_cart = B.k_red and so forth. We can use some
!  vector identities: $\Delta G \times (k+G) = B\Delta G_red \times B(k+G)_red =
!  det(B)(B^{-1,T})\Delta G_red\times(k+G)_red$. det(B) = 1/ucvol, accounting for
!  the second factor in imfac. Finally, as we need
!  $m.(B^{-1,T})\Delta G_red\times(k+G)_red$, we can save time by applying B^{-1,T}
!  to the left on m rather than at each step to $\Delta G_red\times (k+G)_red$ on the right. Thus
!  we need m^T(B^{-1,T}) = B^{-1}m = A^Tm = rprimd^T m.
   mured(1) = rprimd(1,1)*gs_ham%nucdipmom(1,iatom)+&
&   rprimd(2,1)*gs_ham%nucdipmom(2,iatom)+&
&   rprimd(3,1)*gs_ham%nucdipmom(3,iatom)
   mured(2) = rprimd(1,2)*gs_ham%nucdipmom(1,iatom)+&
&   rprimd(2,2)*gs_ham%nucdipmom(2,iatom)+&
&   rprimd(3,2)*gs_ham%nucdipmom(3,iatom)
   mured(3) = rprimd(1,3)*gs_ham%nucdipmom(1,iatom)+&
&   rprimd(2,3)*gs_ham%nucdipmom(2,iatom)+&
&   rprimd(3,3)*gs_ham%nucdipmom(3,iatom)

   do igp=1, gs_ham%npw_k

     phgvecp = cmplx(gs_ham%ph3d_k(1,igp,iatom),-gs_ham%ph3d_k(2,igp,iatom),kind=dpc)

     hggindx = (igp-1)*igp/2
     do ig = 1,igp-1 ! ig = igp term is strictly zero
       hggindx = hggindx+1

       dgvec(1:3) = gs_ham%kg_k(1:3,ig) - gs_ham%kg_k(1:3,igp)

       ldgvec = DOT_PRODUCT(dgvec,MATMUL(gs_ham%gmet,dgvec))

       phgvec = cmplx(gs_ham%ph3d_k(1,ig,iatom),gs_ham%ph3d_k(2,ig,iatom),kind=dpc)
       kplusg = gs_ham%kpt_k(1:3) + gs_ham%kg_k(1:3,ig)
       crossp(1) =  dgvec(2)*kplusg(3) - dgvec(3)*kplusg(2)
       crossp(2) = -dgvec(1)*kplusg(3) + dgvec(3)*kplusg(1)
       crossp(3) =  dgvec(1)*kplusg(2) - dgvec(2)*kplusg(1)

!      cterm = <igp|H|ig>
       cterm = imfac*phgvecp*phgvec*DOT_PRODUCT(mured,crossp)/ldgvec

!      this is a sum because multiple atoms can contribute to each <G'|H|G> element
       hgg(hggindx) = hgg(hggindx) + cterm


     end do  ! end loop over ig

   end do  ! end loop over igp

 end do !end loop over atoms

! load input wavefunction into complex storage
 do igp = 1, gs_ham%npw_k
   inwave(igp) = cmplx(cwavef(1,igp),cwavef(2,igp),kind=dpc)
 end do

! apply hamiltonian hgg to input wavefunction inwave, result in hggc
 call zhpmv('L',cwavedim,cone,hgg,inwave,1,czero,hggc,1)

 do igp=1,gs_ham%npw_k
   ghcnd(1,igp) = dreal(hggc(igp))
   ghcnd(2,igp) = dimag(hggc(igp))
 end do


 ABI_DEALLOCATE(hgg)
 ABI_DEALLOCATE(hggc)
 ABI_DEALLOCATE(inwave)

 ABI_DEALLOCATE(atom_nd)

end subroutine getghcnd
!!***
