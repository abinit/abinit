!{\src2tex{textfont=tt}}
!!****f* ABINIT/qijb_kk
!! NAME
!! qijb_kk
!!
!! FUNCTION
!! Routine which computes PAW onsite part of wavefunction overlap for Bloch
!! functions at two k-points k and k+b. These
!! quantities are used in PAW-based computations of polarization and magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dkvecs(3) :: $\Delta k$ input vector
!!  expibi(2,my_natom,3) :: phase factors at each atomic site for given k offset
!!  gprimd(3,3)=dimensioned primitive translations of reciprocal lattice
!!  lmn2max :: lmnmax*(lmnmax+1)/2
!!  natom=number of atoms in unit cell
!!  ntypat=number of types of atoms in unit cell
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat=typat(natom) list of atom types
!!
!! OUTPUT
!!  calc_qijb(2,lmn2max,natom) :: PAW on-site overlaps of wavefunctions at neighboring
!!                                   k point
!!
!! SIDE EFFECTS
!!
!! NOTES
!! this function computes the on-site data for the PAW version of
!! <u_nk|u_mk+b>, that is, two Bloch vectors at two different k points.
!!
!! PARENTS
!!      initberry,overlap_k1k2_paw
!!
!! CHILDREN
!!      initylmr,sbf8,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine qijb_kk(calc_qijb,dkvecs,expibi,gprimd,lmn2max,natom,ntypat,&
&                   pawang,pawrad,pawtab,typat)

 use m_profiling_abi

 use defs_basis
 use m_errors

 use m_xmpi, only : xmpi_sum
 use m_special_funcs, only : sbf8
 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, simp_gen
 use m_pawtab, only : pawtab_type
 use m_paw_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'qijb_kk'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: lmn2max,natom,ntypat
 type(pawang_type),intent(in) :: pawang
 real(dp),intent(out) :: calc_qijb(2,lmn2max,natom)
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: dkvecs(3),expibi(2,natom),gprimd(3,3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: iatom,ir,isel,itypat
 integer :: klm,kln,klmn,lbess,lbesslm,lmin,lmax,mbess,mesh_size
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: arg,bessg,bnorm,intg,rterm
 complex(dpc) :: cterm,etb,ifac
!arrays
 real(dp) :: bb(3),bbn(3),bcart(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),ylmb(:),sb_out(:)
! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/cone,j_dpc,-cone,-j_dpc/)

! *************************************************************************

 calc_qijb(:,:,:) = zero

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

 ABI_ALLOCATE(sb_out, (pawang%l_size_max))

 do iatom = 1, natom

   itypat = typat(iatom)
   mesh_size = pawtab(itypat)%mesh_size

   ABI_ALLOCATE(j_bessel,(mesh_size,pawang%l_size_max))
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(ylmb,(pawang%l_size_max*pawang%l_size_max))

   !    here is exp(-i b.R) for current atom: recall storage in expibi
   etb = cmplx(expibi(1,iatom),expibi(2,iatom))

   !    note the definition used for the k-dependence of the PAW basis functions:
   !$|\phi_{i,k}\rangle = exp(-i k\cdot r)|\phi_i\rangle
   !    see Umari, Gonze, and Pasquarello, PRB 69,235102 Eq. 23. Thus the k-vector on the
   !    bra side enters as k, while on the ket side it enters as -k.
   bb(:) = -dkvecs(:)

   !    reference bb to cartesian axes
   bcart(1:3)=MATMUL(gprimd(1:3,1:3),bb(1:3))

   !    bbn is b-hat (the unit vector in the b direction)
   bnorm=dsqrt(dot_product(bcart,bcart))
   bbn(:) = bcart(:)/bnorm

   !    as an argument to the bessel function, need 2pi*b*r = 1 so b is re-normed to two_pi
   bnorm = two_pi*bnorm
   do ir=1,mesh_size
     arg=bnorm*pawrad(itypat)%rad(ir)
     call sbf8(pawang%l_size_max,arg,sb_out) ! spherical bessel functions at each mesh point
     j_bessel(ir,:) = sb_out
   end do ! end loop over mesh

   !    compute Y_LM(b) here
   call initylmr(pawang%l_size_max,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,bbn,ylmb(:),ylmgr)

   do klmn = 1, pawtab(itypat)%lmn2_size
     klm =pawtab(itypat)%indklmn(1,klmn)
     kln =pawtab(itypat)%indklmn(2,klmn)
     lmin=pawtab(itypat)%indklmn(3,klmn)
     lmax=pawtab(itypat)%indklmn(4,klmn)
     do lbess = lmin, lmax, 2    ! only possible choices for L s.t. Gaunt integrals
                                  !        will be non-zero
       ifac = il(mod(lbess,4))
       do mbess = -lbess, lbess
         lbesslm = lbess*lbess+lbess+mbess+1
         isel=pawang%gntselect(lbesslm,klm)
         if (isel > 0) then
           bessg = pawang%realgnt(isel)
           ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&           -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&           *j_bessel(1:mesh_size,lbess+1)
           call simp_gen(intg,ff,pawrad(itypat))
           rterm = four_pi*bessg*intg*ylmb(lbesslm)
           cterm = etb*ifac*rterm
           calc_qijb(1,klmn,iatom) = &
&           calc_qijb(1,klmn,iatom) + dreal(cterm)
           calc_qijb(2,klmn,iatom) = &
&           calc_qijb(2,klmn,iatom) + dimag(cterm)

         end if ! end selection on non-zero Gaunt factors
       end do ! end loop on mbess = -lbess, lbess
     end do ! end loop on lmin-lmax bessel l values
   end do ! end loop on lmn2_size klmn basis pairs

   ABI_DEALLOCATE(j_bessel)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ylmb)
 end do ! end loop over atoms

 ABI_DEALLOCATE(sb_out)

 end subroutine qijb_kk
!!***

