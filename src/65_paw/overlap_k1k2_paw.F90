!{\src2tex{textfont=tt}}
!!****f* ABINIT/overlap_k1k2_paw
!! NAME
!! overlap_k1k2_paw
!!
!! FUNCTION
!! compute PAW overlap between two k points,
!! similar to smatrix_k_paw.F90 but more generic
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj_k1 (pawcprj_type) :: cprj for occupied bands at point k1
!!  cprj_k2 :: cprj for occupied bands at point k2
!!  dk(3) :: vector k2 - k1
!!  gprimd(3,3)=dimensioned primitive translations of reciprocal lattice
!!  lmn2max :: lmnmax*(lmnmax+1)/2
!!  lmnsize(ntypat) :: lmnsize for each atom type
!!  mband :: number of bands
!!  natom=number of atoms in unit cell
!!  nspinor :: number of spinors (1 or 2)
!!  ntypat=number of types of atoms in unit cell
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat=typat(natom) list of atom types
!!  xred(natom,3) :: locations of atoms in cell
!!
!! OUTPUT
!! k1k2_paw(2,mband,mband) :: array of the on-site PAW parts of the overlaps between Bloch states at points
!!   k1 and k2, for the various pairs of bands, that is, the on-site part of 
!!   <u_nk1|u_mk2>
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine assumes that the cprj are not explicitly ordered by 
!! atom type.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine overlap_k1k2_paw(cprj_k1,cprj_k2,dk,gprimd,k1k2_paw,lmn2max,lmnsize,mband,&
&                           natom,nspinor,ntypat,pawang,pawrad,pawtab,typat,xred)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_pawang, only : pawang_type
 use m_pawcprj, only : pawcprj_type
 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'overlap_k1k2_paw'
 use interfaces_65_paw, except_this_one => overlap_k1k2_paw
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: lmn2max,mband,natom,nspinor,ntypat
 type(pawang_type),intent(in) :: pawang
 type(pawcprj_type),intent(in) :: cprj_k1(natom,mband),cprj_k2(natom,mband)

!arrays
 integer,intent(in) :: lmnsize(ntypat),typat(natom)
 real(dp),intent(in) :: dk(3),gprimd(3,3),xred(natom,3)
 real(dp),intent(out) :: k1k2_paw(2,mband,mband)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: iatom,iband,ibs,ilmn,ispinor,itypat
 integer :: jband,jbs,jlmn,klmn
 complex(dpc) :: cpk1,cpk2,cterm,paw_onsite

 ! arrays
 real(dp),allocatable :: calc_expibi(:,:),calc_qijb(:,:,:)

! *************************************************************************

!initialize k1k2_paw output variable
 k1k2_paw(:,:,:) = zero

 ! obtain the atomic phase factors for the input k vector shift
 ABI_ALLOCATE(calc_expibi,(2,natom))
 call expibi(calc_expibi,dk,natom,xred)

 ! obtain the onsite PAW terms for the input k vector shift
 ABI_ALLOCATE(calc_qijb,(2,lmn2max,natom))
 call qijb_kk(calc_qijb,dk,calc_expibi,gprimd,lmn2max,natom,ntypat,pawang,pawrad,pawtab,typat)
 ABI_DEALLOCATE(calc_expibi)

 do iatom = 1, natom
    itypat = typat(iatom)

    do ilmn=1,lmnsize(itypat)
       do jlmn=1,lmnsize(itypat)
          klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
          paw_onsite = cmplx(calc_qijb(1,klmn,iatom),calc_qijb(2,klmn,iatom))
          do iband = 1, mband
             do jband = 1, mband
                do ispinor = 1, nspinor
                   ibs = nspinor*(iband-1) + ispinor
                   jbs = nspinor*(jband-1) + ispinor
                   cpk1=cmplx(cprj_k1(iatom,ibs)%cp(1,ilmn),cprj_k1(iatom,ibs)%cp(2,ilmn))
                   cpk2=cmplx(cprj_k2(iatom,jbs)%cp(1,jlmn),cprj_k2(iatom,jbs)%cp(2,jlmn))
                   cterm = conjg(cpk1)*paw_onsite*cpk2
                   k1k2_paw(1,iband,jband) = k1k2_paw(1,iband,jband)+real(cterm)
                   k1k2_paw(2,iband,jband) = k1k2_paw(2,iband,jband)+aimag(cterm)
                end do ! end loop over ispinor
             end do ! end loop over jband
          end do ! end loop over iband
       end do ! end loop over ilmn
    end do ! end loop over jlmn

 end do ! end loop over atoms

 ABI_DEALLOCATE(calc_qijb)

 end subroutine overlap_k1k2_paw
!!***

