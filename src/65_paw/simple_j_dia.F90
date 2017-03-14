!{\src2tex{textfont=tt}}
!!****f* ABINIT/simple_j_dia
!! NAME
!! simple_j_dia
!!
!! FUNCTION
!! simple test current for H atoms
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT

!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine simple_j_dia(jdia,natom,nfft,pawfgrtab)


 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawfgrtab, only : pawfgrtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'simple_j_dia'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft
!arrays
 real(dp),intent(out) :: jdia(3,3,nfft)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,ifgd,ifftsph
 real(dp) :: Bvecsum,nrm,nrmmax,psi2,scale
!arrays
 real(dp) :: axb(3),B(3,3),rvec(3)
 real(dp),allocatable :: Bvec(:,:)

! ************************************************************************

 DBG_ENTER("COLL")

!make sure current field is set to zero everywhere
 jdia(:,:,:) = zero

!define the directions of the external B field
 B = zero
 B(1,1) = one; B(2,2) = one; B(3,3) = one;

!loop over atoms in cell
 scale  = half
 nrm    = zero
 nrmmax = zero
 Bvecsum = zero
 do iatom = 1, natom
   ABI_ALLOCATE(Bvec,(3,pawfgrtab(iatom)%nfgd))
   do ifgd=1, pawfgrtab(iatom)%nfgd
     ifftsph = pawfgrtab(iatom)%ifftsph(ifgd)
     rvec(:) = pawfgrtab(iatom)%rfgd(:,ifgd)
     nrm = sqrt(dot_product(rvec,rvec))
     if (nrm > nrmmax) nrmmax = nrm
     psi2=exp(-two*nrm/scale)/(pi*scale*scale*scale)
     do idir=1, 3
       call acrossb(B(idir,:),rvec,axb)
       jdia(idir,:,ifftsph) = -psi2*axb(:)/(two*Sp_Lt)
     end do ! end loop over idir for jdia
     if (nrm > zero) then
       call acrossb(rvec,jdia(3,:,ifftsph),axb)
       Bvec(:,ifgd) = axb(:)/(Sp_Lt*nrm*nrm*nrm)
     end if
     Bvecsum = Bvecsum + Bvec(3,ifgd)
   end do ! end loop over nfgd points in sphere
   write(std_out,'(i8,f12.8,f12.8)')pawfgrtab(iatom)%nfgd,Bvecsum,Bvecsum*(four_pi*third*nrm**3)/pawfgrtab(iatom)%nfgd
   ABI_DEALLOCATE(Bvec)
 end do ! end loop over atoms in cell

 DBG_EXIT("COLL")

 end subroutine simple_j_dia
!!***
