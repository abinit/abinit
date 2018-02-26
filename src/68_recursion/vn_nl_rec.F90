!{\src2tex{textfont=tt}}
!!****f* ABINIT/vn_nl_rec
!! NAME
!! vn_nl_rec
!! 
!! FUNCTION
!! this routine computes the contribution to the vector vn, during
!! recursion, due to the non-local psp.
!!
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  vn(:,:,:)=the vector on the real-space grid.
!!  inf_ucvol=volume of infinitesimal cell
!!  natom=number of atoms
!!  typat(natom)=the type of psps associated to the atoms
!!  ngfftrec(3)=first 3 components of ngfftrec (truncated box, if different from ngfft) for the real-space grid
!!  nlrec<type(nlpsprec_type)> in recursion_type containing information concerning psp 
!!  projec(ngfftrec(1),ngfftrec(2),ngfftrec(3),lmnmax,natom) is the  vector, on the ngfftrec grid containing 
!!  the non-lacal projector $Y_{lm}(r-R_A)f_{lk}(r-R_A)
!!
!! OUTPUT
!! vn_nl(:,:,:)=the non_local contribution to vn
!!  
!! PARENTS
!!      recursion,recursion_nl
!!
!! CHILDREN
!!      timab
!!
!! NOTES 
!! 
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vn_nl_rec(vn,natom,typat,ngfftrec,&
  &                    inf_ucvol,nlrec,projec)

 use m_profiling_abi
 
 use defs_basis
 use defs_rectypes
 use m_per_cond,only             :per_cond
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vn_nl_rec'
 use interfaces_18_timing
!End of the abilint section

 implicit none
 
!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom
 real(dp),intent(in) :: inf_ucvol
 type(nlpsprec_type),intent(in) :: nlrec
!arrays
 integer,intent(in) :: ngfftrec(3),typat(natom)
 real(dp),intent(in) :: projec(0:,0:,0:,1:,1:)
 real(dp),intent(inout):: vn(0:ngfftrec(1)*ngfftrec(2)*ngfftrec(3)-1)
!Local variables-------------------------------
!scalars
 integer :: iatom,nfftrec
 integer :: jlmn,il,in,jn
 integer :: ipsp,ilmn
 integer :: npsp,lmnmax
 real(dp):: vn_nl_loc 
!arrays
 real(dp):: vn_nl(0:ngfftrec(1)-1,0:ngfftrec(2)-1,0:ngfftrec(3)-1)
 real(dp):: vtempo(0:ngfftrec(1)-1,0:ngfftrec(2)-1,0:ngfftrec(3)-1)
 real(dp):: tsec(2)
! *************************************************************************

 call timab(615,1,tsec)
!--Initialisation
 
 vn_nl = zero
 npsp = nlrec%npsp
 lmnmax = nlrec%lmnmax
 nfftrec = product(ngfftrec)
 vtempo(:,:,:) = reshape(source=vn,shape=ngfftrec(:3))

!--Sum_iatom \int dr1 E(r-r_a,r1-r_a)vn(r1) *infucvol
 do iatom=1,natom !--Loop on atoms
   ipsp = typat(natom)

!  --If psp(typat(iatom)) is local then cycle
   if(all(nlrec%pspinfo(:,ipsp)==0))  cycle 

   
!  write(std_out,*)'lmnmax',nlrec%lmnmax,lmnmax

   do ilmn = 1, lmnmax
     do jlmn = 1,lmnmax
       if(nlrec%indlmn(4,ilmn,ipsp)==nlrec%indlmn(4,jlmn,ipsp)) then
         il = 1+nlrec%indlmn(1,jlmn,ipsp)
         in = nlrec%indlmn(3,ilmn,ipsp)
         jn = nlrec%indlmn(3,jlmn,ipsp)
         vn_nl_loc = ddot(nfftrec,projec(:,:,:,jlmn,iatom),1,vtempo,1)
         vn_nl = vn_nl+projec(:,:,:,ilmn,iatom)*vn_nl_loc*nlrec%mat_exp_psp_nl(in,jn,il,ipsp)
       end if
     end do
   end do
 end do !--End loop on atoms
 vtempo = vtempo + vn_nl*inf_ucvol

 vn = reshape(source=vtempo,shape=(/nfftrec/))

 call timab(615,2,tsec)
end subroutine vn_nl_rec
!!***
