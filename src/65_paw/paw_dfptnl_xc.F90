!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_dfptnl_xc
!! NAME
!! paw_dfptnl_xc
!!
!! FUNCTION
!! Compute a contribution of the third derivative of XC energy of ONE PAW spheres Om_a.
!! It is equal to:
!!   E_at(kxc,rho1,rho2,rho3) = Int_{Om_a} dr kxc(r) * rho1(r) * rho2(r) * rho3(r)
!! where kxc,rho1,rho2 and rho3 are inputs.
!! This routine is similar to m_pawxc.F90:pawxc_dfpt(...) but is implemented independently
!! in order to not overload the original routine.
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!!
!! COPYRIGHT
!! Copyright (C) 2018-2018 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex_1-2-3= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  d3exc1_iat=third-order derivative to compute
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nrad,pawang%angl_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect1-2-3(lm_size)=select the non-zero LM-moments of input density rhor1-2-3
!!  nhat1-2-3(cplex_den*nrad,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor1-2-3(cplex_den*nrad,lm_size,nspden)=first-order change of density
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!
!! OUTPUT
!!  d3exc1_iat = E_at(kxc,rho1,rho2,rho3) (see FUNCTION above)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      paw_dfptnl_energy
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,kxc,lm_size,lmselect1,lmselect2,lmselect3,&
&                 nhat1,nhat2,nhat3,nkxc,nrad,nspden,pawang,pawrad,rhor1,rhor2,rhor3,usexcnhat)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type,simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_dfptnl_xc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_1,cplex_2,cplex_3,ixc,lm_size,nkxc,nrad,nspden,usexcnhat
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect1(lm_size),lmselect2(lm_size),lmselect3(lm_size)
 real(dp),intent(out) :: d3exc1_iat(2)
 real(dp),intent(in) :: kxc(nrad,pawang%angl_size,nkxc)
 real(dp),intent(in) :: nhat1(cplex_1*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat2(cplex_2*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat3(cplex_3*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor1(cplex_1*nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor2(cplex_2*nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor3(cplex_3*nrad,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ilm,ipts,ispden,lm_size_eff,npts
 real(dp) :: d3exc1_int,rho1u,rho1d,rho2u,rho2d,rho3u,rho3d
 character(len=500) :: msg
!arrays
! real(dp) :: tsec(2)
 real(dp),allocatable :: ff(:),rho1arr(:,:),rho2arr(:,:),rho3arr(:,:)

! *************************************************************************

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

 npts=pawang%angl_size
 lm_size_eff=min(lm_size,pawang%ylm_size)

 d3exc1_iat(:) = zero

!Special case: no XC applied
 if (ixc==0) then
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

 ABI_ALLOCATE(rho1arr,(cplex_1*nrad,nspden))
 ABI_ALLOCATE(rho2arr,(cplex_2*nrad,nspden))
 ABI_ALLOCATE(rho3arr,(cplex_3*nrad,nspden))

!Restriction : all cplex must be 1
 if (cplex_1/=1.or.cplex_2/=1.or.cplex_3/=1) then
   msg='All cplex must be one (for the moment...)'
   MSG_BUG(msg)
 end if
!Restriction : nspden must be 1
! if (nkxc>1) then
!   msg='nkxc must be one (<=> nspden=1) (for the moment...)'
!   MSG_BUG(msg)
! end if

 ABI_ALLOCATE(ff,(nrad))

!!----------------------------------------------------------------------
!!----- Loop on the angular part and inits
!!----------------------------------------------------------------------

!Do loop on the angular part (theta,phi)
 do ipts=1,npts

!  Copy the input 1st-order density for this (theta,phi) - PERT1
   rho1arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&         +rhor1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&         +(rhor1(:,ilm,ispden)+nhat1(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  Copy the input 1st-order density for this (theta,phi) - PERT2
   rho2arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect2(ilm)) rho2arr(:,ispden)=rho2arr(:,ispden) &
&         +rhor2(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect2(ilm)) rho2arr(:,ispden)=rho2arr(:,ispden) &
&         +(rhor2(:,ilm,ispden)+nhat2(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  Copy the input 1st-order density for this (theta,phi) - PERT3
   rho3arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect3(ilm)) rho3arr(:,ispden)=rho3arr(:,ispden) &
&         +rhor3(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect3(ilm)) rho3arr(:,ispden)=rho3arr(:,ispden) &
&         +(rhor3(:,ilm,ispden)+nhat3(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store 3nd-order change of XC energy
!  ----------------------------------------------------------------------

   if (cplex_1==1.and.cplex_2==1.and.cplex_3==1) then ! all cplex are 1 :
     if (nspden==1) then
       ff(:)=kxc(:,ipts,1)*rho1arr(:,1)*rho2arr(:,1)*rho3arr(:,1)
     else if (nspden==2) then
       do ii=1,nrad
         rho1u=rho1arr(ii,2)
         rho1d=rho1arr(ii,1)-rho1arr(ii,2)
         rho2u=rho2arr(ii,2)
         rho2d=rho2arr(ii,1)-rho2arr(ii,2)
         rho3u=rho3arr(ii,2)
         rho3d=rho3arr(ii,1)-rho3arr(ii,2)
         ff(ii)=&
!          uuu                                uud
&          kxc(ii,ipts,1)*rho1u*rho2u*rho3u + kxc(ii,ipts,2)*rho1u*rho2u*rho3d + &
!          udu                                udd
&          kxc(ii,ipts,2)*rho1u*rho2d*rho3u + kxc(ii,ipts,3)*rho1u*rho2d*rho3d + &
!          duu                                dud
&          kxc(ii,ipts,2)*rho1d*rho2u*rho3u + kxc(ii,ipts,3)*rho1d*rho2u*rho3d + &
!          ddu                                ddd
&          kxc(ii,ipts,3)*rho1d*rho2d*rho3u + kxc(ii,ipts,4)*rho1d*rho2d*rho3d
       end do
     end if
   end if

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(d3exc1_int,ff,pawrad)
   d3exc1_iat(1)=d3exc1_iat(1)+d3exc1_int*pawang%angwgth(ipts)

!  ----- End of the loop on npts (angular part)
 end do

 d3exc1_iat = d3exc1_iat*four_pi

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(rho1arr)
 ABI_DEALLOCATE(rho2arr)
 ABI_DEALLOCATE(rho3arr)

end subroutine paw_dfptnl_xc
!!***
