!{\src2tex{textfont=tt}}
!!****f* ABINIT/vso_realspace_local
!! NAME
!!   vso_realspace_local
!!
!! FUNCTION
!!
!!  Calculate real space (local - (r,r)) values of the SO part of the
!!   pseudopotential. Reconstructed explicitly in the HGH/GTH case.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group (Mver)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      spin_current
!!
!! CHILDREN
!!      gamma_function,spline,splint,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine vso_realspace_local(dtset,hdr,position_op,psps,vso_realspace)

 use m_profiling_abi

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vso_realspace_local'
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments -------------------------------
   
  type(hdr_type),intent(inout) :: hdr
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type),intent(in) :: psps

  real(dp),intent(in) :: position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3))

  real(dp),intent(out) :: vso_realspace(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),&
      & dtset%nspinor,dtset%nspinor,3)

!Local variables -------------------------

  integer :: i,j,l, lmax,ipsp,iatom, ir1,ir2,ir3
  integer :: rcexponent,irealsp
  integer :: nradgrid,iradgrid

  real(dp) :: gammai, gammaj, relative_position(3), radial_cutoff, norm_rel_pos
  real(dp) :: expfact,lfact, vso_interpol, x,y,z
  real(dp) :: xcart(3,dtset%natom),splint_x(1),splint_y(1)

  real(dp), allocatable :: radial_grid(:)
  real(dp), allocatable :: prefact_ijl(:,:,:,:),tmpvso(:),tmpvso_pp(:)
  real(dp), allocatable :: vso_radial(:,:),vso_radial_pp(:,:),tmp_spline(:)
  real(dp), allocatable :: offdiag_l_fact(:,:,:),kpar_matrix(:,:)

! *********************************************************************

!recalculate xcart (option = 1)
 call xred2xcart(dtset%natom,hdr%rprimd,xcart,hdr%xred)


 lmax = psps%mpsang-1

!content of gth pseudo type:
!These are {rloc, C(1...4)} coefficients for psppar(0, :, :) indices,
!Followed by the h coefficients for psppar(1:2, 1:, :) indices.
!size (0:2, 0:4, npsp)
!potential radius r_l is in psppar(l+1,0,ipsp)
!real(dp), pointer :: psppar(:, :, :)
!The covalence radii for each pseudo (?) size (npsp)
!real(dp), pointer :: radii_cov(:)
!Cut-off radii for core part and long-range part.
!radii_cf(:, 1) is for the long-range cut-off and
!radii_cf(:, 2) is for the core cut-off.
!size (npsp, 2)
!real(dp), pointer :: radii_cf(:, :)
!Spin orbit coefficients in HGH/GTH formats: k11p
!etc... see psp3ini.F90
!dimension = num l channels, 3 coeffs, num psp =
!(1:lmax+1,1:3,npsp)
!real(dp), pointer :: psp_k_par(:, :, :)

!v_SO^l (r,r') = sum_i sum_j sum_m Y_{lm} (\hat{r}) p_i^l (r) k_{ij}^l p_j^l(r') Y^{*}_lm (\hat{r'})
!
!v_SO^l (r,r)  = sum_ij  p_i^l (r) k_{ij}^l p_j^l(r) sum_m Y_{lm} (\hat{r}) Y^{*}_lm (\hat{r})
!= (2l+1)/4\pi sum_ij  p_i^l (r) k_{ij}^l p_j^l(r) (eq B.17 Patrick Rinke thesis)
!p are gaussian projectors (from HGH paper prb 58 3641)
!sum_l v_SO^l (r,r) is a purely radial quantity (function of |r|), so spline it

!maximum distance needed in unit cell
 radial_cutoff = four * maxval(psps%gth_params%psppar(:, 0, :))

!setup radial grid; Should we use a logarithmic grid? The spline functions can
!take it...
 nradgrid = 201 ! this is heuristic
 ABI_ALLOCATE(radial_grid,(nradgrid))
 do iradgrid=1,nradgrid
   radial_grid(iradgrid) = (iradgrid-1)*radial_cutoff/(nradgrid-1)
 end do

!calculate prefactors independent of r
 ABI_ALLOCATE(prefact_ijl,(3,3,0:lmax,psps%npsp))
 ABI_ALLOCATE(offdiag_l_fact,(3,3,0:lmax))
 ABI_ALLOCATE(kpar_matrix,(3,3))

!these factors complete the full 3x3 matrix of k (or h) parameters for the
!HGH pseudos
 offdiag_l_fact = zero
!l=0
 offdiag_l_fact(1,2,0) = -half*sqrt(three/five)
 offdiag_l_fact(1,3,0) = half*sqrt(five/21._dp)
 offdiag_l_fact(2,3,0) = -half*sqrt(100._dp/63._dp)
!l=1
 offdiag_l_fact(1,2,1) = -half*sqrt(five/seven)
 offdiag_l_fact(1,3,1) = sixth*sqrt(35._dp/11._dp)
 offdiag_l_fact(2,3,1) = -sixth*14._dp/sqrt(11._dp)
!l=2
 if (lmax >= 2) then
   offdiag_l_fact(1,2,2) = -half*sqrt(seven/nine)
   offdiag_l_fact(1,3,2) = half*sqrt(63._dp/143._dp)
   offdiag_l_fact(2,3,2) = -half*18._dp /sqrt(143._dp)
 end if
!l=3
 if (lmax >= 3) then
   offdiag_l_fact(1,2,3) = zero
   offdiag_l_fact(1,3,3) = zero
   offdiag_l_fact(2,3,3) = zero
 end if
!get prefactors for evaluation of V_SO: terms that do not depend on r
 prefact_ijl = zero
 do l=0,lmax
!  first the diagonal i=j term
   do i=1,3
     call gamma_function(l+(4._dp*i-1._dp)*0.5_dp, gammai)
     gammai = sqrt(gammai)
     rcexponent = 2*l+2*i+2*i-1
     do ipsp=1,psps%npsp
       prefact_ijl(i,i,l,ipsp) = psps%gth_params%psp_k_par(l+1,i,ipsp) &
&       / ( (psps%gth_params%psppar(l+1,0,ipsp))**(rcexponent) &
&       * gammai * gammai)
     end do
   end do
!  now the off diagonal elements
   do ipsp=1,psps%npsp
     kpar_matrix(1,2) = offdiag_l_fact (1,2,l)* psps%gth_params%psp_k_par(l+1,2,ipsp)
     kpar_matrix(2,1) = kpar_matrix(1,2)
     kpar_matrix(1,3) = offdiag_l_fact (1,3,l)* psps%gth_params%psp_k_par(l+1,3,ipsp)
     kpar_matrix(3,1) = kpar_matrix(1,3)
     kpar_matrix(2,3) = offdiag_l_fact (2,3,l)* psps%gth_params%psp_k_par(l+1,3,ipsp)
     kpar_matrix(3,2) = kpar_matrix(2,3)
   end do

!  for the f case only the 1,1 matrix element is non 0 - it is done above and
!  all these terms are actually 0
   if (l > 2) cycle

   do i=1,3
     call gamma_function(l+(4._dp*i-1._dp)*0.5_dp, gammai)
     gammai = sqrt(gammai)
     do j=1,3
       if (j==i) cycle
       rcexponent = 2*l+2*i+2*j-1
       call gamma_function(l+(4._dp*j-1._dp)*0.5_dp,gammaj)
       gammaj = sqrt(gammaj)
       do ipsp=1,psps%npsp
         prefact_ijl(i,j,l,ipsp) = kpar_matrix(i,j) &
&         / ( (psps%gth_params%psppar(l+1,0,ipsp))**rcexponent &
&         * gammai * gammaj )
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(kpar_matrix)
 ABI_DEALLOCATE(offdiag_l_fact)

 prefact_ijl = prefact_ijl * two

!calculate v_SO on radial grid
! MGNAG Runtime Error: *** Arithmetic exception: Floating invalid operation - aborting
 ABI_ALLOCATE(vso_radial,(nradgrid,psps%npsp))
 vso_radial = zero
 do l=0,lmax
   lfact=(2._dp*l+1._dp)/four/pi
   do iradgrid=1,nradgrid
     norm_rel_pos = radial_grid(iradgrid)
     do ipsp=1,psps%npsp
       expfact = exp(-norm_rel_pos**2 / &
&       (psps%gth_params%psppar(l+1,0,ipsp))**2)

       do i=1,3
         do j=1,3
           rcexponent = 2*l +2*i+2*j-4
           if(prefact_ijl(i,j,l,ipsp)/=0) then !vz_d 0**0
             vso_radial(iradgrid,ipsp) = vso_radial(iradgrid,ipsp) + &
&             prefact_ijl(i,j,l,ipsp)*(norm_rel_pos**rcexponent) * expfact
           end if  !vz_d
         end do ! j
       end do ! i
     end do ! ipsp
   end do ! iradgrid
 end do ! lmax

!spline v_SO(radial coord): get second derivative coefficients 
 ABI_ALLOCATE(vso_radial_pp,(nradgrid,psps%npsp))

 ABI_ALLOCATE(tmp_spline,(nradgrid))
 ABI_ALLOCATE(tmpvso,(nradgrid))
 ABI_ALLOCATE(tmpvso_pp,(nradgrid))
 do ipsp=1,psps%npsp
   tmpvso = vso_radial(:,ipsp)
   call spline( radial_grid, tmpvso, nradgrid, zero, radial_grid(nradgrid), tmpvso_pp )
   vso_radial_pp(:,ipsp) = tmpvso_pp
 end do
 ABI_DEALLOCATE(tmp_spline)
 ABI_DEALLOCATE(tmpvso)
 ABI_DEALLOCATE(tmpvso_pp)

!to optimize this I should precalculate the distances which are actually needed by
!symmetry, or only sum over irreducible points in space and use weights

!for each physical atom present in unit cell
 vso_realspace = zero
 do iatom=1,dtset%natom
!  atom type will be dtset%typat(iatom)
   
!  for each point on grid
   do ir3=1,dtset%ngfft(3)
     do ir2=1,dtset%ngfft(2)
       do ir1=1,dtset%ngfft(1)
         irealsp = ir1 + (ir2-1)*dtset%ngfft(1) + (ir3-1)*dtset%ngfft(2)*dtset%ngfft(1)

!        relative position from atom to point
         relative_position = position_op(:,ir1,ir2,ir3) - xcart(:,iatom)
         x=relative_position(1)
         y=relative_position(2)
         z=relative_position(3)

!        calculate norm^2
         norm_rel_pos = relative_position(1)**2+relative_position(2)**2+relative_position(3)**2

!        if norm^2 is too large, skip this point
         if (norm_rel_pos > radial_cutoff*radial_cutoff) cycle

!        calculate norm
         splint_x(1) = sqrt(norm_rel_pos)

!        spline interpolate vso only depends on position (through pos - atomic position)
         call splint (nradgrid,radial_grid,vso_radial(:,dtset%typat(iatom)),&
&         vso_radial_pp(:,dtset%typat(iatom)),1,splint_x,splint_y)
         vso_interpol=splint_y(1)

!        multiply by vectorial spin factor (S x r)
!        NOTE: this r is taken relative to atom center. It could be that the r operator should
!        applied in an absolute way wrt the origin...
!        
!        Is this correct: accumulated sum over atoms ?
         vso_realspace(1,irealsp,:,:,1) = vso_realspace(1,irealsp,:,:,1) + &
&         vso_interpol * reshape((/y,   zero,zero,-y/),(/2,2/))
         vso_realspace(2,irealsp,:,:,1) = vso_realspace(2,irealsp,:,:,1) + &
&         vso_interpol * reshape((/zero,z,  -z,    zero/),(/2,2/))

         vso_realspace(1,irealsp,:,:,2) = vso_realspace(1,irealsp,:,:,2) + &
&         vso_interpol * reshape((/-x,  z,   z,   x/),(/2,2/))
         vso_realspace(2,irealsp,:,:,2) = vso_realspace(2,irealsp,:,:,2) + &
&         vso_interpol * reshape((/zero,zero,zero,zero/),(/2,2/))

         vso_realspace(1,irealsp,:,:,3) = vso_realspace(1,irealsp,:,:,3) + &
&         vso_interpol * reshape((/zero,-y, -y,   zero/),(/2,2/))
         vso_realspace(2,irealsp,:,:,3) = vso_realspace(2,irealsp,:,:,3) + &
&         vso_interpol * reshape((/zero,-x, -x,    zero/),(/2,2/))


       end do  ! ir3
     end do  ! ir2
   end do  ! ir1
 end do ! iatom

 ABI_DEALLOCATE(prefact_ijl)
 ABI_DEALLOCATE(vso_radial)
 ABI_DEALLOCATE(vso_radial_pp)
 ABI_DEALLOCATE(radial_grid)

end subroutine vso_realspace_local
!!***
