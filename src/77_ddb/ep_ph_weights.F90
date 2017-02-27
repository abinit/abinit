!{\src2tex{textfont=tt}}
!!****f* ABINIT/ep_ph_weights
!!
!! NAME
!! ep_ph_weights
!!
!! FUNCTION
!! This routine calculates the phonon integration weights
!!  for the electron phonon routines, by different methods
!!    1) Gaussian smearing
!!    0) Tetrahedron method
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MVer,BXu)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   phfrq = phonon energies
!!   elphsmear = smearing width for Gaussian method
!!   omega = input phonon energy
!!   gprimd = Reciprocal lattice vectors (dimensionful)
!!   kptrlatt = k-point grid vectors (if divided by determinant of present matrix)
!!   telphint = option for FS integration:
!!      0 Tetrahedron method
!!      1 Gaussian smearing
!!   k_obj%nkpt = number of FS k-points
!!   k_obj%kpt = FS k-points
!!   k_obj%full2full = mapping of FS k-points in full grid under symops
!!
!! OUTPUT
!!   tmp_wtq = integration weights
!!
!! TODO
!!   weights should be recalculated on-the-fly! The present implementation is not flexible!
!!
!! PARENTS
!!      get_tau_k,mka2f,mka2f_tr
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ep_ph_weights(phfrq,elphsmear,omega_min,omega_max,nomega,gprimd,kptrlatt,nbranch,telphint,k_obj,tmp_wtq)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors
 use m_tetrahedron
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ep_ph_weights'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: k_obj
 integer,intent(in) :: nbranch
 real(dp), intent(in) :: elphsmear
 real(dp), intent(in) :: omega_min,omega_max
 real(dp), intent(in) :: gprimd(3,3)
 integer, intent(in) :: kptrlatt(3,3)
 integer, intent(in) :: nomega
 integer, intent(in) :: telphint

! arrays
 real(dp), intent(in) :: phfrq(nbranch,k_obj%nkpt)
 real(dp), intent(out) :: tmp_wtq(nbranch,k_obj%nkpt,nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ikpt, ib1, ibranch
 integer :: ierr, iomega
 real(dp) :: rcvol, max_occ
 real(dp) :: smdeltaprefactor, smdeltafactor, xx, gaussmaxarg
 real(dp) :: domega,omega

! arrays
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp), allocatable :: tweight(:,:), dtweightde(:,:)
 character (len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra

! *************************************************************************

 !write(std_out,*) 'ep_ph : nqpt ', k_obj%nkpt
!===================================
!Set up integration weights for FS
!===================================
 max_occ = one
 gaussmaxarg = sqrt(-log(1.d-100))
 domega = (omega_max - omega_min)/(nomega - 1)

 if (telphint == 0) then

!  =========================
!  Tetrahedron integration
!  =========================

   rlatt(:,:) = kptrlatt(:,:)
   call matr3inv(rlatt,klatt)

   call init_tetra(k_obj%full2full(1,1,:), gprimd,klatt,k_obj%kpt, k_obj%nkpt,&
&   tetrahedra, ierr, errstr)
   ABI_CHECK(ierr==0,errstr)

   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!  do all the omega points for tetrahedron weight calculation

   ABI_ALLOCATE(tweight,(k_obj%nkpt,nomega))
   ABI_ALLOCATE(dtweightde,(k_obj%nkpt,nomega))

   do ibranch = 1,nbranch
     call get_tetra_weight(phfrq(ibranch,:),omega_min,omega_max,&
&     max_occ,nomega,k_obj%nkpt,tetrahedra,bcorr0,&
&     tweight,dtweightde,xmpi_comm_self)

     tmp_wtq(ibranch,:,:) = dtweightde(:,:)*k_obj%nkpt
   end do
   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)

   call destroy_tetra(tetrahedra)

 else if (telphint == 1) then

!  ==============================================================
!  Gaussian or integration:
!  Each kpt contributes a gaussian of integrated weight 1 
!  for each branch. The gaussian being centered at the input energy
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = max_occ*sqrt(piinv)/elphsmear
   smdeltafactor = one/elphsmear

   tmp_wtq = zero
   omega = omega_min
   do iomega = 1, nomega
     omega = omega + domega
     do ikpt=1, k_obj%nkpt
       do ib1=1,nbranch
         xx = smdeltafactor*(phfrq(ib1,ikpt)-omega)
         if (abs(xx) < gaussmaxarg) then
           tmp_wtq(ib1,ikpt,iomega) = exp(-xx*xx)*smdeltaprefactor
         end if
       end do
     end do
   end do
 end if ! if telphint


end subroutine ep_ph_weights
!!***
