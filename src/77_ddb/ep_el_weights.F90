!{\src2tex{textfont=tt}}
!!****f* ABINIT/ep_el_weights
!!
!! NAME
!! ep_el_weights
!!
!! FUNCTION
!! This routine calculates the Fermi Surface integration weights
!!  for the electron phonon routines, by different methods
!!    1) Gaussian smearing
!!    2) Tetrahedron method
!!    3) Window in bands for all k-points
!!    4) Fermi Dirac smearing, follows gaussian with a different smearing function
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   ep_b_min = minimal band to include in FS window integration
!!   ep_b_max = maximal band to include in FS window integration
!!   eigenGS = Ground State eigenvalues
!!   elphsmear = smearing width for Gaussian method
!!   fermie = Fermi level
!!   gprimd = Reciprocal lattice vectors (dimensionful)
!!   irredtoGS = mapping of elph k-points to ground state grid
!!   kptrlatt = k-point grid vectors (if divided by determinant of present matrix)
!!   max_occ = maximal occupancy for a band
!!   minFSband = minimal band included for Fermi Surface integration in Gaussian and Tetrahedron cases
!!   nFSband = number of bands in FS integration
!!   nsppol = number of spin polarizations
!!   telphint = option for FS integration:
!!      0 Tetrahedron method
!!      1 Gaussian smearing
!!      2 Window in bands for all k-points
!!      3 Fermi Dirac smearing
!!   k_obj%nkpt = number of FS k-points
!!   k_obj%kpt = FS k-points
!!   k_obj%full2irr = mapping of FS k-points from full grid to irred points
!!   k_obj%full2full = mapping of FS k-points in full grid under symops
!!
!! OUTPUT
!!
!! TODO
!!   weights should be recalculated on-the-fly! The present implementation is not flexible!
!!
!! PARENTS
!!      get_nv_fs_en,get_tau_k
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ep_el_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, enemin, enemax, nene, gprimd, &
&    irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj, tmp_wtk)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors
 use m_tetrahedron
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ep_el_weights'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(in) :: k_obj
 integer, intent(in) :: ep_b_min
 integer, intent(in) :: ep_b_max
 integer,intent(in) :: nband,nene
 real(dp), intent(in) :: elphsmear
 real(dp), intent(in) :: enemin,enemax
 real(dp), intent(in) :: gprimd(3,3)
 integer, intent(in) :: kptrlatt(3,3)
 real(dp), intent(in) :: max_occ
 integer, intent(in) :: minFSband
 integer, intent(in) :: nFSband
 integer, intent(in) :: nsppol
 integer, intent(in) :: telphint

! arrays
 real(dp), intent(in) :: eigenGS(nband,k_obj%nkptirr,nsppol)
 real(dp), intent(out) :: tmp_wtk(nFSband,k_obj%nkpt,nsppol,nene)
 integer, intent(in) :: irredtoGS(k_obj%nkptirr)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ikpt, ikptgs, ib1, iband
 integer :: ierr, ie, isppol
 real(dp) :: deltaene, rcvol, fermie
 real(dp) :: smdeltaprefactor, smdeltafactor, xx

! arrays
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp), allocatable :: tmp_eigen(:), tweight(:,:), dtweightde(:,:)
 character (len=500) :: message
 character (len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra

! *************************************************************************

 ! Initialize tmp_wtk with zeros
 tmp_wtk = zero

 !write(std_out,*) 'ep_el : nkpt ', k_obj%nkpt
!===================================
!Set up integration weights for FS
!===================================
 deltaene = (enemax-enemin)/dble(nene-1)

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

!  fix small window around fermie for tetrahedron weight calculation
   deltaene = (enemax-enemin)/dble(nene-1)

   ABI_ALLOCATE(tmp_eigen,(k_obj%nkpt))
   ABI_ALLOCATE(tweight,(k_obj%nkpt,nene))
   ABI_ALLOCATE(dtweightde,(k_obj%nkpt,nene))

   do iband = 1,nFSband
!    for each spin pol
     do isppol=1,nsppol
!    For this band get its contribution
       tmp_eigen(:) = zero
       do ikpt=1,k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         tmp_eigen(ikpt) = eigenGS(minFSband+iband-1,ikptgs,isppol)
       end do
!      calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223
       call get_tetra_weight(tmp_eigen,enemin,enemax,&
&       max_occ,nene,k_obj%nkpt,tetrahedra,bcorr0,&
&       tweight,dtweightde,xmpi_comm_self)

       tmp_wtk(iband,:,isppol,:) = dtweightde(:,:)*k_obj%nkpt
     end do
   end do
   ABI_DEALLOCATE(tmp_eigen)
   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)

   call destroy_tetra(tetrahedra)

 else if (telphint == 1) then

!  ==============================================================
!  Gaussian or integration:
!  Each kpt contributes a gaussian of integrated weight 1 
!  for each band. The gaussian being centered at the Fermi level.
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  MJV 18/5/2008 does smdeltaprefactor need to contain max_occ?

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = max_occ*sqrt(piinv)/elphsmear
   smdeltafactor = one/elphsmear

!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
     fermie = enemin
     do ie = 1, nene
       fermie = fermie + deltaene
!      fine grid
       do ikpt=1, k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         do ib1=1,nFSband
           xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
           if (abs(xx) < 40._dp) then
             tmp_wtk(ib1,ikpt,isppol,ie) = exp(-xx*xx)*smdeltaprefactor
           end if
         end do
       end do
     end do
   end do
   

 else if (telphint == 2) then ! range of bands occupied

!  SPPOL eventually be able to specify bands for up and down separately
   fermie = enemin
   do ie = 1, nene
     fermie = fermie + deltaene
     do ikpt=1,k_obj%nkpt
       do ib1=ep_b_min, ep_b_max
!        for the moment both spin channels same
         tmp_wtk(ib1,ikpt,:,ie) = max_occ
       end do
     end do
   end do
   
   write(std_out,*) ' ep_el_weights : DOS is calculated from states in bands ',ep_b_min,' to ',ep_b_max
   
 else if (telphint == 3) then

!  ==============================================================
!  Fermi Dirac integration:
!  Each kpt contributes a Fermi Dirac smearing function of integrated weight 1 
!  for each band. The function being centered at the Fermi level.
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  MJV 18/5/2008 does smdeltaprefactor need to contain max_occ?

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = half*max_occ/elphsmear
   smdeltafactor = one/elphsmear

!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
     fermie = enemin
     do ie = 1, nene
       fermie = fermie + deltaene
!      fine grid
       do ikpt=1, k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         do ib1=1,nFSband
           xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
           tmp_wtk(ib1,ikpt,isppol,ie) = smdeltaprefactor / (one + cosh(xx))
         end do
       end do
     end do
   end do
   

 else 
   write (message,'(a,i0)')" telphint should be between 0 and 3, found: ",telphint
   MSG_BUG(message)
 end if ! if telphint


end subroutine ep_el_weights
!!***
