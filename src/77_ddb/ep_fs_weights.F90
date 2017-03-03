!{\src2tex{textfont=tt}}
!!****f* ABINIT/ep_fs_weights
!!
!! NAME
!! ep_fs_weights
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
!!   k_obj%wtk = integration weights
!!
!! TODO
!!   weights should be recalculated on-the-fly! The present implementation is not flexible!
!!
!! PARENTS
!!      elphon,get_nv_fs_en,get_nv_fs_temp
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ep_fs_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, fermie, gprimd, &
&    irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj)


 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_xmpi
 use m_tetrahedron

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ep_fs_weights'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: k_obj
 integer, intent(in) :: ep_b_min
 integer, intent(in) :: ep_b_max
 integer,intent(in) :: nband
 real(dp), intent(in) :: elphsmear
 real(dp), intent(in) :: fermie
 real(dp), intent(in) :: gprimd(3,3)
 integer, intent(in) :: kptrlatt(3,3)
 real(dp), intent(in) :: max_occ
 integer, intent(in) :: minFSband
 integer, intent(in) :: nFSband
 integer, intent(in) :: nsppol
 integer, intent(in) :: telphint

! arrays
 real(dp), intent(in) :: eigenGS(nband,k_obj%nkptirr,nsppol)
 integer, intent(in) :: irredtoGS(k_obj%nkptirr)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ikpt, ikptgs, ib1, isppol, iband
 integer :: nene, ifermi
 integer :: ierr

 real(dp) :: enemin, enemax, deltaene, rcvol
 real(dp) :: smdeltaprefactor, smdeltafactor, xx

! arrays
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp), allocatable :: tmp_eigen(:), tweight(:,:), dtweightde(:,:)

 character (len=500) :: message
 character (len=80) :: errstr

 type(t_tetrahedron) :: tetrahedra

! *************************************************************************

 write(std_out,*) 'ep_fs : nkpt ', k_obj%nkpt
 write(message, '(a)' ) '- ep_fs_weights  1  = '
 call wrtout(std_out,message,'PERS')

!===================================
!Set up integration weights for FS
!===================================

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

!  just do weights at FS
   nene = 100

!  fix small window around fermie for tetrahedron weight calculation
   deltaene = 2*elphsmear/dble(nene-1)
   ifermi = int(nene/2)
   enemin = fermie - dble(ifermi-1)*deltaene
   enemax = enemin + dble(nene-1)*deltaene

   ABI_ALLOCATE(tmp_eigen,(k_obj%nkpt))
   ABI_ALLOCATE(tweight,(k_obj%nkpt,nene))
   ABI_ALLOCATE(dtweightde,(k_obj%nkpt,nene))

   do iband = 1,nFSband
!    for each spin pol
     do isppol=1,nsppol
!      For this band get its contribution
       tmp_eigen(:) = zero
       do ikpt=1,k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         tmp_eigen(ikpt) = eigenGS(minFSband+iband-1,ikptgs,isppol)
       end do
!      calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223

       call get_tetra_weight(tmp_eigen,enemin,enemax,&
&       max_occ,nene,k_obj%nkpt,tetrahedra,bcorr0,&
&       tweight,dtweightde,xmpi_comm_self)

       k_obj%wtk(iband,:,isppol) = dtweightde(:,ifermi)*k_obj%nkpt
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

   k_obj%wtk = zero
!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
!    fine grid
     do ikpt=1, k_obj%nkpt
       ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
       do ib1=1,nFSband
         xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
         if (abs(xx) < 40._dp) then
           k_obj%wtk(ib1,ikpt,isppol) = exp(-xx*xx)*smdeltaprefactor
         end if
       end do
     end do
   end do


 else if (telphint == 2) then ! range of bands occupied

!  SPPOL eventually be able to specify bands for up and down separately
   k_obj%wtk = zero
   do ikpt=1,k_obj%nkpt
     do ib1=ep_b_min, ep_b_max
!      for the moment both spin channels same
       k_obj%wtk(ib1,ikpt,:) = max_occ
     end do
   end do
   
   write(std_out,*) ' ep_fs_weights : DOS is calculated from states in bands ',ep_b_min,' to ',ep_b_max
   
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

   k_obj%wtk = zero
!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
!    fine grid
     do ikpt=1, k_obj%nkpt
       ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
       do ib1=1,nFSband
         xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
         k_obj%wtk(ib1,ikpt,isppol) = smdeltaprefactor / (one + cosh(xx))
       end do
     end do
   end do
   

 else 
   write (message,'(a,i0)')" telphint should be between 0 and 3, found: ",telphint
   MSG_BUG(message)
 end if ! if telphint


end subroutine ep_fs_weights
!!***
