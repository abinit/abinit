!!****m* ABINIT/m_vcoul/m_cutoff_sphere
!! NAME
!!  m_cutoff_sphere
!!
!! FUNCTION
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


module m_cutoff_sphere

 use defs_basis
 use m_abicore
 use m_geometry,        only : normv

 implicit none

 private

 public ::  cutoff_sphere
!!***

CONTAINS  !========================================================================================
!!***

!!****f* m_cutoff/cutoff_sphere
!! NAME
!! cutoff_sphere
!!
!! FUNCTION
!!  Calculate the Fourier transform of the Coulomb interaction with a spherical cutoff in real-space.
!!
!!   $ v_{cut}(G)= \frac{4\pi}{|q+G|^2} [ 1-cos(|q+G|*R_cut) ] $  (1)
!!
!!  For |q|<small and G=0 we use 2pi.R_cut^2, namely we consider the limit q-->0 of Eq. (1)
!!
!! INPUTS
!!  qpt(3)=q-point where the cutoff Coulomb is required.
!!  ngvec=Number of G vectors
!!  gvec(3,ngvec)=G vectors in reduced coordinates.
!!  gmet(3,3)=Metric in reciprocal space.
!!  rcut=Cutoff radius of the sphere.
!!
!! OUTPUT
!!  vc_cut(ngvec)=Fourier components of the effective Coulomb interaction.
!!
!! SOURCE

subroutine cutoff_sphere(qpt, ngvec, gvec, gmet, rcut, vc_cut)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngvec
 real(dp),intent(in) :: rcut
!arrays
 integer,intent(in) :: gvec(3,ngvec)
 real(dp),intent(in) :: gmet(3,3),qpt(3)
 real(dp),intent(out) :: vc_cut(ngvec)

!Local variables-------------------------------
!scalars
 integer :: ig,igs
 real(dp) :: qpg
 logical :: ltest

!************************************************************************

 !------------------------------------------------------------
 ! Code modification by Bruno Rousseau, Montreal, 06/11/2013
 !------------------------------------------------------------
 !
 ! In order for the code below to run in parallel using MPI
 ! within the gwls code (by Laflamme-Jansen, Cote and Rousseau)
 ! the ABI_CHECK below must be removed.
 !
 ! The ABI_CHECK below will fail if G-vectors are distributed on
 ! many processors, as only the master process has G=(0,0,0).
 !
 ! The test does not seem necessary; IF a process has G=(0,0,0)
 ! (namely, the master process), then it will be the first G vector.
 ! If a process does not have G=(0,0,0) (ie, all the other processes),
 ! then they don't need to worry about the G-> 0 limit.
 !

 ltest = ALL(gvec(:,1) == 0)
 !ABI_CHECK(ltest,'The first G vector should be Gamma')

 igs=1
 if (ltest .and. normv(qpt, gmet, 'G') < tol4) then ! For small q and G=0, use the limit q-->0.
   vc_cut(1)=two_pi*rcut**2
   igs=2
 end if
 do ig=igs,ngvec
   qpg = normv(qpt(:) + gvec(:,ig), gmet, 'G')
   vc_cut(ig) = four_pi* (one-COS(rcut*qpg)) / qpg**2
 end do

end subroutine cutoff_sphere
!!***

end module m_cutoff_sphere
!!***
