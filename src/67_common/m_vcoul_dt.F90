!!****m* ABINIT/m_vcoul_dt
!! NAME
!!  m_vcoul_dt
!!
!! FUNCTION
!!  This module contains the vcoul datatype 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (MG, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

MODULE m_vcoul_dt

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 public
!!***

!!****t* m_vcoul/vcoul_t
!! NAME
!!  vcoul_t

 type,public :: vcoul_t

  ! TODO: Remove it
  integer  :: nfft
  ! Number of points in FFT grid

  integer  :: ng
   ! Number of G-vectors

  integer  :: nqibz
   ! Number of irreducible q-points

  integer  :: nqlwl
   ! Number of small q-points around Gamma

  real(dp) :: alpha(3)
   ! Lenght of the finite surface

  real(dp) :: rcut
   ! Cutoff radius

  real(dp) :: i_sz
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  real(dp) :: i_sz_resid
   ! Residual difference between the i_sz in the sigma self-energy for exchange,
   ! and the i_sz already present in the generalized Kohn-Sham eigenenergies
   ! Initialized to the same value as i_sz

  real(dp) :: hcyl
   ! Length of the finite cylinder along the periodic dimension

  real(dp) :: ucvol
    ! Volume of the unit cell

  character(len=50) :: mode
   ! String defining the cutoff mode, possible values are: sphere,cylinder,surface,crystal

  integer :: pdir(3)
   ! 1 if the system is periodic along this direction

  ! TODO: Remove it
  integer :: ngfft(18)
    ! Information on the FFT grid

  real(dp) :: boxcenter(3)
   ! 1 if the point in inside the cutoff region 0 otherwise
   ! Reduced coordinates of the center of the box (input variable)

  real(dp) :: vcutgeo(3)
    ! For each reduced direction gives the length of the finite system
    ! 0 if the system is infinite along that particular direction
    ! negative value to indicate that a finite size has to be used

  real(dp) :: rprimd(3,3)
    ! Lattice vectors in real space.

  real(dp),allocatable :: qibz(:,:)
   ! qibz(3,nqibz)
   ! q-points in the IBZ.

  real(dp),allocatable :: qlwl(:,:)
   ! qibz(3,nqlwl)
   ! q-points for the treatment of the Coulomb singularity.

  complex(gwpc),allocatable :: vc_sqrt(:,:)
    ! vc_sqrt(ng,nqibz)
    ! Square root of the Coulomb interaction in reciprocal space.
    ! A cut might be applied.

  complex(gwpc),allocatable :: vcqlwl_sqrt(:,:)
    ! vcqs_sqrt(ng,nqlwl)
    ! Square root of the Coulomb term calculated for small q-points

  complex(gwpc),allocatable :: vc_sqrt_resid(:,:)
    ! vc_sqrt_resid(ng,nqibz)
    ! Square root of the residual difference between the Coulomb interaction in the sigma self-energy for exchange,
    ! and the Coulomb interaction already present in the generalized Kohn-Sham eigenenergies (when they come from an hybrid)
    ! Given in reciprocal space. At the call to vcoul_init, it is simply initialized at the value of vc_sqrt(:,:),
    ! and only later modified.
    ! A cut might be applied.

 end type vcoul_t
!!***

contains

end module m_vcoul_dt
!!***
