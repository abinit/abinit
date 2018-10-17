!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_model
!! NAME
!! m_spin_observables
!!
!! FUNCTION
!! This module contains the subroutines to calculate the observables of spin dynamics
!!
!!
!! Datatypes:
!!
!! Subroutines:
!!
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_spin_observables

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  implicit none
  private
  !!***
  type, public :: spin_observable_t
     logical :: calc_thermo_obs, calc_correlation_obs, calc_traj_obs
     integer :: nspins, nsublatt, ntime
     integer, allocatable :: isublatt(:), nspins_sub(:)

     ! Ms_coeff: coefficient to calcualte staggered Mst
     ! Mst_sub: M staggered for sublattice
     ! Mst_sub : norm of Mst_sub
     real(dp), allocatable ::  Ms_coeff(:),  Mst_sub(:, :), Mst_sub_norm(:)
     ! M_total: M total 
     ! M_norm: ||M_total||
     ! Mst_total: staggerd M total
     ! Mst_norm : ||Mst_total||
     real(dp) ::  M_total(3), Mst_total(3), M_total_norm,  Mst_norm_total, Snorm_total
  end type spin_observable_t

  public :: ob_initialize
  public :: ob_finalize
  public :: ob_calc_staggered_M
  public :: ob_calc_thermo_obs
  public :: ob_calc_correlation_obs
  public :: ob_calc_traj_obs
  public :: ob_calc_observables
contains


  subroutine ob_initialize(self, supercell, params)

  use m_spin_terms, only: spin_terms_t
  use m_multibinit_dataset, only: multibinit_dtset_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_initialize'
!End of the abilint section

    class(spin_observable_t) :: self
    type(spin_terms_t) :: supercell
    type(multibinit_dtset_type) :: params
    integer i
    complex(dp) :: i2pi = (0.0, two_pi)

    self%calc_thermo_obs=  (params%spin_calc_thermo_obs ==1)
    self%calc_traj_obs= (params%spin_calc_traj_obs ==1)
    self%calc_correlation_obs=(params%spin_calc_correlation_obs ==1)

    self%nspins=supercell%nspins
    self%nsublatt=maxval(supercell%ispin_prim)

    self%ntime=0

    ABI_ALLOCATE(self%isublatt,(self%nspins) )
    self%isublatt(:)=supercell%ispin_prim(:)

    ABI_ALLOCATE(self%nspins_sub, (self%nsublatt))
    self%nspins_sub(:)=0
    do i =1, self%nspins
       self%nspins_sub(self%isublatt(i)) = self%nspins_sub(self%isublatt(i)) + 1
    end do

    ABI_ALLOCATE(self%Ms_coeff,(self%nspins) )

    self%M_total(:) =0.0
    self%M_total_norm=0.0

    self%Mst_norm_total=0.0

    ABI_ALLOCATE(self%Mst_sub,(3, self%nsublatt) )
    self%Mst_sub(:,:)=0.0
    ABI_ALLOCATE(self%Mst_sub_norm, (self%nsublatt))
    self%Mst_sub_norm=0.0

    do i =1, self%nspins
       self%Ms_coeff(i) = real(exp(i2pi * dot_product(params%spin_qpoint, supercell%Rvec(:,i))))
    end do
  end subroutine ob_initialize

  subroutine ob_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_finalize'
!End of the abilint section

    class(spin_observable_t) :: self
    if (allocated(self%isublatt)) then
       ABI_DEALLOCATE(self%isublatt)
    endif

    if (allocated(self%Ms_coeff)) then
       ABI_DEALLOCATE(self%Ms_coeff)
    endif
    if (allocated(self%Mst_sub)) then
       ABI_DEALLOCATE(self%Mst_sub)
    endif
    if (allocated(self%Mst_sub_norm)) then
       ABI_DEALLOCATE(self%Mst_sub_norm)
    endif

  end subroutine ob_finalize

  subroutine ob_calc_staggered_M(self, S, Snorm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_calc_staggered_M'
!End of the abilint section

    class(spin_observable_t), intent(inout) :: self
    real(dp), intent(in):: S(3,self%nspins), Snorm(self%nspins)
    integer :: i, isub
    self%Mst_sub(:,:)=0.0
    self%M_total(:)=0.0
    do i = 1, self%nspins
       isub=self%isublatt(i)
       self%Mst_sub(:, isub) = self%Mst_sub(:, isub) +  S(:, i)* self%Ms_coeff(i) * Snorm(i)
       self%M_total(:) = self%M_total + S(:, i)*Snorm(i)
    end do

    self%Mst_sub_norm =sqrt(sum(self%Mst_sub**2, dim=1))
    self%Mst_norm_total= sum(self%Mst_sub_norm)
    self%M_total_norm = sqrt(sum(self%M_total**2))
    self%Snorm_total = sum(Snorm)

  end subroutine ob_calc_staggered_M

  subroutine ob_calc_traj_obs(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_calc_traj_obs'
!End of the abilint section

      class(spin_observable_t) :: self
  end subroutine ob_calc_traj_obs

  subroutine ob_calc_thermo_obs(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_calc_thermo_obs'
!End of the abilint section

    class(spin_observable_t) :: self
  end subroutine ob_calc_thermo_obs

  subroutine ob_calc_correlation_obs(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_calc_correlation_obs'
!End of the abilint section

    class(spin_observable_t) :: self
  end subroutine ob_calc_correlation_obs

  subroutine ob_calc_observables(self, S, Snorm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ob_calc_observables'
!End of the abilint section

    class(spin_observable_t) :: self
    real(dp), intent(in) :: S(3,self%nspins), Snorm(self%nspins)
    call ob_calc_staggered_M(self, S, Snorm)
  end subroutine ob_calc_observables

end module m_spin_observables
