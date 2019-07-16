!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_observables
!! NAME
!! m_spin_observables
!!
!! FUNCTION
!! This module contains the subroutines to calculate the observables of spin dynamics
!!
!!
!! Datatypes:
!! spin_observable_t: store data to calculate observables
!! Subroutines:
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
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
  use m_spin_terms, only: spin_terms_t
  use m_multibinit_dataset, only: multibinit_dtset_type

  implicit none

  private
  !!***
  type, public :: spin_observable_t
     logical :: calc_thermo_obs, calc_correlation_obs, calc_traj_obs
     integer :: nspins, nsublatt, ntime, nscell
     real(dp) :: temperature
     integer, allocatable :: isublatt(:), nspins_sub(:)
     real(dp) :: energy
     ! Ms_coeff: coefficient to calcualte staggered Mst
     ! Mst_sub: M staggered for sublattice
     ! Mst_sub : norm of Mst_sub
     real(dp), allocatable :: S(:,:), Snorm(:)
     real(dp), allocatable ::  Ms_coeff(:),  Mst_sub(:, :), Mst_sub_norm(:)
     ! M_total: M total
     ! M_norm: ||M_total||
     ! Mst_total: staggerd M total
     ! Mst_norm : ||Mst_total||
     real(dp) ::  M_total(3), Mst_total(3), M_total_norm,  Mst_norm_total, Snorm_total
     real(dp), allocatable :: Avg_Mst_sub_norm(:)
     real(dp) :: Avg_Mst_norm_total
     real(dp) :: binderU4, chi, Cv

     ! variables for calculate Cv
     real(dp) :: avg_E_t   ! average energy E over time
     real(dp) :: avg_E2_t  ! average E^2 over time
     real(dp) :: avg_m_t   ! average of sum(Mst_norm_total) over t
     real(dp) :: avg_m2_t  ! average of sum(Mst_norm_total**2) over t
     real(dp) :: avg_m4_t  !average of sum(Mst_sub_total_norm**4) over t

  end type spin_observable_t

  public :: ob_initialize
  public :: ob_finalize
  public :: ob_reset
  public :: ob_calc_staggered_M
  public :: ob_calc_thermo_obs
  public :: ob_calc_correlation_obs
  public :: ob_calc_traj_obs
  public :: ob_calc_observables
contains


  subroutine ob_initialize(self, supercell, params)

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


    ABI_ALLOCATE(self%S, (3, self%nspins))
    ABI_ALLOCATE(self%Snorm, (self%nspins))

    ABI_ALLOCATE(self%isublatt,(self%nspins) )
    self%isublatt(:)=supercell%ispin_prim(:)

    ABI_ALLOCATE(self%nspins_sub, (self%nsublatt))
    self%nspins_sub(:)=0
    do i =1, self%nspins
       self%nspins_sub(self%isublatt(i)) = self%nspins_sub(self%isublatt(i)) + 1
    end do

    ABI_ALLOCATE(self%Ms_coeff,(self%nspins) )
    ABI_ALLOCATE(self%Mst_sub,(3, self%nsublatt) )
    ABI_ALLOCATE(self%Mst_sub_norm, (self%nsublatt))

    ABI_ALLOCATE(self%Avg_Mst_sub_norm, (self%nsublatt))

    do i =1, self%nspins
       self%Ms_coeff(i) = real(exp(i2pi * dot_product(params%spin_qpoint, supercell%Rvec(:,i))))
    end do

    call ob_reset(self, params)
  end subroutine ob_initialize

  subroutine ob_reset(self, params)
    ! set values to zeros.
    class(spin_observable_t), intent(inout) :: self
    class(multibinit_dtset_type), optional, intent(in) :: params
    if (present(params)) then
       self%temperature=params%spin_temperature
       self%nscell = product(params%ncell)
    end if
    self%ntime=0
    self%Cv=0.0
    self%binderU4=0.0
    self%M_total(:) =0.0
    self%M_total_norm=0.0
    self%Mst_norm_total=0.0
    self%Mst_sub(:,:)=0.0
    self%Mst_sub_norm(:)=0.0
    self%Avg_Mst_sub_norm(:) =0.0
    self%Avg_Mst_norm_total = 0.0
    self%avg_e_t=0.0
    self%avg_e2_t=0.0

    self%avg_m_t=0.0
    self%avg_m2_t=0.0
    self%avg_m4_t=0.0
  end subroutine ob_reset

  subroutine ob_finalize(self)
    class(spin_observable_t) :: self
    if (allocated(self%isublatt)) then
       ABI_DEALLOCATE(self%isublatt)
    endif

    if (allocated(self%nspins_sub)) then
       ABI_DEALLOCATE(self%nspins_sub)
    endif

    if(allocated(self%S)) then
       ABI_DEALLOCATE(self%S)
    endif

    if(allocated(self%Snorm)) then
       ABI_DEALLOCATE(self%Snorm)
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
    if (allocated(self%Avg_Mst_sub_norm)) then
       ABI_DEALLOCATE(self%Avg_Mst_sub_norm)
    endif

  end subroutine ob_finalize

  subroutine ob_update(self, S, Snorm, energy)
    class(spin_observable_t), intent(inout) :: self
    real(dp), intent(in):: S(3,self%nspins), Snorm(self%nspins), energy
    self%S=S
    self%Snorm=Snorm
    self%energy=energy
  end subroutine ob_update

  subroutine ob_calc_staggered_M(self)

    class(spin_observable_t), intent(inout) :: self
    integer :: i, isub
    self%Mst_sub(:,:)=0.0
    self%M_total(:)=0.0
    do i = 1, self%nspins
       isub=self%isublatt(i)
       self%Mst_sub(:, isub) = self%Mst_sub(:, isub) +  self%S(:, i)* self%Ms_coeff(i) * self%Snorm(i)
       self%M_total(:) = self%M_total + self%S(:, i)*self%Snorm(i)
    end do

    self%Mst_sub_norm(:) =sqrt(sum(self%Mst_sub**2, dim=1))/self%nscell
    self%Mst_norm_total= sum(self%Mst_sub_norm(:))
    !self%M_total_norm = sqrt(sum(self%M_total**2))/self%nscell
    self%Snorm_total = sum(self%Snorm)/self%nscell

    self%avg_Mst_sub_norm(:)=(self%avg_Mst_sub_norm(:)*self%ntime + self%Mst_sub_norm(:))/(self%ntime+1)
    self%avg_Mst_norm_total=(self%avg_Mst_norm_total*self%ntime + self%Mst_norm_total)/(self%ntime+1)
  end subroutine ob_calc_staggered_M

  subroutine ob_calc_traj_obs(self)
    class(spin_observable_t) :: self
    ABI_UNUSED(self%avg_E_t)
  end subroutine ob_calc_traj_obs

  subroutine ob_calc_thermo_obs(self)
    class(spin_observable_t) :: self
    real(dp) :: avgm
    ! Cv
    self%avg_E_t = (self%avg_E_t*self%ntime + self%energy)/(self%ntime+1)
    self%avg_E2_t = (self%avg_E2_t*self%ntime + self%energy**2)/(self%ntime+1)
    if(self%temperature<1d-12) then
       self%Cv=0.0d0
    else
       self%Cv = (self%avg_E2_t-self%avg_E_t**2)/self%temperature**2/kb_SI
    end if


    !
    avgm=self%Mst_norm_total

    self%avg_m_t =  (self%avg_m_t*self%ntime + avgm)/(self%ntime+1)
    self%avg_m2_t = (self%avg_m2_t*self%ntime + avgm**2)/(self%ntime+1)
    self%avg_m4_t = (self%avg_m4_t*self%ntime + avgm**4)/(self%ntime+1)

    self%binderU4 = 1.0-self%avg_m4_t/self%avg_m2_t**2/3.0

    if(self%temperature<1d-12) then
       self%chi=(self%avg_m2_t-self%avg_m_t**2)
    else
       self%chi = (self%avg_m2_t-self%avg_m_t**2)/self%temperature/kb_SI
    endif

  end subroutine ob_calc_thermo_obs

  subroutine ob_calc_correlation_obs(self)

    class(spin_observable_t) :: self
    ABI_UNUSED(self%ntime)
  end subroutine ob_calc_correlation_obs

  subroutine ob_calc_observables(self, S, Snorm, energy)

    class(spin_observable_t) :: self
    real(dp), intent(in) :: S(3,self%nspins), Snorm(self%nspins), energy
    call ob_update(self, S, Snorm, energy)
    call ob_calc_staggered_M(self)
    if(self%calc_traj_obs) then
       call ob_calc_traj_obs(self)
    end if
    if(self%calc_thermo_obs) then
       call ob_calc_thermo_obs(self)
    end if
    if(self%calc_correlation_obs) then
       call ob_calc_correlation_obs(self)
    endif
    self%ntime=self%ntime+1

  end subroutine ob_calc_observables


end module m_spin_observables
