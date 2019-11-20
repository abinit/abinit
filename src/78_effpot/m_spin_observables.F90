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
  use m_spin_potential, only: spin_potential_t
  use m_multibinit_cell, only: mbsupercell_t
  use m_multibinit_dataset, only: multibinit_dtset_type

  implicit none

  private
  !!***

  !-----------------------------------------------------------------------
  !> @brief spin_observale_t: observables for spin dynamics.
  !-----------------------------------------------------------------------
  type, public :: spin_observable_t
     ! Switches for what kind of obs should be calculated
     logical :: calc_thermo_obs    ! theromostatic obs: susceptibility, specific heat...
     logical :: calc_correlation_obs ! correlation function related
     logical :: calc_traj_obs       ! trajectory related (winding number, etc)

     integer :: nspin, nsublatt, ntime, nscell
     ! nspin: number of spin
     ! nsublatt: number of sublattice (currently each spin in primitive cell is one sublattice)
     ! ntime: number of time steps
     ! nscell: number of cell in supercell

     real(dp) :: temperature
     integer, allocatable :: isublatt(:), nspin_sub(:)
     ! isublatt: index of sublattice for each spin
     ! nspin_sub: 

     real(dp) :: energy
     real(dp), allocatable :: S(:,:), Snorm(:)
     real(dp), allocatable ::  Ms_coeff(:),  Mst_sub(:, :), Mst_sub_norm(:)
     ! Ms_coeff: coefficient to calcualte staggered Mst.
     ! Staggerd means a phase factor is multiplied to each spin.
     ! Ms_coeff is that phase factor. Mst_coeff = e ^{iq R}

     ! Mst_sub: M staggered for sublattice: \sum_(i in sublattice)\S_i phase_i
     ! Mst_sub : norm of Mst_sub

     real(dp) ::  M_total(3), Mst_total(3), M_total_norm,  Mst_norm_total, Snorm_total
     ! M_total: M total    \sum M_i where M_i =|S_i|
     ! Mst_total: staggerd M total  |sum M_I| 
     ! Mst_norm : ||Mst_total||
     real(dp), allocatable :: Avg_Mst_sub_norm(:)
     real(dp) :: Avg_Mst_norm_total

     real(dp) :: binderU4, chi, Cv
     ! binderU4: binder U4
     ! chi: susceptibility
     ! Cv: specific heat

     ! variables for calculate Cv
     real(dp) :: avg_E_t   ! average energy E over time
     real(dp) :: avg_E2_t  ! average E^2 over time
     real(dp) :: avg_m_t   ! average of sum(Mst_norm_total) over t
     real(dp) :: avg_m2_t  ! average of sum(Mst_norm_total**2) over t
     real(dp) :: avg_m4_t  !average of sum(Mst_sub_total_norm**4) over t
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: reset
     procedure :: get_staggered_M
     procedure :: get_thermo_obs
     procedure :: get_correlation_obs
     procedure :: get_traj_obs
     procedure :: get_observables

  end type spin_observable_t

contains


  subroutine initialize(self, supercell , params)

    class(spin_observable_t) :: self
    type(mbsupercell_t) :: supercell
    type(multibinit_dtset_type) :: params
    integer i
    complex(dp) :: i2pi = (0.0, two_pi)

    self%calc_thermo_obs=  (params%spin_calc_thermo_obs ==1)
    self%calc_traj_obs= (params%spin_calc_traj_obs ==1)
    self%calc_correlation_obs=(params%spin_calc_correlation_obs ==1)

    self%nspin=supercell%spin%nspin
    self%nsublatt=maxval(supercell%spin%ispin_prim)


    ABI_ALLOCATE(self%S, (3, self%nspin))
    ABI_ALLOCATE(self%Snorm, (self%nspin))

    ABI_ALLOCATE(self%isublatt,(self%nspin) )
    self%isublatt(:)=supercell%spin%ispin_prim(:)

    ABI_ALLOCATE(self%nspin_sub, (self%nsublatt))
    self%nspin_sub(:)=0
    do i =1, self%nspin
       self%nspin_sub(self%isublatt(i)) = self%nspin_sub(self%isublatt(i)) + 1
    end do

    ABI_ALLOCATE(self%Ms_coeff,(self%nspin) )
    ABI_ALLOCATE(self%Mst_sub,(3, self%nsublatt) )
    ABI_ALLOCATE(self%Mst_sub_norm, (self%nsublatt))

    ABI_ALLOCATE(self%Avg_Mst_sub_norm, (self%nsublatt))

    do i =1, self%nspin
       self%Ms_coeff(i) = real(exp(i2pi * dot_product(params%spin_qpoint, supercell%spin%Rvec(:,i))))
    end do

    call reset(self, params)
  end subroutine initialize

  subroutine reset(self, params)
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
  end subroutine reset

  !-----------------------------------------------------------------------
  !> @brief finalize
  !-----------------------------------------------------------------------
  subroutine finalize(self)
    class(spin_observable_t) :: self
    if (allocated(self%isublatt)) then
       ABI_DEALLOCATE(self%isublatt)
    endif

    if (allocated(self%nspin_sub)) then
       ABI_DEALLOCATE(self%nspin_sub)
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

  end subroutine finalize

  !-----------------------------------------------------------------------
  !> @brief set S, Snorm, and energy (after one dynamics step)
  !> @param [in]  S
  !> @param [in] Snorm
  !> @param [in] energy
  !-----------------------------------------------------------------------
  subroutine update(self, S, Snorm, energy)
    class(spin_observable_t), intent(inout) :: self
    real(dp), intent(in):: S(3,self%nspin), Snorm(self%nspin), energy
    self%S=S
    self%Snorm=Snorm
    self%energy=energy
  end subroutine update

  !-----------------------------------------------------------------------
  !> @brief Calculate M in each sublattice
  !>  sum of S * phase factor in each sublattice
  !-----------------------------------------------------------------------
  subroutine get_staggered_M(self)

    class(spin_observable_t), intent(inout) :: self
    integer :: i, isub
    self%Mst_sub(:,:)=0.0
    self%M_total(:)=0.0
    do i = 1, self%nspin
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
  end subroutine get_staggered_M

  !-----------------------------------------------------------------------
  !> @brief calculate observable related to the topology of trajectory
  !> like the winding number
  !-----------------------------------------------------------------------
  subroutine get_traj_obs(self)
    class(spin_observable_t) :: self
    ABI_UNUSED(self%nspin)
  end subroutine get_traj_obs

  !-----------------------------------------------------------------------
  !> @brief calculate thermostatistic observables
  !> Cv, binderU4 and chi
  !-----------------------------------------------------------------------
  subroutine get_thermo_obs(self )
    class(spin_observable_t) :: self
    real(dp) :: avgm
    ABI_UNUSED(self%nspin)
    ! Cv
    self%avg_E_t = (self%avg_E_t*self%ntime + self%energy)/(self%ntime+1)
    self%avg_E2_t = (self%avg_E2_t*self%ntime + self%energy**2)/(self%ntime+1)
    if(self%temperature<1d-12) then
       self%Cv=0.0d0
    else
       self%Cv = (self%avg_E2_t-self%avg_E_t**2)/self%temperature**2
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
       self%chi = (self%avg_m2_t-self%avg_m_t**2)/self%temperature
    endif
  end subroutine get_thermo_obs


  !-----------------------------------------------------------------------
  !> @brief calculate correlation function relatated observables
  !-----------------------------------------------------------------------
  subroutine get_correlation_obs(self)
    class(spin_observable_t) :: self
    ABI_UNUSED(self%nspin)
  end subroutine get_correlation_obs

  !-----------------------------------------------------------------------
  !> @brief calculate all observables from input
  !> @param [in] S
  !> @param [in] Snorm
  !> @param [in] energy
  !-----------------------------------------------------------------------
  subroutine get_observables(self, S, Snorm, energy)

    class(spin_observable_t) :: self
    real(dp), intent(in) :: S(3,self%nspin), Snorm(self%nspin), energy
    call update(self, S, Snorm, energy)
    call get_staggered_M(self)
    if(self%calc_traj_obs) then
       call get_traj_obs(self)
    end if
    if(self%calc_thermo_obs) then
       call get_thermo_obs(self)
    end if
    if(self%calc_correlation_obs) then
       call get_correlation_obs(self)
    endif
    self%ntime=self%ntime+1

  end subroutine get_observables


end module m_spin_observables
