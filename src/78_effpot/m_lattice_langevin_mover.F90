!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_langevin_mover
!! NAME
!! m_lattice_langevin_mover
!!
!! FUNCTION
!! This module contains the langevin  (NVT) lattice mover.
!! It is a translation from the ASE (GPL licenced) Langevin mover python code to fortran.
!! The original code can be found at
!! https://gitlab.com/ase/ase/blob/master/ase/md/langevin.py
!! The method is described in
!! E. V.-Eijnden, and G. Ciccotti, Chem. Phys. Lett. 429, 310 (2006)
!!  https://doi.org/10.1016/j.cplett.2006.07.086
!!
!!
!! Datatypes:
!!
!! * lattice_langevin_mover_t: defines the lattice movers
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
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

module m_lattice_langevin_mover
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_mover, only: lattice_mover_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_random_xoroshiro128plus, only:  rng_t
  use m_hashtable_strval, only: hash_table_t
!!***

  implicit none

  private

  !-------------------------------------------------------------------!
  ! Lattice_langevin_mover_t
  ! 
  !-------------------------------------------------------------------!
  type, public, extends(lattice_mover_t) :: lattice_langevin_mover_t
     ! c1 to c5 are a constants (temperature and mass dependent)
     ! c3 c4 c5 has dimension of natom.
     real(dp) :: c1, c2
     real(dp), allocatable :: c3(:), c4(:), c5(:)
     real(dp) :: fr =1e-4 ! friction, usually 1e-2~1e-4
     ! xi and eta: random numbers of dimension (3, natom)
     real(dp), allocatable :: xi(:,:), eta(:,:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_one_step
     procedure :: update_vars
  end type lattice_langevin_mover_t

contains

  !-------------------------------------------------------------------!
  ! Initialize:
  !  read parameters
  !  allocate memory and call update_vars()
  !-------------------------------------------------------------------!
  subroutine initialize(self,params, supercell, rng)
    class(lattice_langevin_mover_t), intent(inout) :: self
    type(multibinit_dtset_type), target, intent(in):: params
    type(mbsupercell_t), target, intent(in) :: supercell
    type(rng_t), target, intent(in) :: rng
    call self%lattice_mover_t%initialize(params, supercell, rng)

    ! TODO: add friction 
    self%fr = params%latt_friction

    ABI_ALLOCATE(self%c3, (self%natom))
    ABI_ALLOCATE(self%c4, (self%natom))
    ABI_ALLOCATE(self%c5, (self%natom))
    ABI_ALLOCATE(self%xi, (3,self%natom))
    ABI_ALLOCATE(self%eta, (3,self%natom))

    call self%update_vars()
  end subroutine initialize


  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(lattice_langevin_mover_t), intent(inout) :: self
    ABI_DEALLOCATE(self%c3)
    ABI_DEALLOCATE(self%c4)
    ABI_DEALLOCATE(self%c5)
    ABI_DEALLOCATE(self%xi)
    ABI_DEALLOCATE(self%eta)
    call self%lattice_mover_t%finalize()
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! update_vars:
  !  calculate c1 to c5 from dt, masses and temperature
  !  It is called by the initialization function.
  !-------------------------------------------------------------------!
  subroutine update_vars(self)
    class(lattice_langevin_mover_t), intent(inout) :: self
    real(dp) :: dt, T, fr, sigma(self%natom)

    dt=self%dt
    T= self%temperature
    fr=self%fr
    sigma(:) = sqrt(2.0*T*fr/self%masses(:))
    self%c1 = dt / 2.0 - dt * dt * fr / 8.0
    self%c2 = dt * fr / 2.0 - dt * dt * fr * fr / 8.0
    self%c3 = sqrt(dt) * sigma / 2.0 - dt**1.5 * fr * sigma / 8.0
    self%c5 = dt**1.5 * sigma / (2.0 * sqrt(3.0))
    self%c4 = fr / 2. * self%c5
  end subroutine update_vars

  !-------------------------------------------------------------------!
  !> @brief: Run_one_step using a Langevin heat bath.
  !>   effpot: the potential (which do the calculation of E and dE/dvar)
  !> param[in]: effpot
  !> param[in]: (optional) displacement
  !> param[in]: (optional) strain
  !> param[in]: (optional) spin
  !> param[in]: (optional) lwf
  ! NOTE: No need to pass the variable already saved in the mover.
  !     e.g. For spin mover, do NOT pass the spin to it.
  !    The other variables are only required if there is coupling with
  !    the mover variable.
  !-------------------------------------------------------------------!
  subroutine run_one_step(self, effpot,displacement, strain, spin, lwf , energy_table)
    class(lattice_langevin_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    integer :: i
    character(len=40) :: key

    ! do not use displacement and strain because they are stored in the mover.
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)

    self%energy = 0.0
    self%forces(:,:) =0.0
    call effpot%calculate( displacement=self%displacement, strain=self%strain, &
         & spin=spin, lwf=lwf, force=self%forces, stress=self%stress, &
         & energy=self%energy, energy_table=energy_table)
    call self%rng%rand_normal_array(self%xi, 3*self%natom)
    call self%rng%rand_normal_array(self%eta, 3*self%natom)



    ! First half of velocity update
    do i =1, self%natom
       self%current_vcart(:,i) = self%current_vcart(:,i) + &
            & self%c1 * self%forces(:,i) / self%masses(i) - &
            & self%c2 * self%current_vcart(:,i) + &
            & self%c3(i) * self%xi(:, i) - self%c4(i) * self%eta(:,i) 

       self%displacement(:, i) = self%displacement(:, i) &
            & + self%dt * self%current_vcart(:, i) + self%c5( i) *self%eta(:,i)
    end do

    ! second half, update the velocity but not the displacement.
    self%energy=0.0
    self%forces=0.0
    call effpot%calculate( displacement=self%displacement, strain=self%strain, &
         & spin=spin, lwf=lwf, force=self%forces, stress=self%stress, &
         & energy=self%energy, energy_table=energy_table)
    do i =1, self%natom
       self%current_vcart(:,i) = self%current_vcart(:,i) + &
            &self%c1 * self%forces(:,i) / self%masses(i) - &
            & self%c2 * self%current_vcart(:,i) + &
            & self%c3(i) * self%xi(:, i) - self%c4(i) * self%eta(:,i)
    end do

    call self%get_T_and_Ek()
    if (present(energy_table)) then
      key = 'Lattice kinetic energy'
       call energy_table%put(key, self%Ek)
    end if

  end subroutine run_one_step

end module m_lattice_langevin_mover

