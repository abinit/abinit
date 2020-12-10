!!****m* ABINIT/m_spin_hist
!! NAME
!! m_spin_hist
!!
!! FUNCTION
!! This module contains definition the type spin_hist_t
!! and its related routines
!! The observables are also calculated. 
!!
!! Datatypes:
!!
!! * spin_hist_t: history record of spin orientations and amplitudes
!!
!! Subroutines:
!!
!! * init
!! * free
!! * spin_hist_t
!! * get_S
!! * findIndex
!! * set_vars
!! * set_params
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

! TODO hexu:
! sync ihist_latt when with lattice dynamics
! add average, variance, etc (should they be here?)
! structural information and some parameters are no longer 
! used here. They should be removed form this file.

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
module m_spin_hist
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  implicit none

  private
!!***

  !----------------------------------------------------------------------

  !!****t* m_spin_hist/spin_hist_t
  !! NAME
  !! spin_hist_t
  !!
  !! FUNCTION
  !! This type has several vectors, and index scalars to store
  !! a proper history of previous evaluations of forces and
  !! stresses,velocities,positions and energies
  !!
  !! It contains:
  !! * mxhist                  : Maximum size of history
  !! * ihist                   : index of history

  !! natoms : number of atoms
  !! nspin: number of magnetic atoms
  !! * acell(3)         : Acell (acell , rprimd, xred: only initial value kept if there is!!  no lattice dynamics. Other wise for each step, the corresponding lattice step number is kept)
  !! * rprimd(3,3)      : Rprimd
  !! * xred(3,natoms)    : Xred
  !! * index_spin     : the index of atom in spin model, -1 if it is not in the spin model 
  !! * heff(3,nspin,mxhist)   : effective magnetic field (cartesian)
  !! * snorm(nspin, mxhist) : magnetitude of spin.
  !! * S(3,nspin,mxhist)   : spin orientation of atoms (cartesian)
  !! * dSdt(3, nspin, mxhist) : dS/dt (cartesian)
  !! * etot(mxhist)            : Electronic total Energy
  !! * entropy(mxhist)         : Entropy
  !! * itime(mxhist)           : index of spin dynamics step.
  !! * time(mxhist)            : Time (or iteration number for GO)
  !!
  !! * has_latt (whether lattice dynamics is also present)
  !! * ihist_latt(mxhist): the corresponding lattice step. 0 if none.
  !! SOURCE

  type, public :: spin_hist_t
     ! scalars
     ! Index of the last element on all records
     integer :: ihist = 0
     integer :: ihist_prev = -1
     ! Maximun size of the historical records
     integer :: mxhist = 0

     integer :: nspin, nspin_prim
     ! whether lattice dynamics is also present
     integer, allocatable :: ihist_latt(:)
     logical :: has_latt

     ! arrays
     !  placeholders for structure-related parameters. They are not used currently. 
     integer :: natoms
     real(dp) :: acell(3)
     real(dp) :: rprimd(3,3)
     real(dp), allocatable :: xred(:, :)
     integer :: ntypat
     integer, allocatable :: typat(:)
     real(dp), allocatable :: znucl(:)
     integer, allocatable :: spin_index(:)

     ! spin
     !heff(3, nspin, mxhist)
     real(dp), allocatable :: heff(:, :, :)
     !snorm(nspin, mxhist)
     real(dp), allocatable :: snorm(:, :)

     !S(3, nspin, mxhist)
     real(dp), allocatable :: S(:, :, :)
     !dSdt(3, nspin, mxhist)
     ! TODO hexu: is it useful?
     real(dp), allocatable :: dSdt(:, :, :)

     ! etot(mxhist)
     real(dp), allocatable :: etot(:)
     real(dp), allocatable :: entropy(:)
     real(dp), allocatable :: time(:)
     integer, allocatable :: itime(:)

     ! spin_nctime: interval of step for writing to netcdf hist file.
     integer :: spin_nctime
     real(dp) :: spin_temperature

     ! observables
     integer:: calc_thermo_obs, calc_traj_obs, calc_correlation_obs

     real(dp), allocatable :: ms_sub(:,:)   ! staggered M. 
     real(dp), allocatable :: Cv(:) ! specfic heat
     real(dp), allocatable :: binderU4_sub(:,:), binderU4(:)
     real(dp), allocatable :: chi_sub(:, :), chi(:) ! magnetic susceptibility
     real(dp), allocatable :: rcorr(:,:)
     real(dp), allocatable :: sp_corr_func(:,:,:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: reset
     procedure :: set_vars
     procedure :: get_S => get_S
     procedure :: findIndex => findIndex
     procedure :: set_params => set_params
     procedure :: inc1
  end type spin_hist_t
  !!***

  !public :: spinhist2var
  !public :: var2spinhist
  !public :: write_sd_hist
  !public :: read_md_hist
  !public :: get_dims_spinhist

contains


!!****f* m_spin_hist/initialize
!!
!! NAME
!! initialize
!!
!! FUNCTION
!! initialize spin hist
!!
!! INPUTS
!! nspin = number of magnetic atoms
!! mxhist = maximum number of hist steps
!! has_latt = whether spin dynamics in with lattice dynamics
!!
!! OUTPUT
!! hist <type(spin_hist_t)()> = spin hist type
!! PARENTS
!!      m_spin_hist
!!
!! CHILDREN
!!
!! SOURCE

  subroutine initialize(self, nspin, mxhist, has_latt)

    implicit none
    class(spin_hist_t), intent(inout) :: self 
    integer, intent(in) :: nspin, mxhist
    logical, intent(in) :: has_latt
    !integer, optional,  intent(in) :: calc_traj_obs, calc_thermo_obs, calc_correlation_obs

    self%nspin=nspin
    self%ntypat=0
    self%ihist=1
    self%ihist_prev=0
    self%mxhist=mxhist
    self%natoms=0
    self%has_latt=has_latt

    ABI_ALLOCATE(self%heff, (3, nspin, mxhist))
    ABI_ALLOCATE(self%snorm, (nspin, mxhist))
    ABI_ALLOCATE(self%S, (3, nspin, mxhist))
    ABI_ALLOCATE(self%dSdt, (3, nspin, mxhist))

    ABI_ALLOCATE(self%etot, (mxhist))
    ABI_ALLOCATE(self%entropy, (mxhist))
    ABI_ALLOCATE(self%time, (mxhist))
    ABI_ALLOCATE(self%itime, (mxhist))

    ABI_ALLOCATE(self%ihist_latt, (mxhist))

    ! TODO: add observable allocation here.

    self%etot(1) =zero
    self%entropy(1) =zero
    self%time(1) =zero

    !self%acell(:)=zero
    !self%rprimd(:, :)=zero
    !self%xred(:,:) =zero
    self%heff(:,:,:)=zero
    self%S(:,:,:)=zero
    self%dSdt(:,:,:)=zero
    self%snorm(:,:)=zero
  end subroutine initialize
!!***

  subroutine reset(self, array_to_zero)

    implicit none
    class(spin_hist_t), intent(inout) :: self
    logical :: array_to_zero
    self%ntypat=0
    self%ihist=1
    self%ihist_prev=0
    self%natoms=0

    self%etot(1) =zero
    self%entropy(1) =zero
    self%time(1) =zero

    if(array_to_zero) then
       self%heff(:,:,1)=zero
       self%S(:,:,1)=zero
       self%dSdt(:,:,1)=zero
       self%snorm(:,1)=zero
       self%Cv( 1)=zero
       self%sp_corr_func(:, :, 1)=zero
    endif


  end subroutine reset

  !!****f* m_spin_hist/set_atomic_structure
  !!
  !! NAME
  !! set_atomic_structure
  !!
  !! FUNCTION
  !! 
  !! set atomic structure 
  !!
  !! INPUTS
  !! acell(3) = acell
  !! rprimd(3, 3) = 
  !! xred(3, natoms) = positions in reduced coordinates
  !! spin_index(3, natoms) = index of atom in spin hamiltonian
  !! ntypat = number of types of atoms
  !! typat(ntypat)=types of atoms
  !! znucl=z of atoms
  !!
  !! OUTPUT
  !! hist <type(spin_hist_t)()> = spin hist type
  !! PARENTS
!!
  !! CHILDREN
!!      self%inc1
!!
  !! SOURCE
  subroutine set_atomic_structure(self, acell, rprimd, xred, spin_index, ntypat,  typat, znucl)

    class(spin_hist_t), intent(inout) :: self
    real(dp), intent(in) :: acell(3), rprimd(3,3), xred(:,:), znucl(:)
    integer, intent(in):: spin_index(:), ntypat, typat(:)
    integer :: natoms
    natoms=size(typat)
    ABI_ALLOCATE(self%xred, (3, natoms))
    ABI_ALLOCATE(self%spin_index, (natoms))
    ABI_ALLOCATE(self%typat,(ntypat))
    ABI_ALLOCATE(self%znucl, (ntypat))

    self%acell(:)=acell(:)
    self%rprimd(:,:)=rprimd(:,:)
    self%xred(:,:)=xred(:,:)
    self%spin_index(:)=spin_index(:)
    self%ntypat=ntypat
    self%typat(:)=typat(:)
    self%znucl(:)=znucl(:)
  end subroutine set_atomic_structure
  !!***


  !!****f* m_spin_hist/set_params
  !!
  !! NAME
  !! set_params
  !!
  !! FUNCTION
  !! 
  !! set parameters for spin_hist_t
  !!
  !! INPUTS
  !! spin_nctime=number of step between two write to netcdf hist file
  !! spin_temperate= temperature of spin
  !!
  !! OUTPUT
  !! hist <type(spin_hist_t)()> = spin hist type
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  
  subroutine set_params(self, spin_nctime, spin_temperature)

    class(spin_hist_t), intent(inout) :: self
    integer, intent(in) :: spin_nctime
    real(dp), intent(in) :: spin_temperature
    self%spin_nctime= spin_nctime
    self%spin_temperature=spin_temperature
  end subroutine set_params
!!***

  !!****f* m_spin_hist/finalize
  !!
  !! NAME
  !! finalize
  !!
  !! FUNCTION
  !! 
  !! free memory for spin_hist_t
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !! hist <type(spin_hist_t)()> = spin hist type
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine finalize(self)

    class(spin_hist_t) , intent(inout) :: self 

    if (allocated(self%xred)) then
       ABI_DEALLOCATE(self%xred)
    end if
    if (allocated(self%typat)) then
       ABI_DEALLOCATE(self%typat)
    end if
    if (allocated(self%znucl)) then
       ABI_DEALLOCATE(self%znucl)
    end if
    if (allocated(self%spin_index)) then
       ABI_DEALLOCATE(self%spin_index)
    end if
    if (allocated(self%heff)) then
       ABI_DEALLOCATE(self%heff)
    end if
    if (allocated(self%snorm)) then
       ABI_DEALLOCATE(self%snorm)
    end if
    if (allocated(self%S)) then
       ABI_DEALLOCATE(self%S)
    end if
    if (allocated(self%dSdt)) then
       ABI_DEALLOCATE(self%dSdt)
    end if
    if (allocated(self%etot)) then
       ABI_DEALLOCATE(self%etot)
    end if
    if (allocated(self%entropy)) then
       ABI_DEALLOCATE(self%entropy)
    end if
    if (allocated(self%time)) then
       ABI_DEALLOCATE(self%time)
    end if
    if (allocated(self%itime)) then
       ABI_DEALLOCATE(self%itime)
    end if
    if (allocated(self%ihist_latt)) then
       ABI_DEALLOCATE(self%ihist_latt)
    end if

  end subroutine finalize
!!***


  
  !!****f* m_spin_hist/get_S
  !!
  !! NAME
  !! get_S
  !!
  !! FUNCTION
  !! 
  !! get the S for step. step=0 is current. step=-1 is last...
  !!
  !! INPUTS
  !! hist <type(spin_hist_t)()> = spin hist type
  !! step = index of step. current step is 0. last step is -1. 
  !! OUTPUT
  !! S(3, nspin)=spin orientations at step
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function get_S(self, step) result(S)
    class(spin_hist_t), intent(inout) :: self
    integer, intent(in), optional:: step
    real(dp) :: S(3, self%nspin)
    integer :: i, j
    if (.not. present(step)) then
       j=0
    else
       j=step
    end if
    i=self%findIndex(step=j)
    S(:,:)=self%S(:,:,i)
  end function get_S
  !!***

  !!****f* m_spin_hist/inc1
  !!
  !! NAME
  !! inc1
  !!
  !! FUNCTION
  !! 
  !! time counter increase
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!   hist <type(spin_hist_t)()> = spin hist type
  !! PARENTS
!!
  !! CHILDREN
!!      self%inc1
!!
  !! SOURCE
  subroutine inc1(self)

    class(spin_hist_t), intent(inout) :: self
    if(self%ihist_prev ==0 ) then
        self%itime(self%ihist)=1
    else
        self%itime(self%ihist)=self%itime(self%ihist_prev)+1
    endif
    self%ihist_prev=self%ihist
    self%ihist=self%findIndex(1)
  end subroutine inc1
  !!***


  !!***f* m_spin_hist/findIndex
  !!
  !! NAME
  !! get_findIndex
  !!
  !! FUNCTION
  !! get the index of the step in the self%S array
  !! INPUTS
  !!
  !! OUTPUT
  !!   index: the index of the step in the self%S array.
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function findIndex(self, step) result(index)

    class(spin_hist_t), intent(inout) :: self
    integer , intent(in) :: step
    integer :: index
    !Local variables-------------------------------
    !scalars
    integer :: mxhist
    !arrays
    character(len=500) :: msg
    ! *************************************************************

    mxhist = self%mxhist
    if ((mxhist ==1.and.step/=+1).or.&
         &    (mxhist /=1.and.abs(step) >=mxhist)) then
       write(msg,'(a,I0,2a)')' The requested step must be less than ',mxhist,ch10,&
            &                     'Action: increase the number of history store in the hist'
       ABI_BUG(msg)
    end if
    index= mod(self%ihist+step, self%mxhist)+1
  end function findIndex
  !!***


  !!***f* m_spin_hist/set_vars
  !!
  !! NAME
  !! get_set_vars
  !!
  !! FUNCTION
  !! put the data into hist
  !! INPUTS
  !! S(3, nspin)=spin orientation
  !! Snorm(nspin)=spin amplitude
  !! dSdt(3,nspin)= dS/dt
  !! Heff(3, nspin) = effective magnetic field
  !! etot = total energy
  !! entropy = entropy
  !! time = time (note: not index of time)
  !! ihist_latt = index of lattice dynamics step.
  !! inc = whether this step is finished. If true, increment counter.
  !! OUTPUT
  !!   index: the index of the step in the self%S array.
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine set_vars(self, S, Snorm, dSdt, Heff, etot, entropy, time, ihist_latt, inc)

    class(spin_hist_t), intent(inout) ::self
    real(dp), optional, intent(in) :: S(3, self%nspin), Snorm(self%nspin), dSdt(3, self%nspin), &
        &  Heff(3, self%nspin), etot, entropy, time
    integer, optional :: ihist_latt
    logical, intent(in), optional :: inc
    integer :: ihist
    ihist=self%ihist
    if(present(inc)) then
       if (inc) then
          call self%inc1()
       end if
    end if
    if(present(S)) then
       self%S(:, :, ihist)=S(:,:)
    end if
    if(present(Snorm)) then
       self%Snorm(:,  ihist)=Snorm(:)
    endif
    if(present(dSdt)) then
       self%dSdt(:, :, ihist)=dSdt(:,:)
    end if
    if(present(Heff)) then
       self%Heff(:, :, ihist)=Heff(:,:)
    end if
    if(present(etot)) then
       self%etot(ihist)=etot
    end if
    if(present(entropy)) then
       self%entropy(ihist)=entropy
    end if
    if(present(time)) then
       self%time( ihist)=time
    end if
    if(present(ihist_latt)) then
       self%ihist_latt(ihist)=ihist_latt
    endif
  end subroutine set_vars
  !!***

end module m_spin_hist
