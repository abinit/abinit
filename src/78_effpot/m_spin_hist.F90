!{\src2tex{textfont=tt}}
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
!! * spin_hist_t_init
!! * spin_hist_t_free
!! * spin_hist_t
!! * spin_hist_t_get_S
!! * spin_hist_t_findIndex
!! * spin_hist_t_set_vars
!! * spin_hist_t_set_params
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

! TODO hexu:
! sync ihist_latt when with lattice dynamics
! add average , variance, etc (should they be here?)
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
  !! nspins: number of magnetic atoms
  !! * acell(3)         : Acell (acell , rprimd, xred: only initial value kept if there is!!  no lattice dynamics. Other wise for each step, the corresponding lattice step number is kept)
  !! * rprimd(3,3)      : Rprimd
  !! * xred(3,natoms)    : Xred
  !! * index_spin     : the index of atom in spin model, -1 if it is not in the spin model 
  !! * heff(3,nspins,mxhist)   : effective magnetic field (cartesian)
  !! * snorm(nspins, mxhist) : magnetitude of spin.
  !! * S(3,nspins,mxhist)   : spin orientation of atoms (cartesian)
  !! * dSdt(3, nspins, mxhist) : dS/dt (cartesian)
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

     integer :: nspins, nspins_prim
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
     !heff(3, nspins, mxhist)
     real(dp), allocatable :: heff(:, :, :)
     !snorm(nspins, mxhist)
     real(dp), allocatable :: snorm(:, :)

     !S(3, nspins, mxhist)
     real(dp), allocatable :: S(:, :, :)
     !dSdt(3, nspins, mxhist)
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

  end type spin_hist_t
  !!***

  public :: spin_hist_t_init
  public :: spin_hist_t_free
  public :: spin_hist_t_reset
  public :: spin_hist_t_get_S
  public :: spin_hist_t_findIndex
  public :: spin_hist_t_set_vars
  public :: spin_hist_t_set_params
  public :: spin_hist_t_inc
  !public :: spinhist2var
  !public :: var2spinhist
  !public :: write_sd_hist
  !public :: read_md_hist
  !public :: get_dims_spinhist

contains


!!****f* m_spin_hist/spin_hist_t_init
!!
!! NAME
!! spin_hist_t_init
!!
!! FUNCTION
!! initialize spin hist
!!
!! INPUTS
!! nspins = number of magnetic atoms
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

  subroutine spin_hist_t_init(hist, nspins, mxhist, has_latt)

    implicit none
    class(spin_hist_t), intent(inout) :: hist
    integer, intent(in) :: nspins, mxhist
    logical, intent(in) :: has_latt
    !integer, optional,  intent(in) :: calc_traj_obs, calc_thermo_obs, calc_correlation_obs

    hist%nspins=nspins
    hist%ntypat=0
    hist%ihist=1
    hist%ihist_prev=0
    hist%mxhist=mxhist
    hist%natoms=0
    hist%has_latt=has_latt

    !print *, "initialize HIST spin"
    ABI_ALLOCATE(hist%heff, (3, nspins, mxhist))
    ABI_ALLOCATE(hist%snorm, (nspins, mxhist))
    ABI_ALLOCATE(hist%S, (3, nspins, mxhist))
    ABI_ALLOCATE(hist%dSdt, (3, nspins, mxhist))

    ABI_ALLOCATE(hist%etot, (mxhist))
    ABI_ALLOCATE(hist%entropy, (mxhist))
    ABI_ALLOCATE(hist%time, (mxhist))
    ABI_ALLOCATE(hist%itime, (mxhist))

    ABI_ALLOCATE(hist%ihist_latt, (mxhist))



    ! TODO: add observable allocation here.

    hist%etot(1) =zero
    hist%entropy(1) =zero
    hist%time(1) =zero

    !hist%acell(:)=zero
    !hist%rprimd(:, :)=zero
    !hist%xred(:,:) =zero
    hist%heff(:,:,1)=zero
    hist%S(:,:,1)=zero
    hist%dSdt(:,:,1)=zero
    hist%snorm(:,1)=zero
    !print *, "Initialization spin hist finished"
  end subroutine spin_hist_t_init
!!***

  subroutine spin_hist_t_reset(hist, array_to_zero)

    implicit none
    class(spin_hist_t), intent(inout) :: hist
    logical :: array_to_zero
    hist%ntypat=0
    hist%ihist=1
    hist%ihist_prev=0
    hist%natoms=0

    hist%etot(1) =zero
    hist%entropy(1) =zero
    hist%time(1) =zero

    if(array_to_zero) then
       hist%heff(:,:,1)=zero
       hist%S(:,:,1)=zero
       hist%dSdt(:,:,1)=zero
       hist%snorm(:,1)=zero
       hist%Cv( 1)=zero
       hist%sp_corr_func(:, :, 1)=zero
    endif


  end subroutine spin_hist_t_reset

  !!****f* m_spin_hist/spin_hist_t_set_atomic_structure
  !!
  !! NAME
  !! spin_hist_t_set_atomic_structure
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
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_hist_t_set_atomic_structure(hist, acell, rprimd, xred, spin_index, ntypat,  typat, znucl)

    class(spin_hist_t), intent(inout) :: hist
    real(dp), intent(in) :: acell(3), rprimd(3,3), xred(:,:), znucl(:)
    integer, intent(in):: spin_index(:), ntypat, typat(:)
    integer :: natoms
    natoms=size(typat)
    ABI_ALLOCATE(hist%xred, (3, natoms))
    ABI_ALLOCATE(hist%spin_index, (natoms))
    ABI_ALLOCATE(hist%typat,(ntypat))
    ABI_ALLOCATE(hist%znucl, (ntypat))

    hist%acell(:)=acell(:)
    hist%rprimd(:,:)=rprimd(:,:)
    hist%xred(:,:)=xred(:,:)
    hist%spin_index(:)=spin_index(:)
    hist%ntypat=ntypat
    hist%typat(:)=typat(:)
    hist%znucl(:)=znucl(:)
  end subroutine spin_hist_t_set_atomic_structure
  !!***


  !!****f* m_spin_hist/spin_hist_t_set_params
  !!
  !! NAME
  !! spin_hist_t_set_params
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
  
  subroutine spin_hist_t_set_params(hist, spin_nctime, spin_temperature)

    class(spin_hist_t), intent(inout) :: hist
    integer, intent(in) :: spin_nctime
    real(dp), intent(in) :: spin_temperature
    hist%spin_nctime= spin_nctime
    hist%spin_temperature=spin_temperature
  end subroutine spin_hist_t_set_params
!!***

  !!****f* m_spin_hist/spin_hist_t_free
  !!
  !! NAME
  !! spin_hist_t_free
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
  subroutine spin_hist_t_free(hist)

    class(spin_hist_t) , intent(inout) :: hist

    if (allocated(hist%xred)) then
       ABI_DEALLOCATE(hist%xred)
    end if
    if (allocated(hist%typat)) then
       ABI_DEALLOCATE(hist%typat)
    end if
    if (allocated(hist%znucl)) then
       ABI_DEALLOCATE(hist%znucl)
    end if
    if (allocated(hist%spin_index)) then
       ABI_DEALLOCATE(hist%spin_index)
    end if
    if (allocated(hist%heff)) then
       ABI_DEALLOCATE(hist%heff)
    end if
    if (allocated(hist%snorm)) then
       ABI_DEALLOCATE(hist%snorm)
    end if
    if (allocated(hist%S)) then
       ABI_DEALLOCATE(hist%S)
    end if
    if (allocated(hist%dSdt)) then
       ABI_DEALLOCATE(hist%dSdt)
    end if
    if (allocated(hist%etot)) then
       ABI_DEALLOCATE(hist%etot)
    end if
    if (allocated(hist%entropy)) then
       ABI_DEALLOCATE(hist%entropy)
    end if
    if (allocated(hist%time)) then
       ABI_DEALLOCATE(hist%time)
    end if
    if (allocated(hist%itime)) then
       ABI_DEALLOCATE(hist%itime)
    end if
    if (allocated(hist%ihist_latt)) then
       ABI_DEALLOCATE(hist%ihist_latt)
    end if

  end subroutine spin_hist_t_free
!!***


  
  !!****f* m_spin_hist/spin_hist_t_get_S
  !!
  !! NAME
  !! spin_hist_t_get_S
  !!
  !! FUNCTION
  !! 
  !! get the S for step. step=0 is current. step=-1 is last...
  !!
  !! INPUTS
  !! hist <type(spin_hist_t)()> = spin hist type
  !! step = index of step. current step is 0. last step is -1. 
  !! OUTPUT
  !! S(3, nspins)=spin orientations at step
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function spin_hist_t_get_S(hist, step) result(S)

    class(spin_hist_t), intent(inout) :: hist
    integer, intent(in), optional:: step
    real(dp) :: S(3, hist%nspins)
    integer :: i, j
    if (.not. present(step)) then
       j=0
    else
       j=step
    end if
    i=spin_hist_t_findIndex(hist,step=j)
    S(:,:)=hist%S(:,:,i)
  end function spin_hist_t_get_S
  !!***

  !!****f* m_spin_hist/spin_hist_t_inc
  !!
  !! NAME
  !! spin_hist_t_get_inc
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
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_hist_t_inc(hist)

    class(spin_hist_t), intent(inout) :: hist
    if(hist%ihist_prev ==0 ) then
        hist%itime(hist%ihist)=1
    else
        hist%itime(hist%ihist)=hist%itime(hist%ihist_prev)+1
    endif
    hist%ihist_prev=hist%ihist
    hist%ihist=spin_hist_t_findIndex(hist, 1)
  end subroutine spin_hist_t_inc
  !!***


  !!***f* m_spin_hist/spin_hist_t_findIndex
  !!
  !! NAME
  !! spin_hist_t_get_findIndex
  !!
  !! FUNCTION
  !! get the index of the step in the hist%S array
  !! INPUTS
  !!
  !! OUTPUT
  !!   index: the index of the step in the hist%S array.
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function spin_hist_t_findIndex(hist, step) result(index)

    type(spin_hist_t), intent(inout) :: hist
    integer , intent(in) :: step
    integer :: index
    !Local variables-------------------------------
    !scalars
    integer :: mxhist
    !arrays
    character(len=500) :: msg
    ! *************************************************************

    mxhist = hist%mxhist
    if ((mxhist ==1.and.step/=+1).or.&
         &    (mxhist /=1.and.abs(step) >=mxhist)) then
       write(msg,'(a,I0,2a)')' The requested step must be less than ',mxhist,ch10,&
            &                     'Action: increase the number of history store in the hist'
       MSG_BUG(msg)
    end if
    index= mod(hist%ihist+step, hist%mxhist)+1
  end function spin_hist_t_findIndex
  !!***


  !!***f* m_spin_hist/spin_hist_t_set_vars
  !!
  !! NAME
  !! spin_hist_t_get_set_vars
  !!
  !! FUNCTION
  !! put the data into hist
  !! INPUTS
  !! S(3, nspins)=spin orientation
  !! Snorm(nspins)=spin amplitude
  !! dSdt(3,nspins)= dS/dt
  !! Heff(3, nspins) = effective magnetic field
  !! etot = total energy
  !! entropy = entropy
  !! time = time (note: not index of time)
  !! ihist_latt = index of lattice dynamics step.
  !! inc = whether this step is finished. If true, increment counter.
  !! OUTPUT
  !!   index: the index of the step in the hist%S array.
  !! PARENTS
  !!      m_spin_hist
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_hist_t_set_vars(hist, S, Snorm, dSdt, Heff, etot, entropy, time, ihist_latt, inc)

    class(spin_hist_t), intent(inout) :: hist
    real(dp), optional, intent(in) :: S(3, hist%nspins), Snorm(hist%nspins), dSdt(3, hist%nspins), &
        &  Heff(3, hist%nspins), etot, entropy, time
    integer, optional :: ihist_latt
    logical, intent(in), optional :: inc
    integer :: ihist
    ihist=hist%ihist
    !print *, "set spin hist vars: ihist=", ihist
    if(present(inc) .and. inc) then
       call spin_hist_t_inc(hist)
    end if
    if(present(S)) then
       hist%S(:, :, ihist)=S(:,:)
    end if
    if(present(Snorm)) then
       hist%Snorm(:,  ihist)=Snorm(:)
    endif
    if(present(dSdt)) then
       hist%dSdt(:, :, ihist)=dSdt(:,:)
    end if
    if(present(Heff)) then
       hist%Heff(:, :, ihist)=Heff(:,:)
    end if
    if(present(etot)) then
       hist%etot(ihist)=etot
    end if
    if(present(entropy)) then
       hist%entropy(ihist)=entropy
    end if
    if(present(time)) then
       hist%time( ihist)=time
    end if
    if(present(ihist_latt)) then
       hist%ihist_latt(ihist)=ihist_latt
    endif
  end subroutine spin_hist_t_set_vars
  !!***

end module m_spin_hist
