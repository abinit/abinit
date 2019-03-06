!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_manager
!! NAME
!! m_multibinit_manager
!!
!! FUNCTION
!! This module contains the manager type, which is a thin layer above ALL 
!! TODO: the structure of this is yet to be discussed
!!
!!
!! Datatypes:
!!
!! * mb_manager_t
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

module m_multibinit_manager
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_init10, only: init10
  use m_mathfuncs, only: diag
  use m_multibinit_global
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_unitcell, only : unitcell_t
  use m_supercell_maker, only: supercell_maker_t
  use m_multibinit_supercell, only: mb_supercell_t
  use m_primitive_potential_list, only: primitive_potential_list_t
  use m_abstract_potential, only: abstract_potential_t
  use m_potential_list, only: potential_list_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_effpot, only : lattice_effpot_t
  use m_spin_potential, only : spin_potential_t
  use m_lattice_mover, only : lattice_mover_t
  use m_spin_mover, only : spin_mover_t
  use m_spin_lattice_coupling_effpot, only : spin_lattice_coupling_effpot_t
  implicit none
  private

  !!***

!-------------------------------------------------------------------!
! Multibinit manager
!-------------------------------------------------------------------!
  type, public :: mb_manager_t
     character(len=fnlen) :: filenames(17)
     type(multibinit_dtset_type) :: params
     type(supercell_maker_t) :: sc_maker
     type(unitcell_t) :: unitcell
     type(mb_supercell_t) :: supercell
     type(primitive_potential_list_t) :: prim_pots
     type(potential_list_t), pointer :: pots

     type(lattice_mover_t) :: lattice_mover
     type(spin_mover_t) :: spin_mover
     ! type(lwf_mover_t) :: lwf_mover

   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: read_params    ! parse input file
     procedure :: prepare_params ! process the parameters. e.g. convert unit, set some flags, etc.
     procedure :: read_potentials ! read primitve cell and potential
     procedure :: fill_supercell
     procedure :: set_movers
     procedure :: run_dynamics
     procedure :: run
     procedure :: run_all
  end type mb_manager_t

contains
  !-------------------------------------------------------------------!
  ! initialize
  !-------------------------------------------------------------------!
  subroutine initialize(self, filenames)
    class(mb_manager_t), intent(inout) :: self
    character(len=fnlen), intent(inout) :: filenames(17)
    self%filenames=filenames
    call self%read_params()
    call self%prepare_params()
    ! read potentials from
  end subroutine initialize


  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(mb_manager_t), intent(inout) :: self
    !call self%params%finalize()
    call self%sc_maker%finalize()
    call self%unitcell%finalize()
    call self%supercell%finalize()
    call self%prim_pots%finalize()
    call self%pots%finalize()
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! read_params: read parameters from input file
  !-------------------------------------------------------------------!
  subroutine read_params(self)
    class(mb_manager_t), intent(inout) :: self
    character(len=500) :: message
    integer :: filetype,ii,lenstr
    integer :: natom,nph1l,nrpt,ntypat
    integer :: option
    character(len=strlen) :: string
    integer:: ierr
    !To automate a maximum calculation, multibinit reads the number of atoms
    !in the file (ddb or xml). If DDB file is present in input, the ifc calculation
    !will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10
    write(message, '(6a)' )' Read the information in the reference structure in ',ch10,&
         & '-',trim(self%filenames(3)),ch10,' to initialize the multibinit input'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')

    call effective_potential_file_getDimSystem(self%filenames(3),natom,ntypat,nph1l,nrpt)

    !Read the input file, and store the information in a long string of characters
    !strlen from defs_basis module
    option=1
    if (iam_master) then
       call instrng (self%filenames(1),lenstr,option,strlen,string)
       !To make case-insensitive, map characters to upper case:
       call inupper(string(1:lenstr))

       !Check whether the string only contains valid keywords
       call chkvars(string)

    end if

    call xmpi_bcast(string,master, comm, ierr)
    call xmpi_bcast(lenstr,master, comm, ierr)

    !Read the input file
    call invars10(self%params,lenstr,natom,string)

    if (iam_master) then
       !  Echo the inputs to console and main output file
       call outvars_multibinit(self%params,std_out)
       call outvars_multibinit(self%params,ab_out)
    end if
  end subroutine read_params

  !-------------------------------------------------------------------!
  ! prepare_params: after read, something has to be done:
  ! e.g. unit conversion
  !-------------------------------------------------------------------!
  subroutine prepare_params(self)
    class(mb_manager_t), intent(inout) :: self
    ! Kelvin to Hartree (In input file, the spin temperature is in K.
    ! convert to a.u.)
    self%params%spin_temperature = self%params%spin_temperature/Ha_K
    self%params%spin_temperature_start=self%params%spin_temperature_start/Ha_K
    self%params%spin_temperature_end=self%params%spin_temperature_end/Ha_K
  end subroutine prepare_params


  !-------------------------------------------------------------------!
  ! Read potentials from file, if needed by dynamics
  !-------------------------------------------------------------------!
  subroutine read_potentials(self)
    class(mb_manager_t), intent(inout) :: self
    !TODO: unitcell.
    if(params%spin_dynamics>0) then
       call spin_pot%read_from_files(self%params, fnames)
    end if
    if(params%dynamics>0) then
       !TODO: LATT
    endif
    !TODO: LWF
  end subroutine read_potentials

  !-------------------------------------------------------------------!
  ! fill supercell. Both primitive cell and potential
  !-------------------------------------------------------------------!
  subroutine fill_supercell(self)
    class(mb_manager_t), target, intent(inout) :: self
    ! unitcell
    call self%unitcell%fill_supercell(self%sc_maker, self%supercell)
    call self%prim_pots%fill_supercell_list(self%sc_maker,self%pots)
  end subroutine fill_supercell

  !-------------------------------------------------------------------!
  ! Fit lattic model
  !-------------------------------------------------------------------!
  subroutine fit_lattice_model(self)
    class(mb_manager_t), intent(inout) :: self
    !TODO:
  end subroutine fit_lattice_model


  !-------------------------------------------------------------------!
  ! initialize movers needed.
  !-------------------------------------------------------------------!
  subroutine set_movers(self)
    class(mb_manager_t), intent(inout) :: self
    if (self%params%spin_dynamics>0) then
       call self%spin_mover%initialize(self%params, nspin=self%supercell%nspin)
    end if

    if (self%params%dynamics>0) then
       call self%lattice_mover%initialize(self%params, self%filenames)
    end if

    ! TODO: LWF MOVER
  end subroutine set_movers


  !-------------------------------------------------------------------!
  ! Run dynamics
  !-------------------------------------------------------------------!
  subroutine run_dynamics(self)
    class(mb_manager_t), intent(inout) :: self
    call self%prim_pots%initialize()
    call self%sc_maker%initialize(diag(self%params%ncell))
    ! read params
    call self%read_potentials()
    call self%fill_supercell()

    call self%set_movers()
  end subroutine run_dynamics

  !-------------------------------------------------------------------!
  ! Run all jobs
  !-------------------------------------------------------------------!
  subroutine run(self)
    class(mb_manager_t), intent(inout) :: self
    ! if ... fit lattice model
    ! if ... fit lwf model
    ! if ... run dynamics...
    ! if ...
  end subroutine run
 !----------------------------------------------------------------------------------------------------
  ! This does everything from initialize to finalize
  subroutine run_all(self, filenames)
    class(mb_manager_t), intent(inout) :: self
    character(len=fnlen), intent(inout) :: filenames(17)
  end subroutine run_all

end module m_multibinit_manager
