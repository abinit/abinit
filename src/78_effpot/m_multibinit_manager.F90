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
  !use m_abicore
  use m_errors
  use m_xmpi

  use m_init10, only: init10
  use m_mathfuncs, only: diag
  use m_multibinit_global
  use m_multibinit_dataset, only: multibinit_dtset_type, invars10, outvars_multibinit
  use m_unitcell, only : unitcell_t
  use m_supercell_maker, only: supercell_maker_t
  use m_multibinit_supercell, only: mb_supercell_t
  use m_primitive_potential_list, only: primitive_potential_list_t
  use m_primitive_potential, only: primitive_potential_t
  use m_spin_primitive_potential, only: spin_primitive_potential_t
  use m_abstract_potential, only: abstract_potential_t
  use m_potential_list, only: potential_list_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_effpot, only : lattice_effpot_t
  use m_spin_potential, only : spin_potential_t
  use m_lattice_mover, only : lattice_mover_t
  use m_spin_mover, only : spin_mover_t
  ! TODO : should these be moved into spin mover?
  use m_spin_hist, only: spin_hist_t
  use m_spin_ncfile, only: spin_ncfile_t
  use m_spin_observables, only: spin_observable_t

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
     type(potential_list_t) :: pots

     type(lattice_mover_t) :: lattice_mover
     type(spin_mover_t) :: spin_mover
     ! type(lwf_mover_t) :: lwf_mover

     type(spin_hist_t):: spin_hist
     type(spin_ncfile_t) :: spin_ncfile
     type(spin_observable_t) :: spin_ob
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
    call self%spin_mover%finalize()
    call self%lattice_mover%finalize()
    !call self%lwf_mover%finalize()
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! read_params: read parameters from input file
  !-------------------------------------------------------------------!
  subroutine read_params(self)
    use m_fstrings,   only : replace, inupper
    use m_parser, only: instrng
    use m_dtset,      only : chkvars
    use m_effective_potential_file
    use m_effective_potential
    class(mb_manager_t), intent(inout) :: self
    character(len=500) :: message
    integer :: lenstr
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
    type(spin_primitive_potential_t), pointer :: spin_pot
    !class(primitive_potential_t), pointer :: t
    call self%unitcell%initialize()

    ! latt : TODO

    ! spin
    if(self%params%spin_dynamics>0) then
       allocate(spin_pot)
       call spin_pot%initialize()
       call spin_pot%set_atoms(self%unitcell)
       call spin_pot%load_from_files(self%params, self%filenames)
       call self%prim_pots%append(spin_pot)

       ! !DEBUG
       ! print*, "nspin", self%unitcell%spin%nspin
       ! print*, "ms", self%unitcell%spin%ms
       ! print*, "gyro_ratio", self%unitcell%spin%gyro_ratio
       ! print*, "damping_factor", self%unitcell%spin%damping_factor
       ! print*, "positions", self%unitcell%spin%spin_positions
       ! !t=>self%prim_pots%data(1)%obj
       ! select type (t=>self%prim_pots%data(1)%obj)
       !    type is(spin_primitive_potential_t)
       !       print *, "nspin", t%nspin
       !       call t%coeff%print()
       ! end select
    end if
    if(self%params%dynamics>0) then
       !TODO: LATT
    endif

    !LWF : TODO
  end subroutine read_potentials

  !-------------------------------------------------------------------!
  ! fill supercell. Both primitive cell and potential
  !-------------------------------------------------------------------!
  subroutine fill_supercell(self)
    class(mb_manager_t), target, intent(inout) :: self
    class(abstract_potential_t), pointer :: q
    integer :: i
    ! unitcell
    call self%unitcell%fill_supercell(self%sc_maker, self%supercell)
    ! print *, "sc_nspin",  self%supercell%nspin
    ! print *, "sc_ms",  self%supercell%ms
    ! print *, "sc_gyro",  self%supercell%gyro_ratio
    ! print *, "sc_damping",  self%supercell%gilbert_damping
    ! print *, "ispin_prim",  self%supercell%ispin_prim
    ! print *, "sc_rvec",  self%supercell%rvec
    call self%pots%initialize()
    call self%pots%set_supercell(self%supercell)
    call self%prim_pots%fill_supercell_list(self%sc_maker,self%pots)

    do i=1, self%pots%size
       q=>self%pots%data(i)%obj
    end do
    !select type (t=>self%prim_pots%data(1)%obj)
    !    type is(spin_primitive_potential_t)
    !       print *, "fill spin pot supercell", t%nspin
    !       call t%fill_supercell(self%sc_maker, q)
    !end select
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
  ! set_spin_mover: prepare for spin mover
  ! Should these be moved into mover?
  ! including the observables and hist
  !-------------------------------------------------------------------!
  subroutine set_spin_mover(self)
    class(mb_manager_t), intent(inout) :: self
    real(dp):: mfield(3, self%supercell%nspin), damping(self%supercell%nspin)
    integer ::  i
    ! params -> mover
    ! params -> calculator

    call xmpi_bcast(self%params%spin_damping, master, comm, ierr)
    if (self%params%spin_damping >=0) then
       damping(:)= self%params%spin_damping
    end if
    call self%spin_mover%set_langevin_params(temperature=self%params%spin_temperature, &
         & damping=damping, ms=self%supercell%ms, gyro_ratio=self%supercell%gyro_ratio)

    call self%spin_mover%initialize(self%params, self%supercell%nspin )

    if(iam_master) then
       call self%spin_hist%initialize(nspin=self%nspin, &
            &   mxhist=3, has_latt=.False.)
       call self%spin_hist%set_params(spin_nctime=self%params%spin_nctime, &
            &     spin_temperature=self%params%spin_temperature)
       call self%spin_mover%set_hist(self%spin_hist)
    endif

    call self%set_initial_spin()
    call self%spin_mover%set_langevin_params(gyro_ratio=self%spin_calculator%supercell%gyro_ratio, &
         & damping=self%spin_calculator%supercell%gilbert_damping, ms=self%spin_calculator%supercell%ms )

    if(iam_master) then
       call self%spin_ob%initialize(self%spin_calculator%supercell, self%params)
       call spin_model_t_prepare_ncfile(self, self%spin_ncfile, trim(self%out_fname)//'_spinhist.nc')
       call self%spin_ncfile%write_one_step(self%spin_hist)
    endif


  end subroutine set_spin_mover

  !-------------------------------------------------------------------!
  ! Run dynamics
  !-------------------------------------------------------------------!
  subroutine run_dynamics(self)
    class(mb_manager_t), intent(inout) :: self
    call self%prim_pots%initialize()
    call self%sc_maker%initialize(diag(self%params%ncell))

    
    if(iam_master) then
       call self%spin_ob%initialize(self%supercell, self%params)
       call spin_model_t_prepare_ncfile(self, self%spin_ncfile, trim(self%out_fname)//'_spinhist.nc')
       call self%spin_ncfile%write_one_step(self%spin_hist)
    endif



    ! read params
    call self%read_potentials()
    call self%fill_supercell()
    call self%set_movers()

    print *, "Mover initialized"
  end subroutine run_dynamics

  !-------------------------------------------------------------------!
  ! Run all jobs
  !-------------------------------------------------------------------!
  subroutine run(self)
    class(mb_manager_t), intent(inout) :: self
    ! if ... fit lattice model
    ! if ... fit lwf model
    ! if ... run dynamics...
    call self%run_dynamics()
    ! if ...
  end subroutine run


  !-------------------------------------------------------------------!
  !run_all: THE function which does everything
  !         from the very begining to end.
  !-------------------------------------------------------------------!
  subroutine run_all(self, filenames)
    class(mb_manager_t), intent(inout) :: self
    character(len=fnlen), intent(inout) :: filenames(17)
    call self%initialize(filenames)
    call self%run()
    call self%finalize()
  end subroutine run_all

end module m_multibinit_manager
