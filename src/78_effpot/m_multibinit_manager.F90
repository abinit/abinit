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
  use m_multibinit_dataset, only: multibinit_dtset_type, invars10, &
       outvars_multibinit, multibinit_dtset_free
 ! random number generator
  use m_random_xoroshiro128plus, only: rng_t

 ! cells
  use m_supercell_maker, only: supercell_maker_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t

  ! primitive potential
  use m_primitive_potential_list, only: primitive_potential_list_t
  use m_primitive_potential, only: primitive_potential_t

  ! 
  use m_abstract_potential, only: abstract_potential_t
  use m_potential_list, only: potential_list_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_effpot, only : lattice_effpot_t

  ! Spin
  use m_spin_primitive_potential, only: spin_primitive_potential_t
  use m_spin_potential, only : spin_potential_t
  use m_spin_mover, only : spin_mover_t
  use m_mpi_scheduler, only: init_mpi_info
  ! TODO : should these be moved into spin mover?
  use m_spin_ncfile, only: spin_ncfile_t

  ! Lattice harmonic
  use m_lattice_harmonic_primitive_potential, only: lattice_harmonic_primitive_potential_t
  use m_lattice_harmonic_potential, only: lattice_harmonic_potential_t

  ! Lattice movers
  use m_lattice_mover, only: lattice_mover_t
  use m_lattice_langevin_mover, only: lattice_langevin_mover_t
  use m_lattice_verlet_mover, only: lattice_verlet_mover_t
  use m_lattice_berendsen_NVT_mover, only: lattice_berendsen_NVT_mover_t
  use m_lattice_berendsen_NPT_mover, only: lattice_berendsen_NPT_mover_t
  use m_lattice_dummy_mover, only: lattice_dummy_mover_t

  ! Spin lattice coupling
  use m_slc_primitive_potential, only: slc_primitive_potential_t
  use m_slc_potential, only : slc_potential_t

  implicit none
  private

  !!***

  !-------------------------------------------------------------------!
  ! Multibinit manager
  !-------------------------------------------------------------------!
  type, public :: mb_manager_t
     character(len=fnlen) :: filenames(17)
     ! pointer to parameters. it is a pointer because it is initialized outside manager
     type(multibinit_dtset_type), pointer :: params=>null() 
     type(supercell_maker_t) :: sc_maker  ! supercell maker
     type(mbcell_t) :: unitcell         ! unitcell
     type(mbsupercell_t) :: supercell   ! supercell
     type(primitive_potential_list_t) :: prim_pots  ! list of primitive potentials
     type(potential_list_t) :: pots     ! potential list
     ! a polymorphic lattice mover so multiple mover could be used.
     class(lattice_mover_t), pointer :: lattice_mover=> null()
     ! as for the spin, there is only one mover which has several methods
     type(spin_mover_t) :: spin_mover  
     ! type(lwf_mover_t) :: lwf_mover

     ! spin netcdf hist file
     type(spin_ncfile_t) :: spin_ncfile
     type(rng_t) :: rng

     ! TODO: this is temporary. Remove after moving to multibinit_main2
     ! It means the parsing of the params are already done outside the manager.
     logical :: use_external_params=.True.

     logical :: has_displacement = .False.
     logical :: has_strain = .False.
     logical :: has_spin = .False.
     logical :: has_lwf = .False.
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: read_params    ! parse input file
     procedure :: prepare_params ! process the parameters. e.g. convert unit, set some flags, etc.
     procedure :: read_potentials ! read primitve cell and potential
     procedure :: fill_supercell
     procedure :: set_movers
     procedure :: set_lattice_mover
     procedure :: run_spin_dynamics
     procedure :: run_MvT
     procedure :: run_lattice_dynamics
     procedure :: run_spin_latt_dynamics
     procedure :: run_coupled_spin_latt_dynamics
     procedure :: run
     procedure :: run_all
  end type mb_manager_t

contains
  !-------------------------------------------------------------------!
  ! initialize
  !-------------------------------------------------------------------!
  subroutine initialize(self, filenames,params)
    class(mb_manager_t), intent(inout) :: self
    character(len=fnlen), intent(inout) :: filenames(17)
    type(multibinit_dtset_type), target, optional, intent(in) :: params
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    integer :: i
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    self%filenames(:)=filenames(:)
    call xmpi_bcast(self%filenames, master, comm, ierr)

    !TODO: remove params as argument. It is here because the params are read
    ! in the multibinit_main function. Once we use multibinit_main2, remove it.
    if (present(params)) then
       self%params=>params
    else
       self%use_external_params=.False.
       ABI_MALLOC_SCALAR(self%params)
       call self%read_params()
    endif

    ! Initialize the random number generator
    call self%rng%set_seed([111111_dp, 2_dp])
    ! use jump so that each cpu generates independent random numbers.
    if(my_rank>0) then
       do i =1,my_rank
          call self%rng%jump()
       end do
    end if


    if(self%params%spin_dynamics>0) then
       self%has_spin=.True.
    endif
    if(self%params%dynamics >0) then
       self%has_displacement=.True.
       self%has_strain=.True.
    endif

    call self%prepare_params()
    ! read potentials from

    ! lwf

  end subroutine initialize


  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(mb_manager_t), intent(inout) :: self
    call self%sc_maker%finalize()
    call self%unitcell%finalize()
    call self%supercell%finalize()
    call self%prim_pots%finalize()
    call self%pots%finalize()
    call self%spin_mover%finalize()
    ! Note that lattice mover is a pointer.
    ! It might be null if there is no lattice part.
    if (associated(self%lattice_mover)) then
       call self%lattice_mover%finalize()
       nullify(self%lattice_mover)
    end if
    if(.not. self%use_external_params) then
       call multibinit_dtset_free(self%params)
       if (associated(self%params)) then
           ABI_FREE_SCALAR(self%params)
       endif
       !deallocate(self%params)
    endif
    nullify(self%params)
    !call self%lwf_mover%finalize()

    self%has_displacement=.False.
    self%has_strain=.False.
    self%has_spin=.False.
    self%has_lwf=.False.
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! read_params: read parameters from input file
  ! TODO: This function is copied from the implementation before using F03
  ! Some work need to be done to move the initialization of effective
  ! potential out of this function.
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
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    
    natom=0
    nph1l=0
    nrpt=0
    ntypat=0
    call init_mpi_info(master, iam_master, my_rank, comm, nproc)
    !To automate a maximum calculation, multibinit reads the number of atoms
    !in the file (ddb or xml). If DDB file is present in input, the ifc calculation
    !will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10
    if(iam_master) then
       write(message, '(6a)' )' Read the information in the reference structure in ',ch10,&
            & '-',trim(self%filenames(3)),ch10,' to initialize the multibinit input'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
    end if


    !FIXME: This should not be here.
    ! It is only for lattice potential
    if(.False.) then
       call effective_potential_file_getDimSystem(self%filenames(3),natom,ntypat,nph1l,nrpt)
    else
       natom=0
       ntypat=0
       nph1l=0
       nrpt=0
    endif

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
    if(self%has_spin) then
       self%params%spin_temperature = self%params%spin_temperature/Ha_K
       self%params%spin_temperature_start=self%params%spin_temperature_start/Ha_K
       self%params%spin_temperature_end=self%params%spin_temperature_end/Ha_K
    end if
    if(self%has_displacement) then
       self%params%temperature = self%params%temperature/Ha_K
    end if
    
  end subroutine prepare_params


  !-------------------------------------------------------------------!
  ! Read potentials from file, if needed by dynamics
  !-------------------------------------------------------------------!
  subroutine read_potentials(self)
    class(mb_manager_t), intent(inout) :: self
    class(primitive_potential_t), pointer :: spin_pot
    class(primitive_potential_t), pointer :: slc_pot
    class(primitive_potential_t), pointer :: lat_ham_pot

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    call self%unitcell%initialize()

    ! latt : TODO (replace this with full lattice)
    ! only toy harmonic part 
    if(self%params%dynamics>100) then
       ABI_DATATYPE_ALLOCATE_SCALAR(lattice_harmonic_primitive_potential_t, lat_ham_pot)
       select type(lat_ham_pot)
       type is (lattice_harmonic_primitive_potential_t)
          call lat_ham_pot%initialize(self%unitcell)
          call lat_ham_pot%load_from_files(self%params, self%filenames)
          call self%prim_pots%append(lat_ham_pot)
       end select
    end if

    ! spin
    call xmpi_bcast(self%params%spin_dynamics, master, comm, ierr)
    if(self%params%spin_dynamics>0) then
       ! The pointer will be allocated and added to the list
       ! and eventually deallocated by the list%finalize
       ABI_DATATYPE_ALLOCATE_SCALAR(spin_primitive_potential_t, spin_pot)

       ! One may wonder why unitcell does not read data from files
       ! That is because the spin_pot (which has an pointer to unitcell)
       ! read the file and set the spin unitcell.
       select type(spin_pot)
       type is (spin_primitive_potential_t)
          call spin_pot%initialize(self%unitcell)
          call spin_pot%load_from_files(self%params, self%filenames)
          call self%prim_pots%append(spin_pot)
       end select
    end if


    !LWF : TODO

    ! spin-lattice coupling
    if(self%params%slc_coupling>0) then
       ABI_DATATYPE_ALLOCATE_SCALAR(slc_primitive_potential_t, slc_pot)
       select type(slc_pot)
       type is (slc_primitive_potential_t)
          call slc_pot%initialize(self%unitcell)
          call slc_pot%load_from_files(self%params, self%filenames)
          call self%prim_pots%append(slc_pot)
       end select 
    endif
  end subroutine read_potentials

  !-------------------------------------------------------------------!
  ! fill supercell. Both primitive cell and potential
  !-------------------------------------------------------------------!
  subroutine fill_supercell(self)
    class(mb_manager_t), target, intent(inout) :: self

    ! build supercell structure
    !call self%unitcell%fill_supercell(self%sc_maker, self%supercell)
    call self%supercell%from_unitcell(self%sc_maker, self%unitcell)

    ! supercell potential
    call self%pots%initialize()
    call self%pots%set_supercell(self%supercell)
    call self%prim_pots%fill_supercell_list(self%sc_maker, self%params, self%pots)

    ! why do this twice.
    call self%pots%set_supercell(self%supercell)
  end subroutine fill_supercell

  !-------------------------------------------------------------------!
  ! Fit lattic model
  !-------------------------------------------------------------------!
  subroutine fit_lattice_model(self)
    class(mb_manager_t), intent(inout) :: self
    ABI_UNUSED_A(self)
    !TODO:
  end subroutine fit_lattice_model


  !-------------------------------------------------------------------!
  ! initialize movers which are needed.
  !-------------------------------------------------------------------!
  subroutine set_movers(self)
    class(mb_manager_t), intent(inout) :: self
    character(len=fnlen) :: fname
    if (self%params%spin_dynamics>0) then
       fname=trim(self%filenames(2))//"_spinhist_input.nc"
       call self%spin_mover%initialize(params=self%params,&
            & supercell=self%supercell, rng=self%rng, &
            & restart_hist_fname=fname)
    end if


    if (self%params%dynamics>0) then
       call self%set_lattice_mover()
    end if

    ! TODO: LWF MOVER
  end subroutine set_movers


  !-------------------------------------------------------------------!
  !Set_lattice_mover
  !-------------------------------------------------------------------!
  subroutine set_lattice_mover(self)
    class(mb_manager_t), intent(inout) :: self
    select case(self%params%dynamics)
    case (101)  ! Velocity Verlet (NVE)
       ABI_DATATYPE_ALLOCATE_SCALAR(lattice_verlet_mover_t, self%lattice_mover)
    case(102)   ! Langevin (NVT)
       ABI_DATATYPE_ALLOCATE_SCALAR(lattice_langevin_mover_t, self%lattice_mover)
    case(103)   ! Berendsen NVT
       ABI_DATATYPE_ALLOCATE_SCALAR(lattice_berendsen_NVT_mover_t, self%lattice_mover)
    case(104)   ! Berendsen NPT (not yet avaliable)
       ABI_DATATYPE_ALLOCATE_SCALAR(lattice_berendsen_NPT_mover_t, self%lattice_mover)
    case(120)   ! Dummy mover (Do not move atoms, For test only.)
       ABI_DATATYPE_ALLOCATE_SCALAR(lattice_dummy_mover_t, self%lattice_mover)
    end select
    call self%lattice_mover%initialize(params=self%params, supercell=self%supercell, rng=self%rng)
  end subroutine set_lattice_mover

 



  !-------------------------------------------------------------------!
  ! Run dynamics
  !-------------------------------------------------------------------!
  subroutine run_spin_dynamics(self)
    class(mb_manager_t), intent(inout) :: self
    call self%prim_pots%initialize()
    call self%read_potentials()
    call self%sc_maker%initialize(diag(self%params%ncell))
    call self%fill_supercell()
    call self%set_movers()
    call self%spin_mover%set_ncfile_name(self%params, self%filenames(2))
    call self%spin_mover%run_time(self%pots)
    call self%spin_mover%spin_ncfile%close()
  end subroutine run_spin_dynamics


  subroutine run_MvT(self)
    class(mb_manager_t), intent(inout) :: self
    call self%prim_pots%initialize()
    call self%sc_maker%initialize(diag(self%params%ncell))
    call self%read_potentials()
    call self%fill_supercell()
    call self%set_movers()
    call self%spin_mover%run_MvT(self%pots, self%filenames(2))
  end subroutine run_MvT


  !-------------------------------------------------------------------!
  ! Run lattice only dynamics
  !-------------------------------------------------------------------!
  subroutine run_lattice_dynamics(self)
    class(mb_manager_t), intent(inout) :: self
    call self%prim_pots%initialize()
    call self%read_potentials()
    call self%sc_maker%initialize(diag(self%params%ncell))
    call self%fill_supercell()
    call self%set_movers()
    call self%lattice_mover%run_time(self%pots)
  end subroutine run_lattice_dynamics

  !-------------------------------------------------------------------!
  ! Run  spin and lattice dynamics sequentially.
  !-------------------------------------------------------------------!
  subroutine run_spin_latt_dynamics(self)
    class(mb_manager_t), intent(inout) :: self
    call self%prim_pots%initialize()
    call self%read_potentials()
    call self%sc_maker%initialize(diag(self%params%ncell))
    call self%fill_supercell()
    call self%set_movers()
    call self%spin_mover%set_ncfile_name(self%params, self%filenames(2))
    call self%spin_mover%run_time(self%pots)
    call self%lattice_mover%run_time(self%pots)

    call self%spin_mover%spin_ncfile%close()
  end subroutine run_spin_latt_dynamics

  !-------------------------------------------------------------------!
  ! Run coupled lattice spin dynamics
  ! TODO: This is only a prototype. It does not have the proper logic
  ! to decide the time step, etc.
  ! TODO: move this to somewhere else, perhaps a spin-lattice mover file.
  !-------------------------------------------------------------------!
  subroutine run_coupled_spin_latt_dynamics(self)
    class(mb_manager_t), intent(inout) :: self
    integer :: istep
    character(len=90) :: msg

    call self%prim_pots%initialize()
    call self%read_potentials()
    
    call self%sc_maker%initialize(diag(self%params%ncell))
    call self%fill_supercell()
    call self%set_movers()

    call self%spin_mover%set_ncfile_name(self%params, self%filenames(2))
    ! use
    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')
    do istep = 1 , self%params%ntime
       call self%lattice_mover%run_one_step(self%pots, spin=self%spin_mover%Stmp)
       write(msg, "(A13, 4X,  I13)")  "Latt_Iter", istep
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       call self%spin_mover%run_one_step(self%pots, displacement=self%lattice_mover%displacement)
       write(msg, "(A13, 4X,  I13)")  "Spin_Iter", istep
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')
    end do
    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')

    call self%spin_mover%spin_ncfile%close()
  end subroutine run_coupled_spin_latt_dynamics


  !-------------------------------------------------------------------!
  ! Run all jobs
  !-------------------------------------------------------------------!
  subroutine run(self)
    class(mb_manager_t), intent(inout) :: self
    ! if ... fit lattice model
    ! if ... fit lwf model
    ! if ... run dynamics...
    if(self%params%spin_dynamics>0 .and. self%params%dynamics<=0) then
       if (self%params%spin_var_temperature==0) then
          call self%run_spin_dynamics()
       elseif (self%params%spin_var_temperature==1) then
          call self%run_MvT()
       end if
    else if (self%params%dynamics>0 .and. self%params%spin_dynamics<=0) then
       call self%run_lattice_dynamics()

    else if (self%params%dynamics>0 .and. self%params%spin_dynamics>0) then
       !call self%run_spin_latt_dynamics()
       call self%run_coupled_spin_latt_dynamics()
    end if

  end subroutine run



  !-------------------------------------------------------------------!
  !run_all: THE function which does everything
  !         from the very begining to end.
  !-------------------------------------------------------------------!
  subroutine run_all(self, filenames, params)
    class(mb_manager_t), intent(inout) :: self
    character(len=fnlen), intent(inout) :: filenames(17)
    type(multibinit_dtset_type), optional, intent(in) :: params
    call self%initialize(filenames, params=params)
    call self%run()
    call self%finalize()
  end subroutine run_all

end module m_multibinit_manager
