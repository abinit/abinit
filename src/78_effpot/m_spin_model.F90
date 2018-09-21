!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_model
!! NAME
!! m_spin_model
!!
!! FUNCTION
!! This module contains the manager of spin dynamics. It call other
!! related modules to run spin dynamics. Itself does nothing in detail.
!!
!!
!! Datatypes:
!!
!! * spin_model_t
!!
!! Subroutines:
!!
!! * spin_model_t
!! * spin_model_t_run
!! * spin_model_t_initialize
!! * spin_model_t_set_initial_spin
!! * spin_model_t_finalize
!! * spin_model_t_read_xml
!! * spin_model_t_make_supercell
!! * spin_model_t_run_one_step
!! * spin_model_t_run_time
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (hexu)
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

module m_spin_model

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_spin_terms, only: spin_terms_t, spin_terms_t_finalize, spin_terms_t_set_external_hfield
  use m_spin_model_primitive, only: spin_model_primitive_t, &
       & spin_model_primitive_t_initialize, &
       & spin_model_primitive_t_print_terms, &
       & spin_model_primitive_t_finalize, &
       & spin_model_primitive_t_read_xml, &
       & spin_model_primitive_t_make_supercell
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_init, spin_hist_t_get_s, spin_hist_t_free, &
       & spin_hist_t_set_params
  use m_spin_mover, only: spin_mover_t, spin_mover_t_initialize, spin_mover_t_finalize, &
       & spin_mover_t_run_time, spin_mover_t_run_one_step
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_init, spin_ncfile_t_close, spin_ncfile_t_def_sd, &
       & spin_ncfile_t_write_primitive_cell, spin_ncfile_t_write_supercell, spin_ncfile_t_write_parameters, &
       & spin_ncfile_t_write_one_step
  implicit none
!!***


  !!****t* m_spin_model/spin_model_t
  !! NAME
  !! spin_model_t
  !!
  !! FUNCTION
  !! This type contains all the data for spin dynamics.
  !!
  !! It contains:
  !! spin_primitive = information of the primitive cell, which is a map to the xml file.
  !! spin_calculator = the calculator for spin dynamics, include information of supercell.
  !! spin_mover = include information for how to run spin dynamics
  !! params = parameters from input file
  !! spin_ncfile = wrapper for spin hist netcdf file.
  !! nmatom= number of magnetic atoms
  !! SOURCE

  type spin_model_t
     type(spin_model_primitive_t) :: spin_primitive
     type(spin_terms_t) :: spin_calculator
     type(spin_hist_t):: spin_hist
     type(spin_mover_t):: spin_mover
     type(multibinit_dtset_type) :: params
     type(spin_ncfile_t) :: spin_ncfile
     integer :: nmatoms
     !  CONTAINS
     !    procedure :: run => spin_model_t_run
     !    procedure :: initialize=>spin_model_t_initialize
     !    procedure :: set_initial_spin => spin_model_t_set_initial_spin
     !    procedure :: finalize => spin_model_t_finalize
     !    procedure :: read_xml => spin_model_t_read_xml
     !    procedure :: make_supercell => spin_model_t_make_supercell
     !    procedure :: run_one_step => spin_model_t_run_one_step
     !    procedure :: run_time => spin_model_t_run_time
     !    procedure :: run_MvT => spin_model_t_run_MvT
     !    procedure :: run_MvH => spin_model_t_run_MvH
  end type spin_model_t
  !!***
contains

  !!****f* m_spin_model/spin_model_t_run
  !!
  !! NAME
  !!  spin_model_t_run
  !!
  !! FUNCTION
  !!  run the whole spin dynamics.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  currently it just call run_time.
  !!  It will be able to call e.g. run_MvT to do other things.
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    integer :: i
    !call spin_mover_t_run_time(self%spin_mover, self%spin_calculator, self%spin_hist)
    call spin_model_t_run_time(self)
  end subroutine spin_model_t_run
  !!***

  !!****f* m_spin_model/spin_model_t_initialize
  !!
  !! NAME
  !!  spin_model_t_initialize
  !!
  !! FUNCTION
  !! initialize spin dynamics from input
  !!
  !! INPUTS
  !!  xml_fname= file of xml file
  !! params = multibinit dataset from input file
  !! OUTPUT
  !!  self <spin_model_t>
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_initialize(self, xml_fname,  params)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_initialize'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    character(len=*), intent(in) :: xml_fname
    integer :: sc_matrix(3,3)
    type(multibinit_dtset_type), intent(in) :: params

    ! read input
    !call self%spin_primitive%initialize()
    call spin_model_primitive_t_initialize(self%spin_primitive)
    self%params=params
    !call self%read_xml(xml_fname)
    call spin_model_t_read_xml(self, xml_fname)
    !call self%spin_primitive%print_terms()
    call spin_model_primitive_t_print_terms(self%spin_primitive)

    ! make supercell
    sc_matrix(:,:)=0.0_dp
    sc_matrix(1,1)=params%ncell(1)
    sc_matrix(2,2)=params%ncell(2)
    sc_matrix(3,3)=params%ncell(3)
    !call self%make_supercell(sc_matrix)
    call spin_model_t_make_supercell(self, sc_matrix)

    ! set parameters to hamiltonian and mover
    self%nmatoms= self%spin_calculator%nmatoms

    call spin_model_t_set_params(self)


    !TODO hexu: mxhist, has_latt, natom should be input with their true values when lattice part also added
    call spin_hist_t_init(hist=self%spin_hist, nmatom=self%nmatoms, mxhist=3, has_latt=.False.)
    call spin_hist_t_set_params(self%spin_hist, spin_nctime=self%params%spin_nctime, &
            &     spin_temperature=self%params%spin_temperature)
    !TODO
    !call spin_hist_t_set_atomic_structure(self%spin_hist, acell, rprimd, xred, spin_index, ntypat,  typat, znucl)

    !call self%set_initial_spin(mode=1)
    call spin_model_t_set_initial_spin(self, mode=0)

    !call self%spin_mover%initialize(self%nmatoms, dt=params%dtspin, total_time=params%dtspin*params%ntime_spin, temperature=self%params%self)
    call spin_mover_t_initialize(self%spin_mover, self%nmatoms, dt=params%spin_dt, &
         &  total_time=params%spin_dt*params%spin_ntime, temperature=self%params%spin_temperature)

    call spin_ncfile_t_init(self%spin_ncfile, 'spinhist.nc')
    call spin_ncfile_t_def_sd(self%spin_ncfile, self%spin_hist )
    !call spin_ncfile_t_write_primitive_cell(self%spin_ncfile, self%spin_primitive)
    call spin_ncfile_t_write_supercell(self%spin_ncfile, self%spin_calculator)
    call spin_ncfile_t_write_parameters(self%spin_ncfile, self%params)

    call spin_ncfile_t_write_one_step(self%spin_ncfile, self%spin_hist)
  end subroutine spin_model_t_initialize
  !!***

  !!****f* m_spin_model/spin_model_t_finalize
  !!
  !! NAME
  !!  spin_model_t_finalize
  !!
  !! FUNCTION
  !! clean memory for spin dynamics
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_finalize'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    !call self%spin_primitive%finalize()
    call spin_model_primitive_t_finalize(self%spin_primitive)
    call spin_terms_t_finalize(self%spin_calculator)
    call spin_hist_t_free(self%spin_hist)
    call spin_mover_t_finalize(self%spin_mover)
    call spin_ncfile_t_close(self%spin_ncfile)
    ! TODO: finalize others
  end subroutine spin_model_t_finalize
  !!***

  !!****f* m_spin_model/spin_model_t_set_params
  !!
  !! NAME
  !!  spin_model_t_set_params
  !!
  !! FUNCTION
  !! set parameters for spin dynamics from input.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !! Most are already done in initialize. (should be moved to here or remove this function)
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_set_params(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_set_params'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    real(dp):: mfield(3, self%nmatoms)
    integer ::  i
    ! params -> mover
    ! params -> calculator
    do i=1, self%nmatoms
      mfield(:, i)=self%params%spin_mag_field(:)
    end do
    call spin_terms_t_set_external_hfield(self%spin_calculator, mfield)
    ! params -> hist

  end subroutine spin_model_t_set_params
  !!***

  !!****f* m_spin_model/spin_model_t_read_xml
  !!
  !! NAME
  !!  spin_model_t_read_xml
  !!
  !! FUNCTION
  !! read xml file and set primitive cell parameters.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_read_xml(self, xml_fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_read_xml'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    character(len=*), intent(in) :: xml_fname
    !call self%spin_primitive%read_xml(xml_fname)
    call spin_model_primitive_t_read_xml(self%spin_primitive, xml_fname)
  end subroutine spin_model_t_read_xml
  !!***


  !!****f* m_spin_model/spin_model_t_make_supercell
  !!
  !! NAME
  !!  spin_model_t_make_supercell
  !!
  !! FUNCTION
  !!  make supercell and prepare calculator
  !!
  !! INPUTS
  !!  sc_mat(3,3) = supercell matrix
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_make_supercell(self, sc_mat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_make_supercell'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    integer , intent(in):: sc_mat(3, 3)
    !call self%spin_primitive%make_supercell(sc_mat, self%spin_calculator)
    call spin_model_primitive_t_make_supercell(self%spin_primitive, sc_mat, self%spin_calculator)
  end subroutine spin_model_t_make_supercell
  !!***

  !!****f* m_spin_model/spin_model_t_set_initial_spin
  !!
  !! NAME
  !!  spin_model_t_set_initial_spin
  !!
  !! FUNCTION
  !!  set initial spin state
  !!
  !! INPUTS
  !!  mode= 0 : all along z. 1: random
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_set_initial_spin(self, mode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_set_initial_spin'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    integer, intent(in) :: mode
    integer :: i
    real(dp) :: S(3, self%nmatoms)
    character(len=500) :: msg
    if(mode==0) then
       ! set all spin to z direction.
       S(1,:)=0.0d0
       S(2,:)=0.0d0
       S(3,:)=1.0d0
    else if (mode==1) then
       ! randomize S using uniform random number
       ! print *, "Initial spin set to random value"
   write(msg,*) "Initial spin set to random value."
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
       call random_number(S)
       S=S-0.5
       do i=1, self%nmatoms
          S(:,i)=S(:,i)/sqrt(sum(S(:, i)**2))
       end do
    else
      write(msg,*) "Error: Set initial spin: mode should be 0 (FM) or 1 (random)"
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')

    end if

    call spin_hist_t_set_vars(self%spin_hist, S=S, time=0.0_dp, ihist_latt=0, inc=.True.)
    !print *, "initial spin", self%spin_hist%current_S
  end subroutine spin_model_t_set_initial_spin
  !!***


  !!****f* m_spin_model/spin_model_t_run_one_step
  !!
  !! NAME
  !!  spin_model_t_run_one_step
  !!
  !! FUNCTION
  !!  run one spin dynamics step
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  This is not called by run_time. insted run_time call mover%run_time
  !!  (in the spirit that this module should do nothing in detail!)
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run_one_step(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_one_step'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    real(dp) :: S_tmp(3,self%nmatoms), etot
    call spin_mover_t_run_one_step(self%spin_mover, self%spin_calculator, &
         spin_hist_t_get_S(self%spin_hist),S_tmp, etot)
    call spin_hist_t_set_vars(self%spin_hist, S=S_tmp, etot=etot, inc=.False.)
  end subroutine spin_model_t_run_one_step
  !!***

  !!****f* m_spin_model/spin_model_t_run_time
  !!
  !! NAME
  !!  spin_model_t_run_time
  !!
  !! FUNCTION
  !!  run all spin time step
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run_time(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_time'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    call spin_mover_t_run_time(self%spin_mover,self%spin_calculator,self%spin_hist, self%spin_ncfile)
  end subroutine spin_model_t_run_time
  !!***


  !!****f* m_spin_model/spin_model_t_run_MvT
  !!
  !! NAME
  !!  spin_model_t_run_MvT
  !!
  !! FUNCTION
  !!  run spin vs temperature
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  Not yet implemented.
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run_MvT(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_MvT'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    !TODO
    write(std_err, *) "MvH is not yet implemented. "
  end subroutine spin_model_t_run_MvT
  !!***

  !!****f* m_spin_model/spin_model_t_run_MvH
  !!
  !! NAME
  !!  spin_model_t_run_MvT
  !!
  !! FUNCTION
  !!  run spin vs external magnetic field
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  Not yet implemented.
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run_MvH(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_MvH'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    write(std_err, *) "MvH is not yet implemented. "
  end subroutine spin_model_t_run_MvH
  !!***

end module m_spin_model
