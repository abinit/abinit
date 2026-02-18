!!****p* ABINIT/atdep
!! NAME
!! atdep
!!
!! FUNCTION
!! Calculations of phonons using molecular dynamic simulations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2026 ABINIT group (FB,JB,GA)
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

program atdep

  use defs_basis
  use m_abicore
  use m_phonons
  use m_errors
  use m_abi_linalg
  use m_xmpi
  use m_abihist
  use m_io_tools
  use m_argparse

  use m_time,             only : asctime, timein, timab
  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_tdep_abitypes,    only : Qbz_type, tdep_init_crystal, tdep_init_ifc, tdep_init_ddb, tdep_write_ddb, &
&                                tdep_destroy_qbz, tdep_ifc2phi2, tdep_read_ifc, tdep_write_ifc
  use m_tdep_latt,        only : Lattice_type, tdep_make_latt
  use m_tdep_sym,         only : tdep_make_sym, Symmetries_type, tdep_destroy_sym
  use m_tdep_dataset,     only : tdep_read_input, tdep_init_MPIdata, &
&                                tdep_destroy_mpidata, atdep_dataset_type, MPI_enreg_type, tdep_destroy_invar, version_string
  use m_tdep_qpt,         only : tdep_make_qptpath, Qpoints_type, tdep_destroy_qpt
  use m_tdep_sampling,    only : tdep_Sampling_type, tdep_sampling_init_read, tdep_sampling_free,&
&                                tdep_sampling_shift_xred, tdep_sampling_rotate, tdep_MatchIdeal2Average
  use m_tdep_shell,       only : Shell_type, tdep_init_shell2at, tdep_init_shell3at, tdep_init_shell4at, &
&                                tdep_init_shell1at, tdep_destroy_shell
  use m_tdep_solver,      only : tdep_Solver_type, tdep_solver_init, tdep_solver_free, tdep_solver_set_residual_forces, &
&                                tdep_calc_phi1fcoeff, tdep_calc_phi2fcoeff, tdep_calc_phi3fcoeff, tdep_calc_phi4fcoeff, &
&                                tdep_calc_MoorePenrose, tdep_calc_constraints
  use m_tdep_model,       only : tdep_Model_type, tdep_model_init, tdep_model_free, tdep_calc_model,&
&                                Phi2_type, tdep_init_phi2, tdep_destroy_phi2
  use m_tdep_phi2,        only : tdep_calc_phi2, tdep_write_phi2, tdep_calc_ftot2, &
&                                Eigen_type, tdep_init_eigen2nd, tdep_destroy_eigen2nd, tdep_calc_phi1, tdep_write_phi1
  use m_tdep_phi3,        only : tdep_calc_phi3ref, tdep_write_phi3, tdep_calc_ftot3, &
&                                tdep_calc_alpha_gamma, tdep_write_gruneisen
  use m_tdep_phi4,        only : tdep_calc_phi4ref, tdep_write_phi4, tdep_calc_ftot4
  use m_tdep_phdos,       only : tdep_calc_phdos,tdep_calc_elastic,tdep_calc_thermo
  use m_tdep_utils,       only : tdep_check_constraints, tdep_print_Aknowledgments

  implicit none

  integer :: print_mem_report
  integer :: stdout,stdlog
  integer :: iorder
  real(dp) :: rotation(3,3)
  real(dp) :: tcpu, tcpui, twall, twalli
  real(dp) :: tsec(2)
  character(len = 24):: start_datetime
  type(args_t) :: args
  type(atdep_dataset_type) :: Invar
  type(MPI_enreg_type) :: MPIdata
  type(abihist) :: Hist
  type(ifc_type) :: Ifc
  type(ddb_type) :: DDB
  type(crystal_t) :: Crystal
  type(Lattice_type) :: Lattice
  type(Symmetries_type) :: Sym
  type(Qpoints_type) :: Qpt
  type(Qbz_type) :: Qbz
  type(Shell_type) :: Shell1at, Shell2at, Shell3at, Shell4at
  type(tdep_Sampling_type) :: MD
  type(tdep_Solver_type) :: Solver
  type(tdep_Model_type) :: Model
  type(Eigen_type) :: Eigen2nd_MP, Eigen2nd_path
  type(phdos_t) :: PHdos

!******************************************************************

!==========================================================================================
!===================== Initialization & Reading  ==========================================
!==========================================================================================
! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)
! Initialize MPI
 call xmpi_init()

 ! Initialisation of the timing
 call timein(tcpui, twalli)
 start_datetime = asctime()
 call timab(1, 0, tsec)

! Parse command line arguments.
 args = args_parser(); if (args%exit /= 0) goto 100

! Initialize memory profiling if activated at configure time.
! if a full report is desired, set the argument of abimem_init to "2" instead of "0" via the command line.
! note that the file can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(args%abimem_level, limit_mb=args%abimem_limit_mb)
#endif

! Read input values from the input.in input file
 call tdep_read_input(args%input_path,Hist,Invar)
 call tdep_init_MPIdata(Invar,MPIdata)
 call tdep_sampling_init_read(MD,Invar,MPIdata,Hist)
 call abihist_free(Hist)

 if (args%dry_run /= 0) then
   call wrtout(std_out, "Dry run mode. Exiting after have read the input")
   call tdep_sampling_free(MD)
   call tdep_destroy_invar(Invar)
   call tdep_destroy_mpidata(MPIdata)
   goto 100
 end if

! Initialize basic quantities
 print_mem_report = 1
 stdout           = Invar%stdout
 stdlog           = Invar%stdlog

!==========================================================================================
!============== Define the ideal lattice, symmetries and Brillouin zone ===================
!==========================================================================================

!Define all the quantities needed to buid the lattice (rprim*, acell*, brav*...)
 call tdep_make_latt(Invar,Lattice,rotation)

!Compute all the symmetries coming from the bravais lattice
 call tdep_make_sym(Invar,Lattice,MPIdata,Sym)

!Initialize the Brillouin zone and compute the q-points path
 call tdep_make_qptpath(Invar,Lattice,MPIdata,Qpt)

!==========================================================================================
!============== Complete the initialization of the Sampling ===============================
!==========================================================================================

!Shift xred to keep atoms in the same unit cell at each step.
 call tdep_sampling_shift_xred(MD, MPIdata)

!Apply rotation to cartesian forces
 call tdep_sampling_rotate(MD, rotation)

!Map the atoms of the input moledular dynamics to the ideal supercell
 call tdep_MatchIdeal2Average(MD,Invar,Lattice,Sym,MPIdata)

!==========================================================================================
!============== Initialize the Shell1at datatype ==========================================
!==========================================================================================
 call tdep_init_shell1at(Shell1at,Invar,MD,Sym,MPIdata)

!==========================================================================================
!============== Initialize the Shell2at datatype ==========================================
!==========================================================================================
 call tdep_init_shell2at(Shell2at,Invar,MD,Sym,MPIdata)

!==========================================================================================
!============== Initialize the Shell3at datatype ==========================================
!==========================================================================================
 if (Invar%order.ge.3) then
   call tdep_init_shell3at(Shell3at,Invar,MD,Sym,MPIdata)
 end if

!==========================================================================================
!============== Initialize the Shell4at datatype ==========================================
!==========================================================================================
 if (Invar%order==4) then
   call tdep_init_shell4at(Shell4at,Invar,MD,Sym,MPIdata)
 end if

!==========================================================================================
!============== Initialize the TDEP Model datatype ========================================
!==========================================================================================
 call tdep_model_init(Model, Invar, Shell3at, Shell4at)

!==========================================================================================
!============== Initialize the TDEP Solver datatype =======================================
!==========================================================================================
 call tdep_solver_init(Solver, Invar, Shell1at, Shell2at, Shell3at, Shell4at)

!==========================================================================================
!============== Initialize Crystal, DDB, and IFC ABINIT Datatypes =========================
!==========================================================================================
 call tdep_init_crystal(Crystal,Invar,Lattice,Sym)
 call tdep_init_ddb(Crystal,DDB,Invar,Lattice,MPIdata,Qbz)
 call tdep_init_ifc(Crystal,DDB,Ifc,Invar,Lattice,MPIdata,Model%Phi2,MD%Rlatt_scaled,Shell2at,Sym)

!==========================================================================================
!================= Copy the forces into the solver ========================================
!==========================================================================================
 !Remove the supercell contribution (the "LR part") included in total forces
 !before computing the "SR part". The full "LR part" will be added later.
 ! GA: FIXME
 !     Since Ifc%ewald_atfrc has not been initialized, Model%Forces are zero!
 if (Invar%loto) then
   call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,Invar%natom_unitcell,1,Model%Phi2,MD%Rlatt_scaled,Shell2at,Sym)
   call tdep_calc_ftot2(Model,Invar,Model%Phi2%Tot,MD%ucart)
 end if

 call tdep_solver_set_residual_forces(Solver, MD, Model)

 ! GA: I dont think this is even necessary
 !if (Invar%loto) then
 !  Model%Phi1(:) = zero
 !  Model%Phi1Ui(:) = zero
 !  Model%Phi2%SR(:,:) = zero
 !  Model%Phi2%LR(:,:) = zero
 !  Model%Phi2%Tot(:,:) = zero
 !  Model%Phi2UiUj(:) = zero
 !  Model%Forces(:) = zero
 !end if

!==========================================================================================
!================= Build fcoeff and compute constraints ===================================
!==========================================================================================

 if (Invar%readifc.ne.1) then
   call tdep_calc_phi1fcoeff(Solver,Invar,Shell1at,Sym)
   call tdep_calc_phi2fcoeff(Solver,Invar,Shell2at,Sym,MD)
 end if

 if (Invar%order.ge.3) then
   call tdep_calc_phi3fcoeff(Solver,Invar,Shell3at,Sym,MD)
 end if

 if (Invar%order.eq.4) then
   call tdep_calc_phi4fcoeff(Solver,Invar,Shell4at,Sym,MD)
 end if

 call tdep_calc_constraints(Solver,MD%distance,Invar,MPIdata,Sym,&
&                           Shell1at,Shell2at,Shell3at,Shell4at)


!==========================================================================================
!============= Compute the pseudo inverse using the Moore-Penrose method ==================
!==========================================================================================

!=================== If all the Orders are solved simultaneously ==========================

 if (Invar%together.eq.1) then
   write(stdout,*) '############### (Solve simultaneously all the orders) #######################'

   if (Invar%readifc.ne.1) then
     call tdep_calc_MoorePenrose(Solver,0,Invar,MPIdata)
     call tdep_calc_phi1(Solver,Shell1at,Sym,Model%Phi1)
     call tdep_calc_phi2(Solver,Shell2at,Sym,Model%Phi2%SR)
   end if
   call tdep_calc_ftot2(Model,Invar,Model%Phi2%SR,MD%ucart)

   if (Invar%order.ge.3) then
     call tdep_calc_phi3ref(Solver,Shell3at,Model%Phi3)
     call tdep_calc_ftot3(Model,Invar,Shell3at,MD%ucart,Sym)
   end if

   if (Invar%order.ge.4) then
     call tdep_calc_phi4ref(Solver,Shell4at,Model%Phi4)
     call tdep_calc_ftot4(Model,Invar,Shell4at,MD%ucart,Sym)
   end if

!=================== If all the Orders are solved successively ============================

!ATTENTION : Le LR semble enleve a l'ordre 2 mais pas a l'ordre 3 et 4.
!            On repart de MD%Forces tout en bas et pas de Solver%Forces
!            (car les ordres 2 et 3 sont supprimes)
 else if (Invar%together.eq.0) then
   write(stdout,*) '################## (Solve successively each order) ##########################'
   do iorder=1,Invar%order-1
     write(stdout,*) ' For order=',iorder+1
     if (Invar%readifc.eq.1) cycle

     call tdep_calc_MoorePenrose(Solver,iorder,Invar,MPIdata)

     if (iorder.eq.1) then
       if (Invar%readifc.ne.1) then
         call tdep_calc_phi1(Solver,Shell1at,Sym,Model%Phi1)
         call tdep_calc_phi2(Solver,Shell2at,Sym,Model%Phi2%SR)
       end if
       call tdep_calc_ftot2(Model,Invar,Model%Phi2%SR,MD%ucart)

     else if (iorder.eq.2) then
       call tdep_calc_phi3ref(Solver,Shell3at,Model%Phi3)
       call tdep_calc_ftot3(Model,Invar,Shell3at,MD%ucart,Sym)

     else if (iorder.eq.3) then
       call tdep_calc_phi4ref(Solver,Shell4at,Model%Phi4)
       call tdep_calc_ftot4(Model,Invar,Shell4at,MD%ucart,Sym)
     end if

     ! Remove forces computed with the previous order
     call tdep_solver_set_residual_forces(Solver, MD, Model)

   end do ! iorder
 end if

! Add the long-range part of the IFC
 if (Invar%loto) then
   Model%Phi2%Tot = Model%Phi2%LR + Model%Phi2%SR
 end if

! Free some memory
 call tdep_solver_free(Solver)
 call tdep_destroy_shell(Shell1at)

!==========================================================================================
!=================== Write the IFC and check the constraints ==============================
!==========================================================================================
 call tdep_write_phi1(Invar,Model%Phi1)
 call tdep_write_phi2(MD%distance,Invar,MPIdata,Model%Phi2%SR,Shell2at)
 if (Invar%order.ge.3) then
   call tdep_write_phi3(MD%distance,Invar,Model%Phi3,Shell3at,Sym)
 end if
 if (Invar%order.ge.4) then
   call tdep_write_phi4(MD%distance,Invar,Model%Phi4,Shell4at,Sym)
 end if

 call tdep_check_constraints(Model,MD%distance,Invar,Sym,Shell3at,Shell4at)

!==========================================================================================
!===================== Convert Phi2 into IFC object =======================================
!==========================================================================================

 call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,Invar%natom_unitcell,0,&
&                   Model%Phi2,MD%Rlatt_scaled,Shell2at,Sym)

!==========================================================================================
!===================== Compute the phonon spectrum, the DOS, ==============================
!=====================  the dynamical matrix and write them ===============================
!==========================================================================================
 call tdep_init_eigen2nd(Eigen2nd_MP,Invar%natom_unitcell,Qbz%nqbz)
 call tdep_init_eigen2nd(Eigen2nd_path,Invar%natom_unitcell,Qpt%nqpt)

 call tdep_calc_phdos(Crystal,DDB,Eigen2nd_MP,Eigen2nd_path,Ifc,Invar,Lattice,MPIdata,Invar%natom,&
&                     Invar%natom_unitcell,Model%Phi2,PHdos,Qbz,Qpt,MD%Rlatt_scaled,Shell2at,Sym)
 call tdep_destroy_shell(Shell2at)


 ! Create a new DDB with the coarse q-point grid in the IBZ.
 call DDB%free()
 call Ifc%to_ddb(DDB,Crystal)
 call tdep_write_ddb(DDB,Crystal,Invar)

 write(stdout,'(a)') ' See the dij.dat, omega.dat and eigenvectors files'
 write(stdout,'(a)') ' See also the DDB file'

!==========================================================================================
!===================== Compute the elastic constants ======================================
!==========================================================================================
 call tdep_calc_elastic(Model%Phi2%SR,MD%distance,Invar,Lattice)

!==========================================================================================
!=========== Compute U_0, the "free energy" from the model ===============
!==========================================================================================
 call tdep_calc_model(Model,MD,Invar,MPIdata)

!==========================================================================================
!===================== Compute the thermodynamical quantities =============================
!==========================================================================================
 call tdep_calc_thermo(Invar,Lattice,MPIdata,PHdos,Model%U0)
 call PHdos%free()

!==========================================================================================
!===================== CALCULATION OF THE 3rd ORDER =======================================
!==========================================================================================
 if (Invar%order>2) then
   if (MPIdata%iam_master) then
     call tdep_write_gruneisen(MD%distance,Eigen2nd_path,Invar,Model%Phi3,Qpt,MD%Rlatt_cart,Shell3at,Sym)
   end if
   call tdep_calc_alpha_gamma(MD%distance,Eigen2nd_MP,Invar,Lattice,MPIdata,Model%Phi3,Qbz,MD%Rlatt_cart,Shell3at,Sym)

!FB call tdep_calc_lifetime1(Crystal,MD%distance,Eigen2nd_MP,Ifc,Invar,Lattice,Model%Phi3,Qbz,MD%Rlatt_cart,Shell3at,Sym)
 end if

!==========================================================================================
!===================== Free memore ========================================================
!==========================================================================================

 call tdep_destroy_eigen2nd(Eigen2nd_path)
 call tdep_destroy_eigen2nd(Eigen2nd_MP)

 if (Invar%order>2) then
   call tdep_destroy_shell(Shell3at)
   if (Invar%order.eq.4) then
     call tdep_destroy_shell(Shell4at)
   end if
 end if

 call Ifc%free()
 call DDB%free()
 call Crystal%free()
 call tdep_destroy_sym(Sym)
 call tdep_destroy_qbz(Qbz)
 call tdep_destroy_qpt(Qpt)
 call tdep_model_free(Model)
 call tdep_sampling_free(MD)
 call tdep_destroy_invar(Invar)
 call tdep_destroy_mpidata(MPIdata)

!==========================================================================================
!===================== End the calculation ================================================
!==========================================================================================

 call tdep_print_Aknowledgments(stdout)
 call flush_unit(stdout)

 call timein(tcpu, twall)
 tsec(1)=tcpu-tcpui; tsec(2)=twall-twalli


! Write YAML document with the final summary.
 if (MPIdata%iam_master) then
   write(stdlog, "(a)")""
   write(stdlog, "(a)")"--- !FinalSummary"
   write(stdlog, "(a)")"program: atdep"
   write(stdlog, "(2a)")"version: ", trim(version_string)
   write(stdlog, "(2a)")"start_datetime: ", start_datetime
   write(stdlog, "(2a)")"end_datetime: ", asctime()
   write(stdlog, "(a, f13.1)")"overall_cpu_time: ", tsec(1)
   write(stdlog, "(a, f13.1)")"overall_wall_time: ", tsec(2)
   write(stdlog, "(a, i0)")"mpi_procs: ", MPIdata%nproc
   write(stdlog, "(a)")"..."
   call flush_unit(stdlog)
 end if

!Memory analysis
 call abinit_doctor(trim(Invar%output_prefix), print_mem_report=print_mem_report)
 call flush_unit(stdlog)
 close(unit=stdout)
100 call xmpi_end()

 end program atdep
!!***
