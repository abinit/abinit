!!****p* ABINIT/atdep
!! NAME
!! atdep
!!
!! FUNCTION
!! Calculations of phonons using molecular dynamic simulations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (FB,JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!.
!!
!! NOTES
!!  The input files are input.in, xred.dat, fcart.dat and etot.dat
!!  See the examples in the test directory
!!
!! TODO
!!  A lot of things
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_init,abinit_doctor,tdep_build_phijnn,tdep_build_pijn
!!      tdep_calc_alpha_gamma,tdep_calc_constraints,tdep_calc_dij
!!      tdep_calc_elastic,tdep_calc_model,tdep_calc_moorepenrose
!!      tdep_calc_phdos,tdep_calc_phijfcoeff,tdep_calc_pijfcoeff
!!      tdep_calc_psijfcoeff,tdep_calc_psijtot,tdep_calc_thermo
!!      tdep_check_constraints,tdep_destroy_eigen2nd,tdep_destroy_shell
!!      tdep_init_crystal,tdep_init_ddb,tdep_init_eigen2nd,tdep_init_ifc
!!      tdep_init_shell1at,tdep_init_shell2at,tdep_init_shell3at,tdep_make_latt
!!      tdep_make_qptpath,tdep_make_sym,tdep_matchideal2average
!!      tdep_print_aknowledgments,tdep_readecho,tdep_write_dij
!!      tdep_write_gruneisen,tdep_write_yaml,wrtout,xmpi_end,xmpi_init
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
  use m_xmpi
  use m_io_tools
  use m_errors
  use m_argparse

  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_tdep_abitypes,    only : tdep_init_crystal, tdep_init_ifc, tdep_init_ddb, tdep_write_ifc
  use m_tdep_psij,        only : tdep_calc_psijfcoeff, tdep_calc_psijtot, tdep_calc_alpha_gamma, tdep_write_gruneisen
  use m_tdep_phij,        only : tdep_calc_phijfcoeff, tdep_calc_pijfcoeff, tdep_build_phijNN, tdep_calc_dij, tdep_write_dij, &
                                 Eigen_Variables_type, tdep_init_eigen2nd, tdep_destroy_eigen2nd, &
                                 tdep_write_yaml, tdep_build_pijN
  use m_tdep_latt,        only : tdep_make_latt, Lattice_Variables_type
  use m_tdep_sym,         only : tdep_make_sym, Symetries_Variables_type
  use m_tdep_readwrite,   only : tdep_print_Aknowledgments, tdep_ReadEcho, Input_Variables_type
  use m_tdep_utils,       only : Coeff_Moore_type, tdep_calc_MoorePenrose, tdep_MatchIdeal2Average, tdep_calc_model
  use m_tdep_qpt,         only : tdep_make_qptpath, Qpoints_type
  use m_tdep_phdos,       only : tdep_calc_phdos,tdep_calc_elastic,tdep_calc_thermo
  use m_tdep_shell,       only : Shell_Variables_type, tdep_init_shell2at, tdep_init_shell3at, tdep_init_shell1at, &
                                 tdep_destroy_shell
  use m_tdep_constraints, only : tdep_calc_constraints, tdep_check_constraints

  implicit none

  integer :: natom,natom_unitcell,ncoeff1st,ncoeff2nd,ncoeff3rd,ntotcoeff,ntotconst
  integer :: stdout,iqpt,nshell_max
  double precision :: U0,Free_Anh
  double precision, allocatable :: ucart(:,:,:),proj1st(:,:,:),proj2nd(:,:,:),proj3rd(:,:,:)
  double precision, allocatable :: proj_tmp(:,:,:),Forces_TDEP(:)
!FB  double precision, allocatable :: fcoeff(:,:),Phij_coeff(:,:),Forces_MD(:),Phij_NN(:,:)
  double precision, allocatable :: Phij_coeff(:,:),Forces_MD(:),Phij_NN(:,:),Pij_N(:),Pij_coeff(:,:)
  double precision, allocatable :: Psij_coeff(:,:),Psij_ref(:,:,:,:),MP_coeff(:,:)
  double precision, allocatable :: distance(:,:,:),Rlatt_cart(:,:,:),Rlatt4Abi(:,:,:)
  double precision, allocatable :: omega (:),ftot3(:,:)
  double precision :: qpt_cart(3)
  double complex  , allocatable :: dij(:,:),eigenV(:,:)
  type(args_t) :: args
  type(phonon_dos_type) :: PHdos
  type(Input_Variables_type) :: InVar
  type(Lattice_Variables_type) :: Lattice
  type(Symetries_Variables_type) :: Sym
  type(Qpoints_type) :: Qpt
  type(ifc_type) :: Ifc
  type(ddb_type) :: DDB
  type(crystal_t) :: Crystal
  type(Shell_Variables_type) :: Shell1at,Shell2at,Shell3at
  type(Coeff_Moore_type) :: CoeffMoore
  type(Eigen_Variables_type) :: Eigen2nd

!******************************************************************

!==========================================================================================
!===================== Initialization & Reading  ==========================================
!==========================================================================================
 call xmpi_init()

 ! Parse command line arguments.
 args = args_parser(); if (args%exit /= 0) goto 100

 ! Initialize memory profiling if activated at configure time.
 ! if a full report is desired, set the argument of abimem_init to "2" instead of "0" via the command line.
 ! note that the file can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(args%abimem_level, limit_mb=args%abimem_limit_mb)
#endif

! Read input values from the input.in input file
 call tdep_ReadEcho(InVar)

 if (args%dry_run /= 0) then
   call wrtout(std_out, "Dry run mode. Exiting after have read the input")
   goto 100
 end if

! Initialize basic quantities
 natom         =InVar%natom
 natom_unitcell=InVar%natom_unitcell
 stdout        =InVar%stdout
 nshell_max    =500

!==========================================================================================
!============== Define the ideal lattice, symmetries and Brillouin zone ===================
!==========================================================================================
! Define all the quantities needed to buid the lattice (rprim*, acell*, brav*...)
 call tdep_make_latt(InVar,Lattice)

! Compute all the symmetries coming from the bravais lattice
 call tdep_make_sym(Invar,Lattice,Sym)

! Initialize the Brillouin zone and compute the q-points path
 call tdep_make_qptpath(InVar,Lattice,Qpt)

!==========================================================================================
!======== 1/ Determine ideal positions and distances ======================================
!======== 2/ Find the matching between the ideal and average ==============================
!========   (from the MD simulations) positions. ==========================================
!======== 3/ Find the symmetry operation between the reference and image bonds =============
!======== 4/ Write output quantities needed to visualize the neighbouring distances =======
!==========================================================================================
 ABI_MALLOC(Rlatt4Abi ,(3,natom_unitcell,natom))   ; Rlatt4Abi (:,:,:)=0.d0
 ABI_MALLOC(distance,(natom,natom,4))              ; distance(:,:,:)=0.d0
 ABI_MALLOC(Rlatt_cart,(3,natom_unitcell,natom))   ; Rlatt_cart(:,:,:)=0.d0
 ABI_MALLOC(ucart,(3,natom,InVar%nstep))           ; ucart(:,:,:)=0.d0
 ABI_MALLOC(Forces_MD,(3*natom*InVar%nstep))       ; Forces_MD(:)=0.d0

 call tdep_MatchIdeal2Average(distance,Forces_MD,InVar,Lattice,Rlatt_cart,Rlatt4Abi,Sym,ucart)

!==========================================================================================
!============== Initialize Crystal and DDB ABINIT Datatypes ===============================
!==========================================================================================
 call tdep_init_crystal(Crystal,InVar,Lattice,Sym)
 call tdep_init_ddb(DDB,InVar,Lattice)

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=# CALCULATION OF THE 2nd ORDER =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

!==========================================================================================
!============== Initialize the Shell1at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(3,3,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell1at(distance,InVar,3,nshell_max,ncoeff1st,1,proj_tmp,Shell1at,Sym)
 ABI_MALLOC(proj1st  ,(3,3,Shell1at%nshell)) ; proj1st(:,:,:)=0.d0
 proj1st = reshape (proj_tmp, (/ 3,3,Shell1at%nshell /))
 ABI_FREE(proj_tmp)
!Rotational invariances (1st order)
!    constraints = 3
 CoeffMoore%nconst_1st=3
 ntotconst=CoeffMoore%nconst_1st

!==========================================================================================
!============== Initialize the Shell2at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(9,9,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell2at(distance,InVar,9,nshell_max,ncoeff2nd,2,proj_tmp,Shell2at,Sym)
 ABI_MALLOC(proj2nd  ,(9,9,Shell2at%nshell)) ; proj2nd(:,:,:)=0.d0
 proj2nd = reshape (proj_tmp, (/ 9,9,Shell2at%nshell /))
 ABI_FREE(proj_tmp)
!Rotational invariances (2nd order) + Symetry of the Dynamical Matrix + Huang invariances
!    constraints = natom*3**2 + (3*natom_unitcell)**2 + 3**4
 CoeffMoore%nconst_rot2nd = natom_unitcell*9
 CoeffMoore%nconst_dynmat = (3*natom_unitcell)**2
 CoeffMoore%nconst_huang  = 81
 CoeffMoore%nconst_2nd = CoeffMoore%nconst_rot2nd + CoeffMoore%nconst_dynmat + CoeffMoore%nconst_huang
 ntotconst=CoeffMoore%nconst_1st + CoeffMoore%nconst_2nd

!==========================================================================================
!============== Initialize the IFC Abinit datatype ========================================
!==========================================================================================
 ABI_MALLOC(Pij_N  ,(3*natom))         ; Pij_N  (:)  =0.d0
 ABI_MALLOC(Phij_NN,(3*natom,3*natom)) ; Phij_NN(:,:)=0.d0
 call tdep_init_ifc(Crystal,DDB,Ifc,InVar,Lattice,Phij_NN,Rlatt4Abi,Shell2at,Sym)

!==========================================================================================
!============== Initialize the Shell3at datatype (if SC_order==1) =========================
!==========================================================================================
 ncoeff3rd=0
 if (InVar%Order==3) then
   ABI_MALLOC(proj_tmp,(27,27,nshell_max)) ; proj_tmp(:,:,:)=0.d0
   call tdep_init_shell3at(distance,InVar,27,nshell_max,ncoeff3rd,3,proj_tmp,Shell3at,Sym)
   ABI_MALLOC(proj3rd  ,(27,27,Shell3at%nshell)) ; proj3rd(:,:,:)=0.d0
   proj3rd = reshape (proj_tmp, (/ 27,27,Shell3at%nshell /))
   ABI_FREE(proj_tmp)
!  Rotational invariances (3rd order) + acoustic sum rules (3rd order)
!    constraints = natom_unitcell*natom*3**3) + 3permutations*3**3*nshell2at
!FB   CoeffMoore%nconst_rot3rd =   natom_unitcell*natom*3**3
   CoeffMoore%nconst_rot3rd = 3*natom_unitcell*natom*3**4
   CoeffMoore%nconst_asr3rd = 8*natom_unitcell*natom*3**3
   CoeffMoore%nconst_3rd = CoeffMoore%nconst_rot3rd + CoeffMoore%nconst_asr3rd
   ntotconst=CoeffMoore%nconst_1st + CoeffMoore%nconst_2nd + CoeffMoore%nconst_3rd
 end if
 ntotcoeff=ncoeff1st+ncoeff2nd+ncoeff3rd
 CoeffMoore%ntotcoeff=ntotcoeff
 CoeffMoore%ntotconst=ntotconst
 CoeffMoore%ncoeff1st=ncoeff1st
 CoeffMoore%ncoeff2nd=ncoeff2nd
 CoeffMoore%ncoeff3rd=ncoeff3rd
 ABI_MALLOC(CoeffMoore%fcoeff,(3*natom*InVar%nstep+ntotconst,ntotcoeff)); CoeffMoore%fcoeff(:,:)=0.d0

!==========================================================================================
!============= Build fcoeff, needed for the Moore-Penrose method just below ===============
!============================== and compute constraints ===================================
!==========================================================================================
 if (InVar%ReadIFC.ne.1) then
   call tdep_calc_pijfcoeff(CoeffMoore,InVar,proj1st,Shell1at,Sym)
   call tdep_calc_phijfcoeff(CoeffMoore,InVar,proj2nd,Shell2at,Sym,ucart)
 end if
 if (InVar%Order==3) then
   call tdep_calc_psijfcoeff(CoeffMoore,InVar,proj3rd,Shell3at,Sym,ucart)
   CoeffMoore%fcoeff(:,ncoeff2nd+1:ntotcoeff)=CoeffMoore%fcoeff(:,ncoeff2nd+1:ntotcoeff)/2.d0
   call tdep_calc_constraints(CoeffMoore,distance,InVar,Shell1at%nshell,Shell2at%nshell,Shell3at%nshell,Sym,&
&                             proj1st,Shell1at,proj2nd,Shell2at,proj3rd,Shell3at)
 else
   call tdep_calc_constraints(CoeffMoore,distance,InVar,Shell1at%nshell,Shell2at%nshell,Shell3at%nshell,Sym,&
&                             proj1st,Shell1at,proj2nd,Shell2at)
 end if

!==========================================================================================
!============= Compute the pseudo inverse using the Moore-Penrose method ==================
!==========================================================================================
 ABI_MALLOC(MP_coeff,(ntotcoeff,1)); MP_coeff(:,:)=0.d0
 if (InVar%ReadIFC.ne.1) then
   call tdep_calc_MoorePenrose(Forces_MD,CoeffMoore,InVar,MP_coeff)
 end if
 ABI_MALLOC(Pij_coeff ,(ncoeff1st,1)); Pij_coeff (:,:)=MP_coeff(1:ncoeff1st,:)
 ABI_MALLOC(Phij_coeff,(ncoeff2nd,1)); Phij_coeff(:,:)=MP_coeff(ncoeff1st+1:ncoeff1st+ncoeff2nd,:)

!==========================================================================================
!==== Reorganize the IFC coefficients into the whole Pij_N & Phij_NN matrices =============
!=================== and check the constraints ============================================
!==========================================================================================
 if (InVar%ReadIFC.ne.1) then
   call tdep_build_pijN(InVar,ncoeff1st,proj1st,Pij_coeff,Pij_N,Shell1at,Sym)
   call tdep_build_phijNN(distance,InVar,ncoeff2nd,proj2nd,Phij_coeff,Phij_NN,Shell2at,Sym)
   call tdep_check_constraints(distance,InVar,Phij_NN,Pij_N,1)
 end if
 ABI_FREE(proj2nd)
 ABI_FREE(Pij_coeff)
 ABI_FREE(Phij_coeff)

!==========================================================================================
!===================== Compute the phonons density of states ==============================
!==========================================================================================
 call tdep_calc_phdos(Crystal,Ifc,InVar,Lattice,natom,natom_unitcell,Phij_NN,PHdos,Qpt,Rlatt4Abi,Shell2at,Sym)
 call tdep_destroy_shell(natom,2,Shell2at)
 ABI_FREE(Rlatt4Abi)

!==========================================================================================
!============= Compute the dynamical matrix and phonon spectrum ===========================
!================ then print Dij, omega, eigenvectors =====================================
!==========================================================================================
 write(stdout,*)' '
 write(stdout,*) '#############################################################################'
 write(stdout,*) '######################## Dynamical matrix ###################################'
 write(stdout,*) '#############################################################################'
 open(unit=53,file=trim(InVar%output_prefix)//'omega.dat')
 open(unit=52,file=trim(InVar%output_prefix)//'dij.dat')
 open(unit=51,file=trim(InVar%output_prefix)//'eigenvectors.dat')
 ABI_MALLOC(dij   ,(3*InVar%natom_unitcell,3*InVar%natom_unitcell))
 ABI_MALLOC(eigenV,(3*InVar%natom_unitcell,3*InVar%natom_unitcell))
 ABI_MALLOC(omega,(3*InVar%natom_unitcell))
 call tdep_init_eigen2nd(Eigen2nd,InVar%natom_unitcell,Qpt%nqpt)
 do iqpt=1,Qpt%nqpt
   dij(:,:)=zero ; eigenV(:,:)=zero ; omega(:)=zero
   qpt_cart(:)=Qpt%qpt_cart(:,iqpt)
   call tdep_calc_dij (dij,eigenV,iqpt,InVar,Lattice,omega,Phij_NN,qpt_cart,Rlatt_cart)
   call tdep_write_dij(dij,eigenV,iqpt,InVar,Lattice,omega,qpt_cart)
   Eigen2nd%eigenval(:,iqpt)=  omega(:)
   Eigen2nd%eigenvec(:,:,iqpt)=eigenV(:,:)
 end do
 ABI_FREE(dij)
 ABI_FREE(eigenV)
 ABI_FREE(omega)
 close(53)
 close(52)
 close(51)
 call tdep_write_yaml(Eigen2nd,Qpt,InVar%output_prefix)
 write(InVar%stdout,'(a)') ' See the dij.dat, omega.dat and eigenvectors files'
!==========================================================================================
!===================== Compute the elastic constants ======================================
!==========================================================================================
 call tdep_calc_elastic(Phij_NN,distance,InVar,Lattice)

!==========================================================================================
!=========== Compute U_0, the "free energy" and the forces (from the model) ===============
!==========================================================================================
 ABI_MALLOC(Forces_TDEP,(3*InVar%natom*InVar%nstep)); Forces_TDEP(:)=0.d0
 call tdep_calc_model(Free_Anh,Forces_MD,Forces_TDEP,InVar,Phij_NN,Pij_N,ucart,U0)

!==========================================================================================
!===================== Compute the thermodynamical quantities =============================
!==========================================================================================
 call tdep_calc_thermo(Free_Anh,InVar,Lattice,PHdos,U0)

 if (InVar%Order==2) then
   call tdep_print_Aknowledgments(InVar)
   call xmpi_end()
   stop
 end if

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=# CALCULATION OF THE 3rd ORDER =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

!==========================================================================================
!============= Initialize Psij_coeff from Moore-Penrose coefficients ======================
!==========================================================================================
 ABI_MALLOC(Psij_coeff,(ncoeff3rd,1)); Psij_coeff(:,:)=0.d0
 Psij_coeff(:,:)=MP_coeff(ncoeff2nd+1:ntotcoeff,:)
 ABI_FREE(MP_coeff)

!==========================================================================================
!============= Compute all the IFC coefficients of the Psi_ijk matrix =====================
!==========================================================================================
 ABI_MALLOC(Psij_ref,(3,3,3,Shell3at%nshell)) ; Psij_ref(:,:,:,:)=0.d0
 call tdep_calc_psijtot(distance,InVar,ncoeff3rd,proj3rd,Psij_coeff,Psij_ref,Shell3at,Sym)
 ABI_MALLOC(ftot3,(3*InVar%natom,InVar%nstep)); ftot3(:,:)=0.d0
 call tdep_check_constraints(distance,InVar,Phij_NN,Pij_N,Shell3at%nshell,ftot3,Psij_ref,Shell3at,Sym,ucart)
 ABI_FREE(Psij_coeff)
 ABI_FREE(proj3rd)
 call tdep_write_gruneisen(distance,Eigen2nd,InVar,Lattice,Psij_ref,Qpt,Rlatt_cart,Shell3at,Sym)
 call tdep_calc_alpha_gamma(Crystal,distance,DDB,Ifc,InVar,Lattice,Psij_ref,Rlatt_cart,Shell3at,Sym)
!FB stop
 ABI_FREE(Rlatt_cart)
 call tdep_destroy_eigen2nd(Eigen2nd)
 call tdep_destroy_shell(natom,3,Shell3at)
 call tdep_calc_model(Free_Anh,Forces_MD,Forces_TDEP,InVar,Phij_NN,Pij_N,ucart,U0,ftot3)
 ABI_FREE(ftot3)
 ABI_FREE(distance)
 ABI_FREE(Forces_MD)
 ABI_FREE(Forces_TDEP)
 ABI_FREE(Pij_N)
 ABI_FREE(Phij_NN)
 ABI_FREE(Psij_ref)
 ABI_FREE(ucart)
!==========================================================================================
!================= Write the last informations (aknowledgments...)  =======================
!==========================================================================================
 call tdep_print_Aknowledgments(InVar)

 call abinit_doctor("__fftprof")

100 call xmpi_end()

 end program atdep
!!***
