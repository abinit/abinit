!{\src2tex{textfont=tt}}
!!****p* ABINIT/tdep
!! NAME
!! tdep
!!
!! FUNCTION
!! Calculations of phonons using molecular dynamic simulations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FB,JB)
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
!!      tdep_build_phijnn,tdep_calc_dij,tdep_calc_elastic,tdep_calc_model
!!      tdep_calc_moorepenrose,tdep_calc_phdos,tdep_calc_phijfcoeff
!!      tdep_calc_thermo,tdep_destroy_shell,tdep_init_crystal,tdep_init_ddb
!!      tdep_init_eigen2nd,tdep_init_ifc,tdep_init_shell2at,tdep_make_latt
!!      tdep_make_qptpath,tdep_make_sym,tdep_matchideal2average
!!      tdep_print_aknowledgments,tdep_readecho,tdep_write_dij,tdep_write_yaml
!!      xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program tdep

  use defs_basis
  use m_abicore
  use m_phonons
  use m_xmpi,             only : xmpi_init, xmpi_end
  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_dynmat,           only : ftifc_r2q, ftifc_q2r, asrif9
  use m_tdep_abitypes,    only : tdep_init_crystal, tdep_init_ifc, tdep_init_ddb, tdep_write_ifc
  use m_tdep_phij,        only : tdep_calc_phijfcoeff, tdep_build_phijNN, tdep_calc_dij, tdep_write_dij, &
&                                Eigen_Variables_type, tdep_init_eigen2nd, tdep_destroy_eigen2nd, tdep_write_yaml
  use m_tdep_latt,        only : tdep_make_latt, Lattice_Variables_type
  use m_tdep_sym,         only : tdep_make_sym, Symetries_Variables_type
  use m_tdep_readwrite,   only : tdep_print_Aknowledgments, tdep_ReadEcho, Input_Variables_type
  use m_tdep_utils,       only : Coeff_Moore_type, tdep_calc_MoorePenrose, tdep_MatchIdeal2Average, tdep_calc_model
  use m_tdep_qpt,         only : tdep_make_qptpath, Qpoints_type
  use m_tdep_phdos,       only : tdep_calc_phdos,tdep_calc_elastic,tdep_calc_thermo
  use m_tdep_shell,       only : Shell_Variables_type, tdep_init_shell2at, tdep_init_shell3at, tdep_destroy_shell
#ifdef HAVE_NETCDF
  use netcdf
#endif
  use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep'
!End of the abilint section

  implicit none

  integer :: natom,natom_unitcell,ntotcoeff,nshell_max
  integer :: order,stdout,norder,iqpt
  double precision :: U0,DeltaFree_AH2
  double precision, allocatable :: ucart(:,:,:),proj(:,:,:),proj_tmp(:,:,:),Forces_TDEP(:)
!FB  double precision, allocatable :: fcoeff(:,:),Phij_coeff(:,:),Forces_MD(:),Phij_NN(:,:)
  double precision, allocatable :: Phij_coeff(:,:),Forces_MD(:),Phij_NN(:,:)
  double precision, allocatable :: distance(:,:,:),Rlatt_cart(:,:,:),Rlatt4Abi(:,:,:)
  double precision, allocatable :: omega (:)
  double precision :: qpt_cart(3)
  double complex  , allocatable :: dij(:,:),eigenV(:,:)
  type(phonon_dos_type) :: PHdos
  type(Input_Variables_type) :: InVar
  type(Lattice_Variables_type) :: Lattice
  type(Symetries_Variables_type) :: Sym
  type(Qpoints_type) :: Qpt
  type(ifc_type) :: Ifc
  type(ddb_type) :: DDB
  type(crystal_t) :: Crystal
  type(Shell_Variables_type) :: Shell2at
  type(Coeff_Moore_type) :: CoeffMoore
  type(Eigen_Variables_type) :: Eigen2nd

!******************************************************************

!==========================================================================================
!===================== Initialization & Reading  ==========================================
!==========================================================================================
 call xmpi_init()

! Read input values from the input.in input file
 call tdep_ReadEcho(InVar)

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
 ABI_MALLOC(ucart,(3,natom,InVar%nstep))    ; ucart(:,:,:)=0.d0
 ABI_MALLOC(Forces_MD,(3*natom*InVar%nstep)); Forces_MD(:)=0.d0

 write(InVar%stdout,*) "Matching structure"
 call flush_unit(InVar%stdout)
 call tdep_MatchIdeal2Average(distance,Forces_MD,InVar,Lattice,Rlatt_cart,Rlatt4Abi,Sym,ucart)
 call flush_unit(InVar%stdout)

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
 order=2
 norder=3**order
 write(InVar%stdout,*) ' '
 write(InVar%stdout,*) '#############################################################################'
 write(InVar%stdout,*) '############################## SECOND ORDER  ################################'
 write(InVar%stdout,*) '################ Now, find the number of coefficients for ###################'
 write(InVar%stdout,*) '########################## a reference interaction ##########################'
 write(InVar%stdout,*) '#############################################################################'
 call flush_unit(InVar%stdout)
 
!==========================================================================================
!============== Initialize the Shell2at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(norder,norder,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell2at(distance,InVar,norder,nshell_max,ntotcoeff,order,proj_tmp,Shell2at,Sym)
 ABI_MALLOC(proj    ,(norder,norder,Shell2at%nshell)) ; proj(:,:,:)=0.d0
 proj = reshape (proj_tmp, (/ norder,norder,Shell2at%nshell /)) 
 ABI_FREE(proj_tmp)

!==========================================================================================
!============== Initialize the IFC Abinit datatype ========================================
!==========================================================================================
 ABI_MALLOC(Phij_NN,(3*natom,3*natom)) ; Phij_NN(:,:)=0.d0
 call tdep_init_ifc(Crystal,DDB,Ifc,InVar,Lattice,Phij_NN,Rlatt4Abi,Shell2at,Sym)

!==========================================================================================
!============= Build fcoeff, needed for the Moore-Penrose method just below ===============
!==========================================================================================
 ABI_MALLOC(CoeffMoore%fcoeff,(3*natom*InVar%nstep,ntotcoeff)); CoeffMoore%fcoeff(:,:)=0.d0 
 if (InVar%ReadIFC.ne.1) then
   call tdep_calc_phijfcoeff(InVar,ntotcoeff,proj,Shell2at,Sym,ucart,CoeffMoore%fcoeff)
 end if

!==========================================================================================
!============= Compute the pseudo inverse using the Moore-Penrose method ==================
!==========================================================================================
 ABI_MALLOC(Phij_coeff,(ntotcoeff,1)); Phij_coeff(:,:)=0.d0
 if (InVar%ReadIFC.ne.1) then
   call tdep_calc_MoorePenrose(Forces_MD,CoeffMoore,InVar,ntotcoeff,Phij_coeff)
 end if

!==========================================================================================
!============= Reorganize the IFC coefficients into the whole Phij_NN matrix ==============
!==========================================================================================
 if (InVar%ReadIFC.ne.1) then
   call tdep_build_phijNN(distance,InVar,ntotcoeff,proj,Phij_coeff,Phij_NN,Shell2at,Sym)
 end if
 ABI_FREE(Phij_coeff)
 ABI_FREE(proj)

!==========================================================================================
!===================== Compute the phonons density of states ==============================
!==========================================================================================
 call tdep_calc_phdos(Crystal,Ifc,InVar,Lattice,natom,natom_unitcell,Phij_NN,PHdos,Qpt,Rlatt4Abi,Shell2at,Sym)
 ABI_FREE(Rlatt4Abi)
 call tdep_destroy_shell(natom,order,Shell2at)

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
 call tdep_calc_model(DeltaFree_AH2,distance,Forces_MD,Forces_TDEP,InVar,Phij_NN,ucart,U0) 

!==========================================================================================
!===================== Compute the thermodynamical quantities =============================
!==========================================================================================
 call tdep_calc_thermo(DeltaFree_AH2,InVar,PHdos,U0)

 !if (InVar%Order==2) then
 !  call tdep_print_Aknowledgments(InVar)
 !  call xmpi_end()
 !  stop
 !end if  

 ABI_FREE(Rlatt_cart)
 ABI_FREE(distance)
 ABI_FREE(Forces_MD)
 ABI_FREE(Forces_TDEP)
 ABI_FREE(Phij_NN)
 ABI_FREE(ucart)
!==========================================================================================
!================= Write the last informations (aknowledgments...)  =======================
!==========================================================================================
 call tdep_print_Aknowledgments(InVar)
 call xmpi_end()

 end program tdep
!!***
