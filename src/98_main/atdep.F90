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
!!      tdep_calc_phi1,tdep_write_phi1
!!      tdep_calc_phi2,tdep_write_phi2,tdep_calc_phi3ref,tdep_write_phi3
!!      tdep_calc_phi4ref,tdep_write_phi4,tdep_calc_alpha_gamma
!!      tdep_calc_constraints,tdep_calc_elastic,tdep_calc_model
!!      tdep_calc_moorepenrose,tdep_calc_phdos,tdep_calc_phi2fcoeff
!!      tdep_calc_phi3fcoeff,tdep_calc_phi4fcoeff,tdep_calc_thermo,tdep_destroy_shell
!!      tdep_destroy_eigen2nd,tdep_init_crystal,tdep_init_ddb
!!      tdep_init_eigen2nd,tdep_init_ifc,tdep_init_shell2at
!!      tdep_init_shell3at,tdep_init_shell4at,tdep_make_latt,tdep_make_qptpath
!!      tdep_make_sym,tdep_matchideal2average,tdep_print_aknowledgments
!!      tdep_readecho,tdep_write_gruneisen
!!      xmpi_end,xmpi_init
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

  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_tdep_abitypes,    only : Qbz_type, tdep_init_crystal, tdep_init_ifc, tdep_init_ddb, tdep_write_ifc, &
&                                tdep_write_ddb  
  use m_tdep_phi4,        only : tdep_calc_phi4fcoeff, tdep_calc_phi4ref, tdep_write_phi4, tdep_calc_ftot4
  use m_tdep_phi3,        only : tdep_calc_phi3fcoeff, tdep_calc_phi3ref, tdep_write_phi3, tdep_calc_ftot3, &
&                                tdep_calc_alpha_gamma, tdep_write_gruneisen, tdep_calc_lifetime1
  use m_tdep_phi2,        only : tdep_calc_phi2fcoeff, tdep_calc_phi1fcoeff, tdep_calc_phi2, tdep_write_phi2, tdep_calc_ftot2, &
&                                Eigen_Variables_type, tdep_init_eigen2nd, tdep_destroy_eigen2nd, tdep_calc_phi1, tdep_write_phi1
  use m_tdep_latt,        only : tdep_make_latt, Lattice_Variables_type
  use m_tdep_sym,         only : tdep_make_sym, Symetries_Variables_type
  use m_tdep_readwrite,   only : tdep_print_Aknowledgments, tdep_read_input, tdep_distrib_data, tdep_init_MPIdata, tdep_clean_MPI, &
&                                Input_Variables_type, MPI_enreg_type
  use m_tdep_utils,       only : Coeff_Moore_type, tdep_calc_MoorePenrose, tdep_MatchIdeal2Average, tdep_calc_model
  use m_tdep_qpt,         only : tdep_make_qptpath, Qpoints_type
  use m_tdep_phdos,       only : tdep_calc_phdos,tdep_calc_elastic,tdep_calc_thermo
  use m_tdep_shell,       only : Shell_Variables_type, tdep_init_shell2at, tdep_init_shell3at, tdep_init_shell4at, tdep_init_shell1at, &
&                                tdep_destroy_shell
  use m_tdep_constraints, only : tdep_calc_constraints, tdep_check_constraints

  implicit none

  integer :: natom,jatom,natom_unitcell,ncoeff1st,ncoeff2nd,ncoeff3rd,ncoeff4th,ntotcoeff,ntotconst
  integer :: stdout,nshell_max
  double precision :: U0
  double precision, allocatable :: ucart(:,:,:),proj1st(:,:,:),proj2nd(:,:,:),proj3rd(:,:,:),proj4th(:,:,:)
  double precision, allocatable :: proj_tmp(:,:,:),Forces_TDEP(:),Fresid(:)
  double precision, allocatable :: Phi1(:)  ,Phi1_coeff(:,:)
  double precision, allocatable :: Phi2(:,:),Phi2_coeff(:,:)
  double precision, allocatable :: Phi3_ref(:,:,:,:)  ,Phi3_coeff(:,:)
  double precision, allocatable :: Phi4_ref(:,:,:,:,:),Phi4_coeff(:,:)
  double precision, allocatable :: Forces_MD(:),MP_coeff(:,:)
  double precision, allocatable :: distance(:,:,:),Rlatt_cart(:,:,:),Rlatt4Abi(:,:,:)
  double precision, allocatable :: Phi1Ui(:),Phi2UiUj(:),Phi3UiUjUk(:),Phi4UiUjUkUl(:)
  double complex  , allocatable :: Gruneisen(:)
  type(args_t) :: args
  type(phonon_dos_type) :: PHdos
  type(Input_Variables_type) :: InVar
  type(Lattice_Variables_type) :: Lattice
  type(Symetries_Variables_type) :: Sym
  type(Qpoints_type) :: Qpt
  type(Qbz_type) :: Qbz
  type(ifc_type) :: Ifc
  type(ddb_type) :: DDB
  type(crystal_t) :: Crystal
  type(Shell_Variables_type) :: Shell1at,Shell2at,Shell3at,Shell4at
  type(Coeff_Moore_type) :: CoeffMoore
  type(Eigen_Variables_type) :: Eigen2nd_MP
  type(Eigen_Variables_type) :: Eigen2nd_path
  type(MPI_enreg_type) :: MPIdata
  type(abihist) :: Hist

!TEST
  integer, allocatable :: data_tmp(:,:),data_loc(:,:),data_gather(:,:),shft_step(:)
  integer, allocatable :: nstep_acc(:),tab_step(:)
  integer :: ii,istep,ierr,iproc,dim1,dim2,remain,idim1
!TEST

!******************************************************************

!==========================================================================================
!===================== Initialization & Reading  ==========================================
!==========================================================================================
! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)
! Initialize MPI
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
 call tdep_read_input(Hist,InVar)
 call tdep_init_MPIdata(InVar,MPIdata)
 call tdep_distrib_data(Hist,InVar,MPIdata)

!FB TEST
!FB InVar%nstep_tot=10
!FB dim1=1
!FB dim2=1
!FB ABI_MALLOC(data_tmp,(dim1,InVar%nstep_tot*dim2)); data_tmp(:,:)=zero
!FB do istep=1,InVar%nstep_tot
!FB   data_tmp(:,dim2*(istep-1)+1:dim1*istep)=istep
!FB end do
!FB if (MPIdata%iam_master) write(6,*) MPIdata%me_step,data_tmp(:,:)
!FB
!FB InVar%my_nstep=int(InVar%nstep_tot/MPIdata%nproc)
!FB remain=InVar%nstep_tot-MPIdata%nproc*InVar%my_nstep
!FB do ii=1,remain
!FB   if ((ii-1).eq.MPIdata%me_step) InVar%my_nstep=InVar%my_nstep+1
!FB end do
!FB!FB ABI_MALLOC(MPIdata%nstep_all,(MPIdata%nproc)); MPIdata%nstep_all(:)=zero
!FB call xmpi_allgather(InVar%my_nstep,MPIdata%nstep_all,MPIdata%comm_step,ierr)
!FB if (MPIdata%me_step.eq.MPIdata%master) write(6,*) 'Nstep(nproc)=',MPIdata%nstep_all(:)
!FB call flush_unit(6)
!FB call xmpi_barrier(MPIdata%comm_step)
!FB
!FB ABI_MALLOC(nstep_acc,(MPIdata%nproc+1)); nstep_acc(:)=zero
!FB nstep_acc(1)=0
!FB do ii=2,MPIdata%nproc+1
!FB   nstep_acc(ii)=nstep_acc(ii-1)+MPIdata%nstep_all(ii-1)
!FB end do
!FB if (MPIdata%me_step.eq.MPIdata%master) write(6,*) 'NSTEP_ACC=',nstep_acc(:)
!FB call flush_unit(6)
!FB call xmpi_barrier(MPIdata%comm_step)
!FB
!FB ABI_MALLOC(tab_step,(InVar%nstep_tot)); tab_step(:)=zero
!FB! First distrib
!FB do iproc=1,MPIdata%nproc
!FB   do istep=1,InVar%nstep_tot
!FB     if ((istep.gt.nstep_acc(iproc)).and.(istep.le.nstep_acc(iproc+1))) then
!FB       tab_step(istep)=iproc-1
!FB     end if
!FB   end do
!FB end do
!FB
!FB! Second distrib
!FB tab_step(:)=zero
!FB do istep=1,InVar%nstep_tot
!FB   do iproc=1,MPIdata%nproc
!FB     if (mod((istep-1),MPIdata%nproc).eq.(iproc-1)) then
!FB       tab_step(istep)=iproc-1
!FB     end if
!FB   end do
!FB end do
!FB
!FB if (MPIdata%me_step.eq.MPIdata%master) write(6,*) 'TAB_STEP=',tab_step(:)
!FB call flush_unit(6)
!FB call xmpi_barrier(MPIdata%comm_step)
!FB
!FB if (nstep_acc(MPIdata%nproc+1).ne.InVar%nstep_tot) then
!FB   write(6,*) 'STOP : pb in nstep_acc'
!FB   stop
!FB end if
!FB
!FB ii=1
!FB ABI_MALLOC(data_loc,(dim1,InVar%my_nstep*dim2)); data_loc(:,:)=zero
!FB do istep=1,InVar%nstep_tot
!FB   if (tab_step(istep).eq.MPIdata%me_step) then
!FB     data_loc(:,dim2*(ii-1)+1:dim2*ii)=data_tmp(:,dim2*(istep-1)+1:dim2*istep)
!FB     ii=ii+1
!FB   end if  
!FB end do
!FB
!FB do iproc=1,MPIdata%nproc
!FB   if (iproc-1.eq.MPIdata%me_step) write(6,*) 'DATA_LOC =',MPIdata%me_step,data_loc(:,:)
!FB end do
!FB call flush_unit(6)
!FB call xmpi_barrier(MPIdata%comm_step)
!FB
!FB!FB ABI_MALLOC(shft_step,(MPIdata%nproc)); shft_step(:)=zero
!FB!FB shft_step(1)=0
!FB!FB do ii=2,MPIdata%nproc
!FB!FB   shft_step(ii)=shft_step(ii-1)+MPIdata%nstep_all(ii-1)
!FB!FB end do
!FB!FB
!FB!FB ABI_MALLOC(data_gather,(dim1,InVar%nstep_tot*dim2)); data_gather(:,:)=zero
!FB!FB do idim1=1,dim1
!FB!FB   call xmpi_gatherv(data_loc(idim1,:),InVar%my_nstep*dim2,data_gather(idim1,:),MPIdata%nstep_all*dim2,dim2*shft_step,&
!FB!FB&                    MPIdata%master,MPIdata%comm_step,ierr)
!FB!FB end do
!FB!FB
!FB!FB if (MPIdata%me_step.eq.MPIdata%master) write(6,*) "============================================"
!FB!FB if (MPIdata%me_step.eq.MPIdata%master) write(6,*) MPIdata%me_step,data_gather(:,:)
!FB!FB call flush_unit(6)
!FB!!FBFB
!FB 
!FB ABI_MALLOC(data_gather,(dim1,InVar%nstep_tot*dim2)); data_gather(:,:)=zero
!FB call xmpi_allgatherv(data_loc,dim1*InVar%my_nstep,data_gather,dim1*MPIdata%nstep_all,dim1*nstep_acc(1:MPIdata%nproc),MPIdata%comm_step,ierr)
!FB!FB call xmpi_scatterv(data_gather,MPIdata%nstep_all,nstep_acc,data_loc,InVar%my_nstep,&
!FB!FB&                    MPIdata%master,MPIdata%comm_step,ierr)
!FB do iproc=1,MPIdata%nproc
!FB   if (MPIdata%me_step.eq.(iproc-1)) write(6,*) "============================================"
!FB   if (MPIdata%me_step.eq.(iproc-1)) write(6,*) MPIdata%me_step,data_gather(:,:)
!FB   call flush_unit(6)
!FB   call xmpi_barrier(MPIdata%comm_step)
!FB end do
!FB write(6,*) '====== END ======='
!FB call flush_unit(6)
!FB call abi_abort("COLL")
!TEST

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
 call tdep_make_sym(Invar,Lattice,MPIdata,Sym)

! Initialize the Brillouin zone and compute the q-points path 
 call tdep_make_qptpath(InVar,Lattice,MPIdata,Qpt)

!==========================================================================================
!======== 1/ Determine ideal positions and distances ======================================
!======== 2/ Find the matching between the ideal and average ==============================
!========   (from the MD simulations) positions. ==========================================
!======== 3/ Find the symmetry operation between the reference and image bonds ============
!======== 4/ Write output quantities needed to visualize the neighbouring distances =======
!==========================================================================================
 ABI_MALLOC(Rlatt4Abi ,(3,natom_unitcell,natom))   ; Rlatt4Abi (:,:,:)=0.d0
 ABI_MALLOC(distance,(natom,natom,4))              ; distance(:,:,:)=0.d0
 ABI_MALLOC(Rlatt_cart,(3,natom_unitcell,natom))   ; Rlatt_cart(:,:,:)=0.d0
 ABI_MALLOC(ucart,(3,natom,InVar%my_nstep))           ; ucart(:,:,:)=0.d0
 ABI_MALLOC(Forces_MD,(3*natom*InVar%my_nstep))       ; Forces_MD(:)=0.d0

 call tdep_MatchIdeal2Average(distance,Forces_MD,InVar,Lattice,MPIdata,Rlatt_cart,Rlatt4Abi,Sym,ucart)

!==========================================================================================
!============== Initialize Crystal and DDB ABINIT Datatypes ===============================
!==========================================================================================
 call tdep_init_crystal(Crystal,InVar,Lattice,Sym)
 call tdep_init_ddb(Crystal,DDB,InVar,Lattice,MPIdata,Qbz)

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=# CALCULATION OF THE 2nd ORDER =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

!==========================================================================================
!============== Initialize the Shell1at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(3,3,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell1at(distance,InVar,MPIdata,3,nshell_max,ncoeff1st,1,proj_tmp,Shell1at,Sym)
 ABI_MALLOC(proj1st  ,(3,3,Shell1at%nshell)) ; proj1st(:,:,:)=0.d0
 proj1st = reshape (proj_tmp, (/ 3,3,Shell1at%nshell /))
 ABI_FREE(proj_tmp)
!Rotational invariances (1st order)
!    constraints = 3
 CoeffMoore%nconst_1st=3**2

!==========================================================================================
!============== Initialize the Shell2at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(9,9,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell2at(distance,InVar,MPIdata,9,nshell_max,ncoeff2nd,2,proj_tmp,Shell2at,Sym)
 ABI_MALLOC(proj2nd  ,(9,9,Shell2at%nshell)) ; proj2nd(:,:,:)=0.d0
 proj2nd = reshape (proj_tmp, (/ 9,9,Shell2at%nshell /))
 ABI_FREE(proj_tmp)
!Rotational invariances (2nd order) + Symetry of the Dynamical Matrix + Huang invariances
!    constraints = natom*3**2 + (3*natom_unitcell)**2 + 3**4
 CoeffMoore%nconst_rot2nd = 3**3*natom_unitcell
 CoeffMoore%nconst_dynmat = (3*natom_unitcell)**2
 CoeffMoore%nconst_huang  = 3**4
 CoeffMoore%nconst_2nd = CoeffMoore%nconst_rot2nd + CoeffMoore%nconst_dynmat + CoeffMoore%nconst_huang

!==========================================================================================
!============== Initialize the IFC Abinit datatype ========================================
!==========================================================================================
 ABI_MALLOC(Phi1,(3*natom))         ; Phi1(:)  =0.d0
 ABI_MALLOC(Phi2,(3*natom,3*natom)) ; Phi2(:,:)=0.d0
 call tdep_init_ifc(Crystal,DDB,Ifc,InVar,Lattice,MPIdata,Phi2,Rlatt4Abi,Shell2at,Sym)

!==========================================================================================
!============== Initialize the Shell3at datatype ==========================================
!==========================================================================================
 ncoeff3rd=0
 CoeffMoore%nconst_3rd=0
 Shell3at%nshell=1
 if (InVar%Order.ge.3) then
   ABI_MALLOC(proj_tmp,(27,27,nshell_max)) ; proj_tmp(:,:,:)=0.d0
   call tdep_init_shell3at(distance,InVar,MPIdata,27,nshell_max,ncoeff3rd,3,proj_tmp,Shell3at,Sym)
 end if
 ABI_MALLOC(proj3rd  ,(27,27,Shell3at%nshell)) ; proj3rd(:,:,:)=0.d0
 if (InVar%Order.ge.3) then
   proj3rd = reshape (proj_tmp, (/ 27,27,Shell3at%nshell /))
   ABI_FREE(proj_tmp)
!  Rotational invariances (3rd order) + acoustic sum rules (3rd order)
!    constraints = natom_unitcell*natom*3**4 + natom_unitcell*natom*3**3
   CoeffMoore%nconst_rot3rd = 3**4*natom_unitcell*natom
   CoeffMoore%nconst_asr3rd = 3**3*natom_unitcell*natom
   CoeffMoore%nconst_3rd = CoeffMoore%nconst_rot3rd + CoeffMoore%nconst_asr3rd
 end if

!==========================================================================================
!============== Initialize the Shell4at datatype ==========================================
!==========================================================================================
 ncoeff4th=0
 CoeffMoore%nconst_4th=0
 Shell4at%nshell=1
 if (InVar%Order==4) then
   ABI_MALLOC(proj_tmp,(81,81,nshell_max)) ; proj_tmp(:,:,:)=0.d0
   call tdep_init_shell4at(distance,InVar,MPIdata,81,nshell_max,ncoeff4th,4,proj_tmp,Shell4at,Sym)
 end if
 ABI_MALLOC(proj4th  ,(81,81,Shell4at%nshell)) ; proj4th(:,:,:)=0.d0
 if (InVar%Order==4) then
   proj4th = reshape (proj_tmp, (/ 81,81,Shell4at%nshell /)) 
   ABI_FREE(proj_tmp)
!  Rotational invariances (4th order) + acoustic sum rules (4th order)
!    constraints = natom_unitcell*natom**2*3**5 + natom_unitcell*natom**2*3**4
!FB   CoeffMoore%nconst_rot4th = 3**5*natom_unitcell*natom**2
!FB4TH   CoeffMoore%nconst_asr4th = 3**4*natom_unitcell*natom**2
   CoeffMoore%nconst_rot4th = 0
   CoeffMoore%nconst_asr4th = 0
   CoeffMoore%nconst_4th = CoeffMoore%nconst_rot4th + CoeffMoore%nconst_asr4th
 end if
 CoeffMoore%ncoeff1st=ncoeff1st
 CoeffMoore%ncoeff2nd=ncoeff2nd
 CoeffMoore%ncoeff3rd=ncoeff3rd
 CoeffMoore%ncoeff4th=ncoeff4th

!==========================================================================================
!================= Build fcoeff and compute constraints ===================================
!==========================================================================================
 ABI_MALLOC(Forces_TDEP,(3*InVar%natom*InVar%my_nstep)); Forces_TDEP(:)=0.d0 
 ABI_MALLOC(Fresid     ,(3*InVar%natom*InVar%my_nstep)); Fresid(:)=Forces_MD(:)
 ABI_MALLOC(Phi1Ui      ,(InVar%my_nstep)); Phi1Ui      (:)=0.d0 
 ABI_MALLOC(Phi2UiUj    ,(InVar%my_nstep)); Phi2UiUj    (:)=0.d0 
 ABI_MALLOC(Phi3UiUjUk  ,(InVar%my_nstep)); Phi3UiUjUk  (:)=0.d0 
 ABI_MALLOC(Phi4UiUjUkUl,(InVar%my_nstep)); Phi4UiUjUkUl(:)=0.d0 
 ABI_MALLOC(Phi3_ref,(3,3,3,Shell3at%nshell))  ; Phi3_ref(:,:,:,:)=0.d0
 ABI_MALLOC(Phi4_ref,(3,3,3,3,Shell4at%nshell)); Phi4_ref(:,:,:,:,:)=0.d0
 write(6,*) 'Coeff=', CoeffMoore%ncoeff1st,CoeffMoore%ncoeff2nd,CoeffMoore%ncoeff3rd,CoeffMoore%ncoeff4th
 ntotcoeff=CoeffMoore%ncoeff1st  + CoeffMoore%ncoeff2nd  + CoeffMoore%ncoeff3rd  + CoeffMoore%ncoeff4th
 ntotconst=CoeffMoore%nconst_1st + CoeffMoore%nconst_2nd + CoeffMoore%nconst_3rd + CoeffMoore%nconst_4th
 CoeffMoore%ntotcoeff=ntotcoeff
 CoeffMoore%ntotconst=ntotconst
 write(6,*) 'bornes=',3*natom*InVar%my_nstep,ntotcoeff
 ABI_MALLOC(CoeffMoore%fcoeff,(3*natom*InVar%my_nstep,ntotcoeff)); CoeffMoore%fcoeff(:,:)=0.d0 
 if (InVar%ReadIFC.ne.1) then
   call tdep_calc_phi1fcoeff(CoeffMoore,InVar,proj1st,Shell1at,Sym)
   call tdep_calc_phi2fcoeff(CoeffMoore,InVar,proj2nd,Shell2at,Sym,ucart)
 end if 
 if (InVar%Order.ge.3) then
   call tdep_calc_phi3fcoeff(CoeffMoore,InVar,proj3rd,Shell3at,Sym,ucart)
 end if  
 if (InVar%Order.eq.4) then
   call tdep_calc_phi4fcoeff(CoeffMoore,InVar,proj4th,Shell4at,Sym,ucart)
 end if  
 call tdep_calc_constraints(CoeffMoore,distance,InVar,MPIdata,Shell1at%nshell,Shell2at%nshell,&
&                           Shell3at%nshell,Shell4at%nshell,proj1st,proj2nd,proj3rd,proj4th,&
&                           Shell1at,Shell2at,Shell3at,Shell4at,Sym) 

!==========================================================================================
!============= Compute the pseudo inverse using the Moore-Penrose method ==================
!==========================================================================================
 ABI_MALLOC(MP_coeff,(ntotcoeff,1)); MP_coeff(:,:)=0.d0
!==========================================================================================
!=================== If all the Orders are solved simultaneously ==========================
!==========================================================================================
 if (InVar%together.eq.1) then
   write(stdout,*) '############### (Solve simultaneously all the orders) #######################'  
   if (InVar%ReadIFC.ne.1) then
     call tdep_calc_MoorePenrose(CoeffMoore,Forces_MD,0,InVar,MP_coeff,MPIdata)
   end if  
   ABI_MALLOC(Phi1_coeff,(ncoeff1st,1)); Phi1_coeff(:,:)=MP_coeff(1:ncoeff1st,:)
   ABI_MALLOC(Phi2_coeff,(ncoeff2nd,1)); Phi2_coeff(:,:)=MP_coeff(ncoeff1st+1:ncoeff1st+ncoeff2nd,:)
   if (InVar%ReadIFC.ne.1) then
     call tdep_calc_phi1(InVar,ncoeff1st,proj1st,Phi1_coeff,Phi1,Shell1at,Sym)
     call tdep_calc_phi2(InVar,ncoeff2nd,proj2nd,Phi2_coeff,Phi2,Shell2at,Sym)
   end if
   ABI_FREE(proj2nd)
   ABI_FREE(Phi1_coeff)
   ABI_FREE(Phi2_coeff)
   call tdep_calc_ftot2(Forces_TDEP,InVar,Phi1,Phi1Ui,Phi2,Phi2UiUj,ucart) 
   if (InVar%Order.ge.3) then 
     ABI_MALLOC(Phi3_coeff,(ncoeff3rd,1)); Phi3_coeff(:,:)=0.d0
     Phi3_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+1:ncoeff1st+ncoeff2nd+ncoeff3rd,:)
     call tdep_calc_phi3ref(InVar,ncoeff3rd,proj3rd,Phi3_coeff,Phi3_ref,Shell3at)
     ABI_FREE(proj3rd)
     ABI_FREE(Phi3_coeff)
     call tdep_calc_ftot3(Forces_TDEP,InVar,MPIdata,Phi3_ref,Phi3UiUjUk,Shell3at,ucart,Sym) 
   end if
   if (InVar%Order.ge.4) then 
     ABI_MALLOC(Phi4_coeff,(ncoeff4th,1)); Phi4_coeff(:,:)=0.d0
     Phi4_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+ncoeff3rd+1:ntotcoeff,:)
     call tdep_calc_phi4ref(InVar,ncoeff4th,proj4th,Phi4_coeff,Phi4_ref,Shell4at)
     ABI_FREE(proj4th)
     ABI_FREE(Phi4_coeff)
     call tdep_calc_ftot4(Forces_TDEP,InVar,MPIdata,Phi4_ref,Phi4UiUjUkUl,Shell4at,ucart,Sym) 
   end if  

!==========================================================================================
!=================== If all the Orders are solved successively ============================
!==========================================================================================
 else if (InVar%together.eq.0) then
   write(stdout,*) '################## (Solve successively each order) ##########################'  
   do ii=1,InVar%Order-1
     write(stdout,*) ' For Order=',ii+1
     if (InVar%ReadIFC.eq.1) cycle
     call tdep_calc_MoorePenrose(CoeffMoore,Fresid,ii,InVar,MP_coeff,MPIdata)
     if (ii.eq.1) then
       ABI_MALLOC(Phi1_coeff,(ncoeff1st,1)); Phi1_coeff(:,:)=MP_coeff(1:ncoeff1st,:)
       ABI_MALLOC(Phi2_coeff,(ncoeff2nd,1)); Phi2_coeff(:,:)=MP_coeff(ncoeff1st+1:ncoeff1st+ncoeff2nd,:)
       if (InVar%ReadIFC.ne.1) then
         call tdep_calc_phi1(InVar,ncoeff1st,proj1st,Phi1_coeff,Phi1,Shell1at,Sym)
         call tdep_calc_phi2(InVar,ncoeff2nd,proj2nd,Phi2_coeff,Phi2,Shell2at,Sym)
       end if
       ABI_FREE(proj2nd)
       ABI_FREE(Phi1_coeff)
       ABI_FREE(Phi2_coeff)
       call tdep_calc_ftot2(Forces_TDEP,InVar,Phi1,Phi1Ui,Phi2,Phi2UiUj,ucart) 
     else if (ii.eq.2) then 
       ABI_MALLOC(Phi3_coeff,(ncoeff3rd,1)); Phi3_coeff(:,:)=0.d0
       Phi3_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+1:ncoeff1st+ncoeff2nd+ncoeff3rd,:)
       call tdep_calc_phi3ref(InVar,ncoeff3rd,proj3rd,Phi3_coeff,Phi3_ref,Shell3at)
       ABI_FREE(proj3rd)
       ABI_FREE(Phi3_coeff)
       call tdep_calc_ftot3(Forces_TDEP,InVar,MPIdata,Phi3_ref,Phi3UiUjUk,Shell3at,ucart,Sym) 
     else if (ii.eq.3) then 
       ABI_MALLOC(Phi4_coeff,(ncoeff4th,1)); Phi4_coeff(:,:)=0.d0
       Phi4_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+ncoeff3rd+1:ntotcoeff,:)
       call tdep_calc_phi4ref(InVar,ncoeff4th,proj4th,Phi4_coeff,Phi4_ref,Shell4at)
       ABI_FREE(proj4th)
       ABI_FREE(Phi4_coeff)
       call tdep_calc_ftot4(Forces_TDEP,InVar,MPIdata,Phi4_ref,Phi4UiUjUkUl,Shell4at,ucart,Sym) 
     end if  
     Fresid(:)=Forces_MD(:)-Forces_TDEP(:)
   end do  
 end if  
 ABI_FREE(CoeffMoore%const)
 ABI_FREE(CoeffMoore%fcoeff)
 ABI_FREE(Fresid)
 ABI_FREE(MP_coeff)

!==========================================================================================
!=================== Write the IFC and check the constraints ==============================
!==========================================================================================
 call tdep_write_phi1(InVar,Phi1)
 call tdep_write_phi2(distance,InVar,MPIdata,Phi2,Shell2at)
 if (InVar%Order.ge.3) then 
   call tdep_write_phi3(distance,InVar,Phi3_ref,Shell3at,Sym)
 end if  
 if (InVar%Order.ge.4) then 
   call tdep_write_phi4(distance,InVar,Phi4_ref,Shell4at,Sym)
 end if  

 call tdep_check_constraints(distance,InVar,Phi2,Phi1,Shell3at%nshell,Shell4at%nshell,&
&                            Phi3_ref,Phi4_ref,Shell3at,Shell4at,Sym,ucart)
 ABI_FREE(Phi1)
 ABI_FREE(ucart)
   
!==========================================================================================
!===================== Compute the phonon spectrum, the DOS, ==============================
!=====================  the dynamical matrix and write them ===============================
!==========================================================================================
 write(stdout,*)' '
 write(stdout,*) '#############################################################################'
 write(stdout,*) '############## Compute the phonon spectrum, the DOS, ########################'
 write(stdout,*) '##############  the dynamical matrix and write them  ########################'
 write(stdout,*) '#############################################################################'
 call tdep_init_eigen2nd(Eigen2nd_MP,InVar%natom_unitcell,Qbz%nqbz)
!FB call tdep_init_eigen2nd(Eigen2nd_MP,InVar%natom_unitcell,Qbz%nqibz)
 call tdep_init_eigen2nd(Eigen2nd_path,InVar%natom_unitcell,Qpt%nqpt)
 call tdep_calc_phdos(Crystal,DDB,Eigen2nd_MP,Eigen2nd_path,Ifc,InVar,Lattice,MPIdata,natom,&
&                          natom_unitcell,Phi2,PHdos,Qbz,Qpt,Rlatt4abi,Rlatt_cart,Shell2at,Sym)
 call tdep_destroy_shell(natom,2,Shell2at)
 ABI_FREE(Rlatt4Abi)
 write(InVar%stdout,'(a)') ' See the dij.dat, omega.dat and eigenvectors files'
 write(InVar%stdout,'(a)') ' See also the DDB file'

!==========================================================================================
!===================== Compute the elastic constants ======================================
!==========================================================================================
 call tdep_calc_elastic(Phi2,distance,InVar,Lattice)
 ABI_FREE(Phi2)

!==========================================================================================
!=========== Compute U_0, the "free energy" and the forces (from the model) ===============
!==========================================================================================
 call tdep_calc_model(Forces_MD,Forces_TDEP,InVar,MPIdata,Phi1Ui,Phi2UiUj,&
&                     Phi3UiUjUk,Phi4UiUjUkUl,U0) 
 ABI_FREE(Forces_MD)
 ABI_FREE(Forces_TDEP)
 ABI_FREE(Phi1Ui)
 ABI_FREE(Phi2UiUj)
 ABI_FREE(Phi3UiUjUk)
 ABI_FREE(Phi4UiUjUkUl)

!==========================================================================================
!===================== Compute the thermodynamical quantities =============================
!==========================================================================================
 call tdep_calc_thermo(InVar,Lattice,MPIdata,PHdos,U0)

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

 if (MPIdata%iam_master) then
   call tdep_write_gruneisen(Crystal,distance,Eigen2nd_path,Ifc,InVar,Lattice,Phi3_ref,Qpt,Rlatt_cart,Shell3at,Sym)
 end if  
 call tdep_calc_alpha_gamma(distance,DDB,Eigen2nd_MP,InVar,Lattice,MPIdata,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)

!FB Begin Lifetime
!FB call tdep_calc_lifetime1(Crystal,distance,Eigen2nd_MP,Ifc,InVar,Lattice,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)
!FB End Lifetime

 ABI_FREE(Rlatt_cart)
 call tdep_destroy_eigen2nd(Eigen2nd_path)
 call tdep_destroy_eigen2nd(Eigen2nd_MP)
 call tdep_destroy_shell(natom,3,Shell3at)
 ABI_FREE(distance)
 ABI_FREE(Phi3_ref)
 if (InVar%Order.eq.4) then
   ABI_FREE(Phi4_ref)
   call tdep_destroy_shell(natom,4,Shell4at)
 end if  
!==========================================================================================
!================= Write the last informations (aknowledgments...)  =======================
!==========================================================================================
 call tdep_print_Aknowledgments(InVar)

 call tdep_clean_mpi(MPIdata)

!FB if (iam_master) close(stdout)
 call abinit_doctor("__fftprof")

100 call xmpi_end()

 end program atdep
!!***
