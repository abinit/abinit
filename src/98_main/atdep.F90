!!****p* ABINIT/atdep
!! NAME
!! atdep
!!
!! FUNCTION
!! Calculations of phonons using molecular dynamic simulations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2025 ABINIT group (FB,JB)
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
&                                tdep_destroy_qbz, tdep_ifc2phi2
  use m_tdep_phi4,        only : tdep_calc_phi4fcoeff, tdep_calc_phi4ref, tdep_write_phi4, tdep_calc_ftot4
  use m_tdep_phi3,        only : tdep_calc_phi3fcoeff, tdep_calc_phi3ref, tdep_write_phi3, tdep_calc_ftot3, &
&                                tdep_calc_alpha_gamma, tdep_write_gruneisen
  use m_tdep_phi2,        only : tdep_calc_phi2fcoeff, tdep_calc_phi1fcoeff, tdep_calc_phi2, tdep_write_phi2, tdep_calc_ftot2, &
&                                Eigen_type, tdep_init_eigen2nd, tdep_destroy_eigen2nd, tdep_calc_phi1, &
&                                tdep_write_phi1, tdep_init_phi2, tdep_destroy_phi2, Phi2_type
  use m_tdep_latt,        only : Lattice_type, tdep_make_latt, tdep_shift_xred
  use m_tdep_sym,         only : tdep_make_sym, Symetries_type, tdep_destroy_sym
  use m_tdep_readwrite,   only : tdep_print_Aknowledgments, tdep_read_input, tdep_distrib_data, tdep_init_MPIdata, &
&                                tdep_destroy_mpidata, Input_type, MPI_enreg_type, tdep_destroy_invar, version_string
  use m_tdep_utils,       only : Coeff_Moore_type, tdep_calc_MoorePenrose, tdep_MatchIdeal2Average, tdep_calc_model
  use m_tdep_qpt,         only : tdep_make_qptpath, Qpoints_type, tdep_destroy_qpt
  use m_tdep_phdos,       only : tdep_calc_phdos,tdep_calc_elastic,tdep_calc_thermo
  use m_tdep_shell,       only : Shell_type, tdep_init_shell2at, tdep_init_shell3at, tdep_init_shell4at, &
&                                tdep_init_shell1at, tdep_destroy_shell
  use m_tdep_constraints, only : tdep_calc_constraints, tdep_check_constraints

  implicit none

  integer :: natom,natom_unitcell,ncoeff1st,ncoeff2nd,ncoeff3rd,ncoeff4th,ntotcoeff,ntotconst
  integer :: stdout,stdlog,nshell_max,ii,jj,ishell,istep,iatom
  integer :: print_mem_report
  double precision :: U0
  real(dp) :: tcpu, tcpui, twall, twalli
  real(dp) :: tsec(2)
  character(len = 24):: start_datetime
  double precision, allocatable :: ucart(:,:,:),proj1st(:,:,:),proj2nd(:,:,:),proj3rd(:,:,:),proj4th(:,:,:)
  double precision, allocatable :: proj_tmp(:,:,:),Forces_TDEP(:),Fresid(:)
  double precision, allocatable :: Phi1(:)  ,Phi1_coeff(:,:)
  double precision, allocatable :: Phi2_coeff(:,:)
  double precision, allocatable :: Phi3_ref(:,:,:,:)  ,Phi3_coeff(:,:)
  double precision, allocatable :: Phi4_ref(:,:,:,:,:),Phi4_coeff(:,:)
  double precision, allocatable :: Forces_MD(:),MP_coeff(:,:)
  double precision, allocatable :: distance(:,:,:),Rlatt_cart(:,:,:),Rlatt4Abi(:,:,:)
  double precision, allocatable :: Phi1Ui(:),Phi2UiUj(:),Phi3UiUjUk(:),Phi4UiUjUkUl(:)
  type(args_t) :: args
  type(phdos_t) :: PHdos
  type(Phi2_type) :: Phi2
  type(Input_type) :: Invar
  type(Lattice_type) :: Lattice
  type(Symetries_type) :: Sym
  type(Qpoints_type) :: Qpt
  type(Qbz_type) :: Qbz
  type(ifc_type) :: Ifc
  type(ddb_type) :: DDB
  type(crystal_t) :: Crystal
  type(Shell_type) :: Shell1at,Shell2at,Shell3at,Shell4at
  type(Coeff_Moore_type) :: CoeffMoore
  type(Eigen_type) :: Eigen2nd_MP
  type(Eigen_type) :: Eigen2nd_path
  type(MPI_enreg_type) :: MPIdata
  type(abihist) :: Hist

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
 call tdep_distrib_data(Hist,Invar,MPIdata)
 call abihist_free(Hist)

 if (args%dry_run /= 0) then
   call wrtout(std_out, "Dry run mode. Exiting after have read the input")
   goto 100
 end if

! Initialize basic quantities
 print_mem_report = 1
 natom            = Invar%natom
 natom_unitcell   = Invar%natom_unitcell
 stdout           = Invar%stdout
 stdlog           = Invar%stdlog
 nshell_max       = 500

!==========================================================================================
!============== Define the ideal lattice, symmetries and Brillouin zone ===================
!==========================================================================================

! Shift xred to keep atoms in the same unit cell at each step.
 call tdep_shift_xred(Invar,MPIdata)

! Define all the quantities needed to buid the lattice (rprim*, acell*, brav*...)
 call tdep_make_latt(Invar,Lattice)

! Compute all the symmetries coming from the bravais lattice
 call tdep_make_sym(Invar,Lattice,MPIdata,Sym)

! Initialize the Brillouin zone and compute the q-points path
 call tdep_make_qptpath(Invar,Lattice,MPIdata,Qpt)

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
 ABI_MALLOC(ucart,(3,natom,Invar%my_nstep))        ; ucart(:,:,:)=0.d0
 ABI_MALLOC(Forces_MD,(3*natom*Invar%my_nstep))    ; Forces_MD(:)=0.d0

 call tdep_MatchIdeal2Average(distance,Forces_MD,Invar,Lattice,MPIdata,Rlatt_cart,Rlatt4Abi,Sym,ucart)

!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=# CALCULATION OF THE 2nd ORDER =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
!#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

!==========================================================================================
!============== Initialize the Shell1at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(3,3,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell1at(distance,Invar,MPIdata,3,nshell_max,ncoeff1st,1,proj_tmp,Shell1at,Sym)
 ABI_MALLOC(proj1st  ,(3,3,Shell1at%nshell)) ; proj1st(:,:,:)=0.d0
 do ishell=1,Shell1at%nshell
   proj1st(:,:,ishell) = proj_tmp(:,:,ishell)
 end do
 ABI_FREE(proj_tmp)
!Rotational invariances (1st order)
!    constraints = 3
 CoeffMoore%nconst_1st=3**2

!==========================================================================================
!============== Initialize the Shell2at datatype ==========================================
!==========================================================================================
 ABI_MALLOC(proj_tmp,(9,9,nshell_max)) ; proj_tmp(:,:,:)=0.d0
 call tdep_init_shell2at(distance,Invar,MPIdata,9,nshell_max,ncoeff2nd,2,proj_tmp,Shell2at,Sym)
 ABI_MALLOC(proj2nd  ,(9,9,Shell2at%nshell)) ; proj2nd(:,:,:)=0.d0
 do ishell=1,Shell2at%nshell
   proj2nd(:,:,ishell) = proj_tmp(:,:,ishell)
 end do
 ABI_FREE(proj_tmp)
!Rotational invariances (2nd order) + Symetry of the Dynamical Matrix + Huang invariances
!    constraints = natom*3**2 + (3*natom_unitcell)**2 + 3**4
 CoeffMoore%nconst_rot2nd = 3**3*natom_unitcell
 CoeffMoore%nconst_dynmat = (3*natom_unitcell)**2
 CoeffMoore%nconst_huang  = 3**4
 CoeffMoore%nconst_2nd = CoeffMoore%nconst_rot2nd + CoeffMoore%nconst_dynmat + CoeffMoore%nconst_huang

!==========================================================================================
!============== Initialize Phi2 datatype ==================================================
!==========================================================================================
 ABI_MALLOC(Phi1,(3*natom)); Phi1(:)  =0.d0
 call tdep_init_phi2(Phi2,Invar%loto,natom)

!==========================================================================================
!============== Initialize Crystal, DDB, and IFC ABINIT Datatypes =========================
!==========================================================================================
 call tdep_init_crystal(Crystal,Invar,Lattice,Sym)
 call tdep_init_ddb(Crystal,DDB,Invar,Lattice,MPIdata,Qbz)
 call tdep_init_ifc(Crystal,DDB,Ifc,Invar,Lattice,MPIdata,Phi2,Rlatt4Abi,Shell2at,Sym)

!==========================================================================================
!============== Initialize the Shell3at datatype ==========================================
!==========================================================================================
 ncoeff3rd=0
 CoeffMoore%nconst_3rd=0
 Shell3at%nshell=1
 if (Invar%order.ge.3) then
   ABI_MALLOC(proj_tmp,(27,27,nshell_max)) ; proj_tmp(:,:,:)=0.d0
   call tdep_init_shell3at(distance,Invar,MPIdata,27,nshell_max,ncoeff3rd,3,proj_tmp,Shell3at,Sym)
   ABI_MALLOC(proj3rd  ,(27,27,Shell3at%nshell)) ; proj3rd(:,:,:)=0.d0
   do ishell=1,Shell3at%nshell
     proj3rd(:,:,ishell) = proj_tmp(:,:,ishell)
   end do
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
 if (Invar%order==4) then
   ABI_MALLOC(proj_tmp,(81,81,nshell_max)) ; proj_tmp(:,:,:)=0.d0
   call tdep_init_shell4at(distance,Invar,MPIdata,81,nshell_max,ncoeff4th,4,proj_tmp,Shell4at,Sym)
   ABI_MALLOC(proj4th  ,(81,81,Shell4at%nshell)) ; proj4th(:,:,:)=0.d0
   do ishell=1,Shell4at%nshell
     proj4th(:,:,ishell) = proj_tmp(:,:,ishell)
   end do
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
 ABI_MALLOC(Forces_TDEP,(3*Invar%natom*Invar%my_nstep)); Forces_TDEP(:)=0.d0
 ABI_MALLOC(Fresid     ,(3*Invar%natom*Invar%my_nstep)); Fresid(:)=Forces_MD(:)
 ABI_MALLOC(Phi1Ui      ,(Invar%my_nstep)); Phi1Ui      (:)=0.d0
 ABI_MALLOC(Phi2UiUj    ,(Invar%my_nstep)); Phi2UiUj    (:)=0.d0
 ABI_MALLOC(Phi3UiUjUk  ,(Invar%my_nstep)); Phi3UiUjUk  (:)=0.d0
 ABI_MALLOC(Phi4UiUjUkUl,(Invar%my_nstep)); Phi4UiUjUkUl(:)=0.d0
 ntotcoeff=CoeffMoore%ncoeff1st  + CoeffMoore%ncoeff2nd  + CoeffMoore%ncoeff3rd  + CoeffMoore%ncoeff4th
 ntotconst=CoeffMoore%nconst_1st + CoeffMoore%nconst_2nd + CoeffMoore%nconst_3rd + CoeffMoore%nconst_4th
 CoeffMoore%ntotcoeff=ntotcoeff
 CoeffMoore%ntotconst=ntotconst

!LOTO
!Remove the supercell contribution (the "LR part") included in total forces
!before computing the "SR part". The full "LR part" will be added later.
 if (Invar%loto) then
   call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,Invar%natom_unitcell,1,Phi2,Rlatt_cart,Shell2at,Sym)
   call tdep_calc_ftot2(Forces_TDEP,Invar,Phi1,Phi1Ui,Phi2%Tot,Phi2UiUj,ucart)
   Fresid(:)=Fresid(:)-Forces_TDEP(:)
   Phi1       (:)=0.d0
   Phi1Ui     (:)=0.d0
   Phi2%SR  (:,:)=0.d0
   if (Invar%loto) then
     Phi2%LR  (:,:)=0.d0
     Phi2%Tot (:,:)=0.d0
   end if
   Phi2UiUj   (:)=0.d0
   Forces_TDEP(:)=0.d0
 end if
!LOTO

 do istep=1,Invar%my_nstep
    do iatom=1,Invar%natom
       do ii=1,3
          Fresid(ii+3*(iatom-1)+3*Invar%natom*(istep-1))=Fresid(ii+3*(iatom-1)+3*Invar%natom*(istep-1))*&
&                                                        Invar%weights(istep)
       end do
    end do
 end do

 ABI_MALLOC(CoeffMoore%fcoeff,(3*natom*Invar%my_nstep,ntotcoeff)); CoeffMoore%fcoeff(:,:)=0.d0
 if (Invar%readifc.ne.1) then
   call tdep_calc_phi1fcoeff(CoeffMoore,Invar,proj1st,Shell1at,Sym)
   call tdep_calc_phi2fcoeff(CoeffMoore,Invar,proj2nd,Shell2at,Sym,ucart)
 end if
 if (Invar%order.ge.3) then
   ABI_MALLOC(Phi3_ref,(3,3,3,Shell3at%nshell))  ; Phi3_ref(:,:,:,:)=0.d0
   call tdep_calc_phi3fcoeff(CoeffMoore,Invar,proj3rd,Shell3at,Sym,ucart)
 end if
 if (Invar%order.eq.4) then
   ABI_MALLOC(Phi4_ref,(3,3,3,3,Shell4at%nshell)); Phi4_ref(:,:,:,:,:)=0.d0
   call tdep_calc_phi4fcoeff(CoeffMoore,Invar,proj4th,Shell4at,Sym,ucart)
 end if
 if (Invar%order.eq.2) then
 call tdep_calc_constraints(CoeffMoore,distance,Invar,MPIdata,Shell1at%nshell,Shell2at%nshell,&
&                           Shell3at%nshell,Shell4at%nshell,proj1st,proj2nd,&
&                           Shell1at,Shell2at,Sym)
 else if (Invar%order.eq.3) then
 call tdep_calc_constraints(CoeffMoore,distance,Invar,MPIdata,Shell1at%nshell,Shell2at%nshell,&
&                           Shell3at%nshell,Shell4at%nshell,proj1st,proj2nd,&
&                           Shell1at,Shell2at,Sym,proj3rd,Shell3at)
 else if (Invar%order.eq.4) then
 call tdep_calc_constraints(CoeffMoore,distance,Invar,MPIdata,Shell1at%nshell,Shell2at%nshell,&
&                           Shell3at%nshell,Shell4at%nshell,proj1st,proj2nd,&
&                           Shell1at,Shell2at,Sym,proj3rd,Shell3at,proj4th,Shell4at)
 end if

!==========================================================================================
!============= Compute the pseudo inverse using the Moore-Penrose method ==================
!==========================================================================================
 ABI_MALLOC(MP_coeff,(ntotcoeff,1)); MP_coeff(:,:)=0.d0
!==========================================================================================
!=================== If all the Orders are solved simultaneously ==========================
!==========================================================================================
 if (Invar%together.eq.1) then
   write(stdout,*) '############### (Solve simultaneously all the orders) #######################'
   if (Invar%readifc.ne.1) then
     call tdep_calc_MoorePenrose(CoeffMoore,Fresid,0,Invar,MP_coeff,MPIdata)
   end if
   ABI_MALLOC(Phi1_coeff,(ncoeff1st,1)); Phi1_coeff(:,:)=MP_coeff(1:ncoeff1st,:)
   ABI_MALLOC(Phi2_coeff,(ncoeff2nd,1)); Phi2_coeff(:,:)=MP_coeff(ncoeff1st+1:ncoeff1st+ncoeff2nd,:)
   if (Invar%readifc.ne.1) then
     call tdep_calc_phi1(Invar,ncoeff1st,proj1st,Phi1_coeff,Phi1   ,Shell1at,Sym)
     call tdep_calc_phi2(Invar,ncoeff2nd,proj2nd,Phi2_coeff,Phi2%SR,Shell2at,Sym)
   end if
   ABI_FREE(proj1st)
   ABI_FREE(proj2nd)
   ABI_FREE(Phi1_coeff)
   ABI_FREE(Phi2_coeff)
   call tdep_calc_ftot2(Forces_TDEP,Invar,Phi1,Phi1Ui,Phi2%SR,Phi2UiUj,ucart)
   if (Invar%order.ge.3) then
     ABI_MALLOC(Phi3_coeff,(ncoeff3rd,1)); Phi3_coeff(:,:)=0.d0
     Phi3_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+1:ncoeff1st+ncoeff2nd+ncoeff3rd,:)
     call tdep_calc_phi3ref(ncoeff3rd,proj3rd,Phi3_coeff,Phi3_ref,Shell3at)
     ABI_FREE(proj3rd)
     ABI_FREE(Phi3_coeff)
     call tdep_calc_ftot3(Forces_TDEP,Invar,Phi3_ref,Phi3UiUjUk,Shell3at,ucart,Sym)
   end if
   if (Invar%order.ge.4) then
     ABI_MALLOC(Phi4_coeff,(ncoeff4th,1)); Phi4_coeff(:,:)=0.d0
     Phi4_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+ncoeff3rd+1:ntotcoeff,:)
     call tdep_calc_phi4ref(ncoeff4th,proj4th,Phi4_coeff,Phi4_ref,Shell4at)
     ABI_FREE(proj4th)
     ABI_FREE(Phi4_coeff)
     call tdep_calc_ftot4(Forces_TDEP,Invar,Phi4_ref,Phi4UiUjUkUl,Shell4at,ucart,Sym)
   end if

!==========================================================================================
!=================== If all the Orders are solved successively ============================
!==========================================================================================
!ATTENTION : Le LR semble enleve a l'ordre 2 mais pas a l'ordre 3 et 4. On repart de
!            Forces_MD tout en bas et pas de Fresid (car les ordres 2 et 3 sont supprimes au
 else if (Invar%together.eq.0) then
   write(stdout,*) '################## (Solve successively each order) ##########################'
   do ii=1,Invar%order-1
     write(stdout,*) ' For order=',ii+1
     if (Invar%readifc.eq.1) cycle
     call tdep_calc_MoorePenrose(CoeffMoore,Fresid,ii,Invar,MP_coeff,MPIdata)
     if (ii.eq.1) then
       ABI_MALLOC(Phi1_coeff,(ncoeff1st,1)); Phi1_coeff(:,:)=MP_coeff(1:ncoeff1st,:)
       ABI_MALLOC(Phi2_coeff,(ncoeff2nd,1)); Phi2_coeff(:,:)=MP_coeff(ncoeff1st+1:ncoeff1st+ncoeff2nd,:)
       if (Invar%readifc.ne.1) then
         call tdep_calc_phi1(Invar,ncoeff1st,proj1st,Phi1_coeff,Phi1   ,Shell1at,Sym)
         call tdep_calc_phi2(Invar,ncoeff2nd,proj2nd,Phi2_coeff,Phi2%SR,Shell2at,Sym)
       end if
       ABI_FREE(proj1st)
       ABI_FREE(proj2nd)
       ABI_FREE(Phi1_coeff)
       ABI_FREE(Phi2_coeff)
       call tdep_calc_ftot2(Forces_TDEP,Invar,Phi1,Phi1Ui,Phi2%SR,Phi2UiUj,ucart)
     else if (ii.eq.2) then
       ABI_MALLOC(Phi3_coeff,(ncoeff3rd,1)); Phi3_coeff(:,:)=0.d0
       Phi3_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+1:ncoeff1st+ncoeff2nd+ncoeff3rd,:)
       call tdep_calc_phi3ref(ncoeff3rd,proj3rd,Phi3_coeff,Phi3_ref,Shell3at)
       ABI_FREE(proj3rd)
       ABI_FREE(Phi3_coeff)
       call tdep_calc_ftot3(Forces_TDEP,Invar,Phi3_ref,Phi3UiUjUk,Shell3at,ucart,Sym)
     else if (ii.eq.3) then
       ABI_MALLOC(Phi4_coeff,(ncoeff4th,1)); Phi4_coeff(:,:)=0.d0
       Phi4_coeff(:,:)=MP_coeff(ncoeff1st+ncoeff2nd+ncoeff3rd+1:ntotcoeff,:)
       call tdep_calc_phi4ref(ncoeff4th,proj4th,Phi4_coeff,Phi4_ref,Shell4at)
       ABI_FREE(proj4th)
       ABI_FREE(Phi4_coeff)
       call tdep_calc_ftot4(Forces_TDEP,Invar,Phi4_ref,Phi4UiUjUkUl,Shell4at,ucart,Sym)
     end if
     do istep=1,Invar%my_nstep
       do iatom=1,Invar%natom
         do jj=1,3
            Fresid(jj+3*(iatom-1)+3*Invar%natom*(istep-1))=(Forces_MD(jj+3*(iatom-1)+3*Invar%natom*(istep-1)) -&
&                                                           Forces_TDEP(jj+3*(iatom-1)+3*Invar%natom*(istep-1)))*&
&                                                           Invar%weights(istep)
         end do ! ii
       end do ! iatom
     end do ! istep
   end do
 end if
 ABI_FREE(CoeffMoore%const)
 ABI_FREE(CoeffMoore%fcoeff)
 ABI_FREE(Fresid)
 ABI_FREE(MP_coeff)
 ABI_FREE(ucart)
 call tdep_destroy_shell(natom,1,Shell1at)

!LOTO
 if (Invar%loto) then
   Phi2%Tot=Phi2%LR+Phi2%SR
 end if
!LOTO

!==========================================================================================
!=================== Write the IFC and check the constraints ==============================
!==========================================================================================
 call tdep_write_phi1(Invar,Phi1)
 call tdep_write_phi2(distance,Invar,MPIdata,Phi2%SR,Shell2at)
 if (Invar%order.ge.3) then
   call tdep_write_phi3(distance,Invar,Phi3_ref,Shell3at,Sym)
 end if
 if (Invar%order.ge.4) then
   call tdep_write_phi4(distance,Invar,Phi4_ref,Shell4at,Sym)
 end if

 if (Invar%order.eq.2) then
   call tdep_check_constraints(distance,Invar,Phi2%SR,Phi1,Shell3at%nshell,Shell4at%nshell,Sym)
 else if (Invar%order.eq.3) then
   call tdep_check_constraints(distance,Invar,Phi2%SR,Phi1,Shell3at%nshell,Shell4at%nshell,Sym,&
&                              Phi3_ref,Shell3at)
 else if (Invar%order.eq.4) then
   call tdep_check_constraints(distance,Invar,Phi2%SR,Phi1,Shell3at%nshell,Shell4at%nshell,Sym,&
&                              Phi3_ref,Shell3at,Phi4_ref,Shell4at)
 end if
 ABI_FREE(Phi1)

!==========================================================================================
!===================== Compute the phonon spectrum, the DOS, ==============================
!=====================  the dynamical matrix and write them ===============================
!==========================================================================================
 write(stdout,*)' '
 write(stdout,*) '#############################################################################'
 write(stdout,*) '############## Compute the phonon spectrum, the DOS, ########################'
 write(stdout,*) '##############  the dynamical matrix and write them  ########################'
 write(stdout,*) '#############################################################################'
 call tdep_init_eigen2nd(Eigen2nd_MP,Invar%natom_unitcell,Qbz%nqbz)
!FB call tdep_init_eigen2nd(Eigen2nd_MP,Invar%natom_unitcell,Qbz%nqibz)
 call tdep_init_eigen2nd(Eigen2nd_path,Invar%natom_unitcell,Qpt%nqpt)
 call tdep_calc_phdos(Crystal,DDB,Eigen2nd_MP,Eigen2nd_path,Ifc,Invar,Lattice,MPIdata,natom,&
&                          natom_unitcell,Phi2,PHdos,Qbz,Qpt,Rlatt4abi,Shell2at,Sym)
 call tdep_destroy_shell(natom,2,Shell2at)
 ABI_FREE(Rlatt4Abi)

 ! Create a new DDB with the coarse q-point grid in the IBZ.
 call DDB%free()
 call Ifc%to_ddb(DDB,Crystal)
 call tdep_write_ddb(DDB,Crystal,Invar)

 write(stdout,'(a)') ' See the dij.dat, omega.dat and eigenvectors files'
 write(stdout,'(a)') ' See also the DDB file'

!==========================================================================================
!===================== Compute the elastic constants ======================================
!==========================================================================================
 call tdep_calc_elastic(Phi2%SR,distance,Invar,Lattice)
 call tdep_destroy_phi2(Phi2,Invar%loto)

!==========================================================================================
!=========== Compute U_0, the "free energy" and the forces (from the model) ===============
!==========================================================================================
 call tdep_calc_model(Forces_MD,Forces_TDEP,Invar,MPIdata,Phi1Ui,Phi2UiUj,&
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
 call tdep_calc_thermo(Invar,Lattice,MPIdata,PHdos,U0)
 call PHdos%free()


!==========================================================================================
!===================== CALCULATION OF THE 3rd ORDER =======================================
!==========================================================================================
 if (Invar%order>2) then
   if (MPIdata%iam_master) then
     call tdep_write_gruneisen(distance,Eigen2nd_path,Invar,Phi3_ref,Qpt,Rlatt_cart,Shell3at,Sym)
   end if
   call tdep_calc_alpha_gamma(distance,Eigen2nd_MP,Invar,Lattice,MPIdata,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)

!FB Begin Lifetime
!FB call tdep_calc_lifetime1(Crystal,distance,Eigen2nd_MP,Ifc,Invar,Lattice,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)
!FB End Lifetime
 end if

!==========================================================================================
!===================== End the calculation ================================================
!==========================================================================================

!Write acknowledgements
 call tdep_print_Aknowledgments(Invar)
 call flush_unit(stdout)

 call timein(tcpu, twall)
 tsec(1)=tcpu-tcpui; tsec(2)=twall-twalli

!Free memory
 ABI_FREE(distance)
 ABI_FREE(Rlatt_cart)
 call tdep_destroy_eigen2nd(Eigen2nd_path)
 call tdep_destroy_eigen2nd(Eigen2nd_MP)

 if (Invar%order>2) then
   call tdep_destroy_shell(natom,3,Shell3at)
   ABI_FREE(Phi3_ref)
   if (Invar%order.eq.4) then
     call tdep_destroy_shell(natom,4,Shell4at)
     ABI_FREE(Phi4_ref)
   end if
 end if

 call Ifc%free()
 call DDB%free()
 call Crystal%free()
 call tdep_destroy_sym(Sym)
 call tdep_destroy_qbz(Qbz)
 call tdep_destroy_qpt(Qpt)

 if (MPIdata%iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
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

 call tdep_destroy_invar(Invar)
 call tdep_destroy_mpidata(MPIdata)

!Memory analysis
 call abinit_doctor(trim(Invar%output_prefix), print_mem_report=print_mem_report)
 call flush_unit(stdlog)
 close(unit=stdout)
100 call xmpi_end()

 end program atdep
!!***
