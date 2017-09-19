!!****p* ABINIT/tdep
!! NAME
!! tdep
!!
!! FUNCTION
!! Calculations of phonons using molecular dynamic simulations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (FB,JB)
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
!!      Abinit routines, our own routines and lapack/blas ones
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program tdep

  use defs_basis
  use m_phonons
  use m_tdepdos
  use m_ifc
  use m_phij,        only : calc_phij_fcoeff, calc_phij_nbcoeff, build_phij, calc_dij
  use m_xmpi,        only : xmpi_init, xmpi_end
  use m_latt,        only : make_latt, Lattice_Variables_type
  use m_sym,         only : make_sym, Symetries_Variables_type
  use m_readwrite,   only : Aknowledgments, ReadEcho, Input_Variables_type
  use m_utils,       only : calc_MoorePenrose, MatchIdeal2Average, model
  use m_qpt,         only : make_qptpath, Qpoints_type
#ifdef HAVE_NETCDF
  use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep'
!End of the abilint section

  implicit none

  integer :: natom,jatom,natom_unitcell,nshell,ntotcoeff,ncoeff,iatcell,ishell,stdout
  integer, allocatable :: bond_ref(:,:,:),shell(:,:)
  double precision :: U0,DeltaFree_AH2
  double precision, allocatable :: ucart(:,:,:),proj(:,:,:)
  double precision, allocatable :: fcoeff(:,:),Phij_coeff(:,:),fcartij(:),Phij_NN(:,:)
  double precision, allocatable :: distance(:,:,:),Rlatt_cart(:,:,:),Rlatt4dos(:,:,:)
  type(phonon_dos_type) :: PHdos
  type(Input_Variables_type) :: InVar
  type(Lattice_Variables_type) :: Lattice
  type(Symetries_Variables_type) :: Sym
  type(Qpoints_type) :: Qpt
  type(ifc_type) :: Ifc

!==========================================================================================
!===================== Initialization & Reading  ==========================================
!==========================================================================================
  call xmpi_init()

! Read input values from the input.in input file
  call ReadEcho(InVar)

! Initialize basic quantities  
  natom         =InVar%natom
  natom_unitcell=InVar%natom_unitcell
  stdout        =InVar%stdout  

!==========================================================================================
!============== Define the ideal lattice and its symetries ================================
!==========================================================================================
! Define all the quantities needed to buid the lattice (rprim*, acell*, brav*...)
  call make_latt(InVar,Lattice)

! Compute all the symetries coming from the bravais lattice
  call make_sym(Invar,Lattice,Sym)

!==========================================================================================
!======== 1/ Determine ideal positions and distances ======================================
!======== 2/ Find the matching between the ideal and average ==============================
!========   (from the MD simulations) positions. ==========================================
!======== 3/ Find the symetry operation between the reference and image bonds =============
!======== 4/ Write output quantities needed to visualize the neighbouring distances =======
!==========================================================================================
  allocate(Rlatt4dos (3,natom_unitcell,natom))   ; Rlatt4dos (:,:,:)=0.d0
  allocate(distance(natom,natom,4))              ; distance(:,:,:)=0.d0
  allocate(Rlatt_cart(3,natom_unitcell,natom))   ; Rlatt_cart(:,:,:)=0.d0
  allocate(ucart(3,natom,InVar%nstep_remain))    ; ucart(:,:,:)=0.d0
  allocate(fcartij(3*natom*InVar%nstep_remain))  ; fcartij(:)=0.d0
  allocate(bond_ref(natom,natom,3))              ; bond_ref(:,:,:)=0
  call MatchIdeal2Average(bond_ref,distance,fcartij,InVar,Lattice,nshell,Rlatt_cart,Rlatt4dos,Sym,ucart)

!==========================================================================================
!========  Find the number of coefficients of the (3x3) Phij for a given shell ============
!==========================================================================================
  write(stdout,*) ' '
  write(stdout,*) '#############################################################################'
  write(stdout,*) '################# The symetry operations (connecting the ####################'
  write(stdout,*) '################# interactions together) have been found. ###################'
  write(stdout,*) '################ Now, find the number of coefficients for ###################'
  write(stdout,*) '########################## a reference interaction ##########################'
  write(stdout,*) '#############################################################################'
  allocate(shell(nshell,4)); shell(:,:)=0
  allocate(proj(9,9,nshell)) ; proj(:,:,:)=0.d0
  write(stdout,*) 'Number of shells=',nshell
  ntotcoeff=0
  ishell=0
  open(unit=16,file='nbcoeff.dat')
  do iatcell=1,natom_unitcell
    write(stdout,*) '>>>>>>>>>>>>>>>>>>>>>>>>>  For iatcell=',iatcell
    ncoeff=0
    do jatom=1,natom
      if ((bond_ref(iatcell,jatom,1).ne.iatcell).or.(bond_ref(iatcell,jatom,2).ne.jatom)) cycle
      ishell=ishell+1
      write(stdout,*) 'Shell number:',ishell 
      write(stdout,'(a,i5,a,i5,a,f16.10)') '  Between atom',iatcell,' and ',jatom,' the distance is=',distance(iatcell,jatom,1)
      call calc_phij_nbcoeff(distance,iatcell,InVar,ishell,jatom,ncoeff,nshell,proj,Sym)
      shell(ishell,1)=ncoeff
      shell(ishell,2)=ntotcoeff
      shell(ishell,3)=iatcell
      shell(ishell,4)=jatom
      ntotcoeff=ntotcoeff+ncoeff
      write(stdout,*)'  ntotcoeff=',ntotcoeff
      write(stdout,*) '============================================================================'
    end do !jatom
  end do !iatcell
  close(16)
    
!==========================================================================================
!============= Build fcoeff, needed for the Moore-Penrose method just below ===============
!==========================================================================================
  allocate(fcoeff(3*natom*InVar%nstep_remain,ntotcoeff)); fcoeff(:,:)=0.d0 
  call calc_phij_fcoeff(InVar,bond_ref,nshell,ntotcoeff,proj,shell,Sym,ucart,fcoeff)

!==========================================================================================
!============= Compute the pseudo inverse using the Moore-Penrose method ==================
!==========================================================================================
  allocate(Phij_coeff(ntotcoeff,1)); Phij_coeff(:,:)=0.d0
  call calc_MoorePenrose(fcartij,fcoeff,InVar,ntotcoeff,Phij_coeff)
  deallocate(fcoeff)

!==========================================================================================
!============= Reorganize the IFC coefficients into the whole Phij_NN matrix ==============
!==========================================================================================
  allocate(Phij_NN(3*natom,         3*natom)) ; Phij_NN(:,:)=0.d0
  call build_phij(distance,InVar,bond_ref,nshell,ntotcoeff,proj,Phij_coeff,Phij_NN,shell,Sym)
  deallocate(Phij_coeff)
  deallocate(bond_ref)
  deallocate(shell)

!==========================================================================================
!=========== Compute U_0, the "free energy" and the forces (from the model) ===============
!==========================================================================================
  call model(DeltaFree_AH2,fcartij,InVar,Phij_NN,ucart,U0)
  deallocate(fcartij)
  deallocate(ucart)

!==========================================================================================
!===== 1/ Initialize the Brillouin zone and compute the q-points path  Compute  ===========
!===== 2/ Compute the dynamical matrix ====================================================
!==========================================================================================
  call make_qptpath(InVar,Lattice,Qpt)
  call calc_dij(InVar,Lattice,Phij_NN,Qpt,Rlatt_cart)
  deallocate(Rlatt_cart)

!==========================================================================================
!===================== Compute the elastic constants ======================================
!==========================================================================================
  call elastic(Phij_NN,distance,InVar,Lattice)
  deallocate(distance)

!==========================================================================================
!===================== Compute the phonons density of states ==============================
!==========================================================================================
  call make_phdos(Phij_NN(1:3*natom_unitcell,1:3*natom),Ifc,InVar,Lattice,natom,natom_unitcell,PHdos,Qpt,Rlatt4dos,Sym)
  deallocate(Rlatt4dos)
  deallocate(Phij_NN)

!==========================================================================================
!===================== Compute the thermodynamical quantities =============================
!==========================================================================================
  call thermo(DeltaFree_AH2,InVar,PHdos,U0)

!==========================================================================================
!================= Write the last informations (aknowledgments...)  =======================
!==========================================================================================
  call Aknowledgments(InVar)
  call xmpi_end()

end program tdep
!!***
