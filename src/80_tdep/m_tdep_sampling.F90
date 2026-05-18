!!****m* ABINIT/m_tdep_sampling
!! NAME
!!  m_tdep_sampling
!!
!! FUNCTION
!!  This module contains the TDEP Sampling data type
!!  which holds the set of configurations from which the IFC will be fit.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2026 ABINIT group (GA,FB,JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_sampling

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_io_tools
 use m_abihist,          only : abihist
 use m_tdep_dataset,     only : atdep_dataset_type, MPI_enreg_type
 use m_tdep_latt,        only : Lattice_type, tdep_make_inbox
 use m_tdep_sym,         only : Symmetries_type, tdep_SearchS_1at

 implicit none

 type tdep_Sampling_type

   integer :: natom
   ! Number of atoms in the supercell

   integer :: natom_unitcell
   ! Number of atoms in the unitcell

   integer :: my_nstep
   ! Number of MD steps held locally

   integer :: nstep_tot
   ! Total of MD steps

   integer, allocatable :: typat_unitcell(:)
   ! typat_unitcell(natom_unitcell)
   ! Atom type in the unitcell.

   integer, allocatable :: typat(:)
   ! typat(natom)
   ! Atom type in the supercell.

   double precision, allocatable :: xred_unitcell(:,:)
   ! xred_unitcell(3, natom_unitcell)
   ! Reduced equilibrium position of the atoms in the unitcell.

   double precision, allocatable :: xred_ideal(:,:)
   ! xred_ideal(3, natom)
   ! Reduced equilibrium position of the atoms in the supercell.

   double precision, allocatable :: xred(:,:,:)
   ! xred(3, natom, my_nstep)
   ! Reduced position of the atom at each step.

   double precision, allocatable :: fcart(:,:,:)
   ! fcart(3, natom, my_nstep)
   ! Cartesian forces at each step.

   double precision, allocatable :: etot(:)
   ! etot(my_nstep)
   ! Total energy at each step.

   double precision, allocatable :: weights(:)
   ! weights(my_nstep)
   ! The weight of each configuration for the fitting.
   ! By default, these will be 1 / nstep_tot.

   double precision, allocatable :: Rlatt_cart(:,:,:)
   ! Rlatt_cart(3, natom_unitcell, natom)
   ! Cartesian coordinate of the lattice vectors of the unitcell
   ! within the supercell, for each atom.

   double precision, allocatable :: Rlatt_scaled(:,:,:)
   ! Rlatt_scaled(3, natom_unitcell, natom)
   ! Rlatt_cart divided by acell_unitcell.

   double precision, allocatable :: ucart(:,:,:)
   ! ucart(3, natom, my_nstep)
   ! Cartesian displacements of the atoms with respect to their equilibrium
   ! positions at each time step.

   double precision, allocatable :: distance(:,:,:)
   ! distance(natom, natom, 4)
   ! Distance between the ideal positions of the atoms in the supercell,
   ! (norm, and cartesian components).

   double precision, allocatable :: Forces(:)
   ! Forces(3*natom*my_nstep)
   ! The cartesian forces for all configurations, as a flat array.

 end type tdep_Sampling_type

 public :: tdep_sampling_init_read
 public :: tdep_sampling_free
 public :: tdep_sampling_shift_xred
 public :: tdep_sampling_rotate
 public :: tdep_MatchIdeal2Average
 public :: tdep_write_xred_average

contains

!=====================================================================================================

 subroutine tdep_sampling_init_read(MD,Invar,MPIdata,Hist)

  type(tdep_Sampling_type), intent(inout) :: MD
  type(atdep_dataset_type), intent(in) :: Invar
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(abihist), intent(in) :: Hist

  integer :: this_istep,istep,iatom,jstep
  double precision :: tmp1,tmp2,tmp3

  MD%natom = Invar%natom
  MD%natom_unitcell = Invar%natom_unitcell
  MD%my_nstep = Invar%my_nstep
  MD%nstep_tot = Invar%nstep_tot

  ABI_CALLOC(MD%typat_unitcell, (MD%natom_unitcell))
  ABI_CALLOC(MD%typat, (MD%natom))
  ABI_CALLOC(MD%xred_unitcell, (3,MD%natom_unitcell))
  ABI_CALLOC(MD%xred_ideal, (3,MD%natom))
  ABI_CALLOC(MD%xred, (3,MD%natom,MD%my_nstep))
  ABI_CALLOC(MD%fcart, (3,MD%natom,MD%my_nstep))
  ABI_CALLOC(MD%etot, (MD%my_nstep))
  ABI_CALLOC(MD%weights, (MD%my_nstep))
  ABI_CALLOC(MD%distance, (MD%natom,MD%natom,4))
  ABI_CALLOC(MD%Rlatt_scaled, (3,MD%natom_unitcell,MD%natom))
  ABI_CALLOC(MD%Rlatt_cart, (3,MD%natom_unitcell,MD%natom))
  ABI_CALLOC(MD%ucart, (3,MD%natom,MD%my_nstep))
  ABI_CALLOC(MD%Forces, (3*MD%natom*MD%my_nstep))

  MD%typat_unitcell(:) = Invar%typat_unitcell(:)
  MD%xred_unitcell(:,:) = Invar%xred_unitcell(:,:)
  MD%typat(:) = Invar%typat(:)

! Read xred.dat, fcart.dat and etot.dat ASCII files or extract them from the HIST.nc netcdf file.
  write(Invar%stdout,'(a)') ' '
  this_istep=0
  jstep=0
  if (Invar%use_weights.eq.1) then
    open(unit=30,file=trim(Invar%input_prefix)//'_weights.dat')
  else if (Invar%use_weights.eq.0) then
    MD%weights=1.0d0/real(MD%nstep_tot)
  endif
  if (Invar%netcdf) then
    do istep=Invar%nstep_min,Invar%nstep_max
      if (mod(istep-Invar%nstep_min,Invar%slice).ne.0) then
        cycle
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) cycle
        this_istep=this_istep+1
        MD%xred(:,:,this_istep) =Hist%xred (:,:,istep)
        MD%fcart(:,:,this_istep)=Hist%fcart(:,:,istep)
        MD%etot(this_istep)     =Hist%etot     (istep)
      end if
    end do !istep
    write(Invar%stdout,'(a)') ' The positions, forces and energies are extracted from the NetCDF file: HIST.nc'
  else
    open(unit=60,file=trim(Invar%input_prefix)//'_fcart.dat')
    open(unit=50,file=trim(Invar%input_prefix)//'_xred.dat')
    open(unit=40,file=trim(Invar%input_prefix)//'_etot.dat')
    do istep=1,Invar%nstep_min-1
      if (Invar%use_weights.eq.1) then
         read(30,*) tmp1
      endif
      read(40,*) tmp1
      do iatom=1,MD%natom
        read(50,*) tmp1,tmp2,tmp3
        read(60,*) tmp1,tmp2,tmp3
      end do
    end do
    do istep=Invar%nstep_min,Invar%nstep_max
      if (mod(istep-Invar%nstep_min,Invar%slice).ne.0) then
        if (Invar%use_weights.eq.1) then
           read(30,*) tmp1
        endif
        read(40,*) tmp1
        do iatom=1,Invar%natom
          read(50,*) tmp1,tmp2,tmp3
          read(60,*) tmp1,tmp2,tmp3
        end do
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) then
          if (Invar%use_weights.eq.1) then
             read(30,*) tmp1
          endif
          read(40,*) tmp1
          do iatom=1,Invar%natom
            read(50,*) tmp1,tmp2,tmp3
            read(60,*) tmp1,tmp2,tmp3
          end do
        else
          this_istep=this_istep+1
          if (Invar%use_weights.eq.1) then
             read(30,*) MD%weights(this_istep)
          endif
          read(40,*) MD%etot(this_istep)
          do iatom=1,MD%natom
            read(50,*) MD%xred (1,iatom,this_istep),MD%xred (2,iatom,this_istep),MD%xred (3,iatom,this_istep)
            read(60,*) MD%fcart(1,iatom,this_istep),MD%fcart(2,iatom,this_istep),MD%fcart(3,iatom,this_istep)
          end do
        end if !my_step
      end if !slice
    end do !istep
    close(40)
    close(50)
    close(60)
    write(Invar%stdout,'(2a)') ' The positions, forces and energies are extracted from the ASCII files:',&
&                              ' xred.dat, fcart.dat & etot.dat'
  end if !netcdf
  if (Invar%use_weights.eq.1) then
    close(30)
  end if

 end subroutine tdep_sampling_init_read

!=====================================================================================================

 subroutine tdep_sampling_free(MD)

  type(tdep_Sampling_type), intent(inout) :: MD

  ABI_FREE(MD%typat_unitcell)
  ABI_FREE(MD%typat)
  ABI_FREE(MD%xred_unitcell)
  ABI_FREE(MD%xred_ideal)
  ABI_FREE(MD%xred)
  ABI_FREE(MD%fcart)
  ABI_FREE(MD%etot)
  ABI_FREE(MD%weights)
  ABI_FREE(MD%distance)
  ABI_FREE(MD%Rlatt_scaled)
  ABI_FREE(MD%Rlatt_cart)
  ABI_FREE(MD%ucart)
  ABI_FREE(MD%Forces)

 end subroutine tdep_sampling_free

!=====================================================================================================

! Shift xred to keep atoms in the same unit cell at each step.
subroutine tdep_sampling_shift_xred(MD,MPIdata)

  type(tdep_Sampling_type), intent(inout) :: MD
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer :: natom,ii,iatom,istep,ierr
  integer :: shift,shift_max,shift_best
  double precision :: xi, dist, best_dist
  double precision, allocatable :: x0(:,:)

  natom = MD%natom
  ABI_MALLOC(x0,(3,natom))

  ! Communicate xred at the first step
  x0(:,:) = zero
  if (MPIdata%my_step(1)) then
    x0(:,:) = MD%xred(:,:,1)
  end if
  call xmpi_sum(x0,MPIdata%comm_step,ierr)

  ! Shift xred from all steps in the same unitcell as the first step
  shift_max = 1
  do istep=1, MD%my_nstep
    do iatom=1,natom
      do ii=1,3
        best_dist = abs(MD%xred(ii,iatom,istep) - x0(ii,iatom))
        shift_best = 0
        do shift=-shift_max,shift_max
          xi = MD%xred(ii,iatom,istep) + shift
          dist = abs(xi - x0(ii,iatom))
          if (dist < best_dist) then
            best_dist = dist
            shift_best = shift
          end if
        end do
        MD%xred(ii,iatom,istep) = MD%xred(ii,iatom,istep) + shift_best
      end do
    end do
  end do

  ABI_FREE(x0)

end subroutine tdep_sampling_shift_xred

!=====================================================================================================

subroutine tdep_sampling_rotate(MD,rotation_cart)

  type(tdep_Sampling_type), intent(inout) :: MD
  integer :: iatom,istep
  double precision :: rotation_cart(3,3)

! Apply rotation to fcart
  do istep=1,MD%my_nstep
    do iatom=1,MD%natom
      MD%fcart(:,iatom,istep) = MATMUL(rotation_cart, MD%fcart(:,iatom,istep))
    end do
  end do

end subroutine tdep_sampling_rotate

!=====================================================================================================

!!****f* ABINIT/m_tdep_sampling/tdep_MatchIdeal2Average
!! NAME
!!  tdep_MatchIdeal2Average
!!
!! FUNCTION
!! Find the mapping between the atoms in the ideal (equilibrium) supercell,
!! and the atoms of the input moledular dynamics using their average positions.
!! Then compute the atom displacements with respect to the equilibrium positions
!! at each time step of the MD.
!!
!! INPUTS
!!  MD = TDEP Sampling object containing the input positions and forces.
!!  Invar = Input object containing the input variables.
!!  Lattice = Lattice object describing the ideal structure.
!!  Sym = Symetries object describing all the symmetry operations of the crystal.
!!  MPIdata = Info on MPI parallelism.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  The following quantities in MD are computed:
!!
!!  distance = Distance between the ideal positions of the atoms in the supercell,
!!             (norm, and cartesian components).
!!  Forces = Cartesian forces on the atoms at each time steps, as a flat array.
!!  ucart = Cartesian displacements of the atoms with respect to their equilibrium
!!          positions at each time step.
!!  Rlatt_cart = Cartesian coordinate of the lattice vectors of the unitcell
!!               within the supercell, for each atom.
!!               This array seems to have an extra dimension, for algorithmic simplicity.
!!  Rlatt_scaled = Rlatt_cart divided by acell_unitcell.
!!              These are used when reading an IFC file, to compare with the R vectors
!!              that are stored in the file.
!!
!! Some of the reduced positions of the atoms MD%xred are shifted by a supercell
!! lattice vector in order to re-center the crystal.
!!
!! NOTES
!!
!! SOURCE

 subroutine tdep_MatchIdeal2Average(MD,Invar,Lattice,Sym,MPIdata)

  type(tdep_Sampling_type),intent(inout) :: MD
  type(atdep_dataset_type),intent(inout) :: Invar
  type(Lattice_type),intent(in) :: Lattice
  type(Symmetries_type),intent(inout) :: Sym
  type(MPI_enreg_type),intent(in) :: MPIdata

  integer :: ii,jj,kk,max_ijk,iatcell,jatcell,iatom,jatom,eatom,fatom,istep
  integer :: iatom_ref,ierr
  integer :: ndir_match,natom_match
  double precision :: tmp(3),tmp1(3),tmp2(3),Rlatt(3),xred_tmp(3),rprimd_md_tmp(3,3),distance_tmp(3)
  double precision, allocatable :: dist_unitcell(:,:,:),xcart_average(:,:)
  double precision, allocatable :: fcart_tmp(:,:,:),ucart_tmp(:,:,:)
  double precision, allocatable  :: xred_average(:,:)
  double precision, allocatable  :: xred_center(:,:)
  double precision, allocatable  :: Rlatt_red (:,:,:)
  double precision, allocatable  :: xred_ideal(:,:)
! double precision, allocatable  :: distance_average(:,:,:)
  integer, allocatable  :: FromIdeal2Average(:)
  double precision, allocatable  :: xcart(:,:,:)
  double precision, allocatable  :: xcart_ideal(:,:)
  logical :: ok,must_shift,discard_R
  character(len=500) :: msg

  write(Invar%stdout,*)' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '###### Find the matching between ideal and average positions  ###############'
  write(Invar%stdout,*) '#############################################################################'

!==========================================================================================
!======== 1/ Determine ideal positions and distances ======================================
!==========================================================================================
  write(Invar%stdout,*)' Determine ideal positions and distances...'
! Define the bigbox with ideal positions
  ABI_CALLOC(Rlatt_red ,(3,MD%natom_unitcell,MD%natom))
  ABI_CALLOC(xred_ideal,(3,MD%natom))
  max_ijk=20
  iatom=1
  do ii=-max_ijk,max_ijk
    do jj=-max_ijk,max_ijk
      do kk=-max_ijk,max_ijk

        Rlatt(1)=real(ii-1)
        Rlatt(2)=real(jj-1)
        Rlatt(3)=real(kk-1)

        discard_R = .false.
        do iatcell=1,MD%natom_unitcell

          if (discard_R) cycle

!         Compute the reduced positions
          tmp(:) = Rlatt(:) + MD%xred_unitcell(:,iatcell)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,tmp(:),1,0.d0,xred_tmp(:),1)

!         If the first atom of the pattern is in the [0;1[ range then keep all the
!         atoms of the pattern (even if the others are outside the box). Else,
!         none are taken.
          if (iatcell==1) then
            if (minval(xred_tmp(:)).lt.0.d0.or.maxval(xred_tmp(:)).ge.(1.d0-tol12)) then
              discard_R = .true.
              cycle
            end if
          end if

          !GA: Why natom+1 ?
          if (iatom.gt.(MD%natom+1)) then
            ABI_ERROR('The number of atoms found in the bigbox exceeds natom' )
          end if

          xred_ideal(:,iatom) = xred_tmp(:)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,Rlatt(:),1,0.d0,Rlatt_red(:,1,iatom),1)
          iatom = iatom + 1
        end do
      end do
    end do
  end do

  if (iatom.lt.MD%natom+1) then
    ABI_ERROR('The number of atoms found in the big box is smaller than natom')
  end if

! Compute the distances between ideal positions in the SUPERcell
  do eatom=1,MD%natom
    do fatom=1,MD%natom
      tmp(:)=xred_ideal(:,fatom)-xred_ideal(:,eatom)
      call tdep_make_inbox(tmp,1,1d-4)
      rprimd_md_tmp(:,:)=Lattice%rprimd_md(:,:)
      distance_tmp(:)=MD%distance(eatom,fatom,2:4)
      call DGEMV('T',3,3,1.d0,rprimd_md_tmp,3,tmp,1,0.d0,distance_tmp,1)
      MD%distance(eatom,fatom,2:4)=distance_tmp(:)
      do ii=1,3
!       Remove the rounding errors before writing (for non regression testing purposes)
        if (abs(MD%distance(eatom,fatom,ii+1)).lt.tol8) MD%distance(eatom,fatom,ii+1)=zero
        MD%distance(eatom,fatom,1)=MD%distance(eatom,fatom,1)+(MD%distance(eatom,fatom,ii+1))**2
      end do
      MD%distance(eatom,fatom,1)=MD%distance(eatom,fatom,1)**0.5
      MD%distance(eatom,fatom,1)=tol12 * dint(MD%distance(eatom,fatom,1) / tol12)
    end do
  end do

! Compute the distances between ideal positions in the UNITcell
  ABI_MALLOC(dist_unitcell,(MD%natom_unitcell,MD%natom_unitcell,3)); dist_unitcell(:,:,:)=zero
  do iatcell=1,MD%natom_unitcell
    do jatcell=1,MD%natom_unitcell
      tmp(:) = xred_ideal(:,jatcell)-xred_ideal(:,iatcell)
      call tdep_make_inbox(tmp,1,tol8)
      dist_unitcell(iatcell,jatcell,:) = tmp(:)
    end do
  end do

!==========================================================================================
!======== 2/ Find the matching between the ideal and average ==============================
!========   (from the MD simulations) positions. ==========================================
!==========================================================================================
!  NOTE: - xred_center is used to find the matching with the ideal positions
!        - xred_average is used to compute the displacements (from MD trajectories)
!        The difference between those two is that xred_center will be shifted to bring
!        one of the average positions at the origin, for an easier comparison with
!        xred_ideal. Some shifts by a supercell lattice vector will be computed
!        from the difference between xred_center and xred_ideal, and those shifts
!        will be applied to xred_average and xred at all steps.

  write(Invar%stdout,*)' Compute average positions...'
  ABI_CALLOC(xred_average,(3,MD%natom))
  ABI_CALLOC(xred_center,(3,MD%natom))
! Average positions from MD (on nstep steps)
  do istep=1,MD%my_nstep
    do iatom=1,MD%natom
      xred_average(:,iatom)=xred_average(:,iatom)+MD%xred(:,iatom,istep)
    end do
  end do
  call xmpi_sum(xred_average,MPIdata%comm_step,ierr)
  xred_average(:,:) = xred_average(:,:) / real(MD%nstep_tot)

! Search the basis of atoms in the supercell
! in order to find iatom_ref
  write(Invar%stdout,*)' Search the unitcell basis of atoms in the MD trajectory...'
  ok=.false.
  xred_center(:,:) = xred_average(:,:)
  iatcell=1
  do iatom=1,MD%natom
    if (MD%typat(iatom).ne.MD%typat_unitcell(iatcell)) cycle
    natom_match = 0
    do jatom=1,MD%natom

      tmp(:)=xred_center(:,jatom)-xred_center(:,iatom)
      call tdep_make_inbox(tmp,1,Invar%tolinbox)

      do jatcell=1,MD%natom_unitcell
        if (MD%typat(jatom).ne.MD%typat_unitcell(jatcell)) cycle
        ndir_match = 0
        do ii=1,3
          if (abs(tmp(ii)-dist_unitcell(iatcell,jatcell,ii)).le.Invar%tolmotif) then
            ndir_match=ndir_match+1
          end if
        end do
        if (ndir_match==3) then
          natom_match = natom_match + 1
          exit
        end if
      end do
    end do
    if (natom_match.eq.MD%natom_unitcell) then
      iatom_ref = iatom
      ok=.true.
      exit
    else if (natom_match.gt.MD%natom_unitcell) then
      write(msg,'(5a)') 'Too many atoms match the unit cell.',ch10,&
                        'Perhaps the value of tolmotif is too large,',ch10,&
                        'or the value of tolinbox is too small.'
      ABI_ERROR(msg)
    endif
  end do
  if (.not.ok) then
    call tdep_write_xred_average(Invar,MPIdata,Lattice,xred_ideal,xred_center)
    write(msg,'(3a)') 'The basis of atoms written in input.in file does not appear in the MD trajectory.',ch10,&
                      'Perhaps, you can adjust the tolerance (tolmotif).'
    ABI_ERROR(msg)
  end if
  ABI_FREE(dist_unitcell)

  write(Invar%stdout,*)' Compare ideal and average positions using PBC...'
! Modification of xred and Rlatt tabs
! for averaged quantities: xred_center, xred_average, xred
! 1/ The "iatom_ref" atom is put in (0.0;0.0;0.0)
  tmp(:) = xred_center(:,iatom_ref)
  do jatom=1,MD%natom
    xred_center(:,jatom) = xred_center(:,jatom) - tmp(:)
  end do
! 2/ All the atoms are put in the range [-0.5;0.5[ (use of PBC)
  do jatom=1,MD%natom
    tmp(:)=xred_center(:,jatom)
    call tdep_make_inbox(tmp,1,Invar%tolinbox,xred_center(:,jatom))
    call tdep_make_inbox(tmp,1,Invar%tolinbox,xred_average(:,jatom))
    do istep=1,MD%my_nstep
      call tdep_make_inbox(tmp,1,Invar%tolinbox,MD%xred(:,jatom,istep))
    end do
  end do
! Modification of xred and Rlatt tabs
! for ideal quantities: Rlatt_red et xred_ideal
!   1/ The atom 1 is put in (0.0;0.0;0.0)
  tmp1(:)=xred_ideal(:,1)
  tmp2(:)=Rlatt_red(:,1,1)
  do jatom=1,MD%natom
    xred_ideal(:,jatom)=  xred_ideal(:,jatom)  -tmp1(:)
    Rlatt_red (:,1,jatom)=Rlatt_red (:,1,jatom)-tmp2(:)
  end do
! 2/ All the atoms are put in the range [-0.5;0.5[ (use of PBC)
  do jatom=1,MD%natom
    tmp(:)=xred_ideal(:,jatom)
    call tdep_make_inbox(tmp,1,tol8,xred_ideal(:,jatom))
    call tdep_make_inbox(tmp,1,tol8,Rlatt_red(:,1,jatom))
!FB      call tdep_make_inbox(Rlatt_red(:,1,jatom),1,tol8)
  end do

! When the multiplicity equals 1 along one direction, there is some trouble
! To clean!!!!!!!
  do ii=1,3
    if ((Invar%multiplicity(ii,ii).eq.1).and.(Invar%multiplicity(ii,mod(ii  ,3)+1).eq.0)&
&                                 .and.(Invar%multiplicity(ii,mod(ii+1,3)+1).eq.0)) then
      Rlatt_red(ii,1,:)=0.d0
      write(Invar%stdout,*) 'WARNING: multiplicity=1 for ii=',ii
    end if
  end do

! Define Rlatt for all the atoms in the basis (Rlatt_red varies as a function of iatcell)
  if (MD%natom_unitcell.gt.1) then
    do iatcell=2,MD%natom_unitcell
      Rlatt_red(:,iatcell,:)=Rlatt_red(:,1,:)
    end do
  end if
  do iatom=1,MD%natom
    do iatcell=1,MD%natom_unitcell
      tmp(:)=xred_ideal(:,iatom)-xred_ideal(:,iatcell)
      call tdep_make_inbox(tmp,1,tol8,Rlatt_red(:,iatcell,iatom))
    end do
  end do
  if (Invar%debug) then
    do iatcell=1,MD%natom_unitcell
      write(Invar%stdout,*) 'For iatcell=',iatcell
      do jatom=1,MD%natom
        write(Invar%stdout,'(a,i4,a,3(f16.10,1x))') 'For jatom=',jatom,', Rlatt=',Rlatt_red(1:3,iatcell,jatom)
      end do
    end do
  end if

! Matching between Ideal and Average positions: xred_ideal and xred_center
! Then, write them in the xred_average.xyz file.
  write(Invar%stdout,*)' Write the xred_average.xyz file with ideal and average positions...'
  ABI_CALLOC(FromIdeal2Average,(MD%natom))
  do iatom=1,MD%natom
    ok =.false.
    do jatom=1,MD%natom
      if (MD%typat(iatom).ne.MD%typat_unitcell(mod(jatom-1,MD%natom_unitcell)+1)) cycle
      must_shift=.false.
      ndir_match=0
      do ii=1,3
        if (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)).le.Invar%tolmatch) then
          ndir_match=ndir_match+1
        else if ((abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)-1.d0).le.Invar%tolmatch) &
&            .or.(abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)+1.d0).le.Invar%tolmatch)) then
          ndir_match=ndir_match+1
          must_shift=.true.
        endif
      end do
      if (ndir_match==3.and..not.must_shift) then
        FromIdeal2Average(jatom)=iatom
        ok=.true.
        exit
      else if (ndir_match==3.and.must_shift) then
        do ii=1,3
          if (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)-1.d0).le.Invar%tolmatch) then
            xred_center(ii,iatom)=xred_center(ii,iatom)-1d0
            xred_average(ii,iatom)=xred_average(ii,iatom)-1d0
            do istep=1,MD%my_nstep
              MD%xred(ii,iatom,istep)=MD%xred(ii,iatom,istep)-1d0
            end do
            FromIdeal2Average(jatom)=iatom
          else if (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)+1.d0).le.Invar%tolmatch) then
            xred_center(ii,iatom)=xred_center(ii,iatom)+1d0
            xred_average(ii,iatom)=xred_average(ii,iatom)+1d0
            do istep=1,MD%my_nstep
              MD%xred(ii,iatom,istep)=MD%xred(ii,iatom,istep)+1d0
            end do
            FromIdeal2Average(jatom)=iatom
          end if
        end do
        ok=.true.
        exit
      end if
    end do
    if (.not.ok) then
      write(Invar%stdlog,*) 'Problem to find the average position for iatom=',iatom
      write(Invar%stdlog,*) '  Reasons:'
      write(Invar%stdlog,*) '    1/ One atom jump to another equilibrium position'
      write(Invar%stdlog,*) '    2/ The system is no more solid'
      write(Invar%stdlog,*) '    3/ Perhaps, you can adjust the tolerance (tolmatch)'
      write(Invar%stdlog,*) '  xred_center=',(xred_center(ii,iatom),ii=1,3)
      do eatom=1,MD%natom
        write(Invar%stdlog,'(a,1x,3(f10.6,1x))') 'I',xred_ideal (:,eatom)
        write(Invar%stdlog,'(a,1x,3(f10.6,1x))') 'C',xred_center(:,eatom)
      end do
      ABI_ERROR('Problem to find the average position')
    end if
  end do

! WARNING: VERY IMPORTANT: The positions are displayed/sorted
! (and used in the following) according to ideal positions xred_ideal.
  call tdep_write_xred_average(Invar,MPIdata,Lattice,xred_ideal,xred_center,FromIdeal2Average)
  ABI_FREE(xred_center)

!====================================================================================
!====================== END OF REDUCED COORDINATES ==================================
!====================================================================================
! a/ Get cartesian coordinates from reduced ones
! b/ Compute ucart and fcart tabs
! c/ The atoms are sorted according the IDEAL arrangement
!    The correspondance function is contained in: FromIdeal2Average
!    WARNING : Consequently the arrangement of the xcart* tabs is not modified.
  write(Invar%stdout,*)' Compute cartesian coordinates and forces...'
  ABI_MALLOC(xcart        ,(3,MD%natom,MD%my_nstep)); xcart(:,:,:)=0.d0
  ABI_MALLOC(xcart_ideal  ,(3,MD%natom))               ; xcart_ideal(:,:)=0.d0
  ABI_MALLOC(xcart_average,(3,MD%natom))               ; xcart_average(:,:)=0.d0
  ABI_MALLOC(ucart_tmp    ,(3,MD%natom,MD%my_nstep)); ucart_tmp(:,:,:)=0.d0
  do iatom=1,MD%natom
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_md(:,:),3,xred_ideal  (:,iatom),1,0.d0,xcart_ideal  (:,iatom),1)
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_md(:,:),3,xred_average(:,iatom),1,0.d0,xcart_average(:,iatom),1)
    do iatcell=1,MD%natom_unitcell
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_md(:,:),3,Rlatt_red(:,iatcell,iatom),1,0.d0,MD%Rlatt_cart(:,iatcell,iatom),1)
    end do
  end do
  do istep=1,MD%my_nstep
    do iatom=1,MD%natom
      jatom = FromIdeal2Average(iatom)
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_md(:,:),3,MD%xred(:,jatom,istep),&
&                1,0.d0,xcart(:,jatom,istep),1)
      if (Invar%use_ideal_positions.eq.0) then
        ucart_tmp(:,iatom,istep) = xcart(:,jatom,istep) - xcart_average(:,jatom)
      else
        ucart_tmp(:,iatom,istep) = xcart(:,jatom,istep) - xcart_ideal(:,iatom)
      end if
    end do
  end do
  ABI_FREE(xred_average)
  ABI_FREE(xcart)
  ABI_FREE(xcart_ideal)
  ABI_FREE(xcart_average)

! Rearrangement of the fcart tabs in column --> MD%Forces
  ABI_CALLOC(fcart_tmp,(3,MD%natom,MD%my_nstep))
  do istep=1,MD%my_nstep
    do iatom=1,MD%natom
      fcart_tmp(:,iatom,istep) = MD%fcart(:,FromIdeal2Average(iatom),istep)
    end do
  end do
  do istep=1,MD%my_nstep
    do jatom=1,MD%natom
      do ii=1,3
        jj = ii + 3*(jatom-1) + 3*MD%natom*(istep-1)
        MD%Forces(jj) = fcart_tmp(ii,jatom,istep)
        MD%ucart(ii,jatom,istep) = ucart_tmp(ii,jatom,istep)
      enddo
    enddo
  enddo
  ABI_FREE(FromIdeal2Average)
  ABI_FREE(ucart_tmp)
  ABI_FREE(fcart_tmp)

! Define Rlatt_scaled, fulfilling the definition of mkphdos (ABINIT routine)
  do ii=1,3
    rprimd_md_tmp(ii,:) = Lattice%rprimd_md(ii,:) / Lattice%acell_unitcell(ii)
  end do
  do iatom=1,MD%natom
    do iatcell=1,MD%natom_unitcell
      call DGEMV('T',3,3,1.d0,rprimd_md_tmp,3,Rlatt_red(:,iatcell,iatom),1,0.d0,MD%Rlatt_scaled(:,iatcell,iatom),1)
    end do
  end do

! Find the symetry operation between 2 atoms
  call tdep_SearchS_1at(Invar,MPIdata,Sym,xred_ideal)
  MD%xred_ideal(:,:)=xred_ideal(:,:)
  ABI_FREE(xred_ideal)
  ABI_FREE(Rlatt_red)

 end subroutine tdep_MatchIdeal2Average

!====================================================================================================

 subroutine tdep_write_xred_average(Invar,MPIdata,Lattice,&
                                    xred_ideal,xred_center,&
                                    FromIdeal2Average)
  type(atdep_dataset_type), intent(in) :: Invar
  type(MPI_enreg_type),intent(in) :: MPIdata
  type(Lattice_type),intent(in) :: Lattice
  double precision,intent(in) :: xred_ideal(3,Invar%natom)
  double precision,intent(in) :: xred_center(3,Invar%natom)
  integer,intent(in),optional :: FromIdeal2Average(Invar%natom)

  integer :: unt
  !integer :: natom,natom_unitcell
  integer :: iatom,jatom,ii,jj
  !logical :: with_xcart
  integer,allocatable :: ideal2average(:)
  double precision :: rprimd(3,3)
  double precision :: xred_C(3),xred_I(3),xcart_C(3),xcart_I(3)

  if (MPIdata%iam_master) then

    rprimd(:,:) = Lattice%rprimd_md(:,:)

    ABI_MALLOC(ideal2average,(Invar%natom))
    ideal2average(:)=0
    if (present(FromIdeal2Average)) then
      ideal2average(:) = FromIdeal2Average(:)
    else
      do iatom=1,Invar%natom
        ideal2average(iatom) = iatom
      end do
    end if

    unt=31
    open(unit=unt,file=trim(Invar%output_prefix)//'_xred_average.xyz')
    write(unt,'(a,i4)') '# natom = ',Invar%natom
    write(unt,'(a,i4)') '# natom_unitcell = ',Invar%natom_unitcell
    write(unt,'(a,9(f4.1,1x))') '# multiplicity = ',((Lattice%multiplicity(ii,jj),jj=1,3),ii=1,3 )
    write(unt,'(a)') '#'

    write(unt,'(a1,1x,a8,2x,a6,2x,2(a5,30x))') '#', 'position', 'iatom', 'xred ', 'xcart'
    write(unt,'(a)')''

    xred_I = zero
    xred_C = zero
    do iatom=1,Invar%natom
      jatom = ideal2average(iatom)
      xred_I = xred_ideal (:,iatom)
      xred_C = xred_center(:,jatom)

      xcart_I(:)=zero
      xcart_C(:)=zero
      call DGEMV('T',3,3,1.d0,rprimd(:,:),3,xred_I,1,0.d0,xcart_I,1)
      call DGEMV('T',3,3,1.d0,rprimd(:,:),3,xred_C,1,0.d0,xcart_C,1)

      write(unt,'(2x,a6,4x,i6,2x,3(f10.6,1x),2x,3(f10.6,1x))')'Ideal ',iatom,xred_I,xcart_I
      write(unt,'(2x,a6,4x,i6,2x,3(f10.6,1x),2x,3(f10.6,1x))')'Center',jatom,xred_C,xcart_C
      write(unt,'(a)')''

    end do

    close(unt)
    ABI_FREE(ideal2average)
  end if

 end subroutine tdep_write_xred_average

!====================================================================================================

end module m_tdep_sampling
!!***
