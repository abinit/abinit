!!****m* ABINIT/m_results_gs
!! NAME
!!  m_results_gs
!!
!! FUNCTION
!!  This module provides the definition of the results_gs_type
!!  used to store results from GS calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2025 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_results_gs

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_energies
 use m_errors
 use m_yaml
 use m_crystal
 use m_stream_string
 use m_pair_list
 use m_nctk
 use netcdf

 use m_io_tools,      only : file_exists
 use m_fstrings,      only : sjoin
 use m_numeric_tools, only : get_trace
 use m_geometry,      only : stress_voigt_to_mat

 implicit none

 private
!!***

!!****t* m_results_gs/results_gs_type
!! NAME
!! results_gs_type
!!
!! FUNCTION
!! This structured datatype contains the results of a GS calculation :
!! energy and its decomposition, forces and their decompositions, stresses
!! and their decompositions
!!
!! SOURCE

 type, public :: results_gs_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom
   ! The number of atoms for this dataset

  integer :: nspden
   ! The number of spin-density components for this dataset

  integer :: nsppol
   ! The number of spin channels for this dataset

  integer :: ngrvdw
   ! Size of grvdw array
   ! Can be 0 (not allocated) or natom

  integer :: berryopt
   ! Store the Berry phase option to know whether to use pel and pion (0 means no, otherwise yes)

! Real (real(dp)) scalars

  real(dp) :: deltae
   ! change in energy (Hartree)

  real(dp) :: diffor
   ! maximal absolute value of changes in the components of force

  real(dp) :: nelect_extfpmd
   ! Contribution of the Extended FPMD model to the number of electrons for high temperature simulations

! All the energies are in Hartree, obtained "per unit cell".
  type(energies_type) :: energies
!!!  real(dp) :: eei      ! local pseudopotential energy (Hartree)
!!!  real(dp) :: eeig     ! sum of eigenvalue energy (Hartree)
!!!  real(dp) :: eew      ! Ewald energy (Hartree)
!!!  real(dp) :: ehart    ! Hartree part of total energy (Hartree)
!!!  real(dp) :: eii      ! pseudopotential core-core energy
!!!  real(dp) :: ek       ! kinetic energy (Hartree)
!!!  real(dp) :: enefield ! the term of the energy functional that depends
!!!                       ! explicitely on the electric field
!!!                       ! enefield = -ucvol*E*P
!!!  real(dp) :: enl      ! nonlocal pseudopotential energy (Hartree)
  real(dp) :: entropy  ! entropy (Hartree)
!!!  real(dp) :: enxc     ! exchange-correlation energy (Hartree)
!!!  real(dp) :: enxcdc   ! exchange-correlation double-counting energy (Hartree)
!!!  real(dp) :: epaw     ! PAW spherical energy (Hartree)
!!!  real(dp) :: epawdc   ! PAW spherical double-counting energy (Hartree)
  real(dp) :: etotal   ! total energy (Hartree)
                       ! for fixed occupation numbers (occopt==0,1,or 2):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl+PAW_spherical_part
                       ! for varying occupation numbers (occopt>=3):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl - tsmear*entropy +PAW_spherical_part
  real(dp) :: fermie   ! Fermi energy (Hartree)
  real(dp) :: fermih   ! Fermi energy (Hartree) for excited holes in case occopt 9
  real(dp) :: residm   ! maximum value for the residual over all bands, all k points,
                       !   and all spins (Hartree or Hartree**2, to be checked !)
  real(dp) :: res2     ! density/potential residual (squared)
  real(dp) :: vxcavg   ! Average of the exchange-correlation energy. The average
                       ! of the local psp pot and the Hartree pot is set to zero (due
                       ! to the usual problem at G=0 for Coulombic system, so vxcavg
                       ! is also the average of the local part of the Hamiltonian

! Real (real(dp)) arrays

  real(dp), allocatable :: fcart(:,:)
   ! fcart(3,natom)
   ! Cartesian forces (Hartree/Bohr)
   ! Note: unlike gred, this array has been corrected by enforcing
   ! the translational symmetry, namely that the sum of force
   ! on all atoms is zero.

  real(dp), allocatable :: gred(:,:)
   ! gred(3,natom)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp),allocatable :: gaps(:,:)
   ! gaps(3,nsppol)
   ! gaps(1,:) : fundamental gap
   ! gaps(2,:) : optical gap
   ! gaps(3,:) : "status" for each channel:
   !   0.0dp if the gap was not computed (because there are only valence bands);
   !   -1.0dp if the system (or spin-channel) is metallic;
   !   1.0dp if the gap was computed

  real(dp), allocatable :: grchempottn(:,:)
   ! grchempottn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the spatially-varying chemical potential

  real(dp), allocatable :: grcondft(:,:)
   ! grcondft(nspden,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the constrained DFT contribution

  real(dp), allocatable :: gresid(:,:)
   ! gresid(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the residual
   ! of the potential

  real(dp), allocatable :: grewtn(:,:)
   ! grewtn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the Ewald energy

  real(dp), allocatable :: grvdw(:,:)
   ! grvdw(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from
   ! Van der Waals DFT-D dispersion (hartree)
   ! ngrvdw can be 0 or natom

  real(dp), allocatable :: grxc(:,:)
   ! grxc(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the XC energy

  real(dp), allocatable :: intgres(:,:)
   ! intgres(nspden,natom)
   ! Derivative of the total energy with respect to changes of constraints, in constrained DFT.

  real(dp) :: pel(3)
   ! ucvol times the electronic polarization in reduced coordinates

  real(dp) :: pion(3)
   ! ucvol times the ionic polarization in reduced coordinates

  real(dp) :: extfpmd_eshift
   ! Energy shift factor of the Extended FPMD model for high temperature simulations

  real(dp) :: strten(6)
   ! Stress tensor in cartesian coordinates (Hartree/Bohr^3)
   ! 6 unique components of this symmetric 3x3 tensor:
   ! Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).

  real(dp), allocatable :: synlgr(:,:)
   ! synlgr(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the non-local energy
   ! The "sy" prefix refer to the fact that this gradient has been
   ! symmetrized.

 contains

  procedure :: yaml_write => results_gs_yaml_write
    ! Write the most important results in Yaml format.

 end type results_gs_type

!public procedures.
 public :: init_results_gs
 public :: init_results_gs_array
 public :: destroy_results_gs
 public :: destroy_results_gs_array
 public :: copy_results_gs
 public :: results_gs_ncwrite
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_results_gs/init_results_gs
!! NAME
!!  init_results_gs
!!
!! FUNCTION
!!  Init all (or part of) scalars and allocatables in a results_gs datastructure
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  nsppol=number of spin channels for this dataset
!!  only_part= --optional, default=false--
!!            if this flag is activated only the following parts of results_gs
!!            are initalized: all scalars, fcart,gred,strten
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs=<type(results_gs_type)>=results_gs datastructure
!!
!! SOURCE

subroutine init_results_gs(natom,nspden,nsppol,results_gs,only_part)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspden,nsppol
 logical,optional,intent(in) :: only_part
!arrays
 type(results_gs_type),intent(inout) :: results_gs
!Local variables-------------------------------
!scalars
 logical :: full_init

!************************************************************************

 !@results_gs_type

 full_init=.true.;if (present(only_part)) full_init=(.not.only_part)

 results_gs%berryopt=0
 results_gs%natom  =natom
 results_gs%ngrvdw =0
 results_gs%nspden =nspden
 results_gs%nsppol =nsppol

 results_gs%deltae =zero
 results_gs%diffor =zero
 results_gs%entropy=zero
 results_gs%etotal =zero
 results_gs%fermie =zero
 results_gs%fermih =zero
 results_gs%nelect_extfpmd=zero
 results_gs%residm =zero
 results_gs%res2   =zero
 results_gs%extfpmd_eshift=zero
 results_gs%vxcavg =zero

 call energies_init(results_gs%energies)

 results_gs%strten=zero
 ABI_MALLOC(results_gs%fcart,(3,natom))
 results_gs%fcart=zero
 ABI_MALLOC(results_gs%gred,(3,natom))
 results_gs%gred =zero
 ABI_MALLOC(results_gs%gaps,(3,nsppol))
 results_gs%gaps =zero
 ABI_MALLOC(results_gs%intgres,(nspden,natom))
 results_gs%intgres=zero

 if (full_init) then
   results_gs%pel=zero
   results_gs%pion=zero

   ABI_MALLOC(results_gs%grchempottn,(3,natom))
   results_gs%grchempottn=zero
   ABI_MALLOC(results_gs%grcondft,(3,natom))
   results_gs%grcondft=zero
   ABI_MALLOC(results_gs%gresid,(3,natom))
   results_gs%gresid=zero
   ABI_MALLOC(results_gs%grewtn,(3,natom))
   results_gs%grewtn=zero
   ABI_MALLOC(results_gs%grvdw,(3,natom))
   results_gs%grvdw=zero
   ABI_MALLOC(results_gs%grxc,(3,natom))
   results_gs%grxc  =zero
   ABI_MALLOC(results_gs%synlgr,(3,natom))
   results_gs%synlgr=zero
 end if

end subroutine init_results_gs
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/init_results_gs_array
!! NAME
!!  init_results_gs_array
!!
!! FUNCTION
!!  Init all (or part of) scalars and allocatables in a 2D-array of results_gs datastructures
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  nsppol=number of spin channels for this dataset
!!  only_part= --optional, default=false--
!!            if this flag is activated only the following parts of results_gs
!!            are initalized: all scalars, fcart,gred,strten
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure 2Darray
!!
!! SOURCE

subroutine init_results_gs_array(natom,nspden,nsppol,results_gs,only_part)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspden,nsppol
 logical,optional,intent(in) :: only_part
!arrays
 type(results_gs_type),intent(inout) :: results_gs(:,:)
!Local variables-------------------------------
!scalars
 integer :: ii,jj,results_gs_size1,results_gs_size2
 logical :: full_init
!arrays

!************************************************************************

 !@results_gs_type

 results_gs_size1=size(results_gs,1)
 results_gs_size2=size(results_gs,2)
 full_init=.true.;if (present(only_part)) full_init=(.not.only_part)

 if (results_gs_size1>0.and.results_gs_size2>0) then

   do ii=1,results_gs_size2
     do jj=1,results_gs_size1

       results_gs(jj,ii)%berryopt=0
       results_gs(jj,ii)%natom  =natom
       results_gs(jj,ii)%ngrvdw =0
       results_gs(jj,ii)%nspden =nspden
       results_gs(jj,ii)%nsppol =nsppol

       results_gs(jj,ii)%deltae =zero
       results_gs(jj,ii)%diffor =zero
       results_gs(jj,ii)%entropy=zero
       results_gs(jj,ii)%etotal =zero
       results_gs(jj,ii)%fermie =zero
       results_gs(jj,ii)%fermih =zero
       results_gs(jj,ii)%nelect_extfpmd=zero
       results_gs(jj,ii)%residm =zero
       results_gs(jj,ii)%res2   =zero
       results_gs(jj,ii)%extfpmd_eshift=zero
       results_gs(jj,ii)%vxcavg =zero

       call energies_init(results_gs(jj,ii)%energies)

       results_gs(jj,ii)%strten=zero
       ABI_MALLOC(results_gs(jj,ii)%fcart,(3,natom))
       results_gs(jj,ii)%fcart=zero
       ABI_MALLOC(results_gs(jj,ii)%gred,(3,natom))
       results_gs(jj,ii)%gred =zero
       ABI_MALLOC(results_gs(jj,ii)%gaps,(3,nsppol))
       results_gs(jj,ii)%gaps =zero
       ABI_MALLOC(results_gs(jj,ii)%intgres,(nspden,natom))
       results_gs(jj,ii)%intgres =zero

       if (full_init) then
         results_gs(jj,ii)%pel=zero
         results_gs(jj,ii)%pion=zero
         ABI_MALLOC(results_gs(jj,ii)%grchempottn,(3,natom))
         results_gs(jj,ii)%grchempottn=zero
         ABI_MALLOC(results_gs(jj,ii)%grcondft,(3,natom))
         results_gs(jj,ii)%grcondft=zero
         ABI_MALLOC(results_gs(jj,ii)%gresid,(3,natom))
         results_gs(jj,ii)%gresid=zero
         ABI_MALLOC(results_gs(jj,ii)%grewtn,(3,natom))
         results_gs(jj,ii)%grewtn=zero
         ABI_MALLOC(results_gs(jj,ii)%grxc,(3,natom))
         results_gs(jj,ii)%grxc  =zero
         ABI_MALLOC(results_gs(jj,ii)%grvdw,(3,results_gs(jj,ii)%ngrvdw))
         results_gs(jj,ii)%grvdw  =zero
         ABI_MALLOC(results_gs(jj,ii)%synlgr,(3,natom))
         results_gs(jj,ii)%synlgr=zero
       end if

     end do
   end do
 end if

end subroutine init_results_gs_array
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/destroy_results_gs
!! NAME
!!  destroy_results_gs
!!
!! FUNCTION
!!  Clean and destroy a results_gs datastructure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure
!!
!! SOURCE

subroutine destroy_results_gs(results_gs)

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(inout) :: results_gs

!************************************************************************

 !@results_gs_type

 results_gs%natom =0
 results_gs%ngrvdw=0
 results_gs%nspden=0
 results_gs%nsppol=0
 results_gs%berryopt=0

 ABI_SFREE(results_gs%fcart)
 ABI_SFREE(results_gs%gred)
 ABI_SFREE(results_gs%gaps)
 ABI_SFREE(results_gs%grcondft)
 ABI_SFREE(results_gs%gresid)
 ABI_SFREE(results_gs%grewtn)
 ABI_SFREE(results_gs%grchempottn)
 ABI_SFREE(results_gs%grvdw)
 ABI_SFREE(results_gs%grxc)
 ABI_SFREE(results_gs%intgres)
 ABI_SFREE(results_gs%synlgr)

end subroutine destroy_results_gs
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/destroy_results_gs_array
!! NAME
!!  destroy_results_gs_array
!!
!! FUNCTION
!!  Clean and destroy a 2D-array of results_gs datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure 2D-array
!!
!! SOURCE

subroutine destroy_results_gs_array(results_gs)

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(inout) :: results_gs(:,:)
!Local variables-------------------------------
!scalars
 integer :: ii,jj,results_gs_size1,results_gs_size2

!************************************************************************

 !@results_gs_type

 results_gs_size1=size(results_gs,1)
 results_gs_size2=size(results_gs,2)
 if (results_gs_size1>0.and.results_gs_size2>0) then

   do ii=1,results_gs_size2
     do jj=1,results_gs_size1
       results_gs(jj,ii)%natom =0
       results_gs(jj,ii)%ngrvdw=0
       results_gs(jj,ii)%nspden=0
       results_gs(jj,ii)%nsppol=0
       results_gs(jj,ii)%berryopt=0

       ABI_SFREE(results_gs(jj,ii)%fcart)
       ABI_SFREE(results_gs(jj,ii)%gred)
       ABI_SFREE(results_gs(jj,ii)%gaps)
       ABI_SFREE(results_gs(jj,ii)%grchempottn)
       ABI_SFREE(results_gs(jj,ii)%grcondft)
       ABI_SFREE(results_gs(jj,ii)%gresid)
       ABI_SFREE(results_gs(jj,ii)%grewtn)
       ABI_SFREE(results_gs(jj,ii)%grvdw)
       ABI_SFREE(results_gs(jj,ii)%grxc)
       ABI_SFREE(results_gs(jj,ii)%intgres)
       ABI_SFREE(results_gs(jj,ii)%synlgr)
     end do
   end do

 end if

end subroutine destroy_results_gs_array
!!***
!----------------------------------------------------------------------

!!****f* m_results_gs/copy_results_gs
!! NAME
!!  copy_results_gs
!!
!! FUNCTION
!!  Copy a results_gs datastructure into another
!!
!! INPUTS
!!  results_gs_in=<type(results_gs_type)>=input results_gs datastructure
!!
!! OUTPUT
!!  results_gs_out=<type(results_gs_type)>=output results_gs datastructure
!!
!! SOURCE

subroutine copy_results_gs(results_gs_in,results_gs_out)

!Arguments ------------------------------------
!arrays
 class(results_gs_type),intent(in) :: results_gs_in
 type(results_gs_type),intent(inout) :: results_gs_out !vz_i

!Local variables-------------------------------
!scalars
 integer :: natom_in,natom_out,ngrvdw_in,nspden_in,nspden_out,nsppol_in,nsppol_out

!************************************************************************

 !@results_gs_type

 natom_in  =results_gs_in%natom
 natom_out =results_gs_out%natom
 ngrvdw_in =results_gs_in%ngrvdw
 nspden_in =results_gs_in%nspden
 nspden_out=results_gs_out%nspden
 nsppol_in =results_gs_in%nsppol
 nsppol_out=results_gs_out%nsppol

 if (natom_in>natom_out) then
   ABI_SFREE(results_gs_out%fcart)
   ABI_SFREE(results_gs_out%gred)
   ABI_SFREE(results_gs_out%grchempottn)
   ABI_SFREE(results_gs_out%grcondft)
   ABI_SFREE(results_gs_out%gresid)
   ABI_SFREE(results_gs_out%grewtn)
   ABI_SFREE(results_gs_out%grvdw)
   ABI_SFREE(results_gs_out%grxc)
   ABI_SFREE(results_gs_out%intgres)
   ABI_SFREE(results_gs_out%synlgr)

   if (allocated(results_gs_in%fcart))   then
     ABI_MALLOC(results_gs_out%fcart,(3,natom_in))
   end if
   if (allocated(results_gs_in%gred))    then
     ABI_MALLOC(results_gs_out%gred,(3,natom_in))
   end if
   if (allocated(results_gs_in%gresid))  then
     ABI_MALLOC(results_gs_out%gresid,(3,natom_in))
   end if
   if (allocated(results_gs_in%grchempottn))  then
     ABI_MALLOC(results_gs_out%grchempottn,(3,natom_in))
   end if
   if (allocated(results_gs_in%grcondft))  then
     ABI_MALLOC(results_gs_out%grcondft,(3,natom_in))
   end if
   if (allocated(results_gs_in%grewtn))  then
     ABI_MALLOC(results_gs_out%grewtn,(3,natom_in))
   end if
   if (allocated(results_gs_in%grvdw))  then
     ABI_MALLOC(results_gs_out%grvdw,(3,ngrvdw_in))
   end if
   if (allocated(results_gs_in%grxc))    then
     ABI_MALLOC(results_gs_out%grxc,(3,natom_in))
   end if
   if (allocated(results_gs_in%synlgr))  then
     ABI_MALLOC(results_gs_out%synlgr,(3,natom_in))
   end if
 end if

 if (nsppol_in>nsppol_out) then
   ABI_SFREE(results_gs_out%gaps)
   if (allocated(results_gs_in%gaps))    then
     ABI_MALLOC(results_gs_out%gaps,(3,nsppol_in))
   end if
 endif

 if (nspden_in>nspden_out .or. natom_in>natom_out) then
   ABI_SFREE(results_gs_out%intgres)
   if (allocated(results_gs_in%intgres))    then
     ABI_MALLOC(results_gs_out%intgres,(max(nspden_in,nspden_out),max(natom_in,natom_out)))
   end if
 endif


 results_gs_out%natom  =results_gs_in%natom
 results_gs_out%ngrvdw =results_gs_in%ngrvdw
 results_gs_out%nspden =results_gs_in%nspden
 results_gs_out%nsppol =results_gs_in%nsppol
 results_gs_out%berryopt=results_gs_in%berryopt
 results_gs_out%deltae =results_gs_in%deltae
 results_gs_out%diffor =results_gs_in%diffor
 results_gs_out%entropy=results_gs_in%entropy
 results_gs_out%etotal =results_gs_in%etotal
 results_gs_out%fermie =results_gs_in%fermie
 results_gs_out%fermih =results_gs_in%fermih
 results_gs_out%nelect_extfpmd=results_gs_in%nelect_extfpmd
 results_gs_out%residm =results_gs_in%residm
 results_gs_out%res2   =results_gs_in%res2
 results_gs_out%extfpmd_eshift=results_gs_in%extfpmd_eshift
 results_gs_out%vxcavg =results_gs_in%vxcavg

 call energies_copy(results_gs_in%energies,results_gs_out%energies)

 results_gs_out%pel(:)=results_gs_in%pel(:)
 results_gs_out%pion(:)=results_gs_in%pion(:)
 results_gs_out%strten(:)=results_gs_in%strten(:)

 if (allocated(results_gs_in%fcart))  results_gs_out%fcart(:,1:natom_in) =results_gs_in%fcart(:,1:natom_in)
 if (allocated(results_gs_in%gred))   results_gs_out%gred(:,1:natom_in)  =results_gs_in%gred(:,1:natom_in)
 if (allocated(results_gs_in%gaps))   results_gs_out%gaps(:,1:nsppol_in) =results_gs_in%gaps(:,1:nsppol_in)
 if (allocated(results_gs_in%grchempottn))&
&  results_gs_out%grchempottn(:,1:natom_in)=results_gs_in%grchempottn(:,1:natom_in)
 if (allocated(results_gs_in%grcondft)) results_gs_out%grcondft(:,1:natom_in)=results_gs_in%grcondft(:,1:natom_in)
 if (allocated(results_gs_in%gresid)) results_gs_out%gresid(:,1:natom_in)=results_gs_in%gresid(:,1:natom_in)
 if (allocated(results_gs_in%grewtn)) results_gs_out%grewtn(:,1:natom_in)=results_gs_in%grewtn(:,1:natom_in)
 if (allocated(results_gs_in%grxc))   results_gs_out%grxc(:,1:natom_in)  =results_gs_in%grxc(:,1:natom_in)
 if (allocated(results_gs_in%intgres))results_gs_out%intgres(1:nspden_in,1:natom_in)  =results_gs_in%intgres(1:nspden_in,1:natom_in)
 if (allocated(results_gs_in%synlgr)) results_gs_out%synlgr(:,1:natom_in)=results_gs_in%synlgr(:,1:natom_in)
 if (allocated(results_gs_in%grvdw).and.ngrvdw_in>0) then
   results_gs_out%grvdw(:,1:ngrvdw_in)=results_gs_in%grvdw(:,1:ngrvdw_in)
 end if

end subroutine copy_results_gs
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/results_gs_ncwrite
!! NAME
!! results_gs_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!!  ncid=NC file handle
!!  ecut, pawecutdg= Input cutoff energies in Ha.
!!
!! OUTPUT
!!
!! SOURCE

integer function results_gs_ncwrite(res, ncid, ecut, pawecutdg) result(ncerr)

!Arguments ------------------------------------
!scalars
 class(results_gs_type),intent(in) :: res
 integer,intent(in) :: ncid
 real(dp),intent(in) :: ecut,pawecutdg
! *************************************************************************

 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 !FIXME: do not handle k_dependent = 1
 ncerr = nctk_def_dims(ncid, [nctkdim_t("six", 6), nctkdim_t("number_of_cartesian_dimensions", 3), &
   nctkdim_t("number_of_atoms", res%natom), nctkdim_t("number_of_spins", res%nsppol)], defmode=.True.)
 NCF_CHECK(ncerr)

! Define variables.
! scalars passed in input (not belonging to results_gs) as well as scalars defined in results_gs
 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   "ecut", "pawecutdg", "deltae", "diffor", "entropy", "etotal", "fermie", "fermih",&
&  "nelect_extfpmd", "residm", "res2", "extfpmd_eshift"])
 NCF_CHECK(ncerr)

 ! arrays
 !
 ! Note: unlike gred, this array has been corrected by enforcing
 ! the translational symmetry, namely that the sum of force on all atoms is zero.

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('cartesian_forces', "dp", "number_of_cartesian_dimensions, number_of_atoms"),&
   nctkarr_t('cartesian_stress_tensor', "dp", 'six')])
 NCF_CHECK(ncerr)

 ! In case of a Berry phase calculation output the polarization
 if (res%berryopt/=0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('reduced_electronic_polarization', "dp", "number_of_cartesian_dimensions"),&
     nctkarr_t('reduced_ionic_polarization', "dp", "number_of_cartesian_dimensions")])
   NCF_CHECK(ncerr)
 end if

! Write data
! Write variables
 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
&  'ecut', 'pawecutdg', 'deltae', 'diffor', 'entropy', 'etotal', 'fermie', 'fermih',&
&  'nelect_extfpmd', 'residm', 'res2', 'extfpmd_eshift'],&
&  [ecut, pawecutdg, res%deltae, res%diffor, res%entropy, res%etotal, res%fermie, res%fermih,&
&  res%nelect_extfpmd, res%residm, res%res2, res%extfpmd_eshift],&
&  datamode=.True.)
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("cartesian_forces"), res%fcart))
 NCF_CHECK(nf90_put_var(ncid, vid("cartesian_stress_tensor"), res%strten))

 if (res%berryopt/=0) then
   NCF_CHECK(nf90_put_var(ncid, vid("reduced_electronic_polarization"), res%pel))
   NCF_CHECK(nf90_put_var(ncid, vid("reduced_ionic_polarization"), res%pion))
 end if

! Add energies
 call energies_ncwrite(res%energies, ncid)

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function results_gs_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/results_gs_yaml_write
!! NAME
!! results_gs_yaml_write
!!
!! FUNCTION
!!  Write results_gs in yaml format to unit.
!!
!! INPUTS
!!  results <type(results_gs_type)>=miscellaneous information about the system after ground state computation
!!  unit= unit of output file
!!  [cryst]: optional Crystal structure
!!  [info]: optional info for the final document
!!  [occopt]: optional Input variable occopt
!!  [with_conv]: optional True if the convergence dictionary with residuals and diffs should be written.
!!
!! SOURCE
subroutine results_gs_yaml_write(results, unit, cryst, info, occopt, with_conv)

 class(results_gs_type),intent(in) :: results
 integer,intent(in) :: unit
 type(crystal_t),intent(in),optional :: cryst
 integer, intent(in),optional :: occopt
 logical,intent(in),optional :: with_conv
 character(len=*),intent(in),optional :: info

!Local variables-------------------------------
 integer,parameter :: width=10
 integer :: ii
 type(yamldoc_t) :: ydoc
!arrays
 real(dp) :: strten(3,3), abc(3), fnorms(results%natom)
 character(len=2) :: species(results%natom)

!************************************************************************

 if (unit == dev_null) return

 if (present(info)) then
   ydoc = yamldoc_open('ResultsGS', info=info, width=width)
 else
   ydoc = yamldoc_open('ResultsGS', width=width)
 end if

 ! Write lattice parameters
 if (present(cryst)) then
  call ydoc%add_real2d('lattice_vectors', cryst%rprimd, real_fmt="(f11.7)")
  !ori abc = [(sqrt(sum(cryst%rprimd(:, ii) ** 2)), ii=1,3)]
  !replace the implicit loop by an explicit one
  !workaround works with both ifort and ifx on oneapi 2024
  do ii=1,3
    abc(ii) = sqrt(sum(cryst%rprimd(:, ii) ** 2))
  end do
  call ydoc%add_real1d('lattice_lengths', abc, real_fmt="(f10.5)")
  call ydoc%add_real1d('lattice_angles', cryst%angdeg, real_fmt="(f7.3)", comment="degrees, (23, 13, 12)")
  call ydoc%add_real('lattice_volume', cryst%ucvol + tol10, real_fmt="(es15.7)")
 endif

 ! Write convergence degree.
 ! It seems there's a portability problem for residm computed with nstep = 0 and iscf -3
 ! because one may get very small value e.g. 7.91-323. residm with nstep > 0 are OK though
 ! so print zero if residm < tol30 or allow the caller not to write the convergence dict.
 if(present(with_conv))then
   if (with_conv) then
     call ydoc%add_reals( &
      "deltae, res2, residm, diffor", &
       [results%deltae, results%res2, merge(results%residm, zero, results%residm > tol30), results%diffor], &
       real_fmt="(es10.3)", dict_key="convergence")
   else
     call ydoc%set_keys_to_string("deltae, res2, residm, diffor", "null", dict_key="convergence")
   end if
 else
   call ydoc%set_keys_to_string("deltae, res2, residm, diffor", "null", dict_key="convergence")
 endif

 ! Write energies.
 if (present(occopt))then
   if (occopt == 9) then
     call ydoc%add_reals("etotal, entropy, fermie, fermih", [results%etotal, results%entropy, results%fermie, results%fermih])
   else
     call ydoc%add_reals("etotal, entropy, fermie", [results%etotal, results%entropy, results%fermie])
   endif
 else
   call ydoc%add_reals("etotal, entropy, fermie", [results%etotal, results%entropy, results%fermie])
 endif

 ! Cartesian stress tensor and forces.
 call stress_voigt_to_mat(results%strten, strten)
 if (strten(1,1) /= MAGIC_UNDEF) then
   call ydoc%add_real2d('cartesian_stress_tensor', strten, comment="hartree/bohr^3")
   call ydoc%add_real('pressure_GPa', - get_trace(strten) * HaBohr3_GPa / three, real_fmt="(es12.4)")
 else
   call ydoc%set_keys_to_string("cartesian_stress_tensor, pressure_GPa", "null")
 end if

 if(present(cryst))then
   species = [(cryst%symbol_iatom(ii), ii=1,cryst%natom)]
   call ydoc%add_real2d('xred', cryst%xred, slist=species, real_fmt="(es12.4)")
 endif

 if (results%fcart(1,1) /= MAGIC_UNDEF) then
   call ydoc%add_real2d('cartesian_forces', results%fcart, comment="hartree/bohr")
   fnorms = [(sqrt(sum(results%fcart(:, ii) ** 2)), ii=1,results%natom)]
   ! Write force statistics
   call ydoc%add_reals('min, max, mean', &
     values=[minval(fnorms), maxval(fnorms), sum(fnorms) / results%natom], dict_key="force_length_stats")
 else
   ! Set entries to null (python None) to facilitate life to the parsing routines!
   call ydoc%add_string('cartesian_forces', "null")
   call ydoc%set_keys_to_string("min, max, mean", "null", dict_key="force_length_stats")
 end if

 call ydoc%write_and_free(unit)

end subroutine results_gs_yaml_write
!!***

!----------------------------------------------------------------------

END MODULE m_results_gs
!!***
