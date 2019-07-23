!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_results_gs
!! NAME
!!  m_results_gs
!!
!! FUNCTION
!!  This module provides the definition of the results_gs_type
!!  used to store results from GS calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
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
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,      only : file_exists
 use m_fstrings,      only : sjoin
 use m_numeric_tools, only : get_trace
 use defs_abitypes,   only : dataset_type

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
   ! Note : unlike fred, this array has been corrected by enforcing
   ! the translational symmetry, namely that the sum of force
   ! on all atoms is zero.

  real(dp), allocatable :: fred(:,:)
   ! fred(3,natom)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp),allocatable :: gaps(:,:)
   ! gaps(3,nsppol)
   ! gaps(1,:) : fundamental gap
   ! gaps(2,:) : optical gap
   ! gaps(3,:) : "status" for each channel : 0.0dp if the gap was not computed
   !   (because there are only valence bands) ; -1.0dp if the system (or spin-channel) is metallic ; 1.0dp if the
   !   gap was computed

  real(dp), allocatable :: gresid(:,:)
   ! gresid(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the residual
   ! of the potential

  real(dp), allocatable :: grewtn(:,:)
   ! grewtn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the Ewald energy

  real(dp), allocatable :: grchempottn(:,:)
   ! grchempottn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the spatially-varying chemical potential

  real(dp), allocatable :: grvdw(:,:)
   ! grvdw(3,ngrvdw)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from
   ! Van der Waals DFT-D2 dispersion (hartree)
   ! ngrvdw can be 0 or natom

  real(dp), allocatable :: grxc(:,:)
   ! grxc(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the XC energy

  real(dp) :: pel(3)
   ! ucvol times the electronic polarization in reduced coordinates

  real(dp) :: pion(3)
   ! ucvol times the ionic polarization in reduced coordinates

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
!!            are initalized: all scalars, fcart,fred,strten
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs=<type(results_gs_type)>=results_gs datastructure
!!
!! PARENTS
!!      m_results_img,mover_effpot
!!
!! CHILDREN
!!      energies_copy,energies_ncwrite
!!
!! SOURCE

subroutine init_results_gs(natom,nsppol,results_gs,only_part)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsppol
 logical,optional,intent(in) :: only_part
!arrays
 type(results_gs_type),intent(inout) :: results_gs
!Local variables-------------------------------
!scalars
 logical :: full_init

!************************************************************************

 !@results_gs_type

 full_init=.true.;if (present(only_part)) full_init=(.not.only_part)

 results_gs%natom  =natom
 results_gs%ngrvdw =0
 results_gs%nsppol =nsppol
 results_gs%berryopt=zero
 results_gs%deltae =zero
 results_gs%diffor =zero
 results_gs%entropy=zero
 results_gs%etotal =zero
 results_gs%fermie =zero
 results_gs%residm =zero
 results_gs%res2   =zero
 results_gs%vxcavg =zero

 call energies_init(results_gs%energies)

 results_gs%strten=zero
 ABI_ALLOCATE(results_gs%fcart,(3,natom))
 results_gs%fcart=zero
 ABI_ALLOCATE(results_gs%fred,(3,natom))
 results_gs%fred =zero
 ABI_ALLOCATE(results_gs%gaps,(3,nsppol))
 results_gs%gaps =zero
 if (full_init) then
   results_gs%pel=zero
   results_gs%pion=zero
   ABI_ALLOCATE(results_gs%gresid,(3,natom))
   results_gs%gresid=zero
   ABI_ALLOCATE(results_gs%grewtn,(3,natom))
   results_gs%grewtn=zero
   ABI_ALLOCATE(results_gs%grchempottn,(3,natom))
   results_gs%grchempottn=zero
   ABI_ALLOCATE(results_gs%grxc,(3,natom))
   results_gs%grxc  =zero
   ABI_ALLOCATE(results_gs%synlgr,(3,natom))
   results_gs%synlgr=zero
   ABI_ALLOCATE(results_gs%grvdw,(3,results_gs%ngrvdw))
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
!!            are initalized: all scalars, fcart,fred,strten
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure 2Darray
!!
!! PARENTS
!!
!! CHILDREN
!!      energies_copy,energies_ncwrite
!!
!! SOURCE

subroutine init_results_gs_array(natom,nsppol,results_gs,only_part)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsppol
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

       results_gs(jj,ii)%natom  =natom
       results_gs(jj,ii)%ngrvdw =0
       results_gs(jj,ii)%nsppol =nsppol
       results_gs(jj,ii)%berryopt=zero
       results_gs(jj,ii)%deltae =zero
       results_gs(jj,ii)%diffor =zero
       results_gs(jj,ii)%entropy=zero
       results_gs(jj,ii)%etotal =zero
       results_gs(jj,ii)%fermie =zero
       results_gs(jj,ii)%residm =zero
       results_gs(jj,ii)%res2   =zero
       results_gs(jj,ii)%vxcavg =zero

       call energies_init(results_gs(jj,ii)%energies)

       results_gs(jj,ii)%strten=zero
       ABI_ALLOCATE(results_gs(jj,ii)%fcart,(3,natom))
       results_gs(jj,ii)%fcart=zero
       ABI_ALLOCATE(results_gs(jj,ii)%fred,(3,natom))
       results_gs(jj,ii)%fred =zero
       ABI_ALLOCATE(results_gs(jj,ii)%gaps,(3,nsppol))
       results_gs(jj,ii)%gaps =zero
       if (full_init) then
         results_gs(jj,ii)%pel=zero
         results_gs(jj,ii)%pion=zero
         ABI_ALLOCATE(results_gs(jj,ii)%gresid,(3,natom))
         results_gs(jj,ii)%gresid=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grewtn,(3,natom))
         results_gs(jj,ii)%grewtn=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grchempottn,(3,natom))
         results_gs(jj,ii)%grchempottn=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grxc,(3,natom))
         results_gs(jj,ii)%grxc  =zero
         ABI_ALLOCATE(results_gs(jj,ii)%synlgr,(3,natom))
         results_gs(jj,ii)%synlgr=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grvdw,(3,results_gs(jj,ii)%ngrvdw))
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
!! PARENTS
!!      m_results_img,mover_effpot
!!
!! CHILDREN
!!      energies_copy,energies_ncwrite
!!
!! SOURCE

subroutine destroy_results_gs(results_gs)

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(inout) :: results_gs
!Local variables-------------------------------

!************************************************************************

 !@results_gs_type

 results_gs%natom =0
 results_gs%ngrvdw=0
 results_gs%nsppol=0
 results_gs%berryopt=0
 if (allocated(results_gs%fcart))   then
   ABI_DEALLOCATE(results_gs%fcart)
 end if
 if (allocated(results_gs%fred))    then
   ABI_DEALLOCATE(results_gs%fred)
 end if
 if (allocated(results_gs%gaps))    then
   ABI_DEALLOCATE(results_gs%gaps)
 end if
 if (allocated(results_gs%gresid))  then
   ABI_DEALLOCATE(results_gs%gresid)
 end if
 if (allocated(results_gs%grewtn))  then
   ABI_DEALLOCATE(results_gs%grewtn)
 end if
 if (allocated(results_gs%grchempottn))  then
   ABI_DEALLOCATE(results_gs%grchempottn)
 end if
 if (allocated(results_gs%grvdw))  then
   ABI_DEALLOCATE(results_gs%grvdw)
 end if
 if (allocated(results_gs%grxc))    then
   ABI_DEALLOCATE(results_gs%grxc)
 end if
 if (allocated(results_gs%synlgr))  then
   ABI_DEALLOCATE(results_gs%synlgr)
 end if

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
!! PARENTS
!!
!! CHILDREN
!!      energies_copy,energies_ncwrite
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
       results_gs(jj,ii)%nsppol=0
       results_gs(jj,ii)%berryopt=0
       if (allocated(results_gs(jj,ii)%fcart))   then
         ABI_DEALLOCATE(results_gs(jj,ii)%fcart)
       end if
       if (allocated(results_gs(jj,ii)%fred))    then
         ABI_DEALLOCATE(results_gs(jj,ii)%fred)
       end if
       if (allocated(results_gs(jj,ii)%gaps))   then
         ABI_DEALLOCATE(results_gs(jj,ii)%gaps)
       end if
       if (allocated(results_gs(jj,ii)%gresid))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%gresid)
       end if
       if (allocated(results_gs(jj,ii)%grewtn))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%grewtn)
       end if
       if (allocated(results_gs(jj,ii)%grchempottn))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%grchempottn)
       end if
       if (allocated(results_gs(jj,ii)%grvdw))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%grvdw)
       end if
       if (allocated(results_gs(jj,ii)%grxc))    then
         ABI_DEALLOCATE(results_gs(jj,ii)%grxc)
       end if
       if (allocated(results_gs(jj,ii)%synlgr))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%synlgr)
       end if
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
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!      energies_copy,energies_ncwrite
!!
!! SOURCE

subroutine copy_results_gs(results_gs_in,results_gs_out)

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(in) :: results_gs_in
 type(results_gs_type),intent(inout) :: results_gs_out !vz_i

!Local variables-------------------------------
!scalars
 integer :: natom_in,natom_out,ngrvdw_in,nsppol_in,nsppol_out

!************************************************************************

 !@results_gs_type

 natom_in =results_gs_in%natom
 natom_out=results_gs_out%natom
 ngrvdw_in =results_gs_in%ngrvdw
 nsppol_in =results_gs_in%nsppol
 nsppol_out =results_gs_out%nsppol

 if (natom_in>natom_out) then
   if (allocated(results_gs_out%fcart))   then
     ABI_DEALLOCATE(results_gs_out%fcart)
   end if
   if (allocated(results_gs_out%fred))    then
     ABI_DEALLOCATE(results_gs_out%fred)
   end if
   if (allocated(results_gs_out%gresid))  then
     ABI_DEALLOCATE(results_gs_out%gresid)
   end if
   if (allocated(results_gs_out%grewtn))  then
     ABI_DEALLOCATE(results_gs_out%grewtn)
   end if
  if (allocated(results_gs_out%grchempottn))  then
     ABI_DEALLOCATE(results_gs_out%grchempottn)
   end if
   if (allocated(results_gs_out%grvdw))  then
     ABI_DEALLOCATE(results_gs_out%grvdw)
   end if
   if (allocated(results_gs_out%grxc))    then
     ABI_DEALLOCATE(results_gs_out%grxc)
   end if
   if (allocated(results_gs_out%synlgr))  then
     ABI_DEALLOCATE(results_gs_out%synlgr)
   end if

   if (allocated(results_gs_in%fcart))   then
     ABI_ALLOCATE(results_gs_out%fcart,(3,natom_in))
   end if
   if (allocated(results_gs_in%fred))    then
     ABI_ALLOCATE(results_gs_out%fred,(3,natom_in))
   end if
   if (allocated(results_gs_in%gresid))  then
     ABI_ALLOCATE(results_gs_out%gresid,(3,natom_in))
   end if
   if (allocated(results_gs_in%grewtn))  then
     ABI_ALLOCATE(results_gs_out%grewtn,(3,natom_in))
   end if
   if (allocated(results_gs_in%grchempottn))  then
     ABI_ALLOCATE(results_gs_out%grchempottn,(3,natom_in))
   end if
   if (allocated(results_gs_in%grvdw))  then
     ABI_ALLOCATE(results_gs_out%grvdw,(3,ngrvdw_in))
   end if
   if (allocated(results_gs_in%grxc))    then
     ABI_ALLOCATE(results_gs_out%grxc,(3,natom_in))
   end if
   if (allocated(results_gs_in%synlgr))  then
     ABI_ALLOCATE(results_gs_out%synlgr,(3,natom_in))
   end if
 end if

 if (nsppol_in>nsppol_out) then
   if (allocated(results_gs_out%gaps))   then
     ABI_DEALLOCATE(results_gs_out%gaps)
   end if
   if (allocated(results_gs_in%gaps))    then
     ABI_ALLOCATE(results_gs_out%gaps,(3,nsppol_in))
   end if
 endif


 results_gs_out%natom  =results_gs_in%natom
 results_gs_out%ngrvdw =results_gs_in%ngrvdw
 results_gs_out%nsppol =results_gs_in%nsppol
 results_gs_out%berryopt=results_gs_in%berryopt
 results_gs_out%deltae =results_gs_in%deltae
 results_gs_out%diffor =results_gs_in%diffor
 results_gs_out%entropy=results_gs_in%entropy
 results_gs_out%etotal =results_gs_in%etotal
 results_gs_out%fermie =results_gs_in%fermie
 results_gs_out%residm =results_gs_in%residm
 results_gs_out%res2   =results_gs_in%res2
 results_gs_out%vxcavg =results_gs_in%vxcavg

 call energies_copy(results_gs_in%energies,results_gs_out%energies)

 results_gs_out%pel(:)=results_gs_in%pel(:)
 results_gs_out%pion(:)=results_gs_in%pion(:)
 results_gs_out%strten(:)=results_gs_in%strten(:)

 if (allocated(results_gs_in%fcart))  results_gs_out%fcart(:,1:natom_in) =results_gs_in%fcart(:,1:natom_in)
 if (allocated(results_gs_in%fred))   results_gs_out%fred(:,1:natom_in)  =results_gs_in%fred(:,1:natom_in)
 if (allocated(results_gs_in%gaps))   results_gs_out%gaps(:,1:nsppol_in) =results_gs_in%gaps(:,1:nsppol_in)
 if (allocated(results_gs_in%gresid)) results_gs_out%gresid(:,1:natom_in)=results_gs_in%gresid(:,1:natom_in)
 if (allocated(results_gs_in%grewtn)) results_gs_out%grewtn(:,1:natom_in)=results_gs_in%grewtn(:,1:natom_in)
 if (allocated(results_gs_in%grchempottn))&
&  results_gs_out%grchempottn(:,1:natom_in)=results_gs_in%grchempottn(:,1:natom_in)
 if (allocated(results_gs_in%grxc))   results_gs_out%grxc(:,1:natom_in)  =results_gs_in%grxc(:,1:natom_in)
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
!! PARENTS
!!      m_results_gs
!!
!! CHILDREN
!!      energies_ncwrite,results_gs_ncwrite
!!
!! SOURCE

integer function results_gs_ncwrite(res, ncid, ecut, pawecutdg) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 real(dp),intent(in) :: ecut,pawecutdg
 type(results_gs_type),intent(in) :: res

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF

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
   "ecut", "pawecutdg", "deltae", "diffor", "entropy", "etotal", "fermie", "residm", "res2"])
 NCF_CHECK(ncerr)

 ! arrays
 !
 ! Note: unlike fred, this array has been corrected by enforcing
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

! Write data.
! Write variables
 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
&  'ecut', 'pawecutdg', 'deltae', 'diffor', 'entropy', 'etotal', 'fermie', 'residm', 'res2'],&
&  [ecut, pawecutdg, res%deltae, res%diffor, res%entropy, res%etotal, res%fermie, res%residm, res%res2],&
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

#else
 MSG_ERROR("netcdf support is not activated.")
#endif

contains
 integer function vid(vname)

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function results_gs_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/results_gs_yaml_write
!!
!! NAME
!! results_gs_yaml_write
!!
!! FUNCTION
!! Write results_gs in yaml format to unit iout
!!
!! INPUTS
!!  results <type(results_gs_type)>=miscellaneous information about the system after ground state computation
!!  iout= unit of output file
!!  [comment] optional comment for the final document
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine results_gs_yaml_write(results, iout, dtset, cryst, comment)

 class(results_gs_type),intent(in) :: results
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: iout
 character(len=*),intent(in),optional :: comment

!Local variables-------------------------------
 integer,parameter :: width=10
 integer :: ii
 type(yamldoc_t) :: ydoc
 type(pair_list) :: dict
 real(dp) :: strten(3,3), abc(3)

!************************************************************************

 if (present(comment)) then
   ydoc = yamldoc_open('ResultsGS', comment, width=width)
 else
   ydoc = yamldoc_open('ResultsGS', '', width=width)
 end if
 ydoc%use_yaml = dtset%use_yaml

 call ydoc%add_int('natom', results%natom)
 call ydoc%add_int('nsppol', results%nsppol)
 call ydoc%add_int('nspinor', dtset%nspinor)
 call ydoc%add_int('nspden', dtset%nspden)
 call ydoc%add_real("nelect", dtset%nelect)
 call ydoc%add_real("charge", dtset%charge)

 call dict%set('ecut', r=dtset%ecut)
 call dict%set('pawecutdg', r=dtset%pawecutdg)
 call ydoc%add_dict('cutoff_energies', dict)
 call dict%free()

 call dict%set('deltae', r=results%deltae)
 call dict%set('res2', r=results%res2)
 call dict%set('residm', r=results%residm)
 call dict%set('diffor', r=results%diffor)
 call ydoc%add_dict('convergence', dict, multiline_trig=2)
 call dict%free()

 abc(:) = [(sqrt(sum(cryst%rprimd(:, ii) ** 2)), ii=1,3)]
 call ydoc%add_real1d('abc', cryst%angdeg)
 call ydoc%add_real1d('alpha_beta_gamma_angles', cryst%angdeg)
 call ydoc%add_real('etotal', results%etotal)
 call ydoc%add_real('entropy', results%entropy)
 call ydoc%add_real('fermie', results%fermie)

 strten(1,1) = results%strten(1)
 strten(2,2) = results%strten(2)
 strten(3,3) = results%strten(3)
 strten(2,3) = results%strten(4)
 strten(3,2) = results%strten(4)
 strten(1,3) = results%strten(5)
 strten(3,1) = results%strten(5)
 strten(1,2) = results%strten(6)
 strten(2,1) = results%strten(6)

 call ydoc%add_real2d('stress_tensor', strten, tag='CartTensor')
 ! Add results in GPa as well
 !strten = strten * HaBohr3_GPa
 !call ydoc%add_real2d('stress_tensor_GPa', strten, tag='CartTensor')
 !call ydoc%add_real('pressure_GPa', get_trace(strten) / three)

 call ydoc%add_real2d('cartesian_forces', results%fcart, tag='CartForces')
 call ydoc%write_and_free(iout)

end subroutine results_gs_yaml_write
!!***

!----------------------------------------------------------------------

END MODULE m_results_gs
!!***
