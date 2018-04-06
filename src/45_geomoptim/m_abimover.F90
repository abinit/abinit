!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abimover
!! NAME
!! m_abimover
!!
!! FUNCTION
!! This module contains definition the type abimover
!! and its related routines
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (DCA, XG, GMR, SE,  Mver, JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abimover

 use defs_basis
 use m_profiling_abi
 use m_atomdata
 use m_errors

 use m_geometry,  only : acrossb

 implicit none

 private

 ! Helper function
 public :: bond_length

 integer,public, parameter :: mover_BEFORE=0
 integer,public, parameter :: mover_AFTER=1
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/delocint
!! NAME
!! delocint
!!
!! FUNCTION
!! Datatype with the important variables in pred_delocint
!!
!! NOTES
!!   deloc <type(delocint)>=Important variables for pred_delocint
!!   |
!!   | nang     = Number of angles
!!   | nbond    = Number of bonds
!!   | ncart    = Number of cartesian directions (used for constraints)
!!   | ndihed   = Number of dihedrals
!!   | nrshift  = Dimension of rshift
!!   | ninternal= Number of internal coordinates
!!   |            ninternal=nbond+nang+ndihed+ncart
!!   | angs(2,3,nang)  = Indexes to characterize angles
!!   | bonds(2,2,nbond)= For a bond between iatom and jatom
!!   |                   bonds(1,1,nbond) = iatom
!!   |                   bonds(2,1,nbond) = icenter
!!   |                   bonds(1,2,nbond) = jatom
!!   |                   bonds(2,2,nbond) = irshift
!!   | carts(2,ncart)  = Index of total primitive internal, and atom (carts(2,:))
!!   | dihedrals(2,4,ndihed)= Indexes to characterize dihedrals
!!   | rshift(3,nrshift)= Shift in xred that must be done to find
!!   |                    all neighbors of a given atom within a
!!   |                    given number of neighboring shells
!!
!! SOURCE

type,public :: delocint

! scalars
 integer :: nang
 integer :: nbond
 integer :: ncart
 integer :: ndihed
 integer :: nrshift
 integer :: ninternal

! arrays
 integer,allocatable :: angs(:,:,:)
 integer,allocatable :: bonds(:,:,:)
 integer,allocatable :: carts(:,:)
 integer,allocatable :: dihedrals(:,:,:)
 real(dp),allocatable :: rshift(:,:)

end type delocint

 public :: make_prim_internals
 public :: delocint_fin
 public :: calc_prim_int
 !public :: xcart2deloc
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/abimover
!! NAME
!! abimover
!!
!! FUNCTION
!! This datatype has the purpose of store all the data taked
!! usually from dtset needed for the different predictors
!! to update positions, acell, etc.
!!
!! NOTES
!!  At present 32 variables are present in abimover
!!  if a new variable is added in abimover it should
!!  be added also for nullify in abimover_nullify
!!
!! STATS
!!  integer         => 13
!!  real            =>  7
!!  integer array   =>  4
!!  real array      =>  6
!!  character array =>  1
!!  structures      =>  1
!!  TOTAL              34
!!
!! SOURCE

type, public :: abimover

! scalars
! Delay of Permutation (Used by pred_langevin only)
integer,pointer  :: delayperm
! DIIS memory (Used by pred_diisrelax only)
integer,pointer  :: diismemory
! Geometry Optimization Precondition option
integer,pointer  :: goprecon
! include a JELLium SLAB in the cell
integer,pointer  :: jellslab
! Number of ATOMs
integer,pointer  :: natom
! Number of CONstraint EQuations
integer,pointer  :: nconeq
! Use by pred_isothermal only
integer,pointer  :: nnos
! Number of SYMmetry operations
integer,pointer  :: nsym
! Number of Types of atoms
integer,pointer  :: ntypat
! OPTimize the CELL shape and dimensions
integer,pointer  :: optcell
! RESTART Xcart and Fred
integer,pointer  :: restartxf
! Sign of Permutation (Used by pred_langevin only)
integer,pointer  :: signperm
! Ion movement
integer,pointer  :: ionmov

! Use by pred_isothermal only
real(dp),pointer :: bmass
! Delta Time for IONs
real(dp),pointer :: dtion
! Used by pred_langevin only
real(dp),pointer :: friction
! Used by pred_langevin only
real(dp),pointer :: mdwall
! Used by pred_nose only
real(dp),pointer :: noseinert
! STRess PRECONditioner
real(dp),pointer :: strprecon
! VIScosity
real(dp),pointer :: vis

! arrays
! Indices of AToms that are FIXed
integer,pointer  :: iatfix(:,:)         ! iatfix(3,natom)
! SYMmetry in REaL space
integer,pointer  :: symrel(:,:,:)       ! symrel(3,3,nsym)
! TYPe of ATom
integer,pointer  :: typat(:)            ! typat(natom)
! PRTint ATom LIST
integer,pointer  :: prtatlist(:)        ! prtatlist(natom)

! Mass of each atom (NOT IN DTSET)
real(dp),pointer :: amass(:)            ! amass(natom)
! Geometry Optimization Preconditioner PaRaMeters
real(dp),pointer :: goprecprm(:)
! Molecular Dynamics Initial and Final Temperature
real(dp),pointer :: mdtemp(:)           ! mdtemp(2) (initial,final)
! STRess TARGET
real(dp),pointer :: strtarget(:)        ! strtarget(6)
! Use by pred_isothermal only
real(dp),pointer :: qmass(:)
! Z number of each NUCLeus
real(dp),pointer :: znucl(:)            ! znucl(npsp)

! Filename for Hessian matrix
character(len=fnlen), pointer :: fnameabi_hes
! Filename for _HIST file
character(len=fnlen), pointer :: filnam_ds(:)   ! dtfil%filnam_ds(5)

! structure for delocalized internal coordinates
type(delocint) :: deloc

end type abimover

 public :: abimover_ini
 public :: abimover_fin
 public :: abimover_nullify
!!***

!----------------------------------------------------------------------

 type,public :: abimover_specs

   !scalars
   integer           :: ncycle
   integer           :: nhist ! Number of step of history needed in the algorithm
   character(len=8)  :: crit4xml
   character(len=10) :: type4xml
   character(len=60) :: method
   logical :: isFconv ! If the convergence is needed
   logical :: isARused
   logical :: isVused

 end type abimover_specs
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/mttk_type
!! NAME
!! mttk_type
!!
!! FUNCTION
!! For Martyna et al. (TTK) reversible MD integration scheme and related data
!!
!! SOURCE

 type, public :: mttk_type

   real(dp) :: glogv
    !Logarithm of the volume

   real(dp) :: vlogv
    !Derivative of logv

  real(dp) :: gboxg(3,3)
   !Imbalance in pressure (see paper)

  real(dp) :: vboxg(3,3)
   !Velocity of log rprimd (see paper)

  real(dp), allocatable :: glogs(:)
   ! glogs(nnos)
   ! Imbalance of kinetic energy

  real(dp), allocatable :: vlogs(:)
   ! vlogs(nnos)
   ! Velocities of thermostat variables

  real(dp), allocatable :: xlogs(:)
   ! xlogs(nnos)
   ! Positions of thermostat variables

 end type mttk_type

 public :: mttk_ini   ! initialize the object
 public :: mttk_fin   ! Release memory
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/abiforstr
!! NAME
!! abiforstr
!!
!! FUNCTION
!! Store forces, stress and energy, cartesian and reduced forces
!! one scalar for energy and 6 element array for stress
!!
!! NOTES
!!
!! SOURCE

type, public :: abiforstr

  ! scalars
  real(dp) :: etotal
   ! Total energy

  ! arrays
  real(dp),allocatable :: fcart(:,:)
   ! Cartesian forces
  real(dp),allocatable :: fred(:,:)
   ! Reduced forces
  real(dp) :: strten(6)
    ! Stress tensor (Symmetrical 3x3 matrix)

end type abiforstr

public :: abiforstr_ini  ! Initialize the object
public :: abiforstr_fin  ! Free memory
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/ab_xfh_type
!! NAME
!! ab_xfh_type
!!
!! FUNCTION
!! Datatype with the old structure for storing history
!! used in gstate and brdmin,delocint, and others
!!
!! NOTES
!! This is a transitional structure, to bridge between
!! the old code and the new one base on abihist
!!
!! SOURCE

type, public :: ab_xfh_type

 integer :: nxfh,nxfhr,mxfh
   ! mxfh = last dimension of the xfhist array
   ! nxfh = actual number of (x,f) history pairs, see xfhist array

 real(dp),allocatable :: xfhist(:,:,:,:)
   ! xfhist(3,natom+4,2,mxfh) = (x,f) history array, also including rprim and stress

end type ab_xfh_type
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/go_bonds
!! NAME
!! go_bonds
!!
!! FUNCTION
!! Datatype all the information relevant to create
!! bonds between atoms inside and outside the cell
!!
!! SOURCE

type, public ::  go_bonds

!scalar
real(dp) :: tolerance ! To decide if consider bond the atom or not
                      ! 1.0 means that only consider values lower
                      ! than the sum of covalent radius

integer  :: nbonds ! Total number of bonds for the system

!arrays

integer,allocatable :: nbondi(:)    ! Number of bonds for atom i
integer,allocatable :: indexi(:,:)  ! Indices of bonds for atom i
                                ! Positive: Vector from i to j
                                ! Negative: Vector from j to i

real(dp),allocatable :: bond_length(:) ! Bond lengths
real(dp),allocatable :: bond_vect(:,:) ! Unitary vectors for bonds

end type go_bonds

public :: make_bonds_new
public :: bonds_free
public :: print_bonds
!!***

!----------------------------------------------------------------------

!!****t* m_abimover/go_angles
!! NAME
!! go_angles
!!
!! FUNCTION
!! Datatype all the information relevant to create
!! angles between atoms inside and outside the cell
!!
!! NOTES
!!  This type is not used
!!
!! SOURCE

type, public :: go_angles

 !scalar
 integer  :: nangles ! Total number of bonds for the system

 !arrays
 integer,allocatable  :: angle_vertex(:)  ! Indices of the vertex atom
 real(dp),allocatable :: angle_value(:)   ! Value of angle in radians
 real(dp),allocatable :: angle_bonds(:,:) ! Indices of the bonds
 real(dp),allocatable :: angle_vect(:,:)  ! Unitary vector perpendicular to the plane

end type go_angles

!public :: make_angles_new ! This routine is broken and should be tested before use.
!!***

!----------------------------------------------------------------------

contains  !=============================================================
!!***


!!****f* m_abimover/abimover_ini
!! NAME
!! abimover_ini
!!
!! FUNCTION
!! Initializes the abimover structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine abimover_ini(ab_mover,specs,&
& delayperm,diismemory,goprecon,jellslab,&
& natom,nconeq,nnos,nsym,ntypat,&
& optcell,restartxf,signperm,ionmov,&
& bmass,dtion,friction,mdwall,noseinert,strprecon,vis,&
& iatfix,symrel,typat,prtatlist,amass,goprecprm,mdtemp,strtarget,&
& qmass,znucl,fnameabi_hes,filnam_ds)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimover_ini'
!End of the abilint section

type(abimover),intent(inout) :: ab_mover
type(abimover_specs),intent(inout) :: specs
! Delay of Permutation (Used by pred_langevin only)
integer,target, intent(in)  :: delayperm
! DIIS memory (Used by pred_diisrelax only)
integer,target, intent(in)  :: diismemory
! Geometry Optimization Precondition option
integer,target, intent(in)  :: goprecon
! include a JELLium SLAB in the cell
integer,target, intent(in)  :: jellslab
! Number of ATOMs
integer,target, intent(in)  :: natom
! Number of CONstraint EQuations
integer,target, intent(in)  :: nconeq
! Use by pred_isothermal only
integer,target, intent(in)  :: nnos
! Number of SYMmetry operations
integer,target, intent(in)  :: nsym
! Number of Types of atoms
integer,target, intent(in)  :: ntypat
! OPTimize the CELL shape and dimensions
integer,target, intent(in)  :: optcell
! RESTART Xcart and Fred
integer,target, intent(in)  :: restartxf
! Sign of Permutation (Used by pred_langevin only)
integer,target, intent(in)  :: signperm
! Ion movement
integer,target, intent(in)  :: ionmov

! Use by pred_isothermal only
real(dp),target, intent(in) :: bmass
! Delta Time for IONs
real(dp),target, intent(in) :: dtion
! Used by pred_langevin only
real(dp),target, intent(in) :: friction
! Used by pred_langevin only
real(dp),target, intent(in) :: mdwall
! Used by pred_nose only
real(dp),target, intent(in) :: noseinert
! STRess PRECONditioner
real(dp),target, intent(in) :: strprecon
! VIScosity
real(dp),target, intent(in) :: vis

! arrays
! Indices of AToms that are FIXed
integer,target, intent(in)  :: iatfix(:,:)         ! iatfix(3,natom)
! SYMmetry in REaL space
integer,target, intent(in)  :: symrel(:,:,:)       ! symrel(3,3,nsym)
! TYPe of ATom
integer,target, intent(in)  :: typat(:)            ! typat(natom)
! PRTint ATom LIST
integer,target, intent(in)  :: prtatlist(:)        ! prtatlist(natom)

! Mass of each atom (NOT IN DTSET)
real(dp),target, intent(in) :: amass(:)            ! amass(natom)
! Geometry Optimization Preconditioner PaRaMeters
real(dp),target, intent(in) :: goprecprm(:)
! Molecular Dynamics Initial and Final Temperature
real(dp),target, intent(in) :: mdtemp(:)           ! mdtemp(2) (initial,final)
! STRess TARGET
real(dp),target, intent(in) :: strtarget(:)        ! strtarget(6)
! Use by pred_isothermal only
real(dp),target, intent(in) :: qmass(:)
! Z number of each NUCLeus
real(dp),target, intent(in) :: znucl(:)            ! znucl(npsp)

! Filename for Hessian matrix
character(len=fnlen), target, intent(in) :: fnameabi_hes
! Filename for _HIST file
character(len=fnlen), target, intent(in) :: filnam_ds(:)   ! dtfil%filnam_ds(5)


!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays

! ***************************************************************


!write(std_out,*) 'mover 01'
!###########################################################
!### 01. Initialization of indexes and allocations of arrays


!Copy the information from the Dataset (dtset)
!to the ab_mover structure
!Delay of Permutation (Used by pred_langevin only)
 ab_mover%delayperm=>delayperm
!DIIS memory (Use by pred_diisrelax only)
 ab_mover%diismemory=>diismemory
!Geometry Optimization Precondition option
 ab_mover%goprecon=>goprecon
!include a JELLium SLAB in the cell
 ab_mover%jellslab=>jellslab
!Number of ATOMs
 ab_mover%natom=>natom
!Number of CONstraint EQuations
 ab_mover%nconeq=>nconeq
!Use by pred_isothermal only
 ab_mover%nnos=>nnos
!Number of SYMmetry operations
 ab_mover%nsym=>nsym
!Number of Type of atoms
 ab_mover%ntypat=>ntypat
!OPTimize the CELL shape and dimensions
 ab_mover%optcell=>optcell
!RESTART Xcart and Fred
 ab_mover%restartxf=>restartxf
!Sign of Permutation (Used by pred_langevin only)
 ab_mover%signperm=>signperm
!Ion movements
 ab_mover%ionmov=>ionmov

!Use by pred_isothermal only
 ab_mover%bmass=>bmass
!Delta Time for IONs
 ab_mover%dtion=>dtion
!Used by pred_langevin only
 ab_mover%friction=>friction
!Used by pred_langevin only
 ab_mover%mdwall=>mdwall
!NOT DOCUMENTED
 ab_mover%noseinert=>noseinert
!STRess PRECONditioner
 ab_mover%strprecon=>strprecon
!VIScosity
 ab_mover%vis=>vis

!Indices of AToms that are FIXed
 ab_mover%iatfix=>iatfix(:,1:natom)
!SYMmetry in REaL space
 ab_mover%symrel=>symrel
!TYPe of ATom
 ab_mover%typat => typat(1:natom)
!PRTint QTom LIST
 ab_mover%prtatlist=>prtatlist(1:natom)

!Mass of each atom (NOT IN DTSET)
 ab_mover%amass=>amass
!Geometry Optimization Preconditioner PaRaMeters
 ab_mover%goprecprm=>goprecprm
!Molecular Dynamics Temperatures
 ab_mover%mdtemp=>mdtemp
!STRess TARGET
 ab_mover%strtarget=>strtarget
!Use by pred_isothermal only
 ab_mover%qmass=>qmass
!Z number of each NUCLeus
 ab_mover%znucl=>znucl

!Filename for Hessian matrix (NOT IN DTSET)
 ab_mover%fnameabi_hes=>fnameabi_hes
!Filename for _HIST file
 ab_mover%filnam_ds=>filnam_ds


!!DEBUG
!call abimover_print(ab_mover,ab_out)
!!DEBUG

!write(std_out,*) 'mover 02'
!###########################################################
!### 02. Particularities of each predictor

!Default values first
!--------------------

!acell and rprimd are never change except if optcell/=0
 if (ab_mover%optcell/=0)then
   specs%isARused=.TRUE.
 else
   specs%isARused=.FALSE.
 end if

!Velocities are never change except for ionmov=1,6,7,8
 specs%isVused=.FALSE.

!In general convergence is needed
 specs%isFconv=.TRUE.

!specs%ncycle is 1 by default except for ionmov=1,9,14
 specs%ncycle=1

!specs%nhist is -1 by default store all the history except for ionmov=12,13
 specs%nhist=-1

!This is the initialization for ionmov==1
!-----------------------------------------
 select case (ab_mover%ionmov)
 case (1)
   specs%ncycle=4 ! Number of internal cycles for first itime
   specs%isFconv=.FALSE.     ! Convergence is not used for MD
   specs%isVused=.TRUE. ! Velocities are used
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='moldyn'
   specs%crit4xml='none'
!  Name of specs%method
   if (abs(ab_mover%vis)<=1.d-8) then
     specs%method = 'Molecular dynamics without viscosity (vis=0)'
   else
     write(specs%method,'(a,1p,e12.5,a)')&
     'Molecular dynamics with viscosity (vis=',ab_mover%vis,')'
   end if
!  Number of history
   specs%nhist = 6
!  This is the initialization for ionmov==2,3
!  -------------------------------------------
 case (2,3)
!  Values use in XML Output
   specs%type4xml='bfgs'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   if (ab_mover%ionmov==2) then
     specs%method = 'Broyden-Fletcher-Goldfard-Shanno method (forces)'
   else
     specs%method = 'Broyden-Fletcher-Goldfard-Shanno method (forces,Tot energy)'
   end if
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==4,5
!  -------------------------------------------
 case (4,5)
!  Values use in XML Output
   specs%type4xml='simple'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   if (ab_mover%ionmov==4) then
     specs%method = 'Conjugate gradient of potential and ionic degrees of freedom'
   else
     specs%method = 'Simple relaxation of ionic positions'
   end if
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==6
!  ------------------------------------------
 case (6)
   specs%isFconv=.FALSE.     ! Convergence is not used for MD
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
   specs%isVused=.TRUE. ! Velocities are used
!  Values use in XML Output
   specs%type4xml='verlet'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Verlet algorithm for molecular dynamics'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==7
!  ------------------------------------------
 case (7)
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
   specs%isVused=.TRUE. ! Velocities are used
!  Values use in XML Output
   specs%type4xml='verlet'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Verlet algorithm blocking every atom where dot(vel,force)<0'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==8
!  ------------------------------------------
 case (8)
   specs%isVused=.TRUE.
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='nose'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Verlet algorithm with a nose-hoover thermostat'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==9
!  ------------------------------------------
 case (9)
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
   specs%isVused=.TRUE.  ! Velocities are used
   specs%ncycle=3
!  Values use in XML Output
   specs%type4xml='langevin'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Langevin molecular dynamics'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==10
!  -------------------------------------------
 case (10)
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='delocint'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'BFGS with delocalized internal coordinates'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==11
!  -------------------------------------------
 case (11)
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='cg'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Conjugate gradient algorithm'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==12
!  -------------------------------------------
 case (12)
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
   specs%isVused=.TRUE.  ! Velocities are used
!  Values use in XML Output
   specs%isFconv=.FALSE.      ! Convergence is not used for MD
   specs%type4xml='isokin'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Isokinetic ensemble molecular dynamics'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==13
!  -------------------------------------------
 case (13)
!  optcell is allow
   specs%isARused=.TRUE. ! RPRIMD and ACELL may change
   specs%isVused=.TRUE.  ! Velocities are used
   specs%isFconv=.FALSE.      ! Convergence is not used for MD
!  Values use in XML Output
   specs%type4xml='isother'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Isothermal/isenthalpic ensemble molecular dynamics'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==14
!  -------------------------------------------
 case (14)
   specs%ncycle=16
   specs%isFconv=.FALSE.     ! Convergence is not used for MD
   specs%isVused=.TRUE. ! Velocities are used
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='srkna14'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Symplectic algorithm Runge-Kutta-Nystrom SRKNa14'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==20
!  -------------------------------------------
 case (20)
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='diisrelax'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Ionic positions relaxation using DIIS'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==21
!  -------------------------------------------
 case (21)
   specs%isARused=.TRUE.
!  Values use in XML Output
   specs%type4xml='steepdesc'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Steepest descend algorithm'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==22
!  -------------------------------------------
 case (22)
!  Values use in XML Output
   specs%type4xml='lbfgs'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Limited-memory Broyden-Fletcher-Goldfard-Shanno method'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==23
!  -------------------------------------------
 case (23)
   specs%ncycle=2
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
   specs%isVused=.TRUE.  ! Velocities are used
!  Values use in XML Output
   specs%isFconv=.FALSE.      ! Convergence is not used for MD
   specs%type4xml='isokin'
   specs%crit4xml='tolmxf'
!  Name of specs%method
   specs%method = 'Using LOTF Molecular dynamics'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==24
!  -------------------------------------------
 case (24)
   specs%ncycle=1
!  TEMPORARLY optcell is not allow
   specs%isARused=.FALSE.
   specs%isVused=.TRUE.  ! Velocities are used
!  Values use in XML Output
   specs%isFconv=.FALSE.      ! Convergence is not used for MD
   specs%type4xml='velver'
   specs%crit4xml='none'
!  Name of specs%method
   specs%method = 'Symplectic velocity verlet Molecular dynamics'
!  Number of history
   specs%nhist = 3
!  This is the initialization for ionmov==25
!  -------------------------------------------
 case (25)                ! Hybrid Monte Carlo algorithm (fixed lattice vectors)
   specs%ncycle = 12      ! Number of internal cycles (10+2)
   specs%isFconv=.FALSE.  ! Convergence is not used for Monte Carlo
   specs%isVused=.TRUE.   ! Velocities are used for update of atomic positions
!  optcell is not allowed
   specs%isARused=.FALSE.
!  Values use in XML Output
   specs%type4xml='hmc'
   specs%crit4xml='none'
!  Name of specs%method
   specs%method = 'Hybrid Monte Carlo'

case default
   write(msg,"(a,i0)")"Wrong value for ionmov: ",ab_mover%ionmov
 end select
end subroutine abimover_ini
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/delocint_fin
!!
!! NAME
!! delocint_fin
!!
!! FUNCTION
!! destructor function for delocint object
!!
!! INPUT
!! deloc= container object for delocalized internal coordinates
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abimover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine delocint_fin(deloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'delocint_fin'
!End of the abilint section

 type(delocint), intent(inout) :: deloc

 if(allocated(deloc%angs))  then
   ABI_DEALLOCATE(deloc%angs)
 end if
 if(allocated(deloc%bonds))  then
   ABI_DEALLOCATE(deloc%bonds)
 end if
 if(allocated(deloc%carts))  then
   ABI_DEALLOCATE(deloc%carts)
 end if
 if(allocated(deloc%dihedrals))  then
   ABI_DEALLOCATE(deloc%dihedrals)
 end if
 if(allocated(deloc%rshift))  then
   ABI_DEALLOCATE(deloc%rshift)
 end if

end subroutine delocint_fin
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/abimover_nullify
!! NAME
!! abimover_nullify
!!
!! FUNCTION
!! Nullify all the pointers in a ab_mover
!!
!! SIDE EFFECTS
!!  ab_mover <type(abimover)> = The ab_mover to nullify
!!
!! PARENTS
!!      m_abimover,mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! NOTES
!!  At present 32 variables are present in abimover
!!  if a new variable is added in abimover it should
!!  be added also for nullify here
!!
!! STATS
!!  integer         => 13
!!  real            =>  7
!!  integer array   =>  4
!!  real array      =>  6
!!  character array =>  1
!!  structures      =>  1
!!  TOTAL              32
!!
!! SOURCE

subroutine abimover_nullify(ab_mover)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimover_nullify'
!End of the abilint section

 type(abimover),intent(inout) :: ab_mover

! ***************************************************************

! Delay of Permutation (Used by pred_langevin only)
 nullify(ab_mover%delayperm)
! Number of Geometry histories (Use by pred_diisrelax only)
 nullify(ab_mover%diismemory)
! Geometry Optimization precondition option
 nullify(ab_mover%goprecon)
! include a JELLium SLAB in the cell
 nullify(ab_mover%jellslab)
! Number of ATOMs
 nullify(ab_mover%natom)
! Number of CONstraint EQuations
 nullify(ab_mover%nconeq)
! Use by pred_isothermal only
 nullify(ab_mover%nnos)
! Number of SYMmetry operations
 nullify(ab_mover%nsym)
! Number of Types of atoms
 nullify(ab_mover%ntypat)
! OPTimize the CELL shape and dimensions
 nullify(ab_mover%optcell)
! RESTART Xcart and Fred
 nullify(ab_mover%restartxf)
! Sign of Permutation (Used by pred_langevin only)
 nullify(ab_mover%signperm)
! Ion movement
 nullify(ab_mover%ionmov)

! Use by pred_isothermal only
 nullify(ab_mover%bmass)
! Delta Time for IONs
 nullify(ab_mover%dtion)
! Used by pred_langevin only
 nullify(ab_mover%friction)
! Used by pred_langevin only
 nullify(ab_mover%mdwall)
! NOT DOCUMENTED
 nullify(ab_mover%noseinert)
! STRess PRECONditioner
 nullify(ab_mover%strprecon)
! VIScosity
 nullify(ab_mover%vis)

!arrays
! Indices of AToms that are FIXed
 nullify(ab_mover%iatfix)
! SYMmetry in REaL space
 nullify(ab_mover%symrel)
! TYPe of ATom
 nullify(ab_mover%typat)
! TYPe of ATom
 nullify(ab_mover%prtatlist)

! Mass of each atom (NOT IN DTSET)
 nullify(ab_mover%amass)!
! Molecular Dynamics Initial Temperature
 nullify(ab_mover%mdtemp)
! STRess TARGET
 nullify(ab_mover%strtarget)
! Use by pred_isothermal only
 nullify(ab_mover%qmass)
! Z number of each NUCLeus
 nullify(ab_mover%znucl)
! Geometry Optimization Preconditioner PaRaMeters
 nullify(ab_mover%goprecprm)

! Filename for Hessian matrix
 nullify(ab_mover%fnameabi_hes)
! Filename for _HIST file
 nullify(ab_mover%filnam_ds)

end subroutine abimover_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/abimover_fin
!! NAME
!! abimover_fin
!!
!! FUNCTION
!! Deallocate all the pointers in a ab_mover
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ab_mover <type(abimover)> = The ab_mover to destroy
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! NOTES
!!  At present 32 variables are present in abimover
!!  if a new variable is added in abimover it should
!!  be added also for deallocate here
!!
!! STATS
!!  integer         => 13
!!  real            =>  7
!!  integer array   =>  4
!!  real array      =>  6
!!  character array =>  1
!!  structures      =>  1
!!  TOTAL              32
!!
!! SOURCE

subroutine abimover_fin(ab_mover)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimover_fin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(abimover),intent(inout) :: ab_mover

! ***************************************************************

 call delocint_fin(ab_mover%deloc)

! ab_mover is only pointer associated to other data, except possibly amass.
! TODO: check for amass in mover.F90
 call abimover_nullify(ab_mover)

end subroutine abimover_fin
!!***

!----------------------------------------------------------------------

!!****f* defs_mover/abimover_print
!! NAME
!! abimover_print
!!
!! FUNCTION
!! Print all the variables in a ab_mover
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ab_mover <type(abimover)> = The ab_mover to nullify
!!
!! PARENTS
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! NOTES
!!  At present 29 variables are present in abimover
!!  if a new variable is added in abimover it should
!!  be added also for print here
!!
!! SOURCE

subroutine abimover_print(ab_mover,iout)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimover_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: iout
 type(abimover),intent(inout) :: ab_mover

!Local variables-------------------------------
!arrays
character(len=1200) :: message
character(len=110)   :: fmt

! ***************************************************************

 fmt='(a,e12.5,a,a,I5,a,a,I5,a,a,I5,a,a,I5,a,a,I5,a,a,I5,a,a,e12.5,a,a,e12.5,a,a,e12.5,a,a,e12.5,a,a,e12.5,a)'

 write(message,fmt)&
& 'Delta Time for IONs',ab_mover%dtion,ch10, &
& 'include a JELLium SLAB in the cell',ab_mover%jellslab,ch10, &
& 'Number of ATOMs',ab_mover%natom,ch10, &
& 'Number of CONstraint EQuations',ab_mover%nconeq,ch10, &
& 'Number of SYMmetry operations',ab_mover%nsym,ch10, &
& 'OPTimize the CELL shape and dimensions',ab_mover%optcell,ch10, &
& 'RESTART Xcart and Fred',ab_mover%restartxf,ch10, &
& 'Molecular Dynamics Initial Temperature',ab_mover%mdtemp(1),ch10, &
& 'Molecular Dynamics Final Temperature',ab_mover%mdtemp(2),ch10, &
& 'NOT DOCUMENTED',ab_mover%noseinert,ch10, &
& 'STRess PRECONditioner',ab_mover%strprecon,ch10, &
& 'VIScosity',ab_mover%vis,ch10

! ! arrays
! ! Indices of AToms that are FIXed
! integer,  pointer :: iatfix(:,:)
! ! SYMmetry in REaL space
! integer,  pointer :: symrel(:,:,:)
! ! Mass of each atom (NOT IN DTSET)
! real(dp), pointer :: amass(:)
! ! STRess TARGET
! real(dp), pointer :: strtarget(:)
! Filename for Hessian matrix
! character(len=fnlen), pointer :: fnameabi_hes

 write(iout,*) 'CONTENTS of ab_mover'
 write(iout,'(a)') message

end subroutine abimover_print
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/mttk_ini
!!
!! NAME
!! mttk_ini
!!
!! FUNCTION
!! destructor function for mttk object
!!
!! INPUT
!! mttk
!!
!! OUTPUT
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine mttk_ini(mttk_vars,nnos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mttk_ini'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mttk_ini'
!End of the abilint section

 integer,intent(in)  :: nnos
 type(mttk_type), intent(out) :: mttk_vars

 ABI_ALLOCATE(mttk_vars%glogs,(nnos))
 ABI_ALLOCATE(mttk_vars%vlogs,(nnos))
 ABI_ALLOCATE(mttk_vars%xlogs,(nnos))

end subroutine mttk_ini
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/mttk_fin
!!
!! NAME
!! mttk_fin
!!
!! FUNCTION
!! destructor function for mttk object
!!
!! INPUT
!! mttk
!!
!! OUTPUT
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine mttk_fin(mttk_vars)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mttk_fin'
!End of the abilint section

 implicit none

 type(mttk_type), intent(inout) :: mttk_vars

 if(allocated(mttk_vars%glogs))  then
  ABI_DEALLOCATE(mttk_vars%glogs)
 end if
 if(allocated(mttk_vars%vlogs))  then
  ABI_DEALLOCATE(mttk_vars%vlogs)
 end if
 if(allocated(mttk_vars%xlogs))  then
  ABI_DEALLOCATE(mttk_vars%xlogs)
 end if

end subroutine mttk_fin
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/abiforstr_ini
!!
!! NAME
!! abiforstr_ini
!!
!! FUNCTION
!! destructor function for abiforstr object
!!
!! INPUT
!! forstr
!!
!! OUTPUT
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine abiforstr_ini(forstr,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abiforstr_ini'
!End of the abilint section

 integer,intent(in)  :: natom
 type(abiforstr), intent(out) :: forstr

 ABI_ALLOCATE(forstr%fcart,(3,natom))
 ABI_ALLOCATE(forstr%fred,(3,natom))

end subroutine abiforstr_ini
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/abiforstr_fin
!!
!! NAME
!! abiforstr_fin
!!
!! FUNCTION
!! destructor function for abiforstr object
!!
!! INPUT
!! forstr
!!
!! OUTPUT
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine abiforstr_fin(forstr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abiforstr_fin'
!End of the abilint section

 type(abiforstr), intent(inout) :: forstr

 if(allocated(forstr%fcart))  then
    ABI_DEALLOCATE(forstr%fcart)
 end if
 if(allocated(forstr%fred))  then
    ABI_DEALLOCATE(forstr%fred)
 end if

end subroutine abiforstr_fin
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/make_prim_internals
!! NAME
!! make_prim_internals
!!
!! FUNCTION
!!  Determine the bonds, angles and dihedrals for a starting
!!  geometry, based on covalent radii for the atoms.
!!
!! INPUTS
!! natom  = Number of atoms (dtset%natom)
!! icenter= index of the center of the number of shifts
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!!   deloc <type(delocint)>=Important variables for
!!   |                           pred_delocint
!!   |
!!   | nang     = Number of angles
!!   | nbond    = Number of bonds
!!   | ncart    = Number of cartesian directions
!!   |             (used for constraints)
!!   | ndihed   = Number of dihedrals
!!   | nrshift  = Dimension of rshift
!!   | ninternal= Number of internal coordinates
!!   |            ninternal=nbond+nang+ndihed+ncart
!!   |
!!   | angs(2,3,nang)  = Indexes to characterize angles
!!   | bonds(2,2,nbond)= For a bond between iatom and jatom
!!   |                   bonds(1,1,nbond) = iatom
!!   |                   bonds(2,1,nbond) = icenter
!!   |                   bonds(1,2,nbond) = jatom
!!   |                   bonds(2,2,nbond) = irshift
!!   | carts(2,ncart)  = Index of total primitive internal,
!!   |                   and atom (carts(2,:))
!!   | dihedrals(2,4,ndihed)= Indexes to characterize dihedrals
!!   |
!!   | rshift(3,nrshift)= Shift in xred that must be done to find
!!   |                    all neighbors of a given atom within a
!!   |                    given number of neighboring shells
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!   Adds cartesian coordinates if the number of internals with a
!!   given atom is < 4 the chosen coordinate could be optimized
!!   to be less dependent of the internals already incorporated.
!!
!! PARENTS
!!      pred_delocint
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine make_prim_internals(deloc,icenter,natom,ntypat,rprimd,typat,xcart,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_prim_internals'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(delocint),intent(inout) :: deloc
 integer,intent(in) :: icenter,natom,ntypat
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: znucl(:) ! znucl(ntypat) or znucl(npsp) ?

!Local variables ------------------------------
! function
!scalars
 integer :: iang,iatom,ibond,icart,idihed,ii
 real(dp) :: spp
!arrays
 integer :: particip_atom(natom)
 integer,allocatable :: badangles(:)
 real(dp) :: rpt1(3),rpt3(3) !rpt2(3)

!************************************************************************

 particip_atom(:) = 0

 call make_bonds(deloc,natom,ntypat,icenter,rprimd,typat,xcart,znucl)

 do ibond=1,deloc%nbond
   write(std_out,'(a,i4,2(2i5,2x))') 'bond ', ibond, deloc%bonds(:,:,ibond)
   particip_atom(deloc%bonds(1,1,ibond)) = particip_atom(deloc%bonds(1,1,ibond))+1
   particip_atom(deloc%bonds(1,2,ibond)) = particip_atom(deloc%bonds(1,2,ibond))+1
 end do

 call make_angles(deloc,icenter,natom)

 ABI_ALLOCATE(badangles,(deloc%nang))
 badangles(:) = 0
 do iang=1,deloc%nang
   write(std_out,'(a,i4,3(2i5,2x))') 'angle ', iang, deloc%angs(:,:,iang)
   particip_atom(deloc%angs(1,1,iang)) = particip_atom(deloc%angs(1,1,iang))+1
   particip_atom(deloc%angs(1,2,iang)) = particip_atom(deloc%angs(1,2,iang))+1
   particip_atom(deloc%angs(1,3,iang)) = particip_atom(deloc%angs(1,3,iang))+1

!  DEBUG
!  rpt1(:) = xcart(:,deloc%angs(1,1,iang)) &
!  & + deloc%rshift(1,deloc%angs(2,1,iang))*rprimd(:,1) &
!  & + deloc%rshift(2,deloc%angs(2,1,iang))*rprimd(:,2) &
!  & + deloc%rshift(3,deloc%angs(2,1,iang))*rprimd(:,3)
!  rpt2(:) = xcart(:,deloc%angs(1,2,iang)) &
!  & + deloc%rshift(1,deloc%angs(2,2,iang))*rprimd(:,1) &
!  & + deloc%rshift(2,deloc%angs(2,2,iang))*rprimd(:,2) &
!  & + deloc%rshift(3,deloc%angs(2,2,iang))*rprimd(:,3)
!  rpt3(:) = xcart(:,deloc%angs(1,3,iang)) &
!  & + deloc%rshift(1,deloc%angs(2,3,iang))*rprimd(:,1) &
!  & + deloc%rshift(2,deloc%angs(2,3,iang))*rprimd(:,2) &
!  & + deloc%rshift(3,deloc%angs(2,3,iang))*rprimd(:,3)
!  write(std_out,*) rpt1,rpt2,rpt3,bond_length(rpt1,rpt2),bond_length(rpt2,rpt3)
!  ENDDEBUG

!  check if angles are 180 degrees: discard the dihedrals in that case.
   rpt1(:) = xcart(:,deloc%angs(1,1,iang)) &
&   + deloc%rshift(1,deloc%angs(2,1,iang))*rprimd(:,1) &
&   + deloc%rshift(2,deloc%angs(2,1,iang))*rprimd(:,2) &
&   + deloc%rshift(3,deloc%angs(2,1,iang))*rprimd(:,3) &
&   - xcart(:,deloc%angs(1,2,iang)) &
&   - deloc%rshift(1,deloc%angs(2,2,iang))*rprimd(:,1) &
&   - deloc%rshift(2,deloc%angs(2,2,iang))*rprimd(:,2) &
&   - deloc%rshift(3,deloc%angs(2,2,iang))*rprimd(:,3)

   rpt3(:) = xcart(:,deloc%angs(1,3,iang)) &
&   + deloc%rshift(1,deloc%angs(2,3,iang))*rprimd(:,1) &
&   + deloc%rshift(2,deloc%angs(2,3,iang))*rprimd(:,2) &
&   + deloc%rshift(3,deloc%angs(2,3,iang))*rprimd(:,3) &
&   - xcart(:,deloc%angs(1,2,iang)) &
&   - deloc%rshift(1,deloc%angs(2,2,iang))*rprimd(:,1) &
&   - deloc%rshift(2,deloc%angs(2,2,iang))*rprimd(:,2) &
&   - deloc%rshift(3,deloc%angs(2,2,iang))*rprimd(:,3)
   spp = (rpt1(1)*rpt3(1)+rpt1(2)*rpt3(2)+rpt1(3)*rpt3(3))&
&   / sqrt(rpt1(1)*rpt1(1)+rpt1(2)*rpt1(2)+rpt1(3)*rpt1(3)) &
&   / sqrt(rpt3(1)*rpt3(1)+rpt3(2)*rpt3(2)+rpt3(3)*rpt3(3))
   if (abs(abs(spp) - one) < tol6) then
     write(std_out,*) 'make_prim_internals : an angle is too close to 180 degrees:'
     write(std_out,*) '   will discard dihedrals using it '
     badangles(iang) = 1
   end if
 end do

 call make_dihedrals(badangles,deloc,icenter)
 ABI_DEALLOCATE(badangles)

 do idihed=1,deloc%ndihed
   write(std_out,'(a,i4,4(2i5,2x))') 'dihedral ', idihed, deloc%dihedrals(:,:,idihed)
   particip_atom(deloc%dihedrals(1,1,idihed)) = particip_atom(deloc%dihedrals(1,1,idihed))+1
   particip_atom(deloc%dihedrals(1,2,idihed)) = particip_atom(deloc%dihedrals(1,2,idihed))+1
   particip_atom(deloc%dihedrals(1,3,idihed)) = particip_atom(deloc%dihedrals(1,3,idihed))+1
   particip_atom(deloc%dihedrals(1,4,idihed)) = particip_atom(deloc%dihedrals(1,4,idihed))+1

!  do ii=1,4
!  write(std_out,'((3E16.6,2x))') xcart(:,deloc%dihedrals(1,ii,idihed)) + &
!  &  deloc%rshift(1,deloc%dihedrals(2,ii,idihed))*rprimd(:,1)   + &
!  &  deloc%rshift(2,deloc%dihedrals(2,ii,idihed))*rprimd(:,2)   + &
!  &  deloc%rshift(2,deloc%dihedrals(2,ii,idihed))*rprimd(:,3)
!  end do
 end do

 write(std_out,*) 'make_deloc_internals: nbond,nang,ndihed = ', deloc%nbond,deloc%nang,deloc%ndihed

!Check all atoms participate in at least 4 primitives. Otherwise, we should
!probably add cartesian coordinates to the internal ones.
 deloc%ncart = 0
 do iatom=1,natom
   if (particip_atom(iatom) < 4) then
     write(std_out,*) ' make_prim_internals : Warning : atom ', iatom, &
&     ' does not belong to enough primitives to determine its'
     write(std_out,*) ' position uniquely ! instead : ', particip_atom(iatom)
     write(std_out,*) ' Will add cartesian coordinates to set of internals.'
!    write(std_out,*) ' Not done yet.'
!    stop
     deloc%ncart = deloc%ncart + 4-particip_atom(iatom)
   end if
 end do
 if (allocated(deloc%carts)) then
   ABI_FREE(deloc%carts)
 end if
 ABI_ALLOCATE(deloc%carts ,(2,deloc%ncart))
 icart = 0
 do iatom=1,natom
   if (particip_atom(iatom) < 4) then
!    kind of arbitrary : include first few directions for the atom: x, then y then z
     do ii=1,4-particip_atom(iatom)
       icart = icart+1
       deloc%carts(1,icart) = ii
       deloc%carts(2,icart) = iatom
     end do
   end if
 end do

!ninternal=nbond+nang+ndihed
 deloc%ninternal=deloc%nbond+deloc%nang+deloc%ndihed+deloc%ncart

end subroutine make_prim_internals
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/make_angles
!! NAME
!! make_angles
!!
!! FUNCTION
!!  (to be completed)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      m_abimover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine make_angles(deloc,icenter,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_angles'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter,natom
 type(delocint),intent(inout) :: deloc
!arrays

!Local variables-------------------------------
!scalars
 integer :: ia1,ia2,iang,ibond,is1,is2,ishift,ja1,ja2
 integer :: jbond,js1,js2
!arrays
 integer,allocatable :: angs_tmp(:,:,:)

! *************************************************************************

!tentative first allocation: < 6 angles per bond.
 ABI_ALLOCATE(angs_tmp,(2,3,72*natom))

 deloc%nang = 0

 do ibond=1, deloc%nbond
   ia1 = deloc%bonds(1,1,ibond)
   is1 = deloc%bonds(2,1,ibond)
   ia2 = deloc%bonds(1,2,ibond)
   is2 = deloc%bonds(2,2,ibond)
   do jbond=ibond+1,deloc%nbond
     ja1 = deloc%bonds(1,1,jbond)
     ja2 = deloc%bonds(1,2,jbond)
     do ishift=-(icenter-1),+(icenter-1)
       js1 = deloc%bonds(2,1,jbond)+ishift
       js2 = deloc%bonds(2,2,jbond)+ishift

       if      (ia1==ja1 .and. is1==js1) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,2,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,3,deloc%nang) = (/ja2,js2/)

       else if (ia1==ja2 .and. is1==js2) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,2,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,3,deloc%nang) = (/ja1,js1/)

       else if (ia2==ja2 .and. is2==js2) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,2,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,3,deloc%nang) = (/ja1,js1/)

       else if (ia2==ja1 .and. is2==js1) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,2,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,3,deloc%nang) = (/ja2,js2/)

       end if
       if (deloc%nang > 72*natom) then
         MSG_ERROR('too many angles found > 72*natom')
       end if
     end do
   end do ! jbond do
 end do ! ibond

 if (allocated(deloc%angs)) then
   ABI_FREE(deloc%angs)
 end if
 ABI_ALLOCATE(deloc%angs,(2,3,deloc%nang))
 do iang=1,deloc%nang
   deloc%angs(:,:,iang) = angs_tmp(:,:,iang)
 end do
 ABI_DEALLOCATE(angs_tmp)

end subroutine make_angles
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/make_dihedrals
!! NAME
!! make_dihedrals
!!
!! FUNCTION
!!  (to be completed)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      m_abimover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine make_dihedrals(badangles,deloc,icenter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_dihedrals'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter
 type(delocint),intent(inout) :: deloc
!arrays
 integer,intent(in) :: badangles(deloc%nang)

!Local variables-------------------------------
!scalars
 integer :: chkdihed,ia1,ia2,ia3,iang,idihed,is1,is2
 integer :: is3,ishift,ja1,ja2,ja3,jang,js1,js2,js3,maxshift
 integer :: minshift
!arrays
 integer,allocatable :: diheds_tmp(:,:,:)

! *************************************************************************

!tentative first allocation: < 6 dihedrals per angle.
 ABI_ALLOCATE(diheds_tmp,(2,4,6*deloc%nang))

 deloc%ndihed = 0
 diheds_tmp(:,:,:) = 0

 do iang=1,deloc%nang
   if (badangles(iang) == 1) cycle
   ia1 = deloc%angs(1,1,iang)
   is1 = deloc%angs(2,1,iang)
   ia2 = deloc%angs(1,2,iang)
   is2 = deloc%angs(2,2,iang)
   ia3 = deloc%angs(1,3,iang)
   is3 = deloc%angs(2,3,iang)

   do jang=iang+1,deloc%nang
     if (badangles(jang) == 1) cycle
     ja1 = deloc%angs(1,1,jang)
     ja2 = deloc%angs(1,2,jang)
     ja3 = deloc%angs(1,3,jang)
     do ishift=-(icenter-1),(icenter-1)
       js1 = deloc%angs(2,1,jang)+ishift
       js2 = deloc%angs(2,2,jang)+ishift
       js3 = deloc%angs(2,3,jang)+ishift

       chkdihed=0
       if (ia2==ja1 .and. is2==js1) then
         if (ia1==ja2 .and. is1==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia3,is3/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja3,js3/)
           chkdihed=1
         else if (ia3==ja2 .and. is3==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia1,is1/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja3,js3/)
           chkdihed=1
         end if
       else if (ia2==ja3 .and. is2==js3) then
         if (ia1==ja2 .and. is1==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia3,is3/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja1,js1/)
           chkdihed=1
         else if (ia3==ja2 .and. is3==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia1,is1/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja1,js1/)
           chkdihed=1
         end if
       end if
       if (deloc%ndihed > 6*deloc%nang) then
         MSG_ERROR('make_dihedrals : too many dihedrals found > 6*nang')
       end if
       if (chkdihed == 1) then
         if (   diheds_tmp(1,4,deloc%ndihed) == diheds_tmp(1,1,deloc%ndihed) .and.&
&         diheds_tmp(2,4,deloc%ndihed) == diheds_tmp(2,1,deloc%ndihed) ) then
           write(std_out,*) 'make_dihedrals : Bad dihedral was found: atom1 == atom4. Discarding.'
           diheds_tmp(:,:,deloc%ndihed) = 0
           deloc%ndihed = deloc%ndihed-1
         end if
       end if
     end do
   end do
!  end jang do
 end do
!end iang do

 if (allocated(deloc%dihedrals)) then
   ABI_FREE(deloc%dihedrals)
 end if

 ABI_ALLOCATE(deloc%dihedrals,(2,4,deloc%ndihed))
 do idihed=1,deloc%ndihed
   deloc%dihedrals(:,:,idihed) = diheds_tmp(:,:,idihed)

!  minshift = minval(diheds_tmp(2,:,idihed))
!  if (minshift <= 0) then
!  deloc%dihedrals(2,:,idihed) = deloc%dihedrals(2,:,idihed)+minshift+1
!  end if
!  maxshift = maxval(diheds_tmp(2,:,idihed))
!  if (maxshift > deloc%nrshift) then
!  deloc%dihedrals(2,:,idihed) = deloc%dihedrals(2,:,idihed)-maxshift
!  end if
!
   minshift = minval(diheds_tmp(2,:,idihed))
   maxshift = maxval(diheds_tmp(2,:,idihed))
   if (minshift <= 0 .or. maxshift > deloc%nrshift) then
     write(std_out,*) ' make_dihedrals : Error : dihedral extends beyond '
     write(std_out,*) '  first neighboring unit cells ! '
     MSG_ERROR("Aborting now")
   end if
 end do
 ABI_DEALLOCATE(diheds_tmp)

end subroutine make_dihedrals
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/make_bonds
!! NAME
!! make_bonds
!!
!! FUNCTION
!!  (to be completed)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      m_abimover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine make_bonds(deloc,natom,ntypat,icenter,rprimd,typat,xcart,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_bonds'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter
 integer,intent(in) :: natom,ntypat
 type(delocint),intent(inout) :: deloc
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: znucl(:) ! znucl(ntypat) or
                                 ! znucl(npsp) ?
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)

!Local variables ------------------------------
!scalars
 integer :: iatom,ibond,irshift,itypat,jatom
 real(dp) :: bl,bondfudge,rcov1,rcov2
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: bonds_tmp(:,:,:)
 real(dp) :: rcov(ntypat),rpt(3)

!************************************************************************

 do itypat=1,ntypat
   call atomdata_from_znucl(atom,znucl(itypat))
   rcov(itypat) = atom%rcov
 end do

!write(std_out,*) ' rcov =', rcov
!write(std_out,*) ' nrshift =', deloc%nrshift
!write(std_out,*) ' xcart =', xcart
!write(std_out,*) ' natom =',natom

!tentative first allocation: < 12 bonds per atom.
 ABI_ALLOCATE(bonds_tmp,(2,2,12*natom))

 bondfudge = 1.1_dp

 deloc%nbond = 0

 do iatom=1,natom
   rcov1 = rcov(typat(iatom))
   do jatom=iatom+1,natom
     rcov2 = rcov(typat(jatom))
     do irshift=1,deloc%nrshift
       rpt(:) = xcart(:,jatom) &
&       + deloc%rshift(1,irshift)*rprimd(:,1) &
&       + deloc%rshift(2,irshift)*rprimd(:,2) &
&       + deloc%rshift(3,irshift)*rprimd(:,3)
       bl =  bond_length(xcart(:,iatom),rpt)

       !write(std_out,*) ' bl, bondfudge*(rcov1+rcov2) = ',bl, bondfudge*(rcov1+rcov2)

       if (bondfudge*(rcov1+rcov2) - bl > tol6) then
         deloc%nbond = deloc%nbond+1
         if (deloc%nbond > 12*natom) then
           MSG_ERROR('make_bonds: error too many bonds !')
         end if
         bonds_tmp(1,1,deloc%nbond) = iatom
         bonds_tmp(2,1,deloc%nbond) = icenter
         bonds_tmp(1,2,deloc%nbond) = jatom
         bonds_tmp(2,2,deloc%nbond) = irshift

         !write(std_out,*) ' ibond bonds = ', deloc%nbond, bonds_tmp(:,:,deloc%nbond),xcart(:,iatom),rpt
       end if
     end do ! jatom
   end do
 end do ! iatom

 if (allocated(deloc%bonds)) then
   ABI_FREE(deloc%bonds)
 end if

 ABI_ALLOCATE(deloc%bonds,(2,2,deloc%nbond))
 do ibond=1,deloc%nbond
   deloc%bonds(:,:,ibond) = bonds_tmp(:,:,ibond)
 end do
 ABI_DEALLOCATE(bonds_tmp)

!do ibond=1,deloc%nbond
!write(std_out,*) 'bond ', ibond, deloc%bonds(:,:,ibond)
!end do

end subroutine make_bonds
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/calc_prim_int
!! NAME
!! calc_prim_int
!!
!! FUNCTION
!!  calculate values of primitive internal coordinates as a function of
!!  cartesian ones.
!!
!! INPUTS
!! angs= number of angles
!! bonds(2,2,nbond)=for a bond between iatom and jatom
!!              bonds(1,1,nbond) = iatom
!!              bonds(2,1,nbond) = icenter
!!              bonds(1,2,nbond) = jatom
!!              bonds(2,2,nbond) = irshift
!! carts(2,ncart)= index of total primitive internal, and atom (carts(2,:))
!! dihedrals(2,4,ndihed)=indexes to characterize dihedrals
!! dtset <type(dataset_type)>=all input variables for this dataset
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of cartesian coordinates used
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! prim_int(ninternal)=values of primitive internal coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      pred_delocint,xcart2deloc
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine calc_prim_int(deloc,natom,rprimd,xcart,prim_int)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_prim_int'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(inout) :: deloc
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(out) :: prim_int(deloc%ninternal)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,i4,iang,ibond,icart,idihed,iprim,s1,s2,s3,s4
!arrays
 real(dp) :: r1(3),r2(3),r3(3),r4(3)

!************************************************************************

!DEBUG
!write(std_out,*) ' calc_prim_int : enter'
!write(std_out,*) shape(angs),shape(deloc%bonds),shape(dihedrals)
!do ibond=1,deloc%nbond
!do i1=1,2
!write(std_out,'(2I5)') deloc%bonds(:,i1,ibond)
!end do
!end do
!ENDDEBUG

 iprim=1
!first: bond values
 do ibond=1,deloc%nbond
   i1 = deloc%bonds(1,1,ibond)
   s1 = deloc%bonds(2,1,ibond)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%bonds(1,2,ibond)
   s2 = deloc%bonds(2,2,ibond)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   prim_int(iprim) = bond_length(r1,r2)
   iprim=iprim+1
 end do

!second: angle values (ang)
 do iang=1,deloc%nang
   i1 = deloc%angs(1,1,iang)
   s1 = deloc%angs(2,1,iang)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%angs(1,2,iang)
   s2 = deloc%angs(2,2,iang)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   i3 = deloc%angs(1,3,iang)
   s3 = deloc%angs(2,3,iang)
   r3(:) = xcart(:,i3)+deloc%rshift(1,s3)*rprimd(:,1)&
&   +deloc%rshift(2,s3)*rprimd(:,2)&
&   +deloc%rshift(3,s3)*rprimd(:,3)
   prim_int(iprim) = angle_ang(r1,r2,r3)
   iprim=iprim+1
 end do

!third: dihedral values
 do idihed=1,deloc%ndihed
   i1 = deloc%dihedrals(1,1,idihed)
   s1 = deloc%dihedrals(2,1,idihed)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%dihedrals(1,2,idihed)
   s2 = deloc%dihedrals(2,2,idihed)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   i3 = deloc%dihedrals(1,3,idihed)
   s3 = deloc%dihedrals(2,3,idihed)
   r3(:) = xcart(:,i3)+deloc%rshift(1,s3)*rprimd(:,1)&
&   +deloc%rshift(2,s3)*rprimd(:,2)&
&   +deloc%rshift(3,s3)*rprimd(:,3)
   i4 = deloc%dihedrals(1,4,idihed)
   s4 = deloc%dihedrals(2,4,idihed)
   r4(:) = xcart(:,i4)+deloc%rshift(1,s4)*rprimd(:,1)&
&   +deloc%rshift(2,s4)*rprimd(:,2)&
&   +deloc%rshift(3,s4)*rprimd(:,3)
   prim_int(iprim) = angle_dihedral(r1,r2,r3,r4)
   iprim=iprim+1
 end do

 do icart=1,deloc%ncart
   prim_int(iprim) = xcart(deloc%carts(1,icart),deloc%carts(2,icart))
   iprim=iprim+1
 end do

!DEBUG
!write(std_out,*) 'Primitive internal coordinate values:'
!do iprim=1,ninternal
!if (iprim <= deloc%nbond) then
!write(std_out,*) iprim, prim_int(iprim)
!else
!write(std_out,*) iprim, prim_int(iprim), prim_int(iprim)/pi*180.0_dp
!end if
!end do
!ENDDEBUG

end subroutine calc_prim_int
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/bond_length
!! NAME
!! bond_length
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function bond_length(r1,r2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bond_length'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: bond_length
!arrays
 real(dp),intent(in) :: r1(3),r2(3)

!Local variables ------------------------------------
!arrays
 real(dp) :: rpt(3)

!******************************************************************
 rpt(:) = r1(:)-r2(:)
 bond_length = sqrt(rpt(1)**2+rpt(2)**2+rpt(3)**2)

end function bond_length
!!***

!!****f* m_abimover/angle_ang
!! NAME
!! angle_ang
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function angle_ang(r1,r2,r3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'angle_ang'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: angle_ang
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n2
!arrays
 real(dp) :: rpt12(3),rpt32(3)

!******************************************************************
 n1=bond_length(r1,r2)
 n2=bond_length(r3,r2)

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)

 cos_ang = (rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3))/n1/n2

 if (cos_ang > one - epsilon(one)*two) then
   cos_ang = one
 else if(cos_ang < -one + epsilon(one)*two) then
   cos_ang = -one
 end if

 angle_ang=acos(cos_ang)

end function angle_ang
!!***

!!****f* m_abimover/angle_dihedral
!! NAME
!! angle_dihedral
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

 function angle_dihedral(r1,r2,r3,r4)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'angle_dihedral'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: angle_dihedral
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)

!Local variables------------------------------------
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,sin_dihedral
!arrays
 real(dp) :: cp1232(3),cp3432(3),cpcp(3),rpt12(3),rpt32(3),rpt34(3)

!******************************************************************

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)
 rpt34(:) = r3(:)-r4(:)

 call acrossb(rpt12,rpt32,cp1232)
 call acrossb(rpt34,rpt32,cp3432)

!DEBUG
!write(std_out,*) ' cos_dihedral : cp1232 = ', cp1232
!write(std_out,*) ' cos_dihedral : cp3432 = ', cp3432
!ENDDEBUG

 n1 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 n2 = sqrt(cp3432(1)**2+cp3432(2)**2+cp3432(3)**2)

 cos_dihedral = (cp1232(1)*cp3432(1)+cp1232(2)*cp3432(2)+cp1232(3)*cp3432(3))/n1/n2
!we use complementary of standard angle, so
 cos_dihedral = -cos_dihedral

 call acrossb(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
!if (abs(sin_dihedral) > tol12) then
!dih_sign = sin_dihedral/abs(sin_dihedral)
!end if
 if (sin_dihedral < -tol12) then
   dih_sign = -one
 end if

!DEBUG
!write(std_out,'(a,3E20.10)') 'angle_dihedral : cos sin dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

 if (cos_dihedral > one - epsilon(one)*two) then
   cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
   cos_dihedral = -one
 end if

 angle_dihedral = dih_sign*acos(cos_dihedral)

end function angle_dihedral
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/make_bonds_new
!! NAME
!! make_bonds_new
!!
!! FUNCTION
!!  Fill the contents of the bonds structure, that contains
!!  all non redundant bonds that could be generated between
!!  all the atoms in the unitary cell and their adjacent cells
!!
!! INPUTS
!!  natom=  Number of atoms
!!  ntypat= Number of type of atoms
!!  rprimd= Dimensional primitive vectors of the cell
!!  xcart=  Cartesian coordinates of the atoms
!!  znucl=  Z number of the atom
!!
!! OUTPUT
!!  bonds= Structure that store all the information about
!!         bonds created by this routine:
!!         nbonds=  Total number of bonds
!!         nbondi=  Number of bonds for atom i
!!         indexi=  Indeces of bonds for atom i
!!         bond_length=  Distances between atoms i and j (including shift)
!!         bond_vect=    Unitary vector for the bond from i to j
!!         tolerance=    The tolerance is multiplied to the
!!                       adition of covalent radius to decide if a bond is created
!!
!! PARENTS
!!      prec_simple
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine make_bonds_new(bonds,natom,ntypat,rprimd,typat,xcart,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_bonds_new'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom,ntypat
!arrays
integer,intent(in) :: typat(natom)
real(dp),intent(in) :: znucl(ntypat)
real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
type(go_bonds),intent(inout) :: bonds

!Local variables ------------------------------
!scalars
integer :: ii,jj,kk,ibond,irshift
real(dp) :: rcov1,rcov2
real(dp) :: bl
type(go_bonds) :: bonds_tmp
type(atomdata_t) :: atom

!arrays
character(len=2) :: symbol(ntypat)
real(dp) :: amu(ntypat)
integer :: shift(3,13) ! Represent all shift vectors that are not equivalent by central symmetry
! For example (1,1,-1) is equivalent to (-1,-1,1)
! It means that bond between atom i in the original cell and atom j in the
! cell with cordinates (1,1,-1) is equivalent to the bond between atom j in
! the orignal cell and atom i in the cell with coordinates (-1,-1,1)
! The trivial shift (0,0,0) is excluded here
real(dp) :: rcov(ntypat) ! Covalent radius
real(dp) :: rpt(3)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

!write(std_out,*) 'make_bonds 01'
!##########################################################
!### 01. Compute covalent radius

 do ii=1,ntypat
   call atomdata_from_znucl(atom,znucl(ii))
   amu(ii) = atom%amu
   rcov(ii) = atom%rcov
   symbol(ii) = atom%symbol
 end do

!write(std_out,*) 'make_bonds 02'
!##########################################################
!### 02. Fill the 13 posible shift conecting adjacent cells

 shift(:,:)=reshape( (/ 1,0,0,&
& 0, 1, 0,&
& 0, 0, 1,&
& 1, 1, 0,&
& 1,-1, 0,&
& 0, 1, 1,&
& 0, 1,-1,&
& 1, 0, 1,&
& 1, 0,-1,&
& 1, 1, 1,&
& 1,-1, 1,&
& 1, 1,-1,&
& 1,-1,-1 /), (/ 3, 13 /))

!write(std_out,*) 'make_bonds 03'
!##########################################################
!### 03. Initialize the values of bonds

!The total number of bonds could not be predicted without
!compute all the distances, but the extreme case is linking
!all the atoms within all adjacent cells (natom*natom*13)
!plus the all the bonds inside the original cell (natom*(natom-1))

 bonds_tmp%nbonds=0
 bonds_tmp%tolerance=bonds%tolerance
 ibond=0

 ABI_ALLOCATE(bonds_tmp%bond_vect,(3,natom*natom*14-natom))
 ABI_ALLOCATE(bonds_tmp%bond_length,(natom*natom*14-natom))

!indexi contains the indeces to the bonds
 ABI_ALLOCATE(bonds_tmp%indexi,(natom,natom*natom*14-natom))

 ABI_ALLOCATE(bonds_tmp%nbondi,(natom))

 bonds_tmp%indexi(:,:)=0
 bonds_tmp%nbondi(:)=0
 bonds_tmp%bond_vect(:,:)=0.0
 bonds_tmp%bond_length(:)=0.0

!write(std_out,*) 'make_bonds 04'
!##########################################################
!### 04. Compute the bonds inside the original cell
!### shift=(0,0,0)

 do ii=1,natom
   rcov1 = rcov(typat(ii))

   do jj=ii+1,natom
     rcov2 = rcov(typat(jj))

     bl=bond_length(xcart(:,ii),xcart(:,jj))

     if (bonds_tmp%tolerance*(rcov1+rcov2) > bl) then
!      We have a new bond, nbonds starts from
!      0, so it could be used to index the
!      locations of bondij and distij

!      Increase the number of bonds
       bonds_tmp%nbonds= bonds_tmp%nbonds+1

!      The number of bonds for atoms ii and jj
!      needs to raise by one
       bonds_tmp%nbondi(ii)= bonds_tmp%nbondi(ii)+1
       bonds_tmp%nbondi(jj)= bonds_tmp%nbondi(jj)+1

       bonds_tmp%indexi(ii,bonds_tmp%nbondi(ii))=bonds_tmp%nbonds
!      The value for jj is negative to indicate that
!      the vector is from ii to jj
       bonds_tmp%indexi(jj,bonds_tmp%nbondi(jj))=-bonds_tmp%nbonds

!      The unitary vector is always from ii to jj
       bonds_tmp%bond_vect(:,bonds_tmp%nbonds)=(xcart(:,jj)-xcart(:,ii))/bl
       bonds_tmp%bond_length(bonds_tmp%nbonds)=bl

     end if

   end do !! jj
 end do !! ii

!write(std_out,*) 'make_bonds 05'
!##########################################################
!### 05. Compute the bonds outside the original cell
!###     13 shifts considered

!Bonds between identical atoms but in diferent cells are
!allowed

 do ii=1,natom
   rcov1 = rcov(typat(ii))
   do jj=1,natom
     rcov2 = rcov(typat(jj))

     do irshift=1,13

       do kk=1,3
         rpt(kk) = xcart(kk,jj)+&
&         shift(1,irshift)*rprimd(kk,1)+ &
&         shift(2,irshift)*rprimd(kk,2)+ &
&         shift(3,irshift)*rprimd(kk,3)
       end do


       bl =bond_length(xcart(:,ii),rpt)

       if (bonds_tmp%tolerance*(rcov1+rcov2) > bl) then

!        We have a new bond, nbonds starts from
!        0, so it could be used to index the
!        locations of bondij and distij

!        Increase the number of bonds
         bonds_tmp%nbonds= bonds_tmp%nbonds+1

!        The number of bonds for atoms ii and jj
!        needs to raise by one
         bonds_tmp%nbondi(ii)= bonds_tmp%nbondi(ii)+1
         bonds_tmp%indexi(ii,bonds_tmp%nbondi(ii))=bonds_tmp%nbonds

!        The value for jj is negative to indicate that
!        the vector is from ii to jj
         bonds_tmp%nbondi(jj)= bonds_tmp%nbondi(jj)+1
         bonds_tmp%indexi(jj,bonds_tmp%nbondi(jj))=-bonds_tmp%nbonds

!        The unitary vector is always from ii to jj
         bonds_tmp%bond_vect(:,bonds_tmp%nbonds)=(rpt(:)-xcart(:,ii))/bl
         bonds_tmp%bond_length(bonds_tmp%nbonds)=bl

         if (ii==jj) then
           bonds_tmp%nbonds= bonds_tmp%nbonds+1
         end if

       end if

     end do !! irshift

   end do !! jj
 end do !! ii

 call print_bonds(amu,bonds_tmp,natom,ntypat,symbol,typat,znucl)


!write(std_out,*) 'make_bonds 05'
!##########################################################
!### 05. Deallocate all the arrays inside bonds
!###     allocate them with the right size and fill them

 call bonds_free(bonds)

 bonds%nbonds=bonds_tmp%nbonds

 if (bonds%nbonds>0) then
!  Allocate the arrays with exactly the rigth nbonds
   ABI_ALLOCATE(bonds%bond_vect,(3,bonds%nbonds))
   ABI_ALLOCATE(bonds%bond_length,(bonds%nbonds))
   ABI_ALLOCATE(bonds%indexi,(natom,bonds%nbonds))
   ABI_ALLOCATE(bonds%nbondi,(natom))

!  Fill the values
   bonds%bond_vect(:,1:bonds%nbonds)=bonds_tmp%bond_vect(:,1:bonds%nbonds)
   bonds%bond_length(1:bonds%nbonds)=bonds_tmp%bond_length(1:bonds%nbonds)
   bonds%indexi(:,1:bonds%nbonds)=bonds_tmp%indexi(:,1:bonds%nbonds)
   bonds%nbondi(:)=bonds_tmp%nbondi(:)
 end if

 call bonds_free(bonds_tmp)

end subroutine make_bonds_new
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/bonds_free
!! NAME
!! bonds_free
!!
!! FUNCTION
!!  Free memory
!!
!! PARENTS
!!      m_abimover,prec_simple
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine bonds_free(bonds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bonds_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(go_bonds),intent(inout) :: bonds

! *********************************************************************

 if (allocated(bonds%bond_vect))then
   ABI_DEALLOCATE(bonds%bond_vect)
 end if

 if (allocated(bonds%bond_length))then
   ABI_DEALLOCATE(bonds%bond_length)
 end if

 if (allocated(bonds%nbondi))then
   ABI_DEALLOCATE(bonds%nbondi)
 end if

 if (allocated(bonds%indexi))then
   ABI_DEALLOCATE(bonds%indexi)
 end if

end subroutine bonds_free
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/print_bonds
!! NAME
!! print_bonds
!!
!! FUNCTION
!!  Print the bonds
!!
!! INPUTS
!!  natom=  Number of atoms
!!  ntypat= Number of type of atoms
!!  znucl=  Z number of the atom
!!
!! OUTPUT
!!  bonds= Structure that store all the information about
!!         bonds created by this routine:
!!         nbonds=  Total number of bonds
!!         bondij=  Unitary vector along the bond direction
!!         distij=  Distances between atoms i and j (including shift)
!!         listij= Indices of bonds going from i to j
!!         listji= Indices of bonds going from j to i
!!         indexij= Number of bonds between i and j
!!         indexji= Number of bonds between j and i
!!         tolerance
!!
!! PARENTS
!!      m_abimover
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

subroutine print_bonds(amu,bonds,natom,ntypat,symbol,typat,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_bonds'
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: natom,ntypat
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(in) :: amu(ntypat)
 character(len=2),intent(in) :: symbol(ntypat)
 type(go_bonds),intent(in) :: bonds

 !Local variables ------------------------------
 !scalars
 integer :: ii,jj,kk

! *********************************************************************

 write(std_out,'(a)') ch10
 write(std_out,'(a,72a,a)') '---BONDS',('-',kk=1,72),ch10
 write(std_out,'(a,i3)') 'Number of atoms:   ',natom
 write(std_out,'(a,i3)') 'Number of bonds:   ',bonds%nbonds
 write(std_out,'(a,f6.3,a,a)') 'Tolerance of bonds: ',bonds%tolerance,' times the sum of covalent radius',ch10

 do ii=1,natom
   write(std_out,'(a,i3)') 'ATOM number:       ',ii
   write(std_out,'(a,f8.3)') '  Z:              ',znucl(typat(ii))
   write(std_out,'(a,f8.3)') '  Weight:         ',amu(typat(ii))
   write(std_out,'(a,a3)') '  Symbol:          ',symbol(typat(ii))
   write(std_out,'(a,i3)') '  Number of bonds: ',bonds%nbondi(ii)

   do jj=1,bonds%nbondi(ii)
     write(std_out,'(a,i3,a,a,i3,a,3f7.3,a,f7.3)') '    [',jj,']',&
&     '    Index of bond: ',bonds%indexi(ii,jj),&
&     '    Unitary vector: ',bonds%bond_vect(:,abs(bonds%indexi(ii,jj))),&
&     '    Bond length: ',bonds%bond_length(abs(bonds%indexi(ii,jj)))
   end do

 end do

 do ii=1,bonds%nbonds

   write(std_out,'(a,i3)') 'BOND Index=',ii
   write(std_out,'(a,3f8.3)') '    Vector',bonds%bond_vect(:,ii)
   write(std_out,'(a,f8.3)')  '    bond Length',bonds%bond_length(ii)

 end do

end subroutine print_bonds
!!***

!----------------------------------------------------------------------

!!****f* m_abimover/make_angles_new
!! NAME
!! make_angles_new
!!
!! FUNCTION
!!  Fill the contents of the angles structure, that contains
!!  all non redundant angles that could be generated between
!!  all the atoms in the unitary cell and their adjacent cells
!!  An angle is establish when an atom has two or more bonds.
!!  The angles structure contains information about the atoms
!!  involved, the value of the angle in radians, and the unitary
!!  vector perpendicular to the plane of the three atoms that
!!  build the angle.
!!
!! INPUTS
!!  natom=  Number of atoms
!!  ntypat= Number of type of atoms
!!  rprimd= Dimensional primitive vectors of the cell
!!  xcart=  Cartesian coordinates of the atoms
!!  znucl=  Z number of the atom
!!  bonds= Structure that store all the information about
!!         bonds created by this routine:
!!         nbonds=  Total number of bonds
!!         nbondi=  Number of bonds for atom i
!!         indexi=  Indeces of bonds for atom i
!!         bond_length=  Distances between atoms i and j (including shift)
!!         bond_vect=    Unitary vector for the bond from i to j
!!         tolerance=    The tolerance is multiplied to the
!!                       adition of covalent radius to decide if a bond is created
!!
!! OUTPUT
!!  angles=  Structure that store the information about
!!           angles created by this routine
!!          nangles= Total number of angles
!!          angle_vertex= Index of the atom for that angle
!!          angle_value= Value of the angle in radians
!!          angle_bonds=  Indices of the bonds
!!          angle_vect=   Unitary vector perpendicular to the plane of the angle
!!
!! PARENTS
!!
!! CHILDREN
!!      atomdata_from_znucl,bonds_free,print_bonds
!!
!! SOURCE

!This routine has been disables since it's broken
#if 0

subroutine make_angles_new(angles,bonds,natom,ntypat,rprimd,typat,xcart,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_angles_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom,ntypat
!arrays
integer,intent(in) :: typat(natom)
real(dp),intent(in) :: znucl(ntypat)
real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
type(go_bonds),intent(in) :: bonds
type(go_angles),intent(inout) :: angles

!Local variables ------------------------------
!scalars
integer :: ii,jj,kk,iangle
type(atomdata_t) :: atom

!arrays
type(go_bonds) :: bonds_tmp
character(len=2) :: symbol(ntypat)
real(dp) :: amu(ntypat)
integer :: shift(3,13) ! Represent all shift vectors that are not equivalent by central symmetry
! For example (1,1,-1) is equivalent to (-1,-1,1)
! It means that bond between atom i in the original cell and atom j in the
! cell with cordinates (1,1,-1) is equivalent to the bond between atom j in
! the orignal cell and atom i in the cell with coordinates (-1,-1,1)
! The trivial shift (0,0,0) is excluded here
real(dp) :: rcov(ntypat) ! Covalent radius
real(dp) :: rpt(3)

!***************************************************************************
!Beginning of executable session
!***************************************************************************
 MSG_ERROR("This routine is not tested")

!write(std_out,*) 'make_bonds 01'
!##########################################################
!### 01. Compute covalent radius

 do ii=1,ntypat
   call atomdata_from_znucl(atom,znucl(ii))
   amu(ii) = atom%amu
   rcov(ii) = atom%rcov
   symbol(ii) = symbol(ii)
 end do

!write(std_out,*) 'make_bonds 02'
!##########################################################
!### 02. Fill the 13 posible shift conecting adjacent cells

 shift(:,:)=reshape( (/ 1,0,0,&
& 0, 1, 0,&
& 0, 0, 1,&
& 1, 1, 0,&
& 1,-1, 0,&
& 0, 1, 1,&
& 0, 1,-1,&
& 1, 0, 1,&
& 1, 0,-1,&
& 1, 1, 1,&
& 1,-1, 1,&
& 1, 1,-1,&
& 1,-1,-1 /), (/ 3, 13 /))

!write(std_out,*) 'make_bonds 03'
!##########################################################
!### 03. Initialize the values of bonds

!The total number of bonds could not be predicted without
!compute all the distances, but the extreme case is linking
!all the atoms within all adjacent cells (natom*natom*13)
!plus the all the bonds inside the original cell (natom*(natom-1))

 bonds_tmp%nbonds=0
 bonds_tmp%tolerance=bonds%tolerance
 ibond=0

 ABI_ALLOCATE(bonds_tmp%bond_vect,(3,natom*natom*14-natom))
 ABI_ALLOCATE(bonds_tmp%bond_length,(natom*natom*14-natom))

!indexi contains the indeces to the bonds
 ABI_ALLOCATE(bonds_tmp%indexi,(natom,natom*natom*14-natom))

 ABI_ALLOCATE(bonds_tmp%nbondi,(natom))

 bonds_tmp%indexi(:,:)=0
 bonds_tmp%nbondi(:)=0

!write(std_out,*) 'make_bonds 04'
!##########################################################
!### 04. Compute the bonds inside the original cell
!### shift=(0,0,0)

 do ii=1,natom
   rcov1 = rcov(typat(ii))

   do jj=ii+1,natom
     rcov2 = rcov(typat(jj))

     bl=bond_length(xcart(:,ii),xcart(:,jj))

     if (bonds_tmp%tolerance*(rcov1+rcov2) > bl) then
!      We have a new bond, nbonds starts from
!      0, so it could be used to index the
!      locations of bondij and distij

!      Increase the number of bonds
       bonds_tmp%nbonds= bonds_tmp%nbonds+1

!      The number of bonds for atoms ii and jj
!      needs to raise by one
       bonds_tmp%nbondi(ii)= bonds_tmp%nbondi(ii)+1
       bonds_tmp%nbondi(jj)= bonds_tmp%nbondi(jj)+1

       bonds_tmp%indexi(ii,bonds_tmp%nbondi(ii))=bonds_tmp%nbonds
!      The value for jj is negative to indicate that
!      the vector is from ii to jj
       bonds_tmp%indexi(jj,bonds_tmp%nbondi(jj))=-bonds_tmp%nbonds

!      The unitary vector is always from ii to jj
       bonds_tmp%bond_vect(:,bonds_tmp%nbonds)=(xcart(:,jj)-xcart(:,ii))/bl
       bonds_tmp%bond_length(bonds_tmp%nbonds)=bl

     end if

   end do !! jj
 end do !! ii

!write(std_out,*) 'make_bonds 05'
!##########################################################
!### 05. Compute the bonds outside the original cell
!###     13 shifts considered

!Bonds between identical atoms but in diferent cells are
!allowed

 do ii=1,natom
   rcov1 = rcov(typat(ii))
   do jj=1,natom
     rcov2 = rcov(typat(jj))

     do irshift=1,13

       do kk=1,3
         rpt(kk) = xcart(kk,jj)+&
&         shift(1,irshift)*rprimd(kk,1)+ &
&         shift(2,irshift)*rprimd(kk,2)+ &
&         shift(3,irshift)*rprimd(kk,3)
       end do


       bl =bond_length(xcart(:,ii),rpt)

       if (bonds_tmp%tolerance*(rcov1+rcov2) > bl) then

!        We have a new bond, nbonds starts from
!        0, so it could be used to index the
!        locations of bondij and distij

!        Increase the number of bonds
         bonds_tmp%nbonds= bonds_tmp%nbonds+1

!        The number of bonds for atoms ii and jj
!        needs to raise by one
         bonds_tmp%nbondi(ii)= bonds_tmp%nbondi(ii)+1
         bonds_tmp%indexi(ii,bonds_tmp%nbondi(ii))=bonds_tmp%nbonds

!        The value for jj is negative to indicate that
!        the vector is from ii to jj
         bonds_tmp%nbondi(jj)= bonds_tmp%nbondi(jj)+1
         bonds_tmp%indexi(jj,bonds_tmp%nbondi(jj))=-bonds_tmp%nbonds

!        The unitary vector is always from ii to jj
         bonds_tmp%bond_vect(:,bonds_tmp%nbonds)=(rpt(:)-xcart(:,ii))/bl
         bonds_tmp%bond_length(bonds_tmp%nbonds)=bl

         if (ii==jj) then
           bonds_tmp%nbonds= bonds_tmp%nbonds+1
         end if

       end if

     end do !! irshift

   end do !! jj
 end do !! ii

 call print_bonds(amu,bonds_tmp,natom,ntypat,symbol,typat,znucl)


!write(std_out,*) 'make_bonds 05'
!##########################################################
!### 05. Deallocate all the arrays inside bonds
!###     allocate them with the right size and fill them

 call bonds_free(bonds)

 bonds%nbonds=bonds_tmp%nbonds

 if (bonds%nbonds>0) then
!  Allocate the arrays with exactly the rigth nbonds
   ABI_ALLOCATE(bonds%bond_vect,(3,bonds%nbonds))
   ABI_ALLOCATE(bonds%bond_length,(bonds%nbonds))
   ABI_ALLOCATE(bonds%indexi,(natom,bonds%nbonds))
   ABI_ALLOCATE(bonds%nbondi,(natom))

!  Fill the values
   bonds%bond_vect(:,1:bonds%nbonds)=bonds_tmp%bond_vect(:,1:bonds%nbonds)
   bonds%bond_length(1:bonds%nbonds)=bonds_tmp%bond_length(1:bonds%nbonds)
   bonds%indexi(:,1:bonds%nbonds)=bonds_tmp%indexi(:,1:bonds%nbonds)
   bonds%nbondi(:)=bonds_tmp%nbondi(:)
 end if

 call bonds_free(bonds_tmp)

end subroutine make_angles_new
!!***

#endif

!----------------------------------------------------------------------

end module m_abimover
