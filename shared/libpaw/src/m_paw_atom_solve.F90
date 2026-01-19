!!****m* ABINIT/m_paw_atom_solve
!! NAME
!!  m_paw_atom_solve
!!
!! FUNCTION
!! This module provides a modified version of atompaw (created by NAWH, MT, FJ) that is needed to implement relaxed core paw.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2025 ABINIT group (MT,NBrouwer, JBoust)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!!  Several functions and types in this module are declared to be private to avoid
!!  problems with the rest of libpaw while maintaining as much of atompaw code as possible.
!!
!! SOURCE

#include "libpaw.h"

module m_paw_atom_solve

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 USE ieee_arithmetic
 use m_libpaw_libxc
 use m_pawtab
 use m_pawrad
 use m_paw_atomorb
 use m_paw_atom,     only : atompaw_ehnzc,atompaw_dij0,atompaw_kij,atompaw_vhnzc
 use m_paw_numeric
 use m_pawpsp
 use m_libpaw_tools, only : libpaw_get_free_unit

 implicit none
 private





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PRIVATE PARAMETERS IMPORTED FROM ATOMPAW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
 ! Numerics
 logical,private,save :: has_to_print
 real(dp),private,save :: machine_zero,machine_precision,machine_infinity,practical_zero
 real(dp), PRIVATE,save :: minlog,maxlog,minexp,maxexp 
 real(dp), PRIVATE,save :: minlogarg,maxlogarg,minexparg,maxexparg
 real(dp), PARAMETER, PRIVATE :: MaxMix=0.5_dp,seterr=tol11,settoosmall=tol16
 REAL(dp), parameter  :: seterrmg=1.d-9,settoosmallmg=1.d-10
 integer, parameter, private :: MaxIter=2000
 ! Grid
 INTEGER, PARAMETER, PRIVATE :: lineargrid=1  ! r(i)=h*(i-1)
 INTEGER, PARAMETER, PRIVATE :: loggrid=2     ! r(i)=r0*(exp(h*(i-1))-1)
 real(dp),PRIVATE, PARAMETER :: coretailtol=tol12
 real(dp), PARAMETER, PRIVATE ::linrange=50._dp,linh=0.0025_dp,mxgridlin=20001
 real(dp), PARAMETER, PRIVATE ::logrange=80._dp,logh=0.020_dp,mxgridlog=2001
 real(dp), PARAMETER, PRIVATE :: v4logrange=100._dp,lor00=tol5
 ! Constants
 real(dp), parameter :: ifsalpha2=InvFineStruct**2
 real(dp), parameter :: fsalpha2=1._dp/InvFineStruct**2
 ! PS scheme
 INTEGER,PRIVATE,PARAMETER :: BLOECHL=1, VANDERBILT=2, CUSTOM=3, MODRRKJ=7,HFPROJ=8
 INTEGER,PRIVATE,PARAMETER :: BLOECHLPS=0, POLYNOM=1,POLYNOM2=2,RRKJ=3,MARSMAN=6
 INTEGER,PRIVATE,PARAMETER :: VANDERBILTORTHO=0, GRAMSCHMIDTORTHO=1
 INTEGER,PRIVATE,PARAMETER :: SVDORTHO=2, HFORTHO=-13
 INTEGER,PRIVATE,PARAMETER :: MTROULLIER=1, ULTRASOFT=2, BESSEL=3, KERKER_E=4,KERKER_P=5
 INTEGER,PRIVATE,PARAMETER :: HARTREE_FOCK=6, SETVLOC=7, VPSMATCHNC=8,VPSMATCHNNC=9
 INTEGER,PARAMETER,PRIVATE :: norbit_max=50,nbasis_add_max=25
 real(dp),PARAMETER,PRIVATE :: logder_min=-5._dp,logder_max=4.95_dp,logder_pts=200
 real(dp),PARAMETER,PRIVATE :: polynom2_pdeg_def=4
 real(dp),PARAMETER,PRIVATE :: polynom2_qcut_def=10._dp
 real(dp),PARAMETER,PRIVATE :: gausstol_def=tol4
 real(dp),PARAMETER,PRIVATE :: hf_coretol_def=tol4
 INTEGER,PARAMETER,PRIVATE :: PROJECTOR_TYPE_BLOECHL   = 1
 INTEGER,PARAMETER,PRIVATE :: PROJECTOR_TYPE_VANDERBILT= 2
 INTEGER,PARAMETER,PRIVATE :: PROJECTOR_TYPE_CUSTOM    = 3
 INTEGER,PARAMETER,PRIVATE :: PROJECTOR_TYPE_MODRRKJ   = 4
 INTEGER,PARAMETER,PRIVATE :: PROJECTOR_TYPE_HF        = 5
 INTEGER,PARAMETER,PRIVATE :: PROJECTOR_TYPE_MARSMAN   = 6
 INTEGER,PARAMETER,PRIVATE :: PSEUDO_TYPE_BLOECHL      = 1
 INTEGER,PARAMETER,PRIVATE :: PSEUDO_TYPE_POLYNOM      = 2
 INTEGER,PARAMETER,PRIVATE :: PSEUDO_TYPE_POLYNOM2     = 3
 INTEGER,PARAMETER,PRIVATE :: PSEUDO_TYPE_RRKJ         = 4
 INTEGER,PARAMETER,PRIVATE :: PSEUDO_TYPE_BLOECHL_K    = 5
 INTEGER,PARAMETER,PRIVATE :: PSEUDO_TYPE_HF           = 6
 INTEGER,PARAMETER,PRIVATE :: ORTHO_TYPE_GRAMSCHMIDT   = 1
 INTEGER,PARAMETER,PRIVATE :: ORTHO_TYPE_VANDERBILT    = 2
 INTEGER,PARAMETER,PRIVATE :: ORTHO_TYPE_SVD           = 3
 INTEGER,PARAMETER,PRIVATE :: ORTHO_TYPE_HF            = 4
 INTEGER,PARAMETER,PRIVATE :: SHAPEFUNC_TYPE_GAUSSIAN  = 1
 INTEGER,PARAMETER,PRIVATE :: SHAPEFUNC_TYPE_SINC      = 2
 INTEGER,PARAMETER,PRIVATE :: SHAPEFUNC_TYPE_BESSEL    = 3
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_MTROULLIER     = 1
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_ULTRASOFT      = 2
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_BESSEL         = 3
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_SETVLOC        = 4
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_KERKER_EXPF    = 5
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_KERKER_POLY    = 6
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_VPSMATCHNC     = 7
 INTEGER,PARAMETER,PRIVATE :: VLOC_TYPE_VPSMATCHNNC    = 8
 INTEGER,PARAMETER,PRIVATE :: UNKNOWN_TYPE             =-1
 ! XC
 INTEGER, PRIVATE, PARAMETER :: LDA_PW=14
 INTEGER, PRIVATE, PARAMETER :: GGA_PBE=16
 INTEGER, PRIVATE, PARAMETER :: GGA_PBESOL=18
 INTEGER, PRIVATE, PARAMETER :: MGGA_R2SCAN_001=13001
 INTEGER, PRIVATE, PARAMETER :: MGGA_R2SCAN_01=1301
 INTEGER, PRIVATE, PARAMETER :: LIBXC=-1
 REAL(dp), PRIVATE, PARAMETER :: kappa= 0.804_dp
 REAL(dp), PRIVATE, PARAMETER :: muorig = 0.2195149727645171_dp
 REAL(dp), PRIVATE, PARAMETER :: betorig = 0.06672455060314922_dp
 REAL(dp), PRIVATE, PARAMETER :: gamm = 0.03109069086965489503494086371273_dp
 REAL(dp), PRIVATE, PARAMETER :: musol = 0.123456790123456_dp
 REAL(dp), PRIVATE, PARAMETER :: betsol = 0.046_dp
 REAL(dp), PRIVATE, PARAMETER :: AA=0.0310907d0
 REAL(dp), PRIVATE, PARAMETER :: a1=0.21370d0
 REAL(dp), PRIVATE, PARAMETER :: b1=7.59570d0
 REAL(dp), PRIVATE, PARAMETER :: b2=3.58760d0
 REAL(dp), PRIVATE, PARAMETER :: b3=1.63820d0
 REAL(dp), PRIVATE, PARAMETER :: b4=0.49294d0
 REAL(dp), PRIVATE, PARAMETER :: cx0 = 1.d0
 REAL(dp), PRIVATE, PARAMETER :: cx1 = -0.667d0
 REAL(dp), PRIVATE, PARAMETER :: cx2 = -0.4445555d0
 REAL(dp), PRIVATE, PARAMETER :: cx3 = -0.663086601049d0
 REAL(dp), PRIVATE, PARAMETER :: cx4 = 1.451297044490d0
 REAL(dp), PRIVATE, PARAMETER :: cx5 = -0.887998041597d0
 REAL(dp), PRIVATE, PARAMETER :: cx6 = 0.234528941479d0
 REAL(dp), PRIVATE, PARAMETER :: cx7 = -0.023185843322d0
 REAL(dp), PRIVATE, PARAMETER :: SCANc1x = 0.667d0
 REAL(dp), PRIVATE, PARAMETER :: SCANc2x = 0.8d0
 REAL(dp), PRIVATE, PARAMETER :: SCANdx = 1.24d0
 REAL(dp), PRIVATE, PARAMETER :: k0 = 0.174d0
 REAL(dp), PRIVATE, PARAMETER :: k1 = 0.065d0
 REAL(dp), PRIVATE, PARAMETER :: mu = 10.d0/81.d0
 REAL(dp), PRIVATE, PARAMETER :: SCANa1 = 4.9479d0
 REAL(dp) :: eta
 REAL(dp), PRIVATE, PARAMETER :: dp2 = 0.361d0
 REAL(dp) :: C2Ceta
 REAL(dp), PRIVATE, PARAMETER :: cc0 = 1.d0
 REAL(dp), PRIVATE, PARAMETER :: cc1 = -0.64d0
 REAL(dp), PRIVATE, PARAMETER :: cc2 = -0.4352d0
 REAL(dp), PRIVATE, PARAMETER :: cc3 = -1.535685604549d0
 REAL(dp), PRIVATE, PARAMETER :: cc4 = 3.061560252175d0
 REAL(dp), PRIVATE, PARAMETER :: cc5 = -1.915710236206d0
 REAL(dp), PRIVATE, PARAMETER :: cc6 = 0.516884468372d0
 REAL(dp), PRIVATE, PARAMETER :: cc7 = -0.051848879792d0
 REAL(dp), PRIVATE, PARAMETER :: SCANc1c = 0.64d0
 REAL(dp), PRIVATE, PARAMETER :: SCANc2c = 1.5d0
 REAL(dp), PRIVATE, PARAMETER :: SCANdc = 0.7d0
 REAL(dp), PRIVATE, PARAMETER :: b1c = 0.0285764d0
 REAL(dp), PRIVATE, PARAMETER :: b2c = 0.0889d0
 REAL(dp), PRIVATE, PARAMETER :: b3c = 0.125541d0
 REAL(dp), PRIVATE, PARAMETER :: betaMB = 0.066725d0
 REAL(dp), PRIVATE, PARAMETER :: chiinfinity = 0.12802585262625815d0
 REAL(dp), PRIVATE, PARAMETER :: Sgam = 0.031090690869655d0
 REAL(dp), PRIVATE, PARAMETER :: Dfc2=cc1+2*cc2+3*cc3+4*cc4+5*cc5+6*cc6+7*cc7
 ! Parameters for the Perdew-Wang (PRB 45,13244 (1992)) LDA correlation
 REAL(dp), PRIVATE, PARAMETER :: LDAA = 0.03109070d0
 REAL(dp), PRIVATE, PARAMETER :: LDAa1 = 0.21370d0
 REAL(dp), PRIVATE, PARAMETER :: LDAb1 = 7.59570d0
 REAL(dp), PRIVATE, PARAMETER :: LDAb2 = 3.58760d0
 REAL(dp), PRIVATE, PARAMETER :: LDAb3 = 1.63820d0
 REAL(dp), PRIVATE, PARAMETER :: LDAb4 = 0.49294d0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PRIVATE DATA TYPES IMPORTED FROM ATOMPAW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  GridInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type,private :: GridInfo
  integer :: type
  integer :: n
  integer :: ishift
  real(dp) :: h,r0,range
  real(dp), pointer :: r(:) => null()
  real(dp), pointer :: drdu(:) => null()   ! for loggrid -- dr/du
  real(dp), pointer :: pref(:) => null()   ! for loggrid -- r0*exp(u/2)
  real(dp), pointer :: rr02(:) => null()   ! for loggrid -- (r+r0)**2
end type GridInfo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  OrbitInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type,private :: OrbitInfo
  character(132) :: exctype
  integer :: nps, npp, npd ,npf, npg, norbit
  integer :: npsc, nppc, npdc ,npfc, npgc
  INTEGER, POINTER :: np(:) => null()
  INTEGER, POINTER :: l(:) => null()
  INTEGER, POINTER :: kappa(:) => null()
  real(dp), POINTER :: eig(:) => null()
  real(dp), POINTER :: occ(:) => null()
  real(dp), POINTER :: wfn(:,:) => null()
  real(dp), POINTER :: lwfn(:,:) => null()
  real(dp), POINTER :: otau(:,:) => null() ! kinetic energy density for orbital
  real(dp), POINTER :: lqp(:,:) => null()  ! only used for HF
  real(dp), POINTER :: X(:,:) => null()    ! identical to HF%SumY(:,:)
  LOGICAL, POINTER :: iscore(:) => null()
  LOGICAL, POINTER :: issemicore(:) => null()
  real(dp),POINTER :: den(:) => null() ! accumulated over states
  real(dp),POINTER :: tau(:) => null() ! accumulated over states
  real(dp),POINTER :: deltatau(:) => null() !tau-tauW   (tauW==Weizsacker)
  ! LIBPAW specific
  real(dp),POINTER :: coreden(:) => null() 
  real(dp),POINTER :: valeden(:) => null() 
  real(dp) :: qval
  logical :: frozencorecalculation
  logical :: frozenvalecalculation
  logical :: diracrelativistic
  logical :: scalarrelativistic
end type OrbitInfo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PotentialInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type,private :: PotentialInfo
  character(2) :: sym
  integer :: nz     !  nz is nuclear charge
  real(dp) :: zz        !  zz=nz is nuclear charge
  real(dp) :: q,v0,v0p  !  q is total electron charge
  !  v0,v0p are potential value and deriv at r=0
  real(dp) :: Nv0,Nv0p    !  finite nucleus value and deriv at 0
  real(dp) , pointer :: rv(:) => null()
  real(dp) , pointer :: rvn(:) => null()
  real(dp) , pointer :: rvh(:) => null()
  real(dp) , pointer :: rvx(:) => null()
  !  rv(n) is  veff * r
  !  rvh is hartree potential for den
  !  rvn is nuclear potential
  !  rvx is exchange-correlation potential
  real(dp) , pointer :: vtau(:) => null() !for meta-gga
  integer :: finitenucleusmodel
  ! Based on models 2, 3, 4, 5 discussed by Dirk Anrae ,
  !   Physics Reports 336 (2000) 413-525
  !    default is 0 for previous Gaussian model
  !    for finitenucleusmodel<0, finite nucleus is false
  ! LIBPAW specific
  logical :: finitenucleus
  logical :: needvtau
  real(dp),allocatable :: ww(:)
  real(dp), allocatable :: jj(:)
end type PotentialInfo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SCFInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type,private :: SCFInfo
  integer :: iter
  real(dp) :: delta,eone,ekin,estatic,ecoul,eexc,oepcs,etot
  real(dp) :: valekin,valecoul,valeexc,corekin,evale ! used in frozencore only
end type SCFInfo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Anderson_Context
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE,private :: Anderson_Context  !** Anderson Mixing context
  real(dp)    :: NewMix    !** Amount of new vectors to mix, ie beta in paper.
  INTEGER :: Nmax      !** Max number of vectors to keep
  INTEGER :: N         !** Current number of vectors in list
  INTEGER :: Slot      !** New fill Slot
  INTEGER :: VecSize   !** Size of each vector
  INTEGER :: Err_Unit  !** Error unit
  INTEGER :: MaxIter   !** MaxIter
  INTEGER :: CurIter   !** Running iteration index
  real(dp) :: err       !** residue convergence tolerance
  real(dp) :: toosmall  !** solution obviously converged
  real(dp) :: res       !** Running convergence error
  Logical :: writelots
  real(dp), POINTER :: Matrix(:,:)
  real(dp), POINTER :: Gamma(:)  !** Gamma as defined in 7.6
  real(dp), POINTER :: DF(:,:)   !** Delta F
  real(dp), POINTER :: Fprev(:)
  real(dp), POINTER :: DX(:,:)
  real(dp), POINTER :: Xprev(:)
  ! temporary constants and arrays needed for each call to Anderson_Mix
  INTEGER, POINTER :: IPIV(:)
  real(dp),  POINTER :: S(:)
  real(dp),  POINTER :: RWork(:)
  real(dp), POINTER :: U(:,:)
  real(dp), POINTER :: VT(:,:)
  real(dp), POINTER :: Work(:)
  real(dp), POINTER :: DupMatrix(:,:)
  INTEGER          :: Lwork
  INTEGER          :: LRwork
  real(dp)           :: ConditionNo
  real(dp)           :: MachAccur
END TYPE Anderson_Context


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PseudoInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE,private ::  Pseudoinfo
  CHARACTER(132) :: exctype
  INTEGER  :: lmax,irc,irc_shap,irc_vloc,irc_core,coretailpoints,mesh_size
  INTEGER  :: ivale,itau,ivion,mxbase
  CHARACTER(15) :: orthogonalization_scheme
  CHARACTER(132) :: Vloc_description
  CHARACTER(132) :: Proj_description
  CHARACTER(132) :: Comp_description
  LOGICAL :: multi_rc,poscorenhat
  real(dp) :: rc,rc_shap,rc_vloc,rc_core,energyoflmax,gausslength
  real(dp), POINTER :: rcio(:) => null()
  real(dp), POINTER :: vloc(:) => null()
  real(dp), POINTER :: abinitvloc(:) => null()
  real(dp), POINTER :: abinitnohat(:) => null()
  real(dp), POINTER :: rveff(:) => null()
  real(dp), POINTER :: AErefrv(:) => null()
  real(dp), POINTER :: rvx(:) => null()
  real(dp), POINTER :: trvx(:) => null()
  real(dp), POINTER :: Ktvtau(:) => null()
  real(dp), POINTER :: Krveff(:) => null() 
  real(dp), POINTER :: Kunscreen(:) => null() ! Kresse form
  real(dp), POINTER :: projshape(:) => null()
  real(dp), POINTER :: hatshape(:) => null()
  real(dp), POINTER :: hatden(:) => null()
  real(dp), POINTER :: hatpot(:) => null()
  real(dp), POINTER :: den(:) => null()
  real(dp), POINTER :: tden(:) => null()
  real(dp), POINTER :: core(:) => null() 
  real(dp), POINTER :: tcore(:) => null()
  real(dp), POINTER :: nhatv(:) => null()
  real(dp), POINTER :: coretau(:) => null()
  real(dp), POINTER :: tcoretau(:) => null()
  real(dp), POINTER :: valetau(:) => null()
  real(dp), POINTER :: tvaletau(:) => null()
  real(dp), POINTER :: vtau(:) => null()
  real(dp), POINTER :: tvtau(:) => null()
  INTEGER :: nbase,ncoreshell
  INTEGER, POINTER :: np(:) => null()
  INTEGER, POINTER :: l(:) => null()
  INTEGER, POINTER :: nodes(:) => null() 
  INTEGER, POINTER :: kappa(:) => null()
  INTEGER, POINTER :: rng(:) => null()      ! rng particularly of continuum states
  CHARACTER(8), POINTER :: label(:) => null()
  real(dp), POINTER :: phi(:,:) => null()
  real(dp), POINTER :: tphi(:,:) => null()
  real(dp), POINTER :: tp(:,:) => null() ! before orthog
  real(dp), POINTER :: ophi(:,:) => null()
  real(dp), POINTER :: otphi(:,:) => null()
  real(dp), POINTER :: otp(:,:) => null() ! after orthog
  real(dp), POINTER :: Kop(:,:) => null()   ! for storing K|phi>
  real(dp), POINTER :: eig(:) => null()
  real(dp), POINTER :: occ(:) => null() 
  real(dp), POINTER :: ck(:) => null() 
  real(dp), POINTER :: vrc(:) => null()
  real(dp), POINTER :: oij(:,:) => null() 
  real(dp), POINTER :: dij(:,:) => null() 
  real(dp), POINTER :: wij(:,:) => null()
  !********** modified parameters for use with KS and HF
  real(dp), POINTER :: rVf(:) => null() 
  real(dp), POINTER :: rtVf(:) => null()
  real(dp), POINTER :: g(:,:) => null()
  real(dp), POINTER :: Kij(:,:) => null() 
  real(dp), POINTER :: Vfij(:,:) => null()
  real(dp), POINTER :: mLij(:,:,:) => null()
  real(dp), POINTER :: DR(:,:,:,:,:) => null()
  real(dp), POINTER :: DRVC(:,:,:) => null() 
  real(dp), POINTER :: TXVC(:,:) => null()  ! now output for DFT also
  real(dp) :: lambshielding 
  real(dp) :: XCORECORE    ! output for DFT
  INTEGER, POINTER :: valencemap(:) => null()   ! valencemap({occ. states})={basis}
  Type(OrbitInfo), POINTER :: OCCwfn => null()
  Type(OrbitInfo), POINTER :: TOCCwfn => null()
  real(dp) :: tkin,tion,tvale,txc,Ea,Etotal,Eaion,Eaionhat,Eaxc
  real(dp) :: VlocCoef,VlocRad
  !***********for HF only
  real(dp), POINTER :: lmbd(:,:) => null() !(Eq. 72) lmbd({occ. states},{basis states})
  real(dp), POINTER :: DRC(:,:,:,:) => null() 
  real(dp), POINTER :: mLic(:,:,:) => null()
  real(dp), POINTER :: DRCC(:,:,:,:) => null()
  real(dp), POINTER :: DRCjkl(:,:,:,:,:) => null()
  real(dp), POINTER :: mLcc(:,:,:) => null()
  real(dp), POINTER :: Dcj(:,:) => null()
  real(dp) :: coretol
END  TYPE Pseudoinfo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PseudoInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE,private ::  splinesolvinfo
  INTEGER :: ns      ! # spline nodes, not including origin
  REAL(dp) :: r0,h
  real(dp),allocatable :: u(:)
  real(dp),allocatable :: pref(:)
  real(dp),allocatable :: rr1(:)
  real(dp),allocatable :: rr2(:)
  real(dp),allocatable :: srv(:)
  real(dp),allocatable :: svtau(:)
  real(dp),allocatable :: sdvt(:)
  real(dp),allocatable :: soneplusvt(:)
  real(dp),allocatable :: fvtau(:)
  real(dp),allocatable :: fdvtaudr(:)
  real(dp),allocatable :: frvx(:)
  real(dp),allocatable :: fden(:)
  real(dp),allocatable :: ftau(:)
  type(GridInfo) :: Grids
  type(GridInfo) :: Gridf 
END  TYPE splinesolvinfo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PUBLIC DATATYPES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!----------------------------------------------------------------------
!!****t* m_paw_atom_solve/atompaw_type
!! NAME
!! atompaw_type
!!
!! FUNCTION
!! Stores data in atompaw format, including:
!!    -Input that was used to create the original PAW potentials
!!    -Current spherical average of the valence band electron density
!!    -Previous core charge density to be used as starting point
!!    -All data needed for solving the atomic problem
!!    -All data needed for pseudoization
!!
!! COPYRIGHT
!! This module is
!!
!! SOURCE
 type,public :: atompaw_type
  !!! Initialized in read input roughly in this order
  ! Atom
  CHARACTER(2) :: atomic_symbol    ! Atomic symbol
  INTEGER      :: atomic_charge 
  ! Algo
  logical :: scalarrelativistic
  logical :: diracrelativistic
  logical :: usespline 
  INTEGER :: splns=400             ! Spline interpolation grid length
  real(dp) :: splr0=0.1_dp           ! Spline interpolation r0 value 
  logical :: BDsolve
  logical :: HFpostprocess
  logical :: finitenucleus
  integer :: finitenucleusmodel
  ! Grid 
  CHARACTER(10) :: gridkey
  INTEGER :: gridpoints            ! Number of points of the radial grid
  real(dp) :: gridrange             ! Range of the radial grid
  real(dp) :: gridmatch             ! A matching radius in the radial grid
  real(dp) :: minlogderiv
  real(dp) :: maxlogderiv
  integer :: nlogderiv
  ! XC
  CHARACTER(132) :: exctype        ! Exchange-correlation type (string) 
  real(dp) :: temp                 ! temeprature
  LOGICAL :: needvtau              ! TRUE if Calculation is performed with full kinetic energy functional
  logical :: localizedcoreexchange
  LOGICAL :: fixed_zero            ! Flag activating the "fixed zero" exact exchange potential calculation
  INTEGER :: fixed_zero_index      ! Option for "fixed zero" calculation in the exact exchange potential
  type(libxc_functional_type) :: xc_functionals(2)
  ! Electronic config
  INTEGER :: np(5)                 ! Electronic configuration: number of s,p,d,f,g shells
  INTEGER :: norbit                ! Electronic configuration: number of orbitals
  INTEGER :: norbit_mod            ! Electronic configuration: number of orbitals with modified occupations
  INTEGER,ALLOCATABLE :: orbit_mod_n(:)   ! Electronic config.: n number of the modified orbital
  INTEGER,ALLOCATABLE :: orbit_mod_l(:)   ! Electronic config.: l number of the modified orbital
  INTEGER,ALLOCATABLE :: orbit_mod_k(:)   ! Electronic config.: kappa number of the modified orbital
  real(dp),ALLOCATABLE :: orbit_mod_occ(:) ! Electronic config.: occupation of the modified orbital
  LOGICAL,ALLOCATABLE :: orbit_iscore(:)  ! Electronic configuration: TRUE for the core orbitals
  INTEGER :: norbit_val            ! Electronic configuration: number of valence orbitals 
  INTEGER,ALLOCATABLE :: orbit_val_n(:)   ! Electronic config.: n number of the valence orbital
  INTEGER,ALLOCATABLE :: orbit_val_l(:)   ! Electronic config.: l number of the valence orbital
  INTEGER,ALLOCATABLE :: orbit_val_k(:)   ! Electronic config.: kappa number of the valence orbital
  INTEGER :: lmax=-1               ! PAW Basis: maximum l value
  ! Cutoff radii
  real(dp) :: rc=0._dp               ! PAW basis: cut-off radius for the augmentation regions
  real(dp) :: rc_shap=0._dp          ! PAW basis: cut-off radius of the compensation charge shape function
  real(dp) :: rc_vloc=0._dp          ! PAW basis: matching radius for the local potential
  real(dp) :: rc_core=0._dp          ! PAW basis: matching radius for the pseudo-core density
  ! Additional basis functions
  INTEGER :: nbasis                ! PAW basis : number of basis functions
  INTEGER :: nbasis_add            ! PAW basis: number of additional basis functions (unbound states)
  INTEGER,ALLOCATABLE :: basis_add_l(:)      ! PAW basis: l number for the additional basis func.
  INTEGER,ALLOCATABLE :: basis_add_k(:)      ! PAW basis: kappa number for the additional basis func.
  real(dp),ALLOCATABLE :: basis_add_energy(:) ! PAW basis: ref. energy for the additional basis func.
  real(dp),ALLOCATABLE :: basis_func_rc(:)    ! PAW basis: rcut for the additional basis func.
  ! Projectors 
  INTEGER :: projector_type        ! Type of projectors (Bloechl, Vanderbilt,...)
  INTEGER :: pseudo_type           ! Type of pseudization scheme (Bessel,polynom, ...)
  INTEGER :: ortho_type            ! Type of orthogonalization scheme(Gram-Schmidt, ...)
  INTEGER :: pseudo_polynom2_pdeg  ! Polynom2 projectors: degree of the polynom
  real(dp) :: pseudo_polynom2_qcut  ! Polynom2 projectors: q-value for Fourier filtering
  INTEGER :: shapefunc_type           ! Compensation shape function type (sinc2, gaussian, ...)
  real(dp) :: shapefunc_gaussian_param ! Compensation shape function: parameter for gaussian type
  real(dp) :: hf_coretol            ! Tolerance for core density (Hartree-Fock only) 
  LOGICAL :: shapetcore            ! Flag activating building of tcore cancelling a negative compensation charge
  ! Local Psp
  INTEGER :: vloc_type             ! Type of local potential pseudization
  INTEGER :: vloc_l                ! Local potential: l quantum number (MTrouillier, Ultrasoft)
  real(dp) :: vloc_ene              ! Local potential: reference energy (MTrouillier, Ultrasoft)
  real(dp) :: vloc_setvloc_coef     ! "SetVloc" local potential: coefficient
  real(dp) :: vloc_setvloc_rad      ! "SetVloc" local potential: radius
  INTEGER :: vloc_kerker_power(4)  ! "Kerker" locazl potential: polynomial powers
  ! LIBPAW specific
  logical :: frozencorecalculation
  logical :: frozenvalecalculation
  logical :: setupfrozencore
  logical :: gaussianshapefunction
  logical :: besselshapefunction
  logical :: ColleSalvetti
  integer :: itype
  integer :: ixc
  integer :: xclevel
  integer :: npsc,nppc,npdc,npfc,npgc
  integer :: vhtnzc_mode, tpaw_mode, elin_mode
  real(dp) :: electrons
  real(dp) :: pot_ref
  real(dp), allocatable :: pot_refo(:)
  type(GridInfo) :: Grid
  TYPE(OrbitInfo) :: Orbit
  TYPE(PotentialInfo) :: Pot
  TYPE(SCFInfo) :: SCF
  TYPE(Pseudoinfo) :: PAW
  type(splinesolvinfo) :: spline
  logical :: has_to_print
 end type atompaw_type
!!***





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PUBLIC SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!***
 public :: atompaw_solve         ! Solve the atomic problem at fixed valence
 public :: atompaw_init          ! Init the atompaw solver
 public :: atompaw_destroy       ! Destroy the atomic solver
!!***


CONTAINS !===========================================================
!!***


!!****f* m_paw_atom_solve/atompaw_solve
!! NAME
!! atompaw_solve
!!
!! FUNCTION
!! Solve the atomic problem (at fixed valence) and update relevant quantities
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
!! Inspired from atompaw program

subroutine atompaw_solve(atp,pawrad,pawtab,&
  & nval,tnval,mqgrid_vl,qgrid_vl,epsatm,vlspl,&
& zion,update_paw,update_tnc,atm)

! TODO : vtau,dirac
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: mqgrid_vl
 logical,intent(in) :: update_paw,update_tnc
 real(dp), intent(inout) :: epsatm
 real(dp), intent(inout) :: zion
 type(atompaw_type), intent(inout) :: atp
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(inout) :: pawtab
 type(atomorb_type), intent(inout) :: atm
 !arrays
 real(dp), intent(in) :: nval(:),tnval(:)
 real(dp), intent(in) :: qgrid_vl(mqgrid_vl)
 real(dp),intent(out) :: vlspl(mqgrid_vl,2)
!Local variables-------------------------------
!scalars
 integer :: io,ir,icor,io2
 logical :: success,paw_proj
 character(len=500) :: msg
 real(dp) :: insph,norm,potshift,ovl
 real(dp) :: yp1,ypn,ekin,delta_zcore,int_pot
!arrays
 real(dp),ALLOCATABLE:: coredens(:),tcoredens(:),vhnzc_tmp(:)
 real(dp),ALLOCATABLE:: ff(:)
 type(pawrad_type) :: radmesh,vloc_mesh
 real(dp) , allocatable :: nhatc(:),vhatc(:)
 
! *************************************************************************

 ! Put valence density and core occupations to atompaw objects
 LIBPAW_ALLOCATE(coredens,(size(atp%Orbit%coreden)))
 atp%Orbit%valeden=zero
 atp%PAW%den=zero
 atp%PAW%tden=zero
 do io=1,pawtab%mesh_size
   atp%PAW%den(io)=nval(io)
   atp%Orbit%valeden(io)=nval(io)
   atp%PAW%tden(io)=tnval(io)
 enddo
 if(any(atp%PAW%tden<zero)) then
   write(std_out,*) 'atompaw_solve : tden is negative for some r ! Smoothening den instead'
   call smoothpower(atp%Grid,2,atp%PAW%den,atp%PAW%tden,atp%PAW)
 endif
 atp%Orbit%den=atp%Orbit%valeden+atp%Orbit%coreden
 coredens=atp%Orbit%coreden
 delta_zcore=zero
 icor=0
 do io=1,atp%Orbit%norbit
   if(atp%Orbit%iscore(io)) then
     icor=icor+1
     delta_zcore=delta_zcore-atp%Orbit%occ(io)
     atp%Orbit%occ(io)=atm%occ(icor,1)
     delta_zcore=delta_zcore+atp%Orbit%occ(io)
   endif
 enddo
 call pawrad_init(radmesh,atp%Grid%n,pawrad%mesh_type,pawrad%rstep,pawrad%lstep)
 LIBPAW_ALLOCATE(ff,(atp%grid%n))

! ! Solve atomic problem
 if(.not.atm%nc_conv) then
   call SCFatom(atp,.false.)
 endif

 ! Compute new core density
 write(msg,'(a)') 'atompaw_solve: orbital,%out sphere'
 call wrtout(std_out,msg,'COLL') 
 atp%Orbit%coreden=zero
 icor=0
 do io=1,atp%Orbit%norbit
   if(atp%Orbit%iscore(io)) then
     icor=icor+1
     norm=overlap(atp%Grid,atp%Orbit%wfn(1:atp%Grid%n,io),atp%Orbit%wfn(1:atp%Grid%n,io),1,atp%Grid%n)
     insph=overlap(atp%Grid,atp%Orbit%wfn(1:atp%PAW%irc,io),atp%Orbit%wfn(1:atp%PAW%irc,io),1,atp%PAW%irc)
     write(msg,*)io,(one -insph/norm)*100.0_dp
     call wrtout(std_out,msg,'COLL')
     atp%Orbit%coreden=atp%Orbit%coreden+atp%Orbit%occ(io)*(atp%Orbit%wfn(:,io))**2
   endif
 enddo

 ! Compute residue
 if(.not.atm%nc_conv) then
   ff(1:atp%Grid%n)=(atp%Orbit%coreden(1:atp%Grid%n)-coredens(1:atp%Grid%n))**2
   atm%nresid_c=sqrt(integrator(atp%Grid,ff)/atp%grid%r(atp%grid%n))
   ff(1:atp%Grid%n)=atp%Orbit%coreden(1:atp%Grid%n)
   if(integrator(atp%Grid,ff)/=zero) then
     atm%nresid_c=atm%nresid_c/integrator(atp%Grid,ff)
   else
     atm%nresid_c=one
   endif
   write(msg,*) 'atompaw_solve: nc residue',atm%nresid_c
   call wrtout(std_out,msg,'COLL')
 endif

 ! Update core dens
 do ir=2,atp%Grid%n
    coredens(ir)=atp%Orbit%coreden(ir)/(four*pi*atp%Grid%r(ir)**2)
 enddo
 call extrapolate(coredens)
 do ir=1,size(pawtab%coredens)
   pawtab%coredens(ir)=coredens(ir)
 enddo

 ! Update vhnzc
 LIBPAW_ALLOCATE(vhnzc_tmp,(radmesh%mesh_size))
 vhnzc_tmp=zero
 call atompaw_vhnzc(coredens,radmesh,vhnzc_tmp,atm%znucl) 
 do ir=1,size(pawtab%vhnzc)
   pawtab%vhnzc(ir)=vhnzc_tmp(ir)
 enddo
 LIBPAW_DEALLOCATE(vhnzc_tmp)

 ! Update edcc
 ff=zero
 ff(2:atp%Grid%n)=atp%Pot%rv(2:atp%Grid%n)*coredens(2:atp%Grid%n)*atp%Grid%r(2:atp%Grid%n)*four_pi
 CALL extrapolate(ff)
 atm%edcc=integrator(atp%Grid,ff)/two
 LIBPAW_DEALLOCATE(ff)

 ! Update ehnzc
 call atompaw_ehnzc(coredens,radmesh,atm%ehnzc,atm%znucl)

 ! update core wfs, kinetic energy and eigen energies
 atm%ekinc=zero
 atm%eeigc=zero
 icor=0
 do io=1,atp%Orbit%norbit
   if(atp%Orbit%iscore(io)) then
     icor=icor+1
     if(.not.atm%nc_conv) then
       atm%eig(icor,:)=atp%Orbit%eig(io)*half
     endif
     do ir=1,pawtab%mesh_size
       atm%phi(ir,icor,1)=atp%Orbit%wfn(ir,io)
     enddo
!     CALL
!     altkinetic(atp%Grid,atp%Orbit%wfn(:,io),atp%Orbit%eig(io),atp%Pot%rv,x)
     CALL kinetic(atp%Grid,atp%Orbit%wfn(:,io),atp%Orbit%l(io),ekin)
     atm%ekinc=atm%ekinc+ekin/two*atp%Orbit%occ(io)
     atm%eeigc=atm%eeigc+atp%Orbit%eig(io)*half*atp%Orbit%occ(io)
   endif
 enddo

 ! Update tnc
 if(update_tnc) then 
   ! Compute new tnc
   call setcoretail(atp%Grid,atp%Orbit%coreden,atp%PAW,atp%needvtau)
   LIBPAW_ALLOCATE(tcoredens,(atp%Grid%n))
   do ir=2,atp%Grid%n
     tcoredens(ir)=atp%PAW%tcore(ir)/(four*pi*atp%Grid%r(ir)**2)
   enddo
   call extrapolate(tcoredens)
   ! Update tnc
   do ir=1,pawtab%mesh_size
     pawtab%tcoredens(ir,1)=tcoredens(ir)
   enddo
   call pawpsp_cg(pawtab%dncdq0,pawtab%d2ncdq0,mqgrid_vl,qgrid_vl,pawtab%tcorespl(:,1),radmesh,tcoredens,yp1,ypn)
   call paw_spline(qgrid_vl,pawtab%tcorespl(:,1),mqgrid_vl,yp1,ypn,pawtab%tcorespl(:,2))
   LIBPAW_DEALLOCATE(tcoredens)
 endif

 ! Update PAW stuff
 if(update_paw.or.atp%vhtnzc_mode==2) then
   ! Compute potential shift
   if(atp%elin_mode==1) then
     call simp_gen(int_pot,atp%Pot%rv(1:pawtab%mesh_size)*pawrad%rad(1:pawtab%mesh_size)*pawtab%shapefunc(1:pawtab%mesh_size,1),pawrad)
     potshift=int_pot-atp%pot_ref
   else
     LIBPAW_ERROR('Not ready yet')
     potshift=atm%eigshift
   endif
   do io=1,atp%Grid%n
     atp%Pot%rvh(io)=atp%Pot%rvh(io)-potshift*atp%Grid%r(io)
     atp%Pot%rv(io)=atp%Pot%rv(io)-potshift*atp%Grid%r(io)
   enddo
   atp%Pot%v0=atp%Pot%v0-potshift
   call setbasis(atp%Grid,atp%Pot,atp%Orbit,atp%PAW,atp,potshift)  
   call SetPAWOptions2(atp,success)
   ! Update PAW transform
   if(update_paw) then
     if(atp%tpaw_mode>1) then
       do io=1,pawtab%basis_size
         paw_proj=.true.
         do io2=1,atp%Orbit%norbit
           if(atp%PAW%valencemap(io2)==io) then
             if(atp%orbit%issemicore(io2).or.atp%tpaw_mode==3) then
               paw_proj=.false.
             endif
           endif 
         enddo
         if(paw_proj) then
           atp%PAW%ophi(:,io)=zero
           atp%PAW%otphi(:,io)=zero
           do ir=1,pawtab%mesh_size
             atp%PAW%ophi(ir,io)=pawtab%phi(ir,io)
             atp%PAW%otphi(ir,io)=pawtab%tphi(ir,io)
           enddo
           do io2=1,atp%Orbit%norbit
             if(atp%Orbit%iscore(io2).and.atp%orbit%l(io2)==atp%PAW%l(io)) then
               ovl=overlap(atp%grid,atp%PAW%ophi(1:atp%PAW%irc-5,io),atp%orbit%wfn(1:atp%PAW%irc-5,io2),1,atp%PAW%irc-5)
               atp%PAW%ophi(1:atp%PAW%irc-5,io)=atp%PAW%ophi(1:atp%PAW%irc-5,io)-&
&                                               ovl*atp%orbit%wfn(1:atp%PAW%irc-5,io2)
             endif
           enddo
         endif
       enddo
     endif
     do io=1,pawtab%basis_size
       do ir=1,pawtab%mesh_size
         pawtab%phi(ir,io)=atp%PAW%ophi(ir,io)
         pawtab%tphi(ir,io)=atp%PAW%otphi(ir,io)
      enddo
     enddo
     ! Update Kij
     if(.not.allocated(pawtab%kij)) then
       LIBPAW_ALLOCATE(pawtab%kij,(pawtab%lmn2_size))
     endif
     call calc_kij(atp%PAW,atp%Grid,pawtab%kij,pawtab,&
&     atp%scalarrelativistic,atp%needvtau)
   endif
 endif

 ! Update vhtnzc, zion and epsatm
 if(abs(delta_zcore)<tol15*atm%zcore) atm%zcore_conv=.true.
 if(.not.atm%zcore_conv) then
   zion=atm%znucl-atm%zcore
   delta_zcore=atm%zcore-atm%zcore_orig
   write(msg,*) 'atompaw_solve: delta_zcore',delta_zcore
   call wrtout(std_out,msg,'COLL')
   call pawrad_init(vloc_mesh,mesh_size=size(pawtab%vhtnzc),mesh_type=pawrad%mesh_type,&
&   rstep=pawrad%rstep,lstep=pawrad%lstep)
   if(atp%vhtnzc_mode==2) then
     call FindVlocfromVeff(atp%Grid,atp%PAW,atp,potshift)
     if(pawtab%usexcnhat==1) then
       pawtab%vhtnzc(1:size(pawtab%vhtnzc))=half*atp%PAW%abinitvloc(1:size(pawtab%vhtnzc))
     else
       pawtab%vhtnzc(1:size(pawtab%vhtnzc))=half*atp%PAW%abinitnohat(1:size(pawtab%vhtnzc))
     endif
   elseif(atp%vhtnzc_mode==1) then
     ! Compute nhatc=shapefunction*delta_zcore
     LIBPAW_ALLOCATE(nhatc,(size(pawtab%vhtnzc)))
     nhatc=zero
     do ir=1,size(pawtab%shapefunc(:,1))
       nhatc(ir)=delta_zcore*pawtab%shapefunc(ir,1)*vloc_mesh%rad(ir)**2
     enddo
     LIBPAW_ALLOCATE(vhatc,(size(pawtab%vhtnzc)))
     call poisson(nhatc,0,vloc_mesh,vhatc)
     do ir=2,vloc_mesh%mesh_size
       vhatc(ir)=vhatc(ir)/vloc_mesh%rad(ir)
     enddo
     call pawrad_deducer0(vhatc,vloc_mesh%mesh_size,vloc_mesh)
     LIBPAW_DEALLOCATE(nhatc)
     ! Add it to original vhtnzc
     pawtab%vhtnzc=atm%vhtnzc_orig+vhatc
     LIBPAW_DEALLOCATE(vhatc)
   endif
   call pawpsp_lo(epsatm,mqgrid_vl,qgrid_vl,vlspl(:,1),&
&                     vloc_mesh,pawtab%vhtnzc,yp1,ypn,&
&                     zion)
   call  paw_spline(qgrid_vl,vlspl(:,1),mqgrid_vl,yp1,ypn,vlspl(:,2))
   write(msg,*) 'atompaw_solve: nc epsatm',epsatm
   call wrtout(std_out,msg,'COLL')
   call pawrad_free(vloc_mesh)
 endif

 ! update dij0
 call atompaw_dij0(pawtab%indlmn,pawtab%kij,pawtab%lmn_size,coredens,0,pawtab,pawrad,radmesh,&
&                      pawrad,pawtab%vhtnzc,atp%Pot%zz)
 LIBPAW_DEALLOCATE(coredens)

 ! update tcoretau : TODO : tau
 ! update kinetic part of dij0 : TODO : positron
 ! Clean up
 call pawrad_free(radmesh) 
 
end subroutine atompaw_solve
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atom_solve/atompaw_init
!! NAME
!! atompaw_init
!!
!! FUNCTION
!! Initialize an atompaw type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
!! Inspired from SCFatom_init in atompaw

subroutine atompaw_init(pawtab,pawrad,atp,znucl,atm,sctol,elin_mode,vhtnzc_mode,tpaw_mode)
 ! TODO : BDsolve
 implicit none
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: znucl,elin_mode,vhtnzc_mode,tpaw_mode
 real(dp), intent(in) :: sctol
 type(pawtab_type), intent(inout) :: pawtab
 type(pawrad_type), intent(in) :: pawrad
 type(atompaw_type), intent(inout) :: atp
 type(atomorb_type), intent(inout) :: atm
!arrays
!Local variables-------------------------------
!scalars
 CHARACTER(len=500) :: input_file
 character(len=500) :: msg
 REAL(dp)    :: a1,a2,a3,hval,r0
 INTEGER :: ii,jj,icor,ir,io,fnln
 logical :: fmt_xml
 real(dp) :: ekin,insph,norm
 type(pawrad_type) :: radmesh
!arrays
 real(dp), allocatable :: ff(:)

! *************************************************************************

 ! File to read
 input_file=trim(atm%fname)
 fnln=len(trim(atm%fname))
 fmt_xml=.false.
 if (fnln>3) then
    fmt_xml=(input_file(fnln-3:fnln)=='.xml')
 endif
 if(.not.fmt_xml) then
   LIBPAW_ERROR('RCPAW ONLY COMPATIBLE WITH XML COREWF FILE')
 endif

! Initialize global constants
 machine_precision = zero
 a1 = 4._dp/3._dp
 DO WHILE (machine_precision == 0._dp)
   a2 = a1 - 1._dp
   a3 = a2 + a2 + a2
   machine_precision = ABS(a3 - 1._dp)
 ENDDO
 machine_zero= machine_precision**5
 machine_infinity = 1._dp/machine_zero
 practical_zero=machine_precision**2
 minlogarg=machine_precision; minlog=LOG(minlogarg)
 maxlogarg=1._dp/machine_precision; maxlog=LOG(maxlogarg)
 minexparg=LOG(machine_precision);  minexp=0._dp
 maxexparg=-LOG(machine_precision);  maxexp=EXP(maxexparg)
 has_to_print=.false.

 ! Read atp input and initialize
 atp%elin_mode=elin_mode
 atp%vhtnzc_mode=vhtnzc_mode
 atp%tpaw_mode=tpaw_mode
 call input_dataset_read(atp,input_file)
 atp%npsc=0
 atp%nppc=1
 atp%npdc=2
 atp%npfc=3
 atp%npgc=4
 atp%frozenvalecalculation=.false.
 atp%frozencorecalculation=.false.
 atp%setupfrozencore=.false.
 atp%gaussianshapefunction=.false.
 atp%besselshapefunction=.false.
 atp%collesalvetti=.false.
 IF (TRIM(atp%gridkey)=='LINEAR') THEN
   hval=atp%gridmatch/(atp%gridpoints-1)
   CALL InitGrid(atp%Grid,hval,atp%gridrange)
   ELSEIF (TRIM(atp%gridkey)=='LOGGRID') THEN
     hval=logh
     CALL findh(real(atp%atomic_charge,kind=8),atp%gridmatch,atp%gridpoints,hval,r0)
     CALL InitGrid(atp%Grid,pawrad%lstep,atp%gridrange,r0=pawrad%rstep)
  ELSEIF (TRIM(atp%gridkey)=='LOGGRID4') THEN
    hval=logh
    CALL findh_given_r0(real(atp%atomic_charge,kind=8),atp%gridmatch,lor00,&
&                       atp%gridpoints,hval)
    CALL InitGrid(atp%Grid,hval,atp%gridrange,r0=lor00/atp%atomic_charge)
 ENDIF
 call print_check_atompaw_params(atp) 

 ! Init potentials
 CALL InitPot(atp%Pot,atp%Grid%n)
 atp%Pot%sym=atp%atomic_symbol
 atp%Pot%zz=0._dp
 atp%Pot%q=0._dp;
 atp%Pot%v0=0._dp
 atp%Pot%v0p=0._dp
 atp%Pot%Nv0=0
 atp%Pot%Nv0p=0
 atp%Pot%nz=znucl 
 atp%Pot%zz=atp%Pot%nz
 atp%Pot%needvtau=atp%needvtau
 atp%Pot%finitenucleus=atp%finitenucleus
 atp%Pot%finitenucleusmodel=atp%finitenucleusmodel
 CALL Get_Nuclearpotential(atp%Grid,atp%Pot)

 ! Init XC
! atp%temp=temp
 call initexch(atp)

 ! Init orbitals
 ii=maxval(atp%np)
 jj=atp%np(1)
 IF(atp%np(2)>0) jj=jj+atp%np(2)-1
 IF(atp%np(3)>0) jj=jj+atp%np(3)-2
 IF(atp%np(4)>0) jj=jj+atp%np(4)-3
 IF(atp%np(5)>0) jj=jj+atp%np(5)-4
 If (atp%diracrelativistic) jj=jj+jj   !  need more orbitals
 CALL InitOrbit(atp%Orbit,jj,atp%Grid%n,atp%exctype,atp%diracrelativistic,atp%scalarrelativistic,&
&     atp%frozencorecalculation,atp%frozenvalecalculation)
 atp%Orbit%nps=atp%np(1);atp%Orbit%npp=atp%np(2);atp%Orbit%npd=atp%np(3)
 atp%Orbit%npf=atp%np(4);atp%Orbit%npg=atp%np(5)
 CALL Prepare_Orbit(atp,ii,jj)
 atp%Orbit%npsc=atp%npsc;atp%Orbit%nppc=atp%nppc;atp%Orbit%npdc=atp%npdc
 atp%Orbit%npfc=atp%npfc;atp%Orbit%npgc=atp%npgc
 atp%Pot%q=atp%electrons

 ! Init SCF
 CALL InitSCF(atp%SCF)
 IF (atp%needvtau) then
   atp%usespline=.true.
 endif
 if(atp%usespline) CALL initsplinesolver(atp%Grid,atp%splns,atp%splr0,atp%needvtau,atp%spline)

 ! Init PAW
 atp%PAW%irc=FindGridIndex(atp%Grid,atp%rc)
 atp%PAW%irc_shap=FindGridIndex(atp%Grid,atp%rc_shap)
 atp%PAW%irc_vloc=FindGridIndex(atp%Grid,atp%rc_vloc)
 atp%PAW%irc_core=FindGridIndex(atp%Grid,atp%rc_core)
 atp%PAW%rc=atp%Grid%r(atp%PAW%irc)
 atp%PAW%rc_shap=atp%Grid%r(atp%PAW%irc_shap)
 atp%PAW%rc_vloc=atp%Grid%r(atp%PAW%irc_vloc)
 atp%PAW%rc_core=atp%Grid%r(atp%PAW%irc_core)
 atp%PAW%lmax=atp%lmax
 call InitPAW(atp%PAW,atp%Grid,atp%Orbit)

 ! Re-solve (temporary, this should be added to atompaw)
 call SCFatom(atp,.true.)
 call simp_gen(atp%pot_ref,atp%Pot%rv(1:pawtab%mesh_size)*pawrad%rad(1:pawtab%mesh_size)*pawtab%shapefunc(1:pawtab%mesh_size,1),pawrad)
 write(msg,'(a)') 'atompaw_init: orbital, %out of sphere, core, semicore' 
 call wrtout(std_out,msg,'COLL')
 do io=1,atp%Orbit%norbit
     norm=overlap(atp%Grid,atp%Orbit%wfn(1:atp%Grid%n,io),atp%Orbit%wfn(1:atp%Grid%n,io),1,atp%Grid%n)
     insph=overlap(atp%Grid,atp%Orbit%wfn(1:atp%PAW%irc,io),atp%Orbit%wfn(1:atp%PAW%irc,io),1,atp%PAW%irc)
     if(.not.atp%Orbit%iscore(io).and.(one-insph/norm)*100.0_dp<sctol) atp%Orbit%issemicore(io)=.true.
     write(msg,*)io,(one-insph/norm)*100.0_dp,atp%Orbit%iscore(io),atp%Orbit%issemicore(io)
     call wrtout(std_out,msg,'COLL')
 enddo
 !!!

 LIBPAW_ALLOCATE(atp%pot_refo,(pawtab%basis_size))
 LIBPAW_ALLOCATE(ff,(pawtab%mesh_size))
 ff(2:pawtab%mesh_size)=atp%Pot%rv(2:pawtab%mesh_size)/atp%grid%r(2:pawtab%mesh_size)
 call extrapolate(ff) 
 do io=1,pawtab%basis_size
   call simp_gen(atp%pot_refo(io),ff(1:pawtab%mesh_size)*(pawtab%phi(1:pawtab%mesh_size,io))**2,pawrad)
 enddo
 LIBPAW_DEALLOCATE(ff)

 ! Densities
 atp%Orbit%valeden=zero
 atp%Orbit%coreden=zero
 do io=1,atp%Orbit%norbit
   if(atp%Orbit%iscore(io)) then
     atp%Orbit%coreden=atp%Orbit%coreden+atp%Orbit%occ(io)*(atp%Orbit%wfn(:,io))**2
   else
     atp%Orbit%valeden(:)=atp%Orbit%valeden(:)+atp%Orbit%occ(io)*(atp%Orbit%wfn(:,io))**2
   endif
 enddo
 ! Core energies
 atm%edcc=zero
 atm%ekinc=zero
 atm%eeigc=zero
 icor=0
 do io=1,atp%Orbit%norbit
   if(atp%Orbit%iscore(io)) then
     icor=icor+1
     CALL kinetic(atp%Grid,atp%Orbit%wfn(:,io),atp%Orbit%l(io),ekin)
     if(abs(atp%Orbit%eig(io)/two-atm%eig(icor,1))>tol1) then
       LIBPAW_ERROR('Inconsistent RCPAW core files')
     endif
     if(abs(atp%Orbit%occ(io)-atm%occ(icor,1))>tol1) then
       LIBPAW_ERROR('Inconsistent RCPAW core files')
     endif
     atm%ekinc=atm%ekinc+ekin/two*atp%Orbit%occ(io)
     atm%eeigc=atm%eeigc+atp%Orbit%eig(io)*atp%Orbit%occ(io)/two
   endif
 enddo
 LIBPAW_ALLOCATE(ff,(atp%grid%n))
 ff=zero
 ff(2:atp%Grid%n)=atp%Pot%rv(2:atp%Grid%n)*&
&                   atp%Orbit%coreden(2:atp%Grid%n)/atp%Grid%r(2:atp%Grid%n)
 CALL extrapolate(ff)
 atm%edcc=integrator(atp%Grid,ff)/two
 ff(2:atp%Grid%n)=atp%Orbit%coreden(2:atp%Grid%n)/(four*pi*atp%Grid%r(2:atp%Grid%n)**2)
 CALL extrapolate(ff)
 call pawrad_init(radmesh,atp%Grid%n,pawrad%mesh_type,pawrad%rstep,pawrad%lstep)
 call atompaw_ehnzc(ff,radmesh,atm%ehnzc,atm%znucl)
 LIBPAW_DEALLOCATE(ff)
 atm%min_eigv=half*minval(atp%Orbit%eig,mask=.not.atp%Orbit%iscore)
 call pawrad_free(radmesh)

 ! Prepare for frozen val calculation
 atp%frozenvalecalculation=.true.
 atp%Orbit%frozenvalecalculation=.true.
 atp%projector_type=PROJECTOR_TYPE_MARSMAN
 atp%PAW%otp=zero
 do io=1,pawtab%basis_size
   do ir=1,atp%PAW%irc
     atp%PAW%otp(ir,io)=pawtab%tproj(ir,io)
   enddo
 enddo

end subroutine atompaw_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atom_solve/atompaw_destroy
!! NAME
!! atompaw_destroy
!!
!! FUNCTION
!! Destroy an atompaw type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine atompaw_destroy(atp)

 implicit none
!Arguments ------------------------------------
!scalars
 type(atompaw_type), intent(inout) :: atp

! *********************************************************************

 if(allocated(atp%orbit_mod_n)) then 
   LIBPAW_DEALLOCATE(atp%orbit_mod_n)
 endif
 if(allocated(atp%orbit_mod_l)) then
   LIBPAW_DEALLOCATE(atp%orbit_mod_l)
 endif
 if(allocated(atp%orbit_mod_k)) then
   LIBPAW_DEALLOCATE(atp%orbit_mod_k)
 endif
 if(allocated(atp%orbit_mod_occ)) then
   LIBPAW_DEALLOCATE(atp%orbit_mod_occ)
 endif
 if(allocated(atp%orbit_iscore)) then
   LIBPAW_DEALLOCATE(atp%orbit_iscore)
 endif
 if(allocated(atp%orbit_val_n)) then
   LIBPAW_DEALLOCATE(atp%orbit_val_n)
 endif
 if(allocated(atp%orbit_val_l)) then
   LIBPAW_DEALLOCATE(atp%orbit_val_l)
 endif
 if(allocated(atp%orbit_val_k)) then
   LIBPAW_DEALLOCATE(atp%orbit_val_k)
 endif
 if(allocated(atp%basis_add_l)) then
   LIBPAW_DEALLOCATE(atp%basis_add_l)
 endif
 if(allocated(atp%basis_add_k)) then
   LIBPAW_DEALLOCATE(atp%basis_add_k)
 endif
 if(allocated(atp%basis_add_energy)) then
   LIBPAW_DEALLOCATE(atp%basis_add_energy)
 endif
 if(allocated(atp%basis_func_rc)) then
   LIBPAW_DEALLOCATE(atp%basis_func_rc)
 endif
 if(allocated(atp%pot_refo)) then
   LIBPAW_DEALLOCATE(atp%pot_refo)
 endif
 call DestroyGrid(atp%Grid)
 call DestroyOrbit(atp%Orbit)
 call DestroyPot(atp%Pot)
 call DestroyPAW(atp%PAW)
 call deallocatesplinesolver(atp%spline,atp%needvtau)

end subroutine atompaw_destroy
!!***





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BELOW THIS ARE PRIVATE ROUTINES MOSTLY IMPORTED FROM ATOMPAW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 0. Developped for RCPAW 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  makebasis_marsman
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine makebasis_marsman(atp)
 type(atompaw_type),intent(inout) :: atp
 integer :: ns,np,nd,nf,ng,io,l
 real(dp):: v0,v0p
 integer, allocatable :: map(:)
 LIBPAW_ALLOCATE(map,(atp%PAW%nbase))
 ns=0
 np=0
 nd=0
 nf=0
 ng=0
 map=0
 call zeropot(atp%Grid,atp%PAW%rveff,v0,v0p)
 do io=1,atp%PAW%nbase
   l=atp%PAW%l(io)
   if(l==0) then
     ns=ns+1
     map(io)=ns
   endif
   if(l==1) then
     np=np+1
     map(io)=np
   endif
   if(l==2) then
     nd=nd+1
     map(io)=nd
   endif
   if(l==3) then
      nf=nf+1
      map(io)=nf
   endif
   if(l==4) then
     ng=ng+1
     map(io)=ng
   endif
 enddo
 if(ns>0) call marsman_tphi(atp,map,0,ns)
 if(np>0) call marsman_tphi(atp,map,1,np)
 if(nd>0) call marsman_tphi(atp,map,2,nd)
 if(nf>0) call marsman_tphi(atp,map,3,nf)
 if(ng>0) call marsman_tphi(atp,map,4,ng)
 ! renormalize phis
 do io=1,atp%PAW%nbase
   atp%PAW%otphi(:,io)=atp%PAW%tphi(:,io)
   atp%PAW%ophi(:,io)=atp%PAW%phi(:,io)*atp%PAW%otphi(atp%PAW%irc,io)/atp%PAW%phi(atp%PAW%irc,io)
   atp%PAW%Kop(1,io)=zero
   atp%PAW%Kop(2:atp%Grid%n,io)=(atp%PAW%eig(io)-atp%Pot%rv(2:atp%Grid%n)/&
&                                atp%Grid%r(2:atp%Grid%n))*atp%PAW%ophi(2:atp%Grid%n,io)
 enddo
 LIBPAW_DEALLOCATE(map)
end subroutine makebasis_marsman



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  marsman_tphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine marsman_tphi(atp,map,l_in,n)
 integer, intent(in) :: l_in
 integer, intent(in) :: n 
 type(atompaw_type),intent(inout) :: atp
 integer, intent(in) :: map(atp%PAW%nbase)
 integer :: io,nodes,ii,match,irc
 integer :: io2,l2,jj,ioj,kk,l,ir,ir1
 real(dp):: energy,v0,v0p,zeroval
 real(dp), allocatable :: ksi_i0(:,:),ksi_ij(:,:,:)
 real(dp), allocatable :: ksi_i0_0(:),ksi_ij_0(:,:)
 REAL(dp), ALLOCATABLE :: p1(:),p2(:)
 real(dp), allocatable :: AA(:,:), BB(:)
 real(dp) :: dp1(atp%PAW%irc),dp2(atp%PAW%irc)
 real(dp) :: x1,y1,y2,a,b,c,r1,r2
 logical :: need_fit
 irc=atp%PAW%irc
 CALL zeropot(atp%Grid,atp%PAW%rveff,v0,v0p)
 LIBPAW_ALLOCATE(p1,(atp%Grid%n))
 LIBPAW_ALLOCATE(p2,(atp%Grid%n))
 LIBPAW_ALLOCATE(ksi_i0,(atp%Grid%n,n))
 LIBPAW_ALLOCATE(ksi_ij,(atp%Grid%n,n,n))
 LIBPAW_ALLOCATE(ksi_i0_0,(n))
 LIBPAW_ALLOCATE(ksi_ij_0,(n,n))
 ksi_i0=zero
 ksi_ij=zero
 ii=0
 do io=1,atp%PAW%nbase
   l=atp%PAW%l(io)
   if(l==l_in) then
     ii=ii+1
     energy=atp%PAW%eig(io)
     CALL ClassicalTurningPoint(atp%Grid,atp%PAW%rveff,l,energy,match)
     if(match>irc) match=2
  !   !! backward Not used anymore
  !   p2=zero
  !   p2(irc:atp%Grid%n)=atp%PAW%phi(irc:atp%Grid%n,io)
  !   CALL backward_numerov(atp%Grid,l,2,energy,atp%PAW%rveff,p2,nend=irc+1)
     !! forward 
     p1=zero
     p1(2)=wfninit(-0.5_dp*atp%PAW%rveff(1),l,v0,v0p,energy,atp%Grid%r(2))
     zeroval=zero
     IF (l==0) zeroval=atp%PAW%rveff(1)
     IF (l==1) zeroval=2
     CALL forward_numerov(atp%Grid,l,irc+5,energy,atp%PAW%rveff,zeroval,p1,nodes)!,p3val=p1(3))
     ksi_i0(:,ii)=p1(:)
     dp1=zero
     call derivative(atp%grid,p1,dp1,1,irc)
!     p2(2:size(p2))=p2(2:size(p2))*atp%grid%r(2:size(p2))**l_in
!     call extrapolate(p2)
     ksi_i0_0(ii)=dp1(irc)
     jj=0
     do io2=1,atp%PAW%nbase
       l2=atp%PAW%l(io2)
       if(l2==l) then
         jj=jj+1
         match=2
!         p2=zero
!         p2(irc:atp%Grid%n)=atp%PAW%phi(irc:atp%Grid%n,io)
!         CALL backward_numerov(atp%Grid,l,match,energy,atp%PAW%rveff,p2,nend=irc+1,proj=atp%PAW%otp(:,io2))
         p1=zero
         p1(2)=wfninit(-0.5_dp*atp%PAW%rveff(1),l,v0,v0p,energy,atp%Grid%r(2))
         CALL forward_numerov(atp%Grid,l,irc+5,energy,atp%PAW%rveff,zeroval,p1,nodes,proj=atp%PAW%otp(:,io2))!,p3val=p1(3))
         ksi_ij(:,ii,jj)=p1(:)
!         p2(2:size(p2))=p2(2:size(p2))*atp%grid%r(2:size(p2))**l_in
!         call extrapolate(p2)
         dp1=zero
         call derivative(atp%grid,p1,dp1,1,irc)
         ksi_ij_0(ii,jj)=dp1(irc)
       endif
     enddo
   endif
 enddo
 LIBPAW_ALLOCATE(AA,(n+1,n+1))
 LIBPAW_ALLOCATE(BB,(n+1))
 do ioj=1,atp%PAW%nbase
   l=atp%PAW%l(ioj)
   if(l==l_in) then
     jj=map(ioj)
     AA=zero
     BB=zero
     ! get derivative of phi
     dp1=zero
     call derivative(atp%grid,atp%PAW%phi(:,ioj),dp1,1,irc)
     ii=0
     do io=1,atp%PAW%nbase
       l2=atp%PAW%l(io)
       if(l==l2) then
         ii=ii+1
         if(ii==jj) BB(ii)=one
         p2(:)=ksi_i0(:,jj)*atp%PAW%otp(:,io)
         AA(ii,n+1)=integrator(atp%Grid,p2)
         do kk=1,n
           p2(:)=ksi_ij(:,jj,kk)*atp%PAW%otp(:,io)
           AA(ii,kk)=integrator(atp%Grid,p2)
         enddo
       endif
     enddo
     do kk=1,n
       AA(n+1,kk)=ksi_ij_0(jj,kk)-dp1(irc)*ksi_ij(irc,jj,kk)/atp%PAW%phi(irc,ioj)
     enddo
     AA(n+1,n+1)=ksi_i0_0(jj)-dp1(irc)*ksi_i0(irc,jj)/atp%PAW%phi(irc,ioj)
     call SolveAXeqBM(n+1,AA,BB,n+1)
     atp%PAW%tphi(:,ioj)=zero
     do kk=1,n
       atp%PAW%tphi(:,ioj)=atp%PAW%tphi(:,ioj)+BB(kk)*ksi_ij(:,jj,kk)
     enddo
     atp%PAW%tphi(:,ioj)=atp%PAW%tphi(:,ioj)+BB(n+1)*ksi_i0(:,jj)
    ! In some cases, numerical error => need to fit tail of the orbital
    ! Check if low lying orbital
     if(atp%basis_func_rc(ioj)<atp%rc) then
       write(std_out,*) 'TT1',ioj,atp%PAW%eig(ioj)
       ! Check if tphi is close to 0 at irc
       if(abs(atp%PAW%tphi(irc,ioj))/maxval(abs(atp%PAW%tphi(:,ioj)))<tol2) then
          write(std_out,*) 'TT2'
          ! Check if tphi is non-decreasing at irc
          need_fit=.false.
          call derivative(atp%grid,abs(atp%PAW%tphi(:,ioj)),dp2,1,irc)
          ir1=FindGridIndex(atp%Grid,atp%basis_func_rc(ioj))
          do ir=ir1,irc+5
            if(dp2(ir)>zero) then
              need_fit=.true.
              exit
            endif
          enddo
          if(need_fit) then
            x1=atp%PAW%tphi(ir1,ioj)
            call derivative(atp%grid,atp%PAW%tphi(:,ioj),dp2,1,ir1)
            y1=dp2(ir1)/x1
            call derivative(atp%grid,atp%PAW%phi(:,ioj),dp1,1,irc)
            y2=dp1(irc)/atp%PAW%phi(irc,ioj)
            r1=atp%grid%r(ir1)
            write(std_out,*) 'R1',r1
            r2=atp%grid%r(irc)
            c=(y1-y2)/(one/r1-one/r2)
            b=c/r1-y1
            a=x1/(r1**c*exp(-b*r1)) 
            do ir=ir1,irc+5
              atp%PAW%tphi(ir,ioj)=a*atp%grid%r(ir)**c*exp(-b*atp%grid%r(ir))
            enddo
         endif
       endif
     endif
   endif
 enddo
 LIBPAW_DEALLOCATE(AA)
 LIBPAW_DEALLOCATE(BB)
 LIBPAW_DEALLOCATE(p1)
 LIBPAW_DEALLOCATE(p2)
 LIBPAW_DEALLOCATE(ksi_i0)
 LIBPAW_DEALLOCATE(ksi_ij)
 LIBPAW_DEALLOCATE(ksi_i0_0)
 LIBPAW_DEALLOCATE(ksi_ij_0)
end subroutine marsman_tphi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  calc_kij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_kij(PAW,grid,kij,pawtab,scalarrelativistic,needvtau)
 logical, intent(in) :: scalarrelativistic,needvtau
 type(gridinfo), intent(in) :: grid
 type(pseudoinfo), intent(in) :: PAW
 type(pawtab_type), intent(in) :: pawtab
 real(dp), intent(inout) :: kij(pawtab%lmn2_size)
 integer :: nbase,l,ib,jb
 integer :: jlmn,j0lmn,jlm,jln,ilmn,klmn,ilm,iln
 real(dp) :: x, y
 real(dp), allocatable :: kij_atp(:,:)
 nbase=PAW%nbase
 LIBPAW_ALLOCATE(kij_atp,(nbase,nbase))
 DO ib=1,nbase
   l=PAW%l(ib)
   DO jb=1,nbase
     IF (PAW%l(jb)==l) THEN
       If (scalarrelativistic) then
         call altdtij(Grid,PAW,ib,jb,x,needvtau)
       Else
         CALL deltakinetic_ij(Grid,PAW%ophi(:,ib),PAW%ophi(:,jb), &
&              PAW%otphi(:,ib),PAW%otphi(:,jb),l,x,PAW%irc)
       Endif
       if(has_to_print) WRITE(STD_OUT,'(" Kinetic ", 3i5, 1p,3e15.7)') ib,jb,l,x
       Kij_atp(ib,jb)=x
     ENDIF
   ENDDO
 ENDDO
 ! Average equivalent terms
 DO ib=1,nbase
   DO jb=ib,nbase
     IF(jb>ib) THEN
       x=Kij_atp(ib,jb); y=Kij_atp(jb,ib)
       x=0.5_dp*(x+y)
       Kij_atp(ib,jb)=x; Kij_atp(jb,ib)=x
     ENDIF
   ENDDO
 ENDDO
 ! convert from atompaw to libpaw
 kij=zero
 do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=pawtab%indlmn(4,jlmn);jln=pawtab%indlmn(5,jlmn)
   do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ilm=pawtab%indlmn(4,ilmn);iln=pawtab%indlmn(5,ilmn)
     if (ilm==jlm) kij(klmn)=half*Kij_atp(iln,jln) ! Conversion Ry->Ha
   enddo
 enddo
 LIBPAW_DEALLOCATE(kij_atp)
end subroutine calc_kij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Prepare_Orbit
!!   Inspired from SCFatom_Init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Prepare_Orbit(atp,ii,jj)
 integer, intent(in) :: ii,jj
 type(atompaw_type), intent(inout) :: atp
 INTEGER :: icount,id,is,ip,jf,ig,io,l,nfix,kappa,i
 real(dp) :: xocc
 INTEGER, ALLOCATABLE :: nl(:,:)
 icount=0
 if(.not.atp%diracrelativistic) then
   LIBPAW_ALLOCATE(nl,(ii,jj));nl=0
   IF (atp%Orbit%nps.GT.0) THEN
     DO is=1,atp%Orbit%nps
       icount=icount+1
       nl(is,1)=icount
       atp%Orbit%occ(icount)=2._dp
       atp%Orbit%np(icount)=is
       atp%Orbit%l(icount)=0
       if(atp%orbit_iscore(icount)) atp%npsc=is
     ENDDO
   ENDIF
   IF (atp%Orbit%npp.GT.1) THEN
     DO ip=2,atp%Orbit%npp
       icount=icount+1
       nl(ip,2)=icount
       atp%Orbit%occ(icount)=6._dp
       atp%Orbit%np(icount)=ip
       atp%Orbit%l(icount)=1
       if(atp%orbit_iscore(icount)) atp%nppc=ip
     ENDDO
   ENDIF
   IF (atp%Orbit%npd.GT.2) THEN
     DO id=3,atp%Orbit%npd
       icount=icount+1
       nl(id,3)=icount
       atp%Orbit%occ(icount)=10._dp
       atp%Orbit%np(icount)=id
       atp%Orbit%l(icount)=2
       if(atp%orbit_iscore(icount)) atp%npdc=id
     ENDDO
   ENDIF
   IF (atp%Orbit%npf.GT.3) THEN
     DO jf=4,atp%Orbit%npf
       icount=icount+1
       nl(jf,4)=icount
       atp%Orbit%occ(icount)=14._dp
       atp%Orbit%np(icount)=jf
       atp%Orbit%l(icount)=3
       if(atp%orbit_iscore(icount)) atp%npfc=jf
     ENDDO
   ENDIF
   IF(atp%Orbit%npg.GT.4) THEN
     DO ig=5,atp%Orbit%npg
       icount=icount+1
       nl(ig,5)=icount
       atp%Orbit%occ(icount)=18._dp
       atp%Orbit%np(icount)=ig
       atp%Orbit%l(icount)=4
       if(atp%orbit_iscore(icount)) atp%npgc=ig
     ENDDO
   ENDIF
   atp%Orbit%norbit=icount
   if(has_to_print) write(std_out,*)' Below are listed the default occupations '
   if(has_to_print) write(std_out,"(' n  l     occupancy')")
   DO io=1,atp%Orbit%norbit
       if(has_to_print) WRITE(std_out,'(i2,1x,i2,4x,1p,1e15.7)') &
&           atp%Orbit%np(io),atp%Orbit%l(io),atp%Orbit%occ(io)
   ENDDO
   atp%Orbit%iscore=atp%orbit_iscore
   DO io=1,atp%norbit_mod
     l=atp%orbit_mod_l(io)
     ip=atp%orbit_mod_n(io)
     xocc=atp%orbit_mod_occ(io)
     nfix=nl(ip,l+1)
     IF (nfix<=0.OR.nfix>atp%Orbit%norbit) THEN
       LIBPAW_ERROR('error in occupations')
     ENDIF
     atp%Orbit%occ(nfix)=xocc
   END DO
   if(has_to_print) WRITE(STD_OUT,*) ' Corrected occupations are: '
   if(has_to_print) WRITE(STD_OUT,"(' n  l     occupancy')")
   atp%electrons=0._dp
   DO io=1,atp%Orbit%norbit
      if(has_to_print) WRITE(STD_OUT,'(i2,1x,i2,4x,1p,1e15.7)')  &
  &        atp%Orbit%np(io),atp%Orbit%l(io),atp%Orbit%occ(io)
      atp%electrons=atp%electrons+atp%Orbit%occ(io)
      if(.not.atp%Orbit%iscore(io)) atp%Orbit%qval=atp%Orbit%qval+atp%Orbit%occ(io)
   ENDDO
 ENDIF ! scalarrelativistic

 If (atp%diracrelativistic) then 
   i=MAX(atp%Orbit%nps,atp%Orbit%npp,atp%Orbit%npd,atp%Orbit%npf,atp%Orbit%npg)
   LIBPAW_ALLOCATE(nl,(i,-5:5))
   nl=0
   IF (atp%Orbit%nps.GT.0) THEN
     DO is=1,atp%Orbit%nps
       icount=icount+1
       nl(is,-1)=icount
       atp%Orbit%occ(icount)=2._dp
       atp%Orbit%np(icount)=is
       atp%Orbit%l(icount)=0
       atp%Orbit%kappa(icount)=-1
       if(atp%orbit_iscore(icount)) atp%npsc=is
     ENDDO
   ENDIF
   IF (atp%Orbit%npp.GT.1) THEN
     DO ip=2,atp%Orbit%npp
       icount=icount+1
       nl(ip,1)=icount
       atp%Orbit%occ(icount)=2._dp
       atp%Orbit%np(icount)=ip
       atp%Orbit%l(icount)=1
       atp%Orbit%kappa(icount)=1
       if(atp%orbit_iscore(icount)) atp%nppc=ip
     ENDDO
     DO ip=2,atp%Orbit%npp
       icount=icount+1
       nl(ip,-2)=icount
       atp%Orbit%occ(icount)=4._dp
       atp%Orbit%np(icount)=ip
       atp%Orbit%l(icount)=1
       atp%Orbit%kappa(icount)=-2
       if(atp%orbit_iscore(icount)) atp%nppc=ip
     ENDDO
   ENDIF
   IF (atp%Orbit%npd.GT.2) THEN
     DO id=3,atp%Orbit%npd
       icount=icount+1
       nl(id,2)=icount
       atp%Orbit%occ(icount)=4._dp
       atp%Orbit%np(icount)=id
       atp%Orbit%l(icount)=2
       atp%Orbit%kappa(icount)=2
       if(atp%orbit_iscore(icount)) atp%npdc=id
     ENDDO
     DO id=3,atp%Orbit%npd
       icount=icount+1
       nl(id,-3)=icount
       atp%Orbit%occ(icount)=6._dp
       atp%Orbit%np(icount)=id
       atp%Orbit%l(icount)=2
       atp%Orbit%kappa(icount)=-3
       if(atp%orbit_iscore(icount)) atp%npdc=id
     ENDDO
   ENDIF
   IF (atp%Orbit%npf.GT.3) THEN
     DO jf=4,atp%Orbit%npf
       icount=icount+1
       nl(jf,3)=icount
       atp%Orbit%occ(icount)=6._dp
       atp%Orbit%np(icount)=jf
       atp%Orbit%l(icount)=3
       atp%Orbit%kappa(icount)=3
       if(atp%orbit_iscore(icount)) atp%npfc=jf
     ENDDO
     DO jf=4,atp%Orbit%npf
       icount=icount+1
       nl(jf,-4)=icount
       atp%Orbit%occ(icount)=8._dp
       atp%Orbit%np(icount)=jf
       atp%Orbit%l(icount)=3
       atp%Orbit%kappa(icount)=-4
       if(atp%orbit_iscore(icount)) atp%npfc=jf
     ENDDO
   ENDIF
   IF(atp%Orbit%npg.GT.4) THEN
     DO ig=5,atp%Orbit%npg
        icount=icount+1
        nl(ig,4)=icount
        atp%Orbit%occ(icount)=8._dp
        atp%Orbit%np(icount)=ig
        atp%Orbit%l(icount)=4
        atp%Orbit%kappa(icount)=4
        if(atp%orbit_iscore(icount)) atp%npgc=ig
     ENDDO
     DO ig=5,atp%Orbit%npg
       icount=icount+1
       nl(ig,-5)=icount
       atp%Orbit%occ(icount)=10._dp
       atp%Orbit%np(icount)=ig
       atp%Orbit%l(icount)=4
       atp%Orbit%kappa(icount)=-5
       if(atp%orbit_iscore(icount)) atp%npgc=ig
     ENDDO
   ENDIF
   atp%Orbit%norbit=icount
   if(has_to_print) write(std_out,*)' Below are listed the default occupations '
   if(has_to_print) write(std_out,"(' n  l kappa     occupancy')")
   DO io=1,atp%Orbit%norbit
     if(has_to_print) write(std_out,'(i2,1x,i2,3x,i2,4x,1p,1e15.7)')  &
& atp%Orbit%np(io),atp%Orbit%l(io),atp%Orbit%kappa(io),atp%Orbit%occ(io)
   ENDDO
   !Corrected occupations (from input dataset)
   atp%Orbit%iscore=atp%orbit_iscore
   DO io=1,atp%norbit_mod
     l=atp%orbit_mod_l(io)
     ip=atp%orbit_mod_n(io)
     kappa=atp%orbit_mod_k(io)
     xocc=atp%orbit_mod_occ(io)
     nfix=nl(ip,kappa)
     IF (nfix<=0.OR.nfix>atp%Orbit%norbit) THEN
       WRITE(STD_OUT,*) 'error in occupations -- ip,l,kappa,xocc:',&
&       ip,l,kappa,xocc,nfix,atp%Orbit%norbit
       STOP
     ENDIF
     atp%Orbit%occ(nfix)=xocc
    END DO
    if(has_to_print)WRITE(STD_OUT,*) ' Corrected occupations are: '
    if(has_to_print) WRITE(STD_OUT,"(' n  l  kappa   occupancy')")
    atp%electrons=0.d0
    DO io=1,atp%Orbit%norbit
       if(has_to_print) WRITE(STD_OUT,'(i2,1x,i2,3x,i2,4x,1p,1e15.7)')  &
  &         atp%Orbit%np(io),atp%Orbit%l(io),atp%Orbit%kappa(io),atp%Orbit%occ(io)
       atp%electrons=atp%electrons+atp%Orbit%occ(io)
       if(.not.atp%Orbit%iscore(io)) atp%Orbit%qval=atp%Orbit%qval+atp%Orbit%occ(io)
    ENDDO
 ENDIF    !   completed occupations
 LIBPAW_DEALLOCATE(nl)
end subroutine Prepare_Orbit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    subroutine print_atompaw_params(atp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_check_atompaw_params(atp)
 type(atompaw_type),intent(in) :: atp
 integer :: norb,io
 if(atp%finitenucleus) then
   LIBPAW_ERROR('Finitenucleus not implemented')
 endif
 if(atp%HFpostprocess) then
   LIBPAW_ERROR('HFpostprocess not implemented')
 endif
 if(atp%localizedcoreexchange) then
   LIBPAW_ERROR('Localized core exchange not implemented')
 endif
 if(atp%BDsolve) then
   LIBPAW_ERROR('BD solver not implemented')
 endif
 if(atp%fixed_zero) then
   LIBPAW_ERROR('Fixed zero not implemented')
 endif
 if(atp%shapetcore) then
   LIBPAW_ERROR('Shapetcore not implemented')
 endif
 if(atp%ColleSalvetti) then
   LIBPAW_ERROR('ColleSalvetti not implemented')
 endif
 if(.not.atp%ortho_type==ORTHO_TYPE_VANDERBILT) then
   LIBPAW_ERROR('RCPAW only possible with Vanderbilt ortho scheme')
 endif
 if(has_to_print) then
   WRITE(STD_OUT,'(/,3x,a)') "===== Atompaw parameters ====="
  ! Atom
   WRITE(STD_OUT,'(3x,a,a2)') "Atomic symbol : ",atp%atomic_symbol
   WRITE(STD_OUT,'(3x,a,i0)') "Atomic charge : ",atp%atomic_charge
  ! Algo
   WRITE(STD_OUT,'(3x,2a)')     "Scalar-relativistic calculation:",MERGE("YES"," NO",atp%scalarrelativistic)
   WRITE(STD_OUT,'(3x,2a)')     "Dirac-relativistic calculation:",MERGE("YES"," NO",atp%diracrelativistic)
   WRITE(STD_OUT,'(3x,2a)')     "Frozen core calculation:",MERGE("YES"," NO",atp%frozencorecalculation)
   WRITE(STD_OUT,'(3x,2a)')     "Frozen vale calculation:",MERGE("YES"," NO",atp%frozenvalecalculation)
   IF (atp%usespline) THEN
     WRITE(STD_OUT,'(3x,a)')    "    - Use a spline solver"
   END IF
   WRITE(STD_OUT,'(3x,2a)')     "Exchange-correlation functional:",TRIM(atp%exctype)
   WRITE(STD_OUT,'(3x,3a)')     " (mGGA kinetic energy functional:",MERGE("YES"," NO",atp%needvtau),")"
   WRITE(STD_OUT,'(3x,2a)')     "Finite-nucleus calculation:",MERGE("YES"," NO",atp%finitenucleus)
   IF (atp%finitenucleus) THEN
     WRITE(STD_OUT,'(3x,a,i0)') "    - Finite-nucleus model:",atp%finitenucleusmodel
   END IF
   WRITE(STD_OUT,'(3x,2a)')     "Block-Davidson calculation:",MERGE("YES"," NO",atp%BDsolve)
  ! Grid
   WRITE(STD_OUT,'(3x,a,i0)')     "Grid type:",atp%Grid%type
   WRITE(STD_OUT,'(3x,a,i0)')   "Grid size:",atp%Grid%n
   WRITE(STD_OUT,'(3x,a,f7.3)') "Grid maximum value:",atp%Grid%r(atp%Grid%n)
   if(atp%usespline) then
   WRITE(STD_OUT,'(3x,a,f7.3,2x,i0)') "Spline grid r0, ns              :",&
&      atp%splr0,atp%splns
   endif
   ! XC
   WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, post-processing:",MERGE("YES"," NO",atp%HFpostprocess)
   WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, localized coreex.:",MERGE("YES"," NO",atp%localizedcoreexchange)
   WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, fixed zero:",MERGE("YES"," NO",atp%fixed_zero)
   IF (atp%fixed_zero) THEN
     WRITE(STD_OUT,'(3x,a,i0)') "    - HF fixed zero index:",atp%fixed_zero_index
   END IF
   IF (atp%BDsolve.and.atp%gridkey=='LINEAR') THEN
     WRITE(STD_OUT,'(/,3x,a)') "WARNING: BlockDavidson solver works very slowlywith linear grid!"
   END IF
   WRITE(STD_OUT,'(3x,a,5(1x,i0))') "Max. quantum numbers (s,p,d,f,g):",atp%np(1:5)
   WRITE(STD_OUT,'(3x,a,i0)') "Total number of orbitals: ",atp%norbit
   WRITE(STD_OUT,'(3x,a)') "Core and valence orbitals:"
   IF (.NOT.atp%diracrelativistic) WRITE(STD_OUT,'(7x,a)') "n l : type"
   IF (atp%diracrelativistic)      WRITE(STD_OUT,'(7x,a)') "n l kappa :type"
!   io=0
!   DO ll=0,4
!     nn=atp%np(ll+1)
!     IF (nn>0) THEN
!       IF (.NOT.atp%diracrelativistic) THEN
!         DO ii=1+ll,nn
!           io=io+1
!           WRITE(STD_OUT,'(7x,i1,1x,i1,2a)') ii,ll," : ", &
!&            MERGE("CORE   ","VALENCE",atp%orbit_iscore(io))
!         END DO
!       ELSE
!         DO ik=1,nkappa(ll+1)
!           kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
!           DO ii=1+ll,nn
!             io=io+1
!             WRITE(STD_OUT,'(7x,i1,1x,i1,2x,i2,2x,2a)') ii,ll,kk," : ", &
!         &      MERGE("CORE   ","VALENCE",atp%orbit_iscore(io))
!           END DO
!         END DO
!       END IF
!     END IF
!   END DO
   WRITE(STD_OUT,'(3x,a,i0)') "Basis, maximum L : ",atp%lmax
  !cutoff radii
   WRITE(STD_OUT,'(3x,a,f7.4)') "Augmentation region radius : ",atp%rc
   WRITE(STD_OUT,'(3x,a,f7.4)') "Core dens. matching radius : ",atp%rc_core
   WRITE(STD_OUT,'(3x,a,f7.4)') "Local pot. matching radius : ",atp%rc_vloc
   WRITE(STD_OUT,'(3x,a,f7.4)') "Compens. shape func radius : ",atp%rc_shap
  !Additional basis
   WRITE(STD_OUT,'(3x,a,i0)') "Initial number of basis functions:",atp%nbasis-atp%nbasis_add
   WRITE(STD_OUT,'(3x,a,i0)') "Number of additional basis functions:",atp%nbasis_add
   WRITE(STD_OUT,'(3x,a,i0)') "Total number of basis functions:",atp%nbasis
   WRITE(STD_OUT,'(3x,a)') "Additional basis functions:"
   IF (.NOT.atp%diracrelativistic) THEN
     WRITE(STD_OUT,'(7x,a)') "l : energy"
     DO io=1,atp%nbasis_add
       WRITE(STD_OUT,'(7x,i1,a,f7.4)') atp%basis_add_l(io),":",atp%basis_add_energy(io)
     END DO
   ELSE
     WRITE(STD_OUT,'(7x,a)') "l kappa : energy"
     DO io=1,atp%nbasis_add
       WRITE(STD_OUT,'(7x,i1,2x,i2,2x,a,f7.4)') atp%basis_add_l(io), &
&          atp%basis_add_k(io)," : " ,atp%basis_add_energy(io)
     END DO
   END IF
  !Projectors
   WRITE(STD_OUT,'(3x,a)') "Projectors description:"
   IF (atp%projector_type==PROJECTOR_TYPE_BLOECHL) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : BLOECHL"
   IF (atp%projector_type==PROJECTOR_TYPE_VANDERBILT) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : VANDERBILT"
   IF (atp%projector_type==PROJECTOR_TYPE_MODRRKJ) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : MODRRKJ"
   IF (atp%projector_type==PROJECTOR_TYPE_CUSTOM) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : CUSTOM"
   IF (atp%projector_type==PROJECTOR_TYPE_HF) &
&    WRITE(STD_OUT,'(7x,a)') "Type : HARTREE-FOCK"
   IF (atp%projector_type/=PROJECTOR_TYPE_HF) THEN
     IF (atp%pseudo_type==PSEUDO_TYPE_BLOECHL) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : BLOECHL"
     IF (atp%pseudo_type==PSEUDO_TYPE_POLYNOM) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : POLYNOM"
     IF (atp%pseudo_type==PSEUDO_TYPE_RRKJ) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : RRKJ"
     IF (atp%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : BLOECHL KERKER"
     IF (atp%pseudo_type==PSEUDO_TYPE_POLYNOM2) &
&      WRITE(STD_OUT,'(7x,a,i0,a,es9.3)') "Pseudization      : POLYNOM2,pdeg=",&
&       atp%pseudo_polynom2_pdeg,", qcut=",atp%pseudo_polynom2_qcut
     IF (atp%ortho_type==ORTHO_TYPE_GRAMSCHMIDT) &
&      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : GRAM-SCHMIDT"
     IF (atp%ortho_type==ORTHO_TYPE_VANDERBILT) &
&      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : VANDERBILT"
     IF (atp%ortho_type==ORTHO_TYPE_SVD) &
&      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : SVD"
   END IF
   IF (atp%shapefunc_type==SHAPEFUNC_TYPE_GAUSSIAN) &
&    WRITE(STD_OUT,'(3x,a,es9.3)') "Compensation charge shape function: GAUSSIAN,tol=",&
&    atp%shapefunc_gaussian_param
   IF (atp%shapefunc_type==SHAPEFUNC_TYPE_SINC) &
&    WRITE(STD_OUT,'(3x,a)') "Compensation charge shape function : SINC2"
   IF (atp%shapefunc_type==SHAPEFUNC_TYPE_BESSEL) &
&    WRITE(STD_OUT,'(3x,a)') "Compensation charge shape function : BESSEL"
   IF (atp%hf_coretol>0) &
&    WRITE(STD_OUT,'(3x,a,es9.3)') "Core tolerance for Hartree-Fock:",atp%hf_coretol
   WRITE(STD_OUT,'(3x,2a)') "Smooth tcore shape (no negative nhat):",MERGE("YES"," NO",atp%shapetcore)
  !Local Psp
      IF (atp%vloc_type==VLOC_TYPE_MTROULLIER) &
&    WRITE(STD_OUT,'(7x,a,i0,a,f7.4)') "Local pseudopotential type:MTROULLIER,l=",&
&          atp%vloc_l,", energy=",atp%vloc_ene
   IF (atp%vloc_type==VLOC_TYPE_ULTRASOFT) &
&    WRITE(STD_OUT,'(7x,a,i0,a,f7.4)') "Local pseudopotential type:ULTRASOFT,l=",&
&          atp%vloc_l,", energy=",atp%vloc_ene
   IF (atp%vloc_type==VLOC_TYPE_BESSEL) &
&    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : BESSEL"
   IF (atp%vloc_type==VLOC_TYPE_VPSMATCHNC) &
&    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNC"
   IF (atp%vloc_type==VLOC_TYPE_VPSMATCHNNC) &
&    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNNC"
   IF (atp%vloc_type==VLOC_TYPE_SETVLOC) THEN
     WRITE(STD_OUT,'(7x,a,es9.4,a,es9.4)') "Local pseudopotential type:SETVLOC,coef=",&
&          atp%vloc_setvloc_coef,", rad=",atp%vloc_setvloc_rad
     IF (atp%needvtau) THEN
       LIBPAW_ERROR('SETVLOC  option not available for MGGA')
     ENDIF
   ENDIF
   IF (atp%vloc_type==VLOC_TYPE_KERKER_EXPF) &
&    WRITE(STD_OUT,'(7x,a,4(1x,i0))') "Local pseudopotential type : KERKER EXPF,powers=",&
&          atp%vloc_kerker_power(1:4)
   IF (atp%vloc_type==VLOC_TYPE_KERKER_POLY) &
&    WRITE(STD_OUT,'(7x,a,4(1x,i0))') "Local pseudopotential type : KERKER POLY,powers=",&
&          atp%vloc_kerker_power(1:4)
   IF (atp%vloc_type==VLOC_TYPE_MTROULLIER.AND.atp%needvtau) THEN
     WRITE(STD_OUT,'(7x,a)') 'NOTE: MTROULLIER Vloc not available for mGGA!'
     WRITE(STD_OUT,'(7x,a)') '      Calling VPSmatch with norm conservation instead.'
     WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNC"
   END IF
   WRITE(STD_OUT,'(3x,a)') "Matching radius for basis functions:"
   IF (.NOT.atp%diracrelativistic) WRITE(STD_OUT,'(7x,a)') " # - n l : radius"
   IF (atp%diracrelativistic) WRITE(STD_OUT,'(7x,a)') " # - n l kappa :radius"
   norb=0
!   DO ll=0,atp%lmax
!     DO ik=1,MERGE(nkappa(ll+1),1,atp%diracrelativistic)
!       kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
!       IF (.NOT.atp%diracrelativistic) kk=0
!       DO io=1,atp%norbit_val
!         IF (atp%orbit_val_l(io)==ll.AND. &
! &          ((.NOT.atp%diracrelativistic).OR.atp%orbit_val_k(io)==kk))THEN
!           norb=norb+1
!           IF (.NOT.atp%diracrelativistic) &
! &           WRITE(STD_OUT,'(7x,i2,a,i1,1x,i1,a,f7.4)') &
! &           norb," - ",atp%orbit_val_n(io),ll," :",atp%basis_func_rc(norb)
!           IF (atp%diracrelativistic) &
! &           WRITE(STD_OUT,'(7x,i2,a,i1,1x,i1,2x,i2,2x,a,f7.4)') &
! &           norb," - ",atp%orbit_val_n(io),ll,kk," :",atp%basis_func_rc(norb)
!         END IF
!       END DO
!       IF (atp%nbasis_add>0) THEN
!         DO io=1,atp%nbasis_add
!           IF (atp%basis_add_l(io)==ll.AND. &
! &          ((.NOT.atp%diracrelativistic).OR.atp%basis_add_k(io)==kk))THEN
!             norb=norb+1
!             IF (.NOT.atp%diracrelativistic) &
! &             WRITE(STD_OUT,'(7x,i2,a,a1,1x,i1,a,f7.4)') &
! &             norb," - ",".",ll," : ",atp%basis_func_rc(norb)
!             IF (atp%diracrelativistic) &
! &             WRITE(STD_OUT,'(7x,i2,a,a1,1x,i1,2x,i2,2x,a,f7.4)') &
! &             norb," - ",".",ll,kk," : ",atp%basis_func_rc(norb)
!           END IF
!         END DO
!       END IF
!     END DO
!   END DO
 endif
end subroutine print_check_atompaw_params





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. aeatom 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Potential_Init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Potential_Init(Orbit,Pot,Grid)
 IMPLICIT NONE
 TYPE(OrbitInfo),INTENT(IN) :: Orbit
 TYPE(PotentialInfo),INTENT(INOUT) :: Pot
 TYPE(GridInfo),INTENT(IN) :: Grid
 real(dp) :: ecoul,v0
 CALL atompaw_poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v0)
 if(has_to_print) write(std_out,*) 'In Potential_Init', Pot%q,ecoul
END SUBROUTINE Potential_Init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SCFatom
!!    Main atomic SCF routine to calculate the all electron atomic solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SCFatom(atp,initialize)
 implicit none
 logical, intent(in) :: initialize
 type(atompaw_type),target, intent(inout) :: atp
 if(initialize) then
   call Orbit_Init(atp%Orbit,atp%Pot,atp)
 endif
 call Potential_Init(atp%Orbit,atp%Pot,atp%Grid)
 if(atp%Orbit%frozenvalecalculation) atp%Pot%q=atp%electrons
 call LDAGGA_SCF(atp)
end subroutine SCFatom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Orbit_Init
!!   From nuclear charge -- generate hydrogenic-like initial wfns
!!   and densities --
!!   fill AEOrbit%wfn, AEOrbit%eig, and AEOrbit%den and AEOrbit%q
!!   also AEOrbit%otau and AEOrbit%tau
!!   Note that both den and tau need to be divided by 4 \pi r^2
!!   Note that tau is only correct for non-relativistic case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Orbit_Init(Orbit,Pot,atp)
 IMPLICIT NONE
 TYPE(PotentialInfo),INTENT(INOUT) :: Pot
 TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
 type(atompaw_type), intent(in) :: atp
 INTEGER  :: io,ir,l,np,kappa
 real(dp) :: qcal,rescale,zeff,xocc,fpi
 REAL(dp),allocatable :: d(:)
 INTEGER :: initialconfig=0
 IF (initialconfig/=0) STOP 'Error in aeatom -- Orbit_Init already called'
 !  calculate initial charge density from hydrogen-like functions
 !  also initial energies
 zeff=Pot%nz
 If  (.not.atp%diracrelativistic) then
   DO io=1,Orbit%norbit
     if(Orbit%frozenvalecalculation.and.(.not.Orbit%iscore(io))) cycle
     np=Orbit%np(io)
     l=Orbit%l(io)
     xocc=Orbit%occ(io)
     if(.not.Orbit%frozenvalecalculation) Orbit%eig(io)=-(zeff/(np))**2!/2.0!Hartree
     if(has_to_print) WRITE(std_out,*) io,np,l,xocc,Orbit%eig(io)
     DO ir=1,atp%Grid%n
       Orbit%wfn(ir,io)=hwfn(zeff,np,l,atp%Grid%r(ir))
       IF (ABS(Orbit%wfn(ir,io))<machine_zero) Orbit%wfn(ir,io)=0._dp
     ENDDO
     zeff=zeff-0.5_dp*xocc
     zeff=MAX(zeff,1._dp)
   ENDDO
 endif
 If  (atp%diracrelativistic) then
   DO io=1,Orbit%norbit
     if(Orbit%frozenvalecalculation.and.(.not.Orbit%iscore(io))) cycle
     np=Orbit%np(io)
     l=Orbit%l(io)
     kappa=Orbit%kappa(io)
     xocc=Orbit%occ(io)
     if(.not.Orbit%frozenvalecalculation) Orbit%eig(io)=-(zeff/(np))**2!/2.0
     DO ir=1,atp%Grid%n
       call dirachwfn(np,kappa,zeff,atp%Grid%r(ir),Orbit%eig(io) &
&            ,Orbit%wfn(ir,io),Orbit%lwfn(ir,io))
       IF (ABS(Orbit%wfn(ir,io))<machine_zero) Orbit%wfn(ir,io)=0._dp
       IF (ABS(Orbit%lwfn(ir,io))<machine_zero) Orbit%lwfn(ir,io)=0._dp
     ENDDO
     zeff=zeff-0.5_dp*xocc
     zeff=MAX(zeff,1._dp)
   ENDDO
 endif
 ! check charge and rescale
 Orbit%den=0._dp
 Orbit%tau=0._dp
 DO io=1,Orbit%norbit
   if(Orbit%frozenvalecalculation.and.(.not.Orbit%iscore(io))) cycle
   CALL taufromwfn(Orbit%otau(:,io),atp%Grid,Orbit%wfn(:,io),Orbit%l(io),&
&                         energy=Orbit%eig(io),rPot=Pot%rv)
   xocc=Orbit%occ(io)
   DO ir=1,atp%Grid%n
     Orbit%den(ir)=Orbit%den(ir)+(Orbit%wfn(ir,io)**2)*xocc
     Orbit%tau(ir)=Orbit%tau(ir)+xocc*Orbit%otau(ir,io)
     If (atp%diracrelativistic) Orbit%den(ir)=Orbit%den(ir) + &
&                 xocc*((Orbit%lwfn(ir,io))**2)
   ENDDO
 ENDDO
!   Note that kinetic energy density (tau) is in Rydberg units ???
!   Note that kinetic energy density is only correct for non-relativistic
!               formulation
 qcal=integrator(atp%Grid,Orbit%den)
 if(Orbit%frozenvalecalculation) qcal=qcal+atp%Orbit%qval
 if(Orbit%frozenvalecalculation) Orbit%den=Orbit%den+Orbit%valeden
 !rescale density
 rescale=atp%electrons/qcal
 Orbit%den(1:atp%Grid%n)=Orbit%den(1:atp%Grid%n)*rescale
 Orbit%tau(1:atp%Grid%n)=Orbit%tau(1:atp%Grid%n)*rescale
 ! determine difference with tauW (Weizsaker)
 LIBPAW_ALLOCATE(d,(atp%Grid%n)); fpi=4*pi
 d(2:atp%Grid%n)=Orbit%den(2:atp%Grid%n)/(fpi*atp%Grid%r(2:atp%Grid%n)**2)
 call extrapolate(d)
 CALL derivative(atp%Grid,d,Orbit%deltatau)
 Do ir=1,atp%Grid%n
   if (d(ir)>machine_zero) then
     Orbit%deltatau(ir)=0.25_dp*(Orbit%deltatau(ir)**2)/d(ir)
   else
     Orbit%deltatau(ir)=0.0_dp
   endif
 enddo
 d(2:atp%Grid%n)=Orbit%tau(2:atp%Grid%n)/(fpi*atp%Grid%r(2:atp%Grid%n)**2)
 call extrapolate(d)
 Orbit%deltatau=d-Orbit%deltatau
 LIBPAW_DEALLOCATE(d)
END SUBROUTINE Orbit_Init





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. atomdata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Subroutine InitOrbit  -- used in CopyOrbit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitOrbit(Orbit,norbit,n,exctype,diracrelativistic,scalarrelativistic,frozencorecalculation,frozenvalecalculation)
 type (OrbitInfo), intent(inout) :: Orbit
 logical, intent(in) :: diracrelativistic
 logical, intent(in) :: scalarrelativistic
 logical, intent(in) :: frozencorecalculation,frozenvalecalculation
 integer, intent(in) :: n,norbit
 character(*),intent(in) :: exctype
 CALL DestroyOrbit(Orbit)
 Orbit%norbit=norbit;Orbit%exctype=trim(exctype)
 Orbit%nps=0;Orbit%npp=0;Orbit%npd=0;Orbit%npf=0;Orbit%npg=0
 Orbit%npsc=0;Orbit%nppc=0;Orbit%npdc=0;Orbit%npfc=0;Orbit%npgc=0
 Orbit%qval=0._dp
 LIBPAW_POINTER_ALLOCATE(Orbit%np,(norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%l,(norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%eig,(norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%occ,(norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%iscore,(norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%issemicore,(norbit))
 Orbit%iscore=.false.
 Orbit%issemicore=.false.
 Orbit%np=0;Orbit%l=0
 Orbit%eig=0._dp;Orbit%occ=0._dp
 LIBPAW_POINTER_ALLOCATE(Orbit%wfn,(n,norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%otau,(n,norbit))
 LIBPAW_POINTER_ALLOCATE(Orbit%den,(n))
 LIBPAW_POINTER_ALLOCATE(Orbit%tau,(n))
 LIBPAW_POINTER_ALLOCATE(Orbit%coreden,(n))
 LIBPAW_POINTER_ALLOCATE(Orbit%valeden,(n))
 LIBPAW_POINTER_ALLOCATE(Orbit%deltatau,(n))
 Orbit%wfn=0._dp;Orbit%den=0._dp;Orbit%tau=0._dp;Orbit%otau=0._dp
 Orbit%deltatau=0._dp
 Orbit%coreden=0._dp ; Orbit%valeden=0._dp
 Orbit%scalarrelativistic=scalarrelativistic
 Orbit%frozencorecalculation=frozencorecalculation
 Orbit%frozenvalecalculation=frozenvalecalculation
 Orbit%diracrelativistic=diracrelativistic
 If (Orbit%diracrelativistic) then
   LIBPAW_POINTER_ALLOCATE(Orbit%lwfn,(n,norbit))
   LIBPAW_POINTER_ALLOCATE(Orbit%kappa,(norbit))
   Orbit%lwfn=0._dp
   Orbit%kappa=0._dp
 else
   nullify(Orbit%lwfn,Orbit%kappa)
 endif
 if (exctype == "HF".or.exctype == "EXXKLI") then
   LIBPAW_POINTER_ALLOCATE(Orbit%lqp,(norbit,norbit))
   LIBPAW_POINTER_ALLOCATE(Orbit%X,(n,norbit))
 else
   nullify(Orbit%lqp,Orbit%X)
 endif
end subroutine InitOrbit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Subroutine DestroyOrbit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DestroyOrbit(Orbit)
 type (OrbitInfo), intent(inout) :: Orbit
 if (associated(Orbit%np)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%np)
 endif
 if (associated(Orbit%l)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%l)
 endif
 if (associated(Orbit%kappa)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%kappa)
 endif
 if (associated(Orbit%iscore)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%iscore)
 endif
 if (associated(Orbit%issemicore)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%issemicore)
 endif
 if (associated(Orbit%eig)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%eig)
 endif
 if (associated(Orbit%occ)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%occ)
 endif
 if (associated(Orbit%wfn)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%wfn)
 endif
 if (associated(Orbit%otau)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%otau)
 endif
 if (associated(Orbit%lwfn)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%lwfn)
 endif
 if (associated(Orbit%den)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%den)
 endif
 if (associated(Orbit%tau)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%tau)
 endif
 IF (ASSOCIATED(Orbit%deltatau)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%deltatau)
 endif
 if (associated(Orbit%lqp)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%lqp)
 endif
 if (associated(Orbit%X)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%X)
 endif
 if(associated(Orbit%coreden)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%coreden)
 endif
 if(associated(Orbit%valeden)) then
   LIBPAW_POINTER_DEALLOCATE(Orbit%valeden)
 endif
end subroutine DestroyOrbit


!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CopyOrbit(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CopyOrbit(SOrbit,COrbit)
 type(OrbitInfo),intent(inout)::SOrbit
 type(OrbitInfo),intent(inout)::COrbit
 integer::n
 n=size(SOrbit%den,1)
 call InitOrbit(COrbit,SOrbit%norbit,n,SOrbit%exctype,SOrbit%diracrelativistic,SOrbit%scalarrelativistic,&
&                  SOrbit%frozencorecalculation,SOrbit%frozenvalecalculation)
 COrbit%nps=SOrbit%nps
 COrbit%npp=SOrbit%npp
 COrbit%npd=SOrbit%npd
 COrbit%npf=SOrbit%npf
 COrbit%npg=SOrbit%npg
 COrbit%npsc=SOrbit%npsc
 COrbit%nppc=SOrbit%nppc
 COrbit%npdc=SOrbit%npdc
 COrbit%npfc=SOrbit%npfc
 COrbit%npgc=SOrbit%npgc
 COrbit%qval=SOrbit%qval
 COrbit%np(1:SOrbit%norbit)=SOrbit%np(1:SOrbit%norbit)
 COrbit%l(1:SOrbit%norbit)=SOrbit%l(1:SOrbit%norbit)
 COrbit%eig(1:SOrbit%norbit)=SOrbit%eig(1:SOrbit%norbit)
 COrbit%occ(1:SOrbit%norbit)=SOrbit%occ(1:SOrbit%norbit)
 COrbit%wfn(:,1:SOrbit%norbit)=SOrbit%wfn(:,1:SOrbit%norbit)
 COrbit%otau(:,1:SOrbit%norbit)=SOrbit%otau(:,1:SOrbit%norbit)
 COrbit%iscore(1:SOrbit%norbit)=SOrbit%iscore(1:SOrbit%norbit)
 COrbit%issemicore(1:SOrbit%norbit)=SOrbit%issemicore(1:SOrbit%norbit)
 COrbit%den=SOrbit%den
 COrbit%coreden=SOrbit%coreden
 COrbit%valeden=SOrbit%valeden
 COrbit%tau=SOrbit%tau
 COrbit%deltatau=SOrbit%deltatau
 if (SOrbit%diracrelativistic) then
   COrbit%lwfn(:,1:SOrbit%norbit)=SOrbit%lwfn(:,1:SOrbit%norbit)
   COrbit%kappa(1:SOrbit%norbit)=SOrbit%kappa(1:SOrbit%norbit)
 endif
 if (SOrbit%exctype == "HF".or.SOrbit%exctype == "EXXKLI") then
   COrbit%X(:,1:SOrbit%norbit)=SOrbit%X(:,1:SOrbit%norbit)
   COrbit%lqp(1:SOrbit%norbit,1:SOrbit%norbit)=SOrbit%lqp(1:SOrbit%norbit,1:SOrbit%norbit)
 endif
end subroutine CopyOrbit


!!!!!!!!!!!!!!!!!!!!!!!!!
!!  InitPot
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitPot(Pot,n)
 integer, intent(in) :: n
 type (PotentialInfo), intent(inout) :: Pot
 CALL DestroyPot(Pot)
 LIBPAW_POINTER_ALLOCATE(Pot%rv,(n))
 LIBPAW_POINTER_ALLOCATE(Pot%rvn,(n))
 LIBPAW_POINTER_ALLOCATE(Pot%rvh,(n))
 LIBPAW_POINTER_ALLOCATE(Pot%rvx,(n))
 LIBPAW_POINTER_ALLOCATE(Pot%vtau,(n))
 LIBPAW_POINTER_ALLOCATE(Pot%ww,(n))
 LIBPAW_POINTER_ALLOCATE(Pot%jj,(n))
 Pot%rv=0._dp;Pot%rvn=0._dp;Pot%rvh=0._dp;Pot%rvx=0._dp;Pot%vtau=0._dp
 Pot%ww=0._dp;Pot%jj=0._dp
end subroutine InitPot


!!!!!!!!!!!!!!!!!!!!!!!!!
!!  DestroyPot
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DestroyPot(Pot)
 type (PotentialInfo), intent(inout) :: Pot
 if (associated(Pot%rv)) then
   LIBPAW_POINTER_DEALLOCATE(Pot%rv)
 endif
 if (associated(Pot%rvn)) then
   LIBPAW_POINTER_DEALLOCATE(Pot%rvn)
 endif
 if (associated(Pot%rvh)) then
   LIBPAW_POINTER_DEALLOCATE(Pot%rvh)
 endif
 if (associated(Pot%rvx)) then
   LIBPAW_POINTER_DEALLOCATE(Pot%rvx)
 endif
 if (associated(Pot%vtau)) then
   LIBPAW_POINTER_DEALLOCATE(Pot%vtau)
 endif
 if (allocated(Pot%ww)) then
   LIBPAW_DEALLOCATE(Pot%ww)
 endif
 if (allocated(Pot%jj)) then
   LIBPAW_DEALLOCATE(Pot%jj)
 endif
end subroutine DestroyPot


!!!!!!!!!!!!!!!!!!!!!!!!!
!  CopyPot(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CopyPot(SPot,CPot)
 type(PotentialInfo),intent(in) :: SPot
 type(PotentialInfo),intent(inout) :: CPot
 integer :: n
 n=SIZE(SPot%rv,1)
 CALL InitPot(CPot,n)
 CPot%nz=SPot%nz
 CPot%zz=SPot%zz
 CPot%sym=SPot%sym
 CPot%q=SPot%q
 CPot%v0=SPot%v0
 CPot%v0p=SPot%v0p
 CPot%finitenucleusmodel=SPot%finitenucleusmodel
 CPot%finitenucleus=SPot%finitenucleus
 CPot%Nv0=SPot%Nv0
 CPot%Nv0p=SPot%Nv0p
 CPot%rv(1:n)=SPot%rv(1:n)
 CPot%rvn(1:n)=SPot%rvn(1:n)
 CPot%rvh(1:n)=SPot%rvh(1:n)
 CPot%rvx(1:n)=SPot%rvx(1:n)
 CPot%vtau(1:n)=SPot%vtau(1:n)
 CPot%ww(1:n)=SPot%ww(1:n)
 CPot%jj(1:n)=SPot%jj(1:n)
 CPot%needvtau=SPot%needvtau
end subroutine CopyPot


!!!!!!!!!!!!!!!!!!!!!!!!!
!  InitSCF
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitSCF(SCF)
 type(SCFInfo),intent(inout)::SCF
 SCF%iter=0
 SCF%delta=0._dp;SCF%eone=0._dp;SCF%ekin=0._dp;SCF%estatic=0._dp
 SCF%ecoul=0._dp;SCF%eexc=0._dp;SCF%oepcs=0._dp;SCF%etot=0._dp
 SCF%valekin=0._dp;SCF%valecoul=0._dp;SCF%valeexc=0._dp
 SCF%corekin=0._dp;SCF%evale=0._dp
end subroutine InitSCF





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3. ldagga_mod 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!
!  LDAGGA_SCF
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LDAGGA_SCF(atp)
 TYPE(atompaw_type), INTENT(inout) :: atp
 TYPE(Anderson_context):: AC
 INTEGER :: n
 real(dp) :: en1,etxc,eex,x,y
 real(dp), ALLOCATABLE :: arg(:)
 LOGICAL :: success
 character(132):: exctypesave
 n=atp%Grid%n
 LIBPAW_ALLOCATE(arg,(n))
 if(atp%needvtau) then ! first converge LDA
   atp%needvtau=.false.
   atp%Pot%needvtau=.false.
   exctypesave=atp%exctype
   atp%exctype="LDA-PW"
   call initexch(atp)
   CALL exch(atp%Grid,atp%Orbit%den,atp%Pot%rvx,etxc,eex,itype=atp%itype,&
&    tau=atp%Orbit%tau,vtau=atp%Pot%vtau,xc_functionals=atp%xc_functionals)
   arg=atp%Pot%rvh+atp%Pot%rvx   ! iterating only on electronic part of pot
   atp%Pot%rv=atp%Pot%rvh+atp%Pot%rvx-atp%Pot%rvx(1)
   CALL zeropot(atp%Grid,atp%Pot%rv,atp%Pot%v0,atp%Pot%v0p)
   atp%Pot%rv=atp%Pot%rv+atp%Pot%rvn+atp%Pot%rvx(1)
   CALL InitAnderson_dr(AC,6,5,n,0.2d0,1.d3,100,seterr,settoosmall,.true.)
   CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success,atp)
   if(has_to_print) WRITE(STD_OUT,*) 'Anderson Mix with LDA',AC%res ,' iter = ',AC%CurIter
   CALL FreeAnderson(AC)
   if(has_to_print) write(STD_OUT,*) 'Completed initial iteration '
   atp%needvtau=.true.
   atp%Pot%needvtau=.true.
   atp%exctype=exctypesave
   call initexch(atp)
 endif
 CALL exch(atp%Grid,atp%Orbit%den,atp%Pot%rvx,etxc,eex,itype=atp%itype,&
&       needvtau=atp%Pot%needvtau,tau=atp%Orbit%tau,vtau=atp%Pot%vtau,&
&       xc_functionals=atp%xc_functionals)
 atp%Pot%rv=atp%Pot%rvh+atp%Pot%rvx-atp%Pot%rvx(1)
 CALL zeropot(atp%Grid,atp%Pot%rv,atp%Pot%v0,atp%Pot%v0p)
 atp%Pot%rv=atp%Pot%rv+atp%Pot%rvn+atp%Pot%rvx(1)
 atp%SCF%iter=0
 atp%SCF%delta=0
 If (atp%needvtau) then
   x=seterrmg; y=settoosmallmg
 else
   x=seterr; y=settoosmall
 endif 
 CALL InitAnderson_dr(AC,6,5,n,MaxMix,1.d3,MaxIter,x,y,.false.)
 If (atp%needvtau) then
   arg=atp%Orbit%den   ! iterating on density
   CALL DoAndersonMix(AC,arg,en1,DENITERsub,success,atp)
   !!!!   evaluate vxc and vtau on universal grid
   !!!!    because it may be bumpy at intermediate range from spline
   !!!!     evaluation
   CALL exch(atp%Grid,atp%Orbit%den,atp%Pot%rvx,etxc,eex,itype=atp%itype,&
&     tau=atp%Orbit%tau,vtau=atp%Pot%vtau)
   atp%Pot%rv=atp%Pot%rvn+atp%Pot%rvh+atp%Pot%rvx
 else
   arg=atp%Pot%rvh+atp%Pot%rvx
   CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success,atp)
 endif
 atp%SCF%iter=AC%CurIter
 atp%SCF%delta=AC%res
 if(has_to_print) WRITE(STD_OUT,*) 'Anderson Mix ',success,AC%res ,' iter = ',AC%CurIter
   if (AC%res>1._dp) then
     LIBPAW_ERROR('Sadly the program has not converged, Consider trying splineinterp')     
   endif   
 CALL FreeAnderson(AC)
 if(has_to_print) write(std_out,*) 'Finished Anderson Mix', en1 ,' success = ', success
 LIBPAW_DEALLOCATE(arg)
END SUBROUTINE LDAGGA_SCF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  LDAGGASub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LDAGGASub(w,energy,residue,err,success,update,atp)
 real(dp), INTENT(INOUT) :: w(:)
 real(dp), INTENT(OUT) :: energy
 real(dp), INTENT(OUT) :: residue(:)
 real(dp), INTENT(OUT) :: err
 LOGICAL, INTENT(OUT) :: success
 LOGICAL, INTENT(IN) :: update
 type(atompaw_type), intent(inout) :: atp
 INTEGER :: i,j,n,io,nw
 real(dp),ALLOCATABLE :: dum(:)
 real(dp) :: x
 TYPE (OrbitInfo) :: tmpOrbit
 TYPE (PotentialInfo) :: tmpPot
 n=atp%Grid%n
 nw=SIZE(w)
 LIBPAW_ALLOCATE(dum,(nw))
 CALL CopyOrbit(atp%Orbit,tmpOrbit)
 CALL CopyPot(atp%Pot,tmpPot)
 CALL Updatewfn(atp%Grid,tmpPot,tmpOrbit,w,success,atp%BDsolve,atp%usespline,atp%spline,atp%itype)
 if(has_to_print) write(std_out,*) 'completed updatewfn with success ', success
 If (.not.success) then   !  attempt to stablize solution
   !w=w+tmpPot%rvn
   if(has_to_print) write(std_out,*) 'Current eigs', (atp%Orbit%eig(io),io=1,atp%Orbit%norbit)
   j=n
   x=atp%Orbit%eig(1)
   if (atp%Orbit%norbit>1) then
     do io = 2, atp%Orbit%norbit
       if (atp%Orbit%eig(io)<0._dp.and.x<atp%Orbit%eig(io)) &
&        x=atp%Orbit%eig(io)
     enddo
   endif
   x=1._dp/sqrt(abs(x))
   j=FindGridIndex(atp%Grid,x)
   if(has_to_print) write(std_out,*) 'index', x,j,atp%Grid%r(j)
   if (j<10)j=10
   if (j>n-10) j=n-10
   w(j+1)=(-1._dp+w(j+1)/2)
   do i=j+2,n
      w(i)=-2._dp
   enddo
   if(has_to_print) write(std_out,*) 'Reset tmpPot ', j
   if(has_to_print) write(std_out,*) '   Last points '
   if(has_to_print) write(std_out,'(1p,20e15.7)') atp%Grid%r(n),w(n)
   !w=w-tmpPot%rvn
   CALL Updatewfn(atp%Grid,tmpPot,tmpOrbit,w,success,atp%BDsolve,atp%usespline,atp%spline,atp%itype)
   if(has_to_print) write(std_out,*) 'after updatwfn from reset ',success;
 Endif
 IF (.NOT.success) THEN
   if(has_to_print) write(std_out,*) 'Bad luck in Sub'
 ENDIF
 CALL Get_KinCoul(atp%Grid,tmpPot,tmpOrbit,atp%SCF,atp%usespline)
 CALL Get_EXC(atp,tmpPot,tmpOrbit)
 dum(1:n)=tmpPot%rvh(1:n)+tmpPot%rvx(1:n)-w(1:n)
 residue=dum
 err=Dot_Product(residue,residue)
 energy=atp%SCF%etot
 IF (update) THEN
   atp%Pot%rv=w+tmpPot%rvn
   atp%Pot%rvh=tmpPot%rvh
   atp%Pot%rvx=tmpPot%rvx
   if (atp%needvtau) atp%Pot%vtau=tmpPot%vtau
   atp%Orbit%wfn=tmpOrbit%wfn
   If(atp%diracrelativistic) atp%Orbit%lwfn=tmpOrbit%lwfn
   atp%Orbit%eig=tmpOrbit%eig
   atp%Orbit%den=tmpOrbit%den
   atp%Orbit%otau=tmpOrbit%otau
   atp%Orbit%tau=tmpOrbit%tau
   atp%Orbit%deltatau=tmpOrbit%deltatau
 ENDIF
 CALL DestroyPot(tmpPot)
 CALL DestroyOrbit(tmpOrbit)
 LIBPAW_DEALLOCATE (dum)
END SUBROUTINE  LDAGGASub


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  DENITERSub        -- w is the electron density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DENITERSub(w,energy,residue,err,success,update,atp)
  REAL(dp), INTENT(INOUT) :: w(:)
  REAL(dp), INTENT(OUT) :: energy
  REAL(dp), INTENT(OUT) :: residue(:)
  REAL(dp), INTENT(OUT) :: err
  type(atompaw_type), intent(inout) :: atp
  LOGICAL, INTENT(OUT) :: success
  LOGICAL, INTENT(IN) :: update
  INTEGER :: n,nw
  REAL(dp),ALLOCATABLE :: dum(:)
  REAL(dp) :: x,y
  TYPE (OrbitInfo) :: tmpOrbit
  TYPE (PotentialInfo) :: tmpPot
  n=atp%Grid%n
  nw=SIZE(w)
  if(n/=nw) then
    write(std_out,*) 'problem in DENITERsub n,nw', n,nw
    stop
  endif
  LIBPAW_ALLOCATE(dum,(nw))
  CALL CopyOrbit(atp%Orbit,tmpOrbit)
  CALL CopyPot(atp%Pot,tmpPot)
  !  w enters subroutine as next iteration density
  !  It needs to be renormalized 
  !  tmpPot%rvh, tmpPot%rvx,  tmpPot%vtau need to be generated 
  x=integrator(atp%Grid,w)
  if(has_to_print) write(std_out,*) 'In DENITERsub norm(m) adjust ', x,atp%Pot%q
  w=w*tmpPot%q/x
  call  atompaw_poisson(atp%Grid,x,w,tmpPot%rvh,y)
  if(has_to_print) write(std_out,*) 'after poisson  q,ecoul ',x,y
  call exch(atp%Grid,w,tmpPot%rvx,x,y,itype=atp%itype,tau=tmpOrbit%tau,vtau=tmpPot%vtau,xc_functionals=atp%xc_functionals)
  if(has_to_print) write(std_out,*) 'after exch   exvct, eexc ',x,y
  tmpPot%rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
  tmpOrbit%den=w
  CALL Updatewfnwden(atp%Grid,tmpPot,tmpOrbit,success,atp%usespline,atp%spline,atp%itype)
  if(has_to_print) write(std_out,*) 'completed updatewfnwden with success ', success
  IF (.NOT.success) THEN
     WRITE(STD_OUT,*) 'Bad luck in Sub'
  ENDIF
  CALL Get_KinCoul(atp%Grid,tmpPot,tmpOrbit,atp%SCF,atp%usespline)
  if(has_to_print) write(std_out,*)  'Check tau ', integrator(atp%Grid,tmpOrbit%tau)
  CALL Get_EXC(atp,tmpPot,tmpOrbit)
  dum(1:n)=tmpOrbit%den(1:n)-w(1:n)
  residue=dum
  err=Dot_Product(residue,residue)
  if(has_to_print) write(STD_OUT,*) 'in DENITERSub   err ', err
  energy=atp%SCF%etot
  IF (update) THEN
    atp%Pot%rv=tmpPot%rv
    atp%Pot%rvh=tmpPot%rvh
    atp%Pot%rvx=tmpPot%rvx
    if (atp%needvtau) atp%Pot%vtau=tmpPot%vtau
    atp%Orbit%wfn=tmpOrbit%wfn
    If(atp%diracrelativistic)atp%Orbit%lwfn=tmpOrbit%lwfn
    atp%Orbit%eig=tmpOrbit%eig
    atp%Orbit%den=w
    atp%Orbit%otau=tmpOrbit%otau
    atp%Orbit%tau=tmpOrbit%tau
    atp%Orbit%deltatau=tmpOrbit%deltatau
  ENDIF
  CALL DestroyOrbit(tmpOrbit)
  CALL DestroyPot(tmpPot)
  LIBPAW_DEALLOCATE(dum)
END SUBROUTINE  DENITERSub


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get_EXC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_EXC(atp,Pot,Orbit)
 !  program to calculate exc energy and potential
 !     assumes Orbit%den already known
 !     also assume kinetic and coulomb energies
 !       calculated and stored in SCF
 TYPE(Potentialinfo), INTENT(INOUT) :: Pot
 TYPE(Orbitinfo), INTENT(INOUT) :: Orbit
 TYPE(atompaw_type), INTENT(inout) :: atp
 real(dp) :: eex,etot,etxc
 real(dp), ALLOCATABLE :: dum(:)
 INTEGER :: n,fin
 n=atp%Grid%n
 fin=atp%grid%n
 if(atp%frozenvalecalculation) fin=atp%PAW%irc+5
 CALL exch(atp%Grid,Orbit%den,Pot%rvx,etxc,eex,itype=atp%itype,&
& needvtau=Pot%needvtau,tau=Orbit%tau,vtau=Pot%vtau,fin=fin,xc_functionals=atp%xc_functionals)
 atp%SCF%eexc=eex
 etot = atp%SCF%ekin+atp%SCF%estatic+atp%SCF%eexc
 atp%SCF%etot=etot
 if(has_to_print) WRITE(STD_OUT,*) '    Total                    :  ',etot
 LIBPAW_ALLOCATE(dum,(n))
 dum=0
 dum(2:n)=Pot%rvx(2:n)*Orbit%den(2:n)/atp%Grid%r(2:n)
 if (Pot%needvtau) then
   dum=dum+Orbit%tau*Pot%vtau       !Kinetic energy correction    
 endif
 if(has_to_print) then
  WRITE(STD_OUT,*) '    Total   (DC form)        :  ',&
&        atp%SCF%eone-atp%SCF%ecoul+eex-integrator(atp%Grid,dum)
 endif
 LIBPAW_DEALLOCATE(dum)
END SUBROUTINE Get_EXC





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4. excor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE initexch(atp)
!!    choose form of exchange-correlation potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initexch(atp)
 type(atompaw_type), intent(inout) :: atp
 integer :: id(2),i_plus,ii
 character*50 :: xcstrg(2)
 CALL Uppercase(atp%exctype)
 if(has_to_print) WRITE(STD_OUT,*) atp%exctype
 SELECT CASE(TRIM(atp%exctype))
 CASE default
   atp%itype = LIBXC
   i_plus=index(atp%exctype,'+')
   if (i_plus<=0) then
     xcstrg(1)=trim(atp%exctype)
     xcstrg(2)=""
   else
     xcstrg(1)=trim(atp%exctype(1:i_plus-1))
     xcstrg(2)=trim(atp%exctype(i_plus+1:))
   end if
   do ii=1,2
     id(ii)=libxc_functionals_getid(xcstrg(ii))
   enddo
   atp%ixc=-(id(1)*1000+id(2))
   if(has_to_print)WRITE(STD_OUT,*) 'Using Libxc --',TRIM(atp%exctype),atp%ixc
   call libxc_functionals_init(atp%ixc,1,atp%xc_functionals,el_temp=zero)
   if(atp%needvtau.and.(.not.libxc_functionals_ismgga(atp%xc_functionals))) then
     WRITE(STD_OUT,*) 'Problem with XC functional choice -- need mgga form for vtau '   
     WRITE(STD_OUT,*) '    Program stopping '
     stop
   endif    
   write(std_out,*) 'END INITEXCH'
 CASE('LDA-PW')
   atp%itype = LDA_PW
   if(has_to_print) WRITE(STD_OUT,*) 'Perdew-Wang correlation'
 CASE('GGA-PBE')
   atp%itype = GGA_PBE
   if(has_to_print) WRITE(STD_OUT,*) 'Perdew-Burke-Ernzerhof GGA'
 CASE('GGA-PBESOL')
   atp%itype = GGA_PBESOL
   if(has_to_print) WRITE(STD_OUT,*) 'Perdew-Burke-Ernzerhof modified (PBEsol) GGA'
 CASE ('MGGA-R2SCAN-001')
   atp%itype = MGGA_R2SCAN_001    
   if(has_to_print) WRITE(STD_OUT,*) 'R2SCAN MGGA with eta=0.001'
   call r2scaninit(0.001d0)
   atp%needvtau=.true.
 CASE ('MGGA-R2SCAN-01')
   atp%itype = MGGA_R2SCAN_01    
   if(has_to_print) WRITE(STD_OUT,*) 'R2SCAN MGGA with eta=0.01'
   call r2scaninit(0.01d0)
   atp%needvtau=.true.
 !CASE ('HF')
 !  itype = NO_XC
 !  WRITE(STD_OUT,*) 'No XC'
 END SELECT
END SUBROUTINE initexch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Logofterm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Logofterm(term)
 real(dp) :: term, Logofterm
 IF (ABS(term)>machine_precision) THEN
    Logofterm=ddlog(1._dp+term)
 ELSE
    Logofterm=term
 ENDIF
 RETURN
END FUNCTION Logofterm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine radialexcpbe
!   Density(:) input on a uniform radial mesh of Npts
!   Grid%r(:) input mesh points
!   Exc - output integrated exchange correlation energy   -- in Rydberg units
!   vxc(:) -- output exchange correlation potential       -- in Rydberg units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE radialexcpbe(Grid,density,Exc,vxc,mu,beta,fin)
 IMPLICIT NONE
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), intent(in) :: mu,beta
 REAL(dp), INTENT(IN) :: density(:)
 REAL(dp), INTENT(OUT) :: Exc, vxc(:)
 INTEGER, INTENT(IN), OPTIONAL :: fin
 INTEGER :: i,Npts
 REAL(dp),ALLOCATABLE :: gradient(:),gradmag(:),gxc(:),dgxcdr(:),fxc(:)
 REAL(dp) :: dfxcdn,dfxcdgbg
 Npts=Grid%n
 IF (PRESENT(fin)) Npts=fin
 LIBPAW_ALLOCATE(gradient,(Npts))
 LIBPAW_ALLOCATE(gradmag,(Npts))
 LIBPAW_ALLOCATE(gxc,(Npts))
 LIBPAW_ALLOCATE(dgxcdr,(Npts))
 LIBPAW_ALLOCATE(fxc,(Npts))
 CALL derivative(Grid,density(1:Npts),gradient(1:Npts),1,Npts)
 gradmag=ABS(gradient)
 DO i=1,Npts
   CALL  pbefunc(density(i),gradmag(i),fxc(i),dfxcdn,dfxcdgbg,mu,beta)
   vxc(i)=dfxcdn
   gxc(i)=dfxcdgbg*gradient(i)
 ENDDO
 CALL derivative(Grid,gxc(1:Npts),dgxcdr(1:Npts),1,Npts)
 DO i=2,Npts
   fxc(i)=2*fxc(i)*4*pi*(Grid%r(i)**2)  !2* changes from Har to Ryd
   vxc(i)=2*vxc(i)-2*dgxcdr(i)-4*gxc(i)/Grid%r(i)  ! Correction thanks
   ! to Marc Torrent and Francois Jollet
 ENDDO
 fxc(1)=zero
 CALL extrapolate(vxc)
 Exc = integrator(Grid,fxc,1,Npts)
 LIBPAW_DEALLOCATE(gradient)
 LIBPAW_DEALLOCATE(gradmag)
 LIBPAW_DEALLOCATE(gxc)
 LIBPAW_DEALLOCATE(dgxcdr)
 LIBPAW_DEALLOCATE(fxc)
 RETURN
END SUBROUTINE radialexcpbe


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! exch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE exch(Grid,den,rvxc,etxc,eexc,itype,fin,v0,v0p,needvtau,tau,vtau,xc_functionals)
 !  calculate exchange correlation potentials and energies
 !    for density functional theory from electron density
 !  den(n) is electron density * (4*pi*r**2)
 !  rvxc(n) is returned as vxc * r
 !  eexc is the total exchange energy (int(den*exc))
 !  etxc is eexc - int(den*vxc)
 !  icx id of the exchange-correlation functional
 !  xclevel level of the exchange-correlation functional used
 !  fin (optional) is integer range of densities and potentials
 !  v0  (optional) is extrapolated value of vxc for r=0
 !  v0p  (optional) is extrapolated value of dvxc/dr for r=0
 !  needvtau (optional) is logical .true. if mgga
 !  tau(n) (optional) is kinetic energy density * (4*pi*r**2)
 !  vtau(n) (optional) is kinetic energy contribution for mgga
 TYPE (GridInfo), INTENT(IN) :: Grid
 REAL(dp), INTENT(IN) :: den(:)
 REAL(dp), INTENT(INOUT) :: rvxc(:),etxc,eexc
 type(libxc_functional_type),intent(inout),optional :: xc_functionals(2)
 integer, intent(in) :: itype
 INTEGER, INTENT(IN), OPTIONAL :: fin
 REAL(dp), INTENT(OUT), OPTIONAL :: v0,v0p
 LOGICAL, INTENT(OUT), OPTIONAL :: needvtau
 REAL(dp), INTENT(IN),OPTIONAL :: tau(:)
 REAL(dp), INTENT(INOUT),OPTIONAL :: vtau(:)
 REAL(dp), ALLOCATABLE ::  tmpd(:),tmpv(:),dum(:)
 REAL(dp), ALLOCATABLE :: exci(:),dfxcdgbg(:,:),gxc(:),dgxcdr(:)
 REAL(dp), ALLOCATABLE :: grad(:),gradmag(:),dum1(:),sigma(:),tmpvt(:)
 REAL(dp), ALLOCATABLE :: tmpt(:),tmpl(:),dgxcdl(:),dexcdn(:),dexcds(:)
 REAL(dp) :: fpi,beta,mu
 INTEGER :: i,n,order,nspden,ndvxc,nd2vxc
 REAL(dp) :: r,r2,exc,vxc
 n=Grid%n
 IF (PRESENT(fin)) n=fin
 fpi=four*pi
 rvxc=0;etxc=0;eexc=0
 if (PRESENT(v0)) v0=0
 if (PRESENT(v0p)) v0p=0
 if (PRESENT(needvtau)) then
   if (needvtau.and.(PRESENT(vtau).eqv..false.)) then
     if(has_to_print) write(std_out,*) 'exch:  Inconsistency in mgga ', needvtau,PRESENT(vtau)
     LIBPAW_ERROR('exch:  stopping execution ')
    endif
 endif
 if (PRESENT(vtau)) vtau=0._dp
 If (itype==GGA_PBE.or.itype==GGA_PBESOL) then 
   LIBPAW_ALLOCATE(tmpd,(n))
   LIBPAW_ALLOCATE(tmpv,(n))
   tmpd=zero
   if(itype==GGA_PBE) then 
     mu=muorig;beta=betorig
   else
     mu=musol;beta=betsol
   endif
   DO i=2,n
     tmpd(i)=den(i)/(fpi*(Grid%r(i)**2))
   ENDDO
   CALL extrapolate(tmpd)
   IF (PRESENT(fin)) THEN
     CALL radialexcpbe(Grid,tmpd,eexc,tmpv,mu,beta,fin)
   ELSE
     CALL radialexcpbe(Grid,tmpd,eexc,tmpv,mu,beta)
   ENDIF
   IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
     CALL derivative(Grid,tmpv,tmpd,1,15)
     v0=tmpv(1)
     v0p=tmpd(1)
   ENDIF
   DO i=1,n
      rvxc(i)=tmpv(i)*Grid%r(i)
      tmpv(i)=tmpv(i)*den(i)
   ENDDO
   etxc=eexc-integrator(Grid,tmpv(1:n),1,n)
   LIBPAW_DEALLOCATE(tmpd)
   LIBPAW_DEALLOCATE(tmpv)
 ELSE IF (itype==LDA_PW) then !!! ! Perdew-Wang LDA !!!!
   LIBPAW_ALLOCATE(tmpd,(n))
   LIBPAW_ALLOCATE(tmpv,(n))
   LIBPAW_ALLOCATE(dum,(n))
   tmpd=0;tmpv=0;rvxc=0;dum=0
   DO i=2,n
     r=Grid%r(i)
     r2=r*r
     tmpd(i)=den(i)/(fpi*r2)
   ENDDO
   CALL extrapolate(tmpd)
   DO i=1,n
     CALL pwldafunc(tmpd(i),exc,vxc)
     tmpd(i)=den(i)*(exc-vxc)
     tmpv(i)=den(i)*exc
     rvxc(i)=Grid%r(i)*vxc
     IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
       IF (i==1) v0=vxc
       dum(i)=vxc
     ENDIF
   ENDDO
   etxc=integrator(Grid,tmpd(1:n),1,n)
   eexc=integrator(Grid,tmpv(1:n),1,n)
   IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
      CALL derivative(Grid,dum,tmpd,1,15)
      v0p=tmpd(1)
   ENDIF
   LIBPAW_DEALLOCATE(tmpd)
   LIBPAW_DEALLOCATE(tmpv)
   LIBPAW_DEALLOCATE(dum)
 ELSE IF (itype==MGGA_R2SCAN_001.or.itype==MGGA_R2SCAN_01) then !!r2scan 
   LIBPAW_ALLOCATE(tmpd,(n))
   LIBPAW_ALLOCATE(tmpv,(n))
   LIBPAW_ALLOCATE(exci,(n))
   LIBPAW_ALLOCATE(tmpt,(n))
   LIBPAW_ALLOCATE(dum,(n))
   LIBPAW_ALLOCATE(dum1,(n))
   tmpd=0.d0; tmpv=0.d0; exci=0.d0;tmpt=0.d0
   dum1=0.d0
   tmpd(2:n)=den(2:n)/(fpi*(Grid%r(2:n)**2))
   call extrapolate(tmpd)
   LIBPAW_ALLOCATE(grad,(n))
   LIBPAW_ALLOCATE(sigma,(n))
   LIBPAW_ALLOCATE(dexcdn,(n))
   LIBPAW_ALLOCATE(dexcds,(n))
   grad=0.d0;sigma=0.d0;dexcdn=0.d0;dexcds=0.d0
   do i=1,n
      dum1(i)=ddlog(tmpd(i))
   enddo
   call derivative(Grid,dum1,grad,1,n)
   grad(1:n)=grad(1:n)*tmpd(1:n)     ! perhaps more accurate???
   sigma=grad**2
   !   Prepare kinetic energy input tau -- used for most mgga
   tmpt(2:n)=tau(2:n)/(fpi*(Grid%r(2:n)**2))
   call extrapolate(tmpt)
   ! convert to Hartree units
   tmpt=0.5d0*tmpt
   do i=1,n
     call r2scanfun(tmpd(i),grad(i),tmpt(i),&
&             exci(i),vtau(i),dexcdn(i),dexcds(i))
   enddo
   ! convert to Rydberg units
   exci=two*exci
   dexcdn=two*dexcdn
   dexcds=four*dexcds       !extra factor of two due to sigma=grad**2
   rvxc=0.d0
   rvxc=dexcdn
   dum(1:n)=dexcds(1:n)
   call derivative(Grid,dum,dum1,1,n)
   dum(2:n)=two*dum(2:n)/(Grid%r(2:n))
   call extrapolate(dum)
   rvxc=rvxc-dum1-dum
   dum(1:n)=exci(1:n)*fpi*(Grid%r(1:n)**2)
   eexc=integrator(Grid,dum(1:n),1,n)
   dum(1:n)=rvxc(1:n)*den(1:n)
   etxc=eexc-integrator(Grid,dum(1:n),1,n)
   rvxc(1:n)=(dexcdn(1:n)-dum1(1:n))*Grid%r(1:n)-two*dexcds(1:n)
   LIBPAW_DEALLOCATE(grad)
   LIBPAW_DEALLOCATE(sigma)
   LIBPAW_DEALLOCATE(dexcdn)
   LIBPAW_DEALLOCATE(dexcds)
   LIBPAW_DEALLOCATE(tmpd)
   LIBPAW_DEALLOCATE(tmpv)
   LIBPAW_DEALLOCATE(exci)
   LIBPAW_DEALLOCATE(tmpt)
   LIBPAW_DEALLOCATE(dum)
   LIBPAW_DEALLOCATE(dum1)
 ELSE IF (itype==LIBXC) then
   ! Parameters
   if(.not.libxc_functionals_islda(xc_functionals).and..not.libxc_functionals_isgga(xc_functionals)) then
     STOP 
   endif
   order=1
   nspden=1
   ndvxc=0
   nd2vxc=0
   ! Allocations
   LIBPAW_ALLOCATE(tmpd,(n))
   LIBPAW_ALLOCATE(tmpv,(n))
   LIBPAW_ALLOCATE(exci,(n))
   LIBPAW_ALLOCATE(tmpt,(n))
   LIBPAW_ALLOCATE(grad,(n))
   LIBPAW_ALLOCATE(gradmag,(n))
   LIBPAW_ALLOCATE(gxc,(n))
   LIBPAW_ALLOCATE(dgxcdr,(n))
   LIBPAW_ALLOCATE(dfxcdgbg,(n,3))
   LIBPAW_ALLOCATE(tmpl,(n))
   LIBPAW_ALLOCATE(dgxcdl,(n))
   LIBPAW_ALLOCATE(tmpvt,(n))
   LIBPAW_ALLOCATE(dum,(n))
   LIBPAW_ALLOCATE(dum1,(n))
   tmpd=0._dp; tmpv=0._dp; exci=0._dp;tmpt=0._dp
   grad=0.d0;gradmag=0.d0;gxc=0.d0;dgxcdr=0.d0;dfxcdgbg=0.d0
   tmpl=0.d0; dgxcdl=0.d0;tmpvt=0.d0;dum=0.d0;dum1=0.d0
   ! Density
   tmpd(2:n)=den(2:n)/(fpi*(Grid%r(2:n)**2))
   call extrapolate(tmpd)
   ! Grad
   tmpv=0.d0
   do i=1,n
     tmpv(i)=ddlog(tmpd(i))
   enddo
   call derivative(Grid,tmpv,grad,1,n)
   grad(1:n)=grad(1:n)*tmpd(1:n)     !  perhaps more accurate???
   tmpv=0.d0
   gradmag=ABS(grad)*ABS(grad)
   ! Tau
   tmpt(2:n)=tau(2:n)/(fpi*(Grid%r(2:n)**2))
   call extrapolate(tmpt)
   ! Laplacian
   call derivative(Grid,grad,tmpl,1,n)
   tmpl(2:n)=tmpl(2:n)+2.d0*grad(2:n)/Grid%r(2:n)
   call extrapolate(tmpl)
   ! calc
   tmpd=tmpd/two
   gradmag=gradmag/four
   call libxc_functionals_getvxc(ndvxc,nd2vxc,n,nspden,order,tmpd,exci,tmpv,&
&   grho2=gradmag,lrho=tmpl,vxcgr=dfxcdgbg,vxclrho=dgxcdl,tau=tmpt,vxctau=tmpvt,&
&   xc_functionals=xc_functionals)
   tmpd=tmpd*two
   gradmag=gradmag*four
   ! Units 
   exci=two*exci
   tmpv=two*tmpv
   dfxcdgbg=dfxcdgbg*two
   ! Post-process
   gxc(1:n)=dfxcdgbg(1:n,3)*grad(1:n)
   call derivative(Grid,gxc,dgxcdr,1,n)
   tmpv(2:n)=tmpv(2:n)-dgxcdr(2:n)-2.d0*gxc(2:n)/Grid%r(2:n)
   call extrapolate(tmpv)
!   dum=0.d0
!   call derivative(Grid,dgxcdl,dum,1,n)
!   dum1=0.d0
!   call derivative(Grid,dum,dum1,1,n)
!   tmpv(2:n)=tmpv(2:n)+dum1(2:n)+2.d0*dum(2:n)/Grid%r(2:n)
!   call extrapolate(tmpv)
!   if (needvtau) vtau=tmpvt
   do i=1,n
     if (.not.ieee_is_normal(tmpv(i)).or.abs(tmpv(i)).lt.practical_zero) tmpv(i)=0.d0
     if (.not.ieee_is_normal(exci(i)).or.abs(exci(i)).lt.practical_zero)exci(i)=0.d0
     if (.not.ieee_is_normal(vtau(i)).or.abs(vtau(i)).lt.practical_zero)vtau(i)=0.d0
     if (.not.ieee_is_normal(tmpd(i)).or.abs(tmpd(i)).lt.practical_zero)tmpd(i)=0.d0
   enddo
   rvxc=0.d0
   rvxc(1:n)=tmpv(1:n)*Grid%r(1:n)
   exci(1:n)=exci(1:n)*tmpd(1:n)*fpi*Grid%r(1:n)**2
   eexc=integrator(Grid,exci,1,n)
   if (present(v0).and.present(v0p)) then
     call derivative(Grid,tmpv,tmpd,1,15)
     v0=tmpv(1);v0p=tmpd(1)
   endif
   tmpv(1:n)=tmpv(1:n)*den(1:n)
   etxc=eexc-integrator(Grid,tmpv(1:n),1,n)
   if(has_to_print) WRITE(STD_OUT,*) 'etxc,eexc = ',etxc,eexc
   ! DEALLOCATE
   LIBPAW_DEALLOCATE(tmpd)
   LIBPAW_DEALLOCATE(tmpv)
   LIBPAW_DEALLOCATE(exci)
   LIBPAW_DEALLOCATE(tmpt)
   LIBPAW_DEALLOCATE(grad)
   LIBPAW_DEALLOCATE(gradmag)
   LIBPAW_DEALLOCATE(gxc)
   LIBPAW_DEALLOCATE(dgxcdr)
   LIBPAW_DEALLOCATE(dfxcdgbg)
   LIBPAW_DEALLOCATE(tmpl)
   LIBPAW_DEALLOCATE(dgxcdl)
   LIBPAW_DEALLOCATE(tmpvt)
   LIBPAW_DEALLOCATE(dum)
   LIBPAW_DEALLOCATE(dum1)
 else
   WRITE(STD_OUT,*) 'Warning (EXCOR): ', itype,' no results returned !'
   STOP
 END if
END SUBROUTINE exch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! pwldafunc
!!  Subroutine to calculate the LDA exchange correlation functionals
!!  using the form of Perdew and Wang (PRB 45, 13244 (1992)
!!  assuming no spin polarization
!!  Inside this routine, energies are in Hartree units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE pwldafunc(den,exc,vxc)
 real(dp), INTENT(IN) :: den    !density
 real(dp), INTENT(OUT) :: exc,vxc
 real(dp), PARAMETER :: AA=0.0310907_dp
 real(dp), PARAMETER :: a1=0.21370_dp
 real(dp), PARAMETER :: b1=7.59570_dp
 real(dp), PARAMETER :: b2=3.58760_dp
 real(dp), PARAMETER :: b3=1.63820_dp
 real(dp), PARAMETER :: b4=0.49294_dp
 ! Variables depending on den
 real(dp) :: n,kf,rs,ks
 real(dp) :: ex,ec,pprs,decdrs
 real(dp) :: term
 n=den
 IF (n < machine_zero)  THEN
    exc=0._dp; vxc=0._dp
    RETURN
 ENDIF
 kf=(3._dp*(pi**2)*n)**0.3333333333333333333333333333_dp
 rs=(3._dp/(4._dp*pi*n))**0.3333333333333333333333333333_dp
 ks=SQRT(4._dp*kf/pi)
 ex=-3._dp*kf/(4._dp*pi)
 pprs=SQRT(rs)*(b1+b3*rs)+rs*(b2+b4*rs)
 term=Logofterm(1._dp/(2._dp*AA*pprs))
 ec=-2._dp*AA*(1._dp+a1*rs)*term
 exc=ex+ec
 decdrs=-(2._dp*AA*a1)*term &
&        +((1._dp+a1*rs)*((b1+3*b3*rs)/(2._dp*SQRT(rs))+&
&        b2+2*b4*rs))/(pprs*(pprs+1._dp/(2._dp*AA)))
 vxc = (4._dp/3._dp)*ex+ec-(decdrs*rs)/3._dp
 IF ((ABS(exc).GT.1.d65).OR.(ABS(vxc).GT.1.d65)) THEN
    if(has_to_print) WRITE(STD_OUT,*) 'Problem in PW',n,rs,ec
 ENDIF
 exc=2*exc; vxc=2*vxc      ! change to Rydberg units
 RETURN
END SUBROUTINE pwldafunc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to calculate the exchange correlation functionals
!   using the form of Perdew, Burke, and Ernzerhof (PRL 77, 3865 (1996))
!   assuming no spin polarization
!  Inside this routine, energies are in Hartree units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE pbefunc(den,grad,fxc,dfxcdn,dfxcdgbg,mu,beta)
 REAL(dp), INTENT(IN) :: den,grad,mu,beta    !density, magnitude of grad(density)
 REAL(dp), INTENT(OUT) :: fxc,dfxcdn,dfxcdgbg
 REAL(dp) :: n,g,kf,rs,ks,s,t,betabygamm
 REAL(dp) :: ex,ec,Fx,H,A,pprs,ppt,At2,dFds,dHdt,decdrs,dHdrs,dHdA,dAdrs
 REAL(dp) :: term,dHdtbg,dFdsbg
 betabygamm=beta/gamm
 n=den
 IF (n < machine_zero)  THEN
   fxc=0.0_dp; dfxcdn=0.0_dp; dfxcdgbg=0.0_dp
   RETURN
 ENDIF
 g=grad
 IF (g < machine_zero) g=machine_zero
 kf=(3.0_dp*(pi**2)*n)**0.3333333333333333333333333333_dp
 rs=(3.0_dp/(4.0_dp*pi*n))**0.3333333333333333333333333333_dp
 ks=SQRT(4.0_dp*kf/pi)
 s=g/(2.0_dp*kf*n)
 t=g/(2.0_dp*ks*n)
 IF (s*s > machine_infinity .or. t*t > machine_infinity)  THEN
   fxc=0.0_dp; dfxcdn=0.0_dp; dfxcdgbg=0.0_dp
   RETURN
 ENDIF
 ex=-3.0_dp*kf/(4.0_dp*pi)
 pprs=SQRT(rs)*(b1+b3*rs)+rs*(b2+b4*rs)
 term=Logofterm(1.0_dp/(2.0_dp*AA*pprs))
 ec=-2.0_dp*AA*(1.0_dp+a1*rs)*term
 Fx=1.0_dp+kappa -kappa/(1.0_dp+(mu/kappa)*s*s)
 A=Aofec(ec,betabygamm,beta)
 At2=A*t*t
 ppt=(1.0_dp+At2*(1.0_dp+At2))
 H=gamm*Logofterm((betabygamm)*(t*t)*((1._dp+At2)/ppt))
 fxc=n*(ex*Fx+ec+H)
 dFds = (2.0_dp*mu*s)/(1.0_dp+(mu/kappa)*(s**2))**2
 dFdsbg = ((2.0_dp*mu)/(1.0_dp+(mu/kappa)*(s**2))**2)/(2._dp*kf*n)
 dHdt = (2._dp*t*beta*gamm*(1._dp+2._dp*At2))/&
&     ((gamm*ppt+beta*t*t*(1._dp+At2))*ppt)
 dHdtbg = ((2._dp*beta*gamm*(1._dp+ &
&     2._dp*At2))/((gamm*ppt+beta*t*t*(1._dp+At2))*ppt))/(2._dp*ks*n)
 decdrs=-(2._dp*AA*a1)*term &
&     +((1._dp+a1*rs)*((b1+3*b3*rs)/(2._dp*SQRT(rs))+ &
&       b2+2*b4*rs))/(pprs*(pprs+1._dp/(2._dp*AA)))
 dHdA=((2._dp+At2)*(At2*t*t*t*t*beta*gamm))/&
&     ((gamm*ppt+beta*t*t*(1._dp+At2))*ppt)
 dAdrs=-ddexp(-ec/gamm)*A*A*decdrs/beta
 dHdrs=dHdA*dAdrs
 dfxcdn = (4._dp/3._dp)*ex*(Fx-dFds*s)+ec-(decdrs*rs)/3._dp+H-(dHdrs*rs)/3._dp &
&     - (7._dp/6._dp)*dHdt*t
 dfxcdgbg = ex*dFdsbg/(2._dp*kf) + dHdtbg/(2._dp*ks)
 IF ((ABS(fxc).GT.1.d65).OR.(ABS(dfxcdn).GT.1.d65).OR.&
&          (ABS(dfxcdgbg).GT.1.d65)) THEN
    if(has_to_print) WRITE(STD_OUT,*) 'Problem in PBE',n,g,rs,s,t,ec,A,H
 ENDIF
 RETURN
END SUBROUTINE pbefunc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Function Aofec -- needed to take care of behavior for small ec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Aofec(ec,betabygamm,beta)
 REAL(dp) :: ec, Aofec, betabygamm,beta
 IF (ABS(ec)>machine_precision) THEN
   Aofec=betabygamm/(ddexp(-ec/gamm)-1.0_dp)
 ELSEIF (ABS(ec)>machine_zero) THEN
   Aofec=beta/(-ec)
 ELSE
   Aofec=-beta*DSIGN(machine_infinity,ec)
 ENDIF
 RETURN
END FUNCTION Aofec




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5. r2scan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  r2scaninit
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE r2scaninit(etain)
  REAL(dp), INTENT(IN) :: etain
  REAL(dp) :: o
  if(has_to_print) write(std_out,*) 'r2scan calculation with eta ', etain
  eta=etain
  o=cx1+two*cx2+three*cx3+four*cx4+five*cx5+six*cx6+seven*cx7
  if(has_to_print) write(std_out,*) 'r2scan check C2' , o*k0
  C2Ceta = (20.0_dp/27.0_dp + (5*eta)/three)*o*k0
END SUBROUTINE r2scaninit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! r2scan functional
!!  length units == Bohr
!!  energy units -- Hartree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE r2scanfun(rho,grad,tau,exc,vtau,vxcn,vxcs)
  REAL(dp), INTENT(IN) :: rho,grad, tau
  REAL(dp), INTENT(INOUT) :: exc,vtau,vxcn,vxcs
  REAL(dp) :: rr,ss,tt,sigma
  REAL(dp) :: kF,ks,rs,U,s,p,t,W,baralpha
  REAL(dp) :: gx,x,h0x,h1x,fx,ffx,exarg,ddfx,vxtau,ex,vx
  REAL(dp) :: ELDA,ELDA0,fc,ddfc,w1,w0,ginfinity,H0c,H1c,ec0,ec1,y
  REAL(dp) :: beta,srs,ddLDA,difddLDA,gg,corarg,ddy,vctau
  REAL(dp) :: alphadenom,dalphadn,dalphads
  REAL(dp) :: gxp,dpdn,dgxdn,dpds,dgxds,drsdn,dt2ds,dt2dn,t2
  REAL(dp) :: h1xx,xp,dh1xdn,dh1xds,vxn,vxs
  REAL(dp) :: vcn,vcs,dLDA0dn,dH0cdn,dw0dn,dgindn,dw1dn
  REAL(dp) :: dH0cds,dginds,dbetadn,dif2ddLDA,dLDAdn
  REAL(dp) :: pp,dppds,dppdn,dyds,dydn,dUdn,dWds,dWdn
  REAL(dp) :: term1,term4,dddyds,dddydnp,dddydnw,dddydnr,dH1cds
  REAL(dp) :: stf1,stf2,stf3,stf4,stf5,stf6,stf7
  REAL(dp) :: dddydnn,dddydn,dH1cdwn,dH1cdn,dec1dn,dec1ds
  REAL(dp) :: dec0ds,dec0dn,bigdenom
  REAL(dp) :: eex,eec,vtx,vtc,vvxn,vvcn,vvxs,vvcs
  eex=0.d0;eec=0.d0;vtx=0.d0;vtc=0.d0;vvxn=0.d0;vvcn=0.d0;vvxs=0.d0;vvcs=0.d0
  if(rho<(machine_zero**0.333333333333333333333333333333333d0)) return
  sigma=grad*grad
  rr=rho;ss=sigma;tt=tau
  rr=max(machine_zero,rr)
  rr=min(machine_infinity,rr)
  ss=max(machine_zero,ss)
  ss=min(machine_infinity,ss)
  tt=max(machine_zero,tt)
  tt=min(machine_infinity,tt)
  !  Note that all derivatives with respect to sigma are also multiplied
  !      by grad
  kF=3.0936677262801359310d0*(rr**0.33333333333333333333333333333d0)
  ks=1.9846863952198559283d0*(rr**0.16666666666666666666666666667d0)
  rs=0.62035049089940001665d0/(rr**0.33333333333333333333333333333d0)
  drsdn=-rs/(3*rr)
  srs=sqrt(rs)
  U=2.8712340001881918160d0*(rr**1.66666666666666666666666666667d0)
  dUdn=(1.666666666666666666666667d0)*2.8712340001881918160d0*(rr**0.66666666666666666666666666667d0)
  s=sqrt(ss)*0.16162045967399548133d0/(rr**1.33333333333333333333333333333d0)
  p=s*s
  dpdn=-2.666666666666666666666666667d0*p/rr
  dpds=0.026121172985233599568d0*grad/(rr**2.666666666666666666666666666666667d0)
  t=sqrt(ss)*0.25192897034224488769d0/(rr**1.166666666666666666666666666667d0)
  t2=t*t
  dt2ds=0.063468206097703704205d0*grad/(rr**2.3333333333333333333333333333d0)
  dt2dn=-t2*(2.333333333333333333333333333333333333d0)/rr
  W=0.125d0*ss/rr
  dWds=0.125d0*grad/rr
  dWdn=-W/rr
  baralpha=(tt-W)/(U+eta*W)
  !!!!!!exchange part
  gx=1.d0-ddexp(-SCANa1/(p**0.25d0))
  x=(c2ceta*ddexp(-(p**2)/dp2**4)+mu)*p
  xp=c2ceta*ddexp(-(p**2)/dp2**4)*(-2*(p**2)/dp2**4+1.d0)+mu
  h0x=1+k0
  h1x=1.d0+k1-k1/(1+x/k1)
  if(baralpha<1.d-13) then
    fx=ddexp(-SCANc1x*baralpha/(1.d0-baralpha))
  elseif (baralpha.lt.2.5d0) then
    fx=cx0+baralpha*(cx1+baralpha*(cx2+baralpha*(cx3+baralpha*( &
&        cx4+baralpha*(cx5+baralpha*(cx6+baralpha*cx7))))))
  else if (baralpha.ge.2.5d0) then
    fx=-SCANdx*ddexp(SCANc2x/(1.d0-baralpha))
  endif
  if(baralpha<1.d-13) then
    ddfx=-(SCANc1x/((1.d0-baralpha)**2))*ddexp(-SCANc1x*baralpha/(1.d0-baralpha))
  else if (baralpha.lt.2.5d0) then
    ddfx=(cx1+baralpha*(2*cx2+baralpha*(3*cx3+baralpha*( &
&        4*cx4+baralpha*(5*cx5+baralpha*(6*cx6+baralpha*7*cx7))))))
  else if (baralpha.ge.2.5d0) then
    ddfx=-(SCANdx*SCANc2x/((1.d0-baralpha)**2))*ddexp(SCANc2x/(1.d0-baralpha))
  endif
  ffx=(h1x+fx*(h0x-h1x))*gx  
  ex=-0.73855876638202240587d0*(rr**1.3333333333333333333333333333333333d0)
  vx=-0.98474502184269654116d0*(rr**0.3333333333333333333333333333333333d0)
  exarg=ex*ffx
  vxtau=ex*ddfx*((h0x-h1x)*gx)/(U+eta*W)          
  ! density and sigma derivative terms
  alphadenom=(U+eta*W)**2
  dalphads=-(U+eta*tau)*dWds/alphadenom
  dalphadn=-((U+eta*tau)*dWdn+(tt-W)*dUdn)/alphadenom
  gxp=-0.25d0*SCANa1*ddexp(-SCANa1/(p**0.25d0))/(p**1.25d0)
  dgxdn=gxp*dpdn
  dgxds=gxp*dpds
  h1xx=1.d0/((1.d0+x/k1)**2) 
  dh1xdn=h1xx*xp*dpdn
  dh1xds=h1xx*xp*dpds
  vxn=vx*ffx+ex*(ddfx*dalphadn*gx*(h0x-h1x)+dgxdn*(h1x+fx*(h0x-h1x)) &
&      +gx*dh1xdn*(1.d0-fx))          
  vxs=ex*(ddfx*dalphads*gx*(h0x-h1x)+dgxds*(h1x+fx*(h0x-h1x)) &
&      +gx*dh1xds*(1.d0-fx))          
  !!!! correlation part
  ELDA=-two*LDAA*(1.d0 + LDAa1*rs)*ddlog(1.d0 + 0.5d0 &
&   /(LDAA*(sqrt(rs)*(LDAb1 + LDAb3*rs) + rs*(LDAb2 + LDAb4*rs))))
  bigdenom=srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs + LDAb2)
  dLDAdn=((-two*LDAA*LDAa1*ddlog(1.d0 + &
&    1.d0/(two*LDAA*bigdenom))) &
&    + (LDAa1*rs + 1.d0)*((LDAb3*rs + LDAb1)/(two*srs) + srs*LDAb3 +  &
&    two*LDAb4*rs + LDAb2)/((bigdenom**2)+bigdenom/(two*LDAA)))*drsdn
  ELDA0= -b1c/(1.d0 + b2c*sqrt(rs) + b3c*rs)
  dLDA0dn=b1c*(0.5d0*b2c/sqrt(rs)+b3c)/((1.d0 + b2c*sqrt(rs) + b3c*rs)**2) &
&    *drsdn
  ddLDA=ELDA0-ELDA
  beta= betaMB*(1.d0 + 0.1*rs)/(1.d0 + 0.1778d0*rs)
  dbetadn=-(0.0778d0*betaMB/(1.d0 + 0.1778d0*rs)**2)*drsdn
  if(baralpha<1.d-13) then
    fc=ddexp(-SCANc1c*baralpha/(1.d0-baralpha))
  else if (baralpha.lt.2.5d0) then
    fc=cc0+baralpha*(cc1+baralpha*(cc2+baralpha*(cc3+baralpha*( &
&        cc4+baralpha*(cc5+baralpha*(cc6+baralpha*cc7))))))
  else if (baralpha.ge.2.5d0) then
    fc=-SCANdc*ddexp(SCANc2c/(1.d0-baralpha))
  endif
  if(baralpha<1.d-13) then
    ddfc=-(SCANc1c/((1.d0-baralpha)**2))*ddexp(-SCANc1c*baralpha/(1.d0-baralpha))
  else if (baralpha.lt.2.5d0) then
    ddfc=(cc1+baralpha*(2*cc2+baralpha*(3*cc3+baralpha*( &
&        4*cc4+baralpha*(5*cc5+baralpha*(6*cc6+baralpha*7*cc7))))))
  else if (baralpha.ge.2.5d0) then
    ddfc=-(SCANdc*SCANc2c/((1.d0-baralpha)**2))*ddexp(SCANc2c/(1.d0-baralpha))
  endif
  w1= ddexp(-ELDA/Sgam) - 1.d0
  dw1dn=-(ddexp(-ELDA/Sgam)/Sgam)*dLDAdn
  w0= ddexp(-ELDA0/b1c) - 1.d0
  dw0dn=-(ddexp(-ELDA0/b1c)/b1c)*dLDA0dn
  ginfinity= 1.d0/(1.d0 + 4*chiinfinity*(p))**0.25d0
  dgindn=-(chiinfinity*ginfinity/(1.d0 + four*chiinfinity*(p)))*dpdn
  dginds=-(chiinfinity*ginfinity/(1.d0 + four*chiinfinity*(p)))*dpds
  H0c=  b1c*ddlog(1.d0 + w0*(1.d0 - ginfinity))
  dH0cdn=(b1c/(1.d0+w0*(1.d0-ginfinity)))*((1.d0-ginfinity)*dw0dn-w0*dgindn)
  dH0cds=(-b1c*w0*dginds/(1.d0+w0*(1.d0-ginfinity)))
  ec0=  ELDA0 + H0c
  dec0dn=dLDA0dn+dH0cdn
  dec0ds=dH0cds
  difddLDA=b1c*(b2c/(two*srs) + b3c)/(1.d0 + b2c*srs + b3c*rs)**2 + &
&  two*LDAA*LDAa1*ddlog(1.d0 + &
&    1.d0/(2*LDAA*(srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs +  LDAb2)))) &
&    - (LDAa1*rs + 1.d0)*((LDAb3*rs + LDAb1)/(two*srs) + srs*LDAb3 +  &
&    two*LDAb4*rs + LDAb2)/(((srs*(LDAb3*rs + LDAb1) + &  
&    rs*(LDAb4*rs + LDAb2))**2)*(1.d0       +    1.d0 &
&  /(two*LDAA*(srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs + LDAb2)))))
  y=beta*(t2)/(Sgam*w1)
  dyds=beta*dt2ds/(Sgam*w1)
  dydn=(1.d0/(Sgam*w1))*(dbetadn*t2+beta*dt2dn-beta*t2*dw1dn/w1)
  ddy=(Dfc2/(27.0_dp*Sgam*w1))*(20.0_dp*rs*difddLDA-45.0_dp*eta*ddLDA)*p&
&     *ddexp(-(p**2)/(dp2**4))
  gg=1.d0/(1.d0+4*(y-ddy))**0.25d0
  H1c=Sgam*ddlog(1.d0+w1*(1.d0-gg))
  ec1=ELDA+H1c
  pp=p*ddexp(-(p**2)/(dp2**4))
  dppds=((ddexp(-p**2/dp2**4))*(dp2**4 - two*p**2)/dp2**4)*dpds
  dppdn=((ddexp(-p**2/dp2**4))*(dp2**4 - two*p**2)/dp2**4)*dpdn
  term1=1.d0+4*(y-ddy)
  term4=term1**0.25d0
  dddyds=(Dfc2/(27.0_dp*Sgam*w1))*(20.0_dp*rs*difddLDA-45.0_dp*eta*ddLDA)*dppds
  dddydnp=(Dfc2/(27.0_dp*Sgam*w1))*(20.0_dp*rs*difddLDA-45.0_dp*eta*ddLDA)*dppdn
  dddydnw=-(Dfc2/(27.0_dp*Sgam*w1))*(20.0_dp*rs*difddLDA-45.0_dp*eta*ddLDA)*pp*dw1dn/w1
  dddydnr=(Dfc2/(27.0_dp*Sgam*w1))*(20.0_dp*difddLDA)*pp*drsdn
  dH1cds=(Sgam*w1/(term1*(term4*(1.d0+w1)-w1)))*(dyds-dddyds)
  stf1=-b1c*(8.0_dp*srs*rs*(b3c**2) + 9.0_dp*b2c*b3c*rs + 3.0_dp*srs*(b2c**2) + b2c)&
&          /(4.0_dp*rs*srs*((1.d0 + b2c*srs + b3c*rs)**3))
  stf2=srs*(LDAb3*rs + LDAb1) + rs*(LDAb4*rs + LDAb2)
  stf3=1.d0+1.d0/(2*LDAA*stf2)
  stf4=(two*LDAa1*((LDAb3*rs+LDAb1)/(two*srs)+srs*LDAb3+two*LDAb4*rs+LDAb2)) &
&      /((stf2**2)*stf3)          
  stf5=(two*(LDAa1*rs+1.d0)*(((LDAb3*rs+LDAb1)/(two*srs)+srs*LDAb3+  &
&     two*LDAb4*rs + LDAb2)**2))/((stf2**3)*stf3)
  stf6=(LDAa1*rs+1.d0)*(-(LDAb3*rs+LDAb1)/(four*srs*rs)+LDAb3/srs + two*LDAb4) &
&   /((stf2**2)*stf3)          
  stf7=(LDAa1*rs+1.d0)*& 
&  (((LDAb3*rs+LDAb1)/(two*srs)+srs*LDAb3+two*LDAb4*rs+LDAb2)**2) &
&   /(two*LDAA*(stf2**4)*(stf3**2))
  dif2ddLDA=stf1-stf4+stf5-stf6-stf7
  dddydnn=(Dfc2/(27.0_dp*Sgam*w1))*(20.0_dp*rs*dif2ddLDA-45.0_dp*eta*difddLDA)*pp*drsdn
  dddydn=dddydnw+dddydnp+dddydnr+dddydnn
  dH1cdwn=Sgam*((term4-1.d0)/(term4*(w1+1.d0)-w1))*dw1dn
  dH1cdn=dH1cdwn+(Sgam*w1/(term1*(term4*(1.d0+w1)-w1)))*(dydn-dddydn)
  dec1dn=dLDAdn+dH1cdn
  dec1ds=dH1cds 
  corarg=rr*(ec1 + fc*(ec0 - ec1))
  vctau=rr*ddfc*(ec0 - ec1)/(U + eta*W)
  ! density and sigma derivative terms
  vcn=ec1+fc*(ec0-ec1)+ddfc*dalphadn*rr*(ec0-ec1)+dec0dn*rr*fc+dec1dn*rr*(1.d0-fc)
  vcs=ddfc*dalphads*rr*(ec0-ec1)+dec0ds*rr*fc+dec1ds*rr*(1.d0-fc)
  eex=exarg;eec=corarg;               exc=exarg+corarg
  vtx=vxtau;vtc=vctau;                vtau=vxtau+vctau
  vvxn=vxn; vvcn=vcn;                 vxcn=vxn+vcn
  vvxs=vxs; vvcs=vcs;                 vxcs=vxs+vcs
END SUBROUTINE r2scanfun




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5. general_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Updatewfn(Grid,Pot,Orbit,rvin,success)
!!   Given new potential rvin, generate new Orbit%wfn,Orbit%eig,Pot%den
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Updatewfn(Grid,Pot,Orbit,rvin,success,BDsolve,usespline,spline,itype)
! TODO : BDsolve
 logical,intent(in) :: usespline
 integer, intent(in) :: itype
 TYPE (GridInfo), INTENT(INOUT) :: Grid
 TYPE (PotentialInfo), INTENT(INOUT) :: Pot
 TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
 type(splinesolvinfo),intent(inout) :: spline
 logical, intent(in) :: bdsolve
 real(dp), INTENT(IN) ::  rvin(:)
 LOGICAL :: success,calc_s,calc_p,calc_d,calc_f,calc_g
 INTEGER :: icount,n,it,start,ierr,nroot,s1,s2,s2t
 INTEGER :: io,l,jierr,nz,kappa
 real(dp) :: h,emin,zz
 real(dp), ALLOCATABLE :: dum(:)
 LOGICAL :: OK
 n=Grid%n; h=Grid%h;    nz=Pot%nz;   zz=Pot%zz
 success=.TRUE.
 LIBPAW_ALLOCATE(dum,(n))
 Pot%rv=rvin+Pot%rvn(1:n)
 dum=rvin
 if (Pot%needvtau) dum(1)=dum(1)-Pot%rvx(1)
 CALL zeropot(Grid,dum,Pot%v0,Pot%v0p)
 IF (ABS(Pot%v0)> 1.d6) Pot%v0=0
 IF (ABS(Pot%v0p)> 1.d6) Pot%v0p=0
 IF (Pot%finitenucleus) then
         Pot%v0=Pot%v0+Pot%Nv0
         Pot%v0p=Pot%v0p+Pot%Nv0p
 Endif
 if(usespline) call initpotforsplinesolver(Grid,Pot,Orbit%den,Orbit%tau,spline,itype)
 !  solve for bound states of Schroedinger equation
 calc_s=.true.
 calc_p=.true.
 calc_d=.true.
 calc_f=.true.
 calc_g=.true.
 icount=0
 jierr=0
 it=0
 !  s states :
 IF (Orbit%nps.GT.0) THEN
   it=it+1
   emin=-nz*nz-0.1_dp
   l=0
   nroot=Orbit%nps
   start=1;s1=start;s2t=start+nroot-1
   if (Orbit%frozenvalecalculation) then
     nroot=Orbit%npsc
     do io=s1+Orbit%npsc,s2t
       if(Orbit%issemicore(io)) nroot=nroot+1
     enddo
     if(nroot.LT.1) calc_s=.false.
   endif
   s2=s1+nroot-1
   if(calc_s) then       
     IF (Orbit%scalarrelativistic) THEN
       Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&               l,nroot,emin,ierr,OK)
     ELSE IF (Orbit%diracrelativistic) THEN
       kappa=-1     
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&          Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
     ELSE IF (Pot%needvtau) THEN
        Call Boundsplinesolver(Grid,l,nroot, &
&              Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
     ELSE
       if(usespline) then
         Call Boundsplinesolver(Grid,l,nroot, &
&            Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
       else
         CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&                l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       endif
     ENDIF
     IF (.NOT.OK) THEN
         success=.FALSE.
     ENDIF
   endif
 ENDIF
    !  p states :
 IF (Orbit%npp.GT.1) then
   it=it+1
   emin=-nz*nz/4._dp-0.5_dp
   l=1
   nroot=Orbit%npp-1
   s1=s2t+1;s2t=s1+nroot-1
   if (Orbit%frozenvalecalculation) then 
     nroot=Orbit%nppc-1
     do io=s1+Orbit%nppc-1,s2t
       if(Orbit%issemicore(io)) nroot=nroot+1
     enddo
     if(nroot.LT.1) calc_p=.false.
   endif
   s2=s1+nroot-1
   if(calc_p) then
     IF (Orbit%scalarrelativistic) THEN
       Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&           l,nroot,emin,ierr,OK)
     ELSE IF (Orbit%diracrelativistic) THEN
       kappa=1     
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       s1=s2+1;s2=s1+nroot-1
       kappa=-2     
       emin=-nz*nz/4._dp-0.5_dp
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&          Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
     ELSE IF (Pot%needvtau) THEN
       Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
     ELSE
       if(usespline) then
         Call Boundsplinesolver(Grid,l,nroot, &
&            Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
       else
         CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&             l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       endif
     ENDIF
     IF (.NOT.OK) THEN
       success=.FALSE.
     ENDIF
   endif
 ENDIF
 !  d states :
 IF (Orbit%npd.GT.2) THEN
   it=it+1
   emin=-nz*nz/9._dp-0.5_dp
   l=2
   nroot=Orbit%npd-2
   s1=s2t+1;s2t=s1+nroot-1
   if (Orbit%frozenvalecalculation) then
     nroot=Orbit%npdc-2
     do io=s1+Orbit%npdc-2,s2t
       if(Orbit%issemicore(io)) nroot=nroot+1
     enddo
     if(nroot.LT.1) calc_d=.false.
   endif
   s2=s1+nroot-1
   if(calc_d) then
     IF (Orbit%scalarrelativistic) THEN
       Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&              l,nroot,emin,ierr,OK)
     ELSE IF (Orbit%diracrelativistic) THEN
       kappa=2     
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&          Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       kappa=-3     
       s1=s2+1;s2=s1+nroot-1
       emin=-nz*nz/9._dp-0.5_dp
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&          Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
     ELSE IF (Pot%needvtau) THEN
       Call Boundsplinesolver(Grid,l,nroot, &
&          Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
     ELSE!          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
       if(usespline) then
         Call Boundsplinesolver(Grid,l,nroot, &
&            Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
       else
         CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&               l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       endif
     ENDIF
     IF (.NOT.OK) THEN
       success=.FALSE.
     ENDIF
   endif
 ENDIF
 !  f states :
 IF (Orbit%npf.GT.3) THEN
   it=it+1
   emin=-nz*nz/16._dp-0.5_dp
   l=3
   nroot=Orbit%npf-3
   s1=s2t+1;s2t=s1+nroot-1
   if (Orbit%frozenvalecalculation) then
     nroot=Orbit%npfc-3
     do io=s1+Orbit%npfc-3,s2t
       if(Orbit%issemicore(io)) nroot=nroot+1
     enddo
     if(nroot.LT.1) calc_f=.false.
   endif
   s2=s1+nroot-1
   if (calc_f) then
     IF (Orbit%scalarrelativistic) THEN
       Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&            l,nroot,emin,ierr,OK)
     ELSE IF (Orbit%diracrelativistic) THEN
       kappa=3     
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       kappa=-4     
       s1=s2+1;s2=s1+nroot-1
       emin=-nz*nz/16._dp-0.5_dp
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
     ELSE IF (Pot%needvtau) THEN
       Call Boundsplinesolver(Grid,l,nroot, &
&            Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
     ELSE
       if(usespline) then
         Call Boundsplinesolver(Grid,l,nroot, &
&            Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
       else
         CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       endif
     ENDIF
     IF (.NOT.OK) THEN
       success=.FALSE.
     ENDIF
   endif
 ENDIF
 !  g states :
 IF (Orbit%npg.GT.4) THEN
   it=it+1
   emin=-nz*nz/25._dp-0.5_dp
   l=4
   nroot=Orbit%npg-4
   s1=s2t+1;s2t=s1+nroot-1
   if (Orbit%frozenvalecalculation) then
     nroot=Orbit%npgc-4
     do io=s1+Orbit%npgc-4,s2t
       if(Orbit%issemicore(io)) nroot=nroot+1
     enddo
     if(nroot.LT.1) calc_g=.false.
   endif
   s2=s1+nroot-1
   if(calc_g) then
     IF (Orbit%scalarrelativistic) THEN
       Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&            l,nroot,emin,ierr,OK)
     ELSE IF (Orbit%diracrelativistic) THEN
       kappa=4     
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
          kappa=-5     
       s1=s2+1;s2=s1+nroot-1
       emin=-nz*nz/25._dp-0.5_dp
       Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
     ELSE IF (Pot%needvtau) THEN
       Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
     ELSE
       if(usespline) then
         Call Boundsplinesolver(Grid,l,nroot, &
&            Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
       else
         CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       endif
     ENDIF
     IF (.NOT.OK) THEN
       success=.FALSE.
     ENDIF
   endif
 ENDIF
 !Update otau according to wfn
 DO io=1,Orbit%norbit
   IF(Orbit%frozenvalecalculation.and.(.not.Orbit%iscore(io))) cycle
   CALL taufromwfn(Orbit%otau(:,io),Grid,Orbit%wfn(:,io),Orbit%l(io), &
&                     energy=Orbit%eig(io),rPot=Pot%rv)
 ENDDO
 LIBPAW_DEALLOCATE(dum)
END SUBROUTINE Updatewfn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE Updatewfnwden(Grid,Pot,Orbit,w,success)
!!      Given new den and consistent Pot%rv
!!      generate new Orbit%wfn,Orbit%eig,Orbit%otau
!!      Orbit%den=w on input
!!     only splinesolver works for now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Updatewfnwden(Grid,Pot,Orbit,success,usespline,spline,itype)
  integer, intent(in) :: itype
  type(splinesolvinfo),intent(inout) :: spline
  TYPE (GridInfo), INTENT(INOUT) :: Grid
  TYPE (PotentialInfo), INTENT(INOUT) :: Pot
  TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
  LOGICAL,intent(inout) :: success
  LOGICAL,intent(in) :: usespline
  INTEGER :: icount,n,it,start,nroot,s1,s2,s2t
  INTEGER :: io,l,jierr,nz
  LOGICAL :: calc_s,calc_p,calc_d,calc_f,calc_g
  REAL(dp) :: h,emin,zz
  REAL(dp), ALLOCATABLE :: dum(:)
  LOGICAL :: OK
  !!!since programmed only for splinesolver case, check
  If (Orbit%scalarrelativistic.or.Orbit%diracrelativistic) then
      write(std_out,*) 'Program Updatewfnwden not written '    
      write(std_out,*) ' for relativistic solver -- error stop '
      stop
  ENDIF    
  If (.not.usespline) THEN
      write(std_out,*) 'Program Updatewfnwden only written '    
      write(std_out,*) ' for usespline case -- error stop '
      stop
  ENDIF    
  n=Grid%n; h=Grid%h;    nz=Pot%nz;   zz=Pot%zz
  success=.TRUE.
  LIBPAW_ALLOCATE(dum,(n))
  call initpotforsplinesolver(Grid,Pot,Orbit%den,Orbit%tau,spline,itype)
  !  solve for bound states of Schroedinger equation
  icount=0
  jierr=0
  it=0
  !  s states :
  IF (Orbit%nps.GT.0) THEN
    it=it+1
    emin=-nz*nz-0.1_dp
    l=0
    nroot=Orbit%nps
    start=1;s1=start;s2t=start+nroot-1
    if (Orbit%frozenvalecalculation) then
      nroot=Orbit%npsc
      do io=s1+Orbit%npsc,s2t
        if(Orbit%issemicore(io)) nroot=nroot+1
      enddo
      if(nroot.LT.1) calc_s=.false.
    endif
    s2=s1+nroot-1
    if(calc_s) then
      Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
      IF (.NOT.OK) THEN
          success=.FALSE.
      ENDIF
    endif
  ENDIF
      !  p states :
  IF (Orbit%npp.GT.1) then
    it=it+1
    emin=-nz*nz/4._dp-0.5_dp
    l=1
    nroot=Orbit%npp-1
    s1=s2t+1;s2t=s1+nroot-1
    if (Orbit%frozenvalecalculation) then
      nroot=Orbit%nppc-1
      do io=s1+Orbit%nppc-1,s2t
        if(Orbit%issemicore(io)) nroot=nroot+1
      enddo
      if(nroot.LT.1) calc_p=.false.
    endif
    s2=s1+nroot-1
    if(calc_p) then
      Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
      IF (.NOT.OK) THEN
        success=.FALSE.
      ENDIF
    endif
  ENDIF
   !  d states :
  IF (Orbit%npd.GT.2) THEN
    it=it+1
    emin=-nz*nz/9._dp-0.5_dp
    l=2
    nroot=Orbit%npd-2
    s1=s2t+1;s2t=s1+nroot-1
    if (Orbit%frozenvalecalculation) then
      nroot=Orbit%npdc-2
      do io=s1+Orbit%npdc-2,s2t
        if(Orbit%issemicore(io)) nroot=nroot+1
      enddo
      if(nroot.LT.1) calc_d=.false.
    endif
    s2=s1+nroot-1
    if(calc_d) then
      Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
      IF (.NOT.OK) THEN
        success=.FALSE.
      ENDIF
    endif
  ENDIF
   !  f states :
  IF (Orbit%npf.GT.3) THEN
    it=it+1
    emin=-nz*nz/16._dp-0.5_dp
    l=3
    nroot=Orbit%npf-3
    s1=s2t+1;s2t=s1+nroot-1
    if (Orbit%frozenvalecalculation) then
      nroot=Orbit%npfc-3
      do io=s1+Orbit%npfc-3,s2t
        if(Orbit%issemicore(io)) nroot=nroot+1
      enddo
      if(nroot.LT.1) calc_f=.false.
    endif
    s2=s1+nroot-1
    if (calc_f) then
      Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
      IF (.NOT.OK) THEN
        success=.FALSE.
      ENDIF
    endif
  ENDIF
   !  g states :
  IF (Orbit%npg.GT.4) THEN
    it=it+1
    emin=-nz*nz/25._dp-0.5_dp
    l=4
    nroot=Orbit%npg-4
    s1=s2t+1;s2t=s1+nroot-1
    if (Orbit%frozenvalecalculation) then
      nroot=Orbit%npgc-4
      do io=s1+Orbit%npgc-4,s2t
        if(Orbit%issemicore(io)) nroot=nroot+1
      enddo
      if(nroot.LT.1) calc_g=.false.
    endif
    s2=s1+nroot-1
    if(calc_g) then
      Call Boundsplinesolver(Grid,l,nroot, &
&         Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),Orbit%otau(:,s1:s2),OK,spline)
      IF (.NOT.OK) THEN
        success=.FALSE.
      ENDIF
    endif
  ENDIF
  LIBPAW_DEALLOCATE(dum)
END SUBROUTINE Updatewfnwden



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Get_KinCoul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_KinCoul(Grid,Pot,Orbit,SCF,usespline,noalt)
 !  program to calculate Kinetic energy and Coulomb Energies from Orbit%wfn
 !   also update Pot%rvh
 logical,intent(in) :: usespline
 TYPE(GridInfo), INTENT(INOUT) :: Grid
 TYPE(PotentialInfo), INTENT(INOUT) :: Pot
 TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
 TYPE(SCFInfo), INTENT(INOUT) :: SCF
 LOGICAL, OPTIONAL :: noalt
 real(dp) :: ecoul,ekin,eone,h,x,qcal,small,rescale
 real(dp) :: electrons,xocc
 INTEGER :: i,n,io
 real(dp), ALLOCATABLE :: dum(:)
 real(dp) :: small0=tol6,fpi
 INTEGER :: counter=1
 n=Grid%n; h=Grid%h
 small=small0
 LIBPAW_ALLOCATE(dum,(n))
 !update density
 Orbit%den=0._dp;Orbit%tau=0._dp
 DO io=1,Orbit%norbit
   IF(Orbit%frozenvalecalculation.and.(.not.Orbit%iscore(io))) cycle
   IF (Orbit%occ(io).GT.small) THEN
     DO i=1,Grid%n
       IF (ABS(Orbit%wfn(i,io))<machine_zero)Orbit%wfn(i,io)=0
       IF (Orbit%diracrelativistic) then
         IF (ABS(Orbit%lwfn(i,io))<machine_zero)Orbit%lwfn(i,io)=0
       ENDIF
     ENDDO
     if(.not.usespline) then
       CALL taufromwfn(Orbit%otau(:,io),Grid,Orbit%wfn(:,io),Orbit%l(io), &
&                                     energy=Orbit%eig(io),rPot=Pot%rv)
     endif
     xocc=Orbit%occ(io)
     Do i=1,Grid%n
       Orbit%tau(i)=Orbit%tau(i)+xocc*Orbit%otau(i,io)
       Orbit%den(i)=Orbit%den(i)+xocc*(Orbit%wfn(i,io)**2)
       IF (Orbit%diracrelativistic) then
         Orbit%den(i)=Orbit%den(i)+xocc*((Orbit%lwfn(i,io))**2)
       ENDIF
     ENDDO
   ENDIF
 ENDDO
 qcal=integrator(Grid,Orbit%den)
 if(has_to_print) WRITE(STD_OUT,*) 'qcal = ', qcal
 IF(Orbit%frozenvalecalculation) qcal=qcal+Orbit%qval
 IF(Orbit%frozenvalecalculation) Orbit%den=Orbit%den+Orbit%valeden
 electrons=Pot%q
 IF(Orbit%frozenvalecalculation) electrons=qcal 
 rescale=electrons/qcal
 Orbit%den(1:n)=Orbit%den(1:n)*rescale
 Orbit%tau(1:n)=Orbit%tau(1:n)*rescale
 !   Determine difference with tauW (Weizsaker)       
 fpi=4*pi
 dum(2:Grid%n)=Orbit%den(2:Grid%n)/(fpi*Grid%r(2:Grid%n)**2)
 CALL extrapolate(dum)
 CALL derivative(Grid,dum,Orbit%deltatau)
 Do i=1,Grid%n
   if (dum(i)>machine_zero) then
     Orbit%deltatau(i)=0.25_dp*(Orbit%deltatau(i)**2)/dum(i)
   else
     Orbit%deltatau(i)=0.0_dp
   endif
 enddo        
 dum(2:Grid%n)=Orbit%tau(2:Grid%n)/(fpi*Grid%r(2:Grid%n)**2)
 call extrapolate(dum)
 Orbit%deltatau=dum-Orbit%deltatau
 call poisson_marc(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul)   
! call atompaw_poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul)
 dum=zero
 dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/Grid%r(2:n)
 SCF%estatic=integrator(Grid,dum)+ecoul
 if(has_to_print) then
    WRITE(STD_OUT,*) ' n  l     occupancy       energy'
    do io=1,Orbit%norbit
      write(std_out,*) Orbit%np(io), Orbit%l(io),Orbit%occ(io),Orbit%eig(io)
    enddo
 endif
 ekin=0.0_dp; if (Orbit%frozencorecalculation) ekin=SCF%corekin
 eone=0.0_dp
 DO io=1,Orbit%norbit
   if(.not.Orbit%frozencorecalculation &
&        .or.Orbit%frozencorecalculation.and.(.not.Orbit%iscore(io))) then
     eone=eone+Orbit%occ(io)*Orbit%eig(io)
     IF (counter>1.and..not.present(noalt)) THEN
       CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
     ELSE
       x=integrator(Grid,Orbit%otau(:,io))
     ENDIF
     ekin=ekin+Orbit%occ(io)*x
   endif
 ENDDO
 if(has_to_print) write(std_out,*) 'KinCoul check ekin ',ekin,integrator(Grid,Orbit%tau)
 SCF%eone=eone
 SCF%ekin=ekin
 SCF%ecoul=ecoul
 counter=counter+1
 LIBPAW_DEALLOCATE(dum)
END SUBROUTINE Get_KinCoul


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Get_Nuclearpotential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Get_Nuclearpotential(Grid,Pot)
! TODO : finitenucleus
 TYPE(GridInfo), INTENT(INOUT) :: Grid
 TYPE(PotentialInfo), INTENT(INOUT) :: Pot
!  Various finite nuclear models follow the manuscript of Andrae
!   Physics Reports 336 (2000) 413-525
!    finitenucleusmodel 2,3,4,5 correspond to the options
!     described in that paper while finitenucleusmodel 0 corresponds to
!     Gaussian model originally programmed
!     Note that logarithmic grid is reset to be compatible with
!       nuclear model with approximately NN integration points within
!       finite nucleus
 INTEGER :: i
 INTEGER, PARAMETER :: NN=651    ! number of grid points within RR
 real(dp), PARAMETER :: gridrange=100._dp
 real(dp), PARAMETER :: bohr=0.529177249_dp  !Ang/Bohr from Andrae
 IF (.NOT.Pot%finitenucleus) THEN
 !  grid already set
   DO i=1,Grid%n
     Pot%rvn(i)=-2*Pot%nz!*0.5!Hartree
   ENDDO
 ELSE
   STOP
 !  write(std_out,*) 'Finite nucleus model  -- readjusting integration grid'
 !  a=bohr*10._dp**(-5)*(0.57_dp+0.836*   &
 !   &     (-1.168_dp+Pot%nz*(2.163_dp+Pot%nz*0.004467_dp)))
 !      !  From Eqs. A.3 and 51 in Andrae paper
 !  write(std_out,*) 'a parameter calculated to be', a
 !  call destroygrid(Grid)
 !  SELECT CASE(Pot%finitenucleusmodel)
 !     CASE DEFAULT
 !       write(std_out,*) 'Error in finitenucleusmodel',Pot%finitenucleusmodel
 !       write(std_out,*) ' Exiting '
 !       Stop
 !     CASE(0)
 !       write(std_out,*) 'Original Gaussian model'
 !       RR=Pot%nz
 !       RR=2.9*10._dp**(-5)*(RR**0.3333333333333333333333333_dp)
 !       h=log(FLOAT(NN))/(NN-1)
 !       r0=RR/(NN-1)
 !       write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
 !       Call InitGrid(Grid,h,gridrange,r0=r0)
 !       write(std_out,*) 'New Grid ', Grid%n
 !       Call DestroyPot(Pot)
 !       Call InitPot(Pot,Grid%n)
 !       DO i=1,Grid%n
 !         Pot%rvn(i)=-2*Pot%nz*derf(Grid%r(i)/RR)
 !       ENDDO
 !       Pot%Nv0=-2*Pot%nz*sqrt(4._dp/pi)
 !       Pot%Nv0p=0._dp
 !     CASE(2)
 !       write(std_out,*) 'Model 2 -- Breit'
 !       RR=sqrt(2._dp)*a
 !       h=log(FLOAT(NN))/(NN-1)
 !       r0=RR/(NN-1)
 !       write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
 !       Call InitGrid(Grid,h,gridrange,r0=r0)
 !       write(std_out,*) 'New Grid ', Grid%n
 !       Call DestroyPot(Pot)
 !       Call InitPot(Pot,Grid%n)
 !       DO i=1,Grid%n
 !         if (Grid%r(i)<RR) then
 !           Pot%rvn(i)=-2*Pot%nz*Grid%r(i)*(2._dp-Grid%r(i)/RR)/RR
 !         else
 !           Pot%rvn(i)=-2*Pot%nz
 !         endif
 !       ENDDO
 !       Pot%Nv0=-2*Pot%nz*2.0_dp/RR
 !       Pot%Nv0p=2*Pot%nz/(RR**2)
 !     CASE(3)
 !       write(std_out,*) 'Model 3 -- uniform'
 !       RR=sqrt(5._dp/3._dp)*a
 !       h=log(FLOAT(NN))/(NN-1)
 !       r0=RR/(NN-1)
 !       write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
 !       Call InitGrid(Grid,h,gridrange,r0=r0)
 !       write(std_out,*) 'New Grid ', Grid%n
 !       Call DestroyPot(Pot)
 !       Call InitPot(Pot,Grid%n)
 !       DO i=1,Grid%n
 !         if (Grid%r(i)<RR) then
 !           Pot%rvn(i)=-3*Pot%nz*Grid%r(i)*&
 !              &     (1._dp-(Grid%r(i)/RR)**2/3)/RR
 !         else
 !           Pot%rvn(i)=-2*Pot%nz
 !         endif
 !       ENDDO
 !         Pot%Nv0=-3*Pot%nz/RR
 !         Pot%Nv0p=0._dp
 !     CASE(4)
 !       write(std_out,*) 'Model 4 -- exponential'
 !       RR=sqrt(1._dp/12._dp)*a
 !       h=log(FLOAT(NN))/(NN-1)
 !       r0=RR/(NN-1)
 !       write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
 !       Call InitGrid(Grid,h,gridrange,r0=r0)
 !       write(std_out,*) 'New Grid ', Grid%n
 !       Call DestroyPot(Pot)
 !       Call InitPot(Pot,Grid%n)
 !       DO i=1,Grid%n
 !        Pot%rvn(i)=-2*Pot%nz*   &
 !          &  (1._dp-exp(-grid%r(i)/RR)*(1._dp+0.5_dp*Grid%r(i)/RR))
 !       ENDDO
 !       Pot%Nv0=-Pot%nz/RR
 !       Pot%Nv0p=0._dp
 !     CASE(5)
 !       write(std_out,*) 'Model 5 -- Gaussian'
 !       RR=sqrt(2._dp/3._dp)*a
 !       h=log(FLOAT(NN))/(NN-1)
 !       r0=RR/(NN-1)
 !       write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
 !       Call InitGrid(Grid,h,gridrange,r0=r0)
 !       write(std_out,*) 'New Grid ', Grid%n
 !       Call DestroyPot(Pot)
 !       Call InitPot(Pot,Grid%n)
 !       DO i=1,Grid%n
 !         Pot%rvn(i)=-2*Pot%nz*erf(Grid%r(i)/RR)
 !       ENDDO
 !       Pot%Nv0=-2*Pot%nz/RR*(sqrt(4._dp/pi))
 !       Pot%Nv0p=0._dp
 !   END SELECT
  ENDIF
END SUBROUTINE Get_Nuclearpotential





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 6. radialsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine Azeroexpand(Grid,Pot,l,energy)
!!      If finitenucleus==.true. assumes potential is non-singular
!!          at origin and Pot%v0 and Pot%v0p are properly set
!!      Otherwise, assumes nuclear potential is -2*Z/r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Azeroexpand(Grid,Pot,l,energy,qq,gamma,c1,c2,MA,MB,nr)
 Type(GridInfo), INTENT(IN) :: Grid
 Type(PotentialInfo), INTENT(INout) :: Pot
 Integer, INTENT(IN) :: l
 real(dp), INTENT(IN) :: energy
 real(dp), intent(inout) :: qq,gamma,c1,c2,MA,MB
 Integer, optional, INTENT(IN) :: nr
 Integer :: n
 real(dp) :: nz,angm,alpha2,balpha2
 real(dp) :: Tm10,Tm11,T00,Tm21,Tm22,term
 n=Grid%n
 if (present(nr)) n=min(n,nr)
! check for possible ionic charge
 n=Grid%n
 qq=-Pot%rv(n)/2
 if(qq<0.001_dp) qq=0
 qq=zero
 nz=Pot%nz
 Pot%ww=0; Pot%jj=0;
 balpha2=InvFineStruct**2
 alpha2=1._dp/balpha2
 Pot%jj(1:n)=(Grid%r(1:n) + &
&       0.25_dp*alpha2*(energy*Grid%r(1:n)-Pot%rv(1:n)))!hartree
 angm=l*(l+1)
 Pot%ww(2:n)=(Pot%rv(2:n)/Grid%r(2:n)-energy) & !*2.0 &
&     + angm/(Grid%r(2:n)*Pot%jj(2:n))!hartree
 Pot%ww(1)=0
 if (.not.Pot%finitenucleus) then
   gamma=sqrt(angm+1._dp-alpha2*nz**2)
   term=1._dp+0.25_dp*alpha2*(energy-Pot%v0)!*2.0!hartree
   Tm21=2*gamma+1;   Tm22=2*(2*gamma+2)
   !hartree!hartree!hartree
   Tm10=nz*(2._dp+alpha2*(energy-Pot%v0))-(2*balpha2/nz)*term*(gamma-1._dp)
   Tm11=nz*(2._dp+alpha2*(energy-Pot%v0))-(2*balpha2/nz)*term*(gamma)
   T00=-alpha2*nz*Pot%v0p+term*(energy-Pot%v0) + &
&      (Pot%v0p/nz+(4*balpha2**2/(nz*nz))*term**2)*(gamma-1._dp)
   c1=-Tm10/Tm21
   c2=-(Tm11*C1+T00)/Tm22
   MA=0; MB=0
 else  ! version for finite nuclear size
   gamma=l+1._dp
   term=1._dp+0.25_dp*alpha2*(energy-Pot%v0)
   Tm21=2*l+2;      Tm22=2*(2*l+3)
   Tm10=(0.25_dp*alpha2*Pot%v0p/term)*(l)
   Tm11=(0.25_dp*alpha2*Pot%v0p/term)*(l+1)
   T00=(energy-Pot%v0)*term+l*((0.25_dp*alpha2*Pot%v0p/term)**2)
   c1=-Tm10/Tm21
   c2=-(Tm11*C1+T00)/Tm22
   MA=0; MB=0
 endif
end subroutine Azeroexpand


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE wfnsrinit(Grid,l,wfn,lwfn,istart)
!! returns the solution of the scalar relativistic equations near r=0
!!  using power series expansion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wfnsrinit(Grid,l,wfn,lwfn,istart,finitenucleus,gamma,c1,c2,MA,MB,jj)
 Type(GridInfo), INTENT(IN) :: Grid
 INTEGER, INTENT(IN) :: l
 real(dp),intent(in) :: gamma,c1,c2,MA,MB
 real(dp),INTENT(INOUT) :: wfn(:),lwfn(:)
 real(dp),intent(in) :: jj(:)
 INTEGER, INTENT(OUT) :: istart
 logical, intent(in) :: finitenucleus
 real(dp) :: rr,M
 INTEGER :: i
 wfn=0; lwfn=0
 istart=6
 do i=1,istart
   rr=Grid%r(i+1)
   if (.not.finitenucleus) then
     wfn(i+1)=1+rr*(c1+rr*c2)
     lwfn(i+1)=(gamma-1)+rr*(c1*gamma+rr*c2*(gamma+1))
     wfn(i+1)=wfn(i+1)*(rr**gamma)
     lwfn(i+1)=lwfn(i+1)*(rr**gamma)/jj(i+1)
   else   ! finite nucleus case
     M=MA-MB*rr
     wfn(i+1)=(1+rr*(c1+rr*c2))*(rr**(l+1))
     lwfn(i+1)=(l+rr*((l+1)*c1+rr*(l+2)*c2))*(rr**(l+1))/M
   endif
 enddo
End SUBROUTINE wfnsrinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE wfnsrasym(Grid,wfn,lwfn,energy,iend)
!!  returns the solution of the scalar relativistic equations near r=inf
!!  using exp(-x*r) for upper component
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wfnsrasym(Grid,wfn,lwfn,energy,iend,qq,jj)
 ! returns the solution of the scalar relativistic equations near r=inf
 !  using exp(-x*r) for upper component
 Type(GridInfo), INTENT(IN) :: Grid
 real(dp),INTENT(INOUT) :: wfn(:),lwfn(:)
 real(dp), INTENT(IN) :: energy,qq
 real(dp),intent(in) :: jj(:)
 INTEGER, INTENT(OUT) :: iend
 real(dp) :: rr,x,m,qx
 INTEGER :: i,n
 if (energy>0._dp) then
   LIBPAW_ERROR('Error in wfnsrasym -- energy > 0')
 endif
 wfn=0; lwfn=0;
 n=Grid%n
 m=1._dp+0.25_dp*energy/(InvFineStruct**2)!Hartree
 x=sqrt(-m*energy)!Hartree
 qx=qq     !  Possible net ionic charge
 qx=(qx/x)*(1._dp+0.5_dp*energy/(InvFineStruct**2))!Hartree
 iend=5
 do i=n-iend,n
   wfn(i)=exp(-x*(Grid%r(i)-Grid%r(n-iend)))
   if (qx>0._dp) then
     rr=(Grid%r(i)/Grid%r(n-iend))**qx
     wfn(i)=wfn(i)*rr
   endif
   lwfn(i)=-wfn(i)*(x*Grid%r(i)+(1._dp-qx))/jj(i)
 enddo
end subroutine wfnsrasym


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      subroutine unboundsr(Grid,Pot,nr,l,energy,wfn,nodes)
!!  pgm to solve radial scalar relativistic equation for unbound states
!!    at energy 'energy' and at angular momentum l
!!
!!    with potential rv/r, given in uniform linear or log mesh of n points
!!   assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
!!
!!  nz=nuclear charge
!!
!!  Does not use Noumerov algorithm -- but uses coupled first-order
!!       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
!!
!! also returns node == number of nodes for calculated state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unboundsr(Grid,Pot,nr,l,energy,wfn,nodes)
 TYPE(GridInfo), INTENT(IN)  :: Grid
 TYPE(PotentialInfo), INTENT(INout)  :: Pot
 INTEGER, INTENT(IN) :: nr,l
 real(dp), INTENT(IN) :: energy
 real(dp), INTENT(INOUT) :: wfn(:)
 INTEGER, INTENT(INOUT) :: nodes
 INTEGER :: n,istart
 real(dp) :: scale,gamma,c1,c2,MA,MB,qq
 real(dp), allocatable :: lwfn(:),zz(:,:,:),yy(:,:)
 n=Grid%n
 IF (nr > n) THEN
   LIBPAW_ERROR('Error in unboundsr')
 ENDIF
 call Azeroexpand(Grid,Pot,l,energy,qq,gamma,c1,c2,MA,MB,nr)
 LIBPAW_ALLOCATE(lwfn,(nr))
 LIBPAW_ALLOCATE(zz,(2,2,nr))
 LIBPAW_ALLOCATE(yy,(2,nr))
 lwfn=0;zz=0;yy=0;
 call wfnsrinit(Grid,l,wfn,lwfn,istart,Pot%finitenucleus,gamma,c1,c2,MA,MB,Pot%jj)
 call prepareforcfdsol(Grid,1,istart,nr,wfn,lwfn,yy,zz,Pot%ww,Pot%jj)
 call cfdsol(Grid,zz,yy,istart,nr)
 call getwfnfromcfdsol(1,nr,yy,wfn)
 nodes=countnodes(2,nr,wfn)
 ! normalize to unity within integration range
 scale=1._dp/overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
 scale=SIGN(SQRT(scale),wfn(nr-2))
 wfn(1:nr)=wfn(1:nr)*scale
 LIBPAW_DEALLOCATE(lwfn)
 LIBPAW_DEALLOCATE(yy)
 LIBPAW_DEALLOCATE(zz)
END SUBROUTINE unboundsr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE boundsr(Grid,Pot,eig,wfn,l,nroot,emin,ierr,success)
!!    pgm to solve radial scalar relativistic equation for nroot bound state
!!      energies and wavefunctions for angular momentum l
!!      with potential rv/r, given in uniform linear or log mesh of n points
!!    nz=nuclear charge
!!    emin=is estimate of lowest eigenvalue; used if nz=0
!!       otherwise, set to the value of -(nz/(l+1))**2
!!
!!    It is assumed that the wavefunction has np-l-1 nodes, where
!!      np is the principle quantum number-- np=1,2,..nroot
!!
!!    Does not use Noumerov algorithm -- but uses coupled first-order
!!         equations from David Vanderbilt, Marc Torrent, and Francois Jollet
!!
!!    Corrections are also needed for r>n*h, depending on:
!!           e0 (current guess of energy eigenvalue
!!           the extrapolated value of rv == r * v
!!
!!   ierr=an nroot digit number indicating status of each root
!!     a digit of 1 indicates success in converging root
!!                2 indicates near success in converging root
!!                9 indicates that root not found
!!
!!   first check how many roots expected =  ntroot (returned as argument)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE boundsr(Grid,Pot,eig,wfn,l,nroot,emin,ierr,success)
 TYPE(GridInfo), INTENT(IN) :: Grid
 TYPE(PotentialInfo), INTENT(INout) :: Pot
 real(dp), INTENT(INOUT) :: eig(:),wfn(:,:)
 INTEGER, INTENT(IN) :: l,nroot
 INTEGER, INTENT(INOUT) :: ierr
 real(dp), INTENT(INOUT) :: emin
 LOGICAL, INTENT(INOUT) :: success
 real(dp), PARAMETER :: convre=tol10,vlrg=10._dp**30
 INTEGER, PARAMETER :: niter=1000
 real(dp), POINTER :: rv(:)
 real(dp), ALLOCATABLE :: p1(:),p2(:),dd(:)
 INTEGER :: n
 real(dp) :: nz,h,v0,v0p
 real(dp) :: err,convrez,energy,gamma,c1,c2,MA,MB
 real(dp) :: scale,emax,best,rout,qq
 real(dp) :: rin,dele,x
 INTEGER :: iter,i,j,node,match,mxroot,ntroot,ir,iroot
 INTEGER :: ifac,istart,iend
 LOGICAL :: ok
 real(dp), allocatable :: lwfn(:),zz(:,:,:),yy(:,:)
 n=Grid%n
 h=Grid%h
 LIBPAW_ALLOCATE(p1,(n))
 LIBPAW_ALLOCATE(p2,(n))
 LIBPAW_ALLOCATE(dd,(n))
 success=.true.
 LIBPAW_ALLOCATE(lwfn,(n))
 LIBPAW_ALLOCATE(zz,(2,2,n))
 LIBPAW_ALLOCATE(yy,(2,n))
 nz=Pot%nz
 v0=Pot%v0
 v0p=Pot%v0p
 rv=>Pot%rv
 err=n*nz*(h**4)!*0.5!hartree
 convrez=convre
 IF (nz>0.001_dp) convrez=convre*nz
 ierr=0
 if(has_to_print) write(std_out,*) 'z , l = ',nz,l
 ! check how many roots expected by integration outward at
 !   energy = 0
 energy = 0
 call Azeroexpand(Grid,Pot,l,energy,qq,gamma,c1,c2,MA,MB)
 lwfn=0;zz=0;yy=0;
 call wfnsrinit(Grid,l,p1,lwfn,istart,Pot%finitenucleus,gamma,c1,c2,MA,MB,Pot%jj)
 !start outward integration
 call prepareforcfdsol(Grid,1,istart,n,p1,lwfn,yy,zz,Pot%ww,Pot%jj)
 call cfdsoliter(Grid,zz,yy,istart,n)
 call getwfnfromcfdsol(1,n,yy,p1)
 node=countnodes(2,n,p1)
 if(has_to_print) write(std_out,*) ' nodes at e=0  ', node
 mxroot=node+1
 ntroot=node
 IF (mxroot.LT.nroot) THEN
   if(has_to_print) write(std_out,*)'error in boundsr - for l = ',l
   if(has_to_print) write(std_out,*) nroot,' states requested but only',mxroot,' possible'
   DO ir=mxroot+1,nroot
     ierr=ierr+9*(10**(ir-1))
   ENDDO
   success=.false.
 ENDIF
 mxroot=min0(mxroot,nroot)
 IF (nz.EQ.0) energy=-ABS(emin)
 IF (nz.NE.0) energy=-(1.1_dp*(nz/(l+1._dp))**2)!*0.5!Hartree
 emin=energy-err
 emax=0._dp
 DO iroot=1,mxroot
   best=1.d10; dele=1.d10
   energy=emin+err
   IF (energy.LT.emin) energy=emin
   IF (energy.GT.emax) energy=emax
   ok=.FALSE.
   BigIter: DO iter=1,niter
     !  start inward integration
     !  start integration at n
     call Azeroexpand(Grid,Pot,l,energy,qq,gamma,c1,c2,MA,MB)
     ! find classical turning point
     call ClassicalTurningPoint(Grid,Pot%rv,l,energy,match)
     match=max(match,10); match=min(match,n-20)
     call wfnsrasym(Grid,p2,lwfn,energy,iend,qq,Pot%jj)
     call prepareforcfdsol(Grid,n-iend,n,n,p2,lwfn,yy,zz,Pot%ww,Pot%jj)
     call cfdsoliter(Grid,zz,yy,n-iend,match)
     call getwfnfromcfdsol(match,n,yy,p2)
     match=match+6
     rin=Gfirstderiv(Grid,match,p2)/p2(match)
     call wfnsrinit(Grid,l,p1,lwfn,istart,Pot%finitenucleus,gamma,c1,c2,MA,MB,Pot%jj)
     call prepareforcfdsol(Grid,1,istart,n,p1,lwfn,yy,zz,Pot%ww,Pot%jj)
     call cfdsoliter(Grid,zz,yy,istart,match+6)
     call getwfnfromcfdsol(1,match+6,yy,p1)
     node= countnodes(2,match+6,p1)
     rout=Gfirstderiv(Grid,match,p1)/p1(match)
     ! check whether node = (iroot-1)
     !   not enough nodes -- raise energy
     IF (node.LT.iroot-1) THEN
       emin=MAX(emin,energy)-err
       energy=emax-(emax-energy)*ranx()
       ifac=9
       !   too many nodes -- lower energy
     ELSEIF (node.GT.iroot-1) THEN
       IF (energy.LE.emin) THEN
         ierr=ierr+9*(10**(iroot-1))
         if(has_to_print) write(std_out,*) 'boundsr error -- emin too high',l,nz,emin,energy
         IF (energy.LE.emin-tol10) THEN
           STOP
         ENDIF
       ENDIF
       emax=MIN(emax,energy+err)
       energy=emin+(energy-emin)*ranx()
       !   correct number of nodes -- estimate correction
     ELSEIF (node.EQ.iroot-1) THEN
       DO j=1,match
         p1(j)=p1(j)/p1(match)
       ENDDO
       DO j=match,n
         p1(j)=p2(j)/p2(match)
       ENDDO
       scale=1._dp/overlap(Grid,p1,p1)
       dele=(rout-rin)*scale
       x=ABS(dele)
       IF (x.LT.best) THEN
         scale=SQRT(scale)
         p1(1:n)=p1(1:n)*scale
         call filter(n,p1,machine_zero)
         wfn(1:n,iroot)=p1(1:n)
         eig(iroot)=energy
         best=x
       ENDIF
       IF (ABS(dele).LE.convrez) THEN
         ok=.TRUE.
         !  eigenvalue found
         ierr=ierr+10**(iroot-1)
         IF (iroot+1.LE.mxroot) THEN
           emin=energy+err
           emax=0
           energy=(emin+emax)/2
           IF (energy.LT.emin) energy=emin
           IF (energy.GT.emax) energy=emax
           best=1.d10
         ENDIF
         EXIT BigIter
       ENDIF
       IF (ABS(dele).GT.convrez) THEN
         energy=energy+dele!*0.5!hartree
         ! if energy is out of range, pick random energy in correct range
         IF (emin-energy.GT.convrez.OR.energy-emax.GT.convrez)         &
              energy=emin+(emax-emin)*ranx()
         ifac=2
       ENDIF
     ENDIF
   ENDDO BigIter !iter
   IF (.NOT.ok) THEN
     success=.false.
     ierr=ierr+ifac*(10**(iroot-1))
     if(has_to_print) write(std_out,*) 'no convergence in boundsr',iroot,l,dele,energy
     if(has_to_print) write(std_out,*) ' best guess of eig, dele = ',eig(iroot),best
     IF (iroot.LT.mxroot) THEN
       DO ir=iroot+1,mxroot
         ierr=ierr+9*(10**(ir-1))
       ENDDO
     ENDIF
     ! reset wfn with hydrogenic form
     j=iroot+l+1
     wfn(:,iroot)=0
     x=(j)*sqrt(abs(eig(iroot)*2.0))
     do i=2,n
       wfn(i,iroot)=hwfn(x,j,l,Grid%r(i))
     enddo
   ENDIF
 ENDDO !iroot
 LIBPAW_DEALLOCATE(p1)
 LIBPAW_DEALLOCATE(p2)
 LIBPAW_DEALLOCATE(dd)
 LIBPAW_DEALLOCATE(lwfn)
 LIBPAW_DEALLOCATE(yy)
 LIBPAW_DEALLOCATE(zz)
END SUBROUTINE Boundsr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! prepareforcfdsol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prepareforcfdsol(Grid,i1,i2,n,wfn,lwfn,yy,zz,ww,jj)
 Type(gridinfo), INTENT(IN) :: Grid
 INTEGER, INTENT(IN) :: i1,i2,n
 real(dp), INTENT(IN) :: wfn(:),lwfn(:)
 real(dp), INTENT(OUT) :: yy(:,:),zz(:,:,:)
 real(dp),intent(in) :: ww(:),jj(:)
 INTEGER :: i
 yy=0;zz=0
 yy(1,i1:i2)=wfn(i1:i2)
 yy(2,i1:i2)=lwfn(i1:i2)
 do  i=2,n
   zz(1,1,i)=1._dp/Grid%r(i)
   zz(1,2,i)=jj(i)/Grid%r(i)
   zz(2,2,i)=-1._dp/Grid%r(i)
   zz(2,1,i)=ww(i)
 enddo
end subroutine prepareforcfdsol





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 7. anderson_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Anderson_Mix
!! Performs the actual mixing of the input vector with the
!!                history and retuns the result.
!!
!!   AC - Anderson context
!!   X  - Current vector on input and new guess on output
!!   F  - F(X) - X. Nonlinear mixing of input vector
!!
!! Modified to call SVD routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Anderson_Mix(AC, X, F)
 SAVE
 TYPE  (Anderson_Context), INTENT(INOUT) :: AC
 real(dp),                  INTENT(INOUT) :: X(:)
 real(dp),                  INTENT(IN)    :: F(:)
 INTEGER :: i, slot, currentdim , n ,j
 real(dp) :: term
 real(dp)  :: tmp
 !First determine where to store the new correction vectors ***
 AC%slot = AC%slot + 1
 IF (AC%Slot>AC%Nmax) AC%Slot = 1
 IF ((AC%N < 0) .OR. (AC%Nmax == 0)) THEN  !** Simple mixing for 1st time ***
   AC%Xprev = X
   X = X + AC%NewMix*F
 ELSE
   slot = AC%Slot
   AC%DF(:,slot) = F - AC%Fprev   !** Make new DF vector
   AC%DX(:,slot) = X - AC%Xprev   !** Make new DX vector
   currentdim=MIN(AC%N+1,AC%Nmax)
   DO i=1, currentdim              !*** Add row/col to matrix
     term = DOT_PRODUCT(AC%DF(:,i), AC%DF(:,slot))
     AC%Matrix(i,slot) = term
     IF (i /= slot) AC%Matrix(slot,i) = (term)
     AC%Gamma(i) = DOT_PRODUCT(AC%DF(:,i), F)
   END DO
   AC%DupMatrix = AC%Matrix
   n = AC%Nmax;   j= currentdim
   CALL DGESDD('A',j,j,AC%DupMatrix(1,1),n,AC%S(1), &
        AC%U(1,1),n,AC%VT(1,1),n,AC%Work(1),AC%Lwork, AC%IPIV(1),i)
   IF (i /= 0) THEN
     LIBPAW_ERROR('Anderson_Mix: Error in DGESDD.')
   END IF
   AC%Work(1:j) = AC%Gamma(1:j)
   AC%Gamma = 0
   tmp=MAX(ABS(AC%S(1))/AC%ConditionNo,AC%Machaccur)
   DO i=1,j
     IF (ABS(AC%S(i)).GT.tmp) THEN
       AC%Gamma(1:j)=AC%Gamma(1:j)+&
            (AC%VT(i,1:j))*DOT_PRODUCT(AC%U(1:j,i),AC%Work(1:j))/AC%S(i)
     ENDIF
   ENDDO
   AC%Xprev = X
   !*** Now calculate the new vector ***
   X = X + AC%NewMix*F
   DO i=1, currentdim               ! updated vector
     X = X - AC%Gamma(i)*(AC%DX(:,i) + AC%NewMix*AC%DF(:,i))
   END DO
 END IF
 AC%Fprev = F
 AC%N = AC%N + 1
 IF (AC%N > AC%Nmax) AC%N = AC%Nmax
 RETURN
END SUBROUTINE Anderson_Mix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Anderson_ResetMix - Resets the mixing history to None
!!     AC - Anderson context to reset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Anderson_ResetMix(AC)
 TYPE  (Anderson_Context), INTENT(INOUT) :: AC
 AC%N = -1
 AC%Slot = -1
 AC%CurIter=0
 RETURN
END SUBROUTINE Anderson_ResetMix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  FreeAnderson - Frees all the data associated with the AC data structure
!!      AC -Pointer to the Anderson context to free
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FreeAnderson(AC)
 TYPE (Anderson_Context), INTENT(INOUT) :: AC
 IF (ASSOCIATED(AC%Matrix)) then
   LIBPAW_POINTER_DEALLOCATE(AC%Matrix)
 endif
 IF (ASSOCIATED(AC%Gamma)) then
   LIBPAW_POINTER_DEALLOCATE(AC%Gamma)
 endif
 IF (ASSOCIATED(AC%DF)) then
   LIBPAW_POINTER_DEALLOCATE(AC%DF)
 endif
 IF (ASSOCIATED(AC%Fprev)) then
   LIBPAW_POINTER_DEALLOCATE(AC%Fprev)
 endif
 IF (ASSOCIATED(AC%DX)) then
   LIBPAW_POINTER_DEALLOCATE(AC%DX)
 endif
 IF (ASSOCIATED(AC%Xprev)) then
   LIBPAW_POINTER_DEALLOCATE(AC%Xprev)
 endif
 IF (ASSOCIATED(AC%IPIV)) then
   LIBPAW_POINTER_DEALLOCATE(AC%IPIV)
 endif
 IF (ASSOCIATED(AC%S)) then
   LIBPAW_POINTER_DEALLOCATE(AC%S)
 endif
 IF (ASSOCIATED(AC%RWork)) then
   LIBPAW_POINTER_DEALLOCATE(AC%RWork)
 endif
 IF (ASSOCIATED(AC%U)) then
   LIBPAW_POINTER_DEALLOCATE(AC%U)
 endif
 IF (ASSOCIATED(AC%VT)) then
   LIBPAW_POINTER_DEALLOCATE(AC%VT)
 endif
 IF (ASSOCIATED(AC%Work)) then
   LIBPAW_POINTER_DEALLOCATE(AC%Work)
 endif
 IF (ASSOCIATED(AC%DupMatrix)) then
   LIBPAW_POINTER_DEALLOCATE(AC%DupMatrix)
 endif
 RETURN
END SUBROUTINE FreeAnderson


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  InitAnderson_dr - Initializes and Anderson_Context data structure for use
!!   AC       - Anderson context created and returned
!!   Err_Unit - Output error unit
!!   Nmax     - Max number of vectors to keep
!!   VecSize  - Size of each vector
!!   NewMix   - Mixing factor
!!   CondNo   - For matrix inversion
!!   MaxIter  - Maximum number of iterations
!!   err      - minimum residue convergence tolerance
!!   toosmall - result obviously converged
!!   verbose  - if true -- write out results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE InitAnderson_dr(AC,Err_Unit,Nmax,VecSize,NewMix,CondNo,&
&      MaxIter,err,toosmall,verbose)
 TYPE (Anderson_Context), INTENT(INOUT)     :: AC
 INTEGER,                 INTENT(IN)  :: Err_Unit
 INTEGER,                 INTENT(IN)  :: Nmax
 INTEGER,                 INTENT(IN)  :: VecSize
 real(dp),                    INTENT(IN)  :: NewMix
 real(dp),                    INTENT(IN)  :: CondNo
 INTEGER, INTENT(IN) :: MaxIter
 real(dp), INTENT(IN) :: err,toosmall
 LOGICAL, INTENT(IN) :: verbose
 real(dp)    :: a1,a2,a3
 AC%Nmax = Nmax          !*** Store the contants
 AC%VecSize = VecSize
 AC%NewMix = NewMix
 AC%Err_Unit = Err_Unit
 AC%MaxIter = MaxIter
 AC%err = err
 AC%toosmall = toosmall
 AC%writelots=verbose
 AC%N = -1                !** Init the rest of the structure
 AC%Slot = -1
 AC%CurIter=0
 LIBPAW_POINTER_ALLOCATE(AC%Xprev,(VecSize))
 LIBPAW_POINTER_ALLOCATE(AC%Fprev,(VecSize))
 LIBPAW_POINTER_ALLOCATE(AC%DX,(VecSize,Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%DF,(VecSize,Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%Matrix,(Nmax,Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%Gamma,(Nmax))
 AC%Lwork=5*Nmax*Nmax+10*Nmax
 AC%LRwork= 5*Nmax*Nmax+7*Nmax
 AC%ConditionNo= CondNo
 ! Calculate machine accuracy
 AC%Machaccur = 0
 a1 = 4._dp/3._dp
 DO WHILE (AC%Machaccur == 0._dp)
   a2 = a1 - 1._dp
   a3 = a2 + a2 + a2
   AC%Machaccur = ABS(a3 - 1._dp)
 ENDDO
 LIBPAW_POINTER_ALLOCATE(AC%DupMatrix,(Nmax,Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%U,(Nmax, Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%VT,(Nmax,Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%Work,(AC%Lwork))
 LIBPAW_POINTER_ALLOCATE(AC%RWork,(AC%LRWork))
 LIBPAW_POINTER_ALLOCATE(AC%IPIV,(8*Nmax))
 LIBPAW_POINTER_ALLOCATE(AC%S,(Nmax))
 AC%Matrix = 0
 RETURN
END SUBROUTINE InitAnderson_dr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  DoAndersonMix
!!    Note residue can be wout-w   or more general residue that tends --> 0
!!    at convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DoAndersonMix(AC,w,E,Esub,success,atp)
 TYPE (Anderson_Context), INTENT(INOUT) :: AC
 real(dp), INTENT(INOUT) :: E,w(:)
 !     External :: Esub
 LOGICAL, INTENT(OUT) :: success
 type(atompaw_type), intent(inout) :: atp
 real(dp), ALLOCATABLE :: residue(:),tmp(:)
 real(dp) :: err,v1,v2,v3,v4
 INTEGER :: i,n
 real(dp), PARAMETER :: conv1=4.d13,conv2=3.d13,conv3=2.d13,conv4=1.d13
 LOGICAL :: OK
 INTERFACE
   SUBROUTINE Esub(w,energy,residue,err,OK,update,atp)
     USE_DEFS
     import atompaw_type
     real(dp), INTENT(INOUT) :: w(:)
     real(dp), INTENT(OUT) :: energy
     real(dp), INTENT(OUT) :: residue(:)
     real(dp), INTENT(OUT) :: err
     LOGICAL, INTENT(OUT) :: OK
     LOGICAL, INTENT(IN)  :: update
     type(atompaw_type), intent(inout) :: atp
   END SUBROUTINE Esub
 END INTERFACE
 n=SIZE(w);success=.FALSE.
 LIBPAW_ALLOCATE(residue,(n))
 LIBPAW_ALLOCATE(tmp,(n))
 err=1.0d10
 v1=conv1;v2=conv2;v3=conv3;v4=conv4;tmp=0
 DO i=1,AC%MaxIter
   AC%CurIter=i
   CALL  Esub(w,E,residue,err,OK,.TRUE.,atp)
   AC%res=err
   if (err<AC%toosmall) THEN
           If(AC%writelots)&
     write(std_out,&
&       '("AndersonMix converged in ",i5," iterations with err = ",1p,1e15.7)')&
&            i, err
     EXIT
   endif
   CALL shift4(v1,v2,v3,v4,err)
   IF (i>=4.AND.OK) THEN
     IF ((.NOT.(v4.LE.v3.AND.v3.LE.v2 &
&         .AND.v2.LE.v1).AND.v4.LE.AC%err).OR.err<AC%toosmall) THEN
        !  converged result
       success=.TRUE.
        If(AC%writelots)&
       write(std_out,&
&    '("AndersonMix converged in ",i5," iterations with err = ",1p,1e15.7)')&
&         i, err
       EXIT
     ENDIF
   ENDIF
   If(AC%writelots)write(std_out,'("AndersonMixIter ",i7,2x,1p,2e20.12)') i,E,err
   IF (.NOT.OK) THEN
     CALL Anderson_ResetMix(AC)
     IF (i>1) THEN
       w=tmp
       AC%NewMix=MAX(0.00001_dp,AC%NewMix/2)
       IF (AC%NewMix<=0.00001_dp) THEN
         write(std_out,*) 'Sorry -- this is not working '
         STOP
       ENDIF
     ENDIF
   ELSE
     AC%NewMix=MIN(MaxMix,AC%NewMix*2)
   ENDIF
   tmp=w
   CALL Anderson_Mix(AC,w,residue)
 ENDDO
 If (AC%CurIter.ge.AC%MaxIter) then
   if(has_to_print) WRITE(STD_OUT,*) 'Anderson Mix reached MaxIter without success',AC%MaxIter
 Endif
 LIBPAW_DEALLOCATE(residue)
 LIBPAW_DEALLOCATE(tmp)
END SUBROUTINE DoAndersonMix





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 8. global_math
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine shapebes(al,ql,ll,rc)
!!    Find al and ql parameters for a "Bessel" shape function:
!!    Shape(r)=al1.jl(ql1.r)+al2.jl(ql2.r)
!!      such as Shape(r) and 2 derivatives are zero at r=rc
!!              Intg_0_rc[Shape(r).r^(l+2).dr]=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE shapebes(al,ql,ll,rc)
  INTEGER,INTENT(IN) :: ll
  REAL(dp),INTENT(IN) :: rc
  REAL(dp),INTENT(OUT) :: al(2),ql(2)
  INTEGER :: i
  REAL(dp) :: alpha,beta,det,qr,jbes,jbesp,jbespp,amat(2,2),bb(2)
  alpha=1.D0;beta=0.D0
  CALL solvbes(ql,alpha,beta,ll,2)
  ql(1:2)=ql(1:2)/rc
  DO i=1,2
    qr=ql(i)*rc
    CALL jbessel(jbes,jbesp,jbespp,ll,1,qr)
    amat(1,i)=jbesp*ql(i)
    CALL jbessel(jbes,jbesp,jbespp,ll+1,0,qr)
    amat(2,i)=jbes*rc**(ll+2)/ql(i)  !  Intg_0_rc[jl(qr).r^(l+2).dr]
  ENDDO
  bb(1)=0.d0;bb(2)=1.d0
  det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
  al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
  al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det
END SUBROUTINE shapebes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ddexp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ddexp(arg)
 real(dp) :: arg,ddexp
 IF (arg>maxexparg) THEN
   ddexp=maxexp
 ELSE IF (arg<minexparg) THEN
   ddexp=minexp
 ELSE
   ddexp=EXP(arg)
 ENDIF
 RETURN
END FUNCTION ddexp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ddlog
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ddlog(arg)
 real(dp) :: arg,ddlog
 IF (arg>maxlogarg) THEN
   ddlog=maxlog
 ELSE IF (arg<minlogarg) THEN
   ddlog=minlog
 ELSE
   ddlog=LOG(arg)
 ENDIF
 RETURN
END FUNCTION ddlog


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ranx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ranx()
 real(dp) :: ranx
 INTEGER, PARAMETER :: konst=125
 INTEGER  :: m=100001
 m=m*konst
 m=m-2796203*(m/2796203)
 ranx=m/2796203._dp
 RETURN
END FUNCTION ranx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  factorial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION factorial(n)
 real(dp) :: factorial
 INTEGER, INTENT(IN) :: n
 INTEGER :: i
 factorial=one
 IF (n.LT.2) RETURN
 DO i=2,n
   factorial=factorial*i
 ENDDO
END FUNCTION factorial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCTION hwfn(z,np,l,r)
!! function to calculate the radial H wfn for nuclear charge z
!!          (note in this version z is real and need not be integral)
!!                                            principal qn   np
!!                                            orbital qn     l
!!   r*(radial H wfn) is returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION hwfn(z,np,l,r)
 real(dp) :: hwfn
 real(dp), INTENT(IN) :: z,r
 INTEGER, INTENT(IN) :: np,l
 INTEGER :: node,k
 real(dp) :: scale_,rho,pref,term,sum_
 node=np-l-1
 scale_=2._dp*z/np
 rho=scale_*r
 pref=scale_*SQRT(scale_*factorial(np+l)/(2*np*factorial(node)))
 if(rho==zero.and.l==0) then
   term=one/factorial(2*l+1)
 else
   term=(rho**l)/factorial(2*l+1)
 endif
 sum_=term
 IF (node.GT.0) THEN
   DO k=1,node
     term=-term*(node-k+1)*rho/(k*(2*l+1+k))
     sum_=sum_+term
   ENDDO
 ENDIF
 hwfn=r*pref*ddexp(-0.5_dp*rho)*sum_
END FUNCTION hwfn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine dirachwfn(np,kappa,z,r,eig,g,f)
!!   Subroutine to calculate eigenenergy and radial wavefunctions*r
!!      for bound state solutions to the Hydrogenic Dirac equation
!!      for nuclear charge z.   Energy in Rydberg atomic units
!!      np is principal quantum number --
!!           np=abs(kappa), abs(kappa)+1, abs(kappa)+2 ..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dirachwfn(np,kappa,z,r,eig,g,f)
 INTEGER, INTENT(IN) :: np, kappa
 real(dp), INTENT(IN) :: z,r
 real(dp), INTENT(INOUT) :: eig,g,f
 INTEGER :: ak,nr
 real(dp) :: norm,s,rho,ne, term1, term2, term0,x
 ak=abs(kappa)
 nr=np-ak
 s=sqrt(ak**2-(z**2)*(fsalpha2))
 ne=sqrt(np**2-2*nr*(ak-s))
 rho=2*z*r/ne
 norm=gammafunc(2*s+nr+1._dp)/(gammafunc(nr+1._dp)*4*ne*(ne-kappa))
 norm=sqrt(norm*2*z/ne)/gammafunc(2*s+1._dp)
 term1=0._dp
 if(nr>0) term1=nr*kummer(-nr+1,2*s+1._dp,rho)
 term2=(ne-kappa)*kummer(-nr,2*s+1._dp,rho)
 term0=norm*ddexp(-0.5_dp*rho)*(rho**s)
 eig=1._dp + ((z**2)*(fsalpha2))/(np - ak +s)**2
 eig=2*ifsalpha2*(1._dp/sqrt(eig) - 1._dp)
 x=0.5_dp*fsalpha2*eig
 g=sqrt(2._dp+x)*term0*(term2-term1)
 f=-sqrt(-x)*term0*(term2+term1)
end subroutine dirachwfn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine filter(n,func,small)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE filter(n,func,small)
 INTEGER, INTENT(IN) :: n
 real(dp), INTENT(INOUT) :: func(:)
 real(dp), INTENT(IN) :: small
 INTEGER :: i
 DO i=1,n
   IF (ABS(func(i)).LT.small) func(i)=0._dp
 ENDDO
END SUBROUTINE filter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine conthomas(n,o,d,sol)
!  use Thomas's algorithm for inverting matrix
!    Dale U. von Rosenberg, "Methods for the Numerical Solution of
!      Partial Differential Equations,
!         Am. Elsevier Pub., 1969, pg. 113
!    On input, sol contains the RHS of the equation
!    On ouput, sol contains the solution of the equation
!     Equation:  o*sol(i-1)+d*sol(i)+o*sol(i+1) = RHS(i)
!       sol(1)==sol(n+1)==0
!     simplified version for constant tridiagonal terms --
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE conthomas(n,o,d,sol)
 INTEGER, INTENT(IN) :: n
 real(dp), INTENT(IN) :: o,d
 real(dp), INTENT(INOUT) :: sol(:)
 real(dp), ALLOCATABLE :: a(:),b(:)
 real(dp) :: ss2
 INTEGER :: i
 LIBPAW_ALLOCATE(a,(n))
 LIBPAW_ALLOCATE(b,(n))
 a(2)=d
 ss2=o*o
 DO i=3,n
   a(i)=d-ss2/a(i-1)
 ENDDO
 b(2)=sol(2)/d
 DO i=3,n
   b(i)=(sol(i)-o*b(i-1))/a(i)
 ENDDO
 sol(n)=b(n)
 DO i=n-1,2,-1
   sol(i)=b(i)-o*sol(i+1)/a(i)
 ENDDO
 sol(1)=0
 LIBPAW_DEALLOCATE(a)
 LIBPAW_DEALLOCATE(b)
END SUBROUTINE conthomas


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   function kummer(n,b,z)
!     function to return confluent hypergeometric function (Kummer)
!         as defined in Handbook of mathematical functions pg. 504
!         assumes n=0, -1, -2, .. for polymomials of order 0, 1, 2, etc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION kummer(n,b,z)
 real(dp) :: kummer
 INTEGER, INTENT(IN) :: n
 real(dp), INTENT(IN) :: b,z
 INTEGER :: i,k,j,nn
 real(dp) :: num,den,fac,bb
 if (n>0) then
   LIBPAW_ERROR('Error in kummer function -- n>0')
 endif
 kummer=1._dp
 if (n==0) return
 j=-n
 k=1
 nn=n
 bb=b
 num=n
 den=b
 fac=z*nn/(bb*k)
 kummer=kummer+fac
 if (j>1) then
   do i=1,j-1
     nn=nn+1
     bb=bb+1
     k=k+1
     fac=fac*z*nn/(bb*k)
     kummer=kummer+fac
   enddo
 endif
END FUNCTION kummer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gamma function function obtained from netlib.org
!  This routine calculates the GAMMA function for a real argument X
!   Computation is based on an algorithm outlined in reference 1
!   The program uses rational functions that approximate the GAMM
!   function to at least 20 significant decimal digits.  Coefficient
!   for the approximation over the interval (1,2) are unpublished
!   Those for the approximation for X .GE. 12 are from reference 2
!   The accuracy achieved depends on the arithmetic system, th
!   compiler, the intrinsic functions, and proper selection of th
!   machine-dependent constants
!
!******************************************************************
!
! Explanation of machine-dependent constant
!
! beta   - radix for the floating-point representatio
! maxexp - the smallest positive power of beta that overflow
! XBIG   - the largest argument for which GAMMA(X) is representabl
!          in the machine, i.e., the solution to the equatio
!                  GAMMA(XBIG) = beta**maxex
! XINF   - the largest machine representable floating-point number
!          approximately beta**maxex
! EPS    - the smallest positive floating-point number such tha
!          1.0+EPS .GT. 1.
! XMININ - the smallest positive floating-point number such tha
!          1/XMININ is machine representabl
!
!     Approximate values for some important machines are
!
!                            beta       maxexp        XBI
!
! CRAY-1         (S.P.)        2         8191        966.96
! Cyber 180/85
!   under NOS    (S.P.)        2         1070        177.80
! IEEE (IBM/XT
!   SUN, etc.)   (S.P.)        2          128        35.04
! IEEE (IBM/XT
!   SUN, etc.)   (D.P.)        2         1024        171.62
! IBM 3033       (D.P.)       16           63        57.57
! VAX D-Format   (D.P.)        2          127        34.84
! VAX G-Format   (D.P.)        2         1023        171.48
!
!                            XINF         EPS        XMINI
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-246
! Cyber 180/85
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-29
! IEEE (IBM/XT
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-3
! IEEE (IBM/XT
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-30
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-7
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-3
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-30
!
!******************************************************************
! Error return
!
!  The program returns the value XINF for singularities o
!     when overflow would occur.  The computation is believe
!     to be free of underflow and overflow
!
!
!  Intrinsic functions required are
!
!     INT, DBLE, EXP, LOG, REAL, SI
!
!
! References: "An Overview of Software Development for Specia
!              Functions", W. J. Cody, Lecture Notes in Mathematics
!              506, Numerical Analysis Dundee, 1975, G. A. Watso
!              (ed.), Springer Verlag, Berlin, 1976
!
!              Computer Approximations, Hart, Et. Al., Wiley an
!              sons, New York, 1968
!
!  Latest modification: October 12, 1989
!
!  Authors: W. J. Cody and L. StoltZ
!           Applied Mathematics DivisioN
!           Argonne National LaboratorY
!           Argonne, IL 60439
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GAMMAFUNC(X)
 real(dp) :: GAMMAFUNC,X
 INTEGER:: I,N
 real(dp) :: &
&  CONV,EPS,FACT,HALF,ONE,RES,SQRTPI,SUM,TWELVE, &
&  TWO,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
 real(dp) :: C(7),P(8),Q(8)
 LOGICAL :: PARITY
!---------------------------------------------------------------------
!  Mathematical constants
!---------------------------------------------------------------------
 DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/, &
&     SQRTPI/0.9189385332046727417803297D0/
!!!!&     PI/3.1415926535897932384626434D0/   (already defined)
!---------------------------------------------------------------------
!  Machine dependent parameter
!---------------------------------------------------------------------
 DATA XBIG,XMININ,EPS,XINF/171.624D0,2.23D-308,2.22D-16,1.79D308/
!---------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minima
!     approximation over (1,2)
!---------------------------------------------------------------------
 DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1, &
&        -3.79804256470945635097577D+2,6.29331155312818442661052D+2, &
&        8.66966202790413211295064D+2,-3.14512729688483675254357D+4, &
&       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
 DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2, &
&       -1.01515636749021914166146D+3,-3.10777167157231109440444D+3, &
&        2.25381184209801510330112D+4,4.75584627752788110767815D+3, &
&       -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!---------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF)
!---------------------------------------------------------------------
 DATA C/-1.910444077728D-03,8.4171387781295D-04 , &
&        -5.952379913043012D-04,7.93650793500350248D-04 , &
&        -2.777777777777681622553D-03,8.333333333333333331554247D-02 , &
&         5.7083835261D-03/
!---------------------------------------------------------------------
!  Statement functions for conversion between integer and floa
!---------------------------------------------------------------------
 CONV(I) = DBLE(I)
 PARITY = .FALSE.
 FACT = ONE
 N = 0
 Y = X
 IF (Y .LE. ZERO) THEN
!---------------------------------------------------------------------
!  Argument is negative
!---------------------------------------------------------------------
   Y = -X
   Y1 = AINT(Y)
   RES = Y - Y1
   IF (RES .NE. ZERO) THEN
     IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
     FACT = -PI / SIN(PI*RES)
     Y = Y + ONE
   ELSE
     RES = XINF
     GO TO 900
   END IF
 END IF
!---------------------------------------------------------------------
!  Argument is positiv
!---------------------------------------------------------------------
 IF (Y .LT. EPS) THEN
!---------------------------------------------------------------------
!  Argument .LT. EPS
!---------------------------------------------------------------------
   IF (Y .GE. XMININ) THEN
     RES = ONE / Y
   ELSE
     RES = XINF
     GO TO 900
   END IF
 ELSE IF (Y .LT. TWELVE) THEN
   Y1 = Y
   IF (Y .LT. ONE) THEN
!---------------------------------------------------------------------
!  0.0 .LT. argument .LT. 1.
!---------------------------------------------------------------------
     Z = Y
     Y = Y + ONE
   ELSE
!---------------------------------------------------------------------
!  1.0 .LT. argument .LT. 12.0, reduce argument if necessar
!---------------------------------------------------------------------
     N = INT(Y) - 1
     Y = Y - CONV(N)
     Z = Y - ONE
   END IF
!---------------------------------------------------------------------
!  Evaluate approximation for 1.0 .LT. argument .LT. 2.
!---------------------------------------------------------------------
   XNUM = ZERO
   XDEN = ONE
   DO 260 I = 1, 8
     XNUM = (XNUM + P(I)) * Z
     XDEN = XDEN * Z + Q(I)
  260   CONTINUE
   RES = XNUM / XDEN + ONE
   IF (Y1 .LT. Y) THEN
!---------------------------------------------------------------------
!  Adjust result for case  0.0 .LT. argument .LT. 1.
!---------------------------------------------------------------------
     RES = RES / Y1
   ELSE IF (Y1 .GT. Y) THEN
!---------------------------------------------------------------------
!  Adjust result for case  2.0 .LT. argument .LT. 12.
!---------------------------------------------------------------------
     DO 290 I = 1, N
       RES = RES * Y
       Y = Y + ONE
  290  CONTINUE
   END IF
 ELSE
!---------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0
!---------------------------------------------------------------------
   IF (Y .LE. XBIG) THEN
     YSQ = Y * Y
     SUM = C(7)
     DO 350 I = 1,6
       SUM = SUM / YSQ + C(I)
  350 CONTINUE
     SUM = SUM/Y - Y + SQRTPI
     SUM = SUM + (Y-HALF)*LOG(Y)
     RES = EXP(SUM)
   ELSE
     RES = XINF
     GO TO 900
   END IF
 END IF
!---------------------------------------------------------------------
!  Final adjustments and retur
!---------------------------------------------------------------------
 IF (PARITY) RES = -RES
 IF (FACT .NE. ONE) RES = FACT / RES
900 GAMMAFUNC = RES
 RETURN
! ---------- Last line of GAMMA ---------
END FUNCTION GAMMAFUNC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE shift4(v1,v2,v3,v4,NEW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE shift4(v1,v2,v3,v4,NEW)
 real(dp), INTENT(IN) :: NEW
 real(dp), INTENT(INOUT) :: v1,v2,v3,v4
 v1=v2
 v2=v3
 v3=v4
 v4=NEW
 RETURN
END SUBROUTINE shift4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine jbessel(bes,besp,bespp,ll,order,xx)
!    Spherical bessel function and derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE jbessel(bes,besp,bespp,ll,order,xx)
 INTEGER,INTENT(IN) :: ll,order
 real(dp),INTENT(IN) :: xx
 real(dp),INTENT(OUT) :: bes,besp,bespp
 INTEGER,PARAMETER :: imax=40
 real(dp),PARAMETER :: prec=tol15
 INTEGER :: ii,il
 real(dp) :: besp1,fact,factp,factpp,jn,jnp,jnpp,jr,xx2,xxinv
 IF (order>2) STOP "Wrong order in jbessel !"
 IF (ABS(xx)<prec) THEN
   bes=0._dp;IF (ll==0) bes=1._dp
   IF (order>=1) THEN
     besp=0._dp;IF (ll==1) besp=1._dp/3._dp
   ENDIF
   IF (order==2) THEN
     bespp=0._dp
     IF (ll==0) bespp=-1._dp/3._dp
     IF (ll==2) bespp=2._dp/15._dp
   ENDIF
   RETURN
 ENDIF
 xxinv=1._dp/xx
 IF (xx<1._dp) THEN
   xx2=0.5_dp*xx*xx
   fact=1.D0;DO il=1,ll;fact=fact*xx/DBLE(2*il+1);ENDDO
   jn=1.D0;jr=1.D0;ii=0
   DO WHILE(ABS(jr)>=prec.AND.ii<imax)
     ii=ii+1;jr=-jr*xx2/DBLE(ii*(2*(ll+ii)+1))
     jn=jn+jr
   ENDDO
   bes=jn*fact
   IF (ABS(jr)>prec) STOP 'Error: Bessel function did not converge !'
   IF (order>=1) THEN
     factp=fact*xx/DBLE(2*ll+3)
     jnp=1.D0;jr=1.D0;ii=0
     DO WHILE(ABS(jr)>=prec.AND.ii<imax)
       ii=ii+1;jr=-jr*xx2/DBLE(ii*(2*(ll+ii)+3))
       jnp=jnp+jr
     ENDDO
     besp=-jnp*factp+jn*fact*xxinv*DBLE(ll)
     IF (ABS(jr)>prec) STOP 'Error: 1st der. of Bessel function did not converge !'
   ENDIF
   IF (order==2) THEN
     factpp=factp*xx/DBLE(2*ll+5)
     jnpp=1.D0;jr=1.D0;ii=0
     DO WHILE(ABS(jr)>=prec.AND.ii<imax)
       ii=ii+1;jr=-jr*xx2/DBLE(ii*(2*(ll+ii)+5))
       jnpp=jnpp+jr
     ENDDO
     besp1=-jnpp*factpp+jnp*factp*xxinv*DBLE(ll+1)
     IF (ABS(jr)>prec) STOP 'Error: 2nd der. of Bessel function did not converge !'
   ENDIF
 ELSE
   jn =SIN(xx)*xxinv
   jnp=(-COS(xx)+jn)*xxinv
   DO il=2,ll+1
     jr=-jn+DBLE(2*il-1)*jnp*xxinv
     jn=jnp;jnp=jr
   ENDDO
   bes=jn
   IF (order>=1) besp =-jnp+jn *xxinv*DBLE(ll)
   IF (order==2) besp1= jn -jnp*xxinv*DBLE(ll+2)
 ENDIF
 IF (order==2) bespp=-besp1+besp*ll*xxinv-bes*ll*xxinv*xxinv
END SUBROUTINE jbessel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine solvbes(root,alpha,l,nq)
!    Find nq first roots of instrinsic equation:
!                            alpha.jl(Q) + beta.Q.djl/dr(Q) = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solvbes(root,alpha,beta,ll,nq)
 INTEGER,INTENT(IN) :: ll,nq
 real(dp),INTENT(IN) :: alpha,beta
 real(dp),INTENT(OUT) :: root(nq)
 real(dp),PARAMETER :: dh=1.D-1, tol=1.D-14
 INTEGER :: nroot
 real(dp) :: dum,y1,y2,jbes,jbesp,qq,qx,hh
 qq=dh;nroot=0
 DO WHILE (nroot<nq)
   CALL jbessel(jbes,jbesp,dum,ll,1,qq)
   y1=alpha*jbes+beta*qq*jbesp
   qq=qq+dh
   CALL jbessel(jbes,jbesp,dum,ll,1,qq)
   y2=alpha*jbes+beta*qq*jbesp
   DO WHILE (y1*y2>=0.D0)
     qq=qq+dh
     CALL jbessel(jbes,jbesp,dum,ll,1,qq)
     y2=alpha*jbes+beta*qq*jbesp
   ENDDO
   hh=dh;qx=qq
   DO WHILE (hh>tol)
     hh=0.5D0*hh
     IF (y1*y2<0) THEN
       qx=qx-hh
     ELSE
       qx=qx+hh
     ENDIF
     CALL jbessel(jbes,jbesp,dum,ll,1,qx)
     y2=alpha*jbes+beta*qx*jbesp
   ENDDO
   nroot=nroot+1
   root(nroot)=qx
 ENDDO
END SUBROUTINE solvbes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    SUBROUTINE linsol(a,b,kk,la,ra,lb,det)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE linsol(a,b,kk,la,ra,lb,det)
 INTEGER, INTENT(IN) :: kk,la,ra,lb
 real(dp), INTENT(INOUT) :: a(la,ra),b(lb)
 real(dp), OPTIONAL, INTENT(OUT) :: det
 real(dp) :: d,s,r
 INTEGER :: kkm,i,j,k,l,ipo,n,kmo
 d = 1._dp
 if (kk>min(la,ra,lb)) then
   LIBPAW_ERROR('Dimension error in linsol ')
 endif
 kkm=kk-1
 IF (kkm == 0) THEN
   b(1)=b(1)/a(1,1)
 ELSE IF (kkm > 0) THEN
   DO i=1, kkm
     s = 0.0_dp
     l=i
     DO j=i,kk
       r=ABS(a(j,i))
       IF(r >  s) THEN
         s=r
         l=j
       ENDIF
     ENDDO
     IF(l /= i) THEN
       DO j=i,kk
         s=a(i,j)
         a(i,j)=a(l,j)
         a(l,j)=s
       ENDDO
       s=b(i)
       b(i)=b(l)
       b(l)=s
       d = -d
     ENDIF
     IF (a(i,i) /= 0.0_dp) THEN
       ipo=i+1
       DO j=ipo,kk
         IF (a(j,i) /= 0.0_dp) THEN
           s=a(j,i)/a(i,i)
           a(j,i) = 0.0_dp
           DO k=ipo,kk
             a(j,k)=a(j,k)-a(i,k)*s
           ENDDO
           b(j)=b(j)-b(i)*s
         ENDIF
       ENDDO
     ENDIF
   ENDDO
   DO i=1,kk
     d=d*a(i,i)
   ENDDO
   kmo=kk-1
   b(kk)=b(kk)/a(kk,kk)
   DO i=1,kmo
     n=kk-i
     DO j=n,kmo
       b(n)=b(n)-a(n,j+1)*b(j+1)
     ENDDO
     b(n)=b(n)/a(n,n)
   ENDDO
 ENDIF
 !write(std_out,*) 'determinant from linsol ' , d
 IF(ABS(d).LT.tol10.and.has_to_print) then
   WRITE(STD_OUT,*) '**warning from linsol --',&
&     'determinant too small --',d
 endif
 If (present(det)) det=d
END SUBROUTINE linsol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    subroutine SolveAXeqB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SolveAXeqB(n,A,B,conditionNo)
  ! General purpose AX=B solver for A, B real.
  ! On return, B stores X
  INTEGER, INTENT(IN) :: n
  REAL(dp), INTENT(IN) :: A(:,:)
  REAL(dp), INTENT(INOUT) :: B(:)
  REAL(dp), INTENT(IN), OPTIONAL :: conditionNo
  REAL(dp), ALLOCATABLE :: C(:,:),U(:,:),VT(:,:),X(:)
  REAL(dp), ALLOCATABLE :: WORK(:)
  REAL(dp), ALLOCATABLE :: S(:)
  REAL(dp), PARAMETER :: rtol=1.d-9
  INTEGER :: i,LWORK
  REAL(dp) :: xx,tol
  REAL(dp), PARAMETER :: one=1,zero=0
  IF (n == 1) THEN
    B(1)=B(1)/A(1,1)
    RETURN
  ENDIF
  LWORK=MAX(200,n*n)
  LIBPAW_ALLOCATE(C,(n,n))
  LIBPAW_ALLOCATE(X,(n))
  LIBPAW_ALLOCATE(U,(n,n))
  LIBPAW_ALLOCATE(VT,(n,n))
  LIBPAW_ALLOCATE(WORK,(LWORK))
  LIBPAW_ALLOCATE(S,(n))
  tol=rtol
  IF (PRESENT(conditionNo)) tol=1.d0/conditionNo
  C(1:n,1:n)=A(1:n,1:n)
  CALL DGESVD('A','A',n,n,C,n,S,U,n,VT,n,WORK,LWORK,i)
  tol=tol*S(1)
  X=0
  DO i=1,n
    if(has_to_print) write(std_out,*) 'Solver, tol ',i,S(i),tol
    IF (S(i)>tol) THEN
      xx=DOT_PRODUCT(U(1:n,i),B(1:n))/S(i)
      X(1:n)=X(1:n)+xx*(VT(i,1:n))
    ENDIF
  ENDDO
  B=X
  LIBPAW_DEALLOCATE(C)
  LIBPAW_DEALLOCATE(X)
  LIBPAW_DEALLOCATE(U)
  LIBPAW_DEALLOCATE(VT)
  LIBPAW_DEALLOCATE(WORK)
  LIBPAW_DEALLOCATE(S) 
END SUBROUTINE SolveAXeqB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    subroutine SolveAXeqBM(n,A,B,many)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SolveAXeqBM(n,A,B,many)
 integer, intent(in) :: n,many
 real(dp), intent(in) :: A(:,:)
 real(dp), intent(inout) :: B(:)
 integer :: i,LWORK
 real(dp), parameter :: rtol=tol9
 real(dp), parameter :: one=1,zero=0
 real(dp) :: xx
 real(dp), allocatable :: C(:,:),U(:,:),VT(:,:),X(:)
 real(dp), allocatable :: WORK(:)
 real(dp), allocatable :: S(:)
 if (many<1.or.many>n) then
   LIBPAW_ERROR('Error in paw_SolveAXeqBM')
 endif
 if (n == 1) then
   B(1)=B(1)/A(1,1)
   return
 endif
 LWORK=max(200,n*n)
 LIBPAW_ALLOCATE(C,(n,n))
 LIBPAW_ALLOCATE(X,(n))
 LIBPAW_ALLOCATE(U,(n,n))
 LIBPAW_ALLOCATE(VT,(n,n))
 LIBPAW_ALLOCATE(WORK,(LWORK))
 LIBPAW_ALLOCATE(S,(n))
 C(1:n,1:n)=A(1:n,1:n)
 call DGESVD('A','A',n,n,C,n,S,U,n,VT,n,WORK,LWORK,i)
 X=0
 do i=1,n
   if (i<=many) then
     xx=DOT_PRODUCT(U(1:n,i),B(1:n))/S(i)
     X(1:n)=X(1:n)+xx*(VT(i,1:n))
   endif
 enddo
 B=X
 LIBPAW_DEALLOCATE(C)
 LIBPAW_DEALLOCATE(X)
 LIBPAW_DEALLOCATE(U)
 LIBPAW_DEALLOCATE(VT)
 LIBPAW_DEALLOCATE(WORK)
 LIBPAW_DEALLOCATE(S)
end subroutine SolveAXeqBM





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 9. Gridmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!  SUBROUTINE taufromwfn(otau,Grid,wfn,l,energy,rv)
!    input radial wfn and output its kinetic energy density
!     note that total wavefunction is wfn/r * Ylm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
SUBROUTINE taufromwfn(otau,Grid,wfn,l,energy,rPot)
 TYPE(GridInfo), INTENT(IN) :: Grid
 REAL(dp), INTENT(IN) :: wfn(:)
 REAL(dp), INTENT(IN), OPTIONAL :: rPot(:),energy
 INTEGER, INTENT(IN) :: l
 REAL(dp), INTENT(OUT) :: otau(:)
 LOGICAL,PARAMETER :: from_e=.false.
 INTEGER :: n
 REAL(dp) :: fac
 REAL(dp), allocatable:: pbr(:),dpdr(:),d2pdr(:)
 n=Grid%n
 fac=l*(l+1)
!Note: Psi(r) = Yl.Wfn(r)/r
!        Tau(r) given in Rydberg units
!Note: there are several kinetic energy formulas
!      differing by something*Laplacian(rho(r))
 IF (.not.from_e) THEN
!  Standard tau formula:
!  4pir^2* Tau(r) = [r.d/dr(Wfn/r)]^2 + l(l+1) [Wfn/r]^2
   LIBPAW_ALLOCATE(pbr,(n))
   LIBPAW_ALLOCATE(dpdr,(n))
   dpdr=0._dp;pbr=0._dp
   pbr(2:n)=wfn(2:n)/Grid%r(2:n)
   CALL derivative(Grid,pbr,dpdr,2,n)
   otau(2:n)=(Grid%r(2:n)*dpdr(2:n))**2 + fac*pbr(2:n)**2
   CALL extrapolate(otau) ; if (l>0) otau(1)=0._dp
   LIBPAW_DEALLOCATE(pbr)
   LIBPAW_DEALLOCATE(dpdr)
 ELSEIF (from_e) THEN
   IF (.NOT.(PRESENT(energy).AND.PRESENT(rPot))) then
     STOP 'Error in taufromwfn: rPot and energy should be present!'
   END IF
!  For testing purpose:
!  Another formula for the kinetic energy density
!  4pir^2* Tau(r) = [Eigenvalue - Veff(r)]*Wfn^2      
   IF (.TRUE.) THEN
     otau(2:n)=(energy-rPot(2:n)/Grid%r(2:n))*wfn(2:n)**2
     CALL extrapolate(otau) ; if (l>0) otau(1)=0._dp
!  Another one:
!  From energy + correction
   ELSE
     LIBPAW_ALLOCATE(pbr,(n))
     LIBPAW_ALLOCATE(dpdr,(n))
     LIBPAW_ALLOCATE(d2pdr,(n))
     pbr(2:n)=wfn(2:n)/Grid%r(2:n)
     CALL derivative(Grid,pbr,dpdr,2,n)
     otau(2:n)=dpdr(2:n)**2 *Grid%r(2:n)**2
     d2pdr(2:n)=pbr(2:n)**2
     call derivative(Grid,d2pdr,dpdr,2,n)
     call derivative(Grid,dpdr,d2pdr,2,n)
     d2pdr(2:n)=d2pdr(2:n)+2.*dpdr(2:n)/Grid%r(2:n)
     otau(2:n)=0.5_dp*d2pdr(2:n)*Grid%r(2:n)**2 &
&             +(energy-rPot(2:n)/Grid%r(2:n))*wfn(2:n)**2
     otau(2:n)=otau(2:n) + fac*pbr(2:n)**2
     CALL extrapolate(otau) ; if (l>0) otau(1)=0._dp
     LIBPAW_DEALLOCATE(pbr)
     LIBPAW_DEALLOCATE(dpdr)
     LIBPAW_DEALLOCATE(d2pdr)
   ENDIF
 ENDIF
 call filter(n,otau,machine_zero)
END SUBROUTINE taufromwfn  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  deltakinetic_ij(Grid,wfn1,wfn2,twfn1,twfn2,l,ekin,last)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE deltakinetic_ij(Grid,wfn1,wfn2,twfn1,twfn2,l,ekin,last)
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: wfn1(:),wfn2(:),twfn1(:),twfn2(:)
 INTEGER, INTENT(IN) :: l
 real(dp), INTENT(OUT) :: ekin
 INTEGER, INTENT(IN), OPTIONAL :: last
 real(dp), ALLOCATABLE :: dfdr1(:),dfdr2(:),arg1(:),arg2(:)
 real(dp), ALLOCATABLE :: tdfdr1(:),tdfdr2(:),targ1(:),targ2(:)
 INTEGER :: i,n
 n=Grid%n
 if (present(last)) n=last
 LIBPAW_ALLOCATE(dfdr1,(n))
 LIBPAW_ALLOCATE(arg1,(n))
 LIBPAW_ALLOCATE(dfdr2,(n))
 LIBPAW_ALLOCATE(arg2,(n))
 LIBPAW_ALLOCATE(tdfdr1,(n))
 LIBPAW_ALLOCATE(targ1,(n))
 LIBPAW_ALLOCATE(tdfdr2,(n))
 LIBPAW_ALLOCATE(targ2,(n))
 CALL derivative(Grid,wfn1,dfdr1,1,n)
 CALL derivative(Grid,wfn2,dfdr2,1,n)
 CALL derivative(Grid,twfn1,tdfdr1,1,n)
 CALL derivative(Grid,twfn2,tdfdr2,1,n)
 arg1=0; arg2=0; targ1=0; targ2=0
 DO i=2,n
   arg1(i)=wfn1(i)/Grid%r(i)
   arg2(i)=wfn2(i)/Grid%r(i)
   targ1(i)=twfn1(i)/Grid%r(i)
   targ2(i)=twfn2(i)/Grid%r(i)
 ENDDO
 DO i=1,n
   arg1(i)=(dfdr1(i)*dfdr2(i)-tdfdr1(i)*tdfdr2(i))&
&          +(l*(l+1))*(arg1(i)*arg2(i)-targ1(i)*targ2(i))
 ENDDO
 ekin=integrator(Grid,arg1,1,n)
 LIBPAW_DEALLOCATE(dfdr1)
 LIBPAW_DEALLOCATE(arg1)
 LIBPAW_DEALLOCATE(dfdr2)
 LIBPAW_DEALLOCATE(arg2)
 LIBPAW_DEALLOCATE(tdfdr1)
 LIBPAW_DEALLOCATE(targ1)
 LIBPAW_DEALLOCATE(tdfdr2)
 LIBPAW_DEALLOCATE(targ2)
END SUBROUTINE deltakinetic_ij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE atompaw_poisson(Grid,q,den,rv,ecoul,v00)
!!  Use Numerov algorithm to solve poisson equation
!!  den(n) is electron density * (4*pi*r**2)
!!  rv(n) is returned as electrostatic potential * r
!!  ecoul is the coulomb interaction energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE atompaw_poisson(Grid,q,den,rv,ecoul,v00)
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN):: den(:)
 real(dp), INTENT(INOUT) :: rv(:),ecoul,q
 real(dp), OPTIONAL, INTENT(OUT) :: v00
 real(dp), ALLOCATABLE :: a(:),b(:)
 real(dp) :: sd,sl,h,h2
 INTEGER :: i,n
 n=Grid%n
 h=Grid%h
 rv=0._dp
 q=integrator(Grid,den)
 LIBPAW_ALLOCATE(a,(n))
 LIBPAW_ALLOCATE(b,(n))
 IF (Grid%type==lineargrid) THEN
   sd=2
   sl=-1
   a(1)=0
   DO i=2,n
     a(i)=h*den(i)/(6*(i-1))
   ENDDO
   rv(1)=0
   rv(2)=10*a(2)+a(3)
   DO i=3,n-1
     rv(i)=10*a(i)+a(i+1)+a(i-1)
   ENDDO
   rv(n)=10*a(n)+a(n-1)+2*q
 ELSEIF (Grid%type==loggrid) THEN
   sd=2+10*h*h/48
   sl=-1+h*h/48
   a(1)=0
   h2=h*h
   DO i=2,n
     a(i)=h2*Grid%rr02(i)*den(i)/(6*Grid%r(i))/Grid%pref(i)
   ENDDO
   rv(1)=0
   rv(2)=10*a(2)+a(3)
   DO i=3,n-1
     rv(i)=10*a(i)+a(i+1)+a(i-1)
   ENDDO
   !   last term is boundary value at point n+1
   rv(n)=10*a(n)+a(n-1)-2*q*sl/(Grid%pref(n)*EXP(h/2))
 ENDIF
 CALL conthomas(n,sl,sd,rv)
 IF (Grid%type==loggrid) rv=rv*Grid%pref
 !  calculate ecoul
 DO i=2,n
   a(i)=den(i)*rv(i)/Grid%r(i)
 ENDDO
 a(1)=0
 ecoul=integrator(Grid,a)*0.5_dp
 IF (PRESENT(v00)) THEN
   a=0
   a(2:n)=den(2:n)/Grid%r(2:n)
   v00=2*integrator(Grid,a)
 ENDIF
 LIBPAW_DEALLOCATE(a)
 LIBPAW_DEALLOCATE(b)
END SUBROUTINE atompaw_poisson


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE poisson_marc(Grid,q,den,rv,ecoul)
!!  use Numerov algorithm to solve poisson equation
!!  den(n) is electron density * (4*pi*r**2)
!!  rv(n) is returned as electrostatic potential * r
!!  ecoul is the coulomb interation energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE poisson_marc(Grid,q,den,rv,ecoul)
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN):: den(:)
 real(dp), INTENT(INOUT) :: rv(:),ecoul,q
 real(dp), ALLOCATABLE :: aa(:),bb(:),cc(:),dd(:)
 real(dp) :: h
 INTEGER :: i,n,ir,jr
 n=Grid%n
 h=Grid%h
 rv=0._dp
 q=integrator(Grid,den)
 LIBPAW_ALLOCATE(aa,(n))
 LIBPAW_ALLOCATE(bb,(n))
 LIBPAW_ALLOCATE(cc,(n))
 LIBPAW_ALLOCATE(dd,(n))
 DO jr=n,2,-1
   ir=n-jr+1
   aa(ir)=den(jr)*Grid%drdu(jr)
   bb(ir)=den(jr)*Grid%drdu(jr)/Grid%r(jr)
 END DO
 cc=0._dp
 cc(5)=aa(n-4);cc(4)=aa(n-3);cc(3)=aa(n-2);cc(2)=aa(n-1)
 cc(1)=cc(4)+3._dp*(cc(2)-cc(3)) !call extrapolate(Grid,cc)
 aa(n)=cc(1)
 cc(5)=bb(n-4);cc(4)=bb(n-3);cc(3)=bb(n-2);cc(2)=bb(n-1)
 cc(1)=cc(4)+3._dp*(cc(2)-cc(3)) !call extrapolate(Grid,cc)
 bb(n)=cc(1)
 cc(1)=0._dp;dd(1)=0._dp
 DO ir=3,n,2
   cc(ir)  =cc(ir-2)+h/3._dp*(aa(ir-2)+4._dp*aa(ir-1)+aa(ir))
   cc(ir-1)=cc(ir-2)+h/3._dp*(1.25_dp*aa(ir-2)+2.0_dp*aa(ir-1)-0.25_dp*aa(ir))
   dd(ir)  =dd(ir-2)+h/3._dp*(bb(ir-2)+4._dp*bb(ir-1)+bb(ir))
   dd(ir-1)=dd(ir-2)+h/3._dp*(1.25_dp*bb(ir-2)+2._dp*bb(ir-1)-0.25_dp*bb(ir))
 END DO
 IF (MOD(n,2)==0) THEN
   cc(n)=cc(n-1)+h/3._dp*(1.25_dp*aa(n-2)+2._dp*aa(n-1)-0.25_dp*aa(n))
   dd(n)=dd(n-1)+h/3._dp*(1.25_dp*bb(n-2)+2._dp*bb(n-1)-0.25_dp*bb(n))
 END IF
 rv(1)=0._dp
 DO ir=2,n
   jr=n-ir+1
   rv(ir)=2._dp*(dd(jr)*Grid%r(ir)+(cc(n)-cc(jr))) !Ha->Ry
 END DO
 if (n<Grid%n) rv(n+1:Grid%n)=rv(n)
 !  calculate ecoul
 aa(1)=0._dp
 do i=2,n
   aa(i)=den(i)*rv(i)/Grid%r(i)
 end do
 ecoul=0.5_dp*integrator(Grid,aa)
 LIBPAW_DEALLOCATE(aa)
 LIBPAW_DEALLOCATE(bb)
 LIBPAW_DEALLOCATE(cc)
 LIBPAW_DEALLOCATE(dd)
END SUBROUTINE poisson_marc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  zeropot(Grid,rv,v0,v0p)
!!    extrapolate potential to value at r=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE zeropot(Grid,rv,v0,v0p)
 TYPE (GridInfo), INTENT(IN):: Grid
 real(dp), INTENT(IN) :: rv(:)    ! Note: rv(1) corresponds to r=0
 real(dp), INTENT(OUT) :: v0,v0p
 real(dp) :: tmp(15),tmp1(15)
 tmp(2:15)=rv(2:15)/Grid%r(2:15)
 CALL extrapolate(tmp(1:15))
 v0=tmp(1)
 CALL derivative(Grid,tmp(1:15),tmp1(1:15),2,15)
 CALL extrapolate(tmp1(1:15))
 v0p=tmp1(1)
END SUBROUTINE zeropot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  subroutine kinetic(Grid,wfn,l,ekin)
!!       calculates expectation value of kinetic energy for wfn
!!        with orbital angular momentum l
!!        wfn == r*radialwfn in Schroedinger Equation
!!        assumes wfn=(constant)*r^(l+1) at small r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE kinetic(Grid,wfn,l,ekin)
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: wfn(:)
 INTEGER, INTENT(IN) :: l
 real(dp), INTENT(OUT) :: ekin
 real(dp), ALLOCATABLE :: dfdr(:),arg(:)
 INTEGER :: i,n
 n=Grid%n
 LIBPAW_ALLOCATE(dfdr,(n))
 LIBPAW_ALLOCATE(arg,(n))
 CALL derivative(Grid,wfn,dfdr)
 arg=0
 DO i=2,n
   arg(i)=wfn(i)/Grid%r(i)
 ENDDO
 DO i=1,n
   arg(i)=(dfdr(i))**2+(l*(l+1))*(arg(i))**2
 ENDDO
 ekin=integrator(Grid,arg)!*0.5!Hartree
 LIBPAW_DEALLOCATE(dfdr)
 LIBPAW_DEALLOCATE(arg)
END SUBROUTINE kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  subroutine altkinetic(Grid,wfn,energy,rv,ekin)
!!       calculates expectation value of kinetic energy for wfn
!!        with orbital wfn by integrating
!!          int(wfn**2 * (energy-rv/r), r=0..rmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE altkinetic(Grid,wfn,energy,rv,ekin)
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: wfn(:),rv(:),energy
 real(dp), INTENT(OUT) :: ekin
 real(dp), ALLOCATABLE :: arg(:)
 INTEGER :: i,n
 n=Grid%n
 LIBPAW_ALLOCATE(arg,(n))
 arg=0
 DO i=2,n
   arg(i)=(wfn(i)**2)*(energy-rv(i)/Grid%r(i))
 ENDDO
 ekin=integrator(Grid,arg)
 LIBPAW_DEALLOCATE(arg)
END SUBROUTINE altkinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! overint(n,h,f1,icorr)
!!    function to calculate the integral of one vectors f1
!!      using simpsons rule assuming a regular grid with
!!      spacing of h and n total points
!!      icorr: optional parameter: used only when n is even
!!             if icorr<0,  a trapezoidal correction is applied
!!                          at the start of interval
!!             if icorr>=0, a trapezoidal correction is applied
!!                          at the end of interval
!!             default (if missing) is icorr=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION overint(n,h,f1,icorr)
 real(dp) :: overint
 INTEGER, INTENT(IN) :: n
 real(dp), INTENT(IN) :: h,f1(:)
 INTEGER, OPTIONAL :: icorr
 real(dp),PARAMETER :: tol=1.D-14
 INTEGER :: i,j,istart,m
 overint=0
 !Eliminate zeros at end of interval
 i=n;DO WHILE(ABS(f1(i))<machine_zero.AND.i>2);i=i-1;ENDDO
 m=MIN(i+1,n)
 IF (m<=1) THEN
   RETURN
 ELSEIF (m==2) THEN
   overint=(f1(1)+f1(2))*(h/2)   ! Trapezoidal rule
   RETURN
 ENDIF
 istart=1
 IF (PRESENT(icorr)) THEN
   IF (icorr<0.AND.MOD(m,2)==0) istart=2
 ENDIF
 overint=f1(istart)+4*f1(istart+1)+f1(istart+2)
 j=((m-istart)/2)*2+istart
 IF (j>=istart+4) THEN
   DO i=istart+4,j,2
     overint=overint+f1(i-2)+4*f1(i-1)+f1(i)
   ENDDO
 ENDIF
 overint=overint*(h/3)
 IF (m>j) overint=overint+(f1(j)+f1(m))*(h/2)
 IF (istart==2) overint=overint+(f1(1)+f1(2))*(h/2)
 RETURN
END FUNCTION overint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! function integrator(Grid,arg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION integrator(Grid,arg,str,fin)
 real(dp) :: integrator
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: arg(:)
 INTEGER, INTENT(IN), OPTIONAL :: str,fin
 real(dp), ALLOCATABLE :: dum(:)
 INTEGER :: n,i1,i2
 n=Grid%n
 i1=1;i2=n
 IF (PRESENT(str).AND.PRESENT(fin)) THEN
   i1=str; i2=fin; n=i2-i1+1
 ENDIF
 SELECT CASE(Grid%type)
 CASE default
   LIBPAW_ERROR('Error in integrator')
 CASE(lineargrid)
   integrator=overint(n,Grid%h,arg(i1:i2))
 CASE(loggrid)
   LIBPAW_BOUND1_ALLOCATE(dum,BOUNDS(i1,i2))
   dum(i1:i2)=arg(i1:i2)*Grid%drdu(i1:i2)
   integrator=overint(n,Grid%h,dum(i1:i2),-1)
   LIBPAW_DEALLOCATE(dum)
 END SELECT
END FUNCTION integrator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! function FindGridIndex(Grid,rpoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION FindGridIndex(Grid,rpoint)
 INTEGER :: FindGridIndex
 TYPE (GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: rpoint
 real(dp) :: r0
 FindGridIndex=0
 IF (Grid%type==lineargrid) THEN
   FindGridIndex=rpoint/Grid%h+1
   IF (Grid%h*(FindGridIndex-1)<rpoint-tol10) FindGridIndex=FindGridIndex+1
 ELSEIF (Grid%type==loggrid) THEN
   r0=Grid%drdu(1)
   FindGridIndex=LOG(rpoint/r0+1)/Grid%h+1
   IF (r0*EXP(Grid%h*(FindGridIndex-1))<rpoint-tol10) &
&         FindGridIndex=FindGridIndex+1
 ENDIF
END FUNCTION FindGridIndex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ClassicalTurningPoint(Grid,rv,l,energy,turningpoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ClassicalTurningPoint(Grid,rv,l,energy,turningpoint)
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: rv(:)
 INTEGER, INTENT(IN) :: l
 real(dp), INTENT(IN) :: energy
 INTEGER, INTENT(OUT) :: turningpoint
 INTEGER :: i,n
 real(dp), ALLOCATABLE :: v(:)
 n=Grid%n
 LIBPAW_ALLOCATE(v,(n))
 v=0
 v(2:n)=rv(2:n)/Grid%r(2:n)+l*(l+1)/(Grid%r(2:n)**2)!hartree
 turningpoint=n
 DO i=n,2,-1
   IF (v(i)<energy) EXIT
 ENDDO
 turningpoint=i
 turningpoint=MIN(turningpoint,FindGridIndex(Grid,10.0_dp))
 LIBPAW_DEALLOCATE(v)
END SUBROUTINE ClassicalTurningPoint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! getwfnfromcfdsol(start,finish,yy,wfn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getwfnfromcfdsol(start,finish,yy,wfn)
 INTEGER, INTENT(IN) :: start,finish
 real(dp), INTENT(IN) :: yy(:,:)
 real(dp), INTENT(INOUT) :: wfn(:)
 INTEGER :: i
 wfn=0
 do i=start,finish
   wfn(i)=yy(1,i)
 enddo
end subroutine getwfnfromcfdsol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! countnodes(start,finish,wfn,filter)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER function countnodes(start,finish,wfn,filter)
 INTEGER, INTENT(IN) :: start,finish
 real(dp), INTENT(IN) :: wfn(:)
 real(dp), INTENT(IN), OPTIONAL :: filter
 INTEGER :: i,nodes
 nodes=0
 do i=start+1,finish
   if (wfn(i)*wfn(i-1)<0._dp) nodes=nodes+1
   if (PRESENT(filter)) then
     If((abs(wfn(i))+abs(wfn(i-1)))<filter) exit
   endif
 enddo
 countnodes=nodes
end function countnodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! extrapolate(Grid,v)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE extrapolate(v)
 ! extrapolate array v to r=0 at v(1)
 real(dp), INTENT(INOUT) :: v(:)  ! assume v(2),v(3)...  given
 v(1)=5._dp*v(2)-10._dp*v(3)+10._dp*v(4)-5._dp*v(5)+v(6) ! fourth order formula
END SUBROUTINE extrapolate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine derivative(Grid,f,dfdr,begin,bend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE derivative(Grid,f,dfdr,begin,bend)
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: f(:)
 real(dp), INTENT(OUT) :: dfdr(:)
 INTEGER, OPTIONAL, INTENT(IN) :: begin,bend
 INTEGER :: i,n,i1,i2
 i1=1;i2=Grid%n;n=i2-i1+1
 IF (PRESENT(begin).OR.PRESENT(bend)) THEN
   IF (begin>=1.AND.bend<= Grid%n) THEN
     i1=begin;i2=bend;n=i2-i1+1
   ELSE
     LIBPAW_ERROR('Error in derivative')
   ENDIF
 ENDIF
 SELECT CASE(Grid%type)
 CASE default
   LIBPAW_ERROR('Error in derivative')
 CASE(lineargrid)
   CALL nderiv(Grid%h,f(i1:i2),dfdr(i1:i2),n,i)
   IF (i/=0) THEN
     LIBPAW_ERROR('Error in derivative -nderiv problem')
   ENDIF
 CASE(loggrid)
   CALL nderiv(Grid%h,f(i1:i2),dfdr(i1:i2),n,i)
   IF (i/=0) THEN
     LIBPAW_ERROR('Error in derivative -nderiv problem')
   ENDIF
   dfdr(i1:i2)=dfdr(i1:i2)/Grid%drdu(i1:i2)
 END SELECT
END SUBROUTINE derivative


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! nderiv(h,y,z,ndim,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE nderiv(h,y,z,ndim,ierr)
  INTEGER, INTENT(IN) :: ndim
  INTEGER, INTENT(INOUT) :: ierr
  real(dp) , INTENT(IN) :: h,y(:)
  real(dp) , INTENT(INOUT) :: z(:)
  real(dp) :: hh,yy,a,b,c
  INTEGER :: i
  ierr=-1
  IF (ndim.LT.5) RETURN
  !        prepare differentiation loop
  hh=.08333333333333333_dp/h
  yy=y(ndim-4)
  b=hh*(-25._dp*y(1)+48._dp*y(2)-36._dp*y(3)+16._dp*y(4)-3._dp*y(5))
  c=hh*(-3._dp*y(1)-10._dp*y(2)+18._dp*y(3)-6._dp*y(4)+y(5))
  !        start differentiation loop
  DO  i=5,ndim
    a=b
    b=c
    c=hh*(y(i-4)-y(i)+8._dp*(y(i-1)-y(i-3)))
    z(i-4)=a
  ENDDO
  !        end of differentiation loop
  !        normal exit
  a=hh*(-yy+6._dp*y(ndim-3)-18._dp*y(ndim-2)+10._dp*y(ndim-1)          &
&      +3._dp*y(ndim))
  z(ndim)=hh*(3._dp*yy-16._dp*y(ndim-3)+36._dp*y(ndim-2)               &
&      -48._dp*y(ndim-1)+25._dp*y(ndim))
  z(ndim-1)=a
  z(ndim-2)=c
  z(ndim-3)=b
  ierr=0
  RETURN
 END SUBROUTINE nderiv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Second derivative for general grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Gsecondderiv(Grid,index,g)
 real(dp) :: Gsecondderiv
 TYPE (GridInfo), INTENT(IN) :: Grid
 INTEGER, INTENT(IN) :: index
 real(dp), INTENT(IN) :: g(:)
 Gsecondderiv=0
 IF (Grid%type==lineargrid) THEN
   Gsecondderiv=secondderiv(index,g,Grid%h)
 ELSEIF  (Grid%type==loggrid) THEN
   Gsecondderiv=(secondderiv(index,g,Grid%h)&
&            -firstderiv(index,g,Grid%h))/Grid%rr02(index)
 ENDIF
END FUNCTION Gsecondderiv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! First  derivative for general grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Gfirstderiv(Grid,index,g)
 real(dp) :: Gfirstderiv
 TYPE (GridInfo), INTENT(IN) :: Grid
 INTEGER, INTENT(IN) :: index
 real(dp), INTENT(IN) :: g(:)
 Gfirstderiv=0
 IF (Grid%type==lineargrid) THEN
   Gfirstderiv=firstderiv(index,g,Grid%h)
 ELSEIF  (Grid%type==loggrid) THEN
   Gfirstderiv=firstderiv(index,g,Grid%h)/Grid%drdu(index)
 ENDIF
END FUNCTION Gfirstderiv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! finite difference second derivative
!!   based on 5 point formula
!!   Ref. Engeln-Mullges & Uhlig (1996)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION secondderiv(index,f,h)
 real(dp) :: secondderiv
 INTEGER, INTENT(IN) :: index
 real(dp), INTENT(IN) :: f(:),h
 INTEGER :: n
 n=SIZE(f)
 secondderiv=0
 if (index==1.and.n>=5) THEN
   secondderiv=(70*f(1)-208*f(2)+228*f(3)-112*f(4)+22*f(5))/(24*h*h)
 else if (index==2.and.n>=5) THEN
   secondderiv=(22*f(1)-40*f(2)+12*f(3)+8*f(4)-2*f(5))/(24*h*h)
 else if (index>2.and.index<=n-2) THEN
   secondderiv=-(f(index-2)+f(index+2))/12 + &
&       4*(f(index-1)+f(index+1))/3 - 5*f(index)/2
   secondderiv=secondderiv/(h*h)
 else if (index>=5.and.index==n-1)   THEN
   secondderiv=(-2*f(n-4)+8*f(n-3)+12*f(n-2)-40*f(n-1)+22*f(n))/(24*h*h)
 else if (index>=5.and.index==n)   THEN
   secondderiv=(22*f(n-4)-112*f(n-3)+228*f(n-2)-208*f(n-1)+70*f(n))/(24*h*h)
 else
   LIBPAW_ERROR('Error in secondderiv')
 ENDIF
END FUNCTION secondderiv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! finite difference first derivative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION firstderiv(index,f,h)
 real(dp) :: firstderiv
 INTEGER, INTENT(IN) :: index
 real(dp), INTENT(IN) :: f(:),h
 INTEGER :: n
 n=SIZE(f)
 firstderiv=0
 if (index==1.and.n>=5) THEN
   firstderiv=(-25*f(1)+48*f(2)-36*f(3)+16*f(4)-3*f(5))/(12*h)
 else if (index==2.and.n>=5) THEN
   firstderiv=(-3*f(1)-10*f(2)+18*f(3)-6*f(4)+f(5))/(12*h)
 else if (index>2.and.index<=n-2) THEN
   firstderiv=(f(index-2)-8*f(index-1)+8*f(index+1)-f(index+2))/(12*h)
 else if (index>=5.and.index==n-1)   THEN
   firstderiv=(-f(n-4)+6*f(n-3)-18*f(n-2)+10*f(n-1)+3*f(n))/(12*h)
 else if (index>=5.and.index==n)   THEN
   firstderiv=(3*f(n-4)-16*f(n-3)+36*f(n-2)-48*f(n-1)+25*f(n))/(12*h)
 else
   LIBPAW_ERROR('Error in firstderiv')
 ENDIF
END FUNCTION firstderiv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  function to calculate the overlap between two vectors f1 and f2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION overlap(Grid,f1,f2,str,fin)
 real(dp) :: overlap
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: f1(:),f2(:)
 INTEGER, INTENT(IN), OPTIONAL :: str,fin
 real(dp), ALLOCATABLE :: dum(:)
 INTEGER :: n,i1,i2
 n=Grid%n
 i1=1;i2=n
 IF (PRESENT(str).AND.PRESENT(fin)) THEN
   i1=str; i2=fin; n=i2-i1+1
 ENDIF
 LIBPAW_BOUND1_ALLOCATE(dum,BOUNDS(i1,i2))
 dum(1:n)=f1(i1:i2)*f2(i1:i2)
 overlap=integrator(Grid,dum(1:n),1,n)
 LIBPAW_DEALLOCATE(dum)
END FUNCTION overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine from David Vanderbilt's USPS code, modified by Marc
!!     Torrent and Francois Jollet, further modified by NAWH
!!===========================================================================
!!      subroutine cfdsol(zz,yy,jj1,jj2,mesh)
!!===========================================================================
!!     routine for solving coupled first order differential equations
!!
!!      d yy(x,1)
!!      ---------   =  zz(x,1,1) * yy(x,1) + zz(x,1,2) * yy(2,1)
!!         dx
!!
!!      d yy(x,2)
!!      ---------   =  zz(x,2,1) * yy(x,1) + zz(x,2,2) * yy(2,1)
!!         dx
!!
!!
!!     using fifth order predictor corrector algorithm
!!
!!     routine integrates from jj1 to jj2 and can cope with both cases
!!     jj1 < jj2 and jj1 > jj2.  first five starting values of yy must
!!     be provided by the calling program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE cfdsol(Grid,zz,yy,jj1,jj2)
 TYPE(gridinfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN):: zz(:,:,:)
 real(dp), INTENT(INOUT):: yy(:,:)
 INTEGER, INTENT(IN)  :: jj1,jj2
 real(dp):: fa(0:5),fb(0:5),abp(1:5),amc(0:4)
 INTEGER :: isgn,i,j,ip,mesh
 real(dp):: arp,brp
 real(dp), ALLOCATABLE :: tmpz(:,:,:)
 real(dp), PARAMETER :: verylarge=1.d30
 real(dp) :: scale
 mesh=SIZE(yy(2,:))
 IF (SIZE(zz(2,2,:))/=mesh) THEN
   LIBPAW_ERROR('cfdsol error - incompatible arrays')
 ENDIF
 isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
 IF ( isgn .EQ. + 1 ) THEN
   IF ( jj1 .LE. 5 .OR. jj2 .GT. mesh ) THEN
     LIBPAW_ERROR(' ***error in subroutine difsol')
   ENDIF
 ELSEIF ( isgn .EQ. - 1 ) THEN
   IF ( jj1 .GE. ( mesh - 4 ) .OR. jj2 .LT. 1 ) THEN
     LIBPAW_ERROR(' ***error in subroutine difsol')
   ENDIF
 ELSE
   if(has_to_print) write(std_out,*) isgn,jj1,jj2,mesh
 ENDIF
 LIBPAW_ALLOCATE(tmpz,(2,2,mesh))
 tmpz=zz
 DO i=1,2
   DO j=1,2
     tmpz(i,j,:)=tmpz(i,j,:)*Grid%h
     if (Grid%TYPE==loggrid) tmpz(i,j,1:mesh)=tmpz(i,j,1:mesh)*Grid%drdu(1:mesh)
   ENDDO
 ENDDO
 abp(1) = 1901._dp / 720._dp
 abp(2) = -1387._dp / 360._dp
 abp(3) = 109._dp / 30._dp
 abp(4) = -637._dp / 360._dp
 abp(5) = 251._dp / 720._dp
 amc(0) = 251._dp / 720._dp
 amc(1) = 323._dp / 360._dp
 amc(2) = -11._dp / 30._dp
 amc(3) = 53._dp / 360._dp
 amc(4) = -19._dp / 720._dp
 DO j = 1,5
   ip = jj1 - isgn * j
   fa(j) = tmpz(1,1,ip) * yy(1,ip) + tmpz(1,2,ip) * yy(2,ip)
   fb(j) = tmpz(2,1,ip) * yy(1,ip) + tmpz(2,2,ip) * yy(2,ip)
 ENDDO
 DO j = jj1,jj2,isgn
   arp = yy(1,j-isgn)
   brp = yy(2,j-isgn)
   IF (ABS(arp)>verylarge.OR.brp>verylarge) THEN
     scale=1._dp/(ABS(arp)+ABS(brp))
     arp=arp*scale
     brp=brp*scale
     fa(:)=fa(:)*scale; fb(:)=fb(:)*scale
     yy=yy*scale
   ENDIF
   DO  i = 1,5
     arp = arp + DBLE(isgn) * abp(i) * fa(i)
     brp = brp + DBLE(isgn) * abp(i) * fb(i)
   ENDDO
   fa(0) = tmpz(1,1,j) * arp + tmpz(1,2,j) * brp
   fb(0) = tmpz(2,1,j) * arp + tmpz(2,2,j) * brp
   yy(1,j) = yy(1,j-isgn)
   yy(2,j) = yy(2,j-isgn)
   DO  i = 0,4,1
     yy(1,j) = yy(1,j) + DBLE(isgn) * amc(i) * fa(i)
     yy(2,j) = yy(2,j) + DBLE(isgn) * amc(i) * fb(i)
   ENDDO
   DO i = 5,2,-1
     fa(i) = fa(i-1)
     fb(i) = fb(i-1)
   ENDDO
   fa(1) = tmpz(1,1,j) * yy(1,j) + tmpz(1,2,j) * yy(2,j)
   fb(1) = tmpz(2,1,j) * yy(1,j) + tmpz(2,2,j) * yy(2,j)
 ENDDO
 LIBPAW_DEALLOCATE(tmpz)
END SUBROUTINE cfdsol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE cfdsoliter(Grid,zz,yy,jj1,jj2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE cfdsoliter(Grid,zz,yy,jj1,jj2)
 TYPE(gridinfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN):: zz(:,:,:)
 real(dp), INTENT(INOUT):: yy(:,:)
 INTEGER, INTENT(IN)  :: jj1,jj2
 real(dp):: fa(0:5),fb(0:5),abp(1:5),amc(0:4),yprev(2),ycorr(2)
 INTEGER :: isgn,i,j,ip,mesh,k
 INTEGER, PARAMETER :: CORRITER=5
 real(dp):: arp,brp
 real(dp), ALLOCATABLE :: tmpz(:,:,:)
 real(dp), PARAMETER :: verylarge=10._dp, smallenough=tol5
 real(dp) :: scale,small
 mesh=SIZE(yy(2,:))
 IF (SIZE(zz(2,2,:))/=mesh) THEN
  LIBPAW_ERROR('cfdsol error - incompatible arrays')
 ENDIF
 isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
 IF ( isgn .EQ. + 1 ) THEN
   IF ( jj1 .LE. 5 .OR. jj2 .GT. mesh ) THEN
     LIBPAW_ERROR(' ***error in subroutine difsol')
   ENDIF
 ELSEIF ( isgn .EQ. - 1 ) THEN
   IF ( jj1 .GE. ( mesh - 4 ) .OR. jj2 .LT. 1 ) THEN
     LIBPAW_ERROR(' ***error in subroutine difsol')
   ENDIF
 ELSE
   if(has_to_print) WRITE(STD_OUT,*) isgn,jj1,jj2,mesh
 ENDIF
 LIBPAW_ALLOCATE(tmpz,(2,2,mesh))
 tmpz=zz
 DO i=1,2
   DO j=1,2
     tmpz(i,j,:)=tmpz(i,j,:)*Grid%h
     if (Grid%TYPE==loggrid) tmpz(i,j,1:mesh)=tmpz(i,j,1:mesh)*Grid%drdu(1:mesh)
   ENDDO
 ENDDO
 abp(1) = 1901._dp / 720._dp
 abp(2) = -1387._dp / 360._dp
 abp(3) = 109._dp / 30._dp
 abp(4) = -637._dp / 360._dp
 abp(5) = 251._dp / 720._dp
 amc(0) = 251._dp / 720._dp
 amc(1) = 323._dp / 360._dp
 amc(2) = -11._dp / 30._dp
 amc(3) = 53._dp / 360._dp
 amc(4) = -19._dp / 720._dp
 DO j = 1,5
   ip = jj1 - isgn * j
   fa(j) = tmpz(1,1,ip) * yy(1,ip) + tmpz(1,2,ip) * yy(2,ip)
   fb(j) = tmpz(2,1,ip) * yy(1,ip) + tmpz(2,2,ip) * yy(2,ip)
 ENDDO
 DO j = jj1,jj2,isgn
   arp = yy(1,j-isgn)
   brp = yy(2,j-isgn)
   IF (ABS(arp)>verylarge.OR.brp>verylarge) THEN
     scale=1._dp/(ABS(arp)+ABS(brp))
     arp=arp*scale
     brp=brp*scale
     fa(:)=fa(:)*scale; fb(:)=fb(:)*scale
     yy=yy*scale
   ENDIF
   DO  i = 1,5
     arp = arp + DBLE(isgn) * abp(i) * fa(i)
     brp = brp + DBLE(isgn) * abp(i) * fb(i)
   ENDDO
   fa(0) = tmpz(1,1,j) * arp + tmpz(1,2,j) * brp
   fb(0) = tmpz(2,1,j) * arp + tmpz(2,2,j) * brp
   yprev(1) = arp
   yprev(2) = brp
   ycorr(1) = yy(1,j-isgn)
   ycorr(2) = yy(2,j-isgn)
   DO  i = 1,4,1
     ycorr(1) = ycorr(1) + DBLE(isgn) * amc(i) * fa(i)
     ycorr(2) = ycorr(2) + DBLE(isgn) * amc(i) * fb(i)
   ENDDO
   DO k=1,CORRITER
     yy(1,j)=ycorr(1) + DBLE(isgn) * amc(0) * fa(0)
     yy(2,j)=ycorr(2) + DBLE(isgn) * amc(0) * fb(0)
     small=abs(yprev(1))+abs(yprev(2))
     small=(abs(yprev(1)-yy(1,j))+abs(yprev(2)-yy(2,j)))/small
     if(small.le.smallenough) exit
     yprev(1)= yy(1,j)
     yprev(2)= yy(2,j)
     fa(0) = tmpz(1,1,j) * yprev(1) + tmpz(1,2,j) * yprev(2)
     fb(0) = tmpz(2,1,j) * yprev(1) + tmpz(2,2,j) * yprev(2)
   ENDDO 
   DO i = 5,2,-1
     fa(i) = fa(i-1)
     fb(i) = fb(i-1)
   ENDDO
   fa(1) = tmpz(1,1,j) * yy(1,j) + tmpz(1,2,j) * yy(2,j)
   fb(1) = tmpz(2,1,j) * yy(1,j) + tmpz(2,2,j) * yy(2,j)
 ENDDO
 LIBPAW_DEALLOCATE(tmpz)
END SUBROUTINE cfdsoliter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  subroutine findh_given_r0(Z,range,r0,n,hval)
!!    find hval for fixed number of input grid points n in loggrid case
!!    assumes form r(i)=(r0/Z)*(exp(h*(i-1))-1);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE findh_given_r0(Z,range,r0,n,hval)
 INTEGER, INTENT(IN) :: n
 real(dp), INTENT(IN) :: Z,range,r0
 real(dp), INTENT(INOUT) :: hval
 hval=log((Z*range/r0 + 1._dp))/(n-1)
end subroutine findh_given_r0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  subroutine findh(Z,range,n,hval,r0)
!!    find hval for fixed number of input grid points n in loggrid case
!!    assumes form r(i)=(h/Z)*(exp(h*(i-1))-1);   r0=h/Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE findh(Z,range,n,hval,r0)
 real(dp), INTENT(IN) :: Z
 INTEGER, INTENT(IN) :: n
 real(dp), INTENT(IN) :: range
 real(dp), INTENT(INOUT) :: hval,r0
 real(dp) :: h0,dh,f,df
 INTEGER :: i
 INTEGER, parameter :: iter=1000
 real(dp), parameter :: eps=1.e-15
 LOGICAL :: success
 h0=hval
 success=.false.
 do i=1,iter
   f=LOG(Z*range/h0+1._dp)/h0
   df=-f/h0-(Z*range/h0**3)/(Z*range/h0+1._dp)
   dh=(n-1-f)/df
   if (ABS(dh)< eps) then
     success=.true.
     exit
   endif
   if (h0+dh<0._dp) then
     h0=h0/2
   else
     h0=h0+dh
   endif
 enddo
 if (.not.success) then
   if(has_to_print) write(std_out,*) 'Warning in findh -- dh > eps ', dh,h0
 endif
 hval=h0
 r0=hval/Z
end subroutine findh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine initgrid(Grid,h,range,r0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE InitGrid(Grid,h,range,r0,do_not_print,pawrad)
 TYPE (GridInfo), INTENT(INOUT) :: Grid
 real(dp), INTENT(IN) :: range
 real(dp), INTENT(IN) :: h
 real(dp), OPTIONAL, INTENT(IN) :: r0
 LOGICAL, OPTIONAL, INTENT(IN) :: do_not_print
 TYPE(pawrad_type),optional, INTENT(IN) :: pawrad
 INTEGER :: i,n
 LOGICAL :: do_print
 do_print=.false.;if (present(do_not_print)) do_print=.not.do_not_print
 if (present(pawrad)) then
   if (pawrad%mesh_type/=2) then
     LIBPAW_ERROR("Error wrong pawrad grid")
     STOP
   end if
   Grid%n=pawrad%mesh_size
   LIBPAW_POINTER_ALLOCATE(Grid%r,(Grid%n))
   LIBPAW_POINTER_ALLOCATE(Grid%drdu,(Grid%n))
   LIBPAW_POINTER_ALLOCATE(Grid%pref,(Grid%n))
   LIBPAW_POINTER_ALLOCATE(Grid%rr02,(Grid%n))
   Grid%r(:)=pawrad%rad(:)
   Grid%drdu(:)=pawrad%radfact(:)
   Grid%rr02(:)=pawrad%radfact(:)**2
   Grid%pref(:)=(pawrad%radfact(:)*pawrad%rstep)**2
   Grid%r0=pawrad%rstep
   Grid%range=pawrad%rmax
   Grid%h=pawrad%lstep
   Grid%type=loggrid
   Grid%ishift=5
 else
   IF (PRESENT(r0)) THEN
     Grid%h=h
     Grid%r0=r0
     Grid%type=loggrid
     Grid%range=range
     n=LOG(range/r0+1)/h+1
     Grid%ishift=5
     IF (r0*(EXP(h*(n-1))-1)<range-tol5) n=n+1
     Grid%n=n
     if (do_print) write(std_out,*) 'InitGrid: -- logarithmic ',n, h,range,r0
     LIBPAW_POINTER_ALLOCATE(Grid%r,(n))
     LIBPAW_POINTER_ALLOCATE(Grid%drdu,(n))
     LIBPAW_POINTER_ALLOCATE(Grid%pref,(n))
     LIBPAW_POINTER_ALLOCATE(Grid%rr02,(n))
     DO i=1,n
       Grid%r(i)=r0*(EXP(Grid%h*(i-1))-1)
       Grid%drdu(i)=r0*EXP(Grid%h*(i-1))
       Grid%pref(i)=r0*EXP(Grid%h*(i-1)/2._dp)
       Grid%rr02(i)=(Grid%r(i)+r0)**2
     ENDDO
   ELSE
     Grid%h=h
     Grid%r0=0._dp
     Grid%type=lineargrid
     Grid%range=range
     n=range/h+1
     Grid%ishift=25
     IF (h*(n-1)<range-tol5) n=n+1
     Grid%n=n
     if (do_print) write(std_out,*) 'InitGrid: -- linear  ', n,h,range
     LIBPAW_POINTER_ALLOCATE(Grid%r,(n))
     LIBPAW_POINTER_ALLOCATE(Grid%drdu,(n))
     LIBPAW_POINTER_ALLOCATE(Grid%pref,(n))
     LIBPAW_POINTER_ALLOCATE(Grid%rr02,(n))
     DO i=1,n
       Grid%r(i)=(Grid%h*(i-1))
       Grid%drdu(i)=1._dp
       Grid%pref(i)=1._dp
       Grid%rr02(i)=1._dp
     ENDDO
   ENDIF
 endif
END SUBROUTINE InitGrid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine destroygrid(Grid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DestroyGrid(Grid)
 TYPE (GridInfo), INTENT(INOUT) :: Grid
 IF (ASSOCIATED(Grid%r)) then
   LIBPAW_POINTER_DEALLOCATE(Grid%r)
 endif
 IF (ASSOCIATED(Grid%drdu)) then
   LIBPAW_POINTER_DEALLOCATE(Grid%drdu)
 endif
 IF (ASSOCIATED(Grid%pref)) then
   LIBPAW_POINTER_DEALLOCATE(Grid%pref)
 endif
 IF (ASSOCIATED(Grid%rr02)) then
   LIBPAW_POINTER_DEALLOCATE(Grid%rr02)
 endif
END SUBROUTINE DestroyGrid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine nullifygrid(Grid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NullifyGrid(Grid)
 TYPE (GridInfo), INTENT(INOUT) :: Grid
 NULLIFY(Grid%r)
 NULLIFY(Grid%drdu)
 NULLIFY(Grid%pref)
 NULLIFY(Grid%rr02)
END SUBROUTINE NullifyGrid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     SUBROUTINE backward_numerov(Grid,l,match,energy,rv,wfn,wgt,nend,proj)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE backward_numerov(Grid,l,match,energy,rv,wfn,wgt,nend,proj)
 TYPE (GridInfo), INTENT(IN) :: Grid
 INTEGER, INTENT(IN) :: l,match
 real(dp), INTENT(IN) :: energy,rv(:)
 real(dp), INTENT(INOUT) :: wfn(:)    ! on input wfn(n-1) and wfn(n) given
 real(dp), INTENT(IN), OPTIONAL :: wgt(:)
 real(dp), intent(in), optional :: proj(:)
 integer, intent(in), optional :: nend
 real(dp), ALLOCATABLE :: a(:),b(:),p(:),c(:)
 real(dp) :: angm,h,h2,scale
 real(dp), PARAMETER :: vlarg=1.d30
 INTEGER :: i,n
 LOGICAL :: withwgt
 withwgt=.false.
 If (PRESENT(wgt)) withwgt=.true.
 n=Grid%n
 if(present(nend)) n=nend
 LIBPAW_ALLOCATE(a,(n))
 LIBPAW_ALLOCATE(b,(n))
 LIBPAW_ALLOCATE(p,(n))
 LIBPAW_ALLOCATE(c,(n))
 p(n)=wfn(n)
 p(n-1)=wfn(n-1)
 angm=l*(l+1)
 a=0
 b=0
 c=0
 h=Grid%h;    h2=h*h
 DO i=match,n
   if(withwgt) then
     a(i)=rv(i)/Grid%r(i)-energy*wgt(i)+angm/(Grid%r(i)**2)
   else        
     a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
   endif
 ENDDO
 if(present(proj)) b(match:n)=0.1_dp*h2*proj(match:n)
 IF (Grid%type==loggrid) THEN
   p(n-1:n)=p(n-1:n)/Grid%pref(n-1:n)
   a(match:n)=0.25_dp+Grid%rr02(match:n)*a(match:n)
   if(present(proj)) b(match:n)=Grid%rr02(match:n)*b(match:n)/Grid%pref(match:n)
 ENDIF
 if(present(proj)) then
   do i=match+1,n-1
     c(i)=10*b(i)+b(i-1)+b(i+1)
   enddo
 endif
 b(match:n)=2.4_dp+h2*a(match:n)
 a(match:n)=1.2_dp-0.1_dp*h2*a(match:n)
 DO i=n-2,match,-1
   p(i)=(b(i+1)*p(i+1)-a(i+2)*p(i+2)-c(i+1))/a(i)
   !renormalize if necessary
   scale=ABS(p(i))
   IF (scale > vlarg) THEN
     scale=1._dp/scale
     p(i:n)=scale*p(i:n)
   ENDIF
 ENDDO
 wfn(match:n)=p(match:n)
 IF (Grid%type==loggrid) THEN
   wfn(match:n)=wfn(match:n)*Grid%pref(match:n)
 ENDIF
 LIBPAW_DEALLOCATE(a)
 LIBPAW_DEALLOCATE(b)
 LIBPAW_DEALLOCATE(p)
 LIBPAW_DEALLOCATE(c)
END SUBROUTINE backward_numerov


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    SUBROUTINE
!forward_numerov(Grid,l,many,energy,rv,zeroval,wfn,nodes,wgt,proj,p3val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE forward_numerov(Grid,l,many,energy,rv,zeroval,wfn,nodes,wgt,proj,p3val)
 TYPE (GridInfo), INTENT(IN) :: Grid
 INTEGER, INTENT(IN) :: l,many
 real(dp), INTENT(IN) :: energy,zeroval,rv(:)
 real(dp), INTENT(INOUT) :: wfn(:)    ! on input wfn(1) and wfn(2) given
 INTEGER, INTENT(OUT) :: nodes
 real(dp), INTENT(IN), OPTIONAL :: wgt(:)
 real(dp), intent(in), optional :: p3val
 real(dp), intent(in), optional :: proj(:)
 real(dp), ALLOCATABLE :: a(:),b(:),p(:),c(:)
 real(dp) :: xx,angm,h,h2,scale
 real(dp), PARAMETER :: vlarg=1.d30
 INTEGER :: i
 LOGICAL :: withwgt
 withwgt=.false.
 If (PRESENT(wgt)) withwgt=.true.
 LIBPAW_ALLOCATE(a,(many))
 LIBPAW_ALLOCATE(b,(many))
 LIBPAW_ALLOCATE(p,(many))
 LIBPAW_ALLOCATE(c,(many))
 p(1)=wfn(1)
 p(2)=wfn(2)
 xx=zeroval
 angm=l*(l+1)
 a=0
 c=0
 h=Grid%h;    h2=h*h
 DO i=2,many
   if (withwgt) then
     a(i)=rv(i)/Grid%r(i)-energy*wgt(i)+angm/(Grid%r(i)**2)
   else
     a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
   endif  
 ENDDO
 if(present(proj)) b(1:many)=0.1_dp*h2*proj(1:many)
 IF (Grid%type==loggrid) THEN
   p(1:2)=wfn(1:2)/Grid%pref(1:2)
   xx=Grid%rr02(1)*xx/Grid%pref(1)
   a=0.25_dp+Grid%rr02(1:many)*a
   if(present(proj)) b(1:many)=Grid%rr02(1:many)*b(1:many)/Grid%pref(1:many)
 ENDIF
 if(present(proj)) then
   do i=2,many-1
     c(i)=10*b(i)+b(i-1)+b(i+1)
   enddo
 endif   
 b=2.4_dp+h2*a
 a=1.2_dp-0.1_dp*h2*a
 p(3)=(b(2)*p(2)+0.1_dp*h2*xx)/a(3)
 if(present(p3val)) then
   p(3)=p3val
   IF (Grid%type==loggrid) p(3)=p3val/Grid%pref(3)
 endif
 nodes=0
 DO i=4,many
   p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2)-c(i-1))/a(i)
   IF (p(i)*p(i-1) < 0._dp) nodes=nodes+1
   !renormalize if necessary
   scale=ABS(p(i))
   IF (scale > vlarg) THEN
     scale=1._dp/scale
     p(1:i)=scale*p(1:i)
   ENDIF
 ENDDO
 wfn(1:many)=p(1:many)
 IF (Grid%type==loggrid) THEN
   wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
 ENDIF
 LIBPAW_DEALLOCATE(a)
 LIBPAW_DEALLOCATE(b)
 LIBPAW_DEALLOCATE(p)
 LIBPAW_DEALLOCATE(c)
END SUBROUTINE forward_numerov


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine initgridwithn(Grid,type,n,r0,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initgridwithn(Grid,type,n,r0,h)
  TYPE (GridInfo), INTENT(INOUT) :: Grid
  INTEGER, INTENT(IN) :: type,n       !type=1 lin,2 log
  REAL(dp), INTENT(IN) :: r0,h
  INTEGER :: i
  IF (type==loggrid) THEN
    Grid%h=h
    Grid%r0=r0
    Grid%type=loggrid
    Grid%ishift=5
    Grid%n=n
    LIBPAW_ALLOCATE(Grid%r,(n))
    LIBPAW_ALLOCATE(Grid%drdu,(n))
    LIBPAW_ALLOCATE(Grid%pref,(n))
    LIBPAW_ALLOCATE(Grid%rr02,(n))
    DO i=1,n
      Grid%r(i)=r0*(EXP(Grid%h*(i-1))-1)
      Grid%drdu(i)=r0*EXP(Grid%h*(i-1))
      Grid%pref(i)=r0*EXP(Grid%h*(i-1)/2.d0)
      Grid%rr02(i)=(Grid%r(i)+r0)**2
    ENDDO
    Grid%range=Grid%r(n)
    if(has_to_print)WRITE(STD_OUT,*) 'InitGridwithn: -- logarithmic ',n, h,Grid%range,r0
  ELSEIF (type==lineargrid) THEN
    Grid%h=h
    Grid%r0=0.d0
    Grid%type=lineargrid
    Grid%ishift=25
    Grid%n=n
    LIBPAW_ALLOCATE(Grid%r,(n))
    LIBPAW_ALLOCATE(Grid%drdu,(n))
    LIBPAW_ALLOCATE(Grid%pref,(n))
    LIBPAW_ALLOCATE(Grid%rr02,(n))
    DO i=1,n
      Grid%r(i)=(Grid%h*(i-1))
      Grid%drdu(i)=1.d0
      Grid%pref(i)=1.d0
      Grid%rr02(i)=1.d0
    ENDDO
    Grid%range=Grid%r(n)
    if(has_to_print) WRITE(STD_OUT,*) 'InitGridwithn: -- linear ',n, h,Grid%range,r0
  ELSE
    write(std_out,*) 'error in initgridwithn -- type =',type    
    stop
  ENDIF
END SUBROUTINE initgridwithn





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 10. pseudo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE sethat(Grid,PAW,gaussparam,besselopt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE sethat(Grid,PAW,gaussparam,besselopt)
  TYPE(GridInfo), INTENT(IN) :: Grid
  TYPE(PseudoInfo), INTENT(INOUT) :: PAW
  INTEGER,INTENT(IN), OPTIONAL :: besselopt
  REAL(dp),INTENT(IN), OPTIONAL :: gaussparam
  INTEGER :: n,irc,irc_shap,i
  REAL(dp), POINTER :: r(:)
  REAL(dp) :: h,con,rc,rc_shap,selfen,d,dd,jbes1,jbes2,qr
  REAL(dp) :: al(2),ql(2)
  n=Grid%n
  h=Grid%h
  irc=PAW%irc
  rc=PAW%rc
  irc_shap=PAW%irc_shap
  rc_shap=PAW%rc_shap
  r=>Grid%r
  PAW%hatden=0
  PAW%projshape=0
  PAW%hatshape=0
  PAW%projshape(1)=1
  PAW%hatshape(1)=1
  DO i=2,irc-1
    PAW%projshape(i)=(SIN(pi*r(i)/rc)/(pi*r(i)/rc))**2
  ENDDO
  if(present(gaussparam)) then
    d=rc_shap/SQRT(LOG(1.d0/gaussparam))
    PAW%gausslength=d
    DO i=2,irc
      PAW%hatshape(i)=EXP(-(r(i)/d)**2)
    ENDDO
    PAW%irc_shap=PAW%irc
    PAW%rc_shap=PAW%rc
  else if(present(besselopt)) then
    call shapebes(al,ql,0,rc_shap)
    DO i=1,irc_shap-1
      qr=ql(1)*r(i);CALL jbessel(jbes1,d,dd,0,0,qr)
      qr=ql(2)*r(i);CALL jbessel(jbes2,d,dd,0,0,qr)
      PAW%hatshape(i)=al(1)*jbes1+al(2)*jbes2
    ENDDO
  else
    DO i=2,irc_shap-1
      PAW%hatshape(i)=(SIN(pi*r(i)/rc_shap)/(pi*r(i)/rc_shap))**2
    ENDDO
  endif
  PAW%hatden(1:irc)=PAW%hatshape(1:irc)*(r(1:irc)**2)
  !  normalize
  if (.not.present(besselopt)) then
    con=integrator(Grid,PAW%hatden,1,PAW%irc_shap)
    if(has_to_print) WRITE(STD_OUT,*) ' check hatden normalization', con
    PAW%hatden=PAW%hatden/con
  endif
  CALL atompaw_poisson(Grid,con,PAW%hatden,PAW%hatpot,selfen)
  if(has_to_print) WRITE(STD_OUT,*) 'Self energy for L=0 hat density  ', selfen
  if(has_to_print) WRITE(STD_OUT,*) 'hatden charge  ', con
END SUBROUTINE sethat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE FindVlocfromVeff(Grid,Orbit,PAW)
!  Assumes prior call to SUBROUTINE calculate_tvtau
!  which now fills PAW%tden and PAW%ttau
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FindVlocfromVeff(Grid,PAW,atp,potshift)
  ! TODO : needvtau,allocates
  TYPE(GridInfo), INTENT(INOUT) :: Grid
  type(atompaw_type),intent(inout) :: atp
  TYPE(PseudoInfo), INTENT(INOUT) :: PAW
  real(dp),intent(in) :: potshift
  REAL(dp), POINTER  :: r(:)
  REAL(dp) :: h,tq,rat,q00,etxc,eexc
  INTEGER :: n,irc,nbase
  REAL(dp), allocatable :: d(:),v(:),vxB(:),vxK(:)
  REAL(dp), allocatable :: t(:),vt(:),vthat(:)
  CALL FillHat(Grid,PAW,atp%besselshapefunction)
  !if (Vlocalindex == SETVLOC) then
  !  write(std_out,*) 'Vloc == VlocCoef*shapefunc  '
  !  return
  !endif
  n=Grid%n;irc=PAW%irc
  nbase=PAW%nbase
  h=Grid%h ; r=>Grid%r
  irc=max(PAW%irc,PAW%irc_shap,PAW%irc_vloc,PAW%irc_core)
  ! Recalculate den and tau      
  PAW%valetau=0.d0;PAW%tvaletau=0.d0
  LIBPAW_ALLOCATE(d,(n))
  LIBPAW_ALLOCATE(vxB,(n))
  LIBPAW_ALLOCATE(v,(n))
  LIBPAW_ALLOCATE(vxK,(n))
  LIBPAW_ALLOCATE(t,(n))
  LIBPAW_ALLOCATE(vt,(n))
  LIBPAW_ALLOCATE(vthat,(n))      
  d=PAW%den-PAW%tden
  tq=integrator(Grid,d,1,irc)
  if(has_to_print) write(std_out,*) ' Delta Qval = ', tq
  !     Compute VH(tDEN+hatDEN)
  d=PAW%tden+tq*PAW%hatden
  call poisson_marc(Grid,q00,d,v,rat)
  if(has_to_print) write(std_out,*) ' Completed Poisson with q00 = ', q00
  if(has_to_print) write(std_out,*) ' Completed Poisson with v(n) = ', v(n)
  !  Compute Blochl exc      
  d=PAW%tden+PAW%tcore
  t=PAW%tcoretau+PAW%tvaletau
  CALL exch(Grid,d,vxB,etxc,eexc,itype=atp%itype,needvtau=atp%Pot%needvtau,&
&           tau=t,vtau=vt,xc_functionals=atp%xc_functionals)
  !     Compute Kresse   exc Vxc(tcore+tDEN+hatDEN)
  d=PAW%tcore+PAW%tden+tq*PAW%hatden
  t=PAW%tcoretau+PAW%tvaletau
  CALL exch(Grid,d,vxK,etxc,eexc,itype=atp%itype,needvtau=atp%Pot%needvtau,&
&           tau=t,vtau=vthat,xc_functionals=atp%xc_functionals)
  PAW%abinitvloc=zero; PAW%abinitnohat=zero
  PAW%abinitnohat(2:n)=(PAW%rveff(2:n)-v(2:n)-vxB(2:n))/r(2:n)+potshift
  call extrapolate(PAW%abinitnohat)
  PAW%abinitvloc(2:n)=(PAW%rveff(2:n)-v(2:n)-vxK(2:n))/r(2:n)+potshift  
  call extrapolate(PAW%abinitvloc)
  ! Reassess poscorenhat      
  ! check if PAW%tcore+PAW%tden+tq*PAW%hatden is positive      
  PAW%poscorenhat=.true.
  LIBPAW_DEALLOCATE(v)
  LIBPAW_DEALLOCATE(vxB)
  LIBPAW_DEALLOCATE(vxK)
  LIBPAW_DEALLOCATE(d)
  LIBPAW_DEALLOCATE(vt)
  LIBPAW_DEALLOCATE(vthat)
END SUBROUTINE FindVlocfromVeff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE SetPAWOptions2(Grid,Orbit,Pot,success,atp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetPAWOptions2(atp,success)
 LOGICAL , INTENT(OUT) :: success
 type(atompaw_type), intent(inout) :: atp
 INTEGER :: l
 integer :: Projectorindex,Vlocalindex
 real(dp) :: e
 success=.true.
 IF (atp%projector_type==PROJECTOR_TYPE_MARSMAN) Projectorindex=MARSMAN
 IF (atp%vloc_type==VLOC_TYPE_MTROULLIER)  Vlocalindex=MTROULLIER
 IF (atp%vloc_type==VLOC_TYPE_ULTRASOFT)   Vlocalindex=ULTRASOFT
 IF (atp%vloc_type==VLOC_TYPE_BESSEL)      Vlocalindex=BESSEL
 IF (atp%vloc_type==VLOC_TYPE_VPSMATCHNC)  Vlocalindex=VPSMATCHNC
 IF (atp%vloc_type==VLOC_TYPE_VPSMATCHNNC) Vlocalindex=VPSMATCHNNC
 !Store the description of Vloc scheme in a string
 if (Vlocalindex==MTROULLIER) then
   l=atp%vloc_l ; e=atp%vloc_ene
   if(has_to_print) then
     WRITE(atp%PAW%Vloc_description,&
&    '("Vloc: Norm-conserving Troullier-Martins with l=",i1,";e=",1p,1e12.4)')l,e
   endif
 endif
 if (Vlocalindex==ULTRASOFT) then
   l=atp%vloc_l ; e=atp%vloc_ene
   if(has_to_print) then
     WRITE(atp%PAW%Vloc_description,&
&    '("Vloc: Non norm-conserving form with l= ",i1,";e= ",1p,1e12.4)')l,e
   endif
 endif
 if (Vlocalindex==BESSEL) then
   if(has_to_print) then
     WRITE(atp%PAW%Vloc_description,&
&      '("Vloc: truncated form - Vps(r)=A.sin(qr)/r for r<rc")')
   endif
 endif
 !Shape function parameters (from input dataset)
 if (atp%shapefunc_type==SHAPEFUNC_TYPE_GAUSSIAN) then
   atp%gaussianshapefunction=.true.
   CALL sethat(atp%Grid,atp%PAW,gaussparam=atp%shapefunc_gaussian_param)    !Gaussian shape function
 else if (atp%shapefunc_type==SHAPEFUNC_TYPE_BESSEL) then
   atp%besselshapefunction=.true.
   CALL sethat(atp%Grid,atp%PAW,besselopt=0)               ! Bessel shape function
 else
   CALL sethat(atp%Grid,atp%PAW)                          ! sinc^2 shape function
 endif
 !Call the routine computing Vloc - Not mGGA
 IF (.NOT.atp%needvtau) THEN
   IF (Vlocalindex==MTROULLIER.and.Projectorindex/=HFPROJ) THEN
     CALL troullier(atp%Grid,atp%Pot,atp%PAW,l,e,atp%needvtau,atp%scalarrelativistic)
   ENDIF
   IF (Vlocalindex==ULTRASOFT) CALL nonncps(atp%Grid,atp%Pot,atp%PAW,l,e,atp%scalarrelativistic)
   IF (Vlocalindex==BESSEL) CALL besselps(atp%Grid,atp%Pot,atp%PAW)
   call makebasis_marsman(atp)
 ENDIF
 !Call the routine computing Vloc - mGGA case
! IF (atp%needvtau) THEN
!    !All compatibility checks in input_dataset_read routine
!    if(has_to_print) WRITE(STD_OUT,*) 'Sequence of dataset construction modified for MGGA'      
!    if(has_to_print) WRITE(STD_OUT,*) ' Not all possibilites tested carefully yet.... '
!    if(has_to_print) WRITE(STD_OUT,*) ' Some possibilites not yet programmed.... '
!    !Calculate PAW%vtau and PAW%tvtau     
!    CALL calculate_tvtau(atp%Grid,atp%PAW,atp%itype)
!    !Set pseudoptentials     
! !   IF (Vlocalindex==MTROULLIER.and.(TRIM(atp%Orbit%exctype)/='HF')) then
! !     if(has_to_print) WRITE(STD_OUT,*) 'TROULLIER PS not available for MGGA '
! !     if(has_to_print) WRITE(STD_OUT,*) ' calling VPSmatch with norm conservation instead '     
! !     CALL VPSmatch(atp%Grid,atp%Pot,atp%PAW,l,e,.true.,atp%scalarrelativistic)
! !   ENDIF         
!    IF (Vlocalindex==VPSMATCHNNC) CALL VPSmatch(atp%Grid,atp%Pot,atp%PAW,l,e,.false.,atp%scalarrelativistic)
!    IF (Vlocalindex==VPSMATCHNC) CALL VPSmatch(atp%Grid,atp%Pot,atp%PAW,l,e,.true.,atp%scalarrelativistic)
!    IF (Vlocalindex==ULTRASOFT) CALL nonncps(atp%Grid,atp%Pot,atp%PAW,l,e,atp%scalarrelativistic)
!    IF (Vlocalindex==BESSEL) CALL besselps(atp%Grid,atp%Pot,atp%PAW)
!    !Calculate projectors
!    call makebasis_marsman(atp)
! ENDIF
  !Output in summary file
 IF (atp%needvtau) THEN
   if(has_to_print) WRITE(std_out,*) 'Sequence of dataset construction steps modified for mGGA'      
   if(has_to_print) WRITE(std_out,*) 'Only projectors from Vanderbilt scheme available'
 ENDIF
 CALL StoreTOCCWFN(atp%PAW)
END SUBROUTINE SetPAWOptions2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    SUBROUTINE Troullier(Grid,Pot,PAW,l,e,needvtau,scalarrelativistic)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Troullier(Grid,Pot,PAW,l,e,needvtau,scalarrelativistic)
 TYPE(Gridinfo), INTENT(IN) :: Grid
 TYPE(Potentialinfo), INTENT(INout) :: Pot
 TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
 INTEGER,INTENT(IN) :: l
 logical,intent(in) :: needvtau,scalarrelativistic
 real(dp),INTENT(IN) :: e
 real(dp), ALLOCATABLE :: VNC(:)
 real(dp) :: A0,A,B,B0,C,C0,D,F,S
 real(dp) :: Coef(6),Coef0,Coef0old
 real(dp) :: h,rc,delta,x,dpp,ddpp,dddpp,ddddpp
 INTEGER :: i,n,iter,nr,nodes,irc
 INTEGER, PARAMETER :: niter=5000
 real(dp), PARAMETER :: small=tol9
 real(dp), ALLOCATABLE ::  wfn(:),p(:),dum(:),aux(:)
 real(dp), POINTER :: r(:),rv(:)
 n=Grid%n
 h=Grid%h
 r=>Grid%r
 rv=>Pot%rv
 nr=min(PAW%irc_vloc+5,n)
 irc=PAW%irc_vloc
 rc=PAW%rc_vloc
 LIBPAW_ALLOCATE(VNC,(n))
 LIBPAW_ALLOCATE(wfn,(nr))
 LIBPAW_ALLOCATE(p,(nr))
 LIBPAW_ALLOCATE(dum,(nr))
 LIBPAW_ALLOCATE(aux,(nr))
 if (scalarrelativistic) then
   CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
 else if (needvtau) then
   CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
 else
   CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
 endif
 IF (wfn(irc)<0) wfn=-wfn
 dum(1:irc)=(wfn(1:irc)**2)
 S=integrator(Grid,dum(1:irc),1,irc)
 A0=LOG(wfn(irc)/(rc**(l+1)))
 B0=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
 C0=rc*(rv(irc)-rc*e)-B0*(B0+2*l+2)
 D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))-2*B0*C0-2*(l+1)*(C0-B0)
 F=rc*(2*rv(irc)-rc*(2*Gfirstderiv(Grid,irc,rv) &
&     -rc*Gsecondderiv(Grid,irc,rv)))+&
&     4*(l+1)*(C0-B0)-2*(l+1)*D-2*C0**2-2*B0*D
 if(has_to_print) WRITE(STD_OUT,*) 'In troullier -- matching parameters',S,A0,B0,C0,D,F
 delta=1.d10
 iter=0
 Coef0=0
 DO WHILE(delta>small.AND.iter<=niter)
   iter=iter+1
   A=A0-Coef0
   B=B0
   C=C0
   CALL EvaluateTp(l,A,B,C,D,F,coef)
   dum=0
   DO  i=1,irc
     x=(r(i)/rc)**2
     p(i)=x*(Coef(1)+x*(Coef(2)+x*(Coef(3)+&
&         x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
     dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
   ENDDO
   Coef0old=Coef0
   x=integrator(Grid,dum(1:irc),1,irc)
   Coef0=(LOG(S/x))/2
   delta=ABS(Coef0-Coef0old)
 ENDDO
 if(has_to_print) WRITE(STD_OUT,*) '  VNC converged in ', iter,'  iterations'
 if(has_to_print) WRITE(STD_OUT,*) '  Coefficients  -- ', Coef0,Coef(1:6)
 ! Now  calculate VNC
 if (needvtau) then
   aux=0._dp
   call derivative(Grid,PAW%tvtau,aux,1,nr)
 endif     
 VNC=0
 DO  i=2,nr
   x=(r(i)/rc)**2
   p(i)=Coef0+x*(Coef(1)+x*(Coef(2)+&
&       x*(Coef(3)+x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
   dpp=2*r(i)/(rc**2)*(Coef(1)+x*(2*Coef(2)+x*(3*Coef(3)+&
&       x*(4*Coef(4)+x*(5*Coef(5)+x*6*Coef(6))))))
   ddpp=(1/(rc**2))*(2*Coef(1)+x*(12*Coef(2)+x*(30*Coef(3)+&
&       x*(56*Coef(4)+x*(90*Coef(5)+x*132*Coef(6))))))
   dddpp=(r(i)/rc**4)*(24*Coef(2)+x*(120*Coef(3)+x*(336*Coef(4)+&
&       x*(720*Coef(5)+x*1320*Coef(6)))))
   ddddpp=(1/(rc**4)*(24*Coef(2)+x*(360*Coef(3)+x*(1680*Coef(4)+&
&       x*(5040*Coef(5)+x*11880*Coef(6))))))
   IF (i==irc) THEN
     if(has_to_print) WRITE(STD_OUT,*) 'check  dp ', dpp,  B0/rc
     if(has_to_print) WRITE(STD_OUT,*) 'check ddp ', ddpp, C0/rc**2
     if(has_to_print) WRITE(STD_OUT,*) 'check dddp', dddpp, D/rc**3
     if(has_to_print) WRITE(STD_OUT,*) 'check ddddp', ddddpp, F/rc**4
   ENDIF
   if (needvtau) then
     VNC(i)=e+(1._dp+PAW%tvtau(i))*(ddpp+dpp*(dpp+2*(l+1)/r(i))) &
&            +aux(i)*(dpp+l/r(i))                       
   else        
     VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
   endif
     dum(i)=(r(i)**(l+1))*EXP(p(i))
 ENDDO
 x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
 if(has_to_print) WRITE(STD_OUT,*) 'check norm ',x,S
 VNC(irc:n)=rv(irc:n)/r(irc:n)
 PAW%rveff(1:n)=VNC(1:n)*r(1:n)
 LIBPAW_DEALLOCATE(VNC)
 LIBPAW_DEALLOCATE(wfn)
 LIBPAW_DEALLOCATE(p)
 LIBPAW_DEALLOCATE(dum)
 LIBPAW_DEALLOCATE(aux)
END SUBROUTINE troullier


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! StoreTOCCWFN(PAW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE StoreTOCCWFN(PAW)
 TYPE(PseudoInfo), INTENT(INOUT) :: PAW
 INTEGER :: io,ib
 do io=1,PAW%TOCCWFN%norbit
   if (PAW%valencemap(io)>0) then
     ib=PAW%valencemap(io)
     PAW%TOCCWFN%wfn(:,io)=PAW%tphi(:,ib)
   endif
 enddo
END SUBROUTINE StoreTOCCWFN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE setbasis(Grid,Pot,Orbit,PAW,atp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE setbasis(Grid,Pot,Orbit,PAW,atp,potshift)
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), intent(in) :: potshift
 TYPE(PotentialInfo), INTENT(INout) :: Pot
 TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
 TYPE(PseudoInfo), intent(inout) :: PAW
 TYPE(atompaw_type), intent(in) :: atp
 INTEGER :: n,irc,nbase,l,lmax,mxbase,currentnode
 INTEGER :: i,k,io,nbl,nr,nodes,ibasis_add
 real(dp) :: h,rc,energy,rat,range
 real(dp),allocatable :: checkden(:),valden(:)
 real(dp), POINTER  :: r(:)
 n=Grid%n
 h=Grid%h
 r=>Grid%r
 irc=PAW%irc
 nr=MIN(irc+100,n-100)
 rc=PAW%rc
 lmax=PAW%lmax
 LIBPAW_ALLOCATE(checkden,(n))
 LIBPAW_ALLOCATE(valden,(n))
 checkden=0._dp;valden=0._dp
!Check beginning valence density
 DO io=1,Orbit%norbit 
   if (.not.Orbit%iscore(io)) then
     valden=valden+Orbit%occ(io)*(Orbit%wfn(:,io)**2)        
   endif
 ENDDO   
 ! set AErefrv
 PAW%AErefrv=Pot%rv
 PAW%rvx=Pot%rvx
 PAW%exctype=Orbit%exctype
 nbase=PAW%nbase
 mxbase=PAW%mxbase
 ! "filter" occupied states for long-range noise
 DO io=1,Orbit%norbit
   Call Filter(n,Orbit%wfn(:,io),machine_zero)
 ENDDO
 call CopyOrbit(Orbit,PAW%OCCwfn)
 call CopyOrbit(Orbit,PAW%TOCCwfn)
 PAW%valencemap=-13
 IF (Orbit%exctype=='HF') THEN
   range=r(n)
   rat=-1.d30; k=1
   do io=1,Orbit%norbit
     if ((Orbit%occ(io)>tol5).and.Orbit%eig(io)>rat) then
       rat=Orbit%eig(io)
       k=io
     endif
   enddo
   do i=n,irc+1,-1
     if (ABS(Orbit%wfn(i,k))>tol4) then
       range=r(i)
       exit
     endif
   enddo
   if(has_to_print) write(std_out,*) 'range ', k,range!; call flush_unit(std_out)
 ENDIF
 if(has_to_print) WRITE(STD_OUT,*) '  basis functions:'
 if(has_to_print) WRITE(STD_OUT,*)' No.   n     l         energy         occ   '
 nbase=0 ; ibasis_add=1
 DO l=0,lmax
   currentnode=-1
   nbl=0
   DO io=1,Orbit%norbit    ! cycle through all configuration
     IF (Orbit%l(io).EQ.l) THEN
       currentnode=Orbit%np(io)-l-1     
       IF (.NOT.Orbit%iscore(io)) THEN
         nbl=nbl+1
         nbase=nbase+1
         PAW%np(nbase)=Orbit%np(io)
         PAW%l(nbase)=l
         PAW%nodes(nbase)=PAW%np(nbase)-l-1
         if(has_to_print) write(std_out,*) 'l,nbase,node',l,nbase,currentnode
         PAW%eig(nbase)=Orbit%eig(io)
         if(Orbit%issemicore(io)) PAW%eig(nbase)=PAW%eig(nbase)-potshift
         PAW%occ(nbase)=Orbit%occ(io)
         if(Orbit%frozenvalecalculation.and.(.not.Orbit%issemicore(io))) then
           energy=PAW%eig(nbase)
           Orbit%wfn(:,io)=zero
           if (Orbit%scalarrelativistic) then
             CALL unboundsr(Grid,Pot,n,l,energy,Orbit%wfn(:,io),nodes)
           else
             CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,&
&                   nr,l,energy,Orbit%wfn(:,io),nodes)
           endif
         endif
         PAW%phi(:,nbase)=Orbit%wfn(:,io)
         if(Orbit%diracrelativistic) then 
           STOP 'Error -- setbasis subroutine not ready for diracrelativistic!'
         endif  
         PAW%valencemap(io)=nbase
         if(has_to_print) WRITE(STD_OUT,'(3i6,1p,2e15.6)') nbase,PAW%np(nbase),l,&
&                 PAW%eig(nbase),PAW%occ(nbase)!; call flush_unit(std_out)
       ENDIF
     ENDIF
   ENDDO
   generalizedloop: DO
     IF (ibasis_add>atp%nbasis_add) EXIT generalizedloop
     IF (atp%basis_add_l(ibasis_add)/=l) EXIT generalizedloop
     energy=atp%basis_add_energy(ibasis_add)
     IF (energy<0._dp.and.has_to_print) then
       WRITE(STD_OUT,*) 'energy is negative',energy,' -- WARNING WARNING !!!'
     endif
     nbase=nbase+1
     IF (nbase > mxbase ) THEN
       LIBPAW_ERROR('Error in  setbasis -- too many functions ')
     ENDIF
     PAW%l(nbase)=l
     PAW%np(nbase)=999
     PAW%nodes(nbase)=currentnode+1
     currentnode=PAW%nodes(nbase)
     if(has_to_print) write(std_out,*) 'l,nbase,node',l,nbase,currentnode
     PAW%eig(nbase)=energy
     PAW%occ(nbase)=0._dp
     PAW%phi(1:n,nbase)=0._dp
     if (Orbit%scalarrelativistic) then
       CALL unboundsr(Grid,Pot,n,l,energy,PAW%phi(:,nbase),nodes)
     else if (Pot%needvtau) then
       CALL unboundked(Grid,Pot,n,l,energy,PAW%phi(:,nbase),nodes)
     else
       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,&
&                 nr,l,energy,PAW%phi(:,nbase),nodes)
     endif
     rat=MAX(ABS(PAW%phi(irc,nbase)),ABS(PAW%phi(irc+1,nbase)))
     rat=DSIGN(rat,PAW%phi(irc,nbase))
     PAW%phi(1:n,nbase)=PAW%phi(1:n,nbase)/rat
     if(has_to_print) write(std_out,*) 'MAX PHI=',nbase,maxval(PAW%phi(:,nbase))
     if(has_to_print) then
       WRITE(STD_OUT,'(3i6,1p,2e15.6)') nbase,PAW%np(nbase),l,             &
&           PAW%eig(nbase),PAW%occ(nbase)
     endif
     nbl=nbl+1
     ibasis_add=ibasis_add+1
   ENDDO generalizedloop
 ENDDO   ! end lmax loop
 if(has_to_print) WRITE(std_out,*) 'completed phi basis with ',nbase,' functions '
 PAW%nbase=nbase     ! reset nbase
 LIBPAW_DEALLOCATE(checkden)
 LIBPAW_DEALLOCATE(valden)
END SUBROUTINE setbasis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE smoothcore(Grid,orig,PAW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE smoothcore(Grid,orig,PAW)
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: orig(:)
 type(pseudoInfo), intent(inout) :: PAW
 real(dp) :: rc,h,x,y,z,u0,u2,u4
 real(dp), allocatable :: d1(:),d2(:)
 INTEGER :: i,n,irc
 n=Grid%n
 h=Grid%h
 irc=PAW%irc_core
 rc=PAW%rc_core
 if(has_to_print) write(std_out,*) 'In smoothcore ', irc,rc
 LIBPAW_ALLOCATE(d1,(n))
 LIBPAW_ALLOCATE(d2,(n))
 CALL derivative(Grid,orig,d1)
 CALL derivative(Grid,d1,d2)
 x=orig(irc)
 y=d1(irc)*rc
 z=d2(irc)*(rc*rc)
 if(has_to_print) write(std_out,*) 'smoothcore: x,y,z = ', x,y,z
 u0=3*x - 9*y/8 + z/8
 u2=-3*x + 7*y/4 - z/4
 u4=x - 5*y/8 + z/8
 if(has_to_print) write(std_out,*) 'smoothcore: u0,u2,u4 = ', u0,u2,u4
 PAW%tcore=orig
 do i=1,irc
   x=(Grid%r(i)/rc)**2
   PAW%tcore(i)= x*(u0+x*(u2+x*u4))
 enddo
 LIBPAW_DEALLOCATE(d1)
 LIBPAW_DEALLOCATE(d2)
END SUBROUTINE smoothcore


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!! SUBROUTINE smoothpower
!!    program to take array orig (all electron tau or den)
!!        and return smooth polynomial function for 0 \le r \le rc_core
!!        matching 4 points less than and equal to rc_core
!!          power could be 2 or 4, representing the leading power
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
SUBROUTINE smoothpower(Grid,power,orig,smooth,PAW)
  TYPE(GridInfo), INTENT(IN) :: Grid
  type(pseudoInfo), intent(in) :: PAW
  INTEGER, INTENT(IN) :: power
  REAL(dp), INTENT(IN) :: orig(:)
  REAL(dp), INTENT(INOUT) :: smooth(:)
  INTEGER, parameter :: terms=5
  REAL(dp) :: rc,h,x
  REAL(dp) :: aa(terms,terms),Ci(terms)
  INTEGER :: i,j,n,irc
  n=Grid%n
  h=Grid%h
  irc=PAW%irc_core
  rc=Grid%r(irc)
  if(has_to_print) write(std_out,*) 'smoothpower -- ', power,irc,rc
  aa=zero; Ci=zero
  do i=1,terms
    x=Grid%r(irc-terms+i)
    Ci(i)=orig(irc-terms+i)/x**power
    do j=1,terms
      aa(i,j)=x**(2*(j-1))
    enddo
  enddo
  call SolveAXeqBM(terms,aa,Ci,terms-1)
  if(has_to_print) write(std_out,*) 'Completed SolveAXeqB with coefficients'
  if(has_to_print) write(std_out,'(1p,10e15.7)') (Ci(i),i=1,terms)
  smooth=orig
  do i=1,irc-1
    smooth(i)=0
    x=Grid%r(i)
    do j=1,terms
      smooth(i)=smooth(i)+Ci(j)*(x**(power+2*(j-1)))
    enddo
  enddo
END SUBROUTINE smoothpower


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      SUBROUTINE setcoretail(Grid,coreden,PAW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE setcoretail(Grid,coreden,PAW,needvtau)
 type(logical),intent(in) :: needvtau
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: coreden(:)
 type(pseudoInfo), intent(inout) :: PAW
 real(dp) :: rc,h,z
 INTEGER :: i,n,irc
 n=Grid%n
 h=Grid%h
 irc=PAW%irc_core
 rc=PAW%rc_core
 If(.not.needvtau) then
   CALL smoothcore(Grid,coreden,PAW)
 else
   CALL smoothpower(Grid,2,coreden,PAW%tcore,PAW)
 endif
 PAW%core=coreden
 ! Find coretailpoints
 z = integrator(Grid,coreden)
 PAW%coretailpoints=PAW%irc+Grid%ishift    !! coretailpoints should be>=PAW%irc
 do i=PAW%irc+Grid%ishift,n
   if(ABS(z-integrator(Grid,coreden,1,i))<coretailtol) then
     PAW%coretailpoints=i
     exit
   endif
 enddo
 if(has_to_print) write(std_out,*) 'coretailpoints = ',PAW%coretailpoints
END SUBROUTINE setcoretail


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE nonncps(lmax,Grid,Pot)
!!  Creates  screened pseudopotential by inverting Schroedinger
!!    equation from a pseudized radial wave function of the form:
!!        Psi(r) = r**(l+1) * exp (a + b*r**2 + c*r**4 + d*r**6)
!!  No norm-conserving condition is imposed on Psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE nonncps(Grid,Pot,PAW,l,e,scalarrelativistic)
  logical, intent(in) :: scalarrelativistic
  TYPE(Gridinfo), INTENT(IN) :: Grid
  TYPE(Potentialinfo), INTENT(INout) :: Pot
  TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
  INTEGER,INTENT(IN) :: l
  REAL(8),INTENT(IN) :: e
  INTEGER :: i,irc,n,nr,nodes,i1,i2,i3,i4
  REAL(8) :: rc,x,y1,y2,y3,p0,p1,p2,p3,sgn
  REAL(8) :: b(4),c(4),d(4),amat(4,4)
  REAL(8),ALLOCATABLE ::  VNC(:),wfn(:),aux(:)
  REAL(8),POINTER :: r(:),rv(:)
  !Polynomial definitions
  p0(x,y1,y2,y3)=(x-y1)*(x-y2)*(x-y3)
  p1(x,y1,y2,y3)=(x-y2)*(x-y3)+(x-y1)*(x-y3)+(x-y1)*(x-y2)
  p2(x,y1,y2,y3)=two*((x-y1)+(x-y2)+(x-y3))
  p3(x,y1,y2,y3)=six
  n=Grid%n
  r=>Grid%r
  rv=>Pot%rv
  nr=min(PAW%irc_vloc+10,n)
  irc=PAW%irc_vloc
  rc=PAW%rc_vloc
  LIBPAW_ALLOCATE(VNC,(n))
  LIBPAW_ALLOCATE(wfn,(nr))
  LIBPAW_ALLOCATE(aux,(nr))
  if (scalarrelativistic) then
    CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
  else if (Pot%needvtau) then
    CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
  else
    CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
  endif
  IF (wfn(irc)<0) wfn=-wfn
  DO i=2,nr
    wfn(i)=wfn(i)/r(i)**dble(l+1)
  ENDDO
  i1=irc-1;i2=i1+1;i3=i2+1;i4=i3+1
  c(1)=wfn(i1)/p0(r(i1),r(i2),r(i3),r(i4))
  c(2)=wfn(i2)/p0(r(i2),r(i3),r(i4),r(i1))
  c(3)=wfn(i3)/p0(r(i3),r(i4),r(i1),r(i2))
  c(4)=wfn(i4)/p0(r(i4),r(i1),r(i2),r(i3))
  d(1)=c(1)*p0(rc,r(i2),r(i3),r(i4)) + c(2)*p0(rc,r(i3),r(i4),r(i1)) + &
&      c(3)*p0(rc,r(i4),r(i1),r(i2)) + c(4)*p0(rc,r(i1),r(i2),r(i3))
  d(2)=c(1)*p1(rc,r(i2),r(i3),r(i4)) + c(2)*p1(rc,r(i3),r(i4),r(i1)) + &
&      c(3)*p1(rc,r(i4),r(i1),r(i2)) + c(4)*p1(rc,r(i1),r(i2),r(i3))
  d(3)=c(1)*p2(rc,r(i2),r(i3),r(i4)) + c(2)*p2(rc,r(i3),r(i4),r(i1)) + &
&      c(3)*p2(rc,r(i4),r(i1),r(i2)) + c(4)*p2(rc,r(i1),r(i2),r(i3))
  d(4)=c(1)*p3(rc,r(i2),r(i3),r(i4)) + c(2)*p3(rc,r(i3),r(i4),r(i1)) + &
&      c(3)*p3(rc,r(i4),r(i1),r(i2)) + c(4)*p3(rc,r(i1),r(i2),r(i3))
  sgn=d(1)/abs(d(1));d(1:4)=d(1:4)*sgn
  b(1)=log(d(1));b(2:4)=d(2:4)
  amat(1,1)= 1.0d0
  amat(2:4,1)= 0.0d0
  amat(1,2)= rc**2
  amat(2,2)= 2.0d0*d(1)*rc
  amat(3,2)= 2.0d0*d(1)   +2.0d0*d(2)*rc
  amat(4,2)=               4.0d0*d(2)   +2.0d0*d(3)*rc
  amat(1,3)= rc**4
  amat(2,3)=  4.0d0*d(1)*rc**3
  amat(3,3)= 12.0d0*d(1)*rc**2+ 4.0d0*d(2)*rc**3
  amat(4,3)= 24.0d0*d(1)*rc   +24.0d0*d(2)*rc**2+4.0d0*d(3)*rc**3
  amat(1,4)= rc**6
  amat(2,4)=   6.0d0*d(1)*rc**5
  amat(3,4)=  30.0d0*d(1)*rc**4+ 6.0d0*d(2)*rc**5
  amat(4,4)= 120.0d0*d(1)*rc**3+60.0d0*d(2)*rc**4+6.0d0*d(3)*rc**5
  CALL linsol(amat,b,4,4,4,4)
  if (Pot%needvtau) then
    aux=zero
    call derivative(Grid,PAW%tvtau,aux,1,nr)
  endif
  PAW%rveff(1)=0.d0
  DO i=2,irc-1
    c(1)=2.0d0*b(2)*r(i)+ 4.0d0*b(3)*r(i)**3+ 6.0d0*b(4)*r(i)**5
    c(2)=2.0d0*b(2)     +12.0d0*b(3)*r(i)**2+30.0d0*b(4)*r(i)**4
    if (pot%needvtau) then
      PAW%rveff(i)=r(i)*(e+(dble(2*l+2)*c(1)/r(i)+c(1)**2+&
&         c(2))*(1.d0+PAW%tvtau(i))+aux(i)*(c(1)+dble(l)/r(i)))
    else
      PAW%rveff(i)=r(i)*(e+dble(2*l+2)*c(1)/r(i)+c(1)**2+c(2))
    endif
  ENDDO
  PAW%rveff(irc:n)=rv(irc:n)
  LIBPAW_DEALLOCATE(VNC)
  LIBPAW_DEALLOCATE(wfn)
  LIBPAW_DEALLOCATE(aux)
END SUBROUTINE nonncps


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE VPSmatch(lmax,Grid,Pot,NC)
!  Creates  screened norm-conserving pseudopotential similar to the
!    approach of N. Troullier and J. L. Martins, PRB 43, 1993 (1991)
!    Uses p(r)=a0+f(r); f(r)=SUMm(Coef(m)*r^(2*m), where
!          m=1,2..6
!    Psi(r) = r^(l+1)*exp(p(r))
!    Modified for MGGA case and norm conserving condition is optional
!    norm conservation controlled with optional variable NC
!      defaults to no norm conservation
!  Note this program assumes that wfn keeps the same sign for
!    all matching points r(irc-match+1)....r(irc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE VPSmatch(Grid,Pot,PAW,l,e,NC,scalarrelativistic)
  logical,intent(in) :: scalarrelativistic
  TYPE(Gridinfo), INTENT(IN) :: Grid
  TYPE(Potentialinfo), INTENT(INout) :: Pot
  TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
  INTEGER,INTENT(IN) :: l
  REAL(dp),INTENT(IN) :: e
  LOGICAL,INTENT(IN), OPTIONAL :: NC
  REAL(dp), ALLOCATABLE :: VNC(:)
  REAL(dp) :: C0,S
  REAL(dp) :: Coef0,Coef0old
  REAL(dp) :: h,rc,delta,x
  INTEGER :: i,j,n,iter,nr,nodes,irc
  INTEGER, parameter :: match=6
  REAL(dp) :: AAA(match,match),BBB(match)
  REAL(dp) :: AAAA(match-1,match-1),BBBB(match-1)
  INTEGER, PARAMETER :: niter=5000
  REAL(dp), PARAMETER :: small=1.0d-9
  REAL(dp), ALLOCATABLE :: wfn(:),p(:),dum(:),aux(:),Kaux(:),v(:),dp1(:),ddp(:)
  REAL(dp), POINTER :: r(:),rv(:)
  LOGICAL :: normcons
  normcons=.false.
  if(PRESENT(NC)) normcons=NC
  if(has_to_print) write(std_out,*) 'Entering VPSmatch with normcons ', normcons
  n=Grid%n
  h=Grid%h
  r=>Grid%r
  rv=>Pot%rv
  nr=min(PAW%irc_vloc+10,n)
  irc=PAW%irc_vloc
  rc=PAW%rc_vloc
  LIBPAW_ALLOCATE(VNC,(n))
  LIBPAW_ALLOCATE(wfn,(n))
  LIBPAW_ALLOCATE(p,(n))
  LIBPAW_ALLOCATE(dum,(n))
  LIBPAW_ALLOCATE(aux,(n))
  LIBPAW_ALLOCATE(dp1,(n))
  LIBPAW_ALLOCATE(ddp,(n))
  LIBPAW_ALLOCATE(Kaux,(n))
  LIBPAW_ALLOCATE(v,(n))
  wfn=zero
  if (scalarrelativistic) then
    CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
  else if (Pot%needvtau) then
    CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
  else
    CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
  endif
  IF (wfn(irc)<0) wfn=-wfn
  !  first solve non-norm conserving results
  AAA=zero;BBB=zero
  Do i=1,match
    x=r(irc-match+i)
    BBB(i)=log(wfn(irc-match+i)/(x**(l+1)))
    do j=1,match
      AAA(i,j)=x**(2*(j-1))
    enddo
  Enddo
  CALL  SolveAXeqB(match,AAA,BBB,1.d20)
  if(has_to_print) write(std_out,*) 'Returned from SolveAXeqB in VPSmatch with Coefficients '
  if(has_to_print) write(std_out,'(1p,50e16.7)') (BBB(i),i=1,match)
  ! Now  calculate VNC
  if (Pot%needvtau) then
    aux=0.d0
    call derivative(Grid,PAW%tvtau,aux,1,nr)
    Kaux=0.d0
    call derivative(Grid,PAW%Ktvtau,Kaux,1,nr)   ! Kresse form
  endif
  VNC=zero;p=zero;dp1=zero;ddp=zero;dum=zero
  !specific for match=6
  DO  i=2,nr
    x=(r(i))**2
    p(i)=BBB(1)+x*(BBB(2)+x*(BBB(3)+x*(BBB(4)+x*(BBB(5)+x*BBB(6)))))
    dp1(i)=two*r(i)*(BBB(2)+x*(two*BBB(3)+x*(three*BBB(4)+x*(four*BBB(5)+five*x*BBB(6)))))
    ddp(i)=2.0_dp*(BBB(2)+x*(6.0_dp*BBB(3)+x*(15.0_dp*BBB(4)+x*(28.0_dp*BBB(5)+45.0_dp*x*BBB(6)))))
    if (Pot%needvtau) then
      VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddp(i)+ &
&       dp1(i)*(dp1(i)+two*(l+1)/r(i))) &
&       +aux(i)*(dp1(i)+l/r(i))
      v(i)=e+(1.d0+PAW%Ktvtau(i))*(ddp(i)+ &
&       dp1(i)*(dp1(i)+two*(l+1)/r(i))) &
&       +Kaux(i)*(dp1(i)+l/r(i))
    else
      VNC(i)=e+ddp(i)+dp1(i)*(dp1(i)+two*(l+1)/r(i))
      v(i)=e+ddp(i)+dp1(i)*(dp1(i)+two*(l+1)/r(i))
    endif
      dum(i)=(r(i)**(l+1))*EXP(p(i))
  ENDDO
  S=overlap(Grid,wfn(1:irc),wfn(1:irc),1,irc)
  C0=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
  if(has_to_print) WRITE(STD_OUT,*) 'check norm ',C0,S
  VNC(irc:n)=rv(irc:n)/r(irc:n)
  v(irc:n)=rv(irc:n)/r(irc:n)
  PAW%rveff(1:n)=VNC(1:n)*r(1:n)
  PAW%Krveff(1:n)=v(1:n)*r(1:n)
  if(has_to_print) write(std_out,*) 'Completed non-norm-conserving PS '
  if(.not.normcons) return
  !  Iterate to find norm conserving results
  C0=C0/(EXP(2*BBB(1)))
  delta=1.d10
  iter=zero
  Coef0=0.5d0*log(S/C0)
  DO WHILE(delta>small.AND.iter<=niter)
    Coef0old=Coef0
    iter=iter+1
    AAAA=zero;BBBB=zero
    Do i=1,match-1
      x=r(irc-match+1+i)
      BBBB(i)=log(wfn(irc-match+1+i)/(x**(l+1)))-Coef0old
      do j=1,match-1
        AAAA(i,j)=x**(2*(j))
      enddo
    Enddo
    CALL  SolveAXeqB(match-1,AAAA,BBBB,1.d20)
    if(has_to_print) write(std_out,*) 'Returned from SolveAXeqB in VPSmatch with Coefficients'
    if(has_to_print) write(std_out,'(1p,50e16.7)') (BBBB(i),i=1,match-1)
    !specific for match-1=5
    p=zero;dum=zero
    DO  i=2,nr
       x=(r(i))**2
       p(i)=x*(BBBB(1)+x*(BBBB(2)+x*(BBBB(3)+x*(BBBB(4)+x*BBBB(5)))))
       dum(i)=(r(i)**(l+1))*EXP(Coef0+p(i))
    ENDDO
    C0=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
    if(has_to_print) WRITE(std_out,*) 'check norm ',C0,S
    Coef0=0.5d0*log(S/C0)
    delta=ABS(Coef0-Coef0old)
    if(has_to_print) WRITE(std_out,'(" VNC: iter Coef0 delta",i5,1p,2e15.7)') iter,Coef0,delta
  ENDDO
  if(has_to_print) WRITE(std_out,*) '  VNC converged in ', iter,'  iterations'
  if(has_to_print) WRITE(std_out,*) '  Coefficients  -- ', Coef0,BBBB(1:match-1)
  ! Now  calculate VNC
  if (Pot%needvtau) then
    aux=0.d0
    call derivative(Grid,PAW%tvtau,aux,1,nr)
    Kaux=0.d0
    call derivative(Grid,PAW%Ktvtau,Kaux,1,nr)
  endif
  VNC=zero;p=zero;dp1=zero;ddp=zero;v=zero
  DO  i=2,nr
    x=(r(i))**2
    p(i)=x*(BBBB(1)+x*(BBBB(2)+x*(BBBB(3)+x*(BBBB(4)+x*BBBB(5)))))
    dp1(i)=two*r(i)*(BBBB(1)+x*(two*BBBB(2)+x*(three*BBBB(3)+x*(four*BBBB(4)+five*x*BBBB(5)))))
    ddp(i)=2.0_dp*(BBBB(1)+x*(6.0_dp*BBBB(2)+x*(15.0_dp*BBBB(3)+x*(28.0_dp*BBBB(4)+45.0_dp*x*BBBB(5)))))
    if (Pot%needvtau) then
      VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddp(i)+ &
&       dp1(i)*(dp1(i)+two*(l+1)/r(i))) &
&       +aux(i)*(dp1(i)+l/r(i))
      v(i)=e+(1.d0+PAW%Ktvtau(i))*(ddp(i)+ &
&       dp1(i)*(dp1(i)+two*(l+1)/r(i))) &
&       +Kaux(i)*(dp1(i)+l/r(i))
    else
      VNC(i)=e+ddp(i)+dp1(i)*(dp1(i)+two*(l+1)/r(i))
      v(i)=e+ddp(i)+dp1(i)*(dp1(i)+two*(l+1)/r(i))
    endif
    dum(i)=(r(i)**(l+1))*EXP(Coef0+p(i))
  ENDDO
  C0=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
  if(has_to_print) WRITE(std_out,*) 'check norm ',C0,S
  VNC(irc:n)=rv(irc:n)/r(irc:n)
  v(irc:n)=rv(irc:n)/r(irc:n)
  PAW%rveff(1:n)=VNC(1:n)*r(1:n)
  PAW%Krveff(1:n)=v(1:n)*r(1:n)
  if(has_to_print)write(std_out,*) 'Completed norm-conserving PS '
  if (iter.ge.niter) then
    write(std_out,*) 'Failed to converged norm-conserving VPS '
    stop
  endif
  LIBPAW_DEALLOCATE(VNC)
  LIBPAW_DEALLOCATE(wfn)
  LIBPAW_DEALLOCATE(p)
  LIBPAW_DEALLOCATE(dp1)
  LIBPAW_DEALLOCATE(ddp)
  LIBPAW_DEALLOCATE(dum)
  LIBPAW_DEALLOCATE(aux)
  LIBPAW_DEALLOCATE(Kaux)
  LIBPAW_DEALLOCATE(v)
  END SUBROUTINE VPSmatch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculate_tvtau    for MGGA case
!     Assume valence pseudo wavefunctions known
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_tvtau(Grid,PAW,itype)
  TYPE(GridInfo), INTENT(IN) :: Grid
  type(pseudoInfo),intent(inout) :: PAW
  integer,intent(in) :: itype
  INTEGER :: i
  REAL(dp), allocatable :: dp1(:),ddp(:),vxc(:),tvxc(:),locald(:),localtd(:)
  REAL(dp), allocatable :: Ktvxc(:),Kd(:)
  REAL(dp) :: exc,texc,sum,tsum
  REAL(dp), parameter :: small=1.d-5
  LIBPAW_ALLOCATE(dp1,(Grid%n))
  LIBPAW_ALLOCATE(ddp,(Grid%n))
  LIBPAW_ALLOCATE(vxc,(Grid%n))
  LIBPAW_ALLOCATE(tvxc,(Grid%n))
  LIBPAW_ALLOCATE(locald,(Grid%n))
  LIBPAW_ALLOCATE(localtd,(Grid%n))
  LIBPAW_ALLOCATE(Ktvxc,(Grid%n))
  LIBPAW_ALLOCATE(Kd,(Grid%n))
  locald=PAW%core
  localtd=PAW%tcore
  PAW%valetau=0.d0
  PAW%tvaletau=0.d0
  do i=1,PAW%nbase
    if(has_to_print) write(std_out,*) 'tvtau -- ', i,PAW%l(i),PAW%occ(i),PAW%eig(i)
    if (PAW%occ(i).gt.small) then
      locald=locald+PAW%occ(i)*(PAW%phi(:,i)**2)
      localtd=localtd+PAW%occ(i)*(PAW%tphi(:,i)**2)
      CALL taufromwfn(dp1,Grid,PAW%phi(:,i),PAW%l(i))
      CALL taufromwfn(ddp,Grid,PAW%tphi(:,i),PAW%l(i))
      PAW%valetau=PAW%valetau+PAW%occ(i)*dp1
      PAW%tvaletau=PAW%tvaletau+PAW%occ(i)*ddp
    endif
  enddo
  dp1=PAW%coretau+PAW%valetau
  ddp=PAW%tcoretau+PAW%tvaletau
  CALL exch(Grid,locald,vxc,sum,exc,itype,&
&    tau=dp1,vtau=PAW%vtau)
  CALL exch(Grid,localtd,tvxc,tsum,texc,itype,&
&    tau=ddp,vtau=PAW%tvtau)
  if(has_to_print) write(std_out,*) 'tvtau exc texc ', exc, texc
  ! Kresse form    
  Kd=locald-PAW%core-localtd+PAW%tcore
  sum=integrator(Grid,Kd)
  if(has_to_print) write(std_out,*) 'compensation charge in Ktvtau ', sum
  Kd=localtd+sum*PAW%hatden
  CALL exch(Grid,Kd,Ktvxc,tsum,texc,itype,&
&     tau=ddp,vtau=PAW%Ktvtau)
  LIBPAW_DEALLOCATE(dp1)
  LIBPAW_DEALLOCATE(ddp)
  LIBPAW_DEALLOCATE(vxc)
  LIBPAW_DEALLOCATE(tvxc)
  LIBPAW_DEALLOCATE(locald)
  LIBPAW_DEALLOCATE(localtd)
  LIBPAW_DEALLOCATE(Kd)
  LIBPAW_DEALLOCATE(Ktvxc)
END SUBROUTINE calculate_tvtau





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 11. Pseudo_sub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine hatL
!!   Calculates density associated with L component
!!    normalized to unity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE hatL(Grid,PAW,l,dhat,besselshapefunction)
  TYPE(GridInfo), INTENT(IN) :: Grid
  TYPE(PseudoInfo), INTENT(IN) :: PAW
  INTEGER, INTENT(IN) :: l
  logical, intent(in) :: besselshapefunction
  REAL(dp), INTENT(OUT) :: dhat(:)
  INTEGER :: n,irc,i
  REAL(dp), POINTER :: r(:)
  REAL(dp), ALLOCATABLE :: den(:),a(:)
  REAL(dp) :: h,con
  REAL(dp) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)
  n=Grid%n
  h=Grid%h
  r=>Grid%r
  irc=PAW%irc
  LIBPAW_ALLOCATE(den,(n))
  LIBPAW_ALLOCATE(a,(n))
  IF (besselshapefunction) THEN
    CALL shapebes(al,ql,l,PAW%rc_shap)
    DO i=1,PAW%irc_shap
      qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
      qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
      den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
    ENDDO
    IF (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
  ELSE
    DO i=1,n
      den(i)=(r(i)**l)*PAW%hatden(i)
    ENDDO
    a(1:n)=den(1:n)*(r(1:n)**l)
    con=integrator(Grid,a,1,PAW%irc_shap)
    den=den/con
  ENDIF
  dhat=zero
  dhat(1:n)=den(1:n)
  LIBPAW_DEALLOCATE(den)
  LIBPAW_DEALLOCATE(a)
END SUBROUTINE hatL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE FillHat(Grid,PAW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FillHat(Grid,PAW,besselshape)
  logical,intent(in) :: besselshape
  TYPE(GridInfo) , INTENT(IN):: Grid
  TYPE(PseudoInfo), INTENT(INOUT) :: PAW
  INTEGER :: ll,n,l
  ll=MAXVAL(PAW%TOCCWFN%l(:)); ll=MAX(ll,PAW%lmax); ll=2*ll
  n=Grid%n
  LIBPAW_ALLOCATE(PAW%g,(n,ll+1))
  DO l=0,ll
    CALL hatL(Grid,PAW,l,PAW%g(:,l+1),besselshape)
  ENDDO
END SUBROUTINE FillHat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE altdtij
! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
!   for r > rc, f1=t1, f2=t2
! on output: tij is difference kinetic energy matrix element in Rydberg units
!   tij =<f1|T|f2>-<t1|T|t2>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE altdtij(Grid,PAW,ib,ic,tij,needvtau)
 logical, intent(in) :: needvtau
 TYPE(GridInfo), INTENT(IN) :: Grid
 TYPE(PseudoInfo), INTENT(IN) :: PAW
 INTEGER, INTENT(IN) :: ib,ic
 REAL(8), INTENT(OUT) :: tij
 INTEGER :: n,i,l,irc
 REAL(8) :: angm
 REAL(8), POINTER :: r(:)
 REAL(8), ALLOCATABLE :: dum(:),tdel1(:),tdel2(:),aux(:),auxp(:)
 tij=0
 IF (PAW%l(ib)/=PAW%l(ic)) RETURN
 n=Grid%n;  r=>Grid%r;  l=PAW%l(ib);  irc=PAW%irc
 LIBPAW_ALLOCATE(dum,(n))
 LIBPAW_ALLOCATE(tdel1,(n))
 LIBPAW_ALLOCATE(tdel2,(n))
 LIBPAW_ALLOCATE(aux,(n))
 LIBPAW_ALLOCATE(auxp,(n))
 dum=zero
 DO i=2,irc
   dum(i)=PAW%ophi(i,ib)*PAW%Kop(i,ic)        !Corrected 6/6/2023 Thanks to MT
 ENDDO
 CALL derivative(Grid,PAW%otphi(:,ic),tdel1)
 CALL derivative(Grid,tdel1,tdel2)
 aux=1.0_dp;auxp=0.0_dp
 if(needvtau) then
   aux=1._dp+PAW%tvtau
   call derivative(Grid,PAW%tvtau,auxp)
 endif
 angm=l*(l+1)
 DO i=2,irc
    dum(i)=dum(i)+PAW%otphi(i,ib)*(aux(i)*(tdel2(i)-&
&        angm*PAW%otphi(i,ic)/(Grid%r(i)**2))+auxp(i)*&
&        (tdel1(i)-PAW%otphi(i,ic)/Grid%r(i)))
 ENDDO
 tij=integrator(Grid,dum,1,irc)
 LIBPAW_DEALLOCATE(dum)
 LIBPAW_DEALLOCATE(tdel1)
 LIBPAW_DEALLOCATE(tdel2)
 LIBPAW_DEALLOCATE(aux)
 LIBPAW_DEALLOCATE(auxp)
END SUBROUTINE altdtij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE besselps(Grid,Pot,PAW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE besselps(Grid,Pot,PAW)
 TYPE(Gridinfo), INTENT(IN) :: Grid
 TYPE(Potentialinfo), INTENT(IN) :: Pot
 TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
 INTEGER :: i,irc,n
 real(dp) :: rc,alpha,beta,vv,vvp,AA,QQ,xx(1)
 real(dp),POINTER :: r(:),rv(:)
 n=Grid%n
 r=>Grid%r
 rv=>Pot%rv
 irc=PAW%irc_vloc
 rc=PAW%rc_vloc
 vv=rv(irc);vvp=Gfirstderiv(Grid,irc,rv)
 alpha=1.D0-rc*vvp/vv;beta=1.D0
 call solvbes(xx,alpha,beta,0,1);QQ=xx(1)
 AA=vv/sin(QQ);QQ=QQ/rc
 PAW%rveff(1)=0._dp
 PAW%rveff(irc+1:n)=rv(irc+1:n)
 do i=2,irc
   PAW%rveff(i)=AA*sin(QQ*r(i))
 enddo
END SUBROUTINE besselps


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE EvaluateTp
!!   Inverts 5x5 matrix used  by troullier subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EvaluateTp(l,A,B,C,D,F,coef)
 INTEGER, INTENT(IN) :: l
 real(dp), INTENT(IN) :: A,B,C,D,F
 real(dp), INTENT(OUT) ::  coef(6)
 real(dp) :: t(6,6),coef10,old
 real(dp), PARAMETER :: small=1.e-10
 INTEGER :: i,n,iter
 INTEGER, PARAMETER :: niter=1000
 old=-1.e30; Coef10=-1; iter=-1
 DO WHILE (iter < niter .AND. ABS(old-coef10)> small)
   iter=iter+1
   t=0
   Coef(1)=A-Coef10; Coef(2)=B-2*Coef10;  Coef(3)=C-2*Coef10;
   Coef(4)=D;    Coef(5)=F
   Coef(6)=-Coef10**2
   DO i=1,6
     t(1,i)=1
     t(2,i)=2*i
     t(3,i)=2*i*(2*i-1)
     t(4,i)=2*i*(2*i-1)*(2*i-2)
     t(5,i)=2*i*(2*i-1)*(2*i-2)*(2*i-3)
   ENDDO
   t(6,1)=2*Coef10;  t(6,2)=2*l+5
   n=6
   CALL linsol(t,Coef,n,6,6,6)
   old=Coef10; Coef10=Coef10+Coef(1)
   Coef(1)=Coef10
 ENDDO
 IF (iter >= niter) THEN
   LIBPAW_ERROR('Error in EvaluateTP -- no convergence')
 ENDIF
END SUBROUTINE EvaluateTp




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 12. Pseudo_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    SUBROUTINE InitPAW(PAW,Grid,Orbit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE InitPAW(PAW,Grid,Orbit)
 TYPE(GridInfo), INTENT(IN) :: Grid
 TYPE(OrbitInfo), INTENT(IN) :: Orbit
 Type(PseudoInfo), INTENT(INOUT) :: PAW
 INTEGER :: io,l,n,mxbase,nbase
!Initialize logical variables
 PAW%multi_rc=.false.
 PAW%poscorenhat=.true.
 CALL DestroyPAW(PAW)
!Compute initial size of basis
 n=Grid%n
 nbase=0
 DO l=0,PAW%lmax
   DO io=1,Orbit%norbit    ! cycle through all configurations
     IF (Orbit%l(io).EQ.l.AND.(.NOT.Orbit%iscore(io))) THEN
       nbase=nbase+1
     ENDIF
   ENDDO
 ENDDO
 PAW%mxbase=nbase+6*max(1,PAW%lmax)
 mxbase=PAW%mxbase !Estimate excess
 PAW%nbase=nbase
 if(has_to_print) WRITE(STD_OUT,*) 'Found ', nbase,' valence basis functions '
 if(has_to_print) WRITE(STD_OUT,*) 'Allocating for ', mxbase, ' total basis functions'
 LIBPAW_POINTER_ALLOCATE(PAW%projshape,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%hatden,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%hatpot,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%hatshape,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%vloc,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%rveff,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%abinitvloc,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%abinitnohat,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%AErefrv,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%rvx,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%trvx,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%den,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%tden,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%core,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%tcore,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%coretau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%tcoretau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%valetau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%tvaletau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%vtau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%tvtau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%nhatv,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%Ktvtau,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%Krveff,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%Kunscreen,(n))
 PAW%projshape=0._dp;PAW%hatden=0._dp;PAW%hatpot=0._dp
 PAW%hatshape=0._dp;PAW%vloc=0._dp;PAW%rveff=0._dp
 PAW%abinitvloc=0._dp;PAW%abinitnohat=0._dp
 PAW%AErefrv=0._dp;PAW%rvx=0._dp;PAW%trvx=0._dp
 PAW%den=0._dp;PAW%tden=0._dp;PAW%core=0._dp;PAW%tcore=0._dp
 PAW%XCORECORE=0._dp;PAW%nhatv=0._dp
 PAW%coretau=0._dp;PAW%tcoretau=0._dp
 PAW%valetau=0._dp;PAW%tvaletau=0._dp
 PAW%vtau=0._dp;PAW%tvtau=0._dp
 PAW%Ktvtau=0._dp;PAW%Krveff=0._dp;PAW%Kunscreen=0._dp
 LIBPAW_POINTER_ALLOCATE(PAW%phi,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%tphi,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%tp,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%ophi,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%otphi,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%otp,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%np,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%l,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%eig,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%occ,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%ck,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%vrc,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%Kop,(n,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%rng,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%rcio,(mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%nodes,(mxbase))
 PAW%phi=0._dp;PAW%tphi=0._dp;PAW%tp=0._dp
 PAW%ophi=0._dp;PAW%otphi=0._dp;PAW%otp=0._dp
 PAW%eig=0._dp;PAW%occ=0._dp;PAW%vrc=0._dp;PAW%ck=0._dp;PAW%Kop=0._dp
 PAW%rcio=0._dp;PAW%np=0;PAW%l=0
 if(Orbit%diracrelativistic) then
   LIBPAW_POINTER_ALLOCATE(PAW%kappa,(mxbase))
   PAW%kappa=0
 endif        
 PAW%rng=Grid%n
 LIBPAW_POINTER_ALLOCATE(PAW%oij,(mxbase,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%dij,(mxbase,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%wij,(mxbase,mxbase))
 PAW%oij=0._dp;PAW%dij=0._dp;PAW%wij=0._dp
 LIBPAW_POINTER_ALLOCATE(PAW%rVf,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%rtVf,(n))
 LIBPAW_POINTER_ALLOCATE(PAW%Kij,(mxbase,mxbase))
 LIBPAW_POINTER_ALLOCATE(PAW%Vfij,(mxbase,mxbase))
 PAW%rVf=0._dp;PAW%rtVf=0._dp;PAW%Kij=0._dp;PAW%Vfij=0._dp
 IF (Orbit%exctype=='HF') THEN
   LIBPAW_POINTER_ALLOCATE(PAW%lmbd,(Orbit%norbit,mxbase))
   PAW%lmbd=0._dp
 ELSE
   nullify(PAW%lmbd)
 ENDIF
 LIBPAW_POINTER_ALLOCATE(PAW%valencemap,(Orbit%norbit))
 LIBPAW_DATATYPE_ALLOCATE(PAW%OCCwfn,)
 LIBPAW_DATATYPE_ALLOCATE(PAW%TOCCwfn,)
END SUBROUTINE InitPAW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      Subroutine DestroyPAW(PAW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine DestroyPAW(PAW)
 Type(PseudoInfo), INTENT(INOUT) :: PAW
 IF(associated(PAW%Ktvtau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Ktvtau)
 endif
 IF(associated(PAW%Krveff)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Krveff)
 endif
 if(associated(PAW%Kunscreen)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Kunscreen)
 endif
 IF (ASSOCIATED(PAW%rcio)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%rcio)
 endif
 If (ASSOCIATED(PAW%vloc)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%vloc)
 endif
 If (ASSOCIATED(PAW%abinitvloc)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%abinitvloc)
 endif
 If (ASSOCIATED(PAW%abinitnohat)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%abinitnohat)
 endif
 If (ASSOCIATED(PAW%rveff)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%rveff)
 endif
 If (ASSOCIATED(PAW%AErefrv)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%AErefrv)
 endif
 If (ASSOCIATED(PAW%rvx)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%rvx)
 endif
 If (ASSOCIATED(PAW%trvx)) then
  LIBPAW_POINTER_DEALLOCATE(PAW%trvx)
 endif
 If (ASSOCIATED(PAW%projshape)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%projshape)
 endif
 If (ASSOCIATED(PAW%hatshape)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%hatshape)
 endif
 If (ASSOCIATED(PAW%hatden)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%hatden)
 endif
 If (ASSOCIATED(PAW%hatpot)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%hatpot)
 endif
 If (ASSOCIATED(PAW%den)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%den)
 endif
 If (ASSOCIATED(PAW%tden)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tden)
 endif
 If (ASSOCIATED(PAW%core)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%core)
 endif
 If (ASSOCIATED(PAW%tcore)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tcore)
 endif
 If (ASSOCIATED(PAW%coretau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%coretau)
 endif
 If (ASSOCIATED(PAW%tcoretau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tcoretau)
 endif
 If (ASSOCIATED(PAW%valetau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%valetau)
 endif
 If (ASSOCIATED(PAW%tvaletau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tvaletau)
 endif
 If (ASSOCIATED(PAW%vtau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%vtau)
 endif
 If (ASSOCIATED(PAW%tvtau)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tvtau)
 endif
 If (ASSOCIATED(PAW%nhatv)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%nhatv)
 endif
 If (ASSOCIATED(PAW%np)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%np)
 endif
 If (ASSOCIATED(PAW%l)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%l)
 endif
 If (ASSOCIATED(PAW%nodes)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%nodes)
 endif
 If (ASSOCIATED(PAW%kappa)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%kappa)
 endif
 If (ASSOCIATED(PAW%rng)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%rng)
 endif
 If (ASSOCIATED(PAW%label)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%label)
 endif
 If (ASSOCIATED(PAW%phi)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%phi)
 endif
 If (ASSOCIATED(PAW%tphi)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tphi)
 endif
 If (ASSOCIATED(PAW%tp)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%tp)
 endif
 If (ASSOCIATED(PAW%ophi)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%ophi)
 endif
 If (ASSOCIATED(PAW%otphi)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%otphi)
 endif
 If (ASSOCIATED(PAW%otp)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%otp)
 endif
 If (ASSOCIATED(PAW%Kop)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Kop)
 endif
 If (ASSOCIATED(PAW%eig)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%eig)
 endif
 If (ASSOCIATED(PAW%occ)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%occ)
 endif
 If (ASSOCIATED(PAW%ck)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%ck)
 endif
 If (ASSOCIATED(PAW%vrc)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%vrc)
 endif
 If (ASSOCIATED(PAW%oij)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%oij)
 endif
 If (ASSOCIATED(PAW%dij)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%dij)
 endif
 If (ASSOCIATED(PAW%wij)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%wij)
 endif
 If (ASSOCIATED(PAW%rVf)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%rVf)
 endif
 If (ASSOCIATED(PAW%rtVf)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%rtVf)
 endif
 If (ASSOCIATED(PAW%g)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%g)
 endif
 If (ASSOCIATED(PAW%Kij)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Kij)
 endif
 If (ASSOCIATED(PAW%Vfij)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Vfij)
 endif
 If (ASSOCIATED(PAW%mLij)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%mLij)
 endif
 If (ASSOCIATED(PAW%DR)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%DR)
 endif
 If (ASSOCIATED(PAW%DRVC)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%DRVC)
 endif
 If (ASSOCIATED(PAW%TXVC)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%TXVC)
 endif
 If (ASSOCIATED(PAW%valencemap)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%valencemap)
 endif
 If (ASSOCIATED(PAW%lmbd)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%lmbd)
 endif
 If (ASSOCIATED(PAW%DRC)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%DRC)
 endif
 If (ASSOCIATED(PAW%DRCC)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%DRCC)
 endif
 If (ASSOCIATED(PAW%DRCjkl)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%DRCjkl)
 endif
 If (ASSOCIATED(PAW%mLic)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%mLic)
 endif
 If (ASSOCIATED(PAW%mLcc)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%mLcc)
 endif
 If (ASSOCIATED(PAW%Dcj)) then
   LIBPAW_POINTER_DEALLOCATE(PAW%Dcj)
 endif
 If (ASSOCIATED(PAW%OCCwfn)) then
   call DestroyOrbit(PAW%OCCwfn)
   LIBPAW_DATATYPE_DEALLOCATE(PAW%OCCwfn)
 end if
 If (ASSOCIATED(PAW%TOCCwfn)) then
   call DestroyOrbit(PAW%TOCCwfn)
   LIBPAW_DATATYPE_DEALLOCATE(PAW%TOCCwfn)
 end if
End Subroutine DestroyPAW





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 13. Numerov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE BoundNumerov(Grid,rv,v0,v0p,nz,l,nroot,Eig,Psi,BDsolve,success)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BoundNumerov(Grid,rv,v0,v0p,nz,l,nroot,Eig,Psi,BDsolve,success)
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: rv(:),v0,v0p
 INTEGER, INTENT(IN) :: nz,l,nroot
 real(dp), INTENT(INOUT) :: Eig(:), Psi(:,:)
 LOGICAL, INTENT(IN) :: BDsolve
 LOGICAL, INTENT(INOUT) :: success
 INTEGER, PARAMETER :: repeat=4
 INTEGER :: j
 ! TODO : BDsolve
 if(BDsolve) then
   STOP
 endif
 if(has_to_print) write(std_out,*) 'Before newboundsch',l,nroot, Eig(1:nroot)
 CALL newboundsch(Grid,rv,v0,v0p,nz,l,nroot,Eig,Psi,success)
 if(has_to_print) write(std_out,*) 'After newboundsch',l,nroot, Eig(1:nroot)
 ! adjust sign
 Do j=1,nroot
   if (Psi(3,j)<0._dp) Psi(:,j)=-Psi(:,j)
 Enddo
END SUBROUTINE BoundNumerov


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE newboundsch(Grid,rv,v0,v0p,nz,l,nroot,Eig,Psi,ok)
!!  pgm to solve radial schroedinger equation for nroot bound state
!!    energies and wavefunctions for angular momentum l
!!    with potential rv/r
!!
!!   Assymptotic from of wfn:
!!     Psi(r) = const*(r**q/kappa)*EXP(-kappa*r), where eig=-kappa**2
!!  uses Noumerov algorithm
!!
!!  For l=0,1 corrections are needed to approximate wfn(r=0)
!!     These depend upon:
!!         e0 (current guess of energy eigenvalue)
!!         l,nz
!!         v(0) == v0 electronic potential at r=0
!!         v'(0) == v0p derivative of electronic potential at r=0
!!
!!  Corrections are also needed for r>n*h, depending on:
!!         e0 (current guess of energy eigenvalue
!!         the extrapolated value of rv == r * v
!!
!! ierr=an nroot digit number indicating status of each root
!!   a digit of 1 indicates success in converging root
!!              2 indicates near success in converging root
!!              9 indicates that root not found
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE newboundsch(Grid,rv,v0,v0p,nz,l,nroot,Eig,Psi,ok)
 TYPE(GridInfo), INTENT(IN) :: Grid
 real(dp), INTENT(IN) :: rv(:),v0,v0p
 INTEGER, INTENT(IN) :: nz,l,nroot
 real(dp), INTENT(INOUT) :: Eig(:), Psi(:,:)
 LOGICAL, INTENT(OUT) :: ok
 real(dp), PARAMETER :: convre=tol10,vlrg=10._dp**30
 INTEGER, PARAMETER :: niter=300
 real(dp), ALLOCATABLE :: p1(:),p2(:),dd(:)
 INTEGER :: n,ierr
 real(dp) :: h,qq
 real(dp) :: err,convrez,energy,zeroval,zz
 real(dp) :: scale,emin,emax,best,rout,ppp
 real(dp) :: rin,dele,x
 INTEGER :: iter,i,j,node,match,ir,iroot
 INTEGER :: ifac
 ifac=0
 n=Grid%n
 h=Grid%h
 LIBPAW_ALLOCATE(p1,(n))
 LIBPAW_ALLOCATE(p2,(n))
 LIBPAW_ALLOCATE(dd,(n))
 zz=nz
 qq=-rv(n)/2
 IF (qq<0.001_dp) qq=0._dp
 qq=zero
 err=n*nz*(h**4);  if (err<tol6) err=tol6
 convrez=convre
 IF (nz.GT.0) convrez=convre*nz
 ierr=0
 emin=(-REAL((nz)**2)/(l+1)**2-0.5_dp)
 emax=0._dp
 DO iroot=1,nroot
   best=1.d10; dele=1.d10
   energy=Eig(iroot)
   IF (energy.LT.emin) energy=emin
   IF (energy.GT.emax) energy=emax
   ok=.FALSE.
   BigIter: DO iter=1,niter
     !  start inward integration
     !  start integration at n
     ! find classical turning point
     CALL ClassicalTurningPoint(Grid,rv,l,energy,match)
     match=MAX(5,match)
     match=MIN(n-15,match)
     ppp=SQRT(ABS(-energy))
     p2=0
     p2(n)=wfnend(l,energy,Grid%r(n),Grid%r(n),qq)
     p2(n-1)=wfnend(l,energy,Grid%r(n-1),Grid%r(n),qq)
     CALL backward_numerov(Grid,l,match,energy,rv,p2)
     match=match+6
     CALL derivative(Grid,p2,dd,match-5,match+5)
     rin=dd(match)/p2(match)
     !  start outward integration
     !    correct behavior near r=0
     ! initialize p1
     p1=0
     p1(2)=wfninit(-0.5_dp*rv(1),l,v0,v0p,energy,Grid%r(2))
     zeroval=0
     IF (l==0) zeroval=rv(1)
     IF (l==1) zeroval=2
     CALL forward_numerov(Grid,l,match+6,energy,rv,zeroval,p1,node)
     CALL derivative(Grid,p1,dd,match-5,match+5)
     rout=dd(match)/p1(match)
     ! check whether node = (iroot-1)
     !   not enough nodes -- raise energy
     IF (node.LT.iroot-1) THEN
       emin=MAX(emin,energy)-tol5
       energy=emax-(emax-energy)*ranx()
       ifac=9
       !   too many nodes -- lower energy
     ELSEIF (node.GT.iroot-1) THEN
       IF (energy.LT.emin) THEN
         ierr=ierr+9*(10**(iroot-1))
         if(has_to_print) WRITE(STD_OUT,*) 'newboundsch error -- emin too high',l,nz,emin,energy
         RETURN
       ENDIF
       emax=MIN(emax,energy+tol5)
       energy=emin+(energy-emin)*ranx()
       !   correct number of nodes -- estimate correction
     ELSEIF (node.EQ.iroot-1) THEN
       DO j=1,match
         p1(j)=p1(j)/p1(match)
       ENDDO
       DO j=match,n
         p1(j)=p2(j)/p2(match)
       ENDDO
       scale=1._dp/overlap(Grid,p1,p1)
       dele=(rout-rin)*scale
       x=ABS(dele)
       IF (x.LT.best) THEN
         scale=SQRT(scale)
         p1(1:n)=p1(1:n)*scale
         Psi(1:n,iroot)=p1(1:n)
         Eig(iroot)=energy
         best=x
       ENDIF
       IF (ABS(dele).LE.convrez) THEN
         if(has_to_print) WRITE(STD_OUT,*) 'converged iter with dele' , iter,dele
         ok=.TRUE.
         !  eigenvalue found
         ierr=ierr+10**(iroot-1)
         IF (iroot+1.LE.nroot) THEN
           emin=energy+tol5
           emax=0
           energy=(emin+emax)/2
           IF (energy.LT.emin) energy=emin
           IF (energy.GT.emax) energy=emax
           best=1.d10
         ENDIF
         EXIT BigIter
       ENDIF
       IF (ABS(dele).GT.convrez) THEN
         energy=energy+dele
         ! if energy is out of range, pick random energy in correct range
         IF (emin-energy.GT.convrez.OR.energy-emax.GT.convrez)         &
&             energy=emin+(emax-emin)*ranx()
         ifac=2
       ENDIF
     ENDIF
   ENDDO BigIter !iter
   IF (.NOT.ok) THEN
     ierr=ierr+ifac*(10**(iroot-1))
     if(has_to_print) WRITE(STD_OUT,*) 'no convergence in newboundsch',iroot,l,dele,energy
     if(has_to_print) WRITE(STD_OUT,*) ' best guess of eig, dele = ',Eig(iroot),best
     IF (iroot.LT.nroot) THEN
       DO ir=iroot+1,nroot
         ierr=ierr+9*(10**(ir-1))
       ENDDO
     ENDIF
     ! reset wfn with hydrogenic form
     j=iroot+l+1
     Psi(:,iroot)=0
     ppp=(j)*SQRT(ABS(Eig(iroot)))
     DO i=2,n
       Psi(i,iroot)=hwfn(ppp,j,l,Grid%r(i))
     ENDDO
   ENDIF
 ENDDO !iroot
 if(has_to_print) WRITE(STD_OUT,'("finish boundsch with eigenvalues -- ",1p,20e15.7)') &
 &    Eig(1:nroot)
 LIBPAW_DEALLOCATE(p1)
 LIBPAW_DEALLOCATE(p2)
 LIBPAW_DEALLOCATE(dd)
 if(has_to_print) WRITE(STD_OUT,*) 'returning from newboundsch -- ierr=',ierr
END SUBROUTINE newboundsch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    FUNCTION wfninit(nz,l,v0,v0p,energy,r)
!! returns the solution of the Schroedinger equation near r=0
!!  using power series expansion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION wfninit(nz,l,v0,v0p,energy,r)
 real(dp) :: wfninit
 INTEGER, INTENT(IN) :: l
 real(dp), INTENT(IN) :: nz,v0,v0p,energy,r
 real(dp) :: c1,c2,c3
 c1=-REAL(nz)/(l+1._dp)
 c2=((v0-energy)-2*nz*c1)/(4*l+6._dp)
 c3=(v0p+(v0-energy)*c1-2*nz*c2)/(6*l+12._dp)
 wfninit=(r**(l+1))*(1+r*(c1+r*(c2+r*c3)))
END FUNCTION wfninit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCTION wfnend(l,energy,r,rN)
!!
!!  Find asymptotic form of wavefunction
!!    assuming equations has form
!!   (- d^2  +l(l+1) -2q   +b^2 )
!!   (  ---   ------ ---        )  P(r) = 0
!!   (  dr^2    r^2    r        )
!!       where b^2=-energy
!!
!!        P(r) = exp(-b*(r-rN))*r^(q/b)(1+(l*(l+1)-q/b)*(q/b-1)/(2*b)/r + ...)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION wfnend(l,energy,r,rN,qq)
 real(dp) :: wfnend
 INTEGER, INTENT(IN) :: l
 real(dp), INTENT(IN) :: energy,r,rN,qq
 real(dp) :: qbb,b,wfn,cn,term,fac
 INTEGER :: i
 INTEGER, PARAMETER :: last=5
 IF (energy>=0._dp) THEN
   wfnend=zero
   RETURN
 ENDIF
 b=SQRT(-energy)
 qbb=qq/b
 cn=l*(l+1)
 fac=DDEXP(-b*(r-rN))*(r**qbb)
 term=1._dp;   wfn=zero
 DO i=1,last
   wfn=wfn+term
   IF (i<last) THEN
     term=-term*((qbb-i+1)*(qbb-i)-cn)/(2*b*i)/r
   ENDIF
 ENDDO
 wfnend=fac*wfn
END FUNCTION wfnend


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      subroutine unboundsch(Grid,rv,v0,v0p,nr,l,energy,wfn,nodes)
!!  pgm to solve radial schroedinger equation for unbound states
!!    at energy 'energy' and at angular momentum l
!!
!!    with potential rv/r, given in uniform mesh of n points
!!   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
!!                               p((n+1)*h)=0
!!  nz=nuclear charge
!!
!!  uses Noumerov algorithm
!!
!!  For l=0,1 corrections are needed to approximate wfn(r=0)
!!     These depend upon:
!!         e0 (current guess of energy eigenvalue)
!!         l,nz
!!         v(0) == v0 electronic potential at r=0
!!         v'(0) == v0p derivative of electronic potential at r=0
!!
!! also returns node == number of nodes for calculated state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unboundsch(Grid,rv,v0,v0p,nr,l,energy,wfn,nodes)
 TYPE(GridInfo), INTENT(IN)  :: Grid
 real(dp), INTENT(IN) :: rv(:),v0,v0p
 INTEGER, INTENT(IN) :: nr,l
 real(dp), INTENT(IN) :: energy
 real(dp), INTENT(INOUT) :: wfn(:)
 INTEGER, INTENT(INOUT) :: nodes
 INTEGER :: n
 real(dp) :: zeroval,scale
 n=Grid%n
 IF (nr > n) THEN
   LIBPAW_ERROR('Error in unboundsch -- nr > n')
 ENDIF
 ! initialize wfn
 wfn=0
 wfn(2)=wfninit(-0.5_dp*rv(1),l,v0,v0p,energy,Grid%r(2))
 zeroval=0
 if (l==0) zeroval=rv(1)
 if (l==1) zeroval=2
 call forward_numerov(Grid,l,nr,energy,rv,zeroval,wfn,nodes)
 ! normalize to unity within integration range
 scale=1._dp/overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
 scale=SIGN(SQRT(scale),wfn(nr-2))
 wfn(1:nr)=wfn(1:nr)*scale
END SUBROUTINE unboundsch





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 14. tools_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        SUBROUTINE extractword(wordindex,stringin,stringout)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE extractword(wordindex,stringin,stringout)
 INTEGER, INTENT(IN) :: wordindex
 CHARACTER(*), INTENT(IN) :: stringin
 CHARACTER(*), INTENT(OUT) :: stringout
 INTEGER :: i,j,n,str,fin,icount
 stringout=''
 n=LEN(stringin)
 i=INDEX(stringin,'!');IF (i==0) i=n
 j=INDEX(stringin,'#');IF (j==0) j=n
 n=MIN(i,j,n)
 str=1;fin=n
 DO icount=1,MAX(1,wordindex-1)
   DO i=str,n
     IF (stringin(i:i)/=' ') EXIT
   ENDDO
   str=i
   IF (n>str) THEN
     DO i=str+1,n
       IF(stringin(i:i)==' ') EXIT
     ENDDO
     fin=i
   ENDIF
   IF (wordindex>2) THEN
     IF (fin<n) THEN
       str=fin+1
     ELSE
       EXIT
     ENDIF
   ENDIF
 ENDDO
 IF (wordindex>1) THEN
   IF (fin>=n) RETURN
   DO i=fin+1,n
     IF (stringin(i:i)/=' ') EXIT
   ENDDO
   str=i
   IF (n>str) THEN
     DO i=str+1,n
       IF(stringin(i:i)==' ') EXIT
     ENDDO
     fin=i
   ENDIF
 ENDIF
 stringout=stringin(str:fin)
END SUBROUTINE extractword


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     SUBROUTINE UpperCase(str)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UpperCase(str)
 CHARACTER(*), INTENT(INOUT) :: str
 INTEGER  :: i, j, k
 j = LEN(Str)
 DO i=1, j
   k = IACHAR(str(i:i))
   IF ((k>96) .AND. (k<123)) str(i:i) = ACHAR(k-32)
 END DO
 RETURN
END SUBROUTINE UpperCase


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    FUNCTION checkline2(inputline,in1,in2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION checkline2(inputline,in1,in2)
 CHARACTER(*), INTENT(IN) :: inputline
 CHARACTER(*), INTENT(IN) :: in1,in2
 INTEGER :: leninput,len1,len2
 CHARACTER(120) :: inputline_u,in1_u,in2_u
 inputline_u=trim(inputline) ; call UpperCase(inputline_u)
 in1_u=trim(in1) ; call UpperCase(in1_u)
 in2_u=trim(in2) ; call UpperCase(in2_u)
 leninput=len(trim(inputline));len1=len(trim(in1));len2=len(trim(in2))
 checkline2=.false.
 if (leninput==len1) checkline2=(inputline_u(1:len1)==trim(in1))
 if ((.not.checkline2).and.leninput==len2) checkline2=(inputline_u(1:len2)==trim(in2))
 RETURN
END FUNCTION checkline2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  eliminate_comment - Eliminate comment on the right of a line (! or #)
!     line - string to convert (output replaces input)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE eliminate_comment(line)
 CHARACTER(*), INTENT(INOUT) :: line
 INTEGER :: i,i0
 i0=-1 ; i=1
 DO WHILE (i0<0.AND.i<LEN(line))
   i=i+1
   IF (line(i:i)=="!".OR.line(i:i)=="#") i0=i
 END DO
 IF (i0 >1) line=line(1:i0-1)
 IF (i0==1) line=""
END SUBROUTINE eliminate_comment 





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 15. radialked
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine unboundked(Grid,Pot,nr,l,energy,wfn,nodes)
!  pgm to solve radial kinetic energy derivative equations for unbound states
!    at energy 'energy' and at angular momentum l
!    with potential rv/r, given in uniform linear or log mesh of n points
!   assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
!  nz=nuclear chargedd
!  Does not use Noumerov algorithm -- but uses coupled first-order
!       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
! also returns node == number of nodes for calculated state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unboundked(Grid,Pot,nr,l,energy,wfn,nodes)
  TYPE(GridInfo), INTENT(IN)  :: Grid
  TYPE(PotentialInfo), INTENT(IN)  :: Pot
  INTEGER, INTENT(IN) :: nr,l
  REAL(dp), INTENT(IN) :: energy
  REAL(dp), INTENT(INOUT) :: wfn(:)
  INTEGER, INTENT(INOUT) :: nodes
  INTEGER :: n,istart
  REAL(dp) :: oneplusvtau(Grid%n),dvtaudr(Grid%n)
  REAL(dp) :: scale,qq,zxc
  REAL(dp), allocatable :: lwfn(:),zz(:,:,:),yy(:,:)
  n=Grid%n
  IF (nr > n) THEN
    write(std_out,*) 'Error in unboundked -- nr > n', nr,n
    STOP
  ENDIF
  call Set_Pot(Grid,Pot,qq,zxc,oneplusvtau,dvtaudr)
  LIBPAW_ALLOCATE(lwfn,(nr))
  LIBPAW_ALLOCATE(zz,(2,2,nr))
  LIBPAW_ALLOCATE(yy,(2,nr))
  lwfn=0;zz=0;yy=0;
  call wfnkedinit(Grid,l,Pot%nz,wfn,lwfn,istart,zxc,oneplusvtau,dvtaudr)
  call setupforcfdsol(Grid,Pot%rv,1,istart,nr,l,energy,wfn,lwfn,yy,zz,oneplusvtau,dvtaudr)
  call cfdsoliter(Grid,zz,yy,istart,nr)
  call getwfnfromcfdsol(1,nr,yy,wfn)
  nodes=countnodes(2,nr,wfn)
  ! normalize to unity within integration range
  scale=1.d0/overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
  scale=SIGN(SQRT(scale),wfn(nr-2))
  wfn(1:nr)=wfn(1:nr)*scale
  LIBPAW_DEALLOCATE(lwfn)
  LIBPAW_DEALLOCATE(yy)
  LIBPAW_DEALLOCATE(zz)
END SUBROUTINE unboundked


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set_Pot
!    Version with polynomial fitting of vtau for 0<=r<=0.001
!       and corresponding reseting of oneplusvtau and dvtaudr in that range
!       NAWH   4/6/2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Set_Pot(Grid,Pot,qq,zxc,oneplusvtau,dvtaudr)
  real(dp),intent(inout) :: qq,zxc
  real(dp),intent(inout) :: oneplusvtau(:),dvtaudr(:)
  Type(GridInfo), INTENT(IN) :: Grid
  TYPE(PotentialInfo), INTENT(IN) :: Pot
  INTEGER :: n
  REAL(dp),parameter ::  smallr=0.001d0
  INTEGER, parameter :: order=4
  n=Grid%n
  !  check for possible ionic charge
  qq=-Pot%rv(n)/two
  if(qq<0.001d0) qq=zero
  oneplusvtau=0.d0;   dvtaudr=0.d0
  oneplusvtau=1.d0+Pot%vtau
  call derivative(Grid,Pot%vtau,dvtaudr)
  zxc=Pot%rvx(1)
END Subroutine Set_Pot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE wfnkedinit(Grid,l,nz,v0,wfn,lwfn,istart)
!   returns the solution of the modified KS equations near r=0
!   using power series expansion assuming including vtau contributions
!   wfn=P(r)   lwfn=(1+vtau)*dP/dr
!   P(r)~~(r**(l+1))*(1+c1*r)
!   Assumes v(r) ~~ -2*nz/r+zxc/r   for r-->0
!   Assumes vtau(r) -- t0 +t1*r  for r-->0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wfnkedinit(Grid,l,nz,wfn,lwfn,istart,zxc,oneplusvtau,dvtaudr)
  real(dp),intent(in) :: zxc
  real(dp),intent(in) :: oneplusvtau(:),dvtaudr(:)
  Type(GridInfo), INTENT(IN) :: Grid
  INTEGER, INTENT(IN) :: l,nz
  REAL(dp),INTENT(INOUT) :: wfn(:),lwfn(:)
  INTEGER, INTENT(OUT) :: istart
  REAL(dp) :: rr,c1,t1,t0
  INTEGER :: i
  t0=oneplusvtau(1);t1=dvtaudr(1)
  wfn=zero; lwfn=zero
  c1=-(two*nz-zxc+l*t1)/(two*(l+1)*(oneplusvtau(1)))
  istart=6
  do i=1,istart
    rr=Grid%r(i)
    wfn(i)=1+rr*c1
    lwfn(i)=(rr**l)*((l+1)*wfn(i)+rr*(c1))
    lwfn(i)=(t0+t1*rr)*lwfn(i)
    wfn(i)=wfn(i)*(rr**(l+1))
  enddo
End SUBROUTINE wfnkedinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  subroutine setupforcfdsol(Grid,rv,i1,i2,n,l,energy,wfn,lwfn,yy,zz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setupforcfdsol(Grid,rv,i1,i2,n,l,energy,wfn,lwfn,yy,zz,oneplusvtau,dvtaudr)
  Type(gridinfo), INTENT(IN) :: Grid
  real(dp),intent(in) :: oneplusvtau(:),dvtaudr(:)
  INTEGER, INTENT(IN) :: i1,i2,n,l
  REAL(dp), INTENT(IN) :: energy
  REAL(dp), INTENT(IN) :: wfn(:),lwfn(:),rv(:)
  REAL(dp), INTENT(INOUT) :: yy(:,:),zz(:,:,:)
  INTEGER :: i
  REAL(dp) :: x
  x=l*(l+1)
  yy=zero;zz=zero
  yy(1,i1:i2)=wfn(i1:i2)
  yy(2,i1:i2)=lwfn(i1:i2)
  do  i=1,n
    zz(1,2,i)=1.d0/oneplusvtau(i)
    if(i==1) then
      zz(2,1,i)=0.d0
    else
      zz(2,1,i)=oneplusvtau(i)*x/(Grid%r(i)*Grid%r(i))+&
&        dvtaudr(i)/Grid%r(i)+(rv(i)/Grid%r(i)-energy)
    endif
  enddo
end subroutine setupforcfdsol





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 16. splinesolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Program uses two grids -- "universal grid" Grid
!   and modified grid Grids, with fixed Grids%n=401 and ns=400
!     and Grids%r0=0.1 which are found to work well for splinesolver
!   Internally need to omit origin and so ns=Grids%n-1
!   For MGGA case (needvtau=.true.) also need fine linear grid
!      Gridf
!  Setup local private grid  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initsplinesolver(Grid,splns,splr0,needvtau,spline)
  type(splinesolvinfo),intent(inout) :: spline
  logical,intent(in) :: needvtau
  Type(GridInfo), INTENT(IN) :: Grid       !Input universal grid
  INTEGER, INTENT(IN) :: splns  !spline grid
  REAL(dp), INTENT(IN) :: splr0  !spline r0
  INTEGER :: i,nf
  REAL(dp) :: hf
  spline%r0=splr0
  spline%ns=splns
  spline%h=log(1.d0+Grid%r(Grid%n)/spline%r0)/spline%ns
  call initgridwithn(spline%Grids,2,spline%ns+1,spline%r0,spline%h)  !local loggrid
  if(has_to_print) write(std_out,*) 'initsplinesolve ', spline%r0,spline%ns,spline%h
  LIBPAW_ALLOCATE(spline%u,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%pref,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%rr1,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%rr2,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%srv,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%svtau,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%sdvt,(spline%ns+1))
  LIBPAW_ALLOCATE(spline%soneplusvt,(spline%ns+1))
  do i=1,spline%ns+1
    spline%u(i)=(i-1)*spline%h
    spline%pref(i)=exp(0.5d0*spline%u(i))
    spline%rr1(i)=((spline%Grids%r(i)+spline%r0))
    spline%rr2(i)=((spline%Grids%r(i)+spline%r0)**2)
  enddo
  if(needvtau) then    ! set up fine linear grid
    nf=20001
    hf=Grid%r(Grid%n)/(nf-1)
    call initgridwithn(spline%Gridf,1,nf,0.d0,hf)  !local fine linear grid
    LIBPAW_ALLOCATE(spline%fvtau,(nf))
    LIBPAW_ALLOCATE(spline%fdvtaudr,(nf))
    LIBPAW_ALLOCATE(spline%frvx,(nf))
    LIBPAW_ALLOCATE(spline%fden,(nf))
    LIBPAW_ALLOCATE(spline%ftau,(nf))
    spline%fvtau=0.d0;spline%fdvtaudr=0.d0;spline%frvx=0.d0
    spline%fden=0.d0;spline%ftau=0.d0
  endif
  spline%soneplusvt=1.d0;
  spline%svtau=0.d0;spline%sdvt=0.d0
END SUBROUTINE initsplinesolver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! SUBROUTINE initpotforsplinesolver
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initpotforsplinesolver(Grid,Pot,den,tau,spline,itype)
  type(splinesolvinfo),intent(inout) :: spline
  Type(GridInfo), INTENT(IN) :: Grid
  integer,intent(in) :: itype
  Type(Potentialinfo), INTENT(INOUT) :: Pot   !Universal grid
  REAL(dp), INTENT(IN) :: den(:),tau(:)  !Universal grid
  INTEGER :: i,n
  REAL(dp), allocatable :: dum(:),dum1(:)
  REAL(dp) :: etxc,eexc
  n=Grid%n
  if(has_to_print) write(std_out,*) 'initpot ', n,spline%ns
  If (Pot%needvtau) then
    LIBPAW_ALLOCATE(dum,(spline%ns+1))
    LIBPAW_ALLOCATE(dum1,(spline%ns+1))
    call interpfunc(n,Grid%r,den,spline%Gridf%n,spline%Gridf%r,spline%fden)
    call interpfunc(n,Grid%r,tau,spline%Gridf%n,spline%Gridf%r,spline%ftau)
    call exch(spline%Gridf,spline%fden,spline%frvx,etxc,eexc,itype=itype,tau=spline%ftau,vtau=spline%fvtau)
    if(has_to_print) write(std_out,*) 'called exch from splinesolver ', eexc
    call nderiv(spline%Gridf%h,spline%fvtau,spline%fdvtaudr,spline%Gridf%n,i)
    dum=Pot%rvn+Pot%rvh    !  presumably these are smooth
    call interpfunc(n,Grid%r,dum,spline%ns+1,spline%Grids%r,spline%srv)
    call interpfunc(spline%Gridf%n,spline%Gridf%r,spline%frvx,spline%Grids%n,spline%Grids%r,dum1)
    spline%srv=spline%srv+dum1
    call interpfunc(spline%Gridf%n,spline%Gridf%r,spline%fvtau,spline%Grids%n,spline%Grids%r,spline%svtau)
    call interpfunc(spline%Gridf%n,spline%Gridf%r,spline%fdvtaudr,spline%Grids%n,spline%Grids%r,spline%sdvt)
    spline%soneplusvt=1.d0+spline%svtau
    spline%sdvt=spline%sdvt*spline%Grids%drdu      ! needed in algorithm
    LIBPAW_DEALLOCATE(dum)
    LIBPAW_DEALLOCATE(dum1)
  else    !  finegrid not needed        
    call interpfunc(n,Grid%r,Pot%rv,spline%ns+1,spline%Grids%r,spline%srv)
    spline%soneplusvt=1.d0;
    spline%svtau=0.d0;spline%sdvt=0.d0
  endif
END SUBROUTINE initpotforsplinesolver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE Boundsplinesolver(Grid,l,neig,eig,wfn,otau,OK)
!!   Note that although the universal Grid and the local Grids
!!     have the same range, they differ by the number of points and
!!     the r0 parameter.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Boundsplinesolver(Grid,l,neig,eig,wfn,otau,OK,spline)
  type(splinesolvinfo),intent(inout) :: spline
  Type(GridInfo), INTENT(IN) :: Grid
  INTEGER, INTENT(IN) :: l,neig
  REAL(dp), INTENT(INOUT) :: eig(:),wfn(:,:),otau(:,:)
  LOGICAL, INTENT(OUT) :: OK
  real(dp), allocatable :: A(:,:),B(:,:),G(:,:),D(:),DL(:),DU(:),work(:)
  real(dp), allocatable :: vl(:,:),vr(:,:),wr(:),wi(:)
  integer, allocatable :: lut(:)
  REAL(dp), allocatable :: P(:),MP(:)
  REAL(dp), allocatable :: dum(:),dP1(:)
  REAL(dp), allocatable :: F1(:,:),F2(:,:),S(:),E(:),V(:)
  integer :: i,m,n, info,lwork,nu
  real(dp) :: ol,x
  if(has_to_print) write(std_out,*) 'entering boundspline with l,neig ', l, neig
  if(has_to_print) write(std_out,*) 'in splinesolver '
  OK=.false.
  n=spline%ns+1       !Grids%n   within this routine only
  nu=Grid%n     ! Universal grid
  ol=l*(l+1)
  lwork=spline%ns**2
  LIBPAW_ALLOCATE(A,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(B,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(G,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(D,(spline%ns))
  LIBPAW_ALLOCATE(DL,(spline%ns))
  LIBPAW_ALLOCATE(DU,(spline%ns))
  LIBPAW_ALLOCATE(work,(lwork))
  LIBPAW_ALLOCATE(vl,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(vr,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(wr,(spline%ns))
  LIBPAW_ALLOCATE(wi,(spline%ns))
  LIBPAW_ALLOCATE(lut,(spline%ns))
  LIBPAW_ALLOCATE(P,(n))
  LIBPAW_ALLOCATE(MP,(n))
  LIBPAW_ALLOCATE(dum,(nu))
  LIBPAW_ALLOCATE(dP1,(nu))
  LIBPAW_ALLOCATE(F1,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(F2,(spline%ns,spline%ns))
  LIBPAW_ALLOCATE(S,(n))
  LIBPAW_ALLOCATE(E,(n))
  LIBPAW_ALLOCATE(V,(n))
  D(1:spline%ns)=2.d0
  DL(1:spline%ns-1)=0.5d0
  DU(1:spline%ns-1)=0.5d0
  !correct first row assuming wfn=r^(l+1)*(W0+r*W1)
  D(1)=0.5d0*(5.d0+two*l)/(1.d0+l)
  B=0.d0
  do i=1,spline%ns
    B(i,i)=1.d0
  enddo   
  call dgtsv(spline%ns,spline%ns,DL,D,DU,B,spline%ns,info)
  if(has_to_print) write(std_out,*) 'Completed dgtsv with info = '
  A=0.d0
  do i=1,spline%ns
    A(i,i)=-2.d0
  enddo
  do i=1,spline%ns-1
    A(i,i+1)=1.d0
    A(i+1,i)=1.d0
  enddo   
  ! correct first row values
  A(1,1)=-0.5d0*(l+4.d0);  
  A=3.0_dp*A/(spline%h**2)
  G=0.d0
  G=MATMUL(B,A)      !  G stores transformation to find M=G*y
  ! Calculate full matrix
  A=0.0_dp;S=0.0_dp;E=0.0_dp;V=0.0_dp
  ! These arrays  have the full range  1..n
  S=-spline%soneplusvt/spline%rr2
  E=-spline%sdvt/spline%rr2
  V=0.5d0*E-0.25d0*S     !Including only non diverging term
  do i=2,n
    V(i)=V(i)+spline%soneplusvt(i)*ol/(spline%Grids%r(i)**2)&
&           +spline%sdvt(i)/(spline%Grids%r(i)*spline%rr1(i)) &
&           +spline%srv(i)/spline%Grids%r(i)
  enddo
  F1=0.d0;F2=0.d0
  do i=1,spline%ns
    F1(i,i)=S(i+1)-E(i+1)*spline%h/three
    F1(i,i+1)=-E(i+1)*spline%h/six       
    F2(i,i)=V(i+1)-E(i+1)/spline%h
    F2(i,i+1)=E(i+1)/spline%h
  ENDDO   
  A=MATMUL(F1,G)+F2
  call dgeev('N','V',spline%ns,A,spline%ns,wr,wi,vl,spline%ns,vr,spline%ns,work,lwork,info)
  if(has_to_print) write(std_out,*) 'dgeev completed with info = ',info
  if (info/=0) then
    OK=.false.
    return
  endif    
  call real_InsSort(wr,lut,.true.)      
  if(has_to_print) write(std_out,*) 'Results  for l = ', l, (wr(lut(i)),i=1,10)
  if(has_to_print) write(std_out,*)  'enumeration for neig solutions ', neig
  do m=1,neig
    if(has_to_print) write(std_out,'(1p,2e20.8)') wr(lut(m)),wi(lut(m))
    eig(m)=wr(lut(m))
    D=vr(:,lut(m))     !Q
    DL=MATMUL(G,D)     !MQ
    ! now extend grid to r=0   Still Q and MQ
    P=0.0_dp;MP=0.0_dp
    P(2:spline%ns+1)=D(1:spline%ns)
    MP(2:spline%ns+1)=DL(1:spline%ns)
    if(l==0) then
      dum=0      
      do i=2,10      
        dum(i)=P(i)/spline%Grids%r(i)
      enddo   
      call extrapolate(dum) 
      x=1.d0/(S(1)-E(1)*spline%h/three)
      MP(1)=(E(1)*(MP(2)*spline%h/six-P(2))-(spline%sdvt(1)/spline%rr1(1)+spline%srv(1))*dum(1))*x
      if(has_to_print) write(std_out,*) 'MP(1) for l=0',MP(1),dum(1)
    endif
    if(l==1) then
      dum=0.0_dp      
      do i=2,10      
        dum(i)=P(i)/(spline%Grids%r(i)**2)
      enddo   
      call extrapolate(dum) 
      x=1.d0/(S(1)-E(1)*spline%h/three)
      MP(1)=(E(1)*(MP(2)*spline%h/six-P(2))-(two*spline%soneplusvt(1))*dum(1))*x
      if(has_to_print) write(std_out,*) 'MP(1) for l=1',MP(1),dum(1)
    endif
    MP=spline%pref*(MP-0.25d0*P)/spline%rr2
    P=spline%pref*P
    call specialinterp(n,spline%Grids%r,P,MP,Grid%n,Grid%r,wfn(:,m),dP1(:))
    dum=0.d0; dum(2:nu)=wfn(2:nu,m)/Grid%r(2:nu)
    call extrapolate(dum)
    do i=1,nu
      otau(i,m)= (dP1(i)-dum(i))**2+l*(l+1)*(dum(i))**2
    enddo   
    x=overlap(Grid,wfn(:,m),wfn(:,m))
    if(has_to_print) write(std_out,*) 'overlap integral ', x
    x=1.d0/x
    otau(:,m)=x*otau(:,m)
    x=sqrt(x)
    wfn(:,m)=x*wfn(:,m)    ! should be normalized now
  enddo    
  OK=.true.
  LIBPAW_DEALLOCATE(A)
  LIBPAW_DEALLOCATE(B)
  LIBPAW_DEALLOCATE(G)
  LIBPAW_DEALLOCATE(D)
  LIBPAW_DEALLOCATE(DL)
  LIBPAW_DEALLOCATE(DU)
  LIBPAW_DEALLOCATE(work)
  LIBPAW_DEALLOCATE(vl)
  LIBPAW_DEALLOCATE(vr)
  LIBPAW_DEALLOCATE(wr)
  LIBPAW_DEALLOCATE(wi)
  LIBPAW_DEALLOCATE(lut)
  LIBPAW_DEALLOCATE(P)
  LIBPAW_DEALLOCATE(MP)
  LIBPAW_DEALLOCATE(dP1)
  LIBPAW_DEALLOCATE(dum)
  LIBPAW_DEALLOCATE(F1)
  LIBPAW_DEALLOCATE(F2)
  LIBPAW_DEALLOCATE(S)
  LIBPAW_DEALLOCATE(E)
  LIBPAW_DEALLOCATE(V)
end SUBROUTINE Boundsplinesolver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   SUBROUTINE deallocatesplinesolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE deallocatesplinesolver(spline,needvtau)
  type(splinesolvinfo),intent(inout) :: spline
  logical,intent(in) :: needvtau
  if (allocated(spline%u)) then 
    LIBPAW_DEALLOCATE(spline%u)
  endif
  if (allocated(spline%pref)) then
    LIBPAW_DEALLOCATE(spline%pref)
  endif
  if (allocated(spline%rr1)) then
    LIBPAW_DEALLOCATE(spline%rr1)
  endif
  if (allocated(spline%rr2)) then
    LIBPAW_DEALLOCATE(spline%rr2)
  endif
  if (allocated(spline%srv)) then
    LIBPAW_DEALLOCATE(spline%srv)
  endif
  if (allocated(spline%svtau)) then
    LIBPAW_DEALLOCATE(spline%svtau)
  endif
  if (allocated(spline%soneplusvt)) then
    LIBPAW_DEALLOCATE(spline%soneplusvt)
  endif
  call destroygrid(spline%Grids)
  if(needvtau) then
    LIBPAW_DEALLOCATE(spline%fvtau)
    LIBPAW_DEALLOCATE(spline%fdvtaudr)
    LIBPAW_DEALLOCATE(spline%frvx)
    LIBPAW_DEALLOCATE(spline%fden)
    LIBPAW_DEALLOCATE(spline%ftau)
    call destroygrid(spline%Gridf)
  endif
END SUBROUTINE deallocatesplinesolver



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 17. radialdirac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SUBROUTINE boundD(Grid,Pot,eig,wfn,lwfn,kappa,nroot,emin,ierr,success)
!!  pgm to solve radial Dirac relativistic equation for nroot bound state
!!    energies and wavefunctions for spin orbit parameter kappa
!!    with potential rv/r, given in uniform linear or log mesh of n points
!!  nz=nuclear charge
!!  emin=is estimate of lowest eigenvalue; used if nz=0
!!     otherwise, set to the value of -(nz/(l+1))**2
!!  It is assumed that the wavefunction has np-l-1 nodes, where
!!    np is the principal quantum number-- np=1,2,..nroot
!!  Does not use Noumerov algorithm -- but uses coupled first-order
!!       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
!!  Corrections are also needed for r>n*h, depending on:
!!         e0 (current guess of energy eigenvalue
!!         the extrapolated value of rv == r * v
!! ierr=an nroot digit number indicating status of each root
!!   a digit of 1 indicates success in converging root
!!              2 indicates near success in converging root
!!              9 indicates that root not found
!! first check how many roots expected =  ntroot (returned as argument)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE boundD(Grid,Pot,eig,wfn,lwfn,kappa,nroot,emin,ierr,success)
  TYPE(GridInfo), INTENT(IN) :: Grid
  TYPE(PotentialInfo), INTENT(INout) :: Pot
  REAL(dp), INTENT(INOUT) :: eig(:),wfn(:,:),lwfn(:,:)
  INTEGER, INTENT(IN) :: kappa,nroot
  INTEGER, INTENT(INOUT) :: ierr
  REAL(dp), INTENT(INOUT) :: emin
  LOGICAL, INTENT(INOUT) :: success
  REAL(dp), PARAMETER :: convre=1.d-10,vlrg=1.d30
  REAL(dp), PARAMETER :: ftr=0.5d0/InvFineStruct
  INTEGER, PARAMETER :: niter=1000
  REAL(dp), POINTER :: rv(:)
  REAL(dp), ALLOCATABLE :: p1(:),lp1(:),p2(:),lp2(:),dd(:)
  INTEGER :: n
  REAL(dp) :: nz,h,v0,v0p
  REAL(dp) :: err,convrez,energy
  REAL(dp) :: scale,emax,best,rout
  REAL(dp) :: arg,rin,dele,x
  INTEGER :: iter,i,j,node,match,mxroot,ntroot,ir,iroot,l
  INTEGER :: ifac,istart,iend
  LOGICAL :: ok
  Real(dp) :: A0,B0,A1,B1,s
  REAL(dp), allocatable :: zz(:,:,:),yy(:,:)
  n=Grid%n
  h=Grid%h
  if (kappa<0)  l=-kappa-1
  LIBPAW_ALLOCATE(p1,(n))
  LIBPAW_ALLOCATE(lp1,(n))
  LIBPAW_ALLOCATE(p2,(n))
  LIBPAW_ALLOCATE(lp2,(n))
  LIBPAW_ALLOCATE(dd,(n))
  LIBPAW_ALLOCATE(zz,(2,2,n))
  LIBPAW_ALLOCATE(yy,(2,n))
  success=.true.
  nz=Pot%nz
  v0=Pot%v0
  v0p=Pot%v0p
  rv=>Pot%rv
  err=n*nz*(h**4)
  convrez=convre
  IF (nz>0.001d0) convrez=convre*nz
  ierr=0
  if(has_to_print) WRITE(STD_OUT,*) 'z ,kappa, l = ',nz,kappa,l
  ! check how many roots expected by integration outward at
  !   energy = 0
  energy = zero
  call Dzeroexpand(Grid,Pot,kappa,energy,A0,A1,B0,B1,s)
  zz=zero;yy=zero
  call wfnDinit(Grid,p1,lp1,istart,A0,A1,B0,B1,s,Pot%finitenucleus)
  !start outward integration
  call prepareforcfdsolD(Grid,1,istart,n,kappa,p1,lp1,yy,zz,Pot%ww,Pot%jj)
  call cfdsoliter(Grid,zz,yy,istart,n)
  call getwfnfromcfdsolD(1,n,yy,p1,lp1)
  node=countnodes(2,n,p1)
  if(has_to_print) WRITE(STD_OUT,*) ' nodes at e=0  ', node
  mxroot=node+1
  ntroot=node
  IF (mxroot.LT.nroot) THEN
    if(has_to_print) WRITE(STD_OUT,*)'error in boundD - for l = ',l
    if(has_to_print) WRITE(STD_OUT,*) nroot,' states requested but only',mxroot,' possible'
    DO ir=mxroot+1,nroot
      ierr=ierr+9*(10**(ir-1))
    ENDDO
    success=.false.
  ENDIF
  mxroot=min0(mxroot,nroot)
  IF (nz.EQ.0) energy=-ABS(emin)
  IF (nz.NE.0) energy=-1.1d0*(nz/(l+1.d0))**2
  emin=energy-err
  emax=0.d0
  DO iroot=1,mxroot
    best=1.d10; dele=1.d10
    energy=emin+err
    IF (energy.LT.emin) energy=emin
    IF (energy.GT.emax) energy=emax
    ok=.FALSE.
    BigIter: DO iter=1,niter
      !  start inward integration
      !  start integration at n
      call Dzeroexpand(Grid,Pot,kappa,energy,A0,A1,B0,B1,s)
      ! find classical turning point
      call ClassicalTurningPoint(Grid,Pot%rv,l,energy,match)
      match=max(match,10); match=min(match,n-20)
      call wfnDasym(Grid,p2,lp2,energy,iend)
      call prepareforcfdsolD(Grid,n-iend,n,n,kappa,p2,lp2,yy,zz,Pot%ww,Pot%jj)
      call cfdsoliter(Grid,zz,yy,n-iend,match)
      call getwfnfromcfdsolD(match,n,yy,p2,lp2)
      match=match+6
      rin=lp2(match)/p2(match)
      call wfnDinit(Grid,p1,lp1,istart,A0,A1,B0,B1,s,Pot%finitenucleus)
      call prepareforcfdsolD(Grid,1,istart,n,kappa,p1,lp1,yy,zz,Pot%ww,Pot%jj)
      call cfdsoliter(Grid,zz,yy,istart,match+6)
      call getwfnfromcfdsolD(1,match+6,yy,p1,lp1)
      node= countnodes(2,match+6,p1)
      rout=lp1(match)/p1(match)
      ! check whether node = (iroot-1)
      !   not enough nodes -- raise energy
      IF (node.LT.iroot-1) THEN
        emin=MAX(emin,energy)-err
        energy=emax-(emax-energy)*ranx()
        ifac=9
        !   too many nodes -- lower energy
      ELSEIF (node.GT.iroot-1) THEN
        IF (energy.LE.emin) THEN
          ierr=ierr+9*(10**(iroot-1))
          if(has_to_print) WRITE(STD_OUT,*) 'boundD error -- emin too high',l,nz,emin,energy
          IF (energy.LE.emin-1.d-10) THEN
            STOP
          ENDIF
        ENDIF
        emax=MIN(emax,energy+err)
        energy=emin+(energy-emin)*ranx()
        !   correct number of nodes -- estimate correction
      ELSEIF (node.EQ.iroot-1) THEN
        DO j=1,match
          p1(j)=p1(j)/p1(match)
          lp1(j)=ftr*lp1(j)/p1(match)
        ENDDO
        DO j=match,n
          p1(j)=p2(j)/p2(match)
          lp1(j)=ftr*lp2(j)/p2(match)
        ENDDO
        scale=overlap(Grid,p1,p1)+overlap(Grid,lp1,lp1)
        dele=(rout-rin)/scale
        x=ABS(dele)
        IF (x.LT.best) THEN
          scale=1.d0/SQRT(scale)
          p1(1:n)=p1(1:n)*scale
          lp1(1:n)=lp1(1:n)*scale
          call filter(n,p1,machine_zero)
          call filter(n,lp1,machine_zero)
          wfn(1:n,iroot)=p1(1:n)
          lwfn(1:n,iroot)=lp1(1:n)
          eig(iroot)=energy
          best=x
        ENDIF
        IF (ABS(dele).LE.convrez) THEN
          ok=.TRUE.
          !  eigenvalue found
          ierr=ierr+10**(iroot-1)
          IF (iroot+1.LE.mxroot) THEN
            emin=energy+err
            emax=zero
            energy=(emin+emax)/2
            IF (energy.LT.emin) energy=emin
            IF (energy.GT.emax) energy=emax
            best=1.d10
          ENDIF
          EXIT BigIter
        ENDIF
        IF (ABS(dele).GT.convrez) THEN
          energy=energy+dele
          ! if energy is out of range, pick random energy in correct range
          IF (emin-energy.GT.convrez.OR.energy-emax.GT.convrez)         &
&              energy=emin+(emax-emin)*ranx()
          ifac=2
        ENDIF
      ENDIF
    ENDDO BigIter !iter
    IF (.NOT.ok) THEN
      success=.false.     
      ierr=ierr+ifac*(10**(iroot-1))
      if(has_to_print) WRITE(STD_OUT,*) 'no convergence in boundD',iroot,l,dele,energy
      if(has_to_print) WRITE(STD_OUT,*) ' best guess of eig, dele = ',eig(iroot),best
      IF (iroot.LT.mxroot) THEN
        DO ir=iroot+1,mxroot
          ierr=ierr+9*(10**(ir-1))
        ENDDO
      ENDIF
      ! reset wfn with hydrogenic form
      j=iroot+l+1
      wfn(:,iroot)=zero
      lwfn(:,iroot)=zero
      x=(j)*sqrt(abs(eig(iroot)))
      do i=2,n
        call dirachwfn(j,kappa,x,Grid%r(i),arg,wfn(i,iroot),lwfn(i,iroot))
      enddo
    ENDIF
  ENDDO !iroot
  LIBPAW_DEALLOCATE(p1)
  LIBPAW_DEALLOCATE(lp1)
  LIBPAW_DEALLOCATE(p2)
  LIBPAW_DEALLOCATE(lp2)
  LIBPAW_DEALLOCATE(dd)
  LIBPAW_DEALLOCATE(yy)
  LIBPAW_DEALLOCATE(zz)
  if(has_to_print) write(std_out,*) 'returning from boundD -- ierr=',ierr
END SUBROUTINE BoundD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine prepareforcfdsolD(Grid,i1,i2,n,kappa,wfn,lwfn,yy,zz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prepareforcfdsolD(Grid,i1,i2,n,kappa,wfn,lwfn,yy,zz,ww,jj)
  real(dp),intent(in) :: ww(:),jj(:)
  Type(gridinfo), INTENT(IN) :: Grid
  INTEGER, INTENT(IN) :: i1,i2,n,kappa
  REAL(dp), INTENT(IN) :: wfn(:),lwfn(:)
  REAL(dp), INTENT(OUT) :: yy(:,:),zz(:,:,:)
  INTEGER :: i
  yy=zero;zz=zero
  yy(1,i1:i2)=wfn(i1:i2)
  yy(2,i1:i2)=lwfn(i1:i2)
  do  i=2,n
    zz(1,1,i)=-kappa/Grid%r(i)
    zz(1,2,i)=jj(i)
    zz(2,2,i)=kappa/Grid%r(i)
    zz(2,1,i)=-ww(i)
  enddo
end subroutine prepareforcfdsolD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Subroutine Dzeroexpand(Grid,Pot,kappa,energy)
!!      If finitenucleus==.true. assumes potential is non-singular
!!          at origin and Pot%v0 and Pot%v0p are properly set
!!          -- actually not programmed yet
!!      Otherwise, assumes nuclear potential is -2*Z/r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Dzeroexpand(Grid,Pot,kappa,energy,A0,A1,B0,B1,s,nr)
  Real(dp), intent(inout) :: A0,B0,A1,B1,s
  Type(GridInfo), INTENT(IN) :: Grid
  Type(PotentialInfo), INTENT(INout) :: Pot
  Integer, INTENT(IN) :: kappa
  Real(dp), INTENT(IN) :: energy
  Integer, optional, INTENT(IN) :: nr
  Integer :: n
  Real(dp) :: nz,alpha2,balpha2
  Real(dp) :: x,y,z
  n=Grid%n
  if (present(nr)) n=min(n,nr)
  nz=Pot%nz
  Pot%ww=zero; Pot%jj=zero
  balpha2=InvFineStruct**2
  alpha2=1.d0/balpha2
  Pot%ww(2:n)=energy-Pot%rv(2:n)/Grid%r(2:n)
  Pot%jj(2:n)=(1.d0 + 0.25d0*alpha2*Pot%ww(2:n))
  if (.not.Pot%finitenucleus) then
    s=sqrt(kappa*kappa-alpha2*nz**2)
    A0=1.d0
    B0=2.d0*(s+kappa)*balpha2/nz
    z=two*s+one
    x=alpha2*(nz**2)
    y=four*alpha2+energy-Pot%v0
    A1=(four*alpha2*x+y*(s+kappa-two*x))/(two*nz*z)
    y=two*alpha2+energy-Pot%v0
    B1=-(y*(two*(s+kappa)+energy-Pot%v0))/z
  else  ! version for finite nuclear size
    write(std_out,*) 'Dirac case not yet programmed for finite nucleus'
    stop
  endif
end subroutine Dzeroexpand


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE wfnDinit(Grid,kappa,wfn,lwfn,istart)
!! returns the solution of the Dirac relativistic equations near r=0
!!  using power series expansion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE wfnDinit(Grid,wfn,lwfn,istart,A0,A1,B0,B1,s,finitenucleus)
  logical,intent(in) :: finitenucleus
  Real(dp), intent(in) :: A0,B0,A1,B1,s
  Type(GridInfo), INTENT(IN) :: Grid
  REAL(dp),INTENT(INOUT) :: wfn(:),lwfn(:)
  INTEGER, INTENT(OUT) :: istart
  REAL(dp) :: rr
  INTEGER :: i
  wfn=zero; lwfn=zero
  istart=6
  do i=1,istart
    rr=Grid%r(i+1)
    if (.not.finitenucleus) then
      wfn(i+1)=(rr**s)*(A0+A1*rr)
      lwfn(i+1)=(rr**s)*(B0+B1*rr)
    else   ! finite nucleus case
      write(std_out,*) 'Dirac case not programmed for finite nucleus'
      STOP
    endif
  enddo
End SUBROUTINE wfnDinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine wfnDasym(Grid,wfn,lwfn,energy,iend)
!! returns the solution of the Dirac relativistic equations near r=inf
!!  using exp(-x*r) for upper component
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wfnDasym(Grid,wfn,lwfn,energy,iend)
  Type(GridInfo), INTENT(IN) :: Grid
  REAL(dp),INTENT(INOUT) :: wfn(:),lwfn(:)
  REAL(dp), INTENT(IN) :: energy
  INTEGER, INTENT(OUT) :: iend
  REAL(dp) :: rr,x,m
  INTEGER :: i,n
  if (energy>0.d0) then
    write(std_out,*) 'Error in wfnDasym -- energy > 0', energy
    stop
  endif
  wfn=zero; lwfn=zero
  n=Grid%n
  m=1.d0+0.25d0*energy/(InvFineStruct**2)
  x=sqrt(-m*energy)
  rr=energy/x
  iend=5
  do i=n-iend,n
    wfn(i)=exp(-x*(Grid%r(i)-Grid%r(n-iend)))
    lwfn(i)=rr*wfn(i)
  enddo
end subroutine wfnDasym


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine getwfnfromcfdsolD(start,finish,yy,wfn,lwfn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getwfnfromcfdsolD(start,finish,yy,wfn,lwfn)
  INTEGER, INTENT(IN) :: start,finish
  REAL(dp), INTENT(IN) :: yy(:,:)
  REAL(dp), INTENT(INOUT) :: wfn(:),lwfn(:)
  INTEGER :: i
  wfn=0
  do i=start,finish
    wfn(i)=yy(1,i)
    lwfn(i)=yy(2,i)
  enddo
end subroutine getwfnfromcfdsolD





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 18. interpolation_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!    cubic spline interpolation for nin>2
!    linear interpolation for nin=2
!    error end if rout points are not within range of rin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
SUBROUTINE interpfunc(nin,rin,fin,nout,rout,fout)
  INTEGER, INTENT(IN) :: nin,nout
  REAL(dp), INTENT(IN) :: rin(:),fin(:),rout(:)
  REAL(dp), INTENT(INOUT) :: fout(:)
  REAL(dp), ALLOCATABLE :: c(:,:)
  INTEGER :: i,j
  REAL(dp) :: x,s
  LOGICAL :: leftin,lefton,leftout,rightin,righton,rightout
  !  check if grid is within interpolation range
  call checkrange(rin(1),rout(1),leftin,lefton,leftout)
  call checkrange(rout(nout),rin(nin),rightin,righton,rightout)
  if (leftout.or.rightout) then
    write(std_out,*) 'Grid error in interpfunc',rin(1),rout(1),rin(nin),rout(nout)
    stop
  endif
  fout=0
  ! linear interpolation if nin=2
  if (nin==2) then
    s=(fin(2)-fin(1))/(rin(2)-rin(1))
    do i=1,nout
      fout(i)=fin(1)+(rout(i)-rin(1))*s
    enddo
    return
  endif
  LIBPAW_ALLOCATE(c,(4,nin))
  c=zero;
  c(1,1:nin)=fin(1:nin)
  call cubspl(rin,c,nin,0,0)  
  do i=1,nout
    do j=1,nin-1
      call checkrange(rin(j),rout(i),leftin,lefton,leftout)
      call checkrange(rout(i),rin(j+1),rightin,righton,rightout)
       if ((leftin.or.lefton).and.(rightin.or.righton)) then
         x=rout(i)-rin(j)
         fout(i)=c(1,j)+x*(c(2,j)+x*(c(3,j)+x*c(4,j)/3)/2)
         exit
       endif
    enddo
  enddo
  LIBPAW_DEALLOCATE(c)
END SUBROUTINE interpfunc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  cubspl
!! from webpage
!!   http://pages.cs.wisc.edu/~deboor/pgs/cubspl.f
!!   Transformed into fortran 90 and REAL(8)
!!  from  * a practical guide to splines *  by c. de boor
!!     ************************  input  ***************************
!!     n = number of data points. assumed to be .ge. 2.
!!     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
!!        data points. tau is assumed to be strictly increasing.
!!     ibcbeg, ibcend = boundary condition indicators, and
!!     c(2,1), c(2,n) = boundary condition information. specifically,
!!        ibcbeg = 0  means no boundary condition at tau(1) is given.
!!           in this case, the not-a-knot condition is used, i.e. the
!!           jump in the third derivative across tau(2) is forced to
!!           zero, thus the first and the second cubic polynomial pieces
!!           are made to coincide.)
!!        ibcbeg = 1  means that the slope at tau(1) is made to equal
!!           c(2,1), supplied by input.
!!        ibcbeg = 2  means that the second derivative at tau(1) is
!!           made to equal c(2,1), supplied by input.
!!        ibcend = 0, 1, or 2 has analogous meaning concerning the
!!           boundary condition at tau(n), with the additional infor-
!!           mation taken from c(2,n).
!!     ***********************  output  **************************
!!     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
!!        of the cubic interpolating spline with interior knots (or
!!        joints) tau(2), ..., tau(n-1). precisely, in the interval
!!        (tau(i), tau(i+1)), the spline f is given by
!!           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!!        where h = x - tau(i). the function program *ppvalu* may be
!!        used to evaluate f or its derivatives from tau,c, l = n-1,
!!        and k=4.
!!****** a tridiagonal linear system for the unknown slopes s(i) of
!!  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
!!  ination, with s(i) ending up in c(2,i), all i.
!!     c(3,.) and c(4,.) are used initially for temporary storage.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
  integer :: ibcbeg,ibcend,n,   i,j,l,m
  real(dp) :: c(4,n),tau(n),   divdf1,divdf3,dtau,g
  l = n-1
  !compute first differences of tau sequence and store in c(3,.). also,
  !compute first divided difference of data and store in c(4,.).
  do  m=2,n
    c(3,m) = tau(m) - tau(m-1)
    c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
  enddo
  !construct first equation from the boundary condition, of the form
  !             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
  if (ibcbeg-1 <0)                  go to 11
  if (ibcbeg-1==0)                  go to 15
  if (ibcbeg-1 >0)                  go to 16
  11 if (n .gt. 2)                     go to 12
  !     no condition at left end and n = 2.
  c(4,1) = 1.d0
  c(3,1) = 1.d0
  c(2,1) = 2.d0*c(4,2)
  go to 25
 !     not-a-knot condition at left end and n .gt. 2.
  12 c(4,1) = c(3,3)
  c(3,1) = c(3,2) + c(3,3)
  c(2,1) =((c(3,2)+2.d0*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))/c(3,1)
  go to 19
  !     slope prescribed at left end.
  15 c(4,1) = 1.d0
  c(3,1) = 0.d0
  go to 18
  !     second derivative prescribed at left end.
  16 c(4,1) = 2.d0
  c(3,1) = 1.d0
  c(2,1) = 3.d0*c(4,2) - c(3,2)/2.d0*c(2,1)
  18 if(n .eq. 2)                      go to 25
  !  if there are interior knots, generate the corresp. equations and car-
  !  ry out the forward pass of gauss elimination, after which the m-th
  !  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
  19 do m=2,l
       g = -c(3,m+1)/c(4,m-1)
       c(2,m) = g*c(2,m-1) + 3.d0*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
       c(4,m) = g*c(3,m-1) + 2.d0*(c(3,m) + c(3,m+1))
  enddo
  !construct last equation from the second boundary condition, of the form
  !           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
  !     if slope is prescribed at right end, one can go directly to back-
  !     substitution, since c array happens to be set up just right for it
  !     at this point.
  if (ibcend-1 <0)                  go to 21
  if (ibcend-1==0)                  go to 30
  if (ibcend-1 >0)                  go to 24
  21 if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
  !     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
  !     left end point.
  g = c(3,n-1) + c(3,n)
  c(2,n) = ((c(3,n)+2.d0*g)*c(4,n)*c(3,n-1) &
&             + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
  g = -g/c(4,n-1)
  c(4,n) = c(3,n-1)
  go to 29
  !     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
  !     knot at left end point).
  22 c(2,n) = 2.d0*c(4,n)
  c(4,n) = 1.d0
  go to 28
  !     second derivative prescribed at right endpoint.
  24 c(2,n) = 3.d0*c(4,n) + c(3,n)/2.d0*c(2,n)
  c(4,n) = 2.d0
  go to 28
  25 continue
  if (ibcend-1 <0)                  go to 26
  if (ibcend-1==0)                  go to 30
  if (ibcend-1 >0)                  go to 24
  26 if (ibcbeg .gt. 0)                go to 22
  !     not-a-knot at right endpoint and at left endpoint and n = 2.
  c(2,n) = c(4,n)
  go to 30
  28 g = -1.d0/c(4,n-1)
  !complete forward pass of gauss elimination.
  29 c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
!carry out back substitution
   30 j = l
   40    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
   j = j - 1
   if (j .gt. 0)                  go to 40
   !****** generate cubic coefficients in each interval, i.e., the deriv.s
   !  at its left endpoint, from value and slope at its endpoints.
   do  i=2,n
     dtau = c(3,i)
     divdf1 = (c(1,i) - c(1,i-1))/dtau
     divdf3 = c(2,i-1) + c(2,i) - 2.d0*divdf1
     c(3,i-1) = 2.d0*(divdf1 - c(2,i-1) - divdf3)/dtau
     c(4,i-1) = (divdf3/dtau)*(6.d0/dtau)
   enddo
END subroutine cubspl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!  Subroutine to use input from Ahlberg spline interpolation
!     for node points rin(1..nin) and function values yin and
!     second derivative values Min to interpolate to radial grid
!     rout with values yout.    For the range 0 \le r \le rin(1)
!     it is assumed that yout=r**(l+1)(W0+W1*r) where l denotes
!     the angular momentum of the function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine specialinterp(nin,rin,yin,MMin,nout,rout,yout,ypout)
  integer, intent(IN) :: nin,nout
  real(dp), intent(IN) :: rin(:),yin(:),MMin(:),rout(:)
  real(dp), intent(INOUT) :: yout(:),ypout(:)
  real(dp), allocatable :: h(:),C(:,:)
  real(dp) :: x
  integer :: i,j
  LOGICAL :: leftin,lefton,leftout,rightin,righton,rightout
  LIBPAW_ALLOCATE(h,(nin))
  LIBPAW_ALLOCATE(C,(4,nin))
  h=0
  do i=2,nin
    h(i)=rin(i)-rin(i-1)
  enddo  
  !h(1) is hopefully not used
  C=0
  C(1,1:nin)=yin(1:nin)
  C(3,1:nin)=MMin(1:nin)
  do i=1,nin-1
    C(4,i)=(MMin(i+1)-MMin(i))/h(i+1)
    C(2,i)=((yin(i+1)-yin(i))/h(i+1)-(MMin(i+1)+2*MMin(i))*h(i+1)/6)
  enddo   
  yout=0;ypout=0
  !  check if grid is within interpolation range
  call checkrange(rin(1),rout(1),leftin,lefton,leftout)
  call checkrange(rout(nout),rin(nin),rightin,righton,rightout)
  if (leftout.or.rightout) then
    write(std_out,*) 'Grid error in specialint',rin(1),rout(1),rin(nin),rout(nout)
    stop
  endif 
  do i=1,nout
    do j=1,nin-1
      call checkrange(rin(j),rout(i),leftin,lefton,leftout)
      call checkrange(rout(i),rin(j+1),rightin,righton,rightout)
      if ((leftin.or.lefton).and.(rightin.or.righton)) then
        x=rout(i)-rin(j)
        yout(i)=c(1,j)+x*(c(2,j)+x*(c(3,j)+x*c(4,j)/3)/2)
        ypout(i)=c(2,j)+x*(c(3,j)+0.5d0*x*c(4,j))
        exit
      endif
    enddo
  enddo
  LIBPAW_DEALLOCATE(h)
  LIBPAW_DEALLOCATE(C)
end subroutine specialinterp  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!    checkrange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
subroutine  checkrange(rin,rout,inrange,onrange,outofbounds)
  LOGICAL, INTENT(OUT) :: inrange,onrange,outofbounds
  REAL(dp), INTENT(IN) :: rin,rout
  REAL(dp), parameter :: tol=1.d-7
  inrange=.false.;onrange=.false.;outofbounds=.false.
  if (rout>=rin) then
    inrange=.true.
    return
  endif
  if (abs(rout-rin).le.tol) then
    onrange=.true.
    return
  endif
  outofbounds=.true.
end subroutine  checkrange





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 19. search_sort
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!!   SUBROUTINE Real_InsSort(A, LUT, Ascending)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
SUBROUTINE Real_InsSort(A, LUT, Ascending)
  REAL(dp),    INTENT(IN)    :: A(:)
  INTEGER, INTENT(INOUT) :: LUT(:)
  LOGICAL, INTENT(IN)    :: Ascending
  INTEGER :: i, j, k, A_Size
  REAL(dp)    :: CurrVal
  A_Size = SIZE(A)
  DO i=1, A_Size
    LUT(i) = i
  END DO
  DO i = 2, A_Size
    CurrVal = A(LUT(i))
    j = i - 1
    DO WHILE ((CurrVal < A(LUT(j))) .AND. (j>1))
      LUT(j+1) = LUT(j)
      j = j - 1
    END DO
    IF (CurrVal < A(LUT(j))) THEN
      LUT(j+1) = LUT(j)
      j = j - 1
    END IF
    LUT(j+1) = i
  END DO
  IF (.NOT. Ascending) THEN
    j = A_Size / 2
    DO i = 1, j
      k = LUT(i)
      LUT(i) = LUT(A_Size - i + 1)
      LUT(A_Size - i + 1) = k
    END DO
  END IF
  RETURN
END SUBROUTINE Real_InsSort





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 20. input_dataset_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME
!!  input_dataset_read
!!
!! FUNCTION
!!  Initialize an input_dataset datastructure by reading it from
!!  a file. If file is omitted, then read from standard input.
!!  Note: we only read here data used to generate the PAW dataset,
!!    not data used for the post-processing (output, explore, scfpaw, ...)
!!
!! INPUTS (all optionals)
!!  [inputfile]= name of input file to be read
!!  [echofile]= name of a file to echo input file content
!!  [read_global_data]= if TRUE, read global data (atom, XC, grid, ...) -
!Default TRUE
!!  [read_elec_data]= if TRUE, read electronic configuration (orbital &
!occupations) - Default TRUE
!!  [read_coreval_data]= if TRUE, read electronic config (core and valence) -
!Default TRUE
!!  [read_basis_data]= if TRUE, read basis data (radii, pseudo scheme, ...) -
!Default TRUE
!!
!! OUTPUT
!!  [input_dt]= datastructure containing the complete input file.
!!              If omitted, then the global public `input_dataset`
!!              is used.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE input_dataset_read(atp,inputfile,echofile,&
&read_global_data,read_elec_data,read_coreval_data,read_basis_data)
!---- Arguments
 CHARACTER*(*),INTENT(IN) :: inputfile
 CHARACTER*(*),INTENT(IN),OPTIONAL :: echofile
 LOGICAL,INTENT(IN),OPTIONAL :: read_global_data,read_elec_data,&
&                               read_coreval_data,read_basis_data
 TYPE(atompaw_type),INTENT(INOUT) :: atp
!---- Local variables
 INTEGER :: ifunit
 INTEGER,PARAMETER :: ecunit=222
 INTEGER,PARAMETER :: nkappa(5)=(/1,2,2,2,2/)
 INTEGER :: input_unit
 INTEGER :: ii,io,nadd,norb,nval,nbl,nn,ik,kk
 INTEGER :: ilin,ilog,inrl,iscl,ipnt,ifin,iend,ihfpp,ilcex,itau
 INTEGER :: igrid,irelat,ilogder,ilogv4,ibd,idirac,ifixz,ll,nstart
 INTEGER :: ispline,isplr0,isplns
 LOGICAL :: has_to_echo
 LOGICAL :: read_global_data_,read_elec_data_,read_coreval_data_,read_basis_data_
 CHARACTER(200) :: inputline,inputword
 !CHARACTER(128) :: exchangecorrelationandgridline
 CHARACTER(256) :: exchangecorrelationandgridline
 CHARACTER(1) :: CHR
 integer,parameter :: XML_RECL=50000
 character (len=XML_RECL) :: line,readline
 logical :: found
 real(dp) :: x1,x2,xocc
 INTEGER :: basis_add_l(nbasis_add_max)
 INTEGER :: basis_add_k(nbasis_add_max)
 real(dp) :: basis_add_energy(nbasis_add_max)
 INTEGER :: tmp_n(norbit_max),tmp_l(norbit_max),tmp_k(norbit_max)
 real(dp) :: tmp_occ(norbit_max)
 ifunit=libpaw_get_free_unit() 
 input_unit=ifunit
 OPEN(ifunit,file=trim(inputfile),form='formatted',action="read")
!Do we echo input file content?
 has_to_echo=PRESENT(echofile)
 IF (has_to_echo) THEN
   OPEN(ecunit,file=trim(echofile),form='formatted')
 END IF
!Select which components have to be read
 read_global_data_=.true.;if (PRESENT(read_global_data))read_global_data_=read_global_data
 read_elec_data_=.true.;if (PRESENT(read_elec_data))read_elec_data_=read_elec_data
 read_coreval_data_=.true.;if (PRESENT(read_coreval_data))read_coreval_data_=read_coreval_data
 read_basis_data_=.true.;if (PRESENT(read_basis_data))read_basis_data_=read_basis_data
!Print a title
 IF(read_global_data_.OR.read_elec_data_.OR.read_coreval_data_.OR.read_basis_data_)THEN
   if(has_to_print) WRITE(STD_OUT,'(/,3x,a)') "===== READING OF INPUT FILE ====="
 END IF

 found=.false.
 do while (.not.found)
   read(input_unit,'(a)',err=10,end=10) readline
   line=adjustl(readline);goto 20
   10 stop
   20 continue
   if (line(1:22)=='<!-- Program:  atompaw') then
     found=.true.
   end if
 enddo
 if(.not.found) then
   LIBPAW_ERROR('XML COREWF FILE NOT CORRECT')
 endif

!------------------------------------------------------------------
!Start reading of AE data
 IF (read_global_data_) THEN
!------------------------------------------------------------------
!=== 1st line: read atomic symbol, atomic number
   READ(input_unit,'(a)') inputline
   IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
   CALL eliminate_comment(inputline)
   READ(inputline,*) atp%atomic_symbol,atp%atomic_charge
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a,a2)') "Atomic symbol : ",atp%atomic_symbol
     WRITE(STD_OUT,'(3x,a,i0)') "Atomic charge : ",atp%atomic_charge
   END IF
   !------------------------------------------------------------------
   !=== 2nd line: read XC type, grid data, relativistic,point-nucleus,
   !              logderiv data, HF data, Block-Davidson keyword 
   !Read full line
   READ(input_unit,'(a)') exchangecorrelationandgridline
   IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(exchangecorrelationandgridline)
   CALL eliminate_comment(exchangecorrelationandgridline)
   CALL Uppercase(exchangecorrelationandgridline)
   exchangecorrelationandgridline=trim(exchangecorrelationandgridline)
   !Retrieve keyword indexes
   ilin=0;ilin=0;ilog=0;ilogv4=0;inrl=0;iscl=0;ipnt=0;ifin=0
   ihfpp=0;ilcex=0;igrid=0;irelat=0;ilogder=0;ibd=0;idirac=0
   ispline=0;isplr0=0;isplns=0
   ilin=INDEX(exchangecorrelationandgridline,'LINEARGRID')
   ilog=INDEX(exchangecorrelationandgridline,'LOGGRID')
   ilogv4=INDEX(exchangecorrelationandgridline,'LOGGRIDV4')
   ibd=INDEX(exchangecorrelationandgridline,'BDSOLVE')
   inrl=INDEX(exchangecorrelationandgridline,'NONRELATIVISTIC')
   iscl=INDEX(exchangecorrelationandgridline,'SCALARRELATIVISTIC')
   idirac=INDEX(exchangecorrelationandgridline,'DIRACRELATIVISTIC')
   ipnt=INDEX(exchangecorrelationandgridline,'POINT-NUCLEUS')
   ifin=INDEX(exchangecorrelationandgridline,'FINITE-NUCLEUS')
   ilogder=INDEX(exchangecorrelationandgridline,'LOGDERIVRANGE')
   ihfpp=INDEX(exchangecorrelationandgridline,'HFPOSTPROCESS')
   ilcex=INDEX(exchangecorrelationandgridline,'LOCALIZEDCOREEXCHANGE')
   ifixz=INDEX(exchangecorrelationandgridline,'FIXED_ZERO')
   itau=INDEX(exchangecorrelationandgridline,'WTAU')
   ispline=INDEX(exchangecorrelationandgridline,'SPLINEINTERP')
   isplr0=INDEX(exchangecorrelationandgridline,'SPLR0')
   isplns=INDEX(exchangecorrelationandgridline,'SPLNS')
   igrid=max(ilin,ilog)  !This line may need attention....
   irelat=max(inrl,iscl) !This line may need attention....
   !!Treat simple logical variables
   atp%scalarrelativistic=(iscl>0.and.inrl==0)
   atp%diracrelativistic=(idirac>0.and.inrl==0)
   atp%usespline=(itau>0.or.ispline>0.and.inrl==0)
   atp%finitenucleus=(ifin>0.and.ipnt==0)
   atp%BDsolve=(ibd>0)
   atp%HFpostprocess=(ihfpp>0)
   !!Treat finite nucleus option
   atp%finitenucleusmodel=-1 
   IF (atp%finitenucleus) THEN
     READ(exchangecorrelationandgridline(ifin+14:ifin+14),'(a)') CHR
     IF (CHR=="2") atp%finitenucleusmodel=2      
     IF (CHR=="3") atp%finitenucleusmodel=3      
     IF (CHR=="4") atp%finitenucleusmodel=4      
     IF (CHR=="5") atp%finitenucleusmodel=5      
   END IF  
   !Treat possible changes to spline grid
   if (isplr0>0) then
     READ(exchangecorrelationandgridline(isplr0+5:),*) atp%splr0
   end if
   if (isplns>0) then
     READ(exchangecorrelationandgridline(isplns+5:),*) atp%splns
   end if
   !!Treat grid data
   atp%gridkey='LINEAR'
   atp%gridpoints=mxgridlin
   atp%gridrange=linrange
   atp%gridmatch=linrange
   IF (ilog>0.and.ilin==0.and.ilogv4==0) THEN
     atp%gridkey='LOGGRID'
     atp%gridpoints=mxgridlog;
     atp%gridrange=logrange
     atp%gridmatch=logrange
   END IF
   IF (ilog>0.and.ilin==0.and.ilogv4>0) THEN
     atp%gridkey='LOGGRID4'
     atp%gridpoints=mxgridlog;
     atp%gridrange=v4logrange
     atp%gridmatch=v4logrange
   END IF
   IF (igrid>0) THEN
     iend=256
     IF (irelat >igrid.and.irelat-1 <iend) iend=irelat -1
     IF (ilogder>igrid.and.ilogder-1<iend) iend=ilogder-1
     IF (ibd>igrid.and.ibd-1<iend) iend=ibd-1
     inputline=""
     IF (ilog>0.and.ilogv4==0.and.iend>igrid+7) &
&      inputline=TRIM(exchangecorrelationandgridline(igrid+7:iend))
     IF (ilog>0.and.ilogv4>0.and.iend>igrid+9) &
&      inputline=TRIM(exchangecorrelationandgridline(igrid+9:iend))
     IF (ilin>0.and.iend>igrid+10) &
&      inputline=TRIM(exchangecorrelationandgridline(igrid+10:iend))
     IF (inputline/="") THEN
       CALL extractword(1,inputline,inputword);inputword=trim(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) atp%gridpoints
         CALL extractword(2,inputline,inputword);inputword=trim(inputword)
         IF (inputword/="") THEN
           READ(inputword,*) atp%gridrange
           atp%gridmatch=atp%gridrange
           CALL extractword(3,inputline,inputword);inputword=trim(inputword)
           IF (inputword/="") read(inputword,*) atp%gridmatch
         END IF
       END IF
     END IF
     IF (atp%gridpoints<=0) STOP "input_dataset: error -- number of grid points should be >0!"
   END IF
   !Treat logderiv data
   atp%minlogderiv=logder_min
   atp%maxlogderiv=logder_max
   atp%nlogderiv=logder_pts
   IF (ilogder>0) THEN
     iend=256
     IF (igrid >ilogder.and.igrid-1 <iend) iend=igrid -1
     IF (irelat>ilogder.and.irelat-1<iend) iend=irelat-1
     inputline=""
     IF (iend>ilogder+13)inputline=trim(exchangecorrelationandgridline(ilogder+13:iend))
     IF (inputline/="") THEN
       CALL extractword(1,inputline,inputword);inputword=trim(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) atp%minlogderiv
         CALL extractword(2,inputline,inputword);inputword=trim(inputword)
         IF (inputword/="") THEN
           READ(inputword,*) atp%maxlogderiv
           CALL extractword(3,inputline,inputword);inputword=trim(inputword)
           IF (inputword/="") READ(inputword,*) atp%nlogderiv
         END IF
       END IF
     END IF
   END IF
   !Treat XC/HF
   if (itau>0) then     
     READ(unit=exchangecorrelationandgridline(itau+5:),fmt=*) atp%exctype
   else  
     READ(unit=exchangecorrelationandgridline(1:),fmt=*) atp%exctype
   endif  
   atp%needvtau=(itau>0.or.TRIM(atp%exctype)=='MGGA-R2SCAN-001'.or.TRIM(atp%exctype)=='MGGA-R2SCAN-01')
   atp%localizedcoreexchange=(ilcex>0)
   atp%fixed_zero=(ifixz>0) ; atp%fixed_zero_index=-1
   IF (atp%fixed_zero) &
  &   READ(unit=exchangecorrelationandgridline(ifixz+10:),fmt=*)atp%fixed_zero_index
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,2a)')     "Scalar-relativistic calculation:",MERGE("YES"," NO",atp%scalarrelativistic)
     WRITE(STD_OUT,'(3x,2a)')     "Dirac-relativistic calculation:",MERGE("YES"," NO",atp%diracrelativistic)
     IF (atp%usespline) THEN
       WRITE(STD_OUT,'(3x,a)')    "    - Use a spline solver"
     END IF
     WRITE(STD_OUT,'(3x,2a)')     "Exchange-correlation functional:",TRIM(atp%exctype)
     WRITE(STD_OUT,'(3x,3a)')     " (mGGA kinetic energy functional:",MERGE("YES"," NO",atp%needvtau),")"
     WRITE(STD_OUT,'(3x,2a)')     "Finite-nucleus calculation:",MERGE("YES"," NO",atp%finitenucleus)
     IF (atp%finitenucleus) THEN
       WRITE(STD_OUT,'(3x,a,i0)') "    - Finite-nucleus model:",atp%finitenucleusmodel
     END IF
     WRITE(STD_OUT,'(3x,2a)')     "Block-Davidson calculation:",MERGE("YES"," NO",atp%BDsolve)
     WRITE(STD_OUT,'(3x,2a)')     "Grid type:",TRIM(atp%gridkey)
     WRITE(STD_OUT,'(3x,a,i0)')   "Grid size:",atp%gridpoints
     WRITE(STD_OUT,'(3x,a,f7.3)') "Grid maximum value:",atp%gridrange
     WRITE(STD_OUT,'(3x,a,f7.3)') "Grid imposed value:",atp%gridmatch
     if(atp%usespline) then
       WRITE(STD_OUT,'(3x,a,f7.3,2x,i0)') "Spline grid r0, ns              :",&
   &      atp%splr0,atp%splns
     endif
     WRITE(STD_OUT,'(3x,a,i0)')   "Log. derivative, number of pts:",atp%nlogderiv
     WRITE(STD_OUT,'(3x,a,f7.3)') "Log. derivative, min. energy:",atp%minlogderiv
     WRITE(STD_OUT,'(3x,a,f7.3)') "Log. derivative, max. energy:",atp%maxlogderiv
     WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, post-processing:",MERGE("YES"," NO",atp%HFpostprocess)
     WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, localized coreex.:",MERGE("YES"," NO",atp%localizedcoreexchange)
     WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, fixed zero:",MERGE("YES"," NO",atp%fixed_zero)
     IF (atp%fixed_zero) THEN
       WRITE(STD_OUT,'(3x,a,i0)') "    - HF fixed zero index:",atp%fixed_zero_index
     END IF
     IF (atp%BDsolve.and.atp%gridkey=='LINEAR') THEN
       WRITE(STD_OUT,'(/,3x,a)') "WARNING: BlockDavidson solver works very slowlywith linear grid!"
     END IF
   END IF
 !------------------------------------------------------------------
 !End reading of global data. Start reading of electronic configuration data
 ENDIF
 IF (read_elec_data_) THEN
   !------------------------------------------------------------------
   !=== 3rd line and following: electronic configuration of atom
   READ(input_unit,'(a)') inputline
   IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
   CALL eliminate_comment(inputline)
   READ(inputline,*) atp%np(1:5)
   DO ll=1,5
     IF(atp%np(ll)<0) atp%np(ll)=0
   END DO
   atp%norbit=atp%np(1)+max(atp%np(2)-1,0)+max(atp%np(3)-2,0) &
&                              +max(atp%np(4)-3,0)+max(atp%np(5)-4,0)
   IF (atp%diracrelativistic) atp%norbit=2*atp%norbit-atp%np(1)
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a,5(1x,i0))') "Max. quantum numbers(s,p,d,f,g):",atp%np(1:5)
     WRITE(STD_OUT,'(3x,a,i0)') "Total number of orbitals: ",atp%norbit
   END IF
   ! CALL input_dataset_read_occ(dataset%norbit_mod,dataset%orbit_mod_l,&
   !&dataset%orbit_mod_n,dataset%orbit_mod_k,dataset%orbit_mod_occ,&
   !&                   dataset%np,dataset%diracrelativistic,&
   !&                   inputfile_unit=input_unit,echofile_unit=ecunit)
   atp%norbit_mod=0
   kk=0
   DO
     READ(input_unit,'(a)') inputline
     IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
     CALL eliminate_comment(inputline)
     if (.not.atp%diracrelativistic) READ(inputline,*) nn,ll,xocc
     if (atp%diracrelativistic) READ(inputline,*) nn,ll,kk,xocc
     IF (nn<=0) EXIT
     IF (xocc<0._dp.OR.&
&      ((.NOT.atp%diracrelativistic).AND.(xocc>2._dp*(2*ll+1))).OR.&
&      ((     atp%diracrelativistic).AND.(xocc>2._dp*ABS(kk)))) THEN
       LIBPAW_ERROR('input_dataset: error in occupations')
     END IF
     atp%norbit_mod=atp%norbit_mod+1
     if (atp%norbit_mod>norbit_max) stop 'input_dataset_occ: error -- to many occupation lines!'
     tmp_l(atp%norbit_mod)=ll
     tmp_n(atp%norbit_mod)=nn
     tmp_k(atp%norbit_mod)=kk
     tmp_occ(atp%norbit_mod)=xocc
   END DO
   IF(ALLOCATED(atp%orbit_mod_l)) then
     LIBPAW_DEALLOCATE(atp%orbit_mod_l)
   endif
   IF(ALLOCATED(atp%orbit_mod_n)) then
     LIBPAW_DEALLOCATE(atp%orbit_mod_n)
   endif
   IF(ALLOCATED(atp%orbit_mod_k)) then
     LIBPAW_DEALLOCATE(atp%orbit_mod_k)
   endif
   IF(ALLOCATED(atp%orbit_mod_occ)) then
     LIBPAW_DEALLOCATE(atp%orbit_mod_occ)
   endif
   LIBPAW_ALLOCATE(atp%orbit_mod_l,(atp%norbit_mod))
   LIBPAW_ALLOCATE(atp%orbit_mod_n,(atp%norbit_mod))
   LIBPAW_ALLOCATE(atp%orbit_mod_k,(atp%norbit_mod))
   LIBPAW_ALLOCATE(atp%orbit_mod_occ,(atp%norbit_mod))
   atp%orbit_mod_l(1:atp%norbit_mod)=tmp_l(1:atp%norbit_mod)
   atp%orbit_mod_n(1:atp%norbit_mod)=tmp_n(1:atp%norbit_mod)
   atp%orbit_mod_k(1:atp%norbit_mod)=tmp_k(1:atp%norbit_mod)
   atp%orbit_mod_occ(1:atp%norbit_mod)=tmp_occ(1:atp%norbit_mod)
 !------------------------------------------------------------------
 !End reading of electronic data. Start reading of core/valence data
 ENDIF
 IF (read_coreval_data_) THEN
 !------------------------------------------------------------------
 !=== Core and valence states
   !Read core and valence states
   IF (ALLOCATED(atp%orbit_iscore)) then
     LIBPAW_DEALLOCATE(atp%orbit_iscore)
   endif
   LIBPAW_ALLOCATE(atp%orbit_iscore,(atp%norbit))
   DO io=1,atp%norbit
     DO
       READ(input_unit,'(a)') inputline
       CALL eliminate_comment(inputline)
       READ(inputline,*) CHR
       IF (CHR=='c'.OR.CHR=='C'.OR.&
&          CHR=='v'.OR.CHR=='V') THEN
         IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
         EXIT
       ELSE
         LIBPAW_ERROR('Please input c or v!')
       END IF
     END DO
     atp%orbit_iscore(io)=(CHR=='c'.OR.CHR=='C')
   END DO
   !Store valence states
   atp%norbit_val=atp%norbit-COUNT(atp%orbit_iscore(:))
   IF (ALLOCATED(atp%orbit_val_n)) then
     LIBPAW_DEALLOCATE(atp%orbit_val_n)
   endif
   IF (ALLOCATED(atp%orbit_val_l)) then
     LIBPAW_DEALLOCATE(atp%orbit_val_l)
   endif
   IF (ALLOCATED(atp%orbit_val_k)) then
     LIBPAW_DEALLOCATE(atp%orbit_val_k)
   endif
   LIBPAW_ALLOCATE(atp%orbit_val_n,(atp%norbit_val))
   LIBPAW_ALLOCATE(atp%orbit_val_l,(atp%norbit_val))
   LIBPAW_ALLOCATE(atp%orbit_val_k,(atp%norbit_val))
   kk=0
   io=0;nval=0
   DO ll=0,4
     nn=atp%np(ll+1)
     IF (nn>0) THEN
       DO ik=1,MERGE(nkappa(ll+1),1,atp%diracrelativistic)
         kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
         IF (.NOT.atp%diracrelativistic) kk=0
         DO ii=1+ll,nn
           io=io+1
           IF (.NOT.atp%orbit_iscore(io)) THEN
             nval=nval+1
             atp%orbit_val_n(nval)=ii
             atp%orbit_val_l(nval)=ll
             atp%orbit_val_k(nval)=kk    
           END IF  
         END DO
       END DO
     END IF
   END DO
   IF (atp%norbit_val/=nval) STOP 'input_dataset: bug -- wrong nval!'
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a)') "Core and valence orbitals:"
     IF (.NOT.atp%diracrelativistic) WRITE(STD_OUT,'(7x,a)') "n l : type"
     IF (atp%diracrelativistic)      WRITE(STD_OUT,'(7x,a)') "n l kappa :type"
     io=0
     DO ll=0,4
       nn=atp%np(ll+1)
       IF (nn>0) THEN
         IF (.NOT.atp%diracrelativistic) THEN
           DO ii=1+ll,nn
             io=io+1
             WRITE(STD_OUT,'(7x,i1,1x,i1,2a)') ii,ll," : ", &
   &            MERGE("CORE   ","VALENCE",atp%orbit_iscore(io))
           END DO
         ELSE
           DO ik=1,nkappa(ll+1) 
             kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
             DO ii=1+ll,nn
               io=io+1
               WRITE(STD_OUT,'(7x,i1,1x,i1,2x,i2,2x,2a)') ii,ll,kk," : ", &
   &              MERGE("CORE   ","VALENCE",atp%orbit_iscore(io))
             END DO
           END DO
         END IF
       END IF
     END DO
   END IF
 !------------------------------------------------------------------
 !End reading of AE data. Start reading of basis data
 ENDIF
 IF (read_basis_data_) THEN
 !------------------------------------------------------------------
 !=== Maximum L for basis functions
   READ(input_unit,'(a)') inputline
   IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
   CALL eliminate_comment(inputline)
   READ(inputline,*) atp%lmax
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a,i0)') "Basis, maximum L : ",atp%lmax
   END IF
   !------------------------------------------------------------------
   !=== Cut-off radii
   READ(input_unit,'(a)') inputline
   IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
   CALL eliminate_comment(inputline)
   CALL extractword(1,inputline,inputword);inputword=trim(inputword)
   IF (inputword/="") READ(inputword,*) atp%rc
   IF (atp%rc<=tol12) THEN
     LIBPAW_ERROR('input_dataset: error -- rc too small ')
   END IF
   CALL extractword(2,inputline,inputword);inputword=trim(inputword)
   IF (inputword/="") THEN
     READ(inputword,*) atp%rc_shap
     CALL extractword(3,inputline,inputword);inputword=trim(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) atp%rc_vloc
       CALL extractword(4,inputline,inputword);inputword=trim(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) atp%rc_core
       ELSE
         LIBPAW_ERROR('input_dataset: error -- rc(core) is missing!')
       END IF
     ELSE
       LIBPAW_ERROR('input_dataset: error -- rc(Vloc) is missing!')
     END IF
     IF (atp%rc_shap<=tol12.OR.atp%rc_vloc<=tol12.OR.&
&        atp%rc_core<=tol12) THEN
       LIBPAW_ERROR('input_dataset: error -- one rc is too small!')
     END IF
     IF (atp%rc_shap>atp%rc.OR.atp%rc_vloc>atp%rc.OR.&
&        atp%rc_core>atp%rc) THEN
       LIBPAW_ERROR('input_dataset: error -- rc_shape, rc_vloc and rc_core must be <rc!')
     END IF
   ENDIF
   IF(atp%rc_shap==zero) atp%rc_shap=atp%rc
   IF(atp%rc_vloc==zero) atp%rc_vloc=atp%rc
   IF(atp%rc_core==zero) atp%rc_core=atp%rc
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a,f7.4)') "Augmentation region radius : ",atp%rc
     WRITE(STD_OUT,'(3x,a,f7.4)') "Core dens. matching radius : ",atp%rc_core
     WRITE(STD_OUT,'(3x,a,f7.4)') "Local pot. matching radius : ",atp%rc_vloc
     WRITE(STD_OUT,'(3x,a,f7.4)') "Compens. shape func radius : ",atp%rc_shap
   END IF
   !------------------------------------------------------------------
   !=== Additional basis functions
   nstart=0 ; atp%nbasis_add=0 ; basis_add_k(:)=0
   DO ll=0,atp%lmax
     nbl=0
     nadd = MERGE(nkappa(ll+1),1,atp%diracrelativistic)
     IF (atp%np(ll+1)>0) THEN
       nbl=COUNT(.NOT.atp%orbit_iscore(nstart+1:nstart+atp%np(ll+1)-ll))
       nstart=nstart+atp%np(ll+1)-ll
     END IF
     DO
       READ(input_unit,'(a)') inputline
       IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
       CALL eliminate_comment(inputline)
       READ(inputline,*) CHR
       IF (CHR/='y'.AND.CHR/='Y') THEN
         IF (CHR/='n'.AND.CHR/='N') STOP 'input_dataset: error -- Please enter Y or N!'
         EXIT
       END IF
       atp%nbasis_add=atp%nbasis_add+nadd
       IF (atp%nbasis_add>nbasis_add_max) STOP 'Too many additional basis functions!'
       basis_add_l(atp%nbasis_add-nadd+1:atp%nbasis_add)=ll
       IF (atp%diracrelativistic) THEN
         basis_add_k(atp%nbasis_add)=-1
         IF (ll/=0) THEN
           basis_add_k(atp%nbasis_add-1)=ll
           basis_add_k(atp%nbasis_add)=-(ll+1)
         END IF
       END IF
       READ(input_unit,'(a)') inputline
       IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
       CALL eliminate_comment(inputline)
       READ(inputline,*) basis_add_energy(atp%nbasis_add-nadd+1:atp%nbasis_add)
     END DO
   END DO
   IF (ALLOCATED(atp%basis_add_l)) then
     LIBPAW_DEALLOCATE(atp%basis_add_l)
   endif
   IF (ALLOCATED(atp%basis_add_k)) then
     LIBPAW_DEALLOCATE(atp%basis_add_k)
   endif
   IF (ALLOCATED(atp%basis_add_energy)) then
     LIBPAW_DEALLOCATE(atp%basis_add_energy)
   endif
   LIBPAW_ALLOCATE(atp%basis_add_l,(atp%nbasis_add))
   LIBPAW_ALLOCATE(atp%basis_add_k,(atp%nbasis_add))
   LIBPAW_ALLOCATE(atp%basis_add_energy,(atp%nbasis_add))
   IF (atp%nbasis_add>0) THEN
     atp%basis_add_l(1:atp%nbasis_add)=basis_add_l(1:atp%nbasis_add)
     atp%basis_add_k(1:atp%nbasis_add)=basis_add_k(1:atp%nbasis_add)
     atp%basis_add_energy(1:atp%nbasis_add)=basis_add_energy(1:atp%nbasis_add)
   END IF
   atp%nbasis=COUNT(.NOT.atp%orbit_iscore(:))+atp%nbasis_add
   !Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a,i0)') "Initial number of basis functions:",atp%nbasis-atp%nbasis_add
     WRITE(STD_OUT,'(3x,a,i0)') "Number of additional basis functions:",atp%nbasis_add
     WRITE(STD_OUT,'(3x,a,i0)') "Total number of basis functions:",atp%nbasis
     WRITE(STD_OUT,'(3x,a)') "Additional basis functions:"
     IF (.NOT.atp%diracrelativistic) THEN
       WRITE(STD_OUT,'(7x,a)') "l : energy"
       DO io=1,atp%nbasis_add
         WRITE(STD_OUT,'(7x,i1,a,f7.4)') atp%basis_add_l(io)," :",atp%basis_add_energy(io)
       END DO
     ELSE
       WRITE(STD_OUT,'(7x,a)') "l kappa : energy"
       DO io=1,atp%nbasis_add
         WRITE(STD_OUT,'(7x,i1,2x,i2,2x,a,f7.4)') atp%basis_add_l(io), &
   &          atp%basis_add_k(io)," : " ,atp%basis_add_energy(io)
       END DO
     END IF
   END IF
   !------------------------------------------------------------------
   !=== Projectors, compensation charge shape function, core tolerance
   READ(input_unit,'(a)') inputline
   IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
   CALL eliminate_comment(inputline)
   CALL Uppercase(inputline)
   inputline=TRIM(inputline)
   atp%pseudo_type=PSEUDO_TYPE_BLOECHL
   atp%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
   atp%pseudo_polynom2_pdeg=polynom2_pdeg_def
   atp%pseudo_polynom2_qcut=polynom2_qcut_def
   atp%shapefunc_type=SHAPEFUNC_TYPE_SINC
   atp%shapefunc_gaussian_param=gausstol_def
   atp%hf_coretol=hf_coretol_def
   READ(unit=inputline,fmt=*) inputword
   IF (TRIM(inputword)=='BLOECHL'.OR.TRIM(inputword)=='VNCT') THEN
     atp%projector_type=PROJECTOR_TYPE_BLOECHL
     atp%pseudo_type=PSEUDO_TYPE_BLOECHL
     atp%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
   ELSE IF (TRIM(inputword)=='VNCK') THEN
     atp%projector_type=PROJECTOR_TYPE_BLOECHL
     atp%pseudo_type=PSEUDO_TYPE_BLOECHL_K
     atp%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
   ELSE IF (TRIM(inputword)=='VANDERBILT'.OR.TRIM(inputword)=='VNCTV') THEN
     atp%projector_type=PROJECTOR_TYPE_VANDERBILT
     atp%pseudo_type=PSEUDO_TYPE_POLYNOM
     atp%ortho_type=ORTHO_TYPE_VANDERBILT
   ELSE IF(TRIM(inputword)=='MODRRKJ') THEN
     atp%projector_type=PROJECTOR_TYPE_MODRRKJ
     atp%pseudo_type=PSEUDO_TYPE_RRKJ
     atp%ortho_type=ORTHO_TYPE_VANDERBILT
     IF (INDEX(inputline,'VANDERBILTORTHO')>0)atp%ortho_type=ORTHO_TYPE_VANDERBILT
     IF (INDEX(inputline,'GRAMSCHMIDTORTHO')>0)atp%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
     IF (INDEX(inputline,'SVDORTHO')>0) atp%ortho_type=ORTHO_TYPE_SVD
   ELSE IF (TRIM(inputword)=='CUSTOM') THEN
     atp%projector_type=PROJECTOR_TYPE_CUSTOM
     IF (INDEX(inputline,'BLOECHLPS')>0) THEN
       atp%pseudo_type=PSEUDO_TYPE_BLOECHL
       atp%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
     ELSE IF (INDEX(inputline,'POLYNOM2')>0) THEN
       atp%pseudo_type=PSEUDO_TYPE_POLYNOM2
       nstart=INDEX(inputline,'POLYNOM2')
       READ(unit=inputline(nstart+8:),fmt=*,err=111,end=111,iostat=io) &
&           atp%pseudo_polynom2_pdeg,atp%pseudo_polynom2_qcut
111  CONTINUE
   ELSE IF (INDEX(inputline,'POLYNOM')>0) THEN
     atp%pseudo_type=PSEUDO_TYPE_POLYNOM
   ELSE IF (INDEX(inputline,'RRKJ')>0) THEN
     atp%pseudo_type=PSEUDO_TYPE_RRKJ
   END IF
   IF (INDEX(inputline,'VANDERBILTORTHO')>0)atp%ortho_type=ORTHO_TYPE_VANDERBILT
   IF (INDEX(inputline,'GRAMSCHMIDTORTHO')>0)atp%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
 END IF
 IF (TRIM(atp%exctype)=='HF') THEN
   atp%projector_type=PROJECTOR_TYPE_HF
   atp%pseudo_type=PSEUDO_TYPE_HF
   atp%ortho_type=ORTHO_TYPE_HF
   if(has_to_print) WRITE(STD_OUT,'(3x,a)') '>> You are using HF XC type: pseudo and orthogonalization line will be ignored!'
 END IF
 IF ((atp%pseudo_type==PSEUDO_TYPE_BLOECHL.OR. &
&     atp%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&   .AND.atp%ortho_type==ORTHO_TYPE_VANDERBILT) STOP &
&  'input_dataset: error -- Vanderbilt orthogonalization not compatible with Bloechls projector scheme!'
 IF ((atp%pseudo_type==PSEUDO_TYPE_BLOECHL.OR. &
&     atp%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&   .AND.atp%ortho_type==ORTHO_TYPE_VANDERBILT) STOP &
&  'input_dataset: error -- Vanderbilt orthogonalization not compatible with Bloechls projector scheme!'
 IF ((atp%projector_type==PROJECTOR_TYPE_BLOECHL) &
&   .AND.atp%needvtau) STOP &
&   'input_dataset: error -- mGGA not compatible the Bloechl projector scheme!'
 !!!! Hopefully this will never happen             
 IF ((atp%projector_type==PROJECTOR_TYPE_HF) &
&   .AND.atp%needvtau) STOP &
&   'input_dataset: error -- mGGA and Hartree-Fock are not compatible!'
 IF ((atp%pseudo_type==PSEUDO_TYPE_BLOECHL.OR. &
&     atp%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&   .AND.atp%needvtau) STOP &
&   'input_dataset: error -- mGGA not compatible the Bloechl pseudization scheme!'
 IF (INDEX(inputline,'SINC2')>0) THEN
   atp%shapefunc_type=SHAPEFUNC_TYPE_SINC
 ELSE IF (INDEX(inputline,'GAUSSIAN')>0) THEN
   atp%shapefunc_type=SHAPEFUNC_TYPE_GAUSSIAN
   nstart=INDEX(inputline,'GAUSSIAN')
   READ(unit=inputline(nstart+8:),fmt=*,err=222,end=222,iostat=io) &
&       atp%shapefunc_gaussian_param
222 CONTINUE
 ELSE IF (INDEX(inputline,'BESSELSHAPE')>0) THEN
   atp%shapefunc_type=SHAPEFUNC_TYPE_BESSEL
 END IF
 nstart=INDEX(inputline,'CORETOL')
 IF (nstart>0) THEN
   READ(unit=inputline(nstart+7:),fmt=*) atp%hf_coretol
 END IF
 atp%shapetcore=(INDEX(inputline,'SHAPETCORE')>0)
 !Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a)') "Projectors description:"
   IF (atp%projector_type==PROJECTOR_TYPE_BLOECHL) &
 &    WRITE(STD_OUT,'(7x,a)') "Type              : BLOECHL"
   IF (atp%projector_type==PROJECTOR_TYPE_VANDERBILT) &
 &    WRITE(STD_OUT,'(7x,a)') "Type              : VANDERBILT"
   IF (atp%projector_type==PROJECTOR_TYPE_MODRRKJ) &
 &    WRITE(STD_OUT,'(7x,a)') "Type              : MODRRKJ"
   IF (atp%projector_type==PROJECTOR_TYPE_CUSTOM) &
 &    WRITE(STD_OUT,'(7x,a)') "Type              : CUSTOM"
   IF (atp%projector_type==PROJECTOR_TYPE_HF) &
 &    WRITE(STD_OUT,'(7x,a)') "Type : HARTREE-FOCK"
   IF (atp%projector_type/=PROJECTOR_TYPE_HF) THEN
     IF (atp%pseudo_type==PSEUDO_TYPE_BLOECHL) &
 &      WRITE(STD_OUT,'(7x,a)') "Pseudization      : BLOECHL"
     IF (atp%pseudo_type==PSEUDO_TYPE_POLYNOM) &
 &      WRITE(STD_OUT,'(7x,a)') "Pseudization      : POLYNOM"
     IF (atp%pseudo_type==PSEUDO_TYPE_RRKJ) &
 &      WRITE(STD_OUT,'(7x,a)') "Pseudization      : RRKJ"
     IF (atp%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
 &      WRITE(STD_OUT,'(7x,a)') "Pseudization      : BLOECHL KERKER"
     IF (atp%pseudo_type==PSEUDO_TYPE_POLYNOM2) &
 &      WRITE(STD_OUT,'(7x,a,i0,a,es9.3)') "Pseudization      : POLYNOM2,pdeg=",&
 &       atp%pseudo_polynom2_pdeg,", qcut=",atp%pseudo_polynom2_qcut
     IF (atp%ortho_type==ORTHO_TYPE_GRAMSCHMIDT) &
 &      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : GRAM-SCHMIDT"
     IF (atp%ortho_type==ORTHO_TYPE_VANDERBILT) &
 &      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : VANDERBILT"
     IF (atp%ortho_type==ORTHO_TYPE_SVD) &
 &      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : SVD"
   END IF
   IF (atp%shapefunc_type==SHAPEFUNC_TYPE_GAUSSIAN) &
 &    WRITE(STD_OUT,'(3x,a,es9.3)') "Compensation charge shape function : GAUSSIAN, tol=",&
 &    atp%shapefunc_gaussian_param
   IF (atp%shapefunc_type==SHAPEFUNC_TYPE_SINC) &
 &    WRITE(STD_OUT,'(3x,a)') "Compensation charge shape function : SINC2"
   IF (atp%shapefunc_type==SHAPEFUNC_TYPE_BESSEL) &
 &    WRITE(STD_OUT,'(3x,a)') "Compensation charge shape function : BESSEL"
   IF (INDEX(inputline,'CORETOL')>0) &
 &    WRITE(STD_OUT,'(3x,a,es9.3)') "Core tolerance for Hartree-Fock:",atp%hf_coretol
     WRITE(STD_OUT,'(3x,2a)') "Smooth tcore shape (no negative nhat):",MERGE("YES"," NO",atp%shapetcore)
 END IF
 !------------------------------------------------------------------
 !=== Local pseudopotential
 READ(input_unit,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)
 call Uppercase(inputline)
 inputline=TRIM(inputline)
 atp%vloc_type=VLOC_TYPE_MTROULLIER
 atp%vloc_l=-1
 atp%vloc_ene=0._dp
 atp%vloc_setvloc_coef=0._dp
 atp%vloc_setvloc_rad=atp%rc
 atp%vloc_kerker_power(:)=0
 IF (INDEX(inputline,'MTROULLIER')>0) THEN
   atp%vloc_type=VLOC_TYPE_MTROULLIER
 ELSE IF (INDEX(inputline,'ULTRASOFT')>0) THEN
   atp%vloc_type=VLOC_TYPE_ULTRASOFT
 ELSE IF (INDEX(inputline,'BESSEL')>0) THEN
   atp%vloc_type=VLOC_TYPE_BESSEL
 ELSE IF (INDEX(inputline,'VPSMATCHNC')>0) THEN
   atp%vloc_type=VLOC_TYPE_VPSMATCHNC
 ELSE IF (INDEX(inputline,'VPSMATCHNNC')>0) THEN
   atp%vloc_type=VLOC_TYPE_VPSMATCHNNC
 ELSE IF (INDEX(inputline,'SETVLOC')>0) THEN
   atp%vloc_type=VLOC_TYPE_SETVLOC
   nstart=INDEX(inputline,'SETVLOC')
   READ(unit=inputline(nstart+8:),fmt=*,err=333,end=333,iostat=io) x1,x2
   IF (x1<10._dp**3.AND.x1>-10._dp**3) atp%vloc_setvloc_coef=x1
   IF (x2>tol8.AND.x2<atp%rc) atp%vloc_setvloc_rad=x2
333  CONTINUE
 ELSE IF (INDEX(inputline,'KERKER')>0.OR.atp%pseudo_type==PSEUDO_TYPE_BLOECHL_K) THEN
   IF (INDEX(inputline,'EXPF')>0) THEN
     atp%vloc_type=VLOC_TYPE_KERKER_EXPF
     nstart=INDEX(inputline,'EXPF')
   ELSE IF (INDEX(inputline,'POLY')>0) THEN
     atp%vloc_type=VLOC_TYPE_KERKER_POLY
     nstart=INDEX(inputline,'POLY')
   ELSE
     STOP "EXPF or POLY keyword missing!"
   END IF   
   READ(unit=inputline(nstart+5:),fmt=*,err=334,end=334,iostat=io) &
&    atp%vloc_kerker_power(1:4)
334  CONTINUE
 END IF
 IF ((atp%vloc_type==VLOC_TYPE_SETVLOC.OR. &
&     atp%vloc_type==VLOC_TYPE_KERKER_EXPF.OR. &
&     atp%vloc_type==VLOC_TYPE_KERKER_POLY) &
&   .AND.atp%needvtau) STOP &
&   'input_dataset: error -- mGGA not compatible the chosen Vloc scheme!'
 IF (atp%vloc_type==VLOC_TYPE_MTROULLIER.OR. &
&    atp%vloc_type==VLOC_TYPE_VPSMATCHNC.OR. &
&    atp%vloc_type==VLOC_TYPE_VPSMATCHNNC.OR. &
&    atp%vloc_type==VLOC_TYPE_ULTRASOFT) THEN
   READ(unit=inputline,fmt=*,err=444,end=444,iostat=io) atp%vloc_l,atp%vloc_ene
444  CONTINUE
   IF (atp%vloc_l<0.or.atp%vloc_l>10) STOP 'input_dataset: error while reading Vloc parameters!'
 END IF
 IF (atp%vloc_type==VLOC_TYPE_MTROULLIER.AND.atp%needvtau) then
   if(has_to_print) WRITE(STD_OUT,'(7x,a)') 'NOTE: MTROULLIER Vloc not available for mGGA!'
   if(has_to_print) WRITE(STD_OUT,'(7x,a)') '      Calling VPSmatch with norm conservation instead.'
   atp%vloc_type=VLOC_TYPE_VPSMATCHNC
 ENDIF
 !Print read data
 IF (has_to_print) THEN
   IF (atp%vloc_type==VLOC_TYPE_MTROULLIER) &
 &    WRITE(STD_OUT,'(7x,a,i0,a,f7.4)') "Local pseudopotential type : MTROULLIER,l=",&
 &          atp%vloc_l,", energy=",atp%vloc_ene
   IF (atp%vloc_type==VLOC_TYPE_ULTRASOFT) &
 &    WRITE(STD_OUT,'(7x,a,i0,a,f7.4)') "Local pseudopotential type : ULTRASOFT,l=",&
 &          atp%vloc_l,", energy=",atp%vloc_ene
   IF (atp%vloc_type==VLOC_TYPE_BESSEL) &
 &    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : BESSEL"
   IF (atp%vloc_type==VLOC_TYPE_VPSMATCHNC) &
 &    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNC"
   IF (atp%vloc_type==VLOC_TYPE_VPSMATCHNNC) &
 &    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNNC"
   IF (atp%vloc_type==VLOC_TYPE_SETVLOC) THEN
     WRITE(STD_OUT,'(7x,a,es9.4,a,es9.4)') "Local pseudopotential type :SETVLOC, coef=",&
 &          atp%vloc_setvloc_coef,", rad=",atp%vloc_setvloc_rad
     IF (atp%needvtau) THEN
       LIBPAW_ERROR('SETVLOC  option not available for MGGA')
     ENDIF     
   ENDIF  
   IF (atp%vloc_type==VLOC_TYPE_KERKER_EXPF) &
 &    WRITE(STD_OUT,'(7x,a,4(1x,i0))') "Local pseudopotential type : KERKER EXPF,powers=",&
 &          atp%vloc_kerker_power(1:4)
   IF (atp%vloc_type==VLOC_TYPE_KERKER_POLY) &
 &    WRITE(STD_OUT,'(7x,a,4(1x,i0))') "Local pseudopotential type : KERKER POLY,powers=",&
 &          atp%vloc_kerker_power(1:4)
   IF (atp%vloc_type==VLOC_TYPE_MTROULLIER.AND.atp%needvtau) THEN
     WRITE(STD_OUT,'(7x,a)') 'NOTE: MTROULLIER Vloc not available for mGGA!'
     WRITE(STD_OUT,'(7x,a)') '      Calling VPSmatch with norm conservation instead.'
     atp%vloc_type=VLOC_TYPE_VPSMATCHNC
     WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNC"
   END IF
 END IF
 !------------------------------------------------------------------
 !=== Matching radii for the basis functions
 !Not for all choice of projectors
 IF (atp%projector_type==PROJECTOR_TYPE_CUSTOM.OR.&
&    atp%projector_type==PROJECTOR_TYPE_VANDERBILT.OR.&
&    atp%projector_type==PROJECTOR_TYPE_MODRRKJ.OR.&
& atp%projector_type==PROJECTOR_TYPE_HF.AND.atp%vloc_type==VLOC_TYPE_MTROULLIER)THEN
   IF (ALLOCATED(atp%basis_func_rc)) then
     LIBPAW_DEALLOCATE(atp%basis_func_rc)
   endif
   LIBPAW_ALLOCATE(atp%basis_func_rc,(atp%nbasis))
   norb=0
   DO ll=0,atp%lmax
     DO ik=1,MERGE(nkappa(ll+1),1,atp%diracrelativistic)
       kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
       IF (.NOT.atp%diracrelativistic) kk=0
       DO io=1,atp%norbit_val
         IF (atp%orbit_val_l(io)==ll.AND. &
&           ((.NOT.atp%diracrelativistic).OR.atp%orbit_val_k(io)==kk)) THEN     
           norb=norb+1     
           READ(input_unit,'(a)') inputline
           IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
           CALL eliminate_comment(inputline)
           READ(inputline,*) atp%basis_func_rc(norb)
         END IF
       END DO
       IF (atp%nbasis_add>0) THEN
         DO io=1,atp%nbasis_add
           IF (atp%basis_add_l(io)==ll.AND. &
&             ((.NOT.atp%diracrelativistic).OR.atp%basis_add_k(io)==kk)) THEN     
             norb=norb+1     
             READ(input_unit,'(a)') inputline
             IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
             CALL eliminate_comment(inputline)
             READ(inputline,*) atp%basis_func_rc(norb)
           END IF
         END DO
       END IF
     END DO
   END DO
   IF (atp%nbasis/=norb) STOP 'input_dataset: error -- inconsistency in the number of basis functions!'
   !  Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a)') "Matching radius for basis functions:"
     IF (.NOT.atp%diracrelativistic) WRITE(STD_OUT,'(7x,a)') " # - n l : radius"
     IF (atp%diracrelativistic) WRITE(STD_OUT,'(7x,a)') " # - n l kappa : radius"
     norb=0
     DO ll=0,atp%lmax
       DO ik=1,MERGE(nkappa(ll+1),1,atp%diracrelativistic)
         kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
         IF (.NOT.atp%diracrelativistic) kk=0
         DO io=1,atp%norbit_val
           IF (atp%orbit_val_l(io)==ll.AND. &
   &          ((.NOT.atp%diracrelativistic).OR.atp%orbit_val_k(io)==kk))THEN     
             norb=norb+1
           IF (.NOT.atp%diracrelativistic) &
   &           WRITE(STD_OUT,'(7x,i2,a,i1,1x,i1,a,f7.4)') &
   &         norb," - ",atp%orbit_val_n(io),ll," :",atp%basis_func_rc(norb)
           IF (atp%diracrelativistic) &
   &          WRITE(STD_OUT,'(7x,i2,a,i1,1x,i1,2x,i2,2x,a,f7.4)') &
   &          norb," - ",atp%orbit_val_n(io),ll,kk," :",atp%basis_func_rc(norb)
           END IF
         END DO
         IF (atp%nbasis_add>0) THEN
           DO io=1,atp%nbasis_add
             IF (atp%basis_add_l(io)==ll.AND. &
   &          ((.NOT.atp%diracrelativistic).OR.atp%basis_add_k(io)==kk))THEN     
               norb=norb+1
               IF (.NOT.atp%diracrelativistic) &
   &             WRITE(STD_OUT,'(7x,i2,a,a1,1x,i1,a,f7.4)') &
   &             norb," - ",".",ll," : ",atp%basis_func_rc(norb)
               IF (atp%diracrelativistic) &
   &             WRITE(STD_OUT,'(7x,i2,a,a1,1x,i1,2x,i2,2x,a,f7.4)') &
   &             norb," - ",".",ll,kk," : ",atp%basis_func_rc(norb)
             END IF
           END DO
         END IF
       END DO
     END DO
   END IF
 ELSE ! Other projectors
   IF (ALLOCATED(atp%basis_func_rc)) then
     LIBPAW_DEALLOCATE(atp%basis_func_rc)
   endif
   LIBPAW_ALLOCATE(atp%basis_func_rc,(0))
 END IF
 !------------------------------------------------------------------
 !End reading of basis data
 ENDIF
 !Final message
 IF(read_global_data_.OR.read_elec_data_.OR.read_coreval_data_.OR.read_basis_data_)THEN
   if(has_to_print) WRITE(STD_OUT,'(3x,a)') "===== END READING OF INPUT FILE ====="
 END IF
 if(has_to_print) WRITE(STD_OUT,'(2/)')
 !------------------------------------------------------------------
 !Close files
 CLOSE(ifunit)
 IF (has_to_echo) THEN
   CLOSE(ecunit)
 END IF
END SUBROUTINE input_dataset_read


end module m_paw_atom_solve
!!***
