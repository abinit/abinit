!!****m* ABINIT/m_pawtab
!! NAME
!!  m_pawtab
!!
!! FUNCTION
!!  This module contains the definition of the pawtab_type structured datatype,
!!  as well as related functions and methods.
!!  pawtab_type variables define TABulated data for PAW (from pseudopotential)
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_pawtab

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_pawtab/wvlpaw_rholoc_type
!! NAME
!! wvlpaw_rholoc_type
!!
!! FUNCTION
!! Objects for WVL+PAW
!!
!! SOURCE

 type,public :: wvlpaw_rholoc_type

  integer :: msz
! mesh size

  real(dp),allocatable :: d(:,:)
! local rho and derivatives

  real(dp),allocatable :: rad(:)
! radial mesh

 end type wvlpaw_rholoc_type

 public :: wvlpaw_rholoc_free    ! Free memory
 public :: wvlpaw_rholoc_nullify ! Nullify content
!!***

!----------------------------------------------------------------------

!!****t* m_pawtab/wvlpaw_type
!! NAME
!! wvlpaw_type
!!
!! FUNCTION
!! Objects for WVL+PAW
!!
!! SOURCE

 type,public :: wvlpaw_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: npspcode_init_guess
   ! This is for the PAW-WVL case, only for the initial guess

  integer :: ptotgau
   ! total number of complex gaussians
   ! for tproj

  integer,allocatable :: pngau(:)
   ! number of complex gaussians per basis element
   ! for tproj

!Real pointers

  real(dp),allocatable :: parg(:,:)
   !argument of Gaussians

  real(dp),allocatable :: pfac(:,:)
   !factors of Gaussians

!Other scalars

  type(wvlpaw_rholoc_type) :: rholoc
   ! local density
   !   d(:,1): local rho
   !   d(:,2): local rho 2nd-derivative
   !   d(:,3): local pot
   !   d(:,4): local pot 2nd-derivative

 end type wvlpaw_type

 public :: wvlpaw_allocate  ! Allocate memory
 public :: wvlpaw_free   ! Free memory
 public :: wvlpaw_nullify
!!***

!----------------------------------------------------------------------

!!****t* m_pawtab/pawtab_type
!! NAME
!! pawtab_type
!!
!! FUNCTION
!! This structured datatype contains TABulated data for PAW (from pseudopotential)
!! used in PAW calculations.
!!
!! SOURCE

 type,public :: pawtab_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: basis_size
   ! Number of elements for the paw nl basis on the considered atom type

  integer :: has_coretau
   ! Flag controling use of core kinetic enrgy density (AE and pseudo)
   ! if 1, [t]coretau() is allocated.
   ! if 2, [t]coretau() is computed and stored.

  integer :: has_fock
   ! if 1, onsite matrix elements of the core-valence Fock operator are allocated
   ! if 2, onsite matrix elements of the core-valence Fock operator are computed and stored

  integer :: has_kij
   ! if 1, onsite matrix elements of the kinetic operator are allocated
   ! if 2, onsite matrix elements of the kinetic operator are computed and stored

  integer :: has_shapefncg
   ! if 1, the spherical Fourier transforms of the radial shape functions are allocated
   ! if 2, the spherical Fourier transforms of the radial shape functions are computed and stored

  integer :: has_nabla
   ! if 1, onsite matrix elements of the nabla operator are allocated.
   ! if 2, onsite matrix elements of the nabla operator are computed and stored.

  integer :: has_nablaphi
   ! if 1, nablaphi are allocated
   ! if 2, nablaphi are computed for MetaGGA and stored.

  integer :: has_tproj
   ! Flag controling use of projectors in real space (0 if tnval is unknown)
   ! if 1, tproj() is allocated.
   ! if 2, tproj() is computed and stored.

  integer :: has_tvale
   ! Flag controling use of pseudized valence density (0 if tnval is unknown)
   ! if 1, tvalespl() is allocated.
   ! if 2, tvalespl() is computed and stored.

  integer :: has_vhtnzc
   ! if 1, space for vhtnzc is allocated
   ! if 2, vhtnzc has been read from PAW file and stored

  integer :: has_vhnzc
   ! if 1, space for vhnzc is allocated
   ! if 2, vhnzc has been computed and stored

  integer :: has_vminushalf
   ! has_vminushalf=0 ; vminushal is not allocated
   ! has_vminushalf=1 ; vminushal is not allocated and stored

  integer :: has_wvl
   ! if 1, data for wavelets (pawwvl) are allocated
   ! if 2, data for wavelets (pawwvl) are computed and stored

  integer :: ij_proj
   ! Number of (i,j) elements for the orbitals on which U acts (PAW+U only)
   ! on the considered atom type (ij_proj=1 (1 projector), 3 (2 projectors)...)
   ! Also used for local exact-exchange

  integer :: ij_size
   ! Number of (i,j) elements for the symetric paw basis
   ! on the considered atom type (ij_size=basis_size*(basis_size+1)/2)

  integer :: lcut_size
   ! Maximum value of l+1 leading to non zero Gaunt coeffs
   ! modified by dtset%pawlcutd
   ! lcut_size=min(2*l_max,dtset%pawlcutd)+1

  integer :: l_size
   ! Maximum value of l+1 leading to non zero Gaunt coeffs
   ! l_size=2*l_max-1

  integer :: lexexch
   ! lexexch gives l on which local exact-exchange is applied for a given type of atom.

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz
   ! lmnmix_sz=number of klmn=(lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix

  integer :: lpawu
   ! lpawu gives l on which U is applied for a given type of atom.

  integer :: nproju
   ! nproju is the number of projectors for orbitals on which paw+u acts.
   ! Also used for local exact-exchange

  integer :: mesh_size
   ! Dimension of radial mesh for generic arrays contained in this pawtab datastructure
   ! The mesh is usually defined up to the PAW augmentation region boundary
   ! (+ a few additional points). May be different from pawrad%mesh_size

  integer :: core_mesh_size
   ! Dimension of radial mesh for core density

  integer :: coretau_mesh_size
   ! Dimension of radial mesh for core density

  integer :: vminus_mesh_size
   ! Dimension of radial mesh for vminushalf

  integer :: partialwave_mesh_size
   ! Dimension of radial mesh for partial waves (phi, tphi)
   ! May be different from pawrad%mesh_size and pawtab%mesh_size

  integer :: tnvale_mesh_size
   ! Dimension of radial mesh for tnvale

  integer :: mqgrid
   ! Number of points in the reciprocal space grid on which
   ! the radial functions (tcorespl, tvalespl...) are specified
   ! Same as psps%mqgrid_vl

  integer :: mqgrid_shp
   ! Number of points in the reciprocal space grid on which
   ! the radial shape functions (shapefncg) are given

  integer :: usespnorb
   ! usespnorb=0 ; no spin-orbit coupling
   ! usespnorb=1 ; spin-orbit coupling

  integer :: shape_lambda
   ! Lambda parameter in gaussian shapefunction (shape_type=2)

  integer :: shape_type
   ! Radial shape function type
   ! shape_type=-1 ; g(r)=numeric (read from psp file)
   ! shape_type= 1 ; g(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2 if r<=rshp, zero if r>rshp
   ! shape_type= 2 ; g(r)=exp[-(r/sigma)**lambda]
   ! shape_type= 3 ; gl(r)=Alpha(1,l)*jl(q(1,l)*r)+Alpha(2,l)*jl(q(2,l)*r) for each l

  integer :: useexexch
   ! useexexch=0 ; do not use local exact-exchange
   ! useexexch=1 ; use local exact-exchange

  integer :: usepawu
   ! usepawu= 0 ; do not use PAW+U formalism
   ! usepawu= 1 ; use PAW+U formalism (Full localized limit)
   ! usepawu= 2 ; use PAW+U formalism (Around Mean Field)
   ! usepawu= 3 ; use PAW+U formalism (Around Mean Field) - Alternative
   ! usepawu= 4 ; use PAW+U formalism (FLL) without polarization in the XC
   ! usepawu=-1 ; use PAW+U formalism (FLL) - No use of the occupation matrix - Experimental
   ! usepawu=-2 ; use PAW+U formalism (AMF) - No use of the occupation matrix - Experimental
   ! usepawu=10 ; use PAW+U within DMFT
   ! usepawu=14 ; use PAW+U within DMFT without polarization in the XC

  integer :: usepotzero
   ! usepotzero=0 if it is the Kresse-Joubert convention
   ! usepotzero=1 if it is the new convention
   ! usepotzero=2 if it is the PWscf convention

  integer :: usetcore
   ! Flag controling use of pseudized core density (0 if tncore=zero)

  integer :: usexcnhat
   ! 0 if compensation charge density is not included in XC terms
   ! 1 if compensation charge density is included in XC terms

!Real (real(dp)) scalars

  real(dp) :: beta
   ! contains the integral of the difference between vH[nZc] and vH[tnZc]

  real(dp) :: dncdq0
   ! Gives 1/q d(tNcore(q))/dq for q=0
   ! (tNcore(q) = FT of pseudo core density)

  real(dp) :: d2ncdq0
   ! Gives contribution of d2(tNcore(q))/d2q for q=0
   ! \int{(16/15)*pi^5*n(r)*r^6* dr}
   ! (tNcore(q) = FT of pseudo core density)

  real(dp) :: dnvdq0
   ! Gives 1/q d(tNvale(q))/dq for q=0
   ! (tNvale(q) = FT of pseudo valence density)

  real(dp) :: dtaucdq0
   ! Gives 1/q d(tTAUcore(q))/dq for q=0
   ! (tTAUcore(q) = FT of pseudo core kinetic density)

  real(dp) :: ex_cc
   ! Exchange energy for the core-core interaction of the Fock operator

  real(dp) :: exccore
   ! Exchange-correlation energy for the core density

  real(dp) :: exchmix
   ! mixing of exact exchange; default is 0.25 (PBE0)

  real(dp) :: f4of2_sla
   ! Ratio of Slater Integrals F4 and F2

  real(dp) :: f6of2_sla
   ! Ratio of Slater Integrals F6 and F4

  real(dp) :: jpawu
   ! jpawu
   ! Value of J parameter for paw+u for a given type.

  real(dp) :: rpaw
   ! Radius of PAW sphere

  real(dp) :: rshp
   ! Compensation charge radius (if r>rshp, g(r)=zero)

  real(dp) :: rcore
   ! Radius of core corrections (rcore >= rpaw)

  real(dp) :: rcoretau
   ! Radius of kinetic core corrections (rcoretau >= rpaw)

  real(dp) :: shape_sigma
   ! Sigma parameter in gaussian shapefunction (shape_type=2)

  real(dp) :: upawu
   ! upawu
   ! Value of U parameter for paw+u for a given type.

!Objects
  type(wvlpaw_type), pointer :: wvl
   !variable containing objects needed
   !for wvl+paw implementation
   !Warning: it is a pointer; it has to be allocated before use

!Integer arrays

  integer, allocatable :: indklmn(:,:)
   ! indklmn(8,lmn2_size)
   ! Array giving klm, kln, abs(il-jl), (il+jl), ilm and jlm, ilmn and jlmn for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn

  integer, allocatable :: indlmn(:,:)
   ! indlmn(6,lmn_size)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=lmn (if useylm=1)

  integer, allocatable :: klmntomn(:,:)
   ! klmntomn(4,lmn2_size)
   ! Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn
   ! NB: klmntomn is an application and not a bijection

  integer, allocatable :: kmix(:)
   ! kmix(lmnmix_sz)
   ! Indirect array selecting the klmn=(lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix

  integer, allocatable :: lnproju(:)
   ! lnproju(nproju) gives ln (index for phi) for each projectors on which U acts (PAW+U only)
   ! nproju is 1 or 2 and  is the number of projectors for correlated orbitals
   ! Also used for local exact-exchange

  integer, allocatable :: orbitals(:)
   ! gives the l quantum number per basis element

!Real (real(dp)) arrays

  real(dp), allocatable :: coredens(:)
   ! coredens(mesh_size)
   ! Gives the core density of the atom

  real(dp), allocatable :: coretau(:)
   ! coretau(mesh_size)
   ! Gives the kinetic energy density of the atom

  real(dp), allocatable :: dij0(:)
   ! dij0(lmn2_size)
   ! Part of the Dij term (non-local operator) completely
   ! calculated in the atomic data part

  real(dp), allocatable :: dltij(:)
   ! dltij(lmn2_size)
   ! Factor used to compute sums over klmn=(ilmn,jlmn)
   ! ((ilmn,ilmn) term has to be added once)
   ! dltij(klmn)=1 if ilmn=jlmn, else dltij(klmn)=2

  real(dp), allocatable :: dshpfunc(:,:,:)
   ! shapefunc(mesh_size,l_size,4)
   ! Gives the 4 first derivatives of  radial shape function
   ! for each l component; used only if shape_type=-1

  real(dp), allocatable :: eijkl(:,:)
   ! eijkl(lmn2_size,lmn2_size)
   ! Hartree kernel for the on-site terms (E_hartree=Sum_ijkl[rho_ij rho_kl e_ijkl])
   ! Used for Hartree and/or Fock contributions

  real(dp), allocatable :: eijkl_sr(:,:)
   ! eijkl_sr(lmn2_size,lmn2_size)
   ! Screened Hartree kernel for the on-site terms (E_hartree=Sum_ijkl[rho_ij rho_kl e_ijkl_sr])
   ! Used for screened Fock contributions

  real(dp), allocatable :: euijkl(:,:,:,:,:,:)
   ! euijkl(2,2,lmn_size,lmn_size,lmn_size,lmn_size)
   ! PAW+U kernel for the on-site terms ( E_PAW+U = 0.5 * Sum_s1s2 Sum_ijkl [rho_ij^s1 rho_kl^s2 euijkl^s1s2] )
   ! Contrary to eijkl and eijkl_sr, euijkl is not invariant with respect to the permutations i <--> j or k <--> l
   ! However, it is still invariant with respect to the permutation i,k <--> j,l, see pawpuxinit.F90
   ! Also, it depends on two spin indexes
   ! Used for PAW+U contributions

  real(dp), allocatable :: euij_fll(:)
   ! euij_fll(lmn2_size)
   ! Double counting part of the PAW+U kernel in the "fully localized limit".This term is only linear with respect to rho_ij,
   ! while euijkl is quadratic.
   ! Used for PAW+U contributions

  real(dp), allocatable :: ex_cvij(:)
  ! ex_cvij(lmn2_size))
  ! Onsite exact_exchange matrix elements for core-valence interactions of the Fock operator

  real(dp), allocatable :: fk(:,:)
   ! fk(6,4)
   ! Slater integrals used for local exact exchange

  real(dp), allocatable :: gammaij(:)
   ! gammaij(lmn2_size)
   ! background contribution from the densities

  real(dp), allocatable :: gnorm(:)
   ! gnorm(l_size)
   ! Give the the normalization factor of each radial shape function

  real(dp), allocatable :: kij(:)
  ! kij(lmn2_size))
  ! Onsite matrix elements <phi|\kinetic|phj>-<tphi|\kinetic|tphj>

  real(dp), allocatable :: nabla_ij(:,:,:)
   ! nabla_ij(3,lmn_size,lmn_size)
   ! Onsite matrix elements <phi|\nabla|phj>-<tphi|\nabla|tphj>

  real(dp), allocatable :: nablaphi(:,:)
   ! nablaphi(partialwave_mesh_size, basis_size)
   ! store the results of dphi/dr-(1/r)phi

  real(dp), allocatable :: phi(:,:)
   ! phi(partialwave_mesh_size, basis_size)
   ! Gives the paw electron wavefunctions on the radial grid

  real(dp), allocatable :: phiphj(:,:)
   ! phiphj(mesh_size,ij_size)
   ! Useful product Phi(:,i)*Phi(:,j)

  real(dp), allocatable :: phiphjint(:)
   ! phiphjint(ij_proj)
   ! Integration of Phi(:,i)*Phi(:,j) for DFT+U/local exact-exchange occupation matrix

  real(dp), allocatable :: ph0phiint(:)
   ! ph0phjint(ij_proj)
   ! Integration of Phi(:,1)*Phi(:,j) for LDA+DMFT projections

  real(dp), allocatable :: qgrid_shp(:)
   ! qgrid_shp(mqgrid_shp)
   ! Grid of points in reciprocal space on which the shape functions are given

  real(dp), allocatable :: qijl(:,:)
   ! qijl(l_size**2,lmn2_size)
   ! The qijl are the moments of the charge density difference between
   ! the AE and PS partial wave for each channel (i,j). They take part
   ! to the building of the compensation charge

  real(dp), allocatable :: rad_for_spline(:)
   ! rad_for_spline(mesh_size)
   ! Radial mesh used to spline quantities on radial mesh;
   ! Allocated and used only when
   !     shape_type=-1 (numerical shape function)
   !  or usedvloc=1 (use of vloc derivative)

  real(dp), allocatable :: rhoij0(:)
   ! rhoij0(lmn2_size)
   ! Initial guess for rhoij

  real(dp), allocatable :: shape_alpha(:,:)
   ! shape_alpha(2,l_size)
   ! Alpha_i parameters in Bessel shapefunctions (shape_type=3)

  real(dp), allocatable :: shape_q(:,:)
   ! shape_q(2,l_size)
   ! Q_i parameters in Bessel shapefunctions (shape_type=3)

  real(dp), allocatable :: shapefunc(:,:)
   ! shapefunc(mesh_size,l_size)
   ! Gives the normalized radial shape function for each l component

  real(dp), allocatable :: shapefncg(:,:,:)
   ! shapefncg(mqgrid_shp,2,l_size)
   ! Gives the spherical Fourier transform of the radial shape function
   ! for each l component (for each qgrid_shp(i)) + second derivative

  real(dp), allocatable :: sij(:)
   ! sij(lmn2_size)
   ! Nonlocal part of the overlap operator

  real(dp), allocatable :: tcoredens(:,:)
   ! tcoredens(core_mesh_size,1)
   ! Gives the pseudo core density of the atom
   ! In PAW+WVL:
   !  tcoredens(core_mesh_size,2:6)
   !  are the first to the fifth derivatives of the pseudo core density.

 real(dp), allocatable :: tcoretau(:)
   ! tcoretau(coretau_mesh_size)
   ! Gives the pseudo core kinetic energy density of the atom

  real(dp), allocatable :: tcorespl(:,:)
   ! tcorespl(mqgrid,2)
   ! Gives the pseudo core density in reciprocal space on a regular grid

  real(dp), allocatable :: tcoretauspl(:,:)
   ! tcoretauspl(mqgrid,2)
   ! Gives the pseudo kinetic core density in reciprocal space on a regular grid

  real(dp), allocatable :: tnablaphi(:,:)
   ! tphi(partialwave_mesh_size,basis_size)
   ! Gives, on the radial grid, the paw atomic pseudowavefunctions

  real(dp), allocatable :: tphi(:,:)
   ! tphi(partialwave_mesh_size,basis_size)
   ! Gives, on the radial grid, the paw atomic pseudowavefunctions

  real(dp), allocatable :: tphitphj(:,:)
   ! tphitphj(mesh_size,ij_size)
   ! Useful product tPhi(:,i)*tPhi(:,j)

  real(dp), allocatable :: tproj(:,:)
   ! non-local projectors

  real(dp), allocatable :: tvalespl(:,:)
   ! tvalespl(mqgrid,2)
   ! Gives the pseudo valence density in reciprocal space on a regular grid

  real(dp), allocatable :: Vee(:,:,:,:)
   ! PAW+U:
   ! Screened interaction matrix deduced from U and J parameters
   ! computed on the basis of orbitals on which U acts.

  real(dp), allocatable :: Vex(:,:,:,:,:)
   ! Local exact-exchange:
   ! Screened interaction matrix deduced from calculation of Slater integrals
   ! computed on the basis of orbitals on which local exact exchange acts.

  real(dp), allocatable :: vhtnzc(:)
   ! vhtnzc(mesh_size)
   ! Hartree potential for pseudized Zc density, v_H[\tilde{n}_{Zc}]
   ! read in from PAW file

  real(dp), allocatable :: VHnZC(:)
   ! VHnZC(mesh_size)
   ! Hartree potential for Zc density, v_H[n_{Zc}]
   ! constructed from core density in PAW file (see psp7in.F90)

  real(dp), allocatable :: vminushalf(:)
   ! vminushalf(mesh_size)
   ! External potential for LDA minus half calculation
   ! read in from PAW file

  real(dp), allocatable :: zioneff(:)
   ! zioneff(ij_proj)
   ! "Effective charge"*n "seen" at r_paw, deduced from Phi at r_paw, n:
   ! pricipal quantum number
   ! good approximation to model wave function outside PAW-sphere through

 end type pawtab_type

 public :: pawtab_free         ! Free memory
 public :: pawtab_nullify      ! Nullify content
 public :: pawtab_get_lsize    ! Get the max. l for a product of 2 partial waves
 public :: pawtab_set_flags    ! Set the value of the internal flags
 public :: pawtab_print        ! Printout of the object.
 public :: pawtab_bcast        ! MPI broadcast the object

 interface pawtab_nullify
   module procedure pawtab_nullify_0D
   module procedure pawtab_nullify_1D
 end interface pawtab_nullify

 interface pawtab_free
   module procedure pawtab_free_0D
   module procedure pawtab_free_1D
 end interface pawtab_free

 interface pawtab_set_flags
   module procedure pawtab_set_flags_0D
   module procedure pawtab_set_flags_1D
 end interface pawtab_set_flags
!!***

CONTAINS !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_nullify_0D
!! NAME
!!  pawtab_nullify_0D
!!
!! FUNCTION
!!  Nullify pointers and flags in a pawtab structure
!!
!! SIDE EFFECTS
!!  Pawtab<type(pawtab_type)>=PAW arrays tabulated.
!!                            Nullified in output
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_nullify_0D(Pawtab)

!Arguments ------------------------------------
!arrays
 type(Pawtab_type),intent(inout) :: Pawtab

!Local variables-------------------------------

! *************************************************************************

 !@Pawtab_type
 nullify(Pawtab%wvl)

 ! === Reset all flags and sizes ===

!Flags controlling optional arrays
 Pawtab%has_fock=0
 Pawtab%has_kij=0
 Pawtab%has_tproj=0
 Pawtab%has_tvale=0
 Pawtab%has_coretau=0
 Pawtab%has_vhtnzc=0
 Pawtab%has_vhnzc=0
 Pawtab%has_vminushalf=0
 Pawtab%has_nabla=0
 Pawtab%has_nablaphi=0
 Pawtab%has_shapefncg=0
 Pawtab%has_wvl=0

 Pawtab%usetcore=0
 Pawtab%usexcnhat=0
 Pawtab%useexexch=0
 Pawtab%usepawu=0
 Pawtab%usepotzero=0
 Pawtab%usespnorb=0
 Pawtab%mqgrid=0
 Pawtab%mqgrid_shp=0

 Pawtab%basis_size=0
 Pawtab%ij_proj=0
 Pawtab%ij_size=0
 Pawtab%lcut_size=0
 Pawtab%l_size=0
 Pawtab%lexexch=-1
 Pawtab%lmn_size=0
 Pawtab%lmn2_size=0
 Pawtab%lmnmix_sz=0
 Pawtab%lpawu=-1
 Pawtab%nproju=0
 Pawtab%mesh_size=0
 Pawtab%partialwave_mesh_size=0
 Pawtab%core_mesh_size=0
 Pawtab%coretau_mesh_size=0
 Pawtab%vminus_mesh_size=0
 Pawtab%tnvale_mesh_size=0
 Pawtab%shape_type=-10

end subroutine pawtab_nullify_0D
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_nullify_1D
!! NAME
!!  pawtab_nullify_1D
!!
!! FUNCTION
!!  Nullify all pointers in an array of pawtab data structures
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_nullify_1D(Pawtab)

!Arguments ------------------------------------
 type(pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ii,nn

! *************************************************************************

 !@pawtab_type

 nn=size(Pawtab)
 if (nn==0) return

 do ii=1,nn
   call pawtab_nullify_0D(Pawtab(ii))
 end do

end subroutine pawtab_nullify_1D
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_free_0D
!! NAME
!!  pawtab_free_0D
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a pawtab structure
!!
!! SIDE EFFECTS
!!  Pawtab<type(pawtab_type)>=PAW arrays tabulated.
!!  All allocated arrays in Pawtab are deallocated
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_free_0D(Pawtab)

!Arguments ------------------------------------
!arrays
 type(Pawtab_type),intent(inout) :: Pawtab

!Local variables-------------------------------

! *************************************************************************

 !@Pawtab_type

 if (allocated(Pawtab%indklmn))  then
  LIBPAW_DEALLOCATE(Pawtab%indklmn)
 end if
 if (allocated(Pawtab%indlmn))  then
  LIBPAW_DEALLOCATE(Pawtab%indlmn)
 end if
 if (allocated(Pawtab%klmntomn))  then
   LIBPAW_DEALLOCATE(Pawtab%klmntomn)
 end if
 if (allocated(Pawtab%kmix))  then
   LIBPAW_DEALLOCATE(Pawtab%kmix)
 end if
 if (allocated(Pawtab%lnproju))  then
   LIBPAW_DEALLOCATE(Pawtab%lnproju)
 end if
 if (allocated(Pawtab%coredens))  then
   LIBPAW_DEALLOCATE(Pawtab%coredens)
 end if
 if (allocated(Pawtab%coretau))  then
   LIBPAW_DEALLOCATE(Pawtab%coretau)
 end if
 if (allocated(Pawtab%dij0))  then
   LIBPAW_DEALLOCATE(Pawtab%dij0)
 end if
 if (allocated(Pawtab%dltij))  then
   LIBPAW_DEALLOCATE(Pawtab%dltij)
 end if
 if (allocated(Pawtab%dshpfunc))  then
   LIBPAW_DEALLOCATE(Pawtab%dshpfunc)
 end if
 if (allocated(Pawtab%eijkl))  then
   LIBPAW_DEALLOCATE(Pawtab%eijkl)
 end if
 if (allocated(Pawtab%eijkl_sr))  then
   LIBPAW_DEALLOCATE(Pawtab%eijkl_sr)
 end if
 if (allocated(Pawtab%euijkl))  then
   LIBPAW_DEALLOCATE(Pawtab%euijkl)
 end if
 if (allocated(Pawtab%euij_fll))  then
   LIBPAW_DEALLOCATE(Pawtab%euij_fll)
 end if
 if (allocated(Pawtab%fk))  then
   LIBPAW_DEALLOCATE(Pawtab%fk)
 end if
 if (allocated(Pawtab%gammaij))  then
   LIBPAW_DEALLOCATE(Pawtab%gammaij)
 end if
 if (allocated(Pawtab%gnorm))  then
   LIBPAW_DEALLOCATE(Pawtab%gnorm)
 end if
 if (allocated(Pawtab%ex_cvij))  then
   LIBPAW_DEALLOCATE(Pawtab%ex_cvij)
 end if
 if (allocated(Pawtab%kij))  then
   LIBPAW_DEALLOCATE(Pawtab%kij)
 end if
 if (allocated(Pawtab%nabla_ij))  then
   LIBPAW_DEALLOCATE(Pawtab%nabla_ij)
 end if
 if (allocated(Pawtab%nablaphi))  then
   LIBPAW_DEALLOCATE(Pawtab%nablaphi)
 end if
 if (allocated(Pawtab%orbitals)) then
   LIBPAW_DEALLOCATE(Pawtab%orbitals)
 end if
 if (allocated(Pawtab%phi))  then
   LIBPAW_DEALLOCATE(Pawtab%phi)
 end if
 if (allocated(Pawtab%phiphj))  then
   LIBPAW_DEALLOCATE(Pawtab%phiphj)
 end if
 if (allocated(Pawtab%phiphjint))  then
   LIBPAW_DEALLOCATE(Pawtab%phiphjint)
 end if
 if (allocated(Pawtab%ph0phiint))  then
   LIBPAW_DEALLOCATE(Pawtab%ph0phiint)
 end if
 if (allocated(Pawtab%qgrid_shp))  then
   LIBPAW_DEALLOCATE(Pawtab%qgrid_shp)
 end if
 if (allocated(Pawtab%qijl))  then
   LIBPAW_DEALLOCATE(Pawtab%qijl)
 end if
 if (allocated(Pawtab%rad_for_spline))  then
   LIBPAW_DEALLOCATE(Pawtab%rad_for_spline)
 end if
 if (allocated(Pawtab%rhoij0))  then
   LIBPAW_DEALLOCATE(Pawtab%rhoij0)
 end if
 if (allocated(Pawtab%shape_alpha))  then
   LIBPAW_DEALLOCATE(Pawtab%shape_alpha)
 end if
 if (allocated(Pawtab%shape_q))  then
   LIBPAW_DEALLOCATE(Pawtab%shape_q)
 end if
 if (allocated(Pawtab%shapefunc))  then
   LIBPAW_DEALLOCATE(Pawtab%shapefunc)
 end if
 if (allocated(Pawtab%shapefncg))  then
   LIBPAW_DEALLOCATE(Pawtab%shapefncg)
 end if
 if (allocated(Pawtab%sij))  then
   LIBPAW_DEALLOCATE(Pawtab%sij)
 end if
 if (allocated(Pawtab%tcoredens))  then
   LIBPAW_DEALLOCATE(Pawtab%tcoredens)
 end if
 if (allocated(Pawtab%tcoretau))  then
   LIBPAW_DEALLOCATE(Pawtab%tcoretau)
 end if
 if (allocated(Pawtab%tcorespl))  then
   LIBPAW_DEALLOCATE(Pawtab%tcorespl)
 end if
 if (allocated(Pawtab%tcoretauspl))  then
   LIBPAW_DEALLOCATE(Pawtab%tcoretauspl)
 end if
 if (allocated(Pawtab%tnablaphi))  then
   LIBPAW_DEALLOCATE(Pawtab%tphi)
 end if
 if (allocated(Pawtab%tphi))  then
   LIBPAW_DEALLOCATE(Pawtab%tphi)
 end if
 if (allocated(Pawtab%tphitphj))  then
   LIBPAW_DEALLOCATE(Pawtab%tphitphj)
 end if
 if (allocated(Pawtab%tproj)) then
   LIBPAW_DEALLOCATE(Pawtab%tproj)
 end if
 if (allocated(Pawtab%tvalespl))  then
   LIBPAW_DEALLOCATE(Pawtab%tvalespl)
 end if
 if (allocated(Pawtab%vee))  then
   LIBPAW_DEALLOCATE(Pawtab%vee)
 end if
 if (allocated(Pawtab%Vex))  then
   LIBPAW_DEALLOCATE(Pawtab%Vex)
 end if
 if (allocated(Pawtab%vhtnzc))  then
   LIBPAW_DEALLOCATE(Pawtab%vhtnzc)
 end if
 if (allocated(Pawtab%VHnZC))  then
   LIBPAW_DEALLOCATE(Pawtab%VHnZC)
 end if
 if (allocated(Pawtab%vminushalf))  then
   LIBPAW_DEALLOCATE(Pawtab%vminushalf)
 end if
 if (allocated(Pawtab%zioneff))  then
   LIBPAW_DEALLOCATE(Pawtab%zioneff)
 end if

 call wvlpaw_free(Pawtab%wvl)

 ! === Reset all flags and sizes ===

!CAUTION: do not reset these flags
!They are set from input data and must be kept
!Pawtab%has_kij=0
!Pawtab%has_tproj=0
!Pawtab%has_coretau=0
!Pawtab%has_tvale=0
!Pawtab%has_vhtnzc=0
!Pawtab%has_vhnzc=0
!Pawtab%has_vminushalf=0
!Pawtab%has_nabla=0
!Pawtab%has_nablaphi=0
!Pawtab%has_shapefncg=0
!Pawtab%has_wvl=0

 Pawtab%usetcore=0
 Pawtab%usexcnhat=0
 Pawtab%useexexch=0
 Pawtab%usepawu=0
 Pawtab%usepotzero=0
 Pawtab%usespnorb=0
 Pawtab%mqgrid=0
 Pawtab%mqgrid_shp=0

 Pawtab%basis_size=0
 Pawtab%ij_proj=0
 Pawtab%ij_size=0
 Pawtab%lcut_size=0
 Pawtab%l_size=0
 Pawtab%lexexch=-1
 Pawtab%lmn_size=0
 Pawtab%lmn2_size=0
 Pawtab%lmnmix_sz=0
 Pawtab%lpawu=-1
 Pawtab%nproju=0
 Pawtab%mesh_size=0
 Pawtab%partialwave_mesh_size=0
 Pawtab%core_mesh_size=0
 Pawtab%coretau_mesh_size=0
 Pawtab%vminus_mesh_size=0
 Pawtab%tnvale_mesh_size=0
 Pawtab%shape_type=-10

end subroutine pawtab_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_free_1D
!! NAME
!!  pawtab_free_1D
!!
!! FUNCTION
!!  Destroy (deallocate) all pointers in an array of pawtab data structures
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_free_1D(Pawtab)

!Arguments ------------------------------------
 type(pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ii,nn

! *************************************************************************

 !@pawtab_type

 nn=size(Pawtab)
 if (nn==0) return

 do ii=1,nn
   call pawtab_free_0D(Pawtab(ii))
 end do

end subroutine pawtab_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_set_flags_0D
!! NAME
!!  pawtab_set_flags_0D
!!
!! FUNCTION
!!  Set flags controlling optional arrays in a pawtab datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_set_flags_0D(Pawtab,has_coretau,has_fock,has_kij,has_tproj,has_tvale,has_vhnzc,&
&                              has_vhtnzc,has_nabla,has_nablaphi,has_shapefncg,has_vminushalf,has_wvl)

!Arguments ------------------------------------
 integer,intent(in),optional :: has_coretau,has_fock,has_kij,has_tproj,has_tvale
 integer,intent(in),optional :: has_vhnzc,has_vhtnzc,has_vminushalf
 integer,intent(in),optional :: has_nabla,has_nablaphi,has_shapefncg,has_wvl
 type(pawtab_type),intent(inout) :: Pawtab

!Local variables-------------------------------

! *************************************************************************

 !@pawtab_type

 Pawtab%has_fock      =0
 Pawtab%has_kij       =0
 Pawtab%has_tproj     =0
 Pawtab%has_tvale     =0
 Pawtab%has_coretau   =0
 Pawtab%has_vhnzc     =0
 Pawtab%has_vhtnzc    =0
 Pawtab%has_nabla     =0
 Pawtab%has_nablaphi  =0
 Pawtab%has_shapefncg =0
 Pawtab%has_vminushalf=0
 Pawtab%has_wvl       =0
 if (present(has_fock))      Pawtab%has_fock=has_fock
 if (present(has_kij))       Pawtab%has_kij=has_kij
 if (present(has_tproj))     Pawtab%has_tproj=has_tproj
 if (present(has_tvale))     Pawtab%has_tvale=has_tvale
 if (present(has_coretau))   Pawtab%has_coretau=has_coretau
 if (present(has_vhnzc))     Pawtab%has_vhnzc=has_vhnzc
 if (present(has_vhtnzc))    Pawtab%has_vhtnzc=has_vhtnzc
 if (present(has_nabla))     Pawtab%has_nabla=has_nabla
 if (present(has_nablaphi))  Pawtab%has_nablaphi=has_nablaphi
 if (present(has_shapefncg) )Pawtab%has_shapefncg=has_shapefncg
 if (present(has_vminushalf))Pawtab%has_vminushalf=has_vminushalf
 if (present(has_wvl))       Pawtab%has_wvl=has_wvl

end subroutine pawtab_set_flags_0D
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_set_flags_1D
!! NAME
!!  pawtab_set_flags_1D
!!
!! FUNCTION
!!  Set flags controlling optional arrays in an array of pawtab datastructures
!! if (present(has_tvale))    Pawtab%has_tvale=has_tvale

!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_set_flags_1D(Pawtab,has_coretau,has_fock,has_kij,has_tproj,has_tvale,has_vhnzc,&
&                              has_vhtnzc,has_nabla,has_nablaphi,has_shapefncg,has_vminushalf,has_wvl)

!Arguments ------------------------------------
 integer,intent(in),optional :: has_coretau,has_fock,has_kij,has_tproj,has_tvale,has_vhnzc,has_vhtnzc
 integer,intent(in),optional :: has_nabla,has_nablaphi,has_shapefncg,has_vminushalf,has_wvl
 type(pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ii,nn

! *************************************************************************

 !@pawtab_type

 nn=size(Pawtab)
 if (nn==0) return

 do ii=1,nn
   Pawtab(ii)%has_fock      =0
   Pawtab(ii)%has_kij       =0
   Pawtab(ii)%has_tproj     =0
   Pawtab(ii)%has_tvale     =0
   Pawtab(ii)%has_coretau   =0
   Pawtab(ii)%has_vhnzc     =0
   Pawtab(ii)%has_vhtnzc    =0
   Pawtab(ii)%has_nabla     =0
   Pawtab(ii)%has_nablaphi  =0
   Pawtab(ii)%has_shapefncg =0
   Pawtab(ii)%has_vminushalf=0
   Pawtab(ii)%has_wvl       =0
   if (present(has_fock))      Pawtab(ii)%has_fock=has_fock
   if (present(has_kij))       Pawtab(ii)%has_kij=has_kij
   if (present(has_tproj))     Pawtab(ii)%has_tproj=has_tproj
   if (present(has_tvale))     Pawtab(ii)%has_tvale=has_tvale
   if (present(has_coretau))   Pawtab(ii)%has_coretau=has_coretau
   if (present(has_vhnzc))     Pawtab(ii)%has_vhnzc=has_vhnzc
   if (present(has_vhtnzc))    Pawtab(ii)%has_vhtnzc=has_vhtnzc
   if (present(has_nabla))     Pawtab(ii)%has_nabla=has_nabla
   if (present(has_nablaphi))  Pawtab(ii)%has_nablaphi=has_nablaphi
   if (present(has_shapefncg)) Pawtab(ii)%has_shapefncg=has_shapefncg
   if (present(has_vminushalf))Pawtab(ii)%has_vminushalf=has_vminushalf
   if (present(has_wvl))       Pawtab(ii)%has_wvl=has_wvl
 end do

end subroutine pawtab_set_flags_1D
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_print
!! NAME
!! pawtab_print
!!
!! FUNCTION
!!  Print out the content of a pawtab datastructure
!!
!! INPUTS
!!  Pawtab<pawtab_type> Only for PAW, TABulated data initialized at start
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_bethe_salpeter,m_screening_driver,m_sigma_driver,m_wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_print(Pawtab,header,unit,prtvol,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
!arrays
 type(Pawtab_type) :: Pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: ityp,ntypat,my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt   =ab_out ; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(msg,'(6a)')&
&  ' ==================================== ',ch10,&
&  ' ==== Info on PAW TABulated data ==== ',ch10,&
&  ' ==================================== ',ch10
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ntypat=SIZE(Pawtab(:))

 do ityp=1,ntypat

  ! Print out integer values (dimensions)
  write(msg,'(a)')'                                 '
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a)')'  ****************************** '
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4,a)')'  **** Atom type ',ityp,' ****   '
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a)')'  ****************************** '
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Number of (n,l) elements ....................... ',Pawtab(ityp)%basis_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Number of (l,m,n) elements ..................... ',Pawtab(ityp)%lmn_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Number of (i,j) elements (packed form) ......... ',Pawtab(ityp)%ij_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Max L+1 leading to non-zero Gaunt .............. ',Pawtab(ityp)%l_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Max L+1 leading to non-zero Gaunt (pawlcutd) ... ',Pawtab(ityp)%lcut_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  lmn2_size ...................................... ',Pawtab(ityp)%lmn2_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  lmnmix_sz ...................................... ',Pawtab(ityp)%lmnmix_sz
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Size of radial mesh ............................ ',Pawtab(ityp)%mesh_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Size of radial mesh for partial waves........... ',Pawtab(ityp)%partialwave_mesh_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Size of radial mesh for [pseudo] core density... ',Pawtab(ityp)%core_mesh_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Size of radial mesh for [pseudo] kin core density',Pawtab(ityp)%coretau_mesh_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Size of radial mesh for pseudo valence density.. ',Pawtab(ityp)%tnvale_mesh_size
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  No of Q-points for tcorespl/tvalespl/tcoretauspl ',Pawtab(ityp)%mqgrid
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  No of Q-points for the radial shape functions .. ',Pawtab(ityp)%mqgrid_shp
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Radial shape function type ..................... ',Pawtab(ityp)%shape_type
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  shape_lambda ................................... ',Pawtab(ityp)%shape_lambda
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Use pseudized core density ..................... ',Pawtab(ityp)%usetcore
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Option for the use of hat density in XC terms .. ',Pawtab(ityp)%usexcnhat
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Use DFT+U ...................................... ',Pawtab(ityp)%usepawu
  call wrtout(my_unt,msg,my_mode)
  if (Pawtab(ityp)%usepawu/=0) then
    write(msg,'(a,i4)')'  L on which U is applied ........................ ',Pawtab(ityp)%lpawu
    call wrtout(my_unt,msg,my_mode)
  end if
  write(msg,'(a,i4)')'  Use Local Exact exchange ....................... ',Pawtab(ityp)%useexexch
  call wrtout(my_unt,msg,my_mode)
  if (Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,i4)')'  L on which local exact-exchange is applied ..... ',Pawtab(ityp)%lexexch
    call wrtout(my_unt,msg,my_mode)
  end if
  if (Pawtab(ityp)%usepawu/=0.or.Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,i4)')'  Number of (i,j) elements for PAW+U or EXX ..... ',Pawtab(ityp)%ij_proj
    call wrtout(my_unt,msg,my_mode)
    write(msg,'(a,i4)')'  Number of projectors on which U or EXX acts .... ',Pawtab(ityp)%nproju
    call wrtout(my_unt,msg,my_mode)
  end if
  write(msg,'(a,i4)')'  Use potential zero ............................. ',Pawtab(ityp)%usepotzero
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Use spin-orbit coupling ........................ ',Pawtab(ityp)%usespnorb
  call wrtout(my_unt,msg,my_mode)

  ! "Has" flags
  write(msg,'(a,i4)')'  Has Fock  ...................................... ',Pawtab(ityp)%has_fock
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has kij   ...................................... ',Pawtab(ityp)%has_kij
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has tproj ...................................... ',Pawtab(ityp)%has_tproj
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has tvale ...................................... ',Pawtab(ityp)%has_tvale
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has coretau .................................... ',Pawtab(ityp)%has_coretau
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has vhtnzc ..................................... ',Pawtab(ityp)%has_vhtnzc
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has vhnzc ...................................... ',Pawtab(ityp)%has_vhnzc
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has vminushalf ................................. ',Pawtab(ityp)%has_vminushalf
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has nabla ...................................... ',Pawtab(ityp)%has_nabla
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has nablaphi ................................... ',Pawtab(ityp)%has_nablaphi
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has shapefuncg ................................. ',Pawtab(ityp)%has_shapefncg
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,i4)')'  Has wvl ........................................ ',Pawtab(ityp)%has_wvl
  call wrtout(my_unt,msg,my_mode)
  !
  ! Real scalars
  write(msg,'(a,es16.8)')'  beta ............................................',Pawtab(ityp)%beta
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,es16.8)')'  1/q d(tNcore(q))/dq for q=0 .....................',Pawtab(ityp)%dncdq0
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,es16.8)')'  d^2(tNcore(q))/dq^2 for q=0 .....................',Pawtab(ityp)%d2ncdq0
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,es16.8)')'  1/q d(tNvale(q))/dq for q=0 .....................',Pawtab(ityp)%dnvdq0
  call wrtout(my_unt,msg,my_mode)
  if (Pawtab(ityp)%has_coretau/=0) then
    write(msg,'(a,es16.8)')'  1/q d(tTAUcore(q))/dq for q=0 ...................',Pawtab(ityp)%dtaucdq0
    call wrtout(my_unt,msg,my_mode)
  end if
  if (Pawtab(ityp)%has_fock/=0) then
    write(msg,'(a,es16.8)')'  Core-core Fock energy  ..........................',Pawtab(ityp)%ex_cc
    call wrtout(my_unt,msg,my_mode)
  end if
  write(msg,'(a,es16.8)')'  XC energy for the core density ..................',Pawtab(ityp)%exccore
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,es16.8)')'  Radius of the PAW sphere ........................',Pawtab(ityp)%rpaw
  call wrtout(my_unt,msg,my_mode)
  write(msg,'(a,es16.8)')'  Compensation charge radius (if >rshp, g(r)=0) ...',Pawtab(ityp)%rshp !(if r>rshp, g(r)=zero)
  call wrtout(my_unt,msg,my_mode)
  if (Pawtab(ityp)%shape_type==2) then
   write(msg,'(a,es16.8)')'  Sigma parameter in gaussian shape function ......',Pawtab(ityp)%shape_sigma !(shape_type=2)
   call wrtout(my_unt,msg,my_mode)
  end if
  if (Pawtab(ityp)%usepawu/=0) then
   write(msg,'(a,es16.8)')'  Value of the U parameter [eV] ...................',Pawtab(ityp)%upawu*Ha_eV
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,es16.8)')'  Value of the J parameter [eV] ...................',Pawtab(ityp)%jpawu*Ha_eV
   call wrtout(my_unt,msg,my_mode)
  end if
  if (Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,es16.8)')'  Mixing of exact exchange (PBE0) .................',Pawtab(ityp)%exchmix
    call wrtout(my_unt,msg,my_mode)
  end if
 if (associated(Pawtab(ityp)%wvl)) then
   write(msg,'(a,es16.8)')'  WARNING: This Pawtab structure contains WVL data.'
   call wrtout(my_unt,msg,my_mode)
 end if

 end do ! ityp

 ! The other (huge) arrays are not reported..

end subroutine pawtab_print
!!***

!----------------------------------------------------------------------

!!****f* m_pawtap/pawtab_get_lsize
!! NAME
!!  pawtab_get_lsize
!!
!! FUNCTION
!!  From an array of pawtab datastructures, get, for each atom, the value
!!  of "l_size" parameter.
!!  l_size is the maximum value of l accessible by a product of 2 partial waves;
!!  it may be cut by dtset%pawlcutd parameter
!!
!! INPUTS
!!   [mpi_atmtab(:)]=--optional-- indexes of the atoms treated by current proc
!!   natom= number of atoms (may be a local or absolute number of atoms)
!!   typat(:)= list of atom types
!!
!! OUTPUT
!!   l_size_atm(natom)=output array of l_size values (for each atom)
!!
!! PARENTS
!!      m_bethe_salpeter,m_classify_bands,m_d2frnl,m_exc_analyze,m_nonlinear
!!      m_paw_mkaewf,m_paw_mkrho,m_respfn_driver,m_scfcv_core
!!      m_screening_driver,m_sigma_driver,m_wfd,m_wfk_analyze
!!
!! CHILDREN
!!
!! NOTES
!!  This function returns an allocatable integer array which may be allocated
!!  on the fly.
!!
!! SOURCE

subroutine pawtab_get_lsize(Pawtab,l_size_atm,natom,typat, &
&                           mpi_atmtab) ! Optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 integer,intent(in) :: typat(:)
 integer,optional,intent(in) :: mpi_atmtab(:)
 integer,allocatable,intent(inout) :: l_size_atm(:)
 type(pawtab_type),intent(in) :: pawtab(:)

!Local variables-------------------------------
 integer :: ia,ityp,natom_typat
 character(len=100) :: msg

! *************************************************************************

 !@pawtab_type

 natom_typat=count(typat>0)
 if (size(pawtab)<maxval(typat)) then
   msg='error on pawtab size!'
   ABI_BUG(msg)
 end if

 if (.not.allocated(l_size_atm)) then
   LIBPAW_ALLOCATE(l_size_atm,(natom))
 else if (size(l_size_atm)/=natom) then
   LIBPAW_DEALLOCATE(l_size_atm)
   LIBPAW_ALLOCATE(l_size_atm,(natom))
 end if

 if (natom==0) return

 if (natom==natom_typat) then

!First case: sequential mode
   do ia=1,natom
     ityp=typat(ia)
     l_size_atm(ia)=pawtab(ityp)%lcut_size
   end do

 else

!2nd case: parallel mode
   if (.not.present(mpi_atmtab)) then
     msg='optional args error!'
     ABI_BUG(msg)
   end if
   do ia=1,natom
     ityp=typat(mpi_atmtab(ia))
     l_size_atm(ia)=pawtab(ityp)%lcut_size
   end do

 end if

end subroutine pawtab_get_lsize
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_bcast
!! NAME
!! pawtab_bcast
!!
!! FUNCTION
!! Communicate pawtab data to all processors
!!
!! INPUTS
!! comm_mpi= communicator used to broadcast data
!! [only_from_file]= (optional, default=FALSE)
!!    If true, only data obtained at the level of the reading
!!    of the PAW dataset file are broadcasted
!!
!! SIDE EFFECTS
!!  pawtab=<type pawtab_type>=a pawtab datastructure
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_bcast(pawtab,comm_mpi,only_from_file)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm_mpi
 logical,optional,intent(in) :: only_from_file
 type(pawtab_type),intent(inout) :: pawtab

!Local variables-------------------------------
!scalars
 integer :: ierr,ii,me,nn_dpr,nn_dpr_arr,nn_int,nn_int_arr
 integer :: siz_indklmn,siz_indlmn,siz_klmntomn,siz_kmix,siz_lnproju,siz_orbitals
 integer :: siz_coredens,siz_coretau,siz_dij0,siz_dltij,siz_dshpfunc,siz_eijkl,siz_eijkl_sr
 integer :: siz_euijkl,siz_euij_fll,siz_fk,siz_gammaij,siz_gnorm,siz_fock,siz_kij,siz_nabla_ij
 integer :: siz_nablaphi,siz_phi,siz_phiphj,siz_phiphjint,siz_ph0phiint,siz_qgrid_shp,siz_qijl
 integer :: siz_rad_for_spline,siz_rhoij0,siz_shape_alpha,siz_shape_q,siz_shapefunc
 integer :: siz_shapefncg,siz_sij,siz_tcoredens,siz_tcoretau,siz_tcorespl,siz_tcoretauspl
 integer :: siz_tnablaphi,siz_tphi,siz_tphitphj,siz_tproj,siz_tvalespl,siz_vee,siz_vex
 integer :: siz_vhtnzc,siz_vhnzc,siz_vminushalf,siz_zioneff,siz_wvlpaw,siz_wvl_pngau
 integer :: siz_wvl_parg,siz_wvl_pfac,siz_wvl_rholoc_rad,siz_wvl_rholoc_d,sz1,sz2
 logical :: full_broadcast
 character (len=500) :: msg,msg0
!arrays
 integer :: nn(4)
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)

!*************************************************************************

 me=xmpi_comm_rank(comm_mpi)
 full_broadcast=.true.;if (present(only_from_file)) full_broadcast=(.not.only_from_file)

 nn_int=0 ; nn_int_arr=0 ; nn_dpr=0 ; nn_dpr_arr=0

!=========================================================================
!Compute the amount of data to communicate
!=========================================================================

 if (me==0) then
   msg=''

!Integers (read from psp file)
!-------------------------------------------------------------------------
!  basis_size,has_coretau,has_fock,has_kij,has_shapefncg,has_nabla,has_nablaphi,has_tproj
!  has_tvale,has_vhtnzc,has_vhnzc,has_vminushalf,has_wvl,ij_size,l_size,lmn_size,lmn2_size
!  mesh_size,partialwave_mesh_size,core_mesh_size,coretau_mesh_size,vminus_mesh_size
!  tnvale_mesh_size,mqgrid,shape_lambda,shape_type,usetcore,usexcnhat
   nn_int=nn_int+28

!Integers (depending on the parameters of the calculation)
!-------------------------------------------------------------------------
!  ij_proj,lcut_size,lexexch,lmnmix_sz,lpawu,mqgrid_shp,nproju,useexexch,usepawu,usepotzero,usespnorb
   if (full_broadcast) nn_int=nn_int+11

!Reals (read from psp file)
!-------------------------------------------------------------------------
!  beta,dncdq0,d2ncdq0,dnvdq0,dtaucdq0,ex_cc,exccore,rpaw,rshp,rcore,rcoretau,shape_sigma
   nn_dpr=nn_dpr+12

!Reals (depending on the parameters of the calculation)
!-------------------------------------------------------------------------
!  exchmix,f4of2_sla,f6of2_sla,jpawu,upawu
   if (full_broadcast) nn_dpr=nn_dpr+5

!Integers arrays (read from psp file)
!-------------------------------------------------------------------------
   siz_indlmn=0 ; siz_orbitals=0
   nn_int=nn_int+2
   if (allocated(pawtab%indlmn)) then
     siz_indlmn=size(pawtab%indlmn)                  !(6,lmn_size)
     if (siz_indlmn/=6*pawtab%lmn_size) msg=trim(msg)//' indlmn'
     nn_int_arr=nn_int_arr+siz_indlmn
   end if
   if (allocated(pawtab%orbitals)) then
     siz_orbitals=size(pawtab%orbitals)              !(basis_size)
     if (siz_orbitals/=pawtab%basis_size) msg=trim(msg)//' orbitals'
     nn_int_arr=nn_int_arr+siz_orbitals
   end if

!Integers arrays (depending on the parameters of the calculation)
!-------------------------------------------------------------------------
   siz_indklmn=0 ; siz_klmntomn=0 ; siz_kmix=0 ; siz_lnproju=0
   if (full_broadcast) then
     nn_int=nn_int+4
     if (allocated(pawtab%indklmn)) then
       siz_indklmn=size(pawtab%indklmn)             !(6,lmn2_size)
       if (siz_indklmn/=8*pawtab%lmn2_size) msg=trim(msg)//' indklmn'
       nn_int_arr=nn_int_arr+siz_indklmn
     end if
     if (allocated(pawtab%klmntomn)) then
       siz_klmntomn=size(pawtab%klmntomn)           !(4,lmn2_size)
       if (siz_klmntomn/=4*pawtab%lmn2_size) msg=trim(msg)//' klmntomn'
       nn_int_arr=nn_int_arr+siz_klmntomn
     end if
     if (allocated(pawtab%kmix)) then
       siz_kmix=size(pawtab%kmix)                   !(lmnmix_sz)
       if (siz_kmix/=6*pawtab%lmnmix_sz) msg=trim(msg)//' kmix'
       nn_int_arr=nn_int_arr+siz_kmix
     end if
     if (allocated(pawtab%lnproju)) then
       siz_lnproju=size(pawtab%lnproju)             !(nproju)
       if (siz_lnproju/=pawtab%nproju) msg=trim(msg)//' lnproju'
       nn_int_arr=nn_int_arr+siz_lnproju
     end if
   end if ! full_broadcast

!Reals arrays (read from psp file)
!-------------------------------------------------------------------------
   siz_coredens=0 ; siz_coretau=0    ; siz_dij0=0     ; siz_kij=0        ; siz_fock=0
   siz_phi=0      ; siz_nablaphi=0   ;  siz_rhoij0=0   ; siz_shape_alpha=0
   siz_shape_q=0  ; siz_shapefunc=0  ; siz_tcoredens=0; siz_tcoretau=0
   siz_tcorespl=0 ; siz_tcoretauspl=0; siz_tphi=0     ; siz_tnablaphi=0  ; siz_tproj=0
   siz_tvalespl=0 ; siz_vhtnzc=0     ; siz_vhnzc=0    ; siz_vminushalf=0
   nn_int=nn_int+22

   if (allocated(pawtab%coredens)) then
     siz_coredens=size(pawtab%coredens)             !(core_mesh_size)
     if (siz_coredens/=pawtab%core_mesh_size) &
&      msg=trim(msg)//' coredens'
     nn_dpr=nn_dpr+siz_coredens
   end if
   if (allocated(pawtab%coretau)) then
     siz_coretau=size(pawtab%coretau)             !(coretau_mesh_size)
     if (siz_coretau/=pawtab%coretau_mesh_size) &
&      msg=trim(msg)//' coretau'
     nn_dpr=nn_dpr+siz_coretau
   end if
   if (allocated(pawtab%dij0)) then
     siz_dij0=size(pawtab%dij0)                     !(lmn2_size)
     if (siz_dij0/=pawtab%lmn2_size) msg=trim(msg)//' dij0'
     nn_dpr=nn_dpr+siz_dij0
   end if
   if (allocated(pawtab%ex_cvij)) then
     siz_fock=size(pawtab%ex_cvij)                  !(lmn2_size)
     if (siz_fock/=pawtab%lmn2_size) msg=trim(msg)//' fock'
     nn_dpr=nn_dpr+siz_fock
   end if
   if (allocated(pawtab%kij)) then
     siz_kij=size(pawtab%kij)                       !(lmn2_size)
     if (siz_kij/=pawtab%lmn2_size) msg=trim(msg)//' kij'
     nn_dpr=nn_dpr+siz_kij
   end if
   if (allocated(pawtab%nablaphi)) then
     siz_phi=size(pawtab%nablaphi)                  !(partialwave_mesh_size, basis_size)
     if (siz_nablaphi/=pawtab%partialwave_mesh_size*pawtab%basis_size) msg=trim(msg)//' nablaphi'
     nn_dpr=nn_dpr+siz_nablaphi
   end if
   if (allocated(pawtab%phi)) then
     siz_phi=size(pawtab%phi)                       !(partialwave_mesh_size, basis_size)
     if (siz_phi/=pawtab%partialwave_mesh_size*pawtab%basis_size) msg=trim(msg)//' phi'
     nn_dpr=nn_dpr+siz_phi
   end if
   if (allocated(pawtab%rhoij0)) then
     siz_rhoij0=size(pawtab%rhoij0)                 !(lmn2_size)
     if (siz_rhoij0/=pawtab%lmn2_size) msg=trim(msg)//' rhoij0'
     nn_dpr=nn_dpr+siz_rhoij0
   end if
   if (allocated(pawtab%shape_alpha)) then
     siz_shape_alpha=size(pawtab%shape_alpha)       !(2,l_size)
     if (siz_shape_alpha/=pawtab%l_size*2) msg=trim(msg)//' shape_alpha'
     nn_dpr=nn_dpr+siz_shape_alpha
   end if
   if (allocated(pawtab%shape_q)) then
     siz_shape_q=size(pawtab%shape_q)               !(2,l_size)
     if (siz_shape_q/=pawtab%l_size*2) msg=trim(msg)//' shape_q'
     nn_dpr=nn_dpr+siz_shape_q
   end if
   if (allocated(pawtab%shapefunc)) then
     siz_shapefunc=size(pawtab%shapefunc)           !(mesh_size,l_size)
     if (siz_shapefunc/=pawtab%mesh_size*pawtab%l_size) msg=trim(msg)//' shapefunc'
     nn_dpr=nn_dpr+siz_shapefunc
   end if
   if (allocated(pawtab%tcoredens)) then
     siz_tcoredens=size(pawtab%tcoredens)           !(core_mesh_size,1 or 6)
     if (siz_tcoredens/=pawtab%core_mesh_size.and.siz_tcoredens/=6*pawtab%core_mesh_size) &
&      msg=trim(msg)//' tcoredens'
     nn_dpr=nn_dpr+siz_tcoredens
   end if
   if (allocated(pawtab%tcoretau)) then
     siz_tcoretau=size(pawtab%tcoretau)             !(coretau_mesh_size,1)
     if (siz_tcoretau/=pawtab%coretau_mesh_size) &
&      msg=trim(msg)//' tcoretau'
     nn_dpr=nn_dpr+siz_tcoretau
   end if
   if (allocated(pawtab%tcorespl)) then
     siz_tcorespl=size(pawtab%tcorespl)             !(mqgrid,2)
     if (siz_tcorespl/=pawtab%mqgrid*2) msg=trim(msg)//' tcorespl'
     nn_dpr=nn_dpr+siz_tcorespl
   end if
   if (allocated(pawtab%tcoretauspl)) then
     siz_tcoretauspl=size(pawtab%tcoretauspl)       !(mqgrid,2)
     if (siz_tcoretauspl/=pawtab%mqgrid*2.and.siz_tcoretauspl/=0) msg=trim(msg)//' tcoretauspl'
     nn_dpr=nn_dpr+siz_tcoretauspl
   end if
   if (allocated(pawtab%tnablaphi)) then
     siz_tnablaphi=size(pawtab%tnablaphi)           !(partialwave_mesh_size, basis_size)
     if (siz_tnablaphi/=pawtab%partialwave_mesh_size*pawtab%basis_size) msg=trim(msg)//' tnablaphi'
     nn_dpr=nn_dpr+siz_tnablaphi
   end if
   if (allocated(pawtab%tphi)) then
     siz_tphi=size(pawtab%tphi)                     !(partialwave_mesh_size, basis_size)
     if (siz_tphi/=pawtab%partialwave_mesh_size*pawtab%basis_size) msg=trim(msg)//' tphi'
     nn_dpr=nn_dpr+siz_tphi
   end if
   if (allocated(pawtab%tproj)) then
     siz_tproj=size(pawtab%tproj)                   !(???,basis_size)
    if (mod(siz_tproj,pawtab%basis_size)/=0) msg=trim(msg)//' tproj'
     nn_dpr=nn_dpr+siz_tproj
   end if
   if (allocated(pawtab%tvalespl)) then
     siz_tvalespl=size(pawtab%tvalespl)             !(mqgrid or mesh_size or tnvale_mesh_size,2)
     if (siz_tvalespl/=2*pawtab%mqgrid.and.siz_tvalespl/=2*pawtab%mesh_size.and. &
&        siz_tvalespl/=2*pawtab%tnvale_mesh_size) msg=trim(msg)//' tvalespl'
     nn_dpr=nn_dpr+siz_tvalespl
   end if
   if (allocated(pawtab%vhtnzc)) then
     siz_vhtnzc=size(pawtab%vhtnzc)                 !(mesh_size)
     if (siz_vhtnzc/=pawtab%mesh_size) msg=trim(msg)//' vhtnzc'
     nn_dpr=nn_dpr+siz_vhtnzc
   end if
   if (allocated(pawtab%vhnzc)) then
     siz_vhnzc=size(pawtab%vhnzc)                   !(mesh_size)
     if (siz_vhnzc/=pawtab%mesh_size) msg=trim(msg)//' vhnzc'
     nn_dpr=nn_dpr+siz_vhnzc
   end if
   if (allocated(pawtab%vminushalf)) then
     siz_vminushalf=size(pawtab%vminushalf)         !(mesh_size)
     if (siz_vminushalf/=pawtab%vminus_mesh_size) msg=trim(msg)//' vvminushalf'
     nn_dpr=nn_dpr+siz_vminushalf
   end if

!Reals arrays (depending on the parameters of the calculation)
!-------------------------------------------------------------------------
   siz_dltij=0    ; siz_dshpfunc=0
   siz_eijkl=0    ; siz_eijkl_sr=0 ; siz_euijkl=0   ; siz_euij_fll=0
   siz_fk=0       ; siz_gammaij=0  ; siz_gnorm=0    ; siz_nabla_ij=0
   siz_phiphj=0   ; siz_phiphjint=0; siz_ph0phiint=0
   siz_qgrid_shp=0; siz_qijl=0     ; siz_rad_for_spline=0
   siz_shapefncg=0; siz_sij=0      ; siz_tphitphj=0
   siz_vee=0      ; siz_vex=0      ; siz_zioneff=0
   if (full_broadcast) then
     nn_int=nn_int+22
     if (allocated(pawtab%dltij)) then
       siz_dltij=size(pawtab%dltij)                   !(lmn2_size)
       if (siz_dltij/=pawtab%lmn2_size) msg=trim(msg)//' dltij'
       nn_dpr=nn_dpr+siz_dltij
     end if
     if (allocated(pawtab%dshpfunc)) then
       siz_dshpfunc=size(pawtab%dshpfunc)             !(mesh_size,l_size,4)
       if (siz_dshpfunc/=pawtab%mesh_size*pawtab%l_size*4) msg=trim(msg)//' dshpfunc'
       nn_dpr=nn_dpr+siz_dshpfunc
     end if
     if (allocated(pawtab%eijkl)) then
       siz_eijkl=size(pawtab%eijkl)                   !(lmn2_size,lmn2_size)
       if (siz_eijkl/=pawtab%lmn2_size*pawtab%lmn2_size) msg=trim(msg)//' eijkl'
       nn_dpr=nn_dpr+siz_eijkl
     end if
     if (allocated(pawtab%eijkl_sr)) then
       siz_eijkl_sr=size(pawtab%eijkl_sr)             !(lmn2_size,lmn2_size)
       if (siz_eijkl_sr/=pawtab%lmn2_size*pawtab%lmn2_size) msg=trim(msg)//' eijkl_sr'
       nn_dpr=nn_dpr+siz_eijkl_sr
     end if
     if (allocated(pawtab%euijkl)) then
       siz_euijkl=size(pawtab%euijkl)                 !(2,2,lmn_size,lmn_size,lmn_size,lmn_size)
       if (siz_euijkl/=4*pawtab%lmn_size*pawtab%lmn_size*pawtab%lmn_size*pawtab%lmn_size) msg=trim(msg)//' euijkl'
       nn_dpr=nn_dpr+siz_euijkl
     end if
     if (allocated(pawtab%euij_fll)) then
       siz_euij_fll=size(pawtab%euij_fll)             !(2,2,lmn_size,lmn_size,lmn_size,lmn_size)
       if (siz_euij_fll/=4*pawtab%lmn_size*pawtab%lmn_size*pawtab%lmn_size*pawtab%lmn_size) msg=trim(msg)//' euij_fll'
       nn_dpr=nn_dpr+siz_euij_fll
     end if
     if (allocated(pawtab%fk)) then
       siz_fk=size(pawtab%fk)                         !(6,4)
       if (siz_fk/=24) msg=trim(msg)//' fk'
       nn_dpr=nn_dpr+siz_fk
     end if
     if (allocated(pawtab%gammaij)) then
       siz_gammaij=size(pawtab%gammaij)               !(l_size)
       if (siz_gammaij/=pawtab%l_size) msg=trim(msg)//' gammaij'
       nn_dpr=nn_dpr+siz_gammaij
     end if
     if (allocated(pawtab%gnorm)) then
       siz_gnorm=size(pawtab%gnorm)                   !(l_size)
       if (siz_gnorm/=pawtab%l_size) msg=trim(msg)//' gnorm'
       nn_dpr=nn_dpr+siz_gnorm
     end if
     if (allocated(pawtab%nabla_ij)) then
       siz_nabla_ij=size(pawtab%nabla_ij)             !(3,lmn_size,lmn_size)
       if (siz_nabla_ij/=pawtab%lmn_size) msg=trim(msg)//' nabla_ij'
       nn_dpr=nn_dpr+siz_nabla_ij
     end if
     if (allocated(pawtab%phiphj)) then
       siz_phiphj=size(pawtab%phiphj)                 !(mesh_size,ij_size)
       if (siz_phiphj/=pawtab%mesh_size*pawtab%ij_size) msg=trim(msg)//' phiphj'
       nn_dpr=nn_dpr+siz_phiphj
     end if
     if (allocated(pawtab%phiphjint)) then
       siz_phiphjint=size(pawtab%phiphjint)           !(ij_proj)
       if (siz_phiphjint/=pawtab%ij_proj) msg=trim(msg)//' phiphjint'
       nn_dpr=nn_dpr+siz_phiphjint
     end if
     if (allocated(pawtab%ph0phiint)) then
       siz_ph0phiint=size(pawtab%ph0phiint)           !(ij_proj)
       if (siz_ph0phiint/=pawtab%ij_proj) msg=trim(msg)//' ph0phiint'
       nn_dpr=nn_dpr+siz_ph0phiint
     end if
     if (allocated(pawtab%qgrid_shp)) then
       siz_qgrid_shp=size(pawtab%qgrid_shp)           !(mqgrid_shp)
       if (siz_qgrid_shp/=pawtab%mqgrid_shp) msg=trim(msg)//' qgrid_shp'
       nn_dpr=nn_dpr+siz_qgrid_shp
     end if
     if (allocated(pawtab%qijl)) then
       siz_qijl=size(pawtab%qijl)                     !(l_size**2,lmn2_size)
       if (siz_qijl/=pawtab%l_size**2*pawtab%lmn2_size) msg=trim(msg)//' qijl'
       nn_dpr=nn_dpr+siz_qijl
     end if
     if (allocated(pawtab%rad_for_spline)) then
       siz_rad_for_spline=size(pawtab%rad_for_spline) !(mesh_size)
       if (siz_rad_for_spline/=pawtab%mesh_size) msg=trim(msg)//' rad_for_spline'
       nn_dpr=nn_dpr+siz_rad_for_spline
     end if
     if (allocated(pawtab%shapefncg)) then
       siz_shapefncg=size(pawtab%shapefncg)           !(mqgrid_shp,2,l_size)
       if (siz_shapefncg/=2*pawtab%mqgrid_shp*pawtab%l_size) msg=trim(msg)//' shapefncg'
       nn_dpr=nn_dpr+siz_shapefncg
     end if
     if (allocated(pawtab%sij)) then
       siz_sij=size(pawtab%sij)                       !(lmn2_size)
       if (siz_sij/=pawtab%lmn2_size) msg=trim(msg)//' sij'
       nn_dpr=nn_dpr+siz_sij
     end if
     if (allocated(pawtab%tphitphj)) then
       siz_tphitphj=size(pawtab%tphitphj)             !(mesh_size,ij_size)
       if (siz_tphitphj/=pawtab%mesh_size*pawtab%ij_size) msg=trim(msg)//' tphitphj'
       nn_dpr=nn_dpr+siz_tphitphj
     end if
     if (allocated(pawtab%vee)) then
       siz_vee=size(pawtab%vee)                       !(2*lpawu+1,2*lpawu+1,2*lpawu+1,2*lpawu+1)
       if (siz_vee/=(2*pawtab%lpawu+1)**4) msg=trim(msg)//' vee'
       nn_dpr=nn_dpr+siz_vee
     end if
     if (allocated(pawtab%vex)) then
       siz_vex=size(pawtab%vex)                       !(2*lexexch+1,2*lexexch+1,2*lexexch+1,2*lexexch+1,4)
       if (siz_vex/=4*(2*pawtab%lpawu+1)**4) msg=trim(msg)//' vex'
       nn_dpr=nn_dpr+siz_vex
     end if
     if (allocated(pawtab%zioneff)) then
       siz_zioneff=size(pawtab%zioneff)               !(ij_proj)
       if (siz_zioneff/=pawtab%ij_proj) msg=trim(msg)//' zioneff'
       nn_dpr=nn_dpr+siz_zioneff
     end if
   end if ! full_broadcast

!Datastructures (read from psp file)
!-------------------------------------------------------------------------
   siz_wvl_pngau=0 ; siz_wvl_parg=0 ; siz_wvl_pfac=0
   siz_wvl_rholoc_rad=0 ; siz_wvl_rholoc_d=0
   siz_wvlpaw=0
   nn_int=nn_int+1
   if (associated(pawtab%wvl)) then
     siz_wvlpaw=1
     nn_int=nn_int+3
!    wvl%npspcode_init_guess,wvl%ptotgau
     nn_int=nn_int+2
     if (allocated(pawtab%wvl%pngau)) then
       siz_wvl_pngau=size(pawtab%wvl%pngau)         !(basis_size)
       if (siz_wvl_pngau/=pawtab%basis_size) msg=trim(msg)//' wvl_pngau'
       nn_int_arr=nn_int_arr+siz_wvl_pngau
     end if
     if (allocated(pawtab%wvl%parg)) then
       siz_wvl_parg=size(pawtab%wvl%parg)          !(2,ptotgau)
       if (siz_wvl_parg/=2*pawtab%wvl%ptotgau) msg=trim(msg)//' wvl_parg'
       nn_dpr_arr=nn_dpr_arr+siz_wvl_parg
     end if
     if (allocated(pawtab%wvl%pfac)) then
       siz_wvl_pfac=size(pawtab%wvl%pfac )         !(2,ptotgau)
       if (siz_wvl_pfac/=2*pawtab%wvl%ptotgau) msg=trim(msg)//' wvl_pfac'
       nn_dpr_arr=nn_dpr_arr+siz_wvl_pfac
     end if
!    wvl%rholoc%msz
     nn_int=nn_int+3
     if (pawtab%wvl%rholoc%msz>0) then
       if (allocated(pawtab%wvl%rholoc%rad)) then
         siz_wvl_rholoc_rad=size(pawtab%wvl%rholoc%rad) !(msz)
         if (siz_wvl_rholoc_rad/=pawtab%wvl%rholoc%msz) msg=trim(msg)//' wvl_rholoc_rad'
         nn_dpr_arr=nn_dpr_arr+siz_wvl_rholoc_rad
       end if
       if (allocated(pawtab%wvl%rholoc%d)) then
         siz_wvl_rholoc_d=size(pawtab%wvl%rholoc%d)     !(msz,4)
         if (siz_wvl_rholoc_d/=4*pawtab%wvl%rholoc%msz) msg=trim(msg)//' wvl_rholoc_d'
         nn_dpr_arr=nn_dpr_arr+siz_wvl_rholoc_d
       end if
     end if
   end if

!Datastructures (depending on the parameters of the calculation)
!-------------------------------------------------------------------------
!  Nothing

!  Are the sizes OK ?
   if (trim(msg)/='') then
     write(msg0,'(3a)') &
&     'There is a problem with the size of the following array(s):',ch10,trim(msg)
     ABI_BUG(msg0)
   end if

 end if ! me=0

!Broadcast the sizes of buffers
!=========================================================================

 if (me==0) then
   nn(1)=nn_int ; nn(2)=nn_int_arr
   nn(3)=nn_dpr ; nn(4)=nn_dpr_arr
 end if
 call xmpi_bcast(nn,0,comm_mpi,ierr)
 if (me/=0) then
   nn_int=nn(1) ; nn_int_arr=nn(2)
   nn_dpr=nn(3) ; nn_dpr_arr=nn(4)
 end if

!Broadcast all the integer: sizes, integer scalars, integer arrays
!=========================================================================

 LIBPAW_ALLOCATE(list_int,(nn_int+nn_int_arr))

!Fill the buffer of the sender
!-------------------------------------------------------------------------
 if (me==0) then
   ii=1

!First the data read from a psp file
!...................................

!Sizes of arrays (read from psp file)
   list_int(ii)=siz_indlmn  ;ii=ii+1
   list_int(ii)=siz_orbitals  ;ii=ii+1
   list_int(ii)=siz_coredens  ;ii=ii+1
   list_int(ii)=siz_coretau  ;ii=ii+1
   list_int(ii)=siz_dij0  ;ii=ii+1
   list_int(ii)=siz_kij  ;ii=ii+1
   list_int(ii)=siz_fock  ;ii=ii+1
   list_int(ii)=siz_phi  ;ii=ii+1
   list_int(ii)=siz_nablaphi; ii=ii+1
   list_int(ii)=siz_rhoij0  ;ii=ii+1
   list_int(ii)=siz_shape_alpha  ;ii=ii+1
   list_int(ii)=siz_shape_q  ;ii=ii+1
   list_int(ii)=siz_shapefunc  ;ii=ii+1
   list_int(ii)=siz_tcoredens  ;ii=ii+1
   list_int(ii)=siz_tcoretau  ;ii=ii+1
   list_int(ii)=siz_tcorespl  ;ii=ii+1
   list_int(ii)=siz_tcoretauspl  ;ii=ii+1
   list_int(ii)=siz_tphi  ;ii=ii+1
   list_int(ii)=siz_tnablaphi  ;ii=ii+1
   list_int(ii)=siz_tproj  ;ii=ii+1
   list_int(ii)=siz_tvalespl  ;ii=ii+1
   list_int(ii)=siz_vhtnzc  ;ii=ii+1
   list_int(ii)=siz_vhnzc  ;ii=ii+1
   list_int(ii)=siz_vminushalf  ;ii=ii+1
   list_int(ii)=siz_wvlpaw  ;ii=ii+1
!Integers (read from psp file)
   list_int(ii)=pawtab%basis_size  ;ii=ii+1
   list_int(ii)=pawtab%has_fock  ;ii=ii+1
   list_int(ii)=pawtab%has_kij  ;ii=ii+1
   list_int(ii)=pawtab%has_shapefncg  ;ii=ii+1
   list_int(ii)=pawtab%has_nabla  ;ii=ii+1
   list_int(ii)=pawtab%has_nablaphi ; ii=ii+1
   list_int(ii)=pawtab%has_tproj  ;ii=ii+1
   list_int(ii)=pawtab%has_tvale  ;ii=ii+1
   list_int(ii)=pawtab%has_coretau  ;ii=ii+1
   list_int(ii)=pawtab%has_vhtnzc  ;ii=ii+1
   list_int(ii)=pawtab%has_vhnzc  ;ii=ii+1
   list_int(ii)=pawtab%has_vminushalf  ;ii=ii+1
   list_int(ii)=pawtab%has_wvl  ;ii=ii+1
   list_int(ii)=pawtab%ij_size  ;ii=ii+1
   list_int(ii)=pawtab%l_size  ;ii=ii+1
   list_int(ii)=pawtab%lmn_size  ;ii=ii+1
   list_int(ii)=pawtab%lmn2_size  ;ii=ii+1
   list_int(ii)=pawtab%mesh_size  ;ii=ii+1
   list_int(ii)=pawtab%partialwave_mesh_size  ;ii=ii+1
   list_int(ii)=pawtab%core_mesh_size  ;ii=ii+1
   list_int(ii)=pawtab%coretau_mesh_size  ;ii=ii+1
   list_int(ii)=pawtab%vminus_mesh_size  ;ii=ii+1
   list_int(ii)=pawtab%tnvale_mesh_size  ;ii=ii+1
   list_int(ii)=pawtab%mqgrid  ;ii=ii+1
   list_int(ii)=pawtab%shape_lambda  ;ii=ii+1
   list_int(ii)=pawtab%shape_type  ;ii=ii+1
   list_int(ii)=pawtab%usetcore  ;ii=ii+1
   list_int(ii)=pawtab%usexcnhat  ;ii=ii+1
!Integer arrays (read from psp file)
   if (siz_indlmn>0) then
     list_int(ii:ii+siz_indlmn-1)=reshape(pawtab%indlmn,(/siz_indlmn/))
     ii=ii+siz_indlmn
   end if
   if (siz_orbitals>0) then
     list_int(ii:ii+siz_orbitals-1)=pawtab%orbitals(1:siz_orbitals)
     ii=ii+siz_orbitals
   end if
!Integers in datastructures (read from psp file)
   if (siz_wvlpaw==1) then
     list_int(ii)=siz_wvl_pngau  ;ii=ii+1
     list_int(ii)=siz_wvl_parg  ;ii=ii+1
     list_int(ii)=siz_wvl_pfac  ;ii=ii+1
     list_int(ii)=pawtab%wvl%npspcode_init_guess  ;ii=ii+1
     list_int(ii)=pawtab%wvl%ptotgau  ;ii=ii+1
     if (siz_wvl_pngau>0) then
       list_int(ii:ii+siz_wvl_pngau-1)=pawtab%wvl%pngau(1:siz_wvl_pngau)
       ii=ii+siz_wvl_pngau
     end if
     list_int(ii)=siz_wvl_rholoc_rad  ;ii=ii+1
     list_int(ii)=siz_wvl_rholoc_d  ;ii=ii+1
     list_int(ii)=pawtab%wvl%rholoc%msz  ;ii=ii+1
   end if

!Then the data initialized later
!...................................
   if (full_broadcast) then

!Sizes of arrays
     list_int(ii)=siz_indklmn  ;ii=ii+1
     list_int(ii)=siz_klmntomn  ;ii=ii+1
     list_int(ii)=siz_kmix  ;ii=ii+1
     list_int(ii)=siz_lnproju  ;ii=ii+1
     list_int(ii)=siz_dltij  ;ii=ii+1
     list_int(ii)=siz_dshpfunc  ;ii=ii+1
     list_int(ii)=siz_eijkl  ;ii=ii+1
     list_int(ii)=siz_eijkl_sr  ;ii=ii+1
     list_int(ii)=siz_euijkl  ;ii=ii+1
     list_int(ii)=siz_euij_fll  ;ii=ii+1
     list_int(ii)=siz_fk  ;ii=ii+1
     list_int(ii)=siz_gammaij ;ii=ii+1
     list_int(ii)=siz_gnorm  ;ii=ii+1
     list_int(ii)=siz_nabla_ij  ;ii=ii+1
     list_int(ii)=siz_phiphj  ;ii=ii+1
     list_int(ii)=siz_phiphjint  ;ii=ii+1
     list_int(ii)=siz_ph0phiint  ;ii=ii+1
     list_int(ii)=siz_qgrid_shp  ;ii=ii+1
     list_int(ii)=siz_qijl  ;ii=ii+1
     list_int(ii)=siz_rad_for_spline  ;ii=ii+1
     list_int(ii)=siz_shapefncg  ;ii=ii+1
     list_int(ii)=siz_sij  ;ii=ii+1
     list_int(ii)=siz_tphitphj  ;ii=ii+1
     list_int(ii)=siz_vee  ;ii=ii+1
     list_int(ii)=siz_vex  ;ii=ii+1
     list_int(ii)=siz_zioneff  ;ii=ii+1
!Integers
     list_int(ii)=pawtab%ij_proj  ;ii=ii+1
     list_int(ii)=pawtab%lcut_size  ;ii=ii+1
     list_int(ii)=pawtab%lexexch  ;ii=ii+1
     list_int(ii)=pawtab%lmnmix_sz  ;ii=ii+1
     list_int(ii)=pawtab%lpawu  ;ii=ii+1
     list_int(ii)=pawtab%mqgrid_shp  ;ii=ii+1
     list_int(ii)=pawtab%nproju  ;ii=ii+1
     list_int(ii)=pawtab%useexexch  ;ii=ii+1
     list_int(ii)=pawtab%usepawu  ;ii=ii+1
     list_int(ii)=pawtab%usepotzero ;ii=ii+1
     list_int(ii)=pawtab%usespnorb ;ii=ii+1
!Integer arrays
     if (siz_indklmn>0) then
       list_int(ii:ii+siz_indklmn-1)=reshape(pawtab%indklmn,(/siz_indklmn/))
       ii=ii+siz_indklmn
     end if
     if (siz_klmntomn>0) then
       list_int(ii:ii+siz_klmntomn-1)=reshape(pawtab%klmntomn,(/siz_klmntomn/))
       ii=ii+siz_klmntomn
     end if
     if (siz_kmix>0) then
       list_int(ii:ii+siz_kmix-1)=pawtab%kmix(1:siz_kmix)
       ii=ii+siz_kmix
     end if
     if (siz_lnproju>0) then
       list_int(ii:ii+siz_lnproju-1)=pawtab%lnproju(1:siz_lnproju)
       ii=ii+siz_lnproju
     end if
   end if ! full_broadcast
   ii=ii-1

   if (ii/=nn_int+nn_int_arr) then
     msg='the number of loaded integers is not correct!'
     ABI_BUG(msg)
   end if

 end if ! me=0

!Perfom the communication
!-------------------------------------------------------------------------

 call xmpi_bcast(list_int,0,comm_mpi,ierr)

!Fill the receiver from the buffer
!-------------------------------------------------------------------------
 if (me/=0) then
   ii=1

!First the data read from a psp file
!...................................

!Sizes of arrays (read from psp file)
   siz_indlmn=list_int(ii)  ;ii=ii+1
   siz_orbitals=list_int(ii)  ;ii=ii+1
   siz_coredens=list_int(ii)  ;ii=ii+1
   siz_coretau=list_int(ii)  ;ii=ii+1
   siz_dij0=list_int(ii)  ;ii=ii+1
   siz_kij=list_int(ii)  ;ii=ii+1
   siz_fock=list_int(ii)  ;ii=ii+1
   siz_phi=list_int(ii)  ;ii=ii+1
   siz_nablaphi=list_int(ii)  ;ii=ii+1
   siz_rhoij0=list_int(ii)  ;ii=ii+1
   siz_shape_alpha=list_int(ii)  ;ii=ii+1
   siz_shape_q=list_int(ii)  ;ii=ii+1
   siz_shapefunc=list_int(ii)  ;ii=ii+1
   siz_tcoredens=list_int(ii)  ;ii=ii+1
   siz_tcoretau=list_int(ii)  ;ii=ii+1
   siz_tcorespl=list_int(ii)  ;ii=ii+1
   siz_tcoretauspl=list_int(ii)  ;ii=ii+1
   siz_tphi=list_int(ii)  ;ii=ii+1
   siz_tnablaphi=list_int(ii)  ;ii=ii+1
   siz_tproj=list_int(ii)  ;ii=ii+1
   siz_tvalespl=list_int(ii)  ;ii=ii+1
   siz_vhtnzc=list_int(ii)  ;ii=ii+1
   siz_vhnzc=list_int(ii)  ;ii=ii+1
   siz_vminushalf=list_int(ii)  ;ii=ii+1
   siz_wvlpaw=list_int(ii)  ;ii=ii+1
!Integers (read from psp file)
   pawtab%basis_size=list_int(ii)  ;ii=ii+1
   pawtab%has_fock=list_int(ii)  ;ii=ii+1
   pawtab%has_kij=list_int(ii)  ;ii=ii+1
   pawtab%has_shapefncg=list_int(ii)  ;ii=ii+1
   pawtab%has_nabla=list_int(ii)  ;ii=ii+1
   pawtab%has_nablaphi=list_int(ii) ; ii=ii+1
   pawtab%has_tproj=list_int(ii)  ;ii=ii+1
   pawtab%has_tvale=list_int(ii)  ;ii=ii+1
   pawtab%has_coretau=list_int(ii)  ;ii=ii+1
   pawtab%has_vhtnzc=list_int(ii)  ;ii=ii+1
   pawtab%has_vhnzc=list_int(ii)  ;ii=ii+1
   pawtab%has_vminushalf=list_int(ii)  ;ii=ii+1
   pawtab%has_wvl=list_int(ii)  ;ii=ii+1
   pawtab%ij_size=list_int(ii)  ;ii=ii+1
   pawtab%l_size=list_int(ii)  ;ii=ii+1
   pawtab%lmn_size=list_int(ii)  ;ii=ii+1
   pawtab%lmn2_size=list_int(ii)  ;ii=ii+1
   pawtab%mesh_size=list_int(ii)  ;ii=ii+1
   pawtab%partialwave_mesh_size=list_int(ii)  ;ii=ii+1
   pawtab%core_mesh_size=list_int(ii)  ;ii=ii+1
   pawtab%coretau_mesh_size=list_int(ii)  ;ii=ii+1
   pawtab%vminus_mesh_size=list_int(ii)  ;ii=ii+1
   pawtab%tnvale_mesh_size=list_int(ii)  ;ii=ii+1
   pawtab%mqgrid=list_int(ii)  ;ii=ii+1
   pawtab%shape_lambda=list_int(ii)  ;ii=ii+1
   pawtab%shape_type=list_int(ii)  ;ii=ii+1
   pawtab%usetcore=list_int(ii)  ;ii=ii+1
   pawtab%usexcnhat=list_int(ii)  ;ii=ii+1
!Integer arrays (read from psp file)
   if (allocated(pawtab%indlmn)) then
     LIBPAW_DEALLOCATE(pawtab%indlmn)
   end if
   if (siz_indlmn>0) then
     LIBPAW_ALLOCATE(pawtab%indlmn,(6,pawtab%lmn_size))
     pawtab%indlmn=reshape(list_int(ii:ii+siz_indlmn-1),(/6,pawtab%lmn_size/))
     ii=ii+siz_indlmn
   end if
   if (allocated(pawtab%orbitals)) then
     LIBPAW_DEALLOCATE(pawtab%orbitals)
   end if
   if (siz_orbitals>0) then
     LIBPAW_ALLOCATE(pawtab%orbitals,(pawtab%basis_size))
     pawtab%orbitals=list_int(ii:ii+pawtab%basis_size-1)
     ii=ii+siz_orbitals
   end if
!Integers in datastructures (read from psp file)
   if (siz_wvlpaw==1) then
     call wvlpaw_allocate(pawtab%wvl)
     siz_wvl_pngau=list_int(ii)  ;ii=ii+1
     siz_wvl_parg=list_int(ii)  ;ii=ii+1
     siz_wvl_pfac=list_int(ii)  ;ii=ii+1
     pawtab%wvl%npspcode_init_guess=list_int(ii)  ;ii=ii+1
     pawtab%wvl%ptotgau=list_int(ii)  ;ii=ii+1
     if (allocated(pawtab%wvl%pngau)) then
       LIBPAW_DEALLOCATE(pawtab%wvl%pngau)
     end if
     if (siz_wvl_pngau>0) then
       LIBPAW_ALLOCATE(pawtab%wvl%pngau,(pawtab%basis_size))
       pawtab%wvl%pngau=list_int(ii:ii+pawtab%basis_size-1)
       ii=ii+siz_wvl_pngau
     end if
     siz_wvl_rholoc_rad=list_int(ii)  ;ii=ii+1
     siz_wvl_rholoc_d=list_int(ii)  ;ii=ii+1
     pawtab%wvl%rholoc%msz=list_int(ii)  ;ii=ii+1
   end if

!Then the data initialized later
!...................................
   if (full_broadcast) then

!Sizes of arrays
     siz_indklmn=list_int(ii)  ;ii=ii+1
     siz_klmntomn=list_int(ii)  ;ii=ii+1
     siz_kmix=list_int(ii)  ;ii=ii+1
     siz_lnproju=list_int(ii)  ;ii=ii+1
     siz_dltij=list_int(ii)  ;ii=ii+1
     siz_dshpfunc=list_int(ii)  ;ii=ii+1
     siz_eijkl=list_int(ii)  ;ii=ii+1
     siz_eijkl_sr=list_int(ii)  ;ii=ii+1
     siz_euijkl=list_int(ii)  ;ii=ii+1
     siz_euij_fll=list_int(ii)  ;ii=ii+1
     siz_fk=list_int(ii)  ;ii=ii+1
     siz_gammaij=list_int(ii)  ;ii=ii+1
     siz_gnorm=list_int(ii)  ;ii=ii+1
     siz_nabla_ij=list_int(ii)  ;ii=ii+1
     siz_phiphj=list_int(ii)  ;ii=ii+1
     siz_phiphjint=list_int(ii)  ;ii=ii+1
     siz_ph0phiint=list_int(ii)  ;ii=ii+1
     siz_qgrid_shp=list_int(ii)  ;ii=ii+1
     siz_qijl=list_int(ii)  ;ii=ii+1
     siz_rad_for_spline=list_int(ii)  ;ii=ii+1
     siz_shapefncg=list_int(ii)  ;ii=ii+1
     siz_sij=list_int(ii)  ;ii=ii+1
     siz_tphitphj=list_int(ii)  ;ii=ii+1
     siz_vee=list_int(ii)  ;ii=ii+1
     siz_vex=list_int(ii)  ;ii=ii+1
     siz_zioneff=list_int(ii)  ;ii=ii+1
!Integers
     pawtab%ij_proj=list_int(ii)  ;ii=ii+1
     pawtab%lcut_size=list_int(ii)  ;ii=ii+1
     pawtab%lexexch=list_int(ii)  ;ii=ii+1
     pawtab%lmnmix_sz=list_int(ii)  ;ii=ii+1
     pawtab%lpawu=list_int(ii)  ;ii=ii+1
     pawtab%mqgrid_shp=list_int(ii)  ;ii=ii+1
     pawtab%nproju=list_int(ii)  ;ii=ii+1
     pawtab%useexexch=list_int(ii)  ;ii=ii+1
     pawtab%usepawu=list_int(ii)  ;ii=ii+1
     pawtab%usepotzero=list_int(ii) ;ii=ii+1
     pawtab%usespnorb=list_int(ii) ;ii=ii+1
!Integer arrays
     if (allocated(pawtab%indklmn)) then
       LIBPAW_DEALLOCATE(pawtab%indklmn)
     end if
     if (siz_indklmn>0) then
       LIBPAW_ALLOCATE(pawtab%indklmn,(8,pawtab%lmn2_size))
       pawtab%indklmn=reshape(list_int(ii:ii+siz_indklmn-1),(/8,pawtab%lmn2_size/))
       ii=ii+siz_indklmn
     end if
     if (allocated(pawtab%klmntomn)) then
       LIBPAW_DEALLOCATE(pawtab%klmntomn)
     end if
     if (siz_klmntomn>0) then
       LIBPAW_ALLOCATE(pawtab%klmntomn,(4,pawtab%lmn2_size))
       pawtab%klmntomn=reshape(list_int(ii:ii+siz_klmntomn-1),(/4,pawtab%lmn2_size/))
       ii=ii+siz_klmntomn
     end if
     if (allocated(pawtab%kmix)) then
       LIBPAW_DEALLOCATE(pawtab%kmix)
     end if
     if (siz_kmix>0) then
       LIBPAW_ALLOCATE(pawtab%kmix,(pawtab%lmnmix_sz))
       pawtab%kmix=list_int(ii:ii+pawtab%lmnmix_sz-1)
       ii=ii+siz_kmix
     end if
     if (allocated(pawtab%lnproju)) then
       LIBPAW_DEALLOCATE(pawtab%lnproju)
     end if
     if (siz_lnproju>0) then
       LIBPAW_ALLOCATE(pawtab%lnproju,(pawtab%nproju))
       pawtab%lnproju=list_int(ii:ii+pawtab%nproju-1)
       ii=ii+siz_lnproju
     end if
   end if ! full_broadcast
   ii=ii-1

   if (ii/=nn_int+nn_int_arr) then
     msg='the number of broadcasted integers is not correct!'
     ABI_BUG(msg)
   end if

 end if ! me/=0
 LIBPAW_DEALLOCATE(list_int)

!Broadcast all the reals
!=========================================================================

 LIBPAW_ALLOCATE(list_dpr,(nn_dpr+nn_dpr_arr))

!Fill the buffer of the sender
!-------------------------------------------------------------------------
 if (me==0) then
   ii=1

!First the data read from a psp file
!...................................

!Reals (read from psp file)
   list_dpr(ii)=pawtab%beta    ;ii=ii+1
   list_dpr(ii)=pawtab%dncdq0  ;ii=ii+1
   list_dpr(ii)=pawtab%d2ncdq0  ;ii=ii+1
   list_dpr(ii)=pawtab%dnvdq0  ;ii=ii+1
   list_dpr(ii)=pawtab%dtaucdq0  ;ii=ii+1
   list_dpr(ii)=pawtab%ex_cc   ;ii=ii+1
   list_dpr(ii)=pawtab%exccore  ;ii=ii+1
   list_dpr(ii)=pawtab%rpaw  ;ii=ii+1
   list_dpr(ii)=pawtab%rshp  ;ii=ii+1
   list_dpr(ii)=pawtab%rcore  ;ii=ii+1
   list_dpr(ii)=pawtab%rcoretau  ;ii=ii+1
   list_dpr(ii)=pawtab%shape_sigma  ;ii=ii+1
!Reals arrays (read from psp file)
   if (siz_coredens>0) then
     list_dpr(ii:ii+siz_coredens-1)=pawtab%coredens(1:siz_coredens)
     ii=ii+siz_coredens
   end if
   if (siz_coretau>0) then
     list_dpr(ii:ii+siz_coretau-1)=pawtab%coretau(1:siz_coretau)
     ii=ii+siz_coretau
   end if
   if (siz_dij0>0) then
     list_dpr(ii:ii+siz_dij0-1)=pawtab%dij0(1:siz_dij0)
     ii=ii+siz_dij0
   end if
   if (siz_fock>0) then
     list_dpr(ii:ii+siz_fock-1)=pawtab%ex_cvij(1:siz_fock)
     ii=ii+siz_fock
   end if
   if (siz_kij>0) then
     list_dpr(ii:ii+siz_kij-1)=pawtab%kij(1:siz_kij)
     ii=ii+siz_kij
   end if
   if (siz_phi>0) then
     list_dpr(ii:ii+siz_phi-1)=reshape(pawtab%phi,(/siz_phi/))
     ii=ii+siz_phi
   end if
   if (siz_nablaphi>0) then
     list_dpr(ii:ii+siz_nablaphi-1)=reshape(pawtab%nablaphi,(/siz_nablaphi/))
     ii=ii+siz_nablaphi
   end if
   if (siz_rhoij0>0) then
     list_dpr(ii:ii+siz_rhoij0-1)=pawtab%rhoij0(1:siz_rhoij0)
     ii=ii+siz_rhoij0
   end if
   if (siz_shape_alpha>0) then
     list_dpr(ii:ii+siz_shape_alpha-1)=reshape(pawtab%shape_alpha,(/siz_shape_alpha/))
     ii=ii+siz_shape_alpha
   end if
   if (siz_shape_q>0) then
     list_dpr(ii:ii+siz_shape_q-1)=reshape(pawtab%shape_q,(/siz_shape_q/))
     ii=ii+siz_shape_q
   end if
   if (siz_shapefunc>0) then
     list_dpr(ii:ii+siz_shapefunc-1)=reshape(pawtab%shapefunc,(/siz_shapefunc/))
     ii=ii+siz_shapefunc
   end if
   if (siz_tcoredens>0) then
     list_dpr(ii:ii+siz_tcoredens-1)=reshape(pawtab%tcoredens,(/siz_tcoredens/))
     ii=ii+siz_tcoredens
   end if
   if (siz_tcoretau>0) then
     list_dpr(ii:ii+siz_tcoretau-1)=reshape(pawtab%tcoretau,(/siz_tcoretau/))
     ii=ii+siz_tcoretau
   end if
   if (siz_tcorespl>0) then
     list_dpr(ii:ii+siz_tcorespl-1)=reshape(pawtab%tcorespl,(/siz_tcorespl/))
     ii=ii+siz_tcorespl
   end if
   if (siz_tcoretauspl>0) then
     list_dpr(ii:ii+siz_tcoretauspl-1)=reshape(pawtab%tcoretauspl,(/siz_tcoretauspl/))
     ii=ii+siz_tcoretauspl
   end if
   if (siz_tphi>0) then
     list_dpr(ii:ii+siz_tphi-1)=reshape(pawtab%tphi,(/siz_tphi/))
     ii=ii+siz_tphi
   end if
   if (siz_tnablaphi>0) then
     list_dpr(ii:ii+siz_tnablaphi-1)=reshape(pawtab%tnablaphi,(/siz_tnablaphi/))
     ii=ii+siz_tnablaphi
   end if
   if (siz_tproj>0) then
     list_dpr(ii:ii+siz_tproj-1)=reshape(pawtab%tproj,(/siz_tproj/))
     ii=ii+siz_tproj
   end if
   if (siz_tvalespl>0) then
     list_dpr(ii:ii+siz_tvalespl-1)=reshape(pawtab%tvalespl,(/siz_tvalespl/))
     ii=ii+siz_tvalespl
   end if
   if (siz_vhtnzc>0) then
     list_dpr(ii:ii+siz_vhtnzc-1)=pawtab%vhtnzc(1:siz_vhtnzc)
     ii=ii+siz_vhtnzc
   end if
   if (siz_vhnzc>0) then
     list_dpr(ii:ii+siz_vhnzc-1)=pawtab%vhnzc(1:siz_vhnzc)
     ii=ii+siz_vhnzc
   end if
   if (siz_vminushalf>0) then
     list_dpr(ii:ii+siz_vminushalf-1)=pawtab%vminushalf(1:siz_vminushalf)
     ii=ii+siz_vminushalf
   end if
!Reals in datastructures (read from psp file)
   if (siz_wvlpaw==1) then
     if (siz_wvl_parg>0) then
       list_dpr(ii:ii+siz_wvl_parg-1)=reshape(pawtab%wvl%parg,(/siz_wvl_parg/))
       ii=ii+siz_wvl_parg
     end if
     if (siz_wvl_pfac>0) then
       list_dpr(ii:ii+siz_wvl_pfac-1)=reshape(pawtab%wvl%pfac,(/siz_wvl_pfac/))
       ii=ii+siz_wvl_pfac
     end if
     if (siz_wvl_rholoc_rad>0) then
        list_dpr(ii:ii+siz_wvl_rholoc_rad-1)=pawtab%wvl%rholoc%rad(1:siz_wvl_rholoc_rad)
        ii=ii+siz_wvl_rholoc_rad
     end if
     if (siz_wvl_rholoc_d>0) then
        list_dpr(ii:ii+siz_wvl_rholoc_d-1)=reshape(pawtab%wvl%rholoc%d,(/siz_wvl_rholoc_d/))
        ii=ii+siz_wvl_rholoc_d
     end if
   end if

!Then the data initialized later
!...................................
   if (full_broadcast) then

!Reals
     list_dpr(ii)=pawtab%exchmix  ;ii=ii+1
     list_dpr(ii)=pawtab%f4of2_sla  ;ii=ii+1
     list_dpr(ii)=pawtab%f6of2_sla  ;ii=ii+1
     list_dpr(ii)=pawtab%jpawu  ;ii=ii+1
     list_dpr(ii)=pawtab%upawu  ;ii=ii+1
!Reals arrays
     if (siz_dltij>0) then
       list_dpr(ii:ii+siz_dltij-1)=pawtab%dltij(1:siz_dltij)
       ii=ii+siz_dltij
     end if
     if (siz_dshpfunc>0) then
       list_dpr(ii:ii+siz_dshpfunc-1)=reshape(pawtab%dshpfunc,(/siz_dshpfunc/))
       ii=ii+siz_dshpfunc
     end if
     if (siz_eijkl>0) then
       list_dpr(ii:ii+siz_eijkl-1)=reshape(pawtab%eijkl,(/siz_eijkl/))
       ii=ii+siz_eijkl
     end if
     if (siz_eijkl_sr>0) then
       list_dpr(ii:ii+siz_eijkl_sr-1)=reshape(pawtab%eijkl_sr,(/siz_eijkl_sr/))
       ii=ii+siz_eijkl_sr
     end if
     if (siz_euijkl>0) then
       list_dpr(ii:ii+siz_euijkl-1)=reshape(pawtab%euijkl,(/siz_euijkl/))
       ii=ii+siz_euijkl
     end if
     if (siz_euij_fll>0) then
       list_dpr(ii:ii+siz_euij_fll-1)=reshape(pawtab%euij_fll,(/siz_euij_fll/))
       ii=ii+siz_euij_fll
     end if
     if (siz_fk>0) then
       list_dpr(ii:ii+siz_fk-1)=reshape(pawtab%fk,(/siz_fk/))
       ii=ii+siz_fk
     end if
     if (siz_gammaij>0) then
       list_dpr(ii:ii+siz_gammaij-1)=pawtab%gammaij(1:siz_gammaij)
       ii=ii+siz_gammaij
     end if
     if (siz_gnorm>0) then
       list_dpr(ii:ii+siz_gnorm-1)=pawtab%gnorm(1:siz_gnorm)
       ii=ii+siz_gnorm
     end if
     if (siz_nabla_ij>0) then
       list_dpr(ii:ii+siz_nabla_ij-1)=reshape(pawtab%nabla_ij,(/siz_nabla_ij/))
       ii=ii+siz_nabla_ij
     end if
     if (siz_phiphj>0) then
       list_dpr(ii:ii+siz_phiphj-1)=reshape(pawtab%phiphj,(/siz_phiphj/))
       ii=ii+siz_phiphj
     end if
     if (siz_phiphjint>0) then
       list_dpr(ii:ii+siz_phiphjint-1)=pawtab%phiphjint(1:siz_phiphjint)
       ii=ii+siz_phiphjint
     end if
     if (siz_ph0phiint>0) then
       list_dpr(ii:ii+siz_ph0phiint-1)=pawtab%ph0phiint(1:siz_ph0phiint)
       ii=ii+siz_ph0phiint
     end if
     if (siz_qgrid_shp>0) then
       list_dpr(ii:ii+siz_qgrid_shp-1)=pawtab%qgrid_shp(1:siz_qgrid_shp)
       ii=ii+siz_qgrid_shp
     end if
     if (siz_qijl>0) then
       list_dpr(ii:ii+siz_qijl-1)=reshape(pawtab%qijl,(/siz_qijl/))
       ii=ii+siz_qijl
     end if
     if (siz_rad_for_spline>0) then
       list_dpr(ii:ii+siz_rad_for_spline-1)=pawtab%rad_for_spline(1:siz_rad_for_spline)
       ii=ii+siz_rad_for_spline
     end if
     if (siz_shapefncg>0) then
       list_dpr(ii:ii+siz_shapefncg-1)=reshape(pawtab%shapefncg,(/siz_shapefncg/))
       ii=ii+siz_shapefncg
     end if
     if (siz_sij>0) then
       list_dpr(ii:ii+siz_sij-1)=pawtab%sij(1:siz_sij)
       ii=ii+siz_sij
     end if
     if (siz_tphitphj>0) then
       list_dpr(ii:ii+siz_tphitphj-1)=reshape(pawtab%tphitphj,(/siz_tphitphj/))
       ii=ii+siz_tphitphj
     end if
     if (siz_vee>0) then
       list_dpr(ii:ii+siz_vee-1)=reshape(pawtab%vee,(/siz_vee/))
       ii=ii+siz_vee
     end if
     if (siz_vex>0) then
       list_dpr(ii:ii+siz_vex-1)=reshape(pawtab%vex,(/siz_vex/))
       ii=ii+siz_vex
     end if
     if (siz_zioneff>0) then
       list_dpr(ii:ii+siz_zioneff-1)=pawtab%zioneff(1:siz_zioneff)
       ii=ii+siz_zioneff
     end if

   end if ! full_broadcast
   ii=ii-1
   if (ii/=nn_dpr+nn_dpr_arr) then
     msg='the number of loaded reals is not correct!'
     ABI_BUG(msg)
   end if

 end if ! me=0

!Perfom the communication
!-------------------------------------------------------------------------

 call xmpi_bcast(list_dpr,0,comm_mpi,ierr)

!Fill the receiver from the buffer
!-------------------------------------------------------------------------
 if (me/=0) then
   ii=1

!First the data read from a psp file
!...................................

!Reals (read from psp file)
   pawtab%beta=list_dpr(ii)    ;ii=ii+1
   pawtab%dncdq0=list_dpr(ii)  ;ii=ii+1
   pawtab%d2ncdq0=list_dpr(ii)  ;ii=ii+1
   pawtab%dnvdq0=list_dpr(ii)  ;ii=ii+1
   pawtab%dtaucdq0=list_dpr(ii)  ;ii=ii+1
   pawtab%ex_cc=list_dpr(ii)  ;ii=ii+1
   pawtab%exccore=list_dpr(ii)  ;ii=ii+1
   pawtab%rpaw=list_dpr(ii)  ;ii=ii+1
   pawtab%rshp=list_dpr(ii)  ;ii=ii+1
   pawtab%rcore=list_dpr(ii)  ;ii=ii+1
   pawtab%rcoretau=list_dpr(ii)  ;ii=ii+1
   pawtab%shape_sigma=list_dpr(ii)  ;ii=ii+1
!Reals arrays (read from psp file)
   if (allocated(pawtab%coredens)) then
     LIBPAW_DEALLOCATE(pawtab%coredens)
   end if
   if (siz_coredens>0) then
     LIBPAW_ALLOCATE(pawtab%coredens,(pawtab%core_mesh_size))
     pawtab%coredens=list_dpr(ii:ii+pawtab%core_mesh_size-1)
     ii=ii+siz_coredens
   end if
   if (allocated(pawtab%coretau)) then
     LIBPAW_DEALLOCATE(pawtab%coretau)
   end if
   if (siz_coretau>0) then
     LIBPAW_ALLOCATE(pawtab%coretau,(pawtab%coretau_mesh_size))
     pawtab%coretau=list_dpr(ii:ii+pawtab%coretau_mesh_size-1)
     ii=ii+siz_coretau
   end if
   if (allocated(pawtab%dij0)) then
     LIBPAW_DEALLOCATE(pawtab%dij0)
   end if
   if (siz_dij0>0) then
     LIBPAW_ALLOCATE(pawtab%dij0,(pawtab%lmn2_size))
     pawtab%dij0=list_dpr(ii:ii+pawtab%lmn2_size-1)
     ii=ii+siz_dij0
   end if
   if (allocated(pawtab%ex_cvij)) then
     LIBPAW_DEALLOCATE(pawtab%ex_cvij)
   end if
   if (siz_fock>0) then
     LIBPAW_ALLOCATE(pawtab%ex_cvij,(pawtab%lmn2_size))
     pawtab%ex_cvij=list_dpr(ii:ii+pawtab%lmn2_size-1)
     ii=ii+siz_fock
   end if
   if (allocated(pawtab%kij)) then
     LIBPAW_DEALLOCATE(pawtab%kij)
   end if
   if (siz_kij>0) then
     LIBPAW_ALLOCATE(pawtab%kij,(pawtab%lmn2_size))
     pawtab%kij=list_dpr(ii:ii+pawtab%lmn2_size-1)
     ii=ii+siz_kij
   end if
   if (allocated(pawtab%phi)) then
     LIBPAW_DEALLOCATE(pawtab%phi)
   end if
   if (siz_phi>0) then
     LIBPAW_ALLOCATE(pawtab%phi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
     pawtab%phi=reshape(list_dpr(ii:ii+siz_phi-1),(/pawtab%partialwave_mesh_size,pawtab%basis_size/))
     ii=ii+siz_phi
   end if
   if (allocated(pawtab%nablaphi)) then
     LIBPAW_DEALLOCATE(pawtab%nablaphi)
   end if
   if (siz_nablaphi>0) then
     LIBPAW_ALLOCATE(pawtab%nablaphi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
     pawtab%nablaphi=reshape(list_dpr(ii:ii+siz_nablaphi-1),(/pawtab%partialwave_mesh_size,pawtab%basis_size/))
     ii=ii+siz_nablaphi
   end if
   if (allocated(pawtab%rhoij0)) then
     LIBPAW_DEALLOCATE(pawtab%rhoij0)
   end if
   if (siz_rhoij0>0) then
     LIBPAW_ALLOCATE(pawtab%rhoij0,(pawtab%lmn2_size))
     pawtab%rhoij0=list_dpr(ii:ii+pawtab%lmn2_size-1)
     ii=ii+siz_rhoij0
   end if
   if (allocated(pawtab%shape_alpha)) then
     LIBPAW_DEALLOCATE(pawtab%shape_alpha)
   end if
   if (siz_shape_alpha>0) then
     LIBPAW_ALLOCATE(pawtab%shape_alpha,(2,pawtab%l_size))
     pawtab%shape_alpha=reshape(list_dpr(ii:ii+siz_shape_alpha-1),(/2,pawtab%l_size/))
     ii=ii+siz_shape_alpha
   end if
   if (allocated(pawtab%shape_q)) then
     LIBPAW_DEALLOCATE(pawtab%shape_q)
   end if
   if (siz_shape_q>0) then
     LIBPAW_ALLOCATE(pawtab%shape_q,(2,pawtab%l_size))
     pawtab%shape_q=reshape(list_dpr(ii:ii+siz_shape_q-1),(/2,pawtab%l_size/))
     ii=ii+siz_shape_q
   end if
   if (allocated(pawtab%shapefunc)) then
     LIBPAW_DEALLOCATE(pawtab%shapefunc)
   end if
   if (siz_shapefunc>0) then
     LIBPAW_ALLOCATE(pawtab%shapefunc,(pawtab%mesh_size,pawtab%l_size))
     pawtab%shapefunc=reshape(list_dpr(ii:ii+siz_shapefunc-1),(/pawtab%mesh_size,pawtab%l_size/))
     ii=ii+siz_shapefunc
   end if
   if (allocated(pawtab%tcoredens)) then
     LIBPAW_DEALLOCATE(pawtab%tcoredens)
   end if
   if (siz_tcoredens>0) then
     sz2=siz_tcoredens/pawtab%core_mesh_size
     LIBPAW_ALLOCATE(pawtab%tcoredens,(pawtab%core_mesh_size,sz2))
     pawtab%tcoredens=reshape(list_dpr(ii:ii+siz_tcoredens-1),(/pawtab%core_mesh_size,sz2/))
     ii=ii+siz_tcoredens
   end if
   if (allocated(pawtab%tcoretau)) then
     LIBPAW_DEALLOCATE(pawtab%tcoretau)
   end if
   if (siz_tcoretau>0) then
     LIBPAW_ALLOCATE(pawtab%tcoretau,(pawtab%coretau_mesh_size))
     pawtab%tcoretau=list_dpr(ii:ii+siz_tcoretau-1)
     ii=ii+siz_tcoretau
   end if
   if (allocated(pawtab%tcorespl)) then
     LIBPAW_DEALLOCATE(pawtab%tcorespl)
   end if
   if (siz_tcorespl>0) then
     LIBPAW_ALLOCATE(pawtab%tcorespl,(pawtab%mqgrid,2))
     pawtab%tcorespl=reshape(list_dpr(ii:ii+siz_tcorespl-1),(/pawtab%mqgrid,2/))
     ii=ii+siz_tcorespl
   end if
   if (allocated(pawtab%tcoretauspl)) then
     LIBPAW_DEALLOCATE(pawtab%tcoretauspl)
   end if
   if (siz_tcoretauspl>0) then
     LIBPAW_ALLOCATE(pawtab%tcoretauspl,(pawtab%mqgrid,2))
     pawtab%tcoretauspl=reshape(list_dpr(ii:ii+siz_tcoretauspl-1),(/pawtab%mqgrid,2/))
     ii=ii+siz_tcoretauspl
   else
     LIBPAW_ALLOCATE(pawtab%tcoretauspl,(pawtab%mqgrid,0))
   end if
   if (allocated(pawtab%tphi)) then
     LIBPAW_DEALLOCATE(pawtab%tphi)
   end if
   if (siz_tphi>0) then
     LIBPAW_ALLOCATE(pawtab%tphi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
     pawtab%tphi=reshape(list_dpr(ii:ii+siz_tphi-1),(/pawtab%partialwave_mesh_size,pawtab%basis_size/))
     ii=ii+siz_tphi
   end if
  if (allocated(pawtab%tnablaphi)) then
     LIBPAW_DEALLOCATE(pawtab%tnablaphi)
   end if
   if (siz_tnablaphi>0) then
     LIBPAW_ALLOCATE(pawtab%tnablaphi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
     pawtab%tphi=reshape(list_dpr(ii:ii+siz_tnablaphi-1),(/pawtab%partialwave_mesh_size,pawtab%basis_size/))
     ii=ii+siz_tnablaphi
   end if
   if (allocated(pawtab%tproj)) then
     LIBPAW_DEALLOCATE(pawtab%tproj)
   end if
   if (siz_tproj>0) then
     sz1=siz_tproj/pawtab%basis_size
     LIBPAW_ALLOCATE(pawtab%tproj,(sz1,pawtab%basis_size))
     pawtab%tproj=reshape(list_dpr(ii:ii+siz_tproj-1),(/sz1,pawtab%basis_size/))
     ii=ii+siz_tproj
   end if
   if (allocated(pawtab%tvalespl)) then
     LIBPAW_DEALLOCATE(pawtab%tvalespl)
   end if
   if (siz_tvalespl>0) then
     sz1=siz_tvalespl/2
     LIBPAW_ALLOCATE(pawtab%tvalespl,(sz1,2))
     pawtab%tvalespl=reshape(list_dpr(ii:ii+siz_tvalespl-1),(/sz1,2/))
     ii=ii+siz_tvalespl
   end if
   if (allocated(pawtab%vhtnzc)) then
     LIBPAW_DEALLOCATE(pawtab%vhtnzc)
   end if
   if (siz_vhtnzc>0) then
     LIBPAW_ALLOCATE(pawtab%vhtnzc,(pawtab%mesh_size))
     pawtab%vhtnzc=list_dpr(ii:ii+pawtab%mesh_size-1)
     ii=ii+siz_vhtnzc
   end if
   if (allocated(pawtab%vhnzc)) then
     LIBPAW_DEALLOCATE(pawtab%vhnzc)
   end if
   if (siz_vhnzc>0) then
     LIBPAW_ALLOCATE(pawtab%vhnzc,(pawtab%mesh_size))
     pawtab%vhnzc=list_dpr(ii:ii+pawtab%mesh_size-1)
     ii=ii+siz_vhnzc
   end if
   if (allocated(pawtab%vminushalf)) then
     LIBPAW_DEALLOCATE(pawtab%vminushalf)
   end if
   if (siz_vminushalf>0) then
     LIBPAW_ALLOCATE(pawtab%vminushalf,(pawtab%mesh_size))
     pawtab%vminushalf=list_dpr(ii:ii+pawtab%mesh_size-1)
     ii=ii+siz_vminushalf
   end if
!Reals in datastructures (read from psp file)
   if (siz_wvlpaw==1) then
     if (allocated(pawtab%wvl%parg)) then
       LIBPAW_DEALLOCATE(pawtab%wvl%parg)
     end if
     if (siz_wvl_parg>0) then
       LIBPAW_ALLOCATE(pawtab%wvl%parg,(2,pawtab%wvl%ptotgau))
       pawtab%wvl%parg=reshape(list_dpr(ii:ii+siz_wvl_parg-1),(/2,pawtab%wvl%ptotgau/))
       ii=ii+siz_wvl_parg
     end if
     if (allocated(pawtab%wvl%pfac)) then
       LIBPAW_DEALLOCATE(pawtab%wvl%pfac)
     end if
     if (siz_wvl_pfac>0) then
       LIBPAW_ALLOCATE(pawtab%wvl%pfac,(2,pawtab%wvl%ptotgau))
       pawtab%wvl%pfac=reshape(list_dpr(ii:ii+siz_wvl_pfac-1),(/2,pawtab%wvl%ptotgau/))
       ii=ii+siz_wvl_pfac
     end if
     if (allocated(pawtab%wvl%rholoc%rad)) then
        LIBPAW_DEALLOCATE(pawtab%wvl%rholoc%rad)
     end if
     if (siz_wvl_rholoc_rad>0) then
        sz1=pawtab%wvl%rholoc%msz
        LIBPAW_ALLOCATE(pawtab%wvl%rholoc%rad,(sz1))
        pawtab%wvl%rholoc%rad=list_dpr(ii:ii+sz1-1)
        ii=ii+siz_wvl_rholoc_rad
     end if
     if (allocated(pawtab%wvl%rholoc%d)) then
        LIBPAW_DEALLOCATE(pawtab%wvl%rholoc%d)
     end if
     if (siz_wvl_rholoc_d>0) then
        sz1=pawtab%wvl%rholoc%msz
        LIBPAW_ALLOCATE(pawtab%wvl%rholoc%d,(sz1,4))
        pawtab%wvl%rholoc%d=reshape(list_dpr(ii:ii+siz_wvl_rholoc_d-1),(/sz1,4/))
        ii=ii+siz_wvl_rholoc_d
     end if
   end if

!Then the data initialized later
!...................................
   if (full_broadcast) then

!Reals
     pawtab%exchmix=list_dpr(ii)  ;ii=ii+1
     pawtab%f4of2_sla=list_dpr(ii)  ;ii=ii+1
     pawtab%f6of2_sla=list_dpr(ii)  ;ii=ii+1
     pawtab%jpawu=list_dpr(ii)  ;ii=ii+1
     pawtab%upawu=list_dpr(ii)  ;ii=ii+1
!Reals arrays
     if (allocated(pawtab%dltij)) then
       LIBPAW_DEALLOCATE(pawtab%dltij)
     end if
     if (siz_dltij>0) then
       LIBPAW_ALLOCATE(pawtab%dltij,(pawtab%lmn2_size))
       pawtab%dltij=list_dpr(ii:ii+pawtab%lmn2_size-1)
       ii=ii+siz_dltij
     end if
     if (allocated(pawtab%dshpfunc)) then
       LIBPAW_DEALLOCATE(pawtab%dshpfunc)
     end if
     if (siz_dshpfunc>0) then
       LIBPAW_ALLOCATE(pawtab%dshpfunc,(pawtab%mesh_size,pawtab%l_size,4))
       pawtab%dshpfunc=reshape(list_dpr(ii:ii+siz_dshpfunc-1),(/pawtab%mesh_size,pawtab%l_size,4/))
       ii=ii+siz_dshpfunc
     end if
     if (allocated(pawtab%eijkl)) then
       LIBPAW_DEALLOCATE(pawtab%eijkl)
     end if
     if (siz_eijkl>0) then
       LIBPAW_ALLOCATE(pawtab%eijkl,(pawtab%lmn2_size,pawtab%lmn2_size))
       pawtab%eijkl=reshape(list_dpr(ii:ii+siz_eijkl-1),(/pawtab%lmn2_size,pawtab%lmn2_size/))
       ii=ii+siz_eijkl
     end if
     if (allocated(pawtab%eijkl_sr)) then
       LIBPAW_DEALLOCATE(pawtab%eijkl_sr)
     end if
     if (siz_eijkl_sr>0) then
       LIBPAW_ALLOCATE(pawtab%eijkl_sr,(pawtab%lmn2_size,pawtab%lmn2_size))
       pawtab%eijkl_sr=reshape(list_dpr(ii:ii+siz_eijkl_sr-1),(/pawtab%lmn2_size,pawtab%lmn2_size/))
       ii=ii+siz_eijkl_sr
     end if
     if (allocated(pawtab%euijkl)) then
       LIBPAW_DEALLOCATE(pawtab%euijkl)
     end if
     if (siz_euijkl>0) then
       LIBPAW_ALLOCATE(pawtab%euijkl,(2,2,pawtab%lmn_size,pawtab%lmn_size,pawtab%lmn_size,pawtab%lmn_size))
       pawtab%euijkl=reshape(list_dpr(ii:ii+siz_euijkl-1),(/2,2,pawtab%lmn_size,pawtab%lmn_size,pawtab%lmn_size,pawtab%lmn_size/))
       ii=ii+siz_euijkl
     end if
     if (allocated(pawtab%euij_fll)) then
       LIBPAW_DEALLOCATE(pawtab%euij_fll)
     end if
     if (siz_euij_fll>0) then
       LIBPAW_ALLOCATE(pawtab%euij_fll,(pawtab%lmn2_size))
       pawtab%euij_fll=reshape(list_dpr(ii:ii+siz_euij_fll-1),(/pawtab%lmn2_size/))
       ii=ii+siz_euij_fll
     end if
     if (allocated(pawtab%fk)) then
       LIBPAW_DEALLOCATE(pawtab%fk)
     end if
     if (siz_fk>0) then
       LIBPAW_ALLOCATE(pawtab%fk,(6,4))
       pawtab%fk=reshape(list_dpr(ii:ii+siz_fk-1),(/6,4/))
       ii=ii+siz_fk
     end if
     if (allocated(pawtab%gammaij)) then
       LIBPAW_DEALLOCATE(pawtab%gammaij)
     end if
     if (siz_gammaij>0) then
       LIBPAW_ALLOCATE(pawtab%gammaij,(pawtab%l_size))
       pawtab%gammaij=list_dpr(ii:ii+pawtab%l_size-1)
       ii=ii+siz_gammaij
     end if
     if (allocated(pawtab%gnorm)) then
       LIBPAW_DEALLOCATE(pawtab%gnorm)
     end if
     if (siz_gnorm>0) then
       LIBPAW_ALLOCATE(pawtab%gnorm,(pawtab%l_size))
       pawtab%gnorm=list_dpr(ii:ii+pawtab%l_size-1)
       ii=ii+siz_gnorm
     end if
     if (allocated(pawtab%nabla_ij)) then
       LIBPAW_DEALLOCATE(pawtab%nabla_ij)
     end if
     if (siz_nabla_ij>0) then
       LIBPAW_ALLOCATE(pawtab%nabla_ij,(3,pawtab%lmn_size,pawtab%lmn_size))
       pawtab%nabla_ij=reshape(list_dpr(ii:ii+siz_nabla_ij-1),(/3,pawtab%lmn_size,pawtab%lmn_size/))
       ii=ii+siz_nabla_ij
     end if
     if (allocated(pawtab%phiphj)) then
       LIBPAW_DEALLOCATE(pawtab%phiphj)
     end if
     if (siz_phiphj>0) then
       LIBPAW_ALLOCATE(pawtab%phiphj,(pawtab%mesh_size,pawtab%ij_size))
       pawtab%phiphj=reshape(list_dpr(ii:ii+siz_phiphj-1),(/pawtab%mesh_size,pawtab%ij_size/))
       ii=ii+siz_phiphj
     end if
     if (allocated(pawtab%phiphjint)) then
       LIBPAW_DEALLOCATE(pawtab%phiphjint)
     end if
     if (siz_phiphjint>0) then
       LIBPAW_ALLOCATE(pawtab%phiphjint,(pawtab%ij_proj))
       pawtab%phiphjint=list_dpr(ii:ii+pawtab%ij_proj-1)
       ii=ii+siz_phiphjint
     end if
     if (allocated(pawtab%ph0phiint)) then
       LIBPAW_DEALLOCATE(pawtab%ph0phiint)
     end if
     if (siz_ph0phiint>0) then
       LIBPAW_ALLOCATE(pawtab%ph0phiint,(pawtab%ij_proj))
       pawtab%ph0phiint=list_dpr(ii:ii+pawtab%ij_proj-1)
       ii=ii+siz_ph0phiint
     end if
     if (allocated(pawtab%qgrid_shp)) then
       LIBPAW_DEALLOCATE(pawtab%qgrid_shp)
     end if
     if (siz_qgrid_shp>0) then
       LIBPAW_ALLOCATE(pawtab%qgrid_shp,(pawtab%mqgrid_shp))
       pawtab%qgrid_shp=list_dpr(ii:ii+pawtab%mqgrid_shp-1)
       ii=ii+siz_qgrid_shp
     end if
     if (allocated(pawtab%qijl)) then
       LIBPAW_DEALLOCATE(pawtab%qijl)
     end if
     if (siz_qijl>0) then
       LIBPAW_ALLOCATE(pawtab%qijl,(pawtab%l_size**2,pawtab%lmn2_size))
       pawtab%qijl=reshape(list_dpr(ii:ii+siz_qijl-1),(/pawtab%l_size**2,pawtab%lmn2_size/))
       ii=ii+siz_qijl
     end if
     if (allocated(pawtab%rad_for_spline)) then
       LIBPAW_DEALLOCATE(pawtab%rad_for_spline)
     end if
     if (siz_rad_for_spline>0) then
       LIBPAW_ALLOCATE(pawtab%rad_for_spline,(pawtab%mesh_size))
       pawtab%rad_for_spline=list_dpr(ii:ii+pawtab%mesh_size-1)
       ii=ii+siz_rad_for_spline
     end if
     if (allocated(pawtab%shapefncg)) then
       LIBPAW_DEALLOCATE(pawtab%shapefncg)
     end if
     if (siz_shapefncg>0) then
       LIBPAW_ALLOCATE(pawtab%shapefncg,(pawtab%mqgrid_shp,2,pawtab%l_size))
       pawtab%shapefncg=reshape(list_dpr(ii:ii+siz_shapefncg-1),(/pawtab%mqgrid_shp,2,pawtab%l_size/))
       ii=ii+siz_shapefncg
     end if
     if (allocated(pawtab%sij)) then
       LIBPAW_DEALLOCATE(pawtab%sij)
     end if
     if (siz_sij>0) then
       LIBPAW_ALLOCATE(pawtab%sij,(pawtab%lmn2_size))
       pawtab%sij=list_dpr(ii:ii+pawtab%lmn2_size-1)
       ii=ii+siz_sij
     end if
     if (allocated(pawtab%tphitphj)) then
       LIBPAW_DEALLOCATE(pawtab%tphitphj)
     end if
     if (siz_tphitphj>0) then
       LIBPAW_ALLOCATE(pawtab%tphitphj,(pawtab%mesh_size,pawtab%ij_size))
       pawtab%tphitphj=reshape(list_dpr(ii:ii+siz_tphitphj-1),(/pawtab%mesh_size,pawtab%ij_size/))
       ii=ii+siz_tphitphj
     end if
     if (allocated(pawtab%vee)) then
       LIBPAW_DEALLOCATE(pawtab%vee)
     end if
     if (siz_vee>0) then
       sz1=2*pawtab%lpawu+1
       LIBPAW_ALLOCATE(pawtab%vee,(sz1,sz1,sz1,sz1))
       pawtab%vee=reshape(list_dpr(ii:ii+siz_vee-1),(/sz1,sz1,sz1,sz1/))
       ii=ii+siz_vee
     end if
     if (allocated(pawtab%vex)) then
       LIBPAW_DEALLOCATE(pawtab%vex)
     end if
     if (siz_vex>0) then
       sz1=2*pawtab%lexexch+1
       LIBPAW_ALLOCATE(pawtab%vex,(sz1,sz1,sz1,sz1,4))
       pawtab%vex=reshape(list_dpr(ii:ii+siz_vex-1),(/sz1,sz1,sz1,sz1,4/))
       ii=ii+siz_vex
     end if
     if (allocated(pawtab%zioneff)) then
       LIBPAW_DEALLOCATE(pawtab%zioneff)
     end if
     if (siz_zioneff>0) then
       LIBPAW_ALLOCATE(pawtab%zioneff,(pawtab%ij_proj))
       pawtab%zioneff=list_dpr(ii:ii+pawtab%ij_proj-1)
       ii=ii+siz_zioneff
     end if

   end if ! full_broadcast
   ii=ii-1

   if (ii/=nn_dpr+nn_dpr_arr) then
     msg='the number of broadcasted reals is not correct!'
     ABI_BUG(msg)
   end if

 end if ! me/=0
 LIBPAW_DEALLOCATE(list_dpr)

end subroutine pawtab_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_allocate
!! NAME
!!  wvlpaw_allocate
!!
!! FUNCTION
!!  Allocate (if necessary) and nullify content of a wvlpaw pointer
!!
!! SIDE EFFECTS
!!  wvlpaw<type(wvlpaw_type)>=datastructure to be allocated.
!!
!! PARENTS
!!      m_pawpsp,m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_allocate(wvlpaw)

!Arguments ------------------------------------
 type(wvlpaw_type),pointer :: wvlpaw

! *************************************************************************

 !@wvlpaw_type

 if (.not.associated(wvlpaw)) then
   LIBPAW_DATATYPE_ALLOCATE(wvlpaw,)
   call wvlpaw_nullify(wvlpaw)
 end if

 wvlpaw%npspcode_init_guess=10

end subroutine wvlpaw_allocate
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_free
!! NAME
!!  wvlpaw_free
!!
!! FUNCTION
!!  Deallocate arrays and nullify flags in a wvlpaw structure
!!
!! SIDE EFFECTS
!!  wvlpaw<type(wvlpaw_type)>=datastructure to be destroyed.
!!  All allocated arrays are deallocated.
!!
!! PARENTS
!!      m_pawpsp,m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_free(wvlpaw)

!Arguments ------------------------------------
 type(wvlpaw_type),pointer :: wvlpaw

! *************************************************************************

 !@wvlpaw_type

 if (.not.associated(wvlpaw)) return

 if(allocated(wvlpaw%pngau)) then
   LIBPAW_DEALLOCATE(wvlpaw%pngau)
 end if
 if(allocated(wvlpaw%parg)) then
   LIBPAW_DEALLOCATE(wvlpaw%parg)
 end if
 if(allocated(wvlpaw%pfac)) then
   LIBPAW_DEALLOCATE(wvlpaw%pfac)
 end if

 wvlpaw%npspcode_init_guess=0
 wvlpaw%ptotgau=0

 call wvlpaw_rholoc_free(wvlpaw%rholoc)

 LIBPAW_DATATYPE_DEALLOCATE(wvlpaw)

end subroutine wvlpaw_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_nullify
!! NAME
!!  wvlpaw_nullify
!!
!! FUNCTION
!!  Nullify flags in a wvlpaw structure
!!
!! SIDE EFFECTS
!!  wvlpaw=datastructure to be nullified
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_nullify(wvlpaw)

!Arguments ------------------------------------
 type(wvlpaw_type),pointer :: wvlpaw

! *************************************************************************

 !@wvlpaw_type
 if (.not.associated(wvlpaw)) return

 wvlpaw%npspcode_init_guess=0
 wvlpaw%ptotgau=0

 call wvlpaw_rholoc_nullify(wvlpaw%rholoc)

end subroutine wvlpaw_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_rholoc_free
!! NAME
!!  wvlpaw_rholoc_free
!!
!! FUNCTION
!!  Deallocate arrays and nullify flags in a wvlpaw%rholoc structure
!!
!! SIDE EFFECTS
!!  wvlpaw_rholoc<type(wvlpaw_rholoc_type)>=datastructure to be destroyed.
!!  All allocated arrays are deallocated.
!!
!! PARENTS
!!      m_pawpsp,m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_rholoc_free(wvlpaw_rholoc)

!Arguments ------------------------------------
 type(wvlpaw_rholoc_type),intent(inout) :: wvlpaw_rholoc

! *************************************************************************

 !@wvlpaw_rholoc_type

 if(allocated(wvlpaw_rholoc%d)) then
   LIBPAW_DEALLOCATE(wvlpaw_rholoc%d)
 end if
 if(allocated(wvlpaw_rholoc%rad)) then
   LIBPAW_DEALLOCATE(wvlpaw_rholoc%rad)
 end if

 wvlpaw_rholoc%msz=0

end subroutine wvlpaw_rholoc_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_rholoc_nullify
!! NAME
!!  wvlpaw_rholoc_nullify
!!
!! FUNCTION
!!  Nullify flags in a wvlpaw%rholoc structure
!!
!! SIDE EFFECTS
!!  wvlpaw_rholoc<type(wvlpaw_rholoc_type)>=datastructure to be nullified.
!!
!! PARENTS
!!      m_pawpsp,m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_rholoc_nullify(wvlpaw_rholoc)

!Arguments ------------------------------------
 type(wvlpaw_rholoc_type),intent(inout) :: wvlpaw_rholoc

! *************************************************************************

 !@wvlpaw_rholoc_type

 wvlpaw_rholoc%msz=0

end subroutine wvlpaw_rholoc_nullify
!!***

!----------------------------------------------------------------------

END MODULE m_pawtab
!!***
