!!****m* ABINIT/m_hamiltonian
!! NAME
!! m_hamiltonian
!!
!! FUNCTION
!!  This module provides the definition of the gs_hamiltonian_type and rf_hamiltonian_type
!!  datastructures used in the "getghc" and "getgh1c" routines to apply the Hamiltonian (or
!!  its derivative) on a wavefunction.
!!  Methods to initialize or destroy the objects are defined here.
!!
!! TODO
!!  All array pointers in H datatypes should be declared as contiguous for efficient reasons
!!  (well, here performance is critical)
!!  Client code should make sure they always point contiguous targets.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MG, MT)
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

module m_hamiltonian

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_datatypes,      only : pseudopotential_type
 use m_copy,              only : addr_copy
 use m_geometry,          only : metric
 use m_pawtab,            only : pawtab_type
 use m_fftcore,           only : sphereboundary
 use m_pawcprj,           only : pawcprj_getdim
 use m_paw_ij,            only : paw_ij_type
 use m_paral_atom,        only : get_my_atmtab, free_my_atmtab
 use m_electronpositron,  only : electronpositron_type, electronpositron_calctype
 use m_kg,                only : ph1d3d, getph
 use m_fock,              only : fock_common_type, fock_BZ_type, fock_ACE_type, fock_type

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

#if defined HAVE_FC_ISO_C_BINDING
 use iso_c_binding, only : c_ptr,c_loc,c_f_pointer
#endif

 implicit none

 private

 public ::  pawdij2ekb
 public ::  pawdij2e1kb

 ! These constants select how H is applied in reciprocal space
 integer,parameter,public :: KPRIME_H_K=1, K_H_KPRIME=2, K_H_K=3, KPRIME_H_KPRIME=4
!!***

!----------------------------------------------------------------------

!!****t* m_hamiltonian/gs_hamiltonian_type
!! NAME
!! gs_hamiltonian_type
!!
!! FUNCTION
!! This datastructure contains the information about one Hamiltonian,
!! needed in the "getghc" routine, that apply the Hamiltonian on a wavefunction.
!! The Hamiltonian is expressed in reciprocal space:
!!
!!       H_k^prime,k = exp(-i.k^prime.r^prime) H exp(i.k.r)
!!
!! In most cases k = k^prime and the k^prime objects are simply pointers to k objects.
!!
!! SOURCE

 type,public :: gs_hamiltonian_type

! ===== Integer scalars

  integer :: dimekb1
   ! First dimension of Ekb
   ! Same as psps%dimekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb1=lnmax
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom
   !                     dimekb1=cplex_dij*lmnmax*(lmnmax+1)/2

  integer :: dimekb2
   ! Second dimension of Ekb
   ! ->Norm conserving psps: dimekb2=ntypat
   ! ->PAW                 : dimekb2=natom

  integer :: dimekbq
   ! Fourth dimension of Ekb
   ! 2 if Ekb factors contain a exp(-iqR) phase, 1 otherwise

  integer :: istwf_k
   ! option parameter that describes the storage of wfs at k

  integer :: istwf_kp
   ! option parameter that describes the storage of wfs at k^prime

  integer :: lmnmax
   ! Maximum number of different l,m,n components over all types of psps.
   ! same as dtset%lmnmax

  integer :: matblk
   ! dimension of the array ph3d

  integer :: mgfft
   ! maximum size for 1D FFTs (same as dtset%mgfft)

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"
   ! same as psps%mpsang

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! For NC pseudopotentials only:
   !   when mpspso=1, mpssoang=mpsang
   !   when mpspso=2, mpssoang=2*mpsang-1
   ! For PAW: same as mpsang
   ! same as psps%mpssoang

  integer :: natom
   ! The number of atoms for this dataset; same as dtset%natom

  integer :: nfft
   ! number of FFT grid points same as dtset%nfft

  integer :: npw_k
   ! number of plane waves at k
   ! In case of band-FFT parallelism, npw_k is the number of plane waves
   ! processed by current proc

  integer :: npw_fft_k
   ! number of plane waves at k used to apply Hamiltonian when band-FFT
   ! parallelism is activated (i.e. data are distributed in the "FFT" configuration)

  integer :: npw_kp
   ! number of plane waves at k^prime
   ! In case of band-FFT parallelism, npw_kp is the number of plane waves
   ! processed by current proc

  integer :: npw_fft_kp
   ! number of plane waves at k^prime used to apply Hamiltonian when band-FFT
   ! parallelism is activated (i.e. data are distributed in the "FFT" configuration)

  integer :: nspinor
   ! Number of spinorial components

  integer :: nsppol
   ! Total number of spin components (1=non-polarized, 2=polarized)

  integer :: ntypat
   ! Number of types of pseudopotentials same as dtset%ntypat

  integer :: nvloc
   ! Number of components of vloc
   ! usually, nvloc=1, except in the non-collinear magnetism case, where nvloc=4

  integer :: n4,n5,n6
   ! same as ngfft(4:6)

  integer :: use_gpu_cuda
  ! governs wheter we do the hamiltonian calculation on gpu or not

  integer :: usecprj
   ! usecprj= 1 if cprj projected WF are stored in memory
   !        = 0 if they are to be computed on the fly

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! ===== Integer arrays

  integer, allocatable :: atindx(:)
   ! atindx(natom)
   ! index table for atoms (see gstate.f)

  integer, allocatable :: atindx1(:)
   ! atindx1(natom)
   ! index table for atoms, inverse of atindx (see gstate.f)

  integer, allocatable :: dimcprj(:)
   ! dimcprj(natom*usepaw)=dimensions of array cprj
   ! dimcprj(ia)=cprj(ia,:)%nlmn
   ! atoms are ordered by atom-type

  integer, allocatable :: gbound_k(:,:)
   ! gbound_k(2*mgfft+8,2)
   ! G sphere boundary, for each plane wave at k

  integer, allocatable :: indlmn(:,:,:)
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

  integer, allocatable :: nattyp(:)
   ! nattyp(ntypat)
   ! # of atoms of each type

  integer :: ngfft(18)
   ! ngfft(1:3)=integer fft box dimensions
   ! ngfft(4:6)=integer fft box dimensions, might be augmented for CPU speed
   ! ngfft(7)=fftalg
   ! ngfft(8)=fftalg

  integer :: nloalg(3)
   ! governs the choice of the algorithm for non-local operator same as dtset%nloalg

  integer, allocatable :: pspso(:)
   ! pspso(ntypat)
   ! For each type of psp, 1 if no spin-orbit component is taken
   ! into account, 2 if a spin-orbit component is used
   ! Revelant for NC-psps and PAW

  integer, allocatable :: typat(:)
   ! typat(natom)
   ! type of each atom

  ! integer, allocatable :: indpw_k(:,:)
   ! indpw_k(4,npw_fft_k)
   ! array which gives fft box index for given basis sphere

! Integer pointers

  integer, ABI_CONTIGUOUS pointer :: gbound_kp(:,:) => null()
   ! gbound_kp(2*mgfft+8,2)
   ! G sphere boundary, for each plane wave at k^prime

  integer, pointer :: kg_k(:,:) => null()
   ! kg_k(3,npw_fft_k)
   ! G vector coordinates with respect to reciprocal lattice translations
   ! at k

  integer, pointer :: kg_kp(:,:) => null()
   ! kg_kp(3,npw_fft_kp)
   ! G vector coordinates with respect to reciprocal lattice translations
   ! at k^prime

! ===== Real scalars

  real(dp) :: ucvol
   ! unit cell volume (Bohr**3)

! ===== Real arrays

  real(dp), allocatable :: ekb_spin(:,:,:,:,:)
   ! ekb_spin(dimekb1,dimekb2,nspinor**2,dimekbq,my_nsppol)
   ! Contains the values of ekb array for all spins treated by current process
   ! See ekb description ; ekb is pointer to ekb_spin(:,:,:,:,my_isppol)

  real(dp), allocatable :: sij(:,:)
   ! sij(dimekb1,ntypat*usepaw) = overlap matrix for paw calculation

  real(dp) :: gmet(3,3)
   ! reciprocal space metric tensor in Bohr**-2

  real(dp) :: gprimd(3,3)
   ! dimensional reciprocal space primitive translations (Bohr^-1)

  real(dp) :: kpt_k(3)
   ! dimensionless k point coordinates wrt reciprocal lattice vectors

  real(dp) :: kpt_kp(3)
   ! dimensionless k^prime point coordinates wrt reciprocal lattice vectors

  real(dp), allocatable :: nucdipmom(:,:)
   ! nucdipmom(3,natom)
   ! nuclear dipole moments at each atomic position

  real(dp), allocatable :: ph1d(:,:)
   ! ph1d(2,3*(2*mgfft+1)*natom)
   ! 1-dim phase arrays for structure factor (see getph.f).

  real(dp), allocatable :: phkxred(:,:)
   ! phkxred(2,natom)
   ! phase factors exp(2 pi k.xred) at k

! ===== Complex arrays

  complex(dpc), allocatable :: nucdipmom_k(:)
   ! nucdipmom_k(npw_k*(npw_k+1)/2)
   ! nuclear dipole moment Hamiltonian in reciprocal space, stored as
   ! lower triangular part of Hermitian matrix

! ===== Real pointers

  real(dp), ABI_CONTIGUOUS pointer :: ekb(:,:,:,:) => null()
   ! ekb(dimekb1,dimekb2,nspinor**2,dimekbq)
   !  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
   !          for number of basis functions (l,n) (lnmax)
   !          and number of atom types (ntypat)
   !          dimekb1=lnmax ; dimekb2=ntypat ; dimekbq=1
   !  ->PAW : (Real, symmetric) Frozen part of Dij coefficients
   !                            to connect projectors
   !          for number of basis functions (l,m,n) (lmnmax)
   !          and number of atom (natom)
   !          dimekb1=lmnmax*(lmnmax+1)/2 ; dimekb2=natom ; dimekbq=1
   ! ekb is spin dependent in the case of PAW calculations.
   ! For each spin component, ekb points to ekb_spin(:,:,:,:,my_isppol)
   ! dimekbq=2 if Ekb factors contain a exp(-iqR) phase, dimekbq=1 otherwise
   ! About the non-local factors symmetry:
   !   - The lower triangular part of the Dij matrix can be deduced from the upper one
   !     with the following relation: D^s2s1_ji = (D^s1s2_ij)^*
   !     where s1,s2 are spinor components

  real(dp), pointer :: ffnl_k(:,:,:,:) => null()
   ! ffnl_k(npw_fft_k,2,dimffnl_k,ntypat)
   ! nonlocal form factors at k

  real(dp), pointer :: ffnl_kp(:,:,:,:) => null()
   ! ffnl_kp(npw_fft_kp,2,dimffnl_kp,ntypat)
   ! nonlocal form factors at k_prime

  real(dp), pointer :: kinpw_k(:) => null()
   ! kinpw_k(npw_fft_k)
   ! (modified) kinetic energy for each plane wave at k

  real(dp), pointer :: kinpw_kp(:) => null()
   ! kinpw_kp(npw_fft_kp)
   ! (modified) kinetic energy for each plane wave at k^prime

  real(dp), pointer :: kpg_k(:,:) => null()
   ! kpg_k(3,npw_fft_k)
   ! k+G vector coordinates at k

  real(dp), pointer :: kpg_kp(:,:) => null()
   ! kpg_kp(3,npw_fft_kp)
   ! k^prime+G vector coordinates at k^prime

  real(dp), ABI_CONTIGUOUS pointer :: phkpxred(:,:) => null()
   ! phkpxred(2,natom)
   ! phase factors exp(2 pi k^prime.xred) at k^prime

  real(dp), pointer :: ph3d_k(:,:,:) => null()
   ! ph3d_k(2,npw_fft_k,matblk)
   ! 3-dim structure factors, for each atom and plane wave at k

  real(dp), pointer :: ph3d_kp(:,:,:) => null()
   ! ph3d_kp(2,npw_fft_kp,matblk)
   ! 3-dim structure factors, for each atom and plane wave at k^prime

  real(dp), pointer :: vectornd(:,:,:,:,:) => null()
   ! vectornd(n4,n5,n6,nvloc,3)
   ! vector potential of nuclear magnetic dipoles
   ! in real space, on the augmented fft grid

  real(dp), pointer :: vlocal(:,:,:,:) => null()
   ! vlocal(n4,n5,n6,nvloc)
   ! local potential in real space, on the augmented fft grid

  real(dp), pointer :: vxctaulocal(:,:,:,:,:) => null()
   ! vxctaulocal(n4,n5,n6,nvloc,4)
   ! derivative of XC energy density with respect to kinetic energy density,
   ! in real space, on the augmented fft grid

  real(dp), ABI_CONTIGUOUS pointer :: xred(:,:) => null()
   ! xred(3,natom)
   ! reduced coordinates of atoms (dimensionless)

! ===== Structured datatype pointers

  type(fock_common_type), pointer :: fockcommon => null()
   ! common quantities needed to calculate Fock exact exchange

  type(fock_BZ_type), pointer :: fockbz => null()
   ! total brillouin zone quantities needed to calculate Fock exact exchange

  type(fock_ACE_type), pointer :: fockACE_k => null()
   ! ACE quantities needed to calculate Fock exact exchange in the ACE context

 contains
   procedure :: free => destroy_hamiltonian
    ! Free the memory in the GS Hamiltonian

   procedure :: load_spin => load_spin_hamiltonian
    ! Setup of the spin-dependent part of the GS Hamiltonian

   procedure :: load_k => load_k_hamiltonian
    ! Setup of the k-dependent part of the GS Hamiltonian

   procedure :: load_kprime => load_kprime_hamiltonian
    ! Setup of the k^prime-dependent part of the GS Hamiltonian

   !procedure :: copy => copy_hamiltonian

 end type gs_hamiltonian_type

 public :: init_hamiltonian         ! Initialize the GS Hamiltonian
 public :: copy_hamiltonian         ! Copy the object
!!***

!----------------------------------------------------------------------

!!****t* m_hamiltonian/rf_hamiltonian_type
!! NAME
!! rf_hamiltonian_type
!!
!! FUNCTION
!! This datastructure contains few data about one 1st-order Hamiltonian,
!! needed in the "getgh1c" routine, that apply the 1st-order Hamiltonian
!! on a wavefunction.
!!
!! SOURCE

 type,public :: rf_hamiltonian_type

! ===== Integer scalars

  integer :: cplex
   ! if 1, real space 1-order functions on FFT grid are REAL; if 2, COMPLEX

  integer :: dime1kb1
   ! First dimension of E1kb, derivative of Ekb with respect to a perturbation

  integer :: dime1kb2
   ! Second dimension of E1kb, derivative of Ekb with respect to a perturbation
   ! NCPP: dime1kb2=ntypat, PAW: dime1kb2=natom

  integer :: npw_k
   ! number of plane waves at k

  integer :: npw_kp
   ! number of plane waves at k^prime

  integer:: nspinor
   ! Number of spinorial components

  integer :: nsppol
   ! Total number of spin components (1=non-polarized, 2=polarized)

  integer :: nvloc
   ! Number of components of vloc
   ! usually, nvloc=1, except in the non-collinear magnetism case, where nvloc=4

  integer :: n4,n5,n6
   ! same as ngfft(4:6)

! ===== Real arrays

  real(dp), allocatable :: e1kbfr_spin(:,:,:,:,:)
   ! e1kbfr_spin(dimekb1,dimekb2,nspinor**2,cplex,my_nsppol)
   ! Contains the values of e1kbfr array for all spins treated by current process
   ! See e1kbfr description ; e1kbfr is pointer to e1kbfr_spin(:,:,:,:,isppol)

  real(dp), allocatable :: e1kbsc_spin(:,:,:,:,:)
   ! e1kbsc_spin(dimekb1,dimekb2,nspinor**2,cplex,my_nsppol)
   ! Contains the values of e1kbsc array for all spins treated by current process
   ! See e1kbsc description ; e1kbsc is pointer to e1kbsc_spin(:,:,:,:,isppol)

! ===== Real pointers

  real(dp), pointer :: dkinpw_k(:) => null()
   ! dkinpw_k(npw_k)
   ! 1st derivative of the (modified) kinetic energy for each plane wave at k

  real(dp), pointer :: dkinpw_kp(:) => null()
   ! dkinpw_kp(npw_kp)
   ! 1st derivative of the (modified) kinetic energy for each plane wave at k^prime

  real(dp), pointer :: ddkinpw_k(:) => null()
   ! ddkinpw_k(npw_k)
   ! 2nd derivative of the (modified) kinetic energy for each plane wave at k

  real(dp), pointer :: ddkinpw_kp(:) => null()
   ! ddkinpw_kp(npw_kp)
   ! 2nd derivative of the (modified) kinetic energy for each plane wave at k^prime

  real(dp), pointer :: e1kbfr(:,:,:,:) => null()
   ! Frozen part of 1st derivative of ekb for the considered perturbation
   ! (part not depending on VHxc^(1))
   ! e1kbfr(dime1kb1,dime1kb2,nspinor**2,cplex)
   ! For each spin component, e1kbfr points to e1kbfr_spin(:,:,:,:,my_isppol)

  real(dp), ABI_CONTIGUOUS pointer :: e1kbsc(:,:,:,:) => null()
   ! Self-consistent 1st derivative of ekb for the considered perturbation
   ! (part depending only on self-consistent VHxc^(1))
   ! e1kbsc(dime1kb1,dime1kb2,nspinor**2,cplex)
   ! For each spin component, e1kbfr points to e1kbfr_spin(:,:,:,:,my_isppol)

  real(dp), pointer :: vlocal1(:,:,:,:) => null()
   ! vlocal1(cplex*n4,n5,n6,nvloc)
   ! 1st-order local potential in real space, on the augmented fft grid

 contains
   procedure :: free => destroy_rf_hamiltonian
    ! Free the memory in the RF Hamiltonian

   procedure :: load_spin => load_spin_rf_hamiltonian
    ! Setup of the spin-dependent part of the RF Hamiltonian.

   procedure :: load_k => load_k_rf_hamiltonian
    ! Setup of the k-dependent part of the RF Hamiltonian

 end type rf_hamiltonian_type

 public :: init_rf_hamiltonian      ! Initialize the RF Hamiltonian
!!***

CONTAINS  !===========================================================

!----------------------------------------------------------------------

!!****f* m_hamiltonian/destroy_hamiltonian
!! NAME
!!  destroy_hamiltonian
!!
!! FUNCTION
!!  Clean and destroy gs_hamiltonian_type datastructure
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=All dynamic memory defined in the structure is deallocated.
!!
!! PARENTS
!!      d2frnl,dfpt_nselt,dfpt_nstdy,dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho
!!      dfptnl_resp,energy,fock2ACE,forstrnps,gwls_hamiltonian,ks_ddiago,m_gkk
!!      m_io_kss,m_phgamma,m_phpi,m_shirley,m_sigmaph,nonlop_test,vtorho
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine destroy_hamiltonian(Ham)

!Arguments ------------------------------------
!scalars
 class(gs_hamiltonian_type),intent(inout),target :: Ham

! *************************************************************************

 DBG_ENTER("COLL")

!@gs_hamiltonian_type

! Integer Pointers
 if (associated(Ham%gbound_kp,Ham%gbound_k)) then
   nullify(Ham%gbound_kp)
 else if (associated(Ham%gbound_kp)) then
   ABI_DEALLOCATE(Ham%gbound_kp)
 end if

! Integer arrays
 ABI_SFREE(Ham%atindx)
 ABI_SFREE(Ham%atindx1)
 ABI_SFREE(Ham%dimcprj)
 ABI_SFREE(Ham%gbound_k)
 ABI_SFREE(Ham%indlmn)
 ABI_SFREE(Ham%nattyp)
 ABI_SFREE(Ham%pspso)
 ABI_SFREE(Ham%typat)

! Real Pointers
 if (associated(Ham%phkpxred,Ham%phkxred)) then
   nullify(Ham%phkpxred)
 else if (associated(Ham%phkpxred)) then
   ABI_DEALLOCATE(Ham%phkpxred)
 end if
 if (allocated(Ham%phkxred)) then
   ABI_DEALLOCATE(Ham%phkxred)
 end if
 if (associated(Ham%ekb)) nullify(Ham%ekb)
 if (associated(Ham%vectornd)) nullify(Ham%vectornd)
 if (associated(Ham%vlocal)) nullify(Ham%vlocal)
 if (associated(Ham%vxctaulocal)) nullify(Ham%vxctaulocal)
 if (associated(Ham%xred)) nullify(Ham%xred)
 if (associated(Ham%kinpw_k)) nullify(Ham%kinpw_k)
 if (associated(Ham%kinpw_kp)) nullify(Ham%kinpw_kp)
 if (associated(Ham%kg_k)) nullify(Ham%kg_k)
 if (associated(Ham%kg_kp)) nullify(Ham%kg_kp)
 if (associated(Ham%kpg_k)) nullify(Ham%kpg_k)
 if (associated(Ham%kpg_kp)) nullify(Ham%kpg_kp)
 if (associated(Ham%ffnl_k)) nullify(Ham%ffnl_k)
 if (associated(Ham%ffnl_kp)) nullify(Ham%ffnl_kp)
 if (associated(Ham%ph3d_k)) nullify(Ham%ph3d_k)
 if (associated(Ham%ph3d_kp)) nullify(Ham%ph3d_kp)

! Real arrays
 ABI_SFREE(Ham%ekb_spin)
 ABI_SFREE(Ham%sij)
 ABI_SFREE(Ham%nucdipmom)
 ABI_SFREE(Ham%ph1d)

! Complex arrays
 ABI_SFREE(Ham%nucdipmom_k)

! Structured datatype pointers
 if (associated(Ham%fockcommon)) nullify(Ham%fockcommon)
 if (associated(Ham%fockACE_k)) nullify(Ham%fockACE_k)
 if (associated(Ham%fockbz)) nullify(Ham%fockbz)
#if defined HAVE_GPU_CUDA
 if(Ham%use_gpu_cuda==1) then
   call gpu_finalize_ham_data()
 end if
#endif

 DBG_EXIT("COLL")

end subroutine destroy_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_hamiltonian
!! NAME
!!  init_hamiltonian
!!
!! FUNCTION
!!  Creation method for the gs_hamiltonian_type structure.
!!  It allocates memory and initializes all quantities that do not depend on the k-point or spin.
!!
!! INPUTS
!!  [comm_atom]=optional, MPI communicator over atoms
!!  [fockcommon <type(fock_common_type)>]= common quantities to calculate Fock exact exchange
!!  [fockbz <type(fock_BZ_type)>]= quantities to calculate Fock exact exchange in the total BZ
!!  natom=Number of atoms in the unit cell.
!!  nfft=Number of FFT grid points (for this processors).
!!  nspinor=Number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=Number of spin density components.
!!  mgfft=Maximum size for 1D FFTs i.e., MAXVAL(ngfft(1:3))
!!  [mpi_atmtab(:)]=optional, indexes of the atoms treated by current proc
!!  [mpi_spintab(2)]=optional, flags defining the spin(s) treated be current process:
!!                   mpi_spintab(1)=1 if non-polarized or spin-up treated
!!                   mpi_spintab(2)=1 if polarized and spin-dn treated
!!  psps<pseudopotential_type>=structure datatype gathering data on the pseudopotentials.
!!  [electronpositron<electronpositron_type>]=Structured datatype storing data for the
!!    electron-positron two-component DFT (optional).
!!  ngfft(18)=integer array with FFT box dimensions and other information on FFTs, for the FINE rectangular grid.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator
!!  [nucdipmom(3,natom)]= (optional) array of nuclear dipole moments at atomic sites
!!  [ph1d(2,3*(2*mgfft+1)*natom)]=1-dimensions phase arrays for structure factor (see getph.f).
!!              Optional, recalculated inside the routine if not present in input.
!!  rprimd(3,3)=Direct lattice vectors in Bohr.
!!  typat(natom)=Type of each atom.
!!  [usecprj]=flag use only for PAW; 1 if cprj datastructure is allocated
!!  xred(3,natom)=Reduced coordinates of the atoms.
!!  pawtab(ntypat*psps%usepaw)<pawtab_type>=PAW TABulated data initialized at start.
!!  [paw_ij(:) <type(paw_ij_type)>]=optional, paw arrays given on (i,j) channels
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=Structured datatype almost completely initialized:
!!   * Basic variables and dimensions are transfered to the structure.
!!   * All pointers are allocated with correct dimensions.
!!   * Quantities that do not depend on the k-point or spin are initialized.
!!
!! PARENTS
!!      d2frnl,dfpt_nselt,dfpt_nstdy,dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho
!!      dfptnl_resp,energy,fock2ACE,forstrnps,ks_ddiago,m_gkk,m_io_kss
!!      m_phgamma,m_phpi,m_shirley,m_sigmaph,nonlop_test,vtorho
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_hamiltonian(ham,Psps,pawtab,nspinor,nsppol,nspden,natom,typat,&
&                           xred,nfft,mgfft,ngfft,rprimd,nloalg,&
&                           ph1d,usecprj,comm_atom,mpi_atmtab,mpi_spintab,paw_ij,&  ! optional
&                           electronpositron,fock,nucdipmom,use_gpu_cuda)           ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,natom,nspinor,nsppol,nspden,mgfft
 integer,optional,intent(in) :: comm_atom,usecprj,use_gpu_cuda
 class(gs_hamiltonian_type),intent(inout),target :: ham
 type(electronpositron_type),optional,pointer :: electronpositron
 type(fock_type),optional,pointer :: fock
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfft(18),nloalg(3),typat(natom)
 integer,optional,intent(in)  :: mpi_atmtab(:),mpi_spintab(2)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in),target :: xred(3,natom)
 real(dp),optional,intent(in) :: nucdipmom(3,natom),ph1d(2,3*(2*mgfft+1)*natom)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_ij_type),optional,intent(in) :: paw_ij(:)

!Local variables-------------------------------
!scalars
 integer :: my_comm_atom,my_nsppol,itypat,iat,ilmn,indx,isp,cplex_dij,jsp
 real(dp) :: ucvol
!arrays
 integer :: my_spintab(2)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable,target :: ekb_tmp(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 !@gs_hamiltonian_type

!Manage optional parameters
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_spintab=0;my_spintab(1:nsppol)=1;if (present(mpi_spintab)) my_spintab(1:2)=mpi_spintab(1:2)
 my_nsppol=count(my_spintab==1)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ABI_CHECK(mgfft==MAXVAL(ngfft(1:3)),"Wrong mgfft")

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 ABI_ALLOCATE(ham%atindx,(natom))
 ABI_ALLOCATE(ham%atindx1,(natom))
 ABI_ALLOCATE(ham%typat,(natom))
 ABI_ALLOCATE(ham%indlmn,(6,psps%lmnmax,psps%ntypat))
 ABI_ALLOCATE(ham%nattyp,(psps%ntypat))
 ABI_ALLOCATE(ham%nucdipmom,(3,natom))
 ABI_ALLOCATE(ham%ph1d,(2,3*(2*mgfft+1)*natom))
 ABI_ALLOCATE(ham%pspso,(psps%ntypat))

!Initialize most of the Hamiltonian
 indx=1
 do itypat=1,psps%ntypat
   ham%nattyp(itypat)=0
   do iat=1,natom
     if (typat(iat)==itypat) then
       ham%atindx (iat )=indx
       ham%atindx1(indx)=iat
       indx=indx+1
       ham%nattyp(itypat)=ham%nattyp(itypat)+1
     end if
   end do
 end do

 ham%gmet(:,:)  =gmet(:,:)
 ham%gprimd(:,:)=gprimd(:,:)
 ham%indlmn(:,:,:)=psps%indlmn(:,:,:)
 ham%lmnmax     =psps%lmnmax
 ham%mgfft      =mgfft
 ham%mpsang     =psps%mpsang
 ham%mpssoang   =psps%mpssoang
 ham%natom      =natom
 ham%nfft       =nfft
 ham%ngfft(:)   =ngfft(:)
 ham%nloalg(:)  =nloalg(:)
 ham%matblk=min(NLO_MINCAT,maxval(ham%nattyp)); if (nloalg(2)>0) ham%matblk=natom
 ham%nsppol     =nsppol
 ham%nspinor    =nspinor
 ham%ntypat     =psps%ntypat
 ham%typat      =typat(1:natom)
 ham%nvloc=1; if(nspden==4)ham%nvloc=4
 ham%n4         =ngfft(4)
 ham%n5         =ngfft(5)
 ham%n6         =ngfft(6)
 ham%usepaw     =psps%usepaw
 ham%ucvol      =ucvol
 ham%useylm     =psps%useylm
 ham%use_gpu_cuda=0 ; if(PRESENT(use_gpu_cuda)) ham%use_gpu_cuda=use_gpu_cuda

 ham%pspso(:)   =psps%pspso(1:psps%ntypat)
 if (psps%usepaw==1) then
   do itypat=1,psps%ntypat
     ham%pspso(itypat)=1+pawtab(itypat)%usespnorb
   end do
 end if

 if (present(nucdipmom)) then
   ham%nucdipmom(:,:) = nucdipmom(:,:)
 else
   ham%nucdipmom(:,:) = zero
 end if

 ham%xred => xred

 if (present(fock)) then
   if (associated(fock)) then
     ham%fockcommon => fock%fock_common
     if (fock%fock_common%use_ACE==0) ham%fockbz=>fock%fock_BZ
   end if
 end if

 if (present(ph1d)) then
   ham%ph1d(:,:) = ph1d(:,:)
 else ! Recalculate structure factor phases
   call getph(ham%atindx,natom,ngfft(1),ngfft(2),ngfft(3),ham%ph1d,xred)
 end if

 if (ham%usepaw==1) then
   ham%usecprj=0;if (present(usecprj)) ham%usecprj=usecprj
   ABI_ALLOCATE(ham%dimcprj,(natom))
   !Be carefull cprj are ordered by atom type (used in non-local operator)
   call pawcprj_getdim(ham%dimcprj,natom,ham%nattyp,ham%ntypat,ham%typat,pawtab,'O')
 else
   ham%usecprj=0
   ABI_ALLOCATE(ham%dimcprj,(0))
 end if

! ===========================
! ==== Non-local factors ====
! ===========================


 if (ham%usepaw==0) then ! Norm-conserving: use constant Kleimann-Bylander energies.
   ham%dimekb1=psps%dimekb
   ham%dimekb2=psps%ntypat
   ham%dimekbq=1
   ABI_ALLOCATE(ham%ekb_spin,(psps%dimekb,psps%ntypat,nspinor**2,1,1))
   ham%ekb => ham%ekb_spin(:,:,:,:,1)
   ABI_ALLOCATE(ham%sij,(0,0))
   ham%ekb(:,:,1,1)=psps%ekb(:,:)
   if (nspinor==2) then
     ham%ekb(:,:,2,1)=psps%ekb(:,:)
     ham%ekb(:,:,3:4,1)=zero
   end if
   if (PRESENT(electronpositron)) then
     if (electronpositron_calctype(electronpositron)==1) ham%ekb(:,:,:,:)=-ham%ekb(:,:,:,:)
   end if

!  Update enl on GPU (will do it later for PAW)
#if defined HAVE_GPU_CUDA
   if (ham%use_gpu_cuda==1) then
     call gpu_update_ham_data(ham%ekb(:,:,:,1),size(ham%ekb),ham%sij,size(ham%sij),ham%gprimd,size(ham%gprimd))
   end if
#endif

 else ! PAW: store overlap coefficients (spin non dependent) and Dij coefficients (spin dependent)
   cplex_dij=1
   if (present(paw_ij)) then
     if (size(paw_ij)>0) cplex_dij=paw_ij(1)%cplex_dij
   end if
   if ((nspinor==2).or.any(abs(ham%nucdipmom)>tol8)) cplex_dij=2
   ham%dimekb1=psps%dimekb*cplex_dij
   ham%dimekb2=natom
   ham%dimekbq=1
   if (present(paw_ij)) then
     if (size(paw_ij)>0) ham%dimekbq=paw_ij(1)%qphase
   end if
   ABI_ALLOCATE(ham%sij,(ham%dimekb1,psps%ntypat))
   do itypat=1,psps%ntypat
     if (cplex_dij==1) then
       ham%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do ilmn=1,pawtab(itypat)%lmn2_size
         ham%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
         ham%sij(2*ilmn  ,itypat)=zero
       end do
     end if
     if (cplex_dij*pawtab(itypat)%lmn2_size<ham%dimekb1) then
       ham%sij(cplex_dij*pawtab(itypat)%lmn2_size+1:ham%dimekb1,itypat)=zero
     end if
   end do
   !We preload here PAW non-local factors in order to avoid a communication over atoms
   ! inside the loop over spins.
   ABI_ALLOCATE(ham%ekb_spin,(ham%dimekb1,ham%dimekb2,nspinor**2,ham%dimekbq,my_nsppol))
   ham%ekb_spin=zero
   if (present(paw_ij)) then
     if (my_nsppol<ham%nsppol) then
       ABI_ALLOCATE(ekb_tmp,(ham%dimekb1,ham%dimekb2,nspinor**2,ham%dimekbq))
     end if
     jsp=0
     do isp=1,ham%nsppol
       if (my_spintab(isp)==1) then
         jsp=jsp+1 ; ham%ekb => ham%ekb_spin(:,:,:,:,jsp)
       else
         ham%ekb => ekb_tmp
       end if
       if (present(mpi_atmtab)) then
         call pawdij2ekb(ham%ekb,paw_ij,isp,my_comm_atom,mpi_atmtab=mpi_atmtab)
       else
         call pawdij2ekb(ham%ekb,paw_ij,isp,my_comm_atom)
       end if
     end do
     if (my_nsppol<ham%nsppol) then
       ABI_DEALLOCATE(ekb_tmp)
     end if
   end if
   nullify(ham%ekb)
 end if

 DBG_EXIT("COLL")

end subroutine init_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_k_hamiltonian
!! NAME
!!  load_k_hamiltonian
!!
!! FUNCTION
!!  Setup of the k-dependent part of the Hamiltonian H_k_k^prime
!!
!! INPUTS
!!  [compute_gbound]=flag. if true, G sphere boundary is computed here
!!  [compute_ph3d]=flag. if true, 3D structure factors are computed here (only if nloalg(1)>0)
!!  [gbound_k]=G sphere boundary (not compatible with compute_gbound=TRUE)
!!  [ffnl_k]=nonlocal form factors on basis sphere
!!  [istwf_k]=parameter that describes the storage of wfs
!!  [kinpw_k]=(modified) kinetic energy for each plane wave
!!  [kg_k]=planewave reduced coordinates in basis sphere (g vectors)
!!  [kpg_k]=(k+g) vectors in reciprocal space
!!  [kpt_k]=k point coordinates
!!  [npw_k]=number of plane waves (processed by current proc when band-FFT parallelism is on)
!!  [npw_fft_k]=number of plane waves used to apply Hamiltonian (in the "FFT" configuration)
!!  [ph3d_k]=3-dim structure factors, for each atom and plane wave
!!
!! SIDE EFFECTS
!!  ham<gs_hamiltonian_type>=structured datatype completed with k-dependent quantitites.
!!          Quantities at k^prime are set equal to quantities at k.
!!    k-dependent scalars and pointers associated
!!    phkxred=exp(.k.xred) for each atom
!!    [ham%gbound_k]=G sphere boundary, for each plane wave
!!    [ham%ph3d_k]=3-dim structure factors, for each atom and plane wave
!!
!! PARENTS
!!      d2frnl,dfpt_nsteltwf,dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,dfptnl_resp
!!      energy,fock2ACE,forstrnps,getgh1c,gwls_hamiltonian,ks_ddiago,m_io_kss
!!      m_shirley,nonlop_test,vtorho,wfd_vnlpsi
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_k_hamiltonian(ham,ffnl_k,fockACE_k,gbound_k,istwf_k,kinpw_k,&
                              kg_k,kpg_k,kpt_k,nucdipmom_k,npw_k,npw_fft_k,ph3d_k,&
                              compute_gbound,compute_ph3d)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: npw_k,npw_fft_k,istwf_k
 logical,intent(in),optional :: compute_gbound,compute_ph3d
 class(gs_hamiltonian_type),intent(inout),target :: ham
!arrays
 integer,intent(in),optional,target :: gbound_k(:,:),kg_k(:,:)
 real(dp),intent(in),optional :: kpt_k(3)
 real(dp),intent(in),optional,target :: ffnl_k(:,:,:,:),kinpw_k(:),kpg_k(:,:),ph3d_k(:,:,:)
 complex(dpc),intent(in),optional :: nucdipmom_k(:)
 type(fock_ACE_type),intent(in),optional,target :: fockACE_k

!Local variables-------------------------------
!scalars
 integer :: iat,iatom
 logical :: compute_gbound_
 real(dp) :: arg
 !character(len=500) :: msg

! *************************************************************************

 DBG_ENTER("COLL")

!@gs_hamiltonian_type

!k-dependent scalars
 if (present(kpt_k)) then
   ham%kpt_k(:)  = kpt_k(:)
   ham%kpt_kp(:) = kpt_k(:)
 end if
 if (present(istwf_k)) then
   ham%istwf_k  = istwf_k
   ham%istwf_kp = istwf_k
 end if
 if (present(npw_k)) then
   ham%npw_k  = npw_k
   ham%npw_kp = npw_k
 end if
 if (present(npw_fft_k)) then
   ham%npw_fft_k  = npw_fft_k
   ham%npw_fft_kp = npw_fft_k
 else if (present(npw_k)) then
   ham%npw_fft_k  = npw_k
   ham%npw_fft_kp = npw_k
 end if

 ! k-dependend complex quantities
  if (present(nucdipmom_k)) then
   ham%nucdipmom_k(:) = nucdipmom_k(:)
  end if

!Pointers to k-dependent quantitites
 if (present(kinpw_k)) then
   ham%kinpw_k  => kinpw_k
   ham%kinpw_kp => kinpw_k
 end if
 if (present(kg_k)) then
   ham%kg_k  => kg_k
   ham%kg_kp => kg_k
 end if
 if (present(kpg_k)) then
   ham%kpg_k  => kpg_k
   ham%kpg_kp => kpg_k
 end if
 if (present(ffnl_k)) then
   ham%ffnl_k  => ffnl_k
   ham%ffnl_kp => ffnl_k
 end if
 if (present(ph3d_k)) then
   ham%ph3d_k  => ph3d_k
   ham%ph3d_kp => ph3d_k
 end if
 if (present(fockACE_k)) then
   ham%fockACE_k  => fockACE_k
 end if
!Compute exp(i.k.R) for each atom
 if (present(kpt_k)) then
   if (associated(Ham%phkpxred).and.(.not.associated(Ham%phkpxred,Ham%phkxred))) then
     ABI_DEALLOCATE(Ham%phkpxred)
   end if
   if (allocated(ham%phkxred)) then
     ABI_DEALLOCATE(ham%phkxred)
   end if
   ABI_ALLOCATE(ham%phkxred,(2,ham%natom))
   do iat=1,ham%natom
     iatom=ham%atindx(iat)
     arg=two_pi*DOT_PRODUCT(kpt_k,ham%xred(:,iat))
     ham%phkxred(1,iatom)=DCOS(arg)
     ham%phkxred(2,iatom)=DSIN(arg)
   end do
   ham%phkpxred => ham%phkxred
 end if

!Compute or copy G sphere boundary at k+g
 compute_gbound_=.false.;if (present(compute_gbound)) compute_gbound_=compute_gbound
 if (present(gbound_k)) compute_gbound_=.true.
 if (compute_gbound_) then
   if (associated(Ham%gbound_kp,Ham%gbound_k)) then
     nullify(Ham%gbound_kp)
   else if (associated(Ham%gbound_kp)) then
     ABI_DEALLOCATE(Ham%gbound_kp)
   end if
   if (allocated(ham%gbound_k)) then
     ABI_DEALLOCATE(ham%gbound_k)
   end if
 end if
 if (.not.allocated(ham%gbound_k)) then
   ABI_ALLOCATE(ham%gbound_k,(2*ham%mgfft+8,2))
   ham%gbound_k(:,:)=0
   ham%gbound_kp => ham%gbound_k
 end if
 if (compute_gbound_) then
   if (present(gbound_k)) then
     ham%gbound_k(:,:)=gbound_k(:,:)
   else
     if (.not.associated(ham%kg_k)) then
       MSG_BUG('Something is missing for gbound_k computation!')
     end if
     !write(std_out,*)"About to call sphereboundary"
     !write(std_out,*)"size(kg_k), npw_k, mgfft",size(ham%kg_k, dim=2), ham%npw_k, ham%mgfft
     call sphereboundary(ham%gbound_k,ham%istwf_k,ham%kg_k,ham%mgfft,ham%npw_k)
   end if
   ham%gbound_kp => ham%gbound_k
 end if

!Compute 3D structure factors for each atom at k+g
 if (present(compute_ph3d).and.present(ph3d_k)) then
   if (compute_ph3d.and.ham%nloalg(2)>0) then
     if ((.not.allocated(ham%phkxred)).or.(.not.associated(ham%kg_k)).or.&
         (.not.associated(ham%ph3d_k))) then
       MSG_BUG('Something is missing for ph3d_k computation!')
     end if
     call ph1d3d(1,ham%natom,ham%kg_k,ham%matblk,ham%natom,ham%npw_k,ham%ngfft(1),&
                 ham%ngfft(2),ham%ngfft(3),ham%phkxred,ham%ph1d,ham%ph3d_k)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine load_k_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_kprime_hamiltonian
!! NAME
!!  load_kprime_hamiltonian
!!
!! FUNCTION
!!  Setup of the k^prime-dependent part of the Hamiltonian H_k_k^prime
!!
!! INPUTS
!!  [compute_gbound]=flag. if true, G sphere boundary is computed here
!!  [compute_ph3d]=flag. if true, 3D structure factors are computed here (only if nloalg(2)>0)
!!  [gbound_kp]=G sphere boundary (not compatible with compute_gbound=TRUE)
!!  [ffnl_kp]=nonlocal form factors on basis sphere
!!  [istwf_kp]=parameter that describes the storage of wfs
!!  [kinpw_kp]=(modified) kinetic energy for each plane wave
!!  [kg_kp]=planewave reduced coordinates in basis sphere (g vectors)
!!  [kpg_kp]=(k+g) vectors in reciprocal space
!!  [kpt_kp]=k point coordinates
!!  [npw_kp]=number of plane waves (processed by current proc when band-FFT parallelism is on)
!!  [npw_fft_kp]=number of plane waves used to apply Hamiltonian (in the "FFT" configuration)
!!  [ph3d_kp]=3-dim structure factors, for each atom and plane wave
!!
!! SIDE EFFECTS
!!  ham<gs_hamiltonian_type>=structured datatype completed with k^prime-dependent quantitites.
!!    k^prime-dependent scalars and pointers associated
!!    phkpxred=exp(.k^prime.xred) for each atom
!!    [ham%gbound_kp]=G sphere boundary, for each plane wave
!!    [ham%ph3d_kp]=3-dim structure factors at k^prime at k, for each atom and plane wave
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,fock_getghc,getgh1c
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_kprime_hamiltonian(ham,ffnl_kp,gbound_kp,istwf_kp,kinpw_kp,&
                                   kg_kp,kpg_kp,kpt_kp,npw_kp,npw_fft_kp,&
                                   ph3d_kp,compute_gbound,compute_ph3d)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: npw_kp,npw_fft_kp,istwf_kp
 logical,intent(in),optional :: compute_gbound,compute_ph3d
 class(gs_hamiltonian_type),intent(inout),target :: ham
!arrays
 integer,intent(in),optional,target :: gbound_kp(:,:),kg_kp(:,:)
 real(dp),intent(in),optional :: kpt_kp(3)
 real(dp),intent(in),optional,target :: ffnl_kp(:,:,:,:),kinpw_kp(:),kpg_kp(:,:),ph3d_kp(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,iatom
 logical :: compute_gbound_
 real(dp) :: arg
 character(len=100) :: msg

! *************************************************************************

 DBG_ENTER("COLL")

!@gs_hamiltonian_type

!k-dependent scalars
 if (present(kpt_kp))   ham%kpt_kp(:)= kpt_kp(:)
 if (present(istwf_kp)) ham%istwf_kp = istwf_kp
 if (present(npw_kp))   ham%npw_kp   = npw_kp
 if (present(npw_fft_kp)) then
    ham%npw_fft_kp = npw_fft_kp
 else if (present(npw_kp)) then
    ham%npw_fft_kp = npw_kp
 end if

!Pointers to k-dependent quantitites
 if (present(kinpw_kp)) ham%kinpw_kp => kinpw_kp
 if (present(kg_kp))    ham%kg_kp    => kg_kp
 if (present(kpg_kp))   ham%kpg_kp   => kpg_kp
 if (present(ffnl_kp))  ham%ffnl_kp  => ffnl_kp
 if (present(ph3d_kp))  ham%ph3d_kp  => ph3d_kp

!Compute exp(i.k^prime.R) for each atom
 if (present(kpt_kp)) then
   if (associated(ham%phkpxred,ham%phkxred)) then
     nullify(ham%phkpxred)
   else if (associated(ham%phkpxred)) then
     ABI_DEALLOCATE(ham%phkpxred)
   end if
   ABI_ALLOCATE(ham%phkpxred,(2,ham%natom))
   do iat=1,ham%natom
     iatom=ham%atindx(iat)
     arg=two_pi*DOT_PRODUCT(kpt_kp,ham%xred(:,iat))
     ham%phkpxred(1,iatom)=DCOS(arg)
     ham%phkpxred(2,iatom)=DSIN(arg)
   end do
 end if

!Compute or copy G sphere boundary at k^prime+g
 compute_gbound_=.false.
 if (present(kpt_kp).and.present(compute_gbound)) compute_gbound_=compute_gbound
 if (present(gbound_kp)) compute_gbound_=.true.
 if (compute_gbound_) then
   if (associated(ham%gbound_kp,ham%gbound_k)) then
     nullify(ham%gbound_kp)
   else if (associated(ham%gbound_kp)) then
     ABI_DEALLOCATE(ham%gbound_kp)
   end if
   if (present(gbound_kp)) then
     ham%gbound_kp(:,:)=gbound_kp(:,:)
   else
     if (.not.associated(ham%kg_kp)) then
       msg='Something is missing for gbound_kp computation!'
       MSG_BUG(msg)
     end if
     ABI_ALLOCATE(ham%gbound_kp,(2*ham%mgfft+8,2))
     call sphereboundary(ham%gbound_kp,ham%istwf_kp,ham%kg_kp,ham%mgfft,ham%npw_kp)
   end if
 end if

!Compute 3D structure factors for each atom at k^prime+g
 if (present(compute_ph3d).and.present(ph3d_kp)) then
   if (compute_ph3d.and.ham%nloalg(2)>0) then
     if ((.not.associated(ham%phkpxred)).or.(.not.associated(ham%kg_kp)).or.&
&        (.not.associated(ham%ph3d_kp))) then
       msg='Something is missing for ph3d_kp computation!'
       MSG_BUG(msg)
     end if
     call ph1d3d(1,ham%natom,ham%kg_kp,ham%matblk,ham%natom,ham%npw_kp,ham%ngfft(1),&
&                ham%ngfft(2),ham%ngfft(3),ham%phkpxred,ham%ph1d,ham%ph3d_kp)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine load_kprime_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/copy_hamiltonian
!! NAME
!!  copy_hamiltonian
!!
!! INPUTS
!!  gs_hamk_in<gs_hamiltonian_type>=Structured datatype completely initialized,
!!                                  to be copied.
!!
!! FUNCTION
!!  Copy a gs_hamiltonian_type variable (gs_hamk_in) in another (gs_hamk_out).
!!  In contrast to an assignment statement (gs_hamk_out=gs_hamk_in), this
!!  subroutine allocate memory space for the pointers contained in the data
!!  structure (gs_hamk_out) and copy the content of the corresponding memory
!!  space of gs_hamk_in in it. In contrast, the assignment statement would
!!  only associate the pointers of gs_hamk_out to the same memory space than
!!  the corresponding ones in gs_hamk_in. This can cause trouble if one data
!!  structure is destroyed before a reading/writing statement for the other
!!  structure, causing access to unallocated memory space (silently, without
!!  segmentation fault being generated).
!!
!! OUTPUT
!!  gs_hamk_out<gs_hamiltonian_type>=Structured datatype containing separate
!!                                   copies of all data of gs_hamk_in upon exit.
!!
!! PARENTS
!!      gwls_hamiltonian
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine copy_hamiltonian(gs_hamk_out,gs_hamk_in)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(in),target :: gs_hamk_in
 type(gs_hamiltonian_type),intent(out),target :: gs_hamk_out

!Local variables-------------------------------
 integer :: tmp2i(5)
#if defined HAVE_FC_ISO_C_BINDING
 type(C_PTR) :: ham_ptr
#endif

! *************************************************************************

 DBG_ENTER("COLL")

!@gs_hamiltonian_type

 gs_hamk_out%dimekb1 = gs_hamk_in%dimekb1
 gs_hamk_out%dimekb2 = gs_hamk_in%dimekb2
 gs_hamk_out%dimekbq = gs_hamk_in%dimekbq
 gs_hamk_out%istwf_k = gs_hamk_in%istwf_k
 gs_hamk_out%istwf_kp = gs_hamk_in%istwf_kp
 gs_hamk_out%lmnmax = gs_hamk_in%lmnmax
 gs_hamk_out%matblk = gs_hamk_in%matblk
 gs_hamk_out%mgfft = gs_hamk_in%mgfft
 gs_hamk_out%mpsang = gs_hamk_in%mpsang
 gs_hamk_out%mpssoang = gs_hamk_in%mpssoang
 gs_hamk_out%natom = gs_hamk_in%natom
 gs_hamk_out%nfft = gs_hamk_in%nfft
 gs_hamk_out%npw_k = gs_hamk_in%npw_k
 gs_hamk_out%npw_kp = gs_hamk_in%npw_kp
 gs_hamk_out%npw_fft_k = gs_hamk_in%npw_fft_k
 gs_hamk_out%npw_fft_kp = gs_hamk_in%npw_fft_kp
 gs_hamk_out%nspinor = gs_hamk_in%nspinor
 gs_hamk_out%nsppol = gs_hamk_in%nsppol
 gs_hamk_out%ntypat = gs_hamk_in%ntypat
 gs_hamk_out%nvloc = gs_hamk_in%nvloc
 gs_hamk_out%n4 = gs_hamk_in%n4
 gs_hamk_out%n5 = gs_hamk_in%n5
 gs_hamk_out%n6 = gs_hamk_in%n6
 gs_hamk_out%use_gpu_cuda = gs_hamk_in%use_gpu_cuda
 gs_hamk_out%usecprj = gs_hamk_in%usecprj
 gs_hamk_out%usepaw = gs_hamk_in%usepaw
 gs_hamk_out%useylm = gs_hamk_in%useylm
 gs_hamk_out%ngfft = gs_hamk_in%ngfft
 gs_hamk_out%nloalg = gs_hamk_in%nloalg
 gs_hamk_out%ucvol = gs_hamk_in%ucvol
 gs_hamk_out%gmet = gs_hamk_in%gmet
 gs_hamk_out%gprimd = gs_hamk_in%gprimd
 gs_hamk_out%kpt_k = gs_hamk_in%kpt_k
 gs_hamk_out%kpt_kp = gs_hamk_in%kpt_kp

 ABI_ALLOCATE(gs_hamk_out%atindx,(gs_hamk_out%natom))
 gs_hamk_out%atindx = gs_hamk_in%atindx
 ABI_ALLOCATE(gs_hamk_out%atindx1,(gs_hamk_out%natom))
 gs_hamk_out%atindx1 = gs_hamk_in%atindx1
 ABI_ALLOCATE(gs_hamk_out%dimcprj,(gs_hamk_out%natom*gs_hamk_out%usepaw))
 if (gs_hamk_out%usepaw==1) gs_hamk_out%dimcprj = gs_hamk_in%dimcprj
 ABI_ALLOCATE(gs_hamk_out%typat,(gs_hamk_out%natom))
 gs_hamk_out%typat = gs_hamk_in%typat
 ABI_ALLOCATE(gs_hamk_out%gbound_k,(2*gs_hamk_out%mgfft+8,2))
 gs_hamk_out%gbound_k = gs_hamk_in%gbound_k
 ABI_ALLOCATE(gs_hamk_out%indlmn,(6,gs_hamk_out%lmnmax,gs_hamk_out%ntypat))
 gs_hamk_out%indlmn = gs_hamk_in%indlmn
 ABI_ALLOCATE(gs_hamk_out%nattyp,(gs_hamk_out%ntypat))
 gs_hamk_out%nattyp = gs_hamk_in%nattyp
 ABI_ALLOCATE(gs_hamk_out%nucdipmom,(3,gs_hamk_out%natom))
 gs_hamk_out%nucdipmom = gs_hamk_in%nucdipmom
 ABI_ALLOCATE(gs_hamk_out%phkxred,(2,gs_hamk_out%natom))
 gs_hamk_out%phkxred = gs_hamk_in%phkxred
 ABI_ALLOCATE(gs_hamk_out%ph1d,(2,3*(2*gs_hamk_out%mgfft+1)*gs_hamk_out%natom))
 gs_hamk_out%ph1d = gs_hamk_in%ph1d
 ABI_ALLOCATE(gs_hamk_out%pspso,(gs_hamk_out%ntypat))
 gs_hamk_out%pspso = gs_hamk_in%pspso
 tmp2i(1:5)=shape(gs_hamk_in%ekb_spin)
 ABI_ALLOCATE(gs_hamk_out%ekb_spin,(tmp2i(1),tmp2i(2),tmp2i(3),tmp2i(4),tmp2i(5)))
 gs_hamk_out%ekb_spin = gs_hamk_in%ekb_spin
 gs_hamk_out%ekb => gs_hamk_out%ekb_spin(:,:,:,:,1)
 tmp2i(1:2)=shape(gs_hamk_in%sij)
 ABI_ALLOCATE(gs_hamk_out%sij,(tmp2i(1),tmp2i(2)))
 gs_hamk_out%sij = gs_hamk_in%sij

 if (associated(gs_hamk_in%gbound_kp,gs_hamk_in%gbound_k)) then
   gs_hamk_out%gbound_kp => gs_hamk_out%gbound_k
 else
   ABI_ALLOCATE(gs_hamk_out%gbound_kp,(2,gs_hamk_out%natom))
   gs_hamk_out%gbound_kp = gs_hamk_in%gbound_kp
 end if
 if (associated(gs_hamk_in%phkpxred,gs_hamk_in%phkxred)) then
   gs_hamk_out%phkpxred => gs_hamk_out%phkxred
 else
   ABI_ALLOCATE(gs_hamk_out%phkpxred,(2,gs_hamk_out%natom))
   gs_hamk_out%phkpxred = gs_hamk_in%phkpxred
 end if

 call addr_copy(gs_hamk_in%xred,gs_hamk_out%xred)
 call addr_copy(gs_hamk_in%vectornd,gs_hamk_out%vectornd)
 call addr_copy(gs_hamk_in%vlocal,gs_hamk_out%vlocal)
 call addr_copy(gs_hamk_in%vxctaulocal,gs_hamk_out%vxctaulocal)
 call addr_copy(gs_hamk_in%kinpw_k,gs_hamk_out%kinpw_k)
 call addr_copy(gs_hamk_in%kinpw_kp,gs_hamk_out%kinpw_kp)
 call addr_copy(gs_hamk_in%kg_k,gs_hamk_out%kg_k)
 call addr_copy(gs_hamk_in%kg_kp,gs_hamk_out%kg_kp)
 call addr_copy(gs_hamk_in%kpg_k,gs_hamk_out%kpg_k)
 call addr_copy(gs_hamk_in%kpg_kp,gs_hamk_out%kpg_kp)
 call addr_copy(gs_hamk_in%ffnl_k,gs_hamk_out%ffnl_k)
 call addr_copy(gs_hamk_in%ffnl_kp,gs_hamk_out%ffnl_kp)
 call addr_copy(gs_hamk_in%ph3d_k,gs_hamk_out%ph3d_k)
 call addr_copy(gs_hamk_in%ph3d_kp,gs_hamk_out%ph3d_kp)

!For pointers to structured datatypes, have to copy the address
!manually because there is no generic addr_copy function for that
 if (associated(gs_hamk_in%fockcommon)) then
#if defined HAVE_FC_ISO_C_BINDING
   ham_ptr=c_loc(gs_hamk_in%fockcommon)
   call c_f_pointer(ham_ptr,gs_hamk_out%fockcommon)
#else
   gs_hamk_out%fockcommon=transfer(gs_hamk_in%fockcommon,gs_hamk_out%fockcommon)
#endif
 else
   nullify(gs_hamk_out%fockcommon)
 end if
 if (associated(gs_hamk_in%fockbz)) then
#if defined HAVE_FC_ISO_C_BINDING
   ham_ptr=c_loc(gs_hamk_in%fockbz)
   call c_f_pointer(ham_ptr,gs_hamk_out%fockbz)
#else
   gs_hamk_out%fockbz=transfer(gs_hamk_in%fockbz,gs_hamk_out%fockbz)
#endif
 else
   nullify(gs_hamk_out%fockbz)
 end if
 if (associated(gs_hamk_in%fockACE_k)) then
#if defined HAVE_FC_ISO_C_BINDING
   ham_ptr=c_loc(gs_hamk_in%fockACE_k)
   call c_f_pointer(ham_ptr,gs_hamk_out%fockACE_k)
#else
   gs_hamk_out%fockACE_k=transfer(gs_hamk_in%fockACE_k,gs_hamk_out%fockACE_k)
#endif
 else
   nullify(gs_hamk_out%fockACE_k)
 end if

 DBG_EXIT("COLL")

end subroutine copy_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_spin_hamiltonian
!! NAME
!!  load_spin_hamiltonian
!!
!! INPUTS
!!  isppol=index of current spin
!!  [vectornd(n4,n5,n6,nvloc,3)]=optional, vector potential of nuclear magnetic dipoles in real space
!!  [vlocal(n4,n5,n6,nvloc)]=optional, local potential in real space
!!  [vxctaulocal(n4,n5,n6,nvloc,4)]=optional, derivative of XC energy density with respect
!!                                  to kinetic energy density in real space
!!  [with_nonlocal]=optional, true if non-local factors have to be loaded
!!
!! FUNCTION
!!  Setup of the spin-dependent part of the GS Hamiltonian.
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=Structured datatype initialization phase:
!!   * Quantities that depend spin are initialized.
!!
!! PARENTS
!!      d2frnl,dfpt_nstdy,dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho,energy,fock2ACE
!!      forstrnps,ks_ddiago,m_gkk,m_io_kss,m_phgamma,m_phpi,m_shirley,m_sigmaph
!!      nonlop_test,vtorho
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_spin_hamiltonian(Ham,isppol,vectornd,vlocal,vxctaulocal,with_nonlocal)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol
 logical,optional,intent(in) :: with_nonlocal
 class(gs_hamiltonian_type),intent(inout),target :: Ham
!arrays
 real(dp),optional,intent(in),target :: vectornd(:,:,:,:,:)
 real(dp),optional,intent(in),target :: vlocal(:,:,:,:),vxctaulocal(:,:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: jsppol

! *************************************************************************

 DBG_ENTER("COLL")

 !@gs_hamiltonian_type
 if (present(vlocal)) then
   ABI_CHECK(size(vlocal)==Ham%n4*Ham%n5*Ham%n6*Ham%nvloc,"Wrong vlocal")
   Ham%vlocal => vlocal
 end if
 if (present(vxctaulocal)) then
   ABI_CHECK(size(vxctaulocal)==Ham%n4*Ham%n5*Ham%n6*Ham%nvloc*4,"Wrong vxctaulocal")
   Ham%vxctaulocal => vxctaulocal
 end if
 if (present(vectornd)) then
   ABI_CHECK(size(vectornd)==Ham%n4*Ham%n5*Ham%n6*Ham%nvloc*3,"Wrong vectornd")
   Ham%vectornd => vectornd
 end if

 ! Retrieve non-local factors for this spin component
 if (present(with_nonlocal)) then
   if (with_nonlocal) then
     jsppol=min(isppol,size(Ham%ekb_spin,5))
     if (jsppol>0) Ham%ekb => Ham%ekb_spin(:,:,:,:,jsppol)
   end if
 end if

 ! Update enl and sij on GPU
#if defined HAVE_GPU_CUDA
 if (Ham%use_gpu_cuda==1) then
   call gpu_update_ham_data(Ham%ekb(:,:,:,1),size(Ham%ekb),Ham%sij,size(Ham%sij),Ham%gprimd,size(Ham%gprimd))
 end if
#endif

 DBG_EXIT("COLL")

end subroutine load_spin_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/destroy_rf_hamiltonian
!! NAME
!!  destroy_rf_hamiltonian
!!
!! FUNCTION
!!  Clean and destroy rf_hamiltonian_type datastructure
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=All dynamic memory defined in the structure is deallocated.
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,dfpt_vtorho,m_gkk,m_phgamma,m_phpi
!!      m_sigmaph
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine destroy_rf_hamiltonian(rf_Ham)

!Arguments ------------------------------------
!scalars
 class(rf_hamiltonian_type),intent(inout) :: rf_Ham

! *************************************************************************

 DBG_ENTER("COLL")

!@rf_hamiltonian_type

! Real arrays
 ABI_SFREE(rf_Ham%e1kbfr_spin)
 ABI_SFREE(rf_Ham%e1kbsc_spin)

! Real pointers
 if (associated(rf_Ham%dkinpw_k)) nullify(rf_Ham%dkinpw_k)
 if (associated(rf_Ham%dkinpw_kp)) nullify(rf_Ham%dkinpw_kp)
 if (associated(rf_Ham%ddkinpw_k)) nullify(rf_Ham%ddkinpw_k)
 if (associated(rf_Ham%ddkinpw_kp)) nullify(rf_Ham%ddkinpw_kp)
 if (associated(rf_Ham%vlocal1)) nullify(rf_Ham%vlocal1)
 if (associated(rf_Ham%e1kbfr)) nullify(rf_Ham%e1kbfr)
 if (associated(rf_Ham%e1kbsc)) nullify(rf_Ham%e1kbsc)

 DBG_EXIT("COLL")

end subroutine destroy_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_rf_hamiltonian
!! NAME
!!  init_rf_hamiltonian
!!
!! FUNCTION
!!  Creation method for the rf_hamiltonian_type structure.
!!  It allocates memory and initializes all quantities that do not depend on the k-point or spin.
!!
!! INPUTS
!!  [comm_atom]=optional, MPI communicator over atoms
!!  cplex_paw=1 if all on-site PAW quantities are real (GS), 2 if they are complex (RF)
!!  gs_Ham<gs_hamiltonian_type>=Structured datatype containing data for ground-state Hamiltonian at (k+q)
!!  [has_e1kbsc]=optional, true if rf_Ham%e1kbsc has to be initialized.
!!               e1kbsc contains the self-consistent 1st-order PAW Dij coefficients (depending on VHxc^(1))
!!  ipert=index of perturbation
!!  [mpi_atmtab(:)]=optional, indexes of the atoms treated by current proc
!!  [mpi_spintab(2)]=optional, flags defining the spin(s) treated be current process:
!!                    mpi_spintab(1)=1 if non-polarized or spin-up treated
!!                    mpi_spintab(2)=1 if polarized and spin-dn treated
!!  [paw_ij1(:)<paw_ij_type>]=Various 1st-order arrays given on (i,j) (partial waves)
!!                            channels (paw_ij1%dij and paw_ij1%difr only used here).
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=Structured datatype almost completely initialized:
!!   * Basic variables and dimensions are transfered to the structure.
!!   * All pointers are allocated with correct dimensions.
!!   * Quantities that do not depend on the k-point or spin are initialized.
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,dfpt_vtorho,m_gkk,m_phgamma,m_phpi
!!      m_sigmaph
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_rf_hamiltonian(cplex,gs_Ham,ipert,rf_Ham,&
&          comm_atom,mpi_atmtab,mpi_spintab,paw_ij1,has_e1kbsc) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ipert
 integer,intent(in),optional :: comm_atom
 logical,intent(in),optional :: has_e1kbsc
 type(gs_hamiltonian_type),intent(in) :: gs_Ham
 type(rf_hamiltonian_type),intent(inout),target :: rf_Ham
!arrays
 integer,optional,intent(in)  :: mpi_atmtab(:),mpi_spintab(2)
 type(paw_ij_type),optional,intent(in) :: paw_ij1(:)

!Local variables-------------------------------
!scalars
 integer :: cplex_dij1,isp,jsp,my_comm_atom,my_nsppol
 logical :: has_e1kbsc_
!arrays
 integer :: my_spintab(2)
 real(dp),allocatable,target :: e1kb_tmp(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!@rf_hamiltonian_type

!Manage optional parameters
 has_e1kbsc_=.false.;if (present(has_e1kbsc)) has_e1kbsc_=has_e1kbsc
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 my_spintab=0;my_spintab(1:gs_Ham%nsppol)=1;if(present(mpi_spintab)) my_spintab=mpi_spintab
 my_nsppol=count(my_spintab==1)

 rf_Ham%cplex    =cplex

 rf_Ham%n4       =gs_Ham%n4
 rf_Ham%n5       =gs_Ham%n5
 rf_Ham%n6       =gs_Ham%n6
 rf_Ham%nvloc    =gs_Ham%nvloc
 rf_Ham%nsppol   =gs_Ham%nsppol
 rf_Ham%nspinor  =gs_Ham%nspinor

 rf_Ham%dime1kb1=0
 rf_Ham%dime1kb2=gs_Ham%dimekb2
 if (gs_Ham%usepaw==1.and.ipert/=gs_Ham%natom+1.and.ipert/=gs_Ham%natom+10) then
   cplex_dij1=1;if ((gs_Ham%nspinor==2).or.any(abs(gs_Ham%nucdipmom)>tol8)) cplex_dij1=2
   rf_Ham%dime1kb1=cplex_dij1*(gs_Ham%lmnmax*(gs_Ham%lmnmax+1))/2
 end if

  ! Allocate the arrays of the 1st-order Hamiltonian
  ! We preload here 1st-order non-local factors in order to avoid
  ! a communication over atoms inside the loop over spins.
 if (gs_Ham%usepaw==1.and.rf_Ham%dime1kb1>0) then
   if ((ipert>=1.and.ipert<=gs_Ham%natom).or.ipert==gs_Ham%natom+2.or.&
        ipert==gs_Ham%natom+3.or.ipert==gs_Ham%natom+4.or.ipert==gs_Ham%natom+11) then

     ABI_ALLOCATE(rf_Ham%e1kbfr_spin,(rf_Ham%dime1kb1,rf_Ham%dime1kb2,rf_Ham%nspinor**2,cplex,my_nsppol))
     rf_Ham%e1kbfr_spin=zero
     if (has_e1kbsc_) then
       ABI_ALLOCATE(rf_Ham%e1kbsc_spin,(rf_Ham%dime1kb1,rf_Ham%dime1kb2,rf_Ham%nspinor**2,cplex,my_nsppol))
       rf_Ham%e1kbsc_spin=zero
     end if

     if (present(paw_ij1)) then

       if (my_nsppol<rf_Ham%nsppol) then
         ABI_ALLOCATE(e1kb_tmp,(rf_Ham%dime1kb1,rf_Ham%dime1kb2,rf_Ham%nspinor**2,cplex))
       end if

!      === Frozen term
       jsp=0
       do isp=1,rf_Ham%nsppol
         if (my_spintab(isp)==1) then
           jsp=jsp+1 ; rf_Ham%e1kbfr => rf_Ham%e1kbfr_spin(:,:,:,:,jsp)
         else
           rf_Ham%e1kbfr => e1kb_tmp
         end if
         if (present(mpi_atmtab)) then
           call pawdij2e1kb(paw_ij1,isp,my_comm_atom,e1kbfr=rf_Ham%e1kbfr,mpi_atmtab=mpi_atmtab)
         else
           call pawdij2e1kb(paw_ij1,isp,my_comm_atom,e1kbfr=rf_Ham%e1kbfr)
         end if
       end do

!      === Self-consistent term
       if (has_e1kbsc_) then
         jsp=0
         do isp=1,rf_Ham%nsppol
           if (my_spintab(isp)==1) then
             jsp=jsp+1 ; rf_Ham%e1kbsc => rf_Ham%e1kbsc_spin(:,:,:,:,jsp)
           else
             rf_Ham%e1kbsc => e1kb_tmp
           end if
           if (present(mpi_atmtab)) then
             call pawdij2e1kb(paw_ij1,isp,my_comm_atom,e1kbsc=rf_Ham%e1kbsc,mpi_atmtab=mpi_atmtab)
           else
             call pawdij2e1kb(paw_ij1,isp,my_comm_atom,e1kbsc=rf_Ham%e1kbsc)
           end if
         end do
       end if

       if (my_nsppol<rf_Ham%nsppol) then
         ABI_DEALLOCATE(e1kb_tmp)
       end if

     end if
   end if
 end if

 if (.not.allocated(rf_Ham%e1kbfr_spin)) then
   ABI_ALLOCATE(rf_Ham%e1kbfr_spin,(0,0,0,0,0))
 end if
 if (.not.allocated(rf_Ham%e1kbsc_spin)) then
   ABI_ALLOCATE(rf_Ham%e1kbsc_spin,(0,0,0,0,0))
 end if
 nullify(rf_Ham%e1kbfr)
 nullify(rf_Ham%e1kbsc)

 DBG_EXIT("COLL")

end subroutine init_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_spin_rf_hamiltonian
!! NAME
!!  load_spin_rf_hamiltonian
!!
!! FUNCTION
!!  Setup of the spin-dependent part of the 1st- and 2nd- order Hamiltonian.
!!
!! INPUTS
!!  isppol=index of current spin
!!  [vlocal1(cplex*n4,n5,n6,nvloc)]=optional, 1st-order local potential in real space
!!  [with_nonlocal]=optional, true if non-local factors have to be loaded
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=Structured datatype initialization phase:
!!   * Quantities that depend on spin are initialized.
!!
!! PARENTS
!!      dfpt_rhofermi,dfpt_vtorho,m_gkk,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_spin_rf_hamiltonian(rf_Ham,isppol,vlocal1,with_nonlocal)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol
 logical,optional,intent(in) :: with_nonlocal
 class(rf_hamiltonian_type),intent(inout),target :: rf_Ham
!arrays
 real(dp),optional,target,intent(in) :: vlocal1(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: jsppol

! *************************************************************************

 DBG_ENTER("COLL")

!@rf_hamiltonian_type

 if (present(vlocal1)) then
   ABI_CHECK(size(vlocal1)==rf_Ham%cplex*rf_Ham%n4*rf_Ham%n5*rf_Ham%n6*rf_Ham%nvloc,"Wrong vlocal1")
   rf_Ham%vlocal1 => vlocal1
 end if

 ! Retrieve non-local factors for this spin component
 if (present(with_nonlocal)) then
   if (with_nonlocal) then
     if (size(rf_Ham%e1kbfr_spin)>0) then
       jsppol=min(isppol,size(rf_Ham%e1kbfr_spin,5))
       if (jsppol>0) rf_Ham%e1kbfr => rf_Ham%e1kbfr_spin(:,:,:,:,jsppol)
     end if
     if (size(rf_Ham%e1kbsc_spin)>0) then
       jsppol=min(isppol,size(rf_Ham%e1kbsc_spin,5))
       if (jsppol>0) rf_Ham%e1kbsc => rf_Ham%e1kbsc_spin(:,:,:,:,jsppol)
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine load_spin_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_k_rf_hamiltonian
!! NAME
!!  load_k_rf_hamiltonian
!!
!! FUNCTION
!!  Setup of the k-dependent part of the 1st- and 2nd- order Hamiltonian
!!
!! INPUTS
!!  [dkinpw_k]=1st derivative of the (modified) kinetic energy for each plane wave
!!  [ddkinpw_k]=2nd derivative of the (modified) kinetic energy for each plane wave
!!  [npw_k]=number of plane waves
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=structured datatype completed with k-dependent quantitites.
!!          Quantities at k^prime are set equal to quantities at k.
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,getgh1c
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_k_rf_hamiltonian(rf_Ham,dkinpw_k,ddkinpw_k,npw_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: npw_k
 class(rf_hamiltonian_type),intent(inout),target :: rf_Ham
!arrays
 real(dp),intent(in),optional,target :: dkinpw_k(:),ddkinpw_k(:)

!Local variables-------------------------------

! *************************************************************************

 DBG_ENTER("COLL")

!@gs_hamiltonian_type

!k-dependent scalars
 if (present(npw_k)) then
   rf_Ham%npw_k  = npw_k
   rf_Ham%npw_kp = npw_k
 end if

!Pointers to k-dependent quantitites
 if (present(dkinpw_k)) then
   rf_Ham%dkinpw_k  => dkinpw_k
   rf_Ham%dkinpw_kp => dkinpw_k
 end if
 if (present(ddkinpw_k)) then
   rf_Ham%ddkinpw_k  => ddkinpw_k
   rf_Ham%ddkinpw_kp => ddkinpw_k
 end if

 DBG_EXIT("COLL")

end subroutine load_k_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/pawdij2ekb
!! NAME
!!  pawdij2ekb
!!
!! FUNCTION
!!  Transfer PAW Dij (on-site GS Hamiltonian) values
!!  from paw_ij datastructure to ekb array
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hamiltonian
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine pawdij2ekb(ekb,paw_ij,isppol,comm_atom,mpi_atmtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,comm_atom
!arrays
 integer,intent(in),optional,target :: mpi_atmtab(:)
 real(dp),intent(out) :: ekb(:,:,:,:)
 type(paw_ij_type),intent(in) :: paw_ij(:)

!Local variables-------------------------------
!scalars
 integer :: dimdij,dimekb1,dimekb3,dimekb4,iatom,iatom_tot,ierr,ii,isp,ispden,my_natom,natom,qphase
 logical :: my_atmtab_allocated,paral_atom
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

 DBG_ENTER("COLL")

 ekb=zero

!Set up parallelism over atoms
 natom=size(ekb,2); my_natom=size(paw_ij)
 paral_atom=(xmpi_comm_size(comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 call get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 !Retrieve PAW Dij coefficients for this spin component
 if (my_natom>0) then
   if (allocated(paw_ij(1)%dij)) then
     dimekb1=size(ekb,1) ; dimekb3=size(ekb,3) ; dimekb4=size(ekb,4)
     qphase=paw_ij(1)%qphase
     ABI_CHECK(qphase<=dimekb4,'paw_ij%qphase>dimekb4!')
     do ii=1,qphase
       do ispden=1,dimekb3
         isp=isppol; if (dimekb3==4) isp=ispden
         do iatom=1,my_natom
           iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
           dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
           ABI_CHECK(dimdij<=dimekb1,'Size of paw_ij%dij>dimekb1!')
           ekb(1:dimdij,iatom_tot,ispden,ii)=paw_ij(iatom)%dij(1+(ii-1)*dimdij:ii*dimdij,isp)
         end do
       end do
     end do
   end if
 end if

!Communication in case of distribution over atomic sites
 if (paral_atom) then
   call xmpi_sum(ekb,comm_atom,ierr)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawdij2ekb
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/pawdij2e1kb
!! NAME
!!  pawdij2e1kb
!!
!! FUNCTION
!!  Transfer PAW Dij (on-site RF Hamiltonian) values
!!  from paw_ij datastructure to e1kb array
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      d2frnl,dfpt_nstpaw,m_hamiltonian
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine pawdij2e1kb(paw_ij1,isppol,comm_atom,mpi_atmtab,e1kbfr,e1kbsc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,comm_atom
!arrays
 integer,intent(in),optional,target :: mpi_atmtab(:)
 real(dp),optional,intent(out) :: e1kbfr(:,:,:,:),e1kbsc(:,:,:,:)
 type(paw_ij_type),intent(in) :: paw_ij1(:)

!Local variables-------------------------------
!scalars
 integer :: dimdij1,dime1kb1,dime1kb3,dime1kb4,iatom,iatom_tot,ierr,isp,ispden
 integer :: my_natom,natom,qphase
 logical :: my_atmtab_allocated,paral_atom
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if ((.not.present(e1kbfr)).and.(.not.present(e1kbsc))) return
 if (present(e1kbfr)) then
   e1kbfr=zero ; natom=size(e1kbfr,2)
 end if
 if (present(e1kbsc)) then
   e1kbsc=zero ; natom=size(e1kbsc,2)
 end if

!Set up parallelism over atoms
 my_natom=size(paw_ij1) ; paral_atom=(xmpi_comm_size(comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 call get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Retrieve 1st-order PAW Dij coefficients for this spin component (frozen)
 if (my_natom>0.and.present(e1kbfr)) then
   if (allocated(paw_ij1(1)%dijfr)) then
     dime1kb1=size(e1kbfr,1) ; dime1kb3=size(e1kbfr,3) ; dime1kb4=size(e1kbfr,4)
     ABI_CHECK(paw_ij1(1)%qphase==dime1kb4,'BUG in pawdij2e1kb (1)!')
     do ispden=1,dime1kb3
       isp=isppol;if (dime1kb3==4) isp=ispden
       do iatom=1,my_natom
         iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
         qphase=paw_ij1(iatom)%qphase
         dimdij1=paw_ij1(iatom)%cplex_dij*paw_ij1(iatom)%lmn2_size
         ABI_CHECK(dimdij1<=dime1kb1,'BUG: size of paw_ij1%dij>dime1kb1!')
         e1kbfr(1:dimdij1,iatom_tot,ispden,1)=paw_ij1(iatom)%dijfr(1:dimdij1,isp)
         if (qphase==2) e1kbfr(1:dimdij1,iatom_tot,ispden,2)=paw_ij1(iatom)%dijfr(dimdij1+1:2*dimdij1,isp)
       end do
     end do
   end if
 end if

!Retrieve 1st-order PAW Dij coefficients for this spin component (self-consistent)
 if (my_natom>0.and.present(e1kbsc)) then
   if (allocated(paw_ij1(1)%dijfr).and.allocated(paw_ij1(1)%dij)) then
     dime1kb1=size(e1kbsc,1) ; dime1kb3=size(e1kbsc,3) ; dime1kb4=size(e1kbsc,4)
     ABI_CHECK(paw_ij1(1)%qphase==dime1kb4,'BUG in pawdij2e1kb (1)!')
     do ispden=1,dime1kb3
       isp=isppol;if (dime1kb3==4) isp=ispden
       do iatom=1,my_natom
         iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
         qphase=paw_ij1(iatom)%qphase
         dimdij1=paw_ij1(iatom)%cplex_dij*paw_ij1(iatom)%lmn2_size
         ABI_CHECK(dimdij1<=dime1kb1,'BUG: size of paw_ij1%dij>dime1kb1!')
         e1kbsc(1:dimdij1,iatom_tot,ispden,1)=paw_ij1(iatom)%dij  (1:dimdij1,isp) &
&                                            -paw_ij1(iatom)%dijfr(1:dimdij1,isp)
         if (qphase==2) e1kbsc(1:dimdij1,iatom_tot,ispden,2)=paw_ij1(iatom)%dij  (dimdij1+1:2*dimdij1,isp) &
&                                                           -paw_ij1(iatom)%dijfr(dimdij1+1:2*dimdij1,isp)
       end do
     end do
   end if
 end if

 ! Communication in case of distribution over atomic sites
 if (paral_atom) then
   if (present(e1kbfr)) then
     call xmpi_sum(e1kbfr,comm_atom,ierr)
   end if
   if (present(e1kbsc)) then
     call xmpi_sum(e1kbsc,comm_atom,ierr)
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawdij2e1kb
!!***

END MODULE m_hamiltonian
!!***
