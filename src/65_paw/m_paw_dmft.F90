!****m* ABINIT/m_paw_dmft
!! NAME
!!  m_paw_dmft
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_dmft

 use defs_basis
 use m_abicore
 use m_CtqmcInterface
 use m_data4entropyDMFT
 use m_dtset
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_crystal, only : crystal_t
 use m_fstrings, only : int2char4
 use m_geometry, only : symredcart
 use m_io_tools, only : open_file
 use m_mpinfo, only : proc_distrb_cycle
 use m_paw_numeric, only : paw_jbessel_4spline
 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_free,pawrad_init,pawrad_type,simp_gen
 use m_pawtab, only : pawtab_type

 implicit none

 private

 public :: init_dmft
 public :: init_sc_dmft
 public :: construct_nwli_dmft
 public :: destroy_dmft
 public :: destroy_sc_dmft
 public :: print_dmft
 public :: print_sc_dmft
 public :: saveocc_dmft
 public :: readocc_dmft

!!***

!----------------------------------------------------------------------

!!****t* m_paw_dmft/mpi_distrib_dmft_type
!! NAME
!!  mpi_distrib_dmft_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data for the MPI
!!  parallelization over frequencies and kpts in DFT+DMFT.
!!
!! SOURCE

 type, public :: mpi_distrib_dmft_type

   ! Two types of parallelization
   ! Type 1: parallelization over kpt and then frequencies
   ! Type 2: parallelization over frequencies only

   integer :: comm_freq
   ! Frequency communicator (type 1)

   integer :: comm_kpt
   ! Kpt communicator (type 1)

   integer :: me_freq
   ! Rank in comm_freq (type 1)

   integer :: me_kpt
   ! Rank in comm_kpt (type 1)

   integer :: nw
   ! Number of frequencies (either nwlo or nwr, same for both types)

   integer :: shiftk
   ! Shift from kpt index on the current CPU to the physical index (type 1)

   integer, allocatable :: nkpt_mem(:)
   ! Number of kpt handled by each CPU of the kpt communicator (type 1)

   integer, allocatable :: nw_mem(:)
   ! Number of frequencies handled by each CPU of the global communicator (type 2)

   integer, allocatable :: nw_mem_kptparal(:)
   ! Number of frequencies handled by each CPU of the frequency communicator (type 1)

   integer, allocatable :: procb(:)
   ! Rank in comm_kpt of the CPU handling each kpt (type 1)

   integer, allocatable :: procf(:)
   ! Rank in the global communicator of the CPU handling each frequency (type 2)

   integer, allocatable :: proct(:)
   ! Rank in comm_freq of the CPU handling each frequency (type 1)

 end type mpi_distrib_dmft_type
!!***

!----------------------------------------------------------------------

!!****t* m_paw_dmft/paw_dmft_type
!! NAME
!!  paw_dmft_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data for the link
!!  between dmft and paw.
!!  occnd(non-diagonal band occupations for self-consistency), band_in
!!  (say which band are taken into account in the calculation), and the
!   dimensions of these arrays.
!!
!! SOURCE

 type, public :: paw_dmft_type

  integer :: dmft_blockdiag
  ! Block diagonalize Hamiltonian in the local basis

  integer :: dmft_dc
  ! Type of double counting used in DMFT

  integer :: dmft_entropy
  ! = 0: do not compute entropy
  ! >= 1: compute entropy with an integration over self-consistent calculations

  integer :: dmft_iter
  ! Nb of iterations for DMFT self-consistent cycle.

  integer :: dmft_kspectralfunc
  ! =0 Default
  ! =1 Activate calculation of k-resolved spectral function

  integer :: dmft_log_freq
  ! = 0: do not use log frequencies
  ! = 1: use log frequencies

  integer :: dmft_magnfield
  ! = 0: do nothing
  ! = 1: apply a magnetic field Bz via Zeeman Hamiltonian on Kohn-Sham energies
  ! = 2: apply a magnetic field Bz via Zeeman Hamiltonian on local impurity Hamiltonian

  integer :: dmft_nwli
  ! Physical index of the last imaginary frequency (/=dmft_nwlo when dmft_log_freq=1)

  integer :: dmft_nwlo
  ! Number of imaginary frequencies

  integer :: dmft_nwr
  ! Number of real frequencies

  integer :: dmft_prgn
  ! Specify the way of printing the green function.
  !  =1   print green
  !  =2   print self

  integer :: dmft_prt_maxent
  ! =1 to print Maxent files

  integer :: dmft_prtself
  ! =1 to keep self-energy files of all previous iterations

  integer :: dmft_prtwan
  ! =1 to print Wannier functions

  integer :: dmft_read_occnd
  ! Flag to read/write DMFT occupations
  ! =0 : Occupations are written but not read
  ! =1 : Occupations are read from I_DMFTOCCND, with I the root for input file
  ! =2 : Occupations are read from O_DMFTOCCND, with O the root for input file

  integer :: dmft_rslf
  ! Flag to read the self-energy at each iteration
  ! =-1 : Self-energy is set to 0
  ! =0 : Self-energy is set to double counting
  ! =1 : Self-energy is read from file

  integer :: dmft_solv
  ! Choice of solver for DMFT

  integer :: dmft_t2g
  ! Only use t2g orbitals

  integer :: dmft_triqs_compute_integral
  ! Only relevant when dmft_triqs_entropy=1.
  ! =1: Compute the impurity entropy by thermodynamic integration over the impurity models.
  ! =0: Do not compute the impurity entropy. All the other contributions to the free
  ! energy are still computed.

  integer :: dmft_triqs_det_init_size
  ! TRIQS CTQMC: Initial size of the hybridization matrix. If it is too low,
  ! the matrix will be resized very often, which can be slow.

  integer :: dmft_triqs_det_n_operations_before_check
  ! TRIQS CTQMC: Number of operations before check of the hybridization matrix.
  ! If it is low, the matrix will be checked too often, which can be slow.

  integer :: dmft_triqs_entropy
  ! TRIQS CTQMC: Compute the DMFT entropy.

  integer :: dmft_triqs_gaussorder
  ! Order of the Gauss-Legendre quadrature for each subdivision of the thermodynamic integration.

  integer :: dmft_triqs_loc_n_min
  ! TRIQS CTQMC: Only configurations with a number of electrons in
  ! [nlocmin,nlocmax] are taken into account.

  integer :: dmft_triqs_loc_n_max
  ! TRIQS CTQMC: Only configurations with a number of electrons in
  ! [nlocmin,nlocmax] are taken into account.

  integer :: dmft_triqs_nleg
  ! TRIQS CTQMC: Nb of Legendre polynomials used for the
  ! Green's function (Phys. Rev. B 84, 075145) [[cite:Boehnke2011]].

  integer :: dmft_triqs_nsubdivisions
  ! Number of regular subdivisions of the interval [0,U], each of which
  ! containing dmft_triqs_gaussorder points

  integer :: dmft_triqs_read_ctqmcdata
  ! TRIQS CTQMC: Read CTQMC data of the previous iteration

  integer :: dmft_triqs_seed_a
  ! TRIQS CTQMC: The CTQMC seed is seed_a + rank * seed_b.

  integer :: dmft_triqs_seed_b
  ! TRIQS CTQMC: The CTQMC seed is seed_a + rank * seed_b.

  integer :: dmft_triqs_therm_restart
  ! TRIQS CTQMC: Number of thermalization steps when we restart from a previous configuration.

  integer :: dmft_wanorthnorm
  ! =2 orthonormalization of Wannier functions for each k-point
  ! =3 orthonormalization over the sum over k-points

  integer :: dmft_x2my2d
  ! Only use x2my2d orbital

  integer :: dmftbandf
  ! Highest correlated band

  integer :: dmftbandi
  ! Lowest correlated band

  integer :: dmftcheck
  ! Check various part of the implementation

  integer :: dmftctqmc_basis
  ! Basis in which to perform the CTQMC calculation
  ! 0 : Slm basis, 1 : diagonalize local Hamiltonian, 2: diagonalize the density matrix
  ! Only for TRIQS: 3: Ylm, 4: JmJ

  integer :: dmftctqmc_check
  ! ABINIT CTQMC: perform a check on the impurity and/or bath operator
  ! only for debug
  ! 0 : nothing, 1 : impurity, 2 : bath, 3 : both

  integer :: dmftctqmc_correl
  ! ABINIT CTQMC: Gives analysis for CTQMC
  ! 0 : nothing, 1 : activated Correlations.dat

  integer :: dmftctqmc_gmove
  ! ABINIT CTQMC: add global move every dmftctqmc_gmove sweeps
  ! >= 0 ; done inside CT-QMC with warning
  ! == 0 ; no global moves

  integer :: dmftctqmc_grnns
  ! ABINIT CTQMC: compute green function noise for each imaginary time
  ! 0 : nothing, 1 : activated

  integer :: dmftctqmc_localprop
  ! ABINIT CTQMC: local properties calculations
  ! 0 : nothing, 1 : Histogram, 2 : magnetic susceptibility, 3 : charge susceptibility

  integer :: dmftctqmc_meas
  ! ABINIT/TRIQS CTQMC: measurements are done every dmftctqmc_meas step

  integer :: dmftctqmc_mov
  ! ABINIT CTQMC: Gives movie for CTQMC
  ! 0 : nothing, 1 : 1 file Movie_RANK.tex for each cpu

  integer :: dmftctqmc_mrka
  ! ABINIT CTQMC: Write a temporary file Spectra_RANK.dat with the sweep evolution of
  ! the number of electron for each flavor
  ! The measurement is done every dmftctqmc_meas*dmftctqmc_mrka sweep
  ! e.g. : meas=2 mrka=10 -> every 20 sweeps sum_i c+(ti)c(t'i) is measured

  integer :: dmftctqmc_order
  ! ABINIT CTQMC: Gives perturbation order of CTQMC solver
  ! 0 : nothing, >=1 max order evaluated in Perturbation.dat

  integer :: dmftqmc_l
  ! Number of points on the imaginary time grid for G(tau) and Delta(tau)

!  integer :: dmft_mag
!  ! 0 if non magnetic calculation, 1 if magnetic calculation

  integer :: dmftqmc_seed
  ! Seed for CTQMC (only for ABINIT)

  integer :: dmftqmc_therm
  ! Number of thermalization steps for CTQMC (only for ABINIT, and for TRIQS when we don't restart from a previous configuration)

  integer :: gpu_option
  ! Wether to use GPU implementation (expected values: ABI_GPU_DISABLED, ABI_GPU_OPENMP)

  integer :: idmftloop
  ! Current iteration in the DFT+DMFT loop

  integer :: ientropy
  ! activate evaluation of terms for alternative calculation of entropy in DMFT

  integer :: ireadctqmcdata
  ! Internal flag to indicate if an input CTQMC_DATA file must be read

  integer :: ireadself
  ! Internal flag to indicate if an input self file must be read

  integer :: kptopt
  ! Option to generate kpts

  integer :: lchipsiortho
  ! Internal flag
  ! =0 <Chi|Psi> is not orthonormalized
  ! =1 <Chi|Psi> is orthonormalized

  integer :: maxlpawu
  ! Maximal correlated l over all atoms

  integer :: maxmeshsize
  ! Maximal size of the radial mesh over all atoms

  integer :: maxnproju
  ! Maximal number of correlated projectors over all atoms

  integer :: mband
  ! Total number of bands

  integer :: mbandc
  ! Total number of correlated bands

  integer :: mkmem
  ! Number of k-points handled by the current process within the DFT
  ! parallelization scheme

  integer :: myproc
  ! Rank in the global communicator

  integer :: natom
  ! Number of atoms

  integer :: natpawu
  ! Number of correlated atoms

  integer :: nkpt
  ! Number of k-points in the IBZ.

  !integer :: nspden
  ! Number of spin densities

  integer :: nproc
  ! Total number of MPI processes

  integer :: nspinor
  ! Number of spinor components

  integer :: nsppol
  ! Number of spin polarizations

  integer :: nsym
  ! Number of symmetries

  integer :: ntypat
  ! Number of atom types

  integer :: prtdos
  ! Print DOS when >=1

  integer :: prtvol
  ! Flag for different print options

  integer :: spacecomm
  ! MPI_COMM_WORLD

  integer :: unpaw
  ! File number for cprj

  integer :: use_dmft
  ! 1 if non diagonal occupations are used, else 0

  integer :: use_fixed_self
  ! Impose a fixed self-energy during the first use_fixed_self iterations

  integer :: use_sc_dmft
  ! 1 for charge-self consistent calculations

  logical :: dmft_triqs_leg_measure
  ! TRIQS CTQMC: Flag to activate Legendre measurement

  logical :: dmft_triqs_measure_density_matrix
  ! TRIQS CTQMC: Flag to activate the measurement of the density matrix

  logical :: dmft_triqs_move_double
  ! TRIQS CTQMC: Flag to activate the double moves

  logical :: dmft_triqs_move_shift
  ! TRIQS CTQMC: Flag to activate the shift move

  logical :: dmft_triqs_off_diag
  ! TRIQS CTQMC: Flag to sample the off-diagonal elements of the Green's function

  logical :: dmft_triqs_time_invariance
  ! TRIQS CTQMC: Flag to activate the use of time invariance for the sampling
  ! of the density matrix

  logical :: dmft_triqs_use_norm_as_weight
  ! TRIQS CTQMC: Flag to activate the use of the norm of the matrix as weight
  ! instead of the trace

  real(dp) :: dmft_charge_prec
  ! Precision on charge required for determination of fermi level (fermi_green)

  real(dp) :: dmft_fermi_prec
  ! Required precision on Fermi level (fermi_green) during the DMFT SCF cycle, (=> ifermie_cv)
  ! used also for self (new_self)  (=> iself_cv).

  real(dp) :: dmft_fermi_step
  ! When dmft_optim = 0, step increment to find the upper and lower bounds of the Fermi level
  ! When dmft_optim = 1, maximal step size in the Fermi level search

  real(dp) :: dmft_lcpr
  ! Required precision on local correlated charge in order to stop SCF
  ! DMFT cycle (integrate_green) => ichargeloc_cv

  real(dp) :: dmft_magnfield_b
  ! Value of the applied magnetic field in Tesla

  real(dp) :: dmft_mxsf
  ! Mixing coefficient for Self-Energy during the SCF DMFT cycle.

  real(dp) :: dmft_tolfreq
  ! Required precision on local correlated density matrix (depends on
  ! frequency mesh), used in m_dmft/dmft_solve

  real(dp) :: dmft_triqs_det_precision_error
  ! TRIQS CTQMC: Error threshold for the deviation of the determinant when a check is performed.

  real(dp) :: dmft_triqs_det_precision_warning
  ! TRIQS CTQMC: Warning threshold for the deviation of the determinant when a check is performed.

  real(dp) :: dmft_triqs_det_singular_threshold
  ! TRIQS CTQMC: Threshold when checking if the determinant is singular.

  real(dp) :: dmft_triqs_epsilon
  ! TRIQS CTQMC: Threshold for singular values of the kernel matrix for the DLR fit

  real(dp) :: dmft_triqs_imag_threshold
  ! TRIQS CTQMC: Threshold for the imaginary part of Delta(tau)

  real(dp) :: dmft_triqs_lambda
  ! TRIQS CTQMC: Cutoff for the real frequency grid for the DLR fit

  real(dp) :: dmft_triqs_pauli_prob
  ! TRIQS CTQMC: Probability for proposing Pauli-aware insert and remove

  real(dp) :: dmft_triqs_tol_block
  ! TRIQS CTQMC: Off-diagonal elements below this threshold are set to 0

  real(dp) :: dmft_wanrad
  ! Maximal radius for print of the Wannier functions

  real(dp) :: dmftqmc_n
  ! ABINIT CTQMC: Nb of sweeps
  ! TRIQS CTQMC: Nb of measurements

  real(dp) :: e_dc
  ! Double counting energy

  real(dp) :: e_hu
  ! Interaction energy

  real(dp) :: fermie
  ! DMFT Fermi level

  real(dp) :: fermie_dft
  ! DFT Fermi level

  real(dp) :: j_for_s
  ! Variable for evaluation of correlation energy for U=0 in the entropic
  ! calculation

  real(dp) :: nelectval
  ! Number of valence electrons

  real(dp) :: sdmft
  ! DFT+DMFT total entropy

  real(dp) :: simp
  ! DFT+DMFT entropy of the impurity electrons

  real(dp) :: temp
  ! Temperature (Ha)

  real(dp) :: u_for_s
  ! Variable for evaluation of correlation energy for U=0 in the entropic
  ! calculation

  character(len=fnlen) :: filapp
  ! Output file name

  character(len=fnlen) :: filctqmcdatain
  ! Input file name for CTQMC_DATA file

  character(len=fnlen) :: filnamei
  ! Input file name

  character(len=fnlen) :: filselfin
  ! Input file name for self file

  integer, allocatable :: bandc_proc(:)
  ! Proc index (on comm_band) for each correlated band in DMFT (for kgb paral)

  integer, allocatable :: exclude_bands(:)
  ! Gives the bands than are not in the DMFT calculations.

  integer, allocatable :: include_bands(:)
  ! For each bands included in the calculation (1..mbandc), include_bands
  ! gives the index in the full band index  (1...mband)

  integer, allocatable :: lpawu(:)
  ! Correlated l for each atom (set to -1 if not correlated)

  integer, allocatable :: siz_proj(:)
  ! Size of the radial mesh for the DMFT orbital, for each atom type.

  logical, allocatable :: band_in(:)
  ! True for each band included in the calculation

  logical, allocatable :: use_bandc(:)
  ! True for each proc wich has at least one band involved in DMFT non diagonal
  ! occupations on band parallelism

  real(dp), allocatable :: edc(:)
  ! Double counting energy for each atom (only used as a temporary for dmft_dc=8)

  real(dp), allocatable :: edcdc(:)
  ! Integral of Vdc * rho for each atom (only used as a temporary for dmft_dc=8)

  real(dp), allocatable :: eigen_dft(:,:,:)
  ! DFT eigenvalues for each correlated band, k-point, polarization

  real(dp), allocatable :: occnd(:,:,:,:,:)
  ! Non diagonal band-occupation for each k-point, polarization.

  real(dp), allocatable :: omega_lo(:)
  ! Imaginary frequencies

  real(dp), allocatable :: omega_r(:)
  ! Real frequencies

  real(dp), allocatable :: phi_int(:,:)
  ! Integral of <Chi|Phi> for every correlated projector and atom type

  real(dp), allocatable :: phimtphi(:,:,:)
  ! Phi-Phi_tilde for every r,correlated projector and atom type

  real(dp), allocatable :: phimtphi_int(:,:)
  ! Integral of <Chi|Phi-Phi_tilde> for every correlated projector and atom type

  real(dp), allocatable :: symrec_cart(:,:,:)
  ! Symmetries in cartesian coordinates

  real(dp), allocatable :: wgt_wlo(:)
  ! Weight of the imaginary frequencies

  real(dp), allocatable :: ylm(:,:,:,:)
  ! Ylm(k+G) for each G,m,l,k

!  real(dp), allocatable :: phi0phiiint(:)
!  ! non diagonal band-occupation for each k-point, polarisation.

  complex(dpc), allocatable :: bessel(:,:,:,:)
  ! 4*pi*(i**l)*jl(|k+G|r)*r/sqrt(ucvol) for each G,r,atom type and kpt

  complex(dpc), allocatable :: bessel_int(:,:,:)
  ! Integral over r of bessel(ig,:,iat,ikpt)

  complex(dpc), allocatable :: buf_psi(:)
  ! Temporary buffer for the computation of Wannier function

  complex(dpc), allocatable :: chipsi(:,:,:,:,:)
  ! Hermitian product <Chi|Psi> for each flavor, correlated band,
  ! k-point, polarization and atom

  complex(dpc), allocatable :: dpro(:,:,:)
  ! Exp(i(k+G).xred(iatom)) for each G,correlated atom and k

  complex(dpc), allocatable :: jmj2ylm(:,:,:)
  ! Transformation matrix from JmJ to Ylm basis for each lpawu

  complex(dpc), allocatable :: slm2ylm(:,:,:)
  ! Transformation matrix from real to complex harmonics for each lpawu

  complex(dpc), allocatable :: wannier(:,:,:)
  ! Wannier functions for each r,flavor and atom

  complex(dpc), allocatable :: zarot(:,:,:,:)
  !  Coeffs of the transformation of real spherical
  !  harmonics under the symmetry operations symrec.

  integer, ABI_CONTIGUOUS pointer :: dmft_nominal(:) => null()
  ! Only relevant when dmft_dc=7. Nominal occupancies for each atom.

  integer, pointer :: indsym(:,:) => null()
  ! Label of atom into which iatom is sent by the INVERSE of the
  ! symmetry operation symrel(isym)

  integer, pointer :: int_meshsz(:) => null()
  ! PAW integration radius for each atom type

  integer, ABI_CONTIGUOUS pointer :: nband(:) => null()
  ! Number of bands for each k-point and polarization

  integer, ABI_CONTIGUOUS pointer :: npwarr(:) => null()
  ! Number of plane waves on current process for each k-point

  integer, ABI_CONTIGUOUS pointer :: typat(:) => null()
  ! Type of each atom

  real(dp), ABI_CONTIGUOUS pointer :: dmft_shiftself(:) => null()
  ! Initial shift of the self-energy for each atom

  real(dp), ABI_CONTIGUOUS pointer :: eigen(:) => null()
  ! DFT eigenvalues

  real(dp), pointer :: fixed_self(:,:,:,:) => null()
  ! Fixed self-energy (only used when use_fixed_self > 0)

  real(dp), ABI_CONTIGUOUS pointer :: wtk(:) => null()
  ! Weights for each k-point

  type(CtqmcInterface), allocatable :: hybrid(:)

  type(data4entropyDMFT_t) :: forentropyDMFT

  type(pawrad_type), allocatable :: radgrid(:)
  ! Radial grid for each type of atom

  type(mpi_distrib_dmft_type) :: distrib
  ! MPI parallelization for imaginary frequencies

  type(mpi_distrib_dmft_type) :: distrib_r
  ! MPI parallelization for real frequencies

 end type paw_dmft_type
!!***

!----------------------------------------------------------------------

CONTAINS  !========================================================================================
!!***

!!****f* m_paw_dmft/init_sc_dmft
!! NAME
!! init_sc_dmft
!!
!! FUNCTION
!!  Allocate variables used in type paw_dmft_type.
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables for this dataset
!! mpsang = highest angular momentum + 1
!! paw_dmft <type(paw_dmft_type)>= paw+dmft related data
!! gprimd(3,3) = dimensional reciprocal space primitive translations
!! kg(3,mpw*mkmem) = reduced planewave coordinates.
!! mpi_enreg = information about MPI parallelization
!! npwarr(nkpt) = number of planewaves in basis at this k point
!! occ = DFT occupations
!! pawang <type(pawang)>=paw angular mesh and related data
!! pawrad <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rprimd(3,3) = dimensional primitive vectors
!! ucvol = unit cell volume in bohr**3.
!! unpaw = unit number for temporary PAW data
!! use_sc_dmft = for charge self-consistent calculations
!! xred(3,natom) = reduced dimensionless atomic coordinates
!! ylm(mpw*mkmem,mpsang*mpsang*useylm) = real spherical harmonics for each G and k point
!!
!! OUTPUTS
!! paw_dmft = datastructure for dmft
!!
!! SOURCE

subroutine init_sc_dmft(dtset,mpsang,paw_dmft,gprimd,kg,mpi_enreg,npwarr,occ,pawang, &
                      & pawrad,pawtab,rprimd,ucvol,unpaw,use_sc_dmft,xred,ylm)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: mpsang
 integer, optional, intent(in) :: unpaw,use_sc_dmft
 real(dp), optional, intent(in) :: ucvol
!type
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(MPI_type), optional, intent(in) :: mpi_enreg
 type(dataset_type), target, intent(in) :: dtset
 type(pawtab_type), optional, intent(in) :: pawtab(dtset%ntypat)
 type(pawang_type), optional, intent(in) :: pawang
 type(pawrad_type), target, optional, intent(in) :: pawrad(dtset%ntypat)
! arrays
 real(dp), optional, intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 integer, target, optional, intent(in) :: npwarr(dtset%nkpt)
 integer, optional, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 real(dp), optional, intent(in) :: gprimd(3,3),rprimd(3,3),xred(3,dtset%natom)
 real(dp), optional, intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang)
!Local variables ------------------------------------
 integer :: bdtot_index,dmft_dc,dmft_solv,dmftbandi,dmftbandf,fac,i,iatom
 integer :: iatom1,iband,icb,ig,ik,ikg,ikpt,im,im1,indproj,iproj,ir,isppol
 integer :: itypat,jc1,jj,jm,ll,lpawu,lpawu1,maxlpawu,mband,mbandc,mesh_size,mesh_type
 integer :: mkmem,ml1,mm,mpw,ms1,myproc,natom,nband_k,ndim,nkpt,nproc,nproju,npw
 integer :: nspinor,nsppol,nsym,ntypat,off_diag,siz_paw,siz_proj,siz_wan,use_dmft
 logical :: t2g,use_full_chipsi,verif,x2my2d
 real(dp) :: bes,besp,invsqrt2lp1,lstep,norm,onem,rad,rint,rstep,sumwtk,xj,xmj
 integer, parameter :: mt2g(3) = (/1,2,4/)
 integer, allocatable :: ind_msml(:,:)
 logical, allocatable :: lcycle(:),typcycle(:)
 real(dp), allocatable :: rmax(:),kpg(:,:),kpg_norm(:)
 character(len=500) :: dc_string,lda_string,message
!************************************************************************

 mband   = dtset%mband
 mkmem   = dtset%mkmem
 natom   = dtset%natom
 nkpt    = dtset%nkpt
 nspinor = dtset%nspinor
 nsppol  = dtset%nsppol
 nsym    = dtset%nsym
 ntypat  = dtset%ntypat

 dmftbandi = dtset%dmftbandi
 dmftbandf = dtset%dmftbandf
 dmft_dc   = dtset%dmft_dc
 dmft_solv = dtset%dmft_solv
 off_diag  = dtset%dmft_triqs_off_diag
 use_dmft  = abs(dtset%usedmft)
 paw_dmft%use_dmft = use_dmft
 paw_dmft%use_sc_dmft = 0

 paw_dmft%dmftbandf = dmftbandf
 paw_dmft%dmftbandi = dmftbandi
 paw_dmft%mband     = mband
 paw_dmft%mkmem     = mkmem
 paw_dmft%nkpt      = nkpt
 paw_dmft%nsym      = nsym
 paw_dmft%ntypat    = ntypat

 ! Spin related variables
 paw_dmft%nsppol    = nsppol
 paw_dmft%nspinor   = nspinor
 paw_dmft%idmftloop = 0
 paw_dmft%mbandc    = 0
 !paw_dmft%nspden      = nspden

 paw_dmft%dmft_read_occnd = dtset%dmft_read_occnd

 if(use_dmft == 10) then
   ABI_MALLOC(paw_dmft%occnd,(2,mband,mband,nkpt,nsppol*1))
   ABI_MALLOC(paw_dmft%band_in,(mband*1))
   ABI_MALLOC(paw_dmft%include_bands,((dmftbandf-dmftbandi+1)*1))
   ABI_MALLOC(paw_dmft%exclude_bands,(mband*1))
 else
   ABI_MALLOC(paw_dmft%occnd,(2,mband,mband,nkpt,nsppol*use_dmft))
   ABI_MALLOC(paw_dmft%band_in,(mband*use_dmft))
   ABI_MALLOC(paw_dmft%include_bands,((dmftbandf-dmftbandi+1)*use_dmft))
   ABI_MALLOC(paw_dmft%exclude_bands,(mband*use_dmft))
 endif

 if (use_dmft == 0) return

 ! In the case where no kpt is treated by the current CPU (ie sum(isppoltab)=0),
 ! dtset%mkmem is set to nkpt by convention in Abinit. We set it to its true value 0.
 if (sum(mpi_enreg%my_isppoltab(1:nsppol)) == 0) mkmem = 0

 ! Check processors for DMFT
 ! Initialize spaceComm, myproc, and nproc
 !spacecomm=mpi_enreg%comm_cell
 !myproc=mpi_enreg%me_cell
 !nproc=mpi_enreg%nproc_cell
 !spacecomm = mpi_enreg%comm_world
 myproc = mpi_enreg%me
 nproc  = mpi_enreg%nproc
 !print *, " spacecomm,myproc,nproc",spacecomm,myproc,nproc
 paw_dmft%spacecomm = mpi_enreg%comm_world
 paw_dmft%myproc    = myproc
 paw_dmft%nproc     = nproc

 paw_dmft%unpaw = unpaw

 if (dtset%nbandkss == 0) paw_dmft%use_sc_dmft = use_sc_dmft

 ! Do not comment these lines: it guarantees the parallelism in DMFT/QMC will work.
 if (xmpi_comm_size(xmpi_world) /= xmpi_comm_size(mpi_enreg%comm_world)) &
   & ABI_ERROR("Someone changed the k-point parallelism again")

 if (dmft_solv == 0) then
   do itypat=1,ntypat
     if (pawtab(itypat)%lpawu == -1) cycle
     if ((pawtab(itypat)%upawu > tol5) .or. (pawtab(itypat)%jpawu > tol5)) then
       write(message,'(2a,i5,2a,2e15.6)') ch10,&
         & ' option dmft_solv=0 requires upaw=jpaw=0 for species',itypat,ch10,&
         & ' Value of upawu and jpawu are here',pawtab(itypat)%upawu,pawtab(itypat)%jpawu
       ABI_ERROR(message)
     end if
   end do ! itypat
 end if ! dmft_solv=0

! todo_ab: why upaw and jpawu are not zero (on bigmac) if lpawu==-1 ?
! if(paw_dmft%dmft_solv==0.and.&
!& (maxval(abs(pawtab(:)%upawu))>tol5.or.maxval(abs(pawtab(:)%jpawu))>tol5))
!then
!   write(message, '(a,a,2f12.3)' )ch10,&
!&   ' option dmft_solv=0 requires
!upaw=jpaw=0',maxval(abs(pawtab(:)%upawu)),maxval(abs(pawtab(:)%jpawu))
!    ABI_WARNING(message)
! endif

 paw_dmft%dmftcheck = dtset%dmftcheck

 write(message,'(2a,i4)') ch10,'-       ( number of procs used in dmft ) = ',nproc
 call wrtout([std_out,ab_out],message,'COLL')
 write(std_out_default,'(2a,i4)') ch10,'       ( current proc is        ) = ',myproc
  ! write(ab_out_default,'(2a,i3)') ch10,'       ( current proc is        ) =', myproc
 if (myproc == nproc-1) write(std_out_default,'(2a,i4)') ch10,'      ( last proc            ) = ',myproc
  !   write(ab_out_default,'(2a,i3)') ch10,'       ( last proc            ) =', myproc

!#ifdef HAVE_MPI
! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ie)
! write(6,*) "nprocs,nb_procs",nproc,nb_procs
! if(nb_procs/=nproc)  then
!   message = ' Number of procs used in DMFT is erroneously computed '
!   ABI_ERROR(message)
! endif
!#endif

!=============================
!==  Associate pointers
!=============================

 paw_dmft%nband => dtset%nband(:)
 paw_dmft%dmft_shiftself => dtset%dmft_shiftself(:)
 paw_dmft%dmft_nominal => dtset%dmft_nominal(:)
 paw_dmft%npwarr => npwarr(:)

! TODO: Make it work for usedmdft = -1 (interface with Wannier90 needs spinor
! generalization)
 if(nspinor==2.and.dtset%nspden==1.and.use_dmft==10) then
   message = ' nspinor==2 and nspden=1 and usedmft=10 is not implemented yet'
   ABI_ERROR(message)
 endif

 paw_dmft%band_in(:) = .false.
 paw_dmft%occnd(:,:,:,:,:) = zero
 paw_dmft%use_dmft    = use_dmft

 ! if (bandkss/=0) then
 !   paw_dmft%use_sc_dmft = 0
 ! else
 !   paw_dmft%use_sc_dmft = use_sc_dmft
 ! endif
 ! paw_dmft%dmft_read_occnd = dmft_read_occnd
 ! paw_dmft%idmftloop=0
 ! paw_dmft%mbandc  = 0

 icb = 0
 mbandc = 0
 do iband=1,mband
  if (iband >= dmftbandi .and. iband <= dmftbandf) then
   paw_dmft%band_in(iband)=.true.
   mbandc = mbandc + 1
   paw_dmft%include_bands(mbandc) = iband
  else
    icb = icb + 1
    paw_dmft%exclude_bands(icb) = iband
  end if ! band>=bandi and band<=bandf
 end do ! iband
 paw_dmft%mbandc = mbandc

 bdtot_index = 1
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
     do iband=1,nband_k
       paw_dmft%occnd(1,iband,iband,ikpt,isppol) = occ(bdtot_index)
       bdtot_index = bdtot_index + 1
     end do ! iband
   end do ! ikpt
 end do ! isppol

 if (paw_dmft%use_sc_dmft /= 0 .and. mpi_enreg%paral_kgb /= 0) then
   call init_sc_dmft_paralkgb(paw_dmft,mpi_enreg)
 end if

 if (mbandc /= dmftbandf-dmftbandi+1 .and. paw_dmft%use_dmft /= 10) then
   write(message,'(5a)') ' BUG init_sc_dmft',ch10,&
    & '  number of bands in dmft is not correctly computed ',ch10, &
    & '  Action : check the code'
   ABI_BUG(message)
 end if

 if (dmft_solv == 6 .or. dmft_solv == 7) then

   if (off_diag == 1) then
     write(message,'(3a)') "WARNING: You have activated the sampling of the off-diagonal elements ", &
                         & "in TRIQS/CTHYB. Some features are not available, and you will not be ", &
                         & "able to compute an energy."
     ABI_WARNING(message)
   end if
#ifdef HAVE_TRIQS_COMPLEX
   if (off_diag == 0) then
     write(message,'(2a)') "WARNING: You have compiled with the complex version of TRIQS/CTHYB, yet you do not", &
                     & " sample any off-diagonal element. This is a waste of computation time."
     ABI_WARNING(message)
   end if
#else
   if (off_diag == 1) then
     write(message,'(3a)') "WARNING: You have compiled with the real version of TRIQS/CTHYB, yet you have", &
                & " activated the sampling of the off-diagonal elements. Thus their imaginary part will be", &
                & " neglected. You'll have to check that this is a valid approximation."
     ABI_WARNING(message)
   end if
#endif
 end if ! dmft_solv=6 or 7

 write(message,'(7a)') ch10, &
  & ch10," ******************************************", &
  & ch10," DFT+DMFT Method is used", &
  & ch10," ******************************************"
 call wrtout([std_out,ab_out],message,'COLL')

 if (use_dmft /= 10) then
   if (dmft_solv == 0) then
     write(message,'(2a)') ch10,' DMFT check: no solver and U=J=0'
   else if (dmft_solv == 1) then
     write(message,'(2a)') ch10,' DMFT check: static solver'
   else if (dmft_solv == -1) then
     write(message,'(2a)') ch10,' DMFT check: static solver without renormalization of projectors: should recover DFT+U'
   else if (dmft_solv == 2) then
     write(message,'(2a)') ch10,' DMFT uses the Hubbard one solver'
   else if (dmft_solv == 4) then
     write(message,'(2a)') ch10,' DMFT uses the Hirsch Fye solver'
   else if (dmft_solv == 5) then
     write(message,'(2a)') ch10,' DMFT uses the Continuous Time Quantum Monte Carlo solver of ABINIT'
   else if (dmft_solv == 6) then
     write(message,'(2a)') ch10,' DMFT uses the Continuous Time Quantum Monte Carlo solver of TRIQS &
       &(with density density interactions)'
   else if (dmft_solv == 7) then
     write(message,'(2a)') ch10,' DMFT uses the Continuous Time Quantum Monte Carlo solver of TRIQS &
       &(with rotationally invariant interactions)'
   else if (dmft_solv == 9) then
     write(message,'(2a)') ch10,' DMFT uses the python invocation of TRIQS, for which you need to &
       & give your personal script'
   end if ! dmft_solv
 else if(use_dmft == 10) then
   write(message, '(a,a)') ch10,' DMFT uses the python invocation and orbitals constructed using Wannier90 '
 endif
 call wrtout([std_out,ab_out],message,'COLL')

 ! OG: What is all that? Something as moved? START
 if (use_dmft /= 10) then
 if (dmft_dc == 1) then
   dc_string = "Magnetic FLL (Full Localized Limit)"
 else if (dmft_dc == 2) then
   dc_string = "Magnetic AMF (Around Mean Field)"
 else if (dmft_dc == 5) then
   dc_string = "Non-Magnetic FLL (Full Localized Limit)"
 else if (dmft_dc == 6) then
   dc_string = "Non-Magnetic AMF (Around Mean Field)"
 else if (dmft_dc == 7) then
   dc_string = "Non-Magnetic nominal"
 else if (dmft_dc == 8) then
   dc_string = "Non-Magnetic exact"
 end if
 dc_string = trim(dc_string) // " double counting"

 lda_string = "Magnetic DFT, with "
 if (dtset%usepawu == 14) lda_string = "Non " // trim(adjustl(lda_string))
 write(message,'(2(a,1x),a)') ch10,trim(adjustl(lda_string)),trim(adjustl(dc_string))
 call wrtout([std_out,ab_out],message,'COLL')

 if (dtset%dmft_entropy == 0 .and. ((dmft_solv /= 6 .and. dmft_solv /= 7)  &
    & .or. (dtset%dmft_triqs_entropy == 0 .or. &
    & dtset%dmft_triqs_compute_integral == 0))) then
   write(message,'(a,1x,a)') ch10,"Entropy is not computed, only the internal energy is printed"
   call wrtout([std_out,ab_out],message,'COLL')
 end if

 if ((dmft_solv == 6 .or. dmft_solv == 7) .and. off_diag == 1)  then
#ifndef HAVE_TRIQS_COMPLEX
   write(message,'(a,1x,a)') ch10,"The imaginary part of the Green's function is neglected"
   call wrtout([std_out,ab_out],message,'COLL')
#endif
 else if (dmft_solv /= 6 .and. dmft_solv /= 7) then
   write(message,'(a,1x,a)') ch10,"The imaginary part of the Green's function is neglected"
   call wrtout([std_out,ab_out],message,'COLL')
 end if
 if (dmft_solv == 5 .or. ((dmft_solv == 6 .or. dmft_solv == 7) .and. off_diag == 0)) then
   write(message,'(a,1x,a)') ch10,"The off-diagonal elements of the Green's function are neglected"
   call wrtout([std_out,ab_out],message,'COLL')
 end if

!=============================
!==  Define integers and reals
!=============================

 paw_dmft%nelectval = dble(dtset%nelect)

 if (dmft_solv /= 6 .and. dmft_solv /= 7) then
   fac = merge(2,1,nsppol==1.and.nspinor==1)
   paw_dmft%nelectval = dble(dtset%nelect-(dmftbandi-1)*nsppol*fac)
 end if ! not use_all_bands

 paw_dmft%natpawu              = dtset%natpawu
 paw_dmft%natom                = natom
 paw_dmft%temp                 = dtset%tsmear!*unit_e
 paw_dmft%dmft_iter            = dtset%dmft_iter
 paw_dmft%dmft_entropy         = dtset%dmft_entropy
 paw_dmft%dmft_kspectralfunc   = dtset%dmft_kspectralfunc
 paw_dmft%dmft_magnfield       = dtset%dmft_magnfield
 paw_dmft%dmft_magnfield_b     = dtset%dmft_magnfield_b
 paw_dmft%dmft_dc              = dmft_dc
 paw_dmft%dmft_wanorthnorm     = dtset%dmft_wanorthnorm
 paw_dmft%prtvol               = dtset%prtvol
 paw_dmft%prtdos               = dtset%prtdos
 paw_dmft%dmft_tolfreq         = dtset%dmft_tolfreq
 paw_dmft%dmft_lcpr            = dtset%dmft_tollc
 paw_dmft%dmft_charge_prec     = dtset%dmft_charge_prec
 paw_dmft%dmft_fermi_prec      = dtset%dmft_charge_prec * ten
 paw_dmft%dmft_fermi_step      = dtset%dmft_fermi_step
 paw_dmft%dmft_prt_maxent      = dtset%dmft_prt_maxent
 paw_dmft%dmft_prtself         = dtset%dmft_prtself
 paw_dmft%dmft_prtwan          = dtset%dmft_prtwan
 paw_dmft%dmft_wanrad          = dtset%dmft_wanrad
 paw_dmft%dmft_t2g             = dtset%dmft_t2g
 paw_dmft%dmft_x2my2d          = dtset%dmft_x2my2d

 ! for entropy (alternate external calculation)
 paw_dmft%ientropy = 0
 paw_dmft%u_for_s  = 4.1_dp
 paw_dmft%j_for_s  = 0.5_dp

 paw_dmft%kptopt = dtset%kptopt

!=======================
!==  Choose solver
!=======================

 paw_dmft%dmft_solv = merge(2,dmft_solv,dmft_solv==-2)
 paw_dmft%dmft_blockdiag = merge(1,0,dmft_solv==-2)

!  0: DFT, no solver
!  1: DFT+U
! -1: DFT+U but DFT values are not renormalized !
! if((paw_dmft%dmft_solv==0.and.paw_dmft%prtvol>4).or.&
!&   (paw_dmft%dmft_solv>=-1.and.paw_dmft%dmft_solv<=2)) then
!   call wrtout(std_out,message,'COLL')
!   call wrtout(ab_out,message,'COLL')
! endif

!=======================
!==  Frequencies
!=======================

 paw_dmft%dmft_log_freq = merge(0,1,dmft_solv==6.or.dmft_solv==7.or.dmft_solv==9)

 paw_dmft%dmft_nwli = dtset%dmft_nwli
 paw_dmft%dmft_nwlo = merge(dtset%dmft_nwlo,dtset%dmft_nwli,paw_dmft%dmft_log_freq==1)
 paw_dmft%dmft_nwr = 800

 paw_dmft%dmft_rslf = dtset%dmft_rslf
 paw_dmft%dmft_mxsf = dtset%dmft_mxsf

!=======================
!==  CTQMC
!=======================

 paw_dmft%dmftqmc_l     = dtset%dmftqmc_l
 paw_dmft%dmftqmc_n     = dtset%dmftqmc_n
 paw_dmft%dmftqmc_seed  = dtset%dmftqmc_seed
 paw_dmft%dmftqmc_therm = dtset%dmftqmc_therm

 paw_dmft%dmftctqmc_basis  = dtset%dmftctqmc_basis
 paw_dmft%dmftctqmc_check  = dtset%dmftctqmc_check
 paw_dmft%dmftctqmc_correl = dtset%dmftctqmc_correl
 paw_dmft%dmftctqmc_gmove  = dtset%dmftctqmc_gmove
 paw_dmft%dmftctqmc_grnns  = dtset%dmftctqmc_grnns
 paw_dmft%dmftctqmc_meas   = dtset%dmftctqmc_meas
 paw_dmft%dmftctqmc_mrka   = dtset%dmftctqmc_mrka
 paw_dmft%dmftctqmc_mov    = dtset%dmftctqmc_mov
 paw_dmft%dmftctqmc_order  = dtset%dmftctqmc_order
 paw_dmft%dmftctqmc_localprop = dtset%dmftctqmc_localprop

 if (dmft_solv == 5 .or. dmft_solv >= 8) then
   write(message,'(2a,i6)') ch10,&
     & '=> Seed for CT-QMC inside DMFT is dmftqmc_seed=',paw_dmft%dmftqmc_seed
   call wrtout(std_out,message,'COLL')
 end if

!=======================
!==  TRIQS CTQMC
!=======================

 paw_dmft%dmft_triqs_nleg                          = dtset%dmft_triqs_nleg
 paw_dmft%dmft_triqs_therm_restart                 = dtset%dmft_triqs_therm_restart
 paw_dmft%dmft_triqs_det_init_size                 = dtset%dmft_triqs_det_init_size
 paw_dmft%dmft_triqs_det_n_operations_before_check = dtset%dmft_triqs_det_n_operations_before_check
 paw_dmft%dmft_triqs_move_shift                    = (dtset%dmft_triqs_move_shift == 1)
 paw_dmft%dmft_triqs_move_double                   = (dtset%dmft_triqs_move_double == 1)
 paw_dmft%dmft_triqs_loc_n_min                     = dtset%dmft_triqs_loc_n_min
 paw_dmft%dmft_triqs_loc_n_max                     = dtset%dmft_triqs_loc_n_max
 paw_dmft%dmft_triqs_seed_a                        = dtset%dmft_triqs_seed_a
 paw_dmft%dmft_triqs_seed_b                        = dtset%dmft_triqs_seed_b
 paw_dmft%dmft_triqs_measure_density_matrix        = (dtset%dmft_triqs_measure_density_matrix == 1)
 paw_dmft%dmft_triqs_time_invariance               = (dtset%dmft_triqs_time_invariance == 1)
 paw_dmft%dmft_triqs_use_norm_as_weight            = (dtset%dmft_triqs_use_norm_as_weight == 1)
 paw_dmft%dmft_triqs_leg_measure                   = (dtset%dmft_triqs_leg_measure == 1)
 paw_dmft%dmft_triqs_off_diag                      = (off_diag == 1)
 paw_dmft%dmft_triqs_imag_threshold                = dtset%dmft_triqs_imag_threshold
 paw_dmft%dmft_triqs_det_precision_warning         = dtset%dmft_triqs_det_precision_warning
 paw_dmft%dmft_triqs_det_precision_error           = dtset%dmft_triqs_det_precision_error
 paw_dmft%dmft_triqs_det_singular_threshold        = dtset%dmft_triqs_det_singular_threshold
 paw_dmft%dmft_triqs_epsilon                       = dtset%dmft_triqs_epsilon
 paw_dmft%dmft_triqs_lambda                        = dtset%dmft_triqs_wmax / dtset%tsmear
 paw_dmft%dmft_triqs_entropy                       = dtset%dmft_triqs_entropy
 paw_dmft%dmft_triqs_compute_integral              = dtset%dmft_triqs_compute_integral
 paw_dmft%dmft_triqs_gaussorder                    = dtset%dmft_triqs_gaussorder
 paw_dmft%dmft_triqs_nsubdivisions                 = dtset%dmft_triqs_nsubdivisions
 paw_dmft%dmft_triqs_tol_block                     = dtset%dmft_triqs_tol_block
 paw_dmft%dmft_triqs_read_ctqmcdata                = dtset%dmft_triqs_read_ctqmcdata
 paw_dmft%dmft_triqs_pauli_prob                    = dtset%dmft_triqs_pauli_prob

!==============================
!==  Variables for DMFT itself
!==============================

 paw_dmft%wtk => dtset%wtk(:)
 ! In the case where we sample the full BZ, don't overwrite the wtk with 1/nkpt when we use TRIQS
 if (dtset%iscf < 0 .and. (dtset%kptopt < 0 .or. &
   & (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7))) paw_dmft%wtk(:) = one / dble(nkpt)
 sumwtk = sum(paw_dmft%wtk(1:nkpt))
 if (abs(sumwtk-one) > tol11) then
   write(message,'(a,f15.11)') ' sum of k-point is incorrect',sumwtk
   ABI_BUG(message)
 end if

 t2g = (paw_dmft%dmft_t2g == 1)
 x2my2d = (paw_dmft%dmft_x2my2d == 1)

 paw_dmft%typat => dtset%typat(:)

 ABI_MALLOC(paw_dmft%lpawu,(natom))

 maxlpawu = 0
 do iatom=1,natom
   lpawu = pawtab(paw_dmft%typat(iatom))%lpawu
   if (t2g .and. lpawu /= -1) lpawu = 1
   if (x2my2d .and. lpawu /= -1) lpawu = 0
   if (lpawu > maxlpawu) maxlpawu = lpawu
   paw_dmft%lpawu(iatom) = lpawu
 end do ! iatom
 paw_dmft%maxlpawu = maxlpawu
 ndim = 2*maxlpawu + 1

 paw_dmft%maxnproju = 0
 ABI_MALLOC(paw_dmft%siz_proj,(ntypat))
 ABI_MALLOC(lcycle,(maxlpawu+1))

 lcycle(:) = .true.

 do itypat=1,ntypat
   lpawu = pawtab(itypat)%lpawu
   if (lpawu == -1) cycle
   if (t2g) lpawu = 1
   if (x2my2d) lpawu = 0
   paw_dmft%siz_proj(itypat) = size(pawtab(itypat)%proj(:))
   lcycle(lpawu+1) = .false.
   if (pawtab(itypat)%nproju > paw_dmft%maxnproju) paw_dmft%maxnproju = pawtab(itypat)%nproju
 end do ! itypat

 ABI_MALLOC(paw_dmft%slm2ylm,(ndim,ndim,maxlpawu+1))
 ABI_MALLOC(paw_dmft%jmj2ylm,(2*ndim,2*ndim,maxlpawu+1))
 ABI_MALLOC(paw_dmft%zarot,(ndim,ndim,nsym,maxlpawu+1))

 paw_dmft%slm2ylm(:,:,:) = czero
 do lpawu=0,maxlpawu
   if (lcycle(lpawu+1)) cycle
   ndim = 2*lpawu + 1
   do im=1,ndim
     mm = im - lpawu - 1 ; jm = - mm + lpawu + 1
     onem = (-1)**mm
     if (mm > 0) then
       paw_dmft%slm2ylm(im,im,lpawu+1) = cmplx(onem/sqrt2,zero,kind=dp)
       paw_dmft%slm2ylm(jm,im,lpawu+1) = cmplx(one/sqrt2,zero,kind=dp)
     end if
     if (mm == 0) paw_dmft%slm2ylm(im,im,lpawu+1) = cone
     if (mm < 0) then
       paw_dmft%slm2ylm(im,im,lpawu+1) =  cmplx(zero,one/sqrt2,kind=dp)
       paw_dmft%slm2ylm(jm,im,lpawu+1) = -cmplx(zero,onem/sqrt2,kind=dp)
     end if
   end do ! im
 end do ! lpawu

 paw_dmft%jmj2ylm(:,:,:) = czero
 do ll=1,maxlpawu
   if (lcycle(ll+1)) cycle
   ABI_MALLOC(ind_msml,(2,-ll:ll))
   jc1 = 0
   do ms1=1,2
     do ml1=-ll,ll
       jc1 = jc1 + 1
       ind_msml(ms1,ml1) = jc1
     end do ! ml1
   end do ! ms1
   invsqrt2lp1 = one / sqrt(dble(2*ll+1))
   jc1 = 0
   do jj=ll,ll+1
     xj = dble(jj) - half ! xj is in {ll-0.5,ll+0.5}
     do jm=-jj,jj-1
       xmj = dble(jm) + half ! xmj is in {-xj,xj}
       jc1 = jc1 + 1 ! Global index for JMJ
       if (nint(xj+half) == ll+1) then ! if xj=ll+0.5
         if (nint(xmj+half) == ll+1) then
           paw_dmft%jmj2ylm(ind_msml(1,ll),jc1,ll+1) = cone   !  J=L+0.5 and m_J=L+0.5
         else if (nint(xmj-half) == -ll-1) then
           paw_dmft%jmj2ylm(ind_msml(2,-ll),jc1,ll+1) = cone   !  J=L+0.5 and m_J=-L-0.5
         else
           paw_dmft%jmj2ylm(ind_msml(1,nint(xmj-half)),jc1,ll+1) = cmplx(invsqrt2lp1*(sqrt(dble(ll)+xmj+half)),zero,kind=dp)
           paw_dmft%jmj2ylm(ind_msml(2,nint(xmj+half)),jc1,ll+1) = cmplx(invsqrt2lp1*(sqrt(dble(ll)-xmj+half)),zero,kind=dp)
         end if
       end if
       if (nint(xj+half) == ll) then  ! if xj=ll-0.5
         paw_dmft%jmj2ylm(ind_msml(2,nint(xmj+half)),jc1,ll+1) = cmplx(invsqrt2lp1*(sqrt(dble(ll)+xmj+half)),zero,kind=dp)
         paw_dmft%jmj2ylm(ind_msml(1,nint(xmj-half)),jc1,ll+1) = cmplx(-invsqrt2lp1*(sqrt(dble(ll)-xmj+half)),zero,kind=dp)
       end if
     end do ! jm
   end do ! jj
   ABI_FREE(ind_msml)
 end do ! ll

 do lpawu=0,maxlpawu
   if (lcycle(lpawu+1)) cycle
   ndim = 2*lpawu + 1
   if (t2g) then
     do im1=1,ndim
       do im=1,ndim
         paw_dmft%zarot(im,im1,:,lpawu+1) = cmplx(pawang%zarot(mt2g(im),mt2g(im1),3,1:nsym),zero,kind=dp)
       end do ! im
     end do ! im1
   else if (x2my2d) then
     paw_dmft%zarot(1,1,:,lpawu+1) = cmplx(pawang%zarot(5,5,3,1:nsym),zero,kind=dp)
   else
     paw_dmft%zarot(1:ndim,1:ndim,:,lpawu+1) = cmplx(pawang%zarot(1:ndim,1:ndim,lpawu+1,1:nsym),zero,kind=dp)
   end if
 end do ! lpawu

!=======================
! Imaginary frequencies
!=======================
! Set up log frequencies
 if (dtset%ucrpa == 0 .and. paw_dmft%dmft_nwlo > 0) then
   call construct_nwlo_dmft(paw_dmft)
 end if

 if (paw_dmft%dmftcheck == 1 .and. dmft_solv < 4) paw_dmft%dmftqmc_l = 64

!==============
! Radial grid
!==============

 if (paw_dmft%dmft_prtwan == 1) then

   ! Initialize rmax (maximum radius for print of Wannier functions)
   ABI_MALLOC(rmax,(maxlpawu+1))
   rmax(:) = paw_dmft%dmft_wanrad

   if (paw_dmft%dmft_wanrad < 0) then  ! default

     ! Set rmax to half the distance from the current atom
     ! to the nearest atom with the same lpawu

     ! First take half * min(|Ri|)
     rmax(:) = zero
     verif = .true.
     do i=1,3
       norm = norm2(rprimd(1:3,i)) * half
       if (verif .or. norm < rmax(1)) then
         rmax(:) = norm
         verif = .false.
       end if
     end do ! i

     ! Now look at the atoms in the same unit cell
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       do iatom1=1,natom
         lpawu1 = paw_dmft%lpawu(iatom1)
         if (lpawu1 == -1) cycle
         if (lpawu /= lpawu1 .or. iatom == iatom1) cycle
         norm = zero
         do i=1,3
           norm = norm + dot_product(xred(1:3,iatom)-xred(1:3,iatom1),rprimd(i,1:3))**2
         end do
         norm = sqrt(norm) * half
         if (norm < rmax(lpawu+1)) rmax(lpawu+1) = norm
       end do ! iatom1
     end do ! iatom

   end if ! dmft_wanrad < 0

 end if ! prtwan=1

 ! Now build radial grid by extending the PAW mesh up to max(rmax,size(proj))
 ! The mesh inside the PAW sphere is still exactly the same.
 use_full_chipsi = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)
 paw_dmft%int_meshsz => pawrad(:)%int_meshsz

 if (use_full_chipsi) then
   ABI_MALLOC(paw_dmft%phimtphi,(maxval(pawrad(1:ntypat)%int_meshsz),paw_dmft%maxnproju,ntypat))
   ABI_MALLOC(paw_dmft%phimtphi_int,(paw_dmft%maxnproju,ntypat))
 else
   ABI_MALLOC(paw_dmft%phi_int,(paw_dmft%maxnproju,ntypat))
 end if

 ABI_MALLOC(paw_dmft%radgrid,(ntypat))

 paw_dmft%maxmeshsize = 0
 do itypat=1,ntypat
   lpawu = pawtab(itypat)%lpawu
   if (lpawu == -1) cycle
   if (t2g) lpawu = 1
   if (x2my2d) lpawu = 0
   mesh_type = pawrad(itypat)%mesh_type
   lstep = pawrad(itypat)%lstep
   rstep = pawrad(itypat)%rstep
   siz_proj = paw_dmft%siz_proj(itypat)
   mesh_size = siz_proj
   if (paw_dmft%dmft_prtwan == 1) then
     if (mesh_type == 1) mesh_size = int(rmax(lpawu+1)/rstep) + 1
     if (mesh_type == 2) mesh_size = int(log(rmax(lpawu+1)/rstep+1)/lstep) + 1
     if (mesh_type == 3) mesh_size = int(log(rmax(lpawu+1)/rstep)/lstep) + 2
     if (mesh_size < siz_proj) then
       message = "Please set wanrad to a value greater than the radius of your DMFT orbital"
       ABI_ERROR(message)
     end if
   end if ! prtwan=1
   ! mesh_type > 3 cannot be extended outside the PAW sphere while keeping the
   ! mesh inside the sphere unchanged.
   if ((mesh_size /= pawrad(itypat)%mesh_size) .and. mesh_type > 3) then
     message = "mesh_type > 3 is only compatible with dmft_orbital=1 and dmft_prtwan=0"
     ABI_ERROR(message)
   end if
   if (mesh_size > pawrad(itypat)%int_meshsz .and. (.not. use_full_chipsi)) then
     message = "You need to activate use_full_chipsi if you use an orbital that extends outside the PAW sphere"
     ABI_ERROR(message)
   end if
   call pawrad_init(paw_dmft%radgrid(itypat),mesh_size,mesh_type,rstep,lstep)
   if (mesh_size > paw_dmft%maxmeshsize) paw_dmft%maxmeshsize = mesh_size
   siz_paw  = min(mesh_size,pawrad(itypat)%int_meshsz)
   siz_proj = min(siz_proj,pawrad(itypat)%int_meshsz)
   rint   = paw_dmft%radgrid(itypat)%rad(siz_proj)
   nproju = pawtab(itypat)%nproju
   do iproj=1,nproju
     indproj = pawtab(itypat)%lnproju(iproj)
     if (use_full_chipsi) then
       ! Precompute <Chi|Phi-Phi_tilde>
       paw_dmft%phimtphi(1:siz_paw,iproj,itypat) = pawtab(itypat)%phi(1:siz_paw,indproj) - &
                                                 & pawtab(itypat)%tphi(1:siz_paw,indproj)
       call simp_gen(paw_dmft%phimtphi_int(iproj,itypat),pawtab(itypat)%proj(1:siz_proj)* &
                   & paw_dmft%phimtphi(1:siz_proj,iproj,itypat),paw_dmft%radgrid(itypat),r_for_intg=rint)
     else
       ! Precompute <Chi|Phi>
       call simp_gen(paw_dmft%phi_int(iproj,itypat),pawtab(itypat)%proj(1:siz_proj)* &
                   & pawtab(itypat)%phi(1:siz_proj,indproj),paw_dmft%radgrid(itypat),r_for_intg=rint)
     end if ! use_full_chipsi
   end do ! iproj
 end do ! itypat

 if (paw_dmft%dmft_prtwan /= 1 .and. use_full_chipsi) then
   ABI_FREE(paw_dmft%phimtphi)
 end if
 ABI_SFREE(rmax)

!==============
! Plane waves
!==============

 if (use_full_chipsi) then

   ! Compute ylm(k+G),exp(j*(k+G).R(iat)) and Bessel functions
   mpw = dtset%mpw
   ABI_MALLOC(paw_dmft%ylm,(mpw,2*maxlpawu+1,maxlpawu+1,mkmem))
   ABI_MALLOC(paw_dmft%dpro,(mpw,natom,mkmem))
   ABI_MALLOC(paw_dmft%bessel,(mpw,paw_dmft%maxmeshsize,ntypat,mkmem))
   ABI_MALLOC(paw_dmft%bessel_int,(mpw,ntypat,mkmem))
   ABI_MALLOC(typcycle,(ntypat))
   ABI_MALLOC(kpg,(3,mpw))
   ABI_MALLOC(kpg_norm,(mpw))

   ik  = 0 ! kpt index on current CPU
   ikg = 0

   do ikpt=1,nkpt

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,1,-1,mpi_enreg%me_kpt)) cycle

     ik  = ik + 1
     npw = npwarr(ikpt)

     do ig=1,npw

       kpg(:,ig) = dtset%kptns(1:3,ikpt) + dble(kg(1:3,ikg+ig))

       norm = zero
       do i=1,3
         norm = norm + dot_product(kpg(:,ig),gprimd(i,1:3))**2
       end do ! i

       kpg_norm(ig) = sqrt(norm)

     end do ! ig

     lcycle(:)   = .false.
     typcycle(:) = .false.

     do iatom=1,natom

       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       itypat = paw_dmft%typat(iatom)
       ndim   = 2*lpawu + 1

       if (.not. lcycle(lpawu+1)) then  ! if this l has not been visited
         if (t2g) then
           do im=1,ndim
             paw_dmft%ylm(1:npw,im,lpawu+1,ik) = ylm(ikg+1:ikg+npw,4+mt2g(im))
           end do ! im
         else if (x2my2d) then
           paw_dmft%ylm(1:npw,1,lpawu+1,ik) = ylm(ikg+1:ikg+npw,9)
         else
           paw_dmft%ylm(1:npw,1:ndim,lpawu+1,ik) = ylm(ikg+1:ikg+npw,lpawu**2+1:lpawu**2+ndim)
         end if
         lcycle(lpawu+1) = .true.
       end if ! not lcycle

       do ig=1,npw
         paw_dmft%dpro(ig,iatom,ik) = exp(j_dpc*two_pi*dot_product(kpg(:,ig),xred(1:3,iatom)))
       end do ! ig

       if (.not. typcycle(itypat)) then   ! if this type has not been visited
         lpawu1 = lpawu ! physical l
         if (t2g .or. x2my2d) lpawu1 = 2
         siz_proj = paw_dmft%siz_proj(itypat)
         rint = paw_dmft%radgrid(itypat)%rad(siz_proj)
         siz_wan = paw_dmft%radgrid(itypat)%mesh_size
         do ir=1,siz_wan
           rad = paw_dmft%radgrid(itypat)%rad(ir)
           do ig=1,npw
             call paw_jbessel_4spline(bes,besp,lpawu1,0,two_pi*kpg_norm(ig)*rad,tol3)
             ! Multiply by r since we want to compute Psi(r) * r, for radial integration
             paw_dmft%bessel(ig,ir,itypat,ik) = four_pi * bes * rad / sqrt(ucvol)
           end do ! ig
         end do ! ir
         do ig=1,npw
           call simp_gen(bes,pawtab(itypat)%proj(1:siz_proj)*dble(paw_dmft%bessel(ig,1:siz_proj,itypat,ik)), &
                       & paw_dmft%radgrid(itypat),r_for_intg=rint)
           paw_dmft%bessel_int(ig,itypat,ik) = bes * (j_dpc**lpawu1) ! CAREFUL: we multiply by j^l AFTER simp_gen since simp_gen doesn_t handle complex
         end do ! ig
         paw_dmft%bessel(1:npw,1:siz_wan,itypat,ik) = paw_dmft%bessel(1:npw,1:siz_wan,itypat,ik) * (j_dpc**lpawu1)
         typcycle(itypat) = .true.
       end if ! not typcycle

     end do ! iatom

     ikg = ikg + npw

   end do ! ikpt

   ABI_FREE(kpg)
   ABI_FREE(kpg_norm)
   ABI_FREE(typcycle)

   if (paw_dmft%dmft_prtwan /= 1) then
     ABI_FREE(paw_dmft%bessel)
   end if

 end if ! use_full_chipsi

 ABI_FREE(lcycle)

 call init_paral_dmft(paw_dmft,paw_dmft%distrib,paw_dmft%dmft_nwlo)

 ! OG: What is all that? Something as moved? START
 endif

end subroutine init_sc_dmft
!!***

!!****f* m_paw_dmft/init_dmft
!! NAME
!! init_dmft
!!
!! FUNCTION
!!  Allocate variables and setup DFT hamiltonian and related data
!!  (init_sc_dmft has to been called before)
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)> = crystal structure data
!!  dmatpawu = fixed occupation matrix of correlated orbitals
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  fermie_dft = DFT Fermi level
!!  filctqmcdatain = input file name for CTQMC_DATA file
!!  filselfin = input file name for self file
!!  fnamei = input file name
!!  fnametmp_app = header for the output filename
!!  ireadctqmcdata = flag to read CTQMC_DATA input file at first iteration
!!  ireadself = flag to read self input file at first iteration
!!  paw_dmft <type(paw_dmft_type)>= paw+dmft related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SOURCE
!!
!! NOTE
!! The part of the code which deals
!! with the use of logarithmic frequencies
!! is a modification of the GNU GPL
!! code available on http://dmft.rutgers.edu/ and
!! described in the  RMP paper written by
!! G.Kotliar,  S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti.

subroutine init_dmft(cryst_struc,dmatpawu,dtset,fermie_dft,filctqmcdatain,filselfin,fnamei,fnametmp_app,ireadctqmcdata,ireadself,paw_dmft,pawtab)

!Arguments ------------------------------------
 real(dp), intent(in) :: fermie_dft
 type(dataset_type), intent(in) :: dtset
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(crystal_t), target, intent(in) :: cryst_struc
 character(len=fnlen), intent(in) :: filctqmcdatain,filselfin,fnamei,fnametmp_app
 integer, intent(in) :: ireadctqmcdata,ireadself
 real(dp), target, intent(in) :: dmatpawu(:,:,:,:)
 type(pawtab_type), intent(in) :: pawtab(dtset%ntypat)
!Local variables ------------------------------------
 integer :: grid_unt,iatom,ierr,ifreq,ioerr,ir,irot,isym
 integer :: itypat,lpawu,meshsz,nflavor,ngrid,nsym,unt
 real(dp) :: int1,step
 logical :: lexist
 character(len=4) :: tag_at
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
! *********************************************************************

 if (dtset%ucrpa == 0) then
   write(message,'(7a)') ch10,' ====================================', &
                       & ch10,' =====  Start of DMFT calculation', &
                       & ch10,' ====================================',ch10
 else if (dtset%ucrpa > 0) then
   write(message,'(6a)') ch10,' ============================================================', &
                       & ch10,' =====  Initialize construction of Wannier in DMFT routines',&
                       & ch10,' ============================================================'
 end if ! ucrpa
 call wrtout(std_out,message,'COLL')

 if (paw_dmft%dmft_dc == 8) then
   ABI_MALLOC(paw_dmft%edc,(paw_dmft%natom))
   ABI_MALLOC(paw_dmft%edcdc,(paw_dmft%natom))
 end if

 nsym = paw_dmft%nsym

!=======================
!==  Check sym
!=======================

 do isym=1,nsym
   if (dtset%symafm(isym) < 0) then
     message = 'symafm negative is not implemented in DMFT'
     ABI_ERROR(message)
   end if
 end do ! isym

 paw_dmft%nsym = cryst_struc%nsym ! very important to update it here
 nsym = paw_dmft%nsym

 ! TODO: this really should be done in init_sc_dmft
 paw_dmft%indsym => cryst_struc%indsym(4,1:nsym,1:paw_dmft%natom)
 if (paw_dmft%nspinor == 2) then
   ABI_MALLOC(paw_dmft%symrec_cart,(3,3,nsym))
   do irot=1,nsym
     call symredcart(cryst_struc%gprimd(:,:),cryst_struc%rprimd(:,:),&
                   & paw_dmft%symrec_cart(:,:,irot),cryst_struc%symrec(:,:,irot))
   end do ! irot
 end if ! nspinor=2

 paw_dmft%filapp         = fnametmp_app
 paw_dmft%filnamei       = fnamei
 paw_dmft%filselfin      = filselfin
 paw_dmft%filctqmcdatain = filctqmcdatain
 paw_dmft%ireadctqmcdata = ireadctqmcdata
 paw_dmft%ireadself      = ireadself

 ! Write orbital on file
 if (paw_dmft%myproc == 0) then
   do itypat=1,paw_dmft%ntypat
     lpawu = pawtab(itypat)%lpawu
     if (lpawu == -1) cycle
     meshsz = paw_dmft%siz_proj(itypat)

     call simp_gen(int1,pawtab(itypat)%proj(1:meshsz)**2,paw_dmft%radgrid(itypat), &
                 & r_for_intg=paw_dmft%radgrid(itypat)%rad(meshsz))
     int1 = sqrt(int1)

     call int2char4(itypat,tag_at)
     ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
     if (open_file(trim(paw_dmft%filapp)//"_DMFTORBITAL_itypat"//tag_at//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)

     write(unt,'(4a)') "# Correlated normalized radial orbital for DMFT. This", &
        & " is not projected on any energy window (you need to use dmft_prtwan for that).",ch10, &
        & "#       Radius (Bohr)           u_l(r) = R_l * r"

     do ir=1,meshsz
       write(unt,*) paw_dmft%radgrid(itypat)%rad(ir),pawtab(itypat)%proj(ir)/int1
     end do ! ir

     close(unt)
   end do ! itypat
 end if ! myproc=0

!==================
! Real frequencies
!==================

 if (dtset%iscf < 0 .and. paw_dmft%dmft_solv >= 5 .and. paw_dmft%dmft_solv <= 8) then
   tmpfil = trim(paw_dmft%filapp)//'_spectralfunction_realgrid'
   inquire(file=trim(tmpfil),exist=lexist)!,recl=nrecl)
   grid_unt = 2000
   if (.not. lexist) then
     write(message,'(4x,a,i5,3a)') "File number",grid_unt," called ",trim(tmpfil)," does not exist"
     call wrtout(std_out,message,'COLL')
     message = "Cannot continue: the missing file coming from Maxent code is needed"
     ABI_ERROR(message)
   end if ! not lexist

   if (paw_dmft%myproc == 0) then
#ifdef FC_NAG
     open(unit=grid_unt,file=trim(tmpfil),status='unknown',form='formatted',recl=ABI_RECL)
#else
     open(unit=grid_unt,file=trim(tmpfil),status='unknown',form='formatted')
#endif
     rewind(grid_unt)
     write(message,'(3a)') ch10,"  == Read real frequency grid from file ",trim(tmpfil)
     call wrtout(std_out,message,'COLL')
     write(message,'(5x,3a,i4)') 'Opened file : ',trim(tmpfil),' on unit ',grid_unt
     call wrtout(std_out,message,'COLL')
     read(grid_unt,*,iostat=ioerr) ngrid
     ABI_MALLOC(paw_dmft%omega_r,(ngrid))
     do ifreq=1,ngrid
       read(grid_unt,*,iostat=ioerr) paw_dmft%omega_r(ifreq)
     end do ! ifreq
     close(grid_unt)
   end if ! myproc=0
   call xmpi_bcast(ioerr,0,xmpi_world,ierr)
   if (ioerr /= 0) ABI_ERROR("Error when reading grid file")
   call xmpi_bcast(ngrid,0,xmpi_world,ierr)
   if (paw_dmft%myproc /= 0) then
     ABI_MALLOC(paw_dmft%omega_r,(ngrid))
   end if
   call xmpi_bcast(paw_dmft%omega_r(:),0,xmpi_world,ierr)
 else
   ABI_MALLOC(paw_dmft%omega_r,(2*paw_dmft%dmft_nwr))
   ! Set up real frequencies for spectral function in Hubbard one.
   step = 0.00005_dp
   paw_dmft%omega_r(2*paw_dmft%dmft_nwr) = pi * step * (two*dble(paw_dmft%dmft_nwr-1)+one)
   do ifreq=1,2*paw_dmft%dmft_nwr-1
     paw_dmft%omega_r(ifreq) = pi*step*(two*dble(ifreq-1)+one) - paw_dmft%omega_r(2*paw_dmft%dmft_nwr)
  !  write(std_out,*) ifreq,paw_dmft%omega_r(ifreq)
   end do ! ifreq

 end if ! iscf<0 and dmft_solv>=5 and dmft_solv<=8

 call init_paral_dmft(paw_dmft,paw_dmft%distrib_r,size(paw_dmft%omega_r(:)))

 !unit_e=2_dp

! paw_dmft%dmft_mag=0
! do iatom=1,dtset%natom
!   do  ii=1,3
!     if ( dtset(ii,iatom) > 0.001 ) paw_dmft%dmft_mag=1
!   enddo
! enddo

 paw_dmft%gpu_option = dtset%gpu_option
 paw_dmft%fermie_dft = fermie_dft ! in Ha
 paw_dmft%fermie = fermie_dft

!========================
!==  Fixed self as input
!========================
 paw_dmft%use_fixed_self = dtset%usedmatpu
 paw_dmft%fixed_self => dmatpawu(:,:,:,:)

 if (paw_dmft%dmftcheck == -1) then
   message = ' init_dmft: dmftcheck=-1 should not happen here'
   ABI_BUG(message)
 end if

 ABI_MALLOC(paw_dmft%eigen_dft,(paw_dmft%mbandc,paw_dmft%nkpt,paw_dmft%nsppol))
 ABI_MALLOC(paw_dmft%chipsi,(paw_dmft%nspinor*(2*paw_dmft%maxlpawu+1),paw_dmft%mbandc,paw_dmft%nkpt,paw_dmft%nsppol,paw_dmft%natom))

 paw_dmft%lchipsiortho = 0

 !=========================================================
 !== if we use ctqmc impurity solver
 !=========================================================
 ! IMPORTANT : paw_dmft%hybrid is corrupted somewhere in DMFT routines on
 ! tikal_psc and max2_open64. Use a local hybrid in qmc_prep even if not optimal.
 ! Anyway initializing ctqmc here is not good and produce the same result for
 ! dmft_iter=1 which speed up the convergence ...
 ! FIXME : Move this to init_sc_dmft and find bug
 if (paw_dmft%dmft_solv == 5) then ! CTQMC initialisation
 !  write(message,'(a,2x,a,f13.5)') ch10," == Initializing CTQMC"
 !   call wrtout(std_out,message,'COLL')

   ABI_MALLOC(paw_dmft%hybrid,(paw_dmft%natom))
   do iatom=1,paw_dmft%natom
     if (paw_dmft%lpawu(iatom) == -1) cycle
     nflavor = 2 * (2*paw_dmft%lpawu(iatom)+1)
#ifdef HAVE_MPI
     call CtqmcInterface_init(paw_dmft%hybrid(iatom),paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
        & paw_dmft%dmftqmc_therm,paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
        & std_out,paw_dmft%spacecomm,nspinor=paw_dmft%nspinor)
#else
     call CtqmcInterface_init(paw_dmft%hybrid(iatom),paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
        & paw_dmft%dmftqmc_therm,paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
        & std_out,nspinor=paw_dmft%nspinor)
#endif
     call CtqmcInterface_setOpts(paw_dmft%hybrid(iatom),&
                                &  opt_Fk       = 1,&
                                &  opt_order    = paw_dmft%dmftctqmc_order, &
                                &  opt_histo    = paw_dmft%dmftctqmc_localprop,&
                                &  opt_movie    = paw_dmft%dmftctqmc_mov,   &
                                &  opt_analysis = paw_dmft%dmftctqmc_correl,&
                                &  opt_check    = paw_dmft%dmftctqmc_check, &
                                &  opt_noise    = paw_dmft%dmftctqmc_grnns, &
                                &  opt_spectra  = paw_dmft%dmftctqmc_mrka,  &
                                &  opt_gmove    = paw_dmft%dmftctqmc_gmove )
   end do ! iatom
   ! write(message,'(a,2x,a,f13.5)') ch10,&
   !&  " == Initialization CTQMC done"
   !call wrtout(std_out,message,'COLL')
 end if ! dmft_solv=5

!************************************************************************
end subroutine init_dmft
!!***

!!****f* m_paw_dmft/construct_nwli_dmft
!! NAME
!! construct_nwli_dmft
!!
!! FUNCTION
!!  Compute linear frequencies
!!
!! INPUTS
!!  paw_dmft=structure for dmft
!!  nwli=number of linear frequencies
!!
!! OUTPUTS
!!  omegali(1:nwli)=computed frequencies
!!
!! SOURCE
!!

subroutine construct_nwli_dmft(paw_dmft,nwli,omega_li)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: nwli
 real(dp), intent(inout) :: omega_li(:)
!Local variables-------------------------------
 !fortran2003 ?
 !real(dp), allocatable, intent(inout) :: omega_li(:)
 integer :: ifreq
 real(dp) :: factor
 character(len=100) :: message
! *********************************************************************

! if (allocated(omega_li)) then
 if (size(omega_li(:)) /= nwli) then
   write(message,'(2a,i8,a,i8)') ch10,"Number of linear frequencies asked is", &
       & nwli,"whereas dimension of array omega_li is",size(omega_li(:))
   ABI_BUG(message)
!     ABI_FREE(omega_li)
!     ABI_MALLOC(omega_li,(nwli))
!     write(*,*) "RESIZE"
!     call flush(6)
 end if
!     write(*,*) "NOTHING"
!     call flush(6)
! else
!     write(*,*) "ALLOCATE"
!     call flush(6)
!   ABI_MALLOC(omega_li,(nwli))
! endif

! Set up linear frequencies
 factor = pi * paw_dmft%temp
 do ifreq=1,nwli
   omega_li(ifreq) = factor * dble(2*ifreq-1)
   ! (2(ifreq-1)+1 = 2ifreq-1
 end do ! ifreq

end subroutine construct_nwli_dmft
!!***

!!****f* m_paw_dmft/construct_nwlo_dmft
!! NAME
!! construct_nwlo_dmft
!!
!! FUNCTION
!!  Allocate log frequencies if used and compute them as well as their weight
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!
!!
!! SOURCE
!!
!! NOTE
!! The part of the code which deals
!! with the use of logarithmic frequencies
!! is a modification of the GNU GPL
!! code available on http://dmft.rutgers.edu/ and
!! described in the  RMP paper written by
!! G.Kotliar,  S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti.
!!

subroutine construct_nwlo_dmft(paw_dmft)

 use m_splines

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables-------------------------------
 integer :: deltaw,ifreq,ifreq2,myproc,nlin,nproc,nwli
 integer :: nwlo,omegaBegin,omegaEnd,residu,spacecomm
 character(len=10) :: tag
 character(len=500) :: message
 real(dp) :: deltaomega,expfac,omegamaxmin,prefacexp,temp,wl
 complex(dpc):: ybcbeg,ybcend
 integer, allocatable :: select_log(:)
 real(dp), allocatable :: omega_li(:),omega_lo_tmp(:),wgt_wlo(:)
 complex(dpc), allocatable :: splined_li(:),tospline_lo(:),ysplin2_lo(:)
! *********************************************************************

 nwlo = paw_dmft%dmft_nwlo
 nwli = paw_dmft%dmft_nwli
 temp = paw_dmft%temp

!==  Variables for DMFT related to frequencies
! the part of the code which deals
! with the use of logarithmic frequencies
! is a modification of the GNU GPL
! code available on http://dmft.rutgers.edu/ and
! described in the  RMP paper written by
! G.Kotliar, S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti

!========================================
!== construct log. freq.
 if (paw_dmft%dmft_log_freq == 1) then
!=======================================

   ABI_MALLOC(omega_lo_tmp,(nwlo))
   ABI_MALLOC(wgt_wlo,(nwlo))
   !cubic_freq=0
   !omegamaxmin=paw_dmft%omega_li(paw_dmft%dmft_nwli)-paw_dmft%omega_li(paw_dmft%dmftqmc_l+1)
   nlin = paw_dmft%dmftqmc_l ! number of linear frequencies
   if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) nlin = int(nwlo*half)
   omegamaxmin = pi * temp * two * dble(nwli-nlin-1)

   !if(cubic_freq==1) then

    ! if (paw_dmft%dmft_solv .eq. 5 ) then
    !   write(message, '(2a)') ch10, "Warning : Cubish Mesh not tested with CT-QMC"
    !   ABI_WARNING(message)
    ! end if
!  ------------  CUBIC MESH MESH
!    useless
     !nlin=dble(paw_dmft%dmft_nwli)
     !nlog=dble(paw_dmft%dmft_nwlo)
     !AA=((nlin-one)/nlin/(nlog**2-one)-one/(three*nlin))/((nlog**3-one)/(nlog**2-one)-seven/three)
     !BB=(one/nlin - seven*AA)/three
     !CC=-AA-BB
!    AA=((nlin-one)/nlin/(nlog-one)-one/(nlin))/((nlog**2-one)/(nlog-one)-three)
!    BB=(one/nlin - three*AA)
!    CC=-AA-BB
     !write(message, '(a,16x,2(2x,a))') ch10,"  Cubic Mesh Parameters are"
     !call wrtout(std_out,message,'COLL')
     !write(message, '(3x,a,3(2x,e13.5))') "AA,BB,CC",AA,BB,CC
     !call wrtout(std_out,message,'COLL')
     !do ifreq=1,paw_dmft%dmft_nwlo
      ! t1=dble(ifreq)
       !omega_lo_tmp(ifreq)=(AA*t1**3+BB*t1**2+CC)*omegamaxmin+paw_dmft%omega_li(1)
      ! omega_lo_tmp(ifreq)=(AA*t1**3+BB*t1**2+CC)*omegamaxmin+paw_dmft%temp*pi
!       paw_dmft%omega_lo(ifreq)=(AA*t1**2+BB*t1+CC)*omegamaxmin+paw_dmft%omega_li(1)
!     write(69,*) paw_dmft%omega_lo(ifreq),0.5
     !enddo
   !else
   if (paw_dmft%dmft_solv < 4) then
     paw_dmft%dmftqmc_l = 0
     nlin = 0
   end if

!  ------------  LOGARITHMIC MESH
   deltaomega = half
   expfac     = log(omegamaxmin/deltaomega) / (dble(nwlo-nlin-1)*half)
   prefacexp  = omegamaxmin / (exp(expfac*dble(nwlo-nlin-1))-one)
   ABI_MALLOC(select_log,(nwlo))

!  ------------ IMPOSE LINEAR MESH for w < 2*w_n=(2*l-1)pi/beta
!         Check variables (Already done in chkinp if dmft_solv==5)
   if (nlin > nwlo) then
     write(message,'(2a,i6)') ch10, &
       & ' ERROR: dmft_nwlo has to be at least equal to 2xdmftqmc_l :',2*paw_dmft%dmftqmc_l
     ABI_ERROR(message)
   end if
!         End Check

   call construct_nwli_dmft(paw_dmft,nlin,omega_lo_tmp(1:nlin))
   select_log(1:nlin) = (/ (ifreq,ifreq=1,nlin) /)

     !do ifreq=1,paw_dmft%dmftqmc_l
     !  omega_lo_tmp(ifreq)=(two*DBLE(ifreq-1)+one)*pi*paw_dmft%temp
     !  select_log(ifreq)=ifreq
     !enddo

!  ------------ COMPLETE FREQUENCIES WITH LOG MESH
   wl = temp * pi * dble(2*nlin+1)
   do ifreq=1,nwlo-nlin
       !omega_lo_tmp(ifreq+paw_dmft%dmftqmc_l)=prefacexp*(exp(expfac*float(ifreq-1))-one)+paw_dmft%omega_li(paw_dmft%dmftqmc_l+1)
     omega_lo_tmp(ifreq+nlin) = prefacexp*(exp(expfac*dble(ifreq-1))-one) + wl
!    -------- Impose that each frequency of the logarithmic mesh is on a Matsubara frequency
! FIXME : This may be done for all solver, not only for QMCs
     if (paw_dmft%dmft_solv >= 4) then
       ! Compute the index "n" of iwn
       ifreq2 = nint((omega_lo_tmp(ifreq+nlin)/(temp*pi)-one)*half)
       ! Compute freq
       omega_lo_tmp(ifreq+nlin) = (dble(ifreq2)*two+one) * pi * temp

       if ((ifreq2+1) > nwli) then
         write(message,'(2a,i8)') ch10,&
          & ' BUG: init_dmft, dimension of array select_log is about to be overflown',(ifreq2+1)
         ABI_BUG(message)
       end if
       select_log(nlin+ifreq) = ifreq2 + 1
     end if ! dmft_solv>=4
   end do ! ifreq

!  -------- Suppress duplicate frequencies
! FIXME : So this also should be done for all solver and remove useless
! frequencies
  if (paw_dmft%dmft_solv >= 4) then
    ifreq2 = 1
    do ifreq=2,nwlo-1
      if (select_log(ifreq2) == select_log(ifreq)) cycle
      ifreq2 = ifreq2 + 1
      omega_lo_tmp(ifreq2) = omega_lo_tmp(ifreq)
      select_log(ifreq2) = select_log(ifreq)
    end do ! ifreq
    paw_dmft%dmft_nwlo = ifreq2 + 1
    nwlo = paw_dmft%dmft_nwlo
  end if ! dmft_solv>=4

  omega_lo_tmp(1) = temp * pi
  omega_lo_tmp(nwlo) = temp * pi * dble(2*nwli-1)

  !==================================
  !== Construct weight for log. freq.
  !==================================

  ABI_MALLOC(tospline_lo,(nwlo))
  ABI_MALLOC(splined_li,(nwli))
  ABI_MALLOC(ysplin2_lo,(nwlo))
  ABI_MALLOC(omega_li,(nwli))
  call construct_nwli_dmft(paw_dmft,nwli,omega_li(:))

  !Parallelization over frequencies!
  ! ============= Set up =============
  myproc = paw_dmft%myproc
  nproc  = paw_dmft%nproc
  spacecomm = paw_dmft%spacecomm
  deltaw = nwlo / nproc
  residu = nwlo - nproc*deltaw
  if (myproc < nproc-residu) then
    omegaBegin = 1 + myproc*deltaw
    omegaEnd   = (myproc+1) * deltaw
  else
    omegaBegin = 1 + myproc*(deltaw+1) - nproc + residu
    omegaEnd   = omegaBegin + deltaw
  end if

  wgt_wlo(1:nwlo) = zero ! very important for xmpi_sum
  ybcbeg = czero
  ybcend = czero
  ! ============= END Set up =============

  tospline_lo(:) = czero

  do ifreq=omegaBegin,omegaEnd
  ! do ifreq1=1,paw_dmft%dmft_nwlo
    tospline_lo(ifreq) = cone
    ! tospline_lo(ifreq1)=ifreq1**2-ifreq1
    ! enddo
    ! ybcbeg=cmplx(one/tol16**2,zero)
    ! ybcend=cmplx(one/tol16**2,zero)

    !==  spline delta function
    call spline_complex(omega_lo_tmp(:),tospline_lo(:),nwlo, &
                      & ybcbeg,ybcend,ysplin2_lo(:))
    ! do ifreq1=1,paw_dmft%dmft_nwlo
    !  write(6588,*) paw_dmft%omega_lo(ifreq1),ysplin2_lo(ifreq1)
    ! enddo

    call splint_complex(nwlo,omega_lo_tmp(:),tospline_lo(:), &
                      & ysplin2_lo(:),nwli,omega_li(:),splined_li(:))

    tospline_lo(ifreq) = czero

    !==         accumulate weights
    wgt_wlo(ifreq) = sum(dble(splined_li(:)))
    ! do ifreq1=1,paw_dmft%dmft_nwlo
    !  write(6688,*) paw_dmft%omega_lo(ifreq1),tospline_lo(ifreq1)
    ! enddo

    ! do ifreq1=1,paw_dmft%dmft_nwli
    !  write(6788,*) paw_dmft%omega_li(ifreq1),splined_li(ifreq1)

  end do ! ifreq
  ! ============= Gatherall  =============
  call xmpi_sum(wgt_wlo(1:nwlo),spacecomm,residu)
  ! ============= END Gatherall ==========
  ! end parallelisation over frequencies

  ABI_FREE(tospline_lo)
  ABI_FREE(splined_li)
  ABI_FREE(ysplin2_lo)
  ! if(abs(dtset%pawprtvol)>=3) then
  write(message,'(a,18x,2(2x,a21))') ch10,"Log. Freq","weight"
  call wrtout(std_out,message,'COLL')
  do ifreq=1,nwlo
    write(message,'(3x,a9,i6,2(2x,e21.14))') "--ifreq--",ifreq,omega_lo_tmp(ifreq),wgt_wlo(ifreq)
    call wrtout(std_out,message,'COLL')
  end do ! ifreq
  write(message,'(3x,a,i6)') "  Total number of log frequencies is",nwlo
  call wrtout(std_out,message,'COLL')
  ifreq2 = 1
  do ifreq=1,min(30,nwlo)
    write(message,'(3x,a9,i6,2(2x,e21.14))') "--ifreq--",ifreq,omega_li(ifreq)
    call wrtout(std_out,message,'COLL')
    if (select_log(ifreq2) == ifreq) then
      write(message,'(3x,a,i4,2(2x,i5))') "--sel_log",1
      ifreq2 = ifreq + 1
    else
      write(message,'(3x,a,i4,2(2x,i5))') "--sel_log",0
    end if
    call wrtout(std_out,message,'COLL')
  end do ! ifreq
  write(message,'(3x,2a)') "--ifreq--","..."
  call wrtout(std_out,message,'COLL')
  write(message,'(3x,a,i6,2(2x,e13.5))') "--ifreq--",nwli,omega_li(nwli)
  call wrtout(std_out,message,'COLL')
  ! endif
  ABI_FREE(select_log)
  ABI_FREE(omega_li)
  ABI_MALLOC(paw_dmft%omega_lo,(nwlo))
  ABI_MALLOC(paw_dmft%wgt_wlo,(nwlo))
  paw_dmft%omega_lo(1:nwlo) = omega_lo_tmp(1:nwlo)
  paw_dmft%wgt_wlo(1:nwlo) = wgt_wlo(1:nwlo)
  ABI_FREE(omega_lo_tmp)
  ABI_FREE(wgt_wlo)

!=========================================================
!== do not construct log. freq. and use linear frequencies
!=========================================================
 else

   ABI_MALLOC(paw_dmft%omega_lo,(nwlo))
   ABI_MALLOC(paw_dmft%wgt_wlo,(nwlo))
   write(tag,'(i10)') nwlo
   write(message,'(4a)') ch10," Use of ",trim(adjustl(tag))," linear Matsubara frequencies for DMFT calculation"
   call wrtout(std_out,message,'COLL')
   call construct_nwli_dmft(paw_dmft,nwli,paw_dmft%omega_lo(:))
   paw_dmft%wgt_wlo(:) = one
 end if ! dmft_log_freq

 ! Should be check but since type definition does not initialize pointer with
 ! =>null() (fortran95 and later) it produces conditional jump in valgrind
 !if ( associated(paw_dmft%omega_lo) ) then
 !  ABI_FREE(paw_dmft%omega_lo)
 !endif
 !if ( associated(paw_dmft%wgt_wlo) ) then
 !  ABI_FREE(paw_dmft%wgt_wlo)
 !endif

end subroutine construct_nwlo_dmft
!!***

!!****f* m_paw_dmft/destroy_dmft
!! NAME
!! destroy_dmft
!!
!! FUNCTION
!!  deallocate some variables related to paw_dmft
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_dmft(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables-------------------------------
 integer :: iatom
! *********************************************************************

 if (paw_dmft%dmft_solv == 5 .and. allocated(paw_dmft%hybrid)) then
   do iatom=1,size(paw_dmft%hybrid) !paw_dmft%natom
       !if(paw_dmft%lpawu(iatom)/=-1) then
     call ctqmcinterface_finalize(paw_dmft%hybrid(iatom))
       !endif
   end do ! iatom
   ABI_FREE(paw_dmft%hybrid)
 end if
 ABI_SFREE(paw_dmft%chipsi)
 ABI_SFREE(paw_dmft%edc)
 ABI_SFREE(paw_dmft%edcdc)
 ABI_SFREE(paw_dmft%eigen_dft)
 ABI_SFREE(paw_dmft%omega_r)
 ABI_SFREE(paw_dmft%symrec_cart)
 paw_dmft%eigen => null()
 paw_dmft%fixed_self => null()
 paw_dmft%indsym => null()
 call destroy_paral_dmft(paw_dmft%distrib_r)

end subroutine destroy_dmft
!!***

!!****f* m_paw_dmft/destroy_sc_dmft
!! NAME
!! destroy_sc_dmft
!!
!! FUNCTION
!!  deallocate paw_dmft
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_sc_dmft(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
!Local variables-------------------------------
 character(len=500) :: message
! *********************************************************************

 if ((.not. allocated(paw_dmft%occnd) .or. .not. allocated(paw_dmft%band_in) &
   & .or. .not. allocated(paw_dmft%include_bands) .or. .not. allocated(paw_dmft%exclude_bands)) &
   & .and. paw_dmft%use_dmft == 1) then
   write(message,'(3a)') &
     & '  an array is not allocated and is not deallocated with use_dmft==1 ',ch10, &
     & '  Action : check the code'
   ABI_WARNING(message)
 end if
 ABI_SFREE(paw_dmft%occnd)
 ABI_SFREE(paw_dmft%band_in)
 ABI_SFREE(paw_dmft%include_bands)
 ABI_SFREE(paw_dmft%exclude_bands)
 ABI_SFREE(paw_dmft%siz_proj)
 ABI_SFREE(paw_dmft%zarot)
 ABI_SFREE(paw_dmft%phimtphi)
 ABI_SFREE(paw_dmft%phimtphi_int)
 ABI_SFREE(paw_dmft%phi_int)
 ABI_SFREE(paw_dmft%ylm)
 ABI_SFREE(paw_dmft%dpro)
 ABI_SFREE(paw_dmft%bessel)
 ABI_SFREE(paw_dmft%bessel_int)
 ABI_SFREE(paw_dmft%lpawu)
 ABI_SFREE(paw_dmft%omega_lo)
 ABI_SFREE(paw_dmft%wgt_wlo)
 ABI_SFREE(paw_dmft%slm2ylm)
 ABI_SFREE(paw_dmft%jmj2ylm)

 paw_dmft%nband => null()
 paw_dmft%dmft_shiftself => null()
 paw_dmft%dmft_nominal => null()
 paw_dmft%npwarr => null()
!   paw_dmft%wtk is only an explicit pointer =>dtset%wtk
!   if (associated(paw_dmft%wtk)) deallocate(paw_dmft%wtk)
 paw_dmft%wtk => null()
 paw_dmft%typat => null()
 paw_dmft%int_meshsz => null()

 if (allocated(paw_dmft%radgrid)) then
   call pawrad_free(paw_dmft%radgrid(:))
   ABI_FREE(paw_dmft%radgrid)
 end if
 call destroy_sc_dmft_paralkgb(paw_dmft)
 if (paw_dmft%use_dmft == 1) then
   call destroy_paral_dmft(paw_dmft%distrib)
 end if

end subroutine destroy_sc_dmft
!!***

!!****f* m_paw_dmft/print_dmft
!! NAME
!! print_dmft
!!
!! FUNCTION
!!  Print relevant data for DMFT cycle.
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!  pawprtvol=flag for print
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_dmft(paw_dmft,pawprtvol)

!Arguments ------------------------------------
 integer, intent(in) :: pawprtvol
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: iband,ifreq,ikpt,isppol
 character(len=500) :: message
! *********************************************************************

 if (abs(pawprtvol) >= 3 ) then
   write(message,'(4a,3(a,2x,e21.14,a))') &
     & "  -------------------------------------------------",ch10,&
     & "  --- Data for DMFT ",ch10,&
     & "  --- fermie     = ",paw_dmft%fermie    ,ch10,&
     & "  --- fermie_dft = ",paw_dmft%fermie_dft,ch10,&
     & "  --- temp       = ",paw_dmft%temp      ,ch10
   call wrtout(std_out,message,'COLL')
   write(message,'(7(a,15x,i8,a),a,2x,e21.14,3a)') &
     & "  --- natpawu    = ",paw_dmft%natpawu   ,ch10,&
     & "  --- dmft_iter  = ",paw_dmft%dmft_iter ,ch10,&
     & "  --- dmft_solv  = ",paw_dmft%dmft_solv ,ch10,&
     & "  --- dmft_nwlo  = ",paw_dmft%dmft_nwlo ,ch10,&
     & "  --- dmft_nwli  = ",paw_dmft%dmft_nwli ,ch10,&
     & "  --- dmft_dc    = ",paw_dmft%dmft_dc   ,ch10,&
     & "  --- dmftqmc_l  = ",paw_dmft%dmftqmc_l ,ch10,&
     & "  --- dmftqmc_n  = ",paw_dmft%dmftqmc_n ,ch10,&
     & "  -------------------------------------------------",ch10
   call wrtout(std_out,message,'COLL')

!  write(message,'(4a,3(a,2x,f8.3,a),8(a,2x,i8,a),a)') "-----------------------------------------------",ch10,&
!&   "--- Data for DMFT ",ch10,&
!&   "--- paw_dmft%fermie     = ",paw_dmft%fermie    ,ch10,&
!&   "--- paw_dmft%fermie_dft = ",paw_dmft%fermie_dft,ch10,&
!&   "--- paw_dmft%temp       = ",paw_dmft%temp      ,ch10,&
!&   "--- paw_dmft%natpawu    = ",paw_dmft%natpawu   ,ch10,&
!&   "--- paw_dmft%dmft_iter  = ",paw_dmft%dmft_iter ,ch10,&
!&   "--- paw_dmft%dmft_solv  = ",paw_dmft%dmft_solv ,ch10,&
!&   "--- paw_dmft%dmft_nwlo  = ",paw_dmft%dmft_nwlo ,ch10,&
!&   "--- paw_dmft%dmft_nwli  = ",paw_dmft%dmft_nwli ,ch10,&
!&   "--- paw_dmft%dmft_dc    = ",paw_dmft%dmft_dc   ,ch10,&
!&   "--- paw_dmft%dmftqmc_l  = ",paw_dmft%dmftqmc_l ,ch10,&
!&   "--- paw_dmft%dmftqmc_n  = ",paw_dmft%dmftqmc_n ,ch10,&
!&   "-----------------------------------------------"
   if (abs(pawprtvol) > 10) then
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') " DFT Eigenvalues "
     do isppol=1,paw_dmft%nsppol
       write(message,'(a,i4)') "--isppol--",isppol
       call wrtout(std_out,message,'COLL')
       do ikpt=1,paw_dmft%nkpt
         write(message,'(a,i4,2x,f14.5,a)') "  -k-pt--",ikpt,paw_dmft%wtk(ikpt),"(<-weight(k-pt))"
         call wrtout(std_out,message,'COLL')
         do iband=1,paw_dmft%mbandc
           write(message,'(a,i4,f10.5)') "   -iband--",iband,paw_dmft%eigen_dft(iband,ikpt,isppol)
           call wrtout(std_out,message,'COLL')
         end do ! iband
       end do ! ikpt
     end do ! isppol
     write(message,'(3x,a)') "Log. Freq"
     call wrtout(std_out,message,'COLL')
     do ifreq=1,paw_dmft%dmft_nwlo
       write(message,'(3x,a,i4,2(2x,e13.5))') "--ifreq--",ifreq,paw_dmft%omega_lo(ifreq),paw_dmft%wgt_wlo(ifreq)
       call wrtout(std_out,message,'COLL')
     end do ! ifreq
   end if ! abs(pawprtvol)>10
 end if ! abs(pawprtvol)>=3

end subroutine print_dmft
!!***

!!****f* m_paw_dmft/print_sc_dmft
!! NAME
!! print_sc_dmft
!!
!! FUNCTION
!!  Print relevant data for self-consistent
!!  DFT+DMFT cycle.
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!  pawprtvol=flag for print
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_sc_dmft(paw_dmft,pawprtvol)

!Arguments ------------------------------------
 integer, intent(in) :: pawprtvol
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: iband
 character(len=500) :: message
! *********************************************************************

 if (abs(pawprtvol) >= 3) then
   write(message,'(5a,7(a,2x,i5,a),a)') ch10, &
     & " -----------------------------------------------",ch10,&
     & " --- Data for self-consistent DFT+DMFT cycle",ch10,&
     & " --- mband             = ",paw_dmft%mband,ch10,&
     & " --- dmftbandi         = ",paw_dmft%dmftbandi,ch10,&
     & " --- dmftbandf         = ",paw_dmft%dmftbandf,ch10,&
     & " --- nb of corr. bands = ",paw_dmft%mbandc,ch10,&
     & " --- nkpt              = ",paw_dmft%nkpt,ch10,&
     & " --- nsppol            = ",paw_dmft%nsppol,ch10,&
     & " --- usedmft           = ",paw_dmft%use_dmft,ch10,&
     !& " --- use_sc_dmft = ",paw_dmft%use_sc_dmft,ch10,&
     & " -----------------------------------------------"
   call wrtout(std_out,message,'COLL')
   write(message,'(2a)') ch10," Indicating correlated bands"
   call wrtout(std_out,message,'COLL')
   write(message,'(100i5)') (iband,iband=1,min(paw_dmft%mband,100))
   call wrtout(std_out,message,'COLL')
   write(message,'(100L5)') (paw_dmft%band_in(iband),iband=1,min(paw_dmft%mband,100))
   call wrtout(std_out,message,'COLL')
   write(message,'(2a)') ch10," Correlated index   Band index"
   call wrtout(std_out,message,'COLL')
   do iband=1,paw_dmft%mbandc
     write(message,'(5x,i5,10x,i5)') iband,paw_dmft%include_bands(iband)
     call wrtout(std_out,message,'COLL')
   end do ! iband
   if (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7) then
     write(message,'(2a,i4,a)') ch10,&
        & 'The',paw_dmft%mband-paw_dmft%dmftbandf+paw_dmft%dmftbandi-1,&
        & '  Following bands are excluded from the DMFT calculation'
     call wrtout(std_out,message,'COLL')
     write(message,'(100i5)') (paw_dmft%exclude_bands(iband),iband=1,min(paw_dmft%mband-paw_dmft%dmftbandf+paw_dmft%dmftbandi-1,100))
     call wrtout(std_out,message,'COLL')
   end if
   write(message,*)
   call wrtout(std_out,message,'COLL')
 end if ! abs(pawprtvol)>=3

end subroutine print_sc_dmft
!!***

!!****f* m_paw_dmft/saveocc_dmft
!! NAME
!! saveocc_dmft
!!
!! FUNCTION
!!  save occnd on disk
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!
!! OUTPUT
!!
!! SOURCE

subroutine saveocc_dmft(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: ib,ib1,ikpt,is,nband_k,nkpt,unitsaveocc
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
! ********************************************************************

 tmpfil = trim(paw_dmft%filapp)//'_DMFTOCCND'
 if (open_file(tmpfil,message,newunit=unitsaveocc,status='unknown',form='formatted') /= 0) ABI_ERROR(message)

 rewind(unitsaveocc)
 write(message,'(2a)') ch10,"  == Print DFT+DMFT non diagonal occupations on disk"
 call wrtout(std_out,message,'COLL')
 write(message,'(5a,2x,4i5,2a)') "# DFT+DMFT off-diagonal occupations f_{ib,ib1} = <Psi^{dagger}_{ib1}|Psi_{ib}>", &
         & ch10,"# natom,nsppol,mbandc,nkpt",ch10, &
         & "####",paw_dmft%natom,paw_dmft%nsppol,paw_dmft%mbandc,paw_dmft%nkpt,ch10, &
         & "#        isppol      ikpt         ib          ib1         Re                        Imag"

 call wrtout(unitsaveocc,message,'COLL')
 nkpt = paw_dmft%nkpt
 do is=1,paw_dmft%nsppol
   do ikpt=1,nkpt
     nband_k = paw_dmft%nband(ikpt+(is-1)*nkpt)
     do ib=1,nband_k
       if ((.not. paw_dmft%band_in(ib)) .and. (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7)) cycle
       do ib1=1,nband_k
         if ((.not. paw_dmft%band_in(ib1)) .and. (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7)) cycle
         write(unitsaveocc,*) is,ikpt,ib,ib1,paw_dmft%occnd(1,ib,ib1,ikpt,is),&
                                           & paw_dmft%occnd(2,ib,ib1,ikpt,is)
       end do ! ib1
     end do ! ib
   end do ! ikpt
 end do ! is
 write(message,'(3a)') "# end of record",ch10,"####  1234 "
 call wrtout(unitsaveocc,message,'COLL')
 close(unitsaveocc)

end subroutine saveocc_dmft
!!***

!!****f* m_paw_dmft/readocc_dmft
!! NAME
!! readocc_dmft
!!
!! FUNCTION
!!  read occnd on disk
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!  filnam_ds3 = root for filname to read (input)
!!  filnam_ds4 = root for filname to read (output)
!!
!! OUTPUT
!!  paw_dmft: occnd
!!
!! SOURCE

subroutine readocc_dmft(paw_dmft,filnam_ds3,filnam_ds4)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
 character(len=fnlen), intent(in) :: filnam_ds3,filnam_ds4
!Local variables-------------------------------
 integer :: dum1,dum2,dum3,dum4,ib,ib1,ikpt,ioerr,is,nband_k,nkpt,unitsaveocc
 logical :: lexist
 character(len=4) :: chtemp
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
! *********************************************************************

 if (paw_dmft%dmft_read_occnd == 0) return
 if (paw_dmft%dmft_read_occnd == 1) tmpfil = trim(filnam_ds3)//'_DMFTOCCND'
 if (paw_dmft%dmft_read_occnd == 2) tmpfil = trim(filnam_ds4)//'_DMFTOCCND'
 inquire(file=trim(tmpfil),exist=lexist)!,recl=nrecl)
 unitsaveocc = 679
 if (lexist) then
   if (open_file(tmpfil,message,unit=unitsaveocc,status='unknown',form='formatted') /= 0) ABI_ERROR(message)
   rewind(unitsaveocc)
   write(message,'(3a)') ch10,"  == Read DMFT non diagonal occupations on disk"
   call wrtout(std_out,message,'COLL')
   read(unitsaveocc,*)
   read(unitsaveocc,*)
   read(unitsaveocc,*,iostat=ioerr) chtemp,dum1,dum2,dum3,dum4
   read(unitsaveocc,*)
   if (ioerr < 0) write(std_out,*) "read",dum1,dum2,dum3,dum4
   write(message,'(2a,4i4)') ch10,"  == natom, nsppol, nbandc, nkpt read are",dum1,dum2,dum3,dum4
   call wrtout(std_out,message,'COLL')
   nkpt = paw_dmft%nkpt
   do is=1,paw_dmft%nsppol
     do ikpt=1,nkpt
       nband_k = paw_dmft%nband(ikpt+(is-1)*nkpt)
       do ib=1,nband_k
         if ((.not. paw_dmft%band_in(ib)) .and. (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7)) cycle
         do ib1=1,nband_k
           if ((.not. paw_dmft%band_in(ib1)) .and. (paw_dmft%dmft_solv /= 6 .and. paw_dmft%dmft_solv /= 7)) cycle
           read(unitsaveocc,*) dum1,dum2,dum3,dum4,&
               & paw_dmft%occnd(1,ib,ib1,ikpt,is),paw_dmft%occnd(2,ib,ib1,ikpt,is)
         end do ! ib1
       end do ! ib
     end do ! ikpt
   end do ! is
!   write(read,'(3a)') "# end of record",ch10&
!&                ,"####  1234 "
!   call wrtout(unitsaveocc,message,'COLL')
 else
   write(message,'(2a,2x,2a)') ch10,"   File",trim(tmpfil),"is not available"
   call wrtout(std_out,message,'COLL')
   write(message,'(4a)') ch10,"  ==> DMFT Occupations not available for restart", ch10, &
       & "      -> The calculation is started with Fermi Dirac scheme for occupations"
   call wrtout(std_out,message,'COLL')
 end if ! lexist

end subroutine readocc_dmft
!!***

!!****f* m_paw_dmft/init_sc_dmft_paralkgb
!! NAME
!! init_sc_dmft_paralkgb
!!
!! FUNCTION
!!  Init some values used with KGB parallelism in self consistent DMFT
!!  calculation.
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!! mpi_enreg = information about MPI parallelization
!!
!! OUTPUT
!!  paw_dmft: bandc_proc, use_bandc
!!
!! SOURCE

subroutine init_sc_dmft_paralkgb(paw_dmft,mpi_enreg)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(MPI_type), intent(in) :: mpi_enreg
!Local variables-------------------------------
 integer :: ib,ibc,mbandc,nproc,proc
! *********************************************************************

 mbandc = paw_dmft%mbandc
 nproc  = mpi_enreg%nproc_band

 ABI_MALLOC(paw_dmft%bandc_proc,(mbandc))
 ABI_MALLOC(paw_dmft%use_bandc,(nproc))
 paw_dmft%bandc_proc(:) = 0
 paw_dmft%use_bandc(:)  = .false.

 do ibc=1,mbandc
   ib = paw_dmft%include_bands(ibc)
   proc = mod((ib-1)/mpi_enreg%bandpp,nproc)
   paw_dmft%bandc_proc(ibc) = proc
   paw_dmft%use_bandc(proc+1) = .true.
 end do ! ibc

end subroutine init_sc_dmft_paralkgb

!!***

!!****f* m_paw_dmft/destroy_sc_dmft_paralkgb
!! NAME
!! destroy_sc_dmft_paralkgb
!!
!! FUNCTION
!!   deallocate bandc_proc and use_bandc
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_sc_dmft_paralkgb(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(inout) :: paw_dmft
! *********************************************************************

 ABI_SFREE(paw_dmft%bandc_proc)
 ABI_SFREE(paw_dmft%use_bandc)

end subroutine destroy_sc_dmft_paralkgb
!!***

!!****f* m_paw_dmft/init_paral_dmft
!! NAME
!! init_paral_dmft
!!
!! FUNCTION
!!  Initialize MPI distribution
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!  distrib <type(mpi_distrib_dmft_type)> = mpi related data
!!  nfreq=number of frequencies over which to parallelize
!!
!! OUTPUT
!!
!! SOURCE

subroutine init_paral_dmft(paw_dmft,distrib,nfreq)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(mpi_distrib_dmft_type), intent(inout) :: distrib
 integer, intent(in) :: nfreq
!Local variables-------------------------------
 integer :: deltakpt,deltaw,ierr,ifreq,ikpt,irank,myproc,nfreq_proc
 integer :: nkpt,nkpt_proc,nproc,nproc_freq,nproc_kpt,residu,spacecomm
! *********************************************************************

 myproc = paw_dmft%myproc
 nkpt   = paw_dmft%nkpt
 nproc  = paw_dmft%nproc
 spacecomm = paw_dmft%spacecomm

 distrib%nw = nfreq
 distrib%shiftk = 0

 ABI_MALLOC(distrib%nw_mem,(nproc))
 ABI_MALLOC(distrib%procb,(nkpt))
 ABI_MALLOC(distrib%procf,(nfreq))
 ABI_MALLOC(distrib%proct,(nfreq))

 ! First initialize parallelization over frequencies only

 deltaw = nfreq / nproc
 residu = nfreq - deltaw*nproc

 nproc_freq = min(nfreq,nproc)

 if (nproc_freq < nproc) distrib%nw_mem(nproc_freq+1:nproc) = 0
 ifreq = 1
 do irank=0,nproc_freq-1
   nfreq_proc = merge(deltaw+1,deltaw,irank<residu)
   distrib%nw_mem(irank+1) = nfreq_proc
   distrib%procf(ifreq:ifreq+nfreq_proc-1) = irank
   ifreq = ifreq + nfreq_proc
 end do ! irank

 ! Next initialize parallelization on kpt and then frequencies

 nproc_kpt  = min(nkpt,nproc)
 nproc_freq = nproc / nkpt

 ABI_MALLOC(distrib%nkpt_mem,(nproc))
 ABI_MALLOC(distrib%nw_mem_kptparal,(nproc_freq+1))

 if (nproc > nproc_kpt) distrib%nkpt_mem(nproc_kpt+1:nproc) = 0
 distrib%nw_mem_kptparal(nproc_freq+1) = 0

 if (nproc_freq <= 1) then ! parallelization on kpt only

   deltakpt = nkpt / nproc
   residu   = nkpt - deltakpt*nproc

   distrib%nw_mem_kptparal(1) = nfreq
   distrib%proct(:) = 0

   ikpt = 1
   do irank=0,nproc_kpt-1
     nkpt_proc = merge(deltakpt+1,deltakpt,irank<residu)
     distrib%nkpt_mem(irank+1) = nkpt_proc
     distrib%procb(ikpt:ikpt+nkpt_proc-1) = irank
     if (myproc == irank) distrib%shiftk = ikpt - 1
     ikpt = ikpt + nkpt_proc
   end do ! irank

   distrib%comm_kpt = spacecomm
   distrib%me_kpt   = mod(myproc,nproc_kpt)
   distrib%me_freq  = myproc / nproc_kpt

   call xmpi_comm_split(spacecomm,distrib%me_kpt,distrib%me_freq,distrib%comm_freq,ierr)

   distrib%me_kpt = myproc

 else ! parallelization on both kpt and frequencies

   deltaw = nfreq / nproc_freq
   residu = nfreq - deltaw*nproc_freq

   distrib%nkpt_mem(1:nkpt) = 1

   do ikpt=1,nkpt
     distrib%procb(ikpt) = ikpt - 1
   end do ! ikpt

   ifreq = 1
   do irank=0,nproc_freq-1
     nfreq_proc = merge(deltaw+1,deltaw,irank<residu)
     if (nfreq_proc > 0) distrib%proct(ifreq:ifreq+nfreq_proc-1) = irank
     ifreq = ifreq + nfreq_proc
     distrib%nw_mem_kptparal(irank+1) = nfreq_proc
   end do ! irank

   distrib%me_kpt  = myproc / nproc_freq
   distrib%me_freq = mod(myproc,nproc_freq)
   distrib%shiftk  = distrib%me_kpt

   call xmpi_comm_split(spacecomm,distrib%me_freq,distrib%me_kpt,distrib%comm_kpt,ierr)

   if (myproc >= nkpt*nproc_freq) then
     distrib%me_kpt  = myproc - nkpt*nproc_freq
     distrib%me_freq = nproc_freq
   end if

   call xmpi_comm_split(spacecomm,distrib%me_kpt,distrib%me_freq,distrib%comm_freq,ierr)

   if (myproc >= nkpt*nproc_freq) distrib%me_kpt = myproc / nproc_freq

 end if ! nproc_freq<=1

end subroutine init_paral_dmft
!!***

!!****f* m_paw_dmft/destroy_paral_dmft
!! NAME
!! destroy_paral_dmft
!!
!! FUNCTION
!!  Deallocate MPI distribution
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!  distrib <type(mpi_distrib_dmft_type)> = mpi related data
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_paral_dmft(distrib)

!Arguments ------------------------------------
 type(mpi_distrib_dmft_type), intent(inout) :: distrib
!Local variables-------------------------------
! *********************************************************************

 ABI_FREE(distrib%nkpt_mem)
 ABI_FREE(distrib%nw_mem)
 ABI_FREE(distrib%nw_mem_kptparal)
 ABI_FREE(distrib%procb)
 ABI_FREE(distrib%procf)
 ABI_FREE(distrib%proct)

 call xmpi_comm_free(distrib%comm_freq)
 call xmpi_comm_free(distrib%comm_kpt)

end subroutine destroy_paral_dmft
!!***

END MODULE m_paw_dmft
!!***
