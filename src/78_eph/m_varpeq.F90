!!****m* ABINIT/m_varpeq
!! NAME
!!  m_varpeq
!!
!! FUNCTION
!!  Description
!!
!! COPYRIGHT
!!  Copyright (C) 2023-2024 ABINIT group (VV, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_varpeq

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_dtset
 use m_dtfil
 use m_crystal
 use m_ebands
 use m_errors
 use m_krank
 use netcdf
 use m_nctk
 use m_xmpi
 use m_ifc
 use m_wfd

 use defs_datatypes,    only : ebands_t, pseudopotential_type
 use m_fstrings,        only : sjoin, ktoa, ftoa, strcat, ltoa, itoa
 use m_time,            only : cwtime_report, cwtime
 use m_io_tools,        only : file_exists, iomode_from_fname, open_file
 use m_pptools,         only : write_xsf
 use m_geometry,        only : xcart2xred
 use m_kpts,            only : kpts_map, kpts_timrev_from_kptopt, bzlint_t, kptrlatt_from_ngkpt
 use m_fft_mesh,        only : supercell_fft, calc_ceikr
 use m_pawtab,          only : pawtab_type
 use m_gstore,          only : gstore_t, gqk_t
 use m_supercell,       only : supercell_type
 use m_paw_sphharm,     only : ylm_angular_mesh
 use m_fftcore,         only : ngfft_seq !, get_kg
 use m_ephtk,           only : ephtk_get_mpw_gmax
 use m_dynmat,          only : phdispl_from_eigvec
 use m_phonons,         only : pheigvec_rotate


 implicit none

 private

!!***

!----------------------------------------------------------------------

!!****t* m_varpeq/polstate_t
!! NAME
!!  polstate_t
!!
!! FUNCTION
!!  Datatype defining MPI-distributed parameters of polaronic states for for a
!!  given spin index (if collinear magnetism i.e. nsppol 2). Local variables and
!!  arrays start with `my_`, global have `*glob*` in their names. MPI-grid is
!!  inherinted from a gstore%gqk object at initialization.
!!
!! SOURCE

 type, public :: polstate_t

  character(len=fnlen) :: aseed = " "
   ! Specifies the type of initial seed for charge localization A_nk
   ! Possible values: "gau_energy", "gau_length", "random", "even"

  logical :: translate = .false.
   ! Flag controlling treatment of polaronic solution invariant by primitive
   ! translations inside supercell
   ! if .true. and np > 1, the next solution is found to be orthogonal to all
   ! previous states + their translated replicas

  integer :: spin = -1
   ! Spin index

  integer :: ip = -1
   ! Current polaronic state

  integer :: np = -1
   ! Number of polaronic states

  integer :: nkbz = -1
   ! Number of k-points in full BZ

  integer :: nqbz = -1
   ! Number of q-points in full BZ

  real(dp) :: e_frohl
   ! Long-range divergence correction of polaron binding energy due to g(0) avg

  integer :: ngkpt(3)
   ! Number of points in the uniform k-grid defining the electronic subspace

  real(dp) :: gpr_energy(2)
  ! Gaussian parameters for energy-based initialization

  real(dp) :: gpr_length(3)
  ! Gaussian parameters for localization length-based initialization

  logical(dp), allocatable :: has_prev_grad(:)
   ! (np)
   ! Flag indicating if an electroic gradient has been computed at previous step

  real(dp), allocatable :: gradres(:)
   ! (np)
   ! L^2-norm of the electronic gradient for each state

  real(dp), allocatable :: enterms(:,:)
   ! (4, np)
   ! Polaronic energy terms for each state
   ! 1 -> electronic term
   ! 2 -> vibrational term
   ! 3 -> electron-phonon term
   ! 4 -> polaron energy level

  real(dp), allocatable :: eig(:,:)
   ! (gqk%nb, gstore%ebands%nkpt)
   ! Electronic eigenstates participating in the polaron formation
   ! Correspond either to valence or conduction manifold, not both
   ! Band edge is set to 0; valence states are also inverted

  real(dp), allocatable :: my_g0(:)
   ! (gqk%my_npert)
   ! Long-range correction for electron-phonon matrix elements at q=\Gamma

  real(dp), allocatable :: my_qpts(:,:)
   ! (3, gqk%my_nq)
   ! q-points treated by this MPI proc

  real(dp), pointer :: my_kpts(:,:) => null()
   ! (3, gqk%my_nk)
   ! k-points treated by this MPI proc, points to gqk%my_kpts(:,:)

  complex(dp), allocatable :: my_a(:,:,:)
   ! (gqk%nb, gqk%my_nk, np)
   ! Electronic coefficients A_nk for each state treated by this MPI proc

  complex(dp), allocatable :: a_glob(:,:)
   ! (gqk%nb, gqk%glob_nk)
   ! Global array of electronic coefficients A_nk at current state

  complex(dp), allocatable :: my_b(:,:,:)
   ! (gqk%my_npert, gqk%my_nq, np)
   ! Vibrational coefficients B_q\nu for each state treated by this MPI proc

  complex(dp), allocatable :: my_pc(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Preconditioner at current state treated by this MPI proc

  complex(dp), allocatable :: my_grad(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Electronic gradient D_nk at current state treated by this MPI proc
   ! orthogonal to states

  complex(dp), allocatable :: my_prev_grad(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Previous electronic gradient D_nk at current state, orthogonal to states

  complex(dp), allocatable :: my_pcgrad(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Preconditioned gradient at current state, orthogonal to states

  complex(dp), allocatable :: my_prev_pcgrad(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Previous preconditioned gradient at current state, orthogonal to states

  complex(dp), allocatable :: my_pcjgrad(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Preconditioned conjugate gradient at current state treated by this MPI proc
   ! orthogonal to current state & normalized

  complex(dp), allocatable :: my_prev_pcjgrad(:,:)
   ! (gqk%nb, gqk%my_nk)
   ! Previous preconditioned conjugate gradient at current state

  complex(dp), allocatable :: pcjgrad_glob(:,:)
   ! (gqk%nb, gqk%glob_nk)
   ! Global preconditioned conjugate gradient at current state
   ! orthogonal to current state & normalized

  class(gqk_t), pointer :: gqk => null()
   ! Datastructure storing e-ph matrix elements treated by this MPI proc

  type(krank_t) :: krank_kpts
   ! Object used to find k-points in BZ

  type(krank_t) :: krank_qpts
   ! Object used to find q-points in BZ

  contains

    procedure :: setup => polstate_setup
    ! Setup optimization proccess by specifing initial electronic vector A_nk

    procedure :: localize => polstate_localize
    ! Localize polaron at current state. From A_nk calculate:
    ! vibrational coefficients B_q\nu, energy terms, polaron level

    procedure :: get_enel => polstate_get_enel
    ! Calculate and return electronic energy term

    procedure :: get_enph => polstate_get_enph
    ! Caclulate and return vibrational energy term

    procedure :: get_enelph => polstate_get_enelph
    ! Calculate and return electron-phobnon energy term

    procedure :: get_lm_theta => polstate_get_lm_theta
    ! Calculate and return line minimization parameter \theta

    procedure :: calc_grad => polstate_calc_grad
    ! Calculate steepest descent vector

    procedure :: calc_pcjgrad => polstate_calc_pcjgrad
    ! Calculate preconditioned conjugate gradient direction

    procedure :: update_pc => polstate_update_pc
    ! Update preconditioner

    procedure :: update_a => polstate_update_a
    ! Update array of electronic coefficients

    procedure :: ort_to_states => polstate_ort_to_states
    ! Orthogonalize a given vector to a set of polaronic states

    procedure :: calc_b_from_a => polstate_calc_b_from_a
    ! Calculate vibrational coefficients B_q\nu from a known set of electronic
    ! coefficients A_nk

    !procedure :: get_b_from_displ =>_get_b_from_displ

    procedure :: seed_a => polstate_seed_a
    ! Seed an initial vector of electronic coefficients A_nk

    procedure :: load_a => polstate_load_a
    ! Initialize A_nk from an existent vector of electronic coefficients

    procedure :: get_sqnorm => polstate_get_sqnorm
    ! Helper function to compute squared L^2-norm of MPI-distributed array

    procedure :: gather => polstate_gather
    ! Helper function to gather MPI-distributed array into a global one

    procedure :: get_krank_glob => polstate_get_krank_glob
    ! Helper function to calculate global krank objects

    procedure :: free => polstate_free
    ! Free memory

 end type polstate_t
!!***

!----------------------------------------------------------------------

!!****t* m_varpeq/varpeq_t
!! NAME
!!  varpeq_t
!!
!! FUNCTION
!!  Variational Polaron Equations datatype. Stores variational optimization
!!  parameters, polaronic states and provides higher-level procedures required
!!  for variational optimization and related input and output.
!!
!! SOURCE

 type, public :: varpeq_t

   character(len=fnlen) :: pkind = " "
   ! Specifies the kind of polaron
   ! Possible values: "hole", "electron"

   character(len=fnlen) :: aseed = " "
   ! Specifies the type of initial seed for charge localization A_nk
   ! Possible values: "gau_energy", "gau_length", "random", "even"

   logical :: is_complete = .False.
   ! Flag indicating if the datatype is completely or partially initialized
   ! Required to distinguish between newly created and loaded-from-disk datatype

   logical :: restart = .False.
   ! Flag to check if a restart from a *VARPEQ.nc file is needed

   logical :: interp = .False.
   ! Flag to check if an interpolation from a *VARPEQ.nc file is needed

   logical :: ld_flag = .False.
   ! Flag indicating if internal variables have been loaded from source

   logical :: g0_flag = .True.
   ! Flag indicating if avarage of el-ph matrix elements at Gamma is computed

   integer :: ncid = nctk_noid
   ! Netcdf file handle used to save results

   integer :: nstep = -1
   ! Maximum number of iterations for optimization of a single polaronic state

   integer :: nstep_ort = -1
   ! Maximum number of iterations for which orthogonalization to all previous
   ! states is performed

   integer :: nsppol = -1
   ! Number of independent spin polarizations

   integer :: nstates = -1
   ! Number of polaronic states to be optimized for each spin polarization

   integer :: natom3 = -1
   ! 3*gstore%cryst%natom
   ! Number of atomic perturbations, mainly used to dimensionalize arrays

   integer :: max_nk = -1
   ! Maximum number of k-points among between independent spin polarizations

   integer :: max_nq = -1
   ! Maximum number of q-points among between independent spin polarizations

   integer :: max_nb = -1
   ! Maximum number of bands among between independent spin polarizations

   integer :: frohl_ntheta = -1
   ! Number of angular mesh division for spherical integration of long-range
   ! divergence of electron-phonon matrix elements

   real(dp) :: e_frohl
   ! Long-range divergence correction of polaron binding energy

   real(dp) :: tolgrs
   ! L^2 gradient norm tolerance

   integer :: ngkpt(3)
   ! Number of points in the uniform k-grid defining the electronic subspace

   integer, allocatable :: nk_spin(:)
   ! (nsppol)
   ! Number of k-points for each spin polarization

   integer, allocatable :: nq_spin(:)
   ! (nsppol)
   ! Number of q-points for each spin polarization

   integer, allocatable :: nb_spin(:)
   ! (nsppol)
   ! Number of bands for each spin polarization

   integer, allocatable :: brange_spin(:,:)
   ! (2, nsppol)
   ! Number of bands for each spin polarization

   integer, allocatable :: cvflag_spin(:,:)
   ! (nstates, nsppol)
   ! Convergence flags at each state for each spin:
   ! 0 --> calculation is not converged
   ! 1 --> calculation is converged

   integer, allocatable :: nstep2cv_spin(:,:)
   ! (nstates, nsppol)
   ! Number of steps to convergence at each state for each spin

   real(dp), allocatable :: scf_hist_spin(:,:,:,:)
   ! (6, nstep, nstates, nsppol)
   ! SCF optimization history at each state for each spin

   real(dp), allocatable :: kpts_spin(:,:,:)
   ! (3, max_nk, nsppol)
   ! k-points for each spin

   real(dp), allocatable :: qpts_spin(:,:,:)
   ! (3, max_nk, nsppol)
   ! q-points for each spin

   complex(dp), allocatable :: a_spin(:,:,:,:)
   ! (max_nb, max_nk, nstates, nsppol)
   ! Optimized electronic coefficients A_nk at each state for each spin

   complex(dp), allocatable :: b_spin(:,:,:,:)
   ! (natom3, max_nq, nstates, nsppol)
   ! Optimized vibrational coefficients B_q\nu at each state for each spin

   class(gstore_t), pointer :: gstore => null()
   ! Object storing el-ph matrix elements and other related quantities

   type(crystal_t) :: cryst
   ! Object storing information on crystal structure & symmetries

   type(gaps_t) :: gaps
   ! Object used to get information on bandgap

   type(polstate_t), allocatable :: polstate(:)
   ! (nsppol)
   ! Datatype providing data and and lower-level methods for polaronic states
   ! at each spin polarization

 contains

    procedure :: init => varpeq_init
    ! Initialize object

    procedure :: load => varpeq_load
    ! Load the initial electronic vector from a *VARPEQ.nc netcdf file

    procedure :: solve => varpeq_solve
    ! Solve variational polaron equations at each polaronic state for each spin

    procedure :: record => varpeq_record
    ! Record the current SCF optimization data

    procedure :: collect => varpeq_collect
    ! Collect SCF optimization results from each spin

    procedure :: print => varpeq_print
    ! Output SCF optimization final results

    procedure :: ncwrite => varpeq_ncwrite
    ! Save results to a *VARPEQ.nc netcdf file

    procedure :: ncread => varpeq_ncread
    ! Initialize an incomplete object from a *VARPEQ.nc netcdf file

    procedure :: compare => varpeq_compare
    ! Compare basic dimensions with another instance of varpeq_t datatype

    procedure :: calc_fravg => varpeq_calc_fravg
    ! Calculate average Fr\"ohlich (long-range) contribution to the polaron
    ! binding energy & average of electron-phonon matrix elements at q=\Gamma

    procedure :: free => varpeq_free
    ! Free memory

 end type varpeq_t
!!***

 public :: varpeq_run
   ! Main entry point

 public :: varpeq_plot
  ! Compute polaron wavefunctions and atomic displacements in the supercell and write results to XSF files

contains !=====================================================================


!!****f* m_varpeq/varpeq_run
!! NAME
!!  varpeq_run
!!
!! FUNCTION
!!  Higher-level subroutine that solves varitaional polaron equations, produces
!!  neccessary output and writes results to a *VARPEQ.nc file.
!!
!! INPUTS
!!  gstore<gstore_t>=Electron-phonon matrix elements and related quantities.
!!  dtset<dataset_type>=All input variables for this dataset.
!!  dtfil<datafiles_types>=Variables related to files.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_run(gstore, dtset, dtfil)

!Arguments ------------------------------------
 class(gstore_t), intent(in) :: gstore
 type(dataset_type), intent(in) :: dtset
 type(datafiles_type), intent(in) :: dtfil

!Local variables-------------------------------
 type(varpeq_t) :: vpq
!----------------------------------------------------------------------

 call vpq%init(gstore, dtset)
 if (vpq%frohl_ntheta > 0) call vpq%calc_fravg(avg_g0=vpq%g0_flag)
 if (vpq%interp .or. vpq%restart) call vpq%load(dtfil, dtset%varpeq_select)
 call vpq%solve()
 call vpq%print()
 call vpq%ncwrite(dtset, dtfil)
 call vpq%free()

end subroutine varpeq_run
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_free
!! NAME
!!  varpeq_free
!!
!! FUNCTION
!!  Free dynamic memory
!!
!! SOURCE

subroutine varpeq_free(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
 integer :: my_is

!----------------------------------------------------------------------

 ! Free allocatable arrays
 ! integer
 ABI_SFREE(self%nk_spin)
 ABI_SFREE(self%nq_spin)
 ABI_SFREE(self%nb_spin)
 ABI_SFREE(self%brange_spin)
 ABI_SFREE(self%cvflag_spin)
 ABI_SFREE(self%nstep2cv_spin)
 ! real
 ABI_SFREE(self%scf_hist_spin)
 ABI_SFREE(self%kpts_spin)
 ABI_SFREE(self%qpts_spin)
 ! complex
 ABI_SFREE(self%a_spin)
 ABI_SFREE(self%b_spin)

 ! Free local datatypes
 call self%cryst%free()

 ! Close netcdf file
 if (self%ncid /= nctk_noid) then
   NCF_CHECK(nf90_close(self%ncid))
   self%ncid = nctk_noid
 end if

 ! If entry is completely initalized (e.g. from self%init call), free remaining
 ! datatypes and nullify pointers
 if (self%is_complete) then
   call self%gaps%free()
   do my_is=1,self%gstore%my_nspins
     call self%polstate(my_is)%free()
   enddo
   self%gstore => null()
 endif

end subroutine varpeq_free
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_compare
!! NAME
!!  varpeq_compare
!!
!! FUNCTION
!!  Compares the current instance of varpeq datatype with another one (possibly
!!  incomplete).
!!
!! INPUTS
!!  other<varpeq_t>=Varpeq datatype to compare with.
!!  bz_mismatch [optional]=if .True. mismatch between BZ sampling is alowed
!!    (required for comparison prior to an interpolation)
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_compare(self, other, bz_mismatch)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(in) :: self, other
 logical, optional, intent(in) :: bz_mismatch

!Local variables-------------------------------
 integer :: ierr
 real(dp) :: cpu, wall, gflops

!----------------------------------------------------------------------

 ! TODO: provide more flexibility for nstates

 call cwtime(cpu, wall, gflops, "start")

 ierr = 0

 ! Compare basic dimensions
 call check_(self%pkind == other%pkind, "Difference found in pkind.")
 call check_(self%nsppol == other%nsppol, "Difference found in nsppol.")
 !call check_(self%nstates == other%nstates, "Difference found in nstates.")
 call check_(self%cryst%compare(other%cryst) == 0, &
   "Difference found in cryst.")
 call check_(all(self%brange_spin == other%brange_spin), &
   "Difference found in brange_spin.")

 ! If bz_mismatch is not allowed, also compare k/q-grids
 if (present(bz_mismatch)) then
   if (.not. bz_mismatch) then
     call check_(all(abs(self%kpts_spin - other%kpts_spin) < tol6), &
       "Difference found in kpts_spin.")
     call check_(all(abs(self%qpts_spin - other%qpts_spin) < tol6), &
     "Difference found in nq_spin.")
   end if
 endif

 ABI_CHECK(ierr == 0, "Error in varpeq_compare, see previous messages!")

 call cwtime_report(" varpeq_compare", cpu, wall, gflops)

 contains
 subroutine check_(cond, msg)
   logical, intent(in) :: cond
   character(len=*), intent(in) :: msg
   ABI_CHECK_NOSTOP(cond, msg, ierr)
 end subroutine check_

end subroutine varpeq_compare
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_ncread
!! NAME
!!  varpeq_ncread
!!
!! FUNCTION
!!  Reads basic dimensions of varpeq_t datatype from a *VARPEQ.nc netcdf file.
!!
!! INPUTS
!!  path=Path a *VARPEQ.nc file to be read.
!!  comm=MPI communicator.
!!  keep_open [optional]=if .True. keep the nc file handle open for further
!!    reading. Default: .False.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_ncread(self, path, comm, keep_open)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 character(len=fnlen), intent(in) :: path
 integer, intent(in) :: comm
 logical, optional, intent(in) :: keep_open

!Local variables-------------------------------
 integer :: ncid, nsppol, natom3, nstates
 real(dp) :: cpu, wall, gflops
 real(dp), ABI_CONTIGUOUS pointer :: rpt_d5(:,:,:,:,:)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 ABI_CHECK(file_exists(path), sjoin(" varpeq_ncread: cannot find *VARPEQ.nc file", path))

 NCF_CHECK(nctk_open_read(ncid, path, comm))

 ! Read crystal structure
 call self%cryst%ncread(ncid)

 ! Read varpeq dimensions
 NCF_CHECK(nctk_get_dim(ncid, "nsppol", self%nsppol))
 NCF_CHECK(nctk_get_dim(ncid, "nstates", self%nstates))
 NCF_CHECK(nctk_get_dim(ncid, "natom3", self%natom3))
 NCF_CHECK(nctk_get_dim(ncid, "max_nk", self%max_nk))
 NCF_CHECK(nctk_get_dim(ncid, "max_nq", self%max_nq))
 NCF_CHECK(nctk_get_dim(ncid, "max_nb", self%max_nb))
 nsppol = self%nsppol
 nstates = self%nstates
 natom3 = self%natom3

 ! Read data
 ! Static arrays
 NCF_CHECK(nf90_get_var(ncid, vid("varpeq_pkind"), self%pkind))
 NCF_CHECK(nf90_get_var(ncid, vid("ngkpt"), self%ngkpt))

 ! Allocatable arrays
 ABI_MALLOC(self%cvflag_spin, (nstates, nsppol))
 ABI_MALLOC(self%nk_spin, (nsppol))
 ABI_MALLOC(self%nq_spin, (nsppol))
 ABI_MALLOC(self%nb_spin, (nsppol))
 ABI_MALLOC(self%brange_spin, (2, nsppol))
 ABI_MALLOC(self%kpts_spin, (3, self%max_nk, nsppol))
 ABI_MALLOC(self%qpts_spin, (3, self%max_nq, nsppol))
 ABI_MALLOC(self%a_spin, (self%max_nb, self%max_nk, nstates, nsppol))
 ABI_MALLOC(self%b_spin, (natom3, self%max_nq, nstates, nsppol))

 ! integer
 NCF_CHECK(nf90_get_var(ncid, vid("cvflag_spin"), self%cvflag_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("nk_spin"), self%nk_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("nq_spin"), self%nq_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("nb_spin"), self%nb_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("brange_spin"), self%brange_spin))
 ! real
 NCF_CHECK(nf90_get_var(ncid, vid("kpts_spin"), self%kpts_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("qpts_spin"), self%qpts_spin))
 ! complex
 ! TODO: is it possible to encapsulate this trick as some abstraction?
 call c_f_pointer(c_loc(self%a_spin), rpt_d5, &
   [2, self%max_nb, self%max_nk, nstates, nsppol])
 NCF_CHECK(nf90_get_var(ncid, vid("a_spin"), rpt_d5))

 call c_f_pointer(c_loc(self%b_spin), rpt_d5, &
   [2, natom3, self%max_nq, nstates, nsppol])
 NCF_CHECK(nf90_get_var(ncid, vid("b_spin"), rpt_d5))

 if (present(keep_open)) then
   if (keep_open) self%ncid = ncid
 else
   NCF_CHECK(nf90_close(ncid))
   self%ncid = nctk_noid
 end if

 self%is_complete = .False.

 call cwtime_report(" varpeq_ncread", cpu, wall, gflops)

!----------------------------------------------------------------------

 contains
  integer function vid(var_name)
    character(len=*),intent(in) :: var_name
    vid = nctk_idname(ncid, var_name)
  end function vid

end subroutine varpeq_ncread
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_ncwrite
!! NAME
!!  varpeq_ncwrite
!!
!! FUNCTION
!!  Dump varpeq variables in a newly created *VARPEQ.nc netcdf file.
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  dtfil<datafiles_types>=Variables related to files.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_ncwrite(self, dtset, dtfil)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 type(dataset_type), intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil

!Local variables-------------------------------
 character(len=fnlen) :: path
 integer, parameter :: master = 0
 integer :: my_rank, ncid, ncerr
 real(dp) :: cpu, wall, gflops
 real(dp), ABI_CONTIGUOUS pointer :: rpt_d5(:,:,:,:,:)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 my_rank = xmpi_comm_rank(self%gstore%comm)

 ! Create netcdf file (only master works, HDF5 + MPI-IO can be handled after
 ! by reopening the file inside ncwrite_comm)
 path = strcat(dtfil%filnam_ds(4), "_VARPEQ.nc")
 if (my_rank == master) then
   ! Master creates the netcdf file used to store the data.
   NCF_CHECK(nctk_open_create(self%ncid, path, xmpi_comm_self))
   ncid = self%ncid

   ! Write the crystal (TR & invsersion symmetry only) & ebands dataset_type
   NCF_CHECK(self%cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(self%gstore%ebands, ncid))

   ! Add varpeq dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nstep", self%nstep), nctkdim_t("nsppol", self%nsppol), &
     nctkdim_t("nstates", self%nstates), nctkdim_t("natom3", self%natom3), &
     nctkdim_t("max_nk", self%max_nk), nctkdim_t("max_nq", self%max_nq), &
     nctkdim_t("max_nb", self%max_nb)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define scalars
   ! integers
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "nkbz", "nqbz", "frohl_ntheta"])
   NCF_CHECK(ncerr)
   ! real
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "varpeq_tolgrs", "e_frohl"])
   NCF_CHECK(ncerr)

   ! Define arrays with results
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("varpeq_pkind", "c", "fnlen"), &
     nctkarr_t("varpeq_aseed", "c", "fnlen"), &
     nctkarr_t("ngkpt", "int", "three"), &
     nctkarr_t("gstore_ngqpt", "int", "three"), &
     nctkarr_t("nk_spin", "int", "nsppol"), &
     nctkarr_t("nq_spin", "int", "nsppol"), &
     nctkarr_t("nb_spin", "int", "nsppol"), &
     nctkarr_t("brange_spin", "int", "two, nsppol"), &
     nctkarr_t("cvflag_spin", "int", "nstates, nsppol"), &
     nctkarr_t("nstep2cv_spin", "int", "nstates, nsppol"), &
     nctkarr_t("scf_hist_spin", "dp", "six, nstep, nstates, nsppol"), &
     nctkarr_t("kpts_spin", "dp", "three, max_nk, nsppol"), &
     nctkarr_t("qpts_spin", "dp", "three, max_nq, nsppol"), &
     nctkarr_t("cb_min_spin", "dp", "nsppol"), &
     nctkarr_t("vb_max_spin", "dp", "nsppol"), &
     nctkarr_t("a_spin", "dp", "two, max_nb, max_nk, nstates, nsppol"), &
     nctkarr_t("b_spin", "dp", "two, natom3, max_nq, nstates, nsppol") &
   ])
   NCF_CHECK(ncerr)

   ! Write data
   NCF_CHECK(nctk_set_datamode(ncid))
   ! Scalars
   ! integer
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "nkbz", "nqbz", "frohl_ntheta"], &
     [dtset%eph_task, self%gstore%nkbz, self%gstore%nqbz, self%frohl_ntheta])
   NCF_CHECK(ncerr)
   ! real
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "varpeq_tolgrs", "e_frohl"], &
     [self%tolgrs, self%e_frohl])
   NCF_CHECK(ncerr)

   ! Arrays
   ! character
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "varpeq_pkind"), self%pkind))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "varpeq_aseed"), self%aseed))
   ! integer
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cvflag_spin"), self%cvflag_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngkpt"), self%ngkpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gstore_ngqpt"), self%gstore%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nk_spin"), self%nk_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nq_spin"), self%nq_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nb_spin"), self%nb_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "brange_spin"), self%brange_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cvflag_spin"), self%cvflag_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nstep2cv_spin"), self%nstep2cv_spin))
   ! real
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "scf_hist_spin"), self%scf_hist_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kpts_spin"), self%kpts_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpts_spin"), self%qpts_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cb_min_spin"), self%gaps%cb_min))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vb_max_spin"), self%gaps%vb_max))
   ! complex
   call c_f_pointer(c_loc(self%a_spin), rpt_d5, &
     [2, self%max_nb, self%max_nk, self%nstates, self%nsppol])
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "a_spin"), rpt_d5))

   call c_f_pointer(c_loc(self%b_spin), rpt_d5, &
     [2, self%natom3, self%max_nq, self%nstates, self%nsppol])
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "b_spin"), rpt_d5))
 end if ! master

 call xmpi_barrier(self%gstore%comm)
 call cwtime_report(" varpeq: ncwrite", cpu, wall, gflops)

end subroutine varpeq_ncwrite
!!***

!!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_print
!! NAME
!!  varpeq_print
!!
!! FUNCTION
!!  Output SCF optimization results
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_print(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
!scalars
 character(len=5000) :: msg
 integer, parameter :: master = 0
 integer :: my_rank, spin, ip, ii
 real(dp) :: enpol, enel, enph, enelph, eps, grs
!arrays
 integer :: units(2)

!----------------------------------------------------------------------

 my_rank = xmpi_comm_rank(self%gstore%comm)

 units = [std_out, ab_out]
 if (my_rank == master) then
   write(msg, '(a)') " "
   call wrtout(units, msg)
   write(msg, '(a)') "  === Variational Polaron Equations ==="
   call wrtout(units, msg)
   write(msg, '(a)') repeat('-', 80)
   call wrtout(units, msg)

   ! For each spin, print the SCF history at each state
   do spin=1,self%gstore%nsppol
     write(msg, '(a,i1,a,i1)') "  * spin: ", spin, "/", self%gstore%nsppol
     call wrtout(units, msg)

     write(msg,'(a4,a13,2a12,a13,a13,a13)') 'Step', 'E_pol', 'E_el', 'E_ph', &
       'E_elph', 'epsilon', '||gradient||'
     call wrtout(units, msg)

     do ip=1,self%nstates
       write(msg, '(a,i2,a,i2)') "  * pstate: ", ip, "/", self%nstates
       call wrtout(units, msg)

       do ii=1,self%nstep2cv_spin(ip, spin)
         enpol = self%scf_hist_spin(1, ii, ip, spin)
         enel = self%scf_hist_spin(2, ii, ip, spin)
         enph = self%scf_hist_spin(3, ii, ip, spin)
         enelph = self%scf_hist_spin(4, ii, ip, spin)
         eps = self%scf_hist_spin(5, ii, ip, spin)
         grs = self%scf_hist_spin(6, ii, ip, spin)

         write(msg,'(i4,es13.4,2es12.4,es13.4,es13.4,es13.4)') ii, enpol, enel, &
           enph, enelph, eps, grs
         call wrtout(units, msg)
       enddo
       call wrtout(units, repeat('-', 80))
       write(msg, '(a13,es17.8)') "  E_pol (eV):", enpol*Ha_eV
       call wrtout(units, msg)
       write(msg, '(a13,es17.8)') "  eps (eV):", eps*Ha_eV
       call wrtout(units, msg)
       call wrtout(units, repeat('-', 80))
     enddo
     !write(msg, '(a)') repeat('-', 80)
     !call wrtout(units, msg)

   enddo
 endif

end subroutine varpeq_print
!!***

!!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_collect
!! NAME
!!  varpeq_collect
!!
!! FUNCTION
!!  Collect SCF optimization results for each spin
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_collect(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 class(polstate_t), pointer :: polstate
 integer :: ierr
 integer :: my_is, spin, my_ik, ik_glob, ip
 integer :: my_iq, iq_glob, my_pert, pert_glob
 integer :: oc_scf, oc_a, oc_b, oc_k, oc_q

!----------------------------------------------------------------------

 ! Gather the SCF process evolution data
 call xmpi_sum(self%cvflag_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%scf_hist_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%nstep2cv_spin, self%gstore%comm, ierr)

 ! Gather electron/phonon vectors and k/q points
 self%a_spin(:,:,:,:) = zero
 self%b_spin(:,:,:,:) = zero
 self%qpts_spin(:,:,:) = zero
 self%kpts_spin(:,:,:) = zero
 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   gqk => self%gstore%gqk(my_is)
   polstate => self%polstate(spin)

   do ip=1,self%nstates

     ! electronic vector
     do my_ik=1,gqk%my_nk
       ik_glob = gqk%my_kstart + my_ik - 1
       self%a_spin(:, ik_glob, ip, spin) = polstate%my_a(:, my_ik, ip)
     enddo

     ! phonon vector
     do my_iq=1,gqk%my_nq
       iq_glob = gqk%my_qstart + my_iq - 1
       do my_pert=1,gqk%my_npert
         pert_glob = gqk%my_pert_start + my_pert - 1
         self%b_spin(pert_glob, iq_glob, :, spin) = polstate%my_b(my_pert, my_iq, :)
       enddo
     enddo

   enddo

   ! k-points
   do my_ik=1,gqk%my_nk
     ik_glob = gqk%my_kstart + my_ik - 1
     self%kpts_spin(:, ik_glob, spin) = polstate%my_kpts(:, my_ik)
   enddo

   ! q-points
   do my_iq=1,gqk%my_nq
     iq_glob = gqk%my_qstart + my_iq - 1
     self%qpts_spin(:, iq_glob, spin) = polstate%my_qpts(:, my_iq)
   enddo

 enddo
 call xmpi_sum(self%a_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%b_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%kpts_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%qpts_spin, self%gstore%comm, ierr)

 ! Hack to mimic the summation over a non-existing spin commnicator
 ! Divide by the number of times we overcount, as we use the global communicator
 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   gqk => self%gstore%gqk(my_is)
   polstate => self%polstate(spin)

   oc_scf = gqk%comm%nproc
   oc_a = gqk%qpt_pert_comm%nproc
   oc_b = gqk%kpt_comm%nproc
   oc_k = oc_a
   oc_q = oc_b * gqk%pert_comm%nproc

   self%cvflag_spin(:,spin) = self%cvflag_spin(:,spin) / oc_scf
   self%scf_hist_spin(:,:,:,spin) = self%scf_hist_spin(:,:,:,spin) / oc_scf
   self%nstep2cv_spin(:,spin) = self%nstep2cv_spin(:,spin) / oc_scf
   self%a_spin(:,:,:,spin) = self%a_spin(:,:,:,spin) / oc_a
   self%b_spin(:,:,:,spin) = self%b_spin(:,:,:,spin) / oc_b
   self%kpts_spin(:,:,spin) = self%kpts_spin(:,:,spin) / oc_k
   self%qpts_spin(:,:,spin) = self%qpts_spin(:,:,spin) / oc_q
 enddo

end subroutine varpeq_collect
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_load
!! NAME
!!  varpeq_load
!!
!! FUNCTION
!!  Load and (optionally) interpolate the initial electronic vector A_nk from
!!  a *VARPEQ.nc netcdf file. Store result in the self%a_spin variable.
!!
!! INPUTS
!!  dtfil<datafiles_types>=Variables related to files.
!!  pselect=Which state to select for reload/interpolation. Non-positive value
!!    selects all states.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_load(self, dtfil, pselect)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 type(datafiles_type), intent(in) :: dtfil
 integer :: pselect

!Local variables-------------------------------
!scalars
 type(varpeq_t) :: vpq_ld
 type(bzlint_t) :: bzlint
 logical :: single_state
 integer, parameter :: master = 0
 integer :: my_rank, comm, ierr
 integer :: spin, ip, nk, nb, ik, ib
 real(dp) :: cpu, wall, gflops
!arrays
 integer :: units(2)
 real(dp) :: kpt(3)
 real(dp), allocatable :: ak(:), kpts_ld(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: rpt_d2(:,:)
 complex(dp), allocatable, target :: a_ld(:,:)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 units = [std_out, ab_out]
 comm = self%gstore%comm; my_rank = xmpi_comm_rank(comm)

 ! Read A_nk from file. Only the master processor reads, then broadcasts the data
 if (my_rank == master) then
   call vpq_ld%ncread(dtfil%filvarpeqin, xmpi_comm_self, keep_open=.False.)

   ! Consitency check
   call self%compare(vpq_ld, bz_mismatch=self%interp)
   if (pselect > 0) then
     ABI_CHECK(self%nstates == 1, "varpeq_pstates must be 1 if varpeq_select > 0.")
     ABI_CHECK(pselect <= vpq_ld%nstates, "varpeq_select must be <= loaded nstates.")
     single_state = .true.
   else
     ABI_CHECK(self%nstates == vpq_ld%nstates, "Diefference found in nstates.")
     single_state = .false.
   endif

   self%a_spin(:,:,:,:) = zero

   if (self%interp) then ! Interpolation
     call wrtout(units, " - interpolating previous A_nk")
     do spin=1,self%nsppol
       ! Setting basic dimensions & arrays
       nk = vpq_ld%nk_spin(spin); nb = vpq_ld%nb_spin(spin)
       ABI_MALLOC(kpts_ld, (3, nk))
       ABI_MALLOC(a_ld, (nb, nk))
       ABI_MALLOC(ak, (2*nb))
       kpts_ld(:,:) = vpq_ld%kpts_spin(:, 1:nk, spin)

       ! Here, interpolation is performed
       do ip=1,self%nstates

         if (single_state) then
           a_ld(:,:) = vpq_ld%a_spin(:, 1:nk, pselect, spin)
         else
           a_ld(:,:) = vpq_ld%a_spin(:, 1:nk, ip, spin)
         endif

         call c_f_pointer(c_loc(a_ld), rpt_d2, [2*nb, nk])

         call bzlint%init(vpq_ld%ngkpt, 2*nb, nk, kpts_ld, rpt_d2)

         do ik=1,self%nk_spin(spin)
           kpt(:) = self%kpts_spin(:, ik, spin)
           call bzlint%interp(kpt, ak)

           do ib=1,nb
             self%a_spin(ib, ik, ip, spin) = ak(2*ib-1) + j_dpc*ak(2*ib)
           enddo
         enddo
         call bzlint%free()

       enddo
       ABI_FREE(kpts_ld)
       ABI_FREE(a_ld)
       ABI_FREE(ak)
     enddo

   else ! Restart
     call wrtout(units, " - restarting from previous A_nk")
     self%a_spin(:,:,:,:) = vpq_ld%a_spin(:,:,:,:)
   endif

   call vpq_ld%free()
 endif

 call xmpi_bcast(self%a_spin, master, comm, ierr)
 self%ld_flag = .True.

 call cwtime_report(" varpeq: load", cpu, wall, gflops)

end subroutine varpeq_load
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_solve
!! NAME
!!  varpeq_solve
!!
!! FUNCTION
!!  Solve the Variational Polaron Equations for each spin for self%nstates
!!  polaronic states.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_solve(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
 class(polstate_t), pointer :: polstate
 !logical :: ld_flag
 integer :: my_is, spin, ip, ii
 real(dp) :: grad_sqnorm
 real(dp) :: cpu, wall, gflops

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 self%scf_hist_spin(:,:,:,:) = zero

 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   polstate => self%polstate(my_is)

   do ip=1,self%nstates
     ! initialize A_nk at this state, orthogonalize to the previous ones
     ! and normalize
     call polstate%setup(ip, a_src=self%a_spin(:,:,ip,spin), load=self%ld_flag)

     do ii=1,self%nstep
       ! gather A, get B_qnu, get energies
       call polstate%localize(ip)

       ! get bare gradient
       call polstate%calc_grad(ip)

       ! calculate and save the L^2 gradient norm
       grad_sqnorm = polstate%get_sqnorm("grad", ip)
       polstate%gradres(ip) = sqrt(grad_sqnorm)

       ! record the energies & gradient norm to varepq datatype
       call self%record(ii, ip, my_is)

       ! check if gradient norm is lower than convergence threshold
       if (polstate%gradres(ip) < self%tolgrs) then
         self%cvflag_spin(ip, spin) = 1
         exit
       endif

       ! calculate the preconditioner
       call polstate%update_pc(ip)

       ! get preconditioned conjugate gradient direction
       call polstate%calc_pcjgrad(ip, ii, self%nstep_ort)

       ! update a based on line minimization and pcj direction
       call polstate%update_a(ip)

     enddo
   enddo
 enddo

 ! Collect results from each polstate
 call self%collect()

 call cwtime_report(" varpeq: solve", cpu, wall, gflops)

end subroutine varpeq_solve
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_record
!! NAME
!!  varpeq_record
!!
!! FUNCTION
!!  Record variational polaron equations results for current iteration.
!!  Used to propagate the SCF cycle history from self%polstate to the
!!  varpeq datatype itself.
!!
!! INPUTS
!!  iter=Current iteration.
!!  ip=Index of a polaronic state.
!!  my_is=Spin polarization treated by this MPI proc.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_record(self, iter, ip, my_is)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 integer, intent(in) :: iter, my_is
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(polstate_t), pointer :: polstate
 integer :: spin, psign
 real(dp) :: enel, enph, enelph, eps

!----------------------------------------------------------------------

 spin = self%gstore%my_spins(my_is)
 polstate => self%polstate(my_is)

 psign = 1
 if (self%pkind == "hole") psign = -1

 enel = polstate%enterms(1, ip); enph = polstate%enterms(2, ip)
 enelph = polstate%enterms(3, ip); eps = polstate%enterms(4, ip)

 self%scf_hist_spin(1, iter, ip, spin) = psign*(enel + enph + enelph)
 self%scf_hist_spin(2, iter, ip, spin) = psign*enel
 self%scf_hist_spin(3, iter, ip, spin) = psign*enph
 self%scf_hist_spin(4, iter, ip, spin) = psign*enelph
 self%scf_hist_spin(5, iter, ip, spin) = psign*eps
 self%scf_hist_spin(6, iter, ip, spin) = polstate%gradres(ip)
 self%nstep2cv_spin(ip, spin) = iter

end subroutine varpeq_record
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_init
!! NAME
!!  varpeq_init
!!
!! FUNCTION
!!  Initialize the oboject by setting basic variables, allocate dynamic arrays.
!!
!! INPUTS
!!  gstore<gstore_t>=Electron-phonon matrix elements and related quantities.
!!  dtset<dataset_type>=All input variables for this dataset.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_init(self, gstore, dtset)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 type(gstore_t), target, intent(in) :: gstore
 type(dataset_type), intent(in) :: dtset

!Local variables-------------------------------
!scalars
 character(len=5000) :: msg
 class(gqk_t), pointer :: gqk
 class(polstate_t), pointer :: polstate
 integer :: ierr, my_is, spin, bstart, bend, my_iq
 real(dp) :: wtq
 real(dp) :: cpu, wall, gflops

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 ! Consistency check
 if (gstore%check_cplex_qkzone_gmode(2, "bz", "bz", "phonon") /= 0) then
   ABI_ERROR("The gstore object is inconsistent with varpeq. See messages above.")
 end if

 ! Scalars
 ! character
 self%pkind = dtset%varpeq_pkind
 self%aseed = dtset%varpeq_aseed
 ! logical
 self%restart = (dtset%eph_restart /= 0)
 self%interp = (dtset%varpeq_interp /= 0)
 self%g0_flag = (dtset%varpeq_avg_g /= 0)
 ! integer
 self%nstep = dtset%varpeq_nstep
 self%nstep_ort = dtset%varpeq_nstep_ort
 self%nsppol = gstore%nsppol
 self%nstates = dtset%varpeq_nstates
 self%natom3 = gstore%cryst%natom*3
 self%max_nk = maxval(gstore%glob_nk_spin)
 self%max_nq = maxval(gstore%glob_nq_spin)
 self%max_nb = maxval(gstore%brange_spin(2,:) - gstore%brange_spin(1,:)) + 1
 self%frohl_ntheta = dtset%eph_frohl_ntheta
 ! real
 self%tolgrs = dtset%varpeq_tolgrs

 ! Static arrays
 ! integer
 self%ngkpt(:) = dtset%ngkpt(:)

 ! Dynamic arrays
 ! integer
 ABI_MALLOC(self%nk_spin, (gstore%nsppol))
 ABI_MALLOC(self%nq_spin, (gstore%nsppol))
 ABI_MALLOC(self%nb_spin, (gstore%nsppol))
 ABI_MALLOC(self%brange_spin, (2, gstore%nsppol))
 self%nk_spin(:) = gstore%glob_nk_spin(:)
 self%nq_spin(:) = gstore%glob_nq_spin(:)
 self%nb_spin(:) = gstore%brange_spin(2,:) - gstore%brange_spin(1,:) + 1
 self%brange_spin(:,:) = gstore%brange_spin(:,:)

 ABI_MALLOC(self%cvflag_spin, (self%nstates, gstore%nsppol))
 ABI_MALLOC(self%nstep2cv_spin, (self%nstates, gstore%nsppol))
 ABI_MALLOC(self%scf_hist_spin, (6, self%nstep, self%nstates, gstore%nsppol))
 ! real
 ABI_MALLOC(self%kpts_spin, (3, self%max_nk, gstore%nsppol))
 ABI_MALLOC(self%qpts_spin, (3, self%max_nq, gstore%nsppol))
 ! complex
 ABI_MALLOC(self%a_spin, (self%max_nb, self%max_nk, self%nstates, gstore%nsppol))
 ABI_MALLOC(self%b_spin, (self%natom3, self%max_nq, self%nstates, gstore%nsppol))

 ! Datatypes and pointers
 self%gstore => gstore
 call gstore%cryst%copy(self%cryst)
 self%gaps = ebands_get_gaps(gstore%ebands, ierr)

 ! Initialize polaronic states for each spin
 ABI_MALLOC(self%polstate, (gstore%my_nspins))
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   polstate => self%polstate(my_is)

   ! Scalars
   ! character
   polstate%aseed = dtset%varpeq_aseed
   ! logical
   polstate%translate = (dtset%varpeq_translate /= 0)
   ! integers
   polstate%np = dtset%varpeq_nstates
   polstate%nkbz = gstore%nkbz
   polstate%nqbz = gstore%nqbz
   ! real
   polstate%e_frohl = zero

   ! Static arrays
   ! integer
   polstate%ngkpt(:) = dtset%ngkpt(:)
   ! real
   polstate%gpr_energy(:) = dtset%varpeq_gpr_energy(:)
   polstate%gpr_length(:) = dtset%varpeq_gpr_length(:)

   ! Dynamic arrays
   ! logical
   ABI_MALLOC(polstate%has_prev_grad, (dtset%varpeq_nstates))
   polstate%has_prev_grad(:) = .false.
   ! real
   ABI_MALLOC(polstate%gradres, (dtset%varpeq_nstates))
   ABI_MALLOC(polstate%enterms, (4, dtset%varpeq_nstates))

   ABI_MALLOC(polstate%eig, (gqk%nb, gstore%ebands%nkpt))
   msg = sjoin(self%gaps%errmsg_spin(spin), ". VarPEq is incompatible with &
     metals and needs CBM/VBM for electron/hole polaron calculations.")
   ABI_CHECK(self%gaps%ierr(spin) == 0, msg)

   bstart = gstore%brange_spin(1, spin)
   bend = bstart + gqk%nb - 1
   select case(dtset%varpeq_pkind)
   case ("electron")
     polstate%eig = &
       gstore%ebands%eig(bstart:bend, :, spin) - self%gaps%cb_min(spin)
   case ("hole")
     polstate%eig = &
       -(gstore%ebands%eig(bstart:bend, :, spin) - self%gaps%vb_max(spin))
   end select

   ABI_MALLOC(polstate%my_g0, (gqk%my_npert))
   polstate%my_g0(:) = zero

   ABI_MALLOC(polstate%my_qpts, (3, gqk%my_nq))
   do my_iq=1,gqk%my_nq
     call gqk%myqpt(my_iq, gstore, wtq, polstate%my_qpts(:, my_iq))
   enddo

   ! complex
   ABI_MALLOC(polstate%my_a, (gqk%nb, gqk%my_nk, dtset%varpeq_nstates))
   ABI_MALLOC(polstate%a_glob, (gqk%nb, gqk%glob_nk))
   ABI_MALLOC(polstate%my_b, (gqk%my_npert, gqk%my_nq, dtset%varpeq_nstates))
   ABI_MALLOC(polstate%my_pc, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_grad, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_prev_grad, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_pcgrad, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_prev_pcgrad, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_pcjgrad, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_prev_pcjgrad, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%pcjgrad_glob, (gqk%nb, gqk%glob_nk))

   ! Datatypes ans pointers
   polstate%gqk => gqk
   polstate%my_kpts => gqk%my_kpts(:,:)
   polstate%krank_kpts = polstate%get_krank_glob("k", gstore%ebands%kptrlatt)
   polstate%krank_qpts = polstate%get_krank_glob("q", gstore%ebands%kptrlatt)

 enddo

 ! This is needed to fill kpts_spin/qpts_spin before the calculation, as they may
 ! be required by the varpeq_load subroutine
 call xmpi_barrier(gstore%comm)
 call self%collect()

 self%is_complete = .True.

 call cwtime_report(" varpeq: init", cpu, wall, gflops)

end subroutine varpeq_init
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_calc_fravg
!! NAME
!!  varpeq_calc_fravg
!!
!! FUNCTION
!!  Calculate the avarage Fr\"ohlich long-range contribution to the polaron
!!  binding energy at Gamma using spherical integration in the spherical region
!!  arond Gamma-point.
!!
!! INPUTS
!!  avg_g0 [optional]=If .True., avarage electron-phonon matrix elements at
!!    Gamma-point. Defatult: .True.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_calc_fravg(self, avg_g0)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 logical, optional, intent(in) :: avg_g0

!Local variables-------------------------------
!scalars
 type(ifc_type), pointer :: ifc
 class(polstate_t), pointer :: polstate
 class(gqk_t), pointer :: gqk
 integer :: comm, my_rank, nproc, ierr
 integer :: ntheta, nphi, angl_size
 integer :: iang, iatom
 integer :: nu, my_is, my_pert, pert_glob
 real(dp) :: inv_qepsq, wqnu, prefactor
 real(dp) :: cpu, wall, gflops
 complex(dpc) :: cnum
!arrays
 real(dp) :: qpt_cart(3)
 real(dp), allocatable :: phfreq(:), displ_cart(:, :, :, :)
 real(dp), allocatable :: qvers_cart(:, :)
 real(dp), allocatable :: angweight(:)
 real(dp), allocatable :: e_frohl_mode(:)
 complex(dpc) :: cp3(3)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Integrate electron-phonon energy in the spherical region around Gamma
 ABI_MALLOC(phfreq, (self%natom3))
 ABI_MALLOC(e_frohl_mode, (self%natom3))
 ABI_MALLOC(displ_cart, (2, 3, self%cryst%natom, self%natom3))

 ifc => self%gstore%ifc

 ! Create the mesh for spherical integrtion
 ntheta = self%frohl_ntheta
 nphi = 2*ntheta
 call ylm_angular_mesh(ntheta, nphi, angl_size, qvers_cart, angweight)

 ! Integrate the contribution from each phonon mode
 ! This is similar to what is done in the src/78_eph/m_sigmaph.f90 module
 e_frohl_mode(:) = zero
 do iang=1,angl_size
   if (mod(iang, nproc) /= my_rank) cycle

   qpt_cart = qvers_cart(:, iang)
   inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))
   call ifc%fourq(self%cryst, qpt_cart, phfreq, displ_cart, nanaqdir="cart")

   do nu=4,self%natom3
     wqnu = phfreq(nu)

     cp3(:) = zero
     do iatom=1,self%cryst%natom
       cp3(:) = cp3(:) + matmul(ifc%zeff(:, :, iatom), &
         cmplx(displ_cart(1,:,iatom,nu), displ_cart(2,:,iatom,nu), kind=dpc))
     enddo
     cnum = dot_product(qpt_cart, cp3)
     if (abs(cnum) < tol12) cycle

     e_frohl_mode(nu) = e_frohl_mode(nu) + &
       angweight(iang)*(abs(cnum)*inv_qepsq/wqnu)**2
   enddo
 enddo
 call xmpi_sum(e_frohl_mode, comm, ierr)

 prefactor = (eight * pi / self%cryst%ucvol) * &
   (three / (four_pi * self%cryst%ucvol * self%gstore%nqbz))**third

 e_frohl_mode(:) = prefactor*e_frohl_mode(:)
 self%e_frohl = sum(e_frohl_mode)
 ! For an electron polaron, the correction has to be negative
 if (self%pkind == "electron") self%e_frohl = -self%e_frohl

 ! If corresponding flag is provided, average the matrix elements at Gamma
 if (present(avg_g0) .and. avg_g0) then
   do my_is=1,self%gstore%my_nspins
     polstate => self%polstate(my_is)
     gqk => polstate%gqk

     polstate%e_frohl = self%e_frohl
     do my_pert=1,gqk%my_npert
       pert_glob = gqk%my_pert_start + my_pert - 1
       wqnu = phfreq(pert_glob)
       polstate%my_g0(my_pert) = &
         sqrt(polstate%nqbz * wqnu * e_frohl_mode(pert_glob) / two)

     enddo
   enddo
 endif

 ABI_FREE(phfreq)
 ABI_FREE(e_frohl_mode)
 ABI_FREE(displ_cart)
 ABI_FREE(qvers_cart)
 ABI_FREE(angweight)

 call cwtime_report(" varpeq: calc_fravg", cpu, wall, gflops)

end subroutine varpeq_calc_fravg
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_free
!! NAME
!!  polstate_free
!!
!! FUNCTION
!!  Free dynamic memory
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_free(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!----------------------------------------------------------------------

 ! Free allocatable arrays

 ! logical
 ABI_SFREE(self%has_prev_grad)
 ! real
 ABI_SFREE(self%gradres)
 ABI_SFREE(self%enterms)
 ABI_SFREE(self%eig)
 ABI_SFREE(self%my_g0)
 ABI_SFREE(self%my_qpts)
 ! complex
 ABI_SFREE(self%my_a)
 ABI_SFREE(self%a_glob)
 ABI_SFREE(self%my_b)
 ABI_SFREE(self%my_pc)
 ABI_SFREE(self%my_grad)
 ABI_SFREE(self%my_prev_grad)
 ABI_SFREE(self%my_pcjgrad)
 ABI_SFREE(self%my_prev_pcjgrad)
 ABI_SFREE(self%pcjgrad_glob)

 ! Free local datatypes & nullify pointers
 self%my_kpts => null()
 self%gqk => null()

 call self%krank_kpts%free()
 call self%krank_qpts%free()

end subroutine polstate_free
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_setup
!! NAME
!!  polstate_setup
!!
!! FUNCTION
!!  Setup optimization process at a given polaronic state.
!!  This routine specfifies an initial electronic vector A_nk, either by
!!  initializaing it by a pre-determined algorithm or loading from a
!!  *VARPEQ.nc netcdf file.
!!
!! INPUTS
!!  ip=Index of the polaronic state.
!!  a_src(self%gqk%nb, self%gqk%glob_nk) [optional]=Global A_nk coefficients at
!!    this state, which have to be provided if load_src=.True.
!!  load_src [optional]=.True. if A_nk is initialized from an external source,
!!    e.g. loaded from file. Default: .False.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_setup(self, ip, a_src, load)

!Arguments ------------------------------------
 class(polstate_t), target, intent(inout) :: self
 integer, intent(in) :: ip
 logical, optional, intent(in) :: load
 complex(dp), optional, intent(in) :: a_src(self%gqk%nb, self%gqk%glob_nk)

!Local variables-------------------------------
 real(dp) :: a_sqnorm
 !complex(dp) :: prod

!----------------------------------------------------------------------

 if (present(load) .and. load) then
   ABI_CHECK(present(a_src), "polstate_setup: A_nk is expected but not provided")
   call self%load_a(a_src, ip)
 else
   call self%seed_a(self%aseed, ip)
 endif

 ! Orthogonalize current states to the previous ones
 call self%ort_to_states(self%my_a(:,:,ip), 1, ip-1, tr_flag=self%translate)

 ! Normalize A_nk at current polaronic state
 a_sqnorm = self%get_sqnorm("a", ip)
 if (a_sqnorm > tol12) then
   self%my_a(:,:,ip) = sqrt(self%nkbz/a_sqnorm) * self%my_a(:,:,ip)
 endif

end subroutine polstate_setup
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_update_a
!! NAME
!!  polstate_update_a
!!
!! FUNCTION
!!  Update the vector of electronic coefficients A_nk by line minimization in
!!  the pcj direction for current polaronic state.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_update_a(self, ip)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
 real(dp) :: theta

!----------------------------------------------------------------------

 theta = self%get_lm_theta(ip)
 self%my_a(:,:,ip) = &
   cos(theta)*self%my_a(:,:,ip) + sin(theta)*self%my_pcjgrad(:,:)

end subroutine polstate_update_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_update_pc
!! NAME
!!  polstate_update_pc
!!
!! FUNCTION
!!  Update the preconditioner for the present configuration. Changing this
!!  procedure may significantly improve (or worsen) the optimization process.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_update_pc(self, ip)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: my_ik, ik_ibz, ib
 real(dp) :: eps

!----------------------------------------------------------------------

 gqk => self%gqk

 eps = self%enterms(4, ip)
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     self%my_pc(ib, my_ik) = &
       one/abs(self%eig(ib, ik_ibz) - two*abs(self%e_frohl) + abs(eps))
   enddo
 enddo

end subroutine polstate_update_pc
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_ort_to_states
!! NAME
!!  polstate_ort_to_states
!!
!! FUNCTION
!!  Orthognongalize a vector wrt polaronic states using the Gram-Schmidt process.
!!  The orthogonalization is performed for a range of states, specified by arguments.
!!
!! INPUTS
!!  my_v(:,:)=Vetor to be orthogonalized
!!  istart=Index of starting polaronic state
!!  iend=Index of final polaronic state
!!  tr_flag=.True. if orthognoalization must include all states invariant by
!!    translations inside a supercell
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_ort_to_states(self, my_v, istart, iend, tr_flag)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 logical, intent(in) :: tr_flag
 integer, intent(in) :: istart, iend
 complex(dp), intent(inout) :: my_v(self%gqk%nb, self%gqk%my_nk)

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ip, my_ik, vx, vy, vz
 complex(dp) :: phase, proj
 integer :: tr_vec(3), ngkpt_tr(3)
 real(dp) :: kpt(3)
 complex(dp) :: a_tr(self%gqk%nb, self%gqk%my_nk)

!----------------------------------------------------------------------

 gqk => self%gqk

 ! TODO: optimize
 ngkpt_tr(:) = 1
 if (tr_flag) ngkpt_tr(:) = self%ngkpt(:)

 do ip=istart,iend

   do vx=1,ngkpt_tr(1)
     tr_vec(1) = vx - 1
     do vy=1,ngkpt_tr(2)
       tr_vec(2) = vy - 1
       do vz=1,ngkpt_tr(3)
         tr_vec(3) = vz - 1

         do my_ik=1,gqk%my_nk
           kpt(:) = self%my_kpts(:, my_ik)
           phase = exp(j_dpc*sum(kpt(:)*tr_vec(:))*two_pi)
           a_tr(:, my_ik) = phase*self%my_a(:, my_ik, ip)
         enddo

         proj = get_proj_(a_tr)
         my_v(:,:) = my_v(:,:) - proj * a_tr(:,:)

       enddo
     enddo
   enddo

 enddo

 !----------------------------------------------------------------------
 contains

 complex(dp) function get_proj_(my_u) result(proj)

  complex(dp), intent(in) :: my_u(gqk%nb, gqk%my_nk)
  integer :: ierr
  real(dp) :: u_sqnorm
 !----------------------------------------------------------------------
  u_sqnorm = sum(abs(my_u(:,:))**2)
  call xmpi_sum(u_sqnorm, gqk%kpt_comm%value, ierr)

  proj = zero
  if (u_sqnorm > tol12) then
    proj = sum(conjg(my_u(:,:))*my_v(:,:))
    call xmpi_sum(proj, gqk%kpt_comm%value, ierr)
    proj = proj/u_sqnorm
  endif
 end function get_proj_

end subroutine polstate_ort_to_states
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_lm_theta
!! NAME
!!  polstate_get_lm_theta
!!
!! FUNCTION
!!  Calculate line minimization parameter theta that minimizes polaron binding
!!  energy at the next electronic configuration:
!!  A^n = A^(n-1)*cos(theta) + D^(n-1)*sin(theta), where
!!  A^(n-1) and D^(n-1) are current electronic configuration and gradient.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!  theta=Line minimization parameter.
!!
!! SOURCE

real(dp) function polstate_get_lm_theta(self, ip) result(theta)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 logical :: q_gamma
 integer :: ierr
 integer :: my_iq, my_pert
 integer :: my_ik, ik_ibz, ik_forw, ib, jb
 real(dp) :: sqnorm
 real(dp) :: term_sin2, term_sincos, eps
 complex(dp) :: a_from, a_forw, d_from, d_forw
 complex(dp) :: g_forw, g0, b
!arrays
 real(dp) :: kpt(3), qpt(3), kpq(3)
 complex(dp) :: ak(self%gqk%nb), akq(self%gqk%nb)
 complex(dp) :: dk(self%gqk%nb), dkq(self%gqk%nb)
 complex(dp) :: bq(self%gqk%my_npert)

!----------------------------------------------------------------------

 gqk => self%gqk

 ! Orthogonalize pcj direction to the current state and normalize
 call self%ort_to_states(self%my_pcjgrad, ip, ip, tr_flag=.false.)

 sqnorm = self%get_sqnorm('pcjgrad', ip)
 self%my_pcjgrad(:,:) = sqrt(self%nkbz/sqnorm)*self%my_pcjgrad(:,:)

 ! Calculation of theta requires globally available pcj direction
 call self%gather("pcjgrad", ip)

! Scattering-dependent part
 term_sin2 = zero
 term_sincos = zero
 do my_ik=1,gqk%my_nk
   kpt(:) = self%my_kpts(:, my_ik)
   ak(:) = self%my_a(:, my_ik, ip)
   dk(:) = self%my_pcjgrad(:, my_ik)

   do my_iq=1,gqk%my_nq
     qpt(:) = self%my_qpts(:, my_iq)

     ! Find k+q-->k' index in krank_kpts
     kpq(:) = kpt(:) + qpt(:)
     ik_forw = self%krank_kpts%get_index(kpq)
     ! If erange filter was used in gstore, some transitions are not valid
     if (ik_forw == -1) cycle

     ! Check if q=\Gamma
     q_gamma = .false.
     if (all(abs(qpt) < tol6)) q_gamma = .true.

     bq(:) = self%my_b(:, my_iq, ip)
     akq(:) = self%a_glob(:, ik_forw)
     dkq(:) = self%pcjgrad_glob(:, ik_forw)

     do ib=1,gqk%nb
       a_from = ak(ib)
       d_from = dk(ib)

       do jb=1,gqk%nb
         a_forw = akq(jb)
         d_forw = dkq(jb)

         do my_pert=1,gqk%my_npert
           b = bq(my_pert)

           g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
           ! Add long-range correction to matrix elements at Gamma
           g0 = self%my_g0(my_pert)
           if (q_gamma .and. (ib == jb)) then
             g_forw = g_forw + g0
           endif

           term_sin2 = term_sin2 + real(d_from*conjg(b)*g_forw*conjg(d_forw))
           term_sincos = term_sincos + &
             real((a_from*conjg(d_forw) + d_from*conjg(a_forw))*conjg(b)*g_forw)
         enddo
       enddo
     enddo
   enddo
 enddo ! Scattering-dependent part
 call xmpi_sum(term_sin2, gqk%qpt_pert_comm%value, ierr)
 call xmpi_sum(term_sincos, gqk%qpt_pert_comm%value, ierr)
 term_sin2 = -two*term_sin2/self%nqbz
 term_sincos = -two*term_sincos/self%nqbz

 ! Scattering-independent part
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     a_from = self%my_a(ib, my_ik, ip)
     d_from = self%my_pcjgrad(ib, my_ik)

     term_sin2 = term_sin2 + self%eig(ib, ik_ibz)*abs(d_from)**2
     term_sincos = term_sincos + &
       self%eig(ib, ik_ibz)*real(d_from*conjg(a_from) + conjg(d_from)*a_from)
   enddo
 enddo ! Scattering-independent part
 call xmpi_sum(term_sin2, gqk%kpt_comm%value, ierr)
 call xmpi_sum(term_sincos, gqk%kpt_comm%value, ierr)
 term_sin2 = term_sin2/self%nkbz
 term_sincos = term_sincos/self%nkbz

 ! Line-minimization theta
 eps = self%enterms(4, ip)
 theta = half*atan2(-term_sincos, term_sin2 - eps)

end function polstate_get_lm_theta
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_calc_pcjgrad
!! NAME
!!  polstate_calc_pcjgrad
!!
!! FUNCTION
!!  Calculate preconditioned conjugate gradient direction, orthogonal to all
!!  already optimized polaronic states. The procedure is similar to the one
!!  described in [Payne et al, Rev. Mod. Phys, 64, 4, 1045-1097 (1992)].
!!
!! INPUTS
!!   ip=Index of a polaronic state.
!!   ii=Iteration number.
!!   nstep_ort=Iteration number, after which the orthogonality constraint
!!     on all PREVIOUS states is lifted.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_calc_pcjgrad(self, ip, ii, nstep_ort)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip
 integer, intent(in) :: ii
 integer, intent(in) :: nstep_ort

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr
 !real(dp) :: eps
 complex(dp) :: beta, beta_num, beta_den

!----------------------------------------------------------------------

 gqk => self%gqk

 ! Orthogonalize current gradient to all previous bands
 if (ii <= nstep_ort) then
   call self%ort_to_states(self%my_grad, 1, ip-1, tr_flag=self%translate)
 endif

 ! Precondtion vector
 self%my_pcgrad(:,:) = self%my_pc(:,:)*self%my_grad(:,:)
 ! Orthogonalize to all bands
 if (ii <= nstep_ort) then
   call self%ort_to_states(self%my_pcgrad, 1, ip-1, tr_flag=self%translate)
 endif
 call self%ort_to_states(self%my_pcgrad, ip, ip, tr_flag=.false.)

 ! Conjugate gradient direction
 if (self%has_prev_grad(ip)) then
   ! Polak-Ribiere coefficient
   beta_num = &
     sum(conjg(self%my_pcgrad(:,:))*(self%my_grad(:,:) - self%my_prev_grad(:,:)))
   beta_den = sum(conjg(self%my_prev_pcgrad(:,:))*self%my_prev_grad(:,:))
   call xmpi_sum(beta_num, gqk%kpt_comm%value, ierr)
   call xmpi_sum(beta_den, gqk%kpt_comm%value, ierr)
   beta = beta_num / beta_den
   if (abs(aimag(beta)) < tol6) beta = real(beta)

   self%my_pcjgrad(:,:) = self%my_pcgrad(:,:) + beta*self%my_prev_pcjgrad(:,:)
 else
   self%my_pcjgrad(:,:) = self%my_pcgrad(:,:)
 endif

 ! Save previous gradients
 self%my_prev_grad(:,:) = self%my_grad(:,:)
 self%my_prev_pcgrad(:,:) = self%my_pcgrad(:,:)
 self%my_prev_pcjgrad(:,:) = self%my_pcjgrad(:,:)
 self%has_prev_grad(ip) = .true.

end subroutine polstate_calc_pcjgrad
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_calc_grad
!! NAME
!!  polstate_calc_grad
!!
!! FUNCTION
!!  Calculate the electronic gradient D_nk at the current configuration.
!!
!! INPUTS
!!   ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_calc_grad(self, ip)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 logical :: q_gamma
 integer :: ierr
 integer :: my_iq, my_pert
 integer :: my_ik, ik_ibz, ik_forw, ik_back, ib, jb
 real(dp) :: eps
 complex(dp) :: a_forw, a_back
 complex(dp) :: g_forw, g_back
 complex(dp) :: b, g0
!arrays
 real(dp) :: kpt(3), qpt(3), kpq(3), kmq(3)
 complex(dp) :: akq(self%gqk%nb), akmq(self%gqk%nb) ! ak(self%gqk%nb),
 complex(dp) :: bq(self%gqk%my_npert)
 complex(dp), allocatable :: gq_gathered(:,:,:,:)

!----------------------------------------------------------------------

 gqk => self%gqk

 ABI_MALLOC(gq_gathered, (gqk%my_npert, gqk%nb, gqk%nb, gqk%glob_nk))

 ! Scattering-dependent part
 self%my_grad(:, :) = zero
 do my_iq=1,gqk%my_nq
   qpt(:) = self%my_qpts(:, my_iq)
   bq(:) = self%my_b(:, my_iq, ip)

   ! Check if q=\Gamma
   q_gamma = .false.
   if (all(abs(qpt) < tol6)) q_gamma = .true.

   ! For this q, gather all the k-distributed matrix elements
   call gqk%gather("q", my_iq, gq_gathered)

   do my_ik=1,gqk%my_nk
     kpt(:) = self%my_kpts(:, my_ik)

     ! Forward scattering
     ! Find k+q-->k' index in krank_kpts
     kpq(:) = kpt(:) + qpt(:)
     ik_forw = self%krank_kpts%get_index(kpq)

     ! If erange filter was used in gstore, some transitions are not valid
     if (ik_forw /= -1) then
       akq(:) = self%a_glob(:, ik_forw)
       do ib=1,gqk%nb
         do jb=1,gqk%nb
           a_forw = akq(jb)

           do my_pert=1,gqk%my_npert
             b = bq(my_pert)

             g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
             ! Add long-range correction to matrix elements at Gamma
             g0 = self%my_g0(my_pert)
             if (q_gamma .and. (ib == jb)) then
               g_forw = g_forw + g0
             endif

             self%my_grad(ib, my_ik) = &
               self%my_grad(ib, my_ik) + a_forw*b*conjg(g_forw)
           enddo
         enddo
       enddo
     endif ! Forward scattering

     ! Backward scattering
     ! Find k-q-->k' index in krank_kpts
     kmq(:) = kpt(:) - qpt(:)
     ik_back = self%krank_kpts%get_index(kmq)

     ! If erange filter was used in gstore, some transitions are not valid
     if (ik_back /= -1) then
       akmq(:) = self%a_glob(:, ik_back)
       do ib=1,gqk%nb
         do jb=1,gqk%nb
           a_back = akmq(jb)

           do my_pert=1,gqk%my_npert
             b = bq(my_pert)

             g_back = gq_gathered(my_pert, ib, jb, ik_back)
             ! Add long-range correction to matrix elements at Gamma
             g0 = self%my_g0(my_pert)
             if (q_gamma .and. (ib == jb)) then
               g_back = g_back + g0
             endif

             self%my_grad(ib, my_ik) = &
               self%my_grad(ib, my_ik) + a_back*conjg(b)*g_back
           enddo
         enddo
       enddo
     endif ! Backward scattering

   enddo
 enddo
 call xmpi_sum(self%my_grad, gqk%qpt_pert_comm%value, ierr)
 self%my_grad(:, :) = -two/(self%nkbz*self%nqbz) * self%my_grad(:, :)

 ! Scattering-independent part
 eps = self%enterms(4, ip)
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     self%my_grad(ib, my_ik) = self%my_grad(ib, my_ik) + &
       two/self%nkbz * (self%eig(ib, ik_ibz) - eps) * self%my_a(ib, my_ik, ip)
   enddo
 enddo

 ABI_FREE(gq_gathered)

end subroutine polstate_calc_grad
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_localize
!! NAME
!!  polstate_localize
!!
!! FUNCTION
!!  Calculate vibrational coefficients, polaron energy terms and polaron energy
!!  level at a current state.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_localize(self, ip)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!----------------------------------------------------------------------

 ! Calculation of B_qnu requires globally available A_nk
 call self%gather("a", ip)
 call self%calc_b_from_a(ip)

 ! Electronic term
 self%enterms(1, ip) = self%get_enel(ip)
 ! Vibrational term
 self%enterms(2, ip) = self%get_enph(ip)
 ! Electron-phonon term
 self%enterms(3, ip) = self%get_enelph(ip)
 ! Polaron energy level
 self%enterms(4, ip) = self%enterms(1, ip) + self%enterms(3, ip)

end subroutine polstate_localize
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_enelph
!! NAME
!!  polstate_get_enelph
!!
!! FUNCTION
!!  Returns the electron-phonon term of the polaron bidning energy at a specified
!!  polaronic state.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!  enph=Electron-phonon term of the polaron binding energy.
!!
!! SOURCE

real(dp) function polstate_get_enelph(self, ip) result(enelph)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 logical :: q_gamma
 integer :: ierr, my_iq, my_pert, my_ik, ik_forw, ib, jb
 complex(dp) :: a_from, a_forw, g_forw, g0, b
 complex(dp) :: ak(self%gqk%nb), akq(self%gqk%nb), bq(self%gqk%my_npert)
!arrays
 real(dp) :: kpt(3), qpt(3), kpq(3)

!----------------------------------------------------------------------

 gqk => self%gqk

 enelph = zero
 do my_ik=1,gqk%my_nk
   kpt(:) = self%my_kpts(:, my_ik)
   ak(:) = self%my_a(:, my_ik, ip)

   do my_iq=1,gqk%my_nq
     qpt(:) = self%my_qpts(:, my_iq)

     ! Find k+q-->k' index in krank_kpts
     kpq(:) = kpt(:) + qpt(:)
     ik_forw = self%krank_kpts%get_index(kpq)
     ! If erange filter was used in gstore, some transitions are not valid
     if (ik_forw == -1) cycle

     ! Check if q=\Gamma
     q_gamma = .false.
     if (all(abs(qpt) < tol6)) q_gamma = .true.

     akq(:) = self%a_glob(:, ik_forw)
     bq(:) = self%my_b(:, my_iq, ip)

     do ib=1,gqk%nb
       a_from = ak(ib)
       do jb=1,gqk%nb
         a_forw = akq(jb)
         do my_pert=1,gqk%my_npert
           b = bq(my_pert)

           g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
           ! Add long-range correction to matrix elements at Gamma
           g0 = self%my_g0(my_pert)
           if (q_gamma .and. (ib == jb)) then
             g_forw = g_forw + g0
           endif

           enelph = enelph + real(a_from*conjg(b)*g_forw*conjg(a_forw))
         enddo
       enddo
     enddo

   enddo
 enddo
 call xmpi_sum(enelph, gqk%comm%value, ierr)
 enelph = -two*enelph/(self%nkbz*self%nqbz)

 end function polstate_get_enelph
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_enph
!! NAME
!!  polstate_get_enph
!!
!! FUNCTION
!!  Returns the vibrational term of the polaron bidning energy at a specified
!!  polaronic state.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!  enph=Vibrational term of the polaron binding energy.
!!
!! SOURCE

real(dp) function polstate_get_enph(self, ip) result(enph)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr, my_iq, my_pert
!----------------------------------------------------------------------

 gqk => self%gqk

 enph = zero
 do my_iq=1,gqk%my_nq
   do my_pert=1,gqk%my_npert
     enph = enph + gqk%my_wnuq(my_pert, my_iq)*abs(self%my_b(my_pert, my_iq, ip))**2
   enddo
 enddo
 call xmpi_sum(enph, gqk%qpt_pert_comm%value, ierr)
 enph = enph/self%nqbz

end function polstate_get_enph
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_enel
!! NAME
!!  polstate_get_enel
!!
!! FUNCTION
!!  Returns the electronic term of the polaron bidning energy at a specified
!!  polaronic state.
!!
!! INPUTS
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!  enel=Electronic term of the polaron binding energy.
!!
!! SOURCE

real(dp) function polstate_get_enel(self, ip) result(enel)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr, my_ik, ik_ibz, ib

!----------------------------------------------------------------------

 gqk => self%gqk

 enel = zero
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     enel = enel + self%eig(ib, ik_ibz)*abs(self%my_a(ib, my_ik, ip))**2
   enddo
 enddo
 call xmpi_sum(enel, gqk%kpt_comm%value, ierr)
 enel = enel/self%nkbz

end function polstate_get_enel
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_calc_b_from_a
!! NAME
!!  polstate_calc_b_from_a
!!
!! FUNCTION
!!  Calculate vibrational coefficients B_q\nu at the current configuration,
!!  so the energy gradient wrt B_q\nu is 0.
!!
!! INPUTS
!!   ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_calc_b_from_a(self, ip)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 logical :: q_gamma
 integer :: ierr, my_iq, my_pert, my_ik, ik_forw, ib, jb
 real(dp) :: wqnu
 complex(dp) :: a_from, a_forw, g_forw, g0, b_tmp
!arrays
 real(dp) :: qpt(3), kpq(3)
 complex(dp) :: ak(self%gqk%nb), akq(self%gqk%nb)

!----------------------------------------------------------------------

 gqk => self%gqk

 do my_iq=1,gqk%my_nq
   qpt(:) = self%my_qpts(:, my_iq)
   ! Check if q=\Gamma
   q_gamma = .false.
   if (all(abs(qpt) < tol6)) q_gamma = .true.

   do my_pert=1,gqk%my_npert
     wqnu = gqk%my_wnuq(my_pert, my_iq)
     ! Skip acoustic modes at Gamma
     if (wqnu == zero) then
       self%my_b(my_pert, my_iq, ip) = zero
       cycle
     endif
     g0 = self%my_g0(my_pert)

     ! For this q and perurbation, calculate B_q\nu sum
     b_tmp = zero
     do my_ik=1,self%gqk%my_nk
       ! Find k+q-->k' index in krank_kpts
       kpq(:) = qpt(:) + self%my_kpts(:, my_ik)
       ik_forw = self%krank_kpts%get_index(kpq)
       ! If erange filter was used in gstore, some transitions are not valid
       if (ik_forw == -1) cycle

       ak(:) = self%my_a(:, my_ik, ip)
       akq(:) = self%a_glob(:, ik_forw)

       do ib=1,gqk%nb
         a_from = ak(ib)
         do jb=1,gqk%nb
           a_forw = akq(jb)

           g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
           ! Add long-range correction to matrix elements at Gamma
           if (q_gamma .and. (ib == jb)) then
             g_forw = g_forw + g0
           endif

           b_tmp = b_tmp + a_from*g_forw*conjg(a_forw)
         enddo
       enddo
     enddo
     call xmpi_sum(b_tmp, gqk%kpt_comm%value, ierr)

     self%my_b(my_pert, my_iq, ip) = b_tmp/(self%nkbz * wqnu)
   enddo
 enddo

end subroutine polstate_calc_b_from_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_load_a
!! NAME
!!  polstate_load_a
!!
!! FUNCTION
!!  Load the initial vector of electronic coefficients A_nk from source.
!!
!! INPUTS
!!  a_src(self%gqk%nb, self%gqk%glob_nk)=Global A_nk array to be loaded.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_load_a(self, a_src, ip)

!Arguments ------------------------------------
!scalars
 class(polstate_t), intent(inout) :: self
 integer, intent(in) :: ip
!arrays
 complex(dp), intent(in) :: a_src(self%gqk%nb, self%gqk%glob_nk)

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: my_ik, ik_glob, ib
 !complex(dp) :: ank

!----------------------------------------------------------------------

 gqk => self%gqk

 do my_ik=1,gqk%my_nk
   ik_glob = gqk%my_kstart + my_ik - 1
   do ib=1,gqk%nb
     self%my_a(ib, my_ik, ip) = a_src(ib, ik_glob)
   enddo
 enddo

end subroutine polstate_load_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_seed_a
!! NAME
!!  polstate_seed_a
!!
!! FUNCTION
!!  Specify initial vector of electronic coefficients A_nk at runtime.
!!
!! INPUTS
!!  mode=Select the initialization type. Possible options:
!!    "gau_energy" ---> Gaussian shape based on the energy of electronic states;
!!    "gau_length" ---> Gaussian shape based on the polaron localization length;
!!    "even" ---> equal contribution from each electronic state;
!!    "random" ---> random initalization.
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_seed_a(self, mode, ip)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 character(len=*), intent(in) :: mode
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr

!----------------------------------------------------------------------

 gqk => self%gqk

 select case(mode)
 case ("gau_energy")
   call gau_energy_()
 case ("gau_length")
   call gau_length_()
 case ("random")
   call random_()
 case ("even")
   self%my_a(:,:,ip) = one + j_dpc*one
 case default
   ABI_ERROR(sjoin("polstate_seed_a, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 subroutine gau_energy_()
  real(dp) :: mu, sigma
  integer :: my_ik, ik_ibz, ib
  real(dp) :: eig
  mu = self%gpr_energy(1); sigma = self%gpr_energy(2)
  ABI_CHECK(sigma /= 0, "gpr_energy: standard deviation must be non-zero")
  do my_ik=1,gqk%my_nk
    ik_ibz = gqk%my_k2ibz(1, my_ik)
    do ib=1,gqk%nb
      eig = self%eig(ib, ik_ibz)
      self%my_a(ib, my_ik, ip) = exp(-half*(eig - mu)**2/sigma**2)
    enddo
  enddo
 end subroutine gau_energy_

 subroutine gau_length_()
  integer :: my_ik, ib
  real(dp) :: kpt(3)
  do my_ik=1,gqk%my_nk
    kpt(:) = self%my_kpts(:, my_ik)
    do ib=1,gqk%nb
      self%my_a(ib, my_ik, ip) = exp(-sum((kpt(:)*self%gpr_length(:))**2))
    enddo
  enddo
 end subroutine gau_length_

 subroutine random_()
   real(dp), allocatable :: re_rand(:,:), im_rand(:,:)
   ABI_MALLOC(re_rand, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(im_rand, (gqk%nb, gqk%my_nk))
   call random_number(re_rand)
   call random_number(im_rand)
   self%my_a(:,:,ip) = re_rand(:,:) + j_dpc*im_rand(:,:)
   call xmpi_sum(self%my_a(:,:,ip), gqk%qpt_pert_comm%value, ierr)
 end subroutine random_

end subroutine polstate_seed_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_sqnorm
!! NAME
!!  polstate_get_sqnorm
!!
!! FUNCTION
!!  Helper function that calculates squared L^2-norm of an MPI-distributed
!!  array at current state.
!!
!! INPUTS
!!  mode=Select which array to gather. Possible options:
!!    "a" ---> A_nk coefficients for this state;
!!    "b" ---> B_q\nu coefficients for this state;
!!    "grad" ---> current gradient;
!!    "pcjgrad" ---> current preconditioned conjugate gradient direction.
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!   Squared L^2-norm of the specified array.
!!
!! SOURCE

real(dp) function polstate_get_sqnorm(self, mode, ip) result(sqnorm)

!Arguments ------------------------------------
 class(polstate_t), target, intent(inout) :: self
 character(len=*), intent(in) :: mode
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk
 select case(mode)
 case ("a")
   sqnorm = get_sqnorm_(self%my_a(:,:,ip), gqk%kpt_comm%value)
 case ("b")
   sqnorm = get_sqnorm_(self%my_b(:,:,ip), gqk%qpt_comm%value)
 case ("grad")
   sqnorm = get_sqnorm_(self%my_grad, gqk%kpt_comm%value)
 case ("pcjgrad")
   sqnorm = get_sqnorm_(self%my_pcjgrad, gqk%kpt_comm%value)
 case default
   ABI_ERROR(sjoin("polstate_get_sqnorm, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 real(dp) function get_sqnorm_(my_arr, comm) result(sqnorm)
  integer, intent(in) :: comm
  complex(dp), intent(in) :: my_arr(:, :)
  integer :: ierr
  sqnorm = sum(abs(my_arr(:,:))**2)
  call xmpi_sum(sqnorm, comm, ierr)
 end function get_sqnorm_

end function polstate_get_sqnorm
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_gather
!! NAME
!!  polstate_gather
!!
!! FUNCTION
!!  Helper function that gathers a MPI-distributed array from the datatype
!!  into a global one.
!!
!! INPUTS
!!  mode=Select which array to gather. Possible options:
!!    "a" ---> my_a(:,:,ip) -- array of A_nk coefficients for this state;
!!    "pcjgrad" ---> my_pcjgrad(:,:) -- gradient.
!!  ip=Index of a polaronic state.
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_gather(self, mode, ip)

!Arguments ------------------------------------
 class(polstate_t), target, intent(inout) :: self
 character(len=*),intent(in) :: mode
 integer, intent(in) :: ip

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk
 select case(mode)
 case ("a")
   call gather_(self%my_a(:,:,ip), gqk%my_nk, gqk%my_kstart, self%a_glob, &
     gqk%kpt_comm%value)
 case ("pcjgrad")
   call gather_(self%my_pcjgrad, gqk%my_nk, gqk%my_kstart, self%pcjgrad_glob, &
     gqk%kpt_comm%value)
 case default
   ABI_ERROR(sjoin("polstate_gather, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 subroutine gather_(my_arr, my_nk, my_kstart, glob_arr, comm)

  integer, intent(in) :: comm, my_nk, my_kstart
  complex(dp), intent(in) :: my_arr(:, :)
  complex(dp), intent(out) :: glob_arr(:, :)

  integer :: ierr, my_ik, ik_glob

 !----------------------------------------------------------------------

  glob_arr(:, :) = zero
  do my_ik=1,my_nk
    ik_glob = my_ik + my_kstart - 1
    glob_arr(:, ik_glob) = my_arr(:, my_ik)
  enddo
  call xmpi_sum(glob_arr, comm, ierr)

 end subroutine gather_

end subroutine polstate_gather
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_krank_glob
!! NAME
!!  polstate_get_krank_glob
!!
!! FUNCTION
!!  Helper function that computes krank objects for the full BZ q/k-meshes.
!!
!! INPUTS
!!  mode=Select reciprocal space for which krank is computed. Possible options:
!!    "k" ---> k-space;
!!    "q" ---> q-space.
!!  kptrlatt(3,3)=Lattice specifiyng the grid.
!!
!! OUTPUT
!!  krank<krank_t>=Object proiding mapping between reciprocal space points and
!!    their indices.
!!
!! SOURCE

type(krank_t) function polstate_get_krank_glob(self, mode, kptrlatt) &
    result(krank_kpts)

!Arguments ------------------------------------
!scalars
 class(polstate_t), target, intent(in) :: self
 character(len=*), intent(in) :: mode
!arrays
 integer, intent(in) :: kptrlatt(3, 3)

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk

 select case(mode)
 case ("k")
   krank_kpts = get_krank_glob_(self%my_kpts, gqk%my_nk, gqk%my_kstart, &
     gqk%glob_nk, gqk%kpt_comm%value)
 case ("q")
   krank_kpts = get_krank_glob_(self%my_qpts, gqk%my_nq, gqk%my_qstart, &
     gqk%glob_nq, gqk%qpt_comm%value)
 case default
   ABI_ERROR(sjoin("polstate_get_krank_glob, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 type(krank_t) function get_krank_glob_(my_kpts, my_nk, my_kstart, glob_nk, &
     comm) result(krank_kpts)

  integer, intent(in) :: comm, my_nk, my_kstart, glob_nk
  real(dp), intent(in) :: my_kpts(3, my_nk)

  type(krank_t) :: krank_tmp
  integer :: ierr, my_ik, ik_glob
  real(dp) :: kpts(3, glob_nk)

 !----------------------------------------------------------------------

  kpts(:, :) = zero
  do my_ik=1,my_nk
    ik_glob = my_ik + my_kstart - 1
    kpts(:, ik_glob) = my_kpts(:, my_ik)
  enddo
  call xmpi_sum(kpts, comm, ierr)

  krank_tmp = krank_from_kptrlatt(glob_nk, kpts, kptrlatt, &
    compute_invrank=.True.)
  krank_kpts = krank_tmp%copy()
  call krank_tmp%free()

 end function get_krank_glob_

end function polstate_get_krank_glob
!!***

!!****f* m_varpeq/varpeq_plot
!! NAME
!!  varpeq_plot
!!
!! FUNCTION
!! Compute polaron wavefunctions and atomic displacements in the supercell and write results to XSF files
!!
!! INPUTS
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_plot(wfk0_path, ngfft, dtset, dtfil, cryst, ebands, pawtab, psps, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: ngfft(18) !,ngfftf(18)
 type(dataset_type), intent(in) :: dtset
 type(datafiles_type), intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, ndat1 = 1
 integer :: my_rank, ib, irsp, ir, uc_idx, mpw, band, bstart ! nfft, nfftf, mgfft, mgfftf, n1, n2, n3, n4, n5, n6,
 integer :: natom, natom3, nsppol, nspinor, nspden, nkibz, mband, spin, ik, sc_nfft, ebands_timrev, cnt, nproc, num_writes, ii
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_kq_ibz, istwf_k_ibz, nkpg_k, ierr, nk, spinor, spad, nkbz, ncid
 integer :: nqibz, iq, qtimrev, iq_ibz, isym_q, trev_q, qptopt, uc_iat, sc_iat, nu, nqbz, ip ! icell,
 logical :: isirr_k, isirr_q, have_scell_q, use_displaced_scell
 real(dp) :: cpu_all, wall_all, gflops_all, spread
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(varpeq_t) :: vpq
 type(wfd_t) :: wfd
 type(supercell_type), target :: scell_q, scell_k
 type(krank_t) :: krank_ibz, qrank_ibz
 complex(dp) :: a_nk, bstar_qnu, cphase, cphase_tr, c3tmp(3)
!arrays
 integer :: sc_ngfft(18), mapl_k(6), kptrlatt_(3,3), qptrlatt_(3,3)
 integer :: units(2), work_ngfft(18), gmax(3), g0_k(3), mapl_qq(6), g0_q(3), ngqpt(3)
 integer,allocatable :: nband(:,:), wfd_istwfk(:), kg_k(:,:), sc2uc(:), gbound_k(:,:)
 real(dp),parameter :: origin0(3) = zero
 real(dp) :: kk(3), kk_ibz(3), kk_sc(3), qq(3), qq_ibz(3), center_cart(3)
 real(dp),allocatable :: scred(:,:), kpg_k(:,:), ug_k(:,:), work(:,:,:,:), pol_rho(:), qibz(:,:)
 real(dp),allocatable :: pheigvec_qibz(:,:,:,:)
 real(dp),allocatable :: displ_cart_qbz(:,:,:,:), pheigvec_qbz(:,:,:,:)  !displ_red_qbz(:,:,:,:), displ_cart_qibz(:,:,:,:),
 real(dp),allocatable :: phfreqs_ibz(:,:), pheigvec_cart_ibz(:,:,:,:,:) !, pheigvec_cart_qbz(:,:,:,:)
 real(dp),allocatable :: sc_displ_cart_re(:,:,:,:), sc_displ_cart_im(:,:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: xcart_ptr(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 complex(gwpc),allocatable :: ur_k(:), pol_wf(:,:,:), ceikr(:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Read A_nk and B_qnu and other useful tables from file
 call vpq%ncread(dtfil%filvarpeqin, comm, keep_open=.False.)
 !call wrtout(std_out, " Reading done")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkibz = ebands%nkpt; mband = ebands%mband

 if (dtfil%filgstorein == ABI_NOFILE) then
   call wrtout(units, "gstore_filepath is not specified in input. Cannot compute polaron-induced displacements!")
   have_scell_q = .False.

 else
   ! Start by reading ph displacements and frequencies in the IBZ from the gstore file.
   ! First compute displaced supercell then polaron wf so that we can use both when writing the XSF file.
   call wrtout(units, sjoin(" Computing polaron-induced displacements. Reading phonons from: ", dtfil%filgstorein))
   call cwtime(cpu_all, wall_all, gflops_all, "start")
   have_scell_q = .True.

   NCF_CHECK(nctk_open_read(ncid, dtfil%filgstorein, comm))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_nqibz", nqibz))
   !NCF_CHECK(nctk_get_dim(ncid, "gstore_nqbz", nqbz))

   ! TODO: Wrap phstore API?
   ! Encaspulate this part as we're gonna re-use it to deal with hopping
   !call gstore_read_ph_qibz(dtfil%filgstorein, ph, comm)
   !call ph%free()

   ABI_MALLOC(qibz, (3, nqibz))
   ABI_MALLOC(phfreqs_ibz, (natom3, nqibz))
   ABI_MALLOC(pheigvec_cart_ibz, (2, 3, cryst%natom, cryst%natom * 3, nqibz))
   if (nproc > 1) then
     NCF_CHECK(nctk_set_collective(ncid, vid("gstore_qibz")))
     NCF_CHECK(nctk_set_collective(ncid, vid("phfreqs_ibz")))
     NCF_CHECK(nctk_set_collective(ncid, vid("pheigvec_cart_ibz")))
   end if
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qibz"), qibz))
   NCF_CHECK(nf90_get_var(ncid, vid("phfreqs_ibz"), phfreqs_ibz))
   NCF_CHECK(nf90_get_var(ncid, vid("pheigvec_cart_ibz"), pheigvec_cart_ibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_ngqpt"), ngqpt))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qptopt"), qptopt))
   NCF_CHECK(nf90_close(ncid))

   ABI_MALLOC(pheigvec_qibz, (2, 3, cryst%natom, natom3))
   ABI_MALLOC(displ_cart_qbz, (2, 3, cryst%natom, cryst%natom * 3))
   ABI_MALLOC(pheigvec_qbz, (2, 3, cryst%natom, 3*cryst%natom))

   qtimrev = kpts_timrev_from_kptopt(qptopt)
   nqbz = product(ngqpt)
   call kptrlatt_from_ngkpt(ngqpt, qptrlatt_)
   qrank_ibz = krank_from_kptrlatt(nqibz, qibz, qptrlatt_, compute_invrank=.False.)

   call scell_q%init(cryst%natom, qptrlatt_, cryst%rprimd, cryst%typat, cryst%xcart, cryst%znucl, xyz_order="xyz")

   ABI_CALLOC(sc_displ_cart_re, (3, scell_q%natom, vpq%nstates, nsppol))
   ABI_CALLOC(sc_displ_cart_im, (3, scell_q%natom, vpq%nstates, nsppol))

   cnt = 0
   do spin=1,nsppol
     do iq=1,vpq%nq_spin(spin)
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism inside comm.
       qq = vpq%qpts_spin(:, iq, spin)
       ! Note symrec here
       if (kpts_map("symrec", qtimrev, cryst, qrank_ibz, 1, qq, mapl_qq) /= 0) then
         ABI_ERROR("Cannot map qBZ to IBZ!")
       end if
       iq_ibz = mapl_qq(1); isym_q = mapl_qq(2)
       trev_q = mapl_qq(6); g0_q = mapl_qq(3:5)
       ! Don't test if umklapp == 0 because we use the periodic gauge:
       !
       !      phfreq(q+G) = phfreq(q) and eigvec(q) = eigvec(q+G)
       !
       isirr_q = (isym_q == 1 .and. trev_q == 0)
       qq_ibz = qibz(:, iq_ibz)
       !if (all(abs(qq_ibz) < tol6)) cycle
       pheigvec_qibz = pheigvec_cart_ibz(:,:,:,:,iq_ibz)

       if (isirr_q) then
         ! Compute phonon displacements in Cartesian coordinates
         call phdispl_from_eigvec(cryst%natom, cryst%ntypat, cryst%typat, cryst%amu, pheigvec_qibz, displ_cart_qbz)

       else
         ! Rotate phonon eigenvectors from q_ibz to q_bz.
         ! This part is needed to enforce the gauge in the ph eigenvectors, including e(-q) = e(q)^*
         call pheigvec_rotate(cryst, qq_ibz, isym_q, trev_q, pheigvec_qibz, pheigvec_qbz, displ_cart_qbz)
       end if

       ! Phase due to the primitive translation (default: 0)
       cphase_tr = exp(+j_dpc * two_pi * dot_product(qq, dtset%varpeq_trvec))

       do sc_iat=1, scell_q%natom
         uc_iat = scell_q%atom_indexing(sc_iat)
         ! Compute phase e^{iq.R}
         cphase = exp(+j_dpc * two_pi * dot_product(qq, scell_q%uc_indexing(:, sc_iat)))
         do ip=1,vpq%nstates
           ! Summing over ph modes.
           do nu=1,natom3
             bstar_qnu = vpq%b_spin(nu, iq, ip, spin) * cphase_tr
             c3tmp = (displ_cart_qbz(1,:,uc_iat,nu) + j_dpc * displ_cart_qbz(2,:,uc_iat,nu)) * bstar_qnu * cphase
             sc_displ_cart_re(:,sc_iat,ip,spin) = sc_displ_cart_re(:,sc_iat,ip,spin) + real(c3tmp)
             sc_displ_cart_im(:,sc_iat,ip,spin) = sc_displ_cart_im(:,sc_iat,ip,spin) + aimag(c3tmp) ! This to check if the imag part is zero.
           end do
         end do ! ip
       end do ! sc_iat

     end do ! iq
   end do ! spin

   ABI_FREE(qibz)
   ABI_FREE(phfreqs_ibz)
   ABI_FREE(pheigvec_cart_ibz)
   ABI_FREE(displ_cart_qbz)
   ABI_FREE(pheigvec_qbz)
   ABI_FREE(pheigvec_qibz)
   call qrank_ibz%free()

   call xmpi_sum_master(sc_displ_cart_re, master, comm, ierr)
   call xmpi_sum_master(sc_displ_cart_im, master, comm, ierr)
   sc_displ_cart_re = -sqrt2 * sc_displ_cart_re / nqbz
   sc_displ_cart_im = -sqrt2 * sc_displ_cart_im / nqbz

   ! Write polaron-induced displacements in XSF format.
   if (my_rank == master) then
     ! Handle spin-polarized case by writing two XSF files.
     do spin=1,nsppol
       ! Handle muttiple polaronic states for each spin.
       do ip=1,vpq%nstates
         write(msg, "(2(a,i0),a,es16.6)") &
           " For spin: ", spin, ": pstate: ", ip, ": maxval(abs(sc_displ_cart_re)): ", maxval(abs(sc_displ_cart_re(:,:,ip,spin)))
         call wrtout(units, msg)
         write(msg, "(2(a,i0),a,es16.6)") &
           " For spin: ", spin, ": pstate: ", ip, ": maxval(abs(sc_displ_cart_im)): ", maxval(abs(sc_displ_cart_im(:,:,ip,spin)))
         call wrtout(units, msg)

         ! Here we displace the atoms in the supercell for this spin (only master has the correct values)
         scell_q%xcart = scell_q%xcart_ref + sc_displ_cart_re(:,:,ip,spin)

         path = strcat(dtfil%filnam_ds(4), "_pstate_", itoa(ip), "_POLARON_DISPL_VECTORS.xsf")
         if (nsppol == 2) path = strcat(dtfil%filnam_ds(4), strcat("_spin_", itoa(spin)), "_pstate_", itoa(ip), "_POLARON_DISPL_VECTORS.xsf")
         call wrtout(units, sjoin(" Writing displacement vectors to:", path))
         call scell_q%write_xsf(path)
       end do
     end do
   end if

   ABI_FREE(sc_displ_cart_im)
   call cwtime_report(" Computation of polaron-induced displacements completed:", cpu_all, wall_all, gflops_all, &
                      pre_str=ch10, end_str=ch10)
 end if

 call wrtout(std_out, " varpeq_plot: computing polaron wavefunction in real space.", pre_newlines=1)

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 krank_ibz = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

 ! Initialize the wave function descriptor.
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 ! Here we use brange_spin to select the number of bands that should be read and stored in memory.
 ! For the time being, spin and k-points are not MPI-distributed inside comm.
 do spin=1,nsppol
   bstart = vpq%brange_spin(1, spin)
   do ik=1, vpq%nk_spin(spin)
     kk = vpq%kpts_spin(:, ik, spin)
     ! Note symrel option
     if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, 1, kk, mapl_k) /= 0) then
       write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k-point could not be generated from a symmetrical one.",trim(ltoa(kk))
       ABI_ERROR(msg)
     end if
     ik_ibz = mapl_k(1)
     do ib=1,vpq%nb_spin(spin)
       band = bstart + ib - 1; bks_mask(band, ik_ibz, spin) = .True.
     end do
   end do
 end do

 ! Impose istwfk = 1 for all k-points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
               dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft, &
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print([std_out], header="Wavefunctions for varpeq_plot")

 if (dtset%boxcutmin >= two) then
   call wrtout(std_out, " To reduce the size of the FFT mesh and the size of the XSF file, reduce boxcutmin from 2 to e.g. 1.1")
 end if

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call wrtout(std_out, " Computing mpw and gmax needed to allocate workspace array.")
 do spin=1,nsppol
   nk = vpq%nk_spin(spin)
   call ephtk_get_mpw_gmax(nk, vpq%kpts_spin(:, 1:nk, spin), dtset%ecut, cryst%gmet, mpw, gmax, comm, &
                           init_with_zero=spin==1)
 end do

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 call wrtout(std_out, sjoin(" Building supercell from ngkpt:", ltoa(vpq%ngkpt)))
 ! note xyz_order="xyz"
 call kptrlatt_from_ngkpt(vpq%ngkpt, kptrlatt_)
 call scell_k%init(cryst%natom, kptrlatt_, cryst%rprimd, cryst%typat, cryst%xcart, cryst%znucl, xyz_order="xyz")

 nkbz = product(vpq%ngkpt)
 call supercell_fft(vpq%ngkpt, ngfft, sc_nfft, sc_ngfft, sc2uc, scred)
 ABI_FREE(scred)

 call wrtout(std_out, " Computing polaron wavefunction in the real-space supercell...")

 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(gbound_k, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(ur_k, (wfd%nfft*nspinor))
 ABI_MALLOC(ceikr, (sc_nfft*nspinor))
 ABI_CALLOC(pol_wf, (sc_nfft*nspinor, vpq%nstates, nsppol)) ! Init output with zeros.

 cnt = 0
 do spin=1,nsppol
   bstart = vpq%brange_spin(1, spin)
   nk = vpq%nk_spin(spin)
   do ik=1, nk
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism inside comm.

     kk = vpq%kpts_spin(:, ik, spin)
     if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, 1, kk, mapl_k) /= 0) then
       write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k-point could not be generated from a symmetrical one.",trim(ltoa(kk))
       ABI_ERROR(msg)
     end if

     ik_ibz = mapl_k(1); isym_k = mapl_k(2); trev_k = mapl_k(6); g0_k = mapl_k(3:5)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     kk_ibz = ebands%kptns(:, ik_ibz)
     istwf_k_ibz = wfd%istwfk(ik_ibz); npw_kq_ibz = wfd%npwarr(ik_ibz)

     ! Get npw_k, kg_k for this k-point.
     call wfd%get_gvec_gbound(cryst%gmet, dtset%ecut, kk, ik_ibz, isirr_k, dtset%nloalg, &  ! in
                              istwf_k, npw_k, kg_k, nkpg_k, kpg_k, gbound_k)                ! out

     ABI_MALLOC(ug_k, (2, npw_k*nspinor))
     kk_sc = kk * vpq%ngkpt
     call calc_ceikr(kk_sc, sc_ngfft, sc_nfft, nspinor, ceikr)

     ! Phase due to the primitive translation (default: 0)
     cphase_tr = exp(-j_dpc * two_pi * dot_product(kk, dtset%varpeq_trvec))

     do ib=1,vpq%nb_spin(spin)
       band = bstart + ib - 1
       call wfd%rotate_cg(band, ndat1, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                          cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_k)
       !print *, "int_omega dr |u(r)}^2:", sum(abs(ur_k) ** 2) / wfd%nfft

       do ip=1,vpq%nstates
         a_nk = vpq%a_spin(ib, ik, ip, spin)*cphase_tr
         do spinor=1,nspinor
           spad = (spinor - 1) * sc_nfft
           do ir=1,sc_nfft
             uc_idx = sc2uc(ir) + (spinor - 1) * wfd%nfft
             irsp = ir + spad
             pol_wf(irsp, ip, spin) = pol_wf(irsp, ip, spin) + a_nk * ur_k(uc_idx) * ceikr(irsp)
           end do
         end do
       end do ! ip
     end do ! ib

     ABI_FREE(ug_k)
     ABI_FREE(kpg_k)
   end do ! ik
 end do ! spin

 ABI_FREE(kg_k)
 ABI_FREE(work)
 ABI_FREE(ur_k)
 ABI_FREE(gbound_k)
 ABI_FREE(sc2uc)
 ABI_FREE(ceikr)
 call krank_ibz%free()

 ! Collect pol_wf on the master rank who's gonna write the polaron density in XSF format.
 call xmpi_sum_master(pol_wf, master, comm, ierr)

 if (my_rank == master) then
   pol_wf = pol_wf / (nkbz * sqrt(cryst%ucvol))
   ABI_MALLOC(pol_rho, (sc_nfft))

   ! Here decide if we are gonna write the polaron wavefunctions with_diplaced atoms or not.
   use_displaced_scell = all(kptrlatt_ == qptrlatt_) .and. have_scell_q
   num_writes = 1; if (use_displaced_scell) num_writes = 2

   if (nspinor == 1) then
     ! Handle spin-polarized case by writing two XSF files.
     do spin=1,nsppol
       ! Handle multiple polaronic states for each spin.
       do ip=1,vpq%nstates
         write(msg, "(2(a,i0),a,es16.6)")&
           " For spin: ", spin, ": pstate: ", ip, ": 1/N_k \sum_nk |A_nk|^2 = ", sum(abs(vpq%a_spin(:,:,ip,spin))**2) / nkbz
         call wrtout(units, msg)
         pol_rho = abs(pol_wf(:, ip, spin))**2
         write(msg, "(2(a,i0),a,es16.6)")" Polaron density for spin: ", spin, ": pstate: ", ip, " integrates to: ", sum(pol_rho) * cryst%ucvol/wfd%nfft
         call wrtout(units, msg)
         write(msg, "(a,es16.6)")" maxval(abs(aimag(pol_wf))): ", maxval(abs(aimag(pol_wf(:, ip, spin))))
         call wrtout(units, msg)
         call center_and_spread(cryst, vpq%ngkpt, sc_ngfft, pol_rho, center_cart, spread, units)

         do ii=1,num_writes
           if (ii == 1) then
             xcart_ptr => scell_k%xcart
             path = strcat(dtfil%filnam_ds(4), "_pstate_", itoa(ip), "_POLARON.xsf")
             if (nsppol == 2) path = strcat(dtfil%filnam_ds(4), strcat("_spin_", itoa(spin)), "_pstate_", itoa(ip), "_POLARON.xsf")
             call wrtout(units, strcat("- Writing the polaron wavefunction with undisplaced atoms to: ", path))
           else
             path = strcat(dtfil%filnam_ds(4), "_pstate_", itoa(ip), "_POLARON_DISPL.xsf")
             if (nsppol == 2) path = strcat(dtfil%filnam_ds(4), strcat("_spin_", itoa(spin)), "_pstate_", itoa(ip), "_POLARON_DISPL.xsf")
             call wrtout(units, strcat("- Writing the polaron wavefunction with diplaced atoms to: ", path))

             ! Here we displace the atoms in the supercell for this spin (only master has the correct values)
             scell_q%xcart = scell_q%xcart_ref + sc_displ_cart_re(:,:,ip,spin)
             xcart_ptr => scell_q%xcart
           end if

           call write_xsf(path, sc_ngfft(1), sc_ngfft(2), sc_ngfft(3), pol_rho, scell_k%rprimd, origin0, &
                          scell_k%natom, scell_k%ntypat, scell_k%typat, xcart_ptr, scell_k%znucl, 0)
         end do ! ii
       end do ! ip
     end do ! spin

   else
     ! Spinor wavefunctions
     do ip=1,vpq%nstates
       write(msg, "(2(a,i0),a,es16.6)")&
         " For spin: ", spin, ": pstate: ", ip, ": 1/N_k \sum_nk |A_nk|^2 = ", sum(abs(vpq%a_spin(:,:,ip,spin))**2) / nkbz
       call wrtout(units, msg)

       pol_rho(:) = abs(pol_wf(1:sc_nfft, ip, 1)) ** 2
       pol_rho(:) = abs(pol_wf(sc_nfft+1:, ip, 1)) ** 2 + pol_rho(:)
       write(msg, "(2(a,i0),a,es16.6)")" Polaron density for spin: ", spin, ": pstate: ", ip, " integrates to: ", sum(pol_rho) * cryst%ucvol/wfd%nfft
       call wrtout(units, msg)
       write(msg, "(a,es16.6)")" maxval(abs(aimag(pol_wf))): ", maxval(abs(aimag(pol_wf(:, ip, 1))))
       call wrtout(units, msg)
       call center_and_spread(cryst, vpq%ngkpt, sc_ngfft, pol_rho, center_cart, spread, units)

       spin = 1
       do ii=1,num_writes
         if (ii == 1) then
           xcart_ptr => scell_k%xcart
           path = strcat(dtfil%filnam_ds(4), "_pstate_", itoa(ip), "_POLARON.xsf")
           call wrtout(units, strcat("- Writing the polaron wavefunction with undiplaced atoms to: ", path))
         else
           path = strcat(dtfil%filnam_ds(4), "_pstate_", itoa(ip), "_POLARON_DISPL.xsf")
           call wrtout(units, strcat("- Writing the polaron wavefunction with diplaced atoms to: ", path))

           ! Here we displace the atoms in the supercell for this spin (only master has the correct values)
           scell_q%xcart = scell_q%xcart_ref + sc_displ_cart_re(:,:,ip,spin)
           xcart_ptr => scell_q%xcart
         end if

         call write_xsf(path, sc_ngfft(1), sc_ngfft(2), sc_ngfft(3), pol_rho, scell_k%rprimd, origin0, &
                        scell_k%natom, scell_k%ntypat, scell_k%typat, xcart_ptr, scell_k%znucl, 0)
       end do ! ii
     end do ! ip
   end if
   ABI_FREE(pol_rho)
 end if ! master

 call cwtime_report(" Computation of polaron wavefunction completed", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10)

 ABI_SFREE(sc_displ_cart_re)
 ABI_FREE(pol_wf)

 call wfd%free()
 call scell_q%free()
 call scell_k%free()
 call vpq%free()

contains
integer function vid(var_name)
  character(len=*),intent(in) :: var_name
  vid = nctk_idname(ncid, var_name)
end function vid

end subroutine varpeq_plot
!!***

subroutine center_and_spread(prim_cryst, ncells, sc_ngfft, rhor, center_cart, spread, units)

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: prim_cryst
 integer,intent(in) :: ncells(3), sc_ngfft(18)
 real(dp),intent(in) :: rhor(sc_ngfft(1), sc_ngfft(2), sc_ngfft(3))
 !real(dp),intent(in) :: rhor(sc_ngfft(4), sc_ngfft(5), sc_ngfft(6)) ! FIXME
 real(dp),intent(out) :: center_cart(3), spread
 integer,intent(in) :: units(:)

!Local variables-------------------------------
 integer :: i1, i2, i3, nfft, num_cells
 real(dp) :: rr(3), rcart(3), rmr0(3), sc_rprimd(3,3), center_red(3), r2_mean, fwhm
 character(len=500) :: msg
!----------------------------------------------------------------------

 nfft = product(sc_ngfft(1:3))
 sc_rprimd(:,1) = prim_cryst%rprimd(:,1) * ncells(1)
 sc_rprimd(:,2) = prim_cryst%rprimd(:,2) * ncells(2)
 sc_rprimd(:,3) = prim_cryst%rprimd(:,3) * ncells(3)
 num_cells = product(ncells)

 ! Compute center_cart = \int r rhor(r) dr
 center_cart = zero !; r2_mean = zero
 do i3=1,sc_ngfft(3)
   rr(3) = (i3 - one) / sc_ngfft(3)
   do i2=1,sc_ngfft(2)
     rr(2) = (i2 - one) / sc_ngfft(2)
     do i1=1,sc_ngfft(1)
       rr(1) = (i1 - one) / sc_ngfft(1)
       ! Go to cartesian coordinates.
       rcart = matmul(sc_rprimd, rr)
       center_cart = center_cart + rhor(i1, i2, i3) * rcart
       r2_mean = r2_mean + rhor(i1, i2, i3) * dot_product(rcart, rcart)
     end do
   end do
 end do

 !center_cart = matmul(sc_rprimd, center_cart)
 center_cart = center_cart * num_cells * prim_cryst%ucvol / nfft
 r2_mean = r2_mean * num_cells * prim_cryst%ucvol / nfft
 spread = sqrt(r2_mean - dot_product(center_cart, center_cart))
 write(msg, "(a,2(es16.6,a))")" Polaron spread: ", spread, " (Bohr)", spread * Bohr_Ang, " (Ang)"
 call wrtout(units, msg)

 ! Compute = \int (r - center_cart)^2 rhor dr = <r^2> - <r>^2
 spread = zero
 do i3=1,sc_ngfft(3)
   rr(3) = (i3 - one) / sc_ngfft(3)
   do i2=1,sc_ngfft(2)
     rr(2) = (i2 - one) / sc_ngfft(2)
     do i1=1,sc_ngfft(1)
       rr(1) = (i1 - one) / sc_ngfft(1)
       ! Go to cartesian coordinates.
       rmr0 = matmul(sc_rprimd, rr) - center_cart
       spread = spread + rhor(i1, i2, i3) * dot_product(rmr0, rmr0)
     end do
   end do
 end do
 spread = sqrt(spread * num_cells * prim_cryst%ucvol / nfft)

 call xcart2xred(1, prim_cryst%rprimd, center_cart, center_red)
 call wrtout(units, sjoin(" Polaron center in Cartesian coordinates: ", ltoa(center_cart), " (Bohr)"))
 call wrtout(units, sjoin(" Fractional coordinates in termso of the primitive cell:", ltoa(center_red)))
 write(msg, "(a,2(es16.6,a))")" Polaron spread: ", spread, " (Bohr)", spread * Bohr_Ang, " (Ang)"
 call wrtout(units, msg)
 ! full width at half-maximum
 fwhm = (two * sqrt(two * log(two))) * spread
 write(msg, "(a,2(es16.6,a))")" Full width at half-maximum (FWHM) ", fwhm, "(Borh) ", fwhm * Bohr_Ang, " (Ang)"
 call wrtout(units, msg)

end subroutine center_and_spread

end module m_varpeq
!!***
