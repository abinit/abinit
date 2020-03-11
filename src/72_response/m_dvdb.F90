!!****m* ABINIT/m_dvdb
!! NAME
!!  m_dvdb
!!
!! FUNCTION
!!  Objects and methods to extract data from the DVDB file.
!!  The DVDB file is Fortran binary file with a collection of DFPT potentials
!!  associated to the different phonon perturbations (idir, ipert, qpt).
!!  DVDB files are produced with the `mrgdv` utility and used in the EPH code
!!  to compute the matrix elements <k+q| dvscf_{idir, ipert, qpt} |k>.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (MG,GA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_dvdb

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_distribfft
 use m_nctk
 use m_sort
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_ddb
 use m_dtset

 use defs_abitypes,   only : mpi_type
 use m_fstrings,      only : strcat, sjoin, itoa, ktoa, ltoa, ftoa, yesno, endswith
 use m_time,          only : cwtime, cwtime_report, sec2str, timab
 use m_io_tools,      only : open_file, file_exists, delete_file
 use m_numeric_tools, only : wrap2_pmhalf, vdiff_t, vdiff_eval, vdiff_print, l2int
 use m_symtk,         only : mati3inv, matr3inv, littlegroup_q
 use m_geometry,      only : littlegroup_pert, irreducible_set_pert, mkradim, xcart2xred
 use m_dynmat,        only : canat9, get_bigbox_and_weights
 use m_copy,          only : alloc_copy
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_ioarr,         only : read_rhor
 use m_fftcore,       only : ngfft_seq
 use m_fft_mesh,      only : rotate_fft_mesh, times_eigr, times_eikr, ig2gfft, get_gftt, calc_ceikr, calc_eigr
 use m_fft,           only : fourdp, zerosym
 use m_crystal,       only : crystal_t
 use m_kpts,          only : kpts_ibz_from_kptrlatt, listkk
 use m_spacepar,      only : symrhg, setsym
 use m_fourier_interpol,only : fourier_interpol
 use m_pawrhoij,      only : pawrhoij_type

 implicit none

 private
!!***

 ! Version 1: header + vscf1(r)
 ! Version 2: header + vscf1(r) + record with rhog1(G=0)
 integer,public,parameter :: dvdb_last_version = 2

 integer,private,parameter :: DVDB_NOMODE    = 0
 integer,private,parameter :: DVDB_READMODE  = 1
 integer,private,parameter :: DVDB_WRITEMODE = 2

 ! FIXME
 !real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors

 ! Uncomment this line and recompile to use double or single precision cache
 !integer,private,parameter :: QCACHE_KIND = dp
 integer,private,parameter :: QCACHE_KIND = sp

 type, private :: qcache_entry_t

   real(QCACHE_KIND), allocatable :: v1scf(:,:,:,:)
     ! v1scf(cplex, nfftf, nspden, my_npert))
     ! cached potentials (distributed inside comm_pert)

 end type qcache_entry_t

!----------------------------------------------------------------------

!!****t* m_dvdb/qcache_t
!! NAME
!! qcache_t
!!
!! FUNCTION
!!
!! SOURCE

 type, private :: qcache_t

   integer :: maxnq = 0
    ! Max number of q-points in Vscf(q) stored in the cache.
    ! Note this is not the size of the key array that is dimensioned with
    ! dvdb%nqpt or nqibz depending on the type of cache in use.

   integer :: nqpt
    ! Number of q-points in the key array

   !integer :: nspden, my_npert, nfft

   real(dp) :: max_mbsize = zero
    ! Max cache size in megabytes.
    !    < 0 to allocate all q-points.
    !    0 has not effect.
    !    > 0 for cache with automatically computed nqpt points.

   real(dp) :: onepot_mb = zero
    ! Size of 1 DFPT potential in megabytes

   logical :: use_3natom_cache = .False.
    ! True if v1scf_3natom_qibz cache is used.

   integer :: stats(4)
    ! Total number of calls, number of cache hit in v1scf_3natom_qibz cache, hit in key cache, cache misses.

   integer,allocatable :: count_qused(:)
    ! count_qused(nq)
    ! Number of times this q-point has been used in dvdb_readsym_qbz

   integer(i1b),allocatable :: itreatq(:)
    ! itreatq(nq)
    ! Table used to distribute q-points in the IBZ among procs.
    ! 0 if this MPI rank should not store and treat this q-point.
    ! 1 or 2-9

   real(dp), allocatable :: v1scf_3natom_qibz(:,:,:,:)
    ! v1scf_3natom(cplex, nfftf, nspden, 3*natom))
    ! Store all the 3*natom perturbations associated to `stored_iqibz_cplex`
    ! Used to reduce the number of MPI communications required to collect all the potentials
    ! on the MPI rank before computing V(Sq) from V(q).
    ! I assume client code is looping over stars/shells. cplex value is stored in `stored_iqibz_cplex`.

   integer :: stored_iqibz_cplex(2) = huge(1)
    ! The index of the qpoint in the IBZ/DVDB file associated to v1scf_3natom_qibz.

   type(qcache_entry_t), allocatable :: key(:)
    ! key(nq)
    ! array of v1scf potentials (ony a subset is usually allocated)

 contains

    procedure :: free => qcache_free
    ! Release dynamic memory.

    procedure :: report_stats => qcache_report_stats
    ! Print info on q-cache stats and reset counters.

    procedure :: get_mbsize => qcache_get_mbsize
    ! Return the (allocated) size of the cache in Mb.

    procedure :: make_room => qcache_make_room

 end type qcache_t
!!***

!----------------------------------------------------------------------

!!****t* m_dvdb/dvdb_t
!! NAME
!!  dvdb_t
!!
!! FUNCTION
!!  Database of DFPT results. The database contains `numv1` perturbations
!!  and the corresponding first order local potentials in real space on the FFT mesh.
!!  Note that one can have different FFT meshes for the different perturbations.
!!  Provides methods to Fourier interpolate the potentials including the
!!  treatment of long-range behaviour in the FT interpolation in polar semiconductors.
!!
!! NOTES
!!  natom, nspden, nspinor, and usepaw are global variables in the sense that it's not possible to add
!!  new entries to the database if these dimension differ from the global ones.
!!
!! SOURCE

 type,public :: dvdb_t

  integer :: fh
   ! file handler
   ! Fortran unit number if iomode==IO_MODE_FORTRAN
   ! MPI file handler if iomode==IO_MODE_MPI

  integer :: comm
  ! Global MPI communicator used for IO.

  integer :: comm_rpt = xmpi_comm_self
   ! MPI communicator used to distributed R-points.

  integer :: nprocs_rpt = 1
   ! Number of cpus for parallelism over R-points.

  integer :: me_rpt = 0
   ! My rank in comm_rpt.

  integer :: comm_pert = xmpi_comm_self
   ! MPI communicator for parallelism over atomic perturbations.

  integer :: nprocs_pert = 1
   ! Number of cpus for parallelism over atomic perturbations.

  integer :: me_pert = 0
   ! My rank in comm over atomic perturbations.

  integer :: my_npert
   ! Number of atomic perturbations or phonon modes treated by this MPI rank

  integer,allocatable :: my_pinfo(:,:)
    ! my_pinfo(3, my_npert)
    ! my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
    ! my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
    ! my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3

  integer,allocatable :: pert_table(:,:)
    ! pert_table(2, natom3)
    !     pert_table(1, npert): rank of the processor treating this atomic perturbation.
    !     pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.

  integer :: version
  ! File format version read from file.

  integer :: iomode = IO_MODE_FORTRAN
  ! Method used to access the DVDB file:
  !   IO_MODE_FORTRAN for usual Fortran IO routines
  !   IO_MODE_MPI if MPI/IO routines.

  integer :: rw_mode = DVDB_NOMODE
   ! (Read|Write) mode

  integer :: current_fpos
  ! The current position of the file pointer used for sequential access with Fortran-IO

  integer :: numv1
  ! Number of v1 potentials present in file.

  integer :: nqpt
  ! Number of q-points

  integer :: natom
   ! Number of atoms

  integer :: natom3
   ! 3 * natom

  integer :: nspden
   ! Number of spin density components

  integer :: nsppol
   ! Number of spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: usepaw
   ! 1 if PAW calculation, 0 otherwise

  integer :: mpert
   ! Maximum number of perturbations

  integer :: my_nrpt = 0
  ! Number of real space points used for Fourier interpolation treated by this MPI rank.

  integer :: nrtot = 0
  ! Total Number of real space points used for Fourier interpolation.

  integer :: prtvol = 0
   ! Verbosity level

  real(dp) :: qdamp = -one
   ! Exponential damping used in the Fourier transform of the long-range potentials
   ! Use negative value to deactivate damping.

  logical :: debug = .False.
   ! Debug flag

  logical :: has_dielt = .False.
  ! True if the dielectric tensor is available.

  logical :: has_zeff = .False.
  ! True if Born effective charges are available.

  logical :: has_quadrupoles = .False.
  ! True if quadrupoles are available.

  logical :: has_efield = .False.
  ! True if electric field perturbations are available.

  integer :: add_lr = 1
   ! Flag defining the treatment of the long range component in the interpolation of the DFPT potentials.
   ! 0 --> No treatment
   ! 1 --> Remove LR model when building W(R,r). Add it back after W(R,r) --> v(q) Fourier interpolation
   !       This is the standard approach for polar materials.
   ! -1 --> Remove LR model when building W(R,r). DO NOT reintroduce it after Fourier interpolation.
   !       This procedure should be used for homopolar materials with (spurious) non-zero BECS
   !       in order to remove the long range component from the DFPT potentials.

  integer :: symv1 = 0
   ! Flag for the symmetrization of v1 potentials.
   ! 0 --> No symmetrization
   ! 1 --> Symmetrization in real space
   ! 2 --> Call v1phq_complete after interpolation of the potentials in ftinterp_qpt

  type(qcache_t) :: qcache
    ! Cache used to store potentials if Fourier interpolation is not used
    ! (see dvdb_readsym_qbz for the implementation)

  type(qcache_t) :: ft_qcache
    ! Cache used to store potentials if Fourier interpolation is used
    ! (see dvdb_get_ftqbz for the implementation)

  character(len=fnlen) :: path = ABI_NOFILE
   ! File name of the DVDB file.

  real(dp) :: dielt(3, 3) = zero
   ! Dielectric tensor in Cartesian coordinates.
   ! Used to deal with the long-range component in the Fourier interpolation.

  integer,allocatable :: pos_dpq(:,:,:)
   ! pos_dpq(3, mpert, nqpt)
   ! The position of the (idir, ipert, iqpt) potential in the file (in units of POT1 blocks)
   ! 0 if the corresponding entry is not available.

  integer,allocatable :: cplex_v1(:)
  ! cplex_v1(numv1)
  ! The value of cplex for each v1(cplex*nfft, nspden) potential
  ! 2 if the potential is complex, 1 if real (q==Gamma)

  integer,allocatable :: symq_table(:,:,:,:)
  ! symq(4,2,nsym,nqpt)
  ! Table computed by littlegroup_q for all q-points found in the DVDB.
  !   three first numbers define the G vector;
  !   fourth number is zero if the q-vector is not preserved, is 1 otherwise
  !   second index is one without time-reversal symmetry, two with time-reversal symmetry

  integer :: ngfft(18) = -1
   ! Info on the FFT to be used for the potentials.

  integer,allocatable :: iv_pinfoq(:,:)
   !iv_pinfoq(4, numv1)
   !  iv_pinfoq(1, iv1) gives the `idir` index of the iv1 potential
   !  iv_pinfoq(2, iv1) gives the `ipert` index of the iv1 potential
   !  iv_pinfoq(3, iv1) gives `pertcase`=idir + (ipert-1)*3
   !  iv_pinfoq(4, iv1) gives the `iqpt` index of the iv1 potential

  integer,allocatable :: ngfft3_v1(:,:)
   ! ngfft3_v1(3, numv1)
   ! The FFT mesh used for each v1 potential (the one used to store data in the file).

  integer,allocatable :: my_irpt2tot(:)
  ! Mapping my_irpt index to full list of R-points.

  real(dp),allocatable :: qpts(:,:)
   ! qpts(3,nqpt)
   ! List of q-points in reduced coordinates.

  real(dp),allocatable :: my_rpt(:,:)
  ! my_rpt(3, my_nrpt)
  ! Real space points for Fourier interpolation (MPI distributed if nprocs_rpt > 1)

  real(kind=dp),allocatable :: v1scf_rpt(:,:,:,:,:)
  ! DFPT potential in the real space supercell representation.
  ! v1scf_rpt(2, my_nrpt, nfft, nspden, my_npert)

  real(dp),allocatable :: my_wratm(:,:)
  ! my_wratm(my_nrpt, minatom:maxatom)
  ! Weight for the FT associated to the atom and the R vector.

  real(dp),allocatable :: rhog1_g0(:,:)
  ! rhog1_g0(2, numv1)
  ! G=0 component of rhog1. Used to treat the long range component in (polar) semiconductors.
  ! NB: For the time being, this quantity is not used. Long range term is treated with Verdi's model.

  real(dp),allocatable :: zeff(:,:,:)
  ! zeff(3, 3, natom)
  ! Effective charges on each atom, versus electric field and atomic displacement in Cartesian coordinates.
  ! Used to deal with the long-range componenent in the Fourier interpolation.

  real(dp),allocatable :: zeff_raw(:,:,:)
  ! Raw Effective charges i.e. values before enforcing the charge-neutrality condition.

  real(dp),allocatable :: qstar(:,:,:,:)
  ! qstar(3, 3, 3, natom)
  ! dynamical quadrupole in Cartesian coordinates.
  ! First two dimension are associated to the q-point, then atomic perturbation in Cart coords.

  real(dp),allocatable :: v1r_efield(:,:,:)
  ! v1r_efield(nfft, 3, nspden)
  ! First order potentials due to the three different directions of the electric field perturbation.
  ! Potentials are in r-space

  type(crystal_t) :: cryst
  ! Crystalline structure read from the the DVDB file.

  type(hdr_type) :: hdr_ref
  ! Header associated to the first potential in the DVDB. Used to backspace.
  ! Gives the number of Fortran records required to backspace the header
  ! Assume headers with same headform and same basic dimensions e.g. npsp

  type(mpi_type) :: mpi_enreg
  ! Internal object used to call fourdp

 contains

   procedure :: open_read => dvdb_open_read
   ! Open the file in read-only mode.

   procedure :: close => dvdb_close
   ! Close the DVDB file.

   procedure :: free => dvdb_free
   ! Release the memory allocated and close the file.

   procedure :: print => dvdb_print
   ! Print info on object.

   procedure :: findq => dvdb_findq
   ! Returns the index of the q-point.

   procedure :: find_qpts => dvdb_find_qpts
   ! Returns the index of a list of q-points.

   procedure :: set_pert_distrib => dvdb_set_pert_distrib

   procedure :: read_onev1 => dvdb_read_onev1
   ! Read and return the DFPT potential for given (idir, ipert, iqpt).

   procedure :: readsym_allv1 => dvdb_readsym_allv1
   ! Read and return all the 3*natom DFPT potentials (either from file or symmetrized)

   procedure :: readsym_qbz => dvdb_readsym_qbz
   ! Reconstruct the DFPT potential for a q-point in the BZ starting
   ! from its symmetrical image in the IBZ.

   procedure :: qcache_read => dvdb_qcache_read
   ! Allocate internal cache for potentials.
   ! Read potentials from DVDB file and store them in cache (COLLECTIVE routine)

   procedure :: qcache_update_from_file => dvdb_qcache_update_from_file
   ! Read selected potentials and update cache (COLLECTIVE routine)

   procedure :: list_perts => dvdb_list_perts
   ! Check if all the (phonon) perts are available taking into account symmetries.

   procedure :: ftinterp_setup => dvdb_ftinterp_setup
   ! Prepare the internal tables for Fourier interpolation.

   procedure :: get_maxw => dvdb_get_maxw

   procedure :: ftinterp_qpt => dvdb_ftinterp_qpt
   ! Fourier interpolation of potentials for given q-point

   procedure :: get_ftqbz => dvdb_get_ftqbz
   ! Retrieve Fourier interpolated potential for a given q-point in the BZ.
   ! Use cache to reduce number of slow FTs.

   procedure :: ftqcache_build => dvdb_ftqcache_build
   ! This function initializes the internal q-cache from W(R,r)

   procedure :: ftqcache_update_from_ft => dvdb_ftqcache_update_from_ft
   ! This function initializes the internal q-cache from W(R,r)

   procedure :: get_v1r_long_range => dvdb_get_v1r_long_range
   ! Long-range part of the phonon potential

   procedure :: load_ddb => dvdb_load_ddb
   ! Load information about the Born effective charges and dielectric tensor from a DDB file

   procedure :: interpolate_v1scf => dvdb_interpolate_v1scf
   ! Fourier interpolation of the phonon potentials

   procedure :: get_v1scf_rpt => dvdb_get_v1scf_rpt
   ! Fourier transform of the phonon potential from qpt to R

   procedure :: get_v1scf_qpt => dvdb_get_v1scf_qpt
   ! Fourier transform of the phonon potential from R to qpt

   procedure :: load_efield => dvdb_load_efield

   procedure :: interpolate_and_write => dvdb_interpolate_and_write
   ! Interpolate the phonon potentials and write a new DVDB file.

   procedure :: qdownsample => dvdb_qdownsample
   ! Downsample the q-mesh. Produce new DVDB file

   procedure :: write_v1qavg => dvdb_write_v1qavg
   ! Debugging tool used to test the model for the long-range of the potential.

 end type dvdb_t

 public :: dvdb_new                ! Initialize the object.

 ! Utilities
 public :: dvdb_merge_files        ! Merge a list of POT1 files.

 ! debugging tools.
 public :: dvdb_test_v1rsym        ! Check symmetries of the DFPT potentials.
 public :: dvdb_test_v1complete    ! Debugging tool used to test the symmetrization of the DFPT potentials.

 public :: dvdb_test_ftinterp      ! Test Fourier interpolation of DFPT potentials.

!----------------------------------------------------------------------

 !type, public star_t
 !  integer :: npts
 !  real(dp) :: weight
 !  real(dp) :: ibz_point(3)
 !  integer,allocatable :: pt2ibz_map(:,:)
 !  ! (6, npts)
 !  real(dp), allocatable :: points(:,:)
 !  ! points(3, npts)
 !end type kstars_t

 !type, public stars_t
 !  integer :: nstars
 !  type(star_t),allocatable :: star(:)
 !  ! star(nstars)
 !end type stars_t

contains
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_new
!! NAME
!!  dvdb_new
!!
!! FUNCTION
!!  Initialize the object from file. This is a COLLECTIVE procedure that must be called
!!  by each process in the MPI communicator comm.
!!
!! INPUTS
!!   path=DVDB Filename.
!!   comm=MPI communicator.
!!
!! PARENTS
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

type(dvdb_t) function dvdb_new(path, comm) result(new)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,timrev2=2
 integer :: iv1,ii,ierr,unt,fform,nqpt,iq,iq_found,cplex,trev_q
 integer :: idir,ipert,my_rank, nprocs, iatom, pertcase
 real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
 type(hdr_type) :: hdr1
!arrays
 integer,allocatable :: tmp_pos(:,:,:)
 real(dp),allocatable :: tmp_qpts(:,:)
 real(dp) :: tsec(2)

!************************************************************************

 ! Keep track of total time spent.
 call timab(1800, 1, tsec)

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 new%path = path; new%comm = comm; new%iomode = IO_MODE_FORTRAN

 call wrtout(std_out, sjoin("- Analyzing DVDB file: ", path, "..."))
 call cwtime(cpu, wall, gflops, "start")

 ! Master reads the header and builds useful tables
 if (my_rank == master) then

   if (open_file(path, msg, newunit=unt, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   read(unt, err=10, iomsg=msg) new%version
   read(unt, err=10, iomsg=msg) new%numv1

   ! Get important dimensions from the first header and rewind the file.
   call hdr_fort_read(new%hdr_ref, unt, fform)
   if (dvdb_check_fform(fform, "read_dvdb", msg) /= 0) then
     MSG_ERROR(sjoin("While reading:", path, ch10, msg))
   end if
   if (new%debug) call new%hdr_ref%echo(fform, 4, unit=std_out)

   rewind(unt)
   read(unt, err=10, iomsg=msg)
   read(unt, err=10, iomsg=msg)

   ! The code below must be executed by the other procs if MPI.
   new%natom = new%hdr_ref%natom
   new%natom3 = 3 * new%hdr_ref%natom
   new%nspden = new%hdr_ref%nspden
   new%nsppol = new%hdr_ref%nsppol
   new%nspinor = new%hdr_ref%nspinor
   new%usepaw = new%hdr_ref%usepaw
   ABI_CHECK(new%usepaw == 0, "PAW not yet supported")

   ! TODO: Write function to return mpert from natom!
   new%mpert = new%natom + 6

   ABI_MALLOC(tmp_qpts, (3, new%numv1))
   ABI_MALLOC(tmp_pos, (3, new%mpert, new%numv1))
   tmp_pos = 0

   ABI_MALLOC(new%cplex_v1, (new%numv1))
   ABI_MALLOC(new%ngfft3_v1, (3, new%numv1))
   ABI_MALLOC(new%iv_pinfoq, (4, new%numv1))
   ABI_MALLOC(new%rhog1_g0, (2, new%numv1))

   nqpt = 0
   do iv1=1,new%numv1
     call hdr_fort_read(hdr1, unt, fform)
     if (dvdb_check_fform(fform, "read_dvdb", msg) /= 0) then
       MSG_ERROR(sjoin("While reading hdr of v1 potential of index:", itoa(iv1), ch10, msg))
     end if

     ! Save cplex and FFT mesh associated to this perturbation.
     cplex = 2; if (hdr1%qptn(1)**2+hdr1%qptn(2)**2+hdr1%qptn(3)**2<1.d-14) cplex = 1
     new%cplex_v1(iv1) = cplex
     new%ngfft3_v1(:, iv1) = hdr1%ngfft(:3)

     ! Skip the records with v1.
     do ii=1,hdr1%nspden
       read(unt, err=10, iomsg=msg)
     end do
     ! Read rhog1_g0 (if available)
     new%rhog1_g0(:, iv1) = zero
     if (new%version > 1) read(unt, err=10, iomsg=msg) new%rhog1_g0(:, iv1)

     ! Check whether this q-point is already in the list.
     ! Assume qpoints are grouped so invert the iq loop for better performace.
     ! This is gonna be slow if lots of q-points and perturbations are not grouped.
     iq_found = 0
     do iq=nqpt,1,-1
       if (all(abs(hdr1%qptn - tmp_qpts(:,iq)) < tol14)) then
         iq_found = iq; exit
       end if
     end do

     ! pertcase = idir + (ipert-1)*3 where ipert=iatom in the interesting cases
     idir = mod(hdr1%pertcase-1, 3) + 1
     ipert = (hdr1%pertcase - idir) / 3 + 1

     ! Increment nqpt is new q-points and update tmp_pos
     if (iq_found == 0) then
       nqpt = nqpt + 1
       tmp_qpts(:, nqpt) = hdr1%qptn
       iq_found = nqpt
     end if
     tmp_pos(idir, ipert, iq_found) = iv1
     new%iv_pinfoq(:,iv1) = [idir, ipert, hdr1%pertcase, iq_found]

     call hdr1%free()
   end do

   ! Allocate arrays with correct nqpt dimension
   new%nqpt = nqpt
   ABI_MALLOC(new%qpts, (3, nqpt))
   new%qpts = tmp_qpts(:,1:nqpt)
   ABI_FREE(tmp_qpts)

   ABI_MALLOC(new%pos_dpq, (3, new%mpert, nqpt))
   new%pos_dpq = tmp_pos(:, :, 1:nqpt)
   ABI_FREE(tmp_pos)

   close(unt)
 end if

 ! Master broadcasts data.
 if (xmpi_comm_size(comm) > 1) then
   call xmpi_bcast(new%version, master, comm, ierr)
   call xmpi_bcast(new%numv1, master, comm, ierr)
   call xmpi_bcast(new%nqpt, master, comm, ierr)
   call new%hdr_ref%bcast(master, my_rank, comm)

   new%natom = new%hdr_ref%natom
   new%natom3 = 3 * new%hdr_ref%natom
   new%nspden = new%hdr_ref%nspden
   new%nsppol = new%hdr_ref%nsppol
   new%nspinor = new%hdr_ref%nspinor
   new%usepaw = new%hdr_ref%usepaw
   new%mpert = new%natom + 6

   if (my_rank /= master) then
     ABI_MALLOC(new%cplex_v1, (new%numv1))
     ABI_MALLOC(new%ngfft3_v1, (3, new%numv1))
     ABI_MALLOC(new%iv_pinfoq, (4, new%numv1))
     ABI_MALLOC(new%qpts, (3, new%nqpt))
     ABI_MALLOC(new%pos_dpq, (3, new%mpert, new%nqpt))
     ABI_MALLOC(new%rhog1_g0, (2, new%numv1))
   end if

   call xmpi_bcast(new%cplex_v1, master, comm, ierr)
   call xmpi_bcast(new%ngfft3_v1, master, comm, ierr)
   call xmpi_bcast(new%iv_pinfoq, master, comm, ierr)
   call xmpi_bcast(new%qpts, master, comm, ierr)
   call xmpi_bcast(new%pos_dpq, master, comm, ierr)
   call xmpi_bcast(new%rhog1_g0, master, comm, ierr)
 end if

 ! Init crystal_t from the hdr read from file.
 new%cryst = new%hdr_ref%get_crystal(timrev2)
 new%my_npert = new%natom3

 ! Init tables assuming no MPI distribution of perturbations.
 ABI_MALLOC(new%my_pinfo, (3, new%natom3))
 ABI_MALLOC(new%pert_table, (2, new%natom3))
 do iatom=1,new%natom
   do idir=1,3
     pertcase = idir + (iatom-1) * 3
     new%my_pinfo(:, pertcase) = [idir, iatom, pertcase]
     new%pert_table(:, pertcase) = [xmpi_comm_self, pertcase]
   end do
 end do

 ! Init Born effective charges
 ABI_CALLOC(new%zeff, (3, 3, new%natom))
 ABI_CALLOC(new%zeff_raw, (3, 3, new%natom))
 ABI_CALLOC(new%qstar, (3, 3, 3, new%natom))

 ! Internal MPI_type needed for calling fourdp!
 call initmpi_seq(new%mpi_enreg)

 ! Precompute symq_table for all q-points in the DVDB.
 ABI_ICALLOC(new%symq_table, (4, 2, new%cryst%nsym, new%nqpt))
 do iq=1,new%nqpt
   if (mod(iq, nprocs) /= my_rank) cycle ! MPI parallelism
   call littlegroup_q(new%cryst%nsym, new%qpts(:,iq), new%symq_table(:,:,:,iq), &
     new%cryst%symrec, new%cryst%symafm, trev_q, prtvol=0)
 end do
 call xmpi_sum(new%symq_table, comm, ierr)

 call cwtime_report("- dvdb_new", cpu, wall, gflops)
 call timab(1800, 2, tsec)

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(sjoin("Error while reading:", path, ch10, msg))

end function dvdb_new
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_open_read
!! NAME
!!  dvdb_open_read
!!
!! FUNCTION
!!  Open the file in read-only mode.
!!
!! INPUTS
!!   ngfft(18)=Info on the FFT mesh used for the DFPT potentials. Note that ngfft
!!     is the mesh used by the parent. In principle, it can differ from the one
!!     found in the file. In this case a Fourier interpolation is required.
!!   comm=MPI communicator
!!
!! PARENTS
!!      m_dvdb,m_gkk,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_open_read(db, ngfft, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: nprocs, unt, ii
 character(len=500) :: msg
!arrays
 character(len=fnlen) :: pot_paths(3)

!************************************************************************

 if (db%rw_mode /= DVDB_NOMODE) then
   MSG_ERROR("DVDB should be in DVDB_NOMODE when open_read is called.")
 end if
 db%rw_mode = DVDB_READMODE

 nprocs = xmpi_comm_size(comm)

 ! Initialize tables to call fourdp in sequential
 db%ngfft = ngfft
 call init_distribfft_seq(db%mpi_enreg%distribfft, 'c', ngfft(2), ngfft(3), 'all')
 call init_distribfft_seq(db%mpi_enreg%distribfft, 'f', ngfft(2), ngfft(3), 'all')

 ! Open the file.
 select case (db%iomode)
 case (IO_MODE_FORTRAN)
   if (open_file(db%path, msg, newunit=db%fh, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   read(db%fh, err=10, iomsg=msg)
   read(db%fh, err=10, iomsg=msg)
   db%current_fpos = 1

 case (IO_MODE_MPI)
   MSG_ERROR("MPI not coded")

 case default
   MSG_ERROR(sjoin("Unsupported iomode:", itoa(db%iomode)))
 end select

 ! Read potentials induced by electric fields
 ! This requires ngfft so for the time being we call it here.
 ! I should try to add efield perturbations to DVDB but then I also have to handle symmetrization wrt idir!
 if (file_exists("__EFIELD_POTS__")) then
   call wrtout(std_out, " Reading Efield potentials from EFIELD_POTS")
   if (open_file("__EFIELD_POTS__", msg, newunit=unt, form="formatted") /= 0) then
     MSG_ERROR(msg)
   end if
   do ii=1,3
    read(unt, "(a)") pot_paths(ii)
   end do
   close(unt)
   call db%load_efield(pot_paths, comm)
 end if

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(sjoin("Error while reading", db%path, ch10, msg))

end subroutine dvdb_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_close
!! NAME
!!  dvdb_close
!!
!! FUNCTION
!! Close the file
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_close(db)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(inout) :: db

!************************************************************************

 select case (db%iomode)
 case (IO_MODE_FORTRAN)
   close(db%fh)
 case default
   MSG_ERROR(sjoin("Unsupported iomode:", itoa(db%iomode)))
 end select

 db%rw_mode = DVDB_NOMODE

end subroutine dvdb_close
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_free
!! NAME
!!  dvdb_free
!!
!! FUNCTION
!! Close the file and release the memory allocated.
!!
!! PARENTS
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_free(db)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(inout) :: db

!************************************************************************

 ! integer arrays
 ABI_SFREE(db%my_pinfo)
 ABI_SFREE(db%pert_table)
 ABI_SFREE(db%pos_dpq)
 ABI_SFREE(db%cplex_v1)
 ABI_SFREE(db%symq_table)
 ABI_SFREE(db%iv_pinfoq)
 ABI_SFREE(db%ngfft3_v1)
 ABI_SFREE(db%my_irpt2tot)

 ! real arrays
 ABI_SFREE(db%qpts)
 ABI_SFREE(db%my_rpt)
 ABI_SFREE(db%v1scf_rpt)
 ABI_SFREE(db%my_wratm)
 ABI_SFREE(db%rhog1_g0)
 ABI_SFREE(db%zeff)
 ABI_SFREE(db%zeff_raw)
 ABI_SFREE(db%qstar)
 ABI_SFREE(db%v1r_efield)

 ! types
 call db%hdr_ref%free()
 call db%cryst%free()
 call destroy_mpi_enreg(db%mpi_enreg)

 ! Clean cache(s)
 call db%qcache%free()
 call db%ft_qcache%free()

 ! Close the file but only if we have performed IO.
 if (db%rw_mode == DVDB_NOMODE) return
 call db%close()

end subroutine dvdb_free
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_print
!! NAME
!!  dvdb_print
!!
!! FUNCTION
!!  Print info on the object.
!!
!! INPUTS
!! [unit]=the unit number for output
!! [prtvol]=verbosity level
!! [mode_paral]=either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      eph,m_dvdb,m_sigmaph,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_print(db, header, unit, prtvol, mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 class(dvdb_t),intent(in) :: db

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol,iv1,iq,idir,ipert,iatom
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt = std_out; if (present(unit)) my_unt = unit
 my_prtvol = 0   ; if (present(prtvol)) my_prtvol = prtvol
 my_mode = 'COLL'; if (present(mode_paral)) my_mode = mode_paral

 msg=' ==== Info on the dvdb% object ==== '
 if (present(header)) msg=' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(my_unt,"(a)")sjoin(" DVDB version:", itoa(db%version))
 write(my_unt,"(a)")sjoin(" File path:", db%path)
 write(my_unt,"(a)")sjoin(" Number of v1scf potentials:", itoa(db%numv1))
 write(my_unt,"(a)")sjoin(" Number of q-points in DVDB: ", itoa(db%nqpt))
 write(my_unt,"(a)")sjoin("-P Number of CPUs for parallelism over perturbations:", itoa(db%nprocs_pert))
 write(my_unt,"(a)")sjoin("-P Number of perturbations treated by this CPU:", itoa(db%my_npert))
 write(my_unt,"(a)")sjoin(" Option for symmetrization of v1scf(r):", itoa(db%symv1))
 write(my_unt,"(a)")" List of q-points: min(10, nqpt)"
 do iq=1,min(db%nqpt, 10)
   write(my_unt,"(a)")sjoin("[", itoa(iq),"]", ktoa(db%qpts(:,iq)))
 end do
 if (db%nqpt > 10) write(my_unt,"(a)")"..."

 write(my_unt,"(a)")sjoin(" Have dielectric tensor:", yesno(db%has_dielt))
 write(my_unt,"(a)")sjoin(" Have Born effective charges:", yesno(db%has_zeff))
 write(my_unt,"(a)")sjoin(" Have quadrupoles:", yesno(db%has_quadrupoles))
 write(my_unt,"(a)")sjoin(" Have electric field:", yesno(db%has_efield))
 write(my_unt,"(a)")sjoin(" Treatment of long-range part in V1scf:", itoa(db%add_lr))
 write(my_unt,"(a, f6.1)")" qdamp:", db%qdamp

 if (db%has_dielt) then
   write(my_unt, '(a,3(/,3es16.6))') ' Dielectric tensor in Cart coords:', &
     db%dielt(1,1), db%dielt(1,2), db%dielt(1,3), &
     db%dielt(2,1), db%dielt(2,2), db%dielt(2,3), &
     db%dielt(3,1), db%dielt(3,2), db%dielt(3,3)
 end if
 if (db%has_zeff) then
   call print_zeff(my_unt, db%zeff, db%cryst, title=' Born effectives charges in Cart coords:')
   !call print_zeff(my_unt, db%zeff_raw, db%cryst, title=' Born effectives charges before chneut: ')
 end if
 if (db%has_quadrupoles) then
   write(my_unt, '(a)') ' Dynamical Quadrupoles in Cartesian Coordinates: '
   do iatom=1,db%natom
     do idir=1,3
       write(my_unt,'(2(a,i0),3(/,3es16.6))')' iatom: ', iatom, ' idir: ', idir, &
         db%qstar(1,1,idir,iatom), db%qstar(1,2,idir,iatom), db%qstar(1,3,idir,iatom), &
         db%qstar(2,1,idir,iatom), db%qstar(2,2,idir,iatom), db%qstar(2,3,idir,iatom), &
         db%qstar(3,1,idir,iatom), db%qstar(3,2,idir,iatom), db%qstar(3,3,idir,iatom)
     end do
   end do
 end if

 if (my_prtvol > 0) then
   call db%cryst%print(header="Crystal structure in DVDB file")
   write(my_unt,"(a)")"FFT mesh for potentials on file:"
   write(my_unt,"(a)")"q-point, idir, ipert, ngfft(:3)"
   do iv1=1,db%numv1
     idir = db%iv_pinfoq(1, iv1); ipert = db%iv_pinfoq(2, iv1); iq = db%iv_pinfoq(4, iv1)
     write(my_unt,"(a)")sjoin(ktoa(db%qpts(:,iq)), itoa(idir), itoa(ipert), ltoa(db%ngfft3_v1(:,iv1)))
   end do
   write(my_unt, "(a)")"q-point, idir, ipert, rhog1(q,G=0)"
   do iv1=1,db%numv1
     idir = db%iv_pinfoq(1, iv1); ipert = db%iv_pinfoq(2, iv1); iq = db%iv_pinfoq(4, iv1)
     write(my_unt,"(a)")sjoin(ktoa(db%qpts(:,iq)), itoa(idir), itoa(ipert), &
       ftoa(db%rhog1_g0(1, iv1)), ftoa(db%rhog1_g0(2, iv1)))
   end do
 end if

contains

subroutine print_zeff(unt, zeff, cryst, title)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unt
 character(len=*),optional,intent(in) :: title
 type(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: zeff(3,3,cryst%natom)

!Local variables-------------------------------
!scalars
 integer :: iatom

! *************************************************************************

 if (present(title)) then
   write(unt, "(a)")trim(title)
 else
   write(unt, '(a)') ' Born effectives charges in Cartesian coordinates: '
 end if

 do iatom=1,cryst%natom
   write(unt,'(a,i0,1x,2a,3(/,3es16.6),a)')' iatom: ', iatom, ", type: ", cryst%symbol_iatom(iatom), &
     zeff(1,1,iatom), zeff(1,2,iatom), zeff(1,3,iatom), &
     zeff(2,1,iatom), zeff(2,2,iatom), zeff(2,3,iatom), &
     zeff(3,1,iatom), zeff(3,2,iatom), zeff(3,3,iatom), ch10
 end do

 write(unt,'(2a,3(/,3es16.6),a)')ch10,' Fulfillment of charge neutrality, \sum_{atom} Z^*_{ij,atom} = 0', &
   sum(zeff(1,1,:)), sum(zeff(1,2,:)), sum(zeff(1,3,:)), &
   sum(zeff(2,1,:)), sum(zeff(2,2,:)), sum(zeff(2,3,:)), &
   sum(zeff(3,1,:)), sum(zeff(3,2,:)), sum(zeff(3,3,:)), ch10

end subroutine print_zeff

end subroutine dvdb_print
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_pinfo
!! NAME
!!  dvdb_get_pinfo
!!
!! FUNCTION
!!  Return information on the perturbations available for a given q-point index.
!!
!! INPUTS
!!  iqpt=Index of the q-point
!!
!! OUTPUT
!!  nperts=Number of perturbations found.
!!  cplex=2 if potentials are complex, 1 for real
!!  pinfo(3,3*db%mpert)=Array with info on the perturbations present on file
!!     pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
!!     pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
!!     pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_get_pinfo(db, iqpt, cplex, pinfo) result(nperts)

!Arguments ------------------------------------
!scalars
 type(dvdb_t),intent(in) :: db
 integer,intent(in) :: iqpt
 integer,intent(out) :: cplex
!arrays
 integer,intent(out) :: pinfo(3,3*db%mpert)

!Local variables-------------------------------
!scalars
 integer :: idir,ipert,iv1

! *************************************************************************

 ! Get the number of perturbations computed for this iqpt
 pinfo = 0; cplex = 0; nperts = 0
 do ipert=1,db%natom ! selects atomic perturbations only.
    do idir=1,3
      iv1 = db%pos_dpq(idir,ipert,iqpt)
      if (iv1 /= 0) then
        nperts = nperts + 1
        pinfo(:, nperts) = [idir, ipert, idir + (ipert-1)*3]
        if (cplex == 0) cplex = db%cplex_v1(iv1)
        ABI_CHECK(cplex == db%cplex_v1(iv1), "cplex should be constant for given q!")
      end if
    end do
 end do

end function dvdb_get_pinfo
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_read_onev1
!! NAME
!!  dvdb_read_onev1
!!
!! FUNCTION
!!  Read the DFPT potential for the specified (idir, ipert, iqpt).
!!  Note that iqpt is the index in dvdb%qpts. Use dvdb_findq to
!!  get the index from the q-point in reduced coordinates.
!!
!! INPUTS
!!  idir=Direction of the perturbation
!!  ipert=Perturbation type.
!!  iqpt=Index of the q-point in dvdb%qpts
!!  cplex=1 if real, 2 if complex potentials.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT.
!!
!! OUTPUT
!!  ierr=Non-zero if error.
!!  v1scf(cplex*nfft, nspden)=DFT potential associated to (idir, ipert, iqpt).
!!  msg=String with error message if ierr /= 0.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_read_onev1(db, idir, ipert, iqpt, cplex, nfft, ngfft, v1scf, msg) result(ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,iqpt,cplex,nfft
 character(len=*),intent(out) :: msg
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: v1scf(cplex*nfft,db%nspden)

!Local variables-------------------------------
!scalars
 integer,save :: enough = 0
 integer :: iv1,ispden,nfftot_file,nfftot_out,ifft
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: ngfft_in(18),ngfft_out(18)
 real(dp),allocatable :: v1r_file(:,:),v1g_in(:,:),v1g_out(:,:)

! *************************************************************************

 ! Consistency checks
 ierr = 1
 iv1 = db%pos_dpq(idir,ipert,iqpt)

 if (iv1 == 0) then
   write(msg,"(3(a,i0))")"Cannot find idir: ",idir,", ipert: ",ipert,", iqpt:",iqpt
   return
 end if

 if (cplex /= db%cplex_v1(iv1)) then
   write(msg,"(2(a,i0))")"Wrong cplex. Expecting: ",db%cplex_v1(iv1),", received: ",cplex
   return
 end if

 ! Find (idir, ipert, iqpt) and skip the header.
 call dvdb_seek(db, idir, ipert, iqpt)
 ierr = my_hdr_skip(db%fh, idir, ipert, db%qpts(:,iqpt), msg)
 if (ierr /= 0) then
   msg = sjoin("In my_hdr_skip:", msg)
   return
 end if

 ! Read v1 from file.
 nfftot_out = product(ngfft(:3)); nfftot_file = product(db%ngfft3_v1(:3, iv1))

 if (all(ngfft(:3) == db%ngfft3_v1(:3, iv1))) then
   do ispden=1,db%nspden
     read(db%fh, err=10, iomsg=msg) (v1scf(ifft, ispden), ifft=1,cplex*nfftot_file)
   end do
 else
   ! The FFT mesh used in the caller differ from the one found in the DVDB --> Fourier interpolation
   ! TODO: Add linear interpolation as well.
   if (enough == 0) MSG_COMMENT("Performing FFT interpolation of DFPT potentials as input ngfft differs from ngfft_file.")
   enough = enough + 1
   ABI_MALLOC(v1r_file, (cplex*nfftot_file, db%nspden))
   do ispden=1,db%nspden
     read(db%fh, err=10, iomsg=msg) (v1r_file(ifft, ispden), ifft=1,cplex*nfftot_file)
   end do

   ! Call fourier_interpol to get v1scf on ngfft mesh.
   ngfft_in = ngfft; ngfft_out = ngfft
   ngfft_in(1:3) = db%ngfft3_v1(1:3, iv1); ngfft_out(1:3) = ngfft(1:3)
   ngfft_in(4:6) = ngfft_in(1:3); ngfft_out(4:6) = ngfft_out(1:3)
   ngfft_in(9:18) = 0; ngfft_out(9:18) = 0
   ngfft_in(10) = 1; ngfft_out(10) = 1

   call initmpi_seq(MPI_enreg_seq)
   ! Which one is coarse? Note that this part is not very robust and can fail!
   if (ngfft_in(2) * ngfft_in(3) < ngfft_out(2) * ngfft_out(3)) then
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft_in(2),ngfft_in(3),'all')
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft_out(2),ngfft_out(3),'all')
   else
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft_in(2),ngfft_in(3),'all')
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft_out(2),ngfft_out(3),'all')
   end if

   ABI_MALLOC(v1g_in,  (2, nfftot_file))
   ABI_MALLOC(v1g_out, (2, nfftot_out))

   call fourier_interpol(cplex,db%nspden,0,0,nfftot_file,ngfft_in,nfftot_out,ngfft_out,&
     MPI_enreg_seq,v1r_file,v1scf,v1g_in,v1g_out)

   ABI_FREE(v1g_in)
   ABI_FREE(v1g_out)
   ABI_FREE(v1r_file)
   call destroy_mpi_enreg(MPI_enreg_seq)
 end if

 ! Skip record with rhog1_g0 (if present)
 if (db%version > 1) read(db%fh, err=10, iomsg=msg)

 db%current_fpos = db%current_fpos + 1
 !write(std_out, *)"incr current_fpos", db%current_fpos

 return

 ! Handle Fortran IO error
10 continue
 ierr = 1
 msg = sjoin("Error while reading", db%path, ch10, msg)

end function dvdb_read_onev1
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_readsym_allv1
!! NAME
!!  dvdb_readsym_allv1
!!
!! FUNCTION
!!  Read all 3*natom DFPT potentials for the given iqpt (only atomic perturbations).
!!
!!  The routine will:
!!
!!     1) Reconstruct the potentials by symmetry if the DVDB contains less than 3*natom potentials.
!!     2) interpolate the data if the input FFT mesh defined by `ngfft` differs
!!        from the one used to store data in the file.
!!
!!  Note that iqpt is the index in dvdb%qpts. Use dvdb_findq to
!!  get the index from the q-point in reduced coordinates.
!!
!! INPUTS
!!  iqpt=Index of the q-point in dvdb%qpts
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  comm=MPI communicator
!!
!! OUTPUT
!!  cplex=1 if real, 2 if complex.
!!  v1scf(cplex, nfft, nspden, 3*natom)= v1scf potentials on the real-space FFT mesh for the 3*natom perturbations.
!!
!! PARENTS
!!      m_dvdb,m_gkk,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_readsym_allv1(db, iqpt, cplex, nfft, ngfft, v1scf, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqpt,nfft,comm
 integer,intent(out) :: cplex
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp) ABI_ASYNC ,allocatable,intent(out) :: v1scf(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ipc,npc,idir,ipert,pcase,my_rank,nproc,ierr,mu
 character(len=500) :: msg
!arrays
 integer :: pinfo(3,3*db%mpert),pflag(3, db%natom)
 real(dp) :: tsec(2)
 integer,allocatable :: requests(:)

! *************************************************************************

 ! Keep track of total time spent.
 call timab(1805, 1, tsec)

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Get number of perturbations computed for this iqpt as well as cplex.
 npc = dvdb_get_pinfo(db, iqpt, cplex, pinfo)
 ABI_CHECK(npc /= 0, "npc == 0!")

 ABI_MALLOC_OR_DIE(v1scf, (cplex, nfft, db%nspden, 3*db%natom), ierr)

 ! Master read all available perturbations and broadcasts data (non-blocking to overlap IO and MPI)
 ABI_MALLOC(requests, (npc))

 do ipc=1,npc
   idir = pinfo(1,ipc); ipert = pinfo(2,ipc); pcase = pinfo(3, ipc)
   if (my_rank == master) then
     if (db%read_onev1(idir, ipert, iqpt, cplex, nfft, ngfft, v1scf(:,:,:,pcase), msg) /= 0) then
       MSG_ERROR(msg)
     end if
   end if
   if (nproc > 1) call xmpi_ibcast(v1scf(:,:,:,pcase), master, comm, requests(ipc), ierr)
 end do

 if (nproc > 1) call xmpi_waitall(requests, ierr)
 ABI_FREE(requests)

 ! Return if all perts are available.
 if (npc == 3*db%natom) then
   if (db%symv1==1) then
     if (db%debug) write(std_out,*)"Potentials are available but will call v1phq_symmetrize because of symv1"
     do mu=1,db%natom3
       !if (mod(mu, nproc) /= my_rank) cycle ! MPI parallelism.
       idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1
       call v1phq_symmetrize(db%cryst,idir,ipert,db%symq_table(:,:,:,iqpt),ngfft,cplex,nfft,&
         db%nspden,db%nsppol,db%mpi_enreg,v1scf(:,:,:,mu))
       !call MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request)
     end do
   end if
   if (db%debug) write(std_out,*)"All perts available. Returning"
   return
 end if

 ! Perturbation are missing and we have to reconstruct them by symmetry.
 ! This is the common case when DFPT calculations are done for independent perturbations only.
 if (db%debug) then
   write(std_out,*)sjoin("Will use symmetries to recostruct:", itoa(3*db%natom - npc), "perturbations")
 end if

 ! 0 if pert is not available.
 ! 1 if pert is on file.
 ! 2 if pert has been reconstructed by symmetry.
 pflag = 0
 do ipc=1,npc
   pflag(pinfo(1,ipc), pinfo(2,ipc)) = 1
 end do

 call v1phq_complete(db%cryst,db%qpts(:,iqpt),ngfft,cplex,nfft,db%nspden,db%nsppol,db%mpi_enreg,db%symv1,pflag,v1scf)

 call timab(1805, 2, tsec)

end subroutine dvdb_readsym_allv1
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_readsym_qbz
!! NAME
!!  dvdb_readsym_qbz
!!
!! FUNCTION
!! This is the MAIN ENTRY POINT for client code.
!! Reconstruct the DFPT potential for a q-point in the BZ starting
!! from its symmetrical image in the IBZ. Implements caching mechanism to reduce IO.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  qbz(3)=Q-point in BZ.
!!  indq2db(6)=Symmetry mapping qbz --> DVDB qpoint produced by listkk.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  comm=MPI communicator (either xmpi_comm_self or comm for perturbations.
!!
!! OUTPUT
!!  cplex=1 if real, 2 if complex.
!!  v1scf(cplex, nfft, nspden, db%my_npert)= v1scf potentials on the real-space FFT mesh
!!   for the db%my_npert perturbations treated by this MPI rank.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_readsym_qbz(db, cryst, qbz, indq2db, cplex, nfft, ngfft, v1scf, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,comm
 integer,intent(out) :: cplex
 type(crystal_t),intent(in) :: cryst
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: indq2db(6)
 real(dp),intent(in) :: qbz(3)
 real(dp),allocatable,intent(out) :: v1scf(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: db_iqpt, itimrev, isym, npc, ierr, imyp, mu, root
 logical :: isirr_q, incache
!arrays
 integer :: pinfo(3,3*db%mpert), g0q(3), requests(db%natom3)
 real(dp) :: tsec(2)
 real(dp) ABI_ASYNC, allocatable :: work(:,:,:,:), work2(:,:,:,:)

! *************************************************************************
 ABI_UNUSED(qbz(1))

 ! Keep track of total time spent.
 call timab(1802, 1, tsec)

 db_iqpt = indq2db(1)
 db%qcache%count_qused(db_iqpt) = db%qcache%count_qused(db_iqpt) + 1
 db%qcache%stats(1) = db%qcache%stats(1) + 1

 ! IS(q_dvdb) + g0q = q_bz
 isym = indq2db(2); itimrev = indq2db(6) + 1; g0q = indq2db(3:5)
 isirr_q = (isym == 1 .and. itimrev == 1 .and. all(g0q == 0))

 if (db%qcache%use_3natom_cache .and. db%qcache%stored_iqibz_cplex(1) == db_iqpt .and. .not. isirr_q) then
   ! All 3 natom potentials for qibz are in cache. Symmetrize to get Sq without MPI communication.
   db%qcache%stats(2) = db%qcache%stats(2) + 1
   cplex = db%qcache%stored_iqibz_cplex(2)
   ABI_MALLOC(work2, (cplex, nfft, db%nspden, db%natom3))
   call v1phq_rotate(cryst, db%qpts(:, db_iqpt), isym, itimrev, g0q, ngfft, cplex, nfft, &
     db%nspden, db%nsppol, db%mpi_enreg, db%qcache%v1scf_3natom_qibz, work2, db%comm_pert)
   ! Extract my data from work2
   ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
   do imyp=1,db%my_npert
     v1scf(:,:,:,imyp) = work2(:,:,:,db%my_pinfo(3, imyp))
   end do
   ABI_FREE(work2)
   call timab(1802, 2, tsec); return
 end if

 ! Check whether db_iqpt is in cache.
 incache = .False.
 if (db%qcache%maxnq > 0) then
   ! Get number of perturbations computed for this iqpt as well as cplex.
   npc = dvdb_get_pinfo(db, db_iqpt, cplex, pinfo)
   ABI_CHECK(npc /= 0, "npc == 0!")

   ! Size of v1scf in qcache depends on db%my_npert
   if (allocated(db%qcache%key(db_iqpt)%v1scf)) then
      if (size(db%qcache%key(db_iqpt)%v1scf, dim=1) == cplex .and. size(db%qcache%key(db_iqpt)%v1scf, dim=2) == nfft) then
        ! Potential in cache --> copy it in output v1scf.
        ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
        v1scf = real(db%qcache%key(db_iqpt)%v1scf, kind=QCACHE_KIND)
        incache = .True.
        db%qcache%stats(3) = db%qcache%stats(3) + 1
      else
        ! This to handle the unlikely event in which the caller changes ngfft!
        MSG_WARNING("different cplex or nfft!")
      end if
   else
      !call wrtout(std_out, sjoin("Cache miss for db_iqpt. Will read it from file...", itoa(db_iqpt)))
      db%qcache%stats(4) = db%qcache%stats(4) + 1
   end if
 end if

 if (.not. incache) then
   ! Read the dvscf potentials in the IBZ for all 3*natom perturbations.
   ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom)
   call db%readsym_allv1(db_iqpt, cplex, nfft, ngfft, v1scf, comm)

   ! Store all 3 natom potentials for q in IBZ in cache.
   if (db%qcache%use_3natom_cache .and. db%qcache%stored_iqibz_cplex(1) /= db_iqpt .and. isirr_q) then
     if (cplex /= db%qcache%stored_iqibz_cplex(2)) then
       ABI_REMALLOC(db%qcache%v1scf_3natom_qibz, (cplex, nfft, db%nspden, db%natom3))
     end if
     db%qcache%v1scf_3natom_qibz = v1scf
     db%qcache%stored_iqibz_cplex = [db_iqpt, cplex]
   end if
 end if

 if (.not. isirr_q) then
   ! Must rotate db_iqpt to get potential for qpoint in the BZ.
   ! Be careful with the shape of output v1scf because the routine returns db%my_npert potentials.

   if (db%my_npert == db%natom3) then
     ABI_MALLOC(work, (cplex, nfft, db%nspden, db%natom3))
     work = v1scf
     call v1phq_rotate(cryst, db%qpts(:, db_iqpt), isym, itimrev, g0q, ngfft, cplex, nfft, &
       db%nspden, db%nsppol, db%mpi_enreg, work, v1scf, db%comm_pert)
     ABI_FREE(work)

   else
     ! Parallelism over perturbations.
     ABI_MALLOC(work2, (cplex, nfft, db%nspden, db%natom3))

     if (incache) then
       ! Cache is distributed --> have to collect all 3*natom perts inside db%comm_pert.
       ABI_MALLOC(work, (cplex, nfft, db%nspden, db%natom3))

       ! IBCAST is much faster than a naive xmpi_sum.
       call timab(1806, 1, tsec)
       do mu=1,db%natom3
         root = db%pert_table(1, mu)
         if (root == db%me_pert) then
           work(:,:,:,mu) = v1scf(:,:,:,db%pert_table(2, mu))
         end if
         call xmpi_ibcast(work(:,:,:,mu), root, db%comm_pert, requests(mu), ierr)
       end do
       call xmpi_waitall(requests, ierr)
       call timab(1806, 2, tsec)

       ! Store all 3 natom potentials for q in IBZ in cache.
       if (db%qcache%use_3natom_cache .and. db%qcache%stored_iqibz_cplex(1) /= db_iqpt) then
         if (cplex /= db%qcache%stored_iqibz_cplex(2)) then
           ABI_REMALLOC(db%qcache%v1scf_3natom_qibz, (cplex, nfft, db%nspden, db%natom3))
         end if
         db%qcache%v1scf_3natom_qibz = work
         db%qcache%stored_iqibz_cplex = [db_iqpt, cplex]
       end if

       ! Now rotate.
       call v1phq_rotate(cryst, db%qpts(:, db_iqpt), isym, itimrev, g0q, ngfft, cplex, nfft, &
         db%nspden, db%nsppol, db%mpi_enreg, work, work2, db%comm_pert)
       ABI_FREE(work)

     else
       ! All 3 natom have been read in v1scf by dvdb_readsym_allv1
       call v1phq_rotate(cryst, db%qpts(:, db_iqpt), isym, itimrev, g0q, ngfft, cplex, nfft, &
         db%nspden, db%nsppol, db%mpi_enreg, v1scf, work2, db%comm_pert)
     end if

     ! Reallocate v1scf with my_npert and extract data from work2.
     ABI_REMALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
     do imyp=1,db%my_npert
       v1scf(:,:,:,imyp) = work2(:,:,:,db%my_pinfo(3, imyp))
     end do
     ABI_FREE(work2)
   end if

 else
   ! Handle potentials read from file in case of parallelism over perturbations.
   if (.not. incache .and. db%my_npert /= db%natom3) then
     ABI_MALLOC(work, (cplex, nfft, db%nspden, db%my_npert))
     do imyp=1,db%my_npert
       work(:,:,:,imyp) = v1scf(:,:,:,db%my_pinfo(3, imyp))
     end do

     ABI_REMALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
     v1scf = work
     ABI_FREE(work)
   end if
 end if ! not isirr_q

 call timab(1802, 2, tsec)

end subroutine dvdb_readsym_qbz
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/qcache_new
!! NAME
!!  qcache_new
!!
!! FUNCTION
!!  Initialize qcache_t object from dimensions.
!!
!! INPUTS
!!  nqpt=Number of q-points (dvdb%nqpt or nqibz depending on cache type)
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  mbsize: Cache size in megabytes.
!!    < 0 to allocate all q-points.
!!    0 has not effect.
!!    > 0 for cache with automatically computed nqpt points.
!!  natom3, my_npert=Total number of perturbations and number of perturbations treated by this CPU.
!!  nspden=Number of spin density components.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(qcache_t) function qcache_new(nqpt, nfft, ngfft, mbsize, natom3, my_npert, nspden) result(qcache)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpt, nfft, natom3, nspden, my_npert
 real(dp),intent(in) :: mbsize
!arrays
 integer,intent(in) :: ngfft(18)

! *************************************************************************

 qcache%nqpt = nqpt
 ABI_ICALLOC(qcache%count_qused, (nqpt))
 ABI_MALLOC(qcache%key, (nqpt))
 ABI_MALLOC(qcache%itreatq, (nqpt))
 qcache%itreatq = 1
 qcache%stats = 0
 qcache%max_mbsize = mbsize

 qcache%onepot_mb = two * product(ngfft(1:3)) * nspden * QCACHE_KIND * b2Mb
 if (abs(mbsize) < tol3) then
   qcache%maxnq = 0
 else if (mbsize < zero) then
   qcache%maxnq = nqpt
 else
   qcache%maxnq = nint(mbsize / (qcache%onepot_mb * my_npert))
   qcache%maxnq = max(min(qcache%maxnq, nqpt), 1)
 end if

 ! Allocate cache with all the 3*natom perturbations.
 ! Disable it if no parallelism over perturbations.
 qcache%use_3natom_cache = .True.
 if (my_npert == natom3) qcache%use_3natom_cache = .False.
 qcache%stored_iqibz_cplex = huge(1)
 if (qcache%use_3natom_cache) then
   ABI_MALLOC(qcache%v1scf_3natom_qibz, (2, nfft, nspden, natom3))
 end if

 call wrtout(std_out, sjoin(" Using cache for Vscf(q) with MAX input size: ", ftoa(mbsize, fmt="f9.2"), " [Mb]"))
 call wrtout(std_out, sjoin(" Max number of q-points stored in memory: ", itoa(qcache%maxnq)))
 call wrtout(std_out, sjoin(" Use extra cache with 3 natom potentials: ", yesno(qcache%use_3natom_cache)))
 call wrtout(std_out, sjoin(" One DFPT potential requires: ", ftoa(qcache%onepot_mb, fmt="f9.2"), " [Mb]"))
 call wrtout(std_out, sjoin(" QCACHE_KIND: ", itoa(QCACHE_KIND)))

end function qcache_new
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_qcache_read
!! NAME
!!  dvdb_qcache_read
!!
!! FUNCTION
!!  This function initializes the internal q-cache from file.
!!  This is a collective routine that must be called by all procs in comm.
!!
!! INPUTS
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  mbsize: Cache size in megabytes.
!!    < 0 to allocate all q-points.
!!    0 has not effect.
!!    > 0 for cache with automatically computed nqpt points.
!!  qselect_dvdb(%nqpt)=0 to ignore this q-point when reading (global array)
!!  itreatq(%nqpt) = 0 if this q-point won't be treated by this CPU else > 0.
!!    Each CPU calls this routine with its own array.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_qcache_read(db, nfft, ngfft, mbsize, qselect_dvdb, itreatq, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,comm
 class(dvdb_t),intent(inout) :: db
 real(dp),intent(in) :: mbsize
!arrays
 integer,intent(in) :: ngfft(18), qselect_dvdb(db%nqpt)
 integer(i1b),intent(in) :: itreatq(db%nqpt)

!Local variables-------------------------------
!scalars
 integer :: db_iqpt, cplex,  imyp, ipc !, ii
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: v1scf(:,:,:,:)
 real(dp) :: tsec(2)

! *************************************************************************

 call timab(1801, 1, tsec)
 call wrtout(std_out, " Loading Vscf(q) from DVDB into qcache...", do_flush=.True.)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 db%qcache = qcache_new(db%nqpt, nfft, ngfft, mbsize, db%natom3, db%my_npert, db%nspden)
 db%qcache%itreatq(:) = itreatq

 do db_iqpt=1,db%nqpt
   ! Ignore points reported by the oracle
   if (qselect_dvdb(db_iqpt) == 0) cycle

   ! All procs are getting the same q-point. Exit when we reach maxnq
   ! TODO: Rewrite this part.
   !qcnt = 0
   !do ii=1,db%nqpt
   !  if (allocated(db%qcache%key(ii)%v1scf)) qcnt = qcnt + 1
   !end do
   !if (qcnt >= db%qcache%maxnq) exit

   call cwtime(cpu, wall, gflops, "start")

   ! Read all 3*natom potentials inside comm
   call db%readsym_allv1(db_iqpt, cplex, nfft, ngfft, v1scf, comm)

   ! Print progress.
   if (db_iqpt <= 10 .or. mod(db_iqpt, 50) == 0) then
     write(msg,'(2(a,i0),a)') " Reading q-point [",db_iqpt,"/",db%nqpt, "]"
     call cwtime_report(msg, cpu, wall, gflops)
   end if

   ! Transfer to cache taking into account my_npert. Note that IBZ may be distributed.
   if (db%qcache%itreatq(db_iqpt) /= 0) then
     ABI_MALLOC(db%qcache%key(db_iqpt)%v1scf, (cplex, nfft, db%nspden, db%my_npert))
     do imyp=1,db%my_npert
       ipc = db%my_pinfo(3, imyp)
       db%qcache%key(db_iqpt)%v1scf(:,:,:,imyp) = real(v1scf(:,:,:,ipc), kind=QCACHE_KIND)
     end do
   end if

   ABI_FREE(v1scf)
 end do

 call wrtout(std_out, sjoin(" Memory allocated for cache:", ftoa(db%qcache%get_mbsize(), fmt="f8.1"), " [Mb] <<< MEM"))
 call cwtime_report(" DVDB qcache IO + symmetrization", cpu_all, wall_all, gflops_all)
 call timab(1801, 2, tsec)

end subroutine dvdb_qcache_read
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_qcache_update_from_file
!! NAME
!!  dvdb_qcache_update_from_file
!!
!! FUNCTION
!!  Read selected potentials and update the internal q-cache.
!!  This is a collective routine that must be called by all procs in comm.
!!
!! INPUTS
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  ineed_qpt(%nqpt)=1 if this MPI rank requires this q-point.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_qcache_update_from_file(db, nfft, ngfft, ineed_qpt, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,comm
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18), ineed_qpt(db%nqpt)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: db_iqpt, cplex, ierr, imyp, ipc, qcnt
 real(dp) :: cpu_all, wall_all, gflops_all, mbsize, max_mbsize
 character(len=500) :: msg
!arrays
 integer :: qselect(db%nqpt)
 real(dp),allocatable :: v1scf(:,:,:,:)
 real(dp) :: tsec(2)

! *************************************************************************

 if (db%qcache%maxnq == 0) return

 ! Take the union of the q-points inside comm because we need to perform IO inside comm
 qselect = ineed_qpt
 call xmpi_sum(qselect, comm, ierr)
 qcnt = count(qselect > 0)
 if (qcnt == 0) then
   call wrtout(std_out, " All qpts in Vscf(q) already in cache. No need to perform IO. Hurrah!", do_flush=.True.)
   return
 end if

 call wrtout(std_out, sjoin(" Need to update cache. Master node will read ", itoa(qcnt), "q-points. It may be slow..."), &
             do_flush=.True.)
 call timab(1807, 1, tsec)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 if (db%qcache%make_room(ineed_qpt, msg) /= 0) then
   MSG_WARNING(msg)
 end if

 do db_iqpt=1,db%nqpt
   if (qselect(db_iqpt) == 0) cycle

   ! Read all 3*natom potentials inside comm.
   call db%readsym_allv1(db_iqpt, cplex, nfft, ngfft, v1scf, comm)

   ! Transfer to cache taking into account my_npert.
   if (ineed_qpt(db_iqpt) /= 0) then
     ABI_MALLOC(db%qcache%key(db_iqpt)%v1scf, (cplex, nfft, db%nspden, db%my_npert))
     do imyp=1,db%my_npert
       ipc = db%my_pinfo(3, imyp)
       db%qcache%key(db_iqpt)%v1scf(:,:,:,imyp) = real(v1scf(:,:,:,ipc), kind=QCACHE_KIND)
     end do
   end if

   ABI_FREE(v1scf)
 end do

 mbsize = db%ft_qcache%get_mbsize()
 call wrtout(std_out, sjoin(" Memory allocated for cache: ", ftoa(mbsize, fmt="f8.1"), " [Mb] <<< MEM"))
 call xmpi_max(mbsize, max_mbsize, comm, ierr)
 call wrtout(std_out, sjoin(" Max memory inside MPI comm: ", ftoa(max_mbsize, fmt="f8.1"), " [Mb] <<< MEM"))
 call cwtime_report(" dvdb_qcache_update_from_file", cpu_all, wall_all, gflops_all)
 call timab(1807, 2, tsec)

end subroutine dvdb_qcache_update_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/qcache_free
!! NAME
!!  qcache_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine qcache_free(qcache)

!Arguments ------------------------------------
!scalars
 class(qcache_t),intent(inout) :: qcache

!Local variables-------------------------------
!scalars
 integer :: iq

! *************************************************************************

 if (allocated(qcache%key)) then
   do iq=1,size(qcache%key)
     ABI_SFREE(qcache%key(iq)%v1scf)
   end do
   ABI_SFREE(qcache%key)
 end if
 ABI_SFREE(qcache%count_qused)
 ABI_SFREE(qcache%itreatq)
 ABI_SFREE(qcache%v1scf_3natom_qibz)

end subroutine qcache_free
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/qcache_report_stats
!! NAME
!!  qcache_report_stats
!!
!! FUNCTION
!!  Print info on q-cache states and reset counters.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine qcache_report_stats(qcache)

!Arguments ------------------------------------
!scalars
 class(qcache_t),intent(inout) :: qcache

! *************************************************************************

 if (qcache%maxnq == 0) then
   write(std_out, "(a)")" qcache deactivated with maxnq == 0"
 else if (qcache%maxnq > 0 .and. qcache%stats(1) /= 0) then
   write(std_out, "(2a)")ch10, " Qcache stats:"
   write(std_out, "(4x,a,i0)")" Total Number of calls: ", qcache%stats(1)
   write(std_out, "(4x,a,i0,2x,a,f5.1,a)") &
     " Cache hit in v1scf_3natom_qibz: ", qcache%stats(2), "(", (100.0_dp * qcache%stats(2)) / qcache%stats(1), "%)"
   write(std_out, "(4x,a,i0,2x,a,f5.1,a)") &
     " Cache hit in MPI-distributed cache: ", qcache%stats(3), "(", (100.0_dp * qcache%stats(3)) / qcache%stats(1), "%)"
   write(std_out, "(4x,a,i0,2x,a,f5.1,a)") &
     " Cache miss: ", qcache%stats(4), "(", (100.0_dp * qcache%stats(4)) / qcache%stats(1), "%)"
   write(std_out, "(a)")sjoin("     Memory allocated for cache: ", ftoa(qcache%get_mbsize(), fmt="f8.1"), " [Mb] <<< MEM")
 end if
 write(std_out, "(a)")
 qcache%stats = 0

end subroutine qcache_report_stats
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/qcache_get_mbsize
!! NAME
!!  qcache_get_mbsize
!!
!! FUNCTION
!!  Return the (allocated) size of the cache in Mb.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure real(dp) function qcache_get_mbsize(qcache) result(mbsize)

!Arguments ------------------------------------
!scalars
 class(qcache_t),intent(in) :: qcache

!Local variables-------------------------------
!scalars
 integer :: qcnt, iq

! *************************************************************************

 mbsize = zero; if (.not. allocated(qcache%key)) return
 ! Compute cache size.
 qcnt = 0
 do iq=1,size(qcache%key)
   if (allocated(qcache%key(iq)%v1scf)) qcnt = qcnt + 1
 end do
 mbsize = qcache%onepot_mb * qcnt

end function qcache_get_mbsize
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/qcache_qcache_make_room
!! NAME
!!  qcache_qcache_make_room
!!
!! FUNCTION
!!  Deallocate entries in cache taking into account the list of q-points required
!!  in the next iteration. Deallocate entries only if the cache size required to
!!  treat all the points specified by ineed_qpt would become greater than %max_bsize
!!  Return ierr /= 0 and msg if max_bsize constraint cannot be respected.
!!
!! INPUTS
!!  ineed_qpt(%nqpt) = 0 if this q-points is not required.
!!
!! OUTPUTS
!!  msg=Warning message if ierr /= 0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function qcache_make_room(qcache, ineed_qpt, msg) result(ierr)

!Arguments ------------------------------------
!scalars
 class(qcache_t),intent(inout) :: qcache
 character(len=*),intent(out) :: msg
!arrays
 integer,intent(in) :: ineed_qpt(qcache%nqpt)

!Local variables-------------------------------
!scalars
 integer :: iq, nq_remove, qcnt, count_qnew
 real(dp) :: mbsize_now, mbsize_extra

! *************************************************************************

 ierr = 0; msg = ""
 if (qcache%max_mbsize < zero) return
 mbsize_now = qcache%get_mbsize()

 ! Count the number of q-points that are not in cache and the extra memory required to allocate everything.
 count_qnew = 0
 do iq=1,qcache%nqpt
   if (ineed_qpt(iq) > 0 .and. .not. allocated(qcache%key(iq)%v1scf)) count_qnew = count_qnew + 1
 end do
 mbsize_extra = mbsize_now + count_qnew * qcache%onepot_mb

 if (mbsize_extra > qcache%max_mbsize) then
   ! Try to deallocate nq_remove q-points provided they are not needed.
   nq_remove = nint((mbsize_extra - qcache%max_mbsize) / qcache%onepot_mb)
   qcnt = 0
   do iq=1,qcache%nqpt
     if (ineed_qpt(iq) == 0 .and. allocated(qcache%key(iq)%v1scf)) then
       ABI_FREE(qcache%key(iq)%v1scf)
       qcnt = qcnt + 1
       if (qcnt == nq_remove) exit
     end if
   end do
   if (iq == qcache%nqpt + 1) then
     ierr = 1
     msg = "Couldn't decrease cache size below input limit. Continuing anyway but we may go out of memory!"
   end if
 end if

end function qcache_make_room
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/v1phq_complete
!! NAME
!! v1phq_complete
!!
!! FUNCTION
!!  Use the symmetries of the little group of the q-point to reconstruct
!!  the first order potentials starting from an initial irreducible set.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  qpt(3)=q-point in reduced coordinates.
!!  ngfft(18)=Info of FFT grid.
!!  cplex=1 if real potentials (qpt==gamma), 2 if complex
!!  nfft=(effective) number of FFT grid points (for this proc).
!!  nspden=number of spin-density components
!!  nsppol=Number of independent spin polarizations
!!  mpi_enreg=information about MPI parallelization
!!  symv1=If 1, the new potentials are symmetrized using the set of symmetries that leaves the
!!    perturbation invariant.
!!
!! SIDE EFFECTS
!!  pflag(3,natom)= For each atomic perturbation:
!!     0 if pert is not available. 1 if pert is available. 2 if pert has been reconstructed by symmetry.
!!     Initialized by the caller. Changed in output.
!!  v1scf(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials.
!!    in input: filled with the irreducible potentials (corresponding pflag set to 1)
!!    output: Contains full set of perturbations.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1phq_complete(cryst,qpt,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,symv1,pflag,v1scf)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,nsppol
 integer,intent(in) :: symv1
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(inout) :: pflag(3, cryst%natom)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(inout) :: v1scf(cplex*nfft,nspden,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: syuse0=0,rfmeth2=2,tim_fourdp0=0
 integer :: idir,ipert,tsign,isym_eq,itirev_eq,ipert_eq !,itirev
 integer :: pcase,trev_q,idir_eq,pcase_eq,ispden,cnt
 integer :: i1,i2,i3,id1,id2,id3,n1,n2,n3,ind1,ind2,j1,j2,j3,l1,l2,l3,k1,k2,k3,nfftot
 real(dp) :: arg
 logical :: has_phase
 logical,parameter :: debug=.False.
 character(len=500) :: msg
 integer,save :: enough=0
!arrays
 integer :: symrel_eq(3,3),symrec_eq(3,3),g0_qpt(3),l0(3),tsm1g(3) !symm(3,3),
 integer :: symq(4,2,cryst%nsym)
 real(dp) :: phnon1(2),tnon(3)
 real(dp),allocatable :: workg(:,:), workg_eq(:,:),v1g(:,:,:)

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")
 id1 = n1/2+2; id2 = n2/2+2; id3 = n3/2+2

 ABI_MALLOC(v1g, (2,nfft,nspden))
 ABI_MALLOC(workg_eq, (2, nfft))
 ABI_MALLOC(workg, (2, nfft))

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,trev_q,prtvol=0)

pcase_loop: &
 do pcase=1,3*cryst%natom
   idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
   if (pflag(idir, ipert) /= 0) cycle ! This pcase is available

   ! Find symmetry which links to the perturbation requested (pcase)
   call find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt, allow_g0=.true.)
   !if (isym_eq == -1) then
   !  call find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt, allow_g0=.true.)
   !end if

   if (isym_eq == -1) then
     if (debug) write(std_out,*)"Cannot find isym eq for idir, ipert:", idir,ipert
     cycle pcase_loop
   end if

   ! set flag since we will reconstruct pcase from isym_eq.
   pflag(idir, ipert) = 2

   symrel_eq = cryst%symrel(:,:,isym_eq)
   symrec_eq = cryst%symrec(:,:,isym_eq)
   tsign = 3-2*itirev_eq

   ! Phase due to L0 + R^{-1}tau
   l0 = cryst%indsym(1:3,isym_eq,ipert)

   tnon = l0 + matmul(transpose(symrec_eq), cryst%tnons(:,isym_eq))
   has_phase = any(abs(tnon) > tol12)
   ! FIXME
   !ABI_CHECK(.not. has_phase, "has phase must be tested")
   if (has_phase) then
     enough = enough + 1
     if (enough == 1) MSG_WARNING("has phase must be tested")
   end if

   workg = zero

   ! Reconstruct DFPT potential. Final results stored in v1g.
   if (debug) write(std_out,*)"Reconstructing idir:", idir, ", ipert:", ipert
   v1g = zero; cnt = 0
   do idir_eq=1,3
     if (symrec_eq(idir, idir_eq) == 0) cycle
     cnt = cnt + 1
     pcase_eq = idir_eq + (ipert_eq-1)*3
     if (debug) write(std_out,*) "idir_eq: ", idir_eq, ", ipert_eq: ", ipert_eq, ", tsign: ", tsign

     if (pflag(idir_eq, ipert_eq) == 0) then
       write(msg, *)"pflag for idir_eq, ipert_eq", idir_eq, ipert_eq, "cannot be zero"
       MSG_ERROR(msg)
     end if

     !if (pflag(idir_eq, ipert_eq) == 0) then
     !  write(msg, *)"pflag for idir_eq, ipert_eq", idir_eq, ipert_eq, "cannot be zero"
     !  MSG_ERROR(msg)
     !end if

     do ispden=1,nspden
       ! Get symmetric perturbation in G-space in workg_eq array.
       call fourdp(cplex,workg_eq,v1scf(:,ispden,pcase_eq),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp0)
       !call zerosym(workg_eq,cplex,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)

       !call rotate_fqg(itirev_eq,symrec_eq,qpt,tnon,ngfft,nfft,nspden,workg_eq,workg)
       ind1=0
       do i3=1,n3
         ! Get location of G vector (grid point) centered at 0 0 0
         l3 = i3-(i3/id3)*n3-1
         do i2=1,n2
           l2 = i2-(i2/id2)*n2-1
           do i1=1,n1
             ind1=ind1+1

             l1 = i1-(i1/id1)*n1-1

             ! Get rotated G vector Gj for each symmetry element
             ! -- here we use the TRANSPOSE of symrel_eq; assuming symrel_eq expresses
             ! the rotation in real space, the transpose is then appropriate
             ! for G space symmetrization (p. 1172d,e of notes, 2 June 1995).
             j1 = tsign * (symrel_eq(1,1)*l1+symrel_eq(2,1)*l2+symrel_eq(3,1)*l3)
             j2 = tsign * (symrel_eq(1,2)*l1+symrel_eq(2,2)*l2+symrel_eq(3,2)*l3)
             j3 = tsign * (symrel_eq(1,3)*l1+symrel_eq(2,3)*l2+symrel_eq(3,3)*l3)

             ! FIXME :TO BE CLARIFIED:
             ! We are not working on the G-sphere thus SG may be outside
             ! of the box. This check is not done in irrzg!!!
             if ( (j1 > n1/2 .or. j1 < -(n1-1)/2) .or. &
                  (j2 > n2/2 .or. j1 < -(n2-1)/2) .or. &
                  (j3 > n3/2 .or. j3 < -(n3-1)/2) ) then
               !write(std_out,*)"got it"
               workg(:, ind1) = zero; cycle
             end if

             tsm1g = [j1,j2,j3] ! +- S^{-1} G

             ! Map into [0,n-1] and then add 1 for array index in [1,n]
             k1=1+mod(n1+mod(j1,n1),n1)
             k2=1+mod(n2+mod(j2,n2),n2)
             k3=1+mod(n3+mod(j3,n3),n3)

             ! Get linear index of rotated point Gj
             ind2 = k1+n1*((k2-1)+n2*(k3-1))

             if (has_phase) then
               ! compute exp(-2*Pi*I*G dot tau) using original G
               ! NB: this phase is same as that in irrzg and phnons1, and corresponds
               ! to complex conjugate of phase from G to Gj;
               ! we use it immediately below, to go _to_ workg_eq(ind1)
               arg = two_pi * dot_product(qpt + tsm1g, tnon)
               phnon1(1) = cos(arg); phnon1(2) = -sin(arg)

               ! rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
               workg(1, ind1) = phnon1(1) * workg_eq(1, ind2) - phnon1(2) * workg_eq(2, ind2)
               workg(2, ind1) = phnon1(1) * workg_eq(2, ind2) + phnon1(2) * workg_eq(1, ind2)
             else
               workg(1, ind1) = workg_eq(1, ind2)
               workg(2, ind1) = workg_eq(2, ind2)
             end if

             ! Take complex conjugate if time-reversal is used.
             if (tsign == -1) workg(2, ind1) = -workg(2, ind1)
           end do
         end do
       end do

       v1g(:,:,ispden) = v1g(:,:,ispden) + workg * symrec_eq(idir, idir_eq)
     end do ! ispden
   end do ! idir_eq
   !if (debug) write(std_out,*)"Used ",cnt," equivalent perturbations"

   ! Get potential in real space (results in v1scf)
   do ispden=1,nspden
     !call zerosym(v1g(:,:,ispden),cplex,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     call fourdp(cplex,v1g(:,:,ispden),v1scf(:,ispden,pcase),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp0)

     ! IS(q) = q + G0
     ! we want q so we have to multiply by exp(iG0r) in real space.
     if (any(g0_qpt /= 0)) then
       ABI_CHECK(cplex==2, "cplex == 1")
       if (debug) write(std_out,*)"Found not zero g0_qpt", g0_qpt ! for idir: ", idir, ", ipert: ", ipert
       call times_eigr(g0_qpt, ngfft, nfft, 1, v1scf(:,ispden,pcase))
     end if
   end do

   if (symv1==1) then
     if (debug) write(std_out,*)"Calling v1phq_symmetrize"
     call v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1scf(:,:,pcase))
   end if
 end do pcase_loop

 ABI_FREE(v1g)
 ABI_FREE(workg)
 ABI_FREE(workg_eq)

 ! Handle possible error.
 if (any(pflag == 0)) then
   write(std_out,"(2a)")&
     "The following perturbations cannot be recostructed by symmetry for q-point: ",trim(ktoa(qpt))
   do ipert=1,cryst%natom
     do idir=1,3
       if (pflag(idir, ipert) == 0) write(std_out,"(2(a,i0))")"idir= ",idir,", ipert= ",ipert
     end do
   end do
   write(msg,"(5a)")&
     "Cannot recostruct all 3*natom atomic perturbations from file",ch10,&
     "This usually happens when the DVDB does not contain all the independent perturbations for this q-point",ch10,&
     "See above message for further information."
   MSG_ERROR(msg)
 end if

end subroutine v1phq_complete
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/find_symeq
!! NAME
!! find_symeq
!!
!! FUNCTION
!!  Find symmetry which links to the perturbation specified by (idir, ipert)
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  idir=Direction of the perturbation
!!  ipert=Perturbation type.
!!  symq(4,2,nsym)=Table produced by littlegroup_q
!!  pflag(3,natom)= For each atomic perturbation:
!!     0 if pert is not available. 1 if pert is available. 2 if pert can been reconstructed by symmetry.
!!
!! OUTPUT
!!  ipert_eq
!!  isym_eq
!!  itirev_eq
!!  g0_qpt(3)
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt, allow_g0)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert
 integer,intent(out) :: ipert_eq,isym_eq,itirev_eq
 type(crystal_t),intent(in) :: cryst
 logical,optional :: allow_g0
!arrays
 integer,intent(in) :: symq(4,2,cryst%nsym),pflag(3,cryst%natom)
 integer,intent(out) :: g0_qpt(3)

!Local variables-------------------------------
!scalars
 logical :: do_allow_g0
 integer :: isym,idir_eq,ip,itirev

! *************************************************************************

 isym_eq = -1; ipert_eq = -1
 do_allow_g0 = .true.; if (present(allow_g0)) do_allow_g0 = allow_g0

symloop: &
 do itirev=1,2
   itirev_eq = itirev
   do isym=1,cryst%nsym

     ! Check that isym preserves the q-point
     !if (symq(4,itirev,isym) /= 1 .or. any(symq(1:3,itirev,isym) /= 0)) cycle
     if (symq(4,itirev,isym) /= 1) cycle ! .or. any(symq(1:3,itirev,isym) /= 0)) cycle
     if (any(symq(1:3,itirev,isym) /= 0) .and. .not.do_allow_g0) cycle
     g0_qpt = symq(1:3,itirev,isym)

     do ip=1,cryst%natom
       !if (.not. cryst%indsym(4,isym,ip) == ipert) cycle
       if (.not. cryst%indsym(4,isym,ipert) == ip) cycle
       isym_eq = isym; ipert_eq = ip
       do idir_eq=1,3
         if (idir_eq == idir .and. ip == ipert .and. cryst%symrec(idir,idir_eq,isym) /= 0) isym_eq = -1
         if (cryst%symrec(idir,idir_eq,isym) /= 0 .and. pflag(idir_eq, ip) == 0) then
         !if (cryst%symrel(idir,idir_eq,isym) /= 0 .and. pflag(idir_eq, ip) == 0) then
           !if (idir_eq == idir .and. ip == ipert) cycle
           isym_eq = -1
         end if
       end do
       if (isym_eq /= -1) exit symloop
     end do
   end do
 end do symloop

 if (isym_eq == -1) then
   ipert_eq = -1; itirev_eq = -1
 end if

end subroutine find_symeq
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/v1phq_rotate
!! NAME
!! v1phq_rotate
!!
!! FUNCTION
!!  Reconstruct the DFPT potential for a q-point in the BZ starting from its symmetrical image in the IBZ.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  qpt_ibz(3)=q-point in the IBZ in reduced coordinates.
!!  ngfft=array of dimensions for different FFT grids
!!  isym, itimrev, g0q:
!!    qpt_bz = I(itimrev) S(isym) q_ibz + g0q
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  nspden=number of spin-density components
!!  nsppol=Number of independent spin polarizations
!!  mpi_enreg=information about MPI parallelization
!!  v1r_qibz(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space
!!    for the irreducible q-point `qpt_ibz`
!!  comm: MPI communicator to distribute natom3 * nspden FFT calls
!!
!! OUTPUT
!!  v1r_qbz(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space for the q-point in the BZ
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1phq_rotate(cryst,qpt_ibz,isym,itimrev,g0q,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r_qibz,v1r_qbz,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isym,itimrev,cplex,nfft,nspden,nsppol,comm
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: g0q(3),ngfft(18)
 real(dp),intent(in) :: qpt_ibz(3)
 real(dp),intent(inout) :: v1r_qibz(cplex*nfft,nspden,3*cryst%natom)
 real(dp) ABI_ASYNC, intent(out) :: v1r_qbz(cplex*nfft,nspden,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0
 integer,save :: enough=0
 integer :: natom3,mu,ispden,idir,ipert,idir_eq,ipert_eq,mu_eq,cnt,tsign,my_rank,nproc,ierr,root
!arrays
 integer :: symrec_eq(3,3),sm1(3,3),l0(3) !g0_qpt(3), symrel_eq(3,3),
 real(dp) :: tnon(3), tsec(2)
 real(dp) ABI_ASYNC, allocatable :: v1g_qibz(:,:,:),workg(:,:),v1g_mu(:,:)
 integer :: requests(nspden, 3*cryst%natom), requests_v1r_qbz(3*cryst%natom)
 logical :: requests_v1g_qibz_done(nspden, 3*cryst%natom)

! *************************************************************************

 ! Keep track of total time spent.
 call timab(1804, 1, tsec)

 ABI_UNUSED(nsppol)
 ABI_CHECK(cplex == 2, "cplex != 2")

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom3 = 3 * cryst%natom; tsign = 3-2*itimrev

 ! Compute IBZ potentials in G-space. Results stored in v1g_qibz(G)
 ABI_MALLOC(v1g_qibz, (2*nfft, nspden, natom3))
 requests_v1g_qibz_done = .False.
 cnt = 0
 do mu=1,natom3
   do ispden=1,nspden
     cnt = cnt + 1; root = mod(cnt, nproc)
     if (root == my_rank) then ! Non-blocking
       call fourdp(cplex,v1g_qibz(:,ispden,mu),v1r_qibz(:,ispden,mu),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp0)
     end if
     call xmpi_ibcast(v1g_qibz(:,ispden,mu), root, comm, requests(ispden, mu), ierr)
   end do
 end do

 ABI_MALLOC(workg, (2*nfft,nspden))
 ABI_MALLOC(v1g_mu, (2*nfft,nspden))

 ! For each perturbation:
 ! FIXME: This is wrong if natom > 1
 symrec_eq = cryst%symrec(:,:,isym)
 call mati3inv(symrec_eq, sm1); sm1 = transpose(sm1)

 !v1r_qbz = zero
 do mu=1,natom3
   root = mod(mu, nproc)
   ! MPI parallelism.
   if (root == my_rank) then
     idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1

     ! Phase due to L0 + R^{-1}tau
     l0 = cryst%indsym(1:3,isym,ipert)
     tnon = l0 + matmul(transpose(symrec_eq), cryst%tnons(:,isym))
     ! FIXME
     !ABI_CHECK(all(abs(tnon) < tol12), "tnon!")
     if (.not.all(abs(tnon) < tol12)) then
       enough = enough + 1
       if (enough == 1) MSG_WARNING("tnon must be tested!")
     end if

     ipert_eq = cryst%indsym(4,isym,ipert)

     v1g_mu = zero; cnt = 0
     do idir_eq=1,3
       if (symrec_eq(idir, idir_eq) == 0) cycle
       mu_eq = idir_eq + (ipert_eq-1)*3
       cnt = cnt + 1

       ! Wait for request before operating on v1g_qibz
       if (.not. all(requests_v1g_qibz_done(:, mu_eq))) then
         do ispden=1,nspden
           call xmpi_wait(requests(ispden, mu_eq), ierr)
           requests_v1g_qibz_done(ispden, mu_eq) = .True.
         end do
       end if

       ! Rotate in G-space and accumulate in workg
       call rotate_fqg(itimrev,sm1,qpt_ibz,tnon,ngfft,nfft,nspden,v1g_qibz(:,:,mu_eq),workg)
       v1g_mu = v1g_mu + workg * symrec_eq(idir, idir_eq)
     end do ! idir_eq

     ABI_CHECK(cnt /= 0, "cnt should not be zero!")

     ! Transform to real space and take into account a possible shift. Results are stored in v1r_qbz.
     do ispden=1,nspden
       call fourdp(cplex,v1g_mu(:,ispden),v1r_qbz(:,ispden,mu),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp0)
       call times_eigr(-g0q, ngfft, nfft, 1, v1r_qbz(:,ispden,mu))
       !call times_eigr(tsign * g0q, ngfft, nfft, 1, v1r_qbz(:,ispden,mu))
     end do

     !call v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r)
   end if ! root == myrank

   call xmpi_ibcast(v1r_qbz(:,:,mu), root, comm, requests_v1r_qbz(mu), ierr)
 end do ! mu

 ! Relase all requests
 call xmpi_waitall(requests, ierr)
 call xmpi_waitall(requests_v1r_qbz, ierr)

 ABI_FREE(workg)
 ABI_FREE(v1g_mu)
 ABI_FREE(v1g_qibz)

 call timab(1804, 2, tsec)

end subroutine v1phq_rotate
!!***

!!****f* m_dvdb/v1phq_symmetrize
!! NAME
!! v1phq_symmetrize
!!
!! FUNCTION
!!  Enforce spatial-symmetry on the DFPT potential.
!!
!! INPUTS
!! cryst<crystal_t>=crystal structure parameters
!! idir=Direction of the perturbation
!! ipert=Perturbation type.
!! symq(4,2,nsym)= Table computed by littlegroup_q.
!!   three first numbers define the G vector;
!!   fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!   second index is one without time-reversal symmetry, two with time-reversal symmetry
!! ngfft=array of dimensions for different FFT grids
!! cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!! nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!! nspden=number of spin-density components
!! mpi_enreg=information about MPI parallelization
!!
!! SIDE EFFECTS
!!  v1r(cplex*nfft,nspden)=Array with first order potentials in real space. Symmetrized in output.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,cplex,nfft,nspden,nsppol
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: symq(4,2,cryst%nsym),ngfft(18)
 real(dp),intent(inout) :: v1r(cplex*nfft,nspden)

!Local variables-------------------------------
 integer,parameter :: syuse0=0,rfmeth2=2,iscf1=1
 integer :: nsym1,nfftot
!arrays
 integer :: symafm1(cryst%nsym),symrel1(3,3,cryst%nsym),symrc1(3,3,cryst%nsym)
 integer,allocatable :: irrzon1(:,:,:),indsy1(:,:,:)
 real(dp) :: tnons1(3,cryst%nsym)
 real(dp),allocatable :: phnons1(:,:,:),v1g(:,:)

! *************************************************************************

 if (cryst%nsym == 1) return

 nfftot = product(ngfft(1:3))
 ABI_CHECK(nfft == nfftot, "MPI-FFT not coded")

 ! Symmetrize (copied from dfpt_looppert)
 ! Determines the set of symmetries that leaves the perturbation invariant.
 call littlegroup_pert(cryst%gprimd,idir,cryst%indsym,dev_null,ipert,cryst%natom,cryst%nsym,nsym1,rfmeth2,&
   cryst%symafm,symafm1,symq,cryst%symrec,cryst%symrel,symrel1,syuse0,cryst%tnons,tnons1,unit=dev_null)

 ! Set up corresponding symmetry data
 ABI_MALLOC(irrzon1, (nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4)))
 ABI_MALLOC(phnons1, (2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4)))
 ABI_MALLOC(indsy1,(4,nsym1,cryst%natom))

 call setsym(indsy1,irrzon1,iscf1,cryst%natom,nfft,ngfft,nspden,nsppol,&
   nsym1,phnons1,symafm1,symrc1,symrel1,tnons1,cryst%typat,cryst%xred)

 !if (psps%usepaw==1) then
 !  ! Allocate/initialize only zarot in pawang1 datastructure
 !  call pawang_init(pawang1,0,0,pawang%l_max-1,0,0,nsym1,0,0,0,0)
 !  call setsym_ylm(gprimd,pawang1%l_max-1,pawang1%nsym,0,rprimd,symrc1,pawang1%zarot)
 !end if

 ! FIXME Be careful here because symrhg was written for densities!
 ABI_CHECK(nsppol == 1 .and. nspden == 1, "symrhg was written for densities, not for potentials")

 ABI_MALLOC(v1g, (2,nfft))
 call symrhg(cplex,cryst%gprimd,irrzon1,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym1,&
    phnons1,v1g,v1r,cryst%rprimd,symafm1,symrel1,tnons1)

 ABI_FREE(irrzon1)
 ABI_FREE(phnons1)
 ABI_FREE(indsy1)
 ABI_FREE(v1g)

end subroutine v1phq_symmetrize
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/rotate_fqg
!! NAME
!!  rotate_fqg
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_fqg(itirev, symm, qpt, tnon, ngfft, nfft, nspden, infg, outfg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itirev,nfft,nspden
!arrays
 integer,intent(in) :: symm(3,3),ngfft(18)
 real(dp),intent(in) :: qpt(3),tnon(3)
 real(dp),intent(in) :: infg(2,nfft,nspden)
 real(dp),intent(out) :: outfg(2,nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,n1,n2,n3,ind1,ind2,j1,j2,j3,l1,l2,l3,k1,k2,k3,nfftot,isp,tsign
 real(dp) :: arg
 logical :: has_phase
!arrays
 integer :: tsg(3)
 real(dp) :: phnon1(2), tsec(2)

! *************************************************************************

 ! Keep track of total time spent.
 call timab(1803, 1, tsec)

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")
 id1 = n1/2+2; id2 = n2/2+2; id3=n3/2+2

 ABI_CHECK(any(itirev == [1,2]), "Wrong itirev")
 tsign = 3-2*itirev; has_phase = any(abs(tnon) > tol12)

 !outfg = zero

 do isp=1,nspden
   ind1 = 0
   do i3=1,n3
     ! Get location of G vector (grid point) centered at 0 0 0
     l3 = i3-(i3/id3)*n3-1
     do i2=1,n2
       l2 = i2-(i2/id2)*n2-1
       do i1=1,n1
         ind1 = ind1 + 1
         !ind1 = 1 + i1 + (i2-1)*n1 + (i3-1)*n1*n2
         !if (mod(ind1, nprocs) /= my_rank) cycle

         l1 = i1-(i1/id1)*n1-1

         ! Get rotated G vector. IS(G)
         j1 = tsign * (symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3)
         j2 = tsign * (symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3)
         j3 = tsign * (symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3)

         ! FIXME :TO BE CLARIFIED:
         ! We are not working on the G-sphere thus SG may be outside
         ! of the box. This check is not done in irrzg!!!
         if ( (j1 > n1/2 .or. j1 < -(n1-1)/2) .or. &
              (j2 > n2/2 .or. j1 < -(n2-1)/2) .or. &
              (j3 > n3/2 .or. j3 < -(n3-1)/2) ) then
           !write(std_out,*)"outsize box!"
           outfg(:,ind1,isp) = zero
           cycle
         end if

         tsg = [j1,j2,j3] ! +- S^{-1} G

         ! Map into [0,n-1] and then add 1 for array index in [1, n]
         k1=1+mod(n1+mod(j1,n1),n1)
         k2=1+mod(n2+mod(j2,n2),n2)
         k3=1+mod(n3+mod(j3,n3),n3)

         ! Get linear index of rotated point Gj
         ind2 = k1+n1*((k2-1)+n2*(k3-1))

         ! TOD: Here I believe there are lots of cache misses, should perform low-level profiling
         ! OMP perhaps can accelerate this part but mind false sharing...
         if (has_phase) then
           ! compute exp(-2*Pi*I*G dot tau) using original G
           arg = two_pi * dot_product(qpt + tsg, tnon)
           phnon1(1) = cos(arg); phnon1(2) =-sin(arg)

           ! rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
           outfg(1, ind1, isp) = phnon1(1) * infg(1, ind2, isp) - phnon1(2) * infg(2, ind2, isp)
           outfg(2, ind1, isp) = phnon1(1) * infg(2, ind2, isp) + phnon1(2) * infg(1, ind2, isp)
         else
           outfg(1, ind1, isp) = infg(1, ind2, isp)
           outfg(2, ind1, isp) = infg(2, ind2, isp)
         end if

         ! Take complex conjugate if time-reversal is used.
         if (tsign == -1) outfg(2, ind1, isp) = -outfg(2, ind1, isp)
       end do
     end do
   end do
 end do ! isp

 !call xmpi_sum(comm, outfg, ierr)

 call timab(1803, 2, tsec)

end subroutine rotate_fqg
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_ftinterp_setup
!! NAME
!!  dvdb_ftinterp_setup
!!
!! FUNCTION
!!  Precompute the v1scf_rpt array with with the DFPT potential in the supercell
!!  required for the Fourier interpolation
!!  This is a collective routine that should be called by all procs inside db%comm.
!!
!! \begin{equation}
!! 	\label{eq:dfpt_pot_realspace}
!!     W_{\kappa\alpha}(\rr,\RR) = \dfrac{1}{N_\qq} \sum_\qq e^{-i\qq\cdot(\RR - \rr)}\,
!!     \partial_{\kappa\alpha\qq}{v^{\text{scf}}}(\rr),
!! \end{equation}
!!
!! INPUTS
!!  ngqpt(3)=Divisions of the ab-initio q-mesh.
!!  qrefine(3)=Defines intial coarse q-mesh if qrefine /= [1, 1, 1]
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  comm_rpt = MPI communicator used to distribute R-lattice points.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_ftinterp_setup(db, ngqpt, qrefine, nqshift, qshift, nfft, ngfft, method, comm_rpt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,comm_rpt, method
 class(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngqpt(3), qrefine(3), ngfft(18)
 real(dp),intent(in) :: qshift(3,nqshift)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0, qptopt1 = 1
 integer :: iq_ibz,nqibz,iq_bz,nqbz !, timerev_q
 integer :: ii,jj,cplex_qibz,ispden,imyp,irpt,idir,ipert,ipc
 integer :: iqst, itimrev, isym
 integer :: ifft, ierr,  my_rstart, my_rstop, iatom
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all
 logical :: isirr_q
 character(len=500) :: msg
!arrays
 integer :: g0q(3)
 !integer :: symq(4,2,db%cryst%nsym)
 integer,allocatable :: indqq(:,:), iperm(:), nqsts(:), iqs_dvdb(:), all_cell(:,:)
 real(dp) :: qpt_bz(3)
 real(dp),allocatable :: qibz(:,:), qbz(:,:), emiqr(:,:), all_rpt(:,:), all_wghatm(:,:,:)
 real(dp),allocatable :: v1r_qibz(:,:,:,:), v1r_qbz(:,:,:,:), v1r_lr(:,:,:)

! *************************************************************************

 ABI_CHECK(all(qrefine == 1), "qrefine not coded")

 ! Set communicator for R-point parallelism.
 ! Note that client code is responsible for calling the interpolation routine dvdb_get_ftqbz (R -> q)
 ! with all procs inside comm_rpt to avoid MPI deadlocks.
 db%comm_rpt = comm_rpt; db%nprocs_rpt = xmpi_comm_size(db%comm_rpt); db%me_rpt = xmpi_comm_rank(db%comm_rpt)

 if (db%add_lr == 4) then
   call wrtout(std_out, " Skipping construction of W(R,r) because add_lr = 4")
   return
 end if

 ABI_CHECK(.not. allocated(db%v1scf_rpt), " v1scf_rpt is already allocated!")
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 call wrtout(std_out, sjoin(ch10, "Building W(R,r) using ngqpt: ", ltoa(ngqpt), &
   ", with nprocs_rpt:", itoa(db%nprocs_rpt)), do_flush=.True.)
 call wrtout(std_out, " Q-mesh qshift:")
 do ii=1,nqshift
   call wrtout(std_out, ltoa(qshift(:, ii)))
 end do
 call wrtout(std_out, "")
 call wrtout(std_out, " It may take some time depending on the number of MPI procs, ngqpt and nfft points.")
 call wrtout(std_out, " Use boxcutmin < 2.0 (> 1.1) to decrease nfft, reduce memory requirements and speedup the calculation.")

 call prepare_ftinterp(db, method, ngqpt, qptopt1, nqshift, qshift, &
     qibz, qbz, indqq, iperm, nqsts, iqs_dvdb, all_rpt, all_wghatm, db%comm)

 nqibz = size(qibz, dim=2); nqbz = size(qbz, dim=2); db%nrtot = size(all_rpt, dim=2)

 ! Distribute R-points inside comm_rpt.
 call xmpi_split_work(db%nrtot, db%comm_rpt, my_rstart, my_rstop)

 ! Select my_rpoints.
 db%my_nrpt = my_rstop - my_rstart + 1
 ABI_CHECK(db%my_nrpt /= 0, "my_nrpt == 0!")
 ABI_MALLOC(db%my_rpt, (3, db%my_nrpt))
 db%my_rpt = all_rpt(:, my_rstart:my_rstop)
 ABI_MALLOC(db%my_irpt2tot, (db%my_nrpt))
 do irpt=1,db%my_nrpt
   db%my_irpt2tot(irpt) = my_rstart + (irpt - 1)
 end do

 ! Copy weights for the atoms treated by this proc.
 ii = minval(db%my_pinfo(2,:)); jj = maxval(db%my_pinfo(2,:))
 ABI_MALLOC(db%my_wratm, (db%my_nrpt, ii:jj))
 db%my_wratm = one
 if (method == 1) then
   do iatom=1,db%cryst%natom
     if (iatom >= ii .and. iatom <= jj) db%my_wratm(:, iatom) = all_wghatm(iatom, iatom, my_rstart:my_rstop)
   end do
 end if

 write(std_out, "(a, i0)")" Using method for integration weights: ", method
 write(std_out, "(a, i0)")" Total number of R-points in real-space big box: ", db%nrtot
 write(std_out, "(a, i0)")" Number of R-points treated by this MPI rank: ", db%my_nrpt
 write(std_out, "(a, 3(i0, 1x))")" ngfft: ", ngfft(1:3)
 write(std_out, "(a, i0)")" dvdb_add_lr: ", db%add_lr
 ! Allocate potential in the supercell. Memory is distributed over my_nrpt and my_npert
 call wrtout(std_out, sjoin(" Memory required for W(R,r): ", &
    ftoa(two * db%my_nrpt * nfft * db%nspden * db%my_npert * dp * b2Mb, fmt="f8.1"), "[Mb] <<< MEM"))

 ABI_SFREE(all_cell)
 ABI_SFREE(all_wghatm)
 ABI_FREE(all_rpt)

 ABI_MALLOC(emiqr, (2, db%my_nrpt))
 ABI_MALLOC(v1r_qbz, (2, nfft, db%nspden, db%natom3))
 ABI_CALLOC(v1r_lr, (2, nfft, db%my_npert))
 ABI_MALLOC_OR_DIE(db%v1scf_rpt, (2, db%my_nrpt, nfft, db%nspden, db%my_npert), ierr)
 db%v1scf_rpt = zero

 ! TODO: Parallelize this part over q-points using comm_rpt. For the time being only pert parallelism.
 iqst = 0
 do iq_ibz=1,nqibz
   call cwtime(cpu, wall, gflops, "start")
   !
   ! Here all procs get all potentials for this IBZ q-point on the real-space FFT mesh.
   ! This call allocates v1r_qibz(cplex_qibz, nfft, nspden, 3*natom)
   ! Note that here we need all 3*natom perturbations because of v1phq_rotate
   call db%readsym_allv1(iqs_dvdb(iq_ibz), cplex_qibz, nfft, ngfft, v1r_qibz, db%comm)

   ! Reconstruct by symmetry the potentials for the star of this q-point,
   ! perform slow FT and accumulate in v1scf_rpt. Be careful with the gamma point.
   do ii=1,nqsts(iq_ibz)
     iqst = iqst + 1
     !if (mod(ii, nproc) /= my_rank) cycle ! MPI parallelism.
     iq_bz = iperm(iqst)
     ABI_CHECK(iq_ibz == indqq(iq_bz,1), "iq_ibz !/ indqq(1)")
     qpt_bz = qbz(:, iq_bz)
     ! IS(q_ibz) + g0q = q_bz
     isym = indqq(iq_bz, 2); itimrev = indqq(iq_bz, 6) + 1; g0q = indqq(iq_bz, 3:5)
     isirr_q = (isym == 1 .and. itimrev == 1 .and. all(g0q == 0))
     !write(std_out, *)"qbz", trim(ktoa(qpt_bz)), " --> qibz ", trim(ktoa(qibz(:,iq_ibz)))
     !write(std_out, *)"via isym, itimrev, g0q:", isym, itimrev, g0q

     ! Compute long-range part of the coupling potential at qpt_bz.
     if (db%add_lr /= 0) then
       do imyp=1,db%my_npert
         idir = db%my_pinfo(1, imyp); ipert = db%my_pinfo(2, imyp)
         call db%get_v1r_long_range(qpt_bz, idir, ipert, nfft, ngfft, v1r_lr(:,:,imyp))
       end do
     end if

     if (cplex_qibz == 1) then
       ! Gamma point.
       ABI_CHECK(nqsts(iq_ibz) == 1, "cplex_qibz == 1 and nq nqst /= 1 (should be gamma)")
       ABI_CHECK(all(g0q == 0), "gamma point with g0q /= 0")

       if (db%add_lr /= 0) then
         ! Substract the long-range part of the potential.
         do imyp=1,db%my_npert
           ipc = db%my_pinfo(3, imyp)
           do ispden=1,db%nspden
             v1r_qibz(1, :, ispden, ipc) = v1r_qibz(1, :, ispden, ipc) - v1r_lr(1, :, imyp)
           end do
         end do
       end if

       ! Slow FT.
       do imyp=1,db%my_npert
         ipc = db%my_pinfo(3, imyp)
         do ispden=1,db%nspden
           do ifft=1,nfft
             do irpt=1,db%my_nrpt
               db%v1scf_rpt(1, irpt, ifft, ispden, imyp) = db%v1scf_rpt(1, irpt, ifft, ispden, imyp) + &
                                                               v1r_qibz(1, ifft, ispden, ipc)
             end do
           end do
         end do
       end do

     else
       ! q /= Gamma. Get the periodic part of the potential in BZ (v1r_qbz)
       if (isirr_q) then
         v1r_qbz = v1r_qibz
       else
         call v1phq_rotate(db%cryst, qibz(:,iq_ibz), isym, itimrev, g0q, &
           ngfft, cplex_qibz, nfft, db%nspden, db%nsppol, db%mpi_enreg, v1r_qibz, v1r_qbz, xmpi_comm_self)
       end if

       !call littlegroup_q(db%cryst%nsym, qpt_bz, symq, db%cryst%symrec, db%cryst%symafm, timerev_q, prtvol=db%prtvol)
       !do ipc=1,db%natom3
       !  idir = mod(ipc - 1, 3) + 1; ipert = (ipc - idir) / 3 + 1
       !  call v1phq_symmetrize(db%cryst, idir, ipert, symq, ngfft, cplex_qibz, nfft, db%nspden, db%nsppol, &
       !                        db%mpi_enreg, v1r_qbz(:,:,:,ipc))
       !end do

       ! Multiply by e^{iqpt_bz.r}
       call times_eikr(qpt_bz, ngfft, nfft, db%nspden * db%natom3, v1r_qbz)

       if (db%add_lr /= 0) then
         ! Substract the long-range part of the potential.
         do imyp=1,db%my_npert
           ipc = db%my_pinfo(3, imyp)
           do ispden=1,db%nspden
             v1r_qbz(:, :, ispden, ipc) = v1r_qbz(:, :, ispden, ipc) - v1r_lr(:, :, imyp)
           end do
         end do
       end if

       ! Compute FT phases for this qpt_bz.
       call calc_eiqr(-qpt_bz, db%my_nrpt, db%my_rpt, emiqr)

       ! Slow FT.
       do imyp=1,db%my_npert
         ipc = db%my_pinfo(3, imyp)
         do ispden=1,db%nspden
           do ifft=1,nfft
             db%v1scf_rpt(1, :, ifft, ispden, imyp) = db%v1scf_rpt(1, :, ifft, ispden, imyp) &
                + emiqr(1, :) * v1r_qbz(1, ifft, ispden, ipc) &
                - emiqr(2, :) * v1r_qbz(2, ifft, ispden, ipc)

             db%v1scf_rpt(2, :, ifft, ispden, imyp) = db%v1scf_rpt(2, :, ifft, ispden, imyp) &
                + emiqr(1, :) * v1r_qbz(2, ifft, ispden, ipc) &
                + emiqr(2, :) * v1r_qbz(1, ifft, ispden, ipc)
           end do
           !call zgerc(db%my_nrpt, nfft, cone, emiqr, 1, v1r_qbz(:,:,ispden,ipc), 1, &
           !           db%v1scf_rpt(:,:,:,ispden,imyp), db%my_nrpt)

         end do ! ispden
       end do ! imyp
     end if

   end do ! iqst

   write(msg,'(2(a,i0),a)') " IBZ q-point [", iq_ibz, "/", nqibz, "]"
   call cwtime_report(msg, cpu, wall, gflops)
   ABI_FREE(v1r_qibz)
 end do ! iq_ibz

 ABI_CHECK(iqst == nqbz, "iqst /= nqbz")
 call wrtout(std_out, ch10//ch10)

 ABI_FREE(iperm)
 ABI_FREE(emiqr)
 ABI_FREE(qibz)
 ABI_FREE(qbz)
 ABI_FREE(indqq)
 ABI_FREE(iqs_dvdb)
 ABI_FREE(nqsts)
 ABI_FREE(v1r_qbz)
 ABI_FREE(v1r_lr)

 !call xmpi_sum(db%v1scf_rpt, db%comm, ierr)
 db%v1scf_rpt = db%v1scf_rpt / nqbz
 !if (method == 1) db%v1scf_rpt(2,:,:,:,:) = zero
 !call dvdb_enforce_asr(db)

 !do irpt=1,db%my_nrpt
 !  write(std_out,*)"v1scf_rpt:", db%v1scf_rpt(:, irpt, 1:5, 1, 1)
 !end do

 !do imyp=1,db%my_npert
 !  write(std_out, *)" Imaginary part for imyp:", imyp
 !  do ispden=1,db%nspden
 !    write(std_out, *)"min, max, avg:", &
 !       minval(abs(db%v1scf_rpt(2, :, :, ispden, imyp))), maxval(abs(db%v1scf_rpt(2, :, :, ispden, imyp))), &
 !       sum(abs(db%v1scf_rpt(2, :, :, ispden, imyp))) / (nfft * db%my_nrpt)
 !  end do
 !end do

 call cwtime_report(" Construction of W(R,r)", cpu_all, wall_all, gflops_all)

end subroutine dvdb_ftinterp_setup
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_maxw
!! NAME
!!  dvdb_get_maxw
!!
!! FUNCTION
!!
!! INPUTS
!!
!! FUNCTION
!!

subroutine dvdb_get_maxw(db, ngqpt, all_rpt, all_rmod, maxw)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(in) :: db
 integer,intent(in) :: ngqpt(3)
 real(dp),allocatable,intent(out) :: maxw(:,:), all_rpt(:,:), all_rmod(:)

!Local variables-------------------------------
!scalars
 integer :: imyp, ipc, irpt, irpt_tot, ifft, nfft, ispden, ierr, ii
 real(dp) :: phre
!arrays
 integer,allocatable :: iperm_irpt(:)
 real(dp) :: sc_rprimd(3,3)

! *************************************************************************

 ABI_CHECK(allocated(db%v1scf_rpt), "v1scf_rpt is not allocated (call dvdb_ftinterp_setup)")

 nfft = size(db%v1scf_rpt, dim=3)
 ABI_CALLOC(maxw, (db%nrtot, db%natom3))

 ! Need all RPTs to sort output.
 ABI_CALLOC(all_rpt, (3, db%nrtot))
 do irpt=1,db%my_nrpt
   irpt_tot = db%my_irpt2tot(irpt)
   all_rpt(:, irpt_tot) = db%my_rpt(:, irpt)
 end do
 call xmpi_sum(all_rpt, db%comm_rpt, ierr)

 do imyp=1,db%my_npert
   ipc = db%my_pinfo(3, imyp)
   do irpt=1,db%my_nrpt
     irpt_tot = db%my_irpt2tot(irpt)
     phre = zero
     do ispden=1,db%nspden
       do ifft=1,nfft
         phre = max(phre, db%v1scf_rpt(1,irpt,ifft,ispden,imyp) ** 2 + db%v1scf_rpt(2,irpt,ifft,ispden,imyp) ** 2)
       end do
     end do
     maxw(irpt_tot, ipc) = sqrt(phre)
   end do
 end do

 ! Handle parallelism
 if (db%nprocs_rpt /= 1) call xmpi_sum(maxw, db%comm_rpt, ierr)
 if (db%nprocs_pert /= 1) call xmpi_sum(maxw, db%comm_pert, ierr)

 sc_rprimd(:, 1) = ngqpt(1) * db%cryst%rprimd(:, 1)
 sc_rprimd(:, 2) = ngqpt(2) * db%cryst%rprimd(:, 2)
 sc_rprimd(:, 3) = ngqpt(3) * db%cryst%rprimd(:, 3)
 call sort_rpts(db%nrtot, all_rpt, db%cryst%rmet, iperm_irpt, rmod=all_rmod)

 ! Sort output results by |R|.
 all_rpt = all_rpt(:, iperm_irpt(:))
 do ii=1,db%natom3
   maxw(:, ii) = maxw(iperm_irpt(:), ii)
 end do

 ABI_FREE(iperm_irpt)

end subroutine dvdb_get_maxw
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_enforce_asr1
!! NAME
!!  dvdb_enforce_asr1
!!
!! FUNCTION
!!  Enforce acoustic sum rule on the DFPT potentials.
!!
!! INPUTS
!!
!! FUNCTION
!!

subroutine dvdb_enforce_asr1(db)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),target,intent(inout) :: db

!Local variables-------------------------------
!scalars
integer :: imyp, idir, ipert, ifft, nfft, ierr, irpt, ispden
!arrays
 real(dp) :: sumr(2)
 real(dp),allocatable :: asr_work(:,:,:,:)

! *************************************************************************

 ABI_CHECK(allocated(db%v1scf_rpt), "v1scf_rpt is not allocated (call dvdb_ftinterp_setup)")
 nfft = size(db%v1scf_rpt, dim=2)

 ! Enforcing ASR \sum_{R,kappa} v1(R, r, kappa, alpha) = f(r, alpha) = 0 for each r and alpha-direction
 call wrtout(std_out, " Enforcing ASR on W(R,r) i.e. \sum_{R,kappa} v1(R, r, kappa, alpha) = 0")
 ABI_CALLOC(asr_work, (2, nfft, db%nspden, 3))

 do imyp=1,db%my_npert
   idir = db%my_pinfo(1, imyp); ipert = db%my_pinfo(2, imyp)
   do ispden=1,db%nspden
     do ifft=1,nfft
       sumr = zero
       do irpt=1,db%my_nrpt
         sumr(:) = sumr(:) + db%my_wratm(irpt, ipert) * db%v1scf_rpt(:,irpt,ifft,ispden,imyp)
       end do
       asr_work(:,ifft,ispden, idir) = asr_work(:,ifft,ispden, idir) + sumr(:)
     end do
   end do
 end do
 if (db%nprocs_rpt /= 1) call xmpi_sum(asr_work, db%comm_rpt, ierr)
 if (db%nprocs_pert /= 1) call xmpi_sum(asr_work, db%comm_pert, ierr)
 asr_work = asr_work / (db%nrtot * db%natom)

 do idir=1,3
   write(std_out, *)" ASR breaking for idir:" ,idir
   write(std_out, *)"   For Re:  minval, maxval_fft:", minval(abs(asr_work(1,:,1,idir))), maxval(abs(asr_work(1,:,1,idir)))
   write(std_out, *)"   For Im:  minval, maxval_fft:", minval(abs(asr_work(2,:,1,idir))), maxval(abs(asr_work(2,:,1,idir)))
 end do

 ! Remove average value.
 do imyp=1,db%my_npert
   idir = db%my_pinfo(1, imyp)
   do ispden=1,db%nspden
     do ifft=1,nfft
       db%v1scf_rpt(1, :, ifft, ispden, imyp) = db%v1scf_rpt(1, :, ifft, ispden, imyp) - asr_work(1, ifft, ispden, idir)
       db%v1scf_rpt(2, :, ifft, ispden, imyp) = db%v1scf_rpt(2, :, ifft, ispden, imyp) - asr_work(2, ifft, ispden, idir)
     end do
   end do
 end do

 ! Here I change only the R=0 term
 !my_ir0 = -1
 !do irpt=1,db%my_nrpt
 !  if (all(abs(db%my_rpt(:,irpt)) <= tol10)) then
 !    my_ir0 = irpt; exit
 !  end if
 !end do

 !if (my_ir0 /= -1) then
 !  asr_work = asr_work * (nrtot * db%natom)
 !  do imyp=1,db%my_npert
 !    idir = db%my_pinfo(1, imyp)
 !    do ispden=1,db%nspden
 !      asr_work(:,:,ispden, idir) = asr_work(:,:,ispden, idir) - db%v1scf_rpt(:,my_ir0,:,ispden,imyp)
 !    end do
 !  end do
 !  asr_work = asr_work / db%natom

 !  do imyp=1,db%my_npert
 !    idir = db%my_pinfo(1, imyp)
 !    do ispden=1,db%nspden
 !      !write(std_out, *)"imyp ispden:", imyp, ispden, "max: ", maxval(v1r_qbz(:, :, ispden, imyp), dim=2)
 !      db%v1scf_rpt(:, my_ir0, :, ispden, imyp) = - asr_work(:, :, ispden, idir)
 !    end do
 !  end do
 !end if

 ABI_FREE(asr_work)

end subroutine dvdb_enforce_asr1
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/prepare_ftinterp
!! NAME
!!  prepare_ftinterp
!!
!! FUNCTION
!!

subroutine prepare_ftinterp(db, method, ngqpt, qptopt, nqshift, qshift, &
                            qibz, qbz, indqq, iperm, nqsts, iqs_dvdb, all_rpt, all_wghatm, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift, qptopt, method, comm
 class(dvdb_t),target,intent(in) :: db
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qshift(3,nqshift)
 integer,allocatable,intent(out) :: indqq(:,:),nqsts(:),iqs_dvdb(:), iperm(:)
 real(dp),allocatable,intent(out) :: all_rpt(:,:), all_wghatm(:,:,:)
 real(dp),allocatable,intent(out) :: qibz(:,:),qbz(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1, timrev1=1, brav1=1, cutmode2=2
 integer :: iq_ibz,nqibz,iq_bz,nqbz
 integer :: ii,iq_dvdb
 integer :: iqst,nqst,ix,iy,iz,nq1,nq2,nq3,r1,r2,r3, nrtot
 real(dp) :: dksqmax,r_inscribed_sphere
 logical :: found
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: bz2ibz_sort(:),all_cell(:,:)
 real(dp) :: shift(3), rcan(3, db%cryst%natom), trans(3, db%cryst%natom)
 real(dp) :: acell(3), rprim(3,3), gprim(3,3)
 real(dp),allocatable :: wtq(:),all_rcart(:,:)

! *************************************************************************

 cryst => db%cryst

 ! Generate q-mesh: find BZ, IBZ and the corresponding weights from ngqpt.
 nq1 = ngqpt(1); nq2 = ngqpt(2); nq3 = ngqpt(3)
 qptrlatt = 0; qptrlatt(1, 1) = ngqpt(1); qptrlatt(2, 2) = ngqpt(2); qptrlatt(3, 3) = ngqpt(3)

 ABI_CHECK(nqshift == 1, "nshift > 1 not supported")
 ABI_CHECK(all(qshift(:, 1) == zero), "qshift != 0 not supported")

 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt, nqshift, qshift, &
   nqibz, qibz, wtq, nqbz, qbz) ! new_kptrlatt, new_shiftk)

 ABI_CHECK(nqbz == product(ngqpt) * nqshift, "nqbz /= product(ngqpt) * nqshift")

 if (method == 0) then
#if 1
   ! We want a gamma centered q-mesh for the FFT.
   ! TODO: This is not needed anymore because we are gonna use slow FFT to handle a possible distribution of rpts
   ABI_FREE(qbz)
   ABI_MALLOC(qbz, (3, nqbz))
   ii = 0
   do iz=0,nq3-1
     do iy=0,nq2-1
       do ix=0,nq1-1
         ii = ii + 1
         qbz(:, ii) = [ix / dble(nq1), iy / dble(nq2), iz / dble(nq3)]
         call wrap2_pmhalf([ix / dble(nq1), iy / dble(nq2), iz / dble(nq3)], qbz(:,ii), shift)
       end do
     end do
   end do
#endif

   ! Compute real-space points.
   ! Use the following indexing (N means ngfft of the adequate direction)
   ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
   ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig
   nrtot = nqbz
   ABI_CALLOC(all_rpt, (3, nrtot))
   ii = 0
   do iz=1,nq3
     r3 = ig2gfft(iz, nq3)
     do iy=1,nq2
       r2 = ig2gfft(iy, nq2)
       do ix=1,nq1
         r1 = ig2gfft(ix, nq1)
         ii = ii + 1
         all_rpt(:, ii) = [r1, r2, r3]
       end do
     end do
   end do

 else if (method == 1) then
   ! Compute rprim, and gprim
   call mkradim(acell, rprim, cryst%rprimd)
   call canat9(brav1, cryst%natom, rcan, rprim, trans, cryst%xred)
   call matr3inv(rprim, gprim)
   ! TODO: In principle one may have N_R that depends on iatom to minimize memory.
   call get_bigbox_and_weights(brav1, cryst%natom, nqbz, ngqpt, nqshift, qshift, rprim, cryst%rprimd, gprim, rcan, &
                               cutmode2, nrtot, all_rcart, all_cell, all_wghatm, r_inscribed_sphere, comm)
   ABI_FREE(all_cell)
   do ii=1,nrtot
     all_rcart(:,ii) = all_rcart(:,ii) * acell(:)
   end do
   ABI_MALLOC(all_rpt, (3, nrtot))
   call xcart2xred(nrtot, cryst%rprimd, all_rcart, all_rpt)
   ABI_FREE(all_rcart)
 else
   MSG_ERROR(sjoin("Wrong method:", itoa(method)))
 end if

 ! Find correspondence BZ --> IBZ. Note:
 ! q --> -q symmetry is always used for phonons.
 ! we use symrec instead of symrel
 ABI_MALLOC(indqq, (nqbz*sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indqq, qibz, qbz, nqibz, nqbz, cryst%nsym, &
   sppoldbl1, cryst%symafm, cryst%symrec, timrev1, comm, exit_loop=.True., use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

 ! Construct sorted mapping BZ --> IBZ to speedup qbz search below.
 ABI_MALLOC(iperm, (nqbz))
 ABI_MALLOC(bz2ibz_sort, (nqbz))
 iperm = [(ii, ii=1,nqbz)]
 bz2ibz_sort = indqq(:,1)
 call sort_int(nqbz, bz2ibz_sort, iperm)

 ! Reconstruct the IBZ according to what is present in the DVDB.
 ABI_MALLOC(nqsts, (nqibz))
 ABI_MALLOC(iqs_dvdb, (nqibz))
 iqs_dvdb = -1

 iqst = 0
 do iq_ibz=1,nqibz
   ! In each q-point star, count the number of q-points and find the one present in the DVDB.
   nqst = 0
   found = .false.
   do ii=iqst+1,nqbz
     if (bz2ibz_sort(ii) /= iq_ibz) exit
     nqst = nqst + 1
     iq_bz = iperm(ii)
     if (.not. found) then
       iq_dvdb = db%findq(qbz(:,iq_bz))
       if (iq_dvdb /= -1) then
         qibz(:,iq_ibz) = qbz(:,iq_bz)
         iqs_dvdb(iq_ibz) = iq_dvdb
         found = .true.
       end if
     end if
   end do

   ! Check that nqst has been counted properly.
   ABI_CHECK(nqst > 0 .and. bz2ibz_sort(iqst+1) == iq_ibz, "Wrong iqst")
   if (abs(nqst - wtq(iq_ibz) * nqbz) > tol12) then
     write(msg, "(a,i0,a,f5.2)")"Error in q-point star or q-weights. nqst:", nqst, "wtq * nqbz = ", wtq(iq_ibz) * nqbz
     MSG_ERROR(msg)
   end if

   ! Check that the q-point has been found in DVDB.
   if (.not. found) then
     MSG_ERROR(sjoin("Cannot find symmetric q-point of:", ktoa(qibz(:,iq_ibz)), "in DVDB file"))
   end if

   iqst = iqst + nqst
   nqsts(iq_ibz) = nqst
 end do

 ABI_FREE(wtq)
 ABI_FREE(bz2ibz_sort)

 ! Redo the mapping with the new IBZ
 call listkk(dksqmax, cryst%gmet, indqq, qibz, qbz, nqibz, nqbz, cryst%nsym, &
   sppoldbl1, cryst%symafm, cryst%symrec, timrev1, comm, exit_loop=.True., use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

end subroutine prepare_ftinterp
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_ftinterp_qpt
!! NAME
!!  dvdb_ftinterp_qpt
!!
!! FUNCTION
!!  Fourier interpolation of potentials for a given q-point
!!  Internal tables must be prepared in advance by calling `dvdb_ftinterp_setup`.
!!
!!  \begin{equation}
!!  \partial v^{scf}_{\tilde\qq\kappa\alpha}(\rr) \approx \sum_\RR e^{+i\tilde{\qq}\cdot(\RR - \rr)} W_{\kappa\alpha}(\rr,\RR).
!!  \end{equation
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates (arbitrary point in the BZ).
!!  nfft=Number of FFT-points treated by this processors.
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  comm=MPI communicator for R-points.
!!  [add_lr]= If present, use this value for the LR treatment instead of dv%add_lr
!!
!! OUTPUT
!!  ov1r(2*nfft, nspden, my_npert)=Interpolated DFPT potentials at the given q-point (periodic part)
!!
!! PARENTS
!!      m_dvdb,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_ftinterp_qpt(db, qpt, nfft, ngfft, ov1r, comm_rpt, add_lr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft, comm_rpt
 integer,optional,intent(in) :: add_lr
 class(dvdb_t),intent(in) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qpt(3)
 ! SHOULD BE v1scf(cplex, nfft, nspden, db%my_npert)= v1scf potentials on the real-space FFT mesh
 real(dp),intent(out) :: ov1r(2, nfft, db%nspden, db%my_npert)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex2 = 2
 integer :: ir, ispden, ifft, imyp, idir, ipert, timerev_q, ierr, my_add_lr
 real(dp) :: wr, wi, qmod
 !complex(dpc) :: beta
!arrays
 integer :: symq(4,2,db%cryst%nsym), rfdir(3)
 integer,allocatable :: pertsy(:,:), rfpert(:), pflag(:,:)
 real(dp) :: qcart(3)
 real(dp),allocatable :: eiqr(:,:), weiqr(:,:), v1r_lr(:,:,:)

! *************************************************************************

 my_add_lr = db%add_lr; if (present(add_lr)) my_add_lr = add_lr

 qcart = two_pi * matmul(db%cryst%gprimd, qpt)
 qmod = sqrt(dot_product(qcart, qcart))

 if (my_add_lr == 4) then
   ov1r = zero
   do imyp=1,db%my_npert
     idir = db%my_pinfo(1, imyp); ipert = db%my_pinfo(2, imyp)
     call db%get_v1r_long_range(qpt, idir, ipert, nfft, ngfft, ov1r(:, :, 1, imyp))
     ! Remove the phase to get the lattice-periodic part.
     call times_eikr(-qpt, ngfft, nfft, 1, ov1r(:, :, 1, imyp))
     if (db%nspden /= 1) ov1r(:, :, 2, imyp) = ov1r(:, :, 1, imyp)
   end do
   return
 end if

 ABI_CHECK(allocated(db%v1scf_rpt), "v1scf_rpt is not allocated (call dvdb_ftinterp_setup)")

 ! Examine the symmetries of the q-wavevector
 call littlegroup_q(db%cryst%nsym, qpt, symq, db%cryst%symrec, db%cryst%symafm, timerev_q, prtvol=db%prtvol)

 ! Compute e^{iqR} FT phases for this q-point.
 ABI_MALLOC(weiqr, (2, db%my_nrpt))
 ABI_MALLOC(eiqr, (2, db%my_nrpt))
 call calc_eiqr(qpt, db%my_nrpt, db%my_rpt, eiqr)

 ! Compute long-range part of the coupling potential
 if (my_add_lr > 0) then
   ABI_MALLOC(v1r_lr, (2, nfft, db%my_npert))
   do imyp=1,db%my_npert
     idir = db%my_pinfo(1, imyp); ipert = db%my_pinfo(2, imyp)
     call db%get_v1r_long_range(qpt, idir, ipert, nfft, ngfft, v1r_lr(:,:,imyp))
   end do
 end if

 ! Interpolate potentials (results in ov1r)
 ov1r = zero
 do imyp=1,db%my_npert
   idir = db%my_pinfo(1, imyp); ipert = db%my_pinfo(2, imyp)
   weiqr(1,:) = db%my_wratm(:, ipert) * eiqr(1,:)
   weiqr(2,:) = db%my_wratm(:, ipert) * eiqr(2,:)

   do ispden=1,db%nspden

     ! Slow FT.
     do ifft=1,nfft
       do ir=1,db%my_nrpt
         wr = db%v1scf_rpt(1, ir, ifft, ispden, imyp)
         wi = db%v1scf_rpt(2, ir, ifft, ispden, imyp)
         ov1r(1, ifft, ispden, imyp) = ov1r(1, ifft, ispden, imyp) + wr * weiqr(1, ir) - wi * weiqr(2, ir)
         ov1r(2, ifft, ispden, imyp) = ov1r(2, ifft, ispden, imyp) + wr * weiqr(2, ir) + wi * weiqr(1, ir)
       end do
       ! Add the long-range part of the potential
       if (my_add_lr > 0) then
         ov1r(1, ifft, ispden, imyp) = ov1r(1, ifft, ispden, imyp) + v1r_lr(1, ifft, imyp)
         ov1r(2, ifft, ispden, imyp) = ov1r(2, ifft, ispden, imyp) + v1r_lr(2, ifft, imyp)
       end if
     end do ! ifft

     !beta = czero
     !if (my_add_lr > 0) then
     !  beta = cone
     !  ov1r(:, :, ispden, imyp) = v1r_lr(:, :, imyp)
     !end if
     !call ZGEMV("T", db%my_nrpt, nfft, cone, db%v1scf_rpt(1,1,1,ispden,imyp), db%mynrpt, weiqr, 1, &
     !           beta, ov1r(1,1,ispden,imyp), 1)

     ! Remove the phase to get the lattice-periodic part.
     call times_eikr(-qpt, ngfft, nfft, 1, ov1r(:, :, ispden, imyp))

     ! Need to collect results if R-points are distributed (TODO Non-blocking API?)
     if (db%nprocs_rpt > 1) call xmpi_sum(ov1r(:,:,ispden,imyp), comm_rpt, ierr)
   end do ! ispden

   ! Be careful with gamma and cplex!
   if (db%symv1==1) then !(.and. reveiver == -1 .or. receiver == db%comm_rpt%my_rank)
     call v1phq_symmetrize(db%cryst, idir, ipert, symq, ngfft, cplex2, nfft, db%nspden, db%nsppol, &
         db%mpi_enreg, ov1r(:,:,:,imyp))
   end if
 end do ! imyp

 if (db%symv1 == 2) then
   ! Initialize the list of perturbations rfpert and rdfir
   ! WARNING: Only phonon perturbations are considered for the time being.
   ABI_MALLOC(rfpert, (db%mpert))
   rfpert = 0; rfpert(1:db%cryst%natom) = 1; rfdir = 1
   ABI_MALLOC(pertsy, (3,db%mpert))
   ABI_MALLOC(pflag, (3, db%natom))

   ! Determine the symmetrical perturbations. Meaning of pertsy:
   !    0 for non-target perturbations
   !    1 for basis perturbations
   !   -1 for perturbations that can be found from basis perturbations
   call irreducible_set_pert(db%cryst%indsym,db%mpert,db%cryst%natom,db%cryst%nsym,&
       pertsy,rfdir,rfpert,symq,db%cryst%symrec,db%cryst%symrel)

   pflag = 0
   do imyp=1,3*db%cryst%natom
     idir = mod(imyp-1, 3) + 1; ipert = (imyp - idir) / 3 + 1
     if (pertsy(idir, ipert) == 1) pflag(idir,ipert) = 1
   end do

   ! Complete potentials
   call v1phq_complete(db%cryst,qpt,ngfft,cplex2,nfft,db%nspden,db%nsppol,db%mpi_enreg,db%symv1,pflag,ov1r)

   ABI_FREE(pertsy)
   ABI_FREE(rfpert)
   ABI_FREE(pflag)
 endif

 ! Set imaginary part to zero if gamma point.
 if (sum(qpt**2) < tol14) ov1r(2, :, :, :) = zero

 ABI_FREE(eiqr)
 ABI_FREE(weiqr)
 ABI_SFREE(v1r_lr)

end subroutine dvdb_ftinterp_qpt
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_ftqbz
!! NAME
!!  dvdb_get_ftqbz
!!
!! FUNCTION
!!  Fourier interpolation of potentials for a given q-point in the BZ (qbz).
!!  Internal tables must be prepared in advance by calling `dvdb_ftinterp_setup`.
!!  Interpolation is avoided if potential is in ft_qcache
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  qbz(3)= q-point in the BZ in reduced coordinates.
!!  qibz(3) = Image of qbz in the IBZ. Needed to apply symmetries IBZ --> qbz
!!  indq2ibz(6)=Symmetry mapping qbz --> IBZ qpoint produced by listkk.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!
!! OUTPUT
!!  ov1r(2*nfft, nspden, my_npert)=Interpolated DFPT potentials at the given q-point.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_get_ftqbz(db, cryst, qbz, qibz, indq2ibz, cplex, nfft, ngfft, v1scf, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft, comm
 integer,intent(out) :: cplex
 type(crystal_t),intent(in) :: cryst
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18), indq2ibz(6)
 real(dp),intent(in) :: qbz(3), qibz(3)
 real(dp),allocatable,intent(out) :: v1scf(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz, itimrev, isym, ierr, imyp, mu, root
 logical :: isirr_q, incache
!!arrays
 integer :: g0q(3), requests(db%natom3)
 real(dp) :: tsec(2)
 real(dp) ABI_ASYNC, allocatable :: work(:,:,:,:), work2(:,:,:,:)

! *************************************************************************
 ABI_UNUSED(comm)

 ! Keep track of total time spent.
 call timab(1809, 1, tsec)

 iq_ibz = indq2ibz(1)
 db%ft_qcache%count_qused(iq_ibz) = db%ft_qcache%count_qused(iq_ibz) + 1
 db%ft_qcache%stats(1) = db%ft_qcache%stats(1) + 1

 ! IS(q_ibz) + g0q = q_bz
 isym = indq2ibz(2); itimrev = indq2ibz(6) + 1; g0q = indq2ibz(3:5)
 isirr_q = (isym == 1 .and. itimrev == 1 .and. all(g0q == 0))

 if (db%ft_qcache%use_3natom_cache .and. db%ft_qcache%stored_iqibz_cplex(1) == iq_ibz .and. .not. isirr_q) then
   ! All 3 natom potentials for qibz are in cache. Symmetrize to get Sq without MPI communication.
   db%ft_qcache%stats(2) = db%ft_qcache%stats(2) + 1
   cplex = db%ft_qcache%stored_iqibz_cplex(2)
   ABI_MALLOC(work2, (cplex, nfft, db%nspden, db%natom3))
   call v1phq_rotate(cryst, qibz, isym, itimrev, g0q, ngfft, cplex, nfft, &
     db%nspden, db%nsppol, db%mpi_enreg, db%ft_qcache%v1scf_3natom_qibz, work2, db%comm_pert)
   ! Extract my data from work2
   ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
   do imyp=1,db%my_npert
     v1scf(:,:,:,imyp) = work2(:,:,:,db%my_pinfo(3, imyp))
   end do
   ABI_FREE(work2)
   call timab(1809, 2, tsec); return
 end if

 ! TODO: Note that cplex is always set to 2 here
 cplex = 2

 ! Check whether iq_ibz is in cache.
 incache = .False.
 if (db%ft_qcache%maxnq > 0) then
   ! Get number of perturbations computed for this iq_ibz as well as cplex.
   ! Remember that the size of v1scf in qcache depends on db%my_npert
   if (allocated(db%ft_qcache%key(iq_ibz)%v1scf)) then
      if (size(db%ft_qcache%key(iq_ibz)%v1scf, dim=1) == cplex .and. &
          size(db%ft_qcache%key(iq_ibz)%v1scf, dim=2) == nfft) then
        ! Potential in cache --> copy it in output v1scf.
        ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
        v1scf = real(db%ft_qcache%key(iq_ibz)%v1scf, kind=QCACHE_KIND)
        incache = .True.
        db%ft_qcache%stats(3) = db%ft_qcache%stats(3) + 1
      else
        MSG_ERROR("Found different cplex/nfft in cache.")
      end if
   else
      !call wrtout(std_out, sjoin("Cache miss for iq_ibz. Will try to interpolate it...", itoa(iq_ibz)))
      db%ft_qcache%stats(4) = db%ft_qcache%stats(4) + 1
   end if
 end if

 if (.not. incache) then
   !MSG_ERROR("For the time being incache should always be true when FT interpolation is used!")
   ! Interpolate the dvscf potentials directly in the **BZ** for my_npert perturbations.
   ! This is possible only if all procs inside comm_rpt call this routine else deadlock
   ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
   call db%ftinterp_qpt(qbz, nfft, ngfft, v1scf, db%comm_rpt)
   call timab(1809, 2, tsec); return
 end if

 if (.not. isirr_q) then
   ! Must rotate iq_ibz to get potential for qpoint in the BZ.
   ! Be careful with the shape of output v1scf because the routine returns db%my_npert potentials.

   if (db%my_npert == db%natom3) then
     ABI_MALLOC(work, (cplex, nfft, db%nspden, db%natom3))
     work = v1scf
     call v1phq_rotate(cryst, qibz, isym, itimrev, g0q, ngfft, cplex, nfft, &
       db%nspden, db%nsppol, db%mpi_enreg, work, v1scf, db%comm_pert)
     ABI_FREE(work)

   else
     ! Parallelism over perturbations.
     ABI_MALLOC(work2, (cplex, nfft, db%nspden, db%natom3))

     if (incache) then
       ! Cache is distributed --> have to collect all 3 natom perts inside db%comm_pert.
       ABI_MALLOC(work, (cplex, nfft, db%nspden, db%natom3))

       ! IBCAST is much faster than a naive xmpi_sum.
       call timab(1806, 1, tsec)
       do mu=1,db%natom3
         root = db%pert_table(1, mu)
         if (root == db%me_pert) then
           imyp = db%pert_table(2, mu)
           work(:,:,:,mu) = v1scf(:,:,:,imyp)
         end if
         call xmpi_ibcast(work(:,:,:,mu), root, db%comm_pert, requests(mu), ierr)
       end do
       call xmpi_waitall(requests, ierr)
       call timab(1806, 2, tsec)

       call v1phq_rotate(cryst, qibz, isym, itimrev, g0q, ngfft, cplex, nfft, &
         db%nspden, db%nsppol, db%mpi_enreg, work, work2, db%comm_pert)

       ! Store all 3 natom potentials for q in IBZ in cache.
       if (db%ft_qcache%use_3natom_cache .and. db%ft_qcache%stored_iqibz_cplex(1) /= iq_ibz) then
         if (cplex /= db%ft_qcache%stored_iqibz_cplex(2)) then
           ABI_REMALLOC(db%ft_qcache%v1scf_3natom_qibz, (cplex, nfft, db%nspden, db%natom3))
         end if
         db%ft_qcache%stored_iqibz_cplex = [iq_ibz, cplex]
         db%ft_qcache%v1scf_3natom_qibz = work
       end if
       ABI_FREE(work)

     else
       NOT_IMPLEMENTED_ERROR()
       ! All 3 natom have been read in v1scf by dvdb_readsym_allv1
       !call v1phq_rotate(cryst, qibz, isym, itimrev, g0q, ngfft, cplex, nfft, &
       !  db%nspden, db%nsppol, db%mpi_enreg, v1scf, work2, db%comm_pert)
     end if

     ! Reallocate v1scf with my_npert and extract data from work2.
     ABI_REMALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
     do imyp=1,db%my_npert
       v1scf(:,:,:,imyp) = work2(:,:,:,db%my_pinfo(3, imyp))
     end do
     ABI_FREE(work2)
   end if
 end if ! not isirr_q

 call timab(1809, 2, tsec)

end subroutine dvdb_get_ftqbz
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_ftqcache_build
!! NAME
!!  dvdb_ftqcache_build
!!
!! FUNCTION
!!  This function initializes the internal q-cache from W(R,r).
!!  This is a collective routine that must be called by all procs inside comm.
!!
!! INPUTS
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  nqibz=Number of points in the IBZ.
!!  qibz(3, nqibz)=q-points in the IBZ.
!!  mbsize: Cache size in megabytes.
!!    < 0 to allocate all q-points.
!!    0 has not effect.
!!    > 0 for cache with automatically computed nqpt points.
!!  qselect_ibz(nqibz)=0 to ignore this q-point (global array)
!!  itreatq(nqibz) = 0 if this q-point won't be treated by this CPU else > 0.
!!    Each CPU calls this routine with its own array.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_ftqcache_build(db, nfft, ngfft, nqibz, qibz, mbsize, qselect_ibz, itreatq, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft, nqibz, comm
 real(dp),intent(in) :: mbsize
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18), qselect_ibz(nqibz)
 integer(i1b),intent(in) :: itreatq(nqibz)
 real(dp),intent(in) :: qibz(3, nqibz)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz, cplex, ierr, my_cplex, db_iqpt, imyp
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all, my_mbsize, max_mbsize
 logical :: read_from_dvdb
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: v1scf(:,:,:,:), all_v1scf(:,:,:,:)

! *************************************************************************

 call timab(1808, 1, tsec)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 call wrtout(std_out, ch10//" Precomputing Vscf(q) from W(R,r) and building qcache...", do_flush=.True.)
 db%ft_qcache = qcache_new(nqibz, nfft, ngfft, mbsize, db%natom3, db%my_npert, db%nspden)
 db%ft_qcache%itreatq(:) = itreatq

 ! TODO: Note that cplex is always set to 2 here
 cplex = 2
 ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
 read_from_dvdb = .False.

 do iq_ibz=1,nqibz
   ! Ignore points reported by the oracle. We can still recompute them on the fly if needed.
   if (qselect_ibz(iq_ibz) == 0) cycle
   if (.not. read_from_dvdb .and. itreatq(iq_ibz) == 0) cycle

   call cwtime(cpu, wall, gflops, "start")

   ! FIXME: DEBUGGING SECTION
   if (read_from_dvdb) then
     ! Read q-points from file! Incompatible with itreatq cycle and comm_rpt
     db_iqpt = db%findq(qibz(:, iq_ibz))
     if (db_iqpt /= -1) then
       call wrtout(std_out, " DBG: Reading V(q) from DVDB...")
       call db%readsym_allv1(db_iqpt, my_cplex, nfft, ngfft, all_v1scf, comm)
       do imyp=1,db%my_npert
         if (my_cplex == 2) then
           v1scf(:,:,:,imyp) = all_v1scf(:,:,:,db%my_pinfo(3, imyp))
         else
           v1scf(1,:,:,imyp) = all_v1scf(1,:,:,db%my_pinfo(3, imyp))
           v1scf(2,:,:,imyp) = zero
         end if
       end do
       ABI_FREE(all_v1scf)
       write(msg,'(2(a,i0),a)') " Reading q-point [", iq_ibz, "/", nqibz, "]"
       call cwtime_report(msg, cpu, wall, gflops)
     else
       call db%ftinterp_qpt(qibz(:, iq_ibz), nfft, ngfft, v1scf, db%comm_rpt)
     end if

   else
     ! Interpolate my_npert potentials inside comm_rpt
     call db%ftinterp_qpt(qibz(:, iq_ibz), nfft, ngfft, v1scf, db%comm_rpt)
   end if

   ! Points in the IBZ may be distributed to reduce memory.
   if (db%ft_qcache%itreatq(iq_ibz) /= 0) then
     ABI_MALLOC(db%ft_qcache%key(iq_ibz)%v1scf, (cplex, nfft, db%nspden, db%my_npert))
     db%ft_qcache%key(iq_ibz)%v1scf = real(v1scf, kind=QCACHE_KIND)
   end if

   ! Print progress.
   if (iq_ibz <= 50 .or. mod(iq_ibz, 100) == 0) then
     write(msg,'(2(a,i0),a)') " Interpolating q-point [", iq_ibz, "/", nqibz, "]"
     call cwtime_report(msg, cpu, wall, gflops)
   end if
 end do

 ABI_FREE(v1scf)

 ! Compute final cache size.
 my_mbsize = db%ft_qcache%get_mbsize()
 call wrtout(std_out, sjoin(" Memory allocated for cache: ", ftoa(my_mbsize, fmt="f8.1"), " [Mb] <<< MEM"))
 call xmpi_max(my_mbsize, max_mbsize, comm, ierr)
 call wrtout(std_out, sjoin(" Max memory inside MPI comm: ", ftoa(max_mbsize, fmt="f8.1"), " [Mb] <<< MEM"))
 call cwtime_report(" Qcache from W(R,r) + symmetrization", cpu_all, wall_all, gflops_all, end_str=ch10)
 call timab(1808, 2, tsec)

 ! This barrier seems to be needed on lemaitre3. DO NOT REMOVE!
 call xmpi_barrier(comm)

end subroutine dvdb_ftqcache_build
!!***

!!****f* m_dvdb/dvdb_ftqcache_update_from_ft
!! NAME
!!  dvdb_ftqcache_update_from_ft
!!
!! FUNCTION
!!  Interpolate selected potentials and update the internal q-cache.
!!  This is a collective routine that must be called by all procs in comm.
!!
!! INPUTS
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  nqibz=Number of q-points in the IBZ.
!!  qibz(3, nqibz)= q-points in the IBZ.
!!  ineed_qpt(%nqpt)=1 if this MPI rank requires this q-point.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_ftqcache_update_from_ft(db, nfft, ngfft, nqibz, qibz, ineed_qpt, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,comm, nqibz
 class(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18), ineed_qpt(nqibz)
 real(dp),intent(in) :: qibz(3, nqibz)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz, cplex, ierr, qcnt
 real(dp) :: cpu_all, wall_all, gflops_all, mbsize, max_mbsize
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: v1scf(:,:,:,:)
 !real(dp) :: tsec(2)

! *************************************************************************

 if (db%ft_qcache%maxnq == 0) return
 qcnt = count(ineed_qpt /= 0)

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 if (qcnt /= 0) then
   !call timab(1807, 1, tsec)
   call wrtout(std_out, sjoin(" Need to update Vscf(q) cache with: ", itoa(qcnt), "q-points from FT..."), do_flush=.True.)
   if (db%ft_qcache%make_room(ineed_qpt, msg) /= 0) then
     MSG_WARNING(msg)
   end if

   cplex = 2
   ABI_MALLOC(v1scf, (cplex, nfft, db%nspden, db%my_npert))
   do iq_ibz=1,nqibz
     if (ineed_qpt(iq_ibz) == 0) cycle
     ! Interpolate my_npert potentials inside comm_rpt.
     call db%ftinterp_qpt(qibz(:, iq_ibz), nfft, ngfft, v1scf, db%comm_rpt)
     ! Transfer to cache.
     if (ineed_qpt(iq_ibz) /= 0) then
       ABI_MALLOC(db%ft_qcache%key(iq_ibz)%v1scf, (cplex, nfft, db%nspden, db%my_npert))
       db%ft_qcache%key(iq_ibz)%v1scf = real(v1scf, kind=QCACHE_KIND)
     end if
   end do

   ABI_FREE(v1scf)
 end if

 mbsize = db%ft_qcache%get_mbsize()
 call wrtout(std_out, sjoin(" Memory allocated for cache: ", ftoa(mbsize, fmt="f8.1"), " [Mb] <<< MEM"))
 call xmpi_max(mbsize, max_mbsize, comm, ierr)
 call wrtout(std_out, sjoin(" Max memory inside MPI comm: ", ftoa(max_mbsize, fmt="f8.1"), " [Mb] <<< MEM"))
 call cwtime_report(" dvdb_ftqcache_update_from_ft", cpu_all, wall_all, gflops_all)
 !call timab(1807, 2, tsec)

end subroutine dvdb_ftqcache_update_from_ft
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_v1scf_rpt
!! NAME
!!  dvdb_get_v1scf_rpt
!!
!! FUNCTION
!!  Compute the phonon perturbation potential in real space lattice representation.
!!  This routine is meant to replace dvdb_ftinterp_setup
!!  and performs the potential interpolation one perturbation at a time.
!!
!! INPUTS
!!  ngqpt(3)=Divisions of the ab-initio q-mesh.
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nrpt=Number of R-points = number of q-points in the full BZ
!!  nspden=Number of spin densities.
!!  ipert=index of the perturbation to be treated [1,natom3]
!!  comm=MPI communicator
!!
!! OUTPUT
!!  v1scf_rpt(2,nrpt,nfft,nspden)
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_get_v1scf_rpt(db, cryst, ngqpt, nqshift, qshift, nfft, ngfft, &
&                             nrpt, nspden, ipert, v1scf_rpt, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,nrpt,nspden,ipert,comm
 class(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngqpt(3),ngfft(18)
 real(dp),intent(in) :: qshift(3,nqshift)
 real(dp),intent(out) :: v1scf_rpt(2,nrpt,nfft,nspden)
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1,timrev1=1
 integer :: my_qptopt,iq_ibz,nqibz,iq_bz,nqbz
 integer :: ii,iq_dvdb,cplex_qibz,ispden,irpt,idir,iat
 integer :: iqst,nqst,itimrev,tsign,isym,ix,iy,iz,nq1,nq2,nq3,r1,r2,r3
 integer :: nproc,my_rank,ifft,cnt,ierr
 character(len=500) :: msg
 real(dp) :: dksqmax
 logical :: isirr_q, found
!arrays
 integer :: qptrlatt(3,3),g0q(3)
 integer,allocatable :: indqq(:,:),iperm(:),bz2ibz_sort(:),nqsts(:),iqs_dvdb(:)
 real(dp) :: qpt_bz(3),shift(3)
 real(dp) :: cpu, wall, gflops
 real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:),emiqr(:,:)
 real(dp),allocatable :: v1r_qibz(:,:,:,:),v1r_qbz(:,:,:,:), v1r_lr(:,:)

! *************************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 nq1 = ngqpt(1); nq2 = ngqpt(2); nq3 = ngqpt(3); my_qptopt = 1 !; if (present(qptopt)) my_qptopt = qptopt

 ! Generate the q-mesh by finding the IBZ and the corresponding weights.
 qptrlatt = 0
 do ii=1,3
   qptrlatt(ii,ii) = ngqpt(ii)
 end do

 ! Get IBZ and BZ.
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, my_qptopt, nqshift, qshift, &
   nqibz, qibz, wtq, nqbz, qbz) ! new_kptrlatt, new_shiftk)

 !write(std_out,*)"Irreducible q-points:"
 !do iq_ibz=1,nqibz; write(std_out,*)trim(ktoa(qibz(:,iq_ibz))),wtq(iq_ibz)*nqbz; end do
 ABI_CHECK(nqbz == product(ngqpt) * nqshift, "nqbz /= product(ngqpt) * nqshift")
 ABI_CHECK(nqbz == nrpt, "nqbz /= nrpt")

 db%my_nrpt = nqbz

 ABI_CHECK(nspden == db%nspden, "nspden /= db%nspden")

 ! We want a gamma centered q-mesh for FFT.
 ABI_CHECK(nqshift == 1, "nshift > 1 not supported")
 ABI_CHECK(all(qshift(:,1) == zero), "qshift != 0 not supported")

 ABI_FREE(qbz)
 ABI_MALLOC(qbz, (3, nqbz))
 ii = 0
 do iz=0,nq3-1
   do iy=0,nq2-1
     do ix=0,nq1-1
       ii = ii + 1
       qbz(:, ii) = [ix/dble(nq1), iy/dble(nq2), iz/dble(nq3)]
       call wrap2_pmhalf([ix/dble(nq1), iy/dble(nq2), iz/dble(nq3)], qbz(:,ii), shift)
     end do
   end do
 end do

 ! Compute real-space points.
 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig
 ABI_MALLOC(db%my_rpt, (3, db%my_nrpt))
 ii = 0
 do iz=1,nq3
   r3 = ig2gfft(iz,nq3)
   do iy=1,nq2
     r2 = ig2gfft(iy,nq2)
     do ix=1,nq1
       r1 = ig2gfft(ix,nq1)
       ii = ii + 1
       db%my_rpt(:,ii) = [r1, r2, r3]
     end do
   end do
 end do

 ! Find correspondence BZ --> IBZ. Note:
 ! q --> -q symmetry is always used for phonons.
 ! we use symrec instead of symrel
 ABI_MALLOC(indqq, (nqbz*sppoldbl1,6))
 call listkk(dksqmax,cryst%gmet,indqq,qibz,qbz,nqibz,nqbz,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,comm,exit_loop=.True., use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

 ! Construct sorted mapping BZ --> IBZ to speedup qbz search below.
 ABI_MALLOC(iperm, (nqbz))
 ABI_MALLOC(bz2ibz_sort, (nqbz))
 iperm = [(ii, ii=1,nqbz)]
 bz2ibz_sort = indqq(:,1)
 call sort_int(nqbz, bz2ibz_sort, iperm)

 ! Reconstruct the IBZ according to what is present in the DVDB.
 ABI_MALLOC(nqsts, (nqibz))
 ABI_MALLOC(iqs_dvdb, (nqibz))
 ABI_MALLOC(v1r_lr, (2,nfft))

 iqst = 0
 do iq_ibz=1,nqibz
   ! In each q-point star, count the number of q-points and find the one present in DVDB.
   nqst = 0
   found = .false.
   do ii=iqst+1,nqbz
     if (bz2ibz_sort(ii) /= iq_ibz) exit
     nqst = nqst + 1

     iq_bz = iperm(ii)
     if (.not. found) then
       iq_dvdb = db%findq(qbz(:,iq_bz))
       if (iq_dvdb /= -1) then
         qibz(:,iq_ibz) = qbz(:,iq_bz)
         iqs_dvdb(iq_ibz) = iq_dvdb
         found = .true.
       end if
     end if
   end do

   ! Check that nqst has been counted properly.
   ABI_CHECK(nqst > 0 .and. bz2ibz_sort(iqst+1) == iq_ibz, "Wrong iqst")
   if (abs(nqst - wtq(iq_ibz) * nqbz) > tol12) then
     write(std_out,*)nqst, wtq(iq_ibz) * nqbz
     MSG_ERROR("Error in counting q-point star or in the weights.")
   end if

   ! Check that the q-point has been found in DVDB.
   if (.not. found) then
     MSG_ERROR(sjoin("Cannot find symmetric q-point of:", ktoa(qibz(:,iq_ibz)), "in DVDB file"))
   end if
   !write(std_out,*)sjoin("qpt irred:",ktoa(qibz(:,iq_ibz)))

   iqst = iqst + nqst
   nqsts(iq_ibz) = nqst
 end do

 ! Redo the mapping with the new IBZ
 call listkk(dksqmax,cryst%gmet,indqq,qibz,qbz,nqibz,nqbz,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,comm,exit_loop=.True., use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

 ABI_MALLOC(emiqr, (2, db%my_nrpt))
 v1scf_rpt = zero

 ABI_MALLOC_OR_DIE(v1r_qbz, (2, nfft, db%nspden, db%natom3), ierr)
 !v1r_qbz = huge(one)

 iqst = 0
 call cwtime(cpu, wall, gflops, "start")
 do iq_ibz=1,nqibz

   ! Get potentials for this IBZ q-point on the real-space FFT mesh.
   ! This call allocates v1r_qibz(cplex_qibz, nfft, nspden, 3*natom)
   ! TODO: Interface with qcache
   call db%readsym_allv1(iqs_dvdb(iq_ibz), cplex_qibz, nfft, ngfft, v1r_qibz, comm)

   ! Reconstruct by symmetry the potentials for the star of this q-point, perform FT and accumulate
   ! Be careful with the gamma point.
   do ii=1,nqsts(iq_ibz)
     iqst = iqst + 1
     iq_bz = iperm(iqst)
     ABI_CHECK(iq_ibz == indqq(iq_bz,1), "iq_ibz !/ ind qq(1)")
     isym = indqq(iq_bz,2); itimrev = indqq(iq_bz,6) + 1; g0q = indqq(iq_bz,3:5) ! IS(q_ibz) + g0q = q_bz
     tsign = 3-2*itimrev

     qpt_bz = qbz(:, iq_bz)
     !write(std_out,*)"  treating:",trim(ktoa(qpt_bz))
     isirr_q = (isym == 1 .and. itimrev == 1 .and. all(g0q == 0))
     !ABI_CHECK(all(g0q == 0), "g0q /= 0")

     ! Compute long-range part of the coupling potential
     !call cwtime(cpu, wall, gflops, "start")
     v1r_lr = zero; cnt = 0
     if (db%add_lr /= 0) then
       idir = mod(ipert-1, 3) + 1; iat = (ipert - idir) / 3 + 1
       call db%get_v1r_long_range(qpt_bz, idir, iat, nfft, ngfft, v1r_lr)
     end if
     !call cwtime_report(" dvdb_get_v1r_long_range", cpu, wall, gflops)

     if (cplex_qibz == 1) then
       ! Gamma point.
       ABI_CHECK(nqsts(iq_ibz) == 1, "cplex_qibz == 1 and nq nqst /= 1 (should be gamma)")
       ABI_CHECK(all(g0q == 0), "gamma point with g0q /= 0")

       ! Substract the long-range part of the potential
       if (db%add_lr /= 0) then
         do ispden=1,db%nspden
           v1r_qibz(1,:,ispden,ipert) = v1r_qibz(1,:,ispden,ipert) - v1r_lr(1,:)
         end do
       end if

       ! SLOW FT.
       !call cwtime(cpu, wall, gflops, "start")
       cnt = 0
       do ispden=1,db%nspden
         do irpt=1,db%my_nrpt
           ! MPI-parallelism
           cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
           do ifft=1,nfft
             v1scf_rpt(1,irpt,ifft,ispden) = v1scf_rpt(1,irpt,ifft,ispden) + &
                                             v1r_qibz(1, ifft, ispden, ipert)
           end do
         end do
       end do
       !call cwtime_report(" slow fft", cpu, wall, gflops)

     else
       ! q /= Gamma
       ! Get the periodic part of the potential in BZ (v1r_qbz)
       if (isirr_q) then
         !write(std_out,*)sjoin("qpt irred:",ktoa(qpt_bz))
         v1r_qbz = v1r_qibz
       else
         !call cwtime(cpu, wall, gflops, "start")
         call v1phq_rotate(cryst, qibz(:,iq_ibz), isym, itimrev, g0q,&
           ngfft, cplex_qibz, nfft, db%nspden, db%nsppol, db%mpi_enreg, v1r_qibz, v1r_qbz, xmpi_comm_self)
         !call cwtime_report(" rotate", cpu, wall, gflops)
         !v1r_qbz = zero; v1r_qbz = v1r_qibz

         !call times_eigr(-tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
         !call times_eigr(+tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
         !if (itimrev == 2) v1r_qbz(2,:,:,:) = -v1r_qbz(2,:,:,:)
         !call times_eigr(-tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
         !call times_eigr(+tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
       end if

       ! Multiply by e^{iqpt_bz.r}
       call times_eikr(qpt_bz, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)

       ! Substract the long-range part of the potential
       if (db%add_lr /= 0) then
         do ispden=1,db%nspden
           v1r_qbz(1,:,ispden,ipert) = v1r_qbz(1,:,ispden,ipert) - v1r_lr(1,:)
           v1r_qbz(2,:,ispden,ipert) = v1r_qbz(2,:,ispden,ipert) - v1r_lr(2,:)
         end do
       end if

       ! Compute FT phases for this qpt_bz.
       call calc_eiqr(-qpt_bz, db%my_nrpt, db%my_rpt, emiqr)
       !call cwtime_report(" phases", cpu, wall, gflops)

       ! SLOW FT.
       cnt = 0
       do ispden=1,db%nspden
         do ifft=1,nfft
           cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.

           v1scf_rpt(1,:,ifft,ispden) = v1scf_rpt(1,:,ifft,ispden) &
                                   + emiqr(1,:) * v1r_qbz(1, ifft, ispden, ipert) &
                                   - emiqr(2,:) * v1r_qbz(2, ifft, ispden, ipert)

           v1scf_rpt(2,:,ifft,ispden) = v1scf_rpt(2,:,ifft,ispden) &
                                   + emiqr(1,:) * v1r_qbz(2, ifft, ispden, ipert) &
                                   + emiqr(2,:) * v1r_qbz(1, ifft, ispden, ipert)
         end do
       end do
       !call cwtime_report(" slow fft", cpu, wall, gflops)
     end if

   end do ! iqst

   write(msg,'(2(a,i0),a)') " q-point [",iq_ibz,"/",nqibz,"]"
   call cwtime_report(msg, cpu, wall, gflops)

   ABI_FREE(v1r_qibz)
 end do ! iq_ibz
 ABI_CHECK(iqst == nqbz, "iqst /= nqbz")

 v1scf_rpt = v1scf_rpt / nqbz
 call xmpi_sum(v1scf_rpt, comm, ierr)

 ABI_FREE(emiqr)
 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(qbz)
 ABI_FREE(indqq)
 ABI_FREE(iperm)
 ABI_FREE(bz2ibz_sort)
 ABI_FREE(iqs_dvdb)
 ABI_FREE(nqsts)
 ABI_FREE(v1r_qbz)
 ABI_FREE(v1r_lr)

end subroutine dvdb_get_v1scf_rpt
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_v1scf_qpt
!! NAME
!!  dvdb_get_v1scf_qpt
!!
!! FUNCTION
!!  Fourier interpolation of potentials for a given q-point
!!  This routine is meant to replace dvdb_ftinterp_qpt
!!  by performing the interpolation one perturbation at a time.
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nrpt=Number of R-points = number of q-points in the full BZ
!!  nspden=Number of spin densities.
!!  ipert=index of the perturbation to be treated [1,natom3]
!!  v1scf_rpt(2,nrpt,nfft,nspden)=phonon perturbation potential in real space lattice representation.
!!  comm=MPI communicator
!!
!! OUTPUT
!!  v1scf_qpt(2*nfft, nspden)=Interpolated DFPT potentials at the given q-point.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_get_v1scf_qpt(db, cryst, qpt, nfft, ngfft, nrpt, nspden, &
&                             ipert, v1scf_rpt, v1scf_qpt, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nrpt,nspden,ipert,comm
 class(dvdb_t),intent(in) :: db
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(in) :: v1scf_rpt(2,nrpt,nfft,db%nspden)
 real(dp),intent(out) :: v1scf_qpt(2,nfft,db%nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex2=2
 integer :: ir,ispden,ifft,idir,iat,timerev_q,nproc,my_rank,cnt,ierr
 real(dp) :: wr,wi
!arrays
 integer :: symq(4,2,db%cryst%nsym)
 real(dp),allocatable :: eiqr(:,:), v1r_lr(:,:)

! *************************************************************************

 ABI_UNUSED(cryst%natom)
 ABI_UNUSED(nspden)

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ABI_MALLOC(v1r_lr, (2,nfft))

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(db%cryst%nsym, qpt, symq, db%cryst%symrec, db%cryst%symafm, timerev_q, prtvol=db%prtvol)

 ! Compute FT phases for this q-point.
 ABI_MALLOC(eiqr, (2, db%my_nrpt))
 call calc_eiqr(qpt, db%my_nrpt, db%my_rpt, eiqr)

 idir = mod(ipert-1, 3) + 1; iat = (ipert - idir) / 3 + 1

 ! Compute long-range part of the coupling potential
 v1r_lr = zero; cnt = 0
 if (db%add_lr > 0) call db%get_v1r_long_range(qpt, idir, iat, nfft, ngfft, v1r_lr)

 ! TODO: If high-symmetry q-points, one could save flops by FFT interpolating the independent
 ! TODO: Use ZGEMM with MPI
 ! perturbations and then rotate ...
 v1scf_qpt = zero; cnt = 0
 do ispden=1,db%nspden
   do ifft=1,nfft
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI-parallelism

     do ir=1,db%my_nrpt
       wr = v1scf_rpt(1, ir,ifft,ispden)
       wi = v1scf_rpt(2, ir,ifft,ispden)
       v1scf_qpt(1,ifft,ispden) = v1scf_qpt(1,ifft,ispden) + wr*eiqr(1,ir) - wi * eiqr(2,ir)
       v1scf_qpt(2,ifft,ispden) = v1scf_qpt(2,ifft,ispden) + wr*eiqr(2,ir) + wi * eiqr(1,ir)
     end do

     ! Add the long-range part of the potential
     if (db%add_lr > 0) then
       v1scf_qpt(1,ifft,ispden) = v1scf_qpt(1,ifft,ispden) + v1r_lr(1,ifft)
       v1scf_qpt(2,ifft,ispden) = v1scf_qpt(2,ifft,ispden) + v1r_lr(2,ifft)
     end if
     if (db%add_lr == 4) then
       v1scf_qpt(1,ifft,ispden) = v1r_lr(1,ifft)
       v1scf_qpt(2,ifft,ispden) = v1r_lr(2,ifft)
     end if
   end do ! ifft

   call xmpi_sum(v1scf_qpt(:,:,ispden), comm, ierr)

   ! Remove the phase.
   call times_eikr(-qpt, ngfft, nfft, 1, v1scf_qpt(:,:,ispden))
 end do

 ! Be careful with gamma and cplex!
 if (db%symv1==1) then
   call v1phq_symmetrize(db%cryst, idir, iat, symq, ngfft, cplex2, nfft, db%nspden, db%nsppol, db%mpi_enreg, v1scf_qpt)
 end if

 ABI_FREE(eiqr)
 ABI_FREE(v1r_lr)

end subroutine dvdb_get_v1scf_qpt
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_interpolate_v1scf
!! NAME
!!  dvdb_interpolate_v1scf
!!
!! FUNCTION
!!  Interpolate the phonon perturbation potential.
!!  This routine is meant to replace dvdb_ftinterp_setup and dvdb_ftinterp_qpt.
!!  It performs the interpolation one perturbation at a time.
!!
!! INPUTS
!!  ngqpt(3)=Divisions of the ab-initio q-mesh.
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nfftf=Number of fft-points on the fine grid for interpolated potential
!!  ngfftf(18)=information on 3D FFT for interpolated potential
!!  comm=MPI communicator
!!
!! OUTPUT
!!  v1scf(2, nfft, nspden, 3*natom)= v1scf potentials on the real-space FFT mesh for the 3*natom perturbations.
!!
!! PARENTS
!!      m_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_interpolate_v1scf(db, cryst, qpt, ngqpt, nqshift, qshift, &
&                                 nfft, ngfft, nfftf, ngfftf, v1scf, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,nfftf,comm
 class(dvdb_t),target,intent(inout) :: db
!arrays
 real(dp),intent(in) :: qpt(3)
 integer,intent(in) :: ngqpt(3),ngfft(18),ngfftf(18)
 real(dp),intent(in) :: qshift(3,nqshift)
 real(dp),allocatable,intent(out) :: v1scf(:,:,:,:)
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer :: ipert, nqbz, ierr, nproc, my_rank
 !real(dp) :: work_size
!arrays
 real(dp),allocatable :: v1scf_rpt(:,:,:,:)

! *************************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 nqbz = product(ngqpt) * nqshift
 db%my_nrpt = nqbz

 ABI_MALLOC_OR_DIE(v1scf, (2,nfftf,db%nspden,db%natom3), ierr)
 ABI_MALLOC_OR_DIE(v1scf_rpt, (2,db%my_nrpt,nfft,db%nspden), ierr)

 do ipert=1,db%natom3
   write(std_out, "(a,i4,a,i4,a)") " Interpolating potential for perturbation ", ipert, " / ", db%natom3, ch10

   ! FIXME I think this should be ngfftf and not ngfft
   !       Also, other calls to dvdb_ftinterp_setup should use ngfftf.
   call dvdb_get_v1scf_rpt(db, cryst, ngqpt, nqshift, qshift, nfft, ngfft, &
                           db%my_nrpt, db%nspden, ipert, v1scf_rpt, comm)

   call dvdb_get_v1scf_qpt(db, cryst, qpt, nfftf, ngfftf, db%my_nrpt, db%nspden, &
                           ipert, v1scf_rpt, v1scf(:,:,:,ipert), comm)

   ABI_FREE(db%my_rpt)
 end do

 ABI_FREE(v1scf_rpt)

end subroutine dvdb_interpolate_v1scf
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_findq
!! NAME
!!  dvdb_findq
!!
!! FUNCTION
!!  Find the index of the q-point in db%qpts. Non zero umklapp vectors are not allowed.
!!  Returns -1 if not found.
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!  [qtol]=Optional tolerance for q-point comparison.
!!         For each reduced direction the absolute difference between the coordinates must be less that qtol
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function dvdb_findq(db, qpt, qtol) result(iqpt)

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: qtol
 class(dvdb_t),intent(in) :: db
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer :: iq
 real(dp) :: my_qtol

! *************************************************************************

 my_qtol = tol6; if (present(qtol)) my_qtol = qtol

 iqpt = -1
 do iq=1,db%nqpt
   if (all(abs(db%qpts(:, iq) - qpt) < my_qtol)) then
      iqpt = iq; exit
   end if
 end do

end function dvdb_findq
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_find_qpts
!! NAME
!!  dvdb_find_qpts
!!
!! FUNCTION
!!  Find the index of the q-point in db%qpts. Non zero umklapp vectors are not allowed.
!!  Returns -1 if not found.
!!
!! INPUTS
!!  nqpt: Number of q-points
!!  qpt(3,nqpt): q-point in reduced coordinates.
!!  comm: MPI communicator
!!
!! OUTPUT
!!  iq2dvdb(nqpt): index of q-points in dvdb%qpts. Set to -1 if not found
!!  ierr= Number of points **not** found
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_find_qpts(db, nqpt, qpts, iq2dvdb, comm) result(notfound)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(in) :: db
 integer,intent(in) :: nqpt, comm
!arrays
 real(dp),intent(in) :: qpts(3, nqpt)
 integer,intent(out) :: iq2dvdb(nqpt)

!Local variables-------------------------------
!scalars
 integer :: iq, my_rank, nprocs, ierr

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 iq2dvdb = 0
 do iq=1,nqpt
   if (mod(iq, nprocs) /= my_rank) cycle ! MPI parallelism
   iq2dvdb(iq) = db%findq(qpts(:, iq))
 end do

 call xmpi_sum(iq2dvdb, comm, ierr)
 notfound = count(iq2dvdb == -1)

end function dvdb_find_qpts
!!***


!!****f* m_dvdb/dvdb_set_pert_distrib
!! NAME
!!  dvdb_set_pert_distrib
!!
!! FUNCTION
!!  Activate MPI distribution of the 3*natom perturbations.
!!
!! INPUTS
!!  my_npert=Number of perturbations treated by this rank
!!  natom3= 3 * natom
!!  my_pinfo(3, my_npert)
!!     my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
!!     my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
!!     my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3
!!  pert_table(2, natom3)
!!      pert_table(1, npert): rank of the processor treating this atomic perturbation.
!!      pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.
!!  comm_pert=MPI communicator used to distribute the 3*natom perturbations
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_set_pert_distrib(self, my_npert, natom3, my_pinfo, pert_table, comm_pert)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(inout) :: self
 integer,intent(in) :: my_npert, natom3, comm_pert
!arrays
 integer,intent(in) :: my_pinfo(3,my_npert), pert_table(2,natom3)

! *************************************************************************

 self%comm_pert = comm_pert
 self%nprocs_pert = xmpi_comm_size(comm_pert)
 self%me_pert = xmpi_comm_rank(comm_pert)
 self%my_npert = my_npert

 ABI_SFREE(self%my_pinfo)
 ABI_SFREE(self%pert_table)
 call alloc_copy(my_pinfo, self%my_pinfo)
 call alloc_copy(pert_table, self%pert_table)

 if (self%debug) then
   write(std_out, *)"Activating perturbation over perturbations:"
   write(std_out, *)"nprocs_pert: ", self%nprocs_pert
   write(std_out, *)"my_pinfo: ",self%my_pinfo
   write(std_out, *)"pert_table: ",self%pert_table
 end if

end subroutine dvdb_set_pert_distrib
!!***

!!****f* m_dvdb/dvdb_seek
!! NAME
!!  dvdb_seek
!!
!! FUNCTION
!!  Move the internal file pointer so that it points to the
!!  block with (idir, ipert, iqpt). Needed only if dvdb%iomode==IO_MODE_FORTRAN
!!
!! INPUTS
!!   idir,ipert,iqpt = (direction, perturbation, q-point) indices
!!
!! SIDE EFFECTS
!!   db<type(dvdb_t)>: modifies db%current_fpos.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_seek(db, idir, ipert, iqpt)

!Arguments ------------------------------------
 integer,intent(in)  :: idir, ipert, iqpt
 type(dvdb_t),intent(inout) :: db

!Local variables-------------------------------
 integer :: pos_now,pos_wanted,ii,ispden,nn,ierr
 real(dp),parameter :: fake_qpt(3)=zero
 character(len=500) :: msg

! *************************************************************************

 if (db%iomode == IO_MODE_FORTRAN) then
   pos_now = db%current_fpos
   pos_wanted = db%pos_dpq(idir,ipert,iqpt)
   ABI_CHECK(pos_wanted /= 0, "pos_wanted cannot be zero!")

   ! Optimal access.
   if (pos_now == pos_wanted) return

   if (pos_wanted < pos_now) then
     ! Backspace previous records and header
     ! but only if nn <= pos_wanted else rewind file and skip pos_wanted potentials (should be faster)
     nn = pos_now - pos_wanted
     if (nn <= pos_wanted) then
       do ii=1,nn
         !write(std_out, *)"backspacing"
         if (db%version > 1) backspace(unit=db%fh, err=10, iomsg=msg)
         do ispden=1,db%nspden
           backspace(unit=db%fh, err=10, iomsg=msg)
         end do
         ierr = db%hdr_ref%backspace(db%fh, msg)
         if (ierr /= 0) goto 10
       end do
       db%current_fpos = pos_wanted; return
     else
       ! rewind the file and read it from the beginning
       if (dvdb_rewind(db, msg) /= 0) then
         MSG_ERROR(msg)
       end if
       nn = pos_wanted
     end if

   else
     nn = pos_wanted - pos_now + 1
   end if

   do ii=1,nn-1
     !write(std_out,*)"in seek with ii: ",ii,"pos_wanted: ",pos_wanted
     if (my_hdr_skip(db%fh, -1, -1, fake_qpt, msg) /= 0) then
       MSG_ERROR(msg)
     end if
     ! Skip the records with v1.
     do ispden=1,db%nspden
       read(db%fh, err=10, iomsg=msg)
     end do
     ! Skip record with rhog1_g0 (if present)
     if (db%version > 1) read(db%fh, err=10, iomsg=msg)
   end do

   db%current_fpos = pos_wanted

 else
   MSG_ERROR("Should not be called when iomode /= IO_MODE_FORTRAN")
 end if

 return

 ! Handle Fortran IO error
10 continue
 msg = sjoin("Error while reading", db%path, ch10, msg)

end subroutine dvdb_seek
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_rewind
!! NAME
!!  dvdb_rewind
!!
!! FUNCTION
!!   Rewind the file and move to the first header. Needed only if dvdb%iomode==IO_MODE_FORTRAN
!!   Return exit code and error message in msg if ierr != 0
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_rewind(db, msg) result(ierr)

!Arguments ------------------------------------
 type(dvdb_t),intent(inout) :: db
 character(len=*),intent(out) :: msg

! *************************************************************************

 ierr = 0
 if (db%iomode == IO_MODE_FORTRAN) then
   rewind(db%fh, err=10, iomsg=msg)
   read(db%fh, err=10, iomsg=msg)  ! version
   read(db%fh, err=10, iomsg=msg)  ! numv1
   db%current_fpos = 1

 else
   ierr = -1
   msg = "should not be called when iomode /= IO_MODE_FORTRAN"
 end if

 return

 ! Handle Fortran IO error
10 continue
 ierr = 1
 msg = sjoin("Error while reading", db%path, ch10, msg)

end function dvdb_rewind
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/my_hdr_skip
!! NAME
!!  my_hdr_skip
!!
!! FUNCTION
!!  Skip the header without rewinding the file. Return exit code.
!!
!! NOTES
!!  Because hdr_skip rewinds the file and I'm not gonna change that ugly code.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function my_hdr_skip(unit, idir, ipert, qpt, msg) result(ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit,idir,ipert
 real(dp),intent(in) :: qpt(3)
 character(len=500),intent(out) :: msg

!Local variables-------------------------------
 integer :: fform
 type(hdr_type) :: tmp_hdr
!************************************************************************

 ierr = 0; msg = ""
 call hdr_fort_read(tmp_hdr, unit, fform)
 ierr = dvdb_check_fform(fform, "read_dvdb", msg)
 if (ierr /= 0) return

 if (idir /= -1 .and. ipert /= -1) then
   if (idir /= mod(tmp_hdr%pertcase-1, 3) + 1 .or. &
       ipert /= (tmp_hdr%pertcase - idir) / 3 + 1 .or. &
       any(abs(qpt - tmp_hdr%qptn) > tol14)) then
         msg = "Perturbation index on file does not match the one expected by the caller"
         ierr = -1
   end if
 end if

 call tmp_hdr%free()

end function my_hdr_skip
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_list_perts
!! NAME
!!  dvdb_list_perts
!!
!! FUNCTION
!!  Given a q-point mesh, this routine checks if all the (phonon) perturbations
!!  are available taking into account symmetries.
!!
!! INPUTS
!!  ngqpt(3)=Q-mesh divisions. If all(ngqpt == -1), the list of q-points in the DVDB
!!    (i.e. db%qpts) is analyzed instead of the q-points generated from ngqpt.
!!  [unit]=Unit number for output. Default `std_out`.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_list_perts(db, ngqpt, unit)

!Arguments ------------------------------------
 class(dvdb_t),target,intent(in) :: db
 integer,optional,intent(in) :: unit
!arrays
 integer,intent(in) :: ngqpt(3)

!Local variables-------------------------------
!scalars
 integer :: tot_miss,tot_weird,miss_q,idir,ipert,iv1,psy,weird_q,enough
 integer :: iq_ibz,nqibz,iq_file,qptopt,nshiftq,ii,timerev_q,unt,nqbz
 character(len=500) :: msg,ptype,found
 type(crystal_t),pointer :: cryst
!arrays
 integer :: rfdir(3),qptrlatt(3,3)
 integer,allocatable :: pertsy(:,:),symq(:,:,:),rfpert(:)
 real(dp) :: qq(3),shiftq(3,1)
 real(dp),allocatable :: qibz(:,:),wtq(:),qbz(:,:)

! *************************************************************************

 unt = std_out; if (present(unit)) unt = unit
 cryst => db%cryst

 if (all(ngqpt == -1)) then
   ! Will test the q-points in db
   call alloc_copy(db%qpts, qibz)
   nqibz = db%nqpt
 else
   ! Will test the q-points in the IBZ associated to ngqpt hence build IBZ and BZ from ngqpt.
   qptopt = 1; shiftq = zero; nshiftq = 1; qptrlatt = 0
   do ii=1,3
     qptrlatt(ii, ii) = ngqpt(ii)
   end do

   call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt, nshiftq, shiftq, &
     nqibz, qibz, wtq, nqbz, qbz)

   ABI_FREE(qbz)
   ABI_FREE(wtq)
 end if

 ! Initialize the list of perturbations rfpert and rdfir
 ! WARNING: Only phonon perturbations are considered for the time being.
 ABI_MALLOC(rfpert,(db%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1

 ABI_MALLOC(symq, (4,2,cryst%nsym))
 ABI_MALLOC(pertsy, (3,db%mpert))

 ! Loop over the q-points in the IBZ and test whether the q-point is present
 ! and if all the independent perturbations are available.
 !   `tot_miss` is the number of irreducible perturbations not found in the DVDB (critical)
 !   `tot_weird` is the number of redundant perturbations found in the DVDB (not critical)
 enough = 5; if (db%prtvol > 0) enough = nqibz + 1
 tot_miss = 0; tot_weird = 0
 do iq_ibz=1,nqibz
   if (iq_ibz == enough)  then
     call wrtout(unt,' More than 20 q-points with prtvol == 0. Only important messages will be printed...')
   end if
   qq = qibz(:,iq_ibz)
   iq_file = db%findq(qq)

   ! Examine the symmetries of the q wavevector
   call littlegroup_q(cryst%nsym,qq,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=db%prtvol)

   ! Determine the symmetrical perturbations. Meaning of pertsy:
   !    0 for non-target perturbations
   !    1 for basis perturbations
   !   -1 for perturbations that can be found from basis perturbations
   call irreducible_set_pert(cryst%indsym,db%mpert,cryst%natom,cryst%nsym,&
     pertsy,rfdir,rfpert,symq,cryst%symrec,cryst%symrel)

   if (iq_file /= -1) then
     ! This q-point is in the DVDB. Test if all the independent perturbations are available.
     if (iq_ibz <= enough)  then
       call wrtout(unt, sjoin("qpoint:", ktoa(qq), "is present in the DVDB file"))
       call wrtout(unt,' The list of irreducible perturbations for this q vector is:')
     end if
     ii = 0; weird_q = 0; miss_q = 0
     do ipert=1,db%mpert
       do idir=1,3
         psy = pertsy(idir,ipert)
         if (psy == 0) cycle
         iv1 = db%pos_dpq(idir,ipert,iq_file)
         ptype = "independent"; if (psy == -1) ptype = "symmetric"
         found = "Yes"; if (iv1 == 0) found = "No"

         if (psy == 1 .and. iv1 == 0) miss_q = miss_q + 1
         if (psy == -1 .and. iv1 /= 0) weird_q = weird_q + 1

         ii=ii+1
         if (iq_ibz <= enough)  then
           write(msg,'(i5,a,i2,a,i4,4a)')ii,')  idir=',idir,', ipert=',ipert,", type=",trim(ptype),", found=",trim(found)
           call wrtout(unt, msg)
         end if
       end do
     end do

     if (weird_q /= 0) then
       write(msg,"(a,i0,a)")"DVDB is overcomplete. ",weird_q, " perturbation(s) can be reconstructed by symmetry."
       call wrtout(unt, msg)
     end if

     tot_weird = tot_weird + weird_q
     tot_miss = tot_miss + miss_q
     if (miss_q /=0) then
       call wrtout(unt, sjoin("WARNING:", itoa(miss_q), "independent perturbation(s) are missing!."))
     end if

   else
     ! This q-point is not present in dvdb. Print the list of independent perturbations.
     call wrtout(unt, sjoin("qpoint:", ktoa(qq), "is NOT present in the DVDB file"))
     call wrtout(unt,' The list of irreducible perturbations for this q vector is:')
     ii = 0
     do ipert=1,db%mpert
       do idir=1,3
         if (pertsy(idir,ipert) == 1) then
           ii=ii+1
           write(msg,'(i5,a,i2,a,i4,a)')ii,')  idir=',idir,', ipert=',ipert,", type=independent, found=No"
           call wrtout(unt, msg)
           tot_miss = tot_miss + 1
         end if
       end do
     end do
   end if

   if (iq_ibz <= enough) call wrtout(unt," ")
 end do ! iq_ibz

 if (tot_miss /= 0) then
   call wrtout(unt, sjoin(ch10, "There are ",itoa(tot_miss), "independent perturbations missing!"))
 else
   call wrtout(unt, "All the independent perturbations are available")
   if (tot_weird /= 0) then
     call wrtout(unt, "Note however that the DVDB is overcomplete as symmetric perturbations are present.")
   end if
 end if

 ABI_FREE(qibz)
 ABI_FREE(rfpert)
 ABI_FREE(symq)
 ABI_FREE(pertsy)

end subroutine dvdb_list_perts
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_merge_files
!! NAME
!!  dvdb_merge_files
!!
!! FUNCTION
!!  Merge a list of POT1 files.
!!
!! INPUT
!!  nfiles=Number of files to be merged.
!!  dvdb_path=Name of output DVDB file.
!!  prtvol=Verbosity level.
!!
!! SIDE EFFECTS
!!   v1files=List of file names to merge. This list could be changed if POT1 files in netcdf format are found.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_merge_files(nfiles, v1files, dvdb_path, prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfiles,prtvol
 character(len=*),intent(in) :: dvdb_path
 character(len=*),intent(inout) :: v1files(nfiles)

!Local variables-------------------------------
!scalars
! Here I made a mistake because 102 corresponds to GS potentials
! as a consequence DVDB files generated with version <= 8.1.6
! contain list of potentials with fform = 102.
 !integer :: fform_pot=102
 integer :: fform_pot=111
 integer :: ii,jj,fform,ount,cplex,nfft,ifft,ispden,nperts
 integer :: n1,n2,n3,v1_varid,ierr
 logical :: qeq0
 character(len=500) :: msg
 type(hdr_type),pointer :: hdr1
 type(dvdb_t) :: dvdb
!arrays
 integer :: units(nfiles)
 real(dp) :: rhog1_g0(2)
 real(dp),allocatable :: v1(:)
 logical :: has_rhog1_g0(nfiles)
 type(hdr_type),target,allocatable :: hdr1_list(:)

!************************************************************************

 if (file_exists(dvdb_path)) then
   MSG_ERROR(sjoin("Cannot overwrite existing file:", dvdb_path))
 end if

 ! If a file is not found, try the netcdf version and change v1files accordingly.
 do ii=1,nfiles
   if (nctk_try_fort_or_ncfile(v1files(ii), msg) /= 0) then
     MSG_ERROR(msg)
   end if
 end do

 ! Read the headers
 ABI_MALLOC(hdr1_list, (nfiles))
 nperts = size(hdr1_list)

 ! Write dvdb file (we only support fortran binary format)
 if (open_file(dvdb_path, msg, newunit=ount, form="unformatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 write(ount, err=10, iomsg=msg) dvdb_last_version
 write(ount, err=10, iomsg=msg) nperts

 ! Validate headers.
 ! TODO: Should perform consistency check on the headers
 ! rearrange them in blocks of q-points for efficiency reason.
 ! ignore POT1 files that do not correspond to atomic perturbations.
 ! Add POT file from GS run to support Sternheimer in eph_task = 4

 do ii=1,nfiles
   write(std_out,"(a,i0,2a)")"- Reading header of file [",ii,"]: ",trim(v1files(ii))

   if (endswith(v1files(ii), ".nc")) then
#ifdef HAVE_NETCDF
      NCF_CHECK(nctk_open_read(units(ii), v1files(ii), xmpi_comm_self))
      call hdr_ncread(hdr1_list(ii),units(ii),fform)
#endif
   else
     if (open_file(v1files(ii), msg, newunit=units(ii), form="unformatted", action="read", status="old") /= 0) then
       MSG_ERROR(msg)
     end if
     call hdr_fort_read(hdr1_list(ii), units(ii), fform)
   end if

   if (dvdb_check_fform(fform, "merge_dvdb", msg) /= 0) then
     MSG_ERROR(sjoin("While reading:", v1files(ii), msg))
   end if
   if (prtvol > 0) call hdr1_list(ii)%echo(fform, 3, unit=std_out)
   if (hdr1_list(ii)%pertcase == 0) then
     MSG_ERROR(sjoin("Found GS potential:", v1files(ii)))
   end if
   !write(std_out,*)"done", trim(v1files(ii))

   ! Supported fform:
   ! 109  POT1 files without vh1(G=0)
   ! 111  POT1 files with extra record with vh1(G=0) after FFT data.
   has_rhog1_g0(ii) = .True.
   if (fform == 109) has_rhog1_g0(ii) = .False.

   write(std_out,"(a,i0,2a)")"- Merging file [",ii,"]: ",trim(v1files(ii))
   jj = ii
   hdr1 => hdr1_list(jj)
   call hdr1%fort_write(ount, fform_pot, ierr)
   ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

   qeq0 = (hdr1%qptn(1)**2+hdr1%qptn(2)**2+hdr1%qptn(3)**2<1.d-14)
   cplex = 2; if (qeq0) cplex = 1
   nfft = product(hdr1%ngfft(1:3))
   n1 = hdr1%ngfft(1); n2 = hdr1%ngfft(2); n3 = hdr1%ngfft(3)

   ABI_MALLOC(v1, (cplex*nfft))

   if (.not. endswith(v1files(ii), ".nc")) then
      ! Fortran IO
      do ispden=1,hdr1%nspden
        read(units(jj), err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfft)
        write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfft)
      end do
      ! Add rhog1(G=0)
      rhog1_g0 = zero
      if (has_rhog1_g0(jj)) read(units(jj), err=10, iomsg=msg) rhog1_g0
      if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0
   else
#ifdef HAVE_NETCDF
      ! Netcdf IO
      ! netcdf array has shape [cplex, n1, n2, n3, nspden]
      NCF_CHECK(nf90_inq_varid(units(ii), "first_order_potential", v1_varid))
      do ispden=1,hdr1%nspden
        NCF_CHECK(nf90_get_var(units(ii), v1_varid, v1, start=[1,1,1,ispden], count=[cplex, n1, n2, n3, 1]))
        write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfft)
      end do
      ! Add rhog1(G=0)
      rhog1_g0 = zero
      if (has_rhog1_g0(jj)) then
        NCF_CHECK(nf90_get_var(units(ii), nctk_idname(units(ii), "rhog1_g0"), rhog1_g0))
      end if
      if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0
#endif
   end if

   if (.not. endswith(v1files(ii), ".nc")) then
     close(units(ii))
   else
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(units(ii)))
#endif
   end if

   ABI_FREE(v1)
 end do ! nperts

 close(ount)

 do ii=1,size(hdr1_list)
   call hdr1_list(ii)%free()
 end do
 ABI_DT_FREE(hdr1_list)

 write(std_out,"(a,i0,a)")"Merged successfully ", nfiles, " files"

 ! List available perturbations.
 dvdb = dvdb_new(dvdb_path, xmpi_comm_self)
 call dvdb%print()
 call dvdb%list_perts([-1, -1, -1])
 call dvdb%free()

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(sjoin("Error while merging files", ch10, msg))

end subroutine dvdb_merge_files
!!***

!!****f* m_dvdb/calc_eiqr
!! NAME
!!  calc_eiqr
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine calc_eiqr(qpt, nrpt, rpt, eiqr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nrpt
!arrays
 real(dp),intent(in) :: qpt(3),rpt(3,nrpt)
 real(dp),intent(out) :: eiqr(2,nrpt)

!Local variables -------------------------
!scalars
 integer :: ir
 real(dp) :: qr

! *********************************************************************

 do ir=1,nrpt
   qr = two_pi * dot_product(qpt, rpt(:,ir))
   eiqr(1,ir) = cos(qr); eiqr(2,ir) = sin(qr)
 end do

end subroutine calc_eiqr
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_check_fform
!! NAME
!!  dvdb_check_fform
!!
!! FUNCTION
!!  Check the value of fform. Return exit status and error message.
!!
!! INPUTS
!!   fform=Value read from the header
!!   mode="merge_dvdb" to check the value of fform when we are merging POT1 files
!!        "read_dvdb" when we are reading POT1 files from a DVDB file.
!!
!! OUTPUT
!!   errmsg=String with error message if ierr /= 0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_check_fform(fform, mode, errmsg) result(ierr)

!Arguments ------------------------------------
 integer,intent(in) :: fform
 character(len=*),intent(in) :: mode
 character(len=*),intent(out) :: errmsg

! *************************************************************************
 ierr = 0

! Here I made a mistake because 102 corresponds to GS potentials
! as a consequence DVDB files generated with version <= 8.1.6
! contain list of potentials with fform = 102.
 !integer :: fform_pot=102
 !integer :: fform_pot=109
 !integer :: fform_pot=111

 if (fform == 0) then
   errmsg = "fform == 0! Either corrupted/nonexistent file or IO error"
   ierr = 42; return
 end if

 select case (mode)
 case ("merge_dvdb")
    if (all(fform /= [109, 111])) then
      errmsg = sjoin("fform:", itoa(fform), "is not supported in `merge_dvdb` mode")
      ierr = 1; return
    end if

 case ("read_dvdb")
    if (all(fform /= [102, 109, 111])) then
      errmsg = sjoin("fform:", itoa(fform), "is not supported in `read_dvdb` mode")
      ierr = 1; return
    end if

 case default
   errmsg = sjoin("Invalid mode:", mode)
   ierr = -1; return
 end select

end function dvdb_check_fform
!!***

!!****f* m_dvdb/dvdb_test_v1rsym
!! NAME
!!  dvdb_test_v1rsym
!!
!! FUNCTION
!!  Debugging tool used to check whether the DFPT potentials in real space fulfill
!!  the correct symmetries on the real space FFT mesh.
!!
!! INPUTS
!!  db_path=Filename
!!  symv1scf=1 to activate symmetrization of DFPT potentials. 0 to disable it.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_test_v1rsym(db_path, symv1scf, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: db_path
 integer,intent(in) :: symv1scf, comm

!Local variables-------------------------------
!scalars
 integer,parameter :: rfmeth2=2,syuse0=0
 integer :: iqpt,idir,ipert,nsym1,cplex,v1pos
 integer :: isym,nfft,ifft,ifft_rot,ispden
 real(dp) :: max_err,re,im,vre,vim !,pre,pim
 character(len=500) :: msg
 logical :: isok
 type(dvdb_t),target :: db
 type(crystal_t),pointer :: cryst
!arrays
 integer :: ngfft(18)
 integer,allocatable :: symafm1(:),symrel1(:,:,:),irottb(:,:)
 real(dp) :: qpt(3)
 real(dp),allocatable :: tnons1(:,:),v1scf(:,:)

! *************************************************************************

 db = dvdb_new(db_path, comm)
 db%debug = .True.
 db%symv1 = symv1scf
 call db%print()
 !call db%list_perts([-1,-1,-1])

 call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
 nfft = product(ngfft(1:3))
 call db%open_read(ngfft, comm)

 ABI_CHECK(db%nspinor==1, "nspinor == 2 not coded")

 cryst => db%cryst
 ABI_MALLOC(symafm1, (cryst%nsym))
 ABI_MALLOC(symrel1, (3,3,cryst%nsym))
 ABI_MALLOC(tnons1, (3,cryst%nsym))

 do iqpt=1,db%nqpt
   qpt = db%qpts(:, iqpt)
   do ipert=1,db%natom
     do idir=1,3
       v1pos = db%pos_dpq(idir, ipert, iqpt); if (v1pos == 0) cycle

       ! Determines the set of symmetries that leaves the perturbation invariant.
       call littlegroup_pert(cryst%gprimd,idir,cryst%indsym,dev_null,ipert,cryst%natom,cryst%nsym,nsym1,rfmeth2,&
         cryst%symafm,symafm1,db%symq_table(:,:,:,iqpt),cryst%symrec,cryst%symrel,symrel1,syuse0,cryst%tnons,&
         tnons1,unit=dev_null)

       cplex = db%cplex_v1(v1pos)
       ngfft(1:3) = db%ngfft3_v1(:, v1pos)
       nfft = product(ngfft(:3))
       ABI_MALLOC(v1scf, (cplex*nfft, db%nspden))

       if (db%read_onev1(idir, ipert, iqpt, cplex, nfft, ngfft, v1scf, msg) /= 0) then
         MSG_ERROR(msg)
       end if

       ABI_MALLOC(irottb, (nfft,nsym1))
       call rotate_fft_mesh(nsym1,symrel1,tnons1,ngfft,irottb,isok)
       if (.not. isok) then
         MSG_WARNING("Real space FFT mesh is not compatible with symmetries!")
       end if

       max_err = zero
       do isym=1,nsym1
         do ispden=1,db%nspden
           do ifft=1,nfft
             ifft_rot = irottb(ifft, isym)
             !pre =  cos(two_pi * dot_product(qpt, tnons1(:,isym)))
             !pim = -sin(two_pi * dot_product(qpt, tnons1(:,isym)))
             if (cplex == 2) then
               re = v1scf(2*ifft_rot-1, ispden)
               im = v1scf(2*ifft_rot  , ispden)
               !vre = re * pre - im * pim
               !vim = re * pim + im * pre
               vre = re; vim = im

               re = v1scf(2*ifft-1, ispden) - vre
               im = v1scf(2*ifft  , ispden) - vim
             else
               re = v1scf(ifft, ispden) - v1scf(ifft_rot, ispden)
               im = zero
             end if
             !if (sqrt(re**2 + im**2) > tol6) write(std_out,*)"ifft,isym,err: ",ifft,isym,sqrt(re**2 + im**2)
             max_err = max(max_err, sqrt(re**2 + im**2))
           end do
         end do
       end do
       if (nsym1>1) then
         write(std_out,"(3(a,i2),a,i2,a,es16.8)")"For iqpt= ",iqpt,&
         ", idir= ",idir,", ipert= ",ipert,", nsym= ",nsym1,", max_err= ",max_err
       end if

       ABI_FREE(irottb)
       ABI_FREE(v1scf)
     end do
   end do

 end do ! iqpt

 ABI_FREE(symafm1)
 ABI_FREE(symrel1)
 ABI_FREE(tnons1)

 call db%free()

end subroutine dvdb_test_v1rsym
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_test_v1complete
!! NAME
!!  dvdb_test_v1complete
!!
!! FUNCTION
!!  Debugging tool used to test the symmetrization of the DFPT potentials.
!!  Assumes DVDB file containing all 3*natom perturbations (generated with nsym == 1 or
!!  other specialized variables e.g. prepgkk)
!!
!! INPUTS
!!  db_path=Filename of the DVDB file.
!!  symv1scf=1 to activate symmetrization of DFPT potentials. 0 ti disable it.
!!  dump_path=File used to dump potentials (empty string to disable output)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_test_v1complete(dvdb_path, symv1scf, dump_path, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: dvdb_path,dump_path
 integer,intent(in) :: symv1scf, comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iqpt,pcase,idir,ipert,cplex,nfft,ispden,timerev_q,ifft,unt,my_rank, ncid
 integer :: i1,i2,i3,n1,n2,n3,id1,id2,id3,cnt
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
 type(dvdb_t),target :: dvdb
!arrays
 integer :: ngfft(18), rfdir(3)
 integer,allocatable :: pflag(:,:)
 real(dp) :: qpt(3)
 integer,allocatable :: pertsy(:,:),rfpert(:),symq(:,:,:)
 real(dp),allocatable :: file_v1scf(:,:,:,:),symm_v1scf(:,:,:,:), work2(:,:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 dvdb = dvdb_new(dvdb_path, comm)
 dvdb%debug = .false.
 dvdb%symv1 = symv1scf
 call dvdb%print()
 call dvdb%list_perts([-1,-1,-1])

 call ngfft_seq(ngfft, dvdb%ngfft3_v1(:,1))
 nfft = product(ngfft(1:3))
 call dvdb%open_read(ngfft, comm)

 cryst => dvdb%cryst
 call ngfft_seq(ngfft, dvdb%ngfft3_v1(:,1))
 nfft = product(ngfft(:3))

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 id1 = n1/2+2; id2 = n2/2+2; id3 = n3/2+2

 ABI_MALLOC(pflag, (3, dvdb%natom))

 ! Initialize the list of perturbations rfpert and rdfir
 ! WARNING: Only phonon perturbations are considered for the time being.
 ABI_MALLOC(rfpert,(dvdb%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1
 ABI_MALLOC(symq, (4,2,cryst%nsym))
 ABI_MALLOC(pertsy, (3,dvdb%mpert))

 unt = -1; ncid = nctk_noid
 if (len_trim(dump_path) /= 0 .and. my_rank == master) then
   write(std_out,"(a)")sjoin("Will write potentials to:", dump_path)
   if (endswith(dump_path, ".nc")) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_create(ncid, dump_path, xmpi_comm_self))
     NCF_CHECK(dvdb%cryst%ncwrite(ncid))
     ncerr = nctk_def_dims(ncid, [&
       nctkdim_t("two", 2), nctkdim_t("three", 3), nctkdim_t("nfft", nfft), nctkdim_t("nspden", dvdb%nspden), &
       nctkdim_t("natom3", cryst%natom * 3), nctkdim_t("mpert", dvdb%mpert), nctkdim_t("nqpt", dvdb%nqpt)], &
       defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "symv1scf"]))
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("qpts", "dp", "three, nqpt")))
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("origin_v1scf", "dp", "two, nfft, nspden, natom3, nqpt")))
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("recons_v1scf", "dp", "two, nfft, nspden, natom3, nqpt")))
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("pertsy_qpt", "int", "three, mpert, nqpt")))
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("ngfft", "int", "three")))
     NCF_CHECK(nctk_set_datamode(ncid))
     ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
       "symv1scf"], [symv1scf])
     NCF_CHECK(ncerr)
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpts"), dvdb%qpts))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngfft"), ngfft(1:3)))
#endif
   else
     if (open_file(dump_path, msg, newunit=unt, action="write", status="unknown", form="formatted") /= 0) then
       MSG_ERROR(msg)
     end if
   end if
 end if

 ABI_CALLOC(work2, (2, nfft, dvdb%nspden, dvdb%natom3))

 do iqpt=1,dvdb%nqpt
   qpt = dvdb%qpts(:,iqpt)
   ! Examine the symmetries of the q wavevector
   call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dvdb%prtvol)

   ! Determine the symmetrical perturbations. Meaning of pertsy:
   !    0 for non-target perturbations
   !    1 for basis perturbations
   !   -1 for perturbations that can be found from basis perturbations
   call irreducible_set_pert(cryst%indsym,dvdb%mpert,cryst%natom,cryst%nsym,&
     pertsy,rfdir,rfpert,symq,cryst%symrec,cryst%symrel)

   ! Read all potentials (here I assume that all perturbations are available)
   call dvdb%readsym_allv1(iqpt, cplex, nfft, ngfft, file_v1scf, dvdb%comm)

   ! Copy basis perturbations in symm_v1scf and set pflag
   ABI_MALLOC(symm_v1scf, (cplex, nfft, dvdb%nspden, dvdb%natom3))
   symm_v1scf = huge(one); pflag = 0
   do pcase=1,3*dvdb%cryst%natom
     idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
     if (pertsy(idir, ipert) == 1) then
       symm_v1scf(:,:,:,pcase) = file_v1scf(:,:,:,pcase)
       pflag(idir,ipert) = 1
     end if
   end do

   ! Complete potentials
   call v1phq_complete(cryst,qpt,ngfft,cplex,nfft,dvdb%nspden,dvdb%nsppol,dvdb%mpi_enreg,dvdb%symv1,pflag,symm_v1scf)

#ifdef HAVE_NETCDF
   if (ncid /= nctk_noid) then
     work2 = zero
     if (cplex == 1) work2(1,:,:,:) = file_v1scf(1,:,:,:)
     if (cplex == 2) work2 = file_v1scf
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "origin_v1scf"), work2, start=[1,1,1,1,iqpt])
     NCF_CHECK(ncerr)
     if (cplex == 1) work2(1,:,:,:) = symm_v1scf(1,:,:,:)
     if (cplex == 2) work2 = symm_v1scf
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "recons_v1scf"), work2, start=[1,1,1,1,iqpt])
     NCF_CHECK(ncerr)
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "pertsy_qpt"), pertsy, start=[1,1,iqpt]))
   end if
#endif

   ! Compare values.
   do pcase=1,3*cryst%natom
     idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
     if (pflag(idir,ipert) /= 2) cycle
     cnt = cnt+1

     do ispden=1,dvdb%nspden
       !write(std_out,"(5(a,i0),3a,es12.4)")"For cnt: ",cnt ,", iqpt: ", iqpt, ", idir: ", idir, &
       !        ", ipert: ", ipert, ", ispden: ", ispden, ", qpt: ", trim(ktoa(qpt)) ,", max_err: ", &
       !        maxval(abs(file_v1scf(:,:,ispden,pcase) - symm_v1scf(:,:,ispden,pcase)))
       !write(std_out,"(a,es10.3)")" max(abs(f1-f2))", maxval(abs(file_v1scf(:,:,ispden,pcase) - symm_v1scf(:,:,ispden,pcase)))
       call vdiff_print(vdiff_eval(cplex,nfft,file_v1scf(:,:,ispden,pcase),symm_v1scf(:,:,ispden,pcase),cryst%ucvol))

       ! Debug: write potentials to file.
       if (unt /= -1) then
         write(unt,*)"# count:", cnt
         write(unt,*)"# q-point:", trim(ktoa(qpt)), ", iqpt: ", trim(itoa(iqpt))
         write(unt,*)"# idir: ",idir,", ipert: ",ipert,", ispden:", ispden
         write(unt,*)"# file_v1scf, symmetrized_v1scf, diff"
         if (cplex == 1) then
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 ifft = i1+n1*((i2-1)+n2*(i3-1))
                 write(unt,"(3i3,3(es12.4,2x))") &
                   i1,i2,i3, &
                   file_v1scf(1,ifft,ispden,pcase), symm_v1scf(1,ifft,ispden,pcase), &
                   file_v1scf(1,ifft,ispden,pcase) - symm_v1scf(1,ifft,ispden,pcase)
               end do
             end do
           end do
         else
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 ifft = i1+n1*((i2-1)+n2*(i3-1))
                 write(unt, "(3i3,6(es12.4,2x))") &
                   i1,i2,i3, &
                   file_v1scf(1,ifft,ispden,pcase), symm_v1scf(1,ifft,ispden,pcase),  &
                   file_v1scf(1,ifft,ispden,pcase) - symm_v1scf(1,ifft,ispden,pcase), &
                   file_v1scf(2,ifft,ispden,pcase), symm_v1scf(2,ifft,ispden,pcase),  &
                   file_v1scf(2,ifft,ispden,pcase) - symm_v1scf(2,ifft,ispden,pcase)
               end do
             end do
           end do
         end if
         write(unt,*)
         write(unt,*)
       end if

     end do
     !write(std_out,*)""
   end do

   ABI_FREE(symm_v1scf)
   ABI_FREE(file_v1scf)
 end do

 ABI_FREE(work2)
 ABI_FREE(pflag)
 ABI_FREE(rfpert)
 ABI_FREE(symq)
 ABI_FREE(pertsy)

 call dvdb%free()

 if (unt /= -1) close(unt)
#ifdef HAVE_NETCDF
 if (ncid /= nctk_noid) then
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

end subroutine dvdb_test_v1complete
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_write_v1qavg
!! NAME
!!  dvdb_write_v1qavg
!!
!! FUNCTION
!!  This routine computes the average over the unit cell of the periodic part of the DFPT potentials
!!  as a function of the q-point and the corresponding quantity obtained with the model for the LR part.
!!  Results are stored in the V1QAVG netcdf file. Two options are available:
!!
!!  eph_task = -15 --> Use list of q-points found in the DVDB file. Mainly used to plot the average
!!    along a q-path. The procedure required to generate a DVDB with a q-path is rather lengthy
!!    as it requires several phonon calculations with WKQ followed by a merge of the POT files.
!!
!!  eph_task = +15 --> Assume DVDB file with q-mesh (dvdb_ngqpt), use Fourier interpolation
!!    to interpolate potentials along the path specified by ph_qpath and ph_nqpath.
!!
!! INPUTS
!!  dtset<dataset_type>= Input variables.
!!  out_ncpath=Filename for output netcdf file.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_write_v1qavg(dvdb, dtset, out_ncpath)

!Arguments ------------------------------------
 class(dvdb_t),target,intent(inout) :: dvdb
 type(dataset_type),target,intent(in) :: dtset
 character(len=*),intent(in) :: out_ncpath

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: nfft, iq, cplex, ispden, comm_rpt, my_rank, idir, ipert, ipc, imyp
 integer :: n1, n2, n3, unt, this_nqpt, method, interpolated
 integer :: i1, i2, i3, ifft, ig, ngsmall, ii
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
 real(dp) :: gsq_max, g2
 !type(vdiff_t) :: vd_max
 character(len=500) :: msg
 character(len=fnlen) :: dump_path
!arrays
 integer :: ngfft(18)
 integer, allocatable :: gfft(:,:),ig2ifft(:), gsmall(:,:)
 real(dp) :: vals2(2)
 real(dp),pointer :: this_qpts(:,:)
 real(dp),allocatable :: file_v1r(:,:,:,:),long_v1r(:,:,:,:),tmp_v1r(:,:,:,:)
 real(dp),allocatable :: maxw(:,:), all_rpt(:,:), all_rmod(:), workg(:,:), work_gsmall(:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(dvdb%comm)

 call wrtout([std_out, ab_out], &
     " Computing average over the unit cell of the periodic part of the DFPT potentials", newlines=2)
 call dvdb%print(unit=std_out)
 !call dvdb%print(unit=ab_out)

 ! Define FFT mesh
 ngfft = dvdb%ngfft
 nfft = product(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 ! Get list of G-vectors in FFT mesh.
 ABI_MALLOC(gfft, (3, nfft))
 call get_gftt(ngfft, [zero, zero, zero], dvdb%cryst%gmet, gsq_max, gfft)
 ABI_MALLOC(workg, (2, nfft))

 ! Select G-vectors in small sphere (ratio of gsq_max)
 do ii=1,2
   if (ii == 2) then
     ABI_MALLOC(ig2ifft, (ngsmall))
   end if
   ngsmall = 0
   do ig=1,nfft
     ! Don't include (2pi)**2 to be consistent with get_gftt
     g2 = dot_product(gfft(:,ig), matmul(dvdb%cryst%gmet, gfft(:, ig)))
     if (g2 <= gsq_max * 0.01_dp) ngsmall = ngsmall + 1
     if (ii == 2) ig2ifft(ngsmall) = ig
   end do
 end do
 write(std_out, *)"Found ngsmall", ngsmall

 !call ig2fft_sphere(dvdb%cryst%gmet, gfft, ig2ifft)
 ABI_MALLOC(gsmall, (3, ngsmall))
 do ig=1,ngsmall
   gsmall(:, ig) = gfft(:, ig2ifft(ig))
 end do
 ABI_MALLOC(work_gsmall, (2, ngsmall))

 ABI_MALLOC(long_v1r, (2, nfft, dvdb%nspden, dvdb%my_npert))
 ABI_MALLOC(file_v1r, (2, nfft, dvdb%nspden, dvdb%my_npert))

 unt = -1; dump_path = ""
 !dump_path = "V1QAVG.dat"
 if (len_trim(dump_path) /= 0 .and. my_rank == master) then
   if (open_file(dump_path, msg, newunit=unt, action="write", status="unknown", form="formatted") /= 0) then
     MSG_ERROR(msg)
   end if
   write(std_out,"(a)")sjoin(" Will write potentials in text format to:", dump_path)
 end if

 ! Select list of q-points depending on eph_task (either from DVDB file or interpolated)
 if (dtset%eph_task == -15) then
   call wrtout([std_out, ab_out], " Using list of q-points found in the DVDB file")
   this_nqpt = dvdb%nqpt
   this_qpts => dvdb%qpts
   interpolated = 0

 else if (dtset%eph_task == +15) then
   msg = sjoin(" Using list of q-points specified by ph_qpath with ", itoa(dtset%ph_nqpath), "qpoints")
   call wrtout([std_out, ab_out], msg)
   ABI_CHECK(dtset%ph_nqpath > 0, "When eph_task = +15, ph_qpath must be given in input.")
   this_nqpt = dtset%ph_nqpath
   this_qpts => dtset%ph_qpath(:, 1:this_nqpt)
   comm_rpt = xmpi_comm_self
   method = dtset%userid
   call dvdb%ftinterp_setup(dtset%dvdb_ngqpt, [1, 1, 1], 1, dtset%ddb_shiftq, nfft, ngfft, method, comm_rpt)
   interpolated = 1
 else
   MSG_ERROR(sjoin("Invalid value for eph_task:", itoa(dtset%eph_task)))
 end if

 call wrtout([std_out, ab_out], sjoin(ch10, "- Results stored in: ", out_ncpath))
 call wrtout([std_out, ab_out], " Use `abiopen.py out_V1QAVG.nc -e` to visualize results")

#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, out_ncpath, xmpi_comm_self))
   NCF_CHECK(dvdb%cryst%ncwrite(ncid))
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nspden", dvdb%nspden), nctkdim_t("natom", dvdb%natom3 / 3), nctkdim_t("nqpt", this_nqpt), &
     nctkdim_t("natom3", dvdb%natom3), nctkdim_t("ngsmall", ngsmall)], defmode=.True.)
   NCF_CHECK(ncerr)

   if (interpolated == 1) then
     ! Define arrays for Max_r |W(R, r)|
     NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nrpt", dvdb%nrtot), nctkdim_t("nfft", nfft)]))
     ncerr = nctk_def_arrays(ncid, [ &
        nctkarr_t("ngqpt", "int", "three"), nctkarr_t("rpt", "dp", "three, nrpt"), nctkarr_t("rmod", "dp", "nrpt"), &
        nctkarr_t("maxw", "dp", "nrpt, natom3") &
     ])
     NCF_CHECK(ncerr)
   end if

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "symdynmat", "symv1scf", "dvdb_add_lr", "interpolated"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "has_dielt", "has_zeff", "has_quadrupoles", "has_efield", "dvdb_add_lr"])
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"]))
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("v1scf_avg", "dp", "two, nspden, three, natom, nqpt"), &
     nctkarr_t("v1lr_avg", "dp", "two, nspden, three, natom, nqpt"), &
     nctkarr_t("v1scfmlr_avg", "dp", "two, nspden, three, natom, nqpt"), &
     nctkarr_t("v1scfmlr_abs_avg", "dp", "two, nspden, three, natom, nqpt"), &
     nctkarr_t("v1scf_abs_avg", "dp", "two, nspden, three, natom, nqpt"), &
     nctkarr_t("v1lr_abs_avg", "dp", "two, nspden, three, natom, nqpt"), &
     nctkarr_t("gsmall", "int", "three, ngsmall"), &
     nctkarr_t("v1scf_gsmall", "dp", "two, ngsmall, nspden, three, natom, nqpt"), &
     nctkarr_t("v1lr_gsmall", "dp", "two, ngsmall, nspden, three, natom, nqpt"), &
     nctkarr_t("qpoints", "dp", "three, nqpt") &
   ])
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpoints"), this_qpts))
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
       "symdynmat", "symv1scf", "dvdb_add_lr", "interpolated"], &
       [dtset%symdynmat, dvdb%symv1, dtset%dvdb_add_lr, interpolated])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "has_dielt", "has_zeff", "has_quadrupoles", "has_efield"], &
     l2int([dvdb%has_dielt, dvdb%has_zeff, dvdb%has_quadrupoles, dvdb%has_efield]))
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"], [dvdb%qdamp]))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gsmall"), gsmall))
 end if
#endif

 do iq=1,this_nqpt

   if (interpolated == 0) then
     call wrtout(std_out, sjoin(" Treating qpt:", ktoa(this_qpts(:,iq))))

     ! Read data from DVDB file, reconstruct all 3*natom perturbations in tmp_v1r.
     call dvdb%readsym_allv1(dvdb%findq(this_qpts(:, iq)), cplex, nfft, ngfft, tmp_v1r, xmpi_comm_self)

     ! Transfer data to file_v1r taking into account my_npert
     do imyp=1,dvdb%my_npert
       ipc = dvdb%my_pinfo(3, imyp)
       if (cplex == 1) then
         file_v1r(1,:,:,imyp) = tmp_v1r(1,:,:,ipc)
         file_v1r(2,:,:,imyp) = zero
       else
         file_v1r(:,:,:,imyp) = tmp_v1r(:,:,:,ipc)
       end if
     end do
     ABI_FREE(tmp_v1r)

   else
     ! Interpolate my_npert potentials for this q-point.
     call wrtout(std_out, sjoin(" Interpolating qpt:", ktoa(this_qpts(:,iq))))
     call dvdb%ftinterp_qpt(this_qpts(:, iq), nfft, ngfft, file_v1r, comm_rpt) !, add_lr=?)
     cplex = 2
   end if

   ! Compute the periodic part of the LR term (note add_qphase = 0 because we want the periodic part)
   do imyp=1,dvdb%my_npert
     idir = dvdb%my_pinfo(1, imyp); ipert = dvdb%my_pinfo(2, imyp); ipc = dvdb%my_pinfo(3, imyp)
     do ispden=1,min(dvdb%nspden, 2)
       call dvdb%get_v1r_long_range(this_qpts(:,iq), idir, ipert, nfft, ngfft, long_v1r(:,:,ispden,imyp), add_qphase=0)
     end do
   end do

   ! Compute average and write to file.
   if (my_rank /= master) cycle

   do imyp=1,dvdb%my_npert
     idir = dvdb%my_pinfo(1, imyp); ipert = dvdb%my_pinfo(2, imyp); ipc = dvdb%my_pinfo(3, imyp)
     do ispden=1,dvdb%nspden

#ifdef HAVE_NETCDF
       vals2 = sum(file_v1r(:,:,ispden,imyp), dim=2) / nfft
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1scf_avg"), vals2, &
                            start=[1,ispden,idir,ipert,iq], count=[2,1,1,1,1])
       NCF_CHECK(ncerr)

       vals2 = sum(abs(file_v1r(:,:,ispden,imyp)), dim=2) / nfft
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1scf_abs_avg"), vals2, &
                            start=[1,ispden,idir,ipert,iq], count=[2,1,1,1,1])
       NCF_CHECK(ncerr)

       vals2 = sum(long_v1r(:,:,ispden,imyp), dim=2) / nfft
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1lr_avg"), vals2, &
                            start=[1,ispden,idir,ipert,iq], count=[2,1,1,1,1])
       NCF_CHECK(ncerr)
       vals2 = sum(abs(long_v1r(:,:,ispden,imyp)), dim=2) / nfft
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1lr_abs_avg"), vals2, &
                            start=[1,ispden,idir,ipert,iq], count=[2,1,1,1,1])
       NCF_CHECK(ncerr)
       vals2 = sum(file_v1r(:,:,ispden,imyp) - long_v1r(:,:,ispden,imyp), dim=2) / nfft
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1scfmlr_avg"), vals2, &
                            start=[1,ispden,idir,ipert,iq], count=[2,1,1,1,1])
       NCF_CHECK(ncerr)
       vals2 = sum(abs(file_v1r(:,:,ispden,imyp) - long_v1r(:,:,ispden,imyp)), dim=2) / nfft
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1scfmlr_abs_avg"), vals2, &
                            start=[1,ispden,idir,ipert,iq], count=[2,1,1,1,1])
       NCF_CHECK(ncerr)

       ! Compute G-components of DFPT potentials and LR model for G in small-sphere and save results to disk
       call fourdp(2, workg, file_v1r(:,:,ispden,imyp), -1, dvdb%mpi_enreg, nfft, 1, ngfft, 0)
       do ig=1,ngsmall
         work_gsmall(:, ig) = workg(:, ig2ifft(ig))
       end do
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1scf_gsmall"), work_gsmall, &
                            start=[1,1,ispden,idir,ipert,iq], count=[2,ngsmall,1,1,1,1])
       NCF_CHECK(ncerr)

       call fourdp(2, workg, long_v1r(:,:,ispden,imyp), -1, dvdb%mpi_enreg, nfft, 1, ngfft, 0)
       do ig=1,ngsmall
         work_gsmall(:, ig) = workg(:, ig2ifft(ig))
       end do
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1lr_gsmall"), work_gsmall, &
                            start=[1,1,ispden,idir,ipert,iq], count=[2,ngsmall,1,1,1,1])
       NCF_CHECK(ncerr)
#endif

       ! Debugging section.
       !write(std_out, "(a)")"--- !DVDB_LONGRANGE_DIFF"
       !write(std_out,"(3a)")"  qpoint: ", trim(ktoa(this_qpts(:,iq))), ","
       !write(std_out,"(a,i0,a)")"  iq: ", iq, ","
       !write(std_out,"(2(a,i0))")"  idir: ", idir, ", ipert:", ipert
       !write(std_out,"(a,i0,a)")"  ispden: ", ispden, ","
       !call vdiff_print(vdiff_eval(2, nfft, file_v1r(:,:,ispden,imyp), long_v1r(:,:,ispden,imyp), &
       !                 dvdb%cryst%ucvol, vd_max=vd_max))
       !write(std_out,"(a)")"..."

       ! Debug: write potentials to file.
       if (unt /= -1) then
         write(unt,*)"# q-point:", trim(ktoa(this_qpts(:,iq))), ", iq: ", trim(itoa(iq))
         write(unt,*)"# idir: ",idir,", ipert: ",ipert,", ispden:", ispden
         write(unt,*)"# file_v1r, long_v1r, diff"

         if (cplex == 1) then
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 ifft = i1+n1*((i2-1)+n2*(i3-1))
                 write(unt,"(3(i0,1x),3(es12.4,2x))") &
                   i1,i2,i3, &
                   file_v1r(1,ifft,ispden,imyp), long_v1r(1,ifft,ispden,imyp), &
                   file_v1r(1,ifft,ispden,imyp) - long_v1r(1,ifft,ispden,imyp)
               end do
             end do
           end do
         else
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 ifft = i1+n1*((i2-1)+n2*(i3-1))
                 write(unt, "(3(i0,1x),6(es12.4,2x))") &
                   i1,i2,i3, &
                   file_v1r(1,ifft,ispden,imyp), long_v1r(1,ifft,ispden,imyp),  &
                   file_v1r(1,ifft,ispden,imyp) - long_v1r(1,ifft,ispden,imyp), &
                   file_v1r(2,ifft,ispden,imyp), long_v1r(2,ifft,ispden,imyp),  &
                   file_v1r(2,ifft,ispden,imyp) - long_v1r(2,ifft,ispden,imyp)
               end do
             end do
           end do
         end if
         write(unt,*)
         write(unt,*)
       end if

     end do
   end do
   !write(std_out,*)" "
 end do ! iq

 ABI_FREE(long_v1r)
 ABI_FREE(file_v1r)
 ABI_FREE(workg)
 ABI_FREE(gfft)
 ABI_FREE(ig2ifft)
 ABI_FREE(gsmall)
 ABI_FREE(work_gsmall)

 if (interpolated == 1) then
   ! Compute max_r |W(R,r)| and write data to file.
   call dvdb%get_maxw(dtset%dvdb_ngqpt, all_rpt, all_rmod, maxw)
   if (my_rank == master) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt"), dtset%dvdb_ngqpt))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "rpt"), all_rpt))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "rmod"), all_rmod))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "maxw"), maxw))
#endif
   end if
   ABI_FREE(all_rpt)
   ABI_FREE(all_rmod)
   ABI_FREE(maxw)
 end if

 if (my_rank == master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

end subroutine dvdb_write_v1qavg
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_test_ftinterp
!! NAME
!!  dvdb_test_ftinterp
!!
!! FUNCTION
!!  Debugging tool used to test the Fourier interpolation of the DFPT potentials.
!!
!! INPUTS
!!  dvdb_path=Filename
!!  dvdb_ngqpt(3)=Divisions of the Q-mesh reported in the DVDB file
!!  dvdb_add_lr=0 to disable treatment of long-range part in Fourier interpolation.
!!  qdamp=Defines exponential damping in LR potential
!!  ddb_path=Path to DDB file. Used to treat LR part.
!!  prtvol=Verbosity level.
!!  coarse_ngqpt(3)= Coarse q-mesh used to analyze the accuracy of the FT interpolation
!!    Must be divisor of dvdb_ngqpt. Use 0 to disable the test.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_test_ftinterp(dvdb_path, method, symv1, dvdb_ngqpt, dvdb_add_lr, qdamp, &
                              ddb_path, prtvol, coarse_ngqpt, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: dvdb_path, ddb_path
 integer,intent(in) :: comm, prtvol, dvdb_add_lr, qdamp, method, symv1
 integer,intent(in) :: dvdb_ngqpt(3), coarse_ngqpt(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: brav1 = 1, master = 0, chneut2 = 2
 integer :: nfft, iq, cplex, mu, ispden, comm_rpt, iblock_dielt, iblock_dielt_zeff, my_rank,  ierr
 logical :: autotest
 type(dvdb_t) :: dvdb, coarse_dvdb
 type(vdiff_t) :: vd_max
 type(ddb_type) :: ddb
 character(len=fnlen) :: coarse_fname
!arrays
 integer :: ngfft(18), qrefine(3)
 real(dp),allocatable :: file_v1r(:,:,:,:),intp_v1r(:,:,:,:),tmp_v1r(:,:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)
 qrefine = 1

 write(std_out,"(a)")sjoin(" Testing Fourier interpolation of V1(r) with ngqpt:", ltoa(dvdb_ngqpt))
 if (len_trim(ddb_path) > 0) then
   write(std_out,"(a)")sjoin(" Reading Zeff and eps_inf from DDB file:", ddb_path)
   write(std_out,"(a)")sjoin(" dvdb_add_lr set to:", itoa(dvdb_add_lr))
 end if

 dvdb = dvdb_new(dvdb_path, comm)
 dvdb%debug = .False.
 ABI_CHECK(any(symv1 == [0, 1, 2]), sjoin("invalid value of symv1:", itoa(symv1)))
 dvdb%symv1 = symv1
 dvdb%add_lr = dvdb_add_lr
 dvdb%qdamp = qdamp

 !call dvdb%set_pert_distrib(sigma%comm_pert, sigma%my_pinfo, sigma%pert_table)

 iblock_dielt = 0; iblock_dielt_zeff = 0
 if (len_trim(ddb_path) > 0) then
   call dvdb%load_ddb(prtvol, chneut2, comm, ddb_path=ddb_path)
 else
   dvdb%add_lr = 0
   MSG_WARNING("ddb_path was not provided --> Setting dvdb_add_lr to zero")
 end if

 call dvdb%print()

 ! Define FFT mesh for real space representation.
 call ngfft_seq(ngfft, dvdb%ngfft3_v1(:,1))
 nfft = product(ngfft(1:3))
 call dvdb%open_read(ngfft, comm)

 ABI_MALLOC(intp_v1r, (2, nfft, dvdb%nspden, dvdb%natom3))
 ABI_MALLOC(file_v1r, (2, nfft, dvdb%nspden, dvdb%natom3))

 ! Prepare FT interpolation.
 comm_rpt = xmpi_comm_self

 autotest = .True.
 if (autotest) then
   call dvdb%ftinterp_setup(dvdb_ngqpt, qrefine, 1, [zero, zero, zero], nfft, ngfft, method, comm_rpt)

   ! First step: Use FT interpolation to get q-points in the initial ab-initio mesh.
   ! We should get the same result...
   do iq=1,dvdb%nqpt
     ! Read data from DVDB file and store it in file_v1r
     call dvdb%readsym_allv1(dvdb%findq(dvdb%qpts(:,iq)), cplex, nfft, ngfft, tmp_v1r, comm)

     if (cplex == 1) then
       file_v1r(1,:,:,:) = tmp_v1r(1,:,:,:)
       file_v1r(2,:,:,:) = zero
     else
       file_v1r = tmp_v1r
     end if
     ABI_FREE(tmp_v1r)

     ! Interpolate data at the same q-point.
     call dvdb%ftinterp_qpt(dvdb%qpts(:,iq), nfft, ngfft, intp_v1r, dvdb%comm_rpt)

     write(std_out,"(a)")sjoin("=== For q-point:", ktoa(dvdb%qpts(:,iq)), "===")
     do mu=1,dvdb%natom3
       do ispden=1,dvdb%nspden
         write(std_out, "(a)")"--- !DVDB_SELF_DIFF"
         write(std_out,"(3a)")"  qpoint: ", trim(ktoa(dvdb%qpts(:,iq))), ","
         write(std_out,"(a,i0,a)")"  iqpt: ", iq, ","
         write(std_out,"(a,i0,a)")"  iatom3: ", mu, ","
         write(std_out,"(a,i0,a)")"  ispden: ", ispden, ","
         call vdiff_print(vdiff_eval(2, nfft, file_v1r(:,:,ispden,mu), intp_v1r(:,:,ispden,mu), &
                                    dvdb%cryst%ucvol, vd_max=vd_max))
         write(std_out,"(a)")"..."
         !do ifft=1,nfft
         !  write(std_out,*)file_v1r(1,ifft,ispden,mu),intp_v1r(1,ifft,ispden,mu),&
         !  file_v1r(2,ifft,ispden,mu),intp_v1r(2,ifft,ispden,mu)
         !end do
       end do
     end do
     write(std_out,*)" "
   end do ! iq

   write(std_out, "(/, a)")" Max values over q-points and perturbations"
   call vdiff_print(vd_max)
   ABI_FREE(dvdb%v1scf_rpt)
 end if

 ! Now downsample the q-mesh, build real-space representation with coarse q-mesh and
 ! compare with ab-intio values in the initial dvdb.
 if (all(coarse_ngqpt /= 0)) then
   write(std_out, "(/, 2a)")" Downsampling Q-mesh using coarse_ngqpt:", trim(ltoa(coarse_ngqpt))

!Flang compiler complains with empty constructors (this bug should be corrected in future versions)
#if defined FC_LLVM || defined FC_ARM
   vd_max = vdiff_t(zero,zero,zero,zero,zero,zero)
#else
   vd_max = vdiff_t()
#endif

   coarse_fname = strcat(dvdb_path, "_COARSE")
   call dvdb%qdownsample(coarse_fname, coarse_ngqpt, comm)

   coarse_dvdb = dvdb_new(coarse_fname, comm)
   call coarse_dvdb%open_read(ngfft, comm)
   !call coarse_dvdb%set_pert_distrib(sigma%comm_pert, sigma%my_pinfo, sigma%pert_table)

   coarse_dvdb%debug = dvdb%debug
   coarse_dvdb%symv1 = dvdb%symv1
   coarse_dvdb%add_lr = dvdb%add_lr
   coarse_dvdb%has_dielt = dvdb%has_dielt
   coarse_dvdb%has_zeff = dvdb%has_zeff
   coarse_dvdb%has_quadrupoles = dvdb%has_quadrupoles
   coarse_dvdb%has_efield = dvdb%has_efield
   coarse_dvdb%dielt = dvdb%dielt
   coarse_dvdb%zeff = dvdb%zeff
   coarse_dvdb%zeff_raw = dvdb%zeff_raw
   coarse_dvdb%qstar = dvdb%qstar
   coarse_dvdb%qdamp = dvdb%qdamp
   !call coarse_dvdb%print()

   ! Prepare FT interpolation using coarse q-mesh.
   call coarse_dvdb%ftinterp_setup(coarse_ngqpt, qrefine, 1, [zero, zero, zero], nfft, ngfft, method, comm_rpt)

   do iq=1,dvdb%nqpt
     ! Read data from DVDB file and store it in file_v1r
     call dvdb%readsym_allv1(dvdb%findq(dvdb%qpts(:,iq)), cplex, nfft, ngfft, tmp_v1r, comm)

     if (cplex == 1) then
       file_v1r(1,:,:,:) = tmp_v1r(1,:,:,:)
       file_v1r(2,:,:,:) = zero
     else
       file_v1r = tmp_v1r
     end if
     ABI_FREE(tmp_v1r)

     ! Interpolate data at the same q-point using coarse Q-mesh
     call coarse_dvdb%ftinterp_qpt(dvdb%qpts(:,iq), nfft, ngfft, intp_v1r, dvdb%comm_rpt)

     write(std_out,"(a)")sjoin("=== For COARSE q-point:", ktoa(dvdb%qpts(:,iq)), "===")
     do mu=1,dvdb%natom3
       do ispden=1,dvdb%nspden
         write(std_out, "(a)")"--- !DVDB_COARSE_DIFF"
         write(std_out,"(3a)")"  qpoint: ", trim(ktoa(dvdb%qpts(:,iq))), ","
         write(std_out,"(a,i0,a)")"  iqpt: ", iq, ","
         write(std_out,"(a,i0,a)")"  iatom3: ", mu, ","
         write(std_out,"(a,i0,a)")"  ispden: ", ispden, ","
         call vdiff_print(vdiff_eval(2, nfft, file_v1r(:,:,ispden,mu), intp_v1r(:,:,ispden,mu), &
                          dvdb%cryst%ucvol, vd_max=vd_max))
         write(std_out,"(a)")"..."
         !do ifft=1,nfft
         !  write(std_out,*)file_v1r(1,ifft,ispden,mu),intp_v1r(1,ifft,ispden,mu),&
         !  file_v1r(2,ifft,ispden,mu),intp_v1r(2,ifft,ispden,mu)
         !end do
       end do
     end do
     write(std_out,*)" "
   end do ! iq

   write(std_out, "(/, a)")" COARSE DVDB: Max values over q-points and perturbations"
   call vdiff_print(vd_max)
   call coarse_dvdb%free()
   if (my_rank == master) call delete_file(coarse_fname, ierr)
 end if

 ABI_FREE(intp_v1r)
 ABI_FREE(file_v1r)

 call dvdb%free()
 call ddb%free()

end subroutine dvdb_test_ftinterp
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_v1r_long_range
!! NAME
!!  dvdb_get_v1r_long_range
!!
!! FUNCTION
!!  Compute the long-range part of the phonon potential
!!  due to the Born effective charges, PRL 115, 176401 (2015) [[cite:Verdi2015]].
!!
!!    V^L_{iatom,idir}(r) = i (4pi/vol) sum_G (q+G) . Zeff_{iatom,idir}
!!                           e^{i (q+G).(r - tau_{iatom})} / ((q+G) . dielt . (q+G))
!!
!!  where Zeff and dielt are the Born effective charge tensor and the dielectric tensor in cart coords,
!!  tau is the atom position, and vol is the volume of the unit cell.
!!  Note that internally the tensors are stored in Cartesian coordinates while in output we need
!!  the contribution due to the displacement of the iatom-sublattice along the reduced direction idir
!!  hence we need to perform some tensor gymnastics to go from Cart to reduced.
!!
!! INPUTS
!!  db = the DVDB object.
!!  qpt = the q-point in reduced coordinates.
!!  idir = direction index.
!!  iatom = atom index.
!!  nfft = number of fft points.
!!  ngfft(18) = FFT mesh.
!!  [add_qphase]=By default, the routine returns the LR potential with the e^{iqr} phase.
!!    Use add_qphase = 0 to get the lattice-periodic part.
!!
!! OUTPUT
!!  v1r_lr = dipole potential
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_get_v1r_long_range(db, qpt, idir, iatom, nfft, ngfft, v1r_lr, add_qphase)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(in) :: db
 integer,intent(in) :: idir, iatom, nfft
 integer,optional,intent(in) :: add_qphase
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: v1r_lr(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: n1, n2, n3, nfftot, ig, iphase, ii, jj, kk, ll, mm, ifft, ispden
 real(dp) :: fac, qGZ, qGS, denom, denom_inv, qtau, re, im, phre, phim, qg_mod, gsq_max
 real(dp),parameter :: tol_denom = tol8
!arrays
 integer, allocatable :: gfft(:,:)
 real(dp) :: gprimd(3,3), rprimd(3,3), dielt_red(3,3)
 real(dp) :: qG_red(3), qG_cart(3), Zstar(3), Sstar(3,3), tau_red(3)
 real(dp), allocatable :: v1G_lr(:,:), v1G_lr33(:,:,:,:), workr(:,:)

! *************************************************************************

 iphase = 1; if (present(add_qphase)) iphase = add_qphase

 ! Make sure FFT parallelism is not used
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")

 ! Allocate memory
 ABI_MALLOC(gfft, (3, nfft))
 ABI_MALLOC(v1G_lr, (2, nfft))

 ! Reciprocal and real space primitive vectors
 gprimd = db%cryst%gprimd; rprimd = db%cryst%rprimd

 ! Prefactor
 fac = four_pi / db%cryst%ucvol

 ! Transform the Born effective charge tensor from Cartesian to reduced coordinates
 ! and select the relevant direction.
 Zstar = matmul(transpose(gprimd), matmul(db%zeff(:,:,iatom), rprimd(:,idir))) * two_pi

 if (db%has_quadrupoles) then
   ! Transform Qstar from Cartesian to reduced coordinates and select the relevant direction.
   Sstar = zero
   do ii=1,3
     do jj=1,3
       do kk=1,3
         do ll=1,3
           do mm=1,3
             Sstar(ii,jj) = Sstar(ii,jj) + &
                   gprimd(mm,jj) * gprimd(ll,ii) * db%qstar(mm,ll,kk,iatom) * rprimd(kk,idir) * two_pi ** 2
           end do
         end do
       end do
     end do
   end do
 end if

 ! Transform the dielectric tensor from Cartesian to reduced coordinates.
 ! q_cart e_cart q_cart = q_red (G^t e_cart G) q_red
 dielt_red = matmul(transpose(gprimd), matmul(db%dielt, gprimd)) * two_pi ** 2

 ! Atom position
 tau_red = db%cryst%xred(:,iatom)

 ! Get the set of G vectors
 call get_gftt(ngfft, qpt, db%cryst%gmet, gsq_max, gfft)

 ! Compute the long-range potential in G-space due to Z* and Q* (if present)
 v1G_lr = zero
 do ig=1,nfft
   ! (q + G)
   qG_red = qpt + gfft(:,ig)
   qG_cart = two_pi * matmul(db%cryst%gprimd, qG_red)
   qG_mod = sqrt(sum(qG_cart ** 2))
   ! (q + G) . Zeff(:,idir,iatom)
   qGZ = dot_product(qG_red, Zstar)
   ! (q + G) . dielt . (q + G)
   denom = dot_product(qG_red, matmul(dielt_red, qG_red))
   ! Avoid (q + G) = 0
   if (denom < tol_denom) cycle
   denom_inv = one / denom
   ! HM hard cutoff, in this case qdamp takes the meaning of an energy cutoff in Hartree (hardcoded to 1 for the moment)
   !if (half*qG_mod**2 > 1) cycle
   if (db%qdamp > zero) denom_inv = denom_inv * exp(-qG_mod ** 2 / (four * db%qdamp))
   qGS = zero
   if (db%has_quadrupoles) then
     do ii=1,3
       do jj=1,3
         qGS = qGS + qG_red(ii) * qG_red(jj) * Sstar(ii,jj) / two
       end do
     end do
   end if

   ! Phase factor exp(-i (q+G) . tau)
   qtau = - two_pi * dot_product(qG_red, tau_red)
   phre = cos(qtau); phim = sin(qtau)
   !phre = one; phim = zero

   re = +fac * qGS * denom_inv !re = zero
   im = fac * qGZ * denom_inv
   v1G_lr(1,ig) = phre * re - phim * im
   v1G_lr(2,ig) = phim * re + phre * im
 end do

 ! FFT to get the long-range potential in r-space
 call fourdp(2, v1G_lr, v1r_lr, 1, db%mpi_enreg, nfft, 1, ngfft, 0)

 if (db%has_efield) then
   ! Add term due to Electric field.
   ! TODO: Change API to account for ispden/isppol. return nspden LR part
   ! although only the electric field part depends on nsppden.
   !v1r_lr = zero ! Comment this line to have only efield contribution
   ABI_CHECK(db%nspden == 1, "nspden != 1 not coded")
   ispden = 1
   ABI_CALLOC(v1G_lr33, (3, 3, 2, nfft))
   do ig=1,nfft
     !if (ig > 1) cycle
     ! (q + G)
     qG_red = qpt + gfft(:,ig)
     qG_cart = two_pi * matmul(db%cryst%gprimd, qG_red)
     qG_mod = sqrt(sum(qG_cart ** 2))
     ! (q + G) . Zeff(:,idir,iatom)
     !qGZ = dot_product(qG_red, Zstar)
     ! (q + G) . dielt . (q + G)
     denom = dot_product(qG_red, matmul(dielt_red, qG_red))
     ! Avoid (q + G) = 0
     if (denom < tol_denom) cycle
     denom_inv = one / denom
     if (db%qdamp > zero) denom_inv = denom_inv * exp(-qG_mod ** 2 / (four * db%qdamp))
     fac = (four_pi / db%cryst%ucvol) * denom_inv !* qGZ
     ! Phase factor exp(-i (q+G) . tau)
     qtau = - two_pi * dot_product(qG_red, tau_red)
     phre = cos(qtau); phim = sin(qtau)
     !phre = one; phim = zero

     do ii=1,3
       do jj=1,3
         v1G_lr33(ii, jj, 1, ig) = fac * phre * qG_red(ii) * qG_red(jj)
         v1G_lr33(ii, jj, 2, ig) = fac * phim * qG_red(ii) * qG_red(jj)
       end do
     end do
   end do ! ig

   ABI_MALLOC(workr, (2, nfft))
   do ii=1,3
     do jj=1,3
       v1G_lr = v1G_lr33(ii, jj, :, :)
       call fourdp(2, v1G_lr, workr, 1, db%mpi_enreg, nfft, 1, ngfft, 0)
       ! Two pi comes for qpt but we should check whether the gradient wrt E-field is in gprimd or 2pi gprimd coordinates.
       ! MG: Remove two_pi factor because jump discontinuity in the real part for G != 0 are overestimated.
       do ifft=1,nfft
         !v1r_lr(:, ifft) = v1r_lr(:, ifft) - Zstar(ii) * db%v1r_efield(ifft, jj, ispden) * workr(:, ifft) * two_pi
         v1r_lr(:, ifft) = v1r_lr(:, ifft) - Zstar(ii) * db%v1r_efield(ifft, jj, ispden) * workr(:, ifft) ! * two_pi
       end do
     end do
   end do

   ABI_FREE(workr)
   ABI_FREE(v1G_lr33)
 end if

 ! Multiply by exp(i q.r)
 if (iphase == 1) call times_eikr(qpt, ngfft, nfft, 1, v1r_lr)

 ABI_FREE(gfft)
 ABI_FREE(v1G_lr)

end subroutine dvdb_get_v1r_long_range
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_load_ddb
!! NAME
!!  dvdb_load_ddb
!!
!! FUNCTION
!!  Load information about the Born effective charges and dielectric tensor from a DDB file
!!
!! TODO
!!  Use this function in eph driver
!!

subroutine dvdb_load_ddb(dvdb, chneut, prtvol, comm, ddb_path, ddb)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(inout) :: dvdb
 integer,intent(in) :: chneut, prtvol, comm
 character(len=*),optional,intent(in) :: ddb_path
 type(ddb_type),optional,target,intent(in) :: ddb

!Local variables ------------------------------
 integer,parameter :: brav1 = 1, master = 0, natifc0 = 0, rfmeth1 = 1, selectz0 = 0
 integer :: my_rank, iblock_dielt, iblock_dielt_zeff
 logical :: free_ddb
 type(crystal_t) :: cryst_ddb
 type(ddb_type),pointer :: ddb_ptr
 type(ddb_type),target :: this_ddb
!arrays
 integer,allocatable :: dummy_atifc(:)
 real(dp) :: dielt(3,3)
 real(dp),allocatable :: zeff(:,:,:), zeff_raw(:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 if (present(ddb_path)) then
   ! Build ddb object from file. Will release memory before returning.
   ABI_CHECK(.not. present(ddb), "ddb argument cannot be present when ddb_path is used")
   ABI_CALLOC(dummy_atifc, (dvdb%natom))
   call ddb_from_file(this_ddb, ddb_path, brav1, dvdb%natom, natifc0, dummy_atifc, cryst_ddb, comm, prtvol=prtvol)
   ABI_FREE(dummy_atifc)
   call cryst_ddb%free()
   ddb_ptr => this_ddb
   free_ddb = .True.
 else
   ! Point input ddb, won't release memory.
   free_ddb = .False.
   ddb_ptr => ddb
 end if

 ! Get dielectric Tensor
 iblock_dielt = ddb_ptr%get_dielt(rfmeth1, dielt)
 dvdb%dielt = dielt

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 ABI_MALLOC(zeff, (3, 3, dvdb%natom))
 ABI_MALLOC(zeff_raw, (3, 3, dvdb%natom))
 iblock_dielt_zeff = ddb_ptr%get_dielt_zeff(dvdb%cryst, rfmeth1, chneut, selectz0, dielt, zeff, zeff_raw=zeff_raw)

 if (my_rank == master) then
   if (iblock_dielt_zeff == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_path))
     call wrtout(ab_out, "Values initialized with zeros")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_path))
   end if
 end if

 if (iblock_dielt /= 0) then
   dvdb%has_dielt = .True.
   dvdb%dielt = dielt
 end if
 if (iblock_dielt_zeff /= 0) then
   dvdb%has_zeff = .True.; dvdb%zeff = zeff; dvdb%zeff_raw = zeff_raw
 end if
 if (dvdb%has_dielt .and. (dvdb%has_zeff .or. dvdb%has_quadrupoles)) then
   if (dvdb%add_lr == 0)  then
     call wrtout([std_out, ab_out], &
       " WARNING: dvdb_add_lr set to 0. Long-range term won't be substracted in Fourier interpolation.")
   end if
 end if

 ABI_FREE(zeff)
 ABI_FREE(zeff_raw)
 if (free_ddb) call ddb_ptr%free()

end subroutine dvdb_load_ddb
!!***

!!****f* m_dvdb/dvdb_load_efield
!! NAME
!!  dvdb_load_efield
!!
!! FUNCTION
!!  Load first oder derivatives wrt the electric file from files
!!
!! INPUTS
!!  pot_paths=List of strings with paths to POT1 files.
!!  comm=MPI communicator.
!!

subroutine dvdb_load_efield(dvdb, pot_paths, comm)

!Arguments ------------------------------------
!scalars
 class(dvdb_t),intent(inout) :: dvdb
 integer,intent(in) :: comm
 character(len=*),intent(in) :: pot_paths(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: pawread0 = 0, cplex1 = 1
 integer :: ii, idir, ipert, nfft
 type(hdr_type) :: hdr
!arrays
 real(dp),allocatable :: v1e_red(:,:,:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)

! *************************************************************************

 ABI_CHECK(all(dvdb%ngfft /= -1), "dbvd%ngfft must be defined!")

 nfft = product(dvdb%ngfft(1:3))
 ABI_CALLOC(v1e_red, (nfft, dvdb%nspden, 3))

 do ii=1,3
   ! Read DFPT potentials due to E-field.
   ! TODO: Should implement symmetries so that only the irred pots are needed.
   call wrtout(std_out, sjoin("Loading Efield DFPT potential from:", pot_paths(ii)))
   call read_rhor(pot_paths(ii), cplex1, dvdb%nspden, nfft, dvdb%ngfft, pawread0, &
     dvdb%mpi_enreg, v1e_red(:,:,ii), hdr, pawrhoij, comm, allow_interp=.True.)

   ! Consistency check: expecting E-field perturbation.
   idir = mod(hdr%pertcase - 1, 3) + 1; ipert = (hdr%pertcase - idir) / 3 + 1
   ABI_CHECK(all(abs(hdr%qptn) < tol12), sjoin("Expecting Gamma point in E-field pert, got qpt:", ktoa(hdr%qptn)))
   ABI_CHECK(ipert == hdr%natom + 2, sjoin("Expecting E-field perturbation, got ipert:", itoa(ipert)))
   ABI_CHECK(idir == ii, sjoin("Expecting E-field perturbation along idir:", itoa(ii), " got idir:", itoa(idir)))
   call hdr%free()
 end do

 ! Transfer data.
 ABI_MALLOC(dvdb%v1r_efield, (nfft, 3, dvdb%nspden))
 do ii=1,3
   dvdb%v1r_efield(:,ii,:) = v1e_red(:,:,ii)
 end do
 ABI_FREE(v1e_red)

 dvdb%has_efield = .True.

end subroutine dvdb_load_efield
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_interpolate_and_write
!! NAME
!!  dvdb_interpolate_and_write
!!
!! FUNCTION
!!  Interpolate the phonon potential onto a fine q-point grid
!!  and write the data in a new DVDB file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_interpolate_and_write(dvdb, dtset, new_dvdb_fname, ngfft, ngfftf, cryst, &
           ngqpt_coarse, nqshift_coarse, qshift_coarse, comm, custom_qpt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift_coarse, comm
 character(len=*),intent(in) :: new_dvdb_fname
 type(crystal_t),intent(in) :: cryst
 class(dvdb_t),intent(inout) :: dvdb
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18), ngfftf(18)
 integer,intent(in) :: ngqpt_coarse(3)
 real(dp),intent(in) :: qshift_coarse(3,nqshift_coarse)
 real(dp),optional,intent(in) :: custom_qpt(:,:)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: fform_pot=111
 integer :: my_rank,nproc,idir,ipert,iat,ipc,ispden, ierr
 integer :: cplex,db_iqpt,natom,natom3,npc,trev_q,nspden
 integer :: nqbz, nqibz, iq, ifft, nqbz_coarse
 integer :: nperts_read, nperts_interpolate, nperts
 integer :: nqpt_read, nqpt_interpolate
 integer :: nfft,nfftf, dimv1
 integer :: ount, unt, fform
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
 logical :: use_netcdf
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all
 character(len=500) :: msg
 character(len=fnlen) :: tmp_fname
 type(hdr_type) :: hdr_ref
!arrays
 integer :: qptrlatt(3,3), rfdir(3)
 integer :: symq(4,2,cryst%nsym)
 integer,allocatable :: pinfo(:,:),rfpert(:),pertsy(:,:,:),iq_read(:),this_pertsy(:,:)
 real(dp) :: qpt(3), rhog1_g0(2)
 real(dp),allocatable :: v1scf(:,:,:), v1scf_rpt(:,:,:,:),v1(:)
 real(dp),allocatable :: wtq(:),qibz(:,:),qbz(:,:),q_interp(:,:),q_read(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 write(msg, '(2a)') " Interpolation of the electron-phonon coupling potential", ch10
 call wrtout(ab_out, msg, do_flush=.True.); call wrtout(std_out, msg, do_flush=.True.)

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 if (dtset%eph_task == 5 .or. present(custom_qpt)) then
   msg = sjoin(" From coarse q-mesh:", ltoa(ngqpt_coarse), "to:", ltoa(dtset%eph_ngqpt_fine))
   call wrtout([std_out, ab_out], msg)
   ! Setup fine q-point grid in the IBZ
   ! Generate the list of irreducible q-points in the grid
   qptrlatt = 0
   qptrlatt(1,1) = dtset%eph_ngqpt_fine(1); qptrlatt(2,2) = dtset%eph_ngqpt_fine(2); qptrlatt(3,3) = dtset%eph_ngqpt_fine(3)
   call kpts_ibz_from_kptrlatt(cryst, qptrlatt, dtset%qptopt, 1, [zero, zero, zero], nqibz, qibz, wtq, nqbz, qbz)

 else if (dtset%eph_task == -5) then
   msg = sjoin(" Using list of q-points specified by ph_qpath with ", itoa(dtset%ph_nqpath), "qpoints")
   call wrtout([std_out, ab_out], msg)
   ABI_CHECK(dtset%ph_nqpath > 0, "ph_nqpath must be specified when eph_task == -5")
   nqibz = dtset%ph_nqpath
   ABI_MALLOC(qibz, (3, nqibz))
   qibz = dtset%ph_qpath(:, 1:nqibz)
   ABI_CALLOC(wtq, (nqibz))
   nqbz = nqibz
   ABI_MALLOC(qbz, (3, nqbz))
   qbz = qibz

 else
   MSG_ERROR(sjoin("Invalid eph_task", itoa(dtset%eph_task)))
 end if

 if (present(custom_qpt)) then
  ABI_SFREE(qibz)
  ABI_SFREE(wtq)
  ABI_SFREE(qbz)
  nqibz = size(custom_qpt,dim=2)
  ABI_MALLOC(qibz, (3, nqibz))
  qibz = custom_qpt
  ABI_CALLOC(wtq, (nqibz))
  nqbz = nqibz
  ABI_MALLOC(qbz, (3, nqbz))
  qbz = qibz
 end if

 nfft = product(ngfft(1:3)); nfftf = product(ngfftf(1:3))

 ! check that ngqpt_coarse is in DVDB.
 nqbz_coarse = product(ngqpt_coarse) * nqshift_coarse

 ! ==========================================
 ! Prepare the header to write the potentials
 ! ==========================================

 ! Read the first header
 if (my_rank == master) then
   if (open_file(dvdb%path, msg, newunit=unt, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   read(unt, err=10, iomsg=msg) dvdb%version
   read(unt, err=10, iomsg=msg) dvdb%numv1

   call hdr_fort_read(hdr_ref, unt, fform)
   if (dvdb_check_fform(fform, "read_dvdb", msg) /= 0) then
     MSG_ERROR(sjoin("While reading:", dvdb%path, ch10, msg))
   end if
   close(unt)
 end if

 ! Reset the symmetries of the header
 ! One might have disable the symmetries in the response function calculation
 ! that produced the initial set of potentials present in the DVDB.
 ! This is because the symmetry features are not used in all parts
 ! of the response function driver.
 !write(std_out,*)hdr_ref%nsym, cryst%nsym
 !ABI_CHECK(hdr_ref%nsym == cryst%nsym, "Diff nsym")
 ABI_SFREE(hdr_ref%symrel)
 ABI_SFREE(hdr_ref%tnons)
 ABI_SFREE(hdr_ref%symafm)
 hdr_ref%nsym = cryst%nsym
 ABI_MALLOC(hdr_ref%symrel, (3,3,hdr_ref%nsym))
 ABI_MALLOC(hdr_ref%tnons, (3,hdr_ref%nsym))
 ABI_MALLOC(hdr_ref%symafm, (hdr_ref%nsym))

 hdr_ref%symrel(:,:,:) = cryst%symrel(:,:,:)
 hdr_ref%tnons(:,:) = cryst%tnons(:,:)
 hdr_ref%symafm(:) = cryst%symafm(:)
 hdr_ref%ngfft = ngfftf(1:3)

 ! =======================================
 ! Open DVDB and copy important dimensions
 ! =======================================

 call dvdb%open_read(ngfftf, xmpi_comm_self)

 ! Besides perturbations with same q-points won't be contiguous on file --> IO is gonna be inefficient.
 call dvdb%print(prtvol=dtset%prtvol)

 natom = cryst%natom
 natom3 = 3 * natom
 nspden = dvdb%nspden

 ! ==================================================
 ! Sort the q-points to read and those to interpolate
 ! and find the irreducible perturbations
 ! ==================================================
 ABI_MALLOC(iq_read, (nqibz))
 ABI_MALLOC(q_read, (3,nqibz))
 ABI_MALLOC(q_interp, (3,nqibz))
 ABI_MALLOC(pertsy, (nqibz,3,dvdb%mpert))
 ABI_MALLOC(this_pertsy, (3,dvdb%mpert))
 ABI_MALLOC(rfpert, (dvdb%mpert))
 ABI_MALLOC(pinfo, (3,3*dvdb%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1

 pertsy = 0
 nqpt_read = 0
 nperts_read = 0
 nqpt_interpolate = 0
 nperts_interpolate = 0

 do iq=1,nqibz
   qpt = qibz(:,iq)

   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qpt)
   !if (db_iqpt /= 1) db_iqpt = -1

   if (db_iqpt /= -1) then
     if (dvdb%prtvol > 0) call wrtout(std_out, sjoin("Q-point: ",ktoa(qpt)," found in DVDB with index ",itoa(db_iqpt)))
     nqpt_read = nqpt_read + 1
     q_read(:,nqpt_read) = qpt(:)
     iq_read(nqpt_read) = db_iqpt

     ! Count the perturbations
     npc = dvdb_get_pinfo(dvdb, db_iqpt, cplex, pinfo)
     do ipc=1,npc
       idir = pinfo(1,ipc); iat = pinfo(2,ipc); ipert = pinfo(3, ipc)
       if (iat .le. natom) nperts_read = nperts_read + 1
     end do

   else
     if (dvdb%prtvol > 0) call wrtout(std_out, sjoin("Q-point: ",ktoa(qpt), "not found in DVDB. Will interpolate."))
     nqpt_interpolate = nqpt_interpolate + 1
     q_interp(:,nqpt_interpolate) = qpt(:)

     ! Examine the symmetries of the q wavevector
     call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,trev_q,prtvol=0)

     ! Find the list of irreducible perturbations for this q-point.
     call irreducible_set_pert(cryst%indsym,dvdb%mpert,cryst%natom,cryst%nsym,&
         this_pertsy,rfdir,rfpert,symq,cryst%symrec,cryst%symrel)
         pertsy(nqpt_interpolate,:,:) = this_pertsy
     !pertsy = 1

     do iat=1,natom
       do idir=1,3
         ipert = (iat-1) * 3 + idir
         if (pertsy(nqpt_interpolate,idir,iat) == 1) nperts_interpolate = nperts_interpolate + 1
       end do
     end do

   end if
 end do

 call wrtout([std_out, ab_out], sjoin(" Number of q-points found in input DVDB:", itoa(nqpt_read)))
 call wrtout([std_out, ab_out], sjoin(" Number of q-points requiring Fourier interpolation", itoa(nqpt_interpolate)))

 ! =================================================
 ! Open the new DVDB file and write preliminary info
 ! =================================================
 nperts = nperts_read + nperts_interpolate

 if (my_rank == master) then
   if (open_file(new_dvdb_fname, msg, newunit=ount, form="unformatted", action="write", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if
   write(ount, err=10, iomsg=msg) dvdb_last_version
   write(ount, err=10, iomsg=msg) nperts
 end if

 ! =================================================================
 ! Master reads all available perturbations and copy in the new DVDB
 ! =================================================================

 rhog1_g0 = zero

 if (my_rank == master) then
   do iq=1,nqpt_read
     qpt = q_read(:,iq)
     db_iqpt = iq_read(iq)

     ! Read each irreducible perturbation potentials
     npc = dvdb_get_pinfo(dvdb, db_iqpt, cplex, pinfo)
     ABI_CHECK(npc /= 0, "npc == 0!")

     ! These arrays depend on cplex.
     ABI_MALLOC(v1scf, (cplex, nfftf, nspden))
     ABI_MALLOC(v1, (cplex*nfftf))

     do ipc=1,npc
       idir = pinfo(1,ipc); iat = pinfo(2,ipc); ipert = pinfo(3, ipc)
       if (dvdb%read_onev1(idir, iat, db_iqpt, cplex, nfftf, ngfftf, v1scf, msg) /= 0) then
         MSG_ERROR(msg)
       end if

       ! Write header
       hdr_ref%qptn = qpt
       hdr_ref%pertcase = ipert
       call hdr_ref%fort_write(ount, fform_pot, ierr)
       ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

       do ispden=1,nspden
         v1 = reshape(v1scf(:,:,ispden), (/cplex*nfftf/))
         write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfftf)
       end do
       if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0
     end do

     ABI_FREE(v1scf)
     ABI_FREE(v1)
   end do
 end if ! master

 call xmpi_barrier(comm)

 ! ================================================================
 ! Interpolate the potential for q-points not in the original DVDB
 ! ================================================================

 dvdb%my_nrpt = nqbz_coarse
 ABI_MALLOC_OR_DIE(v1scf_rpt, (2, dvdb%my_nrpt, nfftf, dvdb%nspden), ierr)

 cplex = 2
 ABI_MALLOC(v1scf, (cplex,nfftf,nspden))
 ABI_MALLOC(v1, (cplex*nfftf))

 use_netcdf = .False.
#ifdef HAVE_NETCDF
 ! Create temporary netcdf file used to write Fortran file with contiguous perturbations.
 use_netcdf = .True.
 if (my_rank == master) then
   tmp_fname = strcat(new_dvdb_fname, "_TEMPORARY_TRANSFER_FILE.nc")
   dimv1 = cplex * nfftf
   NCF_CHECK(nctk_open_create(ncid, tmp_fname, xmpi_comm_self))
   ncerr = nctk_def_dims(ncid, [&
     nctkdim_t("dimv1", dimv1), nctkdim_t("nspden", nspden), &
     nctkdim_t("natom", natom), nctkdim_t("nqpt_intp", nqpt_interpolate), &
     nctkdim_t("nrpt", dvdb%my_nrpt), nctkdim_t("nfft", nfftf), nctkdim_t("natom3", natom * 3) &
   ])
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("v1", "dp", "dimv1, nspden, three, natom, nqpt_intp")))
   !NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("v1scf_rpt", "dp", "two, nrpt, nfft, nspden, natom3")))
   NCF_CHECK(nctk_set_datamode(ncid))
 end if
#endif

 do iat=1,natom
   do idir=1,3
     ipert = (iat-1) * 3 + idir

     ! Entry set to -1 for perturbations that can be found from basis perturbations.
     if (sum(pertsy(:,idir,iat)) == -nqpt_interpolate) cycle

     call wrtout(std_out, sjoin(" Interpolating perturbation iat, idir = ",itoa(iat), itoa(idir)), do_flush=.True.)
     call cwtime(cpu, wall, gflops, "start")

     ! TODO: This part is slow.
     ! Compute phonon potential in real space lattice representation.
     call dvdb_get_v1scf_rpt(dvdb, cryst, ngqpt_coarse, nqshift_coarse, &
                             qshift_coarse, nfftf, ngfftf, &
                             dvdb%my_nrpt, dvdb%nspden, ipert, v1scf_rpt, comm)

     !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1scf_rpt"), v1scf_rpt, start=[1,1,1,1,ipert]))
     call cwtime_report(" v1scf_rpt built", cpu, wall, gflops)

     do iq=1,nqpt_interpolate
       if (pertsy(iq,idir,iat) == -1) cycle
       qpt = q_interp(:,iq)

       ! Interpolate the phonon potential
       call dvdb_get_v1scf_qpt(dvdb, cryst, qpt, nfftf, ngfftf, dvdb%my_nrpt, &
                               dvdb%nspden, ipert, v1scf_rpt, v1scf, comm)

       !call wrtout(std_out, sjoin("Writing q-point", itoa(iq)))
       if (my_rank == master) then
         if (use_netcdf) then
#ifdef HAVE_NETCDF
           ncerr = nf90_put_var(ncid, nctk_idname(ncid, "v1"), v1scf, &
               start=[1,1,idir,iat,iq], count=[dimv1,nspden,1,1,1])
           NCF_CHECK(ncerr)
#endif
         else
           ! Master writes the file (change also qpt and ipert in hdr%)
           hdr_ref%qptn = qpt
           hdr_ref%pertcase = ipert
           call hdr_ref%fort_write(ount, fform_pot, ierr)
           ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

           do ispden=1,nspden
             v1 = reshape(v1scf(:,:,ispden), [cplex*nfftf])
             write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfftf)
           end do
           if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0
         end if
       end if
     end do

     call cwtime_report(" q-points interpolated and written to new DVDB file.", cpu, wall, gflops)
     ABI_FREE(dvdb%my_rpt)
   end do
 end do

 if (use_netcdf .and. my_rank == master) then
   do iq=1,nqpt_interpolate
     qpt = q_interp(:,iq)
     do iat=1,natom
       do idir=1,3
         if (pertsy(iq,idir,iat) == -1) cycle
         ipert = (iat-1) * 3 + idir
         hdr_ref%qptn = qpt
         hdr_ref%pertcase = ipert
         call hdr_ref%fort_write(ount, fform_pot, ierr)
#ifdef HAVE_NETCDF
         ncerr = nf90_get_var(ncid, nctk_idname(ncid, "v1"), v1scf, &
             start=[1,1,idir,iat,iq], count=[dimv1,nspden,1,1,1])
         NCF_CHECK(ncerr)
#endif
         do ispden=1,nspden
           v1 = reshape(v1scf(:,:,ispden), [cplex*nfftf])
           write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfftf)
         end do
         if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0
       end do
    end do
   end do
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ncid))
#endif
   call delete_file(tmp_fname, ierr)
 end if

 if (my_rank == master) close(ount)

 ! Free memory
 ABI_FREE(v1scf)
 ABI_FREE(v1)
 ABI_FREE(v1scf_rpt)
 ABI_FREE(qbz)
 ABI_FREE(qibz)
 ABI_FREE(q_interp)
 ABI_FREE(q_read)
 ABI_FREE(wtq)
 ABI_FREE(iq_read)
 ABI_FREE(pertsy)
 ABI_FREE(this_pertsy)
 ABI_FREE(rfpert)
 ABI_FREE(pinfo)

 call hdr_ref%free()

 write(msg, '(2a)') "Interpolation of the electron-phonon coupling potential completed", ch10
 call wrtout([std_out, ab_out], msg, do_flush=.True.)

 call cwtime_report(" Overall time:", cpu_all, wall_all, gflops_all)

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(msg)

end subroutine dvdb_interpolate_and_write
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_qdownsample
!! NAME
!!  dvdb_qdownsample
!!
!! FUNCTION
!!  Downsample the q-mesh. Produce new DVDB file
!!
!! INPUTS
!!  new_dvdb_fname=Path of output DVDB
!!  ngqpt(3)=Division of coarse Q-mesh
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_qdownsample(dvdb, new_dvdb_fname, ngqpt, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 class(dvdb_t),intent(inout) :: dvdb
 character(len=*),intent(in) :: new_dvdb_fname
!arrays
 integer,intent(in) :: ngqpt(3)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0, qptopt1=1, timrev1=1
 integer :: fform_pot=111
 integer :: ierr,my_rank,nproc,idir,ipert,iat,ipc,ispden
 integer :: cplex, db_iqpt, npc, nqbz, nqibz, iq, ifft, nperts_read, nfft, ount
 character(len=500) :: msg
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: iq_read(:), pinfo(:,:)
 real(dp) :: rhog1_g0(2)
 real(dp),allocatable :: v1scf(:,:,:), v1(:), wtq(:), qibz(:,:), qbz(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 if (my_rank /= master) goto 20

 nfft = product(dvdb%ngfft(1:3))

 ! =======================
 ! Setup fine q-point grid
 ! =======================
 ! Generate the list of irreducible q-points in the coarse grid
 qptrlatt = 0; qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(dvdb%cryst, qptrlatt, qptopt1, 1, &
                             [zero, zero, zero], nqibz, qibz, wtq, nqbz, qbz)

 ! =======================================
 ! Open DVDB and copy important dimensions
 ! =======================================

 ABI_MALLOC(iq_read, (nqibz))
 ABI_MALLOC(pinfo, (3, 3*dvdb%mpert))
 nperts_read = 0

 do iq=1,nqibz
   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qibz(:, iq))
   ABI_CHECK(db_iqpt /= -1, sjoin("Q-point:", ktoa(qibz(:, iq)), "not found in DVDB!"))
   iq_read(iq) = db_iqpt

   ! Count the number of perturbations.
   npc = dvdb_get_pinfo(dvdb, db_iqpt, cplex, pinfo)
   do ipc=1,npc
     idir = pinfo(1,ipc); iat = pinfo(2,ipc); ipert = pinfo(3, ipc)
     if (iat <= dvdb%cryst%natom) nperts_read = nperts_read + 1
   end do
 end do

 ! =================================================
 ! Open the new DVDB file and write preliminary info
 ! =================================================
 !nperts = nperts_read + nperts_interpolate
 if (open_file(new_dvdb_fname, msg, newunit=ount, form="unformatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 write(ount, err=10, iomsg=msg) dvdb_last_version
 write(ount, err=10, iomsg=msg) nperts_read

 ! Read all perturbations on the coarse Q-mesh and write them to the new DVDB
 rhog1_g0 = zero

 do iq=1,nqibz
   db_iqpt = iq_read(iq)

   ! Read each irreducible perturbation potentials
   npc = dvdb_get_pinfo(dvdb, db_iqpt, cplex, pinfo)
   ABI_CHECK(npc /= 0, "npc == 0!")

   ABI_MALLOC(v1scf, (cplex, nfft, dvdb%nspden))
   ABI_MALLOC(v1, (cplex*nfft))

   do ipc=1,npc
     idir = pinfo(1,ipc); iat = pinfo(2,ipc); ipert = pinfo(3, ipc)
     if (dvdb%read_onev1(idir, iat, db_iqpt, cplex, nfft, dvdb%ngfft, v1scf, msg) /= 0) then
       MSG_ERROR(msg)
     end if

     ! Change the header.
     dvdb%hdr_ref%qptn = qibz(:, iq)
     dvdb%hdr_ref%pertcase = ipert

     ! Write header
     call dvdb%hdr_ref%fort_write(ount, fform_pot, ierr)
     ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

     do ispden=1,dvdb%nspden
       v1 = reshape(v1scf(:,:,ispden), [cplex*nfft])
       write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfft)
     end do
     if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0
   end do

   ABI_FREE(v1scf)
   ABI_FREE(v1)
 end do

 close(ount)

 ! Free memory
 ABI_FREE(qbz)
 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(iq_read)
 ABI_FREE(pinfo)

 write(msg, '(2a)') " Downsampling of the e-ph coupling potential completed", ch10
 call wrtout(std_out, msg, do_flush=.True.)

20 continue
 call xmpi_barrier(comm)

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(msg)

end subroutine dvdb_qdownsample
!!***

end module m_dvdb
!!***
