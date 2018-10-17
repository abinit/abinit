!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dvdb
!! NAME
!!  m_dvdb
!!
!! FUNCTION
!!  Objects and methods to extract data from the DVDB file.
!!  The DVDB file is Fortran binary file with a collection of DFPT potentials
!!  associated to the different phonon perturbations (idir, ipert, qpt).
!!  DVDB files are produced with the `mrgdv` utility and are used in the EPH code
!!  to compute the matrix elements <k+q| dvscf_{idir, ipert, qpt} |k>.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (MG,GA)
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

MODULE m_dvdb

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

 use defs_abitypes,   only : hdr_type, mpi_type, datafiles_type
 use m_fstrings,      only : strcat, sjoin, itoa, ktoa, ltoa, ftoa, yesno, endswith
 use m_io_tools,      only : open_file, file_exists
 use m_numeric_tools, only : wrap2_pmhalf, vdiff_eval, vdiff_print
 use m_symtk,         only : mati3inv, littlegroup_q
 use m_geometry,      only : littlegroup_pert, irreducible_set_pert
 use m_copy,          only : alloc_copy
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_fftcore,       only : ngfft_seq
 use m_fft_mesh,      only : rotate_fft_mesh, times_eigr, times_eikr, ig2gfft, get_gftt, calc_ceikr, calc_eigr
 use m_fft,           only : fourdp
 use m_crystal,       only : crystal_t, crystal_free, crystal_print
 use m_crystal_io,    only : crystal_from_hdr
 use m_kpts,          only : kpts_ibz_from_kptrlatt, listkk
 use m_spacepar,      only : symrhg, setsym
 use m_fourier_interpol,only : fourier_interpol

 implicit none

 private
!!***

 ! Version 1: header + vscf1(r)
 ! Version 2: header + vscf1(r) + record with rhog1(G=0)
 !integer,public,parameter :: dvdb_last_version = 1
 integer,public,parameter :: dvdb_last_version = 2

 integer,private,parameter :: DVDB_NOMODE    = 0
 integer,private,parameter :: DVDB_READMODE  = 1
 integer,private,parameter :: DVDB_WRITEMODE = 2

 ! FIXME
 !real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors

 integer,private,parameter :: FPOS_EOF = -1

!----------------------------------------------------------------------

!!****t* m_dvdb/dvdb_t
!! NAME
!!  dvdb_t
!!
!! FUNCTION
!!  Database of DFPT results. The database contains `numv1` perturbations
!!  and the corresponding first order local potential in real space on the FFT mesh.
!!  Note that one can have different FFT meshes for the different perturbations.
!!
!! NOTES
!!  natom, nspden, nspinor, and usepaw are global variables in the sense that it's not possible to add
!!  new entries to the database if these dimension differ from the global ones.
!!
!! TODO
!!  1) Add option to store potentials in memory to reduce IO pressure
!!  2) handle q_bz = S q_ibz case by symmetrizing the potentials already available in the DVDB.
!!     without performing FT interpolation.
!!  3) Implement treatment of long-range behaviour in FT interpolation in polar semiconductors.
!!
!! SOURCE

 type,public :: dvdb_t

  integer :: fh
   ! file handler
   !  Fortran unit number if iomode==IO_MODE_FORTRAN
   !  MPI file handler if iomode==IO_MODE_MPI

  integer :: comm
  ! MPI communicator used for IO.

  integer :: version
  ! File format version read from file.

  integer :: iomode=IO_MODE_FORTRAN
  ! Method used to access the DVDB file:
  !   IO_MODE_FORTRAN for usual Fortran IO routines
  !   IO_MODE_MPI if MPI/IO routines.

  integer :: rw_mode = DVDB_NOMODE
   ! (Read|Write) mode

  integer :: current_fpos
  ! The current position of the file pointer used for sequential access with Fortran-IO
  !  FPOS_EOF signals the end of file

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

  integer :: nrpt=0
  ! Number of real space points used for Fourier interpolation

  integer :: prtvol=0
   ! Verbosity level

  logical :: debug=.False.
   ! Debug flag

  logical :: has_dielt_zeff=.False.
   ! Does the dvdb have the dielectric tensor and Born effective charges

  logical :: symv1=.False.
   ! Activate symmetrization of v1 potentials.

  character(len=fnlen) :: path = ABI_NOFILE
   ! File name

  real(dp) :: dielt(3,3)
   ! Dielectric tensor

  integer,allocatable :: pos_dpq(:,:,:)
   ! pos_dpq(3, mpert, nqpt)
   ! The position of the (idir, ipert, iqpt) potential in the file (in units for POT1 blocks)
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

  integer :: ngfft(18)=-1
   ! Info on the FFT to be used for the potentials.

  integer,allocatable :: iv_pinfoq(:,:)
   !iv_pinfoq(4, numv1)
   !  iv_pinfoq(1, iv1) gives the `idir` index of the iv1 potential
   !  iv_pinfoq(2, iv1) gives the `ipert` index of the iv1 potential
   !  iv_pinfoq(3, iv1) gives `pertcase`=idir + (ipert-1)*3
   !  iv_pinfoq(4, iv1) gives the `iqpt` index of the iv1 potential

  ! FIXME: (3) or (18)
  integer,allocatable :: ngfft3_v1(:,:)
   ! ngfft3_v1(3, numv1)
   ! The FFT mesh used for each v1 potential (the one used to store data in the file).

  real(dp),allocatable :: qpts(:,:)
   ! qpts(3,nqpt)
   ! List of q-points in reduced coordinates.

  real(dp),allocatable :: rpt(:,:)
  ! rpt(3,nrpt)
  ! Real space points for Fourier interpolation.

  real(dp),allocatable :: v1scf_rpt(:,:,:,:,:)
  ! Wannier representation (real values)
  ! v1scf_rpt(nrpt, nfft , nspden, 3*natom)

  real(dp),allocatable :: rhog1_g0(:,:)
  ! rhog1_g0(2, numv1)
  ! G=0 component of rhog1. Used to treat the long range component
  ! in (polar) semiconductors

  real(dp),allocatable :: zeff(:,:,:)
  ! zeff(3,3,natom)
  ! Effective charge on each atom, versus electric field and atomic displacement.

  type(crystal_t) :: cryst
  ! Crystalline structure read from the the DVDB file.

  type(mpi_type) :: mpi_enreg
  ! Internal object used to call fourdp

 end type dvdb_t

 public :: dvdb_init              ! Initialize the object.
 public :: dvdb_open_read         ! Open the file in read-only mode.
 public :: dvdb_close             ! Close file.
 public :: dvdb_free              ! Release the memory allocated and close the file.
 public :: dvdb_print             ! Release memory.
 public :: dvdb_findq             ! Returns the index of the q-point.
 public :: dvdb_read_onev1        ! Read and return the DFPT potential for given (idir, ipert, iqpt).
 public :: dvdb_readsym_allv1     ! Read and return all the 3*natom DFPT potentials (either from file or symmetrized)
 public :: dvdb_list_perts        ! Check if all the (phonon) perts are available taking into account symmetries.
 public :: dvdb_ftinterp_setup    ! Prepare the internal tables for Fourier interpolation.
 public :: dvdb_ftinterp_qpt      ! Fourier interpolation of potentials for given q-point
 public :: dvdb_merge_files       ! Merge a list of POT1 files.
 public :: dvdb_v1r_long_range    ! Long-range part of the phonon potential
 public :: dvdb_interpolate_v1scf ! Fourier interpolation of the phonon potentials
 public :: dvdb_get_v1scf_rpt     ! Fourier transform of the phonon potential from qpt to R
 public :: dvdb_get_v1scf_qpt     ! Fourier transform of the phonon potential from R to qpt
 public :: dvdb_interpolate_and_write ! Interpolate the phonon potentials and write a new DVDB file.

! Debugging tools.
 public :: dvdb_test_v1rsym       ! Check symmetries of the DFPT potentials.
 public :: dvdb_test_v1complete   ! Debugging tool used to test the symmetrization of the DFPT potentials.
 public :: dvdb_test_ftinterp     ! Test Fourier interpolation of DFPT potentials.
!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_init
!! NAME
!!  dvdb_init
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_init(db, path, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm
 type(dvdb_t),intent(inout) :: db

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,timrev2=2
 integer :: iv1,ii,ierr,unt,fform,nqpt,iq,iq_found,cplex,trev_q
 integer :: idir,ipert,my_rank
 character(len=500) :: msg
 type(hdr_type) :: hdr1,hdr_ref
!arrays
 integer,allocatable :: tmp_pos(:,:,:)
 real(dp),allocatable :: tmp_qpts(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm)
 db%path = path; db%comm = comm; db%iomode = IO_MODE_FORTRAN

 ! Master reads the header and builds useful tables
 if (my_rank == master) then

   if (open_file(path, msg, newunit=unt, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   read(unt, err=10, iomsg=msg) db%version
   read(unt, err=10, iomsg=msg) db%numv1

   ! Get important dimensions from the first header and rewind the file.
   call hdr_fort_read(hdr_ref, unt, fform)
   if (dvdb_check_fform(fform, "read_dvdb", msg) /= 0) then
     MSG_ERROR(sjoin("While reading:", path, ch10, msg))
   end if
   if (db%debug) call hdr_echo(hdr_ref,fform,4,unit=std_out)

   rewind(unt)
   read(unt, err=10, iomsg=msg)
   read(unt, err=10, iomsg=msg)

   ! The code below must be executed by the other procs if MPI.
   db%natom = hdr_ref%natom
   db%natom3 = 3 * hdr_ref%natom
   db%nspden = hdr_ref%nspden
   db%nsppol = hdr_ref%nsppol
   db%nspinor = hdr_ref%nspinor
   db%usepaw = hdr_ref%usepaw
   ABI_CHECK(db%usepaw == 0, "PAW not yet supported")

   ! TODO: Write function that returns mpert from natom!
   db%mpert = db%natom + 6

   ABI_MALLOC(tmp_qpts, (3, db%numv1))
   ABI_MALLOC(tmp_pos, (3, db%mpert, db%numv1))
   tmp_pos = 0

   ABI_MALLOC(db%cplex_v1, (db%numv1))
   ABI_MALLOC(db%ngfft3_v1, (3, db%numv1))
   ABI_MALLOC(db%iv_pinfoq, (4, db%numv1))
   ABI_MALLOC(db%rhog1_g0, (2, db%numv1))

   nqpt = 0
   do iv1=1,db%numv1
     call hdr_fort_read(hdr1, unt, fform)
     if (dvdb_check_fform(fform, "read_dvdb", msg) /= 0) then
       MSG_ERROR(sjoin("While reading hdr of v1 potential of index:", itoa(iv1), ch10, msg))
     end if

     ! Save cplex and FFT mesh associated to this perturbation.
     cplex = 2; if (hdr1%qptn(1)**2+hdr1%qptn(2)**2+hdr1%qptn(3)**2<1.d-14) cplex = 1
     db%cplex_v1(iv1) = cplex
     db%ngfft3_v1(:, iv1) = hdr1%ngfft(:3)

     ! Skip the records with v1.
     do ii=1,hdr1%nspden
       read(unt, err=10, iomsg=msg)
     end do
     ! Read rhog1_g0 (if available)
     db%rhog1_g0(:, iv1) = zero
     if (db%version > 1) read(unt, err=10, iomsg=msg) db%rhog1_g0(:, iv1)

     ! Check whether this q-point is already in the list.
     iq_found = 0
     do iq=1,nqpt
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
     db%iv_pinfoq(:,iv1) = [idir, ipert, hdr1%pertcase, iq_found]

     call hdr_free(hdr1)
   end do

   ! Allocate arrays with correct nqpt dimension
   db%nqpt = nqpt
   ABI_MALLOC(db%qpts, (3, nqpt))
   db%qpts = tmp_qpts(:,1:nqpt)
   ABI_FREE(tmp_qpts)

   ABI_MALLOC(db%pos_dpq, (3, db%mpert, nqpt))
   db%pos_dpq = tmp_pos(:, :, 1:nqpt)
   ABI_FREE(tmp_pos)

   close(unt)
 end if

 ! Master broadcasts data.
 if (xmpi_comm_size(comm) > 1) then
   call xmpi_bcast(db%version, master, comm, ierr)
   call xmpi_bcast(db%numv1, master, comm, ierr)
   call xmpi_bcast(db%nqpt, master, comm, ierr)
   call hdr_bcast(hdr_ref, master, my_rank, comm)

   db%natom = hdr_ref%natom
   db%natom3 = 3 * hdr_ref%natom
   db%nspden = hdr_ref%nspden
   db%nsppol = hdr_ref%nsppol
   db%nspinor = hdr_ref%nspinor
   db%usepaw = hdr_ref%usepaw
   db%mpert = db%natom + 6

   if (my_rank /= master) then
     ABI_MALLOC(db%cplex_v1, (db%numv1))
     ABI_MALLOC(db%ngfft3_v1, (3, db%numv1))
     ABI_MALLOC(db%iv_pinfoq, (4, db%numv1))
     ABI_MALLOC(db%qpts, (3, db%nqpt))
     ABI_MALLOC(db%pos_dpq, (3, db%mpert, db%nqpt))
     ABI_MALLOC(db%rhog1_g0, (2, db%numv1))
   end if

   call xmpi_bcast(db%cplex_v1, master, comm, ierr)
   call xmpi_bcast(db%ngfft3_v1, master, comm, ierr)
   call xmpi_bcast(db%iv_pinfoq, master, comm, ierr)
   call xmpi_bcast(db%qpts, master, comm, ierr)
   call xmpi_bcast(db%pos_dpq, master, comm, ierr)
   call xmpi_bcast(db%rhog1_g0, master, comm, ierr)
 end if

 ! Init crystal_t from the hdr read from file.
 call crystal_from_hdr(db%cryst,hdr_ref,timrev2)
 call hdr_free(hdr_ref)

 ! Init Born effective charges
 ABI_MALLOC(db%zeff, (3, 3, db%natom))

 ! Internal MPI_type needed for calling fourdp!
 call initmpi_seq(db%mpi_enreg)

 ! Precompute symq_table for all q-points in the DVDB.
 ABI_MALLOC(db%symq_table, (4, 2, db%cryst%nsym, db%nqpt))
 do iq=1,db%nqpt
   call littlegroup_q(db%cryst%nsym, db%qpts(:,iq), db%symq_table(:,:,:,iq), &
     db%cryst%symrec, db%cryst%symafm, trev_q, prtvol=0)
 end do

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(sjoin("Error while reading:", path, ch10, msg))

end subroutine dvdb_init
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_open_read(db, ngfft, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_open_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: nprocs
 character(len=500) :: msg

!************************************************************************

 if (db%rw_mode /= DVDB_NOMODE) then
   MSG_ERROR("DVDB should be in DVDB_NOMODE when open_read is called.")
 end if
 db%rw_mode = DVDB_READMODE

 nprocs = xmpi_comm_size(comm)

 ! Initialize tables to call fourdp in sequential
 db%ngfft = ngfft
 call init_distribfft_seq(db%mpi_enreg%distribfft,'c',ngfft(2),ngfft(3),'all')
 call init_distribfft_seq(db%mpi_enreg%distribfft,'f',ngfft(2),ngfft(3),'all')

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_close'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dvdb_t),intent(inout) :: db

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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_free(db)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dvdb_t),intent(inout) :: db

!************************************************************************

 ! integer arrays
 if (allocated(db%pos_dpq)) then
   ABI_FREE(db%pos_dpq)
 end if
 if (allocated(db%cplex_v1)) then
   ABI_FREE(db%cplex_v1)
 end if
 if (allocated(db%symq_table)) then
   ABI_FREE(db%symq_table)
 end if
 if (allocated(db%iv_pinfoq)) then
   ABI_FREE(db%iv_pinfoq)
 end if
 if (allocated(db%ngfft3_v1)) then
   ABI_FREE(db%ngfft3_v1)
 end if

 ! real arrays
 if (allocated(db%qpts)) then
   ABI_FREE(db%qpts)
 end if
 if (allocated(db%rpt)) then
   ABI_FREE(db%rpt)
 end if
 if (allocated(db%v1scf_rpt)) then
   ABI_FREE(db%v1scf_rpt)
 end if
 if (allocated(db%rhog1_g0)) then
   ABI_FREE(db%rhog1_g0)
 end if
 if (allocated(db%zeff)) then
   ABI_FREE(db%zeff)
 end if

 ! types
 call crystal_free(db%cryst)
 call destroy_mpi_enreg(db%mpi_enreg)

 ! Close the file but only if we have performed IO.
 if (db%rw_mode == DVDB_NOMODE) return
 call dvdb_close(db)

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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_print(db, header, unit, prtvol, mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(dvdb_t),intent(in) :: db

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol,iv1,iq,idir,ipert,iatom
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt = std_out; if (PRESENT(unit)) my_unt = unit
 my_prtvol = 0   ; if (PRESENT(prtvol)) my_prtvol = prtvol
 my_mode = 'COLL'; if (PRESENT(mode_paral)) my_mode = mode_paral

 msg=' ==== Info on the dvdb% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,"(a)")sjoin("DVDB version:", itoa(db%version))
 write(std_out,"(a)")sjoin("File path:", db%path)
 write(std_out,"(a)")sjoin("Number of v1scf potentials:", itoa(db%numv1))
 write(std_out,"(a)")sjoin("Number of q-points in DVDB: ", itoa(db%nqpt))
 write(std_out,"(a)")sjoin("Activate symmetrization of v1scf(r):", yesno(db%symv1))
 write(std_out,"(a)")"List of q-points: min(10, nqpt)"
 do iq=1,min(db%nqpt, 10)
   write(std_out,"(a)")sjoin("[", itoa(iq),"]", ktoa(db%qpts(:,iq)))
 end do
 if (db%nqpt > 10) write(std_out,"(a)")"..."

 write(std_out,"(a)")sjoin("Have dielectric tensor and Born effective charges:", yesno(db%has_dielt_zeff))
 if (db%has_dielt_zeff) then
   write(std_out, '(a,3(/,3es16.6))') ' Dielectric Tensor:', &
     db%dielt(1,1), db%dielt(1,2), db%dielt(1,3), &
     db%dielt(2,1), db%dielt(2,2), db%dielt(2,3), &
     db%dielt(3,1), db%dielt(3,2), db%dielt(3,3)
   write(std_out, '(a)') ' Effectives Charges: '
   do iatom=1,db%natom
     write(std_out,'(a,i0,3(/,3es16.6))')' iatom: ',iatom, &
       db%zeff(1,1,iatom), db%zeff(1,2,iatom), db%zeff(1,3,iatom), &
       db%zeff(2,1,iatom), db%zeff(2,2,iatom), db%zeff(2,3,iatom), &
       db%zeff(3,1,iatom), db%zeff(3,2,iatom), db%zeff(3,3,iatom)
   end do
 end if

 if (my_prtvol > 0) then
   call crystal_print(db%cryst, header="Crystal structure in DVDB file")
   write(std_out,"(a)")"FFT mesh for potentials on file:"
   write(std_out,"(a)")"q-point, idir, ipert, ngfft(:3)"
   do iv1=1,db%numv1
     idir = db%iv_pinfoq(1, iv1); ipert = db%iv_pinfoq(2, iv1); iq = db%iv_pinfoq(4, iv1)
     write(std_out,"(a)")sjoin(ktoa(db%qpts(:,iq)), itoa(idir), itoa(ipert), ltoa(db%ngfft3_v1(:,iv1)))
   end do

   write(std_out, "(a)")"q-point, idir, ipert, rhog1(q,G=0)"
   do iv1=1,db%numv1
     idir = db%iv_pinfoq(1, iv1); ipert = db%iv_pinfoq(2, iv1); iq = db%iv_pinfoq(4, iv1)
     write(std_out,"(a)")sjoin(ktoa(db%qpts(:,iq)), itoa(idir), itoa(ipert), &
       ftoa(db%rhog1_g0(1, iv1)), ftoa(db%rhog1_g0(2, iv1)))
   end do
 end if

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_get_pinfo'
!End of the abilint section

 implicit none

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
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_read_onev1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,iqpt,cplex,nfft
 character(len=*),intent(out) :: msg
 type(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: v1scf(cplex*nfft,db%nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: paral_kgb0=0
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
   write(msg,"(3(a,i0))")"Cannot find idir=",idir,", ipert=",ipert,", iqpt=",iqpt
   return
 end if

 if (cplex /= db%cplex_v1(iv1)) then
   write(msg,"(2(a,i0))")"Wrong cplex. Expecting=",db%cplex_v1(iv1),", received= ",cplex
   return
 end if

 ! Find (idir, ipert, iqpt) and skip the header.
 call dvdb_seek(db, idir, ipert, iqpt)
 ierr = my_hdr_skip(db%fh, idir, ipert, db%qpts(:,iqpt), msg)
 if (ierr /= 0) return

 ! Read v1 from file.
 nfftot_out = product(ngfft(:3)); nfftot_file = product(db%ngfft3_v1(:3, iv1))

 if (all(ngfft(:3) == db%ngfft3_v1(:3, iv1))) then
   do ispden=1,db%nspden
     read(db%fh, err=10, iomsg=msg) (v1scf(ifft, ispden), ifft=1,cplex*nfftot_file)
   end do
 else
   ! The FFT mesh used in the caller differ from the one found in the DBDB --> Fourier interpolation
   ! TODO: Add linear interpolation as well.
   !MSG_ERROR("FFT interpolation of DFPT potentials must be tested.")
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
     paral_kgb0,MPI_enreg_seq,v1r_file,v1scf,v1g_in,v1g_out)

   ABI_FREE(v1g_in)
   ABI_FREE(v1g_out)
   ABI_FREE(v1r_file)
   call destroy_mpi_enreg(MPI_enreg_seq)
 end if

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
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_readsym_allv1(db, iqpt, cplex, nfft, ngfft, v1scf, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_readsym_allv1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqpt,nfft,comm
 integer,intent(out) :: cplex
 type(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),allocatable,intent(out) :: v1scf(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ipc,npc,idir,ipert,pcase,my_rank,nproc,ierr,mu
 character(len=500) :: msg
!arrays
 integer :: pinfo(3,3*db%mpert),pflag(3, db%natom)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Get number of perturbations computed for this iqpt as well as cplex.
 npc = dvdb_get_pinfo(db, iqpt, cplex, pinfo)
 ABI_CHECK(npc /= 0, "npc == 0!")

 ABI_STAT_MALLOC(v1scf, (cplex,nfft,db%nspden,3*db%natom), ierr)
 ABI_CHECK(ierr==0, "OOM in v1scf")

 ! Master read all available perturbations and broadcasts data.
 if (my_rank == master) then
   do ipc=1,npc
     idir = pinfo(1,ipc); ipert = pinfo(2,ipc); pcase = pinfo(3, ipc)
     if (dvdb_read_onev1(db, idir, ipert, iqpt, cplex, nfft, ngfft, v1scf(:,:,:,pcase), msg) /= 0) then
       MSG_ERROR(msg)
     end if
   end do
 end if
 call xmpi_bcast(v1scf, master, comm, ierr)

 ! Return if all perts are available.
 if (npc == 3*db%natom) then
   if (db%symv1) then
     if (db%debug) write(std_out,*)"Potentials are available but will call v1phq_symmetrize because of symv1"
     do mu=1,db%natom3
       idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1
       call v1phq_symmetrize(db%cryst,idir,ipert,db%symq_table(:,:,:,iqpt),ngfft,cplex,nfft,&
         db%nspden,db%nsppol,db%mpi_enreg,v1scf(:,:,:,mu))
     end do
   end if
   if (db%debug) write(std_out,*)ABI_FUNC,": All perts available. Returning"
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

end subroutine dvdb_readsym_allv1
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
!!  symv1=If True, the new potentials are symmetrized using the set of symmetries that leaves the
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine v1phq_complete(cryst,qpt,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,symv1,pflag,v1scf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'v1phq_complete'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,nsppol
 logical,intent(in) :: symv1
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
 integer :: symrel_eq(3,3),symrec_eq(3,3),symm(3,3),g0_qpt(3),l0(3),tsm1g(3)
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
   call find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt)
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
   call mati3inv(symrel_eq, symm); symm = transpose(symm)

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
   !write(std_out,*)"Reconstructing idir, ipert",idir,ipert
   v1g = zero; cnt = 0
   do idir_eq=1,3
     if (symrec_eq(idir, idir_eq) == 0) cycle
     cnt = cnt + 1
     pcase_eq = idir_eq + (ipert_eq-1)*3
     if (debug) write(std_out,*)"idir_eq, ipert_eq, tsign",idir_eq, ipert_eq, tsign

     do ispden=1,nspden
       ! Get symmetric perturbation in G-space in workg_eq array.
       call fourdp(cplex,workg_eq,v1scf(:,ispden,pcase_eq),-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)

       !call rotate_fqg(itirev_eq,symrec_eq,qpt,tnon,ngfft,nfft,nspden,workg_eq,workg)
       ind1=0
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ind1=ind1+1

             ! Get location of G vector (grid point) centered at 0 0 0
             l1 = i1-(i1/id1)*n1-1
             l2 = i2-(i2/id2)*n2-1
             l3 = i3-(i3/id3)*n3-1

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
   if (debug) write(std_out,*)"Used ",cnt," equivalent perturbations"

   ! Get potential in real space (results in v1scf)
   do ispden=1,nspden
     call fourdp(cplex,v1g(:,:,ispden),v1scf(:,ispden,pcase),+1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)

     ! IS(q) = q + G0
     ! we want q so we have to multiply by exp(iG0r) in real space.
     if (any(g0_qpt /= 0)) then
       ABI_CHECK(cplex==2, "cplex == 1")
       if (debug) write(std_out,*)"Found not zero g0_qpt for idir, ipert:",idir,ipert
       call times_eigr(g0_qpt, ngfft, nfft, 1, v1scf(:,ispden,pcase))
     end if
   end do

   if (symv1) then
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
     "See message above for further information."
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_symeq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert
 integer,intent(out) :: ipert_eq,isym_eq,itirev_eq
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: symq(4,2,cryst%nsym),pflag(3,cryst%natom)
 integer,intent(out) :: g0_qpt(3)

!Local variables-------------------------------
!scalars
 integer :: isym,idir_eq,ip,itirev

! *************************************************************************

 isym_eq = -1; ipert_eq = -1

symloop: &
 do itirev=1,2
   itirev_eq = itirev
   do isym=1,cryst%nsym

     ! Check that isym preserves the q-point
     !if (symq(4,itirev,isym) /= 1 .or. any(symq(1:3,itirev,isym) /= 0)) cycle
     if (symq(4,itirev,isym) /= 1) cycle ! .or. any(symq(1:3,itirev,isym) /= 0)) cycle
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
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  nspden=number of spin-density components
!!  nsppol=Number of independent spin polarizations
!!  mpi_enreg=information about MPI parallelization
!!  v1r_qibz(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space
!!   for the irreducible q-point `qpt_ibz`
!!
!! OUTPUT
!!  v1r_qbz(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space for the q-point in the BZ
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine v1phq_rotate(cryst,qpt_ibz,isym,itimrev,g0q,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r_qibz,v1r_qbz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'v1phq_rotate'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isym,itimrev,cplex,nfft,nspden,nsppol
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: g0q(3),ngfft(18)
 real(dp),intent(in) :: qpt_ibz(3)
 real(dp),intent(inout) :: v1r_qibz(cplex*nfft,nspden,3*cryst%natom)
 real(dp),intent(out) :: v1r_qbz(cplex*nfft,nspden,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0
 integer,save :: enough=0
 integer :: natom3,mu,ispden,idir,ipert,idir_eq,ipert_eq,mu_eq,cnt,tsign
!arrays
 integer :: symrec_eq(3,3),sm1(3,3),l0(3) !g0_qpt(3), symrel_eq(3,3),
 real(dp) :: tnon(3)
 real(dp),allocatable :: v1g_qibz(:,:,:),workg(:,:),v1g_mu(:,:)

! *************************************************************************

 ABI_UNUSED(nsppol)
 ABI_CHECK(cplex == 2, "cplex != 2")

 natom3 = 3 * cryst%natom; tsign = 3-2*itimrev

 ! Compute IBZ potentials in G-space (results in v1g_qibz)
 ABI_MALLOC(v1g_qibz, (2*nfft,nspden,natom3))
 do mu=1,natom3
   do ispden=1,nspden
     call fourdp(cplex,v1g_qibz(:,ispden,mu),v1r_qibz(:,ispden,mu),-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)
   end do
 end do

 ABI_MALLOC(workg, (2*nfft,nspden))
 ABI_MALLOC(v1g_mu, (2*nfft,nspden))

 ! For each perturbation:
 ! FIXME: This is wrong if natom > 1
 symrec_eq = cryst%symrec(:,:,isym)
 call mati3inv(symrec_eq, sm1); sm1 = transpose(sm1)

 do mu=1,natom3
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

     ! Rotate in G-space and accumulate in workg
     call rotate_fqg(itimrev,sm1,qpt_ibz,tnon,ngfft,nfft,nspden,v1g_qibz(:,:,mu_eq),workg)
     v1g_mu = v1g_mu + workg * symrec_eq(idir, idir_eq)
   end do ! idir_eq
   ABI_CHECK(cnt /= 0, "cnt should not be zero!")

   ! Transform to real space and take into account a possible shift.
   ! (results are stored in v1r_qbz)
   do ispden=1,nspden
     call fourdp(cplex,v1g_mu(:,ispden),v1r_qbz(:,ispden,mu),+1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)
     call times_eigr(-g0q, ngfft, nfft, 1, v1r_qbz(:,ispden,mu))
     !call times_eigr(tsign * g0q, ngfft, nfft, 1, v1r_qbz(:,ispden,mu))
   end do

   !call v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r)
 end do ! mu

 ABI_FREE(workg)
 ABI_FREE(v1g_mu)
 ABI_FREE(v1g_qibz)

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
!!  v1r(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space
!!   Symmetrized in output.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'v1phq_symmetrize'
!End of the abilint section

 implicit none

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
 !  call pawang_init(pawang1,0,pawang%l_max-1,0,nsym1,0,1,0,0,0)
 !  call setsym_ylm(gprimd,pawang1%l_max-1,pawang1%nsym,0,rprimd,symrc1,pawang1%zarot)
 !end if

 ! FIXME Be careful here because symrhg was written for densities!
 ABI_CHECK(nsppol == 1 .and. nspden == 1, "symrhg was written for densities, not for potentials")

 ABI_MALLOC(v1g, (2,nfft))
 call symrhg(cplex,cryst%gprimd,irrzon1,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym1,&
    mpi_enreg%paral_kgb,phnons1,v1g,v1r,cryst%rprimd,symafm1,symrel1)

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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine rotate_fqg(itirev, symm, qpt, tnon, ngfft, nfft, nspden, infg, outfg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_fqg'
!End of the abilint section

 implicit none

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
 real(dp) :: phnon1(2)

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")
 id1 = n1/2+2; id2 = n2/2+2; id3=n3/2+2

 ABI_CHECK(any(itirev == [1,2]), "Wrong itirev")
 tsign = 3-2*itirev; has_phase = any(abs(tnon) > tol12)

 do isp=1,nspden
   ind1 = 0
   do i3=1,n3
     do i2=1,n2
       do i1=1,n1
         ind1 = ind1 + 1

         ! Get location of G vector (grid point) centered at 0 0 0
         l1 = i1-(i1/id1)*n1-1
         l2 = i2-(i2/id2)*n2-1
         l3 = i3-(i3/id3)*n3-1

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
           !write(std_out,*)"got it"
           outfg(:,ind1,isp) = zero
           cycle
         end if

         tsg = [j1,j2,j3] ! +- S^{-1} G

         ! Map into [0,n-1] and then add 1 for array index in [1,n]
         k1=1+mod(n1+mod(j1,n1),n1)
         k2=1+mod(n2+mod(j2,n2),n2)
         k3=1+mod(n3+mod(j3,n3),n3)

         ! Get linear index of rotated point Gj
         ind2 = k1+n1*((k2-1)+n2*(k3-1))

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

end subroutine rotate_fqg
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_ftinterp_setup
!! NAME
!!  dvdb_ftinterp_setup
!!
!! FUNCTION
!!  Prepare the internal tables used for the Fourier interpolation of the DFPT potentials.
!!
!! INPUTS
!!  ngqpt(3)=Divisions of the ab-initio q-mesh.
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  comm=MPI communicator
!!
!! PARENTS
!!      m_dvdb,m_phgamma,m_sigmaph
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_ftinterp_setup(db,ngqpt,nqshift,qshift,nfft,ngfft,comm,cryst_op)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_ftinterp_setup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,comm
 type(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngqpt(3),ngfft(18)
 real(dp),intent(in) :: qshift(3,nqshift)
 type(crystal_t),intent(in),target,optional :: cryst_op

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1,timrev1=1,tim_fourdp0=0
 integer :: my_qptopt,iq_ibz,nqibz,iq_bz,nqbz
 integer :: ii,iq_dvdb,cplex_qibz,ispden,mu,irpt,idir,ipert
 integer :: iqst,nqst,itimrev,tsign,isym,ix,iy,iz,nq1,nq2,nq3,r1,r2,r3
 integer :: nproc,my_rank,ifft,cnt,ierr
 real(dp) :: dksqmax,phre,phim
 logical :: isirr_q, found
 type(crystal_t),pointer :: cryst
 type(mpi_type) :: mpi_enreg_seq
!arrays
 integer :: qptrlatt(3,3),g0q(3),ngfft_qspace(18)
 integer,allocatable :: indqq(:,:),iperm(:),bz2ibz_sort(:),nqsts(:),iqs_dvdb(:)
 real(dp) :: qpt_bz(3),shift(3) !,qpt_ibz(3)
 real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:),emiqr(:,:)
 real(dp),allocatable :: v1r_qibz(:,:,:,:),v1r_qbz(:,:,:,:), v1r_lr(:,:,:) !,all_v1qr(:,:,:,:)

! *************************************************************************

 if (allocated(db%v1scf_rpt)) then
   if (db%debug) then
     call wrtout(std_out, "v1scf_rpt is already computed. Returning")
   end if
   return
 end if

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 if (present(cryst_op)) then
   cryst => cryst_op
 else
   cryst => db%cryst
 end if
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
 db%nrpt = nqbz
 ABI_MALLOC(db%rpt, (3, db%nrpt))
 ii = 0
 do iz=1,nq3
   r3 = ig2gfft(iz,nq3) !* nq3
   do iy=1,nq2
     r2 = ig2gfft(iy,nq2) !* nq2
     do ix=1,nq1
       r1 = ig2gfft(ix,nq1) !* nq1
       ii = ii + 1
       db%rpt(:,ii) = [r1, r2, r3]
     end do
   end do
 end do

 ! Find correspondence BZ --> IBZ. Note:
 !   * q --> -q symmetry is always used for phonons.
 !   * we use symrec instead of symrel
 ABI_MALLOC(indqq, (nqbz*sppoldbl1,6))
 call listkk(dksqmax,cryst%gmet,indqq,qibz,qbz,nqibz,nqbz,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

 ! Construct sorted mapping BZ --> IBZ to speedup qbz search below.
 ABI_MALLOC(iperm, (nqbz))
 ABI_MALLOC(bz2ibz_sort, (nqbz))
 iperm = [(ii, ii=1,nqbz)]
 bz2ibz_sort = indqq(:,1)
 call sort_int(nqbz,bz2ibz_sort,iperm)

 ! Reconstruct the IBZ according to what is present in the DVDB .
 ABI_MALLOC(nqsts, (nqibz))
 ABI_MALLOC(iqs_dvdb, (nqibz))

 ABI_MALLOC(v1r_lr, (2,nfft,db%natom3))

 iqst = 0
 do iq_ibz=1,nqibz

   ! In each q-point star, count the number of q-points
   ! and find the one present in DVDB.
   nqst = 0
   found = .false.
   do ii=iqst+1,nqbz
     if (bz2ibz_sort(ii) /= iq_ibz) exit
     nqst = nqst + 1

     iq_bz = iperm(ii)
     if (.not. found) then
       iq_dvdb = dvdb_findq(db, qbz(:,iq_bz))
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
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)

 ABI_MALLOC(emiqr, (2,db%nrpt))
 ABI_STAT_MALLOC(db%v1scf_rpt,(2,db%nrpt, nfft, db%nspden, db%natom3), ierr)
 ABI_CHECK(ierr==0, 'out of memory in db%v1scf_rpt')
 db%v1scf_rpt = zero

 ABI_MALLOC(v1r_qbz, (2, nfft, db%nspden, db%natom3))
 !v1r_qbz = huge(one)

 iqst = 0
 do iq_ibz=1,nqibz

   ! Get potentials for this IBZ q-point on the real-space FFT mesh.
   ! This call allocates v1r_qibz(cplex_qibz, nfft, nspden, 3*natom)
   call dvdb_readsym_allv1(db, iqs_dvdb(iq_ibz), cplex_qibz, nfft, ngfft, v1r_qibz, comm)

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
     v1r_lr = zero; cnt = 0
     if (db%has_dielt_zeff) then
       do mu=1,db%natom3
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.
         idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1
         call dvdb_v1r_long_range(db,qpt_bz,ipert,idir,nfft,ngfft,v1r_lr(:,:,mu))
       end do
       call xmpi_sum(v1r_lr, comm, ierr)
     end if

     if (cplex_qibz == 1) then
       ! Gamma point.
       ABI_CHECK(nqsts(iq_ibz) == 1, "cplex_qibz == 1 and nq nqst /= 1 (should be gamma)")
       ABI_CHECK(all(g0q == 0), "gamma point with g0q /= 0")

       ! Substract the long-range part of the potential
       do mu=1,db%natom3
         do ispden=1,db%nspden
           v1r_qibz(1, :, ispden, mu) = v1r_qibz(1, :, ispden, mu) - v1r_lr(1, :, mu)
         end do
       end do

       ! Slow FT.
       cnt = 0
       do mu=1,db%natom3
         do ispden=1,db%nspden
           do irpt=1,db%nrpt
             ! MPI-parallelism
             cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
             do ifft=1,nfft
               db%v1scf_rpt(1, irpt, ifft, ispden, mu) = db%v1scf_rpt(1, irpt, ifft, ispden, mu) + &
                                                         v1r_qibz(1, ifft, ispden, mu)
             end do
           end do
         end do
       end do

     else
       ! Get the periodic part of the potential in BZ (v1r_qbz)
       if (isirr_q) then
         !write(std_out,*)sjoin("qpt irred:",ktoa(qpt_bz))
         v1r_qbz = v1r_qibz
       else
         call v1phq_rotate(cryst, qibz(:,iq_ibz), isym, itimrev, g0q,&
           ngfft, cplex_qibz, nfft, db%nspden, db%nsppol, db%mpi_enreg, v1r_qibz, v1r_qbz)
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
       do mu=1,db%natom3
         do ispden=1,db%nspden
           v1r_qbz(1, :, ispden, mu) = v1r_qbz(1, :, ispden, mu) - v1r_lr(1, :, mu)
           v1r_qbz(2, :, ispden, mu) = v1r_qbz(2, :, ispden, mu) - v1r_lr(2, :, mu)
         end do
       end do

       ! Compute FT phases for this qpt_bz.
       call calc_eiqr(-qpt_bz,db%nrpt,db%rpt,emiqr)

       ! Slow FT.
       cnt = 0
       do mu=1,db%natom3

         do ispden=1,db%nspden
           do irpt=1,db%nrpt
             cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.
             phre = emiqr(1,irpt); phim = emiqr(2,irpt)

             do ifft=1,nfft
               db%v1scf_rpt(1,irpt, ifft, ispden, mu) = db%v1scf_rpt(1, irpt, ifft, ispden, mu) + &
                  phre * v1r_qbz(1, ifft, ispden, mu) - phim * v1r_qbz(2, ifft, ispden, mu)

               db%v1scf_rpt(2, irpt, ifft, ispden, mu) = db%v1scf_rpt(2, irpt, ifft, ispden, mu) + &
                  phre * v1r_qbz(2, ifft, ispden, mu) + phim * v1r_qbz(1, ifft, ispden, mu)
             end do

           end do
         end do
       end do
     end if

   end do ! iqst

   ABI_FREE(v1r_qibz)
 end do ! iq_ibz
 ABI_CHECK(iqst == nqbz, "iqst /= nqbz")

 db%v1scf_rpt = db%v1scf_rpt / nqbz
 call xmpi_sum(db%v1scf_rpt, comm, ierr)

 ! Build mpi_type for executing fourdp in sequential.
 call ngfft_seq(ngfft_qspace, [nq1, nq2, nq3])
 call initmpi_seq(mpi_enreg_seq)
 call init_distribfft_seq(mpi_enreg_seq%distribfft,'c',ngfft_qspace(2),ngfft_qspace(3),'all')
 call init_distribfft_seq(mpi_enreg_seq%distribfft,'f',ngfft_qspace(2),ngfft_qspace(3),'all')

 !ABI_STAT_MALLOC(all_v1qr, (nqbz, 2, nfft, db%nspden, db%natom3), ierr)
 !ABI_CHECK(ierr==0, "oom in all_v1qr")
 !all_v1qr = zero

 cnt = 0
 do mu=1,db%natom3
   do ispden=1,db%nspden
     do ifft=1,nfft
       !cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
       !call fourdp(1,all_v1qr(:,:,ifft,ispden,mu),db%v1scf_rpt(:,ifft,ispden,mu),+1,&
       ! mpi_enreg_seq,nqbz,ngfft_qspace,paral_kgb0,tim_fourdp0)
     end do
   end do
 end do

 !ABI_FREE(all_v1qr)
 call destroy_mpi_enreg(mpi_enreg_seq)

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


end subroutine dvdb_ftinterp_setup
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
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  comm=MPI communicator
!!
!! OUTPUT
!!  ov1r(2*nfft, nspden, 3*natom)=Interpolated DFPT potentials at the given q-point.
!!
!! PARENTS
!!      m_dvdb,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_ftinterp_qpt(db, qpt, nfft, ngfft, ov1r, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_ftinterp_qpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,comm
 type(dvdb_t),intent(in) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: ov1r(2,nfft,db%nspden,db%natom3)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex2=2
 integer :: ir,ispden,ifft,mu,idir,ipert,timerev_q,nproc,my_rank,cnt,ierr
 real(dp) :: wr,wi
!arrays
 integer :: symq(4,2,db%cryst%nsym)
 real(dp),allocatable :: eiqr(:,:), v1r_lr(:,:,:)

! *************************************************************************

 ABI_CHECK(allocated(db%v1scf_rpt), "v1scf_rpt is not allocated (call dvdb_ftinterp_setup)")

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ABI_MALLOC(v1r_lr, (2,nfft,db%natom3))

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(db%cryst%nsym,qpt,symq,db%cryst%symrec,db%cryst%symafm,timerev_q,prtvol=db%prtvol)

 ! Compute FT phases for this q-point.
 ABI_MALLOC(eiqr, (2,db%nrpt))
 call calc_eiqr(qpt,db%nrpt,db%rpt,eiqr)

 ! Compute long-range part of the coupling potential
 v1r_lr = zero; cnt = 0
 if (db%has_dielt_zeff) then
   do mu=1,db%natom3
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI-parallelism
     idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1
     call dvdb_v1r_long_range(db,qpt,ipert,idir,nfft,ngfft,v1r_lr(:,:,mu))
   end do
   call xmpi_sum(v1r_lr, comm, ierr)
 end if

 ! TODO: If high-symmetry q-points, one could save flops by FFT interpolating the independent
 ! perturbations and then rotate ...
 ov1r = zero; cnt = 0
 do mu=1,db%natom3
   idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1

   do ispden=1,db%nspden
     do ifft=1,nfft
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI-parallelism

       do ir=1,db%nrpt
         wr = db%v1scf_rpt(1, ir,ifft,ispden,mu)
         wi = db%v1scf_rpt(2, ir,ifft,ispden,mu)
         ov1r(1,ifft,ispden,mu) = ov1r(1,ifft,ispden,mu) + wr*eiqr(1,ir) - wi * eiqr(2,ir)
         ov1r(2,ifft,ispden,mu) = ov1r(2,ifft,ispden,mu) + wr*eiqr(2,ir) + wi * eiqr(1,ir)
       end do

       ! Add the long-range part of the potential
       ov1r(1,ifft,ispden,mu) = ov1r(1,ifft,ispden,mu) + v1r_lr(1,ifft,mu)
       ov1r(2,ifft,ispden,mu) = ov1r(2,ifft,ispden,mu) + v1r_lr(2,ifft,mu)

     end do

     ! Remove the phase.
     call times_eikr(-qpt, ngfft, nfft, 1, ov1r(:,:,ispden,mu))
   end do

   ! Be careful with gamma and cplex!
   if (db%symv1) then
     call v1phq_symmetrize(db%cryst,idir,ipert,symq,ngfft,cplex2,nfft,db%nspden,db%nsppol,db%mpi_enreg,ov1r(:,:,:,mu))
   end if
 end do ! mu

 call xmpi_sum(ov1r, comm, ierr)

 ABI_FREE(eiqr)
 ABI_FREE(v1r_lr)


end subroutine dvdb_ftinterp_qpt
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_v1scf_rpt
!! NAME
!!  dvdb_get_v1scf_rpt
!!
!! FUNCTION
!!  Compute the phonon perturbation potential in real space lattice
!!  representation.
!!  This routine is meant to replace dvdb_ftinterp_setup
!!  and performs the potential interpolation one perturbation at a time.
!!
!! INPUTS
!!  ngqpt(3)=Divisions of the ab-initio q-mesh.
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_get_v1scf_rpt(db, cryst, ngqpt, nqshift, qshift, nfft, ngfft, &
&                             nrpt, nspden, ipert, v1scf_rpt, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_get_v1scf_rpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,nrpt,nspden,ipert,comm
 type(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngqpt(3),ngfft(18)
 real(dp),intent(in) :: qshift(3,nqshift)
 real(dp),intent(out) :: v1scf_rpt(2,nrpt,nfft,nspden)
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1,timrev1=1,tim_fourdp0=0
 integer :: my_qptopt,iq_ibz,nqibz,iq_bz,nqbz
 integer :: ii,iq_dvdb,cplex_qibz,ispden,irpt,idir,iat
 integer :: iqst,nqst,itimrev,tsign,isym,ix,iy,iz,nq1,nq2,nq3,r1,r2,r3
 integer :: nproc,my_rank,ifft,cnt,ierr
 real(dp) :: dksqmax,phre,phim
 logical :: isirr_q, found
!arrays
 integer :: qptrlatt(3,3),g0q(3)
 integer,allocatable :: indqq(:,:),iperm(:),bz2ibz_sort(:),nqsts(:),iqs_dvdb(:)
 real(dp) :: qpt_bz(3),shift(3) !,qpt_ibz(3)
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

 db%nrpt = nqbz

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
 ABI_MALLOC(db%rpt, (3, db%nrpt))
 ii = 0
 do iz=1,nq3
   r3 = ig2gfft(iz,nq3) !* nq3
   do iy=1,nq2
     r2 = ig2gfft(iy,nq2) !* nq2
     do ix=1,nq1
       r1 = ig2gfft(ix,nq1) !* nq1
       ii = ii + 1
       db%rpt(:,ii) = [r1, r2, r3]
     end do
   end do
 end do

 ! Find correspondence BZ --> IBZ. Note:
 !   * q --> -q symmetry is always used for phonons.
 !   * we use symrec instead of symrel
 ABI_MALLOC(indqq, (nqbz*sppoldbl1,6))
 call listkk(dksqmax,cryst%gmet,indqq,qibz,qbz,nqibz,nqbz,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

 ! Construct sorted mapping BZ --> IBZ to speedup qbz search below.
 ABI_MALLOC(iperm, (nqbz))
 ABI_MALLOC(bz2ibz_sort, (nqbz))
 iperm = [(ii, ii=1,nqbz)]
 bz2ibz_sort = indqq(:,1)
 call sort_int(nqbz,bz2ibz_sort,iperm)

 ! Reconstruct the IBZ according to what is present in the DVDB .
 ABI_MALLOC(nqsts, (nqibz))
 ABI_MALLOC(iqs_dvdb, (nqibz))

 ABI_MALLOC(v1r_lr, (2,nfft))

 iqst = 0
 do iq_ibz=1,nqibz

   ! In each q-point star, count the number of q-points
   ! and find the one present in DVDB.
   nqst = 0
   found = .false.
   do ii=iqst+1,nqbz
     if (bz2ibz_sort(ii) /= iq_ibz) exit
     nqst = nqst + 1

     iq_bz = iperm(ii)
     if (.not. found) then
       iq_dvdb = dvdb_findq(db, qbz(:,iq_bz))
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
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)

 ABI_MALLOC(emiqr, (2,db%nrpt))
 v1scf_rpt = zero

 ABI_MALLOC(v1r_qbz, (2, nfft, db%nspden, db%natom3))
 !v1r_qbz = huge(one)

 iqst = 0
 do iq_ibz=1,nqibz

   ! Get potentials for this IBZ q-point on the real-space FFT mesh.
   ! This call allocates v1r_qibz(cplex_qibz, nfft, nspden, 3*natom)
   call dvdb_readsym_allv1(db, iqs_dvdb(iq_ibz), cplex_qibz, nfft, ngfft, v1r_qibz, comm)

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
     v1r_lr = zero; cnt = 0
     if (db%has_dielt_zeff) then
       idir = mod(ipert-1, 3) + 1; iat = (ipert - idir) / 3 + 1
       call dvdb_v1r_long_range(db,qpt_bz,iat,idir,nfft,ngfft,v1r_lr)
     end if

     if (cplex_qibz == 1) then
       ! Gamma point.
       ABI_CHECK(nqsts(iq_ibz) == 1, "cplex_qibz == 1 and nq nqst /= 1 (should be gamma)")
       ABI_CHECK(all(g0q == 0), "gamma point with g0q /= 0")

       ! Substract the long-range part of the potential
       do ispden=1,db%nspden
         v1r_qibz(1,:,ispden,ipert) = v1r_qibz(1,:,ispden,ipert) - v1r_lr(1,:)
       end do

       ! Slow FT.
       cnt = 0
       do ispden=1,db%nspden
         do irpt=1,db%nrpt
           ! MPI-parallelism
           cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
           do ifft=1,nfft
             v1scf_rpt(1,irpt,ifft,ispden) = v1scf_rpt(1,irpt,ifft,ispden) + &
                                             v1r_qibz(1, ifft, ispden, ipert)
           end do
         end do
       end do

     else
       ! Get the periodic part of the potential in BZ (v1r_qbz)
       if (isirr_q) then
         !write(std_out,*)sjoin("qpt irred:",ktoa(qpt_bz))
         v1r_qbz = v1r_qibz
       else
         call v1phq_rotate(cryst, qibz(:,iq_ibz), isym, itimrev, g0q,&
           ngfft, cplex_qibz, nfft, db%nspden, db%nsppol, db%mpi_enreg, v1r_qibz, v1r_qbz)
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
       do ispden=1,db%nspden
         v1r_qbz(1,:,ispden,ipert) = v1r_qbz(1,:,ispden,ipert) - v1r_lr(1,:)
         v1r_qbz(2,:,ispden,ipert) = v1r_qbz(2,:,ispden,ipert) - v1r_lr(2,:)
       end do

       ! Compute FT phases for this qpt_bz.
       call calc_eiqr(-qpt_bz,db%nrpt,db%rpt,emiqr)

       ! Slow FT.
       cnt = 0
       do ispden=1,db%nspden
         do irpt=1,db%nrpt
           cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.
           phre = emiqr(1,irpt); phim = emiqr(2,irpt)

           do ifft=1,nfft
             v1scf_rpt(1,irpt,ifft,ispden) = v1scf_rpt(1,irpt,ifft,ispden) &
             &                      + phre * v1r_qbz(1, ifft, ispden, ipert) &
             &                      - phim * v1r_qbz(2, ifft, ispden, ipert)

             v1scf_rpt(2,irpt,ifft,ispden) = v1scf_rpt(2,irpt,ifft,ispden) &
             &                      + phre * v1r_qbz(2, ifft, ispden, ipert) &
             &                      + phim * v1r_qbz(1, ifft, ispden, ipert)
           end do

         end do
       end do
     end if

   end do ! iqst

   ABI_FREE(v1r_qibz)
 end do ! iq_ibz
 ABI_CHECK(iqst == nqbz, "iqst /= nqbz")

 v1scf_rpt = v1scf_rpt / nqbz
 call xmpi_sum(v1scf_rpt, comm, ierr)

 !! Build mpi_type for executing fourdp in sequential.
 !call ngfft_seq(ngfft_qspace, [nq1, nq2, nq3])
 !call initmpi_seq(mpi_enreg_seq)
 !call init_distribfft_seq(mpi_enreg_seq%distribfft,'c',ngfft_qspace(2),ngfft_qspace(3),'all')
 !call init_distribfft_seq(mpi_enreg_seq%distribfft,'f',ngfft_qspace(2),ngfft_qspace(3),'all')

 !!ABI_STAT_MALLOC(all_v1qr, (nqbz, 2, nfft, db%nspden, db%natom3), ierr)
 !!ABI_CHECK(ierr==0, "oom in all_v1qr")
 !!all_v1qr = zero

 !cnt = 0
 !do mu=1,db%natom3
 !  do ispden=1,db%nspden
 !    do ifft=1,nfft
 !      !cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
 !      !call fourdp(1,all_v1qr(:,:,ifft,ispden,mu),db%v1scf_rpt(:,ifft,ispden,mu),+1,&
 !      ! mpi_enreg_seq,nqbz,ngfft_qspace,paral_kgb0,tim_fourdp0)
 !    end do
 !  end do
 !end do

 !!ABI_FREE(all_v1qr)
 !call destroy_mpi_enreg(mpi_enreg_seq)

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
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nrpt=Number of R-points = number of q-points in the full BZ
!!  nspden=Number of spin densities.
!!  ipert=index of the perturbation to be treated [1,natom3]
!!  v1scf_rpt(2,nrpt,nfft,nspden)=phonon perturbation potential
!!                                in real space lattice representation.
!!  comm=MPI communicator
!!
!! OUTPUT
!!  v1scf_qpt(2*nfft, nspden)=Interpolated DFPT potentials at the given q-point.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_get_v1scf_qpt(db, cryst, qpt, nfft, ngfft, nrpt, nspden, &
&                             ipert, v1scf_rpt, v1scf_qpt, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_get_v1scf_qpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nrpt,nspden,ipert,comm
 type(dvdb_t),intent(in) :: db
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
 call littlegroup_q(db%cryst%nsym,qpt,symq,db%cryst%symrec,db%cryst%symafm,timerev_q,prtvol=db%prtvol)

 ! Compute FT phases for this q-point.
 ABI_MALLOC(eiqr, (2,db%nrpt))
 call calc_eiqr(qpt,db%nrpt,db%rpt,eiqr)

 idir = mod(ipert-1, 3) + 1; iat = (ipert - idir) / 3 + 1

 ! Compute long-range part of the coupling potential
 v1r_lr = zero; cnt = 0
 if (db%has_dielt_zeff) then
   call dvdb_v1r_long_range(db,qpt,iat,idir,nfft,ngfft,v1r_lr)
 end if

 ! TODO: If high-symmetry q-points, one could save flops by FFT interpolating the independent
 ! perturbations and then rotate ...
 v1scf_qpt = zero; cnt = 0
 do ispden=1,db%nspden
   do ifft=1,nfft
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI-parallelism

     do ir=1,db%nrpt
       wr = v1scf_rpt(1, ir,ifft,ispden)
       wi = v1scf_rpt(2, ir,ifft,ispden)
       v1scf_qpt(1,ifft,ispden) = v1scf_qpt(1,ifft,ispden) + wr*eiqr(1,ir) - wi * eiqr(2,ir)
       v1scf_qpt(2,ifft,ispden) = v1scf_qpt(2,ifft,ispden) + wr*eiqr(2,ir) + wi * eiqr(1,ir)
     end do

     ! Add the long-range part of the potential
     v1scf_qpt(1,ifft,ispden) = v1scf_qpt(1,ifft,ispden) + v1r_lr(1,ifft)
     v1scf_qpt(2,ifft,ispden) = v1scf_qpt(2,ifft,ispden) + v1r_lr(2,ifft)

   end do

   ! Remove the phase.
   call times_eikr(-qpt, ngfft, nfft, 1, v1scf_qpt(:,:,ispden))
 end do

 ! Be careful with gamma and cplex!
 if (db%symv1) then
   call v1phq_symmetrize(db%cryst,idir,iat,symq,ngfft,cplex2,nfft,db%nspden,db%nsppol,db%mpi_enreg,v1scf_qpt(:,:,:))
 end if

 call xmpi_sum(v1scf_qpt, comm, ierr)

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
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_interpolate_v1scf(db, cryst, qpt, ngqpt, nqshift, qshift, &
&                                 nfft, ngfft, nfftf, ngfftf, v1scf, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_interpolate_v1scf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,nfftf,comm
 type(dvdb_t),target,intent(inout) :: db
!arrays
 real(dp),intent(in) :: qpt(3)
 integer,intent(in) :: ngqpt(3),ngfft(18),ngfftf(18)
 real(dp),intent(in) :: qshift(3,nqshift)
 real(dp),allocatable,intent(out) :: v1scf(:,:,:,:)
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer :: ipert, nqbz
 integer :: nproc,my_rank
!arrays
 real(dp),allocatable :: v1scf_rpt(:,:,:,:)

! *************************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 nqbz = product(ngqpt) * nqshift
 db%nrpt = nqbz

 ABI_MALLOC(v1scf, (2,nfftf,db%nspden,db%natom3))

 ABI_MALLOC(v1scf_rpt, (2,db%nrpt,nfft,db%nspden))

 do ipert=1,db%natom3

   write(std_out, "(a,i4,a,i4,a)") "Interpolating potential for perturbation ", &
   &                           ipert, " / ", db%natom3, ch10

   ! FIXME I think this should be ngfftf and not ngfft
   !       Also, other calls to dvdb_ftinterp_setup should use ngfftf.
   call dvdb_get_v1scf_rpt(db, cryst, ngqpt, nqshift, qshift, nfft, ngfft, &
   &                       db%nrpt, db%nspden, ipert, v1scf_rpt, comm)

   call dvdb_get_v1scf_qpt(db, cryst, qpt, nfftf, ngfftf, db%nrpt, db%nspden, &
   &                       ipert, v1scf_rpt, v1scf(:,:,:,ipert), comm)

   ABI_FREE(db%rpt)

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_findq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: qtol
 type(dvdb_t),intent(in) :: db
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer :: iq
 real(dp) :: my_qtol

! *************************************************************************

 my_qtol = 0.0001_dp; if (present(qtol)) my_qtol = qtol

 iqpt = -1
 do iq=1,db%nqpt
   if (all(abs(db%qpts(:, iq) - qpt) < my_qtol)) then
      iqpt = iq; exit
   end if
 end do

end function dvdb_findq
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_seek
!! NAME
!!  dvdb_seek
!!
!! FUNCTION
!!   Move the internal file pointer so that it points to the
!!   block with (idir, ipert, iqpt). Needed only if dvdb%iomode==IO_MODE_FORTRAN
!!
!! INPUTS
!!   idir,ipert,iqpt = (direction, perturbation, q-point) indices
!!
!! SIDE EFFECTS
!!   db<type(dvdb_t)> : modifies db%f90_fptr and the internal F90 file pointer.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_seek(db, idir, ipert, iqpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_seek'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: idir, ipert, iqpt
 type(dvdb_t),intent(inout) :: db

!Local variables-------------------------------
 integer :: pos_wanted,ii,ispden
 real(dp),parameter :: fake_qpt(3)=zero
 character(len=500) :: msg

! *************************************************************************

 if (db%iomode == IO_MODE_FORTRAN) then
   ! Numb code: rewind the file and read it from the beginning
   ! TODO: Here I need a routine to backspace the hdr!
   if (dvdb_rewind(db, msg) /= 0) then
     MSG_ERROR(msg)
   end if
   pos_wanted = db%pos_dpq(idir,ipert,iqpt)

   do ii=1,pos_wanted-1
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
 !ierr = 1
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_rewind'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_hdr_skip'
!End of the abilint section

 implicit none

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
         msg = "Perturbation index on file does not match the one expected by caller"
         ierr = -1
   end if
 end if

 call hdr_free(tmp_hdr)

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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_list_perts(db, ngqpt, unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_list_perts'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dvdb_t),target,intent(in) :: db
 integer,optional,intent(in) :: unit
!arrays
 integer,intent(in) :: ngqpt(3)

!Local variables-------------------------------
!scalars
 integer :: tot_miss,tot_weird,miss_q,idir,ipert,iv1,psy,weird_q
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
 tot_miss = 0; tot_weird = 0
 do iq_ibz=1,nqibz
   qq = qibz(:,iq_ibz)
   iq_file = dvdb_findq(db, qq)

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
     call wrtout(unt, sjoin("qpoint:", ktoa(qq), "is present in the DVDB file"))
     call wrtout(unt,' The list of irreducible perturbations for this q vector is:')
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
         write(msg,'(i5,a,i2,a,i4,4a)')ii,')  idir=',idir,', ipert=',ipert,", type=",trim(ptype),", found=",trim(found)
         call wrtout(unt, msg)
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

   call wrtout(unt," ")
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_merge_files(nfiles, v1files, dvdb_path, prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_merge_files'
!End of the abilint section

 implicit none

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
 ABI_DT_MALLOC(hdr1_list, (nfiles))

 nperts = size(hdr1_list)

 ! Write dvdb file (we only support fortran binary format)
 if (open_file(dvdb_path, msg, newunit=ount, form="unformatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 write(ount, err=10, iomsg=msg) dvdb_last_version
 write(ount, err=10, iomsg=msg) nperts

 ! Validate headers.
 ! TODO: Should perform consistency check on the headers, rearrange them in blocks of q-points.
 ! ignore POT1 files that do not correspond to atomic perturbations.

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
   if (prtvol > 0) call hdr_echo(hdr1_list(ii), fform, 3, unit=std_out)
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
   call hdr_fort_write(hdr1, ount, fform_pot, ierr)
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
   call hdr_free(hdr1_list(ii))
 end do
 ABI_DT_FREE(hdr1_list)

 write(std_out,"(a,i0,a)")"Merged successfully ",nfiles," files"

 ! List available perturbations.
 call dvdb_init(dvdb, dvdb_path, xmpi_comm_self)
 call dvdb_print(dvdb)
 call dvdb_list_perts(dvdb, [-1,-1,-1])
 call dvdb_free(dvdb)

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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine calc_eiqr(qpt,nrpt,rpt,eiqr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_eiqr'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_check_fform'
!End of the abilint section

 implicit none

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
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_test_v1rsym(db_path, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_test_v1rsym'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: db_path
 integer,intent(in) :: comm

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

 call dvdb_init(db, db_path, comm)
 db%debug = .True.
 !db%symv1 = .True.
 db%symv1 = .False.
 call dvdb_print(db)

 call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
 nfft = product(ngfft(1:3))
 call dvdb_open_read(db, ngfft, comm)

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

       ! TODO: (3) or (18)
       cplex = db%cplex_v1(v1pos)
       ngfft(1:3) = db%ngfft3_v1(:, v1pos)
       nfft = product(ngfft(:3))
       ABI_MALLOC(v1scf, (cplex*nfft, db%nspden))

       if (dvdb_read_onev1(db, idir, ipert, iqpt, cplex, nfft, ngfft, v1scf, msg) /= 0) then
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
       write(std_out,"(3(a,i0),a,es16.8)")"For iqpt= ",iqpt,", idir= ",idir,", ipert= ",ipert,", max_err= ",max_err

       ABI_FREE(irottb)
       ABI_FREE(v1scf)
     end do
   end do

 end do ! iqpt

 ABI_FREE(symafm1)
 ABI_FREE(symrel1)
 ABI_FREE(tnons1)

 call dvdb_free(db)

end subroutine dvdb_test_v1rsym
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_test_v1complete
!! NAME
!!  dvdb_test_v1complete
!!
!! FUNCTION
!!  Debugging tool used to test the symmetrization of the DFPT potentials.
!!
!! INPUTS
!!  db_path=Filename of the DVDB file.
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_test_v1complete(db_path, dump_path, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_test_v1complete'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: db_path,dump_path
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iqpt,pcase,idir,ipert,cplex,nfft,ispden,timerev_q,ifft,unt,my_rank
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
 type(dvdb_t),target :: db
!arrays
 integer :: ngfft(18),rfdir(3)
 integer,allocatable :: pflag(:,:)
 real(dp) :: qpt(3)
 integer,allocatable :: pertsy(:,:),rfpert(:),symq(:,:,:)
 real(dp),allocatable :: file_v1scf(:,:,:,:),symm_v1scf(:,:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call dvdb_init(db, db_path, comm)
 db%debug = .True.
 db%symv1 = .False. !; db%symv1 = .True.
 call dvdb_print(db)
 call dvdb_list_perts(db, [-1,-1,-1])

 call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
 nfft = product(ngfft(1:3))
 call dvdb_open_read(db, ngfft, comm)

 cryst => db%cryst
 call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
 nfft = product(ngfft(:3))

 ABI_MALLOC(pflag, (3, db%natom))

 ! Initialize the list of perturbations rfpert and rdfir
 ! WARNING: Only phonon perturbations are considered for the time being.
 ABI_MALLOC(rfpert,(db%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1

 ABI_MALLOC(symq, (4,2,cryst%nsym))
 ABI_MALLOC(pertsy, (3,db%mpert))

 unt = -1
 if (len_trim(dump_path) /= 0 .and. my_rank == master) then
   if (open_file(dump_path, msg, newunit=unt, action="write", status="unknown", form="formatted") /= 0) then
     MSG_ERROR(msg)
   end if
   write(std_out,"(a)")sjoin("Will write potentials to:", dump_path)
 end if

 do iqpt=1,db%nqpt
   qpt = db%qpts(:,iqpt)

   ! Examine the symmetries of the q wavevector
   call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=db%prtvol)

   ! Determine the symmetrical perturbations. Meaning of pertsy:
   !    0 for non-target perturbations
   !    1 for basis perturbations
   !   -1 for perturbations that can be found from basis perturbations
   call irreducible_set_pert(cryst%indsym,db%mpert,cryst%natom,cryst%nsym,&
     pertsy,rfdir,rfpert,symq,cryst%symrec,cryst%symrel)

   ! Read all potentials (here I assume that all perturbations are available)
   call dvdb_readsym_allv1(db, iqpt, cplex, nfft, ngfft, file_v1scf, db%comm)

   ! Copy basis perturbations in symm_v1scf and set pflag
   ABI_MALLOC(symm_v1scf, (cplex,nfft,db%nspden,db%natom3))
   symm_v1scf = huge(one); pflag = 0
   do pcase=1,3*db%cryst%natom
     idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
     if (pertsy(idir, ipert) == 1) then
       symm_v1scf(:,:,:,pcase) = file_v1scf(:,:,:,pcase)
       pflag(idir,ipert) = 1
     end if
   end do

   ! Complete potentials
   call v1phq_complete(cryst,qpt,ngfft,cplex,nfft,db%nspden,db%nsppol,db%mpi_enreg,db%symv1,pflag,symm_v1scf)

   ! Compare values.
   do pcase=1,3*cryst%natom
     idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
     if (pflag(idir,ipert) == 1) cycle

     write(std_out,"(2(a,i0),2a)")"For idir: ",idir, ", ipert: ", ipert, ", qpt: ",trim(ktoa(qpt))
     do ispden=1,db%nspden
       !write(std_out,"(a,es10.3)")"  max(abs(f1-f2))", &
       ! maxval(abs(file_v1scf(:,:,ispden,pcase) - symm_v1scf(:,:,ispden,pcase)))
       call vdiff_print(vdiff_eval(cplex,nfft,file_v1scf(:,:,ispden,pcase),symm_v1scf(:,:,ispden,pcase),cryst%ucvol))

       ! Debug: Write potentials to file.
       if (unt /= -1) then
         write(unt,*)"# q-point:", trim(ktoa(qpt))
         write(unt,*)"# idir: ",idir,", ipert: ",ipert,", ispden:", ispden
         write(unt,*)"# file_v1scf, symmetrized_v1scf, diff"
         if (cplex == 1) then
           do ifft=1,nfft
             write(unt,"(3(es12.4,2x))") &
               file_v1scf(1,ifft,ispden,pcase), symm_v1scf(1,ifft,ispden,pcase), &
               file_v1scf(1,ifft,ispden,pcase) - symm_v1scf(1,ifft,ispden,pcase)
           end do
         else
           do ifft=1,nfft
             write(unt, "(6(es12.4,2x))") &
               file_v1scf(1,ifft,ispden,pcase), symm_v1scf(1,ifft,ispden,pcase),   &
               file_v1scf(1,ifft,ispden,pcase) - symm_v1scf(1,ifft,ispden,pcase), &
               file_v1scf(2,ifft,ispden,pcase), symm_v1scf(2,ifft,ispden,pcase),   &
               file_v1scf(2,ifft,ispden,pcase) - symm_v1scf(2,ifft,ispden,pcase)
           end do
         end if
       end if

     end do
     write(std_out,*)""
   end do

   ABI_FREE(symm_v1scf)
   ABI_FREE(file_v1scf)
 end do

 ABI_FREE(pflag)
 ABI_FREE(rfpert)
 ABI_FREE(symq)
 ABI_FREE(pertsy)

 call dvdb_free(db)

 if (unt /= -1) close(unt)

end subroutine dvdb_test_v1complete
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
!!  db_path=Filename
!!  ngqpt(3)=Divisions of the Q-mesh reported in the DVDB file
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_test_ftinterp(db_path, ngqpt, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_test_ftinterp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: db_path
 integer,intent(in) :: comm
 integer,intent(in) :: ngqpt(3)

!Local variables-------------------------------
!scalars
 integer :: nfft,iq,cplex,mu,ispden
 type(dvdb_t) :: db
!arrays
 integer :: ngfft(18)
 real(dp),allocatable :: file_v1r(:,:,:,:),intp_v1r(:,:,:,:),tmp_v1r(:,:,:,:)

! *************************************************************************

 call dvdb_init(db, db_path, comm)
 db%debug = .True.
 db%symv1 = .True.
 db%symv1 = .False.
 call dvdb_print(db)

 ! Define FFT mesh for real space representation.
 call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
 nfft = product(ngfft(1:3))
 call dvdb_open_read(db, ngfft, xmpi_comm_self)

 ABI_MALLOC(intp_v1r, (2,nfft,db%nspden,db%natom3))
 ABI_MALLOC(file_v1r, (2,nfft,db%nspden,db%natom3))

 call dvdb_ftinterp_setup(db,ngqpt,1,[zero,zero,zero],nfft,ngfft,comm)

 do iq=1,db%nqpt
   ! Read data from file
   call dvdb_readsym_allv1(db, dvdb_findq(db, db%qpts(:,iq)), cplex, nfft, ngfft, tmp_v1r, comm)
   if (cplex == 1) then
     file_v1r(1,:,:,:) = tmp_v1r(1,:,:,:)
     file_v1r(2,:,:,:) = zero
   else
     file_v1r = tmp_v1r
   end if
   ABI_FREE(tmp_v1r)

   ! Interpolate data on the same q-point
   call dvdb_ftinterp_qpt(db, db%qpts(:,iq), nfft, ngfft, intp_v1r, comm)

   write(std_out,"(a)")sjoin("=== For q-point:", ktoa(db%qpts(:,iq)), "===")
   do mu=1,db%natom3
     do ispden=1,db%nspden
       write(std_out,"(2(a,i0))")"  iatom3: ", mu, ", ispden: ", ispden
       call vdiff_print(vdiff_eval(2,nfft,file_v1r(:,:,ispden,mu),intp_v1r(:,:,ispden,mu),db%cryst%ucvol))

       !do ifft=1,nfft
       !  write(std_out,*)file_v1r(1,ifft,ispden,mu),intp_v1r(1,ifft,ispden,mu),&
       !  file_v1r(2,ifft,ispden,mu),intp_v1r(2,ifft,ispden,mu)
       !end do
       write(std_out,*)" "
     end do
   end do
   write(std_out,*)" "
 end do ! iq

 ABI_FREE(intp_v1r)
 ABI_FREE(file_v1r)

 call dvdb_free(db)

end subroutine dvdb_test_ftinterp
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_v1r_long_range
!! NAME
!!  dvdb_v1r_long_range
!!
!! FUNCTION
!!  Compute the long-range part of the phonon potential
!!  due to the Born effective charges, PRL 115, 176401 (2015) [[cite:Verdi2015]].
!!
!!    V^L_{iatom,idir}(r) = i (4pi/vol) sum_G (q+G) . Zeff_{iatom,idir}
!!                           e^{i (q + G) . (r - tau_{iatom})} / ((q + G) . dielt . (q + G))
!!
!!  where Zeff and dielt are the Born effective charge tensor and the dielectric tensor,
!!  tau is the atom position, and vol is the volume of the unit cell.
!!
!! INPUTS
!!  db = the DVDB object.
!!  qpt = the q-point in reduced coordinates.
!!  iatom = atom index.
!!  idir = direction index.
!!  nfft = number of fft points.
!!  ngfft(18) = FFT mesh.
!!
!! OUTPUT
!!  v1r_lr = dipole potential
!!
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_v1r_long_range(db,qpt,iatom,idir,nfft,ngfft,v1r_lr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_v1r_long_range'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dvdb_t),intent(in) :: db
 integer,intent(in) :: iatom, idir
 integer,intent(in) :: nfft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: v1r_lr(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: n1, n2, n3, nfftot, ig
 real(dp) :: fac, qGZ, denom, denom_inv, qtau
 real(dp) :: re, im, phre, phim
 real(dp) :: dummy
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer, allocatable :: gfft(:,:)
 real(dp) :: gprimd(3,3), rprimd(3,3)
 real(dp) :: dielt(3,3), zeff(3,3), dielt_red(3,3)
 real(dp) :: qG(3), Zstar(3), tau(3)
 real(dp), allocatable :: v1G_lr(:,:)


! *************************************************************************

 ! Make sure FFT parallelism is not used
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")

 ! Fake MPI_type for sequential fft execution
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')

 ! Allocate memory
 ABI_MALLOC(gfft, (3, nfft))
 ABI_MALLOC(v1G_lr, (2,nfft))

 ! Reciprocal and real space primitive vectors
 gprimd = db%cryst%gprimd
 rprimd = db%cryst%rprimd

 !Prefactor
 fac = four_pi / db%cryst%ucvol

 ! Transform the Born effective charge tensor from Cartesian to reduced coordinates
 ! and select the relevant direction.
 zeff = db%zeff(:,:,iatom)
 Zstar = matmul(transpose(gprimd), matmul(zeff, rprimd(:,idir))) * two_pi

 ! Transform the dielectric tensor from Cartesian to reduced coordinates.
 dielt = db%dielt
 dielt_red = matmul(transpose(gprimd), matmul(dielt, gprimd)) * two_pi ** 2

 ! Atom position
 tau = db%cryst%xred(:,iatom)

 ! Get the set of G vectors
 call get_gftt(ngfft,qpt,db%cryst%gmet,dummy,gfft)

 ! Compute the long-range potential in G-space
 v1G_lr = zero
 do ig=1,nfft

   ! (q + G)
   qG(:) = qpt(:) + gfft(:,ig)

   ! (q + G) . Zeff(:,idir,iatom)
   qGZ = dot_product(qG, Zstar)

   ! (q + G) . dielt . (q + G)
   denom = dot_product(qG, matmul(dielt_red, qG))

   ! Avoid (q+G) = 0
   denom_inv = denom / (denom ** 2 + tol10 ** 2)

   ! Phase factor exp(-i (q+G) . tau)
   qtau = - two_pi * dot_product(qG, tau)
   phre = cos(qtau); phim = sin(qtau)

   re = zero
   im = fac * qGZ * denom_inv

   v1G_lr(1,ig) = phre * re - phim * im
   v1G_lr(2,ig) = phim * re + phre * im

 end do

 ! Free memory
 ABI_FREE(gfft)

 ! FFT to get the long-range potential in r space
 v1r_lr = zero
 call fourdp(2,v1G_lr,v1r_lr,1,MPI_enreg_seq,nfft,ngfft,1,0)

 ! Multiply by  exp(i q . r)
 call times_eikr(qpt,ngfft,nfft,1,v1r_lr)

 ! Free memory
 ABI_FREE(v1G_lr)

 call destroy_mpi_enreg(MPI_enreg_seq)

end subroutine dvdb_v1r_long_range
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
!!      cwtime,dvdb_get_v1scf_qpt,dvdb_get_v1scf_rpt,dvdb_open_read
!!      hdr_fort_read,hdr_fort_write,hdr_free,irreducible_set_pert
!!      kpts_ibz_from_kptrlatt,littlegroup_q,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dvdb_interpolate_and_write(dtfil, ngfft, ngfftf, cryst, dvdb, &
&          ngqpt_coarse, nqshift_coarse, qshift_coarse, &
&          ngqpt, qptopt, mpi_enreg, comm)

 use m_time,           only : cwtime

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_interpolate_and_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 integer,intent(in) :: qptopt,nqshift_coarse
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(dvdb_t),target,intent(inout) :: dvdb
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18), ngfftf(18)
 integer,intent(in) :: ngqpt(3), ngqpt_coarse(3)
 real(dp),intent(in) :: qshift_coarse(3,nqshift_coarse)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: fform_pot=111
 integer :: ierr
 integer :: my_rank,nproc,idir,ipert,iat,ipc,ispden
 integer :: cplex,db_iqpt,natom,natom3,npc,trev_q,nspden
 integer :: nqbz, nqibz, iq, ifft, nqbz_coarse
 integer :: nperts_read, nperts_interpolate, nperts
 integer :: nqpt_read, nqpt_interpolate
 integer :: nfft,nfftf
 integer :: ount, unt, fform
 real(dp) :: cpu,wall,gflops
 logical :: i_am_master
 character(len=500) :: msg
 character(len=fnlen) :: new_ddb_fname
!arrays
 integer :: qptrlatt(3,3)
 integer :: symq(4,2,cryst%nsym)
 integer :: rfdir(3)
 integer,allocatable :: pinfo(:,:),rfpert(:),pertsy(:,:,:),iq_read(:)
 real(dp) :: qpt(3)
 real(dp) :: rhog1_g0(2)
 real(dp),allocatable :: v1scf(:,:,:), v1scf_rpt(:,:,:,:),v1(:)
 real(dp),allocatable :: wtq(:),qibz(:,:),qbz(:,:),q_interp(:,:),q_read(:,:)
 type(hdr_type) :: hdr_ref

!************************************************************************
 ABI_UNUSED(mpi_enreg%nproc)

 call cwtime(cpu,wall,gflops,"start")

 write(msg, '(2a)') "Interpolation of the electron-phonon coupling potential", ch10
 call wrtout(ab_out, msg, "COLL", do_flush=.True.)
 call wrtout(std_out, msg, "COLL", do_flush=.True.)

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 i_am_master = (my_rank == master)

 ! ======================= !
 ! Setup fine q-point grid
 ! ======================= !

 nfftf = product(ngfftf(1:3))
 nfft = product(ngfft(1:3))

 ! Generate the list of irreducible q-points in the grid
 qptrlatt = zero
 qptrlatt(1,1) = ngqpt(1); qptrlatt(2,2) = ngqpt(2); qptrlatt(3,3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt, 1, &
 &                           [zero, zero, zero], nqibz, qibz, wtq, nqbz, qbz)

 nqbz_coarse = product(ngqpt_coarse) * nqshift_coarse

 ! ========================================== !
 ! Prepare the header to write the potentials
 ! ========================================== !

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
 if (allocated(hdr_ref%symrel)) ABI_FREE(hdr_ref%symrel)
 if (allocated(hdr_ref%tnons)) ABI_FREE(hdr_ref%tnons)
 if (allocated(hdr_ref%symafm)) ABI_FREE(hdr_ref%symafm)
 hdr_ref%nsym = cryst%nsym
 ABI_MALLOC(hdr_ref%symrel, (3,3,hdr_ref%nsym))
 ABI_MALLOC(hdr_ref%tnons, (3,hdr_ref%nsym))
 ABI_MALLOC(hdr_ref%symafm, (hdr_ref%nsym))
 hdr_ref%symrel(:,:,:) = cryst%symrel(:,:,:)
 hdr_ref%tnons(:,:) = cryst%tnons(:,:)
 hdr_ref%symafm(:) = cryst%symafm(:)

 hdr_ref%ngfft = ngfftf(1:3)


 ! ======================================= !
 ! Open DVDB and copy important dimensions
 ! ======================================= !

 call dvdb_open_read(dvdb, ngfftf, xmpi_comm_self)

 cplex = 2
 natom = cryst%natom
 natom3 = 3 * natom
 nspden = dvdb%nspden

 ! ================================================== !
 ! Sort the q-points to read and those to interpolate
 ! and find the irreducible perturbations
 ! ================================================== !
 ABI_MALLOC(iq_read, (nqibz))
 ABI_MALLOC(q_read, (3,nqibz))
 ABI_MALLOC(q_interp, (3,nqibz))

 ABI_MALLOC(pertsy, (nqibz,3,dvdb%mpert))
 ABI_MALLOC(rfpert, (dvdb%mpert))
 ABI_MALLOC(pinfo, (3,3*dvdb%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1

 pertsy = zero

 nqpt_read = 0
 nperts_read = 0
 nqpt_interpolate = 0
 nperts_interpolate = 0

 do iq=1,nqibz

   qpt = qibz(:,iq)

   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb_findq(dvdb, qpt)

   if (db_iqpt /= -1) then

     call wrtout(std_out, sjoin("Q-point: ",ktoa(qpt)," found in DVDB with index ",itoa(db_iqpt)))
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

     call wrtout(std_out, sjoin("Q-point: ",ktoa(qpt), "not found in DVDB. Will interpolate."))
     nqpt_interpolate = nqpt_interpolate + 1
     q_interp(:,nqpt_interpolate) = qpt(:)


     ! Examine the symmetries of the q wavevector
     call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,trev_q,prtvol=0)

     ! Find the list of irreducible perturbations for this q-point.
     call irreducible_set_pert(cryst%indsym,dvdb%mpert,cryst%natom,cryst%nsym,&
     &    pertsy(nqpt_interpolate,:,:),rfdir,rfpert,symq,cryst%symrec,cryst%symrel)

     do iat=1,natom
       do idir=1,3
         ipert = (iat-1) * 3 + idir
         if (pertsy(nqpt_interpolate,idir,iat) == 1) then
           nperts_interpolate = nperts_interpolate + 1
         end if

       end do
     end do

   end if
 end do


 ! ================================================= !
 ! Open the new DVDB file and write preliminary info
 ! ================================================= !
 nperts = nperts_read + nperts_interpolate

 if (my_rank == master) then
   new_ddb_fname = strcat(dtfil%filnam_ds(4), '_DVDB')
   if (open_file(new_ddb_fname, msg, newunit=ount, form="unformatted", action="write", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if
   write(ount, err=10, iomsg=msg) dvdb_last_version
   write(ount, err=10, iomsg=msg) nperts
 end if


 ! ================================================================= !
 ! Master reads all available perturbations and copy in the new DVDB
 ! ================================================================= !

 ABI_MALLOC(v1scf, (cplex,nfftf,nspden))
 ABI_MALLOC(v1, (cplex*nfftf))

 rhog1_g0 = zero

 if (my_rank == master) then
   do iq=1,nqpt_read

     qpt = q_read(:,iq)
     db_iqpt = iq_read(iq)

     ! Read each irreducible perturbation potentials
     npc = dvdb_get_pinfo(dvdb, db_iqpt, cplex, pinfo)
     ABI_CHECK(npc /= 0, "npc == 0!")

     do ipc=1,npc
       idir = pinfo(1,ipc); iat = pinfo(2,ipc); ipert = pinfo(3, ipc)
       if (dvdb_read_onev1(dvdb, idir, iat, db_iqpt, cplex, nfftf, ngfftf, v1scf, msg) /= 0) then
         MSG_ERROR(msg)
       end if

       hdr_ref%qptn = qpt
       hdr_ref%pertcase = ipert

       ! Write header
       call hdr_fort_write(hdr_ref, ount, fform_pot, ierr)
       ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

       do ispden=1,nspden
         v1 = reshape(v1scf(:,:,ispden), (/cplex*nfftf/))
         write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfftf)
       end do

       if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0

     end do

   end do
 end if

 call xmpi_barrier(comm)


 ! ================================================================ !
 ! Interpolate the potential for q-points not in the original DVDB
 ! ================================================================ !

 dvdb%nrpt = nqbz_coarse
 ABI_MALLOC(v1scf_rpt, (2,dvdb%nrpt,nfftf,dvdb%nspden))

 do iat=1,natom
   do idir=1,3
     ipert = (iat-1) * 3 + idir

     if (sum(pertsy(:,idir,iat)) == -nqpt_interpolate) cycle

     call wrtout(std_out, sjoin("Interpolating perturbation iat,idir = ",itoa(iat), itoa(idir)))

     ! Compute phonon potential in real space lattice representation
     call dvdb_get_v1scf_rpt(dvdb, cryst, ngqpt_coarse, nqshift_coarse, &
     &                       qshift_coarse, nfftf, ngfftf, &
     &                       dvdb%nrpt, dvdb%nspden, ipert, v1scf_rpt, comm)


     do iq=1,nqpt_interpolate

       if (pertsy(iq,idir,iat) == -1) cycle

       qpt = q_interp(:,iq)

       ! Interpoalte the phonon potential
       call dvdb_get_v1scf_qpt(dvdb, cryst, qpt, nfftf, ngfftf, dvdb%nrpt, &
       &                       dvdb%nspden, ipert, v1scf_rpt, v1scf, comm)

       ! Master writes the file
       if (my_rank == master) then

         hdr_ref%qptn = qpt
         hdr_ref%pertcase = ipert
         ! Write header
         call hdr_fort_write(hdr_ref, ount, fform_pot, ierr)
         ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

         do ispden=1,nspden
           v1 = reshape(v1scf(:,:,ispden), (/cplex*nfftf/))
           write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfftf)
         end do

         if (dvdb_last_version > 1) write(ount, err=10, iomsg=msg) rhog1_g0

       end if

     end do

     ABI_FREE(dvdb%rpt)

   end do
 end do

 if (my_rank == master) close(ount)

 ! ========================================================================== !
 ! Free memory
 ! ========================================================================== !

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
 ABI_FREE(rfpert)
 ABI_FREE(pinfo)

 call hdr_free(hdr_ref)

 call cwtime(cpu,wall,gflops,"stop")

 write(msg, '(2a)') "Interpolation of the electron-phonon coupling potential completed", ch10
 call wrtout(ab_out, msg, "COLL", do_flush=.True.)
 call wrtout(std_out, msg, "COLL", do_flush=.True.)


 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(msg)

end subroutine dvdb_interpolate_and_write
!!***

END MODULE m_dvdb
