!!****m* ABINIT/m_wfk
!! NAME
!!  m_wfk
!!
!! FUNCTION
!!  This module defines the wfk_t object providing a high-level API
!!  to perform common IO operations on the WFK file produced by ABINIT.
!!  The API wraps thee different formats/io-libraries:
!!
!!    1) binary Fortran files with sequential Fortran IO (read, write)
!!    2) binary Fortran files with MPI-IO primitives (C-steam + Fortran records)
!!    3) Netcdf files with parallel IO a.k.a HDF5
!!
!!  and emulate random access when binary Fortran files are used in *read-only* mode.
!!  See notes below for more info.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! NOTES
!!  1) The wfk_t object supports random access also when plain Fortran-IO is used.
!!     One can easily *read* the block of wavefunctions with a given (kpt,spin)
!!     by simply passing the appropriate indices (ik_ibz,spin) to the wfk_read_ routines.
!!     Note however that the same feature is not available in write mode when Fortran IO
!!     is used. In this case indeed one should access the block of wavefunctions
!!     according to their (kpt,spin) indices in order to write the correct record markers
!!     MPI-IO and NETCDF do not have such limitation.
!!
!!  2) MPI-IO read operations are done with file views even for contiguous data of the same type.
!!     I found, indeed, that mixing views with explicit offset calls causes
!!     wrong results unless the file is closed and re-open! Very strange since, according to
!!     the documentation, the two APIs do not interfere and can be mixed. Calls to MPI_FILE_SEEK
!!     to reset the pointer to the start of the file do not solve the problem. Don't know if it's
!!     a feature or a bug (the problem showed up with MPICH2, I haven't tested other MPI libraries)
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_wfk

 use defs_basis
 use m_abicore
 use m_build_info
 use m_errors
 use m_dtset
#ifdef HAVE_MPI2
 use mpi
#endif
 use m_xmpi
 use m_mpiotk
 use m_hdr
 use m_sort
 use m_crystal
 use m_pawtab
 use m_ebands
 use m_pawrhoij
 use m_wffile
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_clib

 use defs_abitypes,  only : MPI_type
 use defs_datatypes, only : pseudopotential_type, ebands_t
 use defs_wvltypes,  only : wvl_internal_type
 use m_geometry,     only : metric
 use m_time,         only : cwtime, cwtime_report, asctime
 use m_fstrings,     only : sjoin, strcat, endswith, itoa, ktoa
 use m_io_tools,     only : get_unit, mvrecord, iomode_from_fname, iomode2str, open_file, close_unit, delete_file, file_exists
 use m_numeric_tools,only : mask2blocks
 use m_cgtk,         only : cgtk_rotate
 use m_fftcore,      only : get_kg, ngfft_seq
 use m_distribfft,   only : init_distribfft_seq
 use m_mpinfo,       only : destroy_mpi_enreg, initmpi_seq
 use m_rwwf,         only : rwwf

 implicit none

 private

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

 integer,private,parameter :: WFK_NOMODE    = 0
 integer,private,parameter :: WFK_READMODE  = 1
 integer,private,parameter :: WFK_WRITEMODE = 2
!!***

!----------------------------------------------------------------------

!!****t* m_wfk/wfk_t
!! NAME
!!  wfk_t
!!
!! FUNCTION
!!  File handler for the WFK file.
!!
!! SOURCE

 type,public :: wfk_t

  integer :: fh
   !  unit number if IO_MODE_FORTRAN
   !  MPI file handler if IO_MODE_MPI
   !  Netcdf file handler if IO_MODE_ETSF

  integer :: iomode
   ! Method used to access the WFK file:
   !   IO_MODE_FORTRAN for usual Fortran IO routines
   !   IO_MODE_MPI if MPI/IO routines.
   !   IO_MODE_ETSF, NetCDF format read/written via etsf-io.

  integer :: mband
  ! Max number of bands stored on file (MAX(Hdr%nband))

  integer :: nkpt
  ! Number of k-points.

  integer :: nsppol
  ! Number of spins

  integer :: nspinor
  ! Number of spinor components.

  integer :: formeig
   ! formeig=format of the eigenvalues
   !    0 => vector of eigenvalues
   !    1 => hermitian matrix of eigenvalues
   ! TODO: this should be reported somewhere in the WFK file, at present is passed to wfk_open

  integer :: fform
   ! File type format of the header

  integer :: rw_mode = WFK_NOMODE
   ! (Read|Write) mode

  character(len=fnlen) :: fname = ABI_NOFILE
   ! File name

  integer :: master
   ! master node of the IO procedure

  integer :: my_rank
   ! index of my processor in the MPI communicator comm

  integer :: nproc
   ! number of processors in comm

  integer :: comm
   ! MPI communicator

  integer :: recn_eof
   ! EOF record number (used for Fortran IO)

  integer(XMPI_OFFSET_KIND) :: offset_eof
  ! EOF offset (used for MPI-IO access)

  logical :: debug=.FALSE.
  !logical :: debug=.TRUE.

  type(hdr_type) :: Hdr
   ! Abinit header.

  integer,allocatable :: nband(:,:)
  ! nband(nkpt,nsppol) = Number of bands at each (k,s)

  integer :: f90_fptr(3) = [0,0,0]
  ! The position of the file pointer used for sequential access with Fortran-IO.
  !  f90_fprt(1) = Index of the k-point associated to the block.
  !  f90_fprt(2) = the spin associated to the block.
  !  f90_fprt(3) = Record Type (see REC_* variables).
  !  [0,0,0] corresponds to the beginning of the file.
  !  FPTR_EOF signals the end of file

  integer,allocatable :: recn_ks(:,:,:)
   ! recn_ks(k,s,1) : record number of  (npw, nspinor, nband_disk)
   ! recn_ks(k,s,2) : record number of the (k+G) vectors.
   ! recn_ks(k,s,3) : record number of the eigenvalues.
   ! recn_ks(k,s,4) : record number of the first wavefunction in the wf coefficients block.

  integer(XMPI_OFFSET_KIND),allocatable :: offset_ks(:,:,:)
   ! offset_ks(k,s,1) : offset of the record: npw, nspinor, nband_disk.
   ! offset_ks(k,s,2) : offset of the Second record: (k+G) vectors.
   ! offset_ks(k,s,3) : offset of the third record eigenvalues.
   ! offset_ks(k,s,4) : offset of the fourth record (wavefunction coefficients).
   !
   ! **********************************************************************
   ! NB: The offset point to the Fortran record marker and not to the data
   ! **********************************************************************

  integer(XMPI_OFFSET_KIND) :: hdr_offset
   ! offset of the header
   ! TODO this should be the output of a hdr method!

  integer(XMPI_OFFSET_KIND) :: chunk_bsize
   ! IO is performed in chunks of max size chunk_bsize [bytes]

  contains

    procedure :: open_write => wfk_open_write
     ! Open the WFK file in write mode.

    procedure :: close => wfk_close
      ! Close the WFK file and release the memory allocated in wfk_t.

    procedure :: print => wfk_print
      ! Print info on the wfk_t object

    procedure :: findk => wfk_findk
      ! Returns the index of the k-point in the WFK file.

    procedure :: compare => wfk_compare
      ! Test two wfk_t objects for consistency.

    procedure :: read_band_block => wfk_read_band_block
      ! Read a contiguous block of bands for a given (kpoint, spin)

    procedure :: read_bks => wfk_read_bks
      ! Read the wavefunction and the eigenvalues for a given (band, k-point, spin)

    procedure :: write_band_block  => wfk_write_band_block
      ! Write a contiguous block of bands for a given (kpoint, spin)

    procedure :: read_bmask => wfk_read_bmask
      ! Read a scattered set of bands for a given (kpoint, spin).

    procedure :: read_eigk => wfk_read_eigk
      ! Read the eigenvalues at a given (kpoint,spin).

    procedure :: write_h1mat => wfk_write_h1mat
      ! Write all the H1 matrix elements.
 end type wfk_t


 public :: wfk_open_read           ! Open the WFK file in read mode.
 public :: wfk_tofullbz            ! Generate a new WFK file with wavefunctions in the full BZ and istwfk==1
                                   ! Mainly used to interface ABINIT with other codes that
                                   ! cannot handle symmetries e.g. lobster
 public :: wfk_nc2fort             ! Convert a netcdf WFK file to a Fortran WFK file.
 public :: wfk_klist2mesh
 public :: wfk_ncdef_dims_vars     ! Define basic dimensions for netcdf file format.
 public :: wfk_read_ebands         ! Read the GS eigenvalues and return ebands_t object.
 public :: wfk_read_eigenvalues    ! Read all the GS eigenvalues stored in the WFK file.
 public :: wfk_read_h1mat          ! Read all the H1 matrix elements.

 ! Profiling tools
 public :: wfk_prof                ! Profiling tool.

 ! Unit tests
 public :: wfk_diff                ! Compare two WFK file for binary equality.
 public :: wfk_create_wfkfile      ! Create a FAKE WFK file.
 public :: wfk_check_wfkfile       ! Read a FAKE WFK file and perform basic tests.

!!***

! Indices associated to the start of the different records of the WFK file.
 integer,private,parameter :: REC_HDR=0
 integer,private,parameter :: REC_NPW=1
 integer,private,parameter :: REC_KG =2
 integer,private,parameter :: REC_EIG=3
 integer,private,parameter :: REC_CG =4
 integer,private,parameter :: REC_NUM=REC_CG

 integer,private,parameter :: FPTR_EOF(3) = [-1,-1,-1]

 integer(XMPI_OFFSET_KIND),private,parameter :: WFK_CHUNK_BSIZE = 2000 * (1024.0_dp**2)
   ! Maximum size (in bytes) of the block of wavefunctions that are (read|written)
   ! in a single MPI-IO call. (Some MPI-IO implementation crashes if we try to
   ! (read|write) a big chunk of data with a single call.

!----------------------------------------------------------------------

!!****t* m_wfk/kvars_t
!! NAME
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: kvars_t
   integer,allocatable  :: kg_k(:,:)
   real(dp),pointer :: occ_k(:)   => null()
   real(dp),pointer :: eig_k(:)   => null()
 end type kvars_t

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_open_read
!! NAME
!!  wfk_open_read
!!
!! FUNCTION
!!  Open the WFK file in read mode.
!!
!! INPUTS
!!  fname = Name of the file
!!  formeig = 0 for GS wavefunctions, 1 for RF wavefunctions.
!!  iomode = access mode
!!  funt = Fortran unit numer. Only used if iomode == IO_MODE_FORTRAN
!!  comm = MPI communicator (used for collective parallel IO)
!!
!! OUTPUT
!!  Wfk<class(wfk_t)> = WFK handler initialized and set in read mode
!!  [Hdr_out]=Copy of the abinit header
!!
!! PARENTS
!!      conducti_nc,d2frnl,dfpt_looppert,dfpt_nstdy,dfpt_nstpaw,fold2Bloch
!!      initwf,ioprof,m_cut3d,m_wfd,m_wfk,mrggkk,optic
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_open_read(Wfk,fname,formeig,iomode,funt,comm,Hdr_out)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,comm,formeig,funt
 character(len=*),intent(in) :: fname
 class(wfk_t),intent(inout) :: Wfk
 type(hdr_type),optional,intent(inout) :: Hdr_out  ! should be intent(out), but psc miscompiles the call!

!Local variables-------------------------------
!scalars
 integer :: ierr,mpierr
 character(len=500) :: msg
#ifdef HAVE_MPI_IO
 integer :: fform !,ncerr
#endif

!************************************************************************

 DBG_ENTER("COLL")

 !Initialize the mandatory data of the Wfk datastructure
 !@wfk_t
 Wfk%rw_mode     = WFK_READMODE
 Wfk%chunk_bsize = WFK_CHUNK_BSIZE

 Wfk%fname     = fname
 Wfk%formeig   = formeig
 Wfk%iomode    = iomode; if (endswith(fname, ".nc")) wfk%iomode = IO_MODE_ETSF
 ! This is to test the different versions.
 !wfk%iomode    = IO_MODE_MPI
 !if (.not. endswith(fname, ".nc") .and. xmpi_comm_size == 1) wfk%iomode == IO_MODE_FORTRAN

 Wfk%comm      = comm
 Wfk%master    = 0
 Wfk%my_rank   = xmpi_comm_rank(comm)
 Wfk%nproc     = xmpi_comm_size(comm)

 ! Reads fform and the Header.
 call hdr_read_from_fname(Wfk%Hdr,fname,Wfk%fform,comm)
 ABI_CHECK(Wfk%fform/=0,"fform ==0")

 if (Wfk%debug) call Wfk%Hdr%echo(Wfk%fform, 4, unit=std_out)

 ! Copy the header if required.
 if (present(Hdr_out)) call hdr_copy(Wfk%Hdr,Hdr_out)

 ! Useful dimensions
 Wfk%mband   = MAXVAL(Wfk%Hdr%nband)
 Wfk%nkpt    = Wfk%Hdr%nkpt
 Wfk%nsppol  = Wfk%Hdr%nsppol
 Wfk%nspinor = Wfk%Hdr%nspinor

 ABI_MALLOC(Wfk%nband, (Wfk%nkpt,Wfk%nsppol))
 Wfk%nband = RESHAPE(Wfk%Hdr%nband, (/Wfk%nkpt,Wfk%nsppol/))

 ierr=0
 select case (wfk%iomode)
 case (IO_MODE_FORTRAN)
   ! All processors see a local Fortran binary file.
   ! Each node opens the file, skip the header and set f90_fptr.
   Wfk%fh = funt
   if (open_file(Wfk%fname,msg,unit=Wfk%fh,form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if

   ! Precompute number of records for Fortran IO.
   call wfk_compute_offsets(Wfk)

   call hdr_skip(Wfk%fh,ierr)
   ABI_CHECK(ierr==0, "hdr_skip returned ierr! /= 0")
   Wfk%f90_fptr = [1,1,REC_NPW]

#ifdef HAVE_MPI_IO
 case (IO_MODE_MPI)
   call MPI_FILE_OPEN(Wfk%comm, Wfk%fname, MPI_MODE_RDONLY, xmpio_info, Wfk%fh, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")
   !call MPI_FILE_SET_VIEW(Wfk%fh,origin,MPI_BYTE,MPI_BYTE,'native',xmpio_info,mpierr)

   call hdr_mpio_skip(Wfk%fh,fform,Wfk%hdr_offset)

   ! Precompute offsets for MPI-IO access
   if (Wfk%hdr_offset > 0) then
     call wfk_compute_offsets(Wfk)
   else
     MSG_ERROR("hdr_offset <=0")
   end if
#endif

#ifdef HAVE_NETCDF
 case (IO_MODE_ETSF)
   NCF_CHECK(nctk_open_read(wfk%fh, wfk%fname, wfk%comm))
#endif

 case default
   MSG_ERROR(sjoin('Wrong or unsupported iomode:', itoa(wfk%iomode)))
 end select

 DBG_EXIT("COLL")

end subroutine wfk_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_open_write
!! NAME
!!  wfk_open_write
!!
!! FUNCTION
!!  Open the WFK file in write mode.
!!
!! INPUTS
!!  fname = Name of the file
!!  formeig = 0 for GS wavefunctions, 1 for RF wavefunctions.
!!  iomode = access mode
!!  funt = Fortran unit numer for  Only used if iomode == IO_MODE_FORTRAN
!!  comm = MPI communicator (used for MPI-IO)
!!  [write_hdr]=True if the header should be written (default)
!!  [write_frm]=True if the fortran record markers should be written (default). Only if Fortran binary file.
!!
!! OUTPUT
!!  Wfk<class(wfk_t)> = WFK handler initialized and set in read mode
!!
!! PARENTS
!!      m_iowf,m_wfd,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_open_write(Wfk,Hdr,fname,formeig,iomode,funt,comm,write_hdr,write_frm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,comm,formeig,funt
 character(len=*),intent(in) :: fname
 logical,optional,intent(in) :: write_hdr,write_frm
 type(hdr_type),intent(in) :: Hdr
 class(wfk_t),intent(out) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: mpierr,ierr
 real(dp) :: cpu,wall,gflops
 logical :: do_write_frm,do_write_hdr
 character(len=500) :: msg
#ifdef HAVE_MPI_IO
 integer :: fform,nfrec,sc_mode
 integer(XMPI_OFFSET_KIND) :: offset
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecords(:)
#endif
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif

!************************************************************************

 DBG_ENTER("COLL")

 do_write_hdr = .TRUE.; if (present(write_hdr)) do_write_hdr = write_hdr
 do_write_frm = .TRUE.; if (present(write_frm)) do_write_frm = write_frm

 !Initialize mandatory data of the Wfk datastructure
 !@wfk_t
 Wfk%rw_mode     = WFK_WRITEMODE
 Wfk%chunk_bsize = WFK_CHUNK_BSIZE

 Wfk%fname     = fname
 Wfk%formeig   = formeig
 Wfk%iomode    = iomode; if (endswith(fname, ".nc")) wfk%iomode = IO_MODE_ETSF
 ! This is to test the different versions.
 !wfk%iomode   = IO_MODE_MPI
 !if (.not. endswith(fname, ".nc") .and. xmpi_comm_size == 1) wfk%iomode == IO_MODE_FORTRAN

 Wfk%comm      = comm
 Wfk%master    = 0
 Wfk%my_rank   = xmpi_comm_rank(comm)
 Wfk%nproc     = xmpi_comm_size(comm)
 Wfk%fform     = 2

 ! Copy the header
 call hdr_copy(Hdr,Wfk%Hdr)

 ! Master writes fform and the Header (write it afterwards if IO_MODE_ETSF)
 if (Wfk%my_rank==Wfk%master .and. do_write_hdr .and. iomode /= IO_MODE_ETSF) then
   call Wfk%Hdr%write_to_fname(Wfk%fname, Wfk%fform)
   if (Wfk%debug) call Wfk%Hdr%echo(Wfk%fform, 4, unit=std_out)
 end if
 call xmpi_barrier(Wfk%comm)

 ! Useful dimensions
 Wfk%mband   = MAXVAL(Wfk%Hdr%nband)
 Wfk%nkpt    = Wfk%Hdr%nkpt
 Wfk%nsppol  = Wfk%Hdr%nsppol
 Wfk%nspinor = Wfk%Hdr%nspinor

 ABI_MALLOC(Wfk%nband, (Wfk%nkpt,Wfk%nsppol))
 Wfk%nband = RESHAPE(Wfk%Hdr%nband, (/Wfk%nkpt,Wfk%nsppol/))

 ierr=0

 select case (wfk%iomode)
 case (IO_MODE_FORTRAN)
   ABI_CHECK(wfk%nproc == 1, "Cannot use Fortran-IO to write WFK file with nprocs > 1")
   Wfk%fh = funt
   if (open_file(Wfk%fname,msg,unit=Wfk%fh,form="unformatted", status="unknown", action="readwrite") /= 0) then
     MSG_ERROR(msg)
   end if

   ! Precompute number of records for Fortran IO.
   call wfk_compute_offsets(Wfk)

   call hdr_skip(Wfk%fh,ierr)
   Wfk%f90_fptr = (/1,1,REC_NPW/)

#ifdef HAVE_MPI_IO
 case (IO_MODE_MPI)

   call cwtime(cpu, wall, gflops, "start")

   ! FIXME: mode flags should be rationalized
   !call MPI_FILE_OPEN(Wfk%comm, Wfk%fname, MPI_MODE_CREATE + MPI_MODE_WRONLY, xmpio_info, Wfk%fh, mpierr)
   !call MPI_FILE_OPEN(Wfk%comm, Wfk%fname, MPI_MODE_CREATE + MPI_MODE_RDWR, xmpio_info, Wfk%fh, mpierr)
   call MPI_FILE_OPEN(Wfk%comm, Wfk%fname,  MPI_MODE_RDWR, xmpio_info, Wfk%fh, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")

   !call MPI_FILE_SET_VIEW(Wfk%fh,origin,MPI_BYTE,MPI_BYTE,'native',xmpio_info,mpierr)

   ! TODO
   !%% call MPI_File_set_size(Wfk%fh, MPI_Offset size, mpierr)
   !ABI_CHECK_MPI(mpierr,"MPI_FILE_SET_SIZE")

   call hdr_mpio_skip(Wfk%fh,fform,Wfk%hdr_offset)
   ABI_CHECK(fform == Wfk%fform,"fform != Wfk%fform")
   !call wfk%Hdr%echo(wfk%fform, 4, unit=std_out)

   ! Precompute offsets for MPI-IO access
   if (Wfk%hdr_offset > 0) then
     call wfk_compute_offsets(Wfk)
   else
     MSG_ERROR("hdr_offset <=0")
   end if
   call cwtime_report(" FILE_OPEN", cpu, wall, gflops)

   ! Write Fortran record markers.
   if (do_write_frm) then
     call cwtime(cpu, wall, gflops, "start")
     call hdr_bsize_frecords(Wfk%Hdr,Wfk%formeig,nfrec,bsize_frecords)

     sc_mode = xmpio_collective
     offset = Wfk%hdr_offset

     if (sc_mode == xmpio_collective) then
       call xmpio_write_frmarkers(Wfk%fh,offset,sc_mode,nfrec,bsize_frecords,ierr)
     else
       ierr = 0
       if (Wfk%my_rank == Wfk%master) then
         call xmpio_write_frmarkers(Wfk%fh,offset,xmpio_single,nfrec,bsize_frecords,ierr)
       end if
     end if
     ABI_CHECK(ierr == 0, "xmpio_write_frmarkers returned ierr!=0")

     !call MPI_FILE_SYNC(Wfk%fh,mpierr)
     !ABI_CHECK_MPI(mpierr,"FILE_SYNC")

     if (Wfk%debug) then
       call xmpio_check_frmarkers(Wfk%fh,offset,sc_mode,nfrec,bsize_frecords,ierr)
       ABI_CHECK(ierr == 0, "xmpio_check_frmarkers returned ierr!=0")
     end if

     ABI_FREE(bsize_frecords)
     call cwtime_report(" write_frmarkers", cpu, wall, gflops)
   end if
#endif

#ifdef HAVE_NETCDF
 CASE (IO_MODE_ETSF)
   !NCF_CHECK(nctk_open_modify(wfk%fh, wfk%fname, wfk%comm))

   if (nctk_has_mpiio) then
#ifdef HAVE_NETCDF_MPI
     ncerr = nf90_create(wfk%fname, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), &
         comm=wfk%comm, info=xmpio_info, ncid=wfk%fh)

     NCF_CHECK_MSG(ncerr, sjoin("nf90_create: ", wfk%fname))
#else
     MSG_ERROR("You should not be here")
#endif
   else
     if (wfk%nproc > 1) then
       MSG_ERROR("Your netcdf library does not support MPI-IO. Cannot write WFK file with nprocs > 1")
     end if

     ncerr = nf90_create(wfk%fname, nf90_write, wfk%fh)
     NCF_CHECK_MSG(ncerr, sjoin("nf90_create: ", wfk%fname))
   end if

   call wfk_ncdef_dims_vars(wfk%fh, hdr, wfk%fform, write_hdr=.True.)
   NCF_CHECK(nctk_def_basedims(wfk%fh))

   ! Switch to data mode.
   NCF_CHECK(nctk_set_datamode(wfk%fh))
#endif

 case default
   MSG_ERROR(sjoin('Wrong/unsupported iomode: ', itoa(wfk%iomode)))
 end select

 DBG_EXIT("COLL")

end subroutine wfk_open_write
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_close
!! NAME
!!  wfk_close
!!
!! FUNCTION
!!  Close the wavefunction file handler and release the memory allocated
!!  Delete the file if `delete` is True. Default: False
!!
!! PARENTS
!!      conducti_nc,d2frnl,dfpt_nstdy,dfpt_nstpaw,dfpt_scfcv,fold2Bloch,initwf
!!      ioprof,m_cut3d,m_iowf,m_wfd,m_wfk,mrggkk,optic
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_close(Wfk, delete)

!Arguments ------------------------------------
!scalars
 class(wfk_t),intent(inout) :: Wfk
 logical,optional,intent(in) :: delete

!Local variables-------------------------------
!scalars
 integer :: ierr
 !character(len=500) :: msg
#ifdef HAVE_MPI_IO
 integer :: mpierr,nfrec
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecords(:)
#endif

! *************************************************************************

 DBG_ENTER("COLL")

 !@wfk_t

 ! Close the file only if it was open.
 if (wfk%rw_mode /= WFK_NOMODE) then
   Wfk%rw_mode = WFK_NOMODE

   select case (Wfk%iomode)
   case (IO_MODE_FORTRAN)
      close(wfk%fh)

#ifdef HAVE_MPI_IO
   case (IO_MODE_MPI)
     call MPI_FILE_CLOSE(Wfk%fh,mpierr)
     ABI_CHECK_MPI(mpierr,"FILE_CLOSE!")

     if (Wfk%debug .and. Wfk%my_rank == Wfk%master) then
       ! Check the fortran records.
       call MPI_FILE_OPEN(xmpi_comm_self, Wfk%fname, MPI_MODE_RDONLY, xmpio_info, Wfk%fh, mpierr)
       ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")
       call hdr_bsize_frecords(Wfk%Hdr,Wfk%formeig,nfrec,bsize_frecords)
       call xmpio_check_frmarkers(Wfk%fh,Wfk%hdr_offset,xmpio_single,nfrec,bsize_frecords,ierr)
       ABI_CHECK(ierr==0,"xmpio_check_frmarkers returned ierr!=0")
       ABI_FREE(bsize_frecords)
       call MPI_FILE_CLOSE(Wfk%fh,mpierr)
       ABI_CHECK_MPI(mpierr,"FILE_CLOSE!")
     end if
#endif

#ifdef HAVE_NETCDF
   case (IO_MODE_ETSF)
     NCF_CHECK(nf90_close(wfk%fh))
#endif

   case default
     MSG_ERROR(sjoin('Wrong/unsupported value of iomode: ', itoa(Wfk%iomode)))
   end select
 end if

 ! Free memory.
 call Wfk%Hdr%free()

 ABI_SFREE(Wfk%nband)
 ABI_SFREE(Wfk%recn_ks)
 ABI_SFREE(Wfk%offset_ks)

 if (present(delete)) then
   if (delete) call delete_file(wfk%fname, ierr)
 end if

 DBG_EXIT("COLL")

end subroutine wfk_close
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_print
!! NAME
!!  wfk_print
!!
!! FUNCTION
!!  Print information on the object.
!!
!! INPUTS
!!  wfk<class(wfk_t)> = WFK handler
!!  [header]=String to be printed as header for additional info.
!!  [unit]=Unit number for output. Defaults to std_out
!!  [prtvol]=Verbosity level
!!
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_print(wfk,unit,header,prtvol)

!Arguments ------------------------------------
!scalars
 class(wfk_t),intent(inout) :: wfk
 integer,optional,intent(in) :: unit,prtvol
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 integer,parameter :: rdwr4=4
 integer :: my_unt,my_prtvol
 character(len=500) :: msg

! *************************************************************************

 my_unt = std_out; if (present(unit)) my_unt = unit
 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 msg=' ==== Info on the wfk_t object ==== '
 if (present(header)) msg=' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(my_unt,msg)
 call wrtout(my_unt, sjoin(" iomode = ",itoa(wfk%iomode)))

 call wfk%hdr%echo(wfk%fform, rdwr4 ,unit=my_unt)

end subroutine wfk_print
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_validate_ks
!! NAME
!!  wfk_validate_ks
!!
!! FUNCTION
!!  Validate the k-point, the spin index and, optionally, the band index.
!!  Return non-zero value if error.
!!
!! INPUTS
!!  wfk<class(wfk_t)> = WFK handler
!!  ik_ibz=k-point index.
!!  spin=Spin index.
!!  [band]=Band index.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function wfk_validate_ks(wfk, ik_ibz, spin, band) result(ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz, spin
 integer,optional,intent(in) :: band
 class(wfk_t),intent(in) :: wfk

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************
 ierr = 0

 if (ik_ibz <= 0 .or. ik_ibz > wfk%nkpt) then
   ierr = ierr + 1
   write(msg, '(2(a,i0))')'ik_ibz = ',ik_ibz,' whereas it should be between 1 and ',wfk%nkpt
   MSG_WARNING(msg)
 end if

 if (spin <= 0 .or. spin > wfk%nsppol) then
   ierr = ierr + 1
   write(msg, '(2(a,i0))')'spin = ',spin,' whereas it should be between 1 and ',wfk%nsppol
   MSG_WARNING(msg)
 end if

 if (present(band)) then
   if (band <=0) then
     ierr = ierr + 1
     MSG_WARNING(sjoin('Negative band index: band = ',itoa(band)))
   end if

   ! Don't touch nband array if wrong indices.
   if (spin > 0 .and. spin <= wfk%nsppol .and. ik_ibz > 0 .and. ik_ibz <= wfk%nkpt) then
      if (band > wfk%nband(ik_ibz, spin)) then
        ierr = ierr + 1
        write(msg, '(2(a,i0))')'band = ',band,' whereas it should be between 1 and ',wfk%nband(ik_ibz,spin)
        MSG_WARNING(msg)
      end if
   end if
 end if

 !if (ierr /= 0) then
 !  MSG_ERROR("Wrong (ik_ibz, spin) args, Aborting now")
 !end if

end function wfk_validate_ks
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_findk
!! NAME
!!  wfk_findk
!!
!! FUNCTION
!!  Find the index of the k-point in the WKF file. umklapp vectors are not allowed.
!!  Return -1 if not found.
!!
!! INPUTS
!!  wfk<class(wfk_t)> = WFK handler initialized and set in read mode
!!  kpt(3)=k-point in reduced coordinates.
!!  [ktol]=Optional tolerance for k-point comparison.
!!         For each reduced direction the absolute difference between the coordinates must be less that ktol
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function wfk_findk(wfk, kpt, ktol) result(ikpt)

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: ktol
 class(wfk_t),intent(in) :: wfk
!arrays
 real(dp),intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 integer :: ik
 real(dp) :: my_ktol

! *************************************************************************

 my_ktol = 0.0001_dp; if (present(ktol)) my_ktol = ktol

 ikpt = -1
 do ik=1,wfk%hdr%nkpt
   if (all(abs(wfk%hdr%kptns(:, ik) - kpt) < my_ktol)) then
     ikpt = ik; exit
   end if
 end do

end function wfk_findk
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_ncdef_dims_vars
!! NAME
!!  wfk_ncdef_dims_vars
!!
!! FUNCTION
!!  Write the heeder, fform as well as the etsf-io dimensions and variables.
!!
!! INPUTS
!!  ncid=netcdf file handler.
!!  hdr<hdr_tyep>=Abinit header
!!  fform=File type format of the header
!!  [write_hdr]=True if the header should be written (default)
!!  [iskss]=True if this is a KSS file (activate kdependent=No)
!!
!! PARENTS
!!      m_io_kss,m_iowf,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_ncdef_dims_vars(ncid, hdr, fform, write_hdr, iskss)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,fform
 type(hdr_type),intent(in) :: hdr
 logical,optional,intent(in) :: write_hdr,iskss

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 character(len=500) :: title,history
 logical :: do_write_hdr,my_iskss
 integer :: ivar,mpw,ncerr
! *************************************************************************

 do_write_hdr = .True.; if (present(write_hdr)) do_write_hdr = write_hdr
 my_iskss = .False.; if (present(iskss)) my_iskss = iskss
 if (do_write_hdr) then
   NCF_CHECK(hdr%ncwrite(ncid, fform, nc_define=.True.))
 end if

 ! Add the etsf header.
 title = sjoin("WFK file generated by Abinit, version: ",  abinit_version)
 if (my_iskss) title = sjoin("KSS file generated by Abinit, version: ",  abinit_version)
 history = sjoin("Generated on: ", asctime())
 NCF_CHECK(nctk_add_etsf_header(ncid, title=title, history=history))

 mpw = MAXVAL(hdr%npwarr)
 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t("real_or_complex_coefficients", 2), nctkdim_t("max_number_of_coefficients", mpw)])
 NCF_CHECK(ncerr)

 ! Define kg_k
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("reduced_coordinates_of_plane_waves", "int", &
     "number_of_reduced_dimensions, max_number_of_coefficients, number_of_kpoints")&
 ])
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_inq_varid(ncid, "reduced_coordinates_of_plane_waves", ivar))
 if (my_iskss) then
   NCF_CHECK(nf90_put_att(ncid, ivar, "k_dependent", "no"))
 else
   NCF_CHECK(nf90_put_att(ncid, ivar, "k_dependent", "yes"))
 end if

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("eigenvalues", "dp", "max_number_of_states, number_of_kpoints, number_of_spins") &
 ])
 NCF_CHECK(ncerr)
 NCF_CHECK(nctk_set_atomic_units(ncid, "eigenvalues"))

 ncerr = nctk_def_arrays(ncid, nctkarr_t("coefficients_of_wavefunctions", "dp", &
   "real_or_complex_coefficients, max_number_of_coefficients, number_of_spinor_components, &
&max_number_of_states, number_of_kpoints, number_of_spins"))
 NCF_CHECK(ncerr)

 !NF90_DEF_VAR_FILL(INTEGER NCID, INTEGER VARID, INTEGER NO_FILL, FILL_VALUE)
 !NCF_CHECK(nf90_inq_varid(ncid, "coefficients_of_wavefunctions", ivar))
 !NCF_CHECK(nf90_def_var_fill(ncid, ivar, 0, -one))

#else
 MSG_ERROR("netcdf not available")
#endif

end subroutine wfk_ncdef_dims_vars
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_compare
!! NAME
!!  wfk_compare
!!
!! FUNCTION
!!  Test two wfk_t objects for consistency. Return non-zero value if test fails.
!!
!! INPUTS
!!  wfk1, wfk2 <class(wfk_t)> = WFK handlers to be compared
!!
!! OUTPUT
!!  ierr
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function wfk_compare(wfk1, wfk2) result(ierr)

!Arguments ------------------------------------
!scalars
 class(wfk_t),intent(in) :: wfk1, wfk2

!Local variables-------------------------------
!scalars
 integer :: restart,restartpaw
 !character(len=500) :: msg

!************************************************************************

 ierr = 0

 ierr=wfk1%hdr%compare(wfk2%hdr)

 ! Test basic dimensions
!if (wfk1%hdr%nsppol /= wfk2%hdr%nsppol) then
!  ierr = ierr + 1; MSG_WARNING("Different nsppol")
!end if
!if (wfk1%hdr%nspinor /= wfk2%hdr%nspinor) then
!  ierr = ierr + 1; MSG_WARNING("Different nspinor")
!end if
!if (wfk1%hdr%nspden /= wfk2%hdr%nspden) then
!  ierr = ierr + 1; MSG_WARNING("Different nspden")
!end if
!if (wfk1%hdr%nkpt /= wfk2%hdr%nkpt) then
!  ierr = ierr + 1; MSG_WARNING("Different nkpt")
!end if
 if (wfk1%formeig /= wfk2%formeig) then
   ierr = ierr + 1; MSG_WARNING("Different formeig")
 end if
!if (wfk1%hdr%usepaw /= wfk2%hdr%usepaw) then
!  ierr = ierr + 1; MSG_WARNING("Different usepaw")
!end if
!if (wfk1%hdr%ntypat /= wfk2%hdr%ntypat) then
!  ierr = ierr + 1; MSG_WARNING("Different ntypat")
!end if
!if (wfk1%hdr%natom /= wfk2%hdr%natom) then
!  ierr = ierr + 1; MSG_WARNING("Different natom")
!end if
 !if (wfk1%hdr%fform /= wfk2%hdr%fform) then
 !  ierr = ierr + 1; MSG_WARNING("Different fform")
 !end if

 ! Return immediately if important dimensions are not equal.
 if (ierr /= 0) return

 ! Test important arrays (rprimd is not tested)
!if (any(wfk1%hdr%typat /= wfk2%hdr%typat)) then
!  ierr = ierr + 1; MSG_WARNING("Different typat")
!end if
!if (any(wfk1%hdr%npwarr /= wfk2%hdr%npwarr)) then
!  ierr = ierr + 1; MSG_WARNING("Different npwarr array")
!end if
 if (any(wfk1%nband /= wfk2%nband)) then
   ierr = ierr + 1; MSG_WARNING("Different nband array")
 end if
!if (any(abs(wfk1%hdr%kptns - wfk2%hdr%kptns) > tol6)) then
!  ierr = ierr + 1; MSG_WARNING("Different kptns array")
!end if

 ! Call hdr_check to get a nice diff of the header but don't check restart and restartpaw.
 call hdr_check(wfk1%fform,wfk2%fform,wfk1%hdr,wfk2%hdr,"PERS",restart,restartpaw)

end function wfk_compare
!!***
!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_band_block
!! NAME
!!  wfk_read_band_block
!!
!! FUNCTION
!!  Read a block of contiguous bands at a given (k-point, spin)
!!
!! INPUTS
!!  Wfk<class(wfk_t)>= WFK file handler object.
!!  band_block(2)=Initial and final band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading in wfk%comm (use it wisely!)
!!
!! OUTPUTS
!!  [kg_k=(:,:)] = G-vectors
!!  [eig_k(:)] = Eigenvectors
!!  [cg_k(:,:)] = Fourier coefficients
!!
!! NOTES
!!  The output arrays eig_k and occ_k contain the *full* set of eigenvalues and occupation
!!  factors stored in the file and are dimensioned with wfk%mband.
!!
!! PARENTS
!!      fold2Bloch,initwf,m_cut3d,m_wfd,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_read_band_block(Wfk,band_block,ik_ibz,spin,sc_mode,kg_k,cg_k,eig_k,occ_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode
 class(wfk_t),intent(inout) :: Wfk
!arrays
 integer,intent(in) :: band_block(2)
 integer,intent(out), DEV_CONTARRD  optional :: kg_k(:,:) !(3,npw_k)
 real(dp),intent(out), DEV_CONTARRD optional :: cg_k(:,:) !(2,npw_k*nspinor*nband)
 real(dp),intent(inout),optional :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),intent(out),optional :: occ_k(Wfk%mband)

!Local variables-------------------------------
!scalars
 integer :: ierr,npw_disk,nspinor_disk,nband_disk,band
 integer :: ipw,my_bcount,npwso,npw_tot,nb_block,base
 integer :: npw_read,nspinor_read,nband_read
 character(len=500) :: msg,errmsg
!arrays
 real(dp),ABI_CONTIGUOUS pointer :: tmp_eigk(:),tmp_occk(:)
#ifdef HAVE_MPI_IO
 integer :: mpierr,bufsz,gkk_type,cgblock_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 integer :: sizes(2),subsizes(2),starts(2),types(2)
#endif
#ifdef HAVE_NETCDF
 integer :: kg_varid,eig_varid,occ_varid,cg_varid,ncerr,h1_varid,idx,ib1,ib2
 real(dp),allocatable :: h1mat(:,:,:)
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfk%rw_mode==WFK_READMODE, "Wfk must be in READMODE")

 if (wfk_validate_ks(wfk, ik_ibz, spin) /= 0) then
   MSG_ERROR("Wrong (ik_ibz, spin) args, Aborting now")
 end if

 ! Look before you leap.
 npw_disk     = Wfk%Hdr%npwarr(ik_ibz)
 nspinor_disk = Wfk%nspinor
 nband_disk   = Wfk%nband(ik_ibz,spin)
 nb_block     = (band_block(2) - band_block(1) + 1)
 ABI_CHECK(nb_block>0,"nband <=0")
 npw_tot      = npw_disk * nspinor_disk * nb_block

 if (present(kg_k)) then
   ABI_CHECK(SIZE(kg_k,DIM=2) >= npw_disk,"kg_k too small")
 end if

 if (present(cg_k)) then
   ABI_CHECK(SIZE(cg_k, DIM=2) >= npw_tot,"cg_k too small")
 end if

 if (present(eig_k)) then
   if (Wfk%formeig==0) then
      ABI_CHECK(SIZE(eig_k) >= nband_disk, "GS eig_k too small")
   else if (Wfk%formeig==1) then
      ABI_CHECK(SIZE(eig_k) >= 2*nband_disk**2, "DFPT eig_k too small")
   else
     MSG_ERROR("formeig != [0,1]")
   end if
 end if

 if (present(occ_k)) then
   ABI_CHECK(SIZE(occ_k) >= nband_disk, "GS eig_k too small")
   if (Wfk%formeig==1) then
     MSG_ERROR("occ_k cannot be used when formeig ==1")
   end if
 end if

 select case (Wfk%iomode)
 case (IO_MODE_FORTRAN)

   ! Rewind the file to have the correct (k,s) block (if needed)
   call wfk_seek(Wfk,ik_ibz,spin)
   !
   ! Read the first record: npw, nspinor, nband_disk
   read(Wfk%fh, err=10, iomsg=errmsg) npw_read, nspinor_read, nband_read

   if (any( [npw_read, nspinor_read, nband_read] /= [npw_disk, nspinor_disk, nband_disk])) then
     write(msg,"(a,6(i0,2x))")"Mismatch between (npw, nspinor, nband) read from WFK and those found in HDR ",&
&      npw_read, nspinor_read, nband_read, npw_disk, nspinor_disk, nband_disk
     MSG_ERROR(msg)
   end if

   ! The second record: (k+G) vectors
   if (present(kg_k)) then
     read(Wfk%fh, err=10, iomsg=errmsg) kg_k(1:3,1:npw_disk)
   else
     read(Wfk%fh, err=10, iomsg=errmsg) ! kg_k(1:3,1:npw_disk)
   end if

   select case (Wfk%formeig)
   case (0)
     ! The third record: eigenvalues and occupation factors.
     ! write(unitwf) (eigen(iband),iband=1,nband_disk),(occ(iband),iband=1,nband_disk)
     if (present(eig_k) .or. present(occ_k)) then

       ABI_MALLOC(tmp_eigk, (nband_disk))
       ABI_MALLOC(tmp_occk, (nband_disk))

       read(Wfk%fh, err=10, iomsg=errmsg) tmp_eigk, tmp_occk

       if (present(eig_k)) eig_k = tmp_eigk
       if (present(occ_k)) occ_k = tmp_occk

       ABI_FREE(tmp_eigk)
       ABI_FREE(tmp_occk)

     else
       read(Wfk%fh, err=10, iomsg=errmsg) ! eig_k(1:nband_disk), occ_k(1:nband_k)
     end if

     ! Fourth record with the wave-functions.
     if (present(cg_k)) then
       npwso = npw_disk*nspinor_disk
       my_bcount = 0
       do band=1,nband_disk
         if (band >= band_block(1) .and. band <= band_block(2)) then
           ipw = my_bcount * npwso
           my_bcount = my_bcount + 1
           read(Wfk%fh, err=10, iomsg=errmsg) cg_k(1:2,ipw+1:ipw+npwso)
         else
           read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
         end if
       end do

     else
       do band=1,nband_disk
         read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       end do
     end if

   case (1)
     npwso = npw_disk*nspinor_disk
     my_bcount = 0
     do band=1,nband_disk

       if (present(eig_k)) then
         ! Read column matrix of size (2*nband_k)
         base = 2*(band-1)*nband_disk
         read(Wfk%fh, err=10, iomsg=errmsg) eig_k(base+1:base+2*nband_disk)
       else
         read(Wfk%fh, err=10, iomsg=errmsg ) ! eig_k(2*nband_disk)
       end if

       if (present(cg_k) .and. (band >= band_block(1) .and. band <= band_block(2)) ) then
         ipw = my_bcount * npwso
         my_bcount = my_bcount + 1
         read(Wfk%fh, err=10, iomsg=errmsg) cg_k(1:2,ipw+1:ipw+npwso)
       else
         read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       end if
     end do

   case default
     MSG_ERROR("formeig != [0,1]")
   end select

   ! Reached the end of the (k,s) block. Update f90_fptr
   if (ik_ibz < Wfk%nkpt) then
     Wfk%f90_fptr = [ik_ibz+1,spin,REC_NPW]
   else
     ABI_CHECK(ik_ibz == wfk%nkpt, "ik_ibz != nkpt")
     if (spin==Wfk%nsppol) then
       Wfk%f90_fptr = FPTR_EOF ! EOF condition
     else
       Wfk%f90_fptr = [1,spin+1,REC_NPW]
     end if
   end if

#ifdef HAVE_MPI_IO
 case (IO_MODE_MPI)
   if (present(kg_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) + xmpio_bsize_frm

     call mpio_read_kg_k(Wfk%fh,my_offset,npw_disk,sc_mode,kg_k,mpierr)
     ABI_CHECK_MPI(mpierr,"reading kg")
   end if

   ! formeig=0 =>  Read both eig and occ in tmp_eigk
   ! formeig=1 =>  Read (nband_k,nband_k) matrix of complex numbers.
   select case (Wfk%formeig)
   case (0)
     if (present(eig_k) .or. present(occ_k)) then
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + xmpio_bsize_frm

       call mpio_read_eigocc_k(Wfk%fh,my_offset,nband_disk,Wfk%formeig,sc_mode,tmp_eigk,mpierr)
       ABI_CHECK_MPI(mpierr,"reading eigocc")

       if (present(eig_k)) eig_k(1:nband_disk) = tmp_eigk(1:nband_disk)
       if (present(occ_k)) occ_k(1:nband_disk) = tmp_eigk(nband_disk+1:)

       ABI_FREE(tmp_eigk)
     end if

     if (present(cg_k)) then
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG)
       sizes     = [npw_disk*nspinor_disk, nband_disk]
       subsizes  = [npw_disk*nspinor_disk, band_block(2)-band_block(1)+1]
       bufsz     = 2 * npw_disk * nspinor_disk * nb_block
       starts    = [1, band_block(1)]

       call mpiotk_read_fsuba_dp2D(Wfk%fh,my_offset,sizes,subsizes,starts,bufsz,cg_k,&
         wfk%chunk_bsize,sc_mode,Wfk%comm,ierr)
       ABI_CHECK(ierr==0,"Fortran record too big")
     end if

   case (1)
     if (present(eig_k)) then
       sizes = [nband_disk, npw_disk*nspinor_disk]
       types = [MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX]

       call xmpio_create_fstripes(nband_disk,sizes,types,gkk_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + my_offpad

       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,gkk_type,'native',xmpio_info,mpierr)
       ABI_CHECK_MPI(mpierr,"")
       call MPI_TYPE_FREE(gkk_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")

       bufsz = nband_disk**2

       if (sc_mode==xmpio_collective) then
         call MPI_FILE_READ_ALL(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_READ(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else
         MSG_ERROR("Wrong sc_mode")
       end if
       ABI_CHECK_MPI(mpierr,"FILE_READ")
     end if

     if (present(cg_k)) then
       types = [MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX]
       sizes = [npw_disk*nspinor_disk, nband_disk]

       call xmpio_create_fstripes(nb_block,sizes,types,cgblock_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

       ! Increment my_offset to account for previous eigen and cg records if band_block(1) != 1
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG) + my_offpad
       my_offset = my_offset + (band_block(1) - 1) * ( &
          (2*npw_disk*wfk%nspinor*xmpi_bsize_dp + 2*xmpio_bsize_frm) + &
          (2*nband_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm) )

       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,cgblock_type,'native',xmpio_info,mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")

       call MPI_TYPE_FREE(cgblock_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")

       bufsz = npw_disk * nspinor_disk * nb_block

       if (sc_mode==xmpio_collective) then
         call MPI_FILE_READ_ALL(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_READ(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else
         MSG_ERROR("Wrong sc_mode")
       end if
       ABI_CHECK_MPI(mpierr,"FILE_READ")
     end if

   case default
     MSG_ERROR("formeig != [0,1]")
   end select
#endif

#ifdef HAVE_NETCDF
 case (IO_MODE_ETSF)
   if (present(kg_k)) then
     ! Read the reduced_coordinates_of_plane_waves for this k point.
     NCF_CHECK(nf90_inq_varid(wfk%fh, "reduced_coordinates_of_plane_waves", kg_varid))
     if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, kg_varid))
     end if
     ncerr = nf90_get_var(wfk%fh, kg_varid, kg_k, start=[1,1,ik_ibz], count=[3,npw_disk,1])
     NCF_CHECK(ncerr)
   end if

   if (Wfk%formeig==0) then
     ! Read eigenvalues and occupations.
     if (present(eig_k)) then
       NCF_CHECK(nf90_inq_varid(wfk%fh, "eigenvalues", eig_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, eig_varid))
       end if
       ncerr = nf90_get_var(wfk%fh, eig_varid, eig_k, start=[1,ik_ibz,spin], count=[nband_disk,1,1])
       NCF_CHECK(ncerr)
     end if

     if (present(occ_k)) then
       NCF_CHECK(nf90_inq_varid(wfk%fh, "occupations", occ_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, occ_varid))
       end if
       ncerr = nf90_get_var(wfk%fh, occ_varid, occ_k, start=[1,ik_ibz,spin], count=[nband_disk,1,1])
       NCF_CHECK_MSG(ncerr, "getting occ_k")
     end if

   else ! formeig == 1

     if (present(eig_k)) then
       ! Read h1 matrix elements. The netcdf array has shape:
       ! [complex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins]
       ABI_MALLOC(h1mat, (2, wfk%mband, wfk%mband))
       NCF_CHECK(nf90_inq_varid(wfk%fh, "h1_matrix_elements", h1_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, h1_varid))
       end if
       ncerr = nf90_get_var(wfk%fh, h1_varid, h1mat, start=[1,1,1,ik_ibz,spin], &
         count=[2, wfk%mband, wfk%mband, 1, 1])
       NCF_CHECK_MSG(ncerr, "getting h1mat_k")

       ! For legacy reasons, I have to pack nband_k**2 elements in the first positions in eig_k
       ! This is important only if nband(:) depends on k i.e. mband != nband_k
       idx=1
       do ib2=1,nband_disk
         do ib1=1,nband_disk
            eig_k(idx:idx+1) = h1mat(:2,ib1,ib2)
            idx=idx+2
          end do
       end do
       ABI_FREE(h1mat)
     end if
   end if

   if (present(cg_k)) then
     ! Read the nb_block bands starting from band_block(1)
     ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
     NCF_CHECK(nf90_inq_varid(wfk%fh, "coefficients_of_wavefunctions", cg_varid))
     if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, cg_varid))
     end if

     ncerr = nf90_get_var(wfk%fh, cg_varid, cg_k, start=[1,1,1,band_block(1),ik_ibz,spin], &
       count=[2,npw_disk,wfk%nspinor,nb_block,1,1])
     NCF_CHECK_MSG(ncerr, "getting cg_k")
  end if
#endif

 case default
   MSG_ERROR(sjoin('Wrong/unsupported iomode: ', itoa(Wfk%iomode)))
 end select

 DBG_EXIT("COLL")

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine wfk_read_band_block
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_bks
!! NAME
!!  wfk_read_bks
!!
!! FUNCTION
!!  Read the wavefunction and, optionally, the *DFPT* matrix elements
!!  for a given (band, k-point, spin).
!!
!! INPUTS
!!  Wfk<class(wfk_t)>=WFK file handler.
!!  band=Band index
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  cg_bks(2,npw_k*nspinor) = Fourier coefficients of the wavefunction
!!  [eig1_bks(2*wfk%mband)] = Matrix elements of the DFPT H1 Hamiltonian at the specified (k, spin).
!!
!! PARENTS
!!      d2frnl,dfpt_nstpaw,dfpt_nstwf,dfpt_vtowfk,rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_read_bks(wfk, band, ik_ibz, spin, sc_mode, cg_bks, eig1_bks)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin,sc_mode
 class(wfk_t),intent(inout) :: wfk
!arrays
 real(dp),DEV_CONTARRD intent(out) :: cg_bks(:,:)
 real(dp),optional,intent(inout) :: eig1_bks(2*wfk%mband)

!Local variables-------------------------------
!scalars
 integer :: start,ib1,nspinor_disk,npw_disk,nband_disk !ierr,
#ifdef HAVE_MPI_IO
 integer :: mpierr,bufsz,gkk_type !,cg_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 integer :: sizes(2),types(2)
#endif
#ifdef HAVE_NETCDF
 integer :: h1_varid,cg_varid,ncerr
#endif
 character(len=500) :: errmsg
!arrays
 real(dp),allocatable :: all_eigk(:)

!************************************************************************

 if (wfk_validate_ks(wfk, ik_ibz, spin, band=band) /= 0) then
   MSG_ERROR("Wrong (ik_ibz, spin, band) args, Aborting now")
 end if
 npw_disk = wfk%Hdr%npwarr(ik_ibz)
 nband_disk = wfk%nband(ik_ibz, spin)
 nspinor_disk = wfk%nspinor
 !cg_bks = one; if (present(eig1_bks)) eig1_bks = zero ; return

 if (.not. present(eig1_bks)) then
   call wfk%read_band_block([band, band], ik_ibz, spin, sc_mode, cg_k=cg_bks)
   return

 else
   ABI_CHECK(wfk%formeig==1, "formeig must be 1 if eig1_bks is present")
   ABI_CHECK(size(cg_bks, dim=2) >= npw_disk*wfk%nspinor,"cg_bks too small")
   ABI_CHECK(size(eig1_bks) >= 2*nband_disk, "eig1_bks too small")

 if (.False.) then
 !if (.True.) then
   ! Due to the API of wfk_read_band_block, we have to read the full set of eigenvalues
   ! and then extract the relevant band
   ! TODO: Should write another routine to avoid doing that.
   ABI_MALLOC(all_eigk, (2*wfk%mband**2))
   call wfk%read_band_block([band, band], ik_ibz, spin, sc_mode, cg_k=cg_bks, eig_k=all_eigk)

   ! Find the index of the slice.
   ! Remember that data in all_eigk does not have a constant stride if nband_disk != mband
   start = (band-1)*2*nband_disk

   !write(std_out,*)size(eig1_bks), nband_disk
   eig1_bks(1:2*nband_disk) = all_eigk(start+1:start+2*nband_disk)
   ABI_FREE(all_eigk)

 else
   ! Improved version
   select case (wfk%iomode)
   case (IO_MODE_FORTRAN)
     ! This code is not optimal because I'm rewinding the (k, s) block
     ! at each call but I'm not gonna spend time on this because I should
     ! refactor a lot of stuff. Use MPI-IO or HDF5!
     call wfk_seek(wfk,ik_ibz,spin)

     ! Read the first record: npw, nspinor, nband_disk
     read(Wfk%fh, err=10, iomsg=errmsg) !npw_read, nspinor_read, nband_read
     ! The second record: (k+G) vectors
     read(wfk%fh, err=10, iomsg=errmsg) ! kg_k(1:3,1:npw_disk)

     do ib1=1,nband_disk
       if (ib1 /= band) then
         read(Wfk%fh, err=10, iomsg=errmsg) ! eig_k(base+1:base+2*nband_disk)
         read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       else
         if (present(eig1_bks)) then
           read(wfk%fh, err=10, iomsg=errmsg) eig1_bks(1:2*nband_disk)
         else
           read(wfk%fh, err=10, iomsg=errmsg)
         end if
         read(wfk%fh, err=10, iomsg=errmsg) cg_bks(:,:npw_disk*wfk%nspinor)
       end if
     end do

     ! Reached the end of the (k,s) block. Update f90_fptr
     call wfk_update_f90ptr(wfk, ik_ibz, spin)

#ifdef HAVE_MPI_IO
   case (IO_MODE_MPI)
     ! MPI-IO operations are done with file views even for contiguous data of the same type.
     ! See NOTES at the beginning of this module.

     sizes = [nband_disk, npw_disk*nspinor_disk]
     types = [MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX]

     call xmpio_create_fstripes(1,sizes,types,gkk_type,my_offpad,mpierr)
     ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

     !call MPI_TYPE_CONTIGUOUS(nband_disk,MPI_DOUBLE_COMPLEX,gkk_type,mpierr)
     !ABI_CHECK_MPI(mpierr,"type_contigous")
     !call MPI_TYPE_COMMIT(gkk_type,mpierr)
     !ABI_CHECK_MPI(mpierr,"mpi_commit")
     !my_offpad = 0

     ! Increment my_offset to account for previous (iband -1) bands.
     my_offset = wfk%offset_ks(ik_ibz,spin,REC_EIG) + my_offpad
     my_offset = my_offset + (band - 1) * ( &
        (2*npw_disk*wfk%nspinor*xmpi_bsize_dp + 2*xmpio_bsize_frm) + &
        (2*nband_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm) )

#if 1
     bufsz = nband_disk
     if (sc_mode==xmpio_collective) then
       call MPI_FILE_READ_AT_ALL(wfk%fh,my_offset,eig1_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     else
       call MPI_FILE_READ_AT(wfk%fh,my_offset,eig1_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     end if

     ! Read the cg_ks(G)
     ! Increment my_offset to account for the previous eigen and cg records of (iband-1) bands.
     my_offset = wfk%offset_ks(ik_ibz,spin,REC_CG) + xmpio_bsize_frm
     my_offset = my_offset + (band - 1) * ( &
        (2*npw_disk*wfk%nspinor*xmpi_bsize_dp + 2*xmpio_bsize_frm) + &
        (2*nband_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm) )

     bufsz = npw_disk * nspinor_disk
     if (sc_mode==xmpio_collective) then
       call MPI_FILE_READ_AT_ALL(wfk%fh,my_offset,cg_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     else
       call MPI_FILE_READ_AT(wfk%fh,my_offset,cg_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     end if

#else
     call MPI_FILE_SET_VIEW(wfk%fh,my_offset,MPI_BYTE,gkk_type,'native',xmpio_info,mpierr)
     ABI_CHECK_MPI(mpierr,"")
     call MPI_TYPE_FREE(gkk_type,mpierr)
     ABI_CHECK_MPI(mpierr,"")

     bufsz = nband_disk
     if (sc_mode==xmpio_collective) then
       call MPI_FILE_READ_ALL(wfk%fh,eig1_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     else if (sc_mode==xmpio_single) then
       call MPI_FILE_READ(wfk%fh,eig1_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     else
       MSG_ERROR("Wrong sc_mode")
     end if
     ABI_CHECK_MPI(mpierr,"FILE_READ")

     ! Read the cg_ks(G)
     ! Increment my_offset to account for the previous eigen and cg records of (iband-1) bands.
     my_offset = wfk%offset_ks(ik_ibz,spin,REC_CG) + xmpio_bsize_frm
     my_offset = my_offset + (band - 1) * ( &
        (2*npw_disk*wfk%nspinor*xmpi_bsize_dp + 2*xmpio_bsize_frm) + &
        (2*nband_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm) )

     call MPI_TYPE_CONTIGUOUS(npw_disk*nspinor_disk,MPI_DOUBLE_COMPLEX,cg_type,mpierr)
     ABI_CHECK_MPI(mpierr,"type_contigous")
     call MPI_TYPE_COMMIT(cg_type,mpierr)
     ABI_CHECK_MPI(mpierr,"mpi_commit")

     call MPI_FILE_SET_VIEW(wfk%fh,my_offset,MPI_BYTE,cg_type,'native',xmpio_info,mpierr)
     ABI_CHECK_MPI(mpierr,"")
     call MPI_TYPE_FREE(cg_type, mpierr)
     ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

     bufsz = npw_disk * nspinor_disk
     if (sc_mode==xmpio_collective) then
       call MPI_FILE_READ_ALL(wfk%fh,cg_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     else if (sc_mode==xmpio_single) then
       call MPI_FILE_READ(wfk%fh,cg_bks,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
     else
       MSG_ERROR("Wrong sc_mode")
     end if
     ABI_CHECK_MPI(mpierr,"FILE_READ")
#endif
#endif

#ifdef HAVE_NETCDF
   case (IO_MODE_ETSF)
     ! Read h1 matrix elements. The netcdf array has shape:
     ! [complex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins]
     NCF_CHECK(nf90_inq_varid(wfk%fh, "h1_matrix_elements", h1_varid))
     if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, h1_varid))
     end if
     ncerr = nf90_get_var(wfk%fh, h1_varid, eig1_bks, start=[1,1,band,ik_ibz,spin], count=[2, nband_disk, 1, 1, 1])
     NCF_CHECK_MSG(ncerr, "getting h1mat_k")

     ! Read the wavefunction.  The coefficients_of_wavefunctions on file have shape:
     ! [cplex, mpw, nspinor, mband, nkpt, nsppol]
     NCF_CHECK(nf90_inq_varid(wfk%fh, "coefficients_of_wavefunctions", cg_varid))
     if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, cg_varid))
     end if

     ncerr = nf90_get_var(wfk%fh, cg_varid, cg_bks, start=[1,1,1,band,ik_ibz,spin], &
       count=[2,npw_disk,wfk%nspinor,1,1,1])
     NCF_CHECK_MSG(ncerr, "getting cg_k")
#endif

  case default
    MSG_ERROR(sjoin('Wrong value for iomode:', itoa(Wfk%iomode)))
  end select
 end if

 end if

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine wfk_read_bks
!!***

!!****f* m_wfk/wfk_write_band_block
!! NAME
!!  wfk_write_band_block
!!
!! FUNCTION
!!  Write a block of contigous bands.
!!
!! INPUTS
!!  Wfk<class(wfk_t)>=
!!  band_block(2)=Initial and final band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  [kg_k=(:,:)] = G-vectors
!!  [cg_k(:,:)]  = Fourier coefficients
!!  [eig_k(:)] = Eigenvectors
!!  [occ_k(:)] = Eigenvectors
!!
!! PARENTS
!!      m_iowf,m_wfd,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_write_band_block(Wfk,band_block,ik_ibz,spin,sc_mode,kg_k,cg_k,eig_k,occ_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode !,mband,rdcg,rdeig,npw_k,nband_k
 class(wfk_t),intent(inout) :: Wfk
!arrays
 integer,intent(in) :: band_block(2)
 integer,intent(in),optional :: kg_k(:,:)  !(3,npw_k)
 real(dp),intent(in),optional :: cg_k(:,:) ! cg_k(2,rdcg*cgsize2) !(2,npw_k*nspinor*nband)
 real(dp),intent(in),optional :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),intent(in),optional :: occ_k(Wfk%mband)

!Local variables-------------------------------
!scalars
 integer :: ierr,npw_disk,nspinor_disk,nband_disk,band
 integer :: ipw,my_bcount,npwso,npw_tot,nb_block,base
 character(len=500) :: errmsg !msg,
!arrays
 real(dp),ABI_CONTIGUOUS pointer :: tmp_eigk(:)
#ifdef HAVE_MPI_IO
 integer :: mpierr,bufsz,recnpw_type,gkk_type,cgblock_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 integer :: sizes(2),subsizes(2),starts(2),dims(3),types(2)
 integer(XMPI_OFFSET_KIND) :: bsize_rec(1)
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecords(:)
#endif
#ifdef HAVE_NETCDF
 integer :: kg_varid,eig_varid,occ_varid,cg_varid,ncerr
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfk%rw_mode==WFK_WRITEMODE, "Wfk must be in WRITEMODE")

 ! Look before you leap.
 npw_disk     = Wfk%Hdr%npwarr(ik_ibz)
 nspinor_disk = Wfk%nspinor
 nband_disk   = Wfk%nband(ik_ibz,spin)
 nb_block     = (band_block(2) - band_block(1) + 1)
 npw_tot      = npw_disk * nspinor_disk * nb_block

 if (PRESENT(kg_k)) then
   ABI_CHECK(SIZE(kg_k,DIM=2) >= npw_disk,"kg_k too small")
 end if

 if (PRESENT(cg_k)) then
   ABI_CHECK(SIZE(cg_k, DIM=2) >= npw_tot,"cg_k too small")
 end if

 if (PRESENT(eig_k)) then
   if (Wfk%formeig==0) then
      ABI_CHECK(SIZE(eig_k) >= nband_disk, "GS eig_k too small")
      ABI_CHECK(PRESENT(occ_k),"both eig_k and occ_k must be present")
   else if (Wfk%formeig==1) then
      ABI_CHECK(SIZE(eig_k) >= 2*nband_disk**2, "DFPT eig_k too small")
   else
     MSG_ERROR("formeig != [0,1]")
   end if
 end if

 if (PRESENT(occ_k)) then
   ABI_CHECK(SIZE(occ_k) >= nband_disk, "GS eig_k too small")
   ABI_CHECK(PRESENT(eig_k),"both eig_k and occ_k must be present")
   ABI_CHECK(Wfk%formeig == 0, "formeig /=0 with occ_k in input!")
 end if

 select case (Wfk%iomode)
 case (IO_MODE_FORTRAN)

   ! Rewind the file to have the correct (k,s) block (if needed)
   call wfk_seek(Wfk,ik_ibz,spin)

   ! Write the first record: npw, nspinor, nband_disk
   write(Wfk%fh, err=10, iomsg=errmsg) npw_disk, nspinor_disk, nband_disk

   ! The second record: (k+G) vectors
   if (PRESENT(kg_k)) then
     write(Wfk%fh, err=10, iomsg=errmsg) kg_k(1:3,1:npw_disk)
   else
     read(Wfk%fh, err=10, iomsg=errmsg) ! kg_k(1:3,1:npw_disk)
   end if

   ! The third record: eigenvalues and occupation factors.
   select case (Wfk%formeig)
   case (0)
     !write(unitwf) (eigen(iband),iband=1,nband_disk),(occ(iband),iband=1,nband_disk)

     if (present(eig_k) .and. present(occ_k)) then
       write(Wfk%fh, err=10, iomsg=errmsg) eig_k, occ_k
     else
       MSG_ERROR("Not coded")
       write(Wfk%fh, err=10, iomsg=errmsg) ! eig_k(1:nband_disk), occ_k(1:nband_k)
     end if

     ! The wave-functions.
     if (present(cg_k)) then
       npwso = npw_disk*nspinor_disk
       my_bcount = 0
       do band=1,nband_disk
         if (band >= band_block(1) .and. band <= band_block(2)) then
           ipw = my_bcount * npwso
           my_bcount = my_bcount + 1
           write(Wfk%fh, err=10, iomsg=errmsg) cg_k(1:2,ipw+1:ipw+npwso)
         else
           MSG_ERROR("Not coded")
           write(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
         end if
       end do

     else
       MSG_ERROR("Not coded")
       do band=1,nband_disk
         write(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       end do
     end if

   case (1)
     ! Write matrix of size (2*nband_k**2)
     ! The wave-functions.
     npwso = npw_disk*nspinor_disk
     my_bcount = 0

     do band=1,nband_disk
       base = 2*(band-1)*nband_disk
       write(Wfk%fh, err=10, iomsg=errmsg) eig_k(base+1:base+2*nband_disk)
       if (band >= band_block(1) .and. band <= band_block(2)) then
         ipw = my_bcount * npwso
         my_bcount = my_bcount + 1
         write(Wfk%fh, err=10, iomsg=errmsg) cg_k(1:2,ipw+1:ipw+npwso)
       else
         MSG_ERROR("Not coded")
         write(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       end if
     end do

   case default
     MSG_ERROR("formeig != [0,1]")
   end select

   ! Reached the end of the (k,s) block. Update f90_fptr
   call wfk_update_f90ptr(wfk, ik_ibz, spin)

#ifdef HAVE_MPI_IO
 case (IO_MODE_MPI)
   my_offset = Wfk%offset_ks(ik_ibz,spin,REC_NPW)

   bsize_rec(1) = 3 * xmpi_bsize_int
   call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,1,bsize_rec,ierr)
   ABI_CHECK(ierr==0,"ierr!=0")

   my_offset = Wfk%offset_ks(ik_ibz,spin,REC_NPW) + xmpio_bsize_frm

   call MPI_TYPE_CONTIGUOUS(3, MPI_INTEGER, recnpw_type, mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   call MPI_TYPE_COMMIT(recnpw_type,mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,recnpw_type,'native',xmpio_info,mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   call MPI_TYPE_FREE(recnpw_type,mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   dims = [npw_disk, nspinor_disk, nband_disk]

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(Wfk%fh,dims,SIZE(dims),MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE(Wfk%fh,dims,SIZE(dims),MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   if (present(kg_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG)

     bsize_rec(1) = 3 * npw_disk * xmpi_bsize_int
     call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,1,bsize_rec,ierr)

     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) + xmpio_bsize_frm

     call mpio_write_kg_k(Wfk%fh,my_offset,npw_disk,sc_mode,kg_k,mpierr)
     ABI_CHECK_MPI(mpierr,"mpio_write_kg_k")
   end if

   if (Wfk%formeig==0) then

     if (present(eig_k) .and. present(occ_k)) then

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG)

       bsize_rec(1) = 2 * nband_disk * xmpi_bsize_dp
       call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,1,bsize_rec,ierr)

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + xmpio_bsize_frm
       !
       ! Write both eig and occ in tmp_eigk
       bufsz = 2*nband_disk
       ABI_MALLOC(tmp_eigk, (bufsz))

       tmp_eigk(1:nband_disk)  = eig_k(1:nband_disk)
       tmp_eigk(nband_disk+1:) = occ_k(1:nband_disk)

       call mpio_write_eigocc_k(Wfk%fh,my_offset,nband_disk,Wfk%formeig,sc_mode,tmp_eigk,mpierr)
       ABI_CHECK_MPI(mpierr,"mpio_write_eigocc_k")

       ABI_FREE(tmp_eigk)
     end if

     if (present(cg_k)) then
       ABI_MALLOC(bsize_frecords, (nb_block))
       bsize_frecords = 2 * npw_disk * nspinor_disk * xmpi_bsize_dp
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG) + (band_block(1)-1) * (bsize_frecords(1) + 2*xmpio_bsize_frm)
       call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,nb_block,bsize_frecords,ierr)
       ABI_CHECK(ierr==0,"ierr!=0")
       ABI_FREE(bsize_frecords)

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG)
       sizes    = (/npw_disk*nspinor_disk, nband_disk/)
       subsizes = (/npw_disk*nspinor_disk, band_block(2)-band_block(1)+1/)
       bufsz = 2 * npw_disk * nspinor_disk * nb_block
       starts = [1, band_block(1)]

       call mpiotk_write_fsuba_dp2D(Wfk%fh,my_offset,sizes,subsizes,starts,bufsz,cg_k,Wfk%chunk_bsize,sc_mode,Wfk%comm,ierr)
       ABI_CHECK(ierr==0,"ierr!=0")
     end if

   else if (Wfk%formeig==1) then

     if (present(eig_k)) then
       types = [MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX]
       sizes = [nband_disk,npw_disk*nspinor_disk]

       call xmpio_create_fstripes(nband_disk,sizes,types,gkk_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + my_offpad

       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,gkk_type,'native',xmpio_info,mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")

       call MPI_TYPE_FREE(gkk_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")

       bufsz = (nband_disk**2)

       if (sc_mode==xmpio_collective) then
         call MPI_FILE_WRITE_ALL(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_WRITE(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else
         MSG_ERROR("Wrong sc_mode")
       end if

       ABI_CHECK_MPI(mpierr,"FILE_WRITE")
     end if

     if (present(cg_k)) then
       ABI_CHECK(band_block(1)==1,"band_block(1) !=1 not coded")

       types = [MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX]
       sizes = [npw_disk*nspinor_disk, nband_disk]

       call xmpio_create_fstripes(nb_block,sizes,types,cgblock_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG) + my_offpad

       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,cgblock_type,'native',xmpio_info,mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")

       call MPI_TYPE_FREE(cgblock_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")

       bufsz = npw_disk * nspinor_disk * nb_block
       if (sc_mode==xmpio_collective) then
         call MPI_FILE_WRITE_ALL(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_WRITE(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else
         MSG_ERROR("Wrong sc_mode")
       end if
       ABI_CHECK_MPI(mpierr,"FILE_WRITE")
     end if

   else
     MSG_ERROR("formeig not in [0,1]")
   end if
#endif

#ifdef HAVE_NETCDF
 case (IO_MODE_ETSF)
   if (present(kg_k)) then
     ! Write the reduced_coordinates_of_plane_waves for this k point.
     NCF_CHECK(nf90_inq_varid(wfk%fh, "reduced_coordinates_of_plane_waves", kg_varid))
     if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, kg_varid))
     end if
     ncerr = nf90_put_var(wfk%fh, kg_varid, kg_k, start=[1,1,ik_ibz], count=[3,npw_disk,1])
     NCF_CHECK_MSG(ncerr, "putting kg_k")
   end if

   ! Write eigenvalues and occupation factors.
   if (Wfk%formeig==0) then

     if (present(eig_k)) then
       NCF_CHECK(nf90_inq_varid(wfk%fh, "eigenvalues", eig_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, eig_varid))
       end if
       ncerr = nf90_put_var(wfk%fh, eig_varid, eig_k, start=[1,ik_ibz,spin], count=[nband_disk,1,1])
       NCF_CHECK_MSG(ncerr, "putting eig_k")
     end if

     if (present(occ_k)) then
       NCF_CHECK(nf90_inq_varid(wfk%fh, "occupations", occ_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, occ_varid))
       end if
       ncerr = nf90_put_var(wfk%fh, occ_varid, occ_k, start=[1,ik_ibz,spin], count=[nband_disk,1,1])
       NCF_CHECK_MSG(ncerr, "putting occ_k")
     end if

   else if (Wfk%formeig==1) then
     if (present(eig_k) .or. present(occ_k)) then
       MSG_ERROR("Don't pass eig_k or occ_k when formeig==1 and ETSF-IO")
     end if

   else
     MSG_ERROR("formeig != [0,1]")
   end if

   if (present(cg_k)) then
     ! Write the nb_block bands starting from band_block(1)
     ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
     NCF_CHECK(nf90_inq_varid(wfk%fh, "coefficients_of_wavefunctions", cg_varid))
    if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, cg_varid))
     end if

     ncerr = nf90_put_var(wfk%fh, cg_varid, cg_k, start=[1,1,1,band_block(1),ik_ibz,spin], &
                          count=[2,npw_disk,wfk%nspinor,nb_block,1,1])
     NCF_CHECK_MSG(ncerr, "putting cg_k")
  end if
#endif

 case default
   MSG_ERROR(sjoin('Wrong value of iomode:', itoa(Wfk%iomode)))
 end select

 DBG_EXIT("COLL")

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine wfk_write_band_block
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_bmask
!! NAME
!!  wfk_read_bmask
!!
!! FUNCTION
!!  Read a set of bands at a given k-point, spin. The bands to be read
!!  are specified by the logical mask `bmask`.
!!
!! INPUTS
!!  Wfk<class(wfk_t)>=
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  [kg_k=(:,:)] = G-vectors
!!  [cg_k(:,:)]  = Fourier coefficients
!!  [eig_k(:)] = Eigenvectors
!!  [occ_k(:)] = Occupation
!!
!! NOTES
!!  The output arrays eig_k and occ_k contain the *full* set of eigenvalues and occupation
!!  factors stored in the file and are dimensioned with wfk%mband.
!!
!! PARENTS
!!      m_wfd,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_read_bmask(Wfk,bmask,ik_ibz,spin,sc_mode,kg_k,cg_k,eig_k,occ_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode
 class(wfk_t),intent(inout) :: Wfk
!arrays
 logical,intent(in) :: bmask(Wfk%mband)
 integer,intent(out), DEV_CONTARRD optional :: kg_k(:,:)  !(3,npw_k)
 real(dp),intent(out), DEV_CONTARRD optional :: cg_k(:,:) !(2,npw_k*nspinor*nband)
 real(dp),intent(out),optional :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),intent(out),optional :: occ_k(Wfk%mband)

!Local variables-------------------------------
!scalars
 integer :: ierr,npw_disk,nspinor_disk,nband_disk,ipw,my_bcount,cnt,npwso,npw_tot,pt1,pt2,band
 integer :: npw_read,nspinor_read,nband_read,nb_tot,ncount,my_bcnt,my_maxb,base,nb
 character(len=500) :: msg,errmsg
!arrays
 real(dp),ABI_CONTIGUOUS pointer :: tmp_eigk(:),tmp_occk(:)
 integer :: mpierr,cgscatter_type,cg_type,method,block,nblocks,nbxblock
 integer :: bstart,bstop,bufsz,ugsz,brest,max_nband
 integer(XMPI_OFFSET_KIND) :: my_offset,base_ofs,my_offpad
 integer :: band_block(2),sizes(2),subsizes(2),starts(2),types(2)
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 real(dp),allocatable :: buffer(:,:)
#ifdef HAVE_NETCDF
 integer :: kg_varid,eig_varid,occ_varid,cg_varid,ncerr
 integer,allocatable :: blocks(:,:)
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfk%rw_mode==WFK_READMODE, "Wfk must be in READMODE")

 !do band=1,wfk%mband
 !  if (.not. bmask(band)) continue
 !  if (wfk_validate_ks(wfk, ik_ibz, spin, band=band) /= 0) then
 !    MSG_ERROR("Wrong (ik_ibz, spin, band) args, Aborting now")
 !  end if
 !end if

 !
 ! Look before you leap.
 npw_disk = Wfk%Hdr%npwarr(ik_ibz)
 nspinor_disk = Wfk%nspinor
 nband_disk = Wfk%nband(ik_ibz,spin)
 nb_tot = COUNT(bmask)
 npw_tot = npw_disk * nspinor_disk * nb_tot

 if (present(kg_k)) then
   ABI_CHECK((SIZE(kg_k,DIM=2) >= npw_disk),"kg_k too small")
 end if

 if (present(cg_k)) then
   ABI_CHECK(SIZE(cg_k, DIM=2) >= npw_tot, "Too small cg_k")
 end if

 if (present(eig_k)) then
   if (Wfk%formeig==0) then
      ABI_CHECK(SIZE(eig_k) >= nband_disk, "GS eig_k too small")
   else if (Wfk%formeig==1) then
      ABI_CHECK(SIZE(eig_k) >= 2*nband_disk**2, "DFPT eig_k too small")
   else
     MSG_ERROR("formeig != [0,1]")
   end if
 end if

 if (present(occ_k)) then
   ABI_CHECK(Wfk%formeig==0,"occ_k with formeig != 0")
   ABI_CHECK(SIZE(occ_k) >= nband_disk, "GS eig_k too small")
 end if

 select case (Wfk%iomode)
 case (IO_MODE_FORTRAN)

   ! Rewind the file to have the correct (k,s) block (if needed)
   call wfk_seek(Wfk,ik_ibz,spin)

   ! Read the first record: npw, nspinor, nband_disk
   read(Wfk%fh, err=10, iomsg=errmsg) npw_read, nspinor_read, nband_read

   if (any([npw_read, nspinor_read, nband_read] /= [npw_disk, nspinor_disk, nband_disk])) then
     write(msg,"(a,6(i0,2x))")"Mismatch between (npw, nspinor, nband) read from WFK and those found in HDR ",&
&      npw_read, nspinor_read, nband_read, npw_disk, nspinor_disk, nband_disk
     MSG_ERROR(msg)
   end if

   ! The second record: (k+G) vectors
   if (present(kg_k)) then
     read(Wfk%fh, err=10, iomsg=errmsg) kg_k(1:3,1:npw_disk)
   else
     read(Wfk%fh, err=10, iomsg=errmsg) ! kg_k(1:3,1:npw_disk)
   end if

   ! The third record: eigenvalues and occupation factors.
   if (Wfk%formeig==0) then

     if (present(eig_k) .or. present(occ_k)) then
       ABI_MALLOC(tmp_eigk, (nband_disk))
       ABI_MALLOC(tmp_occk, (nband_disk))

       read(Wfk%fh, err=10, iomsg=errmsg) tmp_eigk, tmp_occk

       if (present(eig_k)) eig_k = tmp_eigk
       if (present(occ_k)) occ_k = tmp_occk

       ABI_FREE(tmp_eigk)
       ABI_FREE(tmp_occk)

     else
       read(Wfk%fh, err=10, iomsg=errmsg) ! eig_k(1:nband_disk)
     end if

     ! The wave-functions.
     if (present(cg_k)) then
       npwso = npw_disk*nspinor_disk
       my_bcount = 0
       do band=1,nband_disk
         if (bmask(band)) then
           ipw = my_bcount * npwso
           my_bcount = my_bcount + 1
           read(Wfk%fh, err=10, iomsg=errmsg) cg_k(1:2,ipw+1:ipw+npwso)
         else
           read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
         end if
       end do
     else
       do band=1,nband_disk
         read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       end do
     end if

   else if (Wfk%formeig==1) then
     ! Read matrix of size (2*nband_k**2)
     npwso = npw_disk*nspinor_disk
     my_bcount = 0

     do band=1,nband_disk
       base = 2*(band-1)*nband_disk
       if (present(eig_k)) then
         read(Wfk%fh, err=10, iomsg=errmsg) eig_k(base+1:base+2*nband_disk)
       else
         read(Wfk%fh, err=10, iomsg=errmsg) ! eig_k(base+1:base+2*nband_disk)
       end if

       if (bmask(band).and.present(cg_k)) then
         ipw = my_bcount * npwso
         my_bcount = my_bcount + 1
         read(Wfk%fh, err=10, iomsg=errmsg) cg_k(1:2,ipw+1:ipw+npwso)
       else
         read(Wfk%fh, err=10, iomsg=errmsg) ! cg_k(1:2,ipw+1:ipw+npwso)
       end if
     end do

   else
     MSG_ERROR("formeig != [0,1]")
   end if

   ! Reached the end of the (k,s) block. Update f90_fptr
   call wfk_update_f90ptr(wfk, ik_ibz, spin)

#ifdef HAVE_MPI_IO
 case (IO_MODE_MPI)
   if (present(kg_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) + xmpio_bsize_frm

     call mpio_read_kg_k(Wfk%fh,my_offset,npw_disk,sc_mode,kg_k,mpierr)
     ABI_CHECK_MPI(mpierr,"mpio_read_kg_k")
   end if

   ! The third record: eigenvalues and occupation factors.
   if (present(eig_k) .or. present(occ_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + xmpio_bsize_frm
     !
     ! formeig=0 =>  Read both eig and occ in tmp_eigk.
     ! formeig=1 =>  Read (nband_k,nband_k) matrix of complex numbers.
     !
     call mpio_read_eigocc_k(Wfk%fh,my_offset,nband_disk,Wfk%formeig,sc_mode,tmp_eigk,mpierr)
     ABI_CHECK_MPI(mpierr,"mpio_read_eigocc")

     if (Wfk%formeig==0) then
       if (present(eig_k)) eig_k(1:nband_disk) = tmp_eigk(1:nband_disk)
       if (present(occ_k)) occ_k(1:nband_disk) = tmp_eigk(nband_disk+1:)
     else if (Wfk%formeig==1) then
       if (present(eig_k)) eig_k(1:2*nband_disk**2) = tmp_eigk(1:2*nband_disk**2)
     else
       MSG_ERROR("formeig not in [0,1]")
     end if

     ABI_FREE(tmp_eigk)
   end if

   if (present(cg_k)) then
     method = 0

     select case (method)
     case (0)
       ! DATA SIEVING:
       !   read max_nband states in chuncks of nbxblock, then extract my states according to bmask.
       !
       ! MAX number of bands read by the procs in the communicator
       my_maxb = nband_disk
       do band=nband_disk,1,-1
         if (bmask(band)) then
           my_maxb = band
           EXIT
         end if
       end do
       call xmpi_max(my_maxb,max_nband,Wfk%comm,mpierr)
       !max_nband = nband_disk
       !
       ! MPI-IO crashes if we try to read a large number of bands in a single call.
       nbxblock = max_nband
       if ((2*npw_disk*nspinor_disk*nbxblock*xmpi_bsize_dp) > Wfk%chunk_bsize) then
         nbxblock = Wfk%chunk_bsize / (2*npw_disk*nspinor_disk*xmpi_bsize_dp)
         if (nbxblock == 0) nbxblock = 50
       end if
       !nbxblock = 2

       nblocks = max_nband / nbxblock
       brest   = MOD(max_nband, nbxblock)
       if (brest /= 0) nblocks = nblocks + 1
       !write(std_out,*)"in buffered bmask: "nb_tot",nb_tot,"nblocks",nblocks,"nbxblock",nbxblock

       base_ofs = Wfk%offset_ks(ik_ibz,spin,REC_CG)
       sizes = [npw_disk*nspinor_disk, nband_disk]

       my_bcnt = 0  ! index of my band in cg_k
       do block=1,nblocks
         bstart = 1 + (block-1) * nbxblock
         bstop  = bstart + nbxblock - 1
         if (bstop > max_nband) bstop = max_nband
         nb = bstop - bstart + 1
         !
         ! Allocate and read the buffer
         band_block = [bstart, bstop]
         ugsz = npw_disk*nspinor_disk
         bufsz = 2*ugsz*(bstop-bstart+1)
         ABI_MALLOC(buffer, (2,bufsz))
         !write(std_out,*)"  bstart,bstop, ",band_block

         ! Read the cg_ks(G). different versions depending on formeig
         if (wfk%formeig == 0) then
           subsizes = (/npw_disk*nspinor_disk, band_block(2)-band_block(1)+1/)
           starts = [1, bstart]

           call mpiotk_read_fsuba_dp2D(Wfk%fh,base_ofs,sizes,subsizes,starts,&
              bufsz,buffer,Wfk%chunk_bsize,sc_mode,Wfk%comm,ierr)
           ABI_CHECK(ierr==0,"Fortran record too big")

         else if (wfk%formeig == 1) then

           ! Increment my_offset to account for the previous eigen and cg records of (iband-1) bands.
           my_offset = wfk%offset_ks(ik_ibz,spin,REC_CG) + xmpio_bsize_frm
           my_offset = my_offset + (bstart - 1) * ( &
              (2*npw_disk*wfk%nspinor*xmpi_bsize_dp + 2*xmpio_bsize_frm) + &
              (2*nband_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm) )

           sizes = [npw_disk*nspinor_disk, nband_disk]
           types = [MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX]

           call xmpio_create_fstripes(nb,sizes,types,cg_type,my_offpad,mpierr)
           ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

           call MPI_FILE_SET_VIEW(wfk%fh,my_offset,MPI_BYTE,cg_type,'native',xmpio_info,mpierr)
           ABI_CHECK_MPI(mpierr,"")
           call MPI_TYPE_FREE(cg_type,mpierr)
           ABI_CHECK_MPI(mpierr,"")

           if (sc_mode==xmpio_collective) then
             call MPI_FILE_READ_ALL(wfk%fh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
           else if (sc_mode==xmpio_single) then
             call MPI_FILE_READ(wfk%fh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
           else
             MSG_ERROR("Wrong sc_mode")
           end if
           ABI_CHECK_MPI(mpierr,"FILE_READ")
         end if

         ! Extract my bands from buffer.
         do band=bstart,bstop
           if (bmask(band)) then
             my_bcnt = my_bcnt + 1
             pt1 = 1 + (my_bcnt - 1) * ugsz
             pt2 = 1 + (band - bstart) * ugsz
             cg_k(:,pt1:pt1+ugsz-1) = buffer(:,pt2:pt2+ugsz-1)
           end if
         end do

         ABI_FREE(buffer)
       end do

     case (1,2)
       ABI_CHECK(wfk%formeig == 0, "formeig == 1 not coded")
       call MPI_TYPE_CONTIGUOUS(npw_disk*nspinor_disk,MPI_DOUBLE_COMPLEX,cg_type,mpierr)
       ABI_CHECK_MPI(mpierr,"type_contigous")

       if (method==1) then
         ncount = nb_tot
         ABI_MALLOC(block_length,(ncount+2))
         ABI_MALLOC(block_type, (ncount+2))
         ABI_MALLOC(block_displ,(ncount+2))

         block_length(1)=1
         block_displ (1)=0
         block_type  (1)=MPI_LB

         my_bcount = 1
         do band=1,Wfk%mband
           if (bmask(band)) then
             my_bcount = my_bcount + 1
             block_length(my_bcount) = 1
             block_type(my_bcount) = cg_type
             block_displ(my_bcount) = xmpio_bsize_frm + &
&              (band-1) * (2*npw_disk*nspinor_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm)
           end if
         end do

         block_length(ncount+2) = 1
         block_displ (ncount+2) = block_displ(my_bcount)
         block_type  (ncount+2) = MPI_UB

       else if (method==2) then
         ! this file view is not efficient but it's similar to the
         ! one used in wff_readwrite. Let's see if MPI-IO likes it!
         ncount = nb_tot* nspinor_disk * npw_disk

         ABI_MALLOC(block_length,(ncount+2))
         ABI_MALLOC(block_type, (ncount+2))
         ABI_MALLOC(block_displ,(ncount+2))

         block_length(1)=1
         block_displ (1)=0
         block_type  (1)=MPI_LB
         !
         ! The view starts at REC_CG
         cnt = 1
         do band=1,Wfk%mband
           if (bmask(band)) then
             base_ofs =  xmpio_bsize_frm + &
&              (band-1) * (2*npw_disk*nspinor_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm)
             do ipw=1,npw_disk*nspinor_disk
               cnt = cnt + 1
               block_length(cnt) = 1
               block_type(cnt)   = MPI_DOUBLE_COMPLEX
               block_displ(cnt)  = base_ofs + 2*(ipw-1)*xmpi_bsize_dp
             end do
           end if
         end do

         block_length(ncount+2) = 1
         block_displ (ncount+2) = block_displ(cnt)
         block_type  (ncount+2) = MPI_UB
       end if

       call xmpio_type_struct(ncount+2,block_length,block_displ,block_type,cgscatter_type,mpierr)
       ABI_CHECK_MPI(mpierr,"type_struct")

       ABI_FREE(block_length)
       ABI_FREE(block_type)
       ABI_FREE(block_displ)

       call MPI_TYPE_FREE(cg_type, mpierr)
       ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG)

       call MPI_FILE_SET_VIEW(Wfk%fh, my_offset, MPI_BYTE, cgscatter_type, 'native', xmpio_info, mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")

       call MPI_TYPE_FREE(cgscatter_type, mpierr)
       ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

       call MPI_FILE_READ_ALL(Wfk%fh, cg_k, npw_tot, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
       ABI_CHECK_MPI(mpierr,"FILE_READ_ALL")

     case default
       MSG_ERROR("Wrong method")
     end select
   end if
#endif

#ifdef HAVE_NETCDF
 case (IO_MODE_ETSF)
   ABI_CHECK(wfk%formeig == 0, "formeig != 0 not coded")
   !write(std_out,*)"bmask: ",bmask

   ! TODO: extract routines (see other similar calls)
   if (present(kg_k)) then
     ! Read the reduced_coordinates_of_plane_waves for this k point.
     NCF_CHECK(nf90_inq_varid(wfk%fh, "reduced_coordinates_of_plane_waves", kg_varid))
     if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
       NCF_CHECK(nctk_set_collective(wfk%fh, kg_varid))
     end if

     ncerr = nf90_get_var(wfk%fh, kg_varid, kg_k, start=[1,1,ik_ibz], count=[3,npw_disk,1])
     NCF_CHECK(ncerr)
   end if

   ! Read eigenvalues and occupations.
   if (Wfk%formeig==0) then
     if (present(eig_k)) then
       NCF_CHECK(nf90_inq_varid(wfk%fh, "eigenvalues", eig_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, eig_varid))
       end if

       ncerr = nf90_get_var(wfk%fh, eig_varid, eig_k, start=[1,ik_ibz,spin], count=[nband_disk,1,1])
       NCF_CHECK(ncerr)
     end if

     if (present(occ_k)) then
       NCF_CHECK(nf90_inq_varid(wfk%fh, "occupations", occ_varid))
       if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
         NCF_CHECK(nctk_set_collective(wfk%fh, occ_varid))
       end if

       ncerr = nf90_get_var(wfk%fh, occ_varid, occ_k, start=[1,ik_ibz,spin], count=[nband_disk,1,1])
       NCF_CHECK_MSG(ncerr, "getting occ_k")
     end if
   else
     MSG_ERROR("formeig !=0 not compatible with ETSF-IO")
   end if

   if (present(cg_k)) then
     ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
     NCF_CHECK(nf90_inq_varid(wfk%fh, "coefficients_of_wavefunctions", cg_varid))
     ! TODO: Collective
     !if (sc_mode == xmpio_collective .and. wfk%nproc > 1) then
     !  NCF_CHECK(nctk_set_collective(wfk%fh, cg_varid))
     !end if
#if 0
     ! Simple and very inefficient version for debugging.
     ipw = 1
     do band=1,wfk%mband
       if (.not. bmask(band)) cycle
       ncerr = nf90_get_var(wfk%fh, cg_varid, cg_k(:,ipw:), start=[1,1,1,band,ik_ibz,spin], &
         count=[2,npw_disk,wfk%nspinor,1,1,1])
       NCF_CHECK_MSG(ncerr, "getting cg_k block")
       ipw = ipw + wfk%nspinor * npw_disk
     end do
#else
     ! Read bands in blocks defined by bmask.
     ! be careful when in collective mode because processors may call the routine with nblocks==0
     ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
     call mask2blocks(bmask, nblocks, blocks)
     !ABI_CHECK(nblocks /= 0, "nblocks==0")

     ipw = 1
     do block=1,nblocks
       band_block = blocks(:,block)
       nb = band_block(2) - band_block(1) + 1
       ncerr = nf90_get_var(wfk%fh, cg_varid, cg_k(:,ipw:), start=[1,1,1,band_block(1),ik_ibz,spin], &
        count=[2,npw_disk,wfk%nspinor,nb,1,1])
       NCF_CHECK_MSG(ncerr, "getting cg_k block")
       ipw = ipw + wfk%nspinor * npw_disk * nb
     end do
     ABI_FREE(blocks)
#endif

     ! Prototype for collective version.
     !min_band = lfind(bmask); if min_
     !max_band = lfind(bmask, back=.True.) ! n+1
     !! TODO: xmpi_min_max
     !call xmpi_max(my_min_band, min_band, comm_cell, ierr)
     !call xmpi_min(my_min_band, max_band, comm_cell, ierr)
     !nb = max_band - min_band + 1
     !ncalls = nb /

     !NCF_CHECK(nf90_var_par_access(ncid, cg_varid, nf90_collective))
     !do block=1,nblocks
     !  ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
     !  ncerr = nf90_get_var(wfk%fh, cg_varid, cg_k, start=[1,1,1,band_block(1),ik_ibz,spin], &
     !   count=[2,npw_disk,wfk%nspinor,band_block(2)-band_block(1)+1,1,1])
     !  NCF_CHECK_MSG(ncerr, "getting cg_k block")
     !   do band=band_start,band_end
     !     if (.not. bmask(band)) cycle
     !     cg_k(:,:) =
     !   end do
     !end do
   end if
#endif

 case default
   MSG_ERROR(sjoin('Wrong/unsupported value of iomode: ', itoa(wfk%iomode)))
 end select

 DBG_EXIT("COLL")

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine wfk_read_bmask
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_ebands
!! NAME
!!  wfk_read_ebands
!!
!! FUNCTION
!!  Read the GS eigenvalues and return ebands_t object.
!!
!! INPUTS
!!  path=WFK file name
!!  comm=MPI communicator
!!
!! OUTPUTS
!!  ebands<ebands_t>=GS band-structure.
!!  [out_hdr]=Abinit header.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function wfk_read_ebands(path, comm, out_hdr) result(ebands)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm
 type(hdr_type),optional,intent(out) :: out_hdr

!Local variables-------------------------------
!scalars
 type(hdr_type) :: hdr
!arrays
 real(dp),pointer :: eigen(:,:,:)

!************************************************************************

 call wfk_read_eigenvalues(path, eigen, hdr, comm)
 ebands = ebands_from_hdr(hdr,maxval(hdr%nband),eigen)
 if (present(out_hdr)) call hdr_copy(hdr, out_hdr)

 ABI_FREE(eigen)
 call hdr%free()

end function wfk_read_ebands
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_eigk
!! NAME
!!  wfk_read_eigk
!!
!! FUNCTION
!!  Helper function to read all the eigenvalues for a given (k-point,spin)
!!
!! INPUTS
!!  Wfk<class(wfk_t)>= WFK file handler
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  eig_k(1:nband_k) = GS Eigenvalues for the given (k,s)
!!  occ_k(1:nband_k) = Occupation factors for the given (k,s)
!!
!! NOTES
!!  The buffers eig_k and occ_k are dimensions with wfk%mband. The routine
!!  will fill the first nband_k positions with data read from file where
!!  nband_k is the number of bands on file i.e. wfk%nband(ik_ibz,spin)
!!
!! PARENTS
!!      conducti_nc,m_wfk,mrggkk,optic
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_read_eigk(Wfk,ik_ibz,spin,sc_mode,eig_k,occ_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode
 class(wfk_t),intent(inout) :: Wfk
!arrays
 real(dp),intent(out) :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),optional,intent(out) :: occ_k(Wfk%mband)

!Local variables-------------------------------
!scalars
 integer,parameter :: band_block00(2) = [0, 0]

!************************************************************************

 if (present(occ_k)) then
   ABI_CHECK(Wfk%formeig == 0, "formeig !=0")
   call wfk%read_band_block(band_block00,ik_ibz,spin,sc_mode,eig_k=eig_k,occ_k=occ_k)
 else
   call wfk%read_band_block(band_block00,ik_ibz,spin,sc_mode,eig_k=eig_k)
 end if

end subroutine wfk_read_eigk
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_eigenvalues
!! NAME
!!  wfk_read_eigenvalues
!!
!! FUNCTION
!!  Read all the GS eigenvalues stored in the WFK file fname.
!!
!! INPUTS
!!  fname=Name of the file
!!  comm=MPI communicator.
!!
!! OUTPUTS
!!  eigen = In input: nullified pointer
!!          In output: eigen(mband,nkpt,nsppol) contains the GS eigevalues.
!!  Hdr_out<hdr_type>=The header of the file
!!
!! PARENTS
!!      dfpt_looppert,eph,m_wfk,setup_bse,setup_bse_interp,setup_screening
!!      setup_sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_read_eigenvalues(fname, eigen, Hdr_out, comm, occ)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=*),intent(in) :: fname
 type(hdr_type),intent(out) :: Hdr_out
!arrays
!TODO: Replace pointers with allocatable.
 real(dp),pointer :: eigen(:,:,:)
 real(dp),pointer,optional :: occ(:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,formeig0=0
 integer :: ik_ibz,spin,my_rank,ierr,iomode,funt,sc_mode,mband
 real(dp) :: cpu, wall, gflops
 type(wfk_t) :: Wfk

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 my_rank = xmpi_comm_rank(comm)
 iomode = iomode_from_fname(fname)
 !iomode = IO_MODE_FORTRAN

 call wrtout(std_out, sjoin(" Reading eigenvalues from:", fname, ", with iomode:", iomode2str(iomode)))

 if (my_rank == master) then
   ! Open the file.
   sc_mode = xmpio_single
   funt = get_unit()
   call wfk_open_read(Wfk,fname,formeig0,iomode,funt,xmpi_comm_self,Hdr_out=Hdr_out)

   ! Read the eigenvalues and optionally the occupation factors.
   ABI_MALLOC(eigen, (Wfk%mband,Wfk%nkpt,Wfk%nsppol))
   eigen = HUGE(zero)
   if (present(occ)) then
     ABI_MALLOC(occ, (Wfk%mband,Wfk%nkpt,Wfk%nsppol))
     occ = HUGE(zero)
   end if

   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       if (present(occ)) then
         call wfk%read_eigk(ik_ibz,spin,sc_mode,eigen(:,ik_ibz,spin),occ_k=occ(:,ik_ibz,spin))
       else
         call wfk%read_eigk(ik_ibz,spin,sc_mode,eigen(:,ik_ibz,spin))
       end if
     end do
   end do

   ! Close the file.
   call wfk%close()
 end if

 ! Broadcast data
 if (xmpi_comm_size(comm) > 1) then
   call Hdr_out%bcast(master, my_rank, comm)
   mband = MAXVAL(Hdr_out%nband)
   if (my_rank/=master) then
     ABI_MALLOC(eigen, (mband,Hdr_out%nkpt,Hdr_out%nsppol))
     if (present(occ)) then
       ABI_MALLOC(occ, (mband,Hdr_out%nkpt,Hdr_out%nsppol))
     end if
   end if
   call xmpi_bcast(eigen,master,comm,ierr)
   if (present(occ)) call xmpi_bcast(occ,master,comm,ierr)
 end if

 call cwtime_report(" wfk_read_eigenvalues", cpu, wall, gflops)

end subroutine wfk_read_eigenvalues
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_write_h1mat
!! NAME
!!  wfk_write_h1mat
!!
!! FUNCTION
!!  Write all H1 matrix elements in the WFK file fname.
!!
!! INPUTS
!!
!! OUTPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_write_h1mat(Wfk,sc_mode,eigen)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sc_mode
 class(wfk_t),intent(inout) :: Wfk
!arrays
 real(dp),intent(in) :: eigen(2*Wfk%mband**2*Wfk%nkpt*Wfk%nsppol)

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,nband_k,ptr
!arrays
 integer,parameter :: band_block00(2)=[0,0]

!************************************************************************

 ptr=1
 do spin=1,Wfk%nsppol
   do ik_ibz=1,Wfk%nkpt
     nband_k = Wfk%nband(ik_ibz,spin)
     call wfk%write_band_block(band_block00,ik_ibz,spin,sc_mode,eig_k=eigen(ptr:))
     ptr = ptr + 2*nband_k**2
   end do
 end do

end subroutine wfk_write_h1mat
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_h1mat
!! NAME
!!  wfk_read_h1mat
!!
!! FUNCTION
!!  Read all H1 matrix elements in the WFK file fname inside the MPI communicator comm.
!!
!! INPUTS
!!  path=File name
!!  comm=MPI communicator.
!!
!! OUTPUTS
!!  eigen(2*hdr_out%mband**2*hdr_out%nkpt*hdr_out%nsppol)=Array with the matrix elements of H1
!!   packed in the first positions. The array is allocated by the procedure.
!!
!!  Hdr_out<hdr_type>=The header of the file
!!
!! PARENTS
!!      m_ddk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_read_h1mat(fname, eigen, hdr_out, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
 integer,intent(in) :: comm
 type(hdr_type),intent(out) :: Hdr_out
!arrays
 real(dp),allocatable,intent(out) :: eigen(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,formeig1=1
 integer :: spin,ik_ibz,nband_k,ptr,ierr,iomode,mband,my_rank
 type(wfk_t) :: wfk
!arrays
 integer,parameter :: band_block00(2)=[0,0]

!************************************************************************

 my_rank = xmpi_comm_rank(comm)

 if (my_rank==master) then
   ! Open the file.
   iomode = iomode_from_fname(fname)
   call wfk_open_read(wfk, fname, formeig1, iomode, get_unit(), xmpi_comm_self, hdr_out=hdr_out)

   ! Read h1 mat and pack them in the first positions.
   ABI_MALLOC(eigen, (2*wfk%mband**2*wfk%nkpt*wfk%nsppol))

   ptr=1
   do spin=1,wfk%nsppol
     do ik_ibz=1,wfk%nkpt
       nband_k = wfk%nband(ik_ibz,spin)
       call wfk_read_band_block(wfk, band_block00, ik_ibz, spin, xmpio_single, eig_k=eigen(ptr:))
       ptr = ptr + 2*nband_k**2
     end do
   end do

   call wfk%close()
 end if

 ! Broadcast data
 if (xmpi_comm_size(comm) > 1) then
   call hdr_out%bcast(master, my_rank, comm)

   mband = maxval(Hdr_out%nband)
   if (my_rank/=master) then
     ABI_MALLOC(eigen, (2*mband**2*hdr_out%nkpt*hdr_out%nsppol))
   end if
   call xmpi_bcast(eigen,master,comm,ierr)
 end if

end subroutine wfk_read_h1mat
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_rewind
!! NAME
!!  wfk_rewind
!!
!! FUNCTION
!!  Rewind the file, skip the header and modifies Wfk%f90_fptr $
!!  Mainly used for debugging purposes when IO_MODE_FORTRAN is used.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_rewind(wfk)

!Arguments ------------------------------------
 class(wfk_t),intent(inout) :: wfk

!Local variables-------------------------------
 integer :: ierr

! *************************************************************************

 select case (wfk%iomode)
 case (IO_MODE_FORTRAN)
   rewind(wfk%fh)
   call hdr_skip(wfk%fh,ierr)
   ABI_CHECK(ierr==0, "hdr_skip returned ierr! /= 0")
   wfk%f90_fptr = [1,1,REC_NPW]

 case default
   MSG_ERROR("should not be called when wfk%iomode /= IO_MODE_FORTRAN")
 end select

end subroutine wfk_rewind
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_seek
!! NAME
!!  wfk_seek
!!
!! FUNCTION
!!   Move the internal file pointer so that it points to the
!!   block (ik_ibz, spin). Needed only if iomode==IO_MODE_FORTRAN
!!
!! INPUTS
!!   ik_ibz,spin = (k-point,spin) indices
!!
!! SIDE EFFECTS
!!   Wfk<class(wfk_t)> : modifies Wfk%f90_fptr and the internal F90 file pointer.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_seek(Wfk,ik_ibz,spin)

!Arguments ------------------------------------
 integer,intent(in)  :: ik_ibz,spin
 class(wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
 integer :: ierr,ik_fpt,spin_fpt,recn_wanted,recn_fpt,rec_type
 character(len=500) :: msg

! *************************************************************************

 select case (Wfk%iomode)
 case (IO_MODE_FORTRAN)
   !
   ! Find the position inside the file.
   if (ALL(Wfk%f90_fptr==FPTR_EOF)) then ! handle the EOF condition
     if (Wfk%debug) call wrtout(std_out,"EOF condition","PERS")
     recn_fpt = Wfk%recn_eof
   else
     ik_fpt   = Wfk%f90_fptr(1)
     spin_fpt = Wfk%f90_fptr(2)
     rec_type = Wfk%f90_fptr(3)
     recn_fpt = Wfk%recn_ks(ik_fpt,spin_fpt, rec_type)
   end if

   recn_wanted = Wfk%recn_ks(ik_ibz,spin, REC_NPW)

   if (Wfk%debug) then
     write(msg,'(a,3(i0,2x))')"seeking ik_ibz, spin, recn_wanted-recn_fpt: ",ik_ibz,spin,recn_wanted - recn_fpt
     call wrtout(std_out,msg,"PERS")
   end if

   call mvrecord(Wfk%fh, (recn_wanted - recn_fpt) ,ierr)
   ABI_CHECK(ierr==0,"error in mvrecord")

   Wfk%f90_fptr = [ik_ibz, spin, REC_NPW]

 case default
   MSG_ERROR("should not be called when Wfk%iomode /= IO_MODE_FORTRAN")
 end select

end subroutine wfk_seek
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_update_f90ptr
!! NAME
!!  wfk_update_f90ptr
!!
!! FUNCTION
!!  Update wfk%f90_ptr. Used if wfk%iomode == IO_MODE_FORTRAN.
!!
!! INPUTS
!!  ik_ibz=K-point index,
!!  spin=Spin index.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_update_f90ptr(wfk, ik_ibz, spin)

!Arguments ------------------------------------
 class(wfk_t),intent(inout) :: wfk
 integer,intent(in) :: ik_ibz,spin

! *************************************************************************

 if (ik_ibz < wfk%nkpt) then
   wfk%f90_fptr = [ik_ibz+1,spin,REC_NPW]
 else
   ABI_CHECK(ik_ibz == wfk%nkpt, "ik_ibz != nkpt")
   if (spin == wfk%nsppol) then
     wfk%f90_fptr = FPTR_EOF ! EOF condition
   else
     wfk%f90_fptr = [1,spin+1,REC_NPW]
   end if
 end if

end subroutine wfk_update_f90ptr
!!***

!----------------------------------------------------------------------

!!****f* m_wfkfile/wfk_compute_offsets
!! NAME
!!  wfk_compute_offsets
!!
!! FUNCTION
!!  Compute the offsets corresponding to the different sections of the file (G-vectors, eigenvalues, u(G).
!!  Needed only for Fortran-IO or MPI-IO.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_compute_offsets(Wfk)

!Arguments ------------------------------------
 class(wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,npw_k,nband_k,bsize_frm,mpi_type_frm,base !,band
 integer(XMPI_OFFSET_KIND) :: offset

! *************************************************************************

 select case (Wfk%iomode)
 case (IO_MODE_FORTRAN)
   !
   ! Compute record number for Fortran IO
   ABI_MALLOC(Wfk%recn_ks,(Wfk%nkpt,Wfk%nsppol,REC_NUM))

   ! We start to count the number of Fortran records from the end of the Header
   ! Hence recn gives the relative position from the header, it's not an absolute position.
   base = 0
   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       nband_k = Wfk%nband(ik_ibz,spin)
       Wfk%recn_ks(ik_ibz,spin, REC_NPW) = base + 1
       Wfk%recn_ks(ik_ibz,spin, REC_KG)  = base + 2
       Wfk%recn_ks(ik_ibz,spin, REC_EIG) = base + 3
       Wfk%recn_ks(ik_ibz,spin, REC_CG)  = base + 4
       base = Wfk%recn_ks(ik_ibz,spin,REC_CG)
       if (Wfk%formeig==0) then
         base = base + (nband_k-1)
       else if (Wfk%formeig==1) then
         base = base + 2*(nband_k-1)
       else
         MSG_ERROR("formeig != [0,1]")
       end if
     end do
   end do
   !
   ! Save EOF position
   Wfk%recn_eof = base + 1

 case (IO_MODE_MPI)
   !
   ! Compute offsets for MPI-IO.
   ABI_MALLOC(Wfk%offset_ks,(Wfk%nkpt,Wfk%nsppol,REC_NUM))

   bsize_frm    = xmpio_bsize_frm    ! Byte length of the Fortran record marker.
   mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

   ! The offset of the Header. TODO
   ! hdr_offset(Hdr)
   offset = Wfk%hdr_offset

   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       npw_k   = Wfk%Hdr%npwarr(ik_ibz)
       nband_k = Wfk%nband(ik_ibz,spin)
       !
       !---------------------------------------------------------------------------
       ! First record: npw, nspinor, nband_disk
       !---------------------------------------------------------------------------
       Wfk%offset_ks(ik_ibz,spin,REC_NPW) = offset

       if (Wfk%Hdr%headform>=40) then
         ! npw, nspinor, nband_disk
         offset = offset +  3*xmpi_bsize_int + 2*bsize_frm
       else
         MSG_ERROR("Old headforms < 40 are not supported")
       end if
       Wfk%offset_ks(ik_ibz,spin,REC_KG) = offset
       !
       !---------------------------------------------------------------------------
       ! Second record: (k+G) vectors
       ! kg_k(1:3,1:npw_k)
       !---------------------------------------------------------------------------
       offset = offset + 3*npw_k*xmpi_bsize_int + 2*bsize_frm
       Wfk%offset_ks(ik_ibz,spin,REC_EIG) = offset
       !
       !---------------------------------------------------------------------------
       ! Third record: eigenvalues
       !---------------------------------------------------------------------------
       if (Wfk%formeig==0) then
         ! eigen(1:nband_k), occ(1:nband_k)
         offset = offset + 2*nband_k*xmpi_bsize_dp + 2*bsize_frm
         Wfk%offset_ks(ik_ibz,spin,REC_CG) = offset
         !
         ! Wavefunction coefficients
         ! do band=1,nband_k; write(unitwf) cg_k(1:2,npw_k*nspinor); end do
         offset = offset + nband_k * (2*npw_k*Wfk%nspinor*xmpi_bsize_dp + 2*bsize_frm)

       else if (Wfk%formeig==1) then
         ! read(unitwf) eigen(2*nband_k)
         Wfk%offset_ks(ik_ibz,spin,REC_CG) = offset + 2*nband_k*xmpi_bsize_dp + 2*bsize_frm

         offset = offset + &
&          nband_k * (2*npw_k*Wfk%nspinor*xmpi_bsize_dp + 2*bsize_frm) + &
&          nband_k * (2*nband_k*xmpi_bsize_dp + 2*bsize_frm)

       else
         MSG_ERROR("Wrong formeig")
       end if

     end do ! ik_ibz
   end do ! spin
   !
   ! Save EOF offset
   Wfk%offset_eof = offset

   ! Check for possible wraparound errors.
   if (ANY(Wfk%offset_ks <= 0) .or. Wfk%offset_eof < 0) then
     MSG_ERROR("Found negative offset. File too large for MPI-IO!!!")
   end if
 end select

 if (Wfk%debug) call wfk_show_offsets(Wfk)

end subroutine wfk_compute_offsets
!!***

!----------------------------------------------------------------------

!!****f* m_wfkfile/wfk_show_offsets
!! NAME
!!  wfk_show_offsets
!!
!! FUNCTION
!!  Print the offsets.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_show_offsets(Wfk)

!Arguments ------------------------------------
 class(wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz

! *************************************************************************

 select case (Wfk%iomode)

 case (IO_MODE_FORTRAN)
   write(std_out,*)"Record number relative to the header."
   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       write(std_out,"(a,2(i0,2x),a,4(a,i0,a))")                   &
&        "(ik_ibz, spin) ",ik_ibz,spin,ch10,                       &
&        "  recn(REC_NPW): ",Wfk%recn_ks(ik_ibz,spin,REC_NPW),ch10,&
&        "  recn(REC_KG) : ",Wfk%recn_ks(ik_ibz,spin,REC_KG), ch10,&
&        "  recn(REC_EIG): ",Wfk%recn_ks(ik_ibz,spin,REC_EIG),ch10,&
&        "  recn(REC_CG) : ",Wfk%recn_ks(ik_ibz,spin,REC_CG),ch10
     end do
   end do

   write(std_out,"(a,i0)")"EOS position: ",Wfk%recn_eof

 case (IO_MODE_MPI)
   write(std_out,"(a,i0)")"hdr_offset ",Wfk%hdr_offset

   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       write(std_out,"(a,2(i0,2x),a,4(a,i0,a))")                       &
&        "(ik_ibz, spin) ",ik_ibz,spin,ch10,                           &
&        "  offset(REC_NPW): ",Wfk%offset_ks(ik_ibz,spin,REC_NPW),ch10,&
&        "  offset(REC_KG) : ",Wfk%offset_ks(ik_ibz,spin,REC_KG), ch10,&
&        "  offset(REC_EIG): ",Wfk%offset_ks(ik_ibz,spin,REC_EIG),ch10,&
&        "  offset(REC_CG) : ",Wfk%offset_ks(ik_ibz,spin,REC_CG),ch10
     end do ! ik_ibz
   end do ! spin
   !
   ! Write EOF position
   write(std_out,"(a,i0)")"offset_eof ",Wfk%offset_eof
 end select

end subroutine wfk_show_offsets
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_read_kg_k
!! NAME
!!  mpio_read_kg_k
!!
!! FUNCTION
!!  Helper functions to read the G-vectors with MPI-IO
!!
!! INPUTS
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  kg_k=(3,npw_disk) = G-vectors
!!  mpierr=MPI error status (error check is delegated to the caller)
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_read_kg_k(fh,offset,npw_disk,sc_mode,kg_k,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,npw_disk,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 integer,intent(out) :: kg_k(3,npw_disk)

!Local variables-------------------------------
!scalars
 integer :: kg_k_type,ncount,myfh
 integer(XMPI_OFFSET_KIND) :: my_offset

!************************************************************************

 ! Workarounds for XLF
 myfh      = fh
 ncount    = 3*npw_disk
 my_offset = offset

 call MPI_TYPE_CONTIGUOUS(ncount, MPI_INTEGER, kg_k_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_COMMIT(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,kg_k_type,'native',xmpio_info,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 if (sc_mode==xmpio_collective) then
   call MPI_FILE_READ_ALL(myfh,kg_k,ncount,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
 else if (sc_mode==xmpio_single) then
   !call MPI_File_seek(myfh, 0, MPI_SEEK_SET,mpierr)
   call MPI_FILE_READ(myfh,kg_k,ncount,MPI_INTEGER, MPI_STATUS_IGNORE,mpierr)
 else
   MSG_ERROR("Wrong sc_mode")
 end if

 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

end subroutine mpio_read_kg_k
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_write_kg_k
!! NAME
!!  mpio_write_kg_k
!!
!! FUNCTION
!!  Helper function to write the G-vectors with MPI-IO
!!
!! INPUTS
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for writing by current proc.
!!    xmpio_collective ==> for collective write.
!!  kg_k=(3,npw_disk) = G-vectors
!!
!! OUTPUTS
!!  mpierr=MPI error status (error check is delegated to the caller)
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_write_kg_k(fh,offset,npw_disk,sc_mode,kg_k,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,npw_disk,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 integer,intent(in) :: kg_k(3,npw_disk)

!Local variables-------------------------------
!scalars
 integer :: myfh,kg_k_type,ncount
 integer(XMPI_OFFSET_KIND) :: my_offset

!************************************************************************

 DBG_ENTER("COLL")

 ! Workarounds for XLF
 myfh      = fh
 ncount    = 3*npw_disk
 my_offset = offset

 call MPI_TYPE_CONTIGUOUS(ncount, MPI_INTEGER, kg_k_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_COMMIT(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,kg_k_type,'native',xmpio_info,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 if (sc_mode==xmpio_collective) then
   call MPI_FILE_WRITE_ALL(myfh,kg_k,ncount,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
 else if (sc_mode==xmpio_single) then
   call MPI_FILE_WRITE(myfh,kg_k,ncount,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
 else
   MSG_ERROR("Wrong sc_mode")
 end if

 ABI_HANDLE_MPIERR(mpierr)

 DBG_EXIT("COLL")

end subroutine mpio_write_kg_k
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_read_eigocc_k
!! NAME
!!  mpio_read_eigocc_k
!!
!! FUNCTION
!!  Helper functions to read the eigenvalues and the occupations with MPI-IO
!!
!! INPUTS
!!  fh
!!  offset
!!  nband_disk
!!  formeig
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  buffer(:)
!!  mpierr=MPI error status.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_read_eigocc_k(fh,offset,nband_disk,formeig,sc_mode,buffer,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,nband_disk,formeig,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 real(dp),pointer :: buffer(:)

!Local variables-------------------------------
!scalars
 integer :: myfh,bufsz,gkk_type,eneocc_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad !,fmarker
!arrays
 integer :: sizes(2),subsizes(2),starts(2)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 SELECT CASE (formeig)
 CASE (0)
   !
   ! Read both eig and occ in buffer
   bufsz = 2*nband_disk
   my_offset = offset
   ABI_MALLOC(buffer, (bufsz))

   call MPI_TYPE_CONTIGUOUS(bufsz, MPI_DOUBLE_PRECISION, eneocc_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_COMMIT(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_FILE_SET_VIEW(myfh, my_offset, MPI_BYTE, eneocc_type, 'native', xmpio_info, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_READ(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE (1)
   ! Read the (nband_k,nband_k) matrix with the (complex) GKK matrix elements.
   bufsz    = (nband_disk**2)
   sizes    = [nband_disk, nband_disk]
   subsizes = [nband_disk, nband_disk]
   starts   = [1, 1]

   ABI_MALLOC(buffer, (2*bufsz))

   !my_offset = offset - xmpio_bsize_frm
   !call xmpio_read_dp(myfh,my_offset,sc_mode,2*nband_disk,buffer,fmarker,mpierr)
   !write(std_out,*)buffer(1:2*nband_disk)
   !MSG_ERROR("Done")

   call xmpio_create_fsubarray_2D(sizes,subsizes,starts,MPI_DOUBLE_COMPLEX,gkk_type,my_offpad,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   ! TODO: Rationalize the offsets
   my_offset = offset + my_offpad - xmpio_bsize_frm

   call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,gkk_type,'native',xmpio_info,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(gkk_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_READ(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE DEFAULT
   MSG_ERROR("formeig not in [0,1]")
 END SELECT

end subroutine mpio_read_eigocc_k
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_write_eigocc_k
!! NAME
!!  mpio_write_eigocc_k
!!
!! FUNCTION
!!  Helper functions to write the eigenvalues and the occupations with MPI-IO
!!
!! INPUTS
!!  fh
!!  offset
!!  nband_disk
!!  formeig
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for writing  by current proc.
!!    xmpio_collective ==> for collective write.
!!
!! OUTPUTS
!!  buffer(:)
!!  mpierr=MPI error status.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_write_eigocc_k(fh,offset,nband_disk,formeig,sc_mode,buffer,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,nband_disk,formeig,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 real(dp),intent(in) :: buffer(:)

!Local variables-------------------------------
!scalars
 integer :: bufsz,gkk_type,eneocc_type,myfh
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
!arrays
 integer :: sizes(2),subsizes(2),starts(2)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 SELECT CASE (formeig)
 CASE (0)
   !
   ! write both eig and occ in buffer
   my_offset = offset

   bufsz = 2*nband_disk
   ABI_CHECK(SIZE(buffer) >= bufsz, "buffer too small")

   call MPI_TYPE_CONTIGUOUS(bufsz, MPI_DOUBLE_PRECISION, eneocc_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_COMMIT(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,eneocc_type,'native',xmpio_info,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE (1)
   !MSG_ERROR("formeig ==1 with MPI-IO not tested")
   ! write the (nband_k,nband_k) matrix with the (complex) GKK matrix elements.
   bufsz    = (nband_disk**2)
   sizes    = [nband_disk, nband_disk]
   subsizes = [nband_disk, nband_disk]
   starts   = [1, 1]

   ABI_CHECK(SIZE(buffer) >= bufsz, "buffer too small")

   call xmpio_create_fsubarray_2D(sizes,subsizes,starts,MPI_DOUBLE_COMPLEX,gkk_type,my_offpad,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   ! TODO: Rationalize the offsets
   my_offset = offset + my_offpad - xmpio_bsize_frm

   call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,gkk_type,'native',xmpio_info,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(gkk_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE DEFAULT
   MSG_ERROR("formeig not in [0,1]")
 END SELECT

end subroutine mpio_write_eigocc_k
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_tofullbz
!! NAME
!!  wfk_tofullbz
!!
!! FUNCTION
!! Generate a new WFK file with wavefunctions in the full BZ and istwfk==1
!! Mainly used to interface ABINIT with other codes that cannot handle symmetries e.g. lobster
!!
!! INPUTS
!!  in_path = Input WFK file
!!  dtset <dataset_type>=all input variables for this dataset
!!  psps <pseudopotential_type>=all the information about psps
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  out_path = Output WFK file.
!!
!! OUTPUT
!!  Output is written to file out_path
!!
!! NOTES
!!  - This routine should be called by a single processor.
!!  - Only GS WFK files are supported (formeig==0)
!!
!! PARENTS
!!      gstate,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_tofullbz(in_path, dtset, psps, pawtab, out_path)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: in_path,out_path
 type(pseudopotential_type),intent(in) :: psps
 type(dataset_type),intent(in) :: dtset
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0=0,kptopt3=3
 integer :: spin,ikf,ik_ibz,nband_k,mpw_ki,mpw_kf,mband,nspinor,nkfull
 integer :: in_iomode,nsppol,nkibz,out_iomode,isym,itimrev
 integer :: npw_ki,npw_kf,istwf_ki,istwf_kf,ii,jj,iqst,nqst
 real(dp) :: ecut_eff,dksqmax,cpu,wall,gflops
 character(len=500) :: msg
 character(len=fnlen) :: my_inpath
 logical :: isirred_kf
 logical,parameter :: force_istwfk1=.True.
 type(wfk_t),target :: iwfk
 type(wfk_t) :: owfk
 type(crystal_t) :: cryst
 type(hdr_type) :: hdr_kfull
 type(hdr_type),pointer :: ihdr
 type(ebands_t) :: ebands_ibz
 type(ebands_t),target :: ebands_full
 type(wvl_internal_type) :: dummy_wvl
!arrays
 integer :: g0(3),work_ngfft(18),gmax_ki(3),gmax_kf(3),gmax(3)
 integer,allocatable :: bz2ibz(:,:),kg_ki(:,:),kg_kf(:,:),iperm(:),bz2ibz_sort(:)
 real(dp) :: kf(3),kibz(3)
 real(dp),allocatable :: cg_ki(:,:),cg_kf(:,:),eig_ki(:),occ_ki(:),work(:,:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: kfull(:,:)

! *************************************************************************

 if (all(dtset%kptrlatt == 0)) then
   write(msg,"(5a)")&
     "Cannot produce full WFK file because kptrlatt == 0",ch10,&
     "Please use nkgpt and shiftk to define a homogeneous k-mesh.",ch10,&
     "Returning to caller"
   MSG_WARNING(msg)
   return
 end if

 call cwtime(cpu, wall, gflops, "start")
 my_inpath = in_path

 if (nctk_try_fort_or_ncfile(my_inpath, msg) /= 0) then
   MSG_ERROR(msg)
 end if
 call wrtout(std_out, sjoin("Converting:", my_inpath, "to", out_path))

 in_iomode = iomode_from_fname(my_inpath)

 ebands_ibz = wfk_read_ebands(my_inpath, xmpi_comm_self)

 ! Open input file, extract dimensions and allocate workspace arrays.
 call wfk_open_read(iwfk,my_inpath,formeig0,in_iomode,get_unit(),xmpi_comm_self)
 ihdr => iwfk%hdr

 mband = iwfk%mband; mpw_ki = maxval(iwfk%Hdr%npwarr); nkibz = iwfk%nkpt
 nsppol = iwfk%nsppol; nspinor = iwfk%nspinor
 ecut_eff = iwfk%hdr%ecut_eff ! ecut * dilatmx**2

 ABI_MALLOC(kg_ki, (3, mpw_ki))
 ABI_MALLOC(cg_ki, (2, mpw_ki*nspinor*mband))
 ABI_MALLOC(eig_ki, ((2*mband)**iwfk%formeig*mband) )
 ABI_MALLOC(occ_ki, (mband))

 cryst = iwfk%hdr%get_crystal(2)

 ! Build new header for owfk. This is the most delicate part since all the arrays in hdr_full
 ! that depend on k-points must be consistent with kfull and nkfull.
 call ebands_expandk(ebands_ibz, cryst, ecut_eff, force_istwfk1, dksqmax, bz2ibz, ebands_full)

 if (dksqmax > tol12) then
   write(msg, '(3a,es16.6,4a)' )&
   'At least one of the k points could not be generated from a symmetrical one.',ch10,&
   'dksqmax=',dksqmax,ch10,&
   'Action: check your WFK file and k-point input variables',ch10,&
   '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   MSG_ERROR(msg)
 end if

 nkfull = ebands_full%nkpt
 kfull => ebands_full%kptns

 ! Build new header and update pawrhoij.
 call hdr_init_lowlvl(hdr_kfull,ebands_full,psps,pawtab,dummy_wvl,abinit_version,&
   ihdr%pertcase,ihdr%natom,ihdr%nsym,ihdr%nspden,ihdr%ecut,dtset%pawecutdg,ihdr%ecutsm,dtset%dilatmx,&
   ihdr%intxc,ihdr%ixc,ihdr%stmbias,ihdr%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,ihdr%so_psp,&
   ihdr%qptn,cryst%rprimd,cryst%xred,ihdr%symrel,ihdr%tnons,ihdr%symafm,ihdr%typat,ihdr%amu,ihdr%icoulomb,&
   kptopt3,dtset%nelect,dtset%charge,dtset%kptrlatt_orig,dtset%kptrlatt,&
   dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk)

 if (psps%usepaw == 1) call pawrhoij_copy(iwfk%hdr%pawrhoij, hdr_kfull%pawrhoij)

 out_iomode = iomode_from_fname(out_path)
 call owfk%open_write(hdr_kfull, out_path, iwfk%formeig, out_iomode, get_unit(), xmpi_comm_self)
 call hdr_kfull%free()

 ! workspace array for BZ wavefunction block.
 mpw_kf = maxval(ebands_full%npwarr)
 ABI_MALLOC(cg_kf, (2,mpw_kf*nspinor*mband))

 if (out_iomode == IO_MODE_FORTRAN) then
   call wrtout(std_out,"Using (slow) Fortran IO version to generate full WFK file", do_flush=.True.)

   ! Fortran IO does not support random access hence the external loop is on the k-points in the full BZ.
   !
   !   For each point in the BZ:
   !     - Find symmetric k-point in the IBZ and read IBZ wavefunctions from iwfk
   !     - Rotate wavefunctions in G-space and write kbz data
   !
   ! Inefficient since we are reading the same IBZ block several times.
   do spin=1,nsppol
     do ikf=1,nkfull
       ik_ibz = bz2ibz(ikf,1); isym = bz2ibz(ikf,2); itimrev = bz2ibz(ikf,6); g0 = bz2ibz(ikf,3:5) ! IS(k_ibz) + g0 = k_bz
       isirred_kf = (isym == 1 .and. itimrev == 0 .and. all(g0 == 0))

       nband_k = iwfk%nband(ik_ibz,spin)
       kf = kfull(:,ikf)
       kibz = ebands_ibz%kptns(:,ik_ibz)

       istwf_ki = iwfk%hdr%istwfk(ik_ibz)
       istwf_kf = owfk%hdr%istwfk(ikf)
       npw_ki = iwfk%hdr%npwarr(ik_ibz)

       ! Read IBZ data.
       call iwfk%read_band_block([1,nband_k],ik_ibz,spin,xmpio_single,&
         kg_k=kg_ki,cg_k=cg_ki,eig_k=eig_ki,occ_k=occ_ki)

       ! The test on npwarr is needed because we may change istwfk e.g. gamma.
       if (isirred_kf .and. iwfk%hdr%npwarr(ik_ibz) == owfk%hdr%npwarr(ikf)) then

         call owfk%write_band_block([1,nband_k],ikf,spin,xmpio_single,&
           kg_k=kg_ki,cg_k=cg_ki,eig_k=eig_ki,occ_k=occ_ki)

       else
         ! Compute G-sphere centered on kf
         call get_kg(kf,istwf_kf,ecut_eff,cryst%gmet,npw_kf,kg_kf)
         ABI_CHECK(npw_kf == owfk%hdr%npwarr(ikf), "Wrong npw_kf")

         ! FFT box must enclose the two spheres centered on ki and kf
         gmax_ki = maxval(abs(kg_ki(:,1:npw_ki)), dim=2)
         gmax_kf = maxval(abs(kg_kf), dim=2)
         do ii=1,3
           gmax(ii) = max(gmax_ki(ii), gmax_kf(ii))
         end do
         gmax = 2*gmax + 1
         call ngfft_seq(work_ngfft, gmax)
         ABI_CALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))

         ! Rotate nband_k wavefunctions (output in cg_kf)
         call cgtk_rotate(cryst,kibz,isym,itimrev,g0,nspinor,nband_k,&
           npw_ki,kg_ki,npw_kf,kg_kf,istwf_ki,istwf_kf,cg_ki,cg_kf,work_ngfft,work)

         ABI_FREE(work)

         ! Write data
         call owfk%write_band_block([1,nband_k],ikf,spin,xmpio_single,&
           kg_k=kg_kf,cg_k=cg_kf,eig_k=eig_ki,occ_k=occ_ki)

         ABI_FREE(kg_kf)
       end if
     end do
   end do

 else
   !
   ! More efficienct algorithm based on random access IO:
   !   For each point in the IBZ:
   !     - Read wavefunctions from iwfk
   !     - For each k-point in the star of kpt_ibz:
   !        - Rotate wavefunctions in G-space to get the k-point in the full BZ.
   !        - Write kbz data to file.
   if (out_iomode == IO_MODE_MPI) call wrtout(std_out,"Using MPI-IO to generate full WFK file", do_flush=.True.)
   if (out_iomode == IO_MODE_ETSF) call wrtout(std_out,"Using Netcdf-IO to generate full WFK file", do_flush=.True.)

   ! Construct sorted mapping BZ --> IBZ to speedup qbz search below.
   ABI_MALLOC(iperm, (nkfull))
   ABI_MALLOC(bz2ibz_sort, (nkfull))
   iperm = [(ii, ii=1,nkfull)]
   bz2ibz_sort = bz2ibz(:,1)
   call sort_int(nkfull, bz2ibz_sort, iperm)

   do spin=1,nsppol
     iqst = 0
     do ik_ibz=1,iwfk%nkpt
       nband_k = iwfk%nband(ik_ibz,spin)
       kibz = ebands_ibz%kptns(:,ik_ibz)
       istwf_ki = iwfk%hdr%istwfk(ik_ibz)
       npw_ki = iwfk%hdr%npwarr(ik_ibz)

       call iwfk%read_band_block([1,nband_k],ik_ibz,spin,xmpio_single,&
         kg_k=kg_ki,cg_k=cg_ki,eig_k=eig_ki,occ_k=occ_ki)

       ! Find number of symmetric q-ponts associated to ik_ibz
       nqst = 0
       do ii=iqst+1,nkfull
         if (bz2ibz_sort(ii) /= ik_ibz) exit
         nqst = nqst + 1
       end do
       ABI_CHECK(nqst > 0 .and. bz2ibz_sort(iqst+1) == ik_ibz, "Wrong iqst")

       do jj=1,nqst
         iqst = iqst + 1
         ikf = iperm(iqst)
         ABI_CHECK(ik_ibz == bz2ibz(ikf,1), "ik_ibz !/ ind qq(1)")

         isym = bz2ibz(ikf,2); itimrev = bz2ibz(ikf,6); g0 = bz2ibz(ikf,3:5) ! IS(k_ibz) + g0 = k_bz
         isirred_kf = (isym == 1 .and. itimrev == 0 .and. all(g0 == 0))

         kf = kfull(:,ikf)
         istwf_kf = owfk%hdr%istwfk(ikf)

         ! The test on npwarr is needed because we may change istwfk e.g. gamma.
         if (isirred_kf .and. iwfk%hdr%npwarr(ik_ibz) == owfk%hdr%npwarr(ikf)) then

           call owfk%write_band_block([1,nband_k],ikf,spin,xmpio_single,&
             kg_k=kg_ki,cg_k=cg_ki,eig_k=eig_ki,occ_k=occ_ki)

         else
           ! Compute G-sphere centered on kf
           call get_kg(kf,istwf_kf,ecut_eff,cryst%gmet,npw_kf,kg_kf)
           ABI_CHECK(npw_kf == owfk%hdr%npwarr(ikf), "Wrong npw_kf")

           ! FFT box must enclose the two spheres centered on ki and kf
           gmax_ki = maxval(abs(kg_ki(:,1:npw_ki)), dim=2)
           gmax_kf = maxval(abs(kg_kf), dim=2)
           do ii=1,3
             gmax(ii) = max(gmax_ki(ii), gmax_kf(ii))
           end do
           gmax = 2*gmax + 1
           call ngfft_seq(work_ngfft, gmax)
           ABI_CALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))

           ! Rotate nband_k wavefunctions (output in cg_kf)
           call cgtk_rotate(cryst,kibz,isym,itimrev,g0,nspinor,nband_k,&
             npw_ki,kg_ki,npw_kf,kg_kf,istwf_ki,istwf_kf,cg_ki,cg_kf,work_ngfft,work)

           ABI_FREE(work)

           ! Write data
           call owfk%write_band_block([1,nband_k],ikf,spin,xmpio_single,&
             kg_k=kg_kf,cg_k=cg_kf,eig_k=eig_ki,occ_k=occ_ki)

           ABI_FREE(kg_kf)
         end if
       end do
     end do
   end do

   ABI_FREE(iperm)
   ABI_FREE(bz2ibz_sort)
 end if

 call cwtime_report(" FULL_WFK written to file. ", cpu, wall, gflops)

 ABI_FREE(kg_ki)
 ABI_FREE(cg_ki)
 ABI_FREE(eig_ki)
 ABI_FREE(occ_ki)
 ABI_FREE(bz2ibz)
 ABI_FREE(cg_kf)

 call cryst%free()
 call ebands_free(ebands_ibz)
 call ebands_free(ebands_full)
 call iwfk%close()
 call owfk%close()

end subroutine wfk_tofullbz
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_nc2fort
!! NAME
!!  wfk_nc2fort
!!
!! FUNCTION
!!  Convert a netcdf WFK file (nc_path) to a Fortran WFK file (fort_path).
!!
!! PARENTS
!!      ioprof
!!
!! NOTES
!!  - This routine should be called by a single processor.
!!  - Only GS WFK files are supported (formeig==0)
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_nc2fort(nc_path, fort_path)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: nc_path,fort_path

!Local variables-------------------------------
!scalars
 integer :: ik,spin,mband,mpw,nband_k
 type(wfk_t) :: iwfk,owfk
!arrays
 integer,parameter :: formeig0=0
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)

! *************************************************************************

 call wrtout(std_out, sjoin("Converting:", nc_path, "to", fort_path))

 ! Open input file, extract dimensions and allocate workspace arrays.
 call wfk_open_read(iwfk,nc_path,formeig0,IO_MODE_ETSF,get_unit(),xmpi_comm_self)

 mpw = maxval(iwfk%hdr%npwarr); mband = iwfk%mband
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(cg_k, (2, mpw*iwfk%nspinor*mband))
 ABI_MALLOC(eig_k, ((2*mband)**iwfk%formeig*mband) )
 ABI_MALLOC(occ_k, (mband))

 ! Open output file.
 call owfk%open_write(iwfk%hdr,fort_path,formeig0,IO_MODE_FORTRAN,get_unit(),xmpi_comm_self)

 do spin=1,iwfk%nsppol
   do ik=1,iwfk%nkpt
     nband_k = iwfk%nband(ik,spin)

     call iwfk%read_band_block([1,nband_k],ik,spin,xmpio_single,&
       kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

     call owfk%write_band_block([1,nband_k],ik,spin,xmpio_single,&
       kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
   end do
 end do

 ABI_FREE(kg_k)
 ABI_FREE(cg_k)
 ABI_FREE(eig_k)
 ABI_FREE(occ_k)

 call iwfk%close()
 call owfk%close()

end subroutine wfk_nc2fort
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_prof
!! NAME
!!  wfk_prof
!!
!! FUNCTION
!!  Profiling tool for IO routines
!!
!! INPUTS
!!  wfk_fname=Filename
!!  formeig=0 for GS file, 1 for DFPT file
!!  nband=Number of bands to read.
!!  comm=MPI communicator
!!
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_prof(wfk_fname, formeig, nband, comm)

!Arguments ------------------------------------
 integer,intent(in) :: nband,formeig,comm
 character(len=*),intent(in) :: wfk_fname

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr1=1,master=0,optkg1=1,option1=1,tim_rwwf0=0,icg0=0,headform0=0
 integer :: iomode,wfk_unt,ik_ibz,spin,ierr,ii,option,mband
 integer :: npw_disk,nband_disk,mcg,fform,nband_read,sc_mode,my_rank,nproc
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 type(hdr_type) :: Hdr
 type(wfk_t) :: Wfk
 type(wffile_type) :: wff
 type(MPI_type) :: MPI_enreg_seq
!arrays
 !integer,parameter :: io_modes(2) = (/IO_MODE_FORTRAN, IO_MODE_MPI/)
 integer,parameter :: io_modes(1) = (/IO_MODE_MPI/)
 integer :: ngfft(18)
 logical,allocatable :: my_bmask(:)
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: eig_k(:),cg_k(:,:),occ_k(:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 sc_mode = xmpio_collective

 call hdr_read_from_fname(hdr,wfk_fname,fform,comm)

 ! nband_read is the max number of bands we can read from this file.
 nband_read = nband
 if (nband_read <= 0) then
   nband_read = minval(hdr%nband)
   call wrtout(std_out, sjoin("nband == 0 --> Setting nband_read to:",itoa(nband_read)))
 end if
 if (nband_read > minval(hdr%nband)) then
   nband_read = minval(hdr%nband)
   call wrtout(std_out, sjoin("nband > hdr%nband --> Setting nband_read to:",itoa(nband_read)))
 end if

 wfk_unt = get_unit()

 do ii=1,SIZE(io_modes)
   iomode = io_modes(ii)
   !do option=1,3
   do option=1,3,2
     write(std_out,*)"iomode, option",iomode,option
     call cwtime(cpu,wall,gflops,"start")

     select case (option)
     case (1)
       call wfk_open_read(Wfk,wfk_fname,formeig,iomode,wfk_unt,comm)

       do spin=1,Hdr%nsppol
         do ik_ibz=1,Hdr%nkpt
           npw_disk   = Hdr%npwarr(ik_ibz)
           nband_disk = Wfk%nband(ik_ibz,spin)

           mcg = npw_disk*Hdr%nspinor*nband_read

           ABI_MALLOC(eig_k,((2*Wfk%mband)**formeig*Wfk%mband))
           ABI_MALLOC(occ_k,(Wfk%mband))

           ABI_MALLOC(kg_k,(3,npw_disk))
           ABI_MALLOC_OR_DIE(cg_k,(2,mcg), ierr)

           ! Read the block of bands for this (k,s).
           call wfk%read_band_block([1,nband_read],ik_ibz,spin,xmpio_collective,&
             kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

           ABI_FREE(eig_k)
           ABI_FREE(occ_k)
           ABI_FREE(kg_k)
           ABI_FREE(cg_k)
         end do !ik_ibz
       end do !spin

       call wfk%close()

     case (2)
       call wfk_open_read(Wfk,wfk_fname,formeig,iomode,wfk_unt,comm)

       do spin=1,Hdr%nsppol
         do ik_ibz=1,Hdr%nkpt
           npw_disk   = Hdr%npwarr(ik_ibz)
           nband_disk = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)

           ABI_MALLOC(my_bmask,(MAXVAL(Hdr%nband)))
           my_bmask=.False.; my_bmask(1:nband_read) = .True.

           ABI_MALLOC(eig_k,((2*nband_disk)**formeig*nband_disk))
           ABI_MALLOC(kg_k,(3,npw_disk))
           ABI_MALLOC(occ_k,(nband_disk))

           mcg = npw_disk*Hdr%nspinor*COUNT(my_bmask)
           ABI_MALLOC_OR_DIE(cg_k,(2,mcg), ierr)

           call wfk%read_bmask(my_bmask,ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
           !call wfk%read_band_block((/1,nband_read/),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

           ABI_FREE(my_bmask)
           ABI_FREE(eig_k)
           ABI_FREE(occ_k)
           ABI_FREE(kg_k)
           ABI_FREE(cg_k)
         end do !ik_ibz
       end do !spin

       call wfk%close()

     case (3)
       !Fake MPI_type for the sequential part.
       ngfft(1:6) = (/12,12,12,13,13,13/)
       call initmpi_seq(MPI_enreg_seq)
       call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')
       call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft(2),ngfft(3),'all')

       call WffOpen(iomode,comm,wfk_fname,ierr,wff,master,my_rank,wfk_unt) !,spaceComm_mpiio) ! optional argument
       ABI_CHECK(ierr==0,"ierr!=0")

       call Hdr%free()
       call hdr_io(fform,Hdr,1,wff)
       call WffKg(wff,optkg1)

       do spin=1,Hdr%nsppol
         do ik_ibz=1,Hdr%nkpt

           npw_disk   = Hdr%npwarr(ik_ibz)
           nband_disk = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)

           mband = MAXVAL(Hdr%nband)
           mcg = npw_disk*Hdr%nspinor*nband_read

           ABI_MALLOC(eig_k,((2*mband)**formeig*mband))
           ABI_MALLOC(occ_k,(mband))

           ABI_MALLOC(kg_k,(3,optkg1*npw_disk))
           ABI_MALLOC_OR_DIE(cg_k,(2,mcg), ierr)
           !
           ! Read the block of bands for this (k,s).
           call rwwf(cg_k,eig_k,formeig,headform0,icg0,ik_ibz,spin,kg_k,mband,mcg,MPI_enreg_seq,nband_read,&
             nband_disk,npw_disk,Hdr%nspinor,occ_k,option1,optkg1,tim_rwwf0,Wff)

           ABI_FREE(eig_k)
           ABI_FREE(occ_k)
           ABI_FREE(kg_k)
           ABI_FREE(cg_k)

         end do !ik_ibz
       end do !spin

       call WffClose(wff,ierr)
       call destroy_mpi_enreg(MPI_enreg_seq)

     case default
       MSG_ERROR("Wrong method")
     end select

     call cwtime(cpu,wall,gflops,"stop")
     write(msg,'(3(a,i2),2(a,f8.2))')&
       " iomode: ",iomode,", nproc: ",nproc,", option: ",option,", cpu: ",cpu,", wall:",wall
     call wrtout(std_out,msg,"COLL")
     !call cwtime_report(" FULL_WFK written to file. ", cpu, wall, gflops)
   end do
 end do

 call Hdr%free()

end subroutine wfk_prof
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_create_wfkfile
!! NAME
!!  wfk_create_wfkfile
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_create_wfkfile(wfk_fname,Hdr,iomode,formeig,Kvars,cwtimes,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,formeig,comm
 character(len=*),intent(in) :: wfk_fname
!arrays
 real(dp),intent(out) :: cwtimes(2)
 type(hdr_type),intent(in) :: Hdr
 type(kvars_t),target,intent(out) :: Kvars(Hdr%nkpt)

!Local variables-------------------------------
!scalars
 integer :: nkpt,nsppol,nspinor,ierr,sc_mode
 integer :: ik_ibz,spin,funt,nband_k,npw_k,istwfk_k
 real(dp) :: cpu,wall,gflops,ucvol
 type(wfk_t) :: Wfk
!arrays
 integer :: nband(Hdr%nkpt,Hdr%nsppol)
 integer,pointer :: kg_k(:,:)
 real(dp) :: kpoint(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)

!************************************************************************

 cwtimes = zero

 nband   = RESHAPE(Hdr%nband, (/Hdr%nkpt,Hdr%nsppol/) )
 nkpt    = Hdr%nkpt
 nsppol  = Hdr%nsppol
 nspinor = Hdr%nspinor

 call metric(gmet,gprimd,dev_null,rmet,Hdr%rprimd,ucvol)

 ! Generate the G-vectors from input Hdr%ecut.
 do ik_ibz=1,nkpt
   kpoint   = Hdr%kptns(:,ik_ibz)
   istwfk_k = Hdr%istwfk(ik_ibz)
   call get_kg(kpoint,istwfk_k,Hdr%ecut,gmet,npw_k,Kvars(ik_ibz)%kg_k)
   ABI_CHECK(npw_k == Hdr%npwarr(ik_ibz),"npw_k != Hdr%npwarr(ik)")
 end do
 !
 ! Open the file for writing.
 sc_mode = xmpio_collective

 call cwtime(cpu,wall,gflops,"start")
 funt = get_unit()

 call wfk%open_write(Hdr,wfk_fname,formeig,iomode,funt,comm,write_frm=.TRUE.)

 call cwtime(cpu,wall,gflops,"stop")
 cwtimes = cwtimes + (/cpu,wall/)

 do spin=1,nsppol
   do ik_ibz=1,nkpt

     nband_k = nband(ik_ibz,spin)
     npw_k   = Hdr%npwarr(ik_ibz)
     ABI_MALLOC(cg_k, (2,npw_k*nspinor*nband_k))
     ABI_MALLOC(eig_k, ((2*Wfk%mband)**formeig*Wfk%mband) )
     ABI_MALLOC(occ_k, (Wfk%mband))

     kg_k => Kvars(ik_ibz)%kg_k
     !
     ! Fill cg_k, eig_k, occ_k using a deterministic algorithm so that
     ! we can check the correctness of the reading.
     call fill_or_check("fill",Hdr,Kvars(ik_ibz),ik_ibz,spin,formeig,kg_k,cg_k,eig_k,occ_k,ierr)
     ABI_CHECK(ierr==0,"filling")

     call cwtime(cpu,wall,gflops,"start")

     call wfk%write_band_block((/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

     call cwtime(cpu,wall,gflops,"stop")
     cwtimes = cwtimes + (/cpu,wall/)

     ABI_FREE(cg_k)
     ABI_FREE(eig_k)
     ABI_FREE(occ_k)
   end do
 end do

 ! Close the file
 call wfk%close()

end subroutine wfk_create_wfkfile
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_check_wfkfile
!! NAME
!!  wfk_check_wfkfile
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_check_wfkfile(wfk_fname,Hdr,iomode,method,formeig,Kvars,cwtimes,comm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,formeig,comm,method
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: wfk_fname
!arrays
 real(dp),intent(out) :: cwtimes(2)
 type(hdr_type),intent(in) :: Hdr
 type(kvars_t),intent(in) :: Kvars(Hdr%nkpt)

!Local variables-------------------------------
!scalars
 integer :: nkpt,nsppol,nspinor,ik_ibz,spin,funt,nband_k,npw_k,sc_mode
 integer :: my_ierr,restart,restartpaw,is,ik,ntests,test,mband
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 type(wfk_t) :: Wfk
!arrays
 integer :: nband(Hdr%nkpt,Hdr%nsppol),spins(Hdr%nsppol),kindices(Hdr%nkpt)
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)
 logical,allocatable :: bmask(:)

!************************************************************************

 !write(msg,"(3a,i2)")"Checking file: ",TRIM(wfk_fname),", with iomode = ",iomode
 !call wrtout(std_out,msg,"COLL")

 ierr = 0
 cwtimes = zero

 nband   = RESHAPE(Hdr%nband, (/Hdr%nkpt,Hdr%nsppol/) )
 nkpt    = Hdr%nkpt
 nsppol  = Hdr%nsppol
 nspinor = Hdr%nspinor

 ! Open the file for writing.
 call cwtime(cpu,wall,gflops,"start")
 funt = get_unit()

 call wfk_open_read(Wfk,wfk_fname,formeig,iomode,funt,comm)
 mband = Wfk%mband

 call cwtime(cpu,wall,gflops,"stop")
 cwtimes = cwtimes + (/cpu,wall/)

 ntests = 2

 do test=1,ntests
   spins    = (/(spin, spin=1,Hdr%nsppol)/)
   kindices = (/(ik_ibz, ik_ibz=1,Hdr%nkpt)/)

   if (test==2) then ! Reverse the indices
     spins    = (/(spin, spin=Hdr%nsppol,1,-1)/)
     kindices = (/(ik_ibz, ik_ibz=Hdr%nkpt,1,-1)/)
   end if
   !
   do is=1,SIZE(spins)
     spin = spins(is)
     do ik=1,SIZE(kindices)
       ik_ibz = kindices(ik)

       if (Wfk%debug) then
         call hdr_check(Wfk%fform,Wfk%fform,Hdr,Wfk%Hdr,"COLL",restart,restartpaw)
       end if

       nband_k = nband(ik_ibz,spin)
       npw_k   = Hdr%npwarr(ik_ibz)

       ABI_MALLOC(kg_k, (3,npw_k))
       ABI_MALLOC(cg_k, (2,npw_k*nspinor*nband_k))
       ABI_MALLOC(eig_k, ((2*mband)**Wfk%formeig*mband) )
       ABI_MALLOC(occ_k, (mband))
       !
       !sc_mode = xmpio_collective
       sc_mode = xmpio_single
       ABI_MALLOC(bmask, (mband))
       bmask = .FALSE.
       bmask(1:nband_k) = .TRUE.

       call cwtime(cpu,wall,gflops,"start")

       if (method==0) then
         call wfk%read_band_block((/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
       else if (method==1) then
         call wfk%read_bmask(bmask,ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
       else
         MSG_ERROR("Wrong method")
       end if

       !call wfk%read_eigk(ik_ibz,spin,sc_mode,eig_k)
       !write(std_out,*)"eig_k",eig_k

       call cwtime(cpu,wall,gflops,"stop")
       cwtimes = cwtimes + (/cpu,wall/)

       ! Check the correctness of the reading.
       call fill_or_check("check",Hdr,Kvars(ik_ibz),ik_ibz,spin,formeig,kg_k,cg_k,eig_k,occ_k,my_ierr)

       if (my_ierr/=0) then
         write(msg,"(a,i0)")"fill_or_check returned my_ierr: ",my_ierr
         ierr = my_ierr
         MSG_WARNING(msg)
       end if

       ABI_FREE(kg_k)
       ABI_FREE(cg_k)
       ABI_FREE(eig_k)
       ABI_FREE(occ_k)

       ABI_FREE(bmask)
     end do
   end do
   !
 end do ! test

 ! Close the file
 call wfk%close()

end subroutine wfk_check_wfkfile
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/fill_or_check
!! NAME
!!  fill_or_check
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine fill_or_check(task,Hdr,Kvars,ik_ibz,spin,formeig,kg_k,cg_k,eig_k,occ_k,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,formeig
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: task
!arrays
 integer,intent(in) :: kg_k(:,:)
 real(dp),intent(inout) :: cg_k(:,:),eig_k(:),occ_k(:)
 type(hdr_type),intent(in) :: Hdr
 type(kvars_t),intent(in) :: Kvars

!Local variables-------------------------------
!scalars
 integer :: nkpt,nsppol,nspinor,nband_k,npw_k,band,ipw,kspad,ii,base,idx,mpw,eigsz
 character(len=500) :: msg
!arrays
 integer,allocatable :: ref_kg_k(:,:)
 real(dp),allocatable :: ref_eig_k(:),ref_occ_k(:),ref_cg_k(:,:)
!************************************************************************

 ierr = 0

 nkpt    = Hdr%nkpt
 nsppol  = Hdr%nsppol
 nspinor = Hdr%nspinor

 nband_k = Hdr%nband(ik_ibz + (spin-1)*nkpt)
 npw_k   = Hdr%npwarr(ik_ibz)

 ABI_MALLOC(ref_kg_k,(3,npw_k))
 ABI_MALLOC(ref_eig_k,((2*nband_k)**formeig*nband_k))
 ABI_MALLOC(ref_occ_k,(nband_k))
 ABI_MALLOC(ref_cg_k,(2,npw_k*nspinor*nband_k))

 ref_kg_k = Kvars%kg_k

 ! Pad values according to (k,s).
 kspad = (spin-1)*nkpt + (ik_ibz-1) * nband_k

 if (formeig==0) then
   eigsz = nband_k
   do band=1,nband_k
     ref_eig_k(band) = half * (kspad + band)
     ref_occ_k(band) = two  * (kspad + band)
   end do
 else if (formeig==1) then
   eigsz = 2*nband_k**2
   base=0
   do band=1,nband_k
     do ii=1,2*nband_k
       idx = base + ii
       ref_eig_k(idx) = ii*(kspad + band)
     end do
     base = base + 2*nband_k
   end do
 end if

 mpw = npw_k*nspinor*nband_k
 do ipw=1,mpw
   ref_cg_k(1,ipw) =  ipw + kspad
   ref_cg_k(2,ipw) = -ipw + kspad
 end do

 SELECT CASE (task)
 CASE ("fill")
   cg_k(:,1:mpw) = ref_cg_k(:,1:mpw)
   if (formeig==0) then
     eig_k(1:nband_k) = ref_eig_k
     occ_k(1:nband_k) = ref_occ_k
   else
     eig_k(1:2*nband_k**2) = ref_eig_k
   end if

 CASE ("check")

   if (ANY( ABS(cg_k(:,1:mpw) - ref_cg_k) > zero)) then
     ierr = ierr + 1
     MSG_WARNING("Difference in cg_k")
   end if

   if (ANY( ABS(kg_k - ref_kg_k) > zero)) then
     ierr = ierr + 2
     MSG_WARNING("Difference in kg_k")
     !write(std_out,*)"ref_kg_k",ref_kg_k
     !write(std_out,*)"kg_k",kg_k
   end if

   if (ANY( ABS(eig_k(1:eigsz) - ref_eig_k) > zero)) then
     ierr = ierr + 4
     MSG_WARNING("Difference in eig_k")
     !write(std_out,*)"ref_eig_k",ref_eig_k
     !write(std_out,*)"eig_k",eig_k
   end if

   if (formeig==0) then
     if (ANY( ABS(occ_k(1:nband_k) - ref_occ_k) > zero)) then
       ierr = ierr + 8
       MSG_WARNING("occ_k")
       !write(std_out,*)"ref_occ_k",ref_occ_k
       !write(std_out,*)"occ_k",occ_k
     end if
   end if

   write(msg,"(a,3(i0,2x))")" (ik_ibz, spin, ierr) ",ik_ibz,spin,ierr
   if (ierr/=0) then
     MSG_WARNING(TRIM(msg)//": FAILED")
   else
     call wrtout(std_out,TRIM(msg)//": OK","COLL")
   end if

 CASE DEFAULT
   MSG_ERROR("Wrong task")
 END SELECT

 ABI_FREE(ref_kg_k)
 ABI_FREE(ref_eig_k)
 ABI_FREE(ref_occ_k)
 ABI_FREE(ref_cg_k)

end subroutine fill_or_check
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_diff
!! NAME
!!  wfk_diff
!!
!! FUNCTION
!!  Compare two WFK file for binary equality
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_diff(fname1,fname2,formeig,comm,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: formeig,comm
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: fname1,fname2

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode1,iomode2,ik_ibz,spin,mband,nband_k
 integer :: npw_k,mcg,fform1,fform2,sc_mode,my_rank,nproc
 character(len=500) :: msg
 type(hdr_type) :: Hdr1,Hdr2
 type(wfk_t) :: Wfk1,Wfk2
!arrays
 integer,allocatable :: kg1_k(:,:),kg2_k(:,:)
 real(dp),allocatable :: eig1_k(:),cg1_k(:,:),occ1_k(:)
 real(dp),allocatable :: eig2_k(:),cg2_k(:,:),occ2_k(:)

! *************************************************************************

 call wrtout(std_out, "wfk_diff: comparing "//TRIM(fname1)//" "//TRIM(fname2))

 my_rank = xmpi_comm_rank(comm); nproc   = xmpi_comm_size(comm)
 sc_mode = xmpio_collective

 call hdr_read_from_fname(Hdr1,fname1,fform1,comm)
 call hdr_read_from_fname(Hdr2,fname2,fform2,comm)

 ABI_CHECK(fform1==fform2,"fform1 != fform2")
 ABI_CHECK(Hdr1%nsppol==Hdr2%nsppol,"nsppol1 != nsppol2")
 ABI_CHECK(Hdr1%nspinor==Hdr2%nspinor,"nspinor1 != nspinor2")
 ABI_CHECK(Hdr1%nkpt==Hdr2%nkpt,"nkpt1 != nkpt2")
 !call hdr_check(fform,fform0,hdr1,hdr2,"COLL",restart,restartpaw)

 iomode1 = iomode_from_fname(fname1)
 iomode2 = iomode_from_fname(fname1)
 ABI_CHECK(iomode1==iomode2,"iomode1 != iomode2")
 !iomode1 = IO_MODE_FORTRAN
 !iomode2 = IO_MODE_MPI

 call wfk_open_read(Wfk1,fname1,formeig,iomode1,get_unit(),comm)
 call wfk_open_read(Wfk2,fname2,formeig,iomode2,get_unit(),comm)

 if (wfk1%compare(wfk2) /= 0) then
   MSG_ERROR("WFK files are not consistent. See above messages")
 end if

 mband = Wfk1%mband
 ABI_CHECK(mband==Wfk2%mband,"different mband")
 ABI_CHECK(all(Wfk1%nband==Wfk2%nband),"different nband")
 ABI_CHECK(all(Wfk1%hdr%npwarr==Wfk2%hdr%npwarr),"different npwarr")

 ierr = 0
 do spin=1,Wfk1%nsppol
   do ik_ibz=1,Wfk1%nkpt
     npw_k    = Wfk1%Hdr%npwarr(ik_ibz)
     nband_k  = Wfk1%nband(ik_ibz,spin)
     ABI_CHECK(npw_k  ==Wfk2%Hdr%npwarr(ik_ibz),"different npw_k")
     ABI_CHECK(nband_k==Wfk2%nband(ik_ibz,spin),"different nband_k")

     mcg = npw_k*Hdr1%nspinor*nband_k

     ABI_MALLOC(eig1_k,((2*mband)**formeig*mband))
     ABI_MALLOC(occ1_k,(mband))
     ABI_MALLOC(kg1_k,(3,npw_k))
     ABI_MALLOC_OR_DIE(cg1_k,(2,mcg), ierr)

     ABI_MALLOC(eig2_k,((2*mband)**formeig*mband))
     ABI_MALLOC(occ2_k,(mband))
     ABI_MALLOC(kg2_k,(3,npw_k))
     ABI_MALLOC_OR_DIE(cg2_k,(2,mcg), ierr)

     ! Read the block of bands for this (k,s).
     call wfk1%read_band_block((/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg1_k,eig_k=eig1_k,occ_k=occ1_k) !, cg_k=cg1_k,
     call wfk2%read_band_block((/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg2_k,eig_k=eig2_k,occ_k=occ2_k) !, cg_k=cg2_k,

     if (ANY( ABS(kg1_k - kg2_k) > zero)) then
       ierr = ierr + 2
       MSG_WARNING("Difference in kg_k")
       !write(std_out,*)"kg1_k",kg1_k
       !write(std_out,*)"kg2_k",kg2_k
     end if

     if (ANY( ABS(eig1_k - eig2_k) > zero)) then
       ierr = ierr + 4
       MSG_WARNING("Difference in eig_k")
       !write(std_out,*)"eig1_k",eig1_k
       !write(std_out,*)"eig2_k",eig2_k
     end if

     if (formeig==0) then
       if (ANY( ABS(occ1_k - occ2_k) > zero)) then
         ierr = ierr + 8
         MSG_WARNING("occ_k")
         write(std_out,*)"occ1_k",occ1_k
         write(std_out,*)"occ2_k",occ2_k
       end if
     end if

     if (ANY( ABS(cg1_k - cg2_k) > zero)) then
       ierr = ierr + 1
       MSG_WARNING("Difference in cg_k")
     end if

     write(msg,"(a,3(i0,2x))")" (ik_ibz, spin, ierr) ",ik_ibz,spin,ierr
     if (ierr/=0) then
       MSG_WARNING(TRIM(msg)//": FAILED")
     else
       call wrtout(std_out,TRIM(msg)//": OK","COLL")
     end if

     ABI_FREE(eig1_k)
     ABI_FREE(occ1_k)
     ABI_FREE(kg1_k)
     ABI_FREE(cg1_k)

     ABI_FREE(eig2_k)
     ABI_FREE(occ2_k)
     ABI_FREE(kg2_k)
     ABI_FREE(cg2_k)
   end do !ik_ibz
 end do !spin

 call wfk1%close()
 call wfk2%close()

 call Hdr1%free()
 call Hdr2%free()

end subroutine wfk_diff
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_klist2mesh
!! NAME
!!  wfk_klist2mesh
!!
!! FUNCTION
!! This routine receives a WFK file with wavefunctions given on a subset of k-points belonging to a k-mesh and
!! generates a new WFK file with the complete list of k-points in the IBZ by filling the missing k-points with zeros.
!!
!! This routine is mainly used to prepare the computation of electron mobilities
!! whose convergence with the k-point sampling is notoriously slow.
!! Since only the electron/hole states close to the band edges contribute (say ~0.5 eV),
!! one can reduce significantly the computational cost of the NSCF run by computing
!! a WFK file with kptopt == 0 and the explicit list of k-points located inside the pockets
!! instead of computing all the k-points of the dense IBZ.
!! Unfortunately, the EPH code expects a WFK on a k-mesh so we need to "convert" the initial WFK
!! with the list of k-points to a new WFK file with k-points on the dense kmesh.
!!
!! INPUTS
!!  in_wfkpath = Input WFK file with k-point list.
!!  kerange_path = path to KERANGE.nc file.
!!  dtset <dataset_type>=all input variables for this dataset
!!  out_wfkpath = Output WFK file.
!!  comm = MPI communicator.
!!
!! OUTPUT
!!  Output is written to file out_wfkpath.
!!
!! NOTES
!!  Only GS WFK files are supported (formeig==0)
!!
!! PARENTS
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfk_klist2mesh(in_wfkpath, kerange_path, dtset, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: in_wfkpath, kerange_path
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: comm
!arrays

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0 = 0, master = 0
 integer :: spin, ikf, ikin, nband_k, mpw, mband, nspinor, ierr, fine_mband
 integer :: nsppol, iomode, npw_k, ii, my_rank, ncid, fform, fform_kerange
 real(dp) :: cpu, wall, gflops, mae_meV, merr
 character(len=500) :: msg
 character(len=fnlen) :: my_inpath, out_wfkpath
 type(wfk_t),target :: iwfk
 type(wfk_t) :: owfk
 type(crystal_t) :: cryst
 type(hdr_type) :: fine_hdr
 type(hdr_type),pointer :: ihdr
 type(ebands_t) :: iwfk_ebands, fine_ebands
!arrays
 integer,allocatable :: kf2kin(:), kg_k(:,:) !, kshe_mask(:,:,:)
 real(dp),allocatable :: cg_k(:,:), eig_k(:), occ_k(:), fine_eigen(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ! IO section are executed by master only, all other procs wait for the new WFK before returning.
 my_rank = xmpi_comm_rank(comm); if (my_rank /= master) goto 100

 ! Read interpolated ebands and kshe_mask from KERANGE file, build fine_ebands object.
 ! KERANGE is written by sigtk_kpts_in_erange in m_sigtk module.
#ifdef HAVE_NETCDF
 NCF_CHECK(nctk_open_read(ncid, kerange_path, xmpi_comm_self))
 ! Read header associated to the fine k-mesh
 call hdr_ncread(fine_hdr, ncid, fform)
 fform_kerange = fform_from_ext("KERANGE.nc")
 ABI_CHECK(fform == fform_kerange, sjoin("Wrong fform. Got: ", itoa(fform), ", Expecting: ", itoa(fform_kerange)))
 ! Read eigenvalues and kmask
 fine_mband = maxval(fine_hdr%nband)
 _IBM6("This to prevent xlf from miscompiling this code")
 ABI_MALLOC(fine_eigen, (fine_mband, fine_hdr%nkpt, fine_hdr%nsppol))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "eigenvalues"), fine_eigen))
 !NCF_CHECK(nctk_get_dim(ncid, "nkpt_inerange", nkpt_inerage))
 !ABI_MALLOC(kshe_mask, (fine_hdr%nkpt, fine_hdr%nsppol, 2))
 !NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "kshe_mask"), kshe_mask))
 !ABI_MALLOC(krange2ibz, (nkpt_inerange))
 !NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "krange2ibz"), krange2ibz))
 !ABI_FREE(krange2ibz)
 NCF_CHECK(nf90_close(ncid))
 ! Build fine_ebands
 fine_ebands = ebands_from_hdr(fine_hdr, fine_mband, fine_eigen)
 !call ebands_print(fine_ebands, header="SKW interpolated energies", prtvol=dtset%prtvol)
 ABI_FREE(fine_eigen)
#else
 MSG_ERROR("wfk_klist2mesh requires NETCDF support.")
#endif

 if (my_rank == master) then
   write(std_out, "(2a)")ch10, repeat("=", 92)
   !call wrtout([std_out, ab_out], msg)
   write(std_out, "(a)")" Generating new WKF file with dense k-mesh:"
   write(std_out, "(2a)")" Take wavefunctions with k-point list from WFK file: ", trim(in_wfkpath)
   write(std_out, "(2a)")" Take eigenvalues and k-point tables from KERANGE file: ", trim(kerange_path)
   write(std_out, "(a, 9(i0, 1x))")"   fine_kptrlatt: ", fine_hdr%kptrlatt
   do ii=1,fine_hdr%nshiftk
     write(std_out, "(a, 3(f5.2, 1x))")"   fine_shiftk: ", fine_hdr%shiftk(:, ii)
   end do
   write(std_out, "(2a)")repeat("=", 92), ch10
 end if

 ! Open WFK file with k-point list, extract dimensions and allocate workspace arrays.
 my_inpath = in_wfkpath
 if (nctk_try_fort_or_ncfile(my_inpath, msg) /= 0) then
   MSG_ERROR(msg)
 end if
 iwfk_ebands = wfk_read_ebands(my_inpath, xmpi_comm_self)
 !call ebands_print(iwfk_ebands, header="iwfk_ebands", unit=std_out, prtvol=dtset%prtvol)

 iomode = iomode_from_fname(my_inpath)
 call wfk_open_read(iwfk, my_inpath, formeig0, iomode, get_unit(), xmpi_comm_self)

 if (my_rank == master .and. dtset%prtvol > 0) then
   fform = 0
   call iwfk%hdr%echo(fform, 3, unit=std_out, header="Header of iwfk file")
   call fine_hdr%echo(fform, 3, unit=std_out, header="Header of fine_hdr")
 end if

 ihdr => iwfk%hdr
 mband = iwfk%mband; nsppol = iwfk%nsppol; nspinor = iwfk%nspinor

 cryst = iwfk%hdr%get_crystal(2)

 ! Find correspondence fine kmesh --> input WFK and handle possible mismatch
 !TODO: Write specialized routine wrapping listkk to find mapping without O(N2) scaling.
 ABI_MALLOC(kf2kin, (fine_ebands%nkpt))
 kf2kin = -1
 !call kpts_map(iwfk_ebands%nkpt, iwfk_ebands%kptns, fine_ebands%nkpt, fine_ebands%kptns, kf2kin, xmpi_comm_self)
 do ikf=1,fine_ebands%nkpt
   do ii=1,iwfk_ebands%nkpt
     if (all(abs(fine_ebands%kptns(:, ikf) - iwfk_ebands%kptns(:, ii)) < tol12)) then
       kf2kin(ikf) = ii; exit
     end if
   end do
 end do

 if (count(kf2kin /= -1) /= iwfk_ebands%nkpt) then
   write(msg, "(2a, 2(a,i0))")"Something wrong in the computation of fine_mesh --> input_mesh table.",ch10, &
    "Expecting: ", iwfk_ebands%nkpt, " matches, found: ", count(kf2kin /= -1)
   MSG_ERROR(msg)
 end if

 ! Check weights (the list of k-points should be a subset of the kmesh specified by sigma_ngkpt).
 ierr = 0
 do ikf=1,fine_ebands%nkpt
   ikin = kf2kin(ikf)
   if (ikin == -1) cycle
   if (abs(ihdr%wtk(ikin) - fine_ebands%wtk(ikf)) > tol12) then
     ierr = ierr + 1
     if (ierr <= 10) write(std_out, *) "ihdr%wtk:", ihdr%wtk(ikin), "fine_ebands%wtk", fine_ebands%wtk(ikf)
   end if
 end do
 if (ierr /= 0) then
   write(msg, "(3a)") &
     "Mismatch between input k-weights and weigths associated to the fine mesh. ", ch10, &
     "Possible inconsistency between k-mesh defined by sigma_nshiftk and the list of k-points found in file."
   MSG_ERROR(msg)
 end if

 ! Build new header for output WFK. This is the most delicate part since all the arrays in fine_hdr
 ! that depend on k-points must be consistent with the fine k-mesh.
 mae_meV = zero
 do ikf=1,fine_ebands%nkpt
   ikin = kf2kin(ikf)
   if (ikin == -1) then
     ! Set npwarr to 1 if k-point is not in input set to reduce file size.
     fine_ebands%npwarr(ikf) = 1
     fine_hdr%npwarr(ikf) = 1
   else
     fine_ebands%npwarr(ikf) = iwfk_ebands%npwarr(ikin)
     fine_hdr%npwarr(ikf) = iwfk_ebands%npwarr(ikin)

     ! Insert ab-initio eigenvalues in the SKW-interpolated fine k-mesh.
     do spin=1,nsppol
       nband_k = iwfk_ebands%nband(ikin + (spin - 1) * iwfk_ebands%nkpt)
       merr = Ha_meV * maxval(abs(fine_ebands%eig(1:nband_k, ikf, spin) - iwfk_ebands%eig(1:nband_k, ikin, spin)))
       write(std_out, "(a, es12.4,a)")" MERR: ", merr, " (meV)"
       !if merr >
       !write(std_out, *)fine_ebands%eig(1:nband_k, ikf, spin) * Ha_eV
       !write(std_out, *)iwfk_ebands%eig(1:nband_k, ikin, spin) * Ha_eV
       write(std_out, *)Ha_meV * (fine_ebands%eig(1:nband_k, ikf, spin) - iwfk_ebands%eig(1:nband_k, ikin, spin))
       !end if
       mae_meV = max(mae_meV, merr)
       fine_ebands%eig(1:nband_k, ikf, spin) = iwfk_ebands%eig(1:nband_k, ikin, spin)
       fine_ebands%occ(1:nband_k, ikf, spin) = iwfk_ebands%occ(1:nband_k, ikin, spin)
     end do
   end if
 end do
 write(std_out, "(a, es12.4,a)") &
    " Max error between SKW interpolated energies and ab-initio quantities:", mae_meV, " (meV)"

 !if (mae_meV > ten) then
 !  write(msg,"(2a,2(a,es12.4),a)") &
 !    "Large error in SKW interpolation!",ch10," MARE: ",mare, ", MAE: ", mae_meV, " (meV)"
 !  call wrtout(ab_out, msg)
 !  MSG_WARNING(msg)
 !end if

 call ebands_update_occ(fine_ebands, dtset%spinmagntarget, prtvol=dtset%prtvol)
 !call pack_eneocc(nkpt, nsppol, mband, nband, bantot, array3d, vect)
 !fine_hdr%occ = reshape(fine_ebands%occ, fine_ebands%mband (1:nband_k, ikin, spin)
 call ebands_print(fine_ebands, header="fine_ebands", unit=std_out, prtvol=dtset%prtvol)

 out_wfkpath = strcat(in_wfkpath, ".tmp")
 if (iomode == IO_MODE_ETSF) out_wfkpath = strcat(out_wfkpath, ".nc")
 call owfk%open_write(fine_hdr, out_wfkpath, iwfk%formeig, iomode, get_unit(), xmpi_comm_self)

 if (iomode == IO_MODE_ETSF) then
  ! Add crystal structure and ebands if netcdf output.
#ifdef HAVE_NETCDF
   NCF_CHECK(cryst%ncwrite(owfk%fh))
   NCF_CHECK(ebands_ncwrite(fine_ebands, owfk%fh))
#endif
 end if

 call fine_hdr%free()

 ! Allocate workspace arrays for wavefunction block.
 mpw = maxval(fine_ebands%npwarr)
 _IBM6("This to prevent xlf from miscompiling this code")
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(cg_k, (2, mpw * nspinor * mband))
 ABI_MALLOC(eig_k, ((2*mband)**iwfk%formeig * mband) )
 ABI_MALLOC(occ_k, (mband))

 do spin=1,nsppol
   do ikf=1,fine_ebands%nkpt
     ikin = kf2kin(ikf)
     nband_k = owfk%nband(ikf, spin)
     npw_k = owfk%hdr%npwarr(ikf)

     !cg_k = zero
     if (ikin /= -1) then
       ! Read wavefunctions from input WFK file.
       ABI_CHECK(npw_k == iwfk%hdr%npwarr(ikin), "Mismatch in npw_k")
       ABI_CHECK(nband_k == iwfk%nband(ikin, spin), "Mismatch in nband_k")
       !ABI_CHECK(owfk%hdr%istwfk(ikf) == iwfk%hdristwfk(ikin), "Mismatch in istwfk_k")
       call iwfk%read_band_block([1, nband_k], ikin, spin, xmpio_single, kg_k=kg_k, cg_k=cg_k) !, eig_k=eig_k, occ_k=occ_k)
     else
       ! Fill wavefunctions with fake data (npw_k == 1)
       kg_k = 0
       cg_k = zero
     end if

     ! Write (kpt, spin) block
     eig_k(1:nband_k) = fine_ebands%eig(1:nband_k, ikf, spin)
     occ_k(1:nband_k) = fine_ebands%occ(1:nband_k, ikf, spin)
     call owfk%write_band_block([1, nband_k], ikf, spin, xmpio_single, kg_k=kg_k, cg_k=cg_k, eig_k=eig_k, occ_k=occ_k)
   end do
 end do

 ! Free memory
 ABI_FREE(kg_k)
 ABI_FREE(cg_k)
 ABI_FREE(eig_k)
 ABI_FREE(occ_k)
 ABI_FREE(kf2kin)
 !ABI_FREE(kshe_mask)

 call cryst%free()
 call ebands_free(iwfk_ebands)
 call ebands_free(fine_ebands)
 call iwfk%close()
 call owfk%close()

 ! Rename files, keep backup copy of input WFK file.
 ABI_CHECK(clib_rename(my_inpath, strcat(my_inpath, ".bkp")) == 0, "Failed to rename input WFK file.")
 ABI_CHECK(clib_rename(out_wfkpath, my_inpath) == 0, "Failed to rename output WFK file.")

 call cwtime_report(" WFK with fine k-mesh written to file.", cpu, wall, gflops)

 ! All procs wait here.
100 call xmpi_barrier(comm)

end subroutine wfk_klist2mesh
!!***

!----------------------------------------------------------------------

end module m_wfk
!!***
