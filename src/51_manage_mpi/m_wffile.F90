!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_wffile
!! NAME
!!  m_wffile
!!
!! FUNCTION
!!  This module provides the definition of the wffile_type used to WF file data.
!!  As the type contains MPI-dependent fields, it has to be declared in a MPI-managed directory.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MT,MB,MVer,ZL,MD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!! wffile_type : a handler for dealing with the IO of a wavefunction file
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_wffile

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_xmpi
#ifdef HAVE_MPI2
 use mpi
#endif
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,   only : mvrecord, open_file
 use m_fstrings,   only : toupper, endswith, sjoin

 implicit none

 private

#ifdef HAVE_MPI1
include 'mpif.h'
#endif

#define DEV_DEBUG_THIS 0

!public procedures.
 public :: WffOpen
 public :: wffclose
 public :: wffdelete
 public :: wffkg
 public :: wffoffset
 public :: wffreaddatarec  !  Generic subroutines to read data in one record of a wavefunction file
 public :: wffreadnpwrec
 public :: wffreadskiprec
 public :: wffreadwrite_mpio
 public :: wffwritedatarec !  Generic subroutines to write data in one record of a wavefunction file
 public :: wffwritenpwrec
 public :: xderiveread     !  Generic subroutines to read wf files.
 public :: xmpi_read_int2d
 public :: xderivewrite

 public :: getRecordMarkerLength_wffile
 public :: xnullifyOff
 public :: xmoveOff
 public :: xderiveWRecEnd
 public :: xderiveWRecInit
 public :: xderiveRRecEnd
 public :: xderiveRRecInit
#if defined HAVE_MPI_IO
 public :: rwRecordMarker
#endif
 public :: clsopn
 public :: wff_usef90
 public :: xdefineOff
!!***

!Generic interface of the routines wffreaddatarec
 interface wffreaddatarec
   module procedure WffReadDataRec_dp1d
   module procedure WffReadDataRec_dp2d
 end interface wffreaddatarec

 !Generic interface of the routines wffwritedatarec
 interface wffwritedatarec
   module procedure WffWriteDataRec_int2d
   module procedure WffWriteDataRec_dp1d
   module procedure WffWriteDataRec_dp2d
 end interface

 !Generic interface of the routines xderiveread
 interface xderiveread
   module procedure xderiveRead_int             !  read integer value
   module procedure xderiveRead_int1d           !  read integer array 1d
   module procedure xderiveRead_int2d           !  read integer array 2d
   module procedure xderiveRead_dp              !  read double precision value
   module procedure xderiveRead_dp1d            !  read double precision array 1d
   module procedure xderiveRead_dp2d            !  read double precision array 2d
   module procedure xderiveRead_int2d_displ     !  read integer array 2d non-contiguous
   module procedure xderiveRead_dp2d_displ      !  read double precision array 2d non-contiguous
   module procedure xderiveReadVal_char         !  read character string
   module procedure xmpi_read_int2d
 end interface xderiveread

!Generic interface of the routines xderivewrite
 interface xderivewrite
  module procedure xderiveWrite_int             ! write integer value
  module procedure xderiveWrite_int1d           ! write integer array 1d
  module procedure xderiveWrite_int2d           ! write integer array 2d
  module procedure xderiveWrite_dp              ! write double precision value
  module procedure xderiveWrite_dp1d            ! write double precision array 1d
  module procedure xderiveWrite_dp2d            ! write double precision array 2d
  module procedure xderiveWrite_dp2d_seq        ! write double precision array 2d in sequential
  module procedure xderiveWrite_int2d_displ     ! write integer array 2d non contiguous
  module procedure xderiveWrite_dp2d_displ      ! write double precision array 2d non contiguous
  module procedure xderiveWrite_char            ! write character string
 end interface xderivewrite

!!****t* m_wffile/wffile_type
!! NAME
!! wffile_type
!!
!! FUNCTION
!! This structure datatype is a handler for dealing with the IO of a
!! wavefunction file.
!! It contains, among other things, the method of access to the file
!! (standard F90 read/write, or NetCDF call, or MPI IO), the unit number
!! if applicable, the filename, the information on the
!! parallelism, etc ...
!!
!! SOURCE

 type, public :: wffile_type

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar
  integer :: unwff
   ! unwff  unit number of unformatted wavefunction disk file

  integer :: iomode
   ! Method to access the wavefunction file
   !   IO_MODE_FORTRAN for usual Fortran IO routines
   !   IO_MODE_FORTRAN_MASTER if usual Fortran IO routines, but only the master node in the parallel case
   !   IO_MODE_MPI if MPI/IO routines (this access method is only available in parallel)
   !   IO_MODE_NETCDF if NetCDF routines (obsolete, do not use)
   !   IO_MODE_ETSF, NetCDF format read via etsf-io.

  integer :: formwff
   ! formwff=format of the eigenvalues
   !   -1 => not used
   !    0 => vector of eigenvalues
   !    1 => hermitian matrix of eigenvalues

  integer :: headform
   ! headform=format of the header

  integer ::  kgwff
   ! kgwff  if 1 , read or write kg_k ; if 0, do not care about kg_k

! Character
  character(len=fnlen) :: fname
   ! filename (if available)

! In case of MPI parallel use
  integer :: master
   ! index of the processor master of the IO procedure when the WffOpen call is issued

  integer :: me
   ! index of my processor in the spaceComm communicator

  integer :: me_mpiio
   ! index of my processor in the spaceComm_mpiio communicator

  integer :: nproc
   ! number of processors that will have access to the file

  integer :: spaceComm
   ! space communicator for the standard FORTRAN access to the file

  integer :: spaceComm_mpiio
   ! space communicator for the MPI/IO access to the file

! In case of MPI/IO : additional information
  integer :: fhwff
   ! file handle used to access the file with MPI/IO.

  integer(kind=XMPI_OFFSET_KIND) :: nbOct_int,nbOct_dp,nbOct_ch
   ! nbOct_int byte number of int value
   ! nbOct_dp  byte number of dp value
   ! nbOct_ch  byte number of character value

  integer(kind=XMPI_OFFSET_KIND) :: nbOct_recMarker
   ! byte number of Fortran file record markers

  integer(kind=XMPI_OFFSET_KIND) :: lght_recs
   ! length of record

  integer :: marker_mpi_type
   ! MPI Datatype for Fortran record markers

  integer(kind=XMPI_OFFSET_KIND)  :: offwff,off_recs
   ! offwff   offset position of unformatted wavefunction disk file
   ! off_recs offset position of start record
   ! (used in parallel MPI-IO)

  integer :: offset_mpi_type
   ! MPI Datatype for INTEGER(kind=MPI_OFFSET_KIND)

 end type wffile_type


CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/getRecordMarkerLength_wffile
!! NAME
!!  getRecordMarkerLength_wffile
!!
!! FUNCTION
!!  Get the record marker length of the FORTRAN header of a file to access it in MPI/IO.
!!  This routine assumes that the header has been written (and flushed) in the file.
!!
!! SIDE EFFECTS
!!  wff=<type(wffile_type)>=structured info for reading/writing the wavefunctions
!!      only%nbOct_recMarker is changed
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine getRecordMarkerLength_wffile(wff)

!Arguments ------------------------------------
!scalars
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: headform,ierr,ii,iimax
 integer(kind=MPI_OFFSET_KIND)  :: posit,rml
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)
#endif

!************************************************************************

#ifdef DEV_DEBUG_THIS
 return
 ! Already done in WffOpen
#endif

#if defined HAVE_MPI_IO

 if (wff%nbOct_recMarker>0) return

!wff%nbOct_recMarker=4;return
!call flush(wff%unwff)
!call MPI_FILE_SYNC(wff%fhwff,ierr)

!Only master do that
 ierr=0
 if (wff%master==wff%me) then

! Define number of INTEGER types to be tested
#if defined HAVE_FC_INT_QUAD
   iimax=4
#else
   iimax=3
#endif

! Try to read headform
   rml=-1;ii=0
   do while (wff%nbOct_recMarker<=0.and.ii<iimax)
     ii=ii+1
     if (ii==1) rml=4
     if (ii==2) rml=8
     if (ii==3) rml=2
     if (ii==4) rml=16
     posit=rml+6*wff%nbOct_ch
     call MPI_FILE_READ_AT(wff%fhwff,posit,headform,1,MPI_INTEGER,statux,ierr)
     if (ierr==MPI_SUCCESS) then
       if (headform==wff%headform) wff%nbOct_recMarker=rml
     end if
    end do

    if (ierr/=MPI_SUCCESS) then
     MSG_BUG("Header problem")
    end if

   if (ii==iimax.and.wff%nbOct_recMarker<=0) then
!      if (iimax>=4) then
!        write(msg,'(3a)') &
! &        ' Your architecture is not able to handle 16, 8, 4 or 2-bytes FORTRAN file record markers !',ch10,&
! &        ' You cannot use ABINIT and MPI/IO.'
!      else
!        write(msg,'(3a)') &
! &        '  Your architecture is not able to handle 8, 4 or 2-bytes FORTRAN file record markers !',ch10,&
! &        '  You cannot use ABINIT and MPI/IO.'
!      end if
     write(msg,'(13a)') &
&      ' Error during FORTRAN file record marker detection:',ch10,&
&      ' It was not possible to read/write a small file!',ch10,&
&      ' ACTION: check your access permissions to the file system.',ch10,&
&      ' Common sources of this problem:',ch10,&
&      '  - Quota limit exceeded,',ch10,&
&      '  - R/W incorrect permissions,',ch10,&
&      '  - WFK file requested as input (irdwfk=1/getwfk=1) but not existing ...'
     MSG_ERROR(msg)
   else
     write(msg,'(a,i0)') &
&     '  MPI/IO accessing FORTRAN file header: detected record mark length=',wff%nbOct_recMarker
     MSG_COMMENT(msg)
   end if

 end if  ! me=master

!Broadcast record marker length
 if (wff%spaceComm/=MPI_COMM_SELF) then
   call MPI_BCAST(wff%nbOct_recMarker,1,wff%offset_mpi_type,wff%master,wff%spaceComm,ierr)
 end if

!Select MPI datatype for markers
 if (wff%nbOct_recMarker==4) then
   wff%marker_mpi_type=MPI_INTEGER4
 else if (wff%nbOct_recMarker==8) then
   wff%marker_mpi_type=MPI_INTEGER8
#if defined HAVE_FC_INT_QUAD && defined HAVE_MPI_INTEGER16
 else if (wff%nbOct_recMarker==16) then
   wff%marker_mpi_type=MPI_INTEGER16
#endif
 else if (wff%nbOct_recMarker==2) then
   wff%marker_mpi_type=MPI_INTEGER2
 end if

#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine getRecordMarkerLength_wffile
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/rwRecordMarker
!! NAME
!!  rwRecordMarker
!!
!! FUNCTION
!!  Read/Write a record marker in a FORTRAN file at a given file pointer position.
!!  This is needed to access data in a FORTRAN file with MPI/IO.
!!
!! INPUTS
!!  option=1 for reading by current proc
!!         2 for writing by current proc
!!         3 for reading by all procs
!!         4 for writing by all procs
!!  posit= position of the MPI/IO file pointer
!!  wff=<type(wffile_type)>=structured info for reading/writing
!!     Use here only:
!!       wff%fhwff= handle of the MPI/IO file
!!       wff%nbOct_recMarker= length of Fortran record markers
!!
!! OUTPUT
!!  ierr= error code
!!
!! SIDE EFFECTS
!!  posit= position of the MPI/IO file pointer
!!         updated after the reading (with the length of the record)
!!  recordmarker= content of the record marker
!!
!! PARENTS
!!      m_hdr,m_wffile
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine rwRecordMarker(option,posit,recordmarker,wff,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
 integer(kind=MPI_OFFSET_KIND),intent(inout) :: posit,recordmarker
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer(kind=2)  :: delim_record2(1)
 integer(kind=4)  :: delim_record4(1)
 integer(kind=8)  :: delim_record8(1)
#if defined HAVE_FC_INT_QUAD
 integer(kind=16) :: delim_record16(1)
#endif
!character(len=500) :: msg
!arrays
 integer  :: statux(MPI_STATUS_SIZE)

!************************************************************************

 ierr=0

 if (option==1) then
   if (wff%nbOct_recMarker==4) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record4(1)
   else if (wff%nbOct_recMarker==8) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record8(1)
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record16(1)
#endif
   else if (wff%nbOct_recMarker==2) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record2(1)
   else
     MSG_BUG('Wrong record marker length!')
   end if

 else if (option==2) then
   if (wff%nbOct_recMarker==4) then
     delim_record4 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
   else if (wff%nbOct_recMarker==8) then
     delim_record8 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     delim_record16 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
#endif
   else if (wff%nbOct_recMarker==2) then
     delim_record2 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
   else
     MSG_BUG('Wrong record marker length!')
   end if

 else if (option==3) then
   if (wff%nbOct_recMarker==4) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record4(1)
   else if (wff%nbOct_recMarker==8) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record8(1)
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record16(1)
#endif
   else if (wff%nbOct_recMarker==2) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record2(1)
   else
     MSG_BUG('Wrong record marker length !')
   end if

 else if (option==4) then
   if (wff%nbOct_recMarker==4) then
     delim_record4 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
   else if (wff%nbOct_recMarker==8) then
     delim_record8 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     delim_record16 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
#endif
   else if (wff%nbOct_recMarker==2) then
     delim_record2 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
   else
     MSG_BUG('Wrong record marker length!')
   end if

 else
   MSG_BUG('Wrong value for option!')
 end if

 posit = posit + recordmarker + 2*wff%nbOct_recMarker

end subroutine rwRecordMarker
#endif
!!***

!------------------------------------------------------------------------------------

!!****f* m_wffile/xnullifyOff
!! NAME
!!  xnullifyOff
!!
!! FUNCTION
!!  In case of MPI I/O, nullify the offset of a WF file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wff=<type(wffile_type)>=structured info for reading/writing
!!
!! PARENTS
!!      m_wffile
!!
!! CHILDREN
!!
!! SOURCE

subroutine xnullifyOff(wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff

! *************************************************************************

#if defined HAVE_MPI_IO
 wff%offwff    = 0
 wff%off_recs  = 0
 wff%lght_recs = 0
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xnullifyOff
!!***

!------------------------------------------------------------------------------------

!!****f* m_wffile/xmoveOff
!! NAME
!!  xmoveOff
!!
!! FUNCTION
!!  In case of MPI I/O, move the offset of a WF file
!!
!! INPUTS
!!  [n_int] = number if integers to skip
!!  [n_dp]  = number if double precision reals to skip
!!  [n_ch]  = number if characters to skip
!!  [n_mark]= number if record markers to skip
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wff=<type(wffile_type)>=structured info for reading/writing
!!
!! PARENTS
!!      posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmoveOff(wff,n_int,n_dp,n_ch,n_mark)

!Arguments ------------------------------------
 integer,intent(in),optional :: n_int,n_dp,n_ch,n_mark
 type(wffile_type),intent(inout) :: wff

! *************************************************************************

#if defined HAVE_MPI_IO
 if (present(n_int) ) wff%offwff=wff%offwff+n_int *wff%nbOct_int
 if (present(n_dp)  ) wff%offwff=wff%offwff+n_dp  *wff%nbOct_dp
 if (present(n_ch)  ) wff%offwff=wff%offwff+n_ch  *wff%nbOct_ch
 if (present(n_mark)) wff%offwff=wff%offwff+n_mark*wff%nbOct_recMarker
#else
!This section should not be used...
 if (present(n_int) .and.(.false.)) write(std_out,*) n_int
 if (present(n_dp)  .and.(.false.)) write(std_out,*) n_dp
 if (present(n_ch)  .and.(.false.)) write(std_out,*) n_ch
 if (present(n_mark).and.(.false.)) write(std_out,*) n_mark
#endif

end subroutine xmoveOff
!!***

!------------------------------------------------------------------------------------

!!****f* m_wffile/xderiveWRecEnd
!! NAME
!!  xderiveWRecEnd
!!
!! FUNCTION
!!  Writes the first and last wavefunction block marker using MPI/IO
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      m_ioarr,m_wffile,outxfhist,posdoppler,rwwf
!!
!! CHILDREN
!!
!! NOTES
!!  We assume that:
!!    wff%offwff contains the position of the end of the record
!!    wff%off_recs contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveWRecEnd(wff,ierr,me_proc)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in),optional :: me_proc
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: me
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit
!arrays
#endif

! *************************************************************************

 ierr=0

#if defined HAVE_MPI_IO
 me=-1;if (present(me_proc)) me=me_proc
 if (me==-1.or.me==0) then

   delim_record=wff%offwff-wff%off_recs-wff%nbOct_recMarker

!  Write the first word of the record
   posit=wff%off_recs
   call rwRecordMarker(2,posit,delim_record,wff,ierr)

!  Write the last word of the record
   posit=wff%offwff
   call rwRecordMarker(2,posit,delim_record,wff,ierr)

 end if

 wff%offwff = wff%offwff + wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED((/wff%me,me_proc/))

end subroutine xderiveWRecEnd
!!***

!------------------------------------------------------------------------------

!!****f* m_wffile/xderiveWRecInit
!! NAME
!!  xderiveWRecInit
!!
!! FUNCTION
!!  Writes the first wavefunction block marker using MPI/IO.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      m_ioarr,m_wffile,outxfhist,posdoppler,rwwf
!!
!! CHILDREN
!!
!! NOTES
!!  We assume that:
!!    wff%offwff contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveWRecInit(wff,ierr,me_proc)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in),optional :: me_proc
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: me
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit
!arrays
#endif

! *************************************************************************

 ierr=0

#if defined HAVE_MPI_IO
 me=-1;if (present(me_proc)) me=me_proc
 if (me==-1.or.me==0) then

!  Write the first word of the record
   posit=wff%offwff;delim_record=0
   call rwRecordMarker(2,posit,delim_record,wff,ierr)

 end if

 wff%off_recs = wff%offwff
 wff%offwff = wff%offwff + wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED((/wff%me,me_proc/))

end subroutine xderiveWRecInit
!!***

!---------------------------------------------------------------------------------

!!****f* m_wffile/xderiveRRecEnd
!! NAME
!!  xderiveRRecEnd
!!
!! FUNCTION
!!  Initializes the end-of-record offset for MPI/IO.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      m_ioarr,m_wffile,outxfhist,rwwf
!!
!! CHILDREN
!!
!! NOTES
!!  We assume that:
!!    wff%off_recs contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveRRecEnd(wff,ierr)

!Arguments ------------------------------------
 integer,intent(out) ::  ierr
 type(wffile_type),intent(inout) :: wff

! *************************************************************************

 ierr=0
#if defined HAVE_MPI_IO
!Define offset end of record
 wff%offwff = wff%off_recs + wff%lght_recs + 2*wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xderiveRRecEnd
!!***

!-------------------------------------------------------------------------------

!!****f* m_wffile/xderiveRRecInit
!! NAME
!!  xderiveRRecInit
!!
!! FUNCTION
!!  Initializes the record length for MPI/IO.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      m_ioarr,m_wffile,outxfhist,rwwf
!!
!! CHILDREN
!!
!! NOTES
!!  We assume that:
!!    wff%offwff contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveRRecInit(wff,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit
#endif

! *************************************************************************

 ierr=0

#if defined HAVE_MPI_IO
 wff%off_recs = wff%offwff

!Read the length of the record
 posit=wff%off_recs
 call rwRecordMarker(1,posit,delim_record,wff,ierr)

 wff%lght_recs = delim_record
 wff%offwff =  wff%offwff + wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xderiveRRecInit
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/clsopn
!! NAME
!! clsopn
!!
!! FUNCTION
!! Close wavefunction file (provided its access is standard F90 IO), then reopen the same.
!! Uses fortran inquire statement to reopen with same characteristics.
!!
!! INPUTS
!!  wff=number of unit to which on which file is already opened.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine clsopn(wff)

!Arguments ------------------------------------
!scalars
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer :: ios,unit
 logical :: nmd,od
 character(len=11) :: fm
 character(len=500) :: message
 character(len=fnlen) :: filnam

! *************************************************************************

 if ( ANY(wff%iomode==(/IO_MODE_FORTRAN_MASTER,IO_MODE_FORTRAN/) ))then

   unit=wff%unwff
   inquire (unit=unit,iostat=ios,opened=od,name=filnam,form=fm,named=nmd)

!  ios is a status specifier.  If an error condition exists,
!  ios is assigned a processor-dependent value > 0.
   if (ios/=0) then
     write(message, '(/,a,/,a,i8,a,i8,/,a,/,a,/,a)' ) &
&     ' clsopn : ERROR -',&
&     '  Attempt to inquire about unit=',unit,&
&     '  indicates error condition iostat=',ios,&
&     '  May be due to temporary problem with file, disks or network.',&
&     '  Action: check whether there might be some external problem,',&
&     '  then resubmit.'
     MSG_ERROR(message)

!    od is a logical variable which is set to true if the specified
!    unit is connected to a file; otherwise it is set to false.
#if !defined FC_HITACHI
   else if (.not.od) then
     write(message, '(/,a,/,a,i8,/,a,/,a,/,a,/,a)' ) &
&     ' clsopn : ERROR -',&
&     '  Tried to inquire about unit',unit,&
&     '  and found it not connected to a file.',&
&     '  May be due to temporary problem with file, disks or network.',&
&     '  Action: check whether there might be some external problem,',&
&     '  then resubmit.'
     MSG_ERROR(message)
#endif

!    nmd is a logical variable assigned the value true if the file
!    has a name; otherwise false.  A scratch file is not named.
   else if (.not.nmd) then

!    No action for the time being. Possibility to debug.

   else

!    May now close the file and then reopen it
!    (file is already opened according to above checks)

#if defined FC_HITACHI
     if (.not.od) then
       write(message, '(a,i0,/,a,/,a,/,a)' ) &
&       '  Tried to inquire about unit',unit,&
&       '  and found it not connected to a file.',&
&       '  May be due to temporary problem with file, disks or network.',&
&       '  Action: disregard this error and continue the process anyway.'
       MSG_WARNING(message)
     end if
#endif
     close (unit=unit)
     open (unit=unit,file=filnam,form=fm,status='old') !VALGRIND complains filnam is just a few thousand bytes inside a block of 8300

   end if

 else if (wff%iomode == IO_MODE_MPI) then
   call xnullifyOff(wff)
 else if (wff%iomode == IO_MODE_ETSF) then
!  We do nothing, ETSF access already not being sequential.
 end if

end subroutine clsopn
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/wff_usef90
!! NAME
!! wff_usef90
!!
!! FUNCTION
!!  1 if a Fortran file is going to be read by this node, 0 otherwise.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wff_usef90(wff)

!Arguments ------------------------------------
!scalars
 integer :: wff_usef90
 type(wffile_type),intent(in) :: wff

! *************************************************************************

 wff_usef90=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode ==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) wff_usef90=1

end function wff_usef90
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/wff_ireadf90
!! NAME
!! wff_ireadf90
!!
!! FUNCTION
!!  1 if a Fortran file is going to be read by this node, 0 otherwise.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wff_ireadf90(wff)

!Arguments ------------------------------------
!scalars
 integer :: wff_ireadf90
 type(wffile_type),intent(in) :: wff

! *************************************************************************

 wff_ireadf90=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) wff_ireadf90=1

end function wff_ireadf90
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffOpen
!! NAME
!! WffOpen
!!
!! FUNCTION
!! This subroutine opens a Wf file. It might be accessed
!! by different mechanisms (usual F90 IO routines,
!!  MPI I/O, or, in the future, NetCDF). The routine
!! provides a file handler, wff (a data structure containing
!! all needed information).
!!
!! INPUTS
!! iomode=access mode (0 means all procs access using usual F90
!!  routines ; -1 means only the master proc access, using usual
!!  F90 routines ; 1 means MPI I/O; 2 means netcdf I/O)
!! filename=name of the file
!! master=the number of the master proc (only needed in parallel)
!! me=my number (only needed in parallel)
!! spaceComm= the space communicator handler (only needed in MPI parallel I/O)
!! spaceWorld= the space communicator for the whole set of procs
!! unwff=the file unit number
!!
!! OUTPUT
!! ier=error code
!! wff= structured info about the wavefunction file
!!
!! PARENTS
!!      conducti_paw,conducti_paw_core,emispec_paw,inwffil,linear_optics_paw
!!      m_ioarr,m_iowf,m_wfk,optics_paw,optics_paw_core,optics_vloc,posdoppler
!!      uderiv
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffOpen(iomode,spaceComm,filename,ier,wff,master,me,unwff,&
&                  spaceComm_mpiio) ! optional argument

!Arguments ------------------------------------
 integer, intent(in)  :: iomode,spaceComm,master,me,unwff
 integer, intent(in),optional  :: spaceComm_mpiio
 integer, intent(out) :: ier
 character(len=fnlen), intent(in) :: filename
 type(wffile_type), intent(inout) :: wff !vz_i

!Local variables-------------------------------
 character(len=500) :: message
 character(len=fnlen) :: fildata
#ifdef HAVE_MPI_IO
 integer :: isize
#endif

! *************************************************************************

!Initialize the mandatory data of the wff datastructure
 wff%unwff  =unwff
 wff%iomode =iomode; if (endswith(filename, ".nc")) wff%iomode = IO_MODE_ETSF
 if (filename/=wff%fname) wff%fname=filename

!Initialize info useful for parallel use
 wff%nproc    =1
 wff%master   =master
 wff%me       =me
 wff%me_mpiio =0
 wff%spaceComm=spaceComm
 wff%spaceComm_mpiio=xmpi_comm_self

#if defined HAVE_MPI
! This case occurs when wff is connected to a DENSITY file
! abinit_comm_output is generally equal to MPI_COMM_WORLD (except if paral. over images)
  if (spaceComm==MPI_COMM_SELF) wff%spaceComm=abinit_comm_output
! if (spaceComm==MPI_COMM_SELF) wff%spaceComm=MPI_COMM_WORLD
  call MPI_COMM_SIZE(wff%spaceComm,wff%nproc,ier)
! Redefine the default MPIIO communicator if MPI, although MPIIO features should not be used unless
! present(spaceComm_mpiio).and.wff%iomode==1
  wff%spaceComm_mpiio=wff%spaceComm
  wff%me_mpiio=wff%me
#endif

  if (present(spaceComm_mpiio).and.any(wff%iomode==[IO_MODE_MPI, IO_MODE_ETSF])) wff%spaceComm_mpiio=spaceComm_mpiio
#if defined HAVE_MPI
  call MPI_COMM_RANK(wff%spaceComm_mpiio,wff%me_mpiio,ier)
#endif

 ier=0
 if (wff%iomode==IO_MODE_FORTRAN) then !  All processors see a local file
   if (open_file(filename, message, unit=unwff, form="unformatted") /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unwff)

 else if (wff%iomode==IO_MODE_FORTRAN_MASTER) then !  Only the master processor see a local file
   if(master==me)then
     if (open_file(filename, message, unit=unwff, form="unformatted") /= 0) then
       MSG_ERROR(message)
     end if
     rewind(unwff)
   end if

#if defined HAVE_MPI_IO
 else if (wff%iomode==IO_MODE_MPI)then ! In the parallel case, only the master open filename file
   if(master==me)then
     if (open_file(filename, message, unit=unwff, form="unformatted") /= 0) then
       MSG_ERROR(message)
     end if
     rewind(unwff)
   end if
   ! MG: Great! These barriers lead to a deadlock if prtded hence MPI_FILE_OPEN is not called by all the processors!
   !call xmpi_barrier(wff%spaceComm)
   !call xmpi_barrier(wff%spaceComm_mpiio)

   call MPI_FILE_OPEN(wff%spaceComm,filename,MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,wff%fhwff,ier)
   ABI_CHECK_MPI(ier,sjoin("WffOpen:", filename))

!  Define all type values
   call MPI_Type_size(MPI_INTEGER,isize,ier)
   wff%nbOct_int=isize
   call MPI_Type_size(MPI_DOUBLE_PRECISION,isize,ier)
   wff%nbOct_dp=isize
   call MPI_Type_size(MPI_CHARACTER,isize,ier)
   wff%nbOct_ch=isize
   wff%nbOct_recMarker=-1;wff%kgwff=-1;wff%formwff=-1
   wff%offwff=0;wff%off_recs=0;wff%lght_recs=0
   wff%marker_mpi_type=MPI_INTEGER ! Default value

#ifdef DEV_DEBUG_THIS
   wff%nbOct_recMarker=xmpio_bsize_frm
   wff%marker_mpi_type=xmpio_mpi_type_frm
#endif

   if (MPI_OFFSET_KIND==4) then
     wff%offset_mpi_type=MPI_INTEGER4
   else  if (MPI_OFFSET_KIND==8) then
     wff%offset_mpi_type=MPI_INTEGER8
#if defined HAVE_FC_INT_QUAD && defined HAVE_MPI_INTEGER16
   else  if (MPI_OFFSET_KIND==16) then
     wff%offset_mpi_type=MPI_INTEGER16
#endif
   else  if (MPI_OFFSET_KIND==2) then
     wff%offset_mpi_type=MPI_INTEGER2
   end if
#endif

#ifdef HAVE_NETCDF
 else if (wff%iomode==IO_MODE_ETSF)then
   fildata = nctk_ncify(filename)
   NCF_CHECK(nctk_open_modify(wff%unwff, fildata, xmpi_comm_self))
   wff%fname = fildata
   !write(message,'(3A,I0)')'WffOpen: opening ', trim(wff%fname)," on unit ", wff%unwff
   !call wrtout(std_out, message, 'COLL')
#endif
 else
   write(message, '(7a,i0,3a)' )&
&   'For the time being the input variable iomode is restricted ',ch10,&
&   'to 0 (all cases), 1 (in case MPI is enabled),',ch10,&
&   'or 3 (only sequential, and if the NetCDF and ETSF_IO libraries have been enabled).',ch10,&
&   'Its value is iomode= ',wff%iomode,'.',ch10,&
&   'Action: change iomode or use ABINIT in parallel or enable NetCDF and/or ETSF_IO.'
   MSG_ERROR(message)
 end if

end subroutine WffOpen
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffClose
!! NAME
!! WffClose
!!
!! FUNCTION
!! This subroutine closes a Wf file.
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! PARENTS
!!      conducti_paw,conducti_paw_core,dfpt_looppert,dfptnl_loop,emispec_paw
!!      gstate,m_ioarr,m_iowf,m_wfk,nonlinear,optics_paw,optics_paw_core
!!      optics_vloc,posdoppler,respfn,uderiv
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffClose(wff,ier)

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(out) :: ier

! *************************************************************************

 ier=0
 if(wff%iomode==IO_MODE_FORTRAN) then ! All processors see a local file
   close(unit=wff%unwff)

#ifdef HAVE_NETCDF
 else if(wff%iomode == IO_MODE_ETSF)then
   NCF_CHECK(nf90_close(wff%unwff))
#endif

 else if(wff%iomode==IO_MODE_FORTRAN_MASTER)then !  Only the master processor see a local file
   if(wff%master==wff%me) close (unit=wff%unwff)    ! VALGRIND complains buf points to uninitialized bytes

#if defined HAVE_MPI_IO
 else if(wff%iomode==IO_MODE_MPI)then
   call MPI_FILE_CLOSE(wff%fhwff,ier)
   if (wff%master==wff%me ) close(unit=wff%unwff)
   wff%offwff=0;wff%off_recs=0;wff%lght_recs=0
   wff%nbOct_recMarker=-1
   wff%kgwff=-1
#endif

 end if

end subroutine WffClose
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffDelete
!! NAME
!! WffDelete
!!
!! FUNCTION
!! This subroutine closes a Wf file, and delete it.
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffDelete(wff,ier)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ier

! *************************************************************************

 ier=0
 if (wff%iomode==IO_MODE_FORTRAN) then !  All processors see a local file
   close(unit=wff%unwff,status='delete')

 else if (wff%iomode==IO_MODE_FORTRAN_MASTER)then !  Only the master processor see a local file
   if (wff%master==wff%me) close (unit=wff%unwff,status='delete')


 else if (wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   if ( wff%fhwff /= -1 )then
     call MPI_FILE_CLOSE(wff%fhwff,ier)
   end if
   if (wff%master==wff%me ) then
     close(unit=wff%unwff,status='delete')
     wff%fhwff = -1
   end if
   wff%offwff=0;wff%off_recs=0;wff%lght_recs=0
   wff%nbOct_recMarker=-1
   wff%kgwff=-1
#endif
 end if

end subroutine WffDelete
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffKg
!! NAME
!! WffKg
!!
!! FUNCTION
!! Check kgwff to  manage WF file in the MPI/IO case
!!
!! INPUTS
!!  wff <type(wffile_type)> = structured info about the wavefunction file
!!  optkg= if 1 , read or write kg_k ; if 0,do not care about kg_k in rwwf
!!
!! OUTPUT
!!
!! PARENTS
!!      m_iowf,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffKg(wff,optkg)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: optkg

! *********************************************************************

#if defined HAVE_MPI_IO
 if (wff%iomode == IO_MODE_MPI) wff%kgwff=optkg
#else
 ABI_UNUSED((/wff%iomode,optkg/))
#endif

end subroutine WffKg
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/WffOffset
!! NAME
!! WffOffset
!!
!! FUNCTION
!! Tool to manage WF file in the MPI/IO case : broadcast the offset of
!! the first k-point data block
!!
!! INPUTS
!!  wff <type(wffile_type)> = structured info about the wavefunction file
!!  sender = id of the sender
!!  spaceComm = id of the space communicator handler
!!
!! OUTPUT
!!  ier = error code returned by the MPI call
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffOffset(wff,sender,spaceComm,ier)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer          ,intent(inout) :: sender
 integer          ,intent(in)    :: spaceComm
 integer          ,intent(out)   :: ier

!Local variables ------------------------------
#if defined HAVE_MPI_IO
 integer :: icom
 integer(kind=MPI_OFFSET_KIND) :: ima
#endif

! *********************************************************************

#if defined HAVE_MPI_IO
 if (wff%iomode == IO_MODE_MPI) then
   call xmpi_max(sender,icom,spaceComm,ier)
   if (icom>=0)then
     ima=wff%offwff
     call MPI_BCAST(ima,1,wff%offset_mpi_type,icom,spaceComm,ier)
     wff%offwff=ima
   end if
 end if ! iomode
#else
 ier = 0
 ABI_UNUSED((/wff%iomode,sender,spaceComm/))
#endif

end subroutine WffOffset
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffReadDataRec_dp1d
!! NAME
!! WffReadDataRec_dp1d
!!
!! FUNCTION
!! Subroutine to read data in one record of a wavefunction file
!! Handles double precision 1D arrays
!!
!! INPUTS
!! ndp=size of the double precision array to be read
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! dparray=array of double precision numbers
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffReadDataRec_dp1d(dparray,ierr,ndp,wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  ndp
 integer,intent(out) :: ierr
 real(dp),intent(out) :: dparray(ndp)

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   read (wff%unwff,iostat=ierr) dparray(1:ndp)

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   call xderiveRRecInit(wff,ierr)
   call xderiveRead(wff,dparray,ndp,MPI_COMM_SELF,ierr)
   call xderiveRRecEnd(wff,ierr)
#endif
 else
   write(msg,'(a,i0)')"Wrong iomode: ",wff%iomode
   MSG_ERROR(msg)
 end if

end subroutine WffReadDataRec_dp1d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/WffReadDataRec_dp2d
!! NAME
!! WffReadDataRec_dp2d
!!
!! FUNCTION
!! Subroutine to read data in one record of a wavefunction file
!! Handles double precision 2D arrays
!!
!! INPUTS
!! n1,n2=sizes of the double precision array to be read
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! dparray=array of double precision numbers
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffReadDataRec_dp2d(dparray,ierr,n1,n2,wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  n1,n2
 integer,intent(out) :: ierr
 real(dp),intent(out) :: dparray(n1,n2)

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   read (wff%unwff,iostat=ierr) dparray(1:n1,1:n2)

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   call xderiveRRecInit(wff,ierr)
   call xderiveRead(wff,dparray,n1,n2,MPI_COMM_SELF,ierr)
   call xderiveRRecEnd(wff,ierr)
#endif
 else
   write(msg,'(a,i0)')"Wrong iomode: ",wff%iomode
   MSG_ERROR(msg)
 end if

end subroutine WffReadDataRec_dp2d
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffReadNpwRec
!! NAME
!! WffReadNpwRec
!!
!! FUNCTION
!! This subroutine read the npw record of a wavefunction file
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!  wff%access == -1 and wf%master == Wff%me:
!!     read binary data
!!  wff%iomode == 0:
!!     read binary data
!!  wff%iomode == 1:
!!     use MPI/IO routines (MPIO defined)
!!  wff%iomode == 2:
!!     read netcdf format (NETCDF defined)
!! ikpt= the i-th kpoint.
!! isppol= the given spin polarisation element.
!!
!! OUTPUT
!! ierr=error code (iostat integer from read statement)
!! nband_disk=number of bands
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffReadNpwRec(ierr,ikpt,isppol,nband_disk,npw,nspinor,wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in)  :: ikpt, isppol
 integer,intent(out) :: ierr,nband_disk,npw,nspinor

!Local variables-------------------------------
 !character(len=500) :: msg
#if defined HAVE_NETCDF
 integer :: vid
#endif

! *************************************************************************

 ierr=0

 if (wff%iomode == IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me) ) then
   read (wff%unwff,iostat=ierr) npw,nspinor,nband_disk

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   call xderiveRRecInit(wff,ierr)
   call xderiveRead(wff,npw,ierr)
   call xderiveRead(wff,nspinor,ierr)
   call xderiveRead(wff,nband_disk,ierr)
   call xderiveRRecEnd(wff,ierr)
#endif

 else if (wff%iomode == IO_MODE_ETSF) then

#if defined HAVE_NETCDF
   !write(std_out,*)"readnpwrec: ikpt, spin", ikpt, spin
   NCF_CHECK(nctk_get_dim(wff%unwff, "number_of_spinor_components", nspinor))
   vid = nctk_idname(wff%unwff, "number_of_coefficients")
   NCF_CHECK(nf90_get_var(wff%unwff, vid, npw, start=[ikpt]))
   vid = nctk_idname(wff%unwff, "number_of_states")
   NCF_CHECK(nf90_get_var(wff%unwff, vid, nband_disk, start=[ikpt, isppol]))
#endif

 else
   ! MG: I don't understand why we have to use this ugly code!!!!!!!!
   ! Only master knows npw,nspinor,nband_disk in IO_MODE_FORTRAN_MASTE mode
   ! To the person who wrote this stuff:
   ! Have you ever heard about the "IF" statement of Fortran and the typical construct
   !
   !      if (rank==master) call mpifoo_seq()

   MSG_WARNING("Skipping read in WffReadNpwRec. Keep fingers crossed")
   ! MG: Must initialze these values somehow to avoid overflows.
   npw = 0; nspinor = 0; nband_disk = 0
 end if

 !write(std_out,*)"nband_disk,npw,nspinor",nband_disk,npw,nspinor
 ABI_CHECK(ierr==0,"ierr!=0")

end subroutine WffReadNpwRec
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffReadSkipRec
!! NAME
!! WffReadSkipRec
!!
!! FUNCTION
!! This subroutine move forward or backward in a Wf file by nrec records.
!!
!! INPUTS
!! nrec=number of records
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! TODO
!! For the future : one should treat the possible errors of backspace
!!
!! PARENTS
!!      gstate,randac,rwwf
!!
!! CHILDREN
!!
!! SOURCE


subroutine WffReadSkipRec(ierr,nrec,wff)

!Arguments ------------------------------------
 integer,intent(in)  :: nrec
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: irec
 integer(kind=MPI_OFFSET_KIND) :: delim_record,offset
#endif

! *************************************************************************

 ierr=0
 if( wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then

   call mvrecord(wff%unwff,nrec,ierr)
   ABI_CHECK(ierr==0,"error in mvrecord")


 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   if (nrec>0) then ! Move forward nrec records
     do irec=1,nrec
       wff%off_recs = wff%offwff
       call rwRecordMarker(1,wff%offwff,delim_record,wff,ierr)
       wff%lght_recs = delim_record
     end do
   else             ! Move backward -nrec records
     do irec=1,-nrec
       offset = wff%offwff-wff%nbOct_recMarker
       call rwRecordMarker(1,offset,delim_record,wff,ierr)
       wff%lght_recs = delim_record
       wff%offwff = wff%offwff - delim_record - 2*wff%nbOct_recMarker
       wff%off_recs = wff%offwff
     end do
   end if
#endif
 end if ! wff%iomode==0,1 or -1

end subroutine WffReadSkipRec
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffReadWrite_mpio
!! NAME
!! WffReadWrite_mpio
!!
!! FUNCTION
!!  This procedure read or write cg in the file _WFK using MPI_IO
!!  when cg are dispatched amoung commcart communicator
!!
!! INPUTS
!!  wff=struct info for wavefunction
!!  nband_disk=number of bands on disk files to be write
!!  icg=shift to be given to the location of the cg array
!!  mcg=second dimention of cg
!!  mpi_enreg=information about parallelisation
!!  depl_mpi_to_seq=for each proc, index of cg in sequential mode
!!  npwso=npw*nspinor number of plane waves treated by this node.
!!  npwsotot=npwtot*nspinor Total number of planewaves Used to calculate the size of data to be written.
!!  rdwr=1 if reading, 2 if writing
!!
!! OUTPUT
!!  ierr=error status
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions,
!!
!! NOTES
!!  cg is written like the following:
!!    BeginMarker cg ( iband = 1 )  EndMarker
!!    BeginMarker cg ( iband = 2 )  EndMarker
!!    ...
!!    BeginMarker cg( iband = nband_disk ) EndMarker
!!
!!  BeginMarker and EndMarker give the value of the total length of cg for one band
!!
!!  For MPI-IO library the performance is improved by the use a "view" of the file for each proc.

!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffReadWrite_mpio(wff,rdwr,cg,mcg,icg,nband_disk,npwso,npwsotot,depl_mpi_to_seq,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,mcg,nband_disk,npwso,npwsotot,rdwr
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff
!arrays
 integer,intent(in) :: depl_mpi_to_seq(npwso)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
#if defined HAVE_MPI_IO
 integer,parameter :: MAXBAND=500, check_markers=1
 integer :: filetype,iband,ibandmax,ibandmin,iblock,ii,iloc,ipw,jerr,jj
 integer :: nb,nband_block,step,totsize1bandcg,wfftempo
 integer(kind=MPI_OFFSET_KIND) :: delim_record,loc_depl_band,offset,totsize1bandByte
 character(len=500) :: msg
!arrays
 integer,allocatable :: BlockLength(:),BlockType(:),map(:),tempo_map(:)
 integer(kind=MPI_OFFSET_KIND),allocatable :: BlockDepl(:)
 integer(kind=2),allocatable :: bufdelim2(:)
 integer(kind=4),allocatable :: bufdelim4(:)
 integer(kind=8),allocatable :: bufdelim8(:)
 real(dp),allocatable :: buf(:),tempo_buf(:)
#if defined HAVE_FC_INT_QUAD
 integer(kind=16),allocatable :: bufdelim16(:)
#endif
#endif

! *********************************************************************

 ierr=0

#if defined HAVE_MPI_IO
!----------------------------------------------
!! Prepare WF data
!----------------------------------------------
!Init offset of record
 wff%off_recs = wff%offwff

!Total size to be written (in number of bands and in bytes)
 totsize1bandcg=2*npwsotot
!call xmpi_sum(totsize1bandcg,wff%spaceComm_mpiio,ierr)

 totsize1bandByte=totsize1bandcg*wff%nbOct_dp+2*wff%nbOct_recMarker

!Check file size
 offset=wff%offwff+nband_disk*totsize1bandByte
 if (offset>Huge(offset)) then
   msg='File is too large for MPI-IO specifications !'
   MSG_ERROR(msg)
 end if

!Open file
 call MPI_FILE_OPEN(wff%spaceComm_mpiio,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 ABI_CHECK_MPI(ierr, sjoin("MPI_FILE_OPEN:", wff%fname))

!----------------------------------------------------------
!Loop blocks of bands (to decrease offsets inside the file)
!----------------------------------------------------------
 ibandmax=0;ibandmin=1
 ii=huge(check_markers)/totsize1bandByte;step=min(ii,MAXBAND,nband_disk)
 do iblock=1,nband_disk/step+1
   ibandmax=min(ibandmin+step-1,nband_disk)
   nband_block=ibandmax-ibandmin+1
   offset=wff%offwff+(ibandmin-1)*totsize1bandByte

!  ----------------------------------------------
!  Read/Write bands
!  ----------------------------------------------

!  Build map; for better performance, map must be in increasing order
   ABI_MALLOC_OR_DIE(map,(2*npwso*nband_block), ierr)

   ABI_STAT_ALLOCATE(buf,(2*npwso*nband_block), ierr)
   ABI_CHECK(ierr==0, "out of memory in wavefunction buffer. Try to decrease MAXBAND in WffReadWrite_mpio")

   if (rdwr==1) then
!    If reading, only build map
     nb=0;loc_depl_band=0
     ABI_ALLOCATE(tempo_map,(2*npwso))
     do iband=ibandmin,ibandmax
       tempo_map(1:2*npwso)=-1
       jj=1;ipw=(iband-1)*npwso+icg
       do ii=1,npwso
         iloc=loc_depl_band+wff%nbOct_recMarker+2*(depl_mpi_to_seq(ii)-1)*wff%nbOct_dp
         tempo_map(jj  )=iloc              ! Real part
         tempo_map(jj+1)=iloc+wff%nbOct_dp ! Imag part
         jj=jj+2
       end do
       do ii=1,2*npwso ! Now, elimate holes
         if (tempo_map(ii)/=-1) then
           nb=nb+1
           map(nb)=tempo_map(ii)
         end if
       end do
       loc_depl_band=loc_depl_band+totsize1bandByte ! Location in bytes
     end do
   else if (rdwr==2) then
!    If writing, build map and store cg in a buffer
     nb=0;loc_depl_band=0
     ABI_ALLOCATE(tempo_map,(2*npwso))
     ABI_ALLOCATE(tempo_buf,(2*npwso))
     do iband=ibandmin,ibandmax
       tempo_map(1:2*npwso)=-1
       jj=1;ipw=(iband-1)*npwso+icg
       do ii=1,npwso
         iloc=loc_depl_band+wff%nbOct_recMarker+2*(depl_mpi_to_seq(ii)-1)*wff%nbOct_dp
         tempo_map(jj  )=iloc              ! Real part
         tempo_map(jj+1)=iloc+wff%nbOct_dp ! Imag part
         tempo_buf(jj:jj+1)=cg(1:2,ipw+ii)
         jj=jj+2
       end do
       do ii=1,2*npwso ! Now, elimate holes
         if (tempo_map(ii)/=-1) then
           nb=nb+1
           map(nb)=tempo_map(ii)
           buf(nb)=tempo_buf(ii)
         end if
       end do
       loc_depl_band=loc_depl_band+totsize1bandByte ! Location in bytes
     end do
     ABI_DEALLOCATE(tempo_map)
     ABI_DEALLOCATE(tempo_buf)
   end if  ! rdwr

!  Build and commit MPI datatype
   ABI_ALLOCATE(BlockLength,(nb+2))
   ABI_ALLOCATE(BlockDepl,(nb+2))
   ABI_ALLOCATE(BlockType,(nb+2))
   BlockLength(1)=1;BlockDepl(1)=0;BlockType(1)=MPI_LB
   do ii=2,nb+1
     BlockLength(ii)=1
     BlockDepl(ii)=map(ii-1)
     BlockType(ii)=MPI_DOUBLE_PRECISION
   end do
   BlockLength(nb+2)=1;BlockDepl(nb+2)=totsize1bandByte*nband_block;BlockType(nb+2)=MPI_UB
   call xmpio_type_struct(nb+2,BlockLength,BlockDepl,BlockType,filetype,ierr)
   call MPI_TYPE_COMMIT(filetype,ierr)
   ABI_DEALLOCATE(BlockLength)
   ABI_DEALLOCATE(BlockDepl)
   ABI_DEALLOCATE(BlockType)

!  Read/Write data on disk
   call MPI_FILE_SET_VIEW(wfftempo,offset,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
   if (rdwr==1) then
     call MPI_FILE_READ_ALL (wfftempo,buf,nb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
   else
     call MPI_FILE_WRITE_ALL(wfftempo,buf,nb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
   end if

!  In case of reading, retrieve cg
   if (rdwr==1) then
     nb=0;loc_depl_band=0
     ABI_ALLOCATE(tempo_buf,(2*npwso))
     do iband=ibandmin,ibandmax
       do ii=1,2*npwso ! Now, elimate holes
         if (tempo_map(ii)/=-1) then
           nb=nb+1;tempo_buf(ii)=buf(nb)
         end if
       end do
       jj=1;ipw=(iband-1)*npwso+icg
       do ii=1,npwso
         iloc=loc_depl_band+wff%nbOct_recMarker+2*(depl_mpi_to_seq(ii)-1)*wff%nbOct_dp
         cg(1:2,ipw+ii)=tempo_buf(jj:jj+1)
         jj=jj+2
       end do
       loc_depl_band=loc_depl_band+totsize1bandByte ! Location in bytes
     end do
     ABI_DEALLOCATE(tempo_map)
     ABI_DEALLOCATE(tempo_buf)
   end if ! rdwr

!  Free memory
   ABI_DEALLOCATE(map)
   ABI_DEALLOCATE(buf)
   call MPI_TYPE_FREE(filetype,ierr)

!  ----------------------------------------------
!  Check/Write record markers (only master proc)
!  ----------------------------------------------
   if ((rdwr==1.and.check_markers==1).or.(rdwr==2)) then

!    Define view for the file
     nb=2*nband_block
     ABI_ALLOCATE(BlockLength,(nb+2))
     ABI_ALLOCATE(BlockDepl,(nb+2))
     ABI_ALLOCATE(BlockType,(nb+2))
     BlockLength(1)=1;BlockDepl(1)=0;BlockType(1)=MPI_LB
     jj=2
     do ii=1,nband_block
       BlockType(jj:jj+1)  =wff%marker_mpi_type
       BlockLength(jj:jj+1)=1
       BlockDepl(jj  )=(ii-1)*totsize1bandByte
       BlockDepl(jj+1)= ii   *totsize1bandByte-wff%nbOct_recMarker
       jj=jj+2
     end do
     BlockLength(nb+2)=1;BlockDepl(nb+2)=nband_block*totsize1bandByte;BlockType(nb+2)=MPI_UB
     call xmpio_type_struct(nb+2,BlockLength,BlockDepl,BlockType,filetype,ierr)
     call MPI_TYPE_COMMIT(filetype,ierr)
     call MPI_FILE_SET_VIEW(wfftempo,offset,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
     ABI_DEALLOCATE(BlockLength)
     ABI_DEALLOCATE(BlockDepl)
     ABI_DEALLOCATE(BlockType)

!    Read/Write all markers (depend on Fortran marker MPI type)
     if (wff%me_mpiio==0) then
       jerr=0;delim_record=totsize1bandByte-2*wff%nbOct_recMarker
       if (wff%nbOct_recMarker==4) then
         ABI_ALLOCATE(bufdelim4,(nb))
         if (rdwr==2) bufdelim4(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim4,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim4(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim4,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim4)
       else if (wff%nbOct_recMarker==8) then
         ABI_ALLOCATE(bufdelim8,(nb))
         if (rdwr==2) bufdelim8(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim8,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim8(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim8,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim8)
#if defined HAVE_FC_INT_QUAD
       else if (wff%nbOct_recMarker==16) then
         ABI_ALLOCATE(bufdelim16,(nb))
         if (rdwr==2) bufdelim16(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim16,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim16(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim16,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim16)
#endif
       else if (wff%nbOct_recMarker==2) then
         ABI_ALLOCATE(bufdelim2,(nb))
         if (rdwr==2) bufdelim2(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim2,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim2(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim2,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim2)
       end if
       if (rdwr==1.and.jerr==1) then
         write(unit=msg,fmt='(2a)') 'Error when reading record markers of file ',trim(wff%fname)
         MSG_ERROR(msg)
       end if
     end if  ! me_mpiio=0

!    Free memory
     call MPI_TYPE_FREE(filetype,ierr)

   end if ! rdwr

!  -----------------------------------------
!  End loop on blocks of bands
!  -----------------------------------------
   if (ibandmax>=nband_disk) exit
   ibandmin=ibandmax+1
 end do

!-----------------------------------------
!End statements
!-----------------------------------------
!Close file
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Update offset
 wff%offwff=wff%offwff+totsize1bandByte*nband_disk
#endif

#if !defined HAVE_MPI_IO
!Dummy check to avoid warning from compilers.
 ABI_UNUSED((/wff%iomode,rdwr,size(cg),mcg,icg,nband_disk,npwso,depl_mpi_to_seq(1),npwsotot/))
#endif

end subroutine WffReadWrite_mpio
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/WffWriteDataRec_int2d
!! NAME
!! WffWriteDataRec_int2d
!!
!! FUNCTION
!! Subroutine to write data in one record of a wavefunction file
!! Handles integer 2D arrays
!!
!! INPUTS
!! intarray=array of integer numbers
!! n1,n2=sizes of the integer array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffWriteDataRec_int2d(intarray,ierr,n1,n2,wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  n1,n2
 integer,intent(out) :: ierr
 integer,intent(in) :: intarray(n1,n2)

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) intarray(1:n1,1:n2)

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   call xderiveWRecInit(wff,ierr)
   call xderiveWrite(wff,intarray,n1,n2,MPI_COMM_SELF,ierr)
   call xderiveWRecEnd(wff,ierr)
#endif
 else
   write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
   MSG_WARNING(msg)
 end if

end subroutine WffWriteDataRec_int2d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/WffWriteDataRec_dp1d
!! NAME
!! WffWriteDataRec_dp1d
!!
!! FUNCTION
!! Subroutine to write data in one record of a wavefunction file
!! Handles double precision 1D arrays
!!
!! INPUTS
!! dparray=array of double precision numbers
!! ndp=size of the double precision array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine WffWriteDataRec_dp1d(dparray,ierr,ndp,wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  ndp
 integer,intent(out) :: ierr
 real(dp),intent(in) :: dparray(ndp)

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) dparray(1:ndp)

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   call xderiveWRecInit(wff,ierr)
   call xderiveWrite(wff,dparray,ndp,MPI_COMM_SELF,ierr)
   call xderiveWRecEnd(wff,ierr)
#endif
 else
   write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
   MSG_WARNING(msg)
 end if

end subroutine WffWriteDataRec_dp1d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/WffWriteDataRec_dp2d
!! NAME
!! WffWriteDataRec_dp2d
!!
!! FUNCTION
!! Subroutine to write data in one record of a wavefunction file
!! Handles double precision 2D arrays
!!
!! INPUTS
!! dparray=array of double precision numbers
!! n1,n2=sizes of the double precision array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


subroutine WffWriteDataRec_dp2d(dparray,ierr,n1,n2,wff)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  n1,n2
 integer,intent(out) :: ierr
 real(dp),intent(in) :: dparray(n1,n2)

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) dparray(1:n1,1:n2)

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   call xderiveWRecInit(wff,ierr)
   call xderiveWrite(wff,dparray,n1,n2,MPI_COMM_SELF,ierr)
   call xderiveWRecEnd(wff,ierr)
#endif
 else
   write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
   MSG_WARNING(msg)
 end if

end subroutine WffWriteDataRec_dp2d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/WffWriteNpwRec
!! NAME
!! WffWriteNpwRec
!!
!! FUNCTION
!! This subroutine writes the npw record of a wavefunction file
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!! nband_disk=number of bands
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions
!! opt_paral=(optional argument, default=1, only used for MPI-IO)
!!           1: all procs in the communicator write the data
!!           2: only master in the communicator writes the data
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!
!! SOURCE


subroutine WffWriteNpwRec(ierr,nband_disk,npw,nspinor,wff,&
&                         opt_paral) ! optional argument

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: nband_disk,npw,nspinor
 integer,intent(in),optional :: opt_paral
 integer,intent(out) :: ierr

!Local variables-------------------------------
 integer :: opt_paral_
#if defined HAVE_MPI_IO
 integer :: me
#endif
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 opt_paral_=1;if (present(opt_paral)) opt_paral_=opt_paral

 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode ==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) npw,nspinor,nband_disk

 else if(wff%iomode==IO_MODE_MPI)then
#if defined HAVE_MPI_IO
   me=-1;if (opt_paral_==2) me=wff%me_mpiio
   if ((me==-1.and.opt_paral_==1).or.(me==0.and.opt_paral_==2)) then
     call xderiveWRecInit(wff,ierr)
     call xderiveWrite(wff,npw,ierr)
     call xderiveWrite(wff,nspinor,ierr)
     call xderiveWrite(wff,nband_disk,ierr)
     call xderiveWRecEnd(wff,ierr)
   end if
   if (opt_paral_==2.and.wff%spaceComm_mpiio/=MPI_COMM_SELF) then
     call xmpi_barrier(wff%spaceComm_mpiio)
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,0,wff%spaceComm_mpiio,ierr)
   end if
#endif
 else
   write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
   MSG_WARNING(msg)
 end if

end subroutine WffWriteNpwRec
!!***

!!****f* m_wffile/xderiveRead_int
!! NAME
!!  xderiveRead_int
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: integer scalar.
!!
!! INPUTS
!! (none)
!!
!! OUTPUT
!!  xval= data buffer
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_int(wff,xval,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
 integer :: tmparr(1)
#endif
 character(len=500) :: msg

! *********************************************************************

 xval=0; ierr=0

#if defined HAVE_MPI_IO
 call MPI_FILE_READ_AT(wff%fhwff,wff%offwff,tmparr,1,MPI_INTEGER,statux,ierr)
 xval = tmparr(1)
 wff%offwff = wff%offwff + wff%nbOct_int
 RETURN
#endif

 ABI_UNUSED(wff%me)

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)

end subroutine xderiveRead_int
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveRead_int1d
!! NAME
!!  xderiveRead_int1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_int1d(wff,xval,n1,spaceComm,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:)
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,nboct,dispoct,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif
 character(len=500) :: msg

! *********************************************************************

 xval(:)=0 ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_int * n1
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit+dispoct-nboct
   end if
   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1,MPI_INTEGER,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call xmpi_max(dispoct,totoct,spaceComm,ierr)
     !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff = wff%offwff + totoct
 return
#endif

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)

end subroutine xderiveRead_int1d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveRead_int2d
!! NAME
!!  xderiveRead_int2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_int2d(wff,xval,n1,n2,spaceComm,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:,:)
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,dispoct,nboct,posit,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif
 character(len=500) :: msg

! *********************************************************************

 xval(:,:)=0 ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_int * n1 * n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit + dispoct - nboct
   end if
   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call xmpi_max(dispoct,totoct,spaceComm,ierr)
     !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff=wff%offwff + totoct
 return
#endif

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)

end subroutine xderiveRead_int2d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveRead_dp
!! NAME
!!  xderiveRead_dp
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision scalar.
!!
!! INPUTS
!! (none)
!!
!! OUTPUT
!!  xval= data buffer
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_dp(wff,xval,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
 real(dp) :: tmparr(1)
#endif
 character(len=500) :: msg

! *********************************************************************

 xval=zero ; ierr=0
 if(.false.)write(std_out,*)wff%me
#if defined HAVE_MPI_IO
 call MPI_FILE_READ_AT(wff%fhwff,wff%offwff,tmparr,1,MPI_DOUBLE_PRECISION,statux,ierr)
 xval = tmparr(1)
 wff%offwff = wff%offwff + wff%nbOct_dp
 return
#endif

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)

end subroutine xderiveRead_dp
!!***

!----------------------------------------------------------------------


!!****f* ABINIT/xderiveRead_dp1d
!! NAME
!!  xderiveRead_dp1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine xderiveRead_dp1d(wff,xval,n1,spaceComm,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,nboct,dispoct,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif
 character(len=500) :: msg

!*********************************************************************

 xval(:)=zero ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_dp * n1
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit + dispoct - nboct
   end if

   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1,MPI_DOUBLE_PRECISION,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call xmpi_max(dispoct,totoct,spaceComm,ierr)
     !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff=wff%offwff + totoct
 return
#endif

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)

end subroutine xderiveRead_dp1d
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveRead_dp2d
!! NAME
!!  xderiveRead_dp2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_dp2d(wff,xval,n1,n2,spaceComm,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:,:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,dispoct,nboct,posit,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif
 character(len=500) :: msg

! *********************************************************************

 xval(:,:)=zero ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_dp * n1 *n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit + dispoct - nboct
   end if
   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1*n2,MPI_DOUBLE_PRECISION,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call xmpi_max(dispoct,totoct,spaceComm,ierr)
     !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff=wff%offwff + totoct
 return
#endif

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)

end subroutine xderiveRead_dp2d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveRead_int2d_displ
!! NAME
!!  xderiveRead_int2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_int2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(out):: xval(:,:)
 integer,intent(in):: displace(:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: buf_val(:),length1(:),type1(:),val(:)
 integer(kind=MPI_OFFSET_KIND),allocatable :: depl(:),depl1(:),depl_val(:)
#endif
 character(len=500) :: msg

! *********************************************************************

 xval(:,:)=0 ; ierr=0
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm,displace

#if defined HAVE_MPI_IO
 nb=n1*n2
 call xmpi_sum(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     depl_val(ipos)=ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     depl(nbval)=depl_val(i1)
   end if
 end do

!Build MPI datatype for view
 ABI_ALLOCATE(length1,(nbval+2))
 ABI_ALLOCATE(depl1,(nbval+2))
 ABI_ALLOCATE(type1,(nbval+2))
 length1(1)=1;depl1(1)=0;type1(1)=MPI_LB
 do i1=2,nbval+1
   length1(i1) = 1
   depl1(i1)= depl(i1-1)*wff%nbOct_int
   type1(i1)= MPI_INTEGER
 end do
 length1(nbval+2)=1;depl1(nbval+2)=totsize*wff%nbOct_int;type1(nbval+2)=MPI_UB
 call xmpio_type_struct(nbval+2,length1,depl1,type1,filetype,ierr)
 call MPI_TYPE_COMMIT(filetype,ierr)
 ABI_DEALLOCATE(length1)
 ABI_DEALLOCATE(depl1)
 ABI_DEALLOCATE(type1)

!Write data
 call MPI_FILE_OPEN(spaceComm,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 ABI_CHECK_MPI(ierr, sjoin("MPI_FILE_OPEN:", wff%fname))
 call MPI_FILE_SET_VIEW(wfftempo,wff%offwff,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
 call MPI_FILE_READ_ALL(wfftempo,val,nbval,MPI_INTEGER,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Retrieve xval
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     buf_val(i1)=val(nbval)
   end if
 end do
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     xval(i1,i2)=buf_val(ipos)
   end do
 end do

!Update offset
 wff%offwff = wff%offwff + totsize*wff%nbOct_int

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
 return
#endif

 write(msg,'(a,i0,a)')' The value of wff%iomode=',wff%iomode,' is not allowed.'
 MSG_WARNING(msg)


end subroutine xderiveRead_int2d_displ
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveRead_dp2d_displ
!! NAME
!!  xderiveRead_dp2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveRead_dp2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out):: xval(:,:)
 integer,intent(in):: displace(:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: length1(:),type1(:)
 integer(kind=MPI_OFFSET_KIND),allocatable :: depl(:),depl1(:),depl_val(:)
 real(dp), allocatable :: buf_val(:),val(:)
#endif

! *********************************************************************

 xval(:,:)=zero ; ierr=0
 if(.false.)write(std_out,*)wff%me,n1,n2,displace,spaceComm

#if defined HAVE_MPI_IO
 nb=n1*n2
 call xmpi_sum(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     depl_val(ipos)=ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     depl(nbval)=depl_val(i1)
   end if
 end do

!Build MPI datatype for view
 ABI_ALLOCATE(length1,(nbval+2))
 ABI_ALLOCATE(depl1,(nbval+2))
 ABI_ALLOCATE(type1,(nbval+2))
 length1(1)=1;depl1(1)=0;type1(1)=MPI_LB
 do i1=2,nbval+1
   length1(i1) = 1
   depl1(i1)= depl(i1-1)*wff%nbOct_dp
   type1(i1)= MPI_DOUBLE_PRECISION
 end do
 length1(nbval+2)=1;depl1(nbval+2)=totsize*wff%nbOct_dp;type1(nbval+2)=MPI_UB
 call xmpio_type_struct(nbval+2,length1,depl1,type1,filetype,ierr)
 call MPI_TYPE_COMMIT(filetype,ierr)
 ABI_DEALLOCATE(length1)
 ABI_DEALLOCATE(depl1)
 ABI_DEALLOCATE(type1)

!Write data
 call MPI_FILE_OPEN(spaceComm,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 ABI_CHECK_MPI(ierr, sjoin("MPI_FILE_OPEN:", wff%fname))
 call MPI_FILE_SET_VIEW(wfftempo,wff%offwff,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
 call MPI_FILE_READ_ALL(wfftempo,val,nbval,MPI_DOUBLE_PRECISION,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Retrieve xval
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     buf_val(i1)=val(nbval)
   end if
 end do
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     xval(i1,i2)=buf_val(ipos)
   end do
 end do

!Update offset
 wff%offwff = wff%offwff + totsize*wff%nbOct_dp

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveRead_dp2d_displ
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveReadVal_char
!! NAME
!!  xderiveReadVal_char
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: character string.
!!
!! INPUTS
!!  n= number of elements in the array
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveReadVal_char(wff,xval,n,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n
 integer,intent(out) :: ierr
 character(len=*),intent(out) :: xval

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
 character(len=len(xval)) :: tmparr(1)
#endif

! *********************************************************************

 xval=' ' ; ierr=0
 if(.false.)write(std_out,*)wff%me,n

#if defined HAVE_MPI_IO
 call MPI_FILE_READ_AT(wff%fhwff,wff%offwff,tmparr,n,MPI_CHARACTER,statux,ierr)
 xval = tmparr(1)
 wff%offwff = wff%offwff + wff%nbOct_ch * n
#endif

end subroutine xderiveReadVal_char
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xmpi_read_int2d
!! NAME
!!  xmpi_read_int2d
!!
!! FUNCTION
!!  Generic routine to read arrays with MPI I/O.
!!  Target: integer two-dimensional arrays.
!!
!! INPUTS
!!  sc_mode=
!!    xmpio_single     ==> Local reading.
!!    xmpio_collective ==> Collective reading.
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmpi_read_int2d(wff,xval,spaceComm,sc_mode,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spaceComm,sc_mode
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff
!array
 integer,intent(out) :: xval(:,:)

!Local variables-------------------------------
 integer :: n1,n2
#ifdef HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,nboct,posit,totoct
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 n1 = SIZE(xval,DIM=1)
 n2 = SIZE(xval,DIM=2)

#ifdef HAVE_MPI_IO
 nboct = wff%nbOct_int * n1 *n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then

   select case (sc_mode)
     case (xmpio_single)
       call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)

     case (xmpio_collective)
       call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)

     case default
       write(msg,('(a,i0)'))" Wrong value for sc_mode: ",sc_mode
       MSG_ERROR(msg)
   end select

   totoct=nboct
 else
   write(msg,('(a,2(i0,1x))'))" delim_record < nboct: ",delim_record,nboct
   MSG_WARNING(msg)
   ierr=MPI_ERR_UNKNOWN
   totoct=0
 end if
!
!Increment the offset.
 wff%offwff=wff%offwff + totoct
#endif

 RETURN
 ABI_UNUSED(xval(1,1))
 ABI_UNUSED((/wff%me,spaceComm,sc_mode/))

end subroutine xmpi_read_int2d
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveWrite_int
!! NAME
!!  xderiveWrite_int
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: integer scalar.
!!
!! INPUTS
!!  xval= data buffer
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_int(wff,xval,ierr)

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 integer,intent(in):: xval
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: tmparr(1)
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval
#if defined HAVE_MPI_IO
 tmparr(1) = xval
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,tmparr,1,MPI_INTEGER,statux,ierr)
 wff%offwff = wff%offwff+wff%nbOct_int
#endif

end subroutine xderiveWrite_int
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveWrite_int1d
!! NAME
!!  xderiveWrite_int1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_int1d(wff,xval,n1,spaceComm,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: dispoct,nboct,posit,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,n1,spaceComm,xval
#if defined HAVE_MPI_IO
 nboct = n1*wff%nbOct_int
 posit = wff%offwff

!dispoct = sum (nboct, rank=0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit+dispoct-nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1,MPI_INTEGER,statux,ierr)
!gather the bigest offset

 if (spaceComm/=MPI_COMM_SELF) then
   call xmpi_max(dispoct,totoct,spaceComm,ierr)
   !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff+totoct

!Disable old code
#endif

end subroutine xderiveWrite_int1d
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveWrite_int2d
!! NAME
!!  xderiveWrite_int2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_int2d(wff,xval,n1,n2,spaceComm,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: dispoct,nboct,posit,totoct
 integer  :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm,xval
#if defined HAVE_MPI_IO
 nboct = n1*n2*wff%nbOct_int
 posit = wff%offwff

!dispoct = sum(nboct, rank=0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit + dispoct-nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)
!gather the biggest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call xmpi_max(dispoct,totoct,spaceComm,ierr)
   !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff + totoct
#endif

end subroutine xderiveWrite_int2d
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveWrite_dp
!! NAME
!!  xderiveWrite_dp
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision scalar.
!!
!! INPUTS
!!  xval= data buffer
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_dp(wff,xval,ierr)

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
 real(dp) :: tmparr(1)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval
#if defined HAVE_MPI_IO
 tmparr(1) = xval
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,tmparr,1,MPI_DOUBLE_PRECISION,statux,ierr)
 wff%offwff = wff%offwff+wff%nbOct_dp
#endif

end subroutine xderiveWrite_dp
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveWrite_dp1d
!! NAME
!!  xderiveWrite_dp1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_dp1d(wff,xval,n1,spaceComm,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: nboct,dispoct,totoct,posit
 integer  :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,n1,spaceComm,xval
#if defined HAVE_MPI_IO
 nboct = n1*wff%nbOct_dp
 posit = wff%offwff
!dispoct = sum (nboct, rank = 0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit + dispoct - nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1,MPI_DOUBLE_PRECISION,statux,ierr)
!Gather the biggest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call xmpi_max(dispoct,totoct,spaceComm,ierr)
   !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff + totoct
#endif

end subroutine xderiveWrite_dp1d
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveWrite_dp2d
!! NAME
!!  xderiveWrite_dp2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_dp2d(wff,xval,n1,n2,spaceComm,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: nboct,dispoct,totoct,posit
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n1,n2,spaceComm

#if defined HAVE_MPI_IO
 nboct = n1*n2*wff%nbOct_dp
 posit = wff%offwff
!dispoct = sum(nboct, rank=0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit+dispoct-nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1*n2,MPI_DOUBLE_PRECISION,statux,ierr)
 posit = posit+nboct
!gather the biggest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call xmpi_max(dispoct,totoct,spaceComm,ierr)
   !call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff+totoct
#endif

end subroutine xderiveWrite_dp2d
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveWrite_dp2d_seq
!! NAME
!!  xderiveWrite_dp2d_seq
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_dp2d_seq(wff,xval,ierr)

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: n1,n2
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval
#if defined HAVE_MPI_IO
 n1=size(xval,1);n2=size(xval,2)
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,xval,n1*n2,MPI_DOUBLE_PRECISION,statux,ierr)
 wff%offwff = wff%offwff+wff%nbOct_dp*n1*n2
#endif

end subroutine xderiveWrite_dp2d_seq
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveWrite_int2d_displ
!! NAME
!!  xderiveWrite_int2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_int2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: displace(:),xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer, allocatable :: buf_val(:),length1(:),type1(:),val(:)
 integer(kind=MPI_OFFSET_KIND),allocatable :: depl(:),depl1(:),depl_val(:)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n1,n2,spaceComm,displace

#if defined HAVE_MPI_IO
 nb = n1*n2
 call xmpi_sum(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
!Put xval in a buffer at its position
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
!    ipos location of xval(i1,i2) in the array associated with record to be written
     ipos=(displace(i2)-1)*n1 + i1-1
     buf_val(ipos) = xval(i1,i2)
     depl_val(ipos) = ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     val(nbval)=buf_val(i1)
     depl(nbval)=depl_val(i1)
   end if
 end do

!Build MPI datatype for view
 ABI_ALLOCATE(length1,(nbval+2))
 ABI_ALLOCATE(depl1,(nbval+2))
 ABI_ALLOCATE(type1,(nbval+2))
 length1(1)=1;depl1(1)=0;type1(1)=MPI_LB
 do i1=2,nbval+1
   length1(i1) = 1
   depl1(i1)= depl(i1-1)*wff%nbOct_int
   type1(i1)= MPI_INTEGER
 end do
 length1(nbval+2)=1;depl1(nbval+2)=totsize*wff%nbOct_int;type1(nbval+2)=MPI_UB
 call xmpio_type_struct(nbval+2,length1,depl1,type1,filetype,ierr)
 call MPI_TYPE_COMMIT(filetype,ierr)
 ABI_DEALLOCATE(length1)
 ABI_DEALLOCATE(depl1)
 ABI_DEALLOCATE(type1)

!Write data
 call MPI_FILE_OPEN(spaceComm,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 ABI_CHECK_MPI(ierr, sjoin("MPI_FILE_OPEN:", wff%fname))
 call MPI_FILE_SET_VIEW(wfftempo,wff%offwff,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
 call MPI_FILE_WRITE_ALL(wfftempo,val,nbval,MPI_INTEGER,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Update offset
 wff%offwff = wff%offwff + totsize*wff%nbOct_int

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveWrite_int2d_displ
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xderiveWrite_dp2d_displ
!! NAME
!!  xderiveWrite_dp2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_dp2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: displace(:)
 real(dp),intent(in) :: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer, allocatable :: length1(:),type1(:)
 integer(kind=MPI_OFFSET_KIND),allocatable :: depl(:),depl1(:),depl_val(:)
 real(dp),allocatable :: buf_val(:),val(:)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n1,n2,spaceComm,displace

#if defined HAVE_MPI_IO
 nb = n1*n2
 call xmpi_sum(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
!Put xval in a buffer at its position
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
!    ipos location of xval(i1,i2) in the array associated with record to be written
     ipos=(displace(i2)-1)*n1 + i1-1
     buf_val(ipos) = xval(i1,i2)
     depl_val(ipos) = ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     val(nbval)=buf_val(i1)
     depl(nbval)=depl_val(i1)
   end if
 end do

!Build MPI datatype for view
 ABI_ALLOCATE(length1,(nbval+2))
 ABI_ALLOCATE(depl1,(nbval+2))
 ABI_ALLOCATE(type1,(nbval+2))
 length1(1)=1;depl1(1)=0;type1(1)=MPI_LB
 do i1=2,nbval+1
   length1(i1) = 1
   depl1(i1)= depl(i1-1)*wff%nbOct_dp
   type1(i1)= MPI_DOUBLE_PRECISION
 end do
 length1(nbval+2)=1;depl1(nbval+2)=totsize*wff%nbOct_dp;type1(nbval+2)=MPI_UB
 call xmpio_type_struct(nbval+2,length1,depl1,type1,filetype,ierr)
 call MPI_TYPE_COMMIT(filetype,ierr)
 ABI_DEALLOCATE(length1)
 ABI_DEALLOCATE(depl1)
 ABI_DEALLOCATE(type1)

!Write data
 call MPI_FILE_OPEN(spaceComm,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 ABI_CHECK_MPI(ierr, sjoin("MPI_FILE_OPEN:", wff%fname))
 call MPI_FILE_SET_VIEW(wfftempo,wff%offwff,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
 call MPI_FILE_WRITE_ALL(wfftempo,val,nbval,MPI_DOUBLE_PRECISION,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

 wff%offwff = wff%offwff + totsize*wff%nbOct_dp

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveWrite_dp2d_displ
!!***

!----------------------------------------------------------------------


!!****f* m_wffile/xderiveWrite_char
!! NAME
!!  xderiveWrite_char
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: character string.
!!
!! INPUTS
!!  xval= data buffer array
!!  n= number of elements in the string
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xderiveWrite_char(wff,xval,n,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: xval

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n

#if defined HAVE_MPI_IO
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,xval,n,MPI_CHARACTER,statux,ierr)
 wff%offwff = wff%offwff + wff%nbOct_ch * n
#endif

end subroutine xderiveWrite_char
!!***

!!****f* m_wffile/xdefineOff
!! NAME
!!  xdefineOff
!!
!! FUNCTION
!!  In case of MPI I/O, defines the offset for each processor
!!
!! INPUTS
!!  formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of eigenvectors with
!!        random numbers, vector of eigenvalues, occupations are present)
!!   1 => respfn format (initialisation of eigenvectors with 0 s,
!!        hermitian matrix of eigenvalues)
!!  nkpt = number of k points
!!  nspinor = total number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  nband(nkpt*nsppol) = number of bands at each k point, for each polarization
!!  npwarr(nkpt) = number of planewaves at each k point
!!  mpi_enreg <type(MPI_type)> = information about MPI parallelization
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  wff <type(wffile_type)> =
!!
!! PARENTS
!!      uderiv
!!
!! CHILDREN
!!
!! SOURCE

subroutine xdefineOff(formeig,wff,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)

!Arguments ------------------------------------
 integer, intent(in) ::  nsppol,nkpt,nspinor,formeig
 integer, intent(in) ::  nband(nkpt*nsppol),npwarr(nkpt)
 type(wffile_type),intent(inout) :: wff
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: comm,iproc
 integer :: nband_k,npw_k,nproc,me,ipp
 integer :: nbrec,isppol,ikpt,nbint,nbreal,nbd,ippband
 integer :: nrecnpw,nreckg
 integer(kind=XMPI_OFFSET_KIND) :: pos_start
!arrays
 integer(kind=XMPI_OFFSET_KIND),allocatable  :: offproc(:)
#endif

! *************************************************************************
!nbOct_int octet number of int value
!nbOct_dp octet number of dp value
!nbOct_ch octet number of character value
!lght_recs length of record

 if(.false.)write(std_out,*)wff%me,mpi_enreg%nproc,formeig,nband,npwarr,nspinor,nkpt
#if defined HAVE_MPI_IO
 if(wff%iomode==IO_MODE_MPI)then

   comm=mpi_enreg%comm_cell
   me=xmpi_comm_rank(comm)
   nproc=xmpi_comm_size(comm)
   pos_start=wff%offwff

   ABI_ALLOCATE(offproc,(0:nproc))
   offproc = 0
   nbrec =2
   nrecnpw=3+nbrec

   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       npw_k=npwarr(ikpt)
       iproc=mpi_enreg%proc_distrb(ikpt,1,isppol)
       if (mpi_enreg%paralbd==1) iproc=mpi_enreg%proc_distrb(ikpt,1,isppol)
!      record kg
       nreckg=nbrec+ wff%kgwff*3*npw_k

!      Record npw,nspinor,nband, Record kg
       offproc(iproc) = offproc(iproc) + wff%nbOct_int*(nrecnpw+nreckg)

       if (formeig == 0) then
!        Records eigen,occ
         nbint=nbrec
         nbreal =  2 *nband_k
         offproc(iproc) = offproc(iproc) + (wff%nbOct_int*nbint+wff%nbOct_dp*nbreal)

!        Records cg
         offproc(iproc) = offproc(iproc) &
&         + (wff%nbOct_int*nbrec+wff%nbOct_dp*2*npw_k*nspinor)*nband_k

         ippband=iproc
         do nbd=1,nband_k
           ipp=mpi_enreg%proc_distrb(ikpt,nbd,isppol)
           if (ipp /= ippband ) then
             ippband=ipp
             offproc(ippband)=offproc(ippband)+ wff%nbOct_int*(nrecnpw+nreckg)
             offproc(ippband) = offproc(ippband) + (wff%nbOct_int*nbint &
&             +wff%nbOct_dp*nbreal)
             offproc(ippband) = offproc(ippband) + (wff%nbOct_int*nbrec &
&             + wff%nbOct_dp*2*npw_k*nspinor)*nband_k
           end if
         end do
       else if (formeig == 1) then
!        record eigen
         offproc(iproc) = offproc(iproc) + (wff%nbOct_int*2*nbrec  &
&         + wff%nbOct_dp*2*npw_k*nspinor &
&         + wff%nbOct_dp*2*nband_k)*nband_k
         ippband=iproc
         do nbd=1,nband_k
           ipp=mpi_enreg%proc_distrb(ikpt,nbd,isppol)
           if (ipp /= ippband) then
             ippband=ipp
             offproc(ippband)=offproc(ippband)+ wff%nbOct_int*(nrecnpw+nreckg)
             offproc(ippband) = offproc(ippband) + (wff%nbOct_int*2*nbrec  &
&             + wff%nbOct_dp*2*npw_k*nspinor &
&             + wff%nbOct_dp*2*nband_k)*nband_k
           end if
         end do
       end if   ! formeig
     end do ! ikpt

   end do ! isppol

!  pos_start=wff%offwff
!  wff%offwff = pos_start

   if (me/=0)then
     do iproc=0,me-1
       wff%offwff=wff%offwff+offproc(iproc)
     end do
   end if
   ABI_DEALLOCATE(offproc)

 end if ! iomode
#endif

end subroutine xdefineOff
!!***

END MODULE m_wffile
!!***
