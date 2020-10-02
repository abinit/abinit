!!****m* ABINIT/m_mpiotk
!! NAME
!!  m_mpiotk
!!
!! FUNCTION
!!  This module provides helper functions for MPI-IO operations.
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
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_mpiotk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

 private

!public procedures.
#ifdef HAVE_MPI_IO
 public :: mpiotk_read_fsuba_dp2D     ! (individual|collective) read of a 2D sub-matrix
 public :: mpiotk_write_fsuba_dp2D    ! (individual|collective) write of a 2D sub-matrix

 public :: mpiotk_read_fsuba_dpc3D    ! (individual|collective) read of a 3D sub-matrix
 !public :: mpiotk_write_fsuba_dpc3D  ! (individual|collective) write of a 3D sub-matrix

 public :: mpiotk_read_fsuba_dpc4D    ! (individual|collective) read of a 4D sub-matrix
 !public :: mpiotk_write_fsuba_dpc4D  ! (individual|collective) write of a 4D sub-matrix
#else
 public :: no_mpiotk
#endif
!!***

CONTAINS
!!***

#ifndef HAVE_MPI_IO

!----------------------------------------------------------------------

!!****f* m_mpiotk/no_mpiotk
!! NAME
!!  no_mpiotk
!!
!! FUNCTION
!!   Empty placeholder.
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read,mpi_file_read_all,mpi_file_set_view,mpi_type_free
!!      xmpi_max,xmpi_min,xmpio_create_fsubarray_4d
!!
!! SOURCE

subroutine no_mpiotk()

! *************************************************************************

end subroutine no_mpiotk
!!***

#else

!!****f* m_mpiotk/setup_fsuba_dp2D
!! NAME
!!  setup_fsuba_dp2D
!!
!! FUNCTION
!!  Setup tables used in (read|write) a 2D array stored in a Fortran file
!!
!! NOTES
!!  The value of ierr should always be checked by the caller
!!
!! INPUTS
!!  sizes(2)
!!  subsizes(2)
!!  starts(2)
!!  chunk_bsize =
!!  comm = MPI communicator
!!
!! OUTPUTS
!!  my_basead(:)
!!  my_subsizes(:,:)
!!  my_starts(:,:)
!!  ierr=status error. A non zero value indicates that chunk_bsize is smaller that the fortran record
!!    and therefore bufsz has not been read.
!!
!! PARENTS
!!      m_mpiotk
!!
!! CHILDREN
!!      mpi_file_read,mpi_file_read_all,mpi_file_set_view,mpi_type_free
!!      xmpi_max,xmpi_min,xmpio_create_fsubarray_4d
!!
!! SOURCE

subroutine setup_fsuba_dp2D(sizes,subsizes,starts,chunk_bsize,&
&  my_basead,my_subsizes,my_starts,my_ncalls,ncalls,comm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 integer,intent(out) :: ncalls,my_ncalls,ierr
 integer(XMPI_OFFSET_KIND),intent(in) :: chunk_bsize
!arrays
 integer,intent(in) :: sizes(2),subsizes(2),starts(2)
 integer,allocatable,intent(out) :: my_basead(:),my_subsizes(:,:),my_starts(:,:)

!Local variables ------------------------------
!scalars
 integer :: mpierr,ny2read,ny_chunk,icall,yrest
 integer :: size_x,subs_x,start_x,start_y,stop_y
 !character(len=500) :: msg

!************************************************************************

 size_x  = sizes(1)
 subs_x  = subsizes(1)
 start_x = starts(1)
 start_y = starts(2)
 stop_y  = start_y + subsizes(2)-1 ! last column to read
 !
 ! Read rows in blocks of size ny_chunk:
 ! MPI-IO crashes if we try to read data > 2Gb in a single call.
 ny2read = subsizes(2)
 ny_chunk = ny2read
 if ((2*subs_x*ny2read*xmpi_bsize_dp) > chunk_bsize) then
   ny_chunk = chunk_bsize / (2*subs_x*xmpi_bsize_dp)
   !if (ny_chunk == 0) ny_chunk = 50
 end if

 call xmpi_min(ny_chunk,ierr,comm,mpierr)
 if (ierr == 0) then
   ierr = 1
   RETURN
 end if
 ierr = 0
 !
 ! my_ncalls : number of read needed to fill my buffer.
 ! ncalls    : max number of read in comm (needed for collective operations).
 my_ncalls = ny2read / ny_chunk
 yrest = MOD(ny2read, ny_chunk)
 if (yrest /= 0) my_ncalls = my_ncalls + 1

 call xmpi_max(my_ncalls,ncalls,comm,mpierr)
 !
 ! Compute arrays used to define the file view.
 ABI_MALLOC(my_subsizes,(2,ncalls))
 ABI_MALLOC(my_starts,(2,ncalls))
 ABI_MALLOC(my_basead,(ncalls))

 do icall=1,my_ncalls

   if (icall*ny_chunk <= ny2read) then
     my_subsizes(:,icall) = (/subs_x, ny_chunk/)
     my_starts(:,icall) = (/start_x, (icall-1) * ny_chunk + start_y/)
     my_basead(icall) = 1 + 2*(icall-1)*ny_chunk*subs_x ! 2 accounts for real and imag part.
   else
     ! Two cases:
     ! 1) ny2read > ny_chunk and not divisible by ny2read
     ! 2) ny2read < ny_chunk
     my_subsizes(:,icall) = (/subs_x, yrest/)
     if (ny2read >= ny_chunk) then
       my_starts(:,icall) = (/start_x, stop_y-yrest+1/)
       my_basead(icall) = 1 + 2 * (ny2read-yrest) * subs_x ! 2 accounts for real and imag part.
     else
       my_starts(:,icall) = starts
       my_basead(icall) = 1
     end if
   end if
 end do
 !write(std_out,*)" >>>> my_ncalls, ncalls, ny2read, ny_chunk ",my_ncalls,ncalls,ny2read,ny_chunk

end subroutine setup_fsuba_dp2D
!!***

!----------------------------------------------------------------------

!!****f* m_mpiotk/mpiotk_read_fsuba_dp2D
!! NAME
!!  mpiotk_read_fsuba_dp2D
!!
!! FUNCTION
!!  Read a block of contiguous data stored in a 2D matrix.
!!  Data is placed within Fortran records. Target: complex data stored in a real array.
!!
!! NOTES
!!  The value of ierr should always be checked by the caller
!!
!! INPUTS
!!  fh = MPI-IO file handler.
!!  offset =
!!  sizes(2)
!!  subsizes(2)
!!  starts(2)
!!  bufsz = dimension of buffer (takes into accout both real and imaginary part)
!!  chunk_bsize =
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!  comm = MPI communicator
!!
!! OUTPUTS
!!  buffer(bufsz)
!!  ierr=status error. A non zero value indicates that chunk_bsize is smaller that the fortran record
!!    and therefore bufsz has not been read.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      mpi_file_read,mpi_file_read_all,mpi_file_set_view,mpi_type_free
!!      xmpi_max,xmpi_min,xmpio_create_fsubarray_4d
!!
!! SOURCE

subroutine mpiotk_read_fsuba_dp2D(fh,offset,sizes,subsizes,starts,bufsz,buffer,chunk_bsize,sc_mode,comm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,comm,bufsz,sc_mode
 integer,intent(out) :: ierr
 integer(XMPI_OFFSET_KIND),intent(in) :: offset,chunk_bsize
!arrays
 integer,intent(in) :: sizes(2),subsizes(2),starts(2)
 real(dp),intent(out) :: buffer(bufsz)

!Local variables ------------------------------
!scalars
 integer :: mpierr,ptr,ncount,myfh
 integer :: fsub_type,my_ncalls,icall,ncalls,subs_x
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 !character(len=500) :: msg
!arrays
 integer :: call_subsizes(2),call_starts(2)
 integer,allocatable :: my_basead(:),my_subsizes(:,:),my_starts(:,:)
 real(dp),allocatable :: dummy_buf(:,:)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 if (bufsz < 2 * PRODUCT(subsizes) ) then
   MSG_ERROR("bufsz is too small")
 end if

 if (sc_mode==xmpio_single) then
   ! This makes the automatic tests fail (cut3d)
   !MSG_WARNING("comm != xmpi_comm_self")
 else if (sc_mode==xmpio_collective) then
   continue
 else
   MSG_ERROR("Wrong sc_mode")
 end if

 subs_x  = subsizes(1)

 call setup_fsuba_dp2D(sizes,subsizes,starts,chunk_bsize,my_basead,my_subsizes,my_starts,my_ncalls,ncalls,comm,ierr)
 if (ierr/=0) RETURN

 do icall=1,ncalls

   if (icall <= my_ncalls) then
     call_subsizes = my_subsizes(:,icall)
     call_starts   = my_starts(:,icall)
     ptr           = my_basead(icall)
   else
     ! Fake values needed to call read_all collectively.
     call_subsizes = (/subs_x, 1/)
     call_starts   = starts
   end if
   ncount = PRODUCT(call_subsizes)
   !write(std_out,*)"  icall,ptr, ncount, ",icall,ptr,ncount
   !write(std_out,*)"  call_starts",call_starts
   !write(std_out,*)"  call_subsizes",call_subsizes

   ! Create subarry file view.
   call xmpio_create_fsubarray_2D(sizes,call_subsizes,call_starts,MPI_DOUBLE_COMPLEX,fsub_type,my_offpad,mpierr)
   ABI_CHECK_MPI(mpierr,"fsubarray_2D")

   ! Update the offset.
   my_offset = offset + my_offpad

   call MPI_FILE_SET_VIEW(myfh, my_offset, MPI_BYTE, fsub_type, 'native', MPI_INFO_NULL, mpierr)
   ABI_CHECK_MPI(mpierr,"SET_VIEW")

   call MPI_TYPE_FREE(fsub_type, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

   if (sc_mode==xmpio_collective) then
     ! Collective read
     if (icall <= my_ncalls) then
       call MPI_FILE_READ_ALL(myfh, buffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     else
       ABI_MALLOC(dummy_buf,(2,subs_x))
       call MPI_FILE_READ_ALL(myfh, dummy_buf, ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
       ABI_FREE(dummy_buf)
     end if
     ABI_CHECK_MPI(mpierr,"FILE_READ_ALL")

   else
     ! Individual read.
     call MPI_FILE_READ(myfh, buffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     ABI_CHECK_MPI(mpierr,"FILE_READ")
   end if

 end do

 ABI_FREE(my_subsizes)
 ABI_FREE(my_starts)
 ABI_FREE(my_basead)

end subroutine mpiotk_read_fsuba_dp2D
!!***

!----------------------------------------------------------------------

!!****f* m_mpiotk/mpiotk_write_fsuba_dp2D
!! NAME
!!  mpiotk_write_fsuba_dp2D
!!
!! FUNCTION
!!  Write a block of contiguous data stored in a 2D matrix.
!!  Data is placed within Fortran records. Target: complex data stored in a real array.
!!
!! NOTES
!!  The value of ierr should always be checked by the caller
!!
!! INPUTS
!!  fh = MPI-IO file handler.
!!  offset =
!!  sizes(2)
!!  subsizes(2)
!!  starts(2)
!!  bufsz = dimension of buffer (takes into accout both real and imaginary part)
!!  chunk_bsize =
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!  comm = MPI communicator
!!
!! OUTPUTS
!!  buffer(bufsz)
!!  ierr=status error. A non zero value indicates that chunk_bsize is smaller that the fortran record
!!    and therefore bufsz has not been written.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      mpi_file_read,mpi_file_read_all,mpi_file_set_view,mpi_type_free
!!      xmpi_max,xmpi_min,xmpio_create_fsubarray_4d
!!
!! SOURCE

subroutine mpiotk_write_fsuba_dp2D(fh,offset,sizes,subsizes,starts,bufsz,buffer,chunk_bsize,sc_mode,comm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,comm,bufsz,sc_mode
 integer,intent(out) :: ierr
 integer(XMPI_OFFSET_KIND),intent(in) :: offset,chunk_bsize
!arrays
 integer,intent(in) :: sizes(2),subsizes(2),starts(2)
 real(dp),intent(in) :: buffer(bufsz)

!Local variables ------------------------------
!scalars
 integer :: mpierr,ptr,ncount,myfh
 integer :: fsub_type,my_ncalls,icall,ncalls,subs_x
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 !character(len=500) :: msg
!arrays
 integer :: call_subsizes(2),call_starts(2)
 integer,allocatable :: my_basead(:),my_subsizes(:,:),my_starts(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ! Workaround for XLF
 myfh = fh

 if (bufsz < 2 * PRODUCT(subsizes) ) then
   MSG_ERROR("bufsz is too small")
 end if

 if (sc_mode==xmpio_single) then
   call setup_fsuba_dp2D(sizes,subsizes,starts,chunk_bsize,my_basead,my_subsizes,my_starts,my_ncalls,ncalls,xmpi_comm_self,ierr)
 else if (sc_mode==xmpio_collective) then
   call setup_fsuba_dp2D(sizes,subsizes,starts,chunk_bsize,my_basead,my_subsizes,my_starts,my_ncalls,ncalls,comm,ierr)
 else
   MSG_ERROR("Wrong sc_mode")
 end if
 if (ierr/=0) RETURN

 subs_x = subsizes(1)
 do icall=1,ncalls

   if (icall <= my_ncalls) then
     call_subsizes = my_subsizes(:,icall)
     call_starts   = my_starts(:,icall)
     ptr           = my_basead(icall)
   else
     ! Fake values needed to call write_all collectively.
     call_subsizes = (/subs_x, 1/)
     call_starts   = starts
   end if
   ncount = PRODUCT(call_subsizes)
   !write(std_out,*)"  icall,ptr, ncount, ",icall,ptr,ncount
   !write(std_out,*)"  call_starts",call_starts
   !write(std_out,*)"  call_subsizes",call_subsizes

   ! Create subarry file view.
   call xmpio_create_fsubarray_2D(sizes,call_subsizes,call_starts,MPI_DOUBLE_COMPLEX,fsub_type,my_offpad,mpierr)
   ABI_CHECK_MPI(mpierr,"fsubarray_2D")

   ! Update the offset.
   my_offset = offset + my_offpad

   call MPI_FILE_SET_VIEW(myfh, my_offset, MPI_BYTE, fsub_type, 'native', MPI_INFO_NULL, mpierr)
   ABI_CHECK_MPI(mpierr,"SET_VIEW")

   call MPI_TYPE_FREE(fsub_type, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

   if (sc_mode==xmpio_collective) then
     ! Collective write
     if (icall <= my_ncalls) then
       call MPI_FILE_WRITE_ALL(myfh, buffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     else
       ! Re-write my first chunk of data.
       call MPI_FILE_WRITE_ALL(myfh, buffer(1), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     end if
     ABI_CHECK_MPI(mpierr,"FILE_WRITE_ALL")

   else
     ! Individual write.
     call MPI_FILE_WRITE(myfh, buffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     ABI_CHECK_MPI(mpierr,"FILE_WRITE")
   end if

 end do

 ABI_FREE(my_subsizes)
 ABI_FREE(my_starts)
 ABI_FREE(my_basead)

 DBG_EXIT("COLL")

end subroutine mpiotk_write_fsuba_dp2D
!!***

!----------------------------------------------------------------------

!!****f* m_mpiotk/mpiotk_read_fsuba_dpc3D
!! NAME
!!  mpiotk_read_fsuba_dpc3D
!!
!! FUNCTION
!!  Read of a block of contiguous data stored in a 3D matrix.
!!  Data is placed within Fortran records. Target: complex data stored in a complex array.
!!
!! NOTES
!!  The value of ierr should always be checked by the caller
!!
!! INPUTS
!!  fh = MPI-IO file handler.
!!  offset =
!!  sizes(3)
!!  subsizes(3)
!!  starts(3)
!!  bufsz = dimension of cbuffer
!!  chunk_bsize =
!!  comm = MPI communicator
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  cbuffer(bufsz)
!!  ierr=status error. A non zero value indicates that chunk_bsize is smaller that the fortran record
!!    and therefore bufsz has not been read.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      mpi_file_read,mpi_file_read_all,mpi_file_set_view,mpi_type_free
!!      xmpi_max,xmpi_min,xmpio_create_fsubarray_4d
!!
!! SOURCE

subroutine mpiotk_read_fsuba_dpc3D(fh,offset,sizes,subsizes,starts,bufsz,cbuffer,chunk_bsize,sc_mode,comm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,comm,bufsz,sc_mode
 integer,intent(out) :: ierr
 integer(XMPI_OFFSET_KIND),intent(in) :: offset,chunk_bsize
!arrays
 integer,intent(in) :: sizes(3),subsizes(3),starts(3)
 complex(dpc),intent(out) :: cbuffer(bufsz)

!Local variables-------------------------------
!scalars
 integer :: mpierr,nz2read,nz_chunk,ptr,ncount,myfh
 integer :: fsub_type,my_ncalls,icall,zrest,ncalls
 integer :: size_x,size_y,subs_x,subs_y,subs_xy,start_x,start_y,start_z,stop_z
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 !character(len=500) :: msg
!arrays
 integer :: call_subsizes(3),call_starts(3)
 integer,allocatable :: my_basead(:),my_subsizes(:,:),my_starts(:,:)
 complex(dpc),allocatable :: dummy_cbuf(:)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 if (bufsz < PRODUCT(subsizes) ) then
   MSG_ERROR("bufsz is too small")
 end if

 if (sc_mode==xmpio_single) then
  ABI_CHECK(comm==xmpi_comm_self,"comm != xmpi_comm_self")
 else if (sc_mode==xmpio_collective) then
   continue
 else
   MSG_ERROR("Wrong sc_mode")
 end if

 size_x  = sizes(1)
 size_y  = sizes(2)
 subs_x  = subsizes(1)
 subs_y  = subsizes(2)
 subs_xy = subs_x * subs_y
 start_x = starts(1)
 start_y = starts(2)
 start_z = starts(3)
 stop_z  = start_z + subsizes(3)-1 ! last column to read

 ! Read rows in blocks of size nz_chunk:
 ! MPI-IO crashes if we try to read data > 2Gb in a single call.
 nz2read = subsizes(3)
 nz_chunk = nz2read
 if ( (subs_xy*nz2read*xmpi_bsize_dpc) > chunk_bsize) then
   nz_chunk = chunk_bsize / (subs_xy*xmpi_bsize_dpc)
   !if (nz_chunk == 0) nz_chunk = 50
 end if

 call xmpi_min(nz_chunk,ierr,comm,mpierr)
 if (ierr == 0) then
   ierr = 1
   RETURN
 end if
 ierr = 0

 ! my_ncalls : number of read needed to fill my cbuffer.
 ! ncalls    : max number of read in comm (needed for collective operations).
 my_ncalls = nz2read / nz_chunk
 zrest = MOD(nz2read, nz_chunk)
 if (zrest /= 0) my_ncalls = my_ncalls + 1

 call xmpi_max(my_ncalls,ncalls,comm,mpierr)

 ! Compute arrays used to define the file view.
 ABI_MALLOC(my_subsizes,(3,ncalls))
 ABI_MALLOC(my_starts,(3,ncalls))
 ABI_MALLOC(my_basead,(ncalls))

 do icall=1,my_ncalls
   if (icall*nz_chunk <= nz2read) then
     my_subsizes(:,icall) = (/subs_x, subs_y, nz_chunk/)
     my_starts(:,icall) = (/start_x, start_y, (icall-1) * nz_chunk + start_z/)
     my_basead(icall) = 1 + (icall-1)*subs_xy*nz_chunk
   else
     ! Two cases:
     ! 1) nz2read > nz_chunk and not divisible by nz2read
     ! 2) nz2read < nz_chunk
     my_subsizes(:,icall) = (/subs_x, subs_y, zrest/)
     if (nz2read >= nz_chunk) then
       my_starts(:,icall) = (/start_x, start_y, stop_z-zrest+1/)
       my_basead(icall) = 1 + (nz2read-zrest) * subs_xy
     else
       my_starts(:,icall) = starts
       my_basead(icall) = 1
     end if
   end if
 end do
 !write(std_out,*)" >>>> my_ncalls, ncalls, nz2read, nz_chunk ",my_ncalls,ncalls,nz2read,nz_chunk

 do icall=1,ncalls

   if (icall <= my_ncalls) then
     call_subsizes = my_subsizes(:,icall)
     call_starts   = my_starts(:,icall)
     ptr           = my_basead(icall)
   else
     ! Fake values needed to call read_all collectively.
     call_subsizes = (/subs_x, 1, 1/)
     call_starts   = starts
   end if
   ncount = PRODUCT(call_subsizes)
   !write(std_out,*)"  icall,ptr, ncount, ",icall,ptr,ncount
   !write(std_out,*)"  call_starts",call_starts
   !write(std_out,*)"  call_subsizes",call_subsizes

   ! Create subarry file view.
   call xmpio_create_fsubarray_3D(sizes,call_subsizes,call_starts,MPI_DOUBLE_COMPLEX,&
&    fsub_type,my_offpad,mpierr)
   ABI_CHECK_MPI(mpierr,"fsubarray_3D")

   ! Update the offset.
   my_offset = offset + my_offpad

   call MPI_FILE_SET_VIEW(myfh, my_offset, MPI_BYTE, fsub_type, 'native', MPI_INFO_NULL, mpierr)
   ABI_CHECK_MPI(mpierr,"SET_VIEW")

   call MPI_TYPE_FREE(fsub_type, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

   if (sc_mode==xmpio_collective) then
     ! Collective read
     if (icall <= my_ncalls) then
       call MPI_FILE_READ_ALL(myfh, cbuffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     else
       ABI_MALLOC(dummy_cbuf,(subs_x))
       call MPI_FILE_READ_ALL(myfh, dummy_cbuf, ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
       ABI_FREE(dummy_cbuf)
     end if
     ABI_CHECK_MPI(mpierr,"FILE_READ_ALL")

   else
     ! Individual read.
     call MPI_FILE_READ(myfh, cbuffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     ABI_CHECK_MPI(mpierr,"FILE_READ")
   end if

 end do

 ABI_FREE(my_subsizes)
 ABI_FREE(my_starts)
 ABI_FREE(my_basead)

end subroutine mpiotk_read_fsuba_dpc3D
!!***

!----------------------------------------------------------------------

!!****f* m_mpiotk/mpiotk_read_fsuba_dpc4D
!! NAME
!!  mpiotk_read_fsuba_dpc4D
!!
!! FUNCTION
!!  Reading a block of contiguous data stored in a 4D matrix.
!!  Data is placed within Fortran records. Target: complex data stored in a complex array.
!!
!! NOTES
!!  The value of ierr should always be checked by the caller
!!
!! INPUTS
!!  fh = MPI-IO file handler.
!!  offset =
!!  sizes(4)
!!  subsizes(4)
!!  starts(4)
!!  bufsz = dimension of cbuffer
!!  chunk_bsize =
!!  comm = MPI communicator
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  cbuffer(bufsz)
!!  ierr=status error. A non zero value indicates that chunk_bsize is smaller that the fortran record
!!    and therefore bufsz has not been read.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      mpi_file_read,mpi_file_read_all,mpi_file_set_view,mpi_type_free
!!      xmpi_max,xmpi_min,xmpio_create_fsubarray_4d
!!
!! SOURCE

subroutine mpiotk_read_fsuba_dpc4D(fh,offset,sizes,subsizes,starts,bufsz,cbuffer,chunk_bsize,sc_mode,comm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,comm,bufsz,sc_mode
 integer,intent(out) :: ierr
 integer(XMPI_OFFSET_KIND),intent(in) :: offset,chunk_bsize
!arrays
 integer,intent(in) :: sizes(4),subsizes(4),starts(4)
 complex(dpc),intent(out) :: cbuffer(bufsz)

!Local variables-------------------------------
!scalars
 integer :: mpierr,na2read,na_chunk,ptr,ncount,myfh
 integer :: fsub_type,my_ncalls,icall,arest,ncalls
 integer :: size_x,size_y,size_z,subs_x,subs_y,subs_z,subs_xyz,start_x,start_y,start_z,start_a,stop_a
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 !character(len=500) :: msg
!arrays
 integer :: call_subsizes(4),call_starts(4)
 integer,allocatable :: my_basead(:),my_subsizes(:,:),my_starts(:,:)
 complex(dpc),allocatable :: dummy_cbuf(:)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 if (bufsz < PRODUCT(subsizes) ) then
   MSG_ERROR("bufsz is too small")
 end if

 if (sc_mode==xmpio_single) then
  ABI_CHECK(comm==xmpi_comm_self,"comm != xmpi_comm_self")
 else if (sc_mode==xmpio_collective) then
   continue
 else
   MSG_ERROR("Wrong sc_mode")
 end if

 size_x  = sizes(1)
 size_y  = sizes(2)
 size_z  = sizes(3)

 subs_x  = subsizes(1)
 subs_y  = subsizes(2)
 subs_z  = subsizes(3)
 subs_xyz = subs_x * subs_y * subs_z

 start_x = starts(1)
 start_y = starts(2)
 start_z = starts(3)
 start_a = starts(4)
 stop_a  = start_a + subsizes(4)-1 ! last column to read

 ! Read rows in blocks of size na_chunk:
 ! MPI-IO crashes if we try to read data > 2Gb in a single call.
 na2read = subsizes(4)
 na_chunk = na2read
 if ( (subs_xyz*na2read*xmpi_bsize_dpc) > chunk_bsize) then
   na_chunk = chunk_bsize / (subs_xyz*xmpi_bsize_dpc)
 end if

 call xmpi_min(na_chunk,ierr,comm,mpierr)
 if (ierr == 0) then
   ierr = 1
   RETURN
 end if
 ierr = 0

 ! my_ncalls : number of read needed to fill my cbuffer.
 ! ncalls    : max number of read in comm (needed for collective operations).
 my_ncalls = na2read / na_chunk
 arest = MOD(na2read, na_chunk)
 if (arest /= 0) my_ncalls = my_ncalls + 1

 call xmpi_max(my_ncalls,ncalls,comm,mpierr)

 ! Compute arrays used to define the file view.
 ABI_MALLOC(my_subsizes,(4,ncalls))
 ABI_MALLOC(my_starts,(4,ncalls))
 ABI_MALLOC(my_basead,(ncalls))

 do icall=1,my_ncalls

   if (icall*na_chunk <= na2read) then
     my_subsizes(:,icall) = (/subs_x, subs_y, subs_z, na_chunk/)
     my_starts(:,icall) = (/start_x, start_y, start_z, (icall-1) * na_chunk + start_a/)
     my_basead(icall) = 1 + (icall-1)*subs_xyz*na_chunk
   else
     ! Two cases:
     ! 1) na2read > na_chunk and not divisible by na2read
     ! 2) na2read < na_chunk
     my_subsizes(:,icall) = (/subs_x, subs_y, subs_z, arest/)
     if (na2read >= na_chunk) then
       my_starts(:,icall) = (/start_x, start_y, start_z, stop_a-arest+1/)
       my_basead(icall) = 1 + (na2read-arest) * subs_xyz
     else
       my_starts(:,icall) = starts
       my_basead(icall) = 1
     end if
   end if
 end do
 !write(std_out,*)" >>>> my_ncalls, ncalls, na2read, na_chunk ",my_ncalls,ncalls,na2read,na_chunk

 do icall=1,ncalls

   if (icall <= my_ncalls) then
     call_subsizes = my_subsizes(:,icall)
     call_starts   = my_starts(:,icall)
     ptr           = my_basead(icall)
   else
     ! Fake values needed to call read_all collectively.
     call_subsizes = (/subs_x, 1, 1, 1/)
     call_starts   = starts
   end if
   ncount = PRODUCT(call_subsizes)
   !
   !write(std_out,*)"  icall,ptr, ncount, ",icall,ptr,ncount
   !write(std_out,*)"  call_starts",call_starts
   !write(std_out,*)"  call_subsizes",call_subsizes
   !
   ! Create subarry file view.
   call xmpio_create_fsubarray_4D(sizes,call_subsizes,call_starts,MPI_DOUBLE_COMPLEX,&
&    fsub_type,my_offpad,mpierr)
   ABI_CHECK_MPI(mpierr,"fsubarray_4D")

   ! Update the offset.
   my_offset = offset + my_offpad

   call MPI_FILE_SET_VIEW(myfh, my_offset, MPI_BYTE, fsub_type, 'native', MPI_INFO_NULL, mpierr)
   ABI_CHECK_MPI(mpierr,"SET_VIEW")

   call MPI_TYPE_FREE(fsub_type, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

   if (sc_mode==xmpio_collective) then
     ! Collective read
     if (icall <= my_ncalls) then
       call MPI_FILE_READ_ALL(myfh, cbuffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     else
       ABI_MALLOC(dummy_cbuf,(subs_x))
       call MPI_FILE_READ_ALL(myfh, dummy_cbuf, ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
       ABI_FREE(dummy_cbuf)
     end if
     ABI_CHECK_MPI(mpierr,"FILE_READ_ALL")
   else
     ! Individual read.
     call MPI_FILE_READ(myfh, cbuffer(ptr), ncount, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     ABI_CHECK_MPI(mpierr,"FILE_READ")
   end if

 end do

 ABI_FREE(my_subsizes)
 ABI_FREE(my_starts)
 ABI_FREE(my_basead)

end subroutine mpiotk_read_fsuba_dpc4D
!!***
#endif

END MODULE m_mpiotk
!!***
