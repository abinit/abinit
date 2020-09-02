!!****m* ABINIT/m_haydock_io
!! NAME
!! m_haydock_io
!!
!! FUNCTION
!!  This module provides routines to read the Haydock file
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (YG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_haydock_io

 use defs_basis
 use m_abicore
 use m_errors

 use m_io_tools,       only : open_file

 implicit none

 private
!!***

 integer,private,parameter :: CUR_VERSION = 1

!!****t* m_haydock_io/haydock_type
!! NAME
!!  haydock_type
!!
!! FUNCTION
!!  The structure defining the content of the haydock file
!!
!! SOURCE

 type,public :: haydock_type

   integer :: version
   ! Version of the file

   integer :: hsize
   ! Size of the hamiltonian

   integer :: use_coupling
   ! 1 if coupling is used, 0 if not

   integer :: op
   ! Type of the file

   integer :: nq
   ! Number of q-points in the file

   integer :: unt
   ! Fortran unit number.

   real(dp) :: broad
   ! Broadening factor

   complex(dpc) :: factor

   real(dp),allocatable :: qpoints(:,:)
   ! qpoints(3,nq)

   integer,allocatable :: niter(:)
   ! niter(nq)

 end type haydock_type

 public :: open_haydock             ! Init the Haydock file
 public :: read_haydock             ! Reads the data associated with 1 q-point
 public :: read_dim_haydock         ! Read dimensions
 public :: write_dim_haydock        ! Write dimensions.
 public :: skip_dim_haydock         ! Skip the section with dimensions.
 public :: write_haydock            ! Writes data related to a q-point inside the file
 public :: close_haydock            ! Close the Haydock file and release memory
!!***

CONTAINS  !====================================================================
!!***

!!****f* m_haydock_io/open_haydock
!! NAME
!! open_haydock
!!
!! FUNCTION
!! Open the file and initialize haydock file descriptor that will
!! be used for later access
!!
!! INPUTS
!!  filename = Name of the file to be opened
!!
!! OUTPUT
!!  haydock_file = file descriptor for the haydock file
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine open_haydock(filename, haydock_file)

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: filename
 type(haydock_type),intent(out) :: haydock_file
!arrays

!Local variables ------------------------------
!scalars
 character(len=500) :: msg

!************************************************************************

 if (open_file(filename,msg,newunit=haydock_file%unt,form="unformatted") /= 0) then
   MSG_ERROR(msg)
 end if

 haydock_file%version = CUR_VERSION

end subroutine open_haydock
!!***

!----------------------------------------------------------------------

!!****f* m_haydock_io/read_dim_haydock
!! NAME
!! read_dim_haydock
!!
!! FUNCTION
!!  Reads the header dimensions of the haydock file and store them in
!!  the haydock file descriptor
!!
!! INPUT/OUTPUT
!!  haydock_file = haydock file descriptor
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_dim_haydock(haydock_file)

 implicit none

!Arguments ------------------------------------
!scalars
 type(haydock_type),intent(inout) :: haydock_file

!Local variables ------------------------------
!scalars
 integer :: niter_file
 integer :: iq, it
 character(len=500) :: msg
!arrays
 real(dp) :: q_file(3)

!************************************************************************

 read(haydock_file%unt) haydock_file%version

 msg = "The haydock file has been produced by an uncompatible version of abinit"
 ABI_CHECK(haydock_file%version == CUR_VERSION, msg)

 read(haydock_file%unt) haydock_file%hsize,haydock_file%use_coupling,&
&   haydock_file%op,haydock_file%nq,haydock_file%broad

 ABI_MALLOC(haydock_file%qpoints,(3,haydock_file%nq))
 ABI_MALLOC(haydock_file%niter,(haydock_file%nq))

 do iq = 1,haydock_file%nq
   read(haydock_file%unt) q_file(:)
   read(haydock_file%unt) niter_file

   haydock_file%qpoints(:,iq) = q_file(:)
   haydock_file%niter(iq) = niter_file
   ! Skip data for this q.
   do it=1,niter_file
     read(haydock_file%unt) ! it,aa(it),bb(it)
   end do
   read(haydock_file%unt) ! phi_nm1
   read(haydock_file%unt) ! phi_n
   read(haydock_file%unt) ! factor
 end do

end subroutine read_dim_haydock
!!***

!----------------------------------------------------------------------

!!****f* m_haydock_io/write_dim_haydock
!! NAME
!! write_dim_haydock
!!
!! FUNCTION
!!  Writes the basic dimensions stored inside the haydock descriptor inside the file
!!
!! INPUTS
!!  haydock_file = haydock file descriptor
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_dim_haydock(haydock_file)

 implicit none

!Arguments ------------------------------------
 type(haydock_type),intent(in) :: haydock_file

! *************************************************************************

 write(haydock_file%unt) haydock_file%version
 write(haydock_file%unt) haydock_file%hsize,haydock_file%use_coupling, &
&  haydock_file%op,haydock_file%nq,haydock_file%broad

end subroutine write_dim_haydock
!!***

!----------------------------------------------------------------------

!!****f* m_haydock_io/skip_dim_haydock
!! NAME
!! skip_dim_haydock
!!
!! FUNCTION
!!  Skip the part of the file reading basic dimensions contained in the header
!!
!! INPUTS
!!  haydock_file = haydock file descriptor
!!
!! OUTPUT
!!
!! PARENTS
!!      m_haydock_io
!!
!! CHILDREN
!!
!! SOURCE

subroutine skip_dim_haydock(haydock_file)

 implicit none

!Arguments ------------------------------------
!scalars
 type(haydock_type),intent(in) :: haydock_file

! *************************************************************************

 read(haydock_file%unt)
 read(haydock_file%unt)

end subroutine skip_dim_haydock
!!***

!----------------------------------------------------------------------

!!****f* m_haydock_io/read_haydock
!! NAME
!! read_haydock
!!
!! FUNCTION
!!  Reads the data related to one q-point inside the file accessed by file descriptor
!!
!! INPUTS
!!  haydock_file = haydock file descriptor
!!  q = q-point to be searched inside the file
!!
!! OUTPUT
!!  aa = coefficients "a" of the lanczos chain
!!  bb = coefficients "b" of the lanczos chain
!!  phi_n = last vector of the lanczos chain
!!  phi_nm1 = penultimate vector of the lanczos chain
!!  niter = number of iterations done
!!  factor = pre-factor used to obtain the green function
!!
!! NOTES
!!  niter = 0 if the q-point has not been found
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_haydock(haydock_file, q, aa, bb, phi_n, phi_nm1, niter, factor)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: niter
 complex(dpc),intent(out) :: factor
 type(haydock_type),intent(in) :: haydock_file
!arrays
 real(dp),intent(in) :: q(3)
 real(dp),allocatable,intent(out) :: bb(:)
 complex(dpc),allocatable,intent(out) :: aa(:),phi_n(:),phi_nm1(:)

!Local variables ------------------------------
!scalars
 integer :: iq, it, inn
 integer :: niter_file
 logical :: found_q
 character(len=500) :: msg
!arrays
 real(dp) :: q_file(3)

! *************************************************************************

 rewind(haydock_file%unt)
 call skip_dim_haydock(haydock_file)

 do iq = 1,haydock_file%nq
   read(haydock_file%unt) q_file(:)
   read(haydock_file%unt) niter_file

   if ( ALL(ABS(q_file - q) < tol6) ) then
      found_q = .TRUE.; EXIT
   else
      ! Skip data for this q.
      do it=1,niter_file
        read(haydock_file%unt) ! it,aa(it),bb(it)
      end do
      read(haydock_file%unt) ! phi_nm1
      read(haydock_file%unt) ! phi_n
      read(haydock_file%unt) ! factor
   end if
 end do

 if(found_q) then
   niter = niter_file
   ABI_ALLOCATE(aa,(niter))
   ABI_ALLOCATE(bb,(niter))
   do inn=1,niter
     read(haydock_file%unt)it,aa(inn),bb(inn)
     if (inn/=it) then
       write(msg,'(2(a,i0))')" Found it_file: ",it," while it should be: ",inn
       MSG_ERROR(msg)
     end if
   end do
   ABI_MALLOC(phi_nm1,(haydock_file%hsize))
   ABI_MALLOC(phi_n,(haydock_file%hsize))
   read(haydock_file%unt)phi_nm1
   read(haydock_file%unt)phi_n
   read(haydock_file%unt)factor
 else
   niter = 0
 end if

end subroutine read_haydock
!!***

!----------------------------------------------------------------------

!!****f* m_haydock_io/write_haydock
!! NAME
!! write_haydock
!!
!! FUNCTION
!!  Writes data related to a q-point inside the file
!!
!! INPUTS
!!  haydock_file = haydock file descriptor
!!  q = q-point to be searched inside the file
!!  aa = coefficients "a" of the lanczos chain
!!  bb = coefficients "b" of the lanczos chain
!!  phi_n = last vector of the lanczos chain
!!  phi_nm1 = penultimate vector of the lanczos chain
!!  niter = number of iterations done
!!  factor = pre-factor used to obtain the green function
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_haydock(haydock_file, hsize, q, aa, bb, phi_n, phi_nm1, niter, factor)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter,hsize
 complex(dpc),intent(in) :: factor
 type(haydock_type),intent(in) :: haydock_file
!arrays
 real(dp),intent(in) :: q(3)
 real(dp),intent(in) :: bb(niter)
 complex(dpc),intent(in) :: aa(niter),phi_n(hsize),phi_nm1(hsize)

!Local variables -----------------------------
!scalars
 integer :: it

! *************************************************************************

 write(haydock_file%unt) q
 write(haydock_file%unt) niter  ! NB if the previous loop completed inn=niter_max+1
 do it=1,niter        ! if we exited then inn is not incremented by one.
   write(haydock_file%unt)it,aa(it),bb(it)
 end do
 write(haydock_file%unt) phi_nm1
 write(haydock_file%unt) phi_n
 write(haydock_file%unt) factor

end subroutine write_haydock
!!***

!----------------------------------------------------------------------

!!****f* m_haydock_io/close_haydock
!! NAME
!! close_haydock
!!
!! FUNCTION
!!  Closes the haydock file and free the memory related to the file descriptor
!!
!! INPUTS
!!  haydock_file = haydock file descriptor
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine close_haydock(haydock_file)

 implicit none

!Arguments ------------------------------------
 type(haydock_type),intent(inout) :: haydock_file

! *************************************************************************

 close(haydock_file%unt)

 if(allocated(haydock_file%qpoints)) then
   ABI_FREE(haydock_file%qpoints)
 end if

 if(allocated(haydock_file%niter)) then
   ABI_FREE(haydock_file%niter)
 end if

end subroutine close_haydock
!!***

!----------------------------------------------------------------------

END MODULE m_haydock_io
!!***
