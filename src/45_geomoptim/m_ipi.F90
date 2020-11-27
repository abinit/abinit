!!****m* ABINIT/m_ipi
!! NAME
!!  m_ipi
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2020-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!!  From the i-pi user guide available at: https://github.com/i-pi/i-pi/blob/master/doc/manual.tex
!!
!!  The server assumes that 4-byte integers, 8-byte floats and 1-byte characters are used.
!!  The typical communication flow is as follows:
!!
!!  1. a header string "STATUS" is sent by the server to the client that has connected to it;
!!
!!  2. a header string is then returned, giving the status of the client code.
!!     Recognized messages are:
!!
!!  "NEEDINIT": if the client code needs any initialising data, it can be sent here.
!!  The server code will then send a header string “INIT”, followed by an integer
!!  corresponding to the bead index, another integer giving the number of bits
!!  in the initialization string, and finally the initialization string itself.
!!
!!  "READY": sent if the client code is ready to calculate the forces. The server
!!  socket will then send a string "POSDATA", then nine floats for the cell vector
!!  matrix, then another nine floats for the inverse matrix. The server socket
!!  will then send one integer giving the number of atoms, then the position data
!!  as 3 floats for each atom giving the 3 cartesian components of its position.
!!
!!  "HAVEDATA": is sent if the client has finished computing the potential and
!!  forces. The server socket then sends a string "GETFORCE", and the client
!!  socket returns "FORCEREADY". The potential is then returned as a float,
!!  the number of atoms as an integer, then the force data as 3 floats per atom
!!  in the same way as the positions, and the virial as 9 floats in the same
!!  way as the cell vector matrix. Finally, the client may return an arbitrary
!!  string containing additional data that have been obtained by the electronic
!!  structure calculation (atomic charges, dipole moment). The client first
!!  returns an integer specifying the number of characters, and then the string,
!!  which will be output verbatim if this "extra" information is requested in the
!!  output section (see 3.2.2).
!!
!!  3. The server socket waits until the force data for each replica of the system has
!!  been calculated and returned, then the MD can be propagated for one more time
!!  step, and new force requests will be dispatched.
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

module m_ipi

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist
 use m_xmpi
 use m_errors

 use m_fstrings,  only : sjoin, itoa
 use m_fsockets !, only : socket_from_string
 use m_geometry,  only : xcart2xred, det3r !fcart2fred, xred2xcart

 implicit none

 private

 public :: ipi_setup
 public :: ipi_shutdown
 public :: ipi_pred
 public :: ipi_check_initial_consistency
!!***

! =============
! private stuff
! =============

 INTEGER, PARAMETER :: HDRLEN=12
 integer, save :: socket !, socket_initialized = 0
 integer, save :: origin_natom = 0
 real(dp), save :: origin_rprimd(3, 3) = zero
 real(dp), save, allocatable :: origin_xred(:,:)

contains
!!***

!!****f* m_ipi/ipi_setup
!! NAME
!! ipi_setup
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! OUTPUT
!!
!! SOURCE

subroutine ipi_setup(string, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string
 integer,intent(in) :: comm

!local variables
 integer,parameter :: master = 0
 integer :: my_rank, ierr
 character(len=HDRLEN) :: header

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, "Initializing i-pi protocol")
 if (my_rank == master) then
   call socket_from_string(string, socket)
   call readbuffer(socket, header, HDRLEN)
   ABI_CHECK(trim(header) == "STATUS", sjoin("Expecting STATUS header, got:", trim(header)))
   call writebuffer(socket, "READY       ", HDRLEN)
   call readbuffer(socket, header, HDRLEN)
   ABI_CHECK(trim(header) == "POSDATA", sjoin("Expecting POSDATA header, got:", trim(header)))
   call handle_posdata(origin_natom, origin_rprimd, origin_xred)
 end if

 call xmpi_bcast(origin_natom, master, comm, ierr)
 call xmpi_bcast(origin_rprimd, master, comm, ierr)
 if (my_rank /= master) then
   ABI_MALLOC(origin_xred, (3, origin_natom))
 end if
 call xmpi_bcast(origin_xred, master, comm, ierr)
 call wrtout(std_out, "ipi_setup completed")

end subroutine ipi_setup
!!***

subroutine ipi_check_initial_consistency(in_natom, in_rprimd, in_xred, ierr)

 integer,intent(in) :: in_natom
 real(dp),intent(in) :: in_rprimd(3,3)
 real(dp),intent(in) :: in_xred(3,in_natom)
 integer,intent(out) :: ierr

 integer :: ii

 call wrtout(std_out, "ipi mode: Checking whether initial geometry from server agrees with input file")

 ierr = 0
 if (in_natom /= origin_natom) then
   MSG_WARNING(sjoin("in_natom:", sjoin(itoa(in_natom), " != origin_natom", itoa(origin_natom))))
   ierr = ierr + 1
 end if

 if (any(abs(in_rprimd - origin_rprimd) > tol6)) then
   MSG_WARNING("Mismatch between input file and data from socket: in_rprimd and origin_rprimd do not agree within 1e-6")
   write(std_out, "(a)")" in_rprind(:,ii), origin_rprimd(:,ii)"
   do ii=1,3
     write(std_out, *)in_rprimd(:,ii), origin_rprimd(:,ii)
   end do
   ierr = ierr + 1
 end if

 if (.not. allocated(origin_xred)) then
   ierr = ierr + 1
   MSG_WARNING("origin_xred is not allocated!")
 end if

 if (in_natom == origin_natom .and. allocated(origin_xred)) then
   if (any(abs(in_xred - origin_xred) > tol6)) then
     MSG_WARNING("Mismatch between input file and data from socket: in_xred and origin_xred do not agree withing 1e-6")
     ierr = ierr + 1
     write(std_out, "(a)")" in_xred(:,ii), origin_xred(:,ii)"
     do ii=1,3
       write(std_out, *)in_xred(:,ii), origin_xred(:,ii)
     end do
   end if
 end if

 !ABI_CHECK(ierr == 0, "Initial structure sent by i-pi server does not match the one read from file! See messages above")

end subroutine ipi_check_initial_consistency
!!***

subroutine handle_posdata(out_natom, out_rprimd, out_xred)

 integer,intent(out) :: out_natom
 real(dp),intent(out) :: out_rprimd(3,3)
 real(dp),allocatable,intent(out) :: out_xred(:,:)

 REAL*8 :: cellh(3,3), gprimd(3,3), mtxbuffer(9)
 REAL*8, ALLOCATABLE :: combuf(:), xcart(:,:)

 call readbuffer(socket, mtxbuffer, 9)
 out_rprimd = transpose(reshape(mtxbuffer, [3, 3]))
 call readbuffer(socket, mtxbuffer, 9)
 gprimd = transpose(reshape(mtxbuffer, [3, 3]))
 call readbuffer(socket, out_natom)
 ABI_MALLOC(combuf, (3 * out_natom))
 call readbuffer(socket, combuf, 3 * out_natom)
 ABI_MALLOC(xcart, (3, out_natom))
 xcart = reshape(combuf, [3, out_natom])
 ABI_MALLOC(out_xred, (3, out_natom))
 call xcart2xred(out_natom, out_rprimd, xcart, out_xred)
 ABI_FREE(combuf)
 ABI_FREE(xcart)

end subroutine handle_posdata
!!***

!!****f* m_ipi/ipi_shutdown
!! NAME
!! ipi_shutdown
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! OUTPUT
!!
!! SOURCE

subroutine ipi_shutdown()

! *************************************************************************

 origin_natom = 0
 origin_rprimd = zero
 ABI_SFREE(origin_xred)
 !call close_socket()

end subroutine ipi_shutdown
!!***

!!****f* m_ipi/ipi_pred
!! NAME
!! ipi_pred
!!
!! FUNCTION
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the predictor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces,acell, rprimd, stresses
!!
!! OUTPUT
!!
!! SOURCE

subroutine ipi_pred(ab_mover, hist, itime, ntime, zDEBUG, iexit, comm_cell)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime, ntime, iexit, comm_cell
 logical,intent(in) :: zDEBUG
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(inout) :: hist

!local variables
 integer,parameter :: master = 0
 integer :: my_rank, natom, new_natom, ierr, ihist_prev
 real(dp) :: etotal, ucvol
 character(len=HDRLEN) :: header
!arrays
 real(dp) :: acell(3), rprimd(3,3), xred(3,ab_mover%natom), strten(6), sigma(3,3)
 real(dp) :: new_rprimd(3, 3)
 real(dp),allocatable :: new_xred(:,:)

! *************************************************************************

 !print *, "in ipi_pred"

 natom = ab_mover%natom
 ! Obtain the present values from the history
 call hist2var(acell, hist, natom, rprimd, xred, zDEBUG)
 strten(:) = hist%strten(:,hist%ihist)
 etotal = hist%etot(hist%ihist)
 !call fcart2fred(hist%fcart(:,:,hist%ihist), residual, rprimd, natom)

 !call metric(gmet, gprimd, iout, rmet, rprimd, ucvol)
 ucvol = det3r(rprimd)

 if (my_rank == master) then
   call readbuffer(socket, header, HDRLEN)
   ABI_CHECK(trim(header) == "STATUS", sjoin("Expecting STATUS header, got:", trim(header)))
   call writebuffer( socket, "HAVEDATA    ", HDRLEN )
   call readbuffer(socket, header, HDRLEN)
   ABI_CHECK(trim(header) == "GETFORCE", sjoin("Expecting GETFORCE header, got:", trim(header)))

   ! Communicate energy info back to i-pi
   !WRITE(*,*) " @ DRIVER MODE: Returning etotal, forces, stress tensor "
   call writebuffer(socket, "FORCEREADY  ", HDRLEN)
   call writebuffer(socket, etotal)
   call writebuffer(socket, natom)
   call writebuffer_dv(socket, hist%fcart(:,:,hist%ihist), 3 * natom)
   sigma(1,1) = strten(1)
   sigma(2,2) = strten(2)
   sigma(3,3) = strten(3)
   sigma(3,2) = strten(4)
   sigma(2,3) = strten(4)
   sigma(3,1) = strten(5)
   sigma(1,3) = strten(5)
   sigma(2,1) = strten(6)
   sigma(1,2) = strten(6)
   !call stress_voigt_to_mat(strten, sigma)
   !sigma = sigma * ucvol
   ! No need to tranpose as it's symmetric
   ! TODO: Check volume
   call writebuffer_dv(socket, sigma, 9)
   !call writebuffer(socket, combuf, 3 * nat)
   !call writebuffer(socket, reshape( vir, (/9/) ), 9)
   !
   ! Note: i-pi can also receive an arbitrary string, that will be printed
   ! out to the "extra" trajectory file. This is useful if you want to
   ! return additional information, e.g. atomic charges, wannier centres,
   ! etc. one must return the number of characters, then the string. Here
   ! we just send back zero characters.
   !
   call writebuffer(socket, 0)

   ! ================================================
   ! Here master waits for new lattice and positions
   ! ================================================
   call readbuffer(socket, header, HDRLEN)
   ABI_CHECK(trim(header) == "STATUS", sjoin("Expecting STATUS header, got:", trim(header)))
   call writebuffer(socket, "READY       ", HDRLEN)
   call readbuffer(socket, header, HDRLEN)
   ABI_CHECK(trim(header) == "POSDATA", sjoin("Expecting POSDATA header, got:", trim(header)))
   call handle_posdata(new_natom, new_rprimd, new_xred)

   ! Some basic consistency checks. Obviously there are lot of things that can go wrong
   ! if the driver changes the space group or the input file is not consistent with
   ! the logic implemented by the server.
   ABI_CHECK(new_natom == ab_mover%natom, "ipi server shall not change the number of atoms!")
   if (ab_mover%optcell == 0 .and. any(abs(new_rprimd - origin_rprimd) > tol6)) then
     MSG_ERROR("Mismatch between origin_rprimd and data from socket: origin_rprimd and new_rprimd do not agree within 1e-6")
   end if
 end if

 ! ====================
 ! All procs block here
 ! ====================
 call xmpi_barrier(comm_cell)
 call xmpi_bcast(new_rprimd, master, comm_cell, ierr)
 call xmpi_bcast(new_xred, master, comm_cell, ierr)

 ! Update the history with the prediction
 ! Increase indexes
 hist%ihist = abihist_findIndex(hist, +1)
 !call mkradim(acell,rprim,rprimd_symm)

 ! Fill the history with the variables xred, acell, rprimd, vel
 acell = one
 call var2hist(acell, hist, ab_mover%natom, new_rprimd, new_xred, zDEBUG)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%vel(:,:,hist%ihist) = hist%vel(:,:,ihist_prev)

 !if(zDEBUG)then
 !  write (std_out,*) 'residual:'
 !  do kk=1,ab_mover%natom
 !    write (std_out,*) residual(:,kk)
 !  end do
 !  write (std_out,*) 'strten:'
 !  write (std_out,*) strten(1:3),ch10,strten(4:6)
 !  write (std_out,*) 'etotal:'
 !  write (std_out,*) etotal
 !end if

end subroutine ipi_pred
!!***

end module m_ipi
!!***
