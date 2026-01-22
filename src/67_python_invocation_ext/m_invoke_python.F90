!!****m* ABINIT/invocation_python_interface
!! NAME
!!  invocation_python_interface
!!
!! FUNCTION
!! Write function
!!
!! COPYRIGHT
!! Copyright (C) 2006-2026 ABINIT group (OGingras)
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

module m_invoke_python
   use defs_basis
   use m_specialmsg
   use m_errors
   use m_invocation_tools
   use m_io_tools,            only : flush_unit

   implicit none

#include "mpif.h"

   public :: invoke_python_run_script
   contains

   subroutine invoke_python_run_script(dft_iter, rank, filapp_in, comm)
     use iso_c_binding, only: c_char, c_null_char
   ! subroutine Invoke_python_triqs (rank, filapp_in) bind(c)
     integer, intent(in)            :: dft_iter
     integer, intent(in)            :: rank
     character(len=500), intent(in) :: filapp_in
     integer, intent(in)            :: comm

     character(len=500) :: msg, triqs_filename, triqs_pythpath, dft_iter_filename
     integer :: ierr, mrank, msize
     logical :: exists
     character(kind=c_char, len=255) :: f2c_string

     write(msg, '(a)') "  --- Using python invocation ---"
     call wrtout(std_out, msg, 'COLL')

     write(msg, '(a, i3)') "   invoke_python_triqs: rank: ", rank
     call wrtout(std_out, msg, 'COLL')

     write(msg, '(a, a)') "   invoke_python_triqs: filapp_in: ", trim(filapp_in)
     call wrtout(std_out, msg, 'COLL')

     call mpi_comm_size(comm, msize, ierr)
     call mpi_comm_rank(comm, mrank, ierr)

     write(triqs_filename, '(2a)') trim(filapp_in), '_PY_INVOCATION_script.py'
     write(triqs_pythpath, '(3a)') './', trim(filapp_in), '_PY_INVOCATION_python_lib'

     if(mrank == 0) then
        write(msg, '(2a)') "   invoke_python_triqs: python script: ", trim(triqs_filename)
        call wrtout(std_out, msg, 'COLL')
        write(msg, '(2a)') "   invoke_python_triqs: python lib: ", trim(triqs_pythpath)
        call wrtout(std_out, msg, 'COLL')

        write(msg, '(a)') "   invoke_python_triqs: checking python library"
        call wrtout(std_out, msg, 'COLL')
     endif

     inquire( file=triqs_pythpath, exist=exists )
     if (.not.exists) then
        write(msg, '(2a)') '   invoke_python_triqs: ERROR cannot find library at ', trim(triqs_pythpath)
        ABI_ERROR(msg)
     endif

     write(f2c_string, '(a)') trim(triqs_pythpath)//c_null_char
     ierr = init_python_interpreter(f2c_string)
     write(msg, '(a, i3)') "   ierr from initialization: ", ierr
     call wrtout(std_out, msg, 'COLL')
     if (ierr == 1) then
        write(msg, '(a)') '   invoke_python_triqs: ERROR could not initialize properly the python interpreter.'
        ABI_ERROR(msg)
     endif

     write(msg, '(a)') "   invoke_python_triqs: interpreter initialized"
     call wrtout(std_out, msg, 'COLL')

     if (mrank == 0) then
        write(msg, '(a)') "   invoke_python_triqs: reading python script"
        call wrtout(std_out, msg, 'COLL')
     endif

     inquire( file=triqs_filename, exist=exists )
     if (.not. exists) then
        write(msg, '(2a)') '   invoke_python_triqs: ERROR cannot find script at ', trim(triqs_filename)
        ABI_ERROR(msg)
     endif

     write(dft_iter_filename, '(2a)') trim(filapp_in), '.dft_iter'
     open(unit=102, file=dft_iter_filename)
       write(102, "(i4)") dft_iter
     close(102)

     call mpi_barrier(comm, ierr)
     write(f2c_string, '(a)') trim(triqs_filename)//c_null_char

     write(msg, '(a)') "   ######################################"
     call wrtout(std_out, msg, 'COLL')
     write(msg, '(a)') "   ### EXECUTION OF THE PYTHON SCRIPT ###"
     call wrtout(std_out, msg, 'COLL')
     write(msg, '(a)') "   ######################################"
     call wrtout(std_out, msg, 'COLL')
     call flush_unit(std_out)

     ierr = execute_python_file(f2c_string)

     write(msg, '(a)') "   ######################################"
     call wrtout(std_out, msg, 'COLL')
     write(msg, '(a)') "   ### END OF THE EXECUTION OF PYTHON ###"
     call wrtout(std_out, msg, 'COLL')
     write(msg, '(a)') "   ######################################"
     call wrtout(std_out, msg, 'COLL')

     call mpi_barrier(comm, ierr)

     if (ierr == 1) then
        write(msg, '(a)') '   invoke_python_triqs: ERROR could not execute the python file correctly.'
        ABI_ERROR(msg)
     endif

     if (mrank == 0) then
        write(msg, '(a)') "   invoke_python_triqs: script runned"
        call wrtout(std_out, msg, 'COLL')
     endif



     ! ierr = close_python_interpreter()
     ! call mpi_finalize(ierr)
     ! if (ierr == 1) then
     !    write(msg, '(a, i3)') "   invoke_python_triqs: MPI is finalized on rank ", mrank
     !    call wrtout(std_out, msg, 'COLL')
     ! endif

   call flush_unit(std_out)

   end subroutine invoke_python_run_script

end module m_invoke_python

