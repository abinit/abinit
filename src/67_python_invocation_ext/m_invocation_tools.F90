!!****m* ABINIT/invocation_python_interface
!! NAME
!!  invocation_python_interface
!!
!! FUNCTION
!! Write function
!!
!! COPYRIGHT
!! Copyright (C) 2006-2024 ABINIT group (OGingras)
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

module m_invocation_tools

   implicit none

   interface

     subroutine test_python() bind(C, name='test_python')
     end subroutine test_python

     function init_python_interpreter(python_so) result(ierr) bind(C, name='init_python_interpreter')
        use iso_c_binding, only: c_int, c_char
        integer(c_int) :: ierr !< Return value
        character(kind=c_char), intent(in) :: python_so(*) !< Python lib path
     end function init_python_interpreter

     function execute_python_file(python_script) result(ierr) bind(C)
        use iso_c_binding, only: c_int, c_char
        integer(c_int) :: ierr !< Return value
        character(kind=c_char), intent(in) :: python_script(*) !< Python script name
     end function execute_python_file

     function close_python_interpreter() result(ierr) bind(C)
        use iso_c_binding, only: c_int, c_char
        integer(c_int) :: ierr !< Return value
     end function close_python_interpreter

   end interface

end module m_invocation_tools

