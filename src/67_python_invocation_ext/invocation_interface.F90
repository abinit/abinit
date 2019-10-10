
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
! =====================================================================
!    Fortran Module Interface for C++ fonction using ISO_C_BINDING
! =====================================================================
!
! CHILDREN: Ctqmc_triqs_run() ; hellocpp()
!
! =====================================================================

MODULE INVOKE_PYTHON

  implicit none

  interface

! =====================================================================
!              Subroutine Ctqmc_triqs_run() linked to TRIQS
! =====================================================================
!
! Ctqmc_triqs_run(Boolean args, int args, real args, pointer array args)
! 
! =====================================================================

    subroutine Invoke_python_triqs (rank, filapp_in) bind(c)
      use iso_c_binding
      integer, value, intent(in) :: rank
      character(kind=c_char) :: filapp_in(*)
    end subroutine Invoke_python_triqs
 
  end interface


END MODULE INVOKE_PYTHON
