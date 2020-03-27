
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

MODULE TRIQS_CTQMC

  implicit none

  interface

! =====================================================================
!              Subroutine Ctqmc_triqs_run() linked to TRIQS
! =====================================================================
!
! Ctqmc_triqs_run(Boolean args, int args, real args, pointer array args)
! 
! =====================================================================

 
    subroutine Ctqmc_triqs_run (    rot_inv, leg_measure, hist, wrt_files, tot_not,                                    &

                                    nflavor, nfreq, ntau, nl, ncycle, cycle_length, ntherm, verbosity_solver,seed,beta,&

                               &    levels,  u_mat_ij, u_mat_ijkl, fiw_nd, g_iw, gtau, gl, comm                        ) bind( c )
     
      use iso_c_binding

      LOGICAL(Kind=1), VALUE, INTENT(IN) :: rot_inv, leg_measure, hist, wrt_files, tot_not
 
      INTEGER,  VALUE, INTENT(IN)        :: nflavor, nfreq, ntau, nl, ncycle, cycle_length, ntherm, verbosity_solver, seed 

      REAL(Kind=8), VALUE, INTENT(IN)    :: beta

      TYPE(C_PTR), VALUE, INTENT(IN)     :: levels,  u_mat_ij, u_mat_ijkl, fiw_nd, g_iw, gtau, gl

      INTEGER,  VALUE, INTENT(IN)        :: comm
 
    end subroutine Ctqmc_triqs_run
   
  end interface


END MODULE TRIQS_CTQMC
