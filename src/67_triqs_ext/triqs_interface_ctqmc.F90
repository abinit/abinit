
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


     SUBROUTINE Ctqmc_triqs_run(rot_inv,leg_measure,off_diag,move_shift,move_double,measure_density_matrix,time_invariance, &
                              & use_norm_as_weight,loc_n_min,loc_n_max,seed_a,seed_b,nflavor,ntau,nl,ncycle,cycle_length, &
                              & ntherm,ntherm2,det_init_size,det_n_operations_before_check,ntau_delta,nbins_histo, &
                              & rank,nspinor,iatom,ilam,beta,move_global_prob,imag_threshold,det_precision_warning, &
                              & det_precision_error,det_singular_threshold,lam,ftau,gtau,gl,udens,vee,levels,moments_self_1, &
                              & moments_self_2,Eu,occ) bind(c)

      use iso_c_binding

      LOGICAL(Kind=1), VALUE, INTENT(IN) :: rot_inv,leg_measure,off_diag,move_shift,move_double

      LOGICAL(Kind=1), VALUE, INTENT(IN) :: measure_density_matrix,time_invariance,use_norm_as_weight

      INTEGER, VALUE, INTENT(IN) :: loc_n_min,loc_n_max,seed_a,seed_b,nflavor,ntau,nl,ncycle,cycle_length,ntherm,ntherm2

      INTEGER, VALUE, INTENT(IN) :: det_init_size,det_n_operations_before_check,ntau_delta,nbins_histo,rank,nspinor,iatom,ilam

      REAL(Kind=8), VALUE, INTENT(IN) :: beta,move_global_prob,imag_threshold,det_precision_warning,det_precision_error

      REAL(Kind=8), VALUE, INTENT(IN) :: det_singular_threshold,lam

      TYPE(C_PTR), VALUE, INTENT(IN) :: ftau,gtau,gl,udens,vee,levels,moments_self_1,moments_self_2,Eu,occ

    end subroutine Ctqmc_triqs_run

    subroutine build_dlr(wdlr_size,ndlr,wdlr,lam,eps) bind(c)

      use iso_c_binding

      INTEGER, VALUE, INTENT(IN) :: wdlr_size

      TYPE(C_PTR), VALUE, INTENT(IN)  :: ndlr,wdlr

      REAL(Kind=8), VALUE, INTENT(IN) :: lam,eps

    end subroutine build_dlr

  end interface


END MODULE TRIQS_CTQMC
