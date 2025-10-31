
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

     SUBROUTINE Ctqmc_triqs_run(rot_inv,leg_measure,move_shift,move_double,measure_density_matrix,time_invariance,use_norm_as_weight, &
                              & integral,loc_n_min,loc_n_max,seed_a,seed_b,num_orbitals,n_tau,n_l,n_cycles,cycle_length,ntherm, &
                              & ntherm_restart,det_init_size,det_n_operations_before_check,rank,nblocks,read_data,verbo,beta, &
                              & imag_threshold,det_precision_warning,det_precision_error,det_singular_threshold,lam_u,pauli_prob, &
                              & block_list,flavor_list,inner_list,siz_list,ftau,gtau,gl,udens_cmplx,vee_cmplx,levels_cmplx,moments_self_1, &
                              & moments_self_2,occ,eu,fname_data,fname_dataw,fname_histo) bind(c)

      use iso_c_binding

      LOGICAL, VALUE, INTENT(IN) :: rot_inv,leg_measure,move_shift,move_double,measure_density_matrix,time_invariance,use_norm_as_weight

      INTEGER, VALUE, INTENT(IN) :: integral,loc_n_min,loc_n_max,seed_a,seed_b,num_orbitals,n_tau,n_l,n_cycles,cycle_length,ntherm,ntherm_restart

      INTEGER, VALUE, INTENT(IN) :: det_init_size,det_n_operations_before_check,rank,nblocks,read_data,verbo

      REAL(KIND=8), VALUE, INTENT(IN) :: beta,imag_threshold,det_precision_warning,det_precision_error,det_singular_threshold,lam_u,pauli_prob

      TYPE(C_PTR), VALUE, INTENT(IN) :: block_list,flavor_list,inner_list,siz_list,ftau,gtau,gl,udens_cmplx,vee_cmplx,levels_cmplx

      TYPE(C_PTR), VALUE, INTENT(IN) :: moments_self_1,moments_self_2,occ,eu,fname_data,fname_dataw,fname_histo

    end subroutine Ctqmc_triqs_run

    subroutine build_dlr(wdlr_size,ndlr,wdlr,lam,eps) bind(c)

      use iso_c_binding

      INTEGER, VALUE, INTENT(IN) :: wdlr_size

      TYPE(C_PTR), VALUE, INTENT(IN)  :: ndlr,wdlr

      REAL(KIND=8), VALUE, INTENT(IN) :: lam,eps

    end subroutine build_dlr

  end interface


END MODULE TRIQS_CTQMC
