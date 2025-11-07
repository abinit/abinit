#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __TRIQS_CTHYB_QMC_H__
#define __TRIQS_CTHYB_QMC_H__

#include <complex>
#include <triqs_cthyb/config.hpp>

using namespace std;
using namespace triqs_cthyb;

extern "C"{

    void ctqmc_triqs_run( bool rot_inv, bool leg_measure, bool move_shift, bool move_double, bool measure_density_matrix, bool time_invariance,
                          bool use_norm_as_weight, bool debug, int integral, int loc_n_min, int loc_n_max, int seed_a, int seed_b,
			  int num_orbitals, int n_tau, int n_l, int n_cycles_, int cycle_length, int ntherm, int ntherm_restart, int det_init_size,
                          int det_n_operations_before_check, int rank, int nblocks, int read_data, int verbo, double beta,
                          double imag_threshold, double det_precision_warning, double det_precision_error, double det_singular_threshold, double lam_u,
                          double pauli_prob, int *block_list, int *flavor_list, int *inner_list, int *siz_list, complex<double> *ftau,
			  complex<double> *gtau, complex<double> *gl, complex<double> *udens_cmplx, complex<double> *vee_cmplx, complex<double> *levels_cmplx,
                          complex<double> *moments_self_1, complex<double> *moments_self_2, complex<double> *occ, complex<double> *eu,
 			  char *fname_data, char *fname_dataw, char *fname_histo ) ;

    many_body_op_t init_Hamiltonian( h_scalar_t *eps, int nflavor, h_scalar_t *udens, int *block_list,
                                     int *inner_list, double lambda);

    many_body_op_t init_fullHamiltonian( h_scalar_t *eps, int nflavor, h_scalar_t *vee, int *block_list,
                                         int *inner_list, double lambda);

    void build_dlr(int wdlr_size, int *ndlr, double *wdlr, double lam, double eps);
}

#endif

