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

    void ctqmc_triqs_run( bool rot_inv, bool leg_measure, bool orb_off_diag, bool spin_off_diag, bool move_shift, bool move_double,
                          bool measure_density_matrix, bool time_invariance, bool use_norm_as_weight, int loc_n_min, int loc_n_max,
                          int seed_a, int seed_b, int num_orbitals, int n_tau, int n_l, int n_cycles_, int cycle_length,
                          int ntherm, int ntherm2, int det_init_size, int det_n_operations_before_check, int ntau_delta,
                          int nbins_histo, int rank, int nspinor, int iatom, int ilam, double beta, double move_global_prob,
                          double imag_threshold, double det_precision_warning, double det_precision_error,
                          double det_singular_threshold, double lam, complex<double> *ftau, complex<double> *gtau,
                          complex<double> *gl, complex<double> *udens_cmplx, complex<double> *vee_cmplx, complex<double> *levels_cmplx,
                          complex<double> *moments_self_1, complex<double> *moments_self_2, double *eu, double *occ ) ;

    many_body_op_t init_Hamiltonian( h_scalar_t *eps, int nflavor, h_scalar_t *udens, bool orb_off_diag,
                                     bool spin_off_diag, double lambda, std::vector<string> &labels );

    many_body_op_t init_fullHamiltonian( h_scalar_t *eps, int nflavor, h_scalar_t *vee, bool orb_off_diag,
                                         bool spin_off_diag, double lambda, std::vector<string> &labels );

    pair<int,int> convert_indexes(int iflavor, bool off_diag, bool spin_off_diag, int ndim);

    int convert_indexes_back(int iblock, int o, bool orb_off_diag, bool spin_off_diag, int ndim);

    void build_dlr(int wdlr_size, int *ndlr, double *wdlr, double lam, double eps);
}

#endif

