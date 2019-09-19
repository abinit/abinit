#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __TRIQS_CTHYB_QMC_H__ 
#define __TRIQS_CTHYB_QMC_H__

#include <iostream>
#include <complex>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/utility/real_or_complex.hpp>  //ADDED for compatibility with TRIQS 1.4

using namespace std;
//using triqs::utility::many_body_operator; //COMMENTED: Not relevant with TRIQS 1.4
using triqs::operators::many_body_operator_generic; //ADDED Instead

extern "C"{
  
    void ctqmc_triqs_run( bool rot_inv, bool leg_measure, bool hist, bool wrt_files, bool tot_not,
                          int n_orbitals, int n_freq, int n_tau, int n_l, int n_cycles_, int cycle_length, int ntherm, int verbo,int seed,  
                          double beta_,
                          double *epsi, double *umat_ij, double *umat_ijkl, std::complex<double> *delta_iw_ptr, std::complex<double> *g_iw_ptr, double *g_tau, double *gl, int rank );
                          //double *epsi, double *umat_ij, double *umat_ijkl, std::complex<double> *delta_iw_ptr, std::complex<double> *g_iw_ptr, double *g_tau, double *gl, MPI_Fint *MPI_world_ptr );

//COMMENTED: Class name changed in TRIQS 1.4
//    many_body_operator<double> init_Hamiltonian( double *eps, int nflavor, double *U );
//    many_body_operator<double> init_fullHamiltonian( double *eps, int nflavor, double *U );
//    many_body_operator<double> init_fullHamiltonianUpDown( double *eps, int nflavor, double *U );

//ADDED Instead
    many_body_operator_generic<double> init_Hamiltonian( double *eps, int nflavor, double *U );
    many_body_operator_generic<double> init_fullHamiltonian( double *eps, int nflavor, double *U );
    many_body_operator_generic<double> init_fullHamiltonianUpDown( double *eps, int nflavor, double *U );
}

#endif

