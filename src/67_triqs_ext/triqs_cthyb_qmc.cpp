#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <bitset>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <triqs_cthyb/solver_core.hpp>
#include <triqs_cthyb_qmc.hpp>

using namespace h5;
using namespace mpi;
using namespace triqs::gfs;

void ctqmc_triqs_run(bool rot_inv, bool leg_measure, bool off_diag, bool move_shift, bool move_double,
                     bool measure_density_matrix, bool time_invariance, bool use_norm_as_weight,
                     int loc_n_min, int loc_n_max, int seed_a, int seed_b, int num_orbitals, int n_tau,
                     int n_l, int n_cycles, int cycle_length, int ntherm, int ntherm_restart, int det_init_size,
                     int det_n_operations_before_check, int ntau_delta, int nbins_histo, int rank,
                     int nspinor, int iatom, int ilam, int nblocks, double beta, double move_global_prob, double imag_threshold,
                     double det_precision_warning, double det_precision_error, double det_singular_threshold,
                     double lam, int *block_list, int *flavor_list, int *inner_list, int *siz_list, complex<double> *ftau,
                     complex<double> *gtau, complex<double> *gl, complex<double> *udens_cmplx, complex<double> *vee_cmplx,
                     complex<double> *levels_cmplx, complex<double> *moments_self_1, complex<double> *moments_self_2,
                     complex<double> *occ, complex<double> *eu) {

  int verbo = 1;
  int ndim = num_orbitals / 2;
  string lam_fname = "";
  if (ilam != 0) lam_fname = "_ilam_" + to_string(ilam);
  string qmc_data_fname = "qmc_data_iatom_" + to_string(iatom) + lam_fname + ".h5";
  auto comm = MPI_COMM_WORLD;
  int nproc;
  MPI_Comm_size(comm,&nproc);
  int itask = 0;
  int iflavor,iflavor1;

  int therm = ntherm;
  if (exists(qmc_data_fname)) therm = ntherm_restart;

  // Print information about the Solver Parameters
  if (rank == 0 && verbo > 0) {

    cout << endl << "   == Key Input Parameters for the TRIQS CTHYB solver ==" << endl << endl;
    cout << "   Beta                  = " << fixed << beta << endl;
    cout << "   Nflavor               = " << num_orbitals << endl;
    cout << "   Ntau                  = " << n_tau << endl;
    cout << "   Nl                    = " << n_l << endl;
    cout << "   Legendre measurement  = " << leg_measure << endl;
    cout << "   Ncycles               = " << n_cycles << endl;
    cout << "   Cycle length          = " << cycle_length << endl;
    cout << "   Ntherm                = " << therm << endl;
    cout << "   Seed                  = " << seed_a << " + " << seed_b << " * rank" << endl;
    cout << "   Move shift            = " << move_shift << endl;
    cout << "   Move double           = " << move_double << endl;
    cout << "   Meas. density mat.    = " << measure_density_matrix << endl;
    cout << "   Use norm as weight    = " << use_norm_as_weight << endl;
    cout << "   N min                 = " << loc_n_min << endl;
    cout << "   N max                 = " << loc_n_max << endl;
    cout << "   Det init size         = " << det_init_size << endl;
    cout << "   Det N ops. bef. check = " << det_n_operations_before_check << endl;
    cout << "   Global moves prob.    = " << scientific << move_global_prob << endl;
    cout << "   Imaginary Threshold   = " << imag_threshold << endl;
    cout << "   Det Precision warning = " << det_precision_warning << endl;
    cout << "   Det Precision error   = " << det_precision_error << endl;
    cout << "   Det sing. threshold   = " << det_singular_threshold << endl;
    cout << "   Time invariance       = " << time_invariance << endl;
    cout << "   Ntau delta            = " << ntau_delta << endl;
    cout << "   Nbins histo           = " << nbins_histo << endl;
    cout << "   Lambda                = " << lam << endl << endl;
  }

  // Hamiltonian definition
  many_body_op_t H;
  many_body_operator Hint;
  h_scalar_t *levels,*udens,*vee;

#ifdef HAVE_TRIQS_COMPLEX
  levels = levels_cmplx;
  udens  = udens_cmplx;
  vee    = vee_cmplx;
#else
  double levels_real [num_orbitals*num_orbitals] = {};
  double udens_real [num_orbitals*num_orbitals] = {};
  double vee_real [num_orbitals*num_orbitals*num_orbitals*num_orbitals] = {};
  for (int i = 0; i < num_orbitals*num_orbitals; ++i) {
    levels_real[i] = levels_cmplx[i].real();
    udens_real[i]  = udens_cmplx[i].real();
  }
  for (int i = 0; i < num_orbitals*num_orbitals*num_orbitals*num_orbitals; ++i)
    vee_real[i] = vee_cmplx[i].real();
  levels = levels_real;
  udens  = udens_real;
  vee    = vee_real;
#endif

  h_scalar_t levels_zero [num_orbitals*num_orbitals] = {0};

  gf_struct_t gf_struct;

  for (int i = 0; i < nblocks; ++i)
    gf_struct.emplace_back(to_string(i),siz_list[i]);

  if (rank == 0 && verbo > 0) cout << "   == Green's Function Structure Initialized ==" << endl << endl;

  // Init Hamiltonian Basic terms
  if (!rot_inv) {
    if (rank == 0 && verbo > 0) cout << "   == Density-Density Terms Included == " << endl << endl;
    H    = init_Hamiltonian(levels,num_orbitals,udens,block_list,inner_list,lam);
    Hint = init_Hamiltonian(levels_zero,num_orbitals,udens,block_list,inner_list,double(1));
  }
  else {
    if (rank == 0 && verbo > 0) cout << "   == Rotationally Invariant Terms Included ==  " << endl << endl;
    H    = init_fullHamiltonian(levels,num_orbitals,vee,block_list,inner_list,lam);
    Hint = init_fullHamiltonian(levels_zero,num_orbitals,vee,block_list,inner_list,double(1));
  }

  // Construct CTQMC solver with mesh parameters
  solver_core solver({beta,gf_struct,(n_tau-1)/2,n_tau,n_l,ntau_delta,true});
  if (rank == 0 && verbo > 0) cout << "   == Solver Core Initialized ==" << endl << endl;

  for (int iblock = 0; iblock < nblocks; ++iblock) {
    for (int tau = 0; tau < ntau_delta; ++tau) {
      for (int o = 0; o < siz_list[iblock]; ++o) {
        iflavor = flavor_list[o+iblock*num_orbitals];
        for (int oo = 0; oo < siz_list[iblock]; ++oo) {
          iflavor1 = flavor_list[oo+iblock*num_orbitals];
          solver.Delta_tau()[iblock][tau](o,oo) = ftau[iflavor+iflavor1*num_orbitals+tau*num_orbitals*num_orbitals];
        }
      }
    }
  }

  if (rank == 0 && verbo > 0) cout << "   == F(tau) Initialized ==      " << endl << endl;

  // Solver parameters
  if (rank == 0 && verbo > 0) cout << "   == Solver Parametrization ==	" << endl << endl;
  auto paramCTQMC = solve_parameters_t(H,n_cycles);

  many_body_op_t hloc0;
  paramCTQMC.h_loc0 = hloc0;
  paramCTQMC.max_time = -1;
  paramCTQMC.measure_weight_ratio = true;
  paramCTQMC.random_name = "";
  paramCTQMC.random_seed = seed_b * rank + seed_a;
  paramCTQMC.length_cycle = cycle_length;
  paramCTQMC.n_warmup_cycles = therm;
  paramCTQMC.time_invariance = time_invariance;
  paramCTQMC.nbins_histo = nbins_histo;

  if (exists(qmc_data_fname)) {

    // Read initial configuration
    if (rank == 0 && verbo > 0) cout << "   == Reading previous configuration == " << endl;
    std::vector<int> sendcounts(nproc);
    file qmc_data_hfile(qmc_data_fname,'r');
    group grp = qmc_data_hfile;
    group gr  = grp.open_group("config");
    if (rank == 0) {
      int ncpus;
      h5_read(gr,"ncpus",ncpus);
      if (ncpus != nproc) TRIQS_RUNTIME_ERROR << "You are trying to read a configuration file generated with a different number of CPUs";
      h5_read(gr,"size",sendcounts);
    }
    MPI_Bcast(&sendcounts[0],nproc,MPI_INTEGER,0,comm);
    int length = sendcounts[rank];
    std::vector<uint64_t> tau_list(length);
    std::vector<int> block_index_list(length);
    std::vector<int> inner_index_list(length);
    std::vector<int> is_dagger_list(length);
    int displs [nproc] = {0};
    for (int i = 1; i < nproc; ++i) displs[i] = displs[i-1] + sendcounts[i-1];
    int tot_size = displs[nproc-1] + sendcounts[nproc-1];
    std::vector<uint64_t> tau_list_tot(tot_size);
    std::vector<int> block_index_list_tot(tot_size);
    std::vector<int> inner_index_list_tot(tot_size);
    std::vector<int> is_dagger_list_tot(tot_size);
    if (rank == 0) {
      h5_read(gr,"tau",tau_list_tot);
      h5_read(gr,"block",block_index_list_tot);
      h5_read(gr,"inner",inner_index_list_tot);
      h5_read(gr,"dagger",is_dagger_list_tot);
    }
    MPI_Scatterv(&tau_list_tot[0],&sendcounts[0],displs,MPI_UNSIGNED_LONG_LONG,&tau_list[0],length,MPI_UNSIGNED_LONG_LONG,0,comm);
    MPI_Scatterv(&block_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&block_index_list[0],length,MPI_INTEGER,0,comm);
    MPI_Scatterv(&inner_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&inner_index_list[0],length,MPI_INTEGER,0,comm);
    MPI_Scatterv(&is_dagger_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&is_dagger_list[0],length,MPI_INTEGER,0,comm);
    for (int i = 0; i < length; ++i) {
      bool is_dagger = (is_dagger_list[i] == 1 ? true : false);
      auto op = op_desc{block_index_list[i],inner_index_list[i],is_dagger,0}; // linear index will be set later, within TRIQS
      time_pt tau = time_pt(tau_list[i],beta);
      paramCTQMC.initial_configuration.insert(tau,op);
    }

    // Read binned time histogram
    if (rank == 0 && verbo > 0) cout << "   == Reading previous histogram == " << endl;
    gr = grp.open_group("histo");
    int nbins_file;
    h5_read(gr,"nbins",nbins_file);
    if (nbins_file != nbins_histo) TRIQS_RUNTIME_ERROR << "You are trying to read a configuration file generated with a different number of time bins";
    std::vector<double> *hist;
    string tag_move;
    for (int i = 0; i < nblocks; ++i) {
      for (int j = 0; j < 2; ++j) {
        if (j == 0) {
          tag_move = "insert";
          hist = &paramCTQMC.hist_insert[to_string(i)];
        }
        else {
          tag_move = "remove";
          hist = &paramCTQMC.hist_remove[to_string(i)];
        }
        *hist = std::vector<double>(nbins_histo);
        if (rank == 0) h5_read(gr,tag_move+"_"+to_string(i),*hist);
        mpi_broadcast(*hist,comm,0);
      }
    }
    qmc_data_hfile.close();
  }

  paramCTQMC.move_shift = move_shift;
  paramCTQMC.move_double = move_double;
  paramCTQMC.measure_density_matrix = measure_density_matrix;
  paramCTQMC.use_norm_as_weight = use_norm_as_weight;
  paramCTQMC.det_init_size = det_init_size;
  paramCTQMC.det_n_operations_before_check = det_n_operations_before_check;
  paramCTQMC.move_global_prob = move_global_prob;
  paramCTQMC.imag_threshold = imag_threshold;
  paramCTQMC.det_precision_warning = det_precision_warning;
  paramCTQMC.det_precision_error = det_precision_error;
  paramCTQMC.det_singular_threshold = det_singular_threshold;
  paramCTQMC.loc_n_min = loc_n_min;
  paramCTQMC.loc_n_max = loc_n_max;
  paramCTQMC.measure_G_l = leg_measure;

  if (rank == 0 && verbo > 0) cout << "   == Starting Solver [node " << rank << "] ==" << endl << endl;
  solver.solve(paramCTQMC);

  if (rank == 0 && verbo > 0) cout << endl << "   == Reporting ==" << endl << endl;

  if (rank == 0 && verbo > 0) cout << "   == Writing final configuration ==      " << endl << endl;

  // Write all final configurations
  auto config = solver.get_configuration();
  auto tau_0 = time_pt(1,beta);
  int length = config.size();
  std::vector<uint64_t> tau_list(length);
  std::vector<int> block_index_list(length);
  std::vector<int> inner_index_list(length);
  std::vector<int> is_dagger_list(length);
  int count = 0;
  for (auto const &o : config) {
    uint64_t tau = floor_div(o.first,tau_0);   // trick to get the value of the integer representing the time_pt
    auto op = o.second;
    int block_index = op.block_index;
    int inner_index = op.inner_index;
    int is_dagger = (op.dagger ? 1 : 0);
    tau_list[count] = tau;
    block_index_list[count] = block_index;
    inner_index_list[count] = inner_index;
    is_dagger_list[count] = is_dagger;
    ++count;
  }
  std::vector<int> recvcounts(nproc);
  MPI_Allgather(&length,1,MPI_INTEGER,&recvcounts[0],1,MPI_INTEGER,comm);
  int displs [nproc] = {0};
  for (int i = 1; i < nproc; ++i) displs[i] = displs[i-1] + recvcounts[i-1];
  int tot_size = displs[nproc-1] + recvcounts[nproc-1];
  std::vector<uint64_t> tau_list_tot(tot_size);
  std::vector<int> block_index_list_tot(tot_size);
  std::vector<int> inner_index_list_tot(tot_size);
  std::vector<int> is_dagger_list_tot(tot_size);
  MPI_Gatherv(&tau_list[0],length,MPI_UNSIGNED_LONG_LONG,&tau_list_tot[0],&recvcounts[0],displs,MPI_UNSIGNED_LONG_LONG,0,comm);
  MPI_Gatherv(&block_index_list[0],length,MPI_INTEGER,&block_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
  MPI_Gatherv(&inner_index_list[0],length,MPI_INTEGER,&inner_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
  MPI_Gatherv(&is_dagger_list[0],length,MPI_INTEGER,&is_dagger_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);

  if (rank == 0) {

    file qmc_data_hfile(qmc_data_fname,'w');   // very important to open in write mode, so that previous data is erased
    group grp = qmc_data_hfile;
    group gr  = grp.create_group("config");

    h5_write(gr,"ncpus",nproc);
    h5_write(gr,"size",recvcounts);
    h5_write(gr,"tau",tau_list_tot);
    h5_write(gr,"block",block_index_list_tot);
    h5_write(gr,"inner",inner_index_list_tot);
    h5_write(gr,"dagger",is_dagger_list_tot);

    // Write all histograms

    gr = grp.create_group("histo");
    h5_write(gr,"nbins",nbins_histo);

    auto hist_insert = solver.get_weight_ratio_insert();
    auto hist_remove = solver.get_weight_ratio_remove();
    double step = beta / nbins_histo;
    std::vector<double> *hist;
    double p_local = 0.9;
    string tag_move;

    for (int i = 0; i < nblocks; ++i) {
      for (int j = 0; j < 2; j++) {

        std::vector<int> ind_vec; // indices of unsampled bins
        ind_vec.reserve(nbins_histo);

        if (j == 0) {
          hist = &hist_insert[to_string(i)];
          tag_move = "insert";
        }
        else {
          tag_move = "remove";
          hist = &hist_remove[to_string(i)];
        }

        double s = 0.;
        for (int k = 0; k < nbins_histo; ++k) {
          if ((*hist)[k] < 1.e-15)
            ind_vec.push_back(k);
          else
            s += (*hist)[k];
        }
        s *= step;

        double x = s * (1. / p_local - 1.) / (step * double(ind_vec.size()));
        for (auto const &ind : ind_vec) (*hist)[ind] = x;

        h5_write(gr,tag_move+"_"+to_string(i),*hist);
      }
    }

    qmc_data_hfile.close();
  }

  for (int iblock = 0; iblock < nblocks; ++iblock) {
    for (int tau = 0; tau < n_tau; ++tau) {
      for (int o = 0; o < siz_list[iblock]; ++o) {
        iflavor = flavor_list[o+iblock*num_orbitals];
        for (int oo = 0; oo < siz_list[iblock]; ++oo) {
          iflavor1 = flavor_list[oo+iblock*num_orbitals];
          gtau[tau+iflavor*n_tau+iflavor1*num_orbitals*n_tau] = (*solver.G_tau)[iblock].data()(tau,o,oo);
        }
      }
    }
  }

  // Report G(l)
  if (leg_measure) {
    for (int iblock = 0; iblock < nblocks; ++iblock) {
      for (int l = 0; l < n_l; ++l) {
        for (int o = 0; o < siz_list[iblock]; ++o) {
          iflavor = flavor_list[o+iblock*num_orbitals];
          for (int oo = 0; oo < siz_list[iblock]; ++oo) {
            iflavor1 = flavor_list[oo+iblock*num_orbitals];
            gl[l+iflavor*n_l+iflavor1*num_orbitals*n_l] = (*solver.G_l)[iblock].data()(l,o,oo);
          }
        }
      }
    }
  }

  if (measure_density_matrix) {

    auto h_loc_diag = solver.h_loc_diagonalization();
    many_body_operator n_op;
    auto subspaces = h_loc_diag.n_subspaces();
    auto rho = solver.density_matrix();

    complex<double> occ_tmp [num_orbitals] = {0};

    for (int iblock = 0; iblock < nblocks; ++iblock) {
      for (int o = 0; o < siz_list[iblock]; ++o) {
        iflavor = flavor_list[o+iblock*num_orbitals];
        n_op = c_dag(to_string(iblock),o) * c(to_string(iblock),o);
        if (itask == rank) occ_tmp[iflavor] = trace_rho_op(rho,n_op,h_loc_diag);
        itask = (itask+1)%nproc;
      }
    }

    MPI_Allreduce(occ_tmp,occ,num_orbitals,MPI_C_DOUBLE_COMPLEX,MPI_SUM,comm);

    if (rank == 0) {
      ofstream occ_file;
      occ_file.open("Occupation_numbers_triqs");

      for (int sub : range(subspaces)) {
        auto fock_states = h_loc_diag.get_fock_states(sub);
        int sub_dim = h_loc_diag.get_subspace_dim(sub);
        auto unit_mat = h_loc_diag.get_unitary_matrix(sub);
        auto rho_tmp = unit_mat * rho[sub] * dagger(unit_mat);  // transform density matrix to Fock basis
        for (int ind : range(sub_dim)) {
          auto f_state = fock_states[ind];
          auto prob = rho_tmp(ind,ind);
          occ_file << bitset<64>(f_state).to_string().substr(64-num_orbitals) << "\t" << scientific << setprecision(17) << prob << endl;
        }
      }
    }

    if (itask == rank) *eu = trace_rho_op(rho,Hint,h_loc_diag);
    mpi_broadcast(*eu,comm,itask);
    itask = (itask+1)%nproc;

    if (!leg_measure) { // Get moments of the self-energy

      many_body_operator commut,commut2,Sinf_op,S1_op;

      complex<double> mself2 [num_orbitals] = {0};
      auto Sinf_mat = matrix_t(num_orbitals,num_orbitals);
      Sinf_mat = h_scalar_t{0};

      for (int iblock = 0; iblock < nblocks; ++iblock) {
        for (int o = 0; o < siz_list[iblock]; ++o) {
          iflavor = flavor_list[o+iblock*num_orbitals];
          commut = Hint*c(to_string(iblock),o) - c(to_string(iblock),o)*Hint;
          Sinf_op = - commut*c_dag(to_string(iblock),o) - c_dag(to_string(iblock),o)*commut;
          commut2 = c_dag(to_string(iblock),o)*Hint - Hint*c_dag(to_string(iblock),o);
          S1_op = commut2*commut + commut*commut2;
          if (itask == rank) Sinf_mat(iflavor,iflavor) = trace_rho_op(rho,Sinf_op,h_loc_diag);
          itask = (itask+1)%nproc;
          if (itask == rank) mself2[iflavor] = trace_rho_op(rho,S1_op,h_loc_diag);
          itask = (itask+1)%nproc;
        }
      }

      Sinf_mat = all_reduce(Sinf_mat,comm);
      MPI_Allreduce(mself2,moments_self_2,num_orbitals,MPI_C_DOUBLE_COMPLEX,MPI_SUM,comm);

      for (int iflavor = 0; iflavor < num_orbitals; ++iflavor)
        moments_self_1[iflavor] = Sinf_mat(iflavor,iflavor);

      Sinf_mat = Sinf_mat * Sinf_mat;

      for (int iflavor = 0; iflavor < num_orbitals; ++iflavor)
        moments_self_2[iflavor] -= Sinf_mat(iflavor,iflavor);

    }   // not legendre

  }   // measure_density_matrix
  if (rank == 0 && verbo > 0) cout << "   == CTHYB-QMC Process Finished [node "<< rank <<"] =="<< endl << endl;
}

/********************************************************/
/****************** Functions Used **********************/
/********************************************************/

// Build density-density Hamiltonian
many_body_op_t init_Hamiltonian(h_scalar_t *eps, int nflavor, h_scalar_t *udens,
                                int *block_list, int *inner_list, double lambda) {

  many_body_op_t H;
  int iblock,iblock1,o,oo;

  for (int i = 0; i < nflavor; ++i) {
    iblock = block_list[i]; o = inner_list[i];
    for (int j = 0; j < nflavor; ++j) {
      iblock1 = block_list[j]; oo = inner_list[j];
      H += 0.5 * lambda * udens[i+j*nflavor] * n(to_string(iblock),o) * n(to_string(iblock1),oo);
      H += eps[i+j*nflavor] * c_dag(to_string(iblock),o) * c(to_string(iblock1),oo);
    }
  }
  return H;
}

// Build full Hamiltonian
many_body_op_t init_fullHamiltonian(h_scalar_t *eps, int nflavor, h_scalar_t *vee,
                                    int *block_list, int *inner_list, double lambda) {

  many_body_op_t H;
  int iblock,iblock1,iblock2,iblock3,o,oo,ooo,oooo;

  for (int i = 0; i < nflavor; ++i) {
    iblock = block_list[i]; o = inner_list[i];
    for (int j = 0; j < nflavor; ++j) {
      iblock1 = block_list[j]; oo = inner_list[j];
      H += eps[i+j*nflavor] * c_dag(to_string(iblock),o) * c(to_string(iblock1),oo);
      for (int k = 0; k < nflavor; ++k) {
        iblock2 = block_list[k]; ooo = inner_list[k];
        for (int l = 0; l < nflavor; ++l) {
          iblock3 = block_list[l]; oooo = inner_list[l];
          H += 0.5 * lambda * vee[i+j*nflavor+k*nflavor*nflavor+l*nflavor*nflavor*nflavor] * c_dag(to_string(iblock),o) *
               c_dag(to_string(iblock1),oo) * c(to_string(iblock3),oooo) * c(to_string(iblock2),ooo);
        }
      }
    }
  }
  return H;
}

// Build DLR frequencies
void build_dlr(int wdlr_size, int *ndlr, double *wdlr, double lam, double eps) {
  auto omega = cppdlr::build_dlr_rf(lam,eps);
  *ndlr = omega.size();
  if (*ndlr > wdlr_size) return;
  for (int i = 0; i < (*ndlr); ++i) wdlr[i] = omega[i];
}

