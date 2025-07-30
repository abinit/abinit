#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <bitset>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <triqs_cthyb/configuration.hpp>
#include <triqs_cthyb/solver_core.hpp>
#include <triqs_cthyb_qmc.hpp>

using namespace h5;
using namespace mpi;
using namespace triqs::gfs;

void ctqmc_triqs_run(bool rot_inv, bool leg_measure, bool move_shift, bool move_double, bool measure_density_matrix,
                     bool time_invariance, bool use_norm_as_weight, bool compute_entropy, int loc_n_min, int loc_n_max,
                     int seed_a, int seed_b, int num_orbitals, int n_tau, int n_l, int n_cycles, int cycle_length,
                     int ntherm, int ntherm_restart, int det_init_size, int det_n_operations_before_check, int rank,
                     int nblocks, int read_data, int verbo, double beta, double imag_threshold, double det_precision_warning,
                     double det_precision_error, double det_singular_threshold, double lam_u, double pauli_prob, int *block_list,
                     int *flavor_list, int *inner_list, int *siz_list, complex<double> *ftau, complex<double> *gtau,
                     complex<double> *gl, complex<double> *udens_cmplx, complex<double> *vee_cmplx, complex<double> *levels_cmplx,
                     complex<double> *moments_self_1, complex<double> *moments_self_2, complex<double> *occ, complex<double> *eu,
                     char *fname_data, char *fname_dataw, char *fname_histo) {

  int ndim = num_orbitals / 2;
  string qmc_data_fname  = string(fname_data);
  string qmc_data_fnamew = string(fname_dataw);
  string hist_fname = string(fname_histo);
  auto comm = MPI_COMM_WORLD;
  int iflavor, iflavor1, nproc;
  MPI_Comm_size(comm,&nproc);
  int itask = 0, therm = ntherm;

  // Hamiltonian definition
  many_body_op_t H;
  many_body_operator Hint;
  h_scalar_t *levels,*udens,*vee;

#ifdef HAVE_TRIQS_COMPLEX
  levels = levels_cmplx; udens = udens_cmplx; vee = vee_cmplx;
#else
  double levels_real [num_orbitals*num_orbitals] = {};
  double udens_real [num_orbitals*num_orbitals] = {};
  double vee_real [num_orbitals*num_orbitals*num_orbitals*num_orbitals] = {};

  for (int i : range(num_orbitals*num_orbitals)) {
    levels_real[i] = levels_cmplx[i].real();
    udens_real[i]  = udens_cmplx[i].real();
  }

  for (int i : range(num_orbitals*num_orbitals*num_orbitals*num_orbitals))
    vee_real[i] = vee_cmplx[i].real();

  levels = levels_real; udens = udens_real; vee = vee_real;
#endif

  h_scalar_t levels_zero [num_orbitals*num_orbitals] = {0};

  gf_struct_t gf_struct;

  for (int i : range(nblocks)) gf_struct.emplace_back(to_string(i),siz_list[i]);

  // Init Hamiltonian Basic terms
  if (!rot_inv) {
    H    = init_Hamiltonian(levels,num_orbitals,udens,block_list,inner_list,lam_u);
    Hint = init_Hamiltonian(levels_zero,num_orbitals,udens,block_list,inner_list,double(1));
  }
  else {
    H    = init_fullHamiltonian(levels,num_orbitals,vee,block_list,inner_list,lam_u);
    Hint = init_fullHamiltonian(levels_zero,num_orbitals,vee,block_list,inner_list,double(1));
  }

  // Construct CTQMC solver with mesh parameters
  solver_core solver({beta,gf_struct,(n_tau-1)/2,n_tau,n_l,true});

  // Fill Hybridization
  for (int iblock : range(nblocks))
    for (int tau : range(n_tau))
      for (int o : range(siz_list[iblock])) {
        iflavor = flavor_list[o+iblock*num_orbitals];
        for (int oo : range(siz_list[iblock])) {
          iflavor1 = flavor_list[oo+iblock*num_orbitals];
          solver.Delta_tau()[iblock][tau](o,oo) = ftau[iflavor+iflavor1*num_orbitals+tau*num_orbitals*num_orbitals];
        }
      }

  // Solver parameters
  auto paramCTQMC = solve_parameters_t(H,n_cycles);

  many_body_op_t hloc0;
  paramCTQMC.h_loc0 = hloc0;
  paramCTQMC.max_time = -1;
  paramCTQMC.random_name = "";
  paramCTQMC.random_seed = seed_a + rank * seed_b;
  paramCTQMC.length_cycle = cycle_length;
#if defined HAVE_TRIQS_v4_0
  paramCTQMC.time_invariance = time_invariance;
  paramCTQMC.pauli_prob = pauli_prob;
#endif

  int restart = (exists(qmc_data_fname) && read_data == 1 ? 1 : 0);

  file qmc_data_hfile;

#if defined HAVE_TRIQS_v4_0
  if (restart == 1 && rank == 0) {

    // Check if the configuration file can be read
    qmc_data_hfile = {qmc_data_fname,'r'};
    group grp = qmc_data_hfile;
    int err = -1;
    double beta_;
    gf_struct_t gf_struct_;
    int nbins_,nproc_;

    h5_read(grp,"beta",beta_);
    h5_read(grp,"gf_struct",gf_struct_);
    h5_read(grp,"nproc",nproc_);

    if (std::abs(beta-beta_) > 1.e-15) {
      restart = 0;
      cout << endl << "   == You are trying to read a CTQMC_DATA file generated at a different temperature !! File will not be read !" << endl;
    }

    if (gf_struct != gf_struct_) {
      restart = 0;
      cout << endl << "   == You are trying to read a CTQMC_DATA file generated with a different block structure !! File will not be read !" << endl;
    }

    if (nproc_ != nproc) {
      restart = 0;
      cout << endl << "   == You are trying to read a CTQMC_DATA file generated with a different number of CPUs !! File will not be read !" << endl;
    }

    if (restart == 0) qmc_data_hfile.close();

    if (restart == 1) cout << endl << "   == Reading CTQMC data from file " << qmc_data_fname << endl;
  }

  MPI_Bcast(&restart,1,MPI_INTEGER,0,comm);

  if (restart == 1) {

    therm = ntherm_restart;
    std::vector<int> sendcounts(nproc);

    if (rank == 0) {
      group grp = qmc_data_hfile;
      group gr = grp.open_group("config");
      h5_read(gr,"size",sendcounts);
    }

    MPI_Bcast(&sendcounts[0],nproc,MPI_INTEGER,0,comm);
    int length = sendcounts[rank];
    std::vector<uint64_t> tau_list(length);
    std::vector<int> block_index_list(length);
    std::vector<int> inner_index_list(length);
    std::vector<int> is_dagger_list(length);
    std::vector<int> lin_index_list(length);
    int displs [nproc] = {0};
    for (int i : range(1,nproc)) displs[i] = displs[i-1] + sendcounts[i-1];
    int tot_size = displs[nproc-1] + sendcounts[nproc-1];
    std::vector<uint64_t> tau_list_tot(tot_size);
    std::vector<int> block_index_list_tot(tot_size);
    std::vector<int> inner_index_list_tot(tot_size);
    std::vector<int> is_dagger_list_tot(tot_size);
    std::vector<int> lin_index_list_tot(tot_size);

    if (rank == 0) {
      group grp = qmc_data_hfile;
      group gr = grp.open_group("config");
      h5_read(gr,"tau",tau_list_tot);
      h5_read(gr,"block",block_index_list_tot);
      h5_read(gr,"inner",inner_index_list_tot);
      h5_read(gr,"dagger",is_dagger_list_tot);
      h5_read(gr,"lin",lin_index_list_tot);
    }

    MPI_Scatterv(&tau_list_tot[0],&sendcounts[0],displs,MPI_UNSIGNED_LONG_LONG,&tau_list[0],length,MPI_UNSIGNED_LONG_LONG,0,comm);
    MPI_Scatterv(&block_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&block_index_list[0],length,MPI_INTEGER,0,comm);
    MPI_Scatterv(&inner_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&inner_index_list[0],length,MPI_INTEGER,0,comm);
    MPI_Scatterv(&is_dagger_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&is_dagger_list[0],length,MPI_INTEGER,0,comm);
    MPI_Scatterv(&lin_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,&lin_index_list[0],length,MPI_INTEGER,0,comm);

    for (int i : range(length)) {
      bool is_dagger = (is_dagger_list[i] == 1 ? true : false);
      auto op = op_desc{block_index_list[i],inner_index_list[i],is_dagger,lin_index_list[i]};
      time_pt tau = time_pt(tau_list[i],beta);
      paramCTQMC.initial_configuration.insert(tau,op);
    }

    if (rank == 0) qmc_data_hfile.close();
  }
#endif

  paramCTQMC.move_shift = move_shift;
  paramCTQMC.move_double = move_double;
  paramCTQMC.measure_density_matrix = measure_density_matrix;
  paramCTQMC.use_norm_as_weight = use_norm_as_weight;
  paramCTQMC.det_init_size = det_init_size;
  paramCTQMC.det_n_operations_before_check = det_n_operations_before_check;
  paramCTQMC.imag_threshold = imag_threshold;
  paramCTQMC.det_precision_warning = det_precision_warning;
  paramCTQMC.det_precision_error = det_precision_error;
  paramCTQMC.det_singular_threshold = det_singular_threshold;
  paramCTQMC.loc_n_min = loc_n_min;
  paramCTQMC.loc_n_max = loc_n_max;
  paramCTQMC.measure_G_l = leg_measure;
  paramCTQMC.n_warmup_cycles = therm;

  if (compute_entropy) {
    paramCTQMC.measure_G_l = false;
    paramCTQMC.measure_G_tau = false;
  }

  if (rank == 0 && verbo == 1) {

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
    cout << "   Imaginary Threshold   = " << scientific << imag_threshold << endl;
    cout << "   Det Precision warning = " << det_precision_warning << endl;
    cout << "   Det Precision error   = " << det_precision_error << endl;
    cout << "   Det sing. threshold   = " << det_singular_threshold << endl;
#if defined HAVE_TRIQS_v4_0
    cout << "   Time invariance       = " << time_invariance << endl;
    cout << "   Pauli prob            = " << pauli_prob << endl;
#endif
  }

  solver.solve(paramCTQMC);

#if defined HAVE_TRIQS_v4_0
  if (rank == 0) cout << endl << "   == Writing CTQMC data on file " << qmc_data_fnamew << endl;

  // Write all final configurations
  auto config = solver.get_configuration();
  auto tau_0 = time_pt(1,beta);
  int length = config.size();
  std::vector<uint64_t> tau_list(length);
  std::vector<int> block_index_list(length);
  std::vector<int> inner_index_list(length);
  std::vector<int> is_dagger_list(length);
  std::vector<int> lin_index_list(length);
  int count = 0;
  for (auto const &o : config) {
    uint64_t tau = floor_div(o.first,tau_0); // trick to get the value of the integer representing the time_pt
    auto op = o.second;
    int block_index = op.block_index;
    int inner_index = op.inner_index;
    int is_dagger = (op.dagger ? 1 : 0);
    int lin_index = op.linear_index;
    tau_list[count] = tau;
    block_index_list[count] = block_index;
    inner_index_list[count] = inner_index;
    is_dagger_list[count] = is_dagger;
    lin_index_list[count] = lin_index;
    ++count;
  }
  std::vector<int> recvcounts(nproc);
  MPI_Allgather(&length,1,MPI_INTEGER,&recvcounts[0],1,MPI_INTEGER,comm);
  int displs [nproc] = {0};
  for (int i : range(1,nproc)) displs[i] = displs[i-1] + recvcounts[i-1];
  int tot_size = displs[nproc-1] + recvcounts[nproc-1];
  std::vector<uint64_t> tau_list_tot(tot_size);
  std::vector<int> block_index_list_tot(tot_size);
  std::vector<int> inner_index_list_tot(tot_size);
  std::vector<int> is_dagger_list_tot(tot_size);
  std::vector<int> lin_index_list_tot(tot_size);
  MPI_Gatherv(&tau_list[0],length,MPI_UNSIGNED_LONG_LONG,&tau_list_tot[0],&recvcounts[0],displs,MPI_UNSIGNED_LONG_LONG,0,comm);
  MPI_Gatherv(&block_index_list[0],length,MPI_INTEGER,&block_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
  MPI_Gatherv(&inner_index_list[0],length,MPI_INTEGER,&inner_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
  MPI_Gatherv(&is_dagger_list[0],length,MPI_INTEGER,&is_dagger_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
  MPI_Gatherv(&lin_index_list[0],length,MPI_INTEGER,&lin_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);

  if (rank == 0) {

    file qmc_data_hfilew = {qmc_data_fnamew,'w'};   // very important to open in write mode, so that previous data is erased
    group grp = qmc_data_hfilew;

    h5_write(grp,"beta",beta);
    h5_write(grp,"gf_struct",gf_struct);
    h5_write(grp,"nproc",nproc);

    group gr = grp.create_group("config");

    h5_write(gr,"size",recvcounts);
    h5_write(gr,"tau",tau_list_tot);
    h5_write(gr,"block",block_index_list_tot);
    h5_write(gr,"inner",inner_index_list_tot);
    h5_write(gr,"dagger",is_dagger_list_tot);
    h5_write(gr,"lin",lin_index_list_tot);

    qmc_data_hfilew.close();
  }
#endif

  if (!compute_entropy) {

    for (int iblock : range(nblocks))
      for (int tau : range(n_tau))
        for (int o : range(siz_list[iblock])) {
          iflavor = flavor_list[o+iblock*num_orbitals];
          for (int oo : range(siz_list[iblock])) {
            iflavor1 = flavor_list[oo+iblock*num_orbitals];
            gtau[tau+iflavor*n_tau+iflavor1*num_orbitals*n_tau] = (*solver.G_tau)[iblock].data()(tau,o,oo);
          }
        }

    // Report G(l)
    if (leg_measure)
      for (int iblock : range(nblocks))
        for (int l : range(n_l))
          for (int o : range(siz_list[iblock])) {
            iflavor = flavor_list[o+iblock*num_orbitals];
            for (int oo : range(siz_list[iblock])) {
              iflavor1 = flavor_list[oo+iblock*num_orbitals];
              gl[l+iflavor*n_l+iflavor1*num_orbitals*n_l] = (*solver.G_l)[iblock].data()(l,o,oo);
            }
          }
  }

  auto h_loc_diag = solver.h_loc_diagonalization();

  if (measure_density_matrix) {

    many_body_operator n_op;
    auto subspaces = h_loc_diag.n_subspaces();
    auto rho = solver.density_matrix();

    if (!compute_entropy) {
      complex<double> occ_tmp [num_orbitals] = {0};

      for (int iblock : range(nblocks))
        for (int o : range(siz_list[iblock])) {
          iflavor = flavor_list[o+iblock*num_orbitals];
          n_op = c_dag(to_string(iblock),o) * c(to_string(iblock),o);
          if (itask == rank) occ_tmp[iflavor] = trace_rho_op(rho,n_op,h_loc_diag);
          itask = (itask+1)%nproc;
        }

      MPI_Allreduce(occ_tmp,occ,num_orbitals,MPI_C_DOUBLE_COMPLEX,MPI_SUM,comm);
    }

    if (rank == 0 && !compute_entropy) {

      ofstream occ_file;
      occ_file.open(hist_fname);

      occ_file << "# Probability of occupying each Fock state" << endl <<
         "# Regardless of the block structure, leftmost digit always represent the first flavor" << endl <<
         "# (spin up is first), and is set to 1 (0) if the flavor is (not) occupied" << endl <<
         "#    Nb elec    Fock state           Probability" << endl;

      double prob;
      bitset<64> f_stateb;
      std::vector<std::vector<pair<double,bitset<64>>>> histogram(num_orbitals+1);
      for (int i : range(num_orbitals+1)) histogram[i].reserve(1 << i);

      for (int sub : range(subspaces)) {
        auto fock_states = h_loc_diag.get_fock_states(sub);
        int sub_dim = h_loc_diag.get_subspace_dim(sub);
        auto unit_mat = h_loc_diag.get_unitary_matrix(sub);
        auto rho_tmp = unit_mat * rho[sub] * dagger(unit_mat);  // rotate density matrix to Fock basis
        for (int ind : range(sub_dim)) {
          auto f_state = fock_states[ind];
#ifdef HAVE_TRIQS_COMPLEX
          double prob = rho_tmp(ind,ind).real();
#else
          double prob = rho_tmp(ind,ind);
#endif
          f_stateb = bitset<64>(f_state);
          bitset<64> f_stateb_;
          // Put orbitals in the correct order
          int i = 0;
          for (int iblock : range(nblocks))
            for (int o : range(siz_list[iblock])) {
              iflavor = flavor_list[o+iblock*num_orbitals];
              f_stateb_[num_orbitals-1-iflavor] = f_stateb[i];
              i++;
            }
          histogram[f_stateb_.count()].emplace_back(make_pair(prob,f_stateb_));
        }
      }
      for (int i : range(num_orbitals+1)) {
        // Sort by descending probabilities
        sort(histogram[i].begin(),histogram[i].end(), [](auto &left, auto &right) { return left.first > right.first; });
        for (int j : range(histogram[i].size())) {
          tie(prob,f_stateb) = histogram[i][j];
          occ_file << "\t" << to_string(i) << "\t" << f_stateb.to_string().substr(64-num_orbitals) << "\t" << scientific
                   << setprecision(17) << prob << endl;
        }
      }
    }

    if (itask == rank || compute_entropy) *eu = trace_rho_op(rho,Hint,h_loc_diag);
    if (!compute_entropy) {
      mpi_broadcast(*eu,comm,itask);
      itask = (itask+1)%nproc;
    }

    if (!leg_measure && !compute_entropy) { // Get moments of the self-energy

      many_body_operator commut,commut2,Sinf_op,S1_op;

      complex<double> mself2 [num_orbitals] = {0};
      auto Sinf_mat = matrix_t(num_orbitals,num_orbitals);
      Sinf_mat = h_scalar_t{0};

      for (int iblock : range(nblocks))
        for (int o : range(siz_list[iblock])) {
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

      Sinf_mat = all_reduce(Sinf_mat,comm);
      MPI_Allreduce(mself2,moments_self_2,num_orbitals,MPI_C_DOUBLE_COMPLEX,MPI_SUM,comm);

      for (int iflavor : range(num_orbitals)) moments_self_1[iflavor] = Sinf_mat(iflavor,iflavor);

      // Matrix product in case the off-diagonal elements are implemented someday
      Sinf_mat = Sinf_mat * Sinf_mat;

      for (int iflavor : range(num_orbitals)) moments_self_2[iflavor] -= Sinf_mat(iflavor,iflavor);

    }   // not legendre

  }   // measure_density_matrix
}

/********************************************************/
/****************** Functions Used **********************/
/********************************************************/

// Build density-density Hamiltonian
many_body_op_t init_Hamiltonian(h_scalar_t *eps, int nflavor, h_scalar_t *udens,
                                int *block_list, int *inner_list, double lambda) {

  many_body_op_t H;
  int iblock,iblock1,o,oo;

  for (int i : range(nflavor)) {
    iblock = block_list[i]; o = inner_list[i];
    for (int j : range(nflavor)) {
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

  for (int i : range(nflavor)) {
    iblock = block_list[i]; o = inner_list[i];
    for (int j : range(nflavor)) {
      iblock1 = block_list[j]; oo = inner_list[j];
      H += eps[i+j*nflavor] * c_dag(to_string(iblock),o) * c(to_string(iblock1),oo);
      for (int k : range(nflavor)) {
        iblock2 = block_list[k]; ooo = inner_list[k];
        for (int l : range(nflavor)) {
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
  for (int i : range((*ndlr))) wdlr[i] = omega[i];
}

