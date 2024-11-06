#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/*
#include <triqs_cthyb/solver_core.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/real_or_complex.hpp>
#include <triqs/gfs.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <bitset>
#include <iterator>
#include <vector>
#include <h5/h5.hpp>
#include <cppdlr/dlr_build.hpp>
#include "triqs_cthyb_qmc.hpp"
#include <mpi.h>

using namespace std;
using namespace triqs_cthyb;
using namespace triqs::gfs;
using namespace triqs::atom_diag;
using triqs::operators::many_body_operator_generic;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
*/

void ctqmc_triqs_run(bool rot_inv, bool leg_measure, bool off_diag, bool move_shift, bool move_double, 
                     bool measure_density_matrix, bool time_invariance, bool use_norm_as_weight,               
                     int loc_n_min, int loc_n_max, int seed_a, int seed_b, int num_orbitals, int n_tau,
                     int n_l, int n_cycles, int cycle_length, int ntherm, int ntherm2, int det_init_size,
                     int det_n_operations_before_check, int ntau_delta, int nbins_histo, int rank, 
                     int nspinor, int iatom, int ilam, double beta, double move_global_prob, double imag_threshold,                 
                     double det_precision_warning, double det_precision_error, double det_singular_threshold,
                     double lam, complex<double> *ftau, complex<double> *gtau,
                     complex<double> *gl, complex<double> *udens, complex<double> *vee,
                     complex<double> *levels, complex<double> *moments_self_1,
                     complex<double> *moments_self_2, double *Eu, double *occ) {
    
    cout.setf(ios::fixed);

    int verbo = 1;
    int ndim = num_orbitals / 2;
    int tndim = ndim * nspinor;
    string config_fname = "configs_iat_"+to_string(iatom)+"_ilam_"+to_string(ilam)+".h5";
    auto comm = MPI_COMM_WORLD;
    int size;
    MPI_Comm_size(comm,&size);
    
    // Print information about the Solver Parameters 
    if (rank == 0 && verbo > 0) {

      int therm_tmp = ntherm;
      if (exists("configs.h5")) therm_tmp = ntherm2;
        
      cout << endl <<"   == Key Input Parameters for the TRIQS CTHYB solver ==" << endl << endl;
      cout << setw(27) << left << "   Beta                  = " << beta << endl;
      cout << setw(27) << left << "   Nflavor               = " << num_orbitals << endl;
      cout << setw(27) << left << "   Ntau                  = " << n_tau << endl;
      cout << setw(27) << left << "   Nl                    = " << n_l << endl;
      cout << setw(27) << left << "   Legendre measurement  = " << leg_measure << endl;
      cout << setw(27) << left << "   Ncycles               = " << n_cycles << endl;
      cout << setw(27) << left << "   Cycle length          = " << cycle_length << endl;
      cout << setw(27) << left << "   Ntherm                = " << therm_tmp << endl;
      cout << setw(27) << left << "   Seed                  = " << seed_a << " + " << seed_b << " * rank" << endl;
      cout << setw(27) << left << "   Off Diag Elements     = " << off_diag << endl;
      cout << setw(27) << left << "   Move shift            = " << move_shift << endl;
      cout << setw(27) << left << "   Move double           = " << move_double << endl;
      cout << setw(27) << left << "   Meas. density mat.    = " << measure_density_matrix << endl;
      cout << setw(27) << left << "   Use norm as weight    = " << use_norm_as_weight << endl;
      cout << setw(27) << left << "   N min                 = " << loc_n_min << endl;
      cout << setw(27) << left << "   N max                 = " << loc_n_max << endl;
      cout << setw(27) << left << "   Det init size         = " << det_init_size << endl;
      cout << setw(27) << left << "   Det N ops. bef. check = " << det_n_operations_before_check << endl;
      cout << setw(27) << left << "   Global moves prob.    = " << scientific <<  move_global_prob << endl;
      cout << setw(27) << left << "   Imaginary Threshold   = " << scientific << imag_threshold << endl;
      cout << setw(27) << left << "   Det Precision warning = " << scientific << det_precision_warning << endl;
      cout << setw(27) << left << "   Det Precision error   = " << scientific << det_precision_error << endl;
      cout << setw(27) << left << "   Det sing. threshold   = " << det_singular_threshold << endl;
      cout << setw(27) << left << "   Time invariance       = " << time_invariance << endl;
      cout << setw(27) << left << "   Ntau delta            = " << ntau_delta << endl;
      cout << setw(27) << left << "   Nbins histo           = " << nbins_histo << endl;

      //for (int o = 0; o < num_orbitals; ++o)
      //  std::cout << setw(12) << left << "   e- level["<< o+1 <<"] = " << setprecision(17) << epsi[o] << endl; 
    
       // std::cout << endl;
       // std::cout << setw(13) << left  << "   U(i,j) = "  <<endl<< endl<<"\t";
    
        //for(int o = 0; o < num_orbitals; ++o){
        //    for(int oo = 0; oo < num_orbitals; ++oo)
        //    std::cout <<  setprecision(17) << fixed <<umat_ij[o+oo*num_orbitals] << "\t"; //ed != cste => vecteur a parcourir ensuite
         //   std::cout << endl<<"\t";
        //}
    }

    // Hamiltonian definition
    many_body_operator_generic<triqs::utility::real_or_complex> H;
    many_body_operator_generic<triqs::utility::real_or_complex> Hint; 

    double levels_tmp [num_orbitals*num_orbitals] = {0};
    
    // Init Hamiltonian Basic terms
    if (!rot_inv) {
      if (rank == 0 && verbo > 0) std::cout << endl << "   == Density-Density Terms Included ==	" << endl << endl; 
      H = init_Hamiltonian(levels,num_orbitals,udens,off_diag,nspinor,lam);
      Hint = init_Hamiltonian(levels_tmp,num_orbitals,udens,off_diag,nspinor,double(1));
    }
    else { // if(rot_inv){
      if (rank == 0 && verbo > 0) std::cout << "   == Rotationally Invariant Terms Included ==	" << endl << endl;
      H = init_fullHamiltonian(levels,num_orbitals,vee,off_diag,nspinor,lam);
      Hint = init_Hamiltonian(levels_tmp,num_orbitals,vee,off_diag,nspinor,double(1));
    }// else if(rank==0 && verbo>0) std::cout <<endl<<"   == Density-Density Terms Not Included ==	"<<endl<<endl;  

    int nblocks = num_orbitals;
    int siz_block = 1;
    if (off_diag && nspinor == 1) {
      nblocks = 2;
      siz_block = ndim;
    }
    if (off_diag && nspinor == 2) {
      nblocks = 1;
      siz_block = num_orbitals;
    }

    vector<string> labels(nblocks);
    vector<pair<string,long>> gf_struct;
    if (off_diag) {
      if (nspinor == 2) 
        labels[0] = "tot";
      else {
        labels[0] = "up";
        labels[1] = "down";
      }
    }
    else {
      for (int o = 0; o < num_orbitals; ++o)
        labels[o] = to_string(o);
    }
    for (int i = 0; i < nblocks; ++i)
      gf_struct.emplace_back(labels[i],siz_block);

    if (rank == 0 && verbo > 0) cout << "   == Green Function Structure Initialized ==" << endl << endl;

    // Construct CTQMC solver with mesh parameters
    solver_core solver({beta,gf_struct,(n_tau-1)/2,n_tau,n_l,true,ntau_delta});
    if (rank == 0 && verbo > 0) cout << "   == Solver Core Initialized ==" << endl << endl;
    
    for (int iblock = 0; iblock < nblocks; ++iblock) {
      for (int tau = 0; tau < ntau_delta; ++tau) {
        for (int o = 0; o < siz_block; ++o) {  
          int iflavor = iblock;
          if (off_diag) iflavor = o + iblock*ndim;
          for (int oo = 0; oo < siz_block; ++oo) {
            int iflavor1 = iblock;
            if (off_diag) iflavor1 = oo + iblock*ndim;
            solver.Delta_tau()[iblock][tau](o,oo) = ftau[iflavor+iflavor1*tndim+tau*tndim*tndim];
          }
        }
      }
    }

    if (rank == 0 && verbo > 0) cout << "   == F(tau) Initialized ==      " << endl << endl;

    // Solver parameters
    if (rank == 0 && verbo > 0) cout << "   == Solver Parametrization ==	" << endl << endl;
    auto paramCTQMC = solve_parameters_t(H,n_cycles);

    many_body_operator_generic<triqs::utility::real_or_complex> hloc0; 
    paramCTQMC.h_loc0 = hloc0;
    paramCTQMC.max_time = -1;
    paramCTQMC.random_name = "";
    paramCTQMC.random_seed = seed_b * rank + seed_a;
    paramCTQMC.length_cycle = cycle_length;
    paramCTQMC.n_warmup_cycles = ntherm;
    paramCTQMC.time_invariance = time_invariance;
    if (exists(config_fname)) {    
      if (rank == 0) cout << "   == Reading previous configuration == " << endl;
      paramCTQMC.n_warmup_cycles = ntherm2;
      vector<int> sendcounts(size);
      h5::file config_hfile(config_fname,'r');
      h5::group gr = config_hfile;
      if (rank == 0) {
        int ncpus;
        h5_read(gr,"ncpus",ncpus);
        if (ncpus != size) TRIQS_RUNTIME_ERROR << "You are trying to read a file configs.h5 generated with a different number of CPUs";
        h5_read(gr,"size",sendcounts);
      }
      MPI_Bcast(&sendcounts[0],size,MPI_INTEGER,0,comm);
      int length = sendcounts[rank];
      vector<uint64_t> tau_list(length);
      vector<int> block_index_list(length);
      vector<int> inner_index_list(length);
      vector<int> is_dagger_list(length);
      int displs [size] = {0};
      for (int i = 1; i < size; ++i) displs[i] = displs[i-1] + sendcounts[i-1];
      int tot_size = displs[size-1] + sendcounts[size-1];
      vector<uint64_t> tau_list_tot(tot_size);
      vector<int> block_index_list_tot(tot_size);
      vector<int> inner_index_list_tot(tot_size);
      vector<int> is_dagger_list_tot(tot_size);
      if (rank == 0) {
        h5_read(gr,"tau",tau_list_tot);
        h5_read(gr,"block",block_index_list_tot);
        h5_read(gr,"inner",inner_index_list_tot);
        h5_read(gr,"dagger",is_dagger_list_tot);
      }
      configs_hfile.close();
      MPI_Scatterv(&tau_list_tot[0],&sendcounts[0],displs,MPI_UNSIGNED_LONG_LONG,tau_list,length,MPI_UNSIGNED_LONG_LONG,0,comm);
      MPI_Scatterv(&block_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,block_index_list,length,MPI_INTEGER,0,comm);
      MPI_Scatterv(&inner_index_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,inner_index_list,length,MPI_INTEGER,0,comm);
      MPI_Scatterv(&is_dagger_list_tot[0],&sendcounts[0],displs,MPI_INTEGER,is_dagger_list,length,MPI_INTEGER,0,comm);
      for (int i = 0; i < length; ++i) {
        bool is_dagger = (is_dagger_list[i] == 1 ? true : false);
        auto op = op_desc{block_index_list[i],inner_index_list[i],is_dagger,0}; // linear index will be set later, within TRIQS
        time_pt tau = time_pt(tau_list[i],beta);
        paramCTQMC.initial_config.insert({tau,op});
      }
    }
    paramCTQMC.nbins_histo = nbins_histo;
    // Read histograms

    bool verif_histo = true;

    if (rank == 0) {
      for (int i = 0; i < nblocks; ++i) {  // loop over blocks
        if (!verif_histo) break;
        paramCTQMC.hist_insert[labels[i]] = vector<double>(nbins_histo);
        paramCTQMC.hist_remove[labels[i]] = vector<double>(nbins_histo);
        for (int j = 0; j < 2; ++j) {  // loop over moves
          string tag_move;
          if (j == 0)
            tag_move = "insert";
          else
            tag_move = "remove";
          string fname = "histogram_"+tag_move+"_"+tag_block;
          if (!exists(fname)) {
            verif_histo = false;
            break;
          }
          ifstream hist_file(fname);
          for (int k = 0; k < nbins_histo; ++k) {
            double tau,val;
            hist_file >> tau >> val;
            if (j == 0)
              paramCTQMC.hist_insert[labels[i]][k] = val;
            else
              paramCTQMC.hist_remove[labels[i]][k] = val;
          }
        }
      }
    }
    mpi::mpi_broadcast(verif,comm,0);
    paramCTQMC.performance_analysis = !verif_histo;
    if (verif_histo) {
      if (rank == 0) cout << "   == Reading previous histogram == " << endl;
      mpi::mpi_broadcast(paramCTQMC.hist_insert,comm,0);
      mpi::mpi_broadcast(paramCTQMC.hist_remove,comm,0);
    }
    else {
      paramCTQMC.hist_insert = {};
      paramCTQMC.hist_remove = {};
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
    indices_map_t ind_map = {};
    for (int o = 0; o < ndim; ++o) {
      indices_t ind1,ind2;
      if (off_diag) {
        if (nspinor == 1) {
          ind1 = {"up",o};
          ind2 = {"down",o};
        }
        else {
          ind1 = {"tot",o};
          ind2 = {"tot",o+ndim}
        }
      }
      else {
        ind1 = {to_string(o),0};
        ind2 = {to_string(o+ndim),0};
      }
      ind_map.insert({ind1,ind2});
      ind_map.insert({ind2,ind1});
    }

    paramCTQMC.move_global.insert({to_string(0),ind_map});
    paramCTQMC.measure_G_l = leg_measure;
  
    if (rank == 0 && verbo > 0) cout <<"   == Starting Solver [node "<< rank <<"] ==" << endl << endl;
    //Solve!
    solver.solve(paramCTQMC);
  
    if (rank == 0 && verbo > 0) cout <<"   == Reporting ==" << endl << endl;
 
    if (rank == 0 && verbo > 0) cout <<"   == Writing final configuration ==      " << endl << endl;
    // Write all final configurations
    auto config = solver.get_configuration();
    auto tau_0 = time_pt(1,beta);
    int length = config.size();
    vector<uint64_t> tau_list(length);
    vector<int> block_index_list(length);
    vector<int> inner_index_list(length);
    vector<int> is_dagger_list(length);
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
    vector<int> recvcounts(size);
    MPI_Allgather(&length,1,MPI_INTEGER,&recvcounts[0],1,MPI_INTEGER,comm);
    int displs [size] = {0};
    for (int i = 1; i < size; ++i) displs[i] = displs[i-1] + recvcounts[i-1];
    int tot_size = displs[size-1] + recvcounts[size-1];
    vector<uint64_t> tau_list_tot(tot_size);
    vector<int> block_index_list_tot(tot_size);
    vector<int> inner_index_list_tot(tot_size);
    vector<int> is_dagger_list_tot(tot_size);
    MPI_Gatherv(&tau_list[0],length,MPI_UNSIGNED_LONG_LONG,&tau_list_tot[0],&recvcounts[0],displs,MPI_UNSIGNED_LONG_LONG,0,comm);
    MPI_Gatherv(&block_index_list[0],length,MPI_INTEGER,&block_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
    MPI_Gatherv(&inner_index_list[0],length,MPI_INTEGER,&inner_index_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);
    MPI_Gatherv(&is_dagger_list[0],length,MPI_INTEGER,&is_dagger_list_tot[0],&recvcounts[0],displs,MPI_INTEGER,0,comm);

    if (rank == 0) {
    
      h5::file configs_hfile("configs.h5",'w');   // very important to open in write mode, so that previous configs are erased
      h5::group gr = configs_hfile;
      
      h5_write(gr,"ncpus",size);
      h5_write(gr,"size",recvcounts);
      h5_write(gr,"tau",tau_list_tot);
      h5_write(gr,"block",block_index_list_tot);
      h5_write(gr,"inner",inner_index_list_tot);
      h5_write(gr,"dagger",is_dagger_list_tot);
      
      configs_hfile.close();
    }

    // Write all histograms 
    if (!verif_histo) {

      auto histo = solver.get_performance_analysis();
      for (int i = 0; i < 2; ++i) {  // loop over moves
        string tag_move;
        if (i == 0)
          tag_move="insert";
        else
          tag_move="remove";
        for (int j = 0; j < nblocks; ++j) {  // loop over blocks
          auto hist_accept = histo[tag_move+"_length_accepted_"+labels[j]];
          auto hist_prop = histo[tag_move+"_length_proposed_"+labels[j]];
          hist_accept = mpi::all_reduce(hist_accept,comm);
          hist_prop = mpi::all_reduce(hist_prop,comm);
          double s = 0.;
          double step = beta / (nbins_histo - 1.);  // stepsize of the histogram
          for (int k = 0; k < nbins_histo; ++k) {
            hist_accept.data()[k] = max(hist_accept.data()[k] / max(hist_prop.data()[k],1.), 1.e-5);  // prevents division by zero
            double fac = step;
            if (k == 0 || k == nbins_histo-1) fac /= 2.;  // the first and last bin are half-sized
            s += fac * hist_accept.data()[k];
          }
          if (rank == 0) {
            ofstream hist_file;
            hist_file.open("histogram_"+tag_move+"_"+labels[j]);
            for (int k = 0; k < nbins_histo; ++k) {
              hist_file << double(k) * beta / double(nbins_histo - 1) << "\t" << hist_accept.data()[k] / s  << endl;
            }
          }
        }
      }

      for (int iblock = 0; iblock < nblocks; ++iblock) {
        for (int tau = 0; tau < n_tau; ++tau) {
          for (int o = 0; o < siz_block; ++o) {
            int iflavor = iblock;
            if (off_diag) iflavor = o + iblock*ndim;
            for (int oo = 0; oo < siz_block; ++oo) {
              int iflavor1 = iblock;
              if (off_diag) iflavor1 = oo + iblock*ndim;
              gtau[tau+iflavor*n_tau+iflavor1*num_orbitals*n_tau] = (*solver.G_tau)[iblock].data()(tau,o,oo);
            }
          }
        }
      }

    if (rank == 0 && verbo > 0) cout <<"Gtau reported" << endl;
    // Report G(l)
    if (leg_measure) {

      for (int iblock = 0; iblock < nblocks; ++iblock) {
        for (int l = 0; l < n_l; ++l) {
          for (int o = 0; o < siz_block; ++o) {
            int iflavor = iblock;
            if (off_diag) iflavor = o + iblock*ndim;
            for (int oo = 0; oo < siz_block; ++oo) {
              int iflavor1 = iblock;
              if (off_diag) iflavor1 = oo + iblock*ndim;
              gl[l+iflavor*n_l+iflavor1*num_orbitals*n_l] = (*solver.G_l)[iblock].data()(l,o,oo);
            }
          }
        }
      }

    }

    if (rank == 0 && verbo > 0) cout << "Gl reported" << endl;
   
    if (measure_density_matrix) {

      auto h_loc_diag = solver.h_loc_diagonalization();

      many_body_operator_generic<triqs::utility::real_or_complex> N_op;

      auto subspaces = h_loc_diag.n_subspaces();

      auto rho = solver.density_matrix();

      for (int iblock = 0; iblock < nblocks; ++iblock) {
        for (int o = 0; o < siz_block; ++o) {
          int iflavor = iblock;
          if (off_diag) iflavor = o + iblock*ndim;
          N_op = c_dag(labels[iblock],o) * c(labels[iblock],o);
          occ_mat[iflavor] = trace_rho_op(rho,N_op,h_loc_diag);
        }
      }

      if (rank == 0) {
        ofstream occ;
        occ.open("Occupation_numbers_triqs");

        for (int sub : range(subspaces)) {
          auto fock_states = h_loc_diag.get_fock_states(sub);
          int sub_dim = h_loc_diag.get_subspace_dim(sub);
          auto unit_mat = h_loc_diag.get_unitary_matrix(sub);
          auto rho_temp = unit_mat * rho[sub] * dagger(unit_mat);  // transform density matrix to Fock basis
          for (int ind : range(sub_dim)) {
            auto f_state = fock_states[ind];
            auto prob = rho_temp(ind,ind);
            occ << bitset<64>(f_state).to_string().substr(64-num_orbitals) << "\t" << scientific << setprecision(17) << prob << endl;
          }
        }
      }

      Eu = trace_rho_op(rho,Hint,h_loc_diag);

      if (!leg_measure) {     // Get moments of the self-energy

        many_body_operator_generic<triqs::utility::real_or_complex> comm;
        many_body_operator_generic<triqs::utility::real_or_complex> comm2;
        many_body_operator_generic<triqs::utility::real_or_complex> Sinf_op;
        many_body_operator_generic<triqs::utility::real_or_complex> S1_op;

        auto Sinf_mat = matrix_t(siz_block,siz_block);

        for (int iblock = 0; iblock < nblocks; ++iblock) {
          for (int o = 0; o < siz_block; ++o) {
            int iflavor = iblock;
            if (off_diag) iflavor = o + iblock*ndim;
            for (int oo = 0; oo < siz_block; ++oo) {
              int iflavor1 = iblock;
              if (off_diag) iflavor1 = oo + iblock*ndim;
              comm = Hint*c(labels[iblock],o) - c(labels[iblock],o)*Hint;
              Sinf_op = -comm*c_dag(labels[iblock],oo) - c_dag(labels[iblock],oo)*comm;
              comm2 = c_dag(labels[iblock],oo)*Hint - Hint*c_dag(labels[iblock],oo);
              S1_op = comm2*comm + comm*comm2;
              auto Sinf = trace_rho_op(rho,Sinf_op,h_loc_diag);
              auto S1 = trace_rho_op(rho,S1_op,h_loc_diag);
              moments_self_1[iflavor+iflavor1*num_orbitals] = Sinf;
              moments_self_2[iflavor+iflavor1*num_orbitals] = S1;
              Sinf_mat(o,oo) = Sinf;
            }
          }
          Sinf_mat = Sinf_mat * Sinf_mat;
          for (int o = 0; o < siz_block; ++o) {
            int iflavor = iblock;
            if (off_diag) iflavor = o + iblock*ndim;
            for (int oo = 0; oo < siz_block; ++oo) {
              int iflavor1 = iblock;
              if (off_diag) iflavor1 = oo + iblock*ndim;
              moments_self_2[iflavor+iflavor1*num_orbitals] -= Sinf_mat(o,oo);
            }
          }
        }
      }   // not legendre
    }   // measure_density_matrix
    if (rank == 0 && verbo > 0) cout << "   == CTHYB-QMC Process Finished [node "<< rank <<"] =="<< endl << endl;
  }

/********************************************************/
/****************** Functions Used **********************/
/********************************************************/

// Hamiltonian Initialization with precalculate coeff  (J=0 and U(i,j))
many_body_operator_generic<triqs::utility::real_or_complex> init_Hamiltonian(complex<double> *eps,int nflavor,complex<double> *U,
                                                                             bool off_diag,int nspinor,double lambda) {

    many_body_operator_generic<triqs::utility::real_or_complex> H;
    int iblock,iblock1,o,oo;

    for (int i = 0; i < nflavor; ++i) {
      tie(iblock,o) = convert_indexes(i,off_diag,nspinor);
      for (int j = 0; j < nflavor; ++j) {
        tie(iblock1,oo) = convert_indexes(j,off_diag,nspinor);
        H += U[i+j*nflavor] * n(labels[iblock],o) * n(labels[iblock1],oo);
        H += eps[i+j*nflavor] * c_dag(labels[iblock],o) * c(labels[iblock1],oo);
      }
    }
    return 0.5 * lambda * H;
}

many_body_operator_generic<triqs::utility::real_or_complex> init_fullHamiltonian(complex<double> *eps,int nflavor,complex<double> *U,
                                                                                 bool off_diag,int nspinor,double lambda) {

    many_body_operator_generic<triqs::utility::real_or_complex> H;
    int iblock,iblock1,iblock2,iblock3,o,oo,ooo,oooo;

    for (int i = 0; i < nflavor; ++i) {
      tie(iblock,o) = convert_indexes(i,off_diag,nspinor);
      for (int j = 0; j < nflavor; ++j) {
        tie(iblock1,oo) = convert_indexes(j,off_diag,nspinor);
        H += eps[i+j*nflavor] * c_dag(labels[iblock],o) * c(labels[iblock1],oo);
        for (int k = 0; k < nflavor; ++k) {
          tie(iblock2,ooo) = convert_indexes(k,off_diag,nspinor);
          for (int l = 0; l < nflavor; ++l) {
            tie(iblock3,oooo) = convert_indexes(l,off_diag,nspinor);
            H += U[i+j*nflavor+k*nflavor*nflavor+l*nflavor*nflavor*nflavor] * c_dag(labels[iblock],o) *
              c_dag(labels[iblock1],oo) * c(labels[iblock3],oooo) * c(labels[iblock2],ooo);
          }
        }
      }
    }
    return 0.5 * lambda * H;
}

pair<int,int> convert_indexes(int iflavor,bool off_diag,int nspinor) {
  int iblock = iflavor;
  int o = 0;
  if (off_diag && nspinor == 1) {
    iblock = iflavor / ndim;
    o = iflavor - iblock*ndim;
  }
  if (off_diag && nspinor == 2) {
    iblock = 1;
    o = iflavor;
  }
  return make_pair(iblock,o);
}

// Build DLR frequencies
void build_dlr(int wdlr_size,int *ndlr,double *wdlr,double lam,double eps) {
  auto omega = cppdlr::build_dlr_rf(lam,eps);
  *ndlr = omega.size();
  if (*ndlr > wdlr_size) return;
  for (int i = 0; i < (*ndlr); ++i) wdlr[i] = omega[i]; 
}

