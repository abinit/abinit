#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <solver_core.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/real_or_complex.hpp>
#include <triqs/gfs.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#if defined HAVE_TRIQS_v1_4
#include <triqs/mpi/boost.hpp>
#elif defined HAVE_TRIQS_v2_0
#include <triqs/h5.hpp>
#endif

#include "triqs_cthyb_qmc.hpp"

using namespace std;
#if defined HAVE_TRIQS_v2_0
using namespace triqs_cthyb;
#elif defined HAVE_TRIQS_v1_4
using namespace cthyb;
#endif
using namespace triqs::gfs;
using triqs::operators::many_body_operator_generic;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;

#include <mpi.h>

void ctqmc_triqs_run(bool rot_inv, bool leg_measure, bool hist,     /*boolean*/
                     bool wrt_files, bool tot_not,                  /*boolean*/	      
                     int n_orbitals, int n_freq, int n_tau,         /*integer*/
                     int n_l, int n_cycles_, int cycle_length,      /*integer*/
                     int ntherm, int verbo, int seed,               /*integer*/
		             double beta_,                                  /*double*/
                     double *epsi, double *umat_ij,                 /*pointers*/
                     double *umat_ijkl,                             /*pointers*/
                     std::complex<double> *f_iw_ptr,                /*pointers*/
                     std::complex<double> *g_iw_ptr,                /*pointers*/
                     double *g_tau, double *gl,                     /*pointers*/
#if defined HAVE_TRIQS_v1_4
                     MPI_Fint *MPI_world_ptr ){                     /*pointers*/
#else
                     int rank ){                     /*pointers*/
#endif

    cout.setf(ios::fixed);
#if defined HAVE_TRIQS_v1_4
    //Initialize Boost mpi environment
    int rank, nprocs;
    boost::mpi::environment env;
    {
        boost::mpi::communicator c;
        c << MPI_Comm_f2c( *MPI_world_ptr );
        rank=c.rank();

        MPI_Comm_size(c, &nprocs);
        std::cout << "Number of processors: " << nprocs << endl;
    }
#else
    std::cout << "Rank: "<< rank << endl;
#endif

    // Parameters relay from Fortran and default values affectation
    double beta = beta_;                //Temperature inverse 
    int num_orbitals = n_orbitals;
    //double U = U_;                    //Impurity site Interaction term
    double mu = 0.0;                    //Chemical Potential Functionality avail but not used here
    //bool n_n_only = !(rot_inv);       //density-density?
  
    int Nfreq = n_freq;                 //Number of Matsubara frequencies 
    int Ntau = n_tau;                   //Number of imaginary times  
    int Nl = n_l;                       //Number of Legendre measures 30
    int n_cycles = n_cycles_;
    int nspin = 2;                      //always it's Nature
    
    //Print informations about the Solver Parameters 
    if(rank==0 && verbo>0){
        std::cout << endl <<"   == Key Input Parameters for the TRIQS CTHYB solver ==	"<<endl<<endl;
    
        std::cout <<"	Beta               = "<< beta <<" [eV]^(-1) "<<endl;
        std::cout <<"	Mu                 = "<< mu <<" [eV]"<<endl;
        std::cout <<"	Nflavour           = "<< num_orbitals <<endl;
        std::cout <<"	Nfreq              = "<< Nfreq << endl;
        std::cout <<"	Ntau               = "<< Ntau << endl;
        std::cout <<"	Nl                 = "<< Nl << endl;
        std::cout <<"	Ncycles            = "<< n_cycles << endl;
        std::cout <<"	Cycle length       = "<< cycle_length << endl;
        std::cout <<"	Ntherm             = "<< ntherm << endl;
        std::cout <<"	Verbosity (0->3)   = "<< verbo << endl;
        std::cout <<"	Abinit Seed        = "<< seed << endl <<endl;
    
        for(int o = 0; o < num_orbitals; ++o)
	        std::cout << "	e- level["<< o+1 <<"] = "<< setprecision(17) << epsi[o] <<" [eV]"<< endl; //ed != cste => vecteur a parcourir ensuite
    
        std::cout << endl;
        std::cout <<"	U(i,j) [eV] = "  <<endl<< endl<<"\t";
    
        for(int o = 0; o < num_orbitals; ++o){
            for(int oo = 0; oo < num_orbitals; ++oo)
	            std::cout << setprecision(17) << fixed <<umat_ij[o+oo*num_orbitals] << "\t"; //ed != cste => vecteur a parcourir ensuite
            std::cout << endl<<"\t";
        }
    }

    // Hamiltonian definition
    many_body_operator_generic<double> H;
    
    // Spin Orbit general case num_orbitals*2=Nflavour
   
    // Init Hamiltonian Basic terms
    if(!rot_inv){
	    if(rank==0 && verbo>0) std::cout <<endl<<"   == Density-Density Terms Included ==	"<<endl<<endl;  
        H = init_Hamiltonian( epsi, num_orbitals, umat_ij );

        // Include density-density term (J != 0) spin flips and pair hopping => density-density term case
    }else{ // if(rot_inv){
	    if(rank==0 && verbo>0)std::cout <<"   == Rotationnaly Invariant Terms Included ==	"<< endl << endl;
	
	    if(tot_not)
            H = init_fullHamiltonian( epsi, num_orbitals, umat_ijkl );
	    else // up and down separation
	        H = init_fullHamiltonianUpDown( epsi, num_orbitals, umat_ijkl );
    }// else if(rank==0 && verbo>0) std::cout <<endl<<"   == Density-Density Terms Not Included ==	"<<endl<<endl;  

#if defined HAVE_TRIQS_v2_0
    std::vector<std::pair<std::string, indices_type>> gf_struct;

    auto idx_vec = indices_type{};
    for(int o = 0; o < num_orbitals; ++o)
        idx_vec.push_back(o);
  
    if( tot_not )
        gf_struct.emplace_back("tot", idx_vec);
    else{// spin orb case: general case
        gf_struct.emplace_back("up", idx_vec);
        gf_struct.emplace_back("down", idx_vec);
    }
#elif defined HAVE_TRIQS_v1_4
    std::map<std::string, indices_type> gf_struct;

    if( tot_not ){
        std::map<std::string, indices_type> gf_struct_tmp{{"tot", indices_type{}}};

        for(int o = 0; o < num_orbitals; ++o)
            gf_struct_tmp["tot"].push_back(o);
    
        gf_struct = gf_struct_tmp;
    }else{ // spin orb case: general case:
        std::map<std::string, indices_type> gf_struct_tmp{{"up", indices_type{}}, {"down", indices_type{}}};

        for(int o = 0; o < num_orbitals; ++o){
            gf_struct_tmp["up"].push_back(o);
            gf_struct_tmp["down"].push_back(o);
        }

        gf_struct = gf_struct_tmp;
    }
#endif

    if(rank==0 && verbo>0) std::cout <<"   == Green Function Structure Initialized ==	"<< endl << endl;

    // Construct CTQMC solver with mesh parameters
#if defined HAVE_TRIQS_v2_0
    solver_core solver({beta, gf_struct, Nfreq, Ntau, Nl});
#elif defined HAVE_TRIQS_v1_4
    solver_core solver(beta, gf_struct, Nfreq, Ntau, Nl);
#endif
    if(rank==0 && verbo>0) std::cout <<"   == Solver Core Initialized ==	"<< endl << endl;

    //Fill in "hybridation+eps"=F(iw)~Delta(iw) coming from fortran inside delta_iw term RHS
    gf<imfreq> delta_iw = gf<imfreq>{{beta, Fermion, Nfreq}, {num_orbitals,num_orbitals}};
    auto w_mesh = delta_iw.mesh();
 
    for(std::size_t w1 = 0; w1 < w_mesh.size(); ++w1){
        int w_index = w1 - w_mesh.size()/2;
        int w = w_index;
        if(w_index < 0) w = -1*w_index - 1; //Maps negative w to correct positive w

        auto iw = w_mesh.index_to_point(w_index);
        auto cell = delta_iw[w_index];

        for(int o = 0; o < num_orbitals; ++o){
	        for(int oo = 0; oo < num_orbitals; ++oo){
                auto f_iw = f_iw_ptr[o+oo*num_orbitals+w*num_orbitals*num_orbitals];
                if(w_index < 0) f_iw = std::conj(f_iw);
                cell(o,oo) = f_iw;
	        }
	    }
    }
    if(rank==0 && verbo>0) std::cout <<"   == F(iw) Initialized ==	"<< endl << endl;

    triqs::clef::placeholder<0> om_;

    auto g0_iw = gf<imfreq>{{beta, Fermion, Nfreq}, {num_orbitals,num_orbitals}};

    g0_iw(om_) <<  om_  + mu - delta_iw(om_); //calculate G0(iw)^-1 matrix

    if(rank==0 && verbo>0) std::cout <<"   == G0(iw)^-1 Initialized ==	"<< endl << endl;
  
    if(tot_not)
        solver.G0_iw()[0] = triqs::gfs::inverse( g0_iw ); //inverse G0(iw) matrix and affect to solver
    else{
        for (int i = 0; i < nspin; ++i)
            solver.G0_iw()[i] = triqs::gfs::inverse( g0_iw );
    }
  
    if(rank==0 && verbo>0) std::cout <<"   == G0(iw)^-1 Inverted => G0(iw) Constructed ==	"<< endl << endl;
  
    std::cout << setprecision(17) << fixed;

    // Solver parameters
    if(rank==0 && verbo>0) std::cout <<"   == Solver Parametrization ==	"<< endl << endl;
    auto paramCTQMC = solve_parameters_t(H, n_cycles);
    paramCTQMC.max_time = -1;
    paramCTQMC.random_name = "";
    paramCTQMC.random_seed = 123 * rank + 567;// seed;
    paramCTQMC.length_cycle = cycle_length;
    paramCTQMC.n_warmup_cycles = ntherm;
    paramCTQMC.verbosity=verbo;
#if defined HAVE_TRIQS_v2_0
    paramCTQMC.measure_G_l = leg_measure;
#elif defined HAVE_TRIQS_v1_4
    paramCTQMC.measure_g_l = leg_measure;
#endif
    // paramCTQMC.move_double = true; 
    // paramCTQMC.make_histograms=hist;
  
    if(rank==0 && verbo>0) std::cout <<"   == Starting Solver [node "<< rank <<"] ==	"<< endl << endl;
    //Solve!
    solver.solve(paramCTQMC);
  
    if(rank==0 && verbo>0) std::cout <<"   == Reporting ==	"<< endl << endl;
    //Report some stuff
    int compteur = 0;

    std::cout << endl <<" === solver.G0_iw === "<< endl;

    for(int oo = 0; oo < num_orbitals; ++oo){      
        for(int o = 0; o < num_orbitals; ++o){ 
            for(std::size_t w_index = w_mesh.size()/2; w_index < w_mesh.size(); ++w_index){
	            g_iw_ptr[compteur] = solver.G0_iw()[0].data()(w_index,o,oo);
	            compteur++;
	        }
        }
    }
  
    // Report G_tau
    compteur = 0;

    for(int oo = 0; oo < num_orbitals; ++oo){ 
        for(int o = 0; o < num_orbitals; ++o){ 
            for(int tau = 0; tau < n_tau; ++tau){
#if defined HAVE_TRIQS_v2_0
                g_tau[compteur]= (*solver.G_tau)[0].data()(tau,o,oo).real();
#elif defined HAVE_TRIQS_v1_4
                g_tau[compteur]= solver.G_tau()[0].data()(tau,o,oo).real();
#endif
                compteur++;
            }
        }  
    }

    // Report G(l)
    if( leg_measure ){
        compteur = 0;
        for(int oo = 0; oo < num_orbitals; ++oo){ 
            for(int o = 0; o < num_orbitals; ++o){ 
                for(int l = 0; l < n_l; ++l){
#if defined HAVE_TRIQS_v2_0
        	        gl[compteur]= (*solver.G_l)[0].data()(l,o,oo).real();
#elif defined HAVE_TRIQS_v1_4
        	        gl[compteur]= solver.G_l()[0].data()(l,o,oo).real();
#endif
        	        compteur++;
                }
            }
        }
    }

  
    // Write G(tau)
    if( rank==0 && wrt_files ){
    
        ofstream gtau;
        gtau.open("Gtau_triqs.dat");
        double _tau_=0.0;
  
        std::cout <<" gtau "<<endl;
        for(int tau = 0; tau < n_tau; ++tau){
            _tau_=((tau*beta)/n_tau)*27.2113845; //in Harthree
            gtau << fixed << setprecision(17) <<_tau_ << "\t";

            for(int o = 0; o < num_orbitals; ++o){ 	    
	            for(int oo = 0; oo < num_orbitals; ++oo)
#if defined HAVE_TRIQS_v2_0
 	                gtau << fixed<< setprecision(17) << (*solver.G_tau)[0].data()(tau,o,oo) <<"\t";
#elif defined HAVE_TRIQS_v1_4
 	                gtau << fixed<< setprecision(17) << solver.G_tau()[0].data()(tau,o,oo) <<"\t";
#endif
	        }	  
 	        gtau << endl;
	  
        }

        ofstream g_iw;
        g_iw.open("giw");
        for(std::size_t w_index = 0; w_index < w_mesh.size(); ++w_index){
            auto iw = w_mesh.index_to_point(w_index);
            for(int o = 0; o < num_orbitals; ++o){
                for(int oo = 0; oo < num_orbitals; ++oo)
	                g_iw << fixed << setprecision(17) << solver.G0_iw()[0].data()(w_index,o,oo);
            }
            g_iw << endl;
        }
    }
  
    // Write U(i,j)
    if( rank==0 && wrt_files ){
        ofstream umat;
        umat.open("umat");

        for(int o = 0; o < num_orbitals; ++o){ 	    
	        for(int oo = 0; oo < num_orbitals; ++oo)
 	            umat <<"U( "<< o+1 <<" , "<< oo+1 <<" )= "<< fixed << setprecision(17) << umat_ij[o+(oo*num_orbitals)] <<endl;
	    }	  
    } 

    if( leg_measure && rank==0 && wrt_files ){
        ofstream g_l;
        g_l.open("gl");
  
        for(int l = 0; l < n_l; ++l){
            for(int o = 0; o < num_orbitals; ++o){ 
                for(int oo = 0; oo < num_orbitals; ++oo)
#if defined HAVE_TRIQS_v2_0
                    g_l << fixed << setprecision(17) << (*solver.G_l)[0].data()(l,o,oo) <<"\t";	      
#elif defined HAVE_TRIQS_v1_4
                    g_l << fixed << setprecision(17) << solver.G_l()[0].data()(l,o,oo) <<"\t";	      
#endif
            }	  
            g_l << endl;
        }
    }

    if(rank==0 && verbo>0) std::cout <<"   == CTHYB-QMC Process Finished [node "<< rank <<"] ==	"<< endl << endl;
}

/********************************************************/
/****************** Functions Used **********************/
/********************************************************/

// Hamiltonian Initialization with precalculate coeff  (J=0 and U(i,j))
many_body_operator_generic<double> init_Hamiltonian( double *eps, int nflavor, double *U){
    many_body_operator_generic<double> H; double coeff;
  
    for(int i = 0; i < nflavor; ++i)
        for(int j = 0; j < nflavor; ++j){
  
    	    //cas diago principale:
	        if(i==j)
	            H += eps[i] * n("tot",i);
    	    else{
	    	    coeff = double(U[i+j*nflavor])*0.5;
		        H += coeff * n("tot",i) * n("tot",j);
    	    }
        }
    return H;
}

// Hamiltonian Initialization with precalculated coeff U(i,j,k,l) in "tot" notation
many_body_operator_generic<double> init_fullHamiltonian( double *eps, int nflavor, double *U){
    many_body_operator_generic<double> H; double coeff;
  
    for(int i = 0; i < nflavor; ++i)
        for(int j = 0; j < nflavor; ++j){
      
            if(i==j)
	            H += eps[i] * n("tot",i) ; //Levels On diagonal Hamiltonian Matrix
      
            for(int k = 0; k < nflavor; ++k)
                for(int l = 0; l < nflavor; ++l){ 
                    //cas diago principale:
		            coeff = 0.5 * double( U[i+j*nflavor+k*nflavor*nflavor+l*nflavor*nflavor*nflavor] );
		            H += coeff * c_dag("tot",i)*c_dag("tot",j)*c("tot",l)*c("tot",k); //changing l <-> k order index to be in accordance with the Abinit definition
	            }
        }
    return H;
}

// Hamiltonian Initialization with precalculated coeff U(i,j,k,l) in "up" and "down" notation
many_body_operator_generic<double> init_fullHamiltonianUpDown( double *eps, int nflavor, double *U ){
    many_body_operator_generic<double> H; double coeff;
  
    for(int i = 0; i < nflavor; ++i)
        for(int j = 0; j < nflavor; ++j){
      
            if(i==j)
	            H += eps[i] * ( n("up",i) + n("down",i) ); //Levels On diagonal Hamiltonian Matrix
      
            for(int k = 0; k < nflavor; ++k)
                for(int l = 0; l < nflavor; ++l){
                    coeff = 0.5 * double(U[i+j*nflavor+k*nflavor*nflavor+l*nflavor*nflavor*nflavor]);
 
//            for(int s = 0; s < nspin; ++s)
//              for(int s1 = 0; s1 < nspin; ++s1){ 
		            H += coeff * c_dag("up",i)   * c_dag("down",j) * c("up",l)   * c("down",k);
		            H += coeff * c_dag("down",i) * c_dag("up",j)   * c("down",l) * c("up",k);
	        }
    }
    return H;
}
