/* abi_gpu_header.h */

/*
 * Copyright (C) 2008-2020 ABINIT Group
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef GPU_FOUR_HEADER_H
#define GPU_FOUR_HEADER_H


#include "config.h"
#include "cuda_header.h"

#define BLOCK_SIZE 128  // CUDA block size
#define MAX_GRID_SIZE 65535  // CUDA max grid size by dimension

//Interfaces
#ifdef __cplusplus
extern "C" {
#endif

  void abi_cabort();


  //Headers for Fourwf :
  void alloc_gpu_fourwf_(int *ngfft,int *ndat,int *npwin,int *npwout);
  void free_gpu_fourwf_();
  void gpu_fourwf_(int *cplex,double *denpot,double *fofgin,double *fofgout,double *fofr,int *gboundin,int *gboundout,int *istwf_k,
                   int *kg_kin,int *kg_kout,int *mgfft,void *mpi_enreg,int *ndat,int *ngfft,int *npwin,int *npwout,int *n4,int *n5,int *n6,int *option,
                   int *paral_kgb,int * tim_fourwf,double *weight_r,double *weight_i);/*,
                   int *use_ndo,double *fofginb)*/
  //Put wf function from cg to sphere in cfft
  void gpu_sphere_in_(double *cg,double *cfft,int *kg_k,int *npw,int *n1,int *n2,int *n3,int* ndat,int *istwfk,cudaStream_t *compute_stream);

  //Extract wf functions from sphere in cfft to cg
  void gpu_sphere_out_(double *cg,double *cfft,int *kg_k,int *npw,int *n1,int *n2,int *n3,int* ndat,cudaStream_t *compute_stream);

  void gpu_density_accumulation_(double *fofr,double* denpot, double* weight_r,double* weight_i,int* nfft_tot,int *ndat,cudaStream_t *compute_stream);

  void gpu_apply_local_potential_(double2 *fofr,double* denpot,int* nfft_tot,int *ndat,cudaStream_t *compute_stream);


  //Headers for nonlop:
 void alloc_nonlop_gpu_(int *npwin,int *npwout,int *nspinor,
                        int *natom,int *ntypat,int *lmnmax,
                        int *indlmn,int *nattyp,
                        int *atindx1,double *gprimd,
                        int *dimffnlin,int *dimffnlout,
                        int *dimenl1, int *dimenl2 );

  void free_nonlop_gpu_();

  void gpu_update_ham_data_(double *enl,int *size_enl, double *sij,int *size_sij,
                            double *gprimd,int *size_gprimd);

  void gpu_update_ffnl_ph3d_(double *ph3din,int *dimph3din,double *ffnlin,int *dimffnlin);

  void gpu_finalize_ham_data_();

  void gpu_finalize_ffnl_ph3d_();

  void gpu_compute_nl_hamiltonian_(double2 *proj_gpu,double2 *dproj_gpu,
                                   double2 *vectin_gpu,
                                   double2 *vectout_gpu,double2 *svectout_gpu,
                                   double2 *val_ajlmn_gpu,double2 *val_sajlmn_gpu,
                                   double2 *ph3dout_gpu,double *ffnlout_gpu,
                                   double *rdlmn_gpu, double *rdproj_gpu,
                                   double *enlout_gpu, double *d_enlk_gpu,
                                   double *enl_gpu,double *sij_gpu,double *gprimd_gpu,
                                   int *atindx1_gpu,unsigned char *nlmn_gpu,
                                   int *indlmn_gpu,unsigned short int *atoms_gpu,
                                   unsigned char *lmn_gpu, unsigned char *typat_gpu,
                                   int *paw_opt,int *dimffnlout,int *dimenl1,
                                   int *npw,int *nb_projections,int *lmnmax,
                                   int *natoms, int *choice, int *signs,
                                   const double *four_pi_by_ucvol,const double *lambda
                                   );
  
  void gpu_compute_nl_projections_(double2 *proj_gpu,double2 *dproj_gpu,
                                   double2 *vectin_gpu,double2 *ph3din_gpu,
                                   double *ffnlin_gpu,double *kpgin_gpu,
                                   int *indlmn_gpu,unsigned short int *atoms_gpu,
                                   unsigned char *lmn_gpu, unsigned char *typat_gpu,
                                   int *nb_proj_to_compute,int *npw,int *choice,
                                   int *dimffnlin,int *lmnmax,
                                   const char *cplex,const double *pi,const double *ucvol
                                   );
  void gpu_mkkpg_(int *kg_gpu,double *kpg_gpu,double *kpt,int *npw);


#ifdef __cplusplus
}
#endif

#endif
