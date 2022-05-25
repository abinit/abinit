//****f* ABINIT/gpu_apply_invovl_inner
//
// NAME
// gpu_apply_invovl_inner
//
// FUNCTION
// Helper function: iteratively solves the inner system.
// This is a direct equivalent of solve_inner in 66_wfs/m_invovl.F90
//
// INPUTS
//
// PARENTS
//      m_invovl
//
// CHILDREN
//      dsymm,zhemm
//
// SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "abi_gpu_header.h"
#include "cuda_api_error_check.h"

#include "gpu_apply_invovl_inner.h" // definition of struct invovl_kpt_gpu_t

#include <mpi.h>

extern "C" void xmpi_sum_dp_c(double* array, int32_t array_size, MPI_Fint* comm, int32_t* ierr);

//! token used to store memory state.
//! if memory     already allocated : gpu_initialized is 1
//! if memory not already allocated : gpu_initialized is 0 and device memory allocation is needed
static int gpu_initialized = 0;

static double* proj_gpu;
static int* array_nlmntot_pp;


//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_alloc(int32_t proj_dim[3], int32_t nproc_fft)
{
  // when memory allocation already done, deallocate first before re-allocating
  if (gpu_initialized==1)
    gpu_apply_invovl_inner_dealloc();
  else
    gpu_initialized = 1;

  int32_t proj_size = proj_dim[0]*proj_dim[1]*proj_dim[2];
  CHECK_CUDA_ERROR( cudaMalloc((void**)&proj_gpu, proj_size) );

  array_nlmntot_pp = (int *) malloc(nproc_fft);

} // gpu_apply_invovl_inner_alloc

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_dealloc()
{
  gpu_initialized = 0;

  CHECK_CUDA_ERROR( cudaFree(proj_gpu) );

  free(array_nlmntot_pp);

} // gpu_apply_invovl_inner_dealloc

//! apply inverse overlap operator (inner part of the computation)
//! \param[in] invovl data
//! \param[in] proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[in] ndatspinor is the product of ndat and nspinor
//! \param[in] MPI communicator (from fortran)
extern "C" void solve_inner_gpu(invovl_kpt_gpu_t* invovl,
                                double* proj,
                                int32_t proj_dim[3],
                                MPI_Fint* f_comm_fft,
                                int32_t me_fft,
                                int32_t nproc_fft,
                                int32_t paral_kgb)
{

  // projections norm array (for MPI all reduce)
  double* normprojs; // size is ndat*nspinor;

  // retrieve MPI communicator
  MPI_Comm comm_fft = MPI_Comm_f2c(*f_comm_fft);

  // number of projections
  int32_t nprojs = invovl->nprojs;

  // MPI error return value
  int ierr;

  int ibeg, iend, nlmntot_this_proc;

  // TO BE REMOVED
  printf("INSIDE solve_inner_gpu\n");

  // compute normproj by summing over the first two dimensions
  normprojs = (double *) malloc(proj_dim[2]*sizeof(double));
  for (int i=0; i<proj_dim[2]; ++i) {
    double norm = 0.0;
    int slice = proj_dim[0]*proj_dim[1];
    for (int j = i*slice; j<(i+1)*slice; ++j) {
      norm += proj[j]*proj[j];
    }
    normprojs[i] = norm;
  }
  xmpi_sum_dp_c(normprojs, proj_dim[2], f_comm_fft, &ierr);

  // Compute work distribution : split nprojs evenly between the fft processors
  if (paral_kgb == 1) {

    for (int i = 0; i < nproc_fft; ++i) {
      if (i< (nprojs % nproc_fft))
        array_nlmntot_pp[i] = nprojs / nproc_fft + 1;
      else
        array_nlmntot_pp[i] = nprojs / nproc_fft;
    }

    // watchout fortran => C
    ibeg = 0;
    for (int i = 0; i < me_fft; ++i) {
      ibeg += array_nlmntot_pp[i];
    }
    iend = ibeg + array_nlmntot_pp[me_fft];
    nlmntot_this_proc = iend - ibeg;

  } else {
    ibeg = 0;
    iend = nprojs;
    nlmntot_this_proc = nprojs;
  }

  free(normprojs);

} // solve_inner_gpu
