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

//! token used to store memory state.
//! if memory     already allocated : gpu_initialized is 1
//! if memory not already allocated : gpu_initialized is 0 and device memory allocation is needed
static int gpu_initialized = 0;

static double* proj_gpu;

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_alloc(int32_t proj_dim[3])
{
  // when memory allocation already done, deallocate first before re-allocating
  if (gpu_initialized==1)
    gpu_apply_invovl_inner_dealloc();
  else
    gpu_initialized = 1;

  int32_t proj_size = proj_dim[0]*proj_dim[1]*proj_dim[2];
  CHECK_CUDA_ERROR( cudaMalloc((void**)&proj_gpu, proj_size) );

} // gpu_apply_invovl_inner_alloc

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_dealloc()
{
  gpu_initialized = 0;

  CHECK_CUDA_ERROR( cudaFree(proj_gpu) );

} // gpu_apply_invovl_inner_dealloc

//! apply inverse overlap operator (inner part of the computation)
//! \param[in] invovl data
//! \param[in] proj is a 3D array of size (2,nprojs, nspinor*ndat)
//! \param[in] MPI communicator (from fortran)
extern "C" void solve_inner_gpu(invovl_kpt_gpu_t* invovl,
                                double* proj,
                                MPI_Fint f_comm_fft, int32_t nproc_fft)
{
  // retrieve MPI communicator
  MPI_Comm comm_fft = MPI_Comm_f2c(f_comm_fft);

  int32_t nprojs = invovl->nprojs;

  printf("INSIDE solve_inner_gpu\n");

  // TODO PK

} // solve_inner_gpu
