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

//! token used to store memory state.
//! if memory     already allocated : gpu_initialized is 1
//! if memory not already allocated : gpu_initialized is 0 and device memory allocation is needed
static int gpu_initialized = 0;

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_alloc()
{
  // when memory allocation already done, deallocate first before re-allocating
  if (gpu_initialized==1)
    gpu_apply_invovl_inner_dealloc();
  else
    gpu_initialized = 1;

} // gpu_apply_invovl_inner_alloc

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_dealloc()
{
    gpu_initialized = 0;

} // gpu_apply_invovl_inner_dealloc

//! apply inverse overlap operator (inner part of the computation)
extern "C" void solve_inner_gpu(invovl_kpt_gpu_t* invovl)
{

  printf("INSIDE solve_inner_gpu\n");

  // TODO

} // solve_inner_gpu_
