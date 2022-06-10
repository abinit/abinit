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
#include <math.h>

#include "abi_gpu_header.h"
#include "cuda_api_error_check.h"

#include "gpu_apply_invovl_inner.h" // definition of struct invovl_kpt_gpu_t

#include <mpi.h>

extern "C" void xmpi_sum_dp_c(double* array, int32_t array_size, MPI_Fint* comm, int32_t* ierr);

//! token used to store memory state.
//! if memory     already allocated : gpu_initialized is 1
//! if memory not already allocated : gpu_initialized is 0 and device memory allocation is needed
static int gpu_initialized = 0;

//! proj, sm1proj, ptp_sm1proj in GPU memory
static double* proj_gpu;
static double* sm1proj_gpu;
static double* ptp_sm1proj_gpu;
static double* temp_proj;

// hamiltonian data on CPU/GPU

//! number of tuple (l,m,n) by type of atomes (should always be smaller than 255 ?)
static uint8_t *nlmn;

//! same functionality as apply_block (in m_invovl.F90) but adapted to GPU
void apply_block_gpu(int32_t cplx, double* mat, int32_t nprojs, int32_t ndat, double* x, double* y, int32_t block_sliced)
{

} // apply_block_gpu


//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_alloc(int32_t proj_dim[3], int32_t ntypat)
{
  // when memory allocation already done, deallocate first before re-allocating
  if (gpu_initialized==1)
    gpu_apply_invovl_inner_dealloc();
  else
    gpu_initialized = 1;

  int32_t proj_size = proj_dim[0]*proj_dim[1]*proj_dim[2];
  CHECK_CUDA_ERROR( cudaMalloc((void**)&proj_gpu, proj_size) );
  CHECK_CUDA_ERROR( cudaMalloc((void**)&sm1proj_gpu, proj_size) );
  CHECK_CUDA_ERROR( cudaMalloc((void**)&ptp_sm1proj_gpu, proj_size) );
  CHECK_CUDA_ERROR( cudaMalloc((void**)&temp_proj, proj_size) );

  nlmn =   (uint8_t *)  malloc(ntypat * sizeof(uint8_t) );

} // gpu_apply_invovl_inner_alloc

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_dealloc()
{
  gpu_initialized = 0;

  CHECK_CUDA_ERROR( cudaFree(proj_gpu) );
  CHECK_CUDA_ERROR( cudaFree(sm1proj_gpu) );
  CHECK_CUDA_ERROR( cudaFree(ptp_sm1proj_gpu) );
  CHECK_CUDA_ERROR( cudaFree(temp_proj) );

  free(nlmn);

} // gpu_apply_invovl_inner_dealloc


//! initialize arrays nlmn
extern "C" void init_invovl_data(int32_t indlmn_dim[3], int32_t* indlmn)
{

  assertm(indlmn_dim[0] == 6, "indlmn array must have first dimension equal to 6 !");

  const int32_t qmax   = indlmn_dim[0];
  const int32_t lmnmax = indlmn_dim[1];
  const int32_t ntypat = indlmn_dim[2];

  for (int itypat=0; itypat<ntypat; ++itypat) {

    uint8_t count_ilmn = 0;

    for (int ilmn=0; ilmn<lmnmax; ilmn++) {
      if(indlmn[3 + qmax * (ilmn + lmnmax * itypat)] > 0) {
        // if count_ilmn == 255 we have a problem here
        count_ilmn++;
      }
    }

    nlmn[itypat] = count_ilmn;

  }

} // init_invovl_data

//! apply inverse overlap operator (inner part of the computation)
//!
//!
//! \param[in] invovl data
//! \param[in] proj_dim gives the dimension of proj, sm1proj and ptp_sm1proj
//! \param[in] proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[inout] sm1proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[inout] ptp_sm1proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[in] ndatspinor is the product of ndat and nspinor
//! \param[in] block_sliced (can only be 0 or 1), switch used inside apply_block_gpu
extern "C" void solve_inner_gpu(invovl_kpt_gpu_t* invovl,
                                int32_t proj_dim[3],
                                double* proj,
                                double* sm1proj,
                                double* ptp_sm1proj,
                                int32_t nattyp_dim,
                                int32_t* nattyp,
                                int32_t cplx,
                                int32_t block_sliced)
{

  // projections norm array (for MPI all reduce)
  double* normprojs; // size is ndat*nspinor;

  // number of projections
  int32_t nprojs = invovl->nprojs;

  // MPI error return value
  int ierr;

  int ibeg = 0;
  int iend = nprojs;
  int nlmntot_this_proc = nprojs;

  int ndat = proj_dim[2];

  // TO BE REMOVED
  printf("INSIDE solve_inner_gpu\n");

  // compute normproj by summing over the first two dimensions
  normprojs = (double *) malloc(ndat*sizeof(double));

  for (int idat=0; idat<ndat; ++idat) {
    double norm = 0.0;
    int slice = proj_dim[0]*proj_dim[1];
    for (int index = idat*slice; index<(idat+1)*slice; ++index) {
      norm += proj[index]*proj[index];
    }
    normprojs[idat] = norm;
  }

  // first guess for sm1proj
  // mat => inv_s_approx
  // x => proj
  // y => sm1proj
  // compute sm1proj = inv_s_approx * proj
  //apply_block_gpu(cplx, invovl, ndat, proj_gpu, sm1proj_gpu, block_sliced)

  // Iterative refinement
  // TODO use a more efficient iterative algorithm than iterative refinement, use locking

  double maxerr, previous_maxerr;
  const double precision = 1e-16; // maximum relative error. TODO: use tolwfr ?
  double convergence_rate;
  int additional_steps_to_take = -1;

  for (int i=1; i < 30; ++i) {

    // compute resid = proj - (D^-1 + PtP)sm1proj
    // TODO apply_block_gpu(ham, cplx, invovl%inv_sij, nprojs, ndat, sm1proj, resid, block_sliced);
    // TODO temp_proj = sm1proj(:,ibeg:iend,:);

    // compute matrix multiplication : PtPsm1proj(:,:,1) = invovl%gram * temp_proj(:,:,1)
    // TODO abi_xgemm('N', 'N', nprojs, ndat, nlmntot_this_proc, cone, &
    // & invovl%gram_projs(:,:,1), nprojs, &
    // & temp_proj(:,:,1), nlmntot_this_proc, czero, &
    // & PtPsm1proj(:,:,1), nprojs, &
    // & x_cplx=cplx)
    // TODO resid = proj - resid - Ptpsm1proj

    // exit check
    // TODO errs = SUM(SUM(resid**2, 1),1)

    // TODO maxerr = sqrt(MAXVAL(errs/normprojs))
    if(maxerr < precision or additional_steps_to_take == 1) {
      exit(EXIT_FAILURE);
      // We might stall and never get to the specified precision because of machine errors.
      // If we got to 1e-10, extrapolate convergence rate and determine the number of additional
      // steps to take to reach precision
    } else if (maxerr < 1e-10 and additional_steps_to_take == -1) {
      convergence_rate = -log(1e-10) / i;
      additional_steps_to_take = ceil(-log(precision/1e-10)/convergence_rate) + 1;
    } else if (additional_steps_to_take > 0) {
      if (previous_maxerr<maxerr)
        exit(EXIT_FAILURE);
      additional_steps_to_take -= 1;
    }
    previous_maxerr=maxerr;

    // add preconditionned residual
    // TODO call apply_block(ham, cplx, invovl%inv_s_approx, nprojs, ndat, resid, precondresid, block_sliced)
    // TODO sm1proj = sm1proj + precondresid
  }


  free(normprojs);

} // solve_inner_gpu
