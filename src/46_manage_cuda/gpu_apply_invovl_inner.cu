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

#include <gpu_linalg.h>
#include "abi_gpu_header.h"
#include "cuda_api_error_check.h"

#include "gpu_apply_invovl_inner.h" // definition of struct invovl_kpt_gpu_t

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

#include <mpi.h>

extern "C" void xmpi_sum_dp_c(double* array, int32_t array_size, MPI_Fint* comm, int32_t* ierr);

//! token used to store memory state for buffers that need to be uploaded each time
//! apply_block is called.
//!
//! if memory     already allocated : gpu_inner_allocated is 1
//! if memory not already allocated : gpu_inner_allocated is 0 and device memory allocation is needed
static int gpu_inner_allocated = 0;

//! proj, sm1proj, ptp_sm1proj in GPU memory
static double* proj_gpu;
static double* sm1proj_gpu;
static double* ptp_sm1proj_gpu;
static double* temp_proj_gpu;
static double* resid_gpu;
static double* precondresid_gpu;

//! token used to store memory state for buffers that need to be uploaded each time
//! make_block is called.
static int gpu_inv_overlap_matrix_allocated = 0;

//! inverse overlap matrix on GPU
static double* gram_projs_gpu;
static double* inv_sij_gpu;
static double* inv_s_approx_gpu;

// hamiltonian data on CPU/GPU

//! number of tuple (l,m,n) by type of atomes (should always be smaller than 255 ?)
static uint8_t *nlmn;

//! same functionality as apply_block (in m_invovl.F90) but adapted to GPU
//! compute sort of y = mat * x
//!
//! x is rectangular matrix, each column correspond to a wave vector.
//! each wave vector is actually divided in pieces (one piece by atome type), and the piece
//! is recast in a 2d matrix of size (nlmn[itypat],nattyp[itypat]), called the "x matrix"
//! in cublas API.
//!
//! We use the cublas batch API to repeat the matrix multiplication for all wave vectors (ndat).
//!
//! \param[in] cplx (1 => real, 2 => complex data)
//! \param[in] ntypat is the number of types of atoms
//! \param[in] nattyp is the array (of size ntypat) which gives the number of atoms of a given type
//! \param[in] mat is square matrix of size lmnmax (and leading dimension if also lmnmax)
//! \param[in] x is a 2d array of size (nprojs, ndat)
//! \param[inout] y is a 2d array of size (nprojs, ndat)
void apply_block_gpu(int32_t cplx,
                     int32_t ntypat, int32_t* nattyp, int32_t lmnmax,
                     const double* mat, int32_t nprojs, int32_t ndat,
                     const double* x,
                     double* y,
                     int32_t block_sliced)
{

  // TODO use cuda streams to dispatch all the small matrix-matrix multiplications and run
  // concurrently on GPU

  if (block_sliced == 1) {

    for (int idat = 0; idat < ndat; ++idat) {

      int32_t shift = 0;

      for (int itypat=0; itypat < ntypat; ++itypat) {

        if (cplx == 2) {

          const cuDoubleComplex c_one  = make_cuDoubleComplex(1.0, 0.0);
          const cuDoubleComplex c_zero = make_cuDoubleComplex(0.0, 0.0);

          cuDoubleComplex *mat_ptr = (cuDoubleComplex *) mat;
          cuDoubleComplex* x_ptr = (cuDoubleComplex*) x;
          cuDoubleComplex* y_ptr = (cuDoubleComplex*) y;

          // move pointers
          mat_ptr += (lmnmax*lmnmax*itypat);
          x_ptr += shift;
          y_ptr += shift;

          cublasZhemm(cublas_handle, CUBLAS_SIDE_LEFT,CUBLAS_FILL_MODE_UPPER, nlmn[itypat], nattyp[itypat], &c_one,
                      mat_ptr, lmnmax,
                      x_ptr, nlmn[itypat], &c_zero,
                      y_ptr, nlmn[itypat]);

        } else {

          const double one  = 1.0;
          const double zero = 0.0;

          double *mat_ptr = (double *) mat;
          double* x_ptr = (double*) x;
          double* y_ptr = (double*) y;

          // move pointers
          mat_ptr += (lmnmax*lmnmax*itypat);
          x_ptr += shift;
          y_ptr += shift;

          cublasDsymm(cublas_handle, CUBLAS_SIDE_LEFT,CUBLAS_FILL_MODE_UPPER, nlmn[itypat], nattyp[itypat], &one,
                      mat_ptr, lmnmax,
                      x_ptr, nlmn[itypat], &zero,
                      y_ptr, nlmn[itypat]);

        } // end if cplx==2

        shift += nlmn[itypat] * nattyp[itypat];

      } // end for itypat

    } // end for idat

  } else { // block_sliced = 0

    // we regroupe matrix-matrix multiplications of same sizes into a strided batch of GEMM for idat=0:ndata-1
    // use regular matrix multipliy because batched hermitian matrices are not supported
    // see https://developer.nvidia.com/blog/cublas-strided-batched-matrix-multiply/

    if (cplx == 2) {

      const cuDoubleComplex c_one  = make_cuDoubleComplex(1.0, 0.0);
      const cuDoubleComplex c_zero = make_cuDoubleComplex(0.0, 0.0);

      cuDoubleComplex *mat_ptr = (cuDoubleComplex *) mat;
      cuDoubleComplex* x_ptr = (cuDoubleComplex*) x;
      cuDoubleComplex* y_ptr = (cuDoubleComplex*) y;

      for (int itypat=0; itypat < ntypat; ++itypat) {

        // move pointers to the next itypat beginning
        mat_ptr += (lmnmax*lmnmax);
        x_ptr   += nlmn[itypat]*nattyp[itypat];
        y_ptr   += nlmn[itypat]*nattyp[itypat];

        // C = alpha*A*B + beta*C
        // same matrix A for all idat => strideA = 0
        cublasZgemmStridedBatched(cublas_handle,
                                  CUBLAS_OP_N, // op A
                                  CUBLAS_OP_N, // op B
                                  nlmn[itypat], nattyp[itypat], nlmn[itypat], // m,n,k
                                  &c_one,                                     // alpha
                                  mat_ptr, lmnmax, 0,                         // A(m,k), lda, strideA
                                  x_ptr, nlmn[itypat], nprojs,                // B(k,n), ldb, strideB
                                  &c_zero,                                    // beta
                                  y_ptr, nlmn[itypat], nprojs,                // C(m,n), ldc, strideC
                                  ndat                                        // batch count
                                  );


      } // end of itypat

    } else { // cplx = 1

      const double one  = 1.0;
      const double zero = 0.0;

      double *mat_ptr = (double *) mat;
      double* x_ptr = (double*) x;
      double* y_ptr = (double*) y;

      for (int itypat=0; itypat < ntypat; ++itypat) {

        // move pointers to the next itypat beginning
        mat_ptr += (lmnmax*lmnmax);
        x_ptr   += nlmn[itypat]*nattyp[itypat];
        y_ptr   += nlmn[itypat]*nattyp[itypat];

        // C = alpha*A*B + beta*C
        // same matrix A for all idat => strideA = 0
        cublasDgemmStridedBatched(cublas_handle,
                                  CUBLAS_OP_N, // op A
                                  CUBLAS_OP_N, // op B
                                  nlmn[itypat], nattyp[itypat], nlmn[itypat], // m,n,k
                                  &one,                                       // alpha
                                  mat_ptr, lmnmax, 0,                         // A(m,k), lda, strideA
                                  x_ptr, nlmn[itypat], nprojs,                // B(k,n), ldb, strideB
                                  &zero,                                      // beta
                                  y_ptr, nlmn[itypat], nprojs,                // C(m,n), ldc, strideC
                                  ndat                                        // batch count
                                  );

      } // end for itypat

    } // end if cplx

  } // end if block_sliced

} // apply_block_gpu


//! \brief Allocation routine for apply inverse overlap operator.
//! we reallocate only if explicitely requested
//!
//! note:
//! proj_dim[0] = cplx
//! proj_dim[1] = nprojs
//! proj_dim[2] = ndat*nspinor
//!
//! \param[in] proj_dim dimension of proj array
//! \param[in] ntypat of nlmn array
//! \param[in] maximun number of lmn used as a maximum dimension for inverse overlap matrices
//! \param[in] realloc if 1, dealloc before reallocating
extern "C" void gpu_apply_invovl_inner_alloc(int32_t proj_dim[3],
                                             int32_t ntypat,
                                             int32_t realloc)
{
  // when memory allocation already done, and realloc is true, then
  // we deallocate first before re-allocating
  if (gpu_inner_allocated==1 and realloc == 1)
    gpu_apply_invovl_inner_dealloc();

  if (gpu_inner_allocated == 0) {
    //int32_t cplx = proj_dim[0];
    //int32_t nprojs = proj_dim[1];
    //int32_t ndat_nspinor = proj_dim[2];

    int32_t proj_size_in_bytes = proj_dim[0]*proj_dim[1]*proj_dim[2]*sizeof(double);

    CHECK_CUDA_ERROR( cudaMalloc((void**)&proj_gpu, proj_size_in_bytes) );
    CHECK_CUDA_ERROR( cudaMalloc((void**)&sm1proj_gpu, proj_size_in_bytes) );
    CHECK_CUDA_ERROR( cudaMalloc((void**)&ptp_sm1proj_gpu, proj_size_in_bytes) );

    CHECK_CUDA_ERROR( cudaMemset( proj_gpu, 0, proj_size_in_bytes) );
    CHECK_CUDA_ERROR( cudaMemset( sm1proj_gpu, 0, proj_size_in_bytes) );
    CHECK_CUDA_ERROR( cudaMemset( ptp_sm1proj_gpu, 0, proj_size_in_bytes) );

    CHECK_CUDA_ERROR( cudaMalloc((void**)&temp_proj_gpu, proj_size_in_bytes) );
    CHECK_CUDA_ERROR( cudaMalloc((void**)&resid_gpu, proj_size_in_bytes) );
    CHECK_CUDA_ERROR( cudaMalloc((void**)&precondresid_gpu, proj_size_in_bytes) );

    nlmn =   (uint8_t *)  malloc(ntypat * sizeof(uint8_t) );

    gpu_inner_allocated = 1;
  }

} // gpu_apply_invovl_inner_alloc

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_inner_dealloc()
{

  // if GPU data are not already allocated, don't do anything
  if (gpu_inner_allocated == 1) {

    CHECK_CUDA_ERROR( cudaFree(proj_gpu) );
    CHECK_CUDA_ERROR( cudaFree(sm1proj_gpu) );
    CHECK_CUDA_ERROR( cudaFree(ptp_sm1proj_gpu) );

    CHECK_CUDA_ERROR( cudaFree(temp_proj_gpu) );
    CHECK_CUDA_ERROR( cudaFree(resid_gpu) );
    CHECK_CUDA_ERROR( cudaFree(precondresid_gpu) );

    free(nlmn);

    gpu_inner_allocated = 0;
  }

} // gpu_apply_invovl_inner_dealloc

//! \brief Allocation routine for apply inverse overlap operator.
//! we reallocate only if explicitely requested
//!
//! \param[in] complex or real data
//! \param[in] maximun number of lmn used as a maximum dimension for inverse overlap matrices
//! \param[in] number of types of atoms
//! \param[in] realloc if 1, dealloc before reallocating
extern "C" void gpu_apply_invovl_matrix_alloc(int32_t cplx,
                                              int32_t nprojs,
                                              int32_t lmnmax,
                                              int32_t ntypat,
                                              int32_t realloc)
{
  // when memory allocation already done, and realloc is true, then
  // we deallocate first before re-allocating
  if (gpu_inv_overlap_matrix_allocated==1 and realloc == 1)
    gpu_apply_invovl_matrix_dealloc();

  if (gpu_inv_overlap_matrix_allocated == 0) {

    CHECK_CUDA_ERROR( cudaMalloc((void**) &gram_projs_gpu, cplx*nprojs*nprojs*sizeof(double)) );
    CHECK_CUDA_ERROR( cudaMalloc((void**) &inv_sij_gpu, cplx*lmnmax*lmnmax*ntypat*sizeof(double)) );
    CHECK_CUDA_ERROR( cudaMalloc((void**) &inv_s_approx_gpu, cplx*lmnmax*lmnmax*ntypat*sizeof(double)) );

    gpu_inv_overlap_matrix_allocated = 1;
  }

} // gpu_apply_invovl_matrix_alloc

//! \brief Allocation routine for apply inverse overlap operator.
extern "C" void gpu_apply_invovl_matrix_dealloc()
{

  // if GPU data are not already allocated, don't do anything
  if (gpu_inv_overlap_matrix_allocated == 1) {

    CHECK_CUDA_ERROR( cudaFree(gram_projs_gpu)   );
    CHECK_CUDA_ERROR( cudaFree(inv_sij_gpu)      );
    CHECK_CUDA_ERROR( cudaFree(inv_s_approx_gpu) );

    gpu_inv_overlap_matrix_allocated = 0;
  }

} // gpu_apply_invovl_matrix_dealloc


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

//! upload inverse overlap matrices
extern "C" void upload_inverse_overlap(const invovl_kpt_gpu_t invovl,
                                       int32_t cplx,
                                       int32_t nprojs,
                                       int32_t lmnmax,
                                       int32_t ntypat)
{

  CHECK_CUDA_ERROR( cudaMemcpy(gram_projs_gpu, invovl.gram_projs, cplx*nprojs*nprojs*sizeof(double), cudaMemcpyHostToDevice) );
  CHECK_CUDA_ERROR( cudaMemcpy(inv_sij_gpu, invovl.inv_sij, cplx*lmnmax*lmnmax*ntypat*sizeof(double), cudaMemcpyHostToDevice) );
  CHECK_CUDA_ERROR( cudaMemcpy(inv_s_approx_gpu, invovl.inv_s_approx, cplx*lmnmax*lmnmax*ntypat*sizeof(double), cudaMemcpyHostToDevice) );

} // upload_inverse_overlap

//! apply inverse overlap operator (inner part of the computation)
//!
//! note:
//! proj_dim[0] = cplx
//! proj_dim[1] = nprojs
//! proj_dim[2] = ndat*nspinor
//!
//! \param[in] proj_dim gives the dimension of proj, sm1proj and ptp_sm1proj
//! \param[in] proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[inout] sm1proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[inout] ptp_sm1proj is a 3D array of size (cplx, nprojs, nspinor*ndat)
//! \param[in] ndatspinor is the product of ndat and nspinor
//! \param[in] ntypat number of types of atoms
//! \param[in] lmnmax
//! \param[in] cplx (complex or real data)
//! \param[in] block_sliced (can only be 0 or 1), switch used inside apply_block_gpu
extern "C" void solve_inner_gpu(int32_t proj_dim[3],
                                const double* proj,
                                double* sm1proj,
                                double* ptp_sm1proj,
                                int32_t nattyp_dim,
                                int32_t* nattyp,
                                int32_t ntypat,
                                int32_t lmnmax,
                                int32_t cplx,
                                int32_t block_sliced)
{

  // projections norm array (for MPI all reduce)
  double* normprojs; // size is ndat*nspinor;

  // number of projections
  int32_t nprojs = proj_dim[1];

  // MPI error return value
  //int ierr;
  //int ibeg = 0;
  //int iend = nprojs;
  //int nlmntot_this_proc = nprojs;

  int ndat = proj_dim[2];
  int proj_size_in_bytes = proj_dim[0]*proj_dim[1]*proj_dim[2]*sizeof(double);

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

  // upload proj to GPU memory
  CHECK_CUDA_ERROR( cudaMemcpy(proj_gpu, proj, proj_size_in_bytes, cudaMemcpyHostToDevice) );

  // first guess for sm1proj
  // mat => inv_s_approx
  // x => proj
  // y => sm1proj
  // compute sm1proj = inv_s_approx * proj
  apply_block_gpu(cplx, ntypat, nattyp, lmnmax,
                  inv_s_approx_gpu, nprojs, ndat,
                  proj_gpu, sm1proj_gpu,
                  block_sliced);

  // Iterative refinement
  // TODO use a more efficient iterative algorithm than iterative refinement, use locking

  double maxerr, previous_maxerr;
  const double precision = 1e-16; // maximum relative error. TODO: use tolwfr ?
  double convergence_rate;
  int additional_steps_to_take = -1;

  //thrust::device_vector<double> d_temp_proj(temp_proj,temp_proj+(proj_dim[0]*proj_dim[1]*proj_dim[2]));

  const cuDoubleComplex c_one  = make_cuDoubleComplex(1.0, 0.0);
  const cuDoubleComplex c_zero = make_cuDoubleComplex(0.0, 0.0);

  for (int i=1; i < 30; ++i) {

    // compute resid = proj - (D^-1 + PtP)sm1proj
    apply_block_gpu(cplx, ntypat, nattyp, lmnmax,
                    inv_sij_gpu, nprojs, ndat,
                    sm1proj_gpu, resid_gpu,
                    block_sliced);

    // copy temp_proj = sm1proj(:,ibeg:iend,:);
    CHECK_CUDA_ERROR( cudaMemcpy(temp_proj_gpu, sm1proj_gpu, proj_size_in_bytes, cudaMemcpyDeviceToDevice));

    // compute matrix multiplication : PtPsm1proj(:,:,1) = invovl%gram * temp_proj(:,:,1)

    // TODO abi_xgemm('N', 'N', nprojs, ndat, nlmntot_this_proc, cone, &
    // & invovl%gram_projs(:,:,1), nprojs, &
    // & temp_proj(:,:,1), nlmntot_this_proc, czero, &
    // & PtPsm1proj(:,:,1), nprojs, &
    // & x_cplx=cplx)

    cublasZgemm(cublas_handle,
                CUBLAS_OP_N,            // op A
                CUBLAS_OP_N,            // op B
                nprojs, ndat, nprojs,   // m,n,k
                &c_one,                 // alpha
                (cuDoubleComplex*) gram_projs_gpu, nprojs, // A(m,k), lda
                (cuDoubleComplex*) temp_proj_gpu, nprojs,  // B(k,n), ldb
                &c_zero,                // beta
                (cuDoubleComplex*) ptp_sm1proj_gpu, nprojs // C(m,n), ldc
                );

    // TODO resid = proj - resid - Ptpsm1proj

    // exit check
    // TODO errs = SUM(SUM(resid**2, 1),1)

    // TODO maxerr = sqrt(MAXVAL(errs/normprojs))
    if (maxerr < precision or additional_steps_to_take == 1) {
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

  // download sm1proj and ptp_sm1proj to host memory
  CHECK_CUDA_ERROR( cudaMemcpy(sm1proj, sm1proj_gpu, proj_size_in_bytes, cudaMemcpyDeviceToHost) );
  CHECK_CUDA_ERROR( cudaMemcpy(ptp_sm1proj, ptp_sm1proj_gpu, proj_size_in_bytes, cudaMemcpyDeviceToHost) );


} // solve_inner_gpu
