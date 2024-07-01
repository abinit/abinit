#ifndef GPU_APPLY_INVOVL_INNER_KERNEL_H_
#define GPU_APPLY_INVOVL_INNER_KERNEL_H_

#define _CG_ABI_EXPERIMENTAL
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template <class T>
struct SharedMemory {
  __device__ inline operator T *() {
    extern __shared__ int __smem[];
    return (T *)__smem;
  }

  __device__ inline operator const T *() const {
    extern __shared__ int __smem[];
    return (T *)__smem;
  }
};

// specialize for double to avoid unaligned memory
// access compile errors
template <>
struct SharedMemory<double> {
  __device__ inline operator double *() {
    extern __shared__ double __smem_d[];
    return (double *)__smem_d;
  }

  __device__ inline operator const double *() const {
    extern __shared__ double __smem_d[];
    return (double *)__smem_d;
  }
};

//! cuda kernel to compute residue : resid = proj - resid - Ptpsm1proj
//! cuda grid size is fixed.
//!
//! \param[in]  proj       is the projector  matrix (cplx,nprojs,ndat)
//! \param[out] resid      is the residue    matrix (cplx,nprojs,ndat)
//! \param[in] ptp_sm1proj is the ptpsm1proj matrix (cplx,nprojs,ndat)
//! \param[in] proj_size is the total size (= cplx*nprojs*ndat)
__global__
void compute_residue_v1(const double* proj, double* resid, const double* ptp_sm1proj, int32_t proj_size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < proj_size) {

    resid[i] = proj[i] - resid[i] - ptp_sm1proj[i];

    i += gridDim.x * blockDim.x;

  }

} // compute_residue_v1

//! cuda kernel to compute residue : resid = proj - resid - Ptpsm1proj
//! cuda grid size is adapted to data size
//!
//! \param[in]  proj       is the projector  matrix (cplx,nprojs,ndat)
//! \param[out] resid      is the residue    matrix (cplx,nprojs,ndat)
//! \param[in] ptp_sm1proj is the ptpsm1proj matrix (cplx,nprojs,ndat)
//! \param[in] proj_size is the total size (= cplx*nprojs*ndat)
__global__
void compute_residue_v2(const double* proj, double* resid, const double* ptp_sm1proj, int32_t proj_size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < proj_size) {

    resid[i] = proj[i] - resid[i] - ptp_sm1proj[i];

  }

} // compute_residue_v2

//! cuda kernel to compute sm1proj : sm1proj = sm1proj + precondresid
//! cuda grid size is fixed.
//!
//! \param[inout]  sm1proj       is the projector  matrix (cplx,nprojs,ndat)
//! \param[in] precondresid
__global__
void compute_sm1proj_v1(double* sm1proj, const double* precondresid, int32_t size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < size) {

    sm1proj[i] += precondresid[i];

    i += gridDim.x * blockDim.x;

  }

} // compute_sm1proj_v1

//! cuda kernel to compute sm1proj : sm1proj = sm1proj + precondresid
//! cuda grid size is fixed.
//!
//! \param[inout]  sm1proj       is the projector  matrix (cplx,nprojs,ndat)
//! \param[in] precondresid
__global__
void compute_sm1proj_v2(double* sm1proj, const double* precondresid, int32_t size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < size) {

    sm1proj[i] += precondresid[i];

  }

} // compute_sm1proj_v2

//! cuda kernel to compute -sm1proj and - ptp_sm1proj
//! cuda grid size is fixed.
//!
//! \param[inout] sm1proj
//! \param[inout] ptp_sm1proj
__global__
void compute_sm1proj_negate(double* sm1proj, double* ptp_sm1proj, int32_t size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < size) {

    sm1proj[i]     = - sm1proj[i];
    ptp_sm1proj[i] = - ptp_sm1proj[i];

  }

} // compute_sm1proj_v2

//
// custom reduction for computing sums of square residues
//

//! functor square computes the square of a number, i.e. f(x) -> x*x
// template <typename T>
// struct square
// {
//   __host__ __device__
//   T operator()(const T& x) const {
//     return x * x;
//   }
// };

template<typename T>
__device__ T square(const T& x)
{
  return x*x;
}

//! compute sum of squares of errors of residues
//! one value per idat.
//!
//! The following code is adapted from nvidia cuda sample, reduce6 kernel, except that:
//! - we use exactly one block per value of idat, i.e. ndat blocks of threads
//! - each block reads data along the proj dimension
//!
//! \param[in] resid  is the residue    matrix (cplx,nprojs,ndat)
//! \param[out] errs  is the residue error vector (size is ndat)
//! \param[in]  cplx
//! \param[in]  nprojs
//! \param[in]  ndat
//!
template <class T, unsigned int blockSize>
__global__ void compute_residue_errs(const T* resid, T* errs, uint32_t cplx, uint32_t nprojs, uint32_t ndat)
{
  // Handle to thread block group
  cg::thread_block cta = cg::this_thread_block();
  T *sdata = SharedMemory<T>();

  T mySum = 0;

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const unsigned int tid = threadIdx.x;
  unsigned int tid2 = tid;

  const unsigned int idat = blockIdx.x;
  const unsigned int idat_offset = cplx*nprojs*idat;

  //if (threadIdx.x == 0) printf("KKKK %d %d %d\n",cplx,nprojs,ndat);

  while (tid2 < cplx*nprojs) {
    //if (idat==0) printf("tid2 = %d, index=%d\n",tid2,idat_offset+tid2);
    mySum += square(resid[idat_offset+tid2]);
    tid2 += blockDim.x;
  }

  // each thread puts its local sum into shared memory
  sdata[tid] = mySum;
  cg::sync(cta);

  // do reduction in shared mem
  if ((blockSize >= 512) && (tid < 256)) {
    sdata[tid] = mySum = mySum + sdata[tid + 256];
  }

  cg::sync(cta);

  if ((blockSize >= 256) && (tid < 128)) {
    sdata[tid] = mySum = mySum + sdata[tid + 128];
  }

  cg::sync(cta);

  if ((blockSize >= 128) && (tid < 64)) {
    sdata[tid] = mySum = mySum + sdata[tid + 64];
  }

  cg::sync(cta);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);

  if (cta.thread_rank() < 32) {
    // Fetch final intermediate sum from 2nd warp
    if (blockSize >= 64) mySum += sdata[tid + 32];
    // Reduce final warp using shuffle
    for (int offset = tile32.size() / 2; offset > 0; offset /= 2) {
      mySum += tile32.shfl_down(mySum, offset);
    }
  }

  // write result for this block to global mem
  if (cta.thread_rank() == 0) errs[idat] = mySum;
}


#endif // GPU_APPLY_INVOVL_INNER_KERNEL_H_
