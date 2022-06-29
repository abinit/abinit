#ifndef GPU_GEMM_NONLOP_KERNEL_H_
#define GPU_GEMM_NONLOP_KERNEL_H_

//! cuda kernel to extract real or imaginary part of a complex array
//!
//! \param[out] data_out   is an array of size   data_size    (real)
//! \param[in]  data_in    is an array of size 2*data_size    (cplx)
//! \param[in]  data_size  is is the size of the output array (real)
//!
//! \tparam offset is an integer (can only be 0 or 1)
//! 0 means we extract the real part
//! 1 means we extract the imag part
//!
//! \note see kernel_insert_real_imag (reverse operation)
template<int offset>
__global__
void kernel_extract_real_imag(double* data_out, const double* data_in, int32_t data_size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < data_size) {

    data_out[i] = data_in[2*i+offset];

  }

} // kernel_extract_real_imag

//! cuda kernel to insert real or imaginary part in a complex array
//!
//! \param[out] data_out   is an array of size 2*data_size   (cplx)
//! \param[in]  data_in    is an array of size   data_size   (real)
//! \param[in]  data_size  is is the size of the input array (real)
//!
//! \tparam offset is an integer (can only be 0 or 1)
//! 0 means we insert the real part
//! 1 means we insert the imag part
//!
//! \note see kernel_extract_real_imag (reverse operation)
template<int offset>
__global__
void kernel_insert_real_imag(double* data_out, const double* data_in, int32_t data_size)
{

  int32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < data_size) {

    data_out[2*i+offset] = data_in[i];

  }

} // kernel_insert_real_imag

//! cuda kernel for specific use in gemm nonlop
//!
//! \param[in,out] data   is an array of size data_size
//! \param[in]     npw_in is the number of plane waves
//! \param[in]     ndat_nspinor is the product of ndat and nspinor
//! \param[in]     option: 0 means divide by 2, 1 means zero-out
//!
__global__
void kernel_fix_realvec(double* data, int32_t npw_in, int32_t ndat_nspinor, int32_t option)
{

  int32_t idat = blockIdx.x * blockDim.x + threadIdx.x;

  if (idat < ndat_nspinor)
  {
    int32_t index = idat * npw_in;
    data[index] = option == 0 ? data[index]/2 : 0;
  }

} // kernel_fix_realvec

#endif // GPU_GEMM_NONLOP_KERNEL_H_
