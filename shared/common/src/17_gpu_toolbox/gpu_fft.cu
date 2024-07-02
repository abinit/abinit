/* gpu_fft.cu */

/*
 * Copyright (C) 2008-2024 ABINIT Group
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt .
 * for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
 *
 * The main goal of this file is to contain cublas and magma encapsulation routines,
 * that will be callable from fortran routines
 *
 */

#include<assert.h>
#include<gpu_fft.h>

cufftHandle plan_fft;
static cudaStream_t stream_compute;

//! utility function to select eigen type
static cufftType select_cufft_type(const int fftType_int)
{
  switch(fftType_int){
    case(CUFFT_R2C): return CUFFT_R2C;
    case(CUFFT_C2R): return CUFFT_C2R;
    case(CUFFT_C2C): return CUFFT_C2C;
    case(CUFFT_D2Z): return CUFFT_D2Z;
    case(CUFFT_Z2D): return CUFFT_Z2D;
    case(CUFFT_Z2Z): return CUFFT_Z2Z;
    default:
      fprintf(stderr, "Provided wrong enum value for FFT type: %d\n",fftType_int);
      fflush(stderr);
      abi_cabort();
      return CUFFT_R2C;
  }
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_plan_many
 *
 * FUNCTION
 *  Initialize a FFT plan with custom dimension, strided and batch size.
 *
 * INPUTS
 *   rank      Dimensionality of the transform (1, 2, or 3).
 *   n         Array of size rank, describing the size of each dimension,
 *             n[0] being the size of the outermost and n[rank-1] innermost
 *             (contiguous) dimension of a transform.
 *   inembed   Pointer of size rank that indicates the storage dimensions
 *             of the input data in memory.
 *             If set to NULL all other advanced data layout parameters are ignored.
 *   istride   Indicates the distance between two successive input elements
 *             in the least significant (i.e., innermost) dimension
 *   idist     Indicates the distance between the first element of two
 *             consecutive signals in a batch of the input data
 *   onembed   Pointer of size rank that indicates the storage dimensions of
 *             the output data in memory.
 *             If set to NULL all other advanced data layout parameters are ignored.
 *   ostride   Indicates the distance between two successive output elements in
 *             the output array in the least significant (i.e., innermost) dimension
 *   odist     Indicates the distance between the first element of two
 *             consecutive signals in a batch of the output data
 *   type      The transform data type
 *             (e.g., FFT_R2C for single precision real to complex)
 *   batch     Batch size for this transform
 */
/*=========================================================================*/

extern "C"
void gpu_fft_plan_many_cpp(int *rank, int **n, int **inembed,
                       int *istride, int *idist, int **onembed, int *ostride,
                       int *odist, int *fft_type, int *batch){

  assert(CUFFT_Z2Z==0x69 && "cuFFT_Type enum value mismatch !(CUDA update?)");
  assert(CUFFT_FORWARD==-1 && "cuFFT direction enum value mismatch (CUDA update?)");
  assert(CUFFT_INVERSE== 1 && "cuFFT direction enum value mismatch (CUDA update?)");

  cufftType type = select_cufft_type(*fft_type);
  CUDA_API_CHECK(cufftPlanMany(
        &plan_fft,
        *rank,
        *n,
        *inembed,
        *istride,
        *idist,
        *onembed,
        *ostride,
        *odist,
        type,
        *batch));
  CUDA_API_CHECK( cudaStreamCreate(&stream_compute) );
  CUDA_API_CHECK( cufftSetStream(plan_fft,stream_compute) );
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_stream_synchronize
 *
 * FUNCTION
 *  Wait for any FFT operations still running on stream
 */
/*=========================================================================*/

extern "C" void gpu_fft_stream_synchronize_cpp()
{
  CUDA_API_CHECK( cudaStreamSynchronize(stream_compute) );
}


/*=========================================================================*/
// NAME
//  gpu_fft_plan_destroy
//
// FUNCTION
//  Destroy FFT plan
//
/*=========================================================================*/

extern "C"
void gpu_fft_plan_destroy_cpp(void){
  CUDA_API_CHECK(cufftDestroy(plan_fft));
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_exec_z2z
 *
 * FUNCTION
 *  Run a Fast Fourrier Transform on double-complex input and output
 *
 * INPUTS
 *   idata       Pointer to the complex input data (in GPU memory) to transform
 *   odata       Pointer to the complex output data (in GPU memory)
 *   direction   The transform direction: FFT_FORWARD or FFT_INVERSE
 *
 * OUTPUT
 *   odata       Contains the complex Fourier coefficients
 */
/*=========================================================================*/

extern "C"
void gpu_fft_exec_z2z_cpp(void **idata, void **odata, int *direction){

  CUDA_API_CHECK(cufftExecZ2Z(plan_fft, (cufftDoubleComplex*) (*idata),
                 (cufftDoubleComplex*) (*odata), *direction));
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_exec_c2c
 *
 * FUNCTION
 *  Run a Fast Fourrier Transform on float complex input and output
 *
 * INPUTS
 *   idata       Pointer to the complex input data (in GPU memory) to transform
 *   odata       Pointer to the complex output data (in GPU memory)
 *   direction   The transform direction: FFT_FORWARD or FFT_INVERSE
 *
 * OUTPUT
 *   odata       Contains the complex Fourier coefficients
 */
/*=========================================================================*/

extern "C"
void gpu_fft_exec_c2c_cpp(void **idata, void **odata, int *direction){
  CUDA_API_CHECK(cufftExecC2C(plan_fft, (cufftComplex*) *idata,
                 (cufftComplex*) *odata, *direction));
}

