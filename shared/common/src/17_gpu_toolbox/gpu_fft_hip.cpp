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

hipfftHandle plan_fft;
static hipStream_t stream_compute;

//! utility function to select eigen type
static hipfftType select_hipfft_type(const int fftType_int)
{
  switch(fftType_int){
    case(HIPFFT_R2C): return HIPFFT_R2C;
    case(HIPFFT_C2R): return HIPFFT_C2R;
    case(HIPFFT_C2C): return HIPFFT_C2C;
    case(HIPFFT_D2Z): return HIPFFT_D2Z;
    case(HIPFFT_Z2D): return HIPFFT_Z2D;
    case(HIPFFT_Z2Z): return HIPFFT_Z2Z;
    default:
      fprintf(stderr, "Provided wrong enum value for FFT type: %d\n",fftType_int);
      fflush(stderr);
      abi_cabort();
      return HIPFFT_R2C;
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

  assert(HIPFFT_Z2Z==0x69 && "hipFFT_Type enum value mismatch !(HIP update?)");
  assert(HIPFFT_FORWARD==-1 && "hipFFT direction enum value mismatch (HIP update?)");
  assert(HIPFFT_BACKWARD == 1 && "hipFFT direction enum value mismatch (HIP update?)");

  hipfftType type = select_hipfft_type(*fft_type);
  HIP_API_CHECK(hipfftPlanMany(
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
  HIP_API_CHECK( hipStreamCreate(&stream_compute) );
  HIP_API_CHECK( hipfftSetStream(plan_fft,stream_compute) );
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
  HIP_API_CHECK( hipStreamSynchronize(stream_compute) );
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
  HIP_API_CHECK(hipfftDestroy(plan_fft));
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

  HIP_API_CHECK(hipfftExecZ2Z(plan_fft, (hipfftDoubleComplex*) (*idata),
                 (hipfftDoubleComplex*) (*odata), *direction));
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
  HIP_API_CHECK(hipfftExecC2C(plan_fft, (hipfftComplex*) *idata,
                 (hipfftComplex*) *odata, *direction));
}

