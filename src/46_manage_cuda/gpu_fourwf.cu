//{\src2tex{textfont=tt}}
//****f* ABINIT/gpu_fourwf
//
// NAME
// gpu_fourwf
//
// FUNCTION
// Carry out composite Fourier transforms between real and reciprocal (G) space.
// Wavefunctions, contained in a sphere in reciprocal space,
// can be FFT to real space. They can also be FFT from real space
// to a sphere. Also, the density maybe accumulated, and a local
// potential can be applied.
//
// The different options are :
// - option=0 --> reciprocal to real space and output the result.
// - option=1 --> reciprocal to real space and accumulate the density.
// - option=2 --> reciprocal to real space, apply the local potential to the wavefunction
//                in real space and produce the result in reciprocal space.
// - option=3 --> real space to reciprocal space.
//
// This function calls different gpu functions to carry out ffts and operation in reciprocal space
//
// COPYRIGHT
// Copyright (C) 1998-2022 ABINIT group (FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
//
// INPUTS
// cplex= if 1 , denpot is real, if 2 , denpot is complex
//    (cplex=2 only allowed for option=2, and istwf_k=1)
//    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
// fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
//                 (intent(in) but the routine sphere can modify it for another iflag)
// gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
// gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
// istwf_k=option parameter that describes the storage of wfs
// kg_kin(3,npwin)=reduced planewave coordinates, input
// kg_kout(3,npwout)=reduced planewave coordinates, output
// mgfft=maximum size of 1D FFTs
// mpi_enreg=information about MPI parallelization
// ndat=number of FFT to do in //
// ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
// npwin=number of elements in fofgin array (for option 0, 1 and 2)
// npwout=number of elements in fofgout array (for option 2 and 3)
// n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
// option= if 0: do direct FFT
//         if 1: do direct FFT, then sum the density
//         if 2: do direct FFT, multiply by the potential, then do reverse FFT
//         if 3: do reverse FFT only
// paral_kgb=Flag related to the kpoint-band-fft parallelism
// tim_fourwf=timing code of the calling routine (can be set to 0 if not attributed)
// weight_r=weight to be used for the accumulation of the density in real space
//         (needed only when option=1)
// weight_i=weight to be used for the accumulation of the density in real space
//         (needed only when option=1 and (fftalg=4 and fftalgc/=0))
// fofginb(2,npwin)=holds second input wavefunction in G vector basis sphere.
//                 (intent(in) but the routine sphere can modify it for another iflag)
//                 (for non diagonal occupation)
// use_ndo = use non diagonal occupations.

// OUTPUT
//  (see side effects)
//
// SIDE EFFECTS
// Input/Output
// for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
//                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
//                no use of denpot, fofgout and npwout.
// for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
//                denpot(cplex*n4,n5,n6) contains the input density at input,
//                and the updated density at output (accumulated);
//                no use of fofgout and npwout.
// for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
//                denpot(cplex*n4,n5,n6) contains the input local potential;
//                fofgout(2,npwout*ndat) contains the output function;
// for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
//                fofgout(2,npwout*ndat) contains its output Fourier transform;
//                no use of fofgin and npwin.
//
// PARENTS
//      fourwf
//
// CHILDREN
//      gpu_sphere_in, gpu_sphere_out, gpu_apply_local_potential,\
//      FFTEXECC2C,gpu_density_accumulation,gpu_apply_local_potential
//
// SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "stdio.h"
#include "cuda_header.h"
#include "abi_gpu_header.h"
#include "cuda_api_error_check.h"

//STATIC vars to avoid too much cuda overhead
//Cuda environnement vars
static int gpu_initialized=0;
static cudaStream_t stream_cpy;
static cudaStream_t stream_compute;

//FFT plan
static int fft_size=-1;
static int ndat_loc=-1;
static int npw=-1;
static cufftHandle plan_fft;

//GPU ffts buffers
static double *work_gpu;
static double *fofr_gpu;
//static double *denpot_gpu;
static double *weightr_gpu;
static double *weighti_gpu;
//static double *fofgin_gpu;
//static double *fofgout_gpu;
//static int *kg_kin_gpu;
//static int *kg_kout_gpu;


//Transfers buffers in pinned memory for async memcopy between cpu & gpu
//static double *buff_denpot;
static double *buff_weightr;
static double *buff_weighti;
//static double *buff_kgkout;


extern "C" void gpu_fourwf_(int *cplex,
                            double *denpot,
                            double *fofgin,
                            double *fofgout,
                            double *fofr,
                            int *gboundin,
                            int *gboundout,
                            int *istwf_k,
                            int *kg_kin,
                            int *kg_kout,
                            int *mgfft,
                            void *mpi_enreg,
                            int *ndat,
                            int *ngfft,
                            int *npwin,
                            int *npwout,
                            int *n4,
                            int *n5,
                            int *n6,
                            int *option,
                            int *paral_kgb,
                            int * tim_fourwf,
                            double *weight_r,
                            double *weight_i)
/*, int *use_ndo,double *fofginb)*/
{
  // !Arguments ------------------------------------
  //   !scalars
  //   integer,intent(in) :: cplex,istwf_k,mgfft,n4,n5,n6,ndat,npwin,npwout,option,paral_kgb
  //   integer,intent(in) :: tim_fourwf
  //   integer,intent(in),optional :: use_ndo
  //   real(dp),intent(in) :: weight_r,weight_i
  //   type(MPI_type),intent(in) :: mpi_enreg
  //   !arrays
  //   integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
  //   integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
  //   real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin*ndat)
  //   ! real(dp) :: fofginb(2,npwin*ndat)
  //   real(dp),intent(inout),optional :: fofginb(:,:)
  //   real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
  //   real(dp),intent(out) :: fofgout(2,npwout*ndat)

  //Gpu buffers

  //Local integer
  int n1,n2,n3,nfft_tot;


  //*************** CUDA INITIALISATION STAGE ****
  if(!(gpu_initialized)){
    alloc_gpu_fourwf_(ngfft,ndat,npwin,npwout);
  }//end of initialisation


  //***********  GPU ALLOCATIONS  ***********************
  n1=ngfft[0];
  n2=ngfft[1];
  n3=ngfft[2];
  nfft_tot=n1*n2*n3;

  //*********** CHECK some compatibilities **************
  if((n1!=(*n4)) || (n2!=(*n5)) || (n3!=(*n6))){
    printf("FFT SIZE ERROR: \n when gpu mode is on the fft grid must not be augmented\n");
    printf("(n1,n2,n3) = (%d,%d,%d) whereas (n4,n5,n6) = (%d,%d,%d) \n",n1,n2,n3,*n4,*n5,*n6);
    fflush(stdout);
    abi_cabort();
  }

  if((*cplex)==2){
    printf("gpu_fourwf: cplex == 2 is not treated in GPU mode\n");
    fflush(stdout);
    abi_cabort();
  }

  //If fft size has changed, we realloc our buffers
  if((nfft_tot!=fft_size)||(*ndat>ndat_loc)||((*npwin)>npw)||((*npwout)>npw)){
    free_gpu_fourwf_();
    alloc_gpu_fourwf_(ngfft,ndat,npwin,npwout);
  }//end if "fft size changed"

  int deviceId;
  CHECK_CUDA_ERROR( cudaGetDevice(&deviceId) );
  if(*option!=3) {
    //printf("debug prefetch : &fofgin=%p fofgin=%f npwin=%d ndat=%d, deviceId=%d\n",fofgin, *fofgin,  *npwin, *ndat, deviceId);
    CHECK_CUDA_ERROR( cudaMemPrefetchAsync ( fofgin, 2*(*npwin)*(*ndat)*sizeof(double), deviceId) );
  }

  if(*option==2 || *option==3)
    CHECK_CUDA_ERROR( cudaMemPrefetchAsync ( fofgout, 2*(*npwout)*(*ndat)*sizeof(double), deviceId) );


  //memcpy cpu => buffer
  //if(*option == 1 || *option == 2)
  //  memcpy(buff_denpot,denpot,nfft_tot*sizeof(double));
  if(*option==1){
    memcpy(buff_weightr,weight_r,(*ndat)*sizeof(double));
    memcpy(buff_weighti,weight_i,(*ndat)*sizeof(double));
  }
  if (*option == 3){
    CHECK_CUDA_ERROR( cudaMemcpy(fofr_gpu,fofr,2*(*ndat)*nfft_tot*sizeof(double),cudaMemcpyHostToDevice) );
  }

  if(*option!=3){
    //CHECK_CUDA_ERROR( cudaMemcpy(kg_kin_gpu,kg_kin,3*(*npwin)*sizeof(int),cudaMemcpyHostToDevice) );
    //CHECK_CUDA_ERROR( cudaMemcpy(fofgin_gpu,fofgin,2*(*npwin)*(*ndat)*sizeof(double),cudaMemcpyHostToDevice) );

    //We launch async transfert of denpot
    // if(*option == 1 || *option == 2){
    //   CHECK_CUDA_ERROR( cudaMemcpyAsync(denpot_gpu,buff_denpot,nfft_tot*sizeof(double),cudaMemcpyHostToDevice,stream_cpy) );
    // }
    //We launch async transfert of denpot
    if(*option == 1){
      CHECK_CUDA_ERROR( cudaMemcpyAsync(weightr_gpu,buff_weightr,(*ndat)*sizeof(double),cudaMemcpyHostToDevice,stream_cpy) );
      CHECK_CUDA_ERROR( cudaMemcpyAsync(weighti_gpu,buff_weighti,(*ndat)*sizeof(double),cudaMemcpyHostToDevice,stream_cpy) );
    }

    //call preprocessing routine on gpu
    gpu_sphere_in_(fofgin,work_gpu,kg_kin,npwin,&n1,&n2,&n3,ndat,istwf_k,&stream_compute);

    //call backward fourrier transform on gpu work_gpu => fofr_gpu
    CHECK_CUDA_ERROR( FFTEXECC2C(plan_fft,(cucmplx *)work_gpu,(cucmplx *)fofr_gpu,CUFFT_INVERSE) );
  }

  if(*option==0){
    //We copy back fofr
    //CHECK_CUDA_ERROR( cudaStreamSynchronize(stream_compute) );
    CHECK_CUDA_ERROR( cudaDeviceSynchronize() );
    CHECK_CUDA_ERROR( cudaMemcpy(fofr,fofr_gpu,2*(*ndat)*nfft_tot*sizeof(double),cudaMemcpyDeviceToHost) );
  }

  if(*option==1){
    //We finish denpot and weight transferts
    CHECK_CUDA_ERROR( cudaStreamSynchronize(stream_cpy) );

    //call density accumulation routine on gpu
    gpu_density_accumulation_(fofr_gpu,denpot,weightr_gpu,weighti_gpu,&nfft_tot,ndat,&stream_compute);

    // when using managed memory, do a device sync before re-using data on host
    //CHECK_CUDA_ERROR( cudaDeviceSynchronize() );

    //We get denpot back on cpu
    //CHECK_CUDA_ERROR( cudaMemcpy(denpot,denpot_gpu,nfft_tot*sizeof(double),cudaMemcpyDeviceToHost) );
  }

  if(*option==2){
    //We finished denpot transfert
    cudaStreamSynchronize(stream_cpy);
    //call gpu routine to  Apply local potential
    gpu_apply_local_potential_((double2*)fofr_gpu,denpot,&nfft_tot,ndat,&stream_compute);

    // when using managed memory, do a device sync before re-using data on host
    CHECK_CUDA_ERROR( cudaDeviceSynchronize() );
  }

  if(*option==2 || *option==3){

    //call forward fourier transform on gpu: fofr_gpu ==> work_gpu
    CHECK_CUDA_ERROR( FFTEXECC2C(plan_fft,(cucmplx *)fofr_gpu,(cucmplx *)work_gpu,CUFFT_FORWARD) );

    //call post processing  routine on gpu
    //CHECK_CUDA_ERROR( cudaMemcpy(kg_kout_gpu,kg_kout,3*(*npwout)*sizeof(int),cudaMemcpyHostToDevice) );
    gpu_sphere_out_(fofgout,work_gpu,kg_kout,npwout,&n1,&n2,&n3,ndat,&stream_compute);

    // when using managed memory, do a device sync before re-using data on host
    //CHECK_CUDA_ERROR( cudaDeviceSynchronize() );

    //We get fofgout back on cpu
    //CHECK_CUDA_ERROR( cudaMemcpy(fofgout,fofgout_gpu,2*(*npwout)*(*ndat)*sizeof(double),cudaMemcpyDeviceToHost) );
  }

  // when using managed memory, do a device sync before re-using data on host
  CHECK_CUDA_ERROR( cudaDeviceSynchronize() );

}//end subroutine gpu_fourwf



//Memory allocation routine
extern "C" void alloc_gpu_fourwf_(int *ngfft, int *ndat, int *npwin, int *npwout)
{

  //printf("alloc_gpu_fourwf called with (nfft,ndat,npwin,npwout)=(%d,%d,%d,%d)\n",ngfft[0]*ngfft[1]*ngfft[2],*ndat,*npwin,*npwout);
  //fflush(stdout);

  gpu_initialized = 1;

  //Creation des streams cuda
  CHECK_CUDA_ERROR( cudaStreamCreate(&stream_cpy) );
  CHECK_CUDA_ERROR( cudaStreamCreate(&stream_compute) );

  //Size modification if too little memory allocation
  int n1=ngfft[0];
  int n2=ngfft[1];
  int n3=ngfft[2];
  fft_size=n1*n2*n3;
  if (ndat_loc < *ndat)
    ndat_loc=*ndat;
  if(npw < *npwin)
    npw = *npwin;
  if(npw < *npwout)
    npw = *npwout;


  //Initialisation des plans FFT
  int t_fft[3];
  t_fft[0]=n3;
  t_fft[1]=n2;
  t_fft[2]=n1;

  //Creation du plan
  CHECK_CUDA_ERROR( cufftPlanMany(&plan_fft,3,t_fft,NULL,1,0,NULL,1,0,FFT_C2C,*ndat) );

  //Association du plan au stream de calcul
  CHECK_CUDA_ERROR(cufftSetStream( plan_fft,stream_compute) );

  CHECK_CUDA_ERROR( cudaMalloc(&work_gpu,2*ndat_loc*fft_size*sizeof(double)) );
  CHECK_CUDA_ERROR( cudaMalloc(&fofr_gpu,2*ndat_loc*fft_size*sizeof(double)) );
  //CHECK_CUDA_ERROR( cudaMalloc(&denpot_gpu,fft_size*sizeof(double)) );
  CHECK_CUDA_ERROR( cudaMalloc(&weightr_gpu,ndat_loc*sizeof(double)) );
  CHECK_CUDA_ERROR( cudaMalloc(&weighti_gpu,ndat_loc*sizeof(double)) );
  //CHECK_CUDA_ERROR( cudaMalloc(&kg_kin_gpu,3*npw*sizeof(int)) );
  //CHECK_CUDA_ERROR( cudaMalloc(&fofgin_gpu,2*npw*ndat_loc*sizeof(double)) );
  //CHECK_CUDA_ERROR( cudaMalloc(&kg_kout_gpu,3*npw*sizeof(int)) );
  //CHECK_CUDA_ERROR( cudaMalloc(&fofgout_gpu,2*npw*ndat_loc*sizeof(double)) );

  // Allocation des buffers cpu "pinned"
  //CHECK_CUDA_ERROR( cudaMallocHost(&buff_denpot,fft_size*sizeof(double))) ;
  CHECK_CUDA_ERROR( cudaMallocHost(&buff_weightr,ndat_loc*sizeof(double)) );
  CHECK_CUDA_ERROR( cudaMallocHost(&buff_weighti,ndat_loc*sizeof(double)) );

}//End of alloc_gpu_fourwf_


extern "C" void free_gpu_fourwf_(){

  // destroy old fftw plan
  CHECK_CUDA_ERROR( cufftDestroy(plan_fft) );

  CHECK_CUDA_ERROR( cudaFree(work_gpu) );
  CHECK_CUDA_ERROR( cudaFree(fofr_gpu) );
  //CHECK_CUDA_ERROR( cudaFree(denpot_gpu) );
  CHECK_CUDA_ERROR( cudaFree(weightr_gpu) );
  CHECK_CUDA_ERROR( cudaFree(weighti_gpu) );
  //CHECK_CUDA_ERROR( cudaFree(kg_kin_gpu) );
  //CHECK_CUDA_ERROR( cudaFree(fofgin_gpu) );
  //CHECK_CUDA_ERROR( cudaFree(kg_kout_gpu) );
  //CHECK_CUDA_ERROR( cudaFree(fofgout_gpu) );
  //CHECK_CUDA_ERROR( cudaFreeHost(buff_denpot) );
  CHECK_CUDA_ERROR( cudaFreeHost(buff_weightr) );
  CHECK_CUDA_ERROR( cudaFreeHost(buff_weighti) );

  CHECK_CUDA_ERROR( cudaStreamDestroy(stream_cpy) );
  CHECK_CUDA_ERROR( cudaStreamDestroy(stream_compute) );

  gpu_initialized = 0;
}//end of free_gpu_fourwf_()
//***
