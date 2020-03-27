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
// Copyright (C) 1998-2020 ABINIT group (FDahm)
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
static double *denpot_gpu;
static double *weightr_gpu;
static double *weighti_gpu;
static double *fofgin_gpu;
static double *fofgout_gpu;
static int *kg_kin_gpu;
static int *kg_kout_gpu;


//Transfers buffers in pinned memory for async memcopy between cpu & gpu
static double *buff_denpot;
static double *buff_weightr;
static double *buff_weighti;
//static double *buff_kgkout;


extern "C" void gpu_fourwf_(int *cplex,double *denpot,double *fofgin,double *fofgout,double *fofr,int *gboundin,int *gboundout,int *istwf_k,
			    int *kg_kin,int *kg_kout,int *mgfft,void *mpi_enreg,int *ndat,int *ngfft,int *npwin,int *npwout,int *n4,int *n5,int *n6,int *option,
			    int *paral_kgb,int * tim_fourwf,double *weight_r,double *weight_i)/*,
												int *use_ndo,double *fofginb)*/
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




  //Cuda return code
  cufftResult cufft_state;
  cudaError_t cuda_return;

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


  //memcpy cpu => buffer
  if(*option == 1 || *option == 2)
    memcpy(buff_denpot,denpot,nfft_tot*sizeof(double));
  if(*option==1){
    memcpy(buff_weightr,weight_r,(*ndat)*sizeof(double));
    memcpy(buff_weighti,weight_i,(*ndat)*sizeof(double));
  }
  if (*option == 3){
    cuda_return = cudaMemcpy(fofr_gpu,fofr,2*(*ndat)*nfft_tot*sizeof(double),cudaMemcpyHostToDevice);
    if(cuda_return != cudaSuccess){
      printf("ERROR while copying fofr to gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }
  }

  if(*option!=3){
    cuda_return = cudaMemcpy(kg_kin_gpu,kg_kin,3*(*npwin)*sizeof(int),cudaMemcpyHostToDevice);
    cuda_return = cudaMemcpy(fofgin_gpu,fofgin,2*(*npwin)*(*ndat)*sizeof(double),cudaMemcpyHostToDevice);
    if(cuda_return != cudaSuccess){
      printf("ERROR while copying input data to gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }

    //We launch async transfert of denpot
    if(*option == 1 || *option == 2){
      cuda_return = cudaMemcpyAsync(denpot_gpu,buff_denpot,nfft_tot*sizeof(double),cudaMemcpyHostToDevice,stream_cpy);
      if(cuda_return != cudaSuccess){
	printf("ERROR while copying denpot to gpu: %s \n",cudaGetErrorString(cuda_return));
	fflush(stdout);
	abi_cabort();
      }
    }
    //We launch async transfert of denpot
    if(*option == 1){
      cuda_return = cudaMemcpyAsync(weightr_gpu,buff_weightr,(*ndat)*sizeof(double),cudaMemcpyHostToDevice,stream_cpy);
      cuda_return = cudaMemcpyAsync(weighti_gpu,buff_weighti,(*ndat)*sizeof(double),cudaMemcpyHostToDevice,stream_cpy);
      if(cuda_return != cudaSuccess){
	printf("ERROR while copying weight to gpu: %s \n",cudaGetErrorString(cuda_return));
	fflush(stdout);
	abi_cabort();
      }
    }

    //call preprocessing routine on gpu
    gpu_sphere_in_(fofgin_gpu,work_gpu,kg_kin_gpu,npwin,&n1,&n2,&n3,ndat,istwf_k,&stream_compute);

    //call backward fourrier transform on gpu work_gpu => fofr_gpu
    cufft_state = FFTEXECC2C(plan_fft,(cucmplx *)work_gpu,(cucmplx *)fofr_gpu,CUFFT_INVERSE);
    if(cufft_state!=CUFFT_SUCCESS){
      printf("ERROR while fft work_gpu ==> fofr:\n%s\n",cufftGetErrorString(cufft_state));
      //printf("last cudaError was : %s \n",cudaGetErrorString(cudaGetLastError()));
      fflush(stdout);
      abi_cabort();
    }
  }

  if(*option==0){
    //We copy back fofr
    //cuda_return = cudaStreamSynchronize(stream_compute);
    cuda_return = cudaThreadSynchronize();
    if(cuda_return != cudaSuccess){
      printf("ERROR when synchronizing after FFT on gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }
    cuda_return = cudaMemcpy(fofr,fofr_gpu,2*(*ndat)*nfft_tot*sizeof(double),cudaMemcpyDeviceToHost);
    if(cuda_return != cudaSuccess){
      printf("ERROR while copying fofr from gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }
  }

  if(*option==1){
    //We finish denpot and weight transferts
    cuda_return = cudaStreamSynchronize(stream_cpy);
    if(cuda_return != cudaSuccess){
      printf("ERROR while getting denpot and weight on gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }

    //call density accumulation routine on gpu
    gpu_density_accumulation_(fofr_gpu,denpot_gpu,weightr_gpu,weighti_gpu,&nfft_tot,ndat,&stream_compute);

    //We get denpot back on cpu
    cuda_return = cudaMemcpy(denpot,denpot_gpu,nfft_tot*sizeof(double),cudaMemcpyDeviceToHost);
    if(cuda_return != cudaSuccess){
      printf("ERROR while copying denpot from gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }
  }

  if(*option==2){
    //We finished denpot transfert
    cudaStreamSynchronize(stream_cpy);
    //call gpu routine to  Apply local potential
    gpu_apply_local_potential_((double2*)fofr_gpu,denpot_gpu,&nfft_tot,ndat,&stream_compute);
  }

  if(*option==2 || *option==3){

    //call forward fourier transform on gpu: fofr_gpu ==> work_gpu
    cufft_state = FFTEXECC2C(plan_fft,(cucmplx *)fofr_gpu,(cucmplx *)work_gpu,CUFFT_FORWARD);
    if(cufft_state!=CUFFT_SUCCESS){
      printf("ERROR while fft fofr ==> work_gpu:\n%s\n",cufftGetErrorString(cufft_state));
      fflush(stdout);
      abi_cabort();
    }

    //call post processing  routine on gpu
    cuda_return = cudaMemcpy(kg_kout_gpu,kg_kout,3*(*npwout)*sizeof(int),cudaMemcpyHostToDevice);
    gpu_sphere_out_(fofgout_gpu,work_gpu,kg_kout_gpu,npwout,&n1,&n2,&n3,ndat,&stream_compute);

    //We get fofgout back on cpu
    cuda_return = cudaMemcpy(fofgout,fofgout_gpu,2*(*npwout)*(*ndat)*sizeof(double),cudaMemcpyDeviceToHost);
    if(cuda_return != cudaSuccess){
      printf("ERROR while copying fofgout from gpu: %s \n",cudaGetErrorString(cuda_return));
      fflush(stdout);
      abi_cabort();
    }
  }

}//end subroutine gpu_fourwf



//Memory allocation routine
extern "C" void alloc_gpu_fourwf_(int *ngfft,int *ndat,int *npwin,int *npwout){

//   printf("alloc_gpu_fourwf called with (nfft,ndat,npwin,npwout)=(%d,%d,%d,%d)\n",ngfft[0]*ngfft[1]*ngfft[2],*ndat,*npwin,*npwout);
//   fflush(stdout);

  cufftResult cufft_state;
  cudaError_t cuda_return;

  gpu_initialized = 1;

  //Creation des streams cuda
  cuda_return = cudaStreamCreate(&stream_cpy);
  cuda_return = cudaStreamCreate(&stream_compute);
  if(cuda_return != cudaSuccess){
    printf("ERROR: when creating cuda streams:\n%s\n",cudaGetErrorString(cuda_return));
    fflush(stdout);
    abi_cabort();
  }

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
  cufft_state=cufftPlanMany(&plan_fft,3,t_fft,NULL,1,0,NULL,1,0,FFT_C2C,*ndat);
  if(cufft_state!=CUFFT_SUCCESS){
    printf("alloc_gpu_fourwf:\n ERROR: At creation of cufftPlan:\n%s\n",cufftGetErrorString(cufft_state));
    fflush(stdout);
    abi_cabort();
  }

  //Association du plan au stream de calcul
  cufft_state = cufftSetStream(plan_fft,stream_compute);
  if(cufft_state!=CUFFT_SUCCESS){
    printf("alloc_gpu_fourwf:\n ERROR: while associating cufftPlan to a stream:\n%s\n",cufftGetErrorString(cufft_state));
    fflush(stdout);
    abi_cabort();
  }

  cuda_return = cudaMalloc(&work_gpu,2*ndat_loc*fft_size*sizeof(double));
  cuda_return = cudaMalloc(&fofr_gpu,2*ndat_loc*fft_size*sizeof(double));
  cuda_return = cudaMalloc(&denpot_gpu,fft_size*sizeof(double));
  cuda_return = cudaMalloc(&weightr_gpu,ndat_loc*sizeof(double));
  cuda_return = cudaMalloc(&weighti_gpu,ndat_loc*sizeof(double));
  cuda_return = cudaMalloc(&kg_kin_gpu,3*npw*sizeof(int));
  cuda_return = cudaMalloc(&fofgin_gpu,2*npw*ndat_loc*sizeof(double));
  cuda_return = cudaMalloc(&kg_kout_gpu,3*npw*sizeof(int));
  cuda_return = cudaMalloc(&fofgout_gpu,2*npw*ndat_loc*sizeof(double));
  if(cuda_return != cudaSuccess){
    printf("alloc_gpu_fourwf: ERROR while allocating memory on gpu: %s \n",cudaGetErrorString(cuda_return));
    fflush(stdout);
    abi_cabort();
  }


  //Allocation des buffers cpu "pinned"
  cuda_return = cudaMallocHost(&buff_denpot,fft_size*sizeof(double));
  cuda_return = cudaMallocHost(&buff_weightr,ndat_loc*sizeof(double));
  cuda_return = cudaMallocHost(&buff_weighti,ndat_loc*sizeof(double));
  if(cuda_return != cudaSuccess){
    printf("alloc_gpu_fourwf:\n ERROR while allocating pinned memory for transfert to gpu: %s \n",cudaGetErrorString(cuda_return));
    fflush(stdout);
    abi_cabort();
  }
}//End of alloc_gpu_fourwf_


extern "C" void free_gpu_fourwf_(){
  cufftResult cufft_state;
  cudaError_t cuda_return;

  //On detruit l'ancien plan
  cufft_state = cufftDestroy(plan_fft);
  if(cufft_state!=CUFFT_SUCCESS){
    printf("free_gpu_fourwf:\n ERROR: at destruction of cufftPlan \n");
    fflush(stdout);
    abi_cabort();
  }

  cuda_return = cudaFree(work_gpu);
  cuda_return = cudaFree(fofr_gpu);
  cuda_return = cudaFree(denpot_gpu);
  cuda_return = cudaFree(weightr_gpu);
  cuda_return = cudaFree(weighti_gpu);
  cuda_return = cudaFree(kg_kin_gpu);
  cuda_return = cudaFree(fofgin_gpu);
  cuda_return = cudaFree(kg_kout_gpu);
  cuda_return = cudaFree(fofgout_gpu);
  cuda_return = cudaFreeHost(buff_denpot);
  cuda_return = cudaFreeHost(buff_weightr);
  cuda_return = cudaFreeHost(buff_weighti);
  if(cuda_return != cudaSuccess){
    printf("free_gpu_fourwf:\n ERROR while freeing memory on gpu: %s \n",cudaGetErrorString(cuda_return));
    fflush(stdout);
    abi_cabort();
  }

  cuda_return = cudaStreamDestroy(stream_cpy);
  cuda_return = cudaStreamDestroy(stream_compute);
  if(cuda_return != cudaSuccess){
    printf("free_gpu_fourwf:\n ERROR: at destruction of streams: %s \n",cudaGetErrorString(cuda_return));
    fflush(stdout);
    abi_cabort();
  }

  gpu_initialized = 0;
}//end of free_gpu_fourwf_()
//***


