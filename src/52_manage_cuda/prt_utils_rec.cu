/* prt_utils_rec.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "cuda_common.h"
#include "cuda_header.h"
#include "cuda_rec_head.h"

/*=========================================================================*/
/*________________________ GPU_function called by HOST_____________________*/
/*=========================================================================*/

__host__ void prt_dbg_arr(cureal* v_d,
			  size_t size,int num,int pos,
			  char bo[20])
{
#ifdef  HAVE_GPU_CUDA_DEBUG
  cudaDeviceSynchronize();
  cureal* v_h  = NULL;
  v_h = (cureal*) malloc(size);
  cudaMemcpy(v_h,&(v_d[pos]),size,cudaMemcpyDeviceToHost);
  printf("%s ",bo);
  for(int jj=0;jj<num; jj++) printf("%9.5e ",v_h[jj]);
  //for(int jj=0;jj<num; jj++) printf("%9.5f    \n",v_h[jj]);
  printf("\n");
  free(v_h);
#endif
  return;
}
__host__ void prt_dbg_arrc(cucmplx* v_d,
			  size_t size,int num,int pos,
			  char bo[20])
{
#ifdef  HAVE_GPU_CUDA_DEBUG
  cudaDeviceSynchronize();
  cucmplx* v_h = NULL;
  v_h = (cucmplx*) malloc(size);
  cudaMemcpy(v_h,&(v_d[pos]),size,cudaMemcpyDeviceToHost);
  printf("%s ",bo);
  for(int jj=0;jj<num; jj++) printf("%9.5e ",v_h[jj].x);
  printf("\n");
  free(v_h);
#endif 
  return;
}

__host__ void  prt_mem_use(size_t un_pitch,size_t hauteur,size_t size)
{
  size_t largeur = sizeof(cureal)*size;
  size_t clargeur = sizeof(cufftComplex)*size;
  size_t totmem = 3*un_pitch*hauteur+2*clargeur+largeur+(sizeof(cureal)+sizeof(int))*hauteur;
  printf( " ___________________________________________________________________\n");
  printf( " ____________  Allocated Memory on Device for Recursion ____________\n");
  printf( "    Pitched Points            %10d       \n",int( un_pitch/sizeof(cureal)));
  printf( "    Number of Vectors         %10d       \n",(int)(hauteur));
  printf( "    Sizes Real Vectors        %10d  (bytes)\n",(int)largeur);
  printf( "    Sizes Complex Vectors     %10d  (bytes)\n",(int)clargeur);
  printf( "    Size Matrix of vectors    %10d  (bytes)\n",(int)(un_pitch*hauteur));
  printf( "    Allocated memory on GPU   %10.2f (Mbytes)\n",(float)totmem/1048576.);
  printf( " ___________________________________________________________________\n\n");
}



__host__ void starttime(cudaEvent_t* start)
{
#ifdef  HAVE_GPU_CUDA_TM
  cudaEventRecord(*start,0);
#endif
  return;
}

__host__ void calctime(cudaEvent_t* stop,cudaEvent_t start,
		       float* timer,int index)
{
#ifdef  HAVE_GPU_CUDA_TM
  if(index>=DEBUGLEN){
    printf("ERROR!! YOU ARE TESTING TIME FOR TOO MUCH KERNELS\n");
    printf("TRY TO INCREASE IT IN \"cuda_rec_head.cu\"\n");
    exit(EXIT_FAILURE);;
  }
  cudaDeviceSynchronize();
  float bo;
  cudaEventRecord(*stop,0);
  cudaEventSynchronize(*stop);
  cudaEventElapsedTime(&bo,start,*stop);
  timer[index] += bo; 
#endif
  return;
}


__host__ void  prt_device_timing(float* timing,int size)
{
#ifdef  HAVE_GPU_CUDA_TM
  float accoum_time = 0.;
  printf( "___________________________________________________________________\n");
  printf( "________________  Timing on devices  ______________________________\n");
  cudaDeviceSynchronize();
  printf(" Allocation time        :    %10.3f (ms)\n", timing[0]);
  accoum_time += timing[0];
  /*FFT-TIMING*/
  accoum_time += timing[3];
  printf(" kernel execution FFT   :    %10.3f (ms)\n", timing[3]);
  /*KERNELS TIMERS*/
  for(unsigned int ii =1; ii<size; ++ii)
    {
      if(ii==3) continue;
      accoum_time += timing[ii];
      if(timing[ii]>1.e-6){
      printf(" debug execution   %02u   :    %10.3f (ms)\n",ii,timing[ii]);}
    }
  printf(" summed execution time  :    %10.3f (ms)\n", accoum_time) ; 

  /*END TOTAL TIMING*/
  printf( "___________________________________________________________________\n");
#endif
  return;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACE WITH FORTRAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

extern "C"
void prt_mem_use_(int* un_pitch,int* hauteur,int* size ){
  prt_mem_use((size_t) *un_pitch,(size_t) *hauteur,(size_t) *size );
  return;
}
