/* prt_utils_rec.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "stdio.h"

/*=========================================================================*/
/*_________________________TIMING IN CUDA ROUTINES_________________________*/
/*=========================================================================*/
/* This file contains some basic utils from the time measuration in
 * cuda subroutines. A more particular version is contained in 
 * prt_utils_rec.cu (to put together)
*/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACE WITH FORTRAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
extern "C" __host__
void start_cuda_tm_(cudaEvent_t* start,cudaEvent_t* stop)
{
   cudaEventCreate(start);
   cudaEventCreate(stop);
   cudaEventRecord(*start,0);
   return;
};

extern "C" __host__
void stop_cuda_tm_(cudaEvent_t* stop)
{
   cudaEventRecord(*stop,0);
   cudaEventSynchronize(*stop);
   printf("stop %d\n",*stop);
   return;
}

extern "C" __host__
void calc_cuda_time_(cudaEvent_t* stop,cudaEvent_t* start,float* time_ms)
{
#if defined HAVE_GPU_CUDA3
   cudaThreadSynchronize();
#else
   cudaDeviceSynchronize();
#endif   
   *time_ms = 0.;
   stop_cuda_tm_(stop);
   cudaEventElapsedTime(time_ms,*start,*stop);
   printf("stop %d\n",*start);
   printf("stop %d\n",*stop);
   printf("stop %f\n",time_ms);
   cudaEventDestroy(*start);
   cudaEventDestroy(*stop);
   return ;
}
