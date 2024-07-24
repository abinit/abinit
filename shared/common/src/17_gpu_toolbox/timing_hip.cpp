/* timing_hip.cpp */

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "stdio.h"
#include <hip/hip_runtime_api.h>

/*=========================================================================*/
/*_________________________TIMING IN HIP ROUTINES__________________________*/
/*=========================================================================*/
/* This file contains some basic utils from the time measuration in
 * hip subroutines. A more particular version is contained in
 * prt_utils_rec.cu (to put together)
*/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACE WITH FORTRAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
extern "C"
void start_cuda_tm_(hipEvent_t* start,hipEvent_t* stop)
{
   hipEventCreate(start);
   hipEventCreate(stop);
   hipEventRecord(*start,0);
   return;
};

extern "C"
void stop_cuda_tm_(hipEvent_t* stop)
{
   hipEventRecord(*stop,0);
   hipEventSynchronize(*stop);
   printf("stop %d\n",*stop);
   return;
}

extern "C"
void calc_cuda_time_(hipEvent_t* stop,hipEvent_t* start,float* time_ms)
{
   hipDeviceSynchronize();
   *time_ms = 0.;
   stop_cuda_tm_(stop);
   hipEventElapsedTime(time_ms,*start,*stop);
   printf("stop %d\n",*start);
   printf("stop %d\n",*stop);
   printf("stop %f\n",time_ms);
   hipEventDestroy(*start);
   hipEventDestroy(*stop);
   return ;
}
