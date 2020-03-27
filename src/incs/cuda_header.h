/* cuda_header.h */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef CUDA_HEADER_H
#define CUDA_HEADER_H



//#include <cutil.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include "cuda_common.h"
#include <stdio.h>

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_GPU_CUDA_DP
  typedef double cureal;
  typedef cufftDoubleComplex cucmplx;
  #define make_cuComplex(a,b) make_cuDoubleComplex(a,b)
#else
  typedef float cureal;
  typedef cufftComplex cucmplx;
  #define rsqrt(a) rsqrtf(a)
  #define sincos(a,b,c) __sincosf(a,b,c)
  #define make_cuComplex(a,b) make_cuFloatComplex(a,b)
  #define cuCreal(a)     cuCrealf(a)
  #define cuCimag(a)     cuCimagf(a)
  #define cuConj(a)      cuConjf(a)
  #define cuCadd(a,b)    cuCaddf(a,b)
  #define cuCmul(a,b)    cuCmulf(a,b)
  #define cuCdiv(a,b)    cuCdivf(a,b)
  #define cuCabs(a)      cuCabsf(a)
  #define cuCfma(a,b,c)  cuCfmaf(a,b,c)

#endif

typedef unsigned int uint;

void check_err(int);
//Since cufft doesn't have a readable return error code, we added it
static inline char* cufftGetErrorString(cufftResult res)
{
  char* message=(char*)malloc(25*sizeof(char));

  switch(res){
  case CUFFT_SUCCESS :
    sprintf(message,"CUFFT_SUCCESS");
    break;
  case CUFFT_INVALID_PLAN :
    sprintf(message,"CUFFT_INVALID_PLAN");
    break;
  case CUFFT_ALLOC_FAILED :
    sprintf(message,"CUFFT_ALLOC_FAILED");
    break;
  case CUFFT_INVALID_TYPE :
    sprintf(message,"CUFFT_INVALID_TYPE");
    break;
  case CUFFT_INVALID_VALUE :
    sprintf(message,"CUFFT_INVALID_VALUE");
    break;
  case CUFFT_INTERNAL_ERROR :
    sprintf(message,"CUFFT_INTERNAL_ERROR");
    break;
  case CUFFT_EXEC_FAILED :
    sprintf(message,"CUFFT_EXEC_FAILED");
    break;
  case CUFFT_SETUP_FAILED :
    sprintf(message,"CUFFT_SETUP_FAILED");
    break;
  case CUFFT_INVALID_SIZE :
    sprintf(message,"CUFFT_INVALID_SIZE");
    break;
  case CUFFT_UNALIGNED_DATA :
    sprintf(message,"CUFFT_UNALIGNED_DATA");
    break;
  default:
    sprintf(message,"CUFFT_UNKNOWN_ERROR");
  }

  return message;
}

#endif
