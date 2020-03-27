/* cuda_common.h */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */


/** 
 * Language standards requires the existance of pre-defined macros
 * TODO 
 * Microsoft Visual C++ does not define __STDC__, 
 * Sun Workshop 4.2 supports C94 without setting __STDC_VERSION__ to the proper value
**/


#ifndef CUDA_COMMON_H
#define CUDA_COMMON_H

#include "config.h"

#undef HAVE_VARMACROS

#ifdef HAVE_GPU_CUDA_DP 
# define  CUDA_KIND (1.0_dp)
# define  cuzero   0.0
# define  cuone    1.0

# define  FFT_R2C  CUFFT_D2Z
# define  FFT_C2R  CUFFT_Z2D
# define  FFT_C2C  CUFFT_Z2Z

# if defined HAVE_VARMACROS  /**|| (__STDC_VERSION__ >= 199901L)**/
#  define  FFTEXECR2C(args...) cufftExecD2Z(args)
#  define  FFTEXECC2R(args...) cufftExecZ2D(args)
#  define  FFTEXECC2C(args...) cufftExecZ2Z(args)
# else
#  define  FFTEXECR2C cufftExecD2Z
#  define  FFTEXECC2R cufftExecZ2D
#  define  FFTEXECC2C cufftExecZ2Z
# endif /** defined(HAVE_VARMACROS) || (__STDC_VERSION__ >= 199901L) **/


#else
# define  CUDA_KIND (1.0)
# define  cuzero   0.0f
# define  cuone    1.0f

# define  FFT_R2C  CUFFT_R2C
# define  FFT_C2R  CUFFT_C2R
# define  FFT_C2C  CUFFT_C2C

# if defined HAVE_VARMACROS  /**|| (__STDC_VERSION__ >= 199901L)**/
#  define  FFTEXECR2C(args...) cufftExecR2C(args)
#  define  FFTEXECC2R(args...) cufftExecC2R(args)
#  define  FFTEXECC2C(args...) cufftExecC2C(args)
# else
#  define  FFTEXECR2C cufftExecR2C
#  define  FFTEXECC2R cufftExecC2R
#  define  FFTEXECC2C cufftExecC2C
# endif /** defined(HAVE_VARMACROS) || (__STDC_VERSION__ >= 199901L) **/


#endif


#endif
