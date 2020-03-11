/* cuda_rec_head.h */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef CUDA_REC_HEAD_H
#define CUDA_REC_HEAD_H

#include "config.h"
#include "cuda_header.h"


#ifdef  HAVE_GPU_CUDA_TM
        #define DEBUGLEN 13
        #define CUDA_SAFE_CALL(a)  a
        #define CUFFT_SAFE_CALL(a) a

#else 
        #define DEBUGLEN 0
        #define CUDA_SAFE_CALL(a)  a
        #define CUFFT_SAFE_CALL(a) a
        #define CUT_SAFE_CALL(a) a

#endif

//Setted tollerance for recursion
#ifdef HAVE_GPU_CUDA_DP
   #define CUDA_TOL 1.e-12
#else
   #define CUDA_TOL 1.e-07
#endif

#define ACCUM_N 1024
#define BLOCK_X 256  // CUDA block sizes

#define IMUL(a, b) __mul24(a, b)

#define GET_TAB(u,x,y,pitch) *(u + IMUL(y, pitch) + x)


/*-------- Constant Memory declaration -----------*/
__device__ __constant__ int  dc_nfftrec;
__device__ __constant__ int  dc_nptrec;
__device__ __constant__ int  dc_pthsize;
__device__ __constant__ int  dc_cvpthsz;


//interfaces
void starttime(cudaEvent_t*);
void calctime(cudaEvent_t*,cudaEvent_t,float*,int);
void prt_dbg_arr(cureal*,size_t,int,int,char[20]);
void prt_dbg_arrc(cucmplx*,size_t,int,int,char[20]);
void prt_mem_use(size_t,size_t,size_t );
void prt_device_timing(float*,int);
float get_max_mem_dev(int);

/* __host__ __device__ static __inline__ cucmplx mul(cucmplx a,cucmplx b ); */
/* __host__ __device__ static __inline__ cucmplx add(cucmplx a,cucmplx b ); */
/* __host__ __device__ static __inline__ cucmplx div(cucmplx a,cucmplx b ); */

void density_calc(const cureal,const cureal,const cureal, 
		  const int,const int, const int,const int,
		  const int,int*, cureal*,cureal*, cureal*, cureal*, 
		  cucmplx*,cucmplx*,cucmplx*, cucmplx* );

/* //interfaces kernels_rec */
__global__ void cmplxtoreal(cureal*,cucmplx*,int);
__global__ void realtocmplx(cureal*,cucmplx*,int);
__global__ void complex_prod(cucmplx*,cucmplx*,int);
__global__ void complex_prod_tot(cucmplx*,cucmplx*,int,int);

__global__ void setting_un(cureal*,cureal*,cureal*,cureal*,cureal*,int,int,cureal); 
__global__ void set_un_gratio(cureal*,cureal*,cureal*,cureal*,cureal*,int*,int,int,cureal,int);
__global__ void un_x_pot(cucmplx*,cureal*,cureal*,int );
__global__ void vn_x_pot_dv(cucmplx*,cureal*,cureal*,cureal,int );

__global__ void setting_un_cut(cureal*,cureal*,cureal*,cureal*,cureal*,cureal,int);
__global__ void get_loc_potent(cureal*,cureal*,int3,int,int,int );
__global__ void un_x_pot_cut(cucmplx*,cureal*,cureal*,int );
__global__ void vn_x_pot_dv_cut(cucmplx*,cureal*,cureal*,cureal,int );

__global__ void un_invsqrt_scale(cureal*,cureal*,int );
__global__ void oldtonew(cureal*,cureal*,cureal*,cureal*,cureal*,int );
__global__ void scalarProdGPU(cureal*,cureal*,cureal*,cureal); 


__host__ void copytoconstmem(int,int,int,int);

__host__ void find_positions(const int3*,const int3*,int,int,int*,const int*,int );

#endif
