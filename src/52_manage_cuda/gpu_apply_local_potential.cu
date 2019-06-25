//****f* ABINIT/gpu_apply_local_potential
// NAME
// gpu_apply_local_potential
//
// FUNCTION
// Apply local potential to wfs in real space on the GPU
//
//
// COPYRIGHT
// Copyright (C) 1998-2019 ABINIT group (FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
//
// INPUTS
// compute_stream = cuda stream to execute kernels
// denpot(n1,n2,n3) = holds potential to apply
// fofr(n1,n2,n3,ndat) = holds wave functions in real space after fft
// ndat = number of fft to do in //
// nttf_tot = size of the fft_box n1*n2*n3
//
// SIDE EFFECTS
//
// NOTES
//
// PARENTS
//      gpu_fourwf
//
// CHILDREN
//      *none*
//
// SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "stdio.h"
#include "abi_gpu_header.h"


/******************************************************************/
/*******                                                 **********/
/*******          CUDA KERNELS DEFINITIONS               **********/
/*******                                                 **********/
/******************************************************************/

__global__ void kernel_apply_local_potential(double2 *fofr,double* denpot,int nfft_tot){
  int thread_id = threadIdx.x + blockDim.x*blockIdx.x;
  int idat = blockIdx.y;
  for(int id=thread_id;id<nfft_tot;id+=blockDim.x*gridDim.x){
    fofr[id + idat*nfft_tot].x *= denpot[id] ;
    fofr[id + idat*nfft_tot].y *= denpot[id] ;
  }
}



/******************************************************************/
/*******                                                 **********/
/*******          FUNCTIONS TO BE CALLED                 **********/
/*******                                                 **********/
/******************************************************************/

extern "C" void gpu_apply_local_potential_(double2 *fofr,double* denpot,int* nfft_tot,int *ndat,cudaStream_t *compute_stream)
{

  //Arguments ------------------------------------
  //scalars
  // integer intent(int) :: nfft_tot,ndat
  //arrays
  // Complex intent(inout) :: fofr
  // double intent(in) :: denpot
  //Locals
  dim3 grid,bloc;
  // *************************************************************************

  bloc.x = BLOCK_SIZE;
  grid.x = min(( *nfft_tot + bloc.x - 1 )/bloc.x,MAX_GRID_SIZE);
  grid.y = *ndat;

  //Call To Kernel
  kernel_apply_local_potential<<<grid,bloc,0,*compute_stream>>>(fofr,denpot,*nfft_tot);

}//end subroutine gpu_apply_local_potential

//***
