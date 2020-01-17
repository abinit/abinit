//{\src2tex{textfont=tt}}
//****f* ABINIT/gpu_density_accumulation
// NAME
// gpu_density_accumulation
//
// FUNCTION
//  Kernels definitions and calling functions for
//  the accumulation of density with each wave function's contribution
//
// COPYRIGHT
// Copyright (C) 1998-2020 ABINIT group (FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
//
// INPUTS
// compute_stream = cuda stream to execute kernels
// denpot(n1,n2,n3) = holds density to store
// fofr(2,n1,n2,n3,ndat) = holds wave functions in real space after fft
// ndat = number of fft to do in //
// nttf_tot = size of the fft_box n1*n2*n3
// weight_i = weight of imaginary part in contribution
// weight_r = weight of real part in contribution
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


/*******************************************************************/
/**********                                      *******************/
/**********         kernels definitions          *******************/
/**********                                      *******************/
/*******************************************************************/

__global__ void kernel_accumulate_density(double *fofr,double* denpot,double* weight_r,double* weight_i,int nfft_tot,int ndat
					  ){
  int thread_id= threadIdx.x + blockDim.x*blockIdx.x;

  for(int id=thread_id; id <nfft_tot; id+=blockDim.x*gridDim.x){
    //double2 loc=fofr[id];
    //    denpot[id] += weight_r*loc.x*loc.x + weight_i*loc.y*loc.y ;
    int idat;
    double sum=0.;
    for(idat=0;idat<ndat;idat++){
      double loc_x=fofr[2*(id + nfft_tot*idat) ];
      double loc_y=fofr[1+2*(id + nfft_tot*idat)];
      sum += weight_r[idat]*loc_x*loc_x + weight_i[idat]*loc_y*loc_y ;
    }

    denpot[id] += sum ;
  }
}


/*******************************************************************/
/**********                                      *******************/
/**********          calling fonction            *******************/
/**********                                      *******************/
/*******************************************************************/

extern "C" void gpu_density_accumulation_(double *fofr,double* denpot, double* weight_r,
                   double* weight_i,int* nfft_tot,int *ndat,cudaStream_t *compute_stream)
{

  //Arguments ------------------------------------
  //scalars
  // integer intent(int) :: nfft_tot,ndat
  // double intent(int) :: weight_r,weight_i
  //arrays
  // double intent(inout) :: denpot
  // double intent(in) :: fofr
  //Locals
  dim3 grid,bloc;
  // *************************************************************************

  bloc.x = BLOCK_SIZE;
  grid.x = min(( *nfft_tot + bloc.x - 1 )/bloc.x,MAX_GRID_SIZE);

  //Call To Kernel
  kernel_accumulate_density<<<grid,bloc,0,*compute_stream>>>(fofr,denpot,weight_r,weight_i,*nfft_tot,*ndat);

}//end subroutine gpu_density_accumulation

//***
