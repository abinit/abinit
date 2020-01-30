//{\src2tex{textfont=tt}}
//****f* ABINIT/gpu_mkkpg
// NAME
// gpu_mkkpg
//
// FUNCTION
// Compute all (k+G) vectors (in reduced coordinates) for given k point.
// Eventually compute related data
//
// COPYRIGHT
// Copyright (C) 1998-2020 ABINIT group (FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
//
// INPUTS
//  kg(3,npw)=integer coords of planewaves in basis sphere
//  kpt(3)=k point in terms of recip. translations
//  nkpg=second dimension of array kpg
//  npw=number of plane waves in reciprocal space
//
// OUTPUT
//  kpg(npw,3)= (k+G) components
//  === if nkpg==9 ===
//    kpg(npw,4:9)= [(k+G)_a].[(k+G)_b] quantities
//
// PARENTS
//      bloch_interp,ctocprj,dyfnl3,forstrnps,ks_ddiago,m_cprj_bspline
//      nstpaw3,nstwf3,prep_bandfft_tabs,resp3dte,rhofermi3,vtorho,vtorho3
//
// CHILDREN
//      leave_new,wrtout
//
// SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


__global__ void kernel_mkkpg(double *kpg,int *kg,
			     const int npw,
			     const double kpt_0,const double kpt_1,const double kpt_2
			     ){

  int ipw = threadIdx.x + blockIdx.x*blockDim.x;
  if(ipw<npw){
    //-- Compute (k+G) -- for mu = 1 --> 3
    kpg[ipw + 0*npw ] = kpt_0 + kg[0 + 3*ipw];
    kpg[ipw + 1*npw ] = kpt_1 + kg[1 + 3*ipw];
    kpg[ipw + 2*npw ] = kpt_2 + kg[2 + 3*ipw];
  }
}

extern "C" void gpu_mkkpg_(int *kg_gpu,double *kpg_gpu,double *kpt,int *npw){


//Arguments ------------------------------------
//scalars
// integer,intent(in) :: nkpg,npw
//arrays
// integer,intent(in) :: kg(3,npw)
// real(dp),intent(in) :: kpt(3)
// real(dp),intent(out) :: kpg(npw,nkpg)

//Local variables-------------------------------
//integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)

  dim3 grid,block;
  block.x = 64;
  grid.x = ((*npw) + block.x - 1)/block.x;
  kernel_mkkpg<<<grid,block>>>(kpg_gpu,kg_gpu,*npw,kpt[0],kpt[1],kpt[2]);
}
//***
