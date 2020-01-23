//{\src2tex{textfont=tt}}
//****f* ABINIT/gpu_sphere
// NAME
//   gpu_sphere
//
// FUNCTION
// Array cg is defined in sphere with npw points. Insert cg inside box
// of n1*n2*n3 points to define array cfft for fft box.
// corresponds to given element in cg.  rest of cfft is filled with 0 s.
//
// There is also the possibility to apply a symmetry operation,
// as well as to make a shift in reciprocal space, or to multiply
// by a constant factor, in the case iflag=-1.
//
// COPYRIGHT
// Copyright (C) 1998-2020 ABINIT group (FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
//
// INPUTS
// cfft(2,n1,n2,n3) = fft box
// cg(2,npw)= contains values for npw G vectors in basis sphere
// compute_stream= stream cuda for kernels' execution
// istwfk=option parameter that describes the storage of wfs
// kg_k(3,npw)=integer coordinates of G vectors in basis sphere
// n1,n2,n3=physical dimension of the box (cfft)
// ndat=number of FFT to do in //
// npw=number of G vectors in basis at this k point
//
// SIDE EFFECTS
// Input/Output
//
// NOTES
// cg and cfft are assumed to be of type COMPLEX, although this routine treats
// them as real of twice the length to avoid nonstandard complex*16.
// If istwf_k differs from 1, then special storage modes must be taken
// into account, for symmetric wavefunctions coming from k=(0 0 0) or other
// special k points.
//
// TODO
//
//
// PARENTS
//      gpu_fourwf
//
// CHILDREN
//
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

__global__ void kernel_set_zero(double *tab,int n){
  int id= threadIdx.x + blockDim.x*(blockIdx.x + gridDim.x*blockIdx.y);
  if(id<n)
    tab[id]=0.;
}

__global__ void kernel_sphere_in(double *cfft,double *cg, const int *kg_k,
				 const int npw,const int ndat,const int n1,const int n2,const int n3,
				 const int shift_inv1,const int shift_inv2,const int shift_inv3,
				 const int istwfk){

  int ipw,idat,thread_id;
  int i1,i2,i3;
  thread_id  = threadIdx.x + blockIdx.x*blockDim.x;
  idat = blockIdx.y;

  for(ipw=thread_id; ipw<npw; ipw+=blockDim.x*gridDim.x){
    i1=kg_k[ipw*3];//kg_k(1,ipw)
    i2=kg_k[ipw*3 + 1];//kg_k(2,ipw)
    i3=kg_k[ipw*3 + 2];//kg_k(3,ipw)
    if(i1<0)
      i1+=n1;
    if(i2<0)
      i2+=n2;
    if(i3<0)
      i3+=n3;

    //We write cfft(i1,i2,i3)
    //(double2): cfft[i1 + n1*(i2 + n2*(i3 + n3*idat))] = cg[ipw + npw*idat]
    cfft[2*(i1 + n1*(i2 + n2*(i3+n3*idat)))] = cg[2*(ipw + npw*idat)];
    cfft[1+ 2*(i1 + n1*(i2 + n2*(i3+n3*idat)))] = cg[1 + 2*(ipw + npw*idat)];

    if(istwfk > 1){
      int i1inv,i2inv,i3inv;
      i1inv = (shift_inv1 - i1) % n1;
      i2inv = (shift_inv2 - i2) % n2;
      i3inv = (shift_inv3 - i3) % n3;
      //cfft(1,i1inv,i2inv,i3inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
      //cfft(2,i1inv,i2inv,i3inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
      cfft[2*(i1inv + n1*(i2inv + n2*(i3inv+n3*idat)))] = cg[2*(ipw + npw*idat)];
      cfft[1+ 2*(i1inv + n1*(i2inv + n2*(i3inv+n3*idat)))] = -cg[1 + 2*(ipw + npw*idat)];
    }
  }
}

__global__ void kernel_sphere_out(double *cfft,double *cg, int *kg_k,int npw,int ndat,int n1,int n2,int n3,double norm){

  int ig,idat,thread_id;
  int i1,i2,i3;
  thread_id = threadIdx.x + blockIdx.x*blockDim.x;
  idat = blockIdx.y;

  for(ig=thread_id; ig<npw; ig+=blockDim.x*gridDim.x){
    i1=kg_k[ig*3];//kg_k(1,ipw)
    i2=kg_k[ig*3 + 1];//kg_k(2,ipw)
    i3=kg_k[ig*3 + 2];//kg_k(3,ipw)
    if(i1<0)
      i1+=n1;
    if(i2<0)
      i2+=n2;
    if(i3<0)
      i3+=n3;

    //We write cg(ig)
    cg[2*(ig + npw*idat)]     = norm * cfft[2*(i1 + n1*(i2 + n2*(i3+n3*idat)))] ;
    cg[1 + 2*(ig + npw*idat)] = norm * cfft[1+ 2*(i1 + n1*(i2 + n2*(i3+n3*idat)))] ;
  }
}

/******************************************************************/
/*******                                                 **********/
/*******          FUNCTIONS TO BE CALLED                 **********/
/*******                                                 **********/
/******************************************************************/

extern "C" void gpu_sphere_in_(double *cg,double *cfft,int *kg_k,int *npw,int *n1,int *n2,int *n3,int *ndat,int *istwfk,cudaStream_t *compute_stream)
{

  //Arguments ------------------------------------
  //scalars
  //  integer,intent(in) :: istwfk,n1,n2,n3,ndat,npw
  //  cuda_stream_t, intent(in) :: compute_stream
  //arrays
  //  integer,intent(in) :: kg_k(3,npw)
  //  real(dp),intent(in) :: cg(2,npw*ndat)
  //  real(dp),intent(inout) :: cfft(2,n1,n2,n3*ndat)

  // //Local variables-------------------------------
  dim3 grid,bloc;
  int cfft_size=2*(*n1)*(*n2)*(*n3)*(*ndat);
  int shift_inv1,shift_inv2,shift_inv3;
  int istwf_k = *istwfk;
  // *************************************************************************

  //Set all work tab to zero
  bloc.x = BLOCK_SIZE;
  grid.x = min((cfft_size + bloc.x - 1 )/bloc.x,MAX_GRID_SIZE);
  grid.y = (cfft_size + bloc.x*grid.x - 1)/(bloc.x*grid.x);
  kernel_set_zero<<<grid,bloc,0,*compute_stream>>>(cfft,cfft_size);


  //During GPU calculation we do some pre-calculation on symetries
  if((istwf_k==2) || (istwf_k==4) || (istwf_k==6) || (istwf_k==8)){
    shift_inv1 = *n1;
  }
  else{
    shift_inv1 = *n1-1;
  }

  if((istwf_k>=2) && (istwf_k<=5)) {
    shift_inv2 = *n2;
  }
  else{
    shift_inv2 = *n2-1;
  }

  if((istwf_k==2) || (istwf_k==3) || (istwf_k==6) || (istwf_k==7)){
    shift_inv3 = *n3;
  }else{
    shift_inv3 = *n3-1;
  }


  grid.x = min((*npw  + bloc.x - 1 )/bloc.x,MAX_GRID_SIZE);
  grid.y = *ndat;
  //Call kernel to put cg into cfft
  kernel_sphere_in<<<grid,bloc,0,*compute_stream>>>(cfft,cg,kg_k,*npw,*ndat,*n1,*n2,*n3,shift_inv1,shift_inv2,shift_inv3,*istwfk);

}//end subroutine gpu_sphere_in



//{\src2tex{textfont=tt}}
//****f* ABINIT/sphere
// NAME
// gpu_sphere_out
//
// FUNCTION
// Array cg is defined in sphere with npw points. Extract cg from box
// of n1*n2*n3 points defined by array cfft for fft box.
//
// COPYRIGHT
// Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, AR)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
//
// INPUTS
// cfft(2,n1,n2,n3) = fft box
// cg(2,npw)= contains values for npw G vectors in basis sphere
// compute_stream= stream cuda for kernels' execution
// istwfk=option parameter that describes the storage of wfs
// kg_k(3,npw)=integer coordinates of G vectors in basis sphere
// n1,n2,n3=physical dimension of the box (cfft)
// ndat=number of FFT to do in //
// npw=number of G vectors in basis at this k point
//
// NOTES
// cg and cfft are assumed to be of type COMPLEX, although this routine treats
// them as real of twice the length to avoid nonstandard complex*16.
// If istwf_k differs from 1, then special storage modes must be taken
// into account, for symmetric wavefunctions coming from k=(0 0 0) or other
// special k points.
//
// TODO
// Order arguments
//
// PARENTS
//      gpu_fourwf
//
// CHILDREN
//
// SOURCE

extern "C" void gpu_sphere_out_(double *cg,double *cfft,int *kg_k,int *npw,int *n1,int *n2,int *n3,int* ndat,cudaStream_t *compute_stream)
{

   //Arguments ------------------------------------
  //scalars
  //  integer,intent(in) :: istwfk,n1,n2,n3,ndat,npw
  //  cuda_stream_t, intent(in) :: compute_stream
  //arrays
  //  integer,intent(in) :: kg_k(3,npw)
  //  real(dp),intent(in) :: cg(2,npw*ndat)
  //  real(dp),intent(inout) :: cfft(2,n1,n2,n3*ndat)

  // //Local variables-------------------------------
  dim3 grid,bloc;
  int cfft_size=(*n1)*(*n2)*(*n3);
  double norme=1./cfft_size;
  // *************************************************************************

  bloc.x = BLOCK_SIZE;
  grid.x = min((*npw  + bloc.x - 1 )/bloc.x ,MAX_GRID_SIZE);
  grid.y = *ndat;
  //Extract wave functions and appy fft normalisation factor before storing
  kernel_sphere_out<<<grid,bloc,0,*compute_stream>>>(cfft,cg,kg_k,*npw,*ndat,*n1,*n2,*n3,norme);

}//end subroutine gpu_sphere_out

//***
