/* rec_kernel.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef __REC_KERNELS_H__
#define __REC_KERNELS_H__

#include "cuda_rec_head.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* This CUDA module contains some generic functions (kernels) to perform Recursion
   Method DFT on GPU devices.

   Kernels function:
   complex_prod - Complex multiplication of two vector of size M
   scalarProdGPU - computes the scalar prudocts of vectorN pairs of vectors of
   dimension elementN. The results is multiplied for scale factor
   (directly inspired from SDK).
   

   Host function:
   copytoconstmem - Copies in constant memory some variable used by
   all kernels
   find_positions - when gratio>1, it computes the position of in the
   grid where recursion has to apply.
*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*========================================================================*/
__global__ void 
setting_un(cureal* un,    //un(height_max,un_pitch)=delta(i) Initial vector
	   cureal* unold, //unold(height_max,un_pitch)=0 Old initial vector
	   cureal* vn,    //vn(height_max,un_pitch)=0  work vector
	   cureal* an,    //an(height_max)=0 recursion coeff an
	   cureal* bn2,   //bn2(height_max)=0 recursion coeff bn2
	   int pos0,      //first pt to calculate
	   int pos1,      //last pt to calculate
	   cureal scale)  //weight of the dirac delta
{ /* Kernel to initalise the arrays on the Device if gratio==1*/
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y) + threadIdx.y;
  
 if(x<dc_pthsize && y<dc_nptrec)
  {
   int ind = IMUL(y, dc_pthsize)+x;
   *(un+ind) = cuzero;
   *(unold+ind) = cuzero;
   *(vn+ind) = cuzero;
   
   if(x==pos0+y && x<pos1){
    *(un+ind) = scale;
    an[y]  = cuzero;
    bn2[y] = cuone;
   }
  }
}

/*========================================================================*/
__global__ void 
set_un_gratio(cureal* un,    //un(height_max,un_pitch)=delta(i) Initial vector
	      cureal* unold, //unold(height_max,un_pitch)=0 Old initial vector
	      cureal* vn,    //vn(height_max,un_pitch)=0  work vector
	      cureal* an,    //an(height_max)=0 recursion coeff an
	      cureal* bn2,   //bn2(height_max)=0 recursion coeff bn2
	      int* position, //positions where to do recursion
	      int pos0,      //first pt to calculate
	      int maxcoord,  //max of pt to calculate
	      cureal scale,  //weight of the dirac delta
	      int gratio)    //gratio                              
{ /* Kernel to initalise the arrays on the Device when grati>1*/
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y) + threadIdx.y;
 
 if(x<dc_pthsize && y<dc_nptrec)
  {
   int ind = IMUL(y, dc_pthsize)+x;
   *(un+ind) = cuzero;
   *(unold+ind) = cuzero;
   *(vn+ind) = cuzero;
   int coor = y+pos0;

   if(coor<maxcoord){
     if(x==position[coor]){
       *(un+ind) = scale;
       an[y] = cuzero;
       bn2[y] = cuone;
     }}
  }
}

/*========================================================================*/
__global__ void 
setting_un_cut(cureal* un,    //un(height_max,un_pitch)=delta(i) Initial vector
	       cureal* unold, //unold(height_max,un_pitch)=0 Old initial vector
	       cureal* vn,    //vn(height_max,un_pitch)=0  work vector
	       cureal* an,    //an(height_max)=0 recursion coeff an
	       cureal* bn2,   //bn2(height_max)=0 recursion coeff bn2
	       cureal scale,  //weight of the dirac delta
	       int target)    //position of the zero
{ /* Kernel to initalise the arrays on the Device if gratio==1*/
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y) + threadIdx.y;

 if(x<dc_pthsize && y<dc_nptrec)
  {
   int ind = IMUL(y, dc_pthsize)+x;
   *(un+ind) = cuzero;
   *(unold+ind) = cuzero;
   *(vn+ind) = cuzero;
   if(x==target){
    *(un+ind) = scale;
    an[y] = cuzero;
    bn2[y] = cuone;
   }
  }
}

/*========================================================================*/
__global__ void 
get_loc_potent(cureal* pot_gpu,     //potential on the ngfft grid
	       cureal* locpot_gpu, //potential on the ngfftrec grid
	       int3    trasl,      //the origin of ngfftrec grid
	       int     ipt,        //coordinate of the point for potential
	       int     ngfftrec,   //linear dim of the cutted grid
	       int     ngfft)      //linear dim of the grid
{ /* Kernel to re-positioning potential accorting to trasl for recrcut!=0 */
  /*we use in this kernel cubic grid (is very simple to extend)   */       

 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y) + threadIdx.y;
 int mult  = IMUL(ngfftrec,ngfftrec);
 int mult2 = IMUL(ngfft,ngfft);
 int part = x+IMUL(ngfftrec,y);
 int modj = y+trasl.y;
 int modi = x+trasl.x;

 if(modj>=ngfft) modj -= ngfft;
 if(modj<0) modj+=ngfft;
 if(modi>=ngfft)  modi -= ngfft;
 if(modi<0) modi+=ngfft;
  
 int part2 = modi+IMUL(ngfft,modj);
 for(int z = 0;z<ngfftrec;++z){
  
  if(x<ngfftrec && y<ngfftrec){
   for(int z = 0;z<ngfftrec;++z){
    int tot = part+IMUL(z,mult);
    int modk = z+trasl.z; 
    if(modk>=ngfft) modk -= ngfft; 
    if(modk<0) modk+=ngfft; 
    GET_TAB(locpot_gpu,tot,ipt,dc_pthsize) = pot_gpu[part2+IMUL(mult2,modk)];
   }  ; __syncthreads();
  };
 } 
}

/*========================================================================*/
__global__ void 
cmplxtoreal(cureal* a_r,  //out real vector
	    cucmplx* a_c, //in complex vector
	    int size)
{
 /* convert cucmplx vector in a cureal one */
 int id = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 if(id<size)  a_r[id] = a_c[id].x;
}


/*========================================================================*/
__global__ void 
realtocmplx(cureal* a_r,  //in real vector
	    cucmplx* a_c, //out complex vector
	    int size)
{
 /* convert cucmplx vector in a cureal one */
 int id = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 if(id<size) { a_c[id].x = a_r[id]; a_c[id].y = cuzero;}
}


/*========================================================================*/
__global__ void 
complex_prod(cucmplx* cvn,    //inout cvn(cucmplx) 
	     cucmplx* ZT_p,   //in  ZT_p(cucmplx)
	     int size)
{ /* Complex multiplication of two vector of size size*/
 int id = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 if(id<size){
   cvn[id] = cuCmul(ZT_p[id],cvn[id]);
 }
}


/*========================================================================*/
__global__ void 
complex_prod_tot(cucmplx* cvn,  //inout cvn(cucmplx) 
		 cucmplx* ZT_p, //in  ZT_p(cucmplx)
		 int height,    //height<height_max number of pts to calculate
		 int size)
{ /* Complex multiplication of two vector of size size*/
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 if(x<size && y<height){
  int indc = IMUL(y, dc_cvpthsz)+x;
  cvn[indc] = cuCmul(ZT_p[x],cvn[indc]);
 }
}

/*========================================================================*/
__global__ void 
scalarProdGPU(cureal* d_C,   //out d_C(vectorN)
	      cureal* d_A,   //in d_A(vectorN,elementN)
	      cureal* d_B,   //in d_B(vectorN,elementN)
	      cureal scale)  //in scale factor
{
 //Accumulators cache
 __shared__ cureal accumResult[ACCUM_N];

 ////////////////////////////////////////////////////////////////////////////
 // Cycle through every pair of vectors,
 // taking into account that vector counts can be different
 // from total number of thread blocks
 ////////////////////////////////////////////////////////////////////////////
 for(int vec = blockIdx.x; vec < dc_nptrec; vec += gridDim.x){
  int vectorBase = IMUL(dc_pthsize, vec);
  int vectorEnd  = vectorBase + dc_pthsize;

  ////////////////////////////////////////////////////////////////////////
  // Each accumulator cycles through vectors with
  // stride equal to number of total number of accumulators ACCUM_N
  // At this stage ACCUM_N is only preferred be a multiple of warp size
  // to meet memory coalescing alignment constraints.
  ////////////////////////////////////////////////////////////////////////
  for(int iAccum = threadIdx.x; iAccum < ACCUM_N; iAccum += blockDim.x){
   cureal sum = cuzero;
	
   for(int pos = vectorBase + iAccum; pos < vectorEnd; pos += ACCUM_N)
    sum += d_A[pos] * d_B[pos];
	
   accumResult[iAccum] = sum;
  }
  ////////////////////////////////////////////////////////////////////////
  // Perform tree-like reduction of accumulators' results.
  //ACCUM_N has to be power of two at this stage
  ////////////////////////////////////////////////////////////////////////
  for(int stride = ACCUM_N / 2; stride > 0; stride >>= 1){
   __syncthreads();
   for(int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x)
    accumResult[iAccum] += accumResult[stride + iAccum];
  }
  if(threadIdx.x == 0) d_C[vec] = accumResult[0]*scale;
 }
}


/*========================================================================*/
__global__ void 
un_x_pot(cucmplx* cvn,  //out vn(height_max,dc_pthsize)
	 cureal* un,    //in un(height_max,dc_pthsize)
	 cureal* pot,   //in pot(size)
	 int height)    //height<height_max number of pts to calculate
{ /* Kernel to multiply un(size*height) times pot(size) element-element */
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 if(x<dc_nfftrec && y<height) {
  int ind  = IMUL(y, dc_pthsize)+x;
  int indc = IMUL(y, dc_cvpthsz)+x;
  cvn[indc].x = un[ind]*pot[x];
  cvn[indc].y = cuzero;
 }
}


/*========================================================================*/

__global__ void 
un_x_pot_cut(cucmplx* cvn,  //out vn(height_max,dc_pthsize)
	     cureal* un,    //in un(height_max,dc_pthsize)
	     cureal* pot,   //in pot(size)
	     int height)    //height<height_max number of pts to calculate

{ /* Kernel to multiply un(size*height) times pot(size) element-element */
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 if(x<dc_nfftrec && y<height){
  int ind  = IMUL(y, dc_pthsize)+x;
  int indc = IMUL(y, dc_cvpthsz)+x;
  cvn[indc].x = un[ind]*pot[ind];
  cvn[indc].y = cuzero;
 }
}



/*========================================================================*/
__global__ void 
vn_x_pot_dv(cucmplx* cvn,  //in cvn(height_max,dc_pthsize)
	    cureal* vn,    //out vn(height_max,dc_pthsize)		     
	    cureal* pot,   //in pot(size)			     
	    cureal scale,  //scale factor				     
	    int height)    //height<height_max number of pts to calculate    
{ /* Kernel to multiply vn(size*height) times pot(size) and scale element- by-element */
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 volatile cucmplx locc;
 if(x<dc_nfftrec && y<height){
  int ind  = IMUL(y, dc_pthsize)+x;     
  int indc = IMUL(y, dc_cvpthsz)+x;
  locc.x = cvn[indc].x;
  locc.y = cvn[indc].y; //needed for performance
  vn[ind] = locc.x*pot[x]*scale;
 }
}


/*========================================================================*/
__global__ void 
vn_x_pot_dv_cut(cucmplx* cvn,  //in vn(height_max,dc_pthsize)
		cureal* vn,    //out vn(height_max,dc_pthsize)		     
		cureal* pot,   //in pot(size)			     
		cureal scale,  //scale factor				     
		int height)    //height<height_max number of pts to calculate    
{ /* Kernel to multiply vn(size*height) times pot(size) and scale element- by-element */
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 volatile cucmplx locc;
 if(x<dc_nfftrec && y<height){
  int ind  = IMUL(y, dc_pthsize)+x;
  int indc = IMUL(y, dc_cvpthsz)+x;
  locc.x = cvn[indc].x;
  locc.y = cvn[indc].y; //needed for performance
  vn[ind] = locc.x*pot[ind]*scale; 
 }
}


/*========================================================================*/
__global__ void 
un_invsqrt_scale(cureal* un,     //inout un(height_max,dc_pthsize)		       
		 cureal* scale,  //in scale(height_max)			     
		 int height)	 //height<height_max number of pts to calculate    
{ /* Kernel to multiply vn(size*height) times scale(height)  */
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 __shared__ cureal val;
 if(threadIdx.x==0) {
  val = rsqrt(scale[y]);
  if(scale[y] == cuzero ) val = cuzero; 
 }
 __syncthreads();
 if(x<dc_nfftrec && y<height)  
  GET_TAB(un,x,y,dc_pthsize) *= val;
}

/*========================================================================*/
__global__ void 
oldtonew(cureal* un,     //inout un(height_max,dc_pthsize)=delta(i) 
	 cureal* vn,     //in vn(height_max,dc_pthsize)=0  work vector		 
	 cureal* unold,  //inout unold(height_max,dc_pthsize)=0 Old  vector 
	 cureal* an,     //an(height_max)=0 recursion coeff an		 
	 cureal* bn2,    //bn2(height_max)=0 recursion coeff bn2		 
	 int height)     //height<height_max number of pts to calculate	 
 
{/* getting new un from vn,unold,an.bn */
 cureal switchu;
 __shared__ cureal anloc[1], bnloc[1];
 int x = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 int y = IMUL(blockDim.y,blockIdx.y);
 if(threadIdx.x==0) {anloc[0] = an[y]; bnloc[0] = sqrt(bn2[y]);}
 __syncthreads();
 
 if(x<dc_nfftrec && y<height) {
  int ind = IMUL(y, dc_pthsize)+x;
  switchu = *(un+ind);
  *(un+ind) = *(vn+ind) -anloc[0]*switchu-*(unold+ind)*bnloc[0];
  *(unold+ind) = switchu;
 }
}

/*========================================================================*/ 
__host__ void 
find_positions(const int3* pt0,  //first pt to calculate
               const int3* pt1,  //last pt to calculate
	       int delta,        //inital linear point
	       int final,        //final lineat point
	       int* pos_cpu,     //points where compute recursion
	       const int* ngfft, //grid 
	       int gratio)       //ratio between fine and coarse grid
{//intial point 
 int nn = 0;
 for(int kk=pt0->z;kk<1+pt1->z;kk+=gratio){
  int boz = kk*ngfft[1];
  for(int jj=0;jj<ngfft[1];jj+=gratio){
   int boy = (jj+boz)*ngfft[0];
   for(int ii=0;ii<ngfft[0];ii+=gratio){
    int tot = ii+boy;
    if(tot<delta) continue;
    if(tot>final) goto end_of_3loop;
    *(pos_cpu+nn) = tot;
    //printf("nn %d pos %d xx,yy,zz %d %d %d\n",pos_cpu[nn],nn,ii,jj,kk);
    nn++;
   }
  }
 }
  end_of_3loop:
 printf("ipoint local on rec %d\n",nn);
}

/*========================================================================*/
__host__ void 
copytoconstmem(int nfftrec,int nptrec,int pth_size,int cvpthsz)
{/*Host charge in constant memory some constant used in all kernels*/
 cudaMemcpyToSymbol( dc_nfftrec, &nfftrec, sizeof(int),0,cudaMemcpyHostToDevice);
 cudaMemcpyToSymbol( dc_nptrec,  &nptrec,  sizeof(int),0,cudaMemcpyHostToDevice);
 cudaMemcpyToSymbol( dc_pthsize, &pth_size,sizeof(int),0,cudaMemcpyHostToDevice);
 cudaMemcpyToSymbol( dc_cvpthsz, &cvpthsz, sizeof(int),0,cudaMemcpyHostToDevice);
}



/*========================================================================*/

#endif
