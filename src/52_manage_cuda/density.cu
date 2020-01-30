/* density.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "cuda_common.h"
#include "cuda_header.h"
#include "cuda_rec_head.h"
#include "rec_dens_calc.cu"

#define PI 3.1415926535897932384626

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* This CUDA module contains the functions (kernels) to perform
   density computation in the Recursion  Method DFT on GPU devices.

   Kernels function: density_kernel

   Host function: alloc_dens_cuda_, dealloc_dens_cuda_, density_cuda_
*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//- Global variables on GPU
cureal *an_d = NULL;
cureal *bn2_d = NULL;
cureal *rho_d = NULL;


//- Kernel to compute density
__global__ void
density_kernel(int ncolumn,
	       int npitch,
	       size_t nrow,
	       int dim_trott,
	       cureal pi_on_rtrotter,
	       cureal mult,
	       cucmplx cinv2rtrotter,
	       cucmplx coeef_mu,
	       cureal* an_d,
	       cureal* bn2_d,
	       cureal* rho_d
 )
{
 int idx = IMUL(blockDim.x,blockIdx.x) + threadIdx.x;
 //int itrot = IMUL(blockDim.y,blockIdx.y) + threadIdx.y;
 if(idx<nrow){
   rho_d[idx] = mult;
   //if(itrot<=dim_trott){
   for(int itrot=0;itrot<dim_trott+1;itrot++){
     if(itrot==0)   rho_d[idx] = mult;
     __syncthreads();
     cureal arg = pi_on_rtrotter*((cureal)itrot +.5);
     cucmplx zj;
     sincos(arg,&(zj.y),&(zj.x));
     zj = cuCmul(zj,coeef_mu);

     cucmplx N = {0.,0.};
     cucmplx D = {1.,0.};
     cucmplx Nold = {0.,0.};
     cucmplx Dold = {0.,0.};
     cucmplx facrec0 = {1.,0.};

     for(int ii=0;ii<ncolumn;ii++)
     {
       cucmplx an_c = {-an_d[idx*npitch+ii],0.0};
       cucmplx bn_c = {-bn2_d[idx*npitch+ii],0.0};

       if(bn_c.x==0.)  break;

//     Modification by MT - October,24th 2011
//     cucmplx Nnew = cuCaddf(cuCmul(zj,facrec0),cuCaddf(cuCmul(N,cuCaddf(zj,an_c)),cuCmul(Nold,bn_c)));
//     cucmplx Dnew = cuCaddf(cuCmul(D,cuCaddf(zj,an_c)),cuCmul(Dold,bn_c));
       cucmplx Nnew = cuCadd(cuCmul(zj,facrec0),cuCadd(cuCmul(N,cuCadd(zj,an_c)),cuCmul(Nold,bn_c)));
       cucmplx Dnew = cuCadd(cuCmul(D,cuCadd(zj,an_c)),cuCmul(Dold,bn_c));
       Nold = N;
       Dold = D;
       N = Nnew;
       D = Dnew;
       facrec0.x = 0.;
     }
     cucmplx ratio = cuCmul(cuCdiv(N,D),cinv2rtrotter);
     rho_d[idx] -= mult*ratio.x ;
   }
 }
}


//- Allocate recursion coefficients on GPU and copy CPU values on it
extern "C" __host__
void alloc_dens_cuda_(int* ntranche,
		      int* nrec,
		      int* dim_trott,
		      int* npitch,
		      cureal* an_h,
		      cureal* bn2_h

 )
{
 size_t largeur  = (size_t)(*nrec+1)*sizeof(cureal);   //- Size of real vectors
 size_t anpitch = largeur;  //- Pitch to put multi-vectors in a matrix: intial guess
 size_t height = *ntranche;

 cudaMallocPitch((void**) &an_d,&anpitch,largeur,height);
 cudaMemcpy2D(an_d,anpitch,an_h,largeur,largeur,height,cudaMemcpyHostToDevice);

 cudaMallocPitch((void**) &bn2_d,&anpitch,largeur,height);
 cudaMemcpy2D(bn2_d,anpitch,bn2_h,largeur,largeur,height,cudaMemcpyHostToDevice);

 cudaMalloc((void**) &rho_d,height*sizeof(cureal));

 *npitch = (size_t)anpitch/sizeof(cureal);

// prt_dbg_arr(bn2_d,200**npitch*sizeof(cureal),200,0,"lala_a");
// prt_dbg_arr(bn2_d,anpitch,10,(int)(anpitch/sizeof(cureal)),"lala_b");
}

//- Deallacate recursion coefficients on GPU
extern "C" __host__
void dealloc_dens_cuda_()
{
 cudaFree(an_d);
 cudaFree(bn2_d);
 cudaFree(rho_d);
 // printf("dealloc_dens_cuda\n");
}


//- Compute density for recursion
extern "C" __host__
void density_cuda_(int* npitch,       //- Pitch for the an and bn2 on gpu */
		   int* ntranche,     //- Number of point where compute density
		   int* nrec,         //- Max number of recursion
		   int* dim_trott,    //- Trotter dimension 2ptro-1
		   cureal* fermie,    //- Fermi energy
		   cureal* temp,      //- tsmear temperature
		   cureal* rtrotter,  //- Trotter real
		   cureal* infvol,    //- Infinitesimal volume
		   cureal* tollerace, //- tollerance
		   cureal* rho_h      //- density
 )
/* NOTES:
   When launched this function the recursion coefficients have to be
   copyed on the GPU previously.
*/

{
 /*------------------------------------------------------------------------------------*/
 /*---------------------------------INITIALIZATION-------------------------------------*/
 /*------------------------------------------------------------------------------------*/
 size_t height = *ntranche;
 cureal pi_on_rtrotter = PI/(*rtrotter);
 cucmplx cinv2rtrotter = {1./(*rtrotter*2.),0.};

 cufftDoubleComplex coeef_mut = {exp(-*fermie/(2.*(*rtrotter)*(*temp))), 0.};
 cucmplx coeef_mu = {(cureal) coeef_mut.x,(cureal) coeef_mut.y};

 //-Block and grid
 dim3 block(256,1,1);
 dim3 grid((height+block.x-1)/block.x,1,1);

 /* printf("rtrotter %f\n",*rtrotter); */
 /* printf("fermie %f\n",*fermie); */
 /* printf("dim_trott %d\n",*dim_trott);  */
 /* printf("cinv2rtrotter %f %f\n",cinv2rtrotter.x,cinv2rtrotter.y);  */
 /* printf("coeef_mu %f %f\n",coeef_mut.x,coeef_mut.y);   */
 /* printf("1/infvol %f\n",2./(*infvol));  */
 /* printf("nrec+1 %d\n",*nrec+1);  */
 /* printf("*npitch %d\n",*npitch);  */
 /* printf("height %d\n",height); */

 density_kernel <<< grid,block >>> (*nrec+1,*npitch,height,
				    *dim_trott,
				    pi_on_rtrotter,
				    2./(*infvol),
				    cinv2rtrotter,
				    coeef_mu,
				    an_d,bn2_d,
				    rho_d);

 cudaMemcpy(rho_h,rho_d,height*sizeof(cureal),cudaMemcpyDeviceToHost);

// for(int ii=0;ii<height;ii++) printf("%f ",rho_h[ii]);printf("\n");
 return;
}

