/* rec_dens_calc.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* This CUDA module contains the functions (host) to perform Recursion
   calculation of density during recursion method.
   This code have to be traslated on a kernel but it is not prioritary
   for the moment
*/


#include "cuda_common.h"
#include "cuda_header.h"
#include "cuda_rec_head.h"

cucmplx cone  = {1.,0.};
cucmplx czero = {0.,0.};

/* __host__ __device__ static __inline__  */
/* cucmplx mul(cucmplx a,cucmplx b ) */
/* { /\* Scalar multiplication between complexes number *\/ */
/*   cucmplx c; */
/*   c.x = a.x*b.x-a.y*b.y; */
/*   c.y = a.y*b.x+a.x*b.y; */
/*   return c; */
/* } */

/* __host__ __device__ static __inline__  */
/* cucmplx add(cucmplx a,cucmplx b ) */
/* { /\* Scalar addition  between complexes number *\/ */
/*   cucmplx c; */
/*   c.x = a.x+b.x; */
/*   c.y = a.y+b.y; */
/*   return c; */
/* } */

/* __host__ __device__ static __inline__  */
/* cucmplx sub(cucmplx a,cucmplx b ) */
/* { /\* Scalar addition  between complexes number *\/ */
/*   cucmplx c; */
/*   c.x = a.x-b.x; */
/*   c.y = a.y-b.y; */
/*   return c; */
/* } */

/* __host__ __device__ static __inline__  */
/* cucmplx div(cucmplx a,cucmplx b ) */
/* { */
/*   /\* Scalar division between complexes number *\/ */
/*   cucmplx conj; */
/*   cureal norm2 = b.x*b.x+b.y*b.y; */
/*   conj.x = b.x/norm2; conj.y = -b.y/norm2; */
/*   return mul(a,conj); */
/* } */

/* __host__ __device__ static __inline__  */
/* cucmplx scal(cucmplx a,cureal b ) */
/* { /\* Scalar multiplication between complexes number *\/ */
/*   cucmplx c; */
/*   c.x = a.x*b; */
/*   c.y = a.y*b; */
/*   return c; */
/* } */


__host__ void 
density_calc(const cureal betafermie, //-in- product of beta and fermi energy
	     const cureal mult,	      //-in- 2. divited inf_vol	
	     const cureal tolrec,     //-in- Tollerance of recursion
	     const int irec,          //-in- current recursion term
	     const int trotter,       //-in- trotter parameter
	     const int npt,           //-in- number of points on this gpu-cpu
	     const int loctranc,      //-in- number of points to compute in this step
	     const int delta_coor,    //-in- step to a add to find coordinates
	     int* contrec,            //-out- number of converged points
	     cureal* bn2,             //-inout- bn2 coeff
	     cureal* an,              //-in- an coeff
	     cureal* erreur,          //-inout- estimated error 
	     cureal* prod_b2,         //-inout- contains the product of bn2's
	     cucmplx* acc_rho,        //-inout- density
	     cucmplx* ND,             //-inout- numerator of cont.fraction
	     cucmplx* NDold,
	     cucmplx* NDnew)
{
  /* Density calculation */
  cureal POS_INF = 1.0/0.0;

  cureal rtrotter = max((cureal) trotter,.5);
  int loctrott = max(4*trotter,2);
  //  if(trotter==0) {rtrotter = .5;loctrott = 2;}
  cucmplx cmu = {exp((-betafermie/(2.*rtrotter))),0.};
  int loccontrec = 0;

  if(irec==0)
    {
      for( int ipt = 0;ipt<loctranc;ipt++){
	int coord = ipt+irec*npt + delta_coor;
	int accoord = 2*ipt;
	acc_rho[ipt] = czero;
	cureal dd = 1.;

	prod_b2[ipt] = 1.;
	erreur[accoord] =  0.;
	erreur[accoord+1] = 0.;
	for( int it=0;it<4*trotter;it+=2 )
	  {
	    int ifo = it + ipt*loctrott;
	    cureal arg = M_PI*((cureal)it/2. + .5 )/rtrotter;
	    cucmplx zj = {cos(arg) , sin(arg)};
	    cucmplx prod1 = cuCmul(zj,cmu);
	  
	    NDold[ifo] = czero;
	    NDold[ifo+1] = cone;

	    NDnew[ifo] = cone;
	    NDnew[ifo+1] = prod1;
	    NDnew[ifo+1].x -= an[coord]; 
	  
	    ND[ifo] = cone;
	    ND[ifo+1] = NDnew[ifo+1];

	    acc_rho[ipt] = cuCadd(acc_rho[ipt],cuCmul(prod1,cuCdiv(ND[ifo],ND[ifo+1])));
	    prod1 = cuCmul(ND[ifo+1],NDold[ifo+1]);
	    dd = sqrt(prod1.x*prod1.x+prod1.y*prod1.y);
	    erreur[accoord+1] +=  2.*rtrotter/dd;	
	  }	
	cucmplx prod = { -2.*rtrotter,0.};
	acc_rho[ipt] = cuCadd(cone, cuCdiv(acc_rho[ipt],prod));
      }
    }
  else{  
    for( int ipt = 0;ipt<loctranc;ipt++){
      int coord = ipt+irec*npt + delta_coor;
      int accoord = 2*ipt;
      acc_rho[ipt] = czero;
      cureal dd = 1.;

      prod_b2[ipt] *= exp(betafermie/(rtrotter)) * (bn2[coord]);
      erreur[accoord] =  erreur[accoord+1];
      erreur[accoord+1] = 0.;
      if( bn2[coord] != 0. ){
	for( int it=0;it<4*trotter;it+=2 )
	  {
	    int ifo = it + ipt*loctrott;
	    cureal arg = M_PI*((cureal)it/2. + .5 )/rtrotter;
	    cucmplx zj = {cos(arg) , sin(arg)};
	    cucmplx prod1 = cuCmul(zj,cmu);
	  
	    cucmplx an_c = {-an[coord],0.0};
	    cucmplx bn_c = {-bn2[coord],0.0};
	    NDnew[ifo] = cuCadd(cuCmul(cuCadd(prod1,an_c),ND[ifo]),cuCmul(NDold[ifo],bn_c));
	    NDnew[ifo+1] = cuCadd(cuCmul(cuCadd(prod1,an_c),ND[ifo+1]),cuCmul(NDold[ifo+1],bn_c));
	  
	    NDold[ifo] = ND[ifo];
	    NDold[ifo+1] = ND[ifo+1];
	    ND[ifo] = NDnew[ifo];
	    ND[ifo+1] = NDnew[ifo+1];

	    acc_rho[ipt] = cuCadd(acc_rho[ipt],cuCmul(prod1,cuCdiv(ND[ifo],ND[ifo+1])));
	    prod1 = cuCmul(ND[ifo+1],NDold[ifo+1]);
	    dd = sqrt(prod1.x*prod1.x+prod1.y*prod1.y);
	    erreur[accoord+1] += abs(prod_b2[ipt])/dd *2.*rtrotter;
	  }
      }
      cucmplx prod = { -2.*rtrotter,0.};
      acc_rho[ipt] = cuCadd(cone, cuCdiv(acc_rho[ipt],prod));
      
      if(irec>2) {
	if(bn2[coord+npt]<CUDA_TOL ||
	   mult*erreur[accoord+1] == POS_INF ||
	   (mult*erreur[accoord+1] < tolrec && 
	    mult*erreur[accoord] < tolrec))
	  {
	    //bn2[coord+npt] = 0.;
	    erreur[accoord+1] = 0.;
	    loccontrec ++;
	  }
      }

    }
  }

  //printf("irec %d contrec %d dens %e err %e \n",irec,loccontrec,.5*mult*acc_rho[0].x,mult*erreur[1]);
  *contrec = loccontrec;  
  return ;
}




