//////////////////////////////////////////////
//  Routines that are callable from Fortran90
//  and use the levmar Levenberg-Marquardt
//  optimisation library
//////////////////////////////////////////////

#include <stdlib.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_LEVMAR
#include <float.h>
#include <levmar.h>

void dim_screening(double *p, double *y, int m, int n, void *adata)
{
 /* The function to be fitted evaluated at the values x
  * p[m] are the parameters
  * y[n] is the measurement vector (containing the y-values)
  * adata[n] is the optional data supplied (fixed) this will contain
  *  the z-values with adata[0],adata[2],...,adata[2*i] containing the
  *  real coordinate and adata[1],adata[3],...,adata[2*i+1] containing 
  *  the imaginary coordinate (0<=i<n) */

register int ip,in,idx;
double *data;
register double realp,imagp;
int npoles = m/3;

  data = (double *)adata;

  /* There is always at least one pole so do assignment in first loop*/
  for(in=0; in<n; ++in){

    realp = -data[2*in+1]*p[2]+p[1]*p[1]-data[2*in]*data[2*in];
    imagp = data[2*in]*(2.0*data[2*in+1]-p[2]);

    y[in] = -p[0]*imagp/(realp*realp+imagp*imagp);
  }

  if(npoles>1) { /* if there are more poles, add the rest*/
    for(in=0; in<n; ++in){
      for(ip=1; ip<npoles; ++ip){

        realp  = -data[2*in+1]*p[3*ip+2]+p[3*ip+1]*p[3*ip+1]-data[2*in]*data[2*in];
        imagp  = data[2*in]*(2.0*data[2*in+1]-p[3*ip+2]);

        y[in] += -p[3*ip]*imagp/(realp*realp+imagp*imagp);

      }
    }
  }
/* end */
}

void dre_and_im_screening(double *p, double *y, int m, int n, void *adata)
{
 /* The function to be fitted evaluated at the values x
  * p[m] are the parameters
  * y[n] is the measurement vector (containing the y-values)
  * adata[n] is the optional data supplied (fixed) this will contain
  *  the z-values with adata[0],adata[2],...,adata[2*i] containing the
  *  real coordinate and adata[1],adata[3],...,adata[2*i+1] containing 
  *  the imaginary coordinate (0<=i<n) */

register int ip,in,idx;
double *data;
register double realp,imagp;
register int npoles = m/3;

  data = (double *)adata;

  /* There is always at least one pole so do assignment in first loop*/
  for(in=0; in<n/2; ++in){

    realp = -data[2*in+1]*p[2]+p[1]*p[1]-data[2*in]*data[2*in];
    imagp = data[2*in]*(2.0*data[2*in+1]-p[2]);

    y[in    ] = -p[0]*imagp/(realp*realp+imagp*imagp);
    y[n/2+in] = -p[0]*realp/(realp*realp+imagp*imagp);
  }

  if(npoles>1) { /* if there are more poles, add the rest*/
    for(in=0; in<n/2; ++in){
      for(ip=1; ip<npoles; ++ip){

        realp  = -data[2*in+1]*p[3*ip+2]+p[3*ip+1]*p[3*ip+1]-data[2*in]*data[2*in];
        imagp  = data[2*in]*(2.0*data[2*in+1]-p[3*ip+2]);

        y[in    ] += -p[3*ip]*imagp/(realp*realp+imagp*imagp);
        y[n/2+in] += -p[3*ip]*realp/(realp*realp+imagp*imagp);

      }
    }
  }
/* end */
}

void FC_FUNC_(dfit_im_screening,DFIT_IM_SCREENING)
(double *re_zvals, double *im_zvals, double *yvals, \
 int *nvals, int *ncoeffs, double *coeffs, int *prtvol)
{
register int i;
int ret;
/* Transfer integer values */
int c_ncoeffs = ncoeffs[0];
int c_nvals   = nvals[0];
int c_prtvol  = prtvol[0];

double adata[2*c_nvals];
double lower_bounds[c_ncoeffs];
double upper_bounds[c_ncoeffs];

double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  //printf("\n From C - npoles: %i \n",c_nvals/3);
  //printf(  " From C - nvals: %i \n",c_nvals);
  //printf("\n From C - zvals:\n");
  //for(i=0; i<c_nvals; ++i) printf("  (%g,i%g)",re_zvals[i],im_zvals[i]);
  //printf("\n\n From C - yvals:\n");
  //for(i=0; i<c_nvals; ++i) printf("  %g",yvals[i]);

  /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  /* Initialise z-values */
  for(i=0; i<c_nvals; ++i){
    adata[2*i  ] = re_zvals[i];
    adata[2*i+1] = im_zvals[i];
  }
  //printf("\n From C - adata array:\n");
  //for(i=0; i<c_nvals; ++i) printf("  (%g,i%g)",adata[2*i],adata[2*i+1]);
 
  for(i=0; i<c_ncoeffs; i+=3){
    lower_bounds[i  ] = -DBL_MAX;
    lower_bounds[i+1] = 1E-16;
    lower_bounds[i+2] = -DBL_MAX;
    upper_bounds[i  ] = DBL_MAX;
    upper_bounds[i+1] = DBL_MAX;
    upper_bounds[i+2] = -1E-16;
  }
   
  /* invoke the optimisation function */
  //ret=dlevmar_dif(dim_screening, coeffs, yvals, c_ncoeffs, c_nvals, 5000, \
       opts, info, NULL, NULL, (void *)&adata); // without Jacobian

  //for(i=0; i<c_ncoeffs; i+=3){
  //  if (coeffs[i+2]>-1E-16){ 
      /* invoke the optimisation function with box boundaries*/
      ret=dlevmar_bc_dif(dim_screening, coeffs, yvals, c_ncoeffs, c_nvals, \
           lower_bounds, upper_bounds, 3000, opts, info, NULL, NULL, \
           (void *)&adata); // Box boundary conditions without Jacobian
    //}
  //}
  

  if(c_prtvol>9){
    printf("\n From C - Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
    printf("\n From C - Returned parameters:\n");
    for(i=0; i<c_ncoeffs/3; ++i) printf(" Peak %i - osc:%.7g w:%.7g gamma:%.7g\n\n",i,coeffs[3*i],coeffs[3*i+1],coeffs[3*i+2]);
  }
/* end*/
}

void FC_FUNC_(dfit_re_and_im_screening,DFIT_RE_AND_IM_SCREENING)
(double *re_zvals, double *im_zvals, double *im_yvals, \
 double *re_yvals, int *nvals, int *ncoeffs, double *coeffs, int *prtvol)
{
register int i;
int ret;
/* Transfer integer values */
int c_ncoeffs = ncoeffs[0];
int c_nvals   = 2*nvals[0];
int c_prtvol  = prtvol[0];

double adata[c_nvals];
double yvals[c_nvals];
double lower_bounds[c_ncoeffs];
double upper_bounds[c_ncoeffs];

double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

  //printf("\n From C - npoles: %i \n",c_nvals/3);
  //printf(  " From C - nvals: %i \n",c_nvals);
  //printf("\n From C - zvals:\n");
  //for(i=0; i<c_nvals; ++i) printf("  (%g,i%g)",re_zvals[i],im_zvals[i]);
  //printf("\n\n From C - yvals:\n");
  //for(i=0; i<c_nvals; ++i) printf("  %g",yvals[i]);

  /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
  opts[0]=LM_INIT_MU; opts[1]=1E-10; opts[2]=1E-10; opts[3]=1E-15;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  /* Initialise z-values */
  for(i=0; i<c_nvals/2; ++i){
    yvals[i]           = im_yvals[i];
    yvals[c_nvals/2+i] = re_yvals[i];
    adata[2*i  ] = re_zvals[i];
    adata[2*i+1] = im_zvals[i];
  }
  //printf("\n From C - adata array:\n");
  //for(i=0; i<c_nvals; ++i) printf("  (%g,i%g)",adata[2*i],adata[2*i+1]);
 
  for(i=0; i<c_ncoeffs; i+=3){
    lower_bounds[i  ] = -DBL_MAX;
    lower_bounds[i+1] = 1E-16;
    lower_bounds[i+2] = -DBL_MAX;
    upper_bounds[i  ] = DBL_MAX;
    upper_bounds[i+1] = DBL_MAX;
    upper_bounds[i+2] = -1E-16;
  }
   
  /* invoke the optimisation function */
  //ret=dlevmar_dif(dim_screening, coeffs, yvals, c_ncoeffs, c_nvals, 5000, \
       opts, info, NULL, NULL, (void *)&adata); // without Jacobian

  //for(i=0; i<c_ncoeffs; i+=3){
  //  if (coeffs[i+2]>-1E-16){ 
      /* invoke the optimisation function with box boundaries*/
      ret=dlevmar_bc_dif(dre_and_im_screening, coeffs, yvals, c_ncoeffs, c_nvals, \
           lower_bounds, upper_bounds, 1000, opts, info, NULL, NULL, \
           (void *)&adata); // Box boundary conditions without Jacobian
    //}
  //}
  

  if(c_prtvol>9){
    printf("\n From C - Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
    printf("\n From C - Returned parameters:\n");
    for(i=0; i<c_ncoeffs/3; ++i) printf(" Peak %i - osc:%.7g w:%.7g gamma:%.7g\n\n",i,coeffs[3*i],coeffs[3*i+1],coeffs[3*i+2]);
  }
/* end*/
}


#else
/* Dummy routines and error messages*/
#endif

