#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_sf_bessel.h>

#define FC_FUNC(name,NAME) name ## _
#define FC_FUNC_(name,NAME) name ## _

void FC_FUNC_(gsl_f90_sf_bessel_jl,GSL_F90_SF_BESSEL_JL)
    (int *l, double *x, double *y)
{
  *y = gsl_sf_bessel_jl(*l,*x);
}
