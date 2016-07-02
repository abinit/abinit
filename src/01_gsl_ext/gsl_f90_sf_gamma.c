#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_sf_gamma.h>

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

void FC_FUNC_(gsl_f90_sf_gamma,GSL_F90_SF_GAMMA)
    (double *x, double *y)
{
  *y = gsl_sf_gamma(*x);
}
