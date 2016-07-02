#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_sf_bessel.h>

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

void FC_FUNC_(gsl_f90_sf_bessel_j0,GSL_F90_SF_BESSEL_J0)
    (double *x, double *y)
{
  *y = gsl_sf_bessel_J0(*x);
}
