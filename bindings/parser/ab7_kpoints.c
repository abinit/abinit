/** - C source file -
 *
 * Copyright (C) 2009-2018 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include <stdlib.h>

#include <config.h>
#include <ab7_kpoints.fortran.h>
#include "ab7_kpoints.h"

Ab7Error ab7_kpoints_get_irreductible_zone(Ab7Symmetry *sym,
					   int *irrzon, double *phnons,
					   int n1, int n2, int n3,
					   int nsppol, int nspden)
{
  Ab7Error ab7_errno;

  KPT_CALL(get_irreductible_zone,GET_IRREDUCTIBLE_ZONE, sym, irrzon, phnons,
	   &n1, &n2, &n3, &nsppol, &nspden, &ab7_errno);
  return ab7_errno;
}

Ab7Error ab7_kpoints_get_mp_k_grid(Ab7Symmetry *sym, int *nkpt, double **kpt,
				   double **wkpt, const int ngkpt[3],
				   const int nshiftk, const double *shiftk)
{
  Ab7Error ab7_errno;
  double kptrlatt[9], kptrlen;
  double shiftk_[24];
  int nshiftk_, i;

  nshiftk_ = nshiftk;
  for (i = 0; i < nshiftk * 3; i++)
    shiftk_[i] = shiftk[i];
  KPT_CALL(binding_mp_k_1,BINDING_MP_K_1, sym, nkpt, ngkpt, kptrlatt, &kptrlen, &nshiftk_, shiftk_, &ab7_errno);
  if (ab7_errno != AB7_NO_ERROR && ab7_errno != AB7_ERROR_SYM_BRAVAIS_XRED) return ab7_errno;
  *kpt = g_malloc(sizeof(double) * (*nkpt) * 3);
  *wkpt = g_malloc(sizeof(double) * (*nkpt));
  KPT_CALL(binding_mp_k_2,BINDING_MP_K_2, sym, nkpt, *kpt, *wkpt, kptrlatt, &kptrlen, &nshiftk_, shiftk_, &ab7_errno);

  return ab7_errno;
}
Ab7Error ab7_kpoints_get_auto_k_grid(Ab7Symmetry *sym, int *nkpt, double **kpt,
				     double **wkpt, const double kptrlen)
{
  Ab7Error ab7_errno;
  double kptrlatt[9], kptrlen_;
  double shiftk[24];
  int nshiftk;

  kptrlen_ = kptrlen;
  KPT_CALL(binding_auto_k_1,BINDING_AUTO_K_1, sym, nkpt, kptrlatt, &kptrlen_, &nshiftk, shiftk, &ab7_errno);
  if (ab7_errno != AB7_NO_ERROR && ab7_errno != AB7_ERROR_SYM_BRAVAIS_XRED) return ab7_errno;
  *kpt = g_malloc(sizeof(double) * (*nkpt) * 3);
  *wkpt = g_malloc(sizeof(double) * (*nkpt));
  KPT_CALL(binding_auto_k_2,BINDING_AUTO_K_2, sym, nkpt, *kpt, *wkpt, kptrlatt, &kptrlen_, &nshiftk, shiftk, &ab7_errno);

  return ab7_errno;
}
