/** C header file
 *
 * Copyright (C) 2008-2018 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef AB7_KPOINTS
#define AB7_KPOINTS

#include "ab7_base.h"
#include "ab7_symmetry.h"

Ab7Error ab7_kpoints_get_irreductible_zone(Ab7Symmetry *sym,
					   int *irrzon, double *phnons,
					   int n1, int n2, int n3,
					   int nsppol, int nspden);
Ab7Error ab7_kpoints_get_mp_k_grid   (Ab7Symmetry *sym, int *nkpt, double **kpt,
				      double **wkpt, const int ngkpt[3],
				      const int nshiftk, const double *shiftk);
Ab7Error ab7_kpoints_get_auto_k_grid (Ab7Symmetry *sym, int *nkpt, double **kpt,
				      double **wkpt, const double kptrlen);
#endif
