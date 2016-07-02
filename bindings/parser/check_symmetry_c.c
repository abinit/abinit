#include "ab7_symmetry.h"
#include "ab7_kpoints.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, const char **argv)
{
  Ab7Error err;
  Ab7Symmetry *sym;
  double rprimd[3][3] = { {0., 0.5, 0.5}, {1, 0., 0.5}, {1, 1, 0.}};
  int typeAt[2] = {1, 1};
  double xRed[6] = { 0., 0., 0., 0.25, 0.25, 0.25};
  double spinAt[6] = { 0., 0., 1., 0., 0., -1.};
  int i, multiplicity, nSym, spgrp, grpMagn, nkpt;
  Ab7SymmetryMat *syms;
  Ab7SymmetryTrans *trans;
  int *symAfm;
  char *sp;
  double magn[3];
  double *kpt, *wkpt;
  double shiftk[3] = {0.5, 0.5, 0.5};
  int ngkpt[3] = {2, 2, 2};
  
  sym = ab7_symmetry_new();
  if (!sym)
    {
      fprintf(stderr, "Can't create symmetry object.\n");
      return 1;
    }

  if (ab7_symmetry_set_lattice(sym, rprimd) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      fprintf(stderr, "Can't set the lattice of symmetry object.\n");
      return 1;
    }

  if (ab7_symmetry_set_structure(sym, 2, typeAt, xRed) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      fprintf(stderr, "Can't set the structure of symmetry object.\n");
      return 1;
    }

  if (ab7_symmetry_set_spin(sym, spinAt) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      fprintf(stderr, "Can't set the spin of symmetry object.\n");
      return 1;
    }

  err = ab7_symmetry_get_matrices(sym, &nSym, &syms, &trans, &symAfm);
  if (err != AB7_NO_ERROR && err != AB7_ERROR_SYM_BRAVAIS_XRED)
    {
      ab7_symmetry_free(sym);
      fprintf(stderr, "Can't get the matrices of symmetry object.\n");
      return 1;
    }
  fprintf(stdout, "nSym:%3d\n", nSym);
  for (i = 0 ; i < nSym ; i++)
    {
      fprintf(stdout, "sym%4.3d:\n", i + 1);
      fprintf(stdout, "%3d%3d%3d%12.6f%3d\n",
	      syms[i].mat[0][0], syms[i].mat[0][1], syms[i].mat[0][2],
	      trans[i].vect[0], symAfm[i]);
      fprintf(stdout, "%3d%3d%3d%12.6f%3d\n",
	      syms[i].mat[1][0], syms[i].mat[1][1], syms[i].mat[1][2],
	      trans[i].vect[1], symAfm[i]);
      fprintf(stdout, "%3d%3d%3d%12.6f%3d\n",
	      syms[i].mat[2][0], syms[i].mat[2][1], syms[i].mat[2][2],
	      trans[i].vect[2], symAfm[i]);
    }

  if (ab7_symmetry_get_multiplicity(sym, &multiplicity) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      return 1;
    }
  fprintf(stdout, "multiplicity:%3d\n", multiplicity);

  if (ab7_symmetry_get_group(sym, &sp, &spgrp, &grpMagn, magn) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      return 1;
    }
  fprintf(stdout, "space group:%s%6d\n", sp, spgrp);

  if (ab7_kpoints_get_mp_k_grid(sym, &nkpt, &kpt, &wkpt, ngkpt, 1, shiftk) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      return 1;
    }
  fprintf(stdout, "k-points (MP 2x2x2):%3d\n", nkpt);
  for (i = 0 ; i < nkpt ; i++)
    fprintf(stdout, "%3d:%10.6f%10.6f%10.6f%12.6f\n",
	    i, kpt[i * 3 + 0], kpt[i * 3 + 1], kpt[i * 3 + 2], wkpt[i]);
  free(kpt);
  free(wkpt);

  if (ab7_kpoints_get_auto_k_grid(sym, &nkpt, &kpt, &wkpt, 2.) != AB7_NO_ERROR)
    {
      ab7_symmetry_free(sym);
      return 1;
    }
  fprintf(stdout, "k-points (kptrlen = 2):%3d\n", nkpt);
  for (i = 0 ; i < nkpt ; i++)
    fprintf(stdout, "%3d:%10.6f%10.6f%10.6f%12.6f\n",
	    i, kpt[i * 3 + 0], kpt[i * 3 + 1], kpt[i * 3 + 2], wkpt[i]);
  free(kpt);
  free(wkpt);


  ab7_symmetry_free(sym);
  free(syms);
  free(symAfm);
  free(trans);
  free(sp);

  return 0;
}
