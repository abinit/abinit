/** C header file
 *
 * Copyright (C) 2008-2018 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef AB7_SYMMETRY
#define AB7_SYMMETRY

#include "ab7_base.h"

/**
 * Ab7Symmetry:
 *
 * An object to handle a set of symmetries.
 */
typedef int Ab7Symmetry;

/**
 * AB7_MAX_SYMMETRIES:
 *
 * Maximum size of internal arrays.
 */
#define AB7_MAX_SYMMETRIES 384

/**
 * Ab7SymmetryMat:
 *
 * A convenience name for 3x3 integer matrices. These matrices are
 * used to describe the symmetries.
 */
typedef struct Ab7SymmetryMat_ {int mat[3][3];} Ab7SymmetryMat;
/**
 * Ab7SymmetryTrans:
 *
 * A convenience name for 3 double vectors. These vectors are used to
 * represent the non symmorphic translations.
 */
typedef struct Ab7SymmetryTrans_ {double vect[3];} Ab7SymmetryTrans;


/**
 * ab7_symmetry_new:
 *
 * Create a new symmetry object. This object can't be used before
 * setting up the lattice with ab7_symmetry_set_lattice() and the
 * coordinates with ab7_symmetry_set_structure().
 *
 * Returns: a newly created object. Use ab7_symmetry_free() to free it.
 */
Ab7Symmetry* ab7_symmetry_new();
/**
 * ab7_symmetry_free:
 * @sym: a pointer on a #Ab7Symmetry object.
 *
 * Free all memory used by an #Ab7Symmetry object.
 */
void     ab7_symmetry_free          (Ab7Symmetry *sym);


Ab7Error ab7_symmetry_set_tolerance (Ab7Symmetry *sym, double tolsym);
/**
 * ab7_symmetry_set_lattice:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @rprimd: a matrix.
 *
 * Set up the lattice of the system. @rprimd describes the three
 * fundamental vectors of the basis set in dimensioned values. After
 * the lattice has been set, one can access the corresponding Bravais
 * lattice calling ab7_symmetry_get_bravais().
 *
 * Returns: #AB7_NO_ERROR if the lattice has been set successfully.
 */
Ab7Error ab7_symmetry_set_lattice   (Ab7Symmetry *sym, double rprimd[3][3]);
/**
 * ab7_symmetry_set_structure:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @natoms: the number of atoms.
 * @typeAt: the types of atoms. The types range if [1:maxtype].
 * @xRed: the coordinates in reduced coordinates.
 *
 * Set up the structure of the system. After the structure has been
 * set, one can access the symmetry matrices calling ab7_symmetry_get_matrices().
 *
 * Returns: #AB7_NO_ERROR if the structure has been set successfully or
 * #AB7_ERROR_SYM_BRAVAIS_XRED if the structure lattice found from
 * @xRed is less symmetric than the Bravais lattice.
 */
Ab7Error ab7_symmetry_set_structure (Ab7Symmetry *sym, int natoms,
				     int *typeAt, double *xRed);
/**
 * ab7_symmetry_set_spin:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @spinAt: a matrix.
 *
 * Set up the spin for each atom of the
 * system. ab7_symmetry_set_structure() must have been called before.
 *
 * Returns: #AB7_NO_ERROR if the value has been set successfully.
 */
Ab7Error ab7_symmetry_set_spin      (Ab7Symmetry *sym, double *spinAt);
/**
 * ab7_symmetry_set_collinear_spin:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @spinAt: a matrix.
 *
 * Set up the spin for each atom of the
 * system. ab7_symmetry_set_structure() must have been called before.
 *
 * Returns: #AB7_NO_ERROR if the value has been set successfully.
 */
Ab7Error ab7_symmetry_set_collinear_spin(Ab7Symmetry *sym, int *spinAt);
/**
 * ab7_symmetry_set_spin_orbit:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @status: a boolean.
 *
 * Set up if the system uses spin-orbit coupling or not.
 *
 * Returns: #AB7_NO_ERROR if the value has been set successfully.
 */
Ab7Error ab7_symmetry_set_spin_orbit(Ab7Symmetry *sym, gboolean status);
/**
 * ab7_symmetry_set_field:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @field: a vector.
 *
 * Set up if the system is under an electrical field. @field gives the
 * direction of the field.
 *
 * Returns: #AB7_NO_ERROR if the value has been set successfully.
 */
Ab7Error ab7_symmetry_set_field     (Ab7Symmetry *sym, double field[3]);
/**
 * ab7_symmetry_set_jellium:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @jellium: a boolean.
 *
 * Set up if the system is charged and uses a jellium compensation charge.
 *
 * Returns: #AB7_NO_ERROR if the value has been set successfully.
 */
Ab7Error ab7_symmetry_set_jellium   (Ab7Symmetry *sym, gboolean jellium);
/**
 * ab7_symmetry_set_periodicity:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @periodic: three booleans.
 *
 * Set up the system boundary conditions for the X, Y and Z
 * directions. If a direction is set to TRUE, then periodic
 * counditions are applied. Default is fully periodic.
 *
 * Returns: #AB7_NO_ERROR if the value has been set successfully.
 */
Ab7Error ab7_symmetry_set_periodicity(Ab7Symmetry *sym, gboolean periodic[3]);


/**
 * ab7_symmetry_get_n_atoms:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @nAtoms: a location to store an integer.
 *
 * Give the number of atoms defined by ab7_symmetry_set_structure().
 *
 * Returns: #AB7_NO_ERROR if the value has been get successfully.
 */
Ab7Error ab7_symmetry_get_n_atoms         (Ab7Symmetry *sym, int *nAtoms);
/**
 * ab7_symmetry_get_n_sym:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @nSym: a location to store an integer.
 *
 * Give the number of computed symmetries.
 *
 * Returns: #AB7_NO_ERROR if the value has been get successfully.
 */
Ab7Error ab7_symmetry_get_n_sym           (Ab7Symmetry *sym, int *nSym);
/**
 * ab7_symmetry_get_multiplicity:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @multiplicity: a location to store an integer.
 *
 * Give the multiplicity of the cell. If the cell is primitive,
 * @multiplicity is set to 1. In that case, one can call
 * ab7_symmetry_get_group() to have information on the space group.
 *
 * Returns: #AB7_NO_ERROR if the value has been get successfully.
 */
Ab7Error ab7_symmetry_get_multiplicity    (Ab7Symmetry *sym, int *multiplicity);
/**
 * ab7_symmetry_get_bravais:
 * @sym: a pointer on a #Ab7Symmetry object.
 * @bravais: a location for a matrix.
 * @holohedry: a location for an integer.
 * @center: a location for an integer.
 * @nBravSym: a location for an integer.
 * @bravSym: a pointer on #Ab7SymmetryMat.
 *
 * This routine is used to get the Bravais lattice of a system. The
 * @bravais lattice is given in integer value of the reciprocal space
 * vectors. The @holohedry is an integer [1:7] giving the lattice
 * system (1 is triclinic, 2 is monoclinic, 3 is orthorhombic...). The
 * @center is an other integer [-3:3] that specifies the Bravais
 * lattice (the corresponding values are ["F", "F", "I", "P", "A",
 * "B", "C"]). This routine also allocate an array of symmetry
 * matrices, stored in @bravSym. The number of allocated array in set
 * in @nBravSym. Use free() (or g_free()) to deallocate the Bravais
 * symmetry matrices after use.
 *
 * Returns: #AB7_NO_ERROR if the value has been get successfully.
 */
Ab7Error ab7_symmetry_get_bravais         (Ab7Symmetry *sym, int bravais[3][3],
					   int *holohedry, int *center, int *nBravSym,
					   Ab7SymmetryMat **bravSym);
Ab7Error ab7_symmetry_get_matrices        (Ab7Symmetry *sym, int *nSym,
					   Ab7SymmetryMat **syms,
					   Ab7SymmetryTrans **transNon, int **symAfm);
Ab7Error ab7_symmetry_get_group           (Ab7Symmetry *sym,
					   char **spaceGroup, int *spaceGroupId,
				           int *pointGroupMagn, double genAfm[3]);
Ab7Error ab7_symmetry_get_equivalent_atom (Ab7Symmetry *sym, int **equiv,
	                                   int *nSym, int iAtom);
Ab7Error ab7_symmetry_get_type            (Ab7Symmetry *sym, int *type, char **label,
                                           int iSym);
#endif
