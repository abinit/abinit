/* Fortran interface for m_ab7_symmetry.F90. */


#ifndef M_AB7_SYMMETRY_EXPORT_H
#define M_AB7_SYMMETRY_EXPORT_H

#define SYM_NAME(a,A)     ABI_FC_MOD(m_ab7_symmetry,M_AB7_SYMMETRY,	\
				     symmetry_ ## a,SYMMETRY_ ## A)
#define SYM_CALL(a,A,...) SYM_NAME(a,A)(__VA_ARGS__)


void SYM_NAME(new,NEW)            (int *dt);
void SYM_NAME(free,FREE)          (int *sym);

void SYM_NAME(set_tolerance,SET_TOLERANCE)          (int *sym, double *tolsym,
						     unsigned int *ab7_errno);
void SYM_NAME(set_lattice,SET_LATTICE)              (int *sym, double *rprimd,
						     unsigned int *ab7_errno);
void SYM_NAME(set_structure,SET_STRUCTURE)          (int *sym, int *natoms,
						     int *typeAt, double *xRed,
						     unsigned int *ab7_errno);
void SYM_NAME(set_collinear_spin,SET_COLLINEAR_SPIN)(int *sym, int *nAtoms,
						     int *spinAt,
						     unsigned int *ab7_errno);
void SYM_NAME(set_spin,SET_SPIN)                    (int *sym, int *nAtoms,
						     double *spinAt,
						     unsigned int *ab7_errno);
void SYM_NAME(set_spin_orbit,SET_SPIN_ORBIT)        (int *sym, int *status,
						     unsigned int *ab7_errno);
void SYM_NAME(set_field,SET_FIELD)                  (int *sym, double *field,
						     unsigned int *ab7_errno);
void SYM_NAME(set_jellium,SET_JELLIUM)              (int *sym, int *jellium, 
						     unsigned int *ab7_errno);
void SYM_NAME(set_periodicity,SET_PERIODICITY)      (int *sym, int *periodic, 
						     unsigned int *ab7_errno);

void SYM_NAME(get_n_atoms,GET_N_ATOMS)              (int *sym, int *nAtoms,
						     unsigned int *ab7_errno);
void SYM_NAME(get_n_sym,GET_N_SYM)                  (int *sym, int *nSym,
						     unsigned int *ab7_errno);
void SYM_NAME(get_multiplicity,GET_MULTIPLICITY)    (int *sym, int *multiplicity,
						     unsigned int *ab7_errno);
void SYM_NAME(get_bravais,GET_BRAVAIS)              (int *sym, int *bravais,
						     int *holohedry, int *center,
						     int *nBravSym, int *bravSym,
						     unsigned int *ab7_errno);
void SYM_NAME(get_matrices,GET_MATRICES)            (int *sym, int *nSym,
						     int *syms, double *transNon,
						     int *symAfm,
						     unsigned int *ab7_errno);
void SYM_NAME(get_group,GET_GROUP)                  (int *sym, char *spaceGroup,
						     int *spaceGroupId,
						     int *pointGroupMagn,
						     double *genAfm,
						     unsigned int *ab7_errno);
void SYM_NAME(get_equivalent_atom,GET_EQUIVALENT_ATOM) (int *sym, int *equiv,
							int *iAtom,
							unsigned int *ab7_errno);
void SYM_NAME(get_type,GET_TYPE)                    (int *sym, int *iSym,
                                                     char *label, int *type,
                                                     unsigned int *ab7_errno);

#endif
