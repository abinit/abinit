/* Fortran interface for m_ab7_kpoints.F90. */


#ifndef M_AB7_KPOINTS_EXPORT_H
#define M_AB7_KPOINTS_EXPORT_H

#define KPT_NAME(a,A)     ABI_FC_MOD(m_ab7_kpoints,M_AB7_KPOINTS,	\
				     kpoints_ ## a,KPOINTS_ ## A)
#define KPT_CALL(a,A,...) KPT_NAME(a,A)(__VA_ARGS__)

void KPT_NAME(get_irreductible_zone,GET_IRREDUCTIBLE_ZONE) (int *sym, int *irrzon, double *phnons,
					 int *n1, int *n2, int *n3, int *nsppol,
					 int *nspden, unsigned int *ab7_errno);

void KPT_NAME(bindings_get_mp_k_grid_1,BINDINGS_GET_MP_K_GRID_1) (int *sym, int *nkpt,
					    int *ngkpt,
					    double *kptrlatt,
					    double *kptrlen,
					    int *nshiftk,
					    double *shiftk,
					    unsigned int *ab7_errno);
void KPT_NAME(bindings_get_mp_k_grid_2,BINDINGS_GET_MP_K_GRID_2) (int *sym, int *nkpt,
					    double *kpt,
					    double *wkpt,
					    double *kptrlatt,
					    double *kptrlen,
					    int *nshiftk,
					    double *shiftk,
					    unsigned int *ab7_errno);
void KPT_NAME(bindings_get_auto_k_grid_1,BINDINGS_GET_AUTO_K_GRID_1) (int *sym, int *nkpt,
					      double *kptrlatt,
					      double *kptrlen,
					      int *nshiftk,
					      double *shiftk,
					      unsigned int *ab7_errno);
void KPT_NAME(bindings_get_auto_k_grid_2,BINDINGS_GET_AUTO_K_GRID_2) (int *sym, int *nkpt,
					      double *kpt,
					      double *wkpt,
					      double *kptrlatt,
					      double *kptrlen,
					      int *nshiftk,
					      double *shiftk,
					      unsigned int *ab7_errno);

#endif
