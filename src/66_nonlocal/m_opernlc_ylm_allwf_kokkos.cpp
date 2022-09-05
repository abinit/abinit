#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <Kokkos_Core.hpp>
#include <typeinfo>

#include <cstdio> // for printf
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream>
#include <cstdint>

extern "C"
void opernlc_ylm_allwf_kokkos_cpp(int32_t cplex, int32_t cplex_enl, int32_t cplex_fac,
                                  int32_t dimel1, int32_t dimel2, int32_t dimekbq,
                                  int32_t iatm, int32_t itypat,
                                  double* enl,
                                  int32_t natom, int32_t nincat, int32_t nspinor,
                                  int32_t nspinortot, int32_t paw_opt, int32_t nlmn,
                                  double * gx_gpu, double* gxfac_gpu, double* gxfac_sij_gpu,
                                  int32_t shift_spinor, int32_t ndat,
                                  int32_t* atindx1, int32_t* indlmn, double* lambda,
                                  double *sij,
                                  int32_t ibeg, int32_t iend, int32_t nattyp_max)
{
  // TODO

  // step 1 : extract range ibeg::iend from projections, vnl_projections and s_projections

  // step 2 : adapt code from opernlc_ylm_allwf_cpu

} // opernlc_ylm_allwf_kokkos_cpp
