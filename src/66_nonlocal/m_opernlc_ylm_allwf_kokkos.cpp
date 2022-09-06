#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <Kokkos_Core.hpp>
#include <typeinfo>

#include <cstdio> // for printf
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream>
#include <cstdint>
#include <cassert>

// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#ifndef STR
#define STR(x) #x
#endif

#define MY_ASSERT(x)

#ifndef ABI_CHECK
#define ABI_CHECK(expr, msg) if (!(expr)) { printf("Abinit failed: (%s), in file %s, line %d : %s.\n", STR(expr), __FILE__, __LINE__, msg); abort(); }
#endif

using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
using memory_trait = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

//! device views types for array of doubles, unmanaged memory (allocated on fortran side)
using AbiView_r64_1d = Kokkos::View<double*, memory_space, memory_trait >;
using AbiView_r64_2d = Kokkos::View<double**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_r64_3d = Kokkos::View<double***, Kokkos::LayoutLeft,  memory_space, memory_trait >;
using AbiView_r64_4d = Kokkos::View<double****, Kokkos::LayoutLeft, memory_space, memory_trait >;

//! device views types for array of doubles, managed memory
using AbiView_r64_3d_managed = Kokkos::View<double***, Kokkos::LayoutLeft,  memory_space >;


//! device views types for array of integers, unmanaged memory (allocated on fortran side)
using AbiView_i32_1d = Kokkos::View<int32_t*, memory_space, memory_trait >;
using AbiView_i32_2d = Kokkos::View<int32_t**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_i32_3d = Kokkos::View<int32_t***, Kokkos::LayoutLeft,   memory_space, memory_trait >;


extern "C"
void opernlc_ylm_allwf_kokkos_cpp(int32_t cplex, int32_t cplex_enl, int32_t cplex_fac,
                                  int32_t dimenl1, int32_t dimenl2, int32_t dimekbq,
                                  int32_t iatm, int32_t itypat, int32_t ntypat, int32_t nprojs,
                                  int32_t natom, int32_t nincat, int32_t nspinor,
                                  int32_t nspinortot, int32_t paw_opt,
                                  int32_t nlmn, int32_t lmnmax,
                                  double* enl_gpu,
                                  double* gx_gpu,
                                  double* gxfac_gpu,
                                  double* gxfac2_gpu,
                                  double* gxfac_sij_gpu,
                                  int32_t shift_spinor, int32_t ndat,
                                  int32_t* atindx1_gpu,
                                  int32_t* indlmn_gpu,
                                  double* lambda_gpu,
                                  double* sij_typ_gpu,
                                  int32_t shift_proj, int32_t nattyp_max)
{

  // remember correspondances and ranges / sizes
  // gx_gpu        <=>     projections_cpu(:, ibeg:iend, 1:nspinor*ndat) ====> cplex,     nprojs, nspinor*ndat
  // gxfac_gpu     <=> vnl_projections_cpu(:, ibeg:iend, 1:nspinor*ndat) ====> cplex_fac, nprojs, nspinor*ndat
  // gxfac_sij_gpu <=>   s_projections_cpu(:, ibeg:iend, 1:nspinor*ndat) ====> cplex    , nprojs, nspinor*ndat

  // note  (fortran index)
  // ibeg = shift_proj + 1
  // iend = shift_proj + nattyp(itypat)*nlmn

  auto gx          = AbiView_r64_3d(gx_gpu,        cplex,     nprojs, nspinor*ndat);
  auto gxfac_sij   = AbiView_r64_3d(gxfac_sij_gpu, cplex,     nprojs, nspinor*ndat);

  auto gxfac_tmp   = AbiView_r64_3d(gxfac_gpu,     cplex_fac, nprojs, nspinor*ndat);
  auto gxfac_tmp2  = AbiView_r64_3d(gxfac2_gpu,    cplex_fac, nprojs, nspinor*ndat);

  auto indlmn      = AbiView_i32_3d(indlmn_gpu,    6,         lmnmax, ntypat);
  auto atindx1     = AbiView_i32_1d(atindx1_gpu,   natom);

  auto enl         = AbiView_r64_4d(enl_gpu,       dimenl1, dimenl2, nspinortot*nspinortot, dimekbq);

  auto lambda      = AbiView_r64_1d(lambda_gpu,    ndat);
  auto sij_typ     = AbiView_r64_1d(sij_typ_gpu,   ((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2);


  // loop starts at 1 (to mimic fortran)
  for (int iphase=0; iphase < dimekbq; ++iphase)
  {

    if (paw_opt==3)
      continue;

    auto gxfac = iphase == 0 ? gxfac_tmp : gxfac_tmp2;

    // notice : gxfac_gpu and gxfac2_gpu have been initialized (set to zero) in the calling fortran code
    // no need, to redo a reset here

    // zero

    //Accumulate gxfac related to non-local operator (Norm-conserving)
    //-------------------------------------------------------------------
    if (paw_opt==0)
    {

      // Enl is E(Kleinman-Bylander)
      ABI_CHECK(cplex_enl!=2, "BUG: invalid value; cplex_enl=2!");
      ABI_CHECK(cplex_fac==cplex,"BUG: invalid value; cplex_fac/=cplex!");
      {
        using policy = Kokkos::MDRangePolicy< Kokkos::Rank<4> >;
        Kokkos::parallel_for(
          "accum_gxfac_norm_conserving",
          policy({0,0,0,0},{nlmn,nincat,nspinor,ndat}),
          KOKKOS_LAMBDA(const int ilmn, const int ia, const int ispinor, const int idat)
          {
            int ispinor_index = ispinor + shift_spinor;
            int idat_ispinor = ispinor + idat*nspinor;
            int iln = indlmn(5-1, ilmn, itypat);
            int iproj = shift_proj + ilmn + ia*nlmn;

            double enl1 = enl(iln, itypat, ispinor_index, iphase);

            gxfac (0, iproj, idat_ispinor) = enl1 * gx(0, iproj, idat_ispinor);
            if (cplex>1)
              gxfac (1, iproj, idat_ispinor) = enl1 * gx(1, iproj, idat_ispinor);
          });
      }

    } // paw_opt == 0

    // Accumulate gxfac related to nonlocal operator (PAW)
    // -------------------------------------------------------------------
    if (paw_opt==1 or paw_opt==2 or paw_opt==4)
    {
      // !Enl is psp strength Dij or (Dij-lambda.Sij)

      // !  === Diagonal term(s) (up-up, down-down)

      if (cplex_enl==1)
      {
        // !  1-Enl is real

        {
          using policy = Kokkos::MDRangePolicy< Kokkos::Rank<4> >;
          Kokkos::parallel_for(
            "accum_gxfac_cplex_enl_real",
            policy({0,0,0,0},{nlmn,nincat,nspinor,ndat}),
            KOKKOS_LAMBDA(const int jlmn, const int ia, const int ispinor, const int idat)
            {

              int ispinor_index = ispinor + shift_spinor;
              int idat_ispinor = ispinor + idat*nspinor;
              int index_enl = atindx1(iatm+ia)-1; // fortran to c index
              int j0lmn = (jlmn+1)*jlmn/2; // jlmn*(jlmn-1)/2;
              int jjlmn = j0lmn + jlmn;

              double enl1 = enl(jjlmn, index_enl, ispinor_index, iphase);
              if (paw_opt==2) enl1 = enl1 - lambda(idat) * sij_typ(jjlmn);

              int iproj_j = shift_proj + jlmn + ia*nlmn;

              double gxfac_accum[2] = {
                enl1 * gx(0, iproj_j, idat_ispinor),
                enl1 * gx(1, iproj_j, idat_ispinor)
              };

              //double gxi[2] = {0, 0};

              for (int ilmn=0; ilmn<nlmn; ++ilmn)
              {
                int iproj_i = shift_proj + ilmn + ia*nlmn;

                if (ilmn < jlmn) {

                  int ijlmn = j0lmn + ilmn;

                  enl1 = enl(ijlmn, index_enl, ispinor_index, iphase);
                  if (paw_opt==2) enl1 = enl1 - lambda(idat) * sij_typ(ijlmn);

                  gxfac_accum[0] += enl1 * gx(0, iproj_i, idat_ispinor);
                  gxfac_accum[1] += enl1 * gx(1, iproj_i, idat_ispinor);

                } else if (ilmn > jlmn) {

                  if (jlmn<nlmn-1) {

                    int i0lmn = (ilmn+1)*(ilmn)/2;
                    int ijlmn = i0lmn + jlmn;

                    enl1 = enl(ijlmn, index_enl, ispinor_index, iphase);
                    if (paw_opt==2) enl1 = enl1 - lambda(idat) * sij_typ(ijlmn);

                    gxfac_accum[0] += enl1 * gx(0, iproj_i, idat_ispinor);
                    gxfac_accum[1] += enl1 * gx(1, iproj_i, idat_ispinor);
                  }

                }

              } // end for ilmn

              // write back accumulated result s
              gxfac(0, iproj_j, idat_ispinor) += gxfac_accum[0];
              gxfac(1, iproj_j, idat_ispinor) += gxfac_accum[1];

            }); // end parallel_for

        } // end kokkos kernel

      } // end if(cplex_enl == 1)

      else

      {
        // !  2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*

        ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!");

        if (nspinortot == 1) // -------------> NO SPINORS
        {

          {
            using policy = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
            Kokkos::parallel_for(
              "accum_gxfac_cplex_enl_complex",
              policy({0,0,0},{nlmn,nincat,ndat}),
              KOKKOS_LAMBDA(const int jlmn, const int ia, const int idat)
              {

                // ispinor is 1
                // nspinor should be 1 too (maybe an ABI_CHECK could ensure that)
                int idat_ispinor = idat;

                int index_enl = atindx1(iatm+ia)-1; // fortran to c index

                int j0lmn = (jlmn+1)*(jlmn)/2;
                int jjlmn = j0lmn + jlmn;

                double enl1 = enl(2*jjlmn, index_enl, 0, iphase);
                double enl2 = 0;
                if (paw_opt==2) enl1 = enl1 - lambda(idat) * sij_typ(jjlmn);

                int iproj_j = shift_proj + jlmn + ia*nlmn;

                double gxj[2] = {
                  gx(0, iproj_j, idat_ispinor),
                  gx(1, iproj_j, idat_ispinor)
                };

                double gxfac_accum[2] = {
                  enl1 * gxj[0],
                 cplex==2 ? enl1 * gxj[1] : 0
                };

                for (int ilmn=0; ilmn < jlmn-1; ++ilmn)
                {

                  int iproj_i = shift_proj + ilmn + ia*nlmn;

                  if (ilmn < jlmn) {

                    int ijlmn = j0lmn + ilmn;

                    //enl_(1:2) = enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1);
                    enl1 = enl(2*ijlmn  , index_enl, 0, iphase);
                    enl2 = enl(2*ijlmn+1, index_enl, 0, iphase);
                    if (paw_opt==2) enl1 = enl1 - lambda(idat) * sij_typ(ijlmn);

                    gxfac_accum[0] += enl1 * gx(0, iproj_i, idat_ispinor);
                    gxfac_accum[1] -= enl2 * gx(0, iproj_i, idat_ispinor);

                    if (cplex==2) {
                      gxfac_accum[0] += enl2 * gx(1, iproj_i, idat_ispinor);
                      gxfac_accum[1] += enl1 * gx(1, iproj_i, idat_ispinor);
                    }

                  } else if (ilmn > jlmn) {

                    if(jlmn<nlmn-1) {

                      int i0lmn = (ilmn+1)*(ilmn)/2;
                      int ijlmn = i0lmn + jlmn;

                      enl1 = enl(2*ijlmn  , index_enl, 0, iphase);
                      enl2 = enl(2*ijlmn+1, index_enl, 0, iphase);
                      if (paw_opt==2) enl1 = enl1 - lambda(idat) * sij_typ(ijlmn);

                      gxfac_accum[0] += enl1 * gx(0, iproj_i, idat_ispinor);
                      gxfac_accum[1] += enl2 * gx(0, iproj_i, idat_ispinor);

                      if (cplex==2) {
                        gxfac_accum[0] -= enl2 * gx(1, iproj_i, idat_ispinor);
                        gxfac_accum[1] += enl1 * gx(1, iproj_i, idat_ispinor);
                      } //end if cplex == 2

                    } // end if jlmn

                  } // end if ilmn

                } // end for ilmn

                // write back accumulated result s
                gxfac(0, iproj_j, idat_ispinor) += gxfac_accum[0];
                gxfac(1, iproj_j, idat_ispinor) += gxfac_accum[1];

              }); // end kokkos kernel

          } // end parallel_for

        } // end spinortot == 1

        else // nspinortot > 1
        {

          // -------------> SPINORIAL CASE

          // === Diagonal term(s) (up-up, down-down)

          // TODO : PK

        }

      } // cplex_enl==1

    } // end paw_opt

    // !End of loop when a exp(-iqR) phase is present
    // !------------------------------------------- ------------------------

    // !When iphase=0, gxfac and gxfac_tmp point to the same memory space
    // !When iphase=1, we add i.gxfac_tmp2 to gxfac_tmp
    if (iphase==1)
    {
      {
        using policy = Kokkos::MDRangePolicy< Kokkos::Rank<4> >;
        Kokkos::parallel_for(
          "accum_gxfac_iphase1",
          policy({0,0,0,0},{nlmn,nincat,nspinor,ndat}),
          KOKKOS_LAMBDA(const int ilmn, const int ia, const int ispinor, const int idat)
          {
            int idat_ispinor = ispinor + idat*nspinor;

            int iproj = shift_proj + ilmn + ia*nlmn;

            gxfac_tmp(0, iproj, idat_ispinor) -= gxfac_tmp2(1, iproj, idat_ispinor);
            gxfac_tmp(1, iproj, idat_ispinor) += gxfac_tmp2(0, iproj, idat_ispinor);
          });

      } // end kokkos kernel

    } // end if iphase==1

  } // end for iphase

  // !Accumulate gxfac related to overlap (Sij) (PAW)
  // !------------------------------------------- ------------------------
  if (paw_opt==3 or paw_opt==4) // ! Use Sij, overlap contribution
  {

    {
      using policy = Kokkos::MDRangePolicy< Kokkos::Rank<4> >;
      Kokkos::parallel_for(
        "accum_gxfac_sij",
        policy({0,0,0,0},{nlmn,nincat,nspinor,ndat}),
        KOKKOS_LAMBDA(const int jlmn, const int ia, const int ispinor, const int idat)
        {

          int idat_ispinor = ispinor + idat*nspinor;

          int j0lmn = (jlmn+1)*(jlmn)/2;
          int jjlmn = j0lmn + jlmn;
          int jlm = indlmn(4-1, jlmn, itypat);

          int iproj_j = shift_proj + jlmn + ia*nlmn;

          double sijr = sij_typ(jjlmn);
          double gxj[2] = {
            gx(0, iproj_j, idat_ispinor),
            gx(1, iproj_j, idat_ispinor)
          };

          double gxfac_sij_accum[2] = {0, 0};

          gxfac_sij_accum[0] += sijr * gxj[0];
          gxfac_sij_accum[1] += sijr * gxj[1];

          for(int ilmn = 0; ilmn<jlmn; ++ilmn)
          {

            int ilm = indlmn(4-1, ilmn, itypat);
            // if (ilm==jlm) {
            int ijlmn = j0lmn + ilmn;
            sijr = sij_typ(ijlmn);

            int iproj_i = shift_proj + ilmn + ia*nlmn;

            double gxi[2] = {
              gx(0, iproj_i, idat_ispinor),
              gx(1, iproj_i, idat_ispinor)
            };

            gxfac_sij_accum[0] += sijr * gxi[0];
            gxfac_sij_accum[1] += sijr * gxi[1];

            // } // end if ilm==jlm

          } // end for ilmn

          if(jlmn<nlmn-1)
          {

            for(int ilmn = jlmn+1; ilmn<nlmn; ++ilmn)
            {
              int ilm = indlmn(4-1, ilmn, itypat);

              int iproj_i = shift_proj + ilmn + ia*nlmn;

              // if (ilm==jlm) {

              int i0lmn = (ilmn+1)*(ilmn)/2;
              int ijlmn = i0lmn + jlmn;
              sijr = sij_typ(ijlmn);

              double gxi[2] = {
                gx(0, iproj_i, idat_ispinor),
                gx(1, iproj_i, idat_ispinor)
              };

              gxfac_sij_accum[0] += sijr * gxi[0];
              gxfac_sij_accum[1] += sijr * gxi[1];

              // } // end if ilm==jlm

            } // end do ilmn

          } // end if jlmn

          gxfac_sij(0, iproj_j, idat_ispinor) = gxfac_sij_accum[0];
          gxfac_sij(1, iproj_j, idat_ispinor) = gxfac_sij_accum[1];

        }); // end Kokkos kernel

    } // end parallel_for

  } // end if paw_opt

} // opernlc_ylm_allwf_kokkos_cpp
