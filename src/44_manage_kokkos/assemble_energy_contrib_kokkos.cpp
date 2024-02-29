#include <abi_common_kokkos.h>

extern "C" void assemble_energy_contribution_kokkos_cpp(
  double *ghc_ptr,
  double *gsc_ptr,
  const double *kinpw_k2_ptr,
  const double *cwavef_ptr,
  const double *gvnlxc_ptr,
  const int32_t ndat,
  const int32_t my_nspinor,
  const int32_t npw_k2,
  const int32_t sij_opt,
  const bool k1_eq_k2,
  const double hugevalue)
{

  // create data views
  auto ghc      = AbiView_r64_2d      (ghc_ptr,      2, ndat*my_nspinor*npw_k2);
  auto gsc      = AbiView_r64_2d      (gsc_ptr,      2, ndat*my_nspinor*npw_k2);
  auto cwavef   = AbiView_r64_2d_const(cwavef_ptr,   2, ndat*my_nspinor*npw_k2);
  auto kinpw_k2 = AbiView_r64_1d_const(kinpw_k2_ptr,                   npw_k2);
  auto gvnlxc_  = AbiView_r64_2d_const(gvnlxc_ptr,   2, ndat*my_nspinor*npw_k2);

  // perform parallel computation
  {
    using policy = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
    Kokkos::parallel_for(
      "assemble_energy_contribution_kokkos_kernel",
      policy({0, 0, 0}, {npw_k2, my_nspinor, ndat}),
      KOKKOS_LAMBDA(const int ig, const int ispinor, const int idat)
      {

        constexpr int re = 0;
        constexpr int im = 1;
        constexpr double zero = 0.0;
        int igspinor = ig + npw_k2*ispinor + npw_k2*my_nspinor*idat;

        if ( kinpw_k2(ig) < hugevalue )
        {
          ghc(re,igspinor) = k1_eq_k2 == true ?
            ghc(re,igspinor) + kinpw_k2(ig)*cwavef(re,igspinor) + gvnlxc_(re,igspinor) :
            ghc(re,igspinor) + gvnlxc_(re,igspinor);

          ghc(im,igspinor) = k1_eq_k2 == true ?
            ghc(im,igspinor) + kinpw_k2(ig)*cwavef(im,igspinor) + gvnlxc_(im,igspinor) :
            ghc(im,igspinor) + gvnlxc_(im,igspinor);
        }
        else
        {
          ghc(re,igspinor) = zero;
          ghc(im,igspinor) = zero;

          if (sij_opt==1)
          {
            gsc(re,igspinor) = zero;
            gsc(im,igspinor) = zero;
          }
        }

      });
  }

} // assemble_energy_contribution_kokkos_cpp
