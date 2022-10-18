#include "m_xg_kokkos.h"

/**
 * Compute dot product of multiple pair of vectors.
 * Inputs are the column vectors of two 2d matrices x and y.
 * Results is a 1d vector containing of all (x dot y).
 *
 * This routine performs a total of ny dot products.
 *
 * Here we use kokkos hierarchical parallelism with two level parallelism:
 * - the outer level of parallelism maps the team (of threads) to the column index
 * - the inner level of parallelism uses a parallel reduction to compute a dot product
 *
 * An example of use of this functor is also available here:
 * https://github.com/pkestene/kokkos-proj-tmpl/blob/master/src/BatchedDotProduct.cpp
 *
 * \param[in] x_ptr pointer to a (nx,ny) 2d array
 * \param[in] y_ptr pointer to a (nx,ny) 2d array
 * \param[out] res_ptr pointer to the dot products, 1d array of size ny
 */
extern "C" void computeBatchedDotProduct_scalar_kokkos_cpp(
  const double *x_ptr,
  const double *y_ptr,
  double       *res_ptr,
  const int32_t nx,
  const int32_t ny)
{

  // create data views
  auto x   = AbiView_r64_2d_const(x_ptr, nx, ny);
  auto y   = AbiView_r64_2d_const(y_ptr, nx, ny);
  auto res = AbiView_r64_1d      (res_ptr, ny);

  // perform parallel computation on device
  {
    using team_policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
    using member_t = team_policy_t::member_type;

    // number of teams is the number of dot products to compute
    int32_t nbTeams = ny;

    // create a team policy
    const team_policy_t policy(nbTeams, Kokkos::AUTO(), Kokkos::AUTO());

    // define compute lambda
    auto dot_prod_lambda = KOKKOS_LAMBDA (const member_t& member)
    {
      // inside team, compute dot product as a parallel reduce operation
      double dot_prod = 0;
      int32_t j = member.league_rank();

      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, nx),
          [&](const int32_t &i, double &update)
          {
            update += x(i,j) * y(i,j);
          },
          dot_prod);

      // only one thread per team, collect the final reduction result, and write it
      // the output view
      Kokkos::single(Kokkos::PerTeam(member), [&]() { res(j) = dot_prod; });
    };

    Kokkos::parallel_for(
      "compute_dot_products_lambda",
      policy,
      dot_prod_lambda);

  }

} // computeBatchedDotProduct_scalar_kokkos_cpp

/**
 * Same as computeBatchedDotProduct_scalar_kokkos_cpp.
 *
 * TODO : refactor using template parameter to have one implementation for both
 * scalar / complex data type
 */
extern "C" void computeBatchedDotProduct_cplx_kokkos_cpp(
  const cplx_t *x_ptr,
  const cplx_t *y_ptr,
  cplx_t       *res_ptr,
  const int32_t nx,
  const int32_t ny)
{

  // create data views
  auto x   = AbiView_c64_2d_const(x_ptr, nx, ny);
  auto y   = AbiView_c64_2d_const(y_ptr, nx, ny);
  auto res = AbiView_c64_1d      (res_ptr, ny);

  // perform parallel computation on device
  {
    using team_policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
    using member_t = team_policy_t::member_type;

    // number of teams is the number of dot products to compute
    int32_t nbTeams = ny;

    // create a team policy
    const team_policy_t policy(nbTeams, Kokkos::AUTO(), Kokkos::AUTO());

    // define compute lambda
    auto dot_prod_lambda = KOKKOS_LAMBDA (const member_t& member)
    {
      // inside team, compute dot product as a parallel reduce operation
      cplx_t dot_prod = 0;
      int32_t j = member.league_rank();

      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, nx),
          [&](const int32_t &i, cplx_t &update)
          {
            update += Kokkos::conj(x(i,j)) * y(i,j);
          },
          dot_prod);

      // only one thread per team, collect the final reduction result, and write it
      // the output view
      Kokkos::single(Kokkos::PerTeam(member), [&]() { res(j) = dot_prod; });
    };

    Kokkos::parallel_for(
      "compute_dot_products_lambda",
      policy,
      dot_prod_lambda);

    //Kokkos::fence();

  }

} // computeBatchedDotProduct_cplx_kokkos_cpp

// =====================================================================
// =====================================================================
// =====================================================================

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMax_scalar_kokkos_cpp(
  const double *data_ptr,
  const int32_t size,
  double       *res
  )
{
  computeMax_kokkos_cpp<double>(data_ptr, size, res);
} // computeMax_scalar_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMax_complex_kokkos_cpp(
  const cplx_t *data_ptr,
  const int32_t size,
  double       *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  double reduced_value;
  Kokkos::Max<double> max_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMax",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, double& tmp) {
      double val = data_ptr[i].real();
      max_reducer.join(tmp, val);
    }, max_reducer);

  *res = reduced_value;
} // computeMax_complex_kokkos_cpp

// =====================================================================
// =====================================================================
// =====================================================================


// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMin_scalar_kokkos_cpp(
  const double *data_ptr,
  const int32_t size,
  double       *res
  )
{
  computeMin_kokkos_cpp<double>(data_ptr, size, res);
} // computeMin_scalar_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMin_complex_kokkos_cpp(
  const cplx_t *data_ptr,
  const int32_t size,
  double       *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  double reduced_value;
  Kokkos::Min<double> min_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMin",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, double& tmp) {
      double val = data_ptr[i].real();
      min_reducer.join(tmp, val);
    }, min_reducer);

  *res = reduced_value;
} // computeMin_complex_kokkos_cpp

// =====================================================================
// =====================================================================
// =====================================================================

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMaxloc_scalar_kokkos_cpp(
  const double *data_ptr,
  const int32_t size,
  int32_t      *res
  )
{
  computeMaxloc_kokkos_cpp<double>(data_ptr, size, res);
} // computeMaxloc_scalar_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMaxloc_scalar_2d_kokkos_cpp(
  const double *data_ptr,
  const int32_t nx,
  const int32_t ny,
  int32_t      *res
  )
{
  int32_t size = nx*ny;

  int32_t tmp;

  // do a 1d maxloc
  computeMaxloc_kokkos_cpp<double>(data_ptr, size, &tmp);

  // extract x,y coordinates
  // tmp = i + nx*j
  res[1] = tmp/nx;
  res[0] = tmp - nx*res[1];

  // return fortran offset values
  res[0] += 1;
  res[1] += 1;

} // computeMaxloc_scalar_2d_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMaxloc_complex_kokkos_cpp(
  const cplx_t *data_ptr,
  const int32_t size,
  int32_t      *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  using reducer_value_t = typename Kokkos::MinLoc<double, int>::value_type;
  reducer_value_t reduced_value;
  Kokkos::MaxLoc<double,int> maxloc_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMaxloc",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, reducer_value_t& tmp) {
      reducer_value_t val;
      val.val = data_ptr[i].real();
      val.loc = i;
      maxloc_reducer.join(tmp, val);
    }, maxloc_reducer);

  *res = reduced_value.loc;
} // computeMaxloc_complex_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMaxloc_complex_2d_kokkos_cpp(
  const cplx_t *data_ptr,
  const int32_t nx,
  const int32_t ny,
  int32_t      *res
  )
{
  int32_t size = nx*ny;

  int32_t tmp;

  // do a 1d maxloc
  computeMaxloc_complex_kokkos_cpp(data_ptr, size, &tmp);

  // extract x,y coordinates
  // tmp = i + nx*j
  res[1] = tmp/nx;
  res[0] = tmp - nx*res[1];

  // return fortran offset values
  res[0] += 1;
  res[1] += 1;

} // computeMaxloc_complex_2d_kokkos_cpp

// =====================================================================
// =====================================================================
// =====================================================================


// =====================================================================
// =====================================================================
// =====================================================================

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMinloc_scalar_kokkos_cpp(
  const double *data_ptr,
  const int32_t size,
  int32_t      *res
  )
{
  computeMinloc_kokkos_cpp<double>(data_ptr, size, res);
} // computeMinloc_scalar_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMinloc_scalar_2d_kokkos_cpp(
  const double *data_ptr,
  const int32_t nx,
  const int32_t ny,
  int32_t      *res
  )
{
  int32_t size = nx*ny;

  int32_t tmp;

  // do a 1d minloc
  computeMinloc_kokkos_cpp<double>(data_ptr, size, &tmp);

  // extract x,y coordinates
  // tmp = i + nx*j
  res[1] = tmp/nx;
  res[0] = tmp - nx*res[1];

  // return fortran offset values
  res[0] += 1;
  res[1] += 1;

} // computeMinloc_scalar_2d_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMinloc_complex_kokkos_cpp(
  const cplx_t *data_ptr,
  const int32_t size,
  int32_t      *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  using reducer_value_t = typename Kokkos::MinLoc<double, int>::value_type;
  reducer_value_t reduced_value;
  Kokkos::MinLoc<double,int> minloc_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMinloc",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, reducer_value_t& tmp) {
      reducer_value_t val;
      val.val = data_ptr[i].real();
      val.loc = i;
      minloc_reducer.join(tmp, val);
    }, minloc_reducer);

  *res = reduced_value.loc;
} // computeMinloc_complex_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeMinloc_complex_2d_kokkos_cpp(
  const cplx_t *data_ptr,
  const int32_t nx,
  const int32_t ny,
  int32_t      *res
  )
{
  int32_t size = nx*ny;

  int32_t tmp;

  // do a 1d minloc
  computeMinloc_complex_kokkos_cpp(data_ptr, size, &tmp);

  // extract x,y coordinates
  // tmp = i + nx*j
  res[1] = tmp/nx;
  res[0] = tmp - nx*res[1];

  // return fortran offset values
  res[0] += 1;
  res[1] += 1;

} // computeMinloc_complex_2d_kokkos_cpp

// =====================================================================
// =====================================================================
// =====================================================================

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeColwiseDivision_scalar_kokkos_cpp(
  const double *x_ptr,
  const double *y_ptr,
  const int32_t size,
  double       *res_ptr
  )
{

  computeColwiseDivision_kokkos_cpp<double>(x_ptr, y_ptr, size, res_ptr);

} // computeColwiseDivision_scalar_kokkos_cpp

// =====================================================================
//! C wrapper for iso_c_binding
extern "C"
void computeColwiseDivision_complex_kokkos_cpp(
  const cplx_t *x_ptr,
  const cplx_t *y_ptr,
  const int32_t size,
  cplx_t       *res_ptr
  )
{

  computeColwiseDivision_kokkos_cpp<cplx_t>(x_ptr, y_ptr, size, res_ptr);

} // computeColwiseDivision_complex_kokkos_cpp

// ================================================================
extern "C" void compute_colwiseCymax_scalar_kokkos_cpp(
  double        *A_ptr,
  const double  *da_ptr,
  const double  *B_ptr,
  const double  *W_ptr,
  const int32_t rows,
  const int32_t cols,
  const int32_t ldim)
{

  compute_colwiseCymax_kokkos_cpp<double>(A_ptr, da_ptr, B_ptr, W_ptr, rows, cols, ldim);

} // compute_colwiseCymax_scalar_kokkos_cpp

// ================================================================
extern "C" void compute_colwiseCymax_cplx_kokkos_cpp(
  cplx_t        *A_ptr,
  const double  *da_ptr,
  const cplx_t  *B_ptr,
  const cplx_t  *W_ptr,
  const int32_t  rows,
  const int32_t  cols,
  const int32_t  ldim)
{

  compute_colwiseCymax_kokkos_cpp<cplx_t>(A_ptr, da_ptr, B_ptr, W_ptr, rows, cols, ldim);

} // compute_colwiseCymax_cplx_kokkos_cpp
