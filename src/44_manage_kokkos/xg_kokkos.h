#ifndef M_XG_KOKKOS_H
#define M_XG_KOKKOS_H

#include <abi_common_kokkos.h>

// =====================================================================
/**
 * Kokkos/c++ max reduction on device.
 *
 * template parameter value_t can be double or Kokkos::complex<double>
 *
 * \note array size is stored on a 32 bit integer (<= 2**31, i.e. 2 10^9)
 * we may need to use 64 bit consistently here and in m_xg.F90
 *
 */
template<typename value_t>
void computeMax_kokkos_cpp(
  const value_t *data_ptr,
  const int32_t  size,
  value_t       *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  value_t reduced_value;
  Kokkos::Max<value_t> max_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMax",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, value_t& tmp) {
      value_t val = data_ptr[i];
      max_reducer.join(tmp, val);
    }, max_reducer);

  *res = reduced_value;
} // computeMax_kokkos_cpp

// =====================================================================
/**
 * Kokkos/c++ min reduction on device.
 *
 * template parameter value_t can be double or Kokkos::complex<double>
 */
template<typename value_t>
void computeMin_kokkos_cpp(
  const value_t *data_ptr,
  const int32_t  size,
  value_t       *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  value_t reduced_value;
  Kokkos::Min<value_t> min_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMin",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, value_t& tmp) {
      value_t val = data_ptr[i];
      min_reducer.join(tmp, val);
    }, min_reducer);

  *res = reduced_value;
} // computeMin_kokkos_cpp

// =====================================================================
/**
 * Kokkos/c++ maxloc reduction on device.
 *
 * template parameter value_t can be double or Kokkos::complex<double>
 */
template<typename value_t>
void computeMaxloc_kokkos_cpp(
  const value_t *data_ptr,
  const int32_t  size,
  int32_t       *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  using reducer_value_t = typename Kokkos::MinLoc<value_t, int>::value_type;
  reducer_value_t reduced_value;
  Kokkos::MaxLoc<value_t,int> maxloc_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMaxloc",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, reducer_value_t& tmp) {
      reducer_value_t val;
      val.val = data_ptr[i];
      val.loc = i;
      maxloc_reducer.join(tmp, val);
    }, maxloc_reducer);

  *res = reduced_value.loc;
} // computeMaxloc_kokkos_cpp

// =====================================================================
/**
 * Kokkos/c++ minloc reduction on device.
 *
 * template parameter value_t can be double or Kokkos::complex<double>
 */
template<typename value_t>
void computeMinloc_kokkos_cpp(
  const value_t *data_ptr,
  const int32_t  size,
  int32_t       *res
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  using reducer_value_t = typename Kokkos::MinLoc<value_t, int>::value_type;
  reducer_value_t reduced_value;
  Kokkos::MinLoc<value_t,int> minloc_reducer(reduced_value);

  Kokkos::parallel_reduce(
    "computeMinloc",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i, reducer_value_t& tmp) {
      reducer_value_t val;
      val.val = data_ptr[i];
      val.loc = i;
      minloc_reducer.join(tmp, val);
    }, minloc_reducer);

  *res = reduced_value.loc;
} // computeMinloc_kokkos_cpp

// =====================================================================
/**
 * Kokkos/c++ equivalent to xgBlock_colwiseDivision.
 *
 * template parameter value_t can be double or Kokkos::complex<double>
 */
template<typename value_t>
void computeColwiseDivision_kokkos_cpp(
  const value_t *x_ptr,
  const value_t *y_ptr,
  const int32_t  size,
  value_t       *res_ptr
  )
{
  using policy_t = Kokkos::RangePolicy<>;

  Kokkos::parallel_for(
    "computeColwiseDivision",
    Kokkos::RangePolicy<>(0,size),
    KOKKOS_LAMBDA(const int32_t i)
    {
      res_ptr[i] = x_ptr[i] / y_ptr[i];
    });

} // computeColwiseDivision_kokkos_cpp

// =====================================================================
/**
 * Compute colwiseCymax (see fortran subroutine xgBlock_colwiseCymax in
 * m_xg.F90).
 *
 * \param[in,out] A_ptr pointer to a 2d array (xgBlock vecR)
 * \param[in    ] da_ptr pointer
 * \param[in    ] B_ptr
 * \param[in    ] W_ptr
 * \param[in    ] rows
 * \param[in    ] cols
 * \param[in    ] ldim
 */
template<typename value_t>
void compute_colwiseCymax_kokkos_cpp(
  value_t       *A_ptr,
  const double  *da_ptr,
  const value_t *B_ptr,
  const value_t *W_ptr,
  const int32_t  rows,
  const int32_t  cols,
  const int32_t  ldim)
{

  // create data views
  auto A   = AbiView_2d<value_t>      (A_ptr,  ldim, cols);
  auto da  = AbiView_1d_const<double> (da_ptr, cols);
  auto B   = AbiView_2d_const<value_t>(B_ptr,  ldim, cols);
  auto W   = AbiView_2d_const<value_t>(W_ptr,  ldim, cols);

  // perform parallel computation on device
  using policy = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;

  Kokkos::parallel_for(
    "colwiseCymax_kokkos_kernel",
    policy({0,0},{rows,cols}),
    KOKKOS_LAMBDA(const int32_t i, const int32_t j)
    {
      A(i,j) = - da(j) * B(i,j) + W(i,j);
    });

} // compute_colwiseCymax_kokkos_cpp

// =====================================================================
/**
 * Compute colwiseMul (see fortran subroutine xgBlock_colwiseMulR / xgBlock_colwiseMulC in
 * m_xg.F90).
 *
 * \param[in,out] data_ptr pointer to an array (xgBlock vecR / vecC)
 * \param[in    ] vec_ptr pointer
 * \param[in    ] shift
 * \param[in    ] rows
 * \param[in    ] cols
 * \param[in    ] ldim
 * \param[in    ] vec_size
 */
template<typename value_t, typename value2_t>
void compute_colwiseMul_kokkos_cpp(
  value_t        *data_ptr,
  const value2_t *vec_ptr,
  const int32_t  shift,
  const int32_t  rows,
  const int32_t  cols,
  const int32_t  ldim,
  const int32_t  vec_size)
{

  // create data views
  auto data = AbiView_2d<value_t>       (data_ptr, ldim, cols);
  auto vec  = AbiView_1d_const<value2_t>(vec_ptr, vec_size);

  // perform parallel computation on device
  using policy = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;

  Kokkos::parallel_for(
    "colwiseMul_kokkos_kernel",
    policy({shift,0}, {std::min(rows, shift+vec_size),cols}),
    KOKKOS_LAMBDA(const int32_t i, const int32_t j)
    {
      data(i,j) *= vec(i-shift);
    });

} // compute_colwiseMul_kokkos_cpp

#endif // M_XG_KOKKOS_H
