#include <abi_common_kokkos.h>

extern "C" void add_array_kokkos_cpp(
  double *array1_ptr,
  double *array2_ptr,
  const int32_t array_size)
{

  Kokkos::parallel_for(
    "add_array_kokkos_kernel",
    Kokkos::RangePolicy<>(0,array_size),
    KOKKOS_LAMBDA(const int32_t i) {
      array1_ptr[i] += array2_ptr[i];
    });

} // add_array_kokkos_cpp
