#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <Kokkos_Core.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stddef.h>

extern "C" {

  void c_abinit_kokkos_print_config()
  {

    std::cout << "##########################\n";
    std::cout << "KOKKOS CONFIG             \n";
    std::cout << "##########################\n";

    std::ostringstream msg;
    std::cout << "Kokkos configuration" << std::endl;
    if ( Kokkos::hwloc::available() ) {
      msg << "hwloc( NUMA[" << Kokkos::hwloc::get_available_numa_count()
          << "] x CORE["    << Kokkos::hwloc::get_available_cores_per_numa()
          << "] x HT["      << Kokkos::hwloc::get_available_threads_per_core()
          << "] )"
          << std::endl ;
    }
    //Kokkos::print_configuration( msg );
    std::cout << msg.str() << "##########################\n" << std::flush;

  } // c_abinit_kokkos_print_config

  void c_kokkos_initialize(int *argc, char **argv) {
    Kokkos::initialize(*argc, argv);
  }

  void c_kokkos_initialize_without_args() {
    Kokkos::initialize();
  }

  void c_kokkos_finalize() {
    Kokkos::finalize();
  }

  bool c_kokkos_is_initialized() {
    return Kokkos::is_initialized();
  }

} // extern "C"
