
#
# SCALAPACK : NETLIB or MKL ?
#
# default behavior is to detect Netlib scalapack, unless cmake var USE_MKL_SCALAPACK is TRUE
#

if (ABINIT_SCALAPACK_FLAVOR MATCHES "NETLIB")

  #
  # SCALAPACK/Netlib
  #

  # first try the MODULE mode (using pkg-config)
  # if not found try CONFIG mode
  # output: target library scalapack

  # Notes:
  # default scalapack library names are scalapack,
  # scalapack-openmpi and scalapack-mpich
  # e.g. on Ubuntu scalapack possible names are scalapack-openmpi,
  # scalapack-mpich
  #
  # if your scalapack has a different name, you can set variable
  # scalapack_name on the cmake configure command line,
  # and it will be check first in pkg_search_modules

  pkg_search_module(ScalapackPkg QUIET IMPORTED_TARGET ${scalapack_name} scalapack scalapack-openmpi scalapack-mpich)

  # On ubuntu 20.04, we must use pkg-config because the cmake config file is ill-formed
  # so here the logic is to first try to detect netlib/scalapack throught pkg-config
  # if not found, then try with cmake-config
  set(scalapack_found_using_pkg_config FALSE)
  set(scalapack_found_using_cmake_target FALSE)

  if(ScalapackPkg_FOUND)

    message("scalapack/netlib found via pkg-config")
    set(scalapack_FOUND TRUE)
    set(scalapack_found_using_pkg_config TRUE)
    add_library(abinit::scalapack ALIAS PkgConfig::ScalapackPkg)
    set(USING_SCALAPACK_NETLIB TRUE)

  else(ScalapackPkg_FOUND)

    message("scalapack/netlib not found via pkg-config")

    # if we are here, it means scalapack/netlib was not found by the MODULE mode (pkgconfig)
    # let try CONFIG mode
    find_package(scalapack CONFIG)
    if(scalapack_FOUND)
      message("scalapack/netlib found via cmake config")
      set(scalapack_found_using_cmake_target TRUE)
    else()
      message("scalapack/netlib not found. Please adjust variable scalapack_DIR / scalapack_ROOT")
    endif()

  endif(ScalapackPkg_FOUND)

  if(ScalapackPkg_FOUND OR scalapack_FOUND)
    set(HAVE_LINALG_SCALAPACK 1)
  endif()

elseif(ABINIT_SCALAPACK_FLAVOR MATCHES "MKL")

  #
  # SCALAPACK / MKL
  #

  #
  # Note : OneAPI MKL (version >= 2022.1.0) provides macro MKLConfig.cmake
  #

  #set(CMAKE_FIND_DEBUG_MODE TRUE)
  find_package(MKL)
  #set(CMAKE_FIND_DEBUG_MODE FALSE)

  if(MKL_FOUND)

    include(mkl_setup)

    add_library(abinit::scalapack ALIAS mkl::${MKL_SCALAPACK_FLAVOR})

    set(HAVE_LINALG_SCALAPACK 1)
    set(HAVE_LINALG_MKL 1)

  else(MKL_FOUND)

    message("MKL not found ! Please provide environment variable MKL_ROOT.")

  endif(MKL_FOUND)

endif()
