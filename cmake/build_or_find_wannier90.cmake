# recent version of wannier has pkg-config support for detection
# see PR : https://github.com/wannier-developers/wannier90/pull/406
#
# 1. If user set ABINIT_WANNIER90_BUILD, we download wannier sources and build them
# 2. Don't build wannier, but check if WANNIER_ROOT is defined and use find_package for setup
# 3. just detect libwannier using pkg-config
#

#
# Does abinit builds libwannier (https://github.com/wannier-developers/wannier90) ?
#
option(ABINIT_WANNIER90_BUILD "Turn ON if you want to build libwannier90 (default: OFF)" OFF)
option(ABINIT_WANNIER90_BUILD_FORCE "Enforce building libwannier90 ? (default: OFF)" OFF)

#
# Option to enable / disable wannier90 detection
#
option(ABINIT_WANNIER90_WANTED "Turn OFF if you don't want libwannier90 (default: ON)" ON)

if(ABINIT_WANNIER90_WANTED)

  # check if user requested a build of libwannier90
  # use carefully, it may strongly increase build time
  if(ABINIT_WANNIER90_BUILD)

    message("[abinit / wannier] Building wannier from source")

    set(WANNIER90_EXTERNAL wannier90_external)

    #set(WANNIER90_SRC_DIR      ${PROJECT_BINARY_DIR}/external/wannier90)
    set(WANNIER90_INSTALL_DIR  ${PROJECT_BINARY_DIR}/external/Install/${WANNIER90_EXTERNAL})
    set(WANNIER90_INC_DIR      ${WANNIER90_INSTALL_DIR}/include)
    set(WANNIER90_LIB_DIR      ${WANNIER90_INSTALL_DIR}/lib)

    if (ABINIT_WANNIER90_BUILD_FORCE)
      set(ABINIT_WANNIER90_BUILD_FORCE_BOOL True)
    else()
      set(ABINIT_WANNIER90_BUILD_FORCE_BOOL False)
    endif()

    include (ExternalProject)

    find_program(MAKE_EXE NAMES gmake nmake make)
    find_package(Git REQUIRED)

    set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/external)

    macro(select_make_inc)
      message("[abinit / wannier] build wannier90, run select_make_inc macro")
      if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(wannier_make_inc make.inc.gfort)
      elseif(CMAKE_Fortran_COMPILER_ID MATCHES "G95")
        set(wannier_make_inc make.inc.g95)
      elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
        set(wannier_make_inc make.inc.ifort)
      elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
        set(wannier_make_inc make.inc.nag)
      elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR
          CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        set(wannier_make_inc make.inc.pgf90)
      else()
        message(FATAL_ERROR "[abinit / wannier] CompilerID not supported")
      endif()
    endmacro()
    select_make_inc()

    # TODO : provide a URL with a tarball, in case git is down
    ExternalProject_Add (${WANNIER90_EXTERNAL}
      GIT_REPOSITORY https://github.com/wannier-developers/wannier90.git
      GIT_TAG d141f9f84dcd3ac54729b9e5874dabd451684237
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy config/${wannier_make_inc} make.inc
      BUILD_COMMAND ${MAKE_EXE} all pkgconfig PREFIX=${WANNIER90_INSTALL_DIR}
      BUILD_IN_SOURCE True
      INSTALL_COMMAND ${MAKE_EXE} install PREFIX=${WANNIER90_INSTALL_DIR}
      LOG_CONFIGURE 1
      LOG_BUILD 1
      LOG_INSTALL 1
      BUILD_ALWAYS 0
      )

    # create alias abinit::wannier through find_package
    set(WANNIER_ROOT ${WANNIER90_INSTALL_DIR})
    find_package(WANNIER)

    if(WANNIER_FOUND)
      message("[abinit / wannier] Wannier found via find_package")
      set(ABINIT_WANNIER90_FOUND True)
      set(HAVE_WANNIER90 1)
      # alias abinit::wannier is already defined inside find_package(WANNIER)
    else()
      message("[abinit / wannier] we shouldn't be here. We just build wannier library, and find_package should have succeeded !")
    endif()

    set(ABINIT_WANNIER90_BUILTIN TRUE)

  elseif (DEFINED WANNIER_ROOT OR DEFINED ENV{WANNIER_ROOT})

    find_package(WANNIER)

    if(WANNIER_FOUND)
      message("[abinit / wannier] Wannier found via find_package")
      set(ABINIT_WANNIER90_FOUND True)
      set(HAVE_WANNIER90 1)
      # alias abinit::wannier is already defined inside find_package(WANNIER)
    endif()

  else()

    # regular detection through PKG_CONFIG_PATH
    pkg_check_modules(ABINIT_WANNIER90 QUIET IMPORTED_TARGET wannier)
    if(ABINIT_WANNIER90_FOUND)
      add_library(abinit::wannier ALIAS PkgConfig::ABINIT_WANNIER90)
      set(HAVE_WANNIER90 1)
      message("[abinit / wannier] wannier FOUND via pkg-config")
    endif()

  endif()

else(ABINIT_WANNIER90_WANTED)

  message("[abinit / wannier] libwannier90 is not wanted")

endif(ABINIT_WANNIER90_WANTED)
