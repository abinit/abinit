# -*- Autoconf -*-
#
# Copyright (C) 2006-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Tricks for external libraries
#



# ABI_TRICKS_ATOMPAW(FC_VENDOR,FC_VERSION)
# ----------------------------------------
#
# Applies tricks and workarounds to have the AtomPAW library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_ATOMPAW],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  AC_MSG_NOTICE([applying AtomPAW tricks (vendor: $1, version: $2)])

  dnl Linear algebra
  tmpflags_atompaw='--with-linalg-libs="$(lib_linalg_libs)"'
  CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} ${tmpflags_atompaw}"

  dnl LibXC
  CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} --enable-libxc"
  tmpflags_atompaw='--with-libxc-incs="$(lib_libxc_incs)"'
  CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} ${tmpflags_atompaw}"
  tmpflags_atompaw='--with-libxc-libs="$(lib_libxc_libs)"'
  CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} ${tmpflags_atompaw}"

  dnl Force static build (shared libraries fail to build)
  CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} --enable-static --disable-shared"

  unset tmpflags_atompaw

  dnl Display changes
  AC_MSG_NOTICE([CFGFLAGS_ATOMPAW = ${CFGFLAGS_ATOMPAW}])
]) # ABI_TRICKS_ATOMPAW



# ABI_TRICKS_BIGDFT(FC_VENDOR,FC_VERSION)
# ---------------------------------------
#
# Applies tricks and workarounds to have the BigDFT library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_BIGDFT],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  AC_MSG_NOTICE([applying BigDFT tricks (vendor: $1, version: $2)])

  tmpflags_libxc='--disable-internal-libxc --with-libxc-incs="$(lib_libxc_incs)" --with-libxc-libs="$(lib_libxc_libs)"'
  tmpflags_libyaml='--enable-internal-libyaml --disable-shared'
  tmpflags_options='--without-archives --with-moduledir="$(fallbacks_instdir)/include"'
  tmpflags_bigdft='--disable-binaries --disable-bindings --enable-libbigdft'
  CFGFLAGS_BIGDFT="${CFGFLAGS_BIGDFT} ${tmpflags_bigdft} ${tmpflags_options} ${tmpflags_libyaml} ${tmpflags_libxc} --program-suffix='-abinit'"

  unset tmpflags_libxc
  unset tmpflags_libyaml
  unset tmpflags_options
  unset tmpflags_bigdft

  dnl Display changes
  AC_MSG_NOTICE([CFGFLAGS_BIGDFT = ${CFGFLAGS_BIGDFT}])
]) # ABI_TRICKS_BIGDFT



# ABI_TRICKS_ETSF_IO(FC_VENDOR,FC_VERSION)
# ----------------------------------------
#
# Applies tricks and workarounds to have the ETSF I/O library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_ETSF_IO],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  AC_MSG_NOTICE([applying ETSF_IO tricks (vendor: $1, version: $2)])

  tmpflags_etsf_io='--with-netcdf-incs="$(lib_netcdf_incs)"'
  CFGFLAGS_ETSF_IO="${CFGFLAGS_ETSF_IO} ${tmpflags_etsf_io}"
  tmpflags_etsf_io='--with-netcdf-libs="$(lib_netcdf_libs)"'
  CFGFLAGS_ETSF_IO="${CFGFLAGS_ETSF_IO} ${tmpflags_etsf_io}"
  tmpflags_etsf_io='--with-moduledir="$(fallbacks_instdir)/include"'
  CFGFLAGS_ETSF_IO="${CFGFLAGS_ETSF_IO} ${tmpflags_etsf_io}"

  case "$1" in
    ibm)
      FCFLAGS_ETSF_IO="${FCFLAGS_ETSF_IO} -qsuffix=cpp=f90:f=f"
      ;;
  esac

  unset tmpflags_etsf_io

  dnl Display changes
  AC_MSG_NOTICE([CFGFLAGS_ETSF_IO = ${CFGFLAGS_ETSF_IO}])
  AC_MSG_NOTICE([FCFLAGS_ETSF_IO  = ${FCFLAGS_ETSF_IO}])
]) # ABI_TRICKS_ETSF_IO



# ABI_TRICKS_FFTW()
# -----------------
#
# Applies tricks and workarounds to have the FFTW library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_FFTW],
[
  AC_MSG_NOTICE([applying FFTW tricks (vendor: $1, version: $2)])
]) # ABI_TRICKS_FFTW



# ABI_TRICKS_LIBXC(FC_VENDOR,FC_VERSION)
# --------------------------------------
#
# Applies tricks and workarounds to have the LIBXC library correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_LIBXC],
[
  AC_MSG_NOTICE([applying LIBXC tricks])

  dnl Enable the build of Fortran modules
  CFGFLAGS_LIBXC="--enable-fortran --enable-static --disable-shared"

  dnl Enforce C99 programming-style
  case "$1" in

    ibm)
      if test "${ac_cv_prog_cc_c99}" != "no"; then
        CFLAGS_LIBXC="${CFLAGS_LIBXC} ${ac_cv_prog_cc_c99}"
      fi
      ;;

  esac

  dnl Display changes
  AC_MSG_NOTICE([CFGFLAGS_LIBXC = ${CFGFLAGS_LIBXC}])
  AC_MSG_NOTICE([CFLAGS_LIBXC   = ${CFLAGS_LIBXC}])
]) # ABI_TRICKS_LIBXC



# ABI_TRICKS_LINALG(FC_VENDOR,FC_VERSION)
# ---------------------------------------
#
# Applies tricks and workarounds to have the optimized linear algebra
# libraries correctly linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_LINALG],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  AC_MSG_NOTICE([applying linear algebra tricks (vendor: $1, version: $2)])

]) # ABI_TRICKS_LINALG



# ABI_TRICKS_NETCDF(FC_VENDOR,FC_VERSION)
# ---------------------------------------
#
# Applies tricks and workarounds to have the optimized linear algebra
# libraries correctly linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_NETCDF],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  AC_MSG_NOTICE([applying NetCDF tricks (vendor: $1, version: $2)])

  CFGFLAGS_NETCDF="${CFGFLAGS_NETCDF} --disable-cxx --disable-cxx-4 --disable-dap --disable-hdf4 --disable-netcdf4 --enable-fortran"
  CFGFLAGS_NETCDF="${CFGFLAGS_NETCDF} --enable-static --disable-shared"
  CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DNDEBUG"

  case "$1" in

    g95)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -Df2cFortran"
      ;;

    gnu)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
      ;;

    ibm)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DIBMR2Fortran"
      FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -WF,-DIBMR2Fortran,-DNDEBUG"
      ;;

    intel)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
      case "$2" in
        9.0|9.1)
          FCFLAGS_NETCDF="${FCFLAGS_NETCDF} -mp"
          ;;
      esac
      ;;

    pathscale)
      case "$2" in
        1.0|4.0|5.0)
          CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
          ;;
        *)
          CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -Df2cFortran"
          ;;
      esac
      ;;

    open64)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -Df2cFortran -DF2CSTYLE"
      ;;

    pgi)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
      ;;

    sun)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DsunFortran"
      ;;

    *)
      CPPFLAGS_NETCDF="${CPPFLAGS_NETCDF} -DpgiFortran"
      ;;
  esac

  FCFLAGS_NETCDF="${CPPFLAGS_NETCDF} ${FCFLAGS_NETCDF}"

  dnl Display changes
  AC_MSG_NOTICE([CFGFLAGS_NETCDF = ${CFGFLAGS_NETCDF}])
  AC_MSG_NOTICE([CPPFLAGS_NETCDF = ${CPPFLAGS_NETCDF}])
  AC_MSG_NOTICE([FCFLAGS_NETCDF  = ${FCFLAGS_NETCDF}])
]) # ABI_TRICKS_NETCDF



# ABI_TRICKS_WANNIER90(FC_VENDOR,FC_VERSION)
# ------------------------------------------
#
# Applies tricks and workarounds to have the Wannier90 bindings correctly
# linked to the binaries.
#
AC_DEFUN([ABI_TRICKS_WANNIER90],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  AC_MSG_NOTICE([applying Wannier90 tricks (vendor: $1, version: $2)])

  dnl Disable shared libraries, because linear algebra fallback is
  dnl only providing static libraries.
  CFGFLAGS_WANNIER90="${CFGFLAGS_WANNIER90} --disable-shared --enable-static"

  dnl Linear algebra
  tmplibs_wannier90='$(lib_linalg_libs)'
  LIBS_WANNIER90="${tmplibs_wannier90} ${LIBS_WANNIER90}"

  case "$1" in

    intel)
      case "${target_cpu}" in
        ia64)
          # Do nothing
          ;;
        *)
          case "$2" in
            9.0|9.1|10.0|10.1|11.0)
            LIBS_WANNIER90="${LIBS_WANNIER90} -lsvml"
            ;;
          esac
          ;;
      esac
      ;;

  esac

  unset tmplibs_wannier90

  dnl Display changes
  AC_MSG_NOTICE([LIBS_WANNIER90 = ${LIBS_WANNIER90}])
]) # ABI_TRICKS_WANNIER90
