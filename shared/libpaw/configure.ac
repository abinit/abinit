# Copyright (C) 2017 Yann Pouillon <devops@materialsevolution.es>

#
# Process this file with autoconf to produce a configure script.
#

# Enhance help message (must be set before anything else)
m4_divert_text([HELP_END], [
More information on how to configure the build of LibPAW and to contribute
is available in the README file from the top source directory of LibPAW.])

# Init Autoconf
AC_PREREQ([2.69])
AC_INIT([LibPAW], [0.5.0],
  [https://gitlab.abinit.org/builders/libpaw/issues],
  [libpaw],
  [https://gitlab.abinit.org/builders/libpaw])
AC_CONFIG_AUX_DIR([config/gnu])
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR([src/m_paw_atom.F90])

# Init Automake (assume GNU Make for now)
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.13 parallel-tests color-tests])
AM_CONFIG_HEADER([config.h])

# Generate timestamp
paw_timestamp=`date '+%Y%m%dT%H%M%S%z'`
AC_SUBST(paw_timestamp)

                    # ------------------------------------ #

#
# LibPAW-specific preprocessing macros
#

AH_BOTTOM([
/* The following was previously provided through libpaw.h */

/*
 * This file is part of the libPAW library.
 * It has to be customized according to the host code.
 * For the time being there are 2 known host codes:
 * ABINIT (www.abinit.org) and BigDFT (bigdft.org).
 */

/*
 * Copyright (C) 2014-2019 Marc Torrent
 *
 * This file is part of the LibPAW software package. For license information,
 * please see the COPYING file in the top-level directory of the LibPAW source
 * distribution.
 */

/* =============================
 * ========= ABINIT ============
 * ============================= */
#if defined HAVE_LIBPAW_ABINIT

/* ABINIT specific macros */
#  include "abi_common.h"

/* Constants and defs */
#  define USE_DEFS use defs_basis

/* MPI wrappers */
#  define USE_MPI_WRAPPERS use m_xmpi

/* Messages, errors */
/* Other macros already defined in abi_common.h */
#  define USE_MSG_HANDLING use m_errors, only : msg_hndl, netcdf_check; use m_abicore
#  undef HAVE_YAML

/* Allocation/deallocation with memory profiling */
#  define USE_MEMORY_PROFILING use m_profiling_abi
/* Use this to allocate/deallocate basic-type arrays with sizes */
#  define LIBPAW_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#  define LIBPAW_DEALLOCATE(ARR) ABI_DEALLOCATE(ARR)
/* Use this to allocate/deallocate basic-type pointers with sizes */
#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#  define LIBPAW_POINTER_DEALLOCATE(ARR) ABI_DEALLOCATE(ARR)
/* Use this to allocate/deallocate user-defined-type arrays with sizes */
#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) ABI_DATATYPE_ALLOCATE(ARR,SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) ABI_DATATYPE_DEALLOCATE(ARR)
/* Use this to allocate basic-type arrays with explicit bounds */
#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) ABI_ALLOCATE(ARR,(BND1))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) ABI_ALLOCATE(ARR,(BND1,BND2))
#  define BOUNDS(LBND,UBND) LBND : UBND

/* =============================
 * ========= BIGDFT ============
 * ============================= */
#elif defined HAVE_LIBPAW_BIGDFT

/* Constants and defs */
#  define USE_DEFS use m_libpaw_defs

/* MPI wrappers */
#  define USE_MPI_WRAPPERS use m_libpaw_mpi

/* Messages, errors */
#  define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#  define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")
/*BigDFT should accept long lines...*/
/*#define MSG_ERROR(msg) call libpaw_msg_hndl(msg,"ERROR","PERS",__FILE__,__LINE__)*/
#  define HAVE_YAML 1

/* Allocation/deallocation with memory profiling */
#  define USE_MEMORY_PROFILING use dynamic_memory
/* Use this to allocate/deallocate basic-type arrays with sizes */
#  define LIBPAW_ALLOCATE(ARR,SIZE) ARR=f_malloc(to_array SIZE )
#  define LIBPAW_DEALLOCATE(ARR) call f_free(ARR)
/* Use this to allocate/deallocate basic-type pointers with sizes */
#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ARR=f_malloc_ptr(to_array SIZE )
#  define LIBPAW_POINTER_DEALLOCATE(ARR) call f_free_ptr(ARR)
/* Use this to allocate/deallocate user-defined-type arrays or pointers with sizes */
#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate user-defined-type arrays with explicit bounds */
#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) ARR=f_malloc((/ BND1 /))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) ARR=f_malloc((/ BND1 , BND2 /))
#  define BOUNDS(LBND,UBND) LBND .to. UBND

/* =============================
 * ========= DEFAULT ===========
 * ============================= */
#else

/* Constants and defs */
#  define USE_DEFS use m_libpaw_defs

/* MPI wrappers */
#  define USE_MPI_WRAPPERS use m_libpaw_mpi

/* Messages, errors */
#  define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#  define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")
#  undef  HAVE_YAML

/* Allocation/deallocation */
#  define USE_MEMORY_PROFILING
/* Use this to allocate/deallocate basic-type arrays with sizes */
#  define LIBPAW_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate/deallocate basic-type pointers with sizes */
#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_POINTER_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate/deallocate user-defined-type arrays with sizes */
#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate user-defined-type arrays with explicit bounds */
#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) allocate(ARR(BND1))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) allocate(ARR(BND1,BND2))
#  define BOUNDS(LBND,UBND) LBND : UBND

/* =============================
 * =========== END =============
 * ============================= */
#endif



/* =============================
 * ===== COMMON DEFINITIONS ====
 * ============================= */

/* F2008 support */
#  define LIBPAW_CONTIGUOUS
#  define LIBPAW_ISO_C_BINDING

/* LibXC support */
#  if defined HAVE_LIBXC
#    define LIBPAW_HAVE_LIBXC HAVE_LIBXC
#  else
#    undef LIBPAW_HAVE_LIBXC
#  endif

/* Netcdf support */
#  if defined HAVE_NETCDF
#    define LIBPAW_HAVE_NETCDF HAVE_NETCDF
#  else
#    undef LIBPAW_HAVE_NETCDF
#  endif

/* FoX support */
#  if defined HAVE_FOX
#    define LIBPAW_HAVE_FOX HAVE_FOX
#  else
#    undef LIBPAW_HAVE_FOX
#  endif

/* Error handlers for netcdf; libpaw_netcdf_check is defined in m_libpaw_tools */
#ifndef NCF_CHECK
#  define NCF_CHECK(ncerr) if (ncerr/=nf90_noerr) call libpaw_netcdf_check(ncerr, "ncf_check")
#endif
])

                    # ------------------------------------ #

#
# System utilities
#

# Check for common programs and features
AC_PROG_AWK
AC_PROG_GREP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_SED

# Workaround for the wrong path to install-sh on Mac systems
AX_PROG_MKDIR_P

                    # ------------------------------------ #

#
# Default settings
#

# Note for developers: you may edit these parameters to tune the behaviour
# of the build system. Please do it with care.

# Optional features
paw_debug_enable_def="no"
paw_multiconfig_enable_def="no"

# Internal variables
LIBPAW_ARCH_SUFFIX=""
LIBPAW_MCFG_MODPFX=""

                    # ------------------------------------ #

#
# Check consistency of defaults
#

# Default triggers for optional features must be yes or no, and not empty
for paw_cfg_default in debug multiconfig; do
  tmp_default=`eval echo \$\{paw_${paw_cfg_default}_enable_def\}`
  if test "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: paw_${paw_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Clean-up
unset tmp_default

                    # ------------------------------------ #

#
# Optional features
#

# Debugging
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
    [Enable verbose function output (default: ${paw_debug_enable_def})])],
  [paw_debug_enable="${enableval}"; paw_debug_type="yon"],
  [paw_debug_enable="${paw_debug_enable_def}"; paw_debug_type="def"])
AC_SUBST(enable_debug)

# Concurrent installation of multiple configurations
AC_ARG_ENABLE([multiconfig],
  [AS_HELP_STRING([--enable-multiconfig],
    [Change library name depending on configuration, a la FFTW, when multiple build flavors are installed concurrently on the system (default: ${paw_multiconfig_enable_def})])],
  [paw_multiconfig_enable="${enableval}"; paw_multiconfig_type="yon"],
  [paw_multiconfig_enable="${paw_multiconfig_enable_def}"; paw_multiconfig_type="def"])
AC_SUBST(enable_multiconfig)

                    # ------------------------------------ #

#
# External dependencies
#

# Init Steredeg core features
SD_CORE_INIT

# Init LibXC
SD_LIBXC_INIT([auto no-fortran])

# Init MPI
SD_MPI_INIT([auto])

                    # ------------------------------------ #

#
# Environment variables
#

# Information about where LibPAW is embedded
AC_ARG_VAR([LIBPAW_EMBED],
  [External program in which to optionally embed LibPAW (abinit or bigdft)])

                    # ------------------------------------ #

#
# Check option consistency
#

# All --enable-* options must be yes or no
for paw_cfg_option in debug; do
  tmp_option=`eval echo \$\{enable_${paw_cfg_option}\}`
  if test "${tmp_option}" != "" -a \
          "${tmp_option}" != "no" -a \
          "${tmp_option}" != "yes"; then
    AC_MSG_ERROR([--enable-${paw_cfg_option} must be "yes" or "no"])
  fi
done
unset tmp_option

                    # ------------------------------------ #

#
# Substitute build-system variables
#

# Defaults
AC_SUBST(paw_debug_enable_def)
AC_SUBST(paw_multiconfig_enable_def)

# Triggers
AC_SUBST(paw_debug_enable)
AC_SUBST(paw_multiconfig_enable)

# Initialization types
AC_SUBST(paw_debug_type)
AC_SUBST(paw_multiconfig_type)

# Internal variables
AC_SUBST(LIBPAW_ARCH_SUFFIX)
AC_SUBST(LIBPAW_MCFG_MODPFX)

                    # ------------------------------------ #

#
# C language support
#

# Init C compiler and preprocessor
if test "${CC}" != "" -a ! -x "${CC}"; then
  paw_cc_probe=`echo "${CC}" | sed -e 's/ .*//'`
  if test ! -x "${paw_cc_probe}"; then
    AC_PATH_PROG([paw_cc_path],[${paw_cc_probe}])
    if test "${paw_cc_path}" = ""; then
      AC_MSG_ERROR([could not run C compiler "${CC}"])
    fi
  fi
fi
AC_PROG_CC
AC_PROG_CPP

                    # ------------------------------------ #

#
# C++ language support
#

# Init C++ compiler
AC_PROG_CXX

                    # ------------------------------------ #

#
# Fortran language support
#

# Init Fortran compiler
if test "${FC}" != "" -a ! -x "${FC}"; then
  paw_fc_probe=`echo "${FC}" | sed -e 's/ .*//'`
  if test ! -x "${paw_fc_probe}"; then
    AC_PATH_PROG([paw_fc_path], [${paw_fc_probe}])
    if test "${paw_fc_path}" = ""; then
      AC_MSG_ERROR([could not run Fortran compiler "${FC}"])
    fi
  fi
fi
AC_PROG_FC

# Enforce strict file extensions
paw_fc_src_ok="unknown"
AC_FC_SRCEXT([F90], [paw_fc_src_ok="yes"], [paw_fc_src_ok="no"])
if test "${paw_fc_src_ok}" != "yes"; then
  AC_MSG_WARN([Fortran file extension could not be changed])
  AC_MSG_WARN([Fortran tests may fail])
fi

# Look for the iso_c_binding module
paw_fc_iso_c_binding="unknown"
AC_MSG_CHECKING([for the ISO_C_BINDING Fortran module])
AC_LANG_PUSH([Fortran])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([], [[use iso_c_binding]])],
  [paw_fc_iso_c_binding="yes"],
  [paw_fc_iso_c_binding="no"])
AC_LANG_POP([Fortran])
AC_MSG_RESULT([${paw_fc_iso_c_binding}])
if test "${paw_fc_iso_c_binding}" != "yes"; then
  AC_MSG_ERROR([the Fortran compiler does not provide the iso_c_binding module.
                This probably happens because this compiler is an antiquity.
                Please use a reasonably recent Fortran compiler.])
fi

# Fortran compiler peculiarities
AX_F90_MODULE_EXTENSION
AX_F90_MODULE_CASE
FMODEXT="${ax_cv_f90_modext}"
AC_SUBST(FMODEXT)

# Inform Automake
AM_CONDITIONAL([FORTRAN_MOD_UPPERCASE],
  [test "${ax_cv_f90_mod_uppercase}" = "yes"])

                    # ------------------------------------ #

#
# Libtool configuration
#

# Init Libtool (must be done once compilers are fully set)
LT_INIT
LT_PREREQ([2.4])
LT_LANG([Fortran])
LTOBJEXT="lo"
AC_SUBST(LTOBJEXT)

                    # ------------------------------------ #

#
# Multicore architecture support
#

# Look for MPI
SD_MPI_DETECT

                    # ------------------------------------ #

#
# External dependencies
#

# LibXC
SD_LIBXC_DETECT

                    # ------------------------------------ #

#
# Concurrent installations support
#

# Set include directory for Fortran modules and
# architecture suffix for library name
if test "${paw_multiconfig_enable}" = "yes"; then
  LIBPAW_MCFG_MODPFX="/libpaw"
  if test "${sd_mpi_enable}" = "yes"; then
    LIBPAW_ARCH_SUFFIX="_mpi"
  fi
fi

                    # ------------------------------------ #

#
# Build parameters for embedded configurations
#

# Check if we are embedded in ABINIT
paw_abinit_common_incs="-I${ac_top_srcdir}/../../common/src/incs -I${ac_top_srcdir}/../../common/src/mods"
if test -s "${ac_top_srcdir}/../common/src/10_defs/defs_basis.F90" -o \
           "${LIBPAW_EMBED}" = "abinit"; then
  FCFLAGS="${FCFLAGS} ${paw_abinit_common_incs}"
fi

                    # ------------------------------------ #

#
# Build-system configuration checkpoint
#

# Write down YAML configuration before any final check
AC_OUTPUT([config/data/libpaw-config.yml])

                    # ------------------------------------ #

#
# Adjust build parameters according to configure options
#

# Debugging support
if test "${paw_debug_enable}" = "yes"; then
  AC_DEFINE([LIBPAW_DEBUG], 1, [Define to 1 to enable internal debugging.])
fi

# FIXME: Remove MPI1/MPI2 distinction
if test "${sd_mpi_ok}" = "yes"; then
  AC_DEFINE([HAVE_MPI2], 1, [Define to 1 if you have modern MPI support.])
fi

# FIXME: Unsupported CPP options
#
# FC_G95
# FC_NAG
# HAVE_AVX_SAFE_MODE
# HAVE_FC_FLUSH
# HAVE_ISO_C_BINDING
# HAVE MPI2_INPLACE
# HAVE_OS_WINDOWS
# HAVE_YAML
# LIBPAW_HAVE_FOX
# LIBPAW_HAVE_LIBXC
# LIBPAW_HAVE_NETCDF
# LIBPAW_ISO_C_BINDING

                    # ------------------------------------ #

# Output configuration
AC_CONFIG_FILES([
  Makefile
  config/data/libpaw.pc
  doc/Makefile
  libxc/Makefile
  src/Makefile])
AC_OUTPUT

# Report configuration
if test -z "${LIBPAW_EMBED}"; then
  tmp_embed_show="no"
else
  tmp_embed_show="${LIBPAW_EMBED}"
fi
if test -z "${LIBPAW_MCFG_MODPFX}"; then
  tmp_moddir_show="${includedir}"
else
  tmp_moddir_show="${includedir}${LIBPAW_MCFG_MODPFX}${LIBPAW_ARCH_SUFFIX}"
fi
AC_MSG_NOTICE([
Final build parameters
----------------------

    TSTAMP       = ${paw_timestamp}

    DEBUG        = ${paw_debug_enable} (init: ${paw_debug_type})
    MULTICONFIG  = ${paw_multiconfig_enable} (init: ${paw_multiconfig_type})

    MPI          = ${sd_mpi_enable} (init: ${sd_mpi_init})

    CPP          = ${CPP}
    CPPFLAGS     = ${CPPFLAGS}
    CC           = ${CC}
    CFLAGS       = ${CFLAGS}
    FC           = ${FC}
    FCFLAGS      = ${FCFLAGS}
    LDFLAGS      = ${LDFLAGS}
    LIBS         = ${LIBS}

    Embedded     : ${tmp_embed_show}
    Module dir   : ${tmp_moddir_show}
    Library name : libpaw${LIBPAW_ARCH_SUFFIX}.(a|so)
])
if test ! -z "${PAW_EMBED}"; then
  AC_MSG_NOTICE([
    *** Please remember to build the ABINIT Common Library before LibPAW ***
])
fi
