/* libpaw.h */

/*
 * This file is part of the libPAW library.
 * It has to be customized according to the host code.
 * For the time being there are 2 known host codes:
 * ABINIT (www.abinit.org) and BigDFT (bigdft.org).
 */

/*
 * Copyright (C) 2014-2020 ABINIT Group (MT)
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 */

/* config.h should contain all preprocessing directives */
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

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
#  undef  HAVE_YAML

/* Allocation/deallocation with memory profiling */
#  define USE_MEMORY_PROFILING use m_profiling_abi
/* Use this to allocate/deallocate basic-type arrays with sizes */
#  define LIBPAW_ALLOCATE(ARR,SIZE) ABI_MALLOC(ARR,SIZE)
#  define LIBPAW_DEALLOCATE(ARR) ABI_FREE(ARR)
/* Use this to allocate/deallocate basic-type pointers with sizes */
#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ABI_MALLOC(ARR,SIZE)
#  define LIBPAW_POINTER_DEALLOCATE(ARR) ABI_FREE(ARR)
/* Use this to allocate/deallocate user-defined-type arrays with sizes */
#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) ABI_MALLOC(ARR,SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) ABI_FREE(ARR)
/* Use this to allocate basic-type arrays with explicit bounds */
#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) ABI_MALLOC(ARR,(BND1))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) ABI_MALLOC(ARR,(BND1,BND2))
#  define BOUNDS(LBND,UBND) LBND : UBND

/* libXC support */
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
#  undef LIBPAW_HAVE_FOX

/* F2008 support */
#  define LIBPAW_CONTIGUOUS ABI_CONTIGUOUS
#  define LIBPAW_ISO_C_BINDING 1

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
#  define ABI_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define ABI_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define ABI_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define ABI_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")
/*BigDFT should accept long lines...*/
/*#define ABI_ERROR(msg) call libpaw_msg_hndl(msg,"ERROR","PERS",__FILE__,__LINE__)*/
#  define HAVE_YAML

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

/* libXC support */
#  if (defined HAVE_LIBXC)
#    define LIBPAW_HAVE_LIBXC HAVE_LIBXC
#  else
#    undef LIBPAW_HAVE_LIBXC
#  endif

/* Netcdf support */
#  undef LIBPAW_HAVE_NETCDF

/* FoX support */
#  undef LIBPAW_HAVE_FOX

/* F2008 support */
#  define LIBPAW_CONTIGUOUS
#  define LIBPAW_ISO_C_BINDING 1

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
#  define ABI_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define ABI_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define ABI_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define ABI_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")
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

/* libXC support */
#  undef LIBPAW_HAVE_LIBXC

/* Netcdf support */
#  undef LIBPAW_HAVE_NETCDF

/* FoX support */
#  undef LIBPAW_HAVE_FOX

/* F2008 support */
#  define LIBPAW_CONTIGUOUS
#  define LIBPAW_ISO_C_BINDING 1

/* =============================
 * =========== END =============
 * ============================= */
#endif



/* =============================
 * ===== COMMON DEFINITIONS ====
 * ============================= */

/* Error handlers for netcdf; libpaw_netcdf_check is defined in m_libpaw_tools */
#ifndef NCF_CHECK
#  define NCF_CHECK(ncerr) if (ncerr/=nf90_noerr) call libpaw_netcdf_check(ncerr, "ncf_check")
#endif
