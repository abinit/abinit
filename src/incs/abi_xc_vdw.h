/* abi_xc_vdw.h */

/*
 * Copyright (C) 2012-2020 ABINIT Group (Yann Pouillon)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef _ABINIT_XC_VDW_H
#define _ABINIT_XC_VDW_H

#include "abi_common.h"

/* Error handler for VDWXC NetCDF calls */
#ifdef HAVE_FC_LONG_LINES
#define NETCDF_VDWXC_CHECK(ncerr) \
  if ( ncerr /= NF90_NOERR ) then NEWLINE \
    call vdw_df_netcdf_ioerr(ncerr,__FILE__,__LINE__) NEWLINE \
  end if
#else
#define NETCDF_VDWXC_CHECK(ncerr) \
  if ( ncerr /= NF90_NOERR ) then NEWLINE \
    call vdw_df_netcdf_ioerr(ncerr) NEWLINE \
  end if
#endif

/* Debug messages */
#ifdef DEBUG_VERBOSE

#  ifdef HAVE_FC_LONG_LINES
#    define VDWXC_DBG_ENTER(mode,funcname) \
       call vdw_df_write_func(ABI_FUNC,mode) NEWLINE \
       call sentinel(1,mode,__FILE__,funcname,__LINE__)
#    define VDWXC_DBG_EXIT(mode,funcname)  \
       call vdw_df_write_func(ABI_FUNC,mode) NEWLINE \
       call sentinel(2,mode,__FILE__,funcname,__LINE__)
#  else
#    define VDWXC_DBG_ENTER(mode,funcname) \
       call vdw_df_write_func(ABI_FUNC,mode) NEWLINE \
       call sentinel(1,mode)
#    define VDWXC_DBG_EXIT(mode,funcname) \
       call vdw_df_write_func(ABI_FUNC,mode) NEWLINE \
       call sentinel(2,mode)
#  endif

#else

#  define VDWXC_DBG_ENTER(mode,funcname)
#  define VDWXC_DBG_EXIT(mode,funcname)

#endif

#endif /* _ABINIT_XC_VDW_H */
