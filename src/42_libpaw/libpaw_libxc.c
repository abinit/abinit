/*
 * Copyright (C) 2015-2018 ABINIT group (MT)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

/* ===============================================================
 * Set of C functions interfacing the LibXC library.
 * (see http://www.tddft.org/programs/Libxc)
 * ===============================================================
 */

#include "libpaw.h"
#include <stdlib.h>

#if defined LIBPAW_HAVE_LIBXC

#include "xc.h"
#include "xc_version.h"
/* if version before 4 get config file*/
#if ( XC_MAJOR_VERSION < 4 )
#include "xc_config.h"
#else
#  define FLOAT double
#endif

/* ===============================================================
 * Get the SINGLE_PRECISION constant
 * ===============================================================
 */
void libpaw_xc_get_singleprecision_constant(int *xc_cst_singleprecision)
{
 if (sizeof(FLOAT)<sizeof(double))
  {*xc_cst_singleprecision = 1;}
 else
  {*xc_cst_singleprecision = 0;}
}

/* ===============================================================
 * Get the FAMILY constants
 * ===============================================================
 */
void libpaw_xc_get_family_constants(int *xc_cst_family_unknown,
                                   int *xc_cst_family_lda,
                                   int *xc_cst_family_gga,
                                   int *xc_cst_family_mgga,
                                   int *xc_cst_family_lca,
                                   int *xc_cst_family_oep,
                                   int *xc_cst_family_hyb_gga,
                                   int *xc_cst_family_hyb_mgga)
{
 *xc_cst_family_unknown  = XC_FAMILY_UNKNOWN;
 *xc_cst_family_lda      = XC_FAMILY_LDA;
 *xc_cst_family_gga      = XC_FAMILY_GGA;
 *xc_cst_family_mgga     = XC_FAMILY_MGGA;
 *xc_cst_family_lca      = XC_FAMILY_LCA;
 *xc_cst_family_oep      = XC_FAMILY_OEP;
 *xc_cst_family_hyb_gga  = XC_FAMILY_HYB_GGA;
 *xc_cst_family_hyb_mgga = XC_FAMILY_HYB_MGGA;
}

/* ===============================================================
 * Get the FLAGS constants
 * ===============================================================
 */
void libpaw_xc_get_flags_constants(int *xc_cst_flags_have_exc,
                                   int *xc_cst_flags_have_vxc,
                                   int *xc_cst_flags_have_fxc,
                                   int *xc_cst_flags_have_kxc,
                                   int *xc_cst_flags_have_lxc)
{
 *xc_cst_flags_have_exc  = XC_FLAGS_HAVE_EXC;
 *xc_cst_flags_have_vxc  = XC_FLAGS_HAVE_VXC;
 *xc_cst_flags_have_fxc  = XC_FLAGS_HAVE_FXC;
 *xc_cst_flags_have_kxc  = XC_FLAGS_HAVE_KXC;
 *xc_cst_flags_have_lxc  = XC_FLAGS_HAVE_LXC;
}

/* ===============================================================
 * Get the KIND constants
 * ===============================================================
 */
void libpaw_xc_get_kind_constants(int *xc_cst_exchange,
                                  int *xc_cst_correlation,
                                  int *xc_cst_exchange_correlation,
                                  int *xc_cst_kinetic)
{
 *xc_cst_exchange              = XC_EXCHANGE;
 *xc_cst_correlation           = XC_CORRELATION;
 *xc_cst_exchange_correlation  = XC_EXCHANGE_CORRELATION;
 *xc_cst_kinetic               = XC_KINETIC;
}

/* ===============================================================
 * Allocate/free xc_func_type pointer
 * ===============================================================
 */
XC(func_type) * libpaw_xc_func_type_malloc()
 {return (XC(func_type) *) malloc(sizeof(XC(func_type)));}

void libpaw_xc_func_type_free(XC(func_type) **xc_func)
 {free(*xc_func);*xc_func = NULL;}

/* ===============================================================
 * Get properties from a xc_func_info_type pointer
 *     These accessors where not provided before libXC v3
 * ===============================================================
 */
#if defined XC_MICRO_VERSION
 /* libXC v3.0 and later */
char const *libpaw_xc_get_info_name(XC(func_type) *xc_func)
 {return xc_func_info_get_name(xc_func->info);}
int libpaw_xc_get_info_flags(XC(func_type) *xc_func)
 {return xc_func_info_get_flags(xc_func->info);}
int libpaw_xc_get_info_kind(XC(func_type) *xc_func)
 {return xc_func_info_get_kind(xc_func->info);}
char const *libpaw_xc_get_info_refs(XC(func_type) *xc_func, const int *number)
 {if (*number>=0&&*number<=4)
#if ( XC_MAJOR_VERSION < 4 ) 
   {if (xc_func_info_get_ref(xc_func->info,*number) != NULL)
    {return xc_func_info_get_ref(xc_func->info,*number);}}
#else
   {if (xc_func_info_get_references(xc_func->info,*number) != NULL)
    {return xc_func_info_get_references(xc_func->info,*number);}}
#endif
  else {return NULL;}}
#else
 /* libXC before v3.0 */
char const *libpaw_xc_get_info_name(XC(func_type) *xc_func)
 {return xc_func->info->name;}
int libpaw_xc_get_info_flags(XC(func_type) *xc_func)
 {return xc_func->info->flags;}
int libpaw_xc_get_info_kind(XC(func_type) *xc_func)
 {return xc_func->info->kind;}
char const *libpaw_xc_get_info_refs(XC(func_type) *xc_func, const int *number)
 {if (*number==0) {return xc_func->info->refs;} else {return NULL;}}
#endif

#endif
