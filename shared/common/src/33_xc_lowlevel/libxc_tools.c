/*
 * Copyright (C) 2015-2020 ABINIT group (MT)
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>

#if defined HAVE_LIBXC

#include "xc.h"
#include "xc_version.h"
#include "xc_funcs.h"

/* if version before 4 get config file */
#if ( XC_MAJOR_VERSION < 4 )
#  include "xc_config.h"
#else
#  define FLOAT double
#endif

/* ===============================================================
 * Get the SINGLE_PRECISION constant
 * ===============================================================
 */
void xc_get_singleprecision_constant(int *xc_cst_singleprecision)
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
void xc_get_family_constants(int *xc_cst_family_unknown,
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
void xc_get_flags_constants(int *xc_cst_flags_have_exc,
                            int *xc_cst_flags_have_vxc,
                            int *xc_cst_flags_have_fxc,
                            int *xc_cst_flags_have_kxc,
                            int *xc_cst_flags_have_lxc,
                            int *xc_cst_flags_needs_laplacian)
{
 *xc_cst_flags_have_exc  = XC_FLAGS_HAVE_EXC;
 *xc_cst_flags_have_vxc  = XC_FLAGS_HAVE_VXC;
 *xc_cst_flags_have_fxc  = XC_FLAGS_HAVE_FXC;
 *xc_cst_flags_have_kxc  = XC_FLAGS_HAVE_KXC;
 *xc_cst_flags_have_lxc  = XC_FLAGS_HAVE_LXC;
#if ( XC_MAJOR_VERSION > 3 )
 *xc_cst_flags_needs_laplacian  = XC_FLAGS_NEEDS_LAPLACIAN;
#else
 *xc_cst_flags_needs_laplacian  = 1;
#endif
}

/* ===============================================================
 * Get the KIND constants
 * ===============================================================
 */
void xc_get_kind_constants(int *xc_cst_exchange,
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
XC(func_type) * xc_func_type_malloc()
 {return (XC(func_type) *) malloc(sizeof(XC(func_type)));}

void xc_func_type_free(XC(func_type) **xc_func)
 {free(*xc_func);*xc_func = NULL;}

/* ===============================================================
 * Get properties from a xc_func_info_type pointer
 *     These accessors where not provided before libXC v3
 * ===============================================================
 */
#if ( XC_MAJOR_VERSION > 3 )
/* ==== libXC v4.0 and later ==== */
char const *xc_get_info_name(XC(func_type) *xc_func)
 {return xc_func_info_get_name(xc_func->info);}
int xc_get_info_flags(XC(func_type) *xc_func)
 {return xc_func_info_get_flags(xc_func->info);}
int xc_get_info_kind(XC(func_type) *xc_func)
 {return xc_func_info_get_kind(xc_func->info);}
char const *xc_get_info_refs(XC(func_type) *xc_func, const int *number)
 {if (*number>=0&&*number<XC_MAX_REFERENCES)
   {if (xc_func_info_get_references(xc_func->info,*number) != NULL)
    {return xc_func_info_get_references(xc_func->info,*number)->ref;}}
  else {return NULL;}
  return NULL;}

#elif ( XC_MAJOR_VERSION > 2 )
/* ==== libXC v3.0 ==== */
char const *xc_get_info_name(XC(func_type) *xc_func)
 {return xc_func_info_get_name(xc_func->info);}
int xc_get_info_flags(XC(func_type) *xc_func)
 {return xc_func_info_get_flags(xc_func->info);}
int xc_get_info_kind(XC(func_type) *xc_func)
 {return xc_func_info_get_kind(xc_func->info);}
char const *xc_get_info_refs(XC(func_type) *xc_func, const int *number)
 {if (*number>=0&&*number<=4)
   {if (xc_func_info_get_ref(xc_func->info,*number) != NULL)
    {return xc_func_info_get_ref(xc_func->info,*number);}}
  else {return NULL;}
  return NULL;}

#else
/* ==== Before libXC v3.0 ==== */
char const *xc_get_info_name(XC(func_type) *xc_func)
 {return xc_func->info->name;}
int xc_get_info_flags(XC(func_type) *xc_func)
 {return xc_func->info->flags;}
int xc_get_info_kind(XC(func_type) *xc_func)
 {return xc_func->info->kind;}
char const *xc_get_info_refs(XC(func_type) *xc_func, const int *number)
 {if (*number==0) {return xc_func->info->refs;} else {return NULL;}
  return NULL;}
#endif

/* ===============================================================
 * Wrapper to xc_func_set_ext_params for backward compatibility
 *    Allows to change the parameters of a XC functional
 * ===============================================================
 */
void xc_func_set_params(XC(func_type) *xc_func, double *ext_params, int n_ext_params)
#if ( XC_MAJOR_VERSION > 3 ) 
/* ==== libXC v4.0 and later ==== */
 {/* set_ext_params function is missing for PBE0 */  
  if (xc_func->info->number == XC_HYB_GGA_XC_PBEH && n_ext_params == 1)
   {xc_func->cam_alpha=ext_params[0];xc_func->mix_coef[0]=1.0-ext_params[0];}

  else if (n_ext_params == xc_func->info->n_ext_params)
   {XC(func_set_ext_params)(xc_func, ext_params);}

#else
/* ==== Before libXC v4.0 ==== */
 {if (xc_func->info->number == XC_LDA_C_XALPHA && n_ext_params == 1)
   {XC(lda_c_xalpha_set_params)(xc_func, *ext_params);}
  else if (xc_func->info->number == XC_MGGA_X_TB09 && n_ext_params == 1)
   {XC(mgga_x_tb09_set_params)(xc_func, *ext_params);}
#if ( XC_MAJOR_VERSION > 2 || ( XC_MAJOR_VERSION > 1 && XC_MINOR_VERSION > 0 ) ) 
  else if (xc_func->info->number == XC_HYB_GGA_XC_PBEH && n_ext_params == 1)
   {XC(hyb_gga_xc_pbeh_set_params)(xc_func, *ext_params);}
  else if (xc_func->info->number == XC_HYB_GGA_XC_HSE03 && n_ext_params == 3)
   {XC(hyb_gga_xc_hse_set_params)(xc_func, ext_params[0], ext_params[2]);
    xc_func->cam_omega=ext_params[1];}
  else if (xc_func->info->number == XC_HYB_GGA_XC_HSE06 && n_ext_params == 3)
   {XC(hyb_gga_xc_hse_set_params)(xc_func, ext_params[0], ext_params[2]);
    xc_func->cam_omega=ext_params[1];}
#else
  else if (xc_func->info->number == XC_HYB_GGA_XC_HSE03 && n_ext_params == 3)
   {XC(hyb_gga_xc_hse_set_params)(xc_func, ext_params[2]);
    xc_func->cam_omega=ext_params[1];}
  else if (xc_func->info->number == XC_HYB_GGA_XC_HSE06 && n_ext_params == 3)
   {XC(hyb_gga_xc_hse_set_params)(xc_func, ext_params[2]);
    xc_func->cam_omega=ext_params[1];}
#endif
#endif
  else
   {fprintf(stderr, "BUG: invalid entry in set_params!\n");abort();}
 }

/* ===============================================================
 * Wrapper to xc_func_set_dens_threshold for backward compatibility
 *    Allows to change the zero-density threshold of a XC functional
 *    Only available from libXC v4
 * ===============================================================
 */
void xc_func_set_density_threshold(XC(func_type) *xc_func, double *dens_threshold)
#if ( XC_MAJOR_VERSION > 3 ) 
/* ==== libXC v4.0 and later ==== */
   {XC(func_set_dens_threshold)(xc_func, *dens_threshold);}
#else
   {fprintf(stderr, "WARNING: setting density threshold not available for libXC<4.0!\n");}
#endif


#endif
