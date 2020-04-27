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
#if ( XC_MAJOR_VERSION > 5 ) 
/* ==== libXC v6.0 and later ==== */
 *xc_cst_family_hyb_gga  = -11;
 *xc_cst_family_hyb_mgga = -11;
#else
/* ==== Before libXC v6.0 ==== */
 *xc_cst_family_hyb_gga  = XC_FAMILY_HYB_GGA;
 *xc_cst_family_hyb_mgga = XC_FAMILY_HYB_MGGA;
#endif
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
 * Get the HYBRID constants
 * ===============================================================
 */
void xc_get_hybrid_constants(int *xc_cst_hyb_none,
							 int *xc_cst_hyb_fock,
							 int *xc_cst_hyb_pt2,
							 int *xc_cst_hyb_erf_sr,
							 int *xc_cst_hyb_yukawa_sr,
							 int *xc_cst_hyb_gaussian_sr,
							 int *xc_cst_hyb_semilocal,
							 int *xc_cst_hyb_hybrid,
							 int *xc_cst_hyb_cam,
							 int *xc_cst_hyb_camy,
							 int *xc_cst_hyb_camg,
							 int *xc_cst_hyb_double_hybrid,
							 int *xc_cst_hyb_mixture)
{
#if ( XC_MAJOR_VERSION > 5 ) 
/* ==== libXC v6.0 and later ==== */
 *xc_cst_hyb_none          = XC_HYB_NONE;
 *xc_cst_hyb_fock          = XC_HYB_FOCK;
 *xc_cst_hyb_pt2           = XC_HYB_PT2;
 *xc_cst_hyb_erf_sr        = XC_HYB_ERF_SR;
 *xc_cst_hyb_yukawa_sr     = XC_HYB_YUKAWA_SR;
 *xc_cst_hyb_gaussian_sr   = XC_HYB_GAUSSIAN_SR;
 *xc_cst_hyb_semilocal     = XC_HYB_SEMILOCAL;
 *xc_cst_hyb_hybrid        = XC_HYB_HYBRID;
 *xc_cst_hyb_cam           = XC_HYB_CAM;
 *xc_cst_hyb_camy          = XC_HYB_CAMY;
 *xc_cst_hyb_camg          = XC_HYB_CAMG;
 *xc_cst_hyb_double_hybrid = XC_HYB_DOUBLE_HYBRID;
 *xc_cst_hyb_mixture       = XC_HYB_MIXTURE;
#else
/* ==== Before libXC v6.0 ==== */
 *xc_cst_hyb_none      = -11; *xc_cst_hyb_fock          = -11;
 *xc_cst_hyb_pt2       = -11; *xc_cst_hyb_erf_sr        = -11;
 *xc_cst_hyb_yukawa_sr = -11; *xc_cst_hyb_gaussian_sr   = -11;
 *xc_cst_hyb_semilocal = -11; *xc_cst_hyb_hybrid        = -11;
 *xc_cst_hyb_cam       = -11; *xc_cst_hyb_camy          = -11;
 *xc_cst_hyb_camg      = -11; *xc_cst_hyb_double_hybrid = -11;
 *xc_cst_hyb_mixture   = -11;
#endif
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
 * Wrappers to the LDA/GGA/MGGA functionals
 * ===============================================================
 */
/* ---------------------------------------------------------------
   ----- LDA ----- */
void xc_get_lda(const XC(func_type) *xc_func, int np, const double *rho, 
        double *zk, double *vrho, double *v2rho2, double *v3rho3)
#if ( XC_MAJOR_VERSION > 4 ) 
/* ==== libXC v5.0 and later ==== */
 {xc_lda(xc_func, np, rho, zk, vrho, v2rho2, v3rho3, NULL);}
#else
/* ==== Before libXC v5.0 ==== */
 {xc_lda(xc_func, np, rho, zk, vrho, v2rho2, v3rho3);}
#endif
/* ---------------------------------------------------------------
   ----- GGA ----- */
void xc_get_gga(const XC(func_type) *xc_func, int np,
        const double *rho, const double *sigma,
        double *zk, double *vrho, double *vsigma,
        double *v2rho2, double *v2rhosigma, double *v2sigma2,
        double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
#if ( XC_MAJOR_VERSION > 4 ) 
/* ==== libXC v5.0 and later ==== */
 {xc_gga(xc_func, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
         NULL, NULL, NULL, NULL, NULL);}
#else
/* ==== Before libXC v5.0 ==== */
 {xc_gga(xc_func, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2,
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);}
#endif
/* ---------------------------------------------------------------
   ----- meta-GGA ----- */
void xc_get_mgga(const XC(func_type) *xc_func, int np,
        const double *rho, const double *sigma, const double *lapl, const double *tau,
        double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
        double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
        double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
        double *v2lapltau, double *v2tau2) 
#if ( XC_MAJOR_VERSION > 4 ) 
/* ==== libXC v5.0 and later ==== */
 {xc_mgga(xc_func, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2,
          v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);}
#else
/* ==== Before libXC v5.0 ==== */
 {xc_mgga(xc_func, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, 
		  v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau,
		  v2sigmalapl, v2sigmatau, v2lapltau);}
#endif

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
#if ( XC_MAJOR_VERSION > 4 ) 
/* ==== libXC v5.0 and later ==== */
 {if (n_ext_params == xc_func->info->ext_params.n)
   {XC(func_set_ext_params)(xc_func, ext_params);}
#elif ( XC_MAJOR_VERSION > 3 ) 
/* ==== libXC v4.0 ==== */
 {if (xc_func->info->number == XC_HYB_GGA_XC_PBEH && n_ext_params == 1)
   /* set_ext_params function is missing for PBE0 */  
   {xc_func->cam_alpha=ext_params[0];xc_func->mix_coef[0]=1.0-ext_params[0];}
  else if (xc_func->info->number == XC_MGGA_X_TB09 && n_ext_params >= 1)
   /* XC_MGGA_X_TB09 has only one parameter */
   {XC(func_set_ext_params)(xc_func, ext_params);}
  else if (n_ext_params == xc_func->info->n_ext_params)
   {XC(func_set_ext_params)(xc_func, ext_params);}

#else
/* ==== Before libXC v4.0 ==== */
 {if (xc_func->info->number == XC_LDA_C_XALPHA && n_ext_params == 1)
   {XC(lda_c_xalpha_set_params)(xc_func, *ext_params);}
  else if (xc_func->info->number == XC_MGGA_X_TB09 && n_ext_params >= 1)
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

/* ===============================================================
 * Missing function:
 *  Return 1 if the functional is hybrid, from its id
 * ===============================================================
 */
int xc_func_is_hybrid_from_id(int func_id)
#if ( XC_MAJOR_VERSION > 5 ) 
/* ==== libXC v6.0 and later ==== */
 {xc_func_type func; int result=0;
  if(xc_func_init(&func,func_id,XC_UNPOLARIZED)==0)
    {if (func.hyb_number_terms>0) {result=1;}}
  xc_func_end(&func);
  return result;
 }
#else
/* ==== Before libXC v6.0 ==== */
 {int family; family=xc_family_from_id(func_id, NULL, NULL);
  if (family==XC_FAMILY_HYB_GGA || family==XC_FAMILY_HYB_MGGA)
   {return 1;}
  else
   {return 0;}
 }
#endif

#endif
