!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_wvltypes
!! NAME
!! defs_wvltypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatypes for the
!! wavelet part of the ABINIT package.
!!
!! List of datatypes :
!! * wvl_projectors_type : structure to store projectors for wavelets calculations.
!! * wvl_wf_type : structure to store wavefunctions for wavelets calculations.
!! * wvl_data : container for all required wavelets data.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2016 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_wvltypes

 use defs_basis
 use m_profiling_abi

#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : atoms_data, locreg_descriptors, local_zone_descriptors, &
      & DFT_PSP_projectors, GPU_pointers, rho_descriptors, &
      & DFT_local_fields, DFT_wavefunction, energy_terms, &
      & rholoc_objects,gaussian_basis,paw_objects,xc_info,coulomb_operator
#endif

 implicit none
!!***

 !!****t* defs_wvltypes/wvl_projectors_type
 !! NAME
 !! wvl_projectors_type
 !!
 !! FUNCTION
 !! This type constains all physical data associated to
 !! projectors of a given system (Z and atomic positions).
 !!
 !! SOURCE
 type wvl_projectors_type

#if defined HAVE_DFT_BIGDFT
    type(DFT_PSP_projectors) :: nlpsp
    !object containing the projectors as a sum of Gaussian functions
    type(gaussian_basis),dimension(:),allocatable :: G
#else
    integer :: nlpsp
#endif

 end type wvl_projectors_type
!!***

 !!****t* defs_wvltypes/wvl_wf_type
 !! NAME
 !! wvl_wf_type
 !!
 !! FUNCTION
 !! This type constains all physical data associated to
 !! wavefunctions of a given system (Z and atomic positions).
 !!
 !! SOURCE
 type wvl_wf_type
#if defined HAVE_DFT_BIGDFT
    type(DFT_wavefunction) :: ks

    ! Some GPU internals.
    type(GPU_pointers) :: GPU
#else
    ! To avoid having an empty type.
    integer :: ks
#endif
 end type wvl_wf_type
!!***

!!****t* defs_wvltypes/wvl_internal_type
!! NAME
!! wvl_internal_type
!!
!! FUNCTION
!! This type is a gathering for all internal variables wavelets required. It is
!! included in the datatypes strutcture.
!!
!! NOTES
!! If you modify this datatype, the copy, init/free routines are in wvl_utils.F90.
!!
!! SOURCE
type wvl_internal_type
  real(dp) :: h(3)
  ! The hgrid values in each direction, after the application of the
  ! boundary conditions. In free boundary conditions, the three values are equal.
  real(dp) :: shift(3)
  ! Shift applied by BigDFT on the atomic position (in cartesian coordinates).

  integer,allocatable :: npspcode_paw_init_guess(:)
  ! This is for the initial guess in PAW pseudos, this 
  ! is only relevant for ABINIT and PAW

#if defined HAVE_DFT_BIGDFT
  type(locreg_descriptors) :: Glr
  ! Contains the description of the global localisation region.

  type(atoms_data) :: atoms
  ! A copy of the current dtset values.

  type(rholoc_objects) :: rholoc
  ! Information for local density

  type(paw_objects) :: paw
  ! Objects for PAW 
#endif

  character(len=4) :: exctxpar
  ! parallelisation scheme of the exact exchange operator
  !   BC (Blocking Collective)
  !   OP2P (Overlap Point-to-Point)  
  real(dp) :: exctxfac
  ! Exact exchange coefficient to mix.  
end type wvl_internal_type
!!***

!!****t* defs_wvltypes/wvl_denspot_type
!! NAME
!! wvl_denspot_type
!!
!! FUNCTION
!! This type contains the 
!! BigDFT information for 
!! Densities and potentials
!!
!! SOURCE
type wvl_denspot_type
#if defined HAVE_DFT_BIGDFT
   ! Densities and potentials, and related metadata, needed for their creation/application
   type(DFT_local_fields) :: denspot
#endif
   ! This is a pointer on the symmetry object present in atoms structure of BigDFT.
   integer :: symObj
end type wvl_denspot_type
!!***

!!****t* defs_wvltypes/wvl_energy_terms
!! NAME
!! wvl_energy_terms
!!
!! FUNCTION
!! This type contains the 
!! BigDFT information for energies
!!
!! SOURCE
type wvl_energy_terms
#if defined HAVE_DFT_BIGDFT
   !object for energies
   type(energy_terms) :: energs
#else
   integer :: energs
#endif
end type wvl_energy_terms
!!***

#ifndef HAVE_DFT_BIGDFT
!!****t* defs_wvltypes/coulomb_operator
!! NAME
!! coulomb_operator
!!
!! FUNCTION
!! This type contains the 
!! BigDFT information for kernels
!!
!! SOURCE
type coulomb_operator
   !this is a fake type when not compiled with BigDFT.
   integer :: co
end type coulomb_operator
!!***
#endif


!!****t* defs_wvltypes/wvl_data
!! NAME
!! wvl_data
!!
!! FUNCTION
!! This type is a container to limit the number of arguments in
!! ABINIT, it should not contains attributes other than types
!! created for wavelets.
!!
!! SOURCE
type wvl_data
   ! The data associated to projectors (computed from pseudo)
   type(wvl_projectors_type) :: projectors
   ! The data associated to the wavefunctions
   type(wvl_wf_type) :: wfs
   ! The pointers to internal BigDFT data structures
   type(wvl_internal_type) :: descr
   ! Densities and potentials
   type(wvl_denspot_type) :: den
   ! Energies
   type(wvl_energy_terms) :: e
end type wvl_data
!!***
end module defs_wvltypes
