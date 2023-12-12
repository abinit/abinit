!!****m* ABINIT/m_manage_kokkos
!! NAME
!!  m_manage_kokkos
!!
!! FUNCTION
!! This module provides iso_c_binding wrappers to kernels written using kokkos.
!! computational routines located in m_xg
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2022 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_manage_kokkos

 use, intrinsic :: iso_c_binding

 implicit none

 interface

   ! ========================================================================
   ! ========================================================================

   subroutine opernlc_ylm_allwf_kokkos(cplex, cplex_enl, cplex_fac, &
     &                                 dimenl1, dimenl2, dimekbq, &
     &                                 iatm, itypat, ntypat, nprojs, &
     &                                 natom, nincat, nspinor, &
     &                                 nspinortot, paw_opt, &
     &                                 nlmn, lmnmax, &
     &                                 enl_gpu, &
     &                                 gx_gpu, &
     &                                 gxfac_gpu, &
     &                                 gxfac2_gpu, &
     &                                 gxfac_sij_gpu, &
     &                                 shift_spinor, ndat, &
     &                                 atindx1_gpu, &
     &                                 indlmn_gpu, &
     &                                 lambda_gpu, &
     &                                 sij_typ_gpu, &
     &                                 shift_proj, &
     &                                 nattyp_max) &
     & bind(c, name='opernlc_ylm_allwf_kokkos_cpp')
     use, intrinsic :: iso_c_binding
     implicit none
     integer(kind=c_int32_t), value, intent(in)    :: cplex, cplex_enl, cplex_fac
     integer(kind=c_int32_t), value, intent(in)    :: dimenl1, dimenl2, dimekbq
     integer(kind=c_int32_t), value, intent(in)    :: iatm, itypat, ntypat, nprojs
     integer(kind=c_int32_t), value, intent(in)    :: natom, nincat, nspinor
     integer(kind=c_int32_t), value, intent(in)    :: nspinortot, paw_opt
     integer(kind=c_int32_t), value, intent(in)    :: nlmn, lmnmax
     type(c_ptr),             value                :: enl_gpu ! (dimenl1, dimenl2, nspinortot**2, dimekbq)
     type(c_ptr),             value                :: gx_gpu ! (cplex,nlmn,nincat,nspinor*ndat)
     type(c_ptr),             value                :: gxfac_gpu ! (cplex_fac,nlmn,nincat,nspinor*ndat)
     type(c_ptr),             value                :: gxfac2_gpu ! (cplex_fac,nlmn,nincat,nspinor*ndat)
     type(c_ptr),             value                :: gxfac_sij_gpu !(cplex,nlmn,nincat,nspinor*ndat*(paw_opt/3))
     integer(kind=c_int32_t), value, intent(in)    :: shift_spinor, ndat
     type(c_ptr),             value                :: atindx1_gpu ! (natom)
     type(c_ptr),             value                :: indlmn_gpu  ! (6,nlmn)
     type(c_ptr),             value                :: lambda_gpu ! (ndat)
     type(c_ptr),             value                :: sij_typ_gpu ! (((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
     integer(kind=c_int32_t), value, intent(in)    :: shift_proj, nattyp_max
   end subroutine opernlc_ylm_allwf_kokkos

   ! ========================================================================
   ! ========================================================================

   !> add arrays on GPU, array already on device (managed memory)
   subroutine add_array_kokkos(array1_ptr, array2_ptr, array_size) &
     & bind(c, name='add_array_kokkos_cpp')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr)            , value             :: array1_ptr
     type(c_ptr)            , value             :: array2_ptr
     integer(kind=c_int32_t), value, intent(in) :: array_size
   end subroutine add_array_kokkos

   ! ========================================================================
   ! ========================================================================

   !> assemble energy contributions into ghc / gsc array
   subroutine assemble_energy_contribution_kokkos(ghc_ptr, &
     &                                            gsc_ptr, &
     &                                            kinpw_k2_ptr, &
     &                                            cwavef_ptr, &
     &                                            gvnlxc_ptr, &
     &                                            ndat, &
     &                                            my_nspinor, &
     &                                            npw_k2, &
     &                                            sij_opt, &
     &                                            k1_eq_k2, &
     &                                            hugevalue) &
     & bind(c, name='assemble_energy_contribution_kokkos_cpp')
     use, intrinsic :: iso_c_binding
     implicit none
     type(c_ptr)            , value             :: ghc_ptr
     type(c_ptr)            , value             :: gsc_ptr
     type(c_ptr)            , value             :: kinpw_k2_ptr
     type(c_ptr)            , value             :: cwavef_ptr
     type(c_ptr)            , value             :: gvnlxc_ptr
     integer(kind=c_int32_t), value, intent(in) :: ndat
     integer(kind=c_int32_t), value, intent(in) :: my_nspinor
     integer(kind=c_int32_t), value, intent(in) :: npw_k2
     integer(kind=c_int32_t), value, intent(in) :: sij_opt
     logical(kind=c_bool),    value, intent(in) :: k1_eq_k2
     real(kind=c_double),     value, intent(in) :: hugevalue
   end subroutine assemble_energy_contribution_kokkos


 end interface

contains
  !!***

end module m_manage_kokkos
!!***

