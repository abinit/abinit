!{\src2tex{textfont=tt}}
!!****f* ABINIT/nullify_wvl_data
!! NAME
!!  nullify_wvl_data
!!
!! FUNCTION
!!  Nullify all wvl pointers
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2016 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      nullify_atoms_data,nullify_dft_local_fields,nullify_diis_objects
!!      nullify_gaussian_basis,nullify_gpu_pointers
!!      nullify_local_zone_descriptors,nullify_locreg_descriptors
!!      nullify_orbitals_data,nullify_paw_objects,nullify_rholoc_objects
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nullify_wvl_data(wvl)

 use m_profiling_abi
 use m_errors
 use defs_basis
 use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API, only : &
& nullify_gaussian_basis,nullify_local_zone_descriptors,&
& nullify_orbitals_data,comms_cubic,comms_linear_null, &
& nullify_diis_objects,p2pComms_null,&
& nullify_GPU_pointers,nullify_locreg_descriptors,&
& nullify_atoms_data,nullify_rholoc_objects,nullify_paw_objects,&
& nullify_DFT_local_fields, DFT_PSP_projectors_null
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_wvl_data'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wvl_data) , intent(inout)  :: wvl

!Local variables-------------------------------
!character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
#if defined HAVE_BIGDFT

 DBG_ENTER("COLL")

!1)   wvl_projectors_type: projectors
 wvl%projectors%nlpsp=DFT_PSP_projectors_null()

!2)   wvl_wf_type: wfs
!2.1) DFT_wavefunction: ks
 nullify(wvl%wfs%ks%psi)
 nullify(wvl%wfs%ks%hpsi)
 nullify(wvl%wfs%ks%psit)
 nullify(wvl%wfs%ks%spsi)
 nullify(wvl%wfs%ks%gaucoeffs)
 nullify(wvl%wfs%ks%confdatarr)
 nullify(wvl%wfs%ks%oldpsis)
 nullify(wvl%wfs%ks%coeff)
!2.1.1) gaussian_basis: gbd 
 call nullify_gaussian_basis(wvl%wfs%ks%gbd)
!2.1.2) local_zone_descriptors: Lzd
 call nullify_local_zone_descriptors(wvl%wfs%ks%Lzd)
!2.1.3) orbitals_data: orbs
 call nullify_orbitals_data(wvl%wfs%ks%orbs)
!2.1.4) communications_arrays: comms
 wvl%wfs%ks%comms=comms_cubic_null()
!2.1.5) diis_objects: diis
 call nullify_diis_objects(wvl%wfs%ks%diis)
!2.1.6) p2pComms: comgp
 wvl%wfs%ks%comgp = p2pComms_null()
!2.1.7) collective_comms: collcom
 wvl%wfs%ks%collcom = comms_linear_null()
!2.1.8) collective_comms: collcom_sr
 wvl%wfs%ks%collcom_sr = comms_linear_null()

!2.2) GPU_pointers: GPU
 call nullify_GPU_pointers(wvl%wfs%GPU)
 wvl%wfs%GPU%OCLconv = .false.

!3) wvl_internal_type: descr
!3.1) locreg_descriptors: Glr
 call nullify_locreg_descriptors(wvl%descr%Glr)
!3.2) atoms_data: atoms
 call nullify_atoms_data(wvl%descr%atoms)
!3.3) rholoc_objects: rholoc
 call nullify_rholoc_objects(wvl%descr%rholoc)
!3.4) paw_objects: paw 
 call nullify_paw_objects(wvl%descr%paw)

!4) wvl_denspot_type: den
!4.1) DFT_local_fields: denspot
 call nullify_DFT_local_fields(wvl%den%denspot)

!5) wvl_energy_terms: e
!there are no pointers here

 DBG_EXIT("COLL")

#else
 ABI_UNUSED(wvl%den%symObj)
#endif

#if defined HAVE_BIGDFT
 contains
!!***

!!****f* nullify_wvl_data/comms_cubic_null
!! NAME
!!  comms_cubic_null
!!
!! FUNCTION
!!  Nullify comms
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function comms_cubic_null() result(comms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'comms_cubic_null'
!End of the abilint section

  implicit none
  type(comms_cubic) :: comms
! *************************************************************************
   nullify(comms%ncntd)
   nullify(comms%ncntt)
   nullify(comms%ndspld)
   nullify(comms%ndsplt)
   nullify(comms%nvctr_par)
 end function comms_cubic_null
!!***
#endif

end subroutine nullify_wvl_data
!!***

