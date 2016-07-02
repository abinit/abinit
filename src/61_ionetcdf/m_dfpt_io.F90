!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfpt_io
!! NAME
!! m_dfpt_io
!!
!! FUNCTION
!! Module containing the methods used to do IO on DFPT results. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (DW)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_dfpt_io

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_nctk
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

 implicit none

 private 
!!***

 public :: elast_ncwrite    ! Dump the elastic and compliance tensors in a netcdf file.

CONTAINS

!!****f* m_dfpt_io/elast_ncwrite 
!! NAME
!!  elast_ncwrite
!!
!! FUNCTION
!! Writes the elastic constants tensors (clamped-ion, relaxed-ion with and
!! without stress corrections) to a netCDF file.
!!
!! INPUTS
!! compl=relaxed-ion compliance tensor(without stress correction) (6*6) in
!!       Voigt notation
!! compl_clamped=clamped-ion compliance tensor(without stress correction) (6*6)
!!               in Voigt notation
!! compl_stress=relaxed-ion compliance tensor(with stress correction) (6*6) in
!!              Voigt notation
!! elast=relaxed-ion elastic tensor(without stress correction) (6*6) in
!!       Voigt notation
!! elast_clamped=clamped-ion elastic tensor(without stress correction) (6*6)
!!               in Voigt notation
!! elast_stress=relaxed-ion elastic tensor(with stress correction) (6*6) in
!!              Voigt notation
!! ncid=NC file handle (open in the caller)
!!
!! OUTPUT
!! Only writing
!!
!! NOTES
!! Units are GPa for elastic constants and GPa^-1 for compliance constants
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE


subroutine elast_ncwrite(compl,compl_clamped,compl_stress,elast,elast_clamped,elast_stress,ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elast_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid 
!arrays
 real(dp),intent(in) :: compl(6,6), compl_clamped(6,6),compl_stress(6,6)
 real(dp),intent(in) :: elast(6,6), elast_clamped(6,6),elast_stress(6,6)

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_NETCDF
 integer :: ncerr

! *************************************************************************

! Define dimensions
 NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))

!arrays
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('compliance_constants_relaxed_ion', "dp", 'six, six'), &
   nctkarr_t('compliance_constants_clamped_ion', "dp", 'six, six'), &
   nctkarr_t('compliance_constants_relaxed_ion_stress_corrected', "dp", "six, six"), &
   nctkarr_t('elastic_constants_relaxed_ion', "dp", 'six, six'), &
   nctkarr_t('elastic_constants_clamped_ion', "dp", 'six, six'), &
   nctkarr_t('elastic_constants_relaxed_ion_stress_corrected', "dp", 'six, six')])
 NCF_CHECK(ncerr)

 ! Write variables. 
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('compliance_constants_relaxed_ion'), compl))
 NCF_CHECK(nf90_put_var(ncid, vid('compliance_constants_clamped_ion'), compl_clamped))
 NCF_CHECK(nf90_put_var(ncid, vid('compliance_constants_relaxed_ion_stress_corrected'), compl_stress))
 NCF_CHECK(nf90_put_var(ncid, vid('elastic_constants_relaxed_ion'), elast))
 NCF_CHECK(nf90_put_var(ncid, vid('elastic_constants_clamped_ion'), elast_clamped))
 NCF_CHECK(nf90_put_var(ncid, vid('elastic_constants_relaxed_ion_stress_corrected'), elast_stress))

#else
 MSG_ERROR("Netcdf support not enabled")
 ABI_UNUSED((/ncid/))
#endif

contains
 integer function vid(vname) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine elast_ncwrite
!!***

end MODULE m_dfpt_io
!!***
