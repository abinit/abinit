!!****m* ABINIT/m_abi_etsf
!! NAME
!! m_abi_etsf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (DCA,YP,MJV,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  Remove this module. At present it's used only in m_iowf.F90
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abi_etsf

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_dtset
 use m_xmpi
 use m_errors
 use m_atomdata
 use m_nctk
 use iso_c_binding

 use m_fstrings,     only : endswith
 use defs_datatypes, only : pseudopotential_type

 implicit none

 private

 public :: abi_etsf_init

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_init
!! NAME
!! abi_etsf_init
!!
!! FUNCTION
!!  Create a NetCDF file following the ETSF file format specifications.
!!  It declares the dimensions and set-up the different variables, according
!!  to the itype argument.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filapp = character string giving the root to form the name of the GEO file
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!  kdep= .true. if the data for the array sizes are dependant on k points.
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs <type(wvl_projector_type)>=wavefunctions information for wavelets.
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!
!! SOURCE

subroutine abi_etsf_init(dtset,filapp,itype,kdep,lmn_size,psps,wfs)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 logical,intent(in) :: kdep
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),optional,intent(in) :: wfs
!arrays
 integer,intent(in) :: lmn_size(psps%npsp)

!Local variables-------------------------------

! *************************************************************************

 ABI_UNUSED(dtset%ntypat)
 ABI_UNUSED(filapp)
 ABI_UNUSED(itype)
 ABI_UNUSED(kdep)
 ABI_UNUSED(lmn_size)
 if(.false. .and. present(wfs))write(std_out,*)' OK'

end subroutine abi_etsf_init
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/ini_wf_etsf
!! NAME
!! ini_wf_etsf
!!
!! FUNCTION
!! Do initialization of additional dimensions and variables in wavefunction files in ETSF format.
!!
!! INPUTS
!!  usewvl=1 if wavelets are used, 0 otherwise
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  npsp=number of pseudopotentials.
!!  ntypat=number of types of atoms in cell.
!!  ncid=the unit of the open NetCDF file.
!!
!! SIDE EFFECTS
!!  New dimensions and variables are added to the initial NetCDF file.
!!
!! PARENTS
!!      m_abi_etsf
!!
!! CHILDREN
!!
!! SOURCE

subroutine ini_wf_etsf(ncid,usewvl,lmn_size,npsp,ntypat)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,usewvl,npsp,ntypat
!arrays
 integer,intent(in) :: lmn_size(npsp)

!Local variables-------------------------------

! *************************************************************************

!if ETSF_IO is undefined, do nothing
 ABI_UNUSED(ncid)
 ABI_UNUSED(ntypat)
 ABI_UNUSED(lmn_size)
 ABI_UNUSED(usewvl)

end subroutine ini_wf_etsf
!!***

!----------------------------------------------------------------------

END MODULE m_abi_etsf
!!***
