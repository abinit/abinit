!{\src2tex{textfont=tt}}
!!****f* ABINIT/create_nc_file
!! NAME
!! create_nc_file
!!
!! FUNCTION
!! Create an NetCDF file including a dimension one definition
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! TODO:
!!  Remove
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine create_nc_file (filename,ncid)

 use defs_basis
 use m_profiling_abi
 use m_errors
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,      only : sjoin

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'create_nc_file'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ncid
!arrays
character(len=*),intent(in) :: filename

!Local variables-------------------------------
#if defined HAVE_NETCDF
integer :: one_id
integer :: ncerr
#endif

! *************************************************************************

 ncid = 0
#if defined HAVE_NETCDF
!Create the NetCDF file
 ncerr=nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=ncid)
 NCF_CHECK_MSG(ncerr, sjoin('Error while creating:', filename))
 ncerr=nf90_def_dim(ncid,'one',1,one_id)
 NCF_CHECK_MSG(ncerr,'nf90_def_dim')
#endif

 end subroutine create_nc_file
!!***
