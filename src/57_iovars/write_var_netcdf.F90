!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_var_netcdf
!!
!! NAME
!! write_var_netcdf
!!
!! FUNCTION
!! Write the history into a netcdf dataset
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! arr_int
!! arr_real
!! marr
!! narr
!! typevar
!! varname
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      prttagm,prttagm_images
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine write_var_netcdf(arr_int,arr_real,marr,narr,ncid,typevar,varname)

 use defs_basis
 use m_profiling_abi
 use m_errors 
#if defined HAVE_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_var_netcdf'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: narr,marr,ncid
 character(len=*),intent(in) :: varname
 character(len=3),intent(in) :: typevar
!arrays
 integer,intent(in) :: arr_int(marr)
 real(dp),intent(in) :: arr_real(marr)

!Local variables-------------------------------
!scalars
 integer :: var_id,var_type,vardim_id,ncerr
 !character(len=500) :: msg

! *************************************************************************

 !write(std_out,*)"about to write varname: ",trim(varname)

#if defined HAVE_NETCDF
 if (ncid>0) then
!  ### Put the file in definition mode
   ncerr=nf90_redef(ncid)
   if (ncerr/=NF90_NOERR.and.ncerr/=NF90_EINDEFINE) then
     NCF_CHECK_MSG(ncerr,'nf90_redef')
   end if
!  ### Define the dimensions
   if (narr==1)then
     ncerr=nf90_inq_dimid(ncid,'one',vardim_id)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
   else
     ncerr=nf90_def_dim(ncid,trim(varname),narr,vardim_id)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')
   end if
!  ### Define the variables
   if (typevar=='INT') then
     var_type=NF90_INT
   else if (typevar=='DPR') then
     var_type=NF90_DOUBLE
   end if
   ncerr=nf90_def_var(ncid, trim(varname), var_type, vardim_id, var_id)
   NCF_CHECK_MSG(ncerr,'nf90_def_var')
!  ### Put the file in data mode
   ncerr=nf90_enddef(ncid)
   if (ncerr/=NF90_NOERR.and.ncerr/=NF90_ENOTINDEFINE) then
     NCF_CHECK_MSG(ncerr,'nf90_enddef')
   end if
!  ### Write variables into the dataset
   if (typevar=='INT') then
     ncerr=nf90_put_var(ncid,var_id,arr_int,start=(/1/),count=(/narr/))
   else if (typevar=='DPR') then
     ncerr=nf90_put_var(ncid,var_id,arr_real,start=(/1/),count=(/narr/))
   end if
   NCF_CHECK_MSG(ncerr,'nf90_put_var')
 end if
#endif

end subroutine write_var_netcdf
!!***
