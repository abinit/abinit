!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_eig
!!
!! NAME
!! write_eig
!!
!! FUNCTION
!! Write the eigenvalues band by band and k point by k point
!! in a NetCDF file format
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
!! filname = Filename of the file where the history will be stored
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      clnup1
!!
!! CHILDREN
!!      ab_define_var
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine write_eig(eigen,filename,kptns,mband,nband,nkpt,nsppol)

 use defs_basis
 use m_profiling_abi
 use m_errors
#ifdef HAVE_NETCDF
 use netcdf
 use m_nctk,         only : ab_define_var
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_eig'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer,intent(in) :: nkpt,nsppol,mband
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ncerr,ncid,ii
 integer :: xyz_id,nkpt_id,mband_id,nsppol_id
 integer :: eig_id,kpt_id,nbk_id,nbk
 integer :: ikpt,isppol,nband_k,band_index
 real(dp):: convrt
!arrays
 integer :: dimEIG(3),dimKPT(2),dimNBK(2)
 integer :: count2(2),start2(2)
 integer :: count3(3),start3(3)
 real(dp):: band(mband)

! *********************************************************************

#if defined HAVE_NETCDF

 convrt=1.0_dp

!1. Create netCDF file
 ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER, ncid=ncid)
 NCF_CHECK_MSG(ncerr," create netcdf EIG file")

!2. Define dimensions
 ncerr = nf90_def_dim(ncid,"xyz",3,xyz_id)
 NCF_CHECK_MSG(ncerr," define dimension xyz")

 ncerr = nf90_def_dim(ncid,"mband",mband,mband_id)
 NCF_CHECK_MSG(ncerr," define dimension mband")

 ncerr = nf90_def_dim(ncid,"nkpt",nkpt,nkpt_id)
 NCF_CHECK_MSG(ncerr," define dimension nkpt")

 ncerr = nf90_def_dim(ncid,"nsppol",nsppol,nsppol_id)
 NCF_CHECK_MSG(ncerr," define dimension nsppol")

!Dimensions for EIGENVALUES
 dimEIG = (/ mband_id, nkpt_id, nsppol_id /)
!Dimensions for kpoint positions
 dimKPT = (/ xyz_id, nkpt_id /)
!Dimensions for number kpoints per band and spin
!dimNBK = (/ nkpt_id, nsppol_id /)
 dimNBK = (/ nkpt_id, nsppol_id /)

!3. Define variables

 call ab_define_var(ncid, dimEIG, eig_id, NF90_DOUBLE,&
& "Eigenvalues",&
& "Values of eigenvalues",&
& "Hartree")
 call ab_define_var(ncid, dimKPT, kpt_id, NF90_DOUBLE,"Kptns",&
& "Positions of K-points in reciprocal space",&
& "Dimensionless")
 call ab_define_var(ncid, dimNBK, nbk_id, NF90_INT,"NBandK",&
& "Number of bands per kpoint and Spin",&
& "Dimensionless")

!4. End define mode
 ncerr = nf90_enddef(ncid)
 NCF_CHECK_MSG(ncerr," end define mode")

!5 Write kpoint positions
 do ikpt=1,nkpt
   start2 = (/ 1, ikpt /)
   count2 = (/ 3, 1 /)
   ncerr = nf90_put_var(ncid, kpt_id,&
&   kptns(1:3,ikpt),&
&   start = start2,&
&   count = count2)
   NCF_CHECK_MSG(ncerr," write variable kptns")
 end do


!6 Write eigenvalues
 band_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     start3 = (/ 1, ikpt, isppol /)
     count3 = (/ mband, 1, 1 /)
     band(:)=zero
     do ii=1,nband_k
       band(ii)=eigen(band_index+ii)
     end do
     ncerr = nf90_put_var(ncid, eig_id,&
&     band,&
&     start = start3,&
&     count = count3)
     NCF_CHECK_MSG(ncerr," write variable band")

     band_index=band_index+nband_k
   end do
 end do

!6 Write Number of bands per kpoint and Spin

 do isppol=1,nsppol
   do ikpt=1,nkpt
     start2 = (/ ikpt, 1 /)
     count2 = (/ 1, 1 /)
     nbk=nband(ikpt+(isppol-1)*nkpt)
     ncerr = nf90_put_var(ncid, nbk_id,&
&     nbk,&
&     start = start2)
     NCF_CHECK_MSG(ncerr," write variable nband")
   end do
 end do

!7 Close file

 ncerr = nf90_close(ncid)
 NCF_CHECK_MSG(ncerr," close netcdf EIG file")
#endif

end subroutine write_eig
!!***
