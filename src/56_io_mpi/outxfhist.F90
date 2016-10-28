!{\src2tex{textfont=tt}}
!!****f* ABINIT/outxfhist
!! NAME
!! outxfhist
!!
!! FUNCTION
!!  read/write xfhist
!!
!! COPYRIGHT
!! Copyright (C) 2003-2016 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  option =
!!   1: write
!!   2: read only nxfh
!!   3: read xfhist
!!  response =
!!   0: GS wavefunctions
!!   1: RF wavefunctions
!!  natom = number of atoms in unit cell
!!  mxfh = last dimension of the xfhist array
!!
!! OUTPUT
!!  ios = error code returned by read operations
!!
!! SIDE EFFECTS
!!  nxfh = actual number of (x,f) history pairs, see xfhist array
!!  wff2 = structured info for wavefunctions
!!  xfhist(3,natom+4,2,ab_xfh%mxfh) = (x,f) history array, also including
!!   rprim and stress
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      xderiveread,xderiverrecend,xderiverrecinit,xderivewrecend
!!      xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outxfhist(ab_xfh,natom,option,wff2,ios)

 use defs_basis
 use m_profiling_abi
 use m_abimover
 use m_xmpi
 use m_wffile
 use m_errors
#if defined HAVE_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outxfhist'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer          ,intent(in)    :: natom,option
 integer          ,intent(out)   :: ios
 type(wffile_type),intent(inout)    :: wff2
 type(ab_xfh_type),intent(inout) :: ab_xfh

!Local variables-------------------------------
 integer :: ierr,ixfh,ncid_hdr,spaceComm,xfdim2
 real(dp),allocatable :: xfhist_tmp(:)
 character(len=500) :: message
!no_abirules
#if defined HAVE_NETCDF
 integer :: ncerr
 integer :: nxfh_id, mxfh_id, xfdim2_id, dim2inout_id, dimr3_id,xfhist_id
 integer :: nxfh_tmp,mxfh_tmp,xfdim2_tmp,dim2inout_tmp
#endif


! *************************************************************************

!DEBUG
!write(std_out,*)'outxfhist  : enter, option = ', option
!ENDDEBUG
 ncid_hdr = wff2%unwff
 xfdim2 = natom+4

 ios = 0

!### (Option=1) Write out content of all iterations
!#####################################################################
 if ( option == 1 ) then

!  Write the (x,f) history
   if (wff2%iomode == IO_MODE_FORTRAN) then
     write(unit=wff2%unwff)ab_xfh%nxfh
     do ixfh=1,ab_xfh%nxfh
       write(unit=wff2%unwff)ab_xfh%xfhist(:,:,:,ixfh)
     end do

   else if (wff2%iomode == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'iomode == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     MSG_ERROR(message)

     write(unit=wff2%unwff)ab_xfh%nxfh
     do ixfh=1,ab_xfh%nxfh
       write(unit=wff2%unwff)ab_xfh%xfhist(:,:,:,ixfh)
     end do

!    insert mpi broadcast here

   else if(wff2%iomode==IO_MODE_MPI)then
     ABI_ALLOCATE(xfhist_tmp,(3*(natom+4)*2))
     spaceComm=xmpi_comm_self
     call xderiveWRecInit(wff2,ierr)
     call xderiveWrite(wff2,ab_xfh%nxfh,ierr)
     call xderiveWRecEnd(wff2,ierr)
     do ixfh=1,ab_xfh%nxfh
       xfhist_tmp(:)=reshape(ab_xfh%xfhist(:,:,:,ixfh),(/3*(natom+4)*2/))
       call xderiveWRecInit(wff2,ierr)
       call xderiveWrite(wff2,xfhist_tmp,3*(natom+4)*2,spaceComm,ierr)
       call xderiveWRecEnd(wff2,ierr)
     end do
     ABI_DEALLOCATE(xfhist_tmp)

#if defined HAVE_NETCDF
   else if (wff2%iomode == IO_MODE_NETCDF) then
!    check if nxfh and xfhist are defined
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)

     if (ncerr /= NF90_NOERR) then
!      need to define everything
       ncerr = nf90_redef (ncid=ncid_hdr)
       NCF_CHECK_MSG(ncerr," outxfhist : going to define mode ")

       ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim2inout",len=2,dimid=dim2inout_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define dim2inout")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="mxfh",len=ab_xfh%mxfh,dimid=mxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define mxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="nxfh",len=ab_xfh%nxfh,dimid=nxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define nxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="xfdim2",len=xfdim2,dimid=xfdim2_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define xfdim2")

       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimr3",dimid=dimr3_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire dimr3")

!      ab_xfh%xfhist(3,natom+4,2,ab_xfh%mxfh)
       ncerr = nf90_def_var(ncid=ncid_hdr,name="xfhist",xtype=NF90_DOUBLE,&
&       dimids=(/dimr3_id,xfdim2_id,dim2inout_id,mxfh_id/),varid=xfhist_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define xfhist")

!      End define mode and go to data mode
       ncerr = nf90_enddef(ncid=ncid_hdr)
       NCF_CHECK_MSG(ncerr," outxfhist : enddef call ")
     else
!      check that the dimensions are correct
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&       len=nxfh_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get nxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="xfdim2",dimid=xfdim2_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire xfdim2")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=xfdim2_id,&
&       len=xfdim2_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get xfdim2")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mxfh",dimid=mxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire mxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mxfh_id,&
&       len=mxfh_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get mxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dim2inout",dimid=dim2inout_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire dim2inout")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=dim2inout_id,&
&       len=dim2inout_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get dim2inout")

       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xfhist",varid=xfhist_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire xfhist")

       if (mxfh_tmp /= ab_xfh%mxfh .or. dim2inout_tmp /= 2 .or. xfdim2_tmp /= xfdim2) then
         write (message,"(A)") 'outxfhist : ERROR xfhist has bad dimensions in NetCDF file. Can not re-write it.'
         MSG_ERROR(message)
       end if

     end if

!    Now fill the data
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=xfhist_id,values=ab_xfh%xfhist)
     NCF_CHECK_MSG(ncerr," outxfhist : fill xfhist")

!    end NETCDF definition ifdef
#endif
   end if  ! end iomode if

!  ### (Option=2) Read in number of iterations
!  #####################################################################
 else if ( option == 2 ) then

   if (wff2%iomode == IO_MODE_FORTRAN) then
     read(unit=wff2%unwff,iostat=ios)ab_xfh%nxfh

   else if (wff2%iomode == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'iomode == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     MSG_ERROR(message)

     read(unit=wff2%unwff,iostat=ios)ab_xfh%nxfh

   else if (wff2%iomode == IO_MODE_MPI) then
     call xderiveRRecInit(wff2,ierr)
     call xderiveRead(wff2,ab_xfh%nxfh,ierr)
     call xderiveRRecEnd(wff2,ierr)

#if defined HAVE_NETCDF
   else if (wff2%iomode == IO_MODE_NETCDF) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
     NCF_CHECK_MSG(ncerr," outxfhist : inquire nxfh")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&     len=ab_xfh%nxfh)
     NCF_CHECK_MSG(ncerr,"  outxfhist : get nxfh")
#endif
   end if

!  ### (Option=3) Read in iteration content
!  #####################################################################
 else if ( option == 3 ) then
   if (wff2%iomode == IO_MODE_FORTRAN) then
     do ixfh=1,ab_xfh%nxfhr
       read(unit=wff2%unwff,iostat=ios)ab_xfh%xfhist(:,:,:,ixfh)
     end do
   else if (wff2%iomode == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'iomode == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     MSG_ERROR(message)

     do ixfh=1,ab_xfh%nxfhr
       read(unit=wff2%unwff,iostat=ios)ab_xfh%xfhist(:,:,:,ixfh)
     end do

   else if (wff2%iomode == IO_MODE_MPI) then
     ABI_ALLOCATE(xfhist_tmp,(3*(natom+4)*2))
     spaceComm=xmpi_comm_self
     do ixfh=1,ab_xfh%nxfhr
       call xderiveRRecInit(wff2,ierr)
       call xderiveRead(wff2,xfhist_tmp,3*(natom+4)*2,spaceComm,ierr)
       call xderiveRRecEnd(wff2,ierr)
       xfhist_tmp(:)=xfhist_tmp(:)
     end do
     ABI_DEALLOCATE(xfhist_tmp)
   end if

!  FIXME: should this be inside the if not mpi as above for options 1 and 2?
!  it is placed here because the netcdf read is a single operation
#if defined HAVE_NETCDF
   if (wff2%iomode == IO_MODE_NETCDF) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
     NCF_CHECK_MSG(ncerr," outxfhist : inquire nxfh")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&     len=ab_xfh%nxfhr)
     NCF_CHECK_MSG(ncerr,"  outxfhist : get nxfh")

     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=xfhist_id,name="xfhist")
     NCF_CHECK_MSG(ncerr," outxfhist : inquire xfhist")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=xfhist_id,values=ab_xfh%xfhist,&
&     start=(/1,1,1,1/),count=(/3,natom+4,2,ab_xfh%nxfhr/))
     NCF_CHECK_MSG(ncerr," outxfhist : read xfhist")
   end if
#endif

 else
!  write(std_out,*)' outxfhist : option ', option , ' not available '
   write(message, "(A,A,A,A,I3,A)") ch10, "outxfhist: ERROR -", ch10, &
&   "option ", option, " not available."
   MSG_ERROR(message)
 end if

end subroutine outxfhist
!!***
