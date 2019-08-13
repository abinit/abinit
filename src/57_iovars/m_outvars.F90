!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_outvars
!! NAME
!!  m_outvars
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR)
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

module m_outvars

 use defs_basis
 use defs_abitypes
 use m_results_out
 use m_dtset
 use m_abicore
 use m_errors
 use m_xomp
 use m_xmpi
#if defined HAVE_NETCDF
 use netcdf
#endif
 use m_outvar_a_h
 use m_outvar_i_n
 use m_outvar_o_z

 use m_nctk,      only : create_nc_file

 implicit none

 private
!!***

 public :: outvars
!!***

contains
!!***

!!****f* ABINIT/outvars
!! NAME
!! outvars
!!
!! FUNCTION
!! Echo variables for the ABINIT code.
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  mxvals=maximum size of some arrays along all datasets, including:
!!         ga_n_rules =maximal value of input ga_n_rules for all the datasets
!!         gw_nqlwl   =maximal value of input gw_nqlwl for all the datasets
!!         lpawu      =maximal value of input lpawu for all the datasets
!!         mband      =maximum number of bands
!!         natom      =maximal value of input natom for all the datasets
!!         natpawu    =maximal value of number of atoms on which +U is applied for all the datasets
!!         natsph     =maximal value of input natsph for all the datasets
!!         natvshift  =maximal value of input natvshift for all the datasets
!!         nconeq     =maximal value of input nconeq for all the datasets
!!         nimage     =maximal value of input nimage for all the datasets
!!         nimfrqs    =maximal value of input cd_customnimfrqs for all the datasets
!!         nkpt       =maximal value of input nkpt for all the datasets
!!         nkptgw     =maximal value of input nkptgw for all the datasets
!!         nkpthf     =maximal value of input nkpthf for all the datasets
!!         nnos       =maximal value of input nnos for all the datasets
!!         nqptdm     =maximal value of input nqptdm for all the datasets
!!         nspinor    =maximal value of input nspinor for all the datasets
!!         nsppol     =maximal value of input nsppol for all the datasets
!!         nsym       =maximum number of symmetries
!!         ntypat     =maximum number of type of atoms
!!         nzchempot  =maximal value of input nzchempot for all the datasets
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!   one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  timopt=input variable to modulate the timing
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      create_nc_file,outvar_a_h,outvar_i_n,outvar_o_z,wrtout
!!
!! SOURCE

subroutine outvars(choice,dmatpuflag,dtsets,filnam4,iout,&
&  mxvals,ndtset,ndtset_alloc,npsp,results_out,timopt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,dmatpuflag,iout
 integer,intent(in) :: ndtset,ndtset_alloc,npsp,timopt
 type(ab_dimensions),intent(in) :: mxvals
 character(len=*),intent(in) :: filnam4
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: first,idtset,iimage,kptopt
 integer :: marr,mu,ncerr
 integer :: nshiftk
 integer :: prtvol_glob,max_nthreads
 integer :: rfddk,rfelfd,rfphon,rfstrs,rfuser,rfmagn,rf2_dkdk,rf2_dkde
 integer :: ncid=0 ! Variables for NetCDF output
 character(len=500) :: message
 character(len=4) :: stringimage
 type(ab_dimensions) :: multivals
!arrays
 integer,allocatable :: jdtset_(:),response_(:)
 character(len=8),allocatable :: strimg(:)

! *************************************************************************

!Set up a 'global' prtvol value
 prtvol_glob=1
 if(sum((dtsets(:)%prtvol)**2)==0)prtvol_glob=0

!###########################################################
!### 00. Echo of selected default values

 if(choice==1)then

   max_nthreads = xomp_get_max_threads()
#ifndef HAVE_OPENMP
   max_nthreads = 0 ! this value signals that OMP is not enabled in ABINIT.
#endif

   write(iout, '(10a)' )&
&   '--------------------------------------------------------------------------------',ch10,&
&   '------------- Echo of variables that govern the present computation ------------',ch10,&
&   '--------------------------------------------------------------------------------',ch10,&
&   '-',ch10,&
&   '- outvars: echo of selected default values                                      '
   write(iout, '(3(a,i3),2a)' )&
&   '-   iomode0 =',dtsets(0)%iomode,' , fftalg0 =',dtsets(0)%ngfft(7),' , wfoptalg0 =',dtsets(0)%wfoptalg,ch10,&
&   '-'
   write(iout, '(3a,(a,i5),2a)' )&
&   '- outvars: echo of global parameters not present in the input file              ',ch10,&
&   '- ',' max_nthreads =',max_nthreads,ch10,&
&   '-'
 end if

!write(std_out,*) 'outvar 01'
!###########################################################
!### 01. First line indicating outvars

 if(choice==1)then
   write(iout, '(a)' )&
&   ' -outvars: echo values of preprocessed input variables --------'
 else
   write(iout, '(a)' )&
&   ' -outvars: echo values of variables after computation  --------'
 end if

!###########################################################
!### 02. Open NetCDF file for export variables

#ifdef HAVE_NETCDF
 ! Enable netcdf output only if the number of datasets is small.
 ! otherwise v6[34] crashes with errmess:
 !    nf90_def_dim - NetCDF library returned:   NetCDF: NC_MAX_DIMS exceeded
 ! because we keep on creating dimensions in write_var_netcdf.
 ! one should use groups for this kind of operations!!

 ncid = 0
 if (ndtset_alloc  < 10) then
   if (iout==std_out)then
     write(iout,*) ch10,' These variables are accessible in NetCDF format (',trim(filnam4)//'_OUT.nc',')',ch10
   end if
   call create_nc_file(trim(filnam4)//"_OUT.nc",ncid)

   if (dtsets(1)%prtvol==-2) then
     if (ncid>0)then
       ncid=-ncid
     else
       ncid=-1
     end if
   end if
 else
   MSG_WARNING("output of OUT.nc has been disabled. Too many datasets")
 end if
#endif
 !ncid = 0

!###########################################################
!##1 03. Set up dimensions : determine whether these are different for different datasets.

 multivals%ga_n_rules=0
 multivals%gw_nqlwl=0
 multivals%mband=0
 multivals%natom=0
 multivals%natpawu=0
 multivals%natsph=0
 multivals%natvshift=0
 multivals%nberry=0
 multivals%nbandhf=0
 multivals%nconeq=0
 multivals%nfreqsp=0
 multivals%nimage=0
 multivals%nimfrqs=0
 multivals%nkpt=0
 multivals%nkptgw=0
 multivals%nkpthf=0
 multivals%nnos=0
 multivals%nqptdm=0
 multivals%nshiftk=0
 multivals%nsp=0
 multivals%nspinor=0
 multivals%nsppol=0
 multivals%nsym=0
 multivals%ntypat=0
 multivals%ntypalch=0
 multivals%nzchempot=0

 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%ga_n_rules/=dtsets(idtset)%ga_n_rules) multivals%ga_n_rules =1
     if(dtsets(1)%gw_nqlwl /=dtsets(idtset)%gw_nqlwl ) multivals%gw_nqlwl =1
     if(dtsets(1)%mband    /=dtsets(idtset)%mband    ) multivals%mband    =1
     if(dtsets(1)%natom    /=dtsets(idtset)%natom    ) multivals%natom    =1
     if(dtsets(1)%natpawu  /=dtsets(idtset)%natpawu  ) multivals%natpawu  =1
     if(dtsets(1)%natsph   /=dtsets(idtset)%natsph   ) multivals%natsph   =1
     if(dtsets(1)%natvshift/=dtsets(idtset)%natvshift) multivals%natvshift=1
     if(dtsets(1)%nberry   /=dtsets(idtset)%nberry   ) multivals%nberry   =1
     if(dtsets(1)%nbandhf  /=dtsets(idtset)%nbandhf  ) multivals%nbandhf  =1
     if(dtsets(1)%nconeq   /=dtsets(idtset)%nconeq   ) multivals%nconeq   =1
     if(dtsets(1)%nfreqsp  /=dtsets(idtset)%nfreqsp  ) multivals%nfreqsp  =1
     if(dtsets(1)%nimage   /=dtsets(idtset)%nimage   ) multivals%nimage   =1
     if(dtsets(1)%cd_customnimfrqs  /=dtsets(idtset)%cd_customnimfrqs  ) multivals%nimfrqs  =1
     if(dtsets(1)%nkpt     /=dtsets(idtset)%nkpt     ) multivals%nkpt     =1
     if(dtsets(1)%nkptgw   /=dtsets(idtset)%nkptgw   ) multivals%nkptgw   =1
     if(dtsets(1)%nkpthf*dtsets(1)%usefock /=dtsets(idtset)%nkpthf*dtsets(idtset)%usefock) multivals%nkpthf=1
     if(dtsets(1)%nnos     /=dtsets(idtset)%nnos     ) multivals%nnos     =1
     if(dtsets(1)%nqptdm   /=dtsets(idtset)%nqptdm   ) multivals%nqptdm   =1
     if(dtsets(1)%nsppol*dtsets(1)%nspinor/=dtsets(idtset)%nsppol*dtsets(idtset)%nspinor) multivals%nsp=1
     if(dtsets(1)%nsppol   /=dtsets(idtset)%nsppol   ) multivals%nsppol   =1
     if(dtsets(1)%nspinor  /=dtsets(idtset)%nspinor  ) multivals%nspinor  =1
     if(dtsets(1)%nsym     /=dtsets(idtset)%nsym     ) multivals%nsym     =1
     if(dtsets(1)%ntypat   /=dtsets(idtset)%ntypat   ) multivals%ntypat   =1
     if(dtsets(1)%ntypalch /=dtsets(idtset)%ntypalch ) multivals%ntypalch =1
     if(dtsets(1)%nzchempot/=dtsets(idtset)%nzchempot) multivals%nzchempot=1
   end do
 end if

!DEBUG
! write(std_out,*)' outvars : multivals%nkpthf =',multivals%nkpthf
! write(std_out,*)' outvars : dtsets(1:ndtset_alloc)%nkpthf =',dtsets(1:ndtset_alloc)%nkpthf
!ENDDEBUG

 nshiftk=1
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   first=0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>=1)then
       if(first==0)then
         first=1
         nshiftk=dtsets(idtset)%nshiftk
       else
         if(nshiftk/=dtsets(idtset)%nshiftk)multivals%nshiftk=1
       end if
     end if
   end do
 end if

!###########################################################
!### 04. Determine whether each dataset is (or not) a response calculation
!## (should use optdriver, isn't it ?)

 ABI_ALLOCATE(response_,(ndtset_alloc))
 response_(:)=0
 do idtset=1,ndtset_alloc
   rfddk=dtsets(idtset)%rfddk
   rfelfd=dtsets(idtset)%rfelfd
   rfphon=dtsets(idtset)%rfphon
   rfstrs=dtsets(idtset)%rfstrs
   rfuser=dtsets(idtset)%rfuser
   rfmagn=dtsets(idtset)%rfmagn
   rf2_dkdk=dtsets(idtset)%rf2_dkdk
   rf2_dkde=dtsets(idtset)%rf2_dkde
   if(rfddk/=0 .or. rfelfd/=0 .or. rfphon/=0 .or. rfstrs/=0 .or. &
&   rfuser/=0 .or. rf2_dkdk/=0 .or. rf2_dkde/=0 .or. rfmagn/=0)then
     response_(idtset)=1
   end if
 end do

!###########################################################
!### 05. Determine size of work arrays

 marr=max(3*mxvals%natom,&
& mxvals%natsph,&
& mxvals%natvshift*mxvals%nsppol*mxvals%natom,&
& 3*mxvals%nberry,&
& mxvals%nimage,&
& 3*mxvals%nkptgw,&
& 3*mxvals%nkpthf,&
& mxvals%nkpt*mxvals%nsppol*mxvals%mband,&
& 3*mxvals%nkpt,npsp,&
& 3*mxvals%nqptdm,&
& mxvals%ntypat,&
& 9*mxvals%nsym,3*8,&
& 3*mxvals%natom*mxvals%nconeq,&
& mxvals%nnos,&
& 3*mxvals%nqptdm,&
& 3*mxvals%nzchempot*mxvals%ntypat,&
& 3*mxvals%gw_nqlwl,&
& (2*mxvals%lpawu+1)**2*max(mxvals%nsppol,mxvals%nspinor)*mxvals%natpawu*dmatpuflag,&
& 30 ) ! used by ga_rules TODO : replace with mxvals% ga_n_rules

!###########################################################
!### 06. Initialize strimg

 ABI_ALLOCATE(strimg,(mxvals%nimage))
 do iimage=1,mxvals%nimage
   if(iimage<10)then
     write(stringimage,'(i1)')iimage
   else if(iimage<100)then
     write(stringimage,'(i2)')iimage
   else if(iimage<1000)then
     write(stringimage,'(i3)')iimage
   else if(iimage<10000)then
     write(stringimage,'(i4)')iimage
   end if
   strimg(iimage)='_'//trim(stringimage)//'img'
 end do
 strimg(1)=''

!###########################################################
!### 07. Initialize jdtset_

 ABI_ALLOCATE(jdtset_,(0:ndtset_alloc))
 jdtset_(0:ndtset_alloc)=dtsets(0:ndtset_alloc)%jdtset


!###########################################################
!### 08. Print variables, for different ranges of names

 call outvar_a_h(choice,dmatpuflag,dtsets,iout,jdtset_,marr,multivals,mxvals,&
& ncid,ndtset,ndtset_alloc,results_out,strimg)

 call outvar_i_n(dtsets,iout,jdtset_,marr,multivals,mxvals,&
& ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,response_,results_out,strimg)

 call outvar_o_z(choice,dtsets,iout,&
& jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,&
& results_out,strimg,timopt)


!###########################################################
!## Deallocations and cleaning

 ABI_DEALLOCATE(jdtset_)
 ABI_DEALLOCATE(response_)
 ABI_DEALLOCATE(strimg)

 write(message,'(a,80a)')ch10,('=',mu=1,80)
 call wrtout(iout,message,'COLL')

#ifdef HAVE_NETCDF
 if (ncid /= 0) then
   ncerr=nf90_close(abs(ncid))
   if (ncerr/=nf90_NoErr) then
     message='Netcdf Error while closing the OUT.nc file: '//trim(nf90_strerror(ncerr))
     MSG_ERROR(message)
   end if
 end if
#endif
 if (.false.) write(std_out,*) ncerr

!**************************************************************************

end subroutine outvars
!!***

end module m_outvars
!!***
