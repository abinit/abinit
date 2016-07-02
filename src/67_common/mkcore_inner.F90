!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkcore_inner
!! NAME
!!  mkcore_inner
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2016 ABINIT group (T. Rangel)
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
!!      mkcore_paw,mkcore_wvl
!!
!! CHILDREN
!!      sort_dp,splint,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkcore_inner(corfra,core_mesh,dyfrx2,&
& grxc1,grxc2,grxc3,ifftsph,msz,natom,ncmax,nfft,nfgd,nfgd_r0,nspden,n3xccc,option,pawtab,&
& rmet,rprimd,rr,rred,rshpm1,strdia,ucvol,vxc,xccc3d)
    
 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_splines
 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore_inner'
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: msz,natom,ncmax,nfft,nfgd,nfgd_r0,nspden,n3xccc,option
 real(dp), intent(in)  :: ucvol
 real(dp), intent(out) :: grxc1,grxc2,grxc3
 !arrays
 integer , intent(in)  :: ifftsph(ncmax)
 real(dp), intent(in)  :: rmet(3,3)
 real(dp), intent(in)  :: rr(ncmax)
 real(dp),intent(inout) :: corfra(3,3)
 real(dp),intent(inout) :: dyfrx2(3,3,natom)
 real(dp),intent(in)   :: rshpm1,rprimd(3,3),vxc(nfft,nspden)
 real(dp),intent(inout)  :: strdia,xccc3d(n3xccc)
 type(pawrad_type),intent(in)::core_mesh
 type(pawtab_type),intent(in)::pawtab
 real(dp)::raux(nfgd),raux2(nfgd),raux3(nfgd) 
 real(dp),pointer::rred(:,:) !associated if usefull
!Local variables-------------------------------
 integer :: iatom,ifgd,ii,jj,mu,nu
 real(dp) :: factor
 real(dp) :: rr_tmp(nfgd)
 character(len=500) :: message                   ! to be uncommented, if needed
 !arrays
 integer :: iindex(nfgd)
 real(dp) ::rcart(3)
 
! *************************************************************************
 if(nspden >1) then
   write(message, '(a)')'mkcore_inner: this is not yet generalized to npsden>1'
   MSG_ERROR(message)
 end if
!if(.not. associated(raux2) .and. option >2) then
!write(message, '(a)')'mkcore_inner: raux2 is not associated and option >2'
!MSG_ERROR(message)
!end if
 if(.not. associated(rred) .and. option >2) then
   write(message, '(a)')'mkcore_inner: rred is not associated and option >2'
   MSG_ERROR(message)
 end if
!if(.not. associated(raux3) .and. option == 4) then
!write(message, '(a)')'mkcore_inner: raux3 is not associated and option == 4'
!MSG_ERROR(message)
!end if

 
!Start from options ne 1, since they all need the first derivative
!of the charge density. Moreover, option 2 is more frequent, according
!to routine "mkcore"

 rr_tmp=rr(1:nfgd)

 if(option .ne. 1) then
!  Evaluate spline fit of 1st der of core charge density
!  from tcoredens(:,2) and tcoredens(:,4)
   do ii=1,nfgd
     iindex(ii)=ii
   end do
   call sort_dp(nfgd,rr_tmp,iindex(1:nfgd),tol16)
   call splint(msz,core_mesh%rad,&
&   pawtab%tcoredens(:,2),pawtab%tcoredens(:,4),&
&   nfgd,rr_tmp,raux(1:nfgd))

 else

!  Option==1
   do ii=1,nfgd
     iindex(ii)=ii
   end do
   call sort_dp(nfgd,rr_tmp,iindex(1:nfgd),tol16)
!  Evaluate spline fit of core charge density
!  from tcoredens(:,1) and tcoredens(:,3)
   call splint(msz,core_mesh%rad,&
&   pawtab%tcoredens(:,1),pawtab%tcoredens(:,3),&
&   nfgd,rr_tmp,raux(1:nfgd))
!  Accumulate contributions to core density on the entire cell
   do ii=1,nfgd
     jj=ifftsph(iindex(ii))
     xccc3d(jj)=xccc3d(jj)+raux(ii)
   end do
!  DEBUG
!  do ii=1,msz
!  write(itmp,'(2(f15.7,1x))')core_mesh%rad(ii),pawtab%tcoredens(ii,1)
!  end do
!  do ii=1,nfgd
!  write(501,'(4f15.7)')rr(ii),rr(iindex(ii)),rr_tmp(ii),rr_tmp(iindex(ii))
!  write(500,'(3(f15.7,1x))')rr(ii),xccc3d(ifftsph(ii))
!  end do
!  END DEBUG
 end if
!
 if (option==2) then    
!  Accumulate contributions to Exc grandients
   raux2(1:nfgd)=vxc(ifftsph(1:nfgd),1)*raux(iindex(1:nfgd))
   raux2(1:nfgd)=raux2(1:nfgd)/rr_tmp(iindex(1:nfgd))
   do ifgd=1,nfgd
     grxc1=grxc1*rred(1,ifgd)*raux2(ifgd)
     grxc2=grxc2*rred(2,ifgd)*raux2(ifgd)
     grxc3=grxc3*rred(3,ifgd)*raux2(ifgd)
   end do
 elseif(option==3) then
!  Accumulate contributions to stress tensor
!  in reduced coordinates
   factor=rshpm1/real(nfft,dp)
   raux2(1:nfgd)=vxc(ifftsph(1:nfgd),1)*raux(iindex(1:nfgd))
   raux2(1:nfgd)=raux2(1:nfgd)*factor/rr_tmp(iindex(1:nfgd))
!  Write out the 6 symmetric components
   do ifgd=1,nfgd
     corfra(1,1)=corfra(1,1)+raux2(ifgd)*rred(1,ifgd)**2
     corfra(2,2)=corfra(2,2)+raux2(ifgd)*rred(2,ifgd)**2
     corfra(3,3)=corfra(3,3)+raux2(ifgd)*rred(3,ifgd)**2
     corfra(3,2)=corfra(3,2)+raux2(ifgd)*rred(3,ifgd)*rred(2,ifgd)
     corfra(3,1)=corfra(3,1)+raux2(ifgd)*rred(3,ifgd)*rred(1,ifgd)
     corfra(2,1)=corfra(2,1)+raux2(ifgd)*rred(2,ifgd)*rred(1,ifgd)
!    (the above still needs to be transformed to cartesian coords)
   end do

!  Also compute a diagonal term
!  Evaluate spline fit of core charge density
!  from tcoredens(:,1) and tcoredens(:3)
   call splint(msz,core_mesh%rad,&
&   pawtab%tcoredens(:,1),pawtab%tcoredens(:,3),&
&   nfgd,rr_tmp,raux(1:nfgd))
   raux2(1:nfgd)=vxc(ifftsph(1:nfgd),1)*raux(iindex(1:nfgd))
   do ifgd=1,nfgd
     strdia=strdia+raux2(ifgd)
   end do

 else if (option==4) then
!  Compute frozen-wf contribution to Dynamical matrix
   if(nfgd>0) then
!    Accumulate contributions to dynamical matrix
     factor=ucvol/real(nfft,dp)*rshpm1
     raux2(1:nfgd)=factor*vxc(ifftsph(1:nfgd),1)/rr_tmp(iindex(1:nfgd))
     raux2(1:nfgd)=raux2(1:nfgd)*raux(iindex(1:nfgd))

!    Evaluate spline fit of 2nd der of core charge density
!    from tcoredens(:,3) and tcoredens(:5)
     call splint(msz,core_mesh%rad,&
&     pawtab%tcoredens(:,3),pawtab%tcoredens(:,5),&
&     nfgd,rr,raux(1:nfgd))
     
     raux3(1:nfgd)=raux2(1:nfgd)*raux(iindex(1:nfgd))*rshpm1/rr_tmp(iindex(1:nfgd))
     
     do ifgd=1,nfgd
       call xred2xcart(1,rprimd,rcart,rred(:,ifgd)) 
       do mu=1,3
         do nu=1,3
           dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)&
&           +(raux3(ifgd)-raux2(ifgd)/rr_tmp(iindex(ifgd))**2)&
&           *rcart(mu)*rcart(nu)&
&           +raux2(ifgd)*rmet(mu,nu)
         end do
       end do
     end do
   end if !nfgd>0

!  Calculate contributions from |r-R|=0
   if(nfgd_r0>0) then
!    Evaluate spline fit of 2nd der of core charge density
!    from tcoredens(:,3) and tcoredens(:5)
     rr_tmp(1)=tol10
     call splint(msz,core_mesh%rad,&
&     pawtab%tcoredens(:,3),pawtab%tcoredens(:,5),&
&     1,rr_tmp(1),raux(1))

     factor=ucvol/real(nfft,dp)*raux(1)*rshpm1**2
     raux2(1:nfgd_r0)=factor*vxc(ifftsph(ncmax:ncmax-nfgd_r0+1),1)
     do ifgd=1,nfgd_r0
       do mu=1,3
         do nu=1,3
           dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)+raux2(ifgd)*rmet(mu,nu)
         end do
       end do
     end do
   end if !nfgd_r0>0
 end if !option==4


end subroutine mkcore_inner
!!***
