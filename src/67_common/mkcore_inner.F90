!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkcore_inner
!! NAME
!!  mkcore_inner
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (T. Rangel)
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
!!      paw_splint,paw_splint_der,sort_dp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkcore_inner(corfra,core_mesh,dyfrx2,grxc1,grxc2,grxc3,ifftsph,msz,natom,ncmax,nfft,&
&          nfgd,nfgd_r0,nspden,n3xccc,option,pawtab,rmet,rprimd,rr,strdia,vxc,xccc3d,&
&          rred) ! optional argument

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_sort,   only : sort_dp
 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type
 use m_paw_numeric, only : paw_splint,paw_splint_der

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore_inner'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer ,intent(in) :: msz,natom,ncmax,nfft,nfgd,nfgd_r0,nspden,n3xccc,option
 real(dp),intent(out) :: grxc1,grxc2,grxc3,strdia
!arrays
 integer,intent(in) :: ifftsph(ncmax)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),rr(ncmax),vxc(nfft,nspden)
 real(dp),intent(in),optional :: rred(:,:)
 real(dp),intent(inout) :: corfra(3,3),dyfrx2(3,3,natom),xccc3d(n3xccc)
 type(pawrad_type),intent(in) :: core_mesh
 type(pawtab_type),intent(in) :: pawtab

!Local variables-------------------------------
!scalars
 integer :: iatom,ifgd,ii,jj,mu,nu
 character(len=500) :: message
 real(dp) :: term,term2
!arrays
 integer :: iindex(nfgd)
 real(dp) :: tcore(nfgd),dtcore(nfgd),rr_tmp(nfgd)
 real(dp),allocatable :: d2tcore(:)

! *************************************************************************

!Checks
 if(nspden >1) then
   write(message, '(a)')'mkcore_inner: this is not yet generalized to npsden>1'
   MSG_ERROR(message)
 end if
 if (present(rred)) then
   if (option>1.and.size(rred)/=3*ncmax) then
     write(message, '(a)')'mkcore_inner: incorrect size for rred'
     MSG_BUG(message)
   end if
 else if (option>1) then
   write(message, '(a)')'mkcore_inner: rred is not present and option >1'
   MSG_BUG(message)
 end if

!Retrieve values of pseudo core density (and derivative)
 rr_tmp=rr(1:nfgd)
 iindex(1:nfgd)=(/(ii,ii=1,nfgd)/)
 call sort_dp(nfgd,rr_tmp,iindex(1:nfgd),tol16)
 if (option==1.or.option==3) then
   call paw_splint(msz,core_mesh%rad,pawtab%tcoredens(:,1),pawtab%tcoredens(:,3),&
&   nfgd,rr_tmp,tcore)
 end if
 if (option>=2) then
   call paw_splint_der(msz,core_mesh%rad,pawtab%tcoredens(:,1),pawtab%tcoredens(:,3),&
&   nfgd,rr_tmp,dtcore)
 end if

!Accumulate contributions to core density on the entire cell
 if (option==1) then
   do ifgd=1,nfgd
     ii=iindex(ifgd);jj=ifftsph(ii)
     xccc3d(jj)=xccc3d(jj)+tcore(ifgd)
   end do

!Accumulate contributions to Exc gradients
 else if (option==2) then
   do ifgd=1,nfgd
     ii=iindex(ifgd);jj=ifftsph(ii)
     term=vxc(jj,1)*dtcore(ifgd)/rr_tmp(ifgd)
     grxc1=grxc1-rred(1,ii)*term
     grxc2=grxc2-rred(2,ii)*term
     grxc3=grxc3-rred(3,ii)*term
   end do

!Accumulate contributions to stress tensor
 else if (option==3) then
!  Write out the 6 symmetric components
   do ifgd=1,nfgd
     ii=iindex(ifgd);jj=ifftsph(ii)
     term=vxc(jj,1)*dtcore(ifgd)/rr_tmp(ifgd)
     corfra(1,1)=corfra(1,1)+term*rred(1,ifgd)**2
     corfra(2,2)=corfra(2,2)+term*rred(2,ifgd)**2
     corfra(3,3)=corfra(3,3)+term*rred(3,ifgd)**2
     corfra(3,2)=corfra(3,2)+term*rred(3,ifgd)*rred(2,ifgd)
     corfra(3,1)=corfra(3,1)+term*rred(3,ifgd)*rred(1,ifgd)
     corfra(2,1)=corfra(2,1)+term*rred(2,ifgd)*rred(1,ifgd)
   end do
!  Also compute diagonal term
   do ifgd=1,nfgd
     ii=iindex(ifgd);jj=ifftsph(ii)
     strdia=strdia+vxc(jj,1)*tcore(ii)
   end do

!Compute frozen-wf contribution to Dynamical matrix
 else if (option==4) then
   ABI_ALLOCATE(d2tcore,(nfgd))
   if(nfgd>0) then
!    Evaluate spline fit of 2nd der of pseudo core density
     call paw_splint(msz,core_mesh%rad,pawtab%tcoredens(:,3),pawtab%tcoredens(:,5),&
&     nfgd,rr_tmp,d2tcore)
     do ifgd=1,nfgd
       ii=iindex(ifgd);jj=ifftsph(ii)
       term=vxc(jj,1)*dtcore(ifgd)/rr_tmp(ifgd)
       term2=term*d2tcore(ifgd)/rr_tmp(ifgd)
       do mu=1,3
         do nu=1,3
           dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)&
&           +(term2-term/rr_tmp(iindex(ifgd))**2)&
&           *rred(mu,iatom)*rred(nu,iatom)+term*rmet(mu,nu)
         end do
       end do
     end do
   end if
!  Contributions from |r-R|=0
   if (nfgd_r0>0) then
     rr_tmp(1)=tol10
     call paw_splint(msz,core_mesh%rad,pawtab%tcoredens(:,3),pawtab%tcoredens(:,5),&
&     1,rr_tmp,d2tcore(1))
     ifgd=1
     ii=iindex(ifgd);jj=ifftsph(ii)
     term=vxc(jj,1)*d2tcore(ifgd)
     do mu=1,3
       do nu=1,3
         dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)+term*rmet(mu,nu)
       end do
     end do
   end if
   ABI_DEALLOCATE(d2tcore)

 end if !option

end subroutine mkcore_inner
!!***
