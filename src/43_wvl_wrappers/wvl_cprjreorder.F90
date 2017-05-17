!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_cprjreorder
!! NAME
!!  wvl_cprjreorder
!!
!! FUNCTION
!! Change the order of a wvl-cprj datastructure
!!   From unsorted cprj to atom-sorted cprj (atm_indx=atindx)
!!   From atom-sorted cprj to unsorted cprj (atm_indx=atindx1)
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2017 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atm_indx(natom)=index table for atoms
!!   From unsorted wvl%paw%cprj to atom-sorted wvl%paw%cprj (atm_indx=atindx)
!!   From atom-sorted wvl%paw%cprj to unsorted wvl%paw%cprj (atm_indx=atindx1)
!!
!! OUTPUT
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cprj_clean,cprj_paw_alloc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_cprjreorder(wvl,atm_indx)
    
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
#if defined HAVE_BIGDFT
 use BigDFT_API,only : cprj_objects,cprj_paw_alloc,cprj_clean
 use dynamic_memory
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_cprjreorder'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: atm_indx(:)
 type(wvl_internal_type),intent(inout),target :: wvl

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: iexit,ii,jj,kk,n1atindx,n1cprj,n2cprj,ncpgr
 character(len=100) :: msg
!arrays
 integer,allocatable :: nlmn(:)
 type(cprj_objects),pointer :: cprj(:,:)
 type(cprj_objects),allocatable :: cprj_tmp(:,:)
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
 cprj => wvl%paw%cprj

 n1cprj=size(cprj,dim=1);n2cprj=size(cprj,dim=2)
 n1atindx=size(atm_indx,dim=1)
 if (n1cprj==0.or.n2cprj==0.or.n1atindx<=1) return
 if (n1cprj/=n1atindx) then
   msg='wrong sizes!'
   MSG_BUG(msg)
 end if

!Nothing to do when the atoms are already sorted
 iexit=1;ii=0
 do while (iexit==1.and.ii<n1atindx)
   ii=ii+1
   if (atm_indx(ii)/=ii) iexit=0
 end do
 if (iexit==1) return

 ABI_ALLOCATE(nlmn,(n1cprj))
 do ii=1,n1cprj
   nlmn(ii)=cprj(ii,1)%nlmn
 end do
 ncpgr=cprj(1,1)%ncpgr

 ABI_DATATYPE_ALLOCATE(cprj_tmp,(n1cprj,n2cprj))
 call cprj_paw_alloc(cprj_tmp,ncpgr,nlmn) 
 do jj=1,n2cprj
   do ii=1,n1cprj
     cprj_tmp(ii,jj)%nlmn=nlmn(ii)
     cprj_tmp(ii,jj)%ncpgr=ncpgr
     cprj_tmp(ii,jj)%cp(:,:)=cprj(ii,jj)%cp(:,:)
     if (ncpgr>0) cprj_tmp(ii,jj)%dcp(:,:,:)=cprj(ii,jj)%dcp(:,:,:)
   end do
 end do

 call cprj_clean(cprj)

 do jj=1,n2cprj
   do ii=1,n1cprj
     kk=atm_indx(ii)
     cprj(kk,jj)%nlmn=nlmn(ii)
     cprj(kk,jj)%ncpgr=ncpgr
     cprj(kk,jj)%cp=f_malloc_ptr((/2,nlmn(ii)/),id='cprj%cp')
     cprj(kk,jj)%cp(:,:)=cprj_tmp(ii,jj)%cp(:,:)
     if (ncpgr>0) then
       cprj(kk,jj)%dcp=f_malloc_ptr((/2,ncpgr,nlmn(ii)/),id='cprj%dcp')
       cprj(kk,jj)%dcp(:,:,:)=cprj_tmp(kk,jj)%dcp(:,:,:)
     end if
   end do
 end do

 call cprj_clean(cprj_tmp)
 ABI_DATATYPE_DEALLOCATE(cprj_tmp)
 ABI_DEALLOCATE(nlmn)

#else
 if (.false.) write(std_out,*) atm_indx(1),wvl%h(1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_cprjreorder
!!***
