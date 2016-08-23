!{\src2tex{textfont=tt}}
!!****f* ABINIT/prttagm_images
!!
!! NAME
!! prttagm_images
!!
!! FUNCTION
!! Extension to prttagm to include the printing of
!! images information, in those cases the same variable
!! is printed several times for each dataset 
!!
!! Cases where images information are relevant includes
!! xcart, xred, acell, fcart.
!!
!! INPUT
!! (see prttagm.F90)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar_a_h,outvar_i_n,outvar_o_z
!!
!! CHILDREN
!!      appdig,prttagm,write_var_netcdf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prttagm_images(dprarr_images,iout,jdtset_,length,&
& marr,narrm,ncid,ndtset_alloc,token,typevarphys,&
& mxnimage,nimage,ndtset,prtimg,strimg,firstchar,forceprint)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prttagm_images'
 use interfaces_32_util
 use interfaces_57_iovars, except_this_one => prttagm_images
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,length,marr,ndtset_alloc,ncid
 integer,intent(in) :: mxnimage,ndtset
 integer,intent(in),optional :: forceprint
 character(len=*),intent(in) :: token
 character(len=3),intent(in) :: typevarphys
 character(len=1),intent(in),optional :: firstchar
!arrays
 integer,intent(in) :: prtimg(mxnimage,0:ndtset_alloc)
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 integer,intent(in) :: nimage(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxnimage)
 integer,intent(in) :: narrm(0:ndtset_alloc)
 real(dp),intent(in) :: dprarr_images(marr,mxnimage,0:ndtset_alloc)

!Local variables-------------------------------
 integer :: iarr,idtset,iimage,jdtset,multi_narr,narr
 integer :: intarr_images(marr,mxnimage,0:ndtset_alloc)
 integer,allocatable :: intarr(:,:)
 real(dp), allocatable :: dprarr(:,:)
 logical :: print_out,print_netcdf,test_multiimages
 character(len=1) :: first_column
 character(len=4) :: appen
 character(len=16) :: keywd
 character(len=50) :: full_format
 character(len=*), parameter :: format_1  ='",a16,t22,'
 character(len=*), parameter :: format_1a ='",a16,a,t22,'
 character(len=*), parameter :: format_2  ='",t22,'
 character(len=*), parameter :: long_dpr  ='3es18.10)'
!character(len=*), parameter :: format01160 ="(1x,a16,1x,(t22,3es18.10)) "
!character(len=*), parameter :: format01160a="(1x,a16,a,1x,(t22,3es18.10)) "

! *************************************************************************

 test_multiimages=.false.
 do idtset=1,ndtset_alloc
   if(nimage(idtset)>1)then
     do iarr=1,narrm(idtset)
       if(sum(abs( dprarr_images(iarr,2:nimage(idtset),idtset)- &
&       dprarr_images(iarr,1              ,idtset)))>tol12)then
         test_multiimages=.true.
       end if
     end do
   end if
 end do

 if(nimage(0)==0)test_multiimages=.true.

!DEBUG
!if(trim(token)=='vel')then
!write(ab_out,*)' test_multiimages=',test_multiimages
!endif
!ENDDEBUG

 if(.not.test_multiimages)then

   narr=narrm(1)
   ABI_ALLOCATE(intarr,(marr,0:ndtset_alloc))
   ABI_ALLOCATE(dprarr,(marr,0:ndtset_alloc))
   do idtset=0,ndtset_alloc
     dprarr(1:narrm(idtset),idtset)=dprarr_images(1:narrm(idtset),1,idtset)

!    DEBUG
!    if(trim(token)=='vel')then
!    write(ab_out,*)' idtset,narrm(idtset),dprarr(1:narrm(idtset),idtset)=',&
!    &    idtset,narrm(idtset),dprarr(1:narrm(idtset),idtset)
!    endif
!    ENDDEBUG

   end do
   multi_narr=0
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(narrm(1)/=narrm(idtset))multi_narr=1
     end do
   end if
!  if(narrm(0)==0)multi_narr=1
!  DEBUG
!  if(trim(token)=='fcart')then
!  write(std_out,*)' will call prttagm with fcart '
!  write(std_out,*)' narrm(0:ndtset_alloc)=',narrm(0:ndtset_alloc)
!  write(std_out,*)' multi_narr=',multi_narr
!  do idtset=0,ndtset_alloc
!  write(std_out,*)' dprarr_images(1:narrm(idtset),1,idtset)=',dprarr_images(1:narrm(idtset),1,idtset)
!  enddo
!  endif
!  ENDDEBUG
   if (present(firstchar).and.present(forceprint)) then
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
&     narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr,&
&     firstchar=firstchar,forceprint=forceprint)
   else if (present(firstchar)) then
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
&     narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr,&
&     firstchar=firstchar)
   else if (present(forceprint)) then
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
&     narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr,&
&     forceprint=forceprint)
   else
     call prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,&
&     narrm,ncid,ndtset_alloc,token,typevarphys,multi_narr)
   end if
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)

 else

   first_column=' ';if (present(firstchar)) first_column=firstchar

   do idtset=1,ndtset_alloc

     if (narrm(idtset)>0)then
       do iimage=1,nimage(idtset)

         print_out=.true.
         if (prtimg(iimage,idtset)==0) print_out=.false.
         if (nimage(0)>=nimage(idtset)) then
           if (sum(abs(dprarr_images(1:narrm(idtset),iimage,idtset) &
&           -dprarr_images(1:narrm(idtset),iimage,0)))<tol12) print_out=.false.
         end if
         print_netcdf=print_out
         
         if (present(forceprint)) then
           if (forceprint==1.or.forceprint==3) print_out=.true.
           if (forceprint==1.or.forceprint==2) print_netcdf=.true.
         end if

         if (print_out.or.print_netcdf.or.(ncid<0))then
           keywd=token//trim(strimg(iimage))

           if(ndtset>0)then
             jdtset=jdtset_(idtset)
             call appdig(jdtset,'',appen)
             if (print_out) then
               full_format='("'//first_column//trim(format_1a)//'("'// &
&               first_column//trim(format_2)//trim(long_dpr)//")"
               write(iout,full_format) &
&               trim(keywd),appen,dprarr_images(1:narrm(idtset),iimage,idtset)
             end if        
#ifdef HAVE_TRIO_NETCDF
             if (print_netcdf) then
               call write_var_netcdf(intarr_images(1:narrm(idtset),iimage,idtset),&
&               dprarr_images(1:narrm(idtset),iimage,idtset),&
&               marr,narrm(idtset),ncid,'DPR',trim(keywd)//appen)
             end if
#endif
           else

             if (print_out) then
               full_format='("'//first_column//trim(format_1)//'("'// &
&               first_column//trim(format_2)//trim(long_dpr)//")"
               write(iout,full_format) &
&               trim(keywd),dprarr_images(1:narrm(idtset),iimage,idtset)
             end if
#ifdef HAVE_TRIO_NETCDF
             if (print_netcdf) then
               call write_var_netcdf(intarr_images(1:narrm(idtset),iimage,idtset),&
&               dprarr_images(1:narrm(idtset),iimage,idtset),&
&               marr,narrm(idtset),abs(ncid),'DPR',trim(keywd))
             end if
#endif
             
           end if
         end if
       end do
     end if
   end do

 end if

end subroutine prttagm_images
!!***
