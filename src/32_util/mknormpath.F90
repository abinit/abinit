!{\src2tex{textfont=tt}}
!!****f* ABINIT/mknormpath
!! NAME
!! mknormpath
!!
!! FUNCTION
!! Please do not use this  routine, use make_normpath instead.
!! mknormpath should be removed
!!
!!  This simple routine generates a normalized path that can be used to plot a band 
!!  structures in an easy way. For normalized path we mean a path where the number 
!!  of division on each segment is proportional to the length of the segment itself. 
!!  To generate the above mentioned path, the subroutine must be called twice. 
!!  The first call reports the total number of divisions in the normalized path, dimension
!!  that is required to correctly allocate the array.  
!!  The second call calculates the reduced coordinates of the circuit.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2016 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nbounds=number of points defining the path
!! ndiv_small=number of points to be used to sample the smallest
!!  segment defined by bounds(:,1:nbounds)
!! bounds(3,nbounds)=points defining the path
!! gmet(3,3)=metric 
!!
!! OUTPUT
!! ndiv(nbounds-1)= number of divisions for each segment
!! npt_tot=total number of points sampled along the circuit
!! path(3,npt_tot)= normalized path in reciprocal space
!!
!! TODO 
!!  Do not use this routine, it is obsolete and should be replaced by make_path in m_bz_mesh.
!!
!! PARENTS
!!      inkpts
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mknormpath(nbounds,bounds,gmet,ndiv_small,ndiv,npt_tot,path)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mknormpath'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !F95 construct, interface required but we can call mknormpath once
 !real(dp),pointer :: path(:,:) 
!scalars
 integer,intent(in) :: nbounds,ndiv_small
 integer,intent(inout) :: npt_tot
!arrays
 integer,intent(inout) :: ndiv(nbounds-1)
 real(dp),intent(in) :: bounds(3,nbounds),gmet(3,3)
 real(dp),intent(out),optional :: path(3,npt_tot)

!Local variables-------------------------------
!scalars
 integer :: idx,ii,jp
 real(dp) :: fct
 character(len=500) :: message
!arrays
 real(dp) :: dd(3),lng(nbounds-1)

! *************************************************************************
 
 if (ndiv_small<=0) then
   write(message,'(3a,i0)')&
&   'The argument ndiv_small should be a positive number,',ch10,&
&   'however, ndiv_small=',ndiv_small
   MSG_ERROR(message)
 end if

 do ii=1,nbounds-1 
   dd(:)=bounds(:,ii+1)-bounds(:,ii)
   lng(ii)= sqrt( dd(1)*gmet(1,1)*dd(1)+ &    
&   dd(2)*gmet(2,2)*dd(2)+ &
&   dd(3)*gmet(3,3)*dd(3)+ &
&   2.0d0*(dd(1)*gmet(1,2)*dd(2)+ &
&   dd(1)*gmet(1,3)*dd(3)+ &
&   dd(2)*gmet(2,3)*dd(3)) &
&   )
 end do
 write(std_out,*)lng
 fct=minval(lng)

!Avoid division by zero if k(:,i+1)=k(:,i)
 if (abs(fct)<tol6) then 
   write(message,'(3a)')&
&   'found two consecutive points in the path which are equal',ch10,&
&   'This is not allowed, please modify the path in your input file'
   MSG_ERROR(message)
 end if

 fct=fct/ndiv_small
 ndiv(:)=nint(lng(:)/fct) 
!The 1 stand for the first point
 npt_tot=sum(ndiv)+1

!allocate(path(3,npt_tot)
 if (.not.present(path)) then 
   write(message,'(2a,i8)')ch10,&
&   ' mknormpath : total number of points on the path : ',npt_tot
   call wrtout(std_out,message,'COLL')
   write(message,'(2a)')ch10,' Number of divisions for each segment of the normalized path : '
   call wrtout(std_out,message,'COLL') 
   do ii=1,nbounds-1
     write(message,'(2(3f8.5,a),i5,a)')&
     bounds(:,ii),' ==> ',bounds(:,ii+1),' ( ndiv : ',ndiv(ii),' )' 
     call wrtout(std_out,message,'COLL')
   end do 
   write(message,'(a)')ch10
   call wrtout(std_out,message,'COLL') 
 else 
   write(message,'(2a)')ch10,' Normalized Path : '
   call wrtout(std_out,message,'COLL')
   idx=1
   do ii=1,nbounds-1
     do jp=1,ndiv(ii)
       path(:,idx)=bounds(:,ii)+(jp-1)*(path(:,ii+1)-path(:,ii))/ndiv(ii)
       write(message,'(i4,4x,3(f8.5,1x))')idx,path(:,idx)
       call wrtout(std_out,message,'COLL')
       idx=idx+1
     end do 
   end do 
 end if

end subroutine mknormpath
!!***
