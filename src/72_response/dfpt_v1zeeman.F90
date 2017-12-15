!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_v1zeeman
!! NAME
!!  dfpt_v1zeeman
!!
!! FUNCTION
!!  Calculate 1st order Zeeman potential = -vec{\sigma}.\vec{b}, where 
!!  sigma is the vector of Pauli matrices and \vec{b} is the unit 
!!  vector indicating the perturbing field direction.
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nspden = number of density matrix components
!!  nfft   = numbder of fft grid points
!!  cplex  = complex or real density matrix
!!  idir   = direction of the perturbing field in Cartesian frame
!!           1: along x
!!           2: along y
!!           3: along z
!!           4: identity matrix at each fft point is returned (for density-density response)
!!
!! OUTPUT
!!  v1zeeman(nfft*cplex,nspden)= 1st order Zeeman potential, or Identity matrix (electrostatic potential) for idir=4
!!
!! SIDE EFFECTS
!!
!!  None
!!
!! NOTES
!!
!!  The definition of components of the potential matrix differ depending on cplex
!!  for nspden=4:
!!  For cplex=1, the potential is defined as (V_upup,V_dndn,Re[V_updn],Im[V_updn])
!!  For cplex=2, the definition is (V_upup,V_dndn,V_updn,i.V_updn)
!!
!! PARENTS
!!
!!  dfpt_rhotov
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_v1zeeman(nspden,nfft,cplex,idir,v1zeeman)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_v1zeeman'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)    :: idir,nfft,cplex,nspden
 real(dp), intent(inout) :: v1zeeman(cplex*nfft,nspden)                       

!Local variables-------------------------------
 integer :: ifft                                
!character(len=500) :: msg                   
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

 select case(cplex)
 case(1)
   if (nspden==4) then
     if(idir==3)then       ! Zeeman field along the 3rd axis (z)   
       v1zeeman(:,1)=-0.5d0
       v1zeeman(:,2)=+0.5d0
       v1zeeman(:,3)= 0.0d0
       v1zeeman(:,4)= 0.0d0
     else if(idir==2)then  ! Zeeman field along the 2nd axis (y)
       v1zeeman(:,1)= 0.0d0
       v1zeeman(:,2)= 0.0d0
       v1zeeman(:,3)= 0.0d0
       v1zeeman(:,4)=+0.5d0
     else                  ! Zeeman field along the 1st axis (x)
       v1zeeman(:,1)= 0.0d0
       v1zeeman(:,2)= 0.0d0
       v1zeeman(:,3)=-0.5d0
       v1zeeman(:,4)= 0.0d0
     end if
   else if (nspden==2) then
     v1zeeman(:,1)=-0.5e0
     v1zeeman(:,2)= 0.5e0
   else
     v1zeeman(:,1)= 0.0e0
   end if
 case(2)
   if (nspden==2) then
     do ifft=1,nfft
       v1zeeman(2*ifft-1,1)  =-0.5e0
       v1zeeman(2*ifft  ,1)  = 0.0e0
       v1zeeman(2*ifft-1,2)  = 0.5e0
       v1zeeman(2*ifft  ,2)  = 0.0e0
     enddo
   else if (nspden==4) then
     select case(idir)
     case(1) !along x, v1=-sigma_x
       do ifft=1,nfft
         v1zeeman(2*ifft-1,1)= 0.0e0 !Re[V^11]
         v1zeeman(2*ifft  ,1)= 0.0e0 !Im[V^11]
         v1zeeman(2*ifft-1,2)= 0.0e0 !Re[V^22]
         v1zeeman(2*ifft  ,2)= 0.0e0 !Im[V^22]
         v1zeeman(2*ifft-1,3)=-0.5e0 !Re[V^12]
         v1zeeman(2*ifft  ,3)= 0.0e0 !Im[V^12]
         v1zeeman(2*ifft-1,4)= 0.0e0 !Re[i.V^21]=Im[V^12]
         v1zeeman(2*ifft  ,4)=-0.5e0 !Im[i.V^21]=Re[V^12]
       enddo
     case(2) !along y, v1 = -sigma_y
       do ifft=1,nfft
         v1zeeman(2*ifft-1,1)= 0.0e0 !Re[V^11]
         v1zeeman(2*ifft  ,1)= 0.0e0 !Im[V^11]
         v1zeeman(2*ifft-1,2)= 0.0e0 !Re[V^22]
         v1zeeman(2*ifft  ,2)= 0.0e0 !Im[V^22]
         v1zeeman(2*ifft-1,3)= 0.0e0 !Re[V^12]
         v1zeeman(2*ifft  ,3)=+0.5e0 !Im[V^12]
         v1zeeman(2*ifft-1,4)=+0.5e0 !Re[i.V^21]=Im[V^12]
         v1zeeman(2*ifft  ,4)= 0.0e0 !Im[i.V^21]=Re[V^12]
       enddo
     case(3)
       do ifft=1,nfft
         v1zeeman(2*ifft-1,1)=-0.5e0 !Re[V^11]
         v1zeeman(2*ifft  ,1)= 0.0e0 !Im[V^11]
         v1zeeman(2*ifft-1,2)= 0.5e0 !Re[V^22]
         v1zeeman(2*ifft  ,2)= 0.0e0 !Im[V^22]
         v1zeeman(2*ifft-1,3)= 0.0e0 !Re[V^12]
         v1zeeman(2*ifft  ,3)= 0.0e0 !Im[V^12]
         v1zeeman(2*ifft-1,4)= 0.0e0 !Re[i.V^21]
         v1zeeman(2*ifft  ,4)= 0.0e0 !Im[i.V^21]
       enddo
     end select
   endif
 end select !cplex

end subroutine dfpt_v1zeeman
!!***
