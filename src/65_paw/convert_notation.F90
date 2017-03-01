!{\src2tex{textfont=tt}}
!!****f* ABINIT/convert_notation
!! NAME
!!  convert_notation
!!
!! FUNCTION
!!  notation: Convert notation index of the elastic tensor to use the  
!!            second derivative of gylm
!!            for example :  to calculate d2(gylm)/d(eps_32)
!!                          - 
!!                          - voigt notation       => 32   
!!                          - normal notation      => 3 3 2 2
!!                          - notation for gylmgr2 => 32 32 32 32 => 4 4 4   
!!            The calculation of the second derivative of gylm wrt strains requires  four
!!            derivative (see the formula)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2017 ABINIT group (AM)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  eps_alpha 
!!  eps_beta  
!!  eps_delta 
!!  eps_gamma 
!!
!! OUTPUT
!!  mu4(4) = array with index for the second derivative of gylm
!!
!! SIDE EFFECTS
!!  return the four index in mu4 for the calculation of the second derivative of gylm
!!
!!
!! NOTES
!!
!! PARENTS
!!      pawgrnl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine convert_notation(mu4,eps_alpha,eps_beta,eps_gamma,eps_delta)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'convert_notation'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalar
 integer,intent(in)  :: eps_alpha,eps_beta
 integer,optional,intent(in)  :: eps_gamma,eps_delta
 !array
 integer,intent(inout) :: mu4(4)

!Local variables-------------------------------
 integer :: eps1,eps2,i,j,k      
 integer,allocatable :: mu_temp(:)

! *************************************************************************
 
!DBG_ENTER("COLL")

 ABI_ALLOCATE(mu_temp,(4))
 if (present(eps_gamma).and.present(eps_delta)) then
   mu_temp(1)=eps_alpha
   mu_temp(2)=eps_beta
   mu_temp(3)=eps_gamma
   mu_temp(4)=eps_delta
 else
   mu_temp(1)=eps_alpha
   mu_temp(2)=eps_beta
   mu_temp(3)= 0
   mu_temp(4)= 0
 end if
 k=1
 do i=1,2
   eps1=mu_temp(i)
   do j=1,2
     eps2=mu_temp(2+j)
     if(eps1==eps2) then
       if(eps1==1) mu4(k)=1;
       if(eps1==2) mu4(k)=2;
       if(eps1==3) mu4(k)=3;
     else
       if((eps1==3.and.eps2==2).or.(eps1==2.and.eps2==3)) mu4(k)=4;
       if((eps1==3.and.eps2==1).or.(eps1==1.and.eps2==3)) mu4(k)=5;
       if((eps1==1.and.eps2==2).or.(eps1==2.and.eps2==1)) mu4(k)=6;
     end if
     k=k+1
   end do
 end do

 ABI_DEALLOCATE(mu_temp)

!DBG_EXIT("COLL")

end subroutine convert_notation
!!***
