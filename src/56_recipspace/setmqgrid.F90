!{\src2tex{textfont=tt}}
!!****f* ABINIT/setmqgrid
!! NAME
!!  setmqgrid
!!
!! FUNCTION
!!  Sets the number of points needed to represent the pseudopotentials in
!!  reciprocal space for a specified resolution.
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2016 ABINIT group (DW)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ecut=cutoff energy for the wavefunctions
!!  ecutdg=cutoff energy for the fine grid in case usepaw==1
!!  gprimd=primitive translation vectors for reciprocal space
!!  nptsgvec=number of points along the smallest primitive translation vector
!!    of the reciprocal space
!!  usepaw=1 if PAW is used, 0 otherwise
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_psps,memory_eval
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine setmqgrid(mqgrid,mqgriddg,ecut,ecutdg,gprimd,nptsgvec,usepaw)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setmqgrid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(inout)  :: mqgrid,mqgriddg
 integer , intent(in)  :: nptsgvec,usepaw
 real(dp), intent(in) :: ecut,ecutdg
 real(dp), intent(in) :: gprimd(3,3)

!Local variables-------------------------------
 integer :: mqgrid2,mqgriddg2
 real(dp) :: gmax,gmaxdg,gvecnorm
 character(len=500) :: message  
 
! *************************************************************************
 
 gvecnorm=sqrt(min(dot_product(gprimd(:,1),gprimd(:,1)), &
& dot_product(gprimd(:,2),gprimd(:,2)), &
& dot_product(gprimd(:,3),gprimd(:,3))))
 gmax=one/(sqrt2*pi)*sqrt(ecut)

 if (mqgrid == 0) then
   mqgrid2=ceiling(gmax/gvecnorm*nptsgvec)
   mqgrid=max(mqgrid2,3001)
   write(message, '(5a,i0,a)' )&
&   'The number of points "mqgrid" in reciprocal space used for the',ch10,&
&   'description of the pseudopotentials has been set automatically',ch10,&
&   'by abinit to: ',mqgrid,'.'
   !MSG_COMMENT(message)
 else
   mqgrid2=ceiling(gmax/gvecnorm*nptsgvec)
   if (mqgrid2>mqgrid) then
     write(message, '(3a,i8,3a,i8,3a)' )&
&     'The number of points "mqgrid" in reciprocal space used for the',ch10,&
&     'description of the pseudopotentials is : ',mqgrid,'.',ch10,&
&     'It would be better to increase it to at least ',mqgrid2,', or',ch10,&
&     'let abinit choose it automatically by setting mqgrid = 0.'
     MSG_WARNING(message)
   end if
 end if

 if (usepaw==1) then
   gmaxdg=one/(sqrt2*pi)*sqrt(ecutdg)
   if (mqgriddg == 0) then
     mqgriddg2=ceiling(gmaxdg/gvecnorm*nptsgvec)
     mqgriddg=max(mqgriddg2,3001)
     write(message, '(5a,i0,a)' )&
&     'The number of points "mqgriddg" in reciprocal space used for the',ch10,&
&     'description of the pseudopotentials has been set automatically',ch10,&
&     'by abinit to: ',mqgriddg,'.'
     !MSG_COMMENT(message)
   else
     mqgriddg2=ceiling(gmax/gvecnorm*nptsgvec)
     if (mqgriddg2>mqgriddg) then
       write(message, '(3a,i8,3a,i8,3a)' )&
&       'The number of points "mqgriddg" in reciprocal space used for the',ch10,&
&       'description of the pseudopotentials (fine grid) is :',mqgriddg,'.',ch10,&
&       'It would be better to increase it to at least ',mqgriddg2,', or',ch10,&
&       'let abinit choose it automatically by setting mqgrid = 0.'
       MSG_WARNING(message)
     end if
   end if
 end if

end subroutine setmqgrid
!!***
