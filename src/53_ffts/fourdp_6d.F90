!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourdp_6d
!! NAME
!! fourdp_6d
!!
!! FUNCTION
!!     calculate a 6-dimensional Fast Fourier Transform
!!
!!     isign=-1 : A(G1,G2) = Sum(r1,r2) A(r1,r2) exp(-iG1.r1) exp(+iG2.r2)
!!                                                  ^            ^
!!     isign=+1 : A(r1,r2) = Sum(G1,G2) A(G1,G2) exp(+iG1.r1) exp(-iG2.r2)
!!                                                  ^            ^
!!     isign=-1 and isign=1 form a transform/inverse-transform pair: calling
!!     one then the other will take you back to the original function,
!!     multiplied by a factor of (nl*nm*nn)**2.
!!     ------------------------------------------------------------------
!!
!!     input:
!!      a: A(r1,r2) [overwritten]
!!     output:
!!      a: A(G1,G2) in the format IGFFT
!!     ------------------------------------------------------------------
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_kxc
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fourdp_6d(cplex,matrix,isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
    
 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourdp_6d'
 use interfaces_53_ffts, except_this_one => fourdp_6d
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 complex(gwpc),intent(inout) :: matrix(nfft,nfft)

!Local variables-------------------------------
!scalars
 !integer,parameter :: cplex=2
 integer :: i1,i2,i3,ifft
 integer :: n1,n2,n3
!arrays
 real(dp),allocatable :: fofg(:,:),fofr(:)
 
! *************************************************************************

!TODO check normalization factor, it is better if we use the GW conventions.
 n1 = ngfft(1)
 n2 = ngfft(2)
 n3 = ngfft(3)

 ABI_MALLOC(fofg,(2,nfft))
 ABI_MALLOC(fofr,(cplex*nfft))

 do i3=0,n3-1
   do i2=0,n2-1
     do i1=0,n1-1

       ifft=1+i1+i2*n1+i3*n1*n2
       if (isign==1) then ! G1 -> r1 transform for each G2 to form A(r1,G2)
         fofg(1,:)=REAL (matrix(:,ifft))
         fofg(2,:)=AIMAG(matrix(:,ifft))
       else if (isign==-1) then ! r1 -> G1 transform for each r2 to form A(G1,r2)
         fofr(1:nfft)       =REAL (matrix(:,ifft))
         fofr(nfft+1:2*nfft)=AIMAG(matrix(:,ifft))
       else 
         MSG_ERROR("Wrong isign")
       end if

       call fourdp(cplex,fofg,fofr,isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

       if (isign==1) then ! Save A(r1,G2)
         matrix(:,ifft)=CMPLX(fofr(1:nfft),fofr(nfft+1:2*nfft))
       else if (isign==-1) then ! Save A(G1,r2)
         matrix(:,ifft)=CMPLX(fofg(1,:),fofg(2,:))
       end if

     end do
   end do
 end do

 do i3=0,n3-1
   do i2=0,n2-1
     do i1=0,n1-1

       ifft=1+i1+i2*n1+i3*n1*n2
       if (isign==1) then ! Do the G2 -> r2 transform of A(r1,G2) to get A(r1,r2)
         fofr(1:nfft       )=REAL (matrix(ifft,:))
         fofr(nfft+1:2*nfft)=AIMAG(matrix(ifft,:))
       else if (isign==-1) then
!        Do the r2 -> G2 transform of A(G1,r2) to get A(G1,G2)
         fofg(1,:)=REAL (matrix(ifft,:))
         fofg(2,:)=AIMAG(matrix(ifft,:))
       end if

       call fourdp(2,fofg,fofr,-isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

       if (isign==1) then
         matrix(ifft,:)=CMPLX(fofg(1,:),fofg(2,:))
       else if (isign==-1) then
         matrix(ifft,:)=CMPLX(fofr(1:nfft),fofr(nfft+1:2*nfft))
       end if

     end do
   end do
 end do

 ABI_FREE(fofg)
 ABI_FREE(fofr)

end subroutine fourdp_6d
!!***
