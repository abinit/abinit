!{\src2tex{textfont=tt}}
!!****f* ABINIT/green_kernel
!! NAME
!! green_kernel
!! 
!! FUNCTION
!! this routine compute the fourrier transform of the Green kernel for the 
!! recursion method  
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  inf_rmet=define the  infinitesimal metric : rprimd*(transpose(rprimd)) divided by the number of discretisation point
!!  inf_ucvol=volume of infinitesimal cell
!!  mult=variance of the Gaussian (=rtrotter/beta)
!!  mpi_enreg=information about MPI parallelization
!!  ngfft=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nfft=total number of fft grid points
!!  debug_rec=debugging variable
!! 
!! OUTPUT
!!  ZT_p=fourier transforme of the Green kernel
!!  
!! PARENTS
!!      first_rec
!!
!! CHILDREN
!!      fourdp,timab,wrtout
!!
!! NOTES 
!!  at this time :
!!       - need a rectangular box
!! 
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine green_kernel(ZT_p,inf_rmet,inf_ucvol,mult,mpi_enreg,ngfft,nfft)
 
 use defs_basis
 use defs_abitypes
 use m_profiling_abi

 use m_time,        only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'green_kernel'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none
 
!Arguments -------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: inf_ucvol,mult
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: inf_rmet(3,3)
 real(dp),intent(out) :: ZT_p(1:2,0:nfft-1)
 
!Local variables-------------------------------
!scalars
 integer,parameter :: n_green_max=5
 integer :: ii,isign,jj,kk,n_green,xx,yy,zz
 real(dp) :: acc, norme
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: T_p(:)

! *************************************************************************
 
 call timab(603,1,tsec)
 
 norme = (mult/pi)**(onehalf)

 ABI_ALLOCATE(T_p,(0:nfft-1))
 
!n_green should be better chosen for non rectangular cell
 do xx=1, n_green_max
   n_green = xx
   if(exp(-mult*dsq_green(xx*ngfft(1),0,0,inf_rmet))<tol14 &
&   .and. exp(-mult*dsq_green(0,xx*ngfft(2),0,inf_rmet))<tol14 &
&   .and. exp(-mult*dsq_green(0,0,xx*ngfft(3),inf_rmet))<tol14 ) exit
 end do
 
 acc = zero
 T_p = zero
 do kk = 0,ngfft(3)-1
   do jj = 0,ngfft(2)-1
     do ii = 0,ngfft(1)-1
       
       do xx=-n_green,n_green-1
         do yy=-n_green,n_green-1
           do zz=-n_green,n_green-1
             
             T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)+ & 
&             exp(-mult*dsq_green(ii+xx*ngfft(1),jj+yy*ngfft(2),kk+zz*ngfft(3),inf_rmet))
             
           end do
         end do
       end do
       
       T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = norme*T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)
       acc = acc + inf_ucvol* T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)
       
     end do
   end do
 end do
 
 T_p(:)= (one/acc)*T_p(:)

!if(debug_rec)then
 write(msg,'(a,d12.3,2(2a,i8),2(2a,3d12.3),2a,d16.6)')&
& ' on the boundary    ', exp(-mult*dsq_green(ngfft(1),0,0,inf_rmet)),ch10, &
& ' no zero            ', count(T_p>tol14),ch10, &
& ' n_green            ', n_green,ch10, &
& ' erreur_n_green     ', exp(-mult*dsq_green(n_green*ngfft(1),0,0,inf_rmet)), &
& exp(-mult*dsq_green(0,n_green*ngfft(2),0,inf_rmet)), &
& exp(-mult*dsq_green(0,0,n_green*ngfft(3),inf_rmet)),ch10,&
& ' erreur_troncat     ', T_p(ngfft(1)/2),  &
& T_p(ngfft(1)*(ngfft(2)/2)), &
& T_P(ngfft(1)*ngfft(2)*(ngfft(3)/2)),ch10, &
& ' erreurT_p          ',abs(acc-1.d0)
 call wrtout(std_out,msg,'COLL')
!endif

 
 isign = -1
 call fourdp(1,ZT_p,T_p,isign,mpi_enreg,nfft,ngfft,1,0)
 
 ABI_DEALLOCATE(T_p)

 ZT_p(:,:) = real(nfft,dp)*ZT_p
 

 call timab(603,2,tsec)

 contains

   function dsq_green(ii,jj,kk,inf_rmet)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dsq_green'
!End of the abilint section

   real(dp) :: dsq_green
   integer,intent(in) :: ii,jj,kk
   real(dp),intent(in) :: inf_rmet(3,3)
   dsq_green= inf_rmet(1,1)*dble(ii**2)&
&   +inf_rmet(2,2)*dble(jj**2)&
&   +inf_rmet(3,3)*dble(kk**2)&
&   +two*(inf_rmet(1,2)*dble(ii*jj)&
&   +inf_rmet(2,3)*dble(jj*kk)&
&   +inf_rmet(3,1)*dble(kk*ii))
 end function dsq_green
 
end subroutine green_kernel
!!***
