!{\src2tex{textfont=tt}}
!!****f* ABINIT/msig
!! NAME
!! msig
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (SMazevet,VR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  fcti(npti)=  conductivity, as calculated in conducti
!!  npti= number of points to calculate conductivity
!!  xi(npti)= energies where the conductivity is calculated
!!
!! OUTPUT
!!   no output, only files
!!
!! NOTES
!!     this program calculates the imaginary part of the conductivity (principal value)
!!     +derived optical properties. 
!!     the calculation is performed on the same grid as the initial input
!!     to calculate the principal value, a trapezoidale integration +taylor expansion to 
!!     third order is used (W.J. Thomson computer in physics vol 12 p94 1998)
!!    two input files are needed inppv.dat (parameters) and sigma.dat (energy,sigma_1)
!!     two output files ppsigma.dat (energy,sigma_1,sigma_2,epsilon_1,epsilon_2)
!!                      abs.dat     (energy,nomega,komega,romega,absomega)
!!     march 2002 s.mazevet
!!
!! PARENTS
!!      conducti_nc,conducti_paw
!!
!! CHILDREN
!!      intrpl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine msig(fcti,npti,xi,filnam_out_sig)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools, only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'msig'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: npti
!arrays
 real(dp),intent(in) :: fcti(npti),xi(npti)
 character(len=fnlen),intent(in) :: filnam_out_sig

!Local variables-------------------------------
!scalars
 integer,parameter :: npt=10000
 integer :: ii,ip,npt1,npt2,eps_unt,abs_unt
 real(dp),parameter :: del=0.001_dp,ohmtosec=9.d11
 real(dp) :: dx,dx1,dx2,eps1,eps2,idel,komega,pole,refl,sigma2,xsum
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: abso(:),fct(:),fct1(:),fct2(:),fct3(:),fct4(:),fct5(:)
 real(dp),allocatable :: fctii(:),fp(:),fpp(:),fppp(:),nomega(:),ppsig(:)
 real(dp),allocatable :: x1(:),x2(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

 write(std_out,'(2a)')ch10,'Calculate the principal value and related optical properties'
 write(std_out,'(a)')'following W.J. Thomson computer in physics vol 12 p94 1998 for '
 write(std_out,'(a)')'the principal value. S. Mazevet'
 write(std_out,'(a)')'OPTIONS'
 write(std_out,'(a)')'use default number of integration pts: npt=10000'
 write(std_out,'(a)')'Use default value for delta interval: del=1e-3'

 if (open_file(trim(filnam_out_sig)//'_eps',msg, newunit=eps_unt,status='replace',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(eps_unt,'(a)')'#energy (eV),sigma_1(Ohm-1cm-1),sigma_2(Ohm-1cm-1),epsilon_1,epsilon_2'

 if (open_file(trim(filnam_out_sig)//'_abs', msg, newunit=abs_unt, status='replace',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(abs_unt,'(a)')'#energy(eV),nomega,komega,refl.,abso.(cm-1)'

 ABI_ALLOCATE(fct,(npt))
 ABI_ALLOCATE(fct2,(npt))
 ABI_ALLOCATE(fct3,(npt))
 ABI_ALLOCATE(fct4,(npt))
 ABI_ALLOCATE(fct5,(npt))
 ABI_ALLOCATE(fp,(npt))
 ABI_ALLOCATE(fpp,(npt))
 ABI_ALLOCATE(fppp,(npt))
 ABI_ALLOCATE(x1,(npt))
 ABI_ALLOCATE(x2,(npt))
 ABI_ALLOCATE(fct1,(npt))
 ABI_ALLOCATE(ppsig,(npt))
 ABI_ALLOCATE(fctii,(npt))
 ABI_ALLOCATE(abso,(npt))
 ABI_ALLOCATE(nomega,(npt))

 if (npti > npt) then
   write (std_out,*) 'msig: input npti is too large for hard coded npt array size = ', npt
   MSG_ERROR("Aborting now")
 end if

!loop on the initial energy grid      
 do ip=1,npti

!  adjust the interval before and after the pole to reflect range/npt interval
   xsum=zero
   dx=(xi(npti)-xi(1))/dble(npt-1)
   pole=xi(ip)
   npt1=int((pole-del)/dx)
   dx1=zero
   if(npt1/=1) dx1=(pole-del)/(npt1-1)
   npt2=int((xi(npti)-pole-del)/dx)
   dx2=(xi(npti)-pole-del)/(npt2-1)

!  for the moment skip the pp calculation when the pole if too close to the end of the range
   if(npt1<=1.or.npt2<=1) then

     xsum=zero
     ppsig(ip)=zero

   else

!    define the fct for which the pp calculation is needed using xi^2-pole^2 factorization
     fctii(1:npti) = zero
     fctii(1:npti)=fcti(1:npti)*pole/(xi(1:npti)+pole)

!    define the grid on each side of the pole x1 before x2 after
     do ii=1,npt1
       x1(ii)=dx1*dble(ii-1)
     end do
     do ii=1,npt2
       x2(ii)=pole+del+dx2*dble(ii-1)
     end do

!    interpolate the initial fct fii on the new grids x1 and x2 (cubic spline)
!    write(std_out,*) npti,npt1

!    MJV 6/12/2008:
!    for each use of fctii should ensure that npt1 npt2 etc... are less than
!    npt=len(fctii)
     call intrpl(npti,xi,fctii,npt1,x1,fct4,fct1,fct5,1)
     call intrpl(npti,xi,fctii,npt2,x2,fct3,fct2,fct5,1)

!    calculate the two integrals from 0-->pole-lamda and pole+lamda--> end range
!    trapezoidal integration
     do ii=1,npt1
       fct1(ii)=fct4(ii)/(x1(ii)-pole)
     end do
     do ii=1,npt2
       fct2(ii)=fct3(ii)/(x2(ii)-pole)
     end do

     do ii=2,npt1-1
       xsum=xsum+fct1(ii)*dx1
     end do
     do ii=2,npt2-1
       xsum=xsum+fct2(ii)*dx2
     end do
     xsum=xsum+half*(fct1(1)+fct1(npt1))*dx1+half*(fct2(1)+fct2(npt2))*dx2

!    calculate the first and third derivative at the pole and add the taylor expansion
     call intrpl(npti,xi,fctii,npti,xi,fct3,fct4,fct5,1)
     call intrpl(npti,xi,fct4,1,(/pole/),fp,fpp,fppp,1)

     idel=two*fp(1)*(del)+fppp(1)*(del**3)/nine
     xsum=xsum+idel

   end if

!  calculate the derivated optical quantities and output the value 
   sigma2=(-two/pi)*xsum
   eps1=one-(four_pi*sigma2/(pole))
   eps2=four*fcti(ip)*pi/(pole)

!  A special treatment of the case where eps2 is very small compared to eps1 is needed
   if(eps2**2 > eps1**2 * tol12)then
     nomega(ip)=sqrt(half*(eps1 + sqrt(eps1**2 + eps2**2)))
     komega=sqrt(half*(-eps1 + sqrt(eps1**2 + eps2**2)))
     abso(ip)=four_pi*fcti(ip)*ohmtosec*Ohmcm/nomega(ip)/(Sp_Lt_SI*100._dp)
   else if(eps1>zero)then
     nomega(ip)=sqrt(half*(eps1 + sqrt(eps1**2 + eps2**2)))
     komega=half*abs(eps2/sqrt(eps1))
     abso(ip)=four_pi*fcti(ip)*ohmtosec*Ohmcm/nomega(ip)/(Sp_Lt_SI*100._dp)
   else if(eps1<zero)then
     nomega(ip)=half*abs(eps2/sqrt(-eps1))
     komega=sqrt(half*(-eps1 + sqrt(eps1**2 + eps2**2)))
     abso(ip)=two*sqrt(-eps1)*pole*ohmtosec*Ohmcm/(Sp_Lt_SI*100._dp)
   end if

   refl=((one-nomega(ip))**2 + komega**2)/ &
&   ((one+nomega(ip))**2 + komega**2)

   write(eps_unt,'(5e18.10)') Ha_eV*pole,fcti(ip)*Ohmcm,sigma2*Ohmcm,eps1,eps2
   write(abs_unt,'(5e18.10)') Ha_eV*pole,nomega(ip),komega,refl,abso(ip) 

 end do

 close(eps_unt)
 close(abs_unt) 

 ABI_DEALLOCATE(fct)
 ABI_DEALLOCATE(x1)
 ABI_DEALLOCATE(x2)
 ABI_DEALLOCATE(fct2)
 ABI_DEALLOCATE(fct3)
 ABI_DEALLOCATE(fct4)
 ABI_DEALLOCATE(fct5)
 ABI_DEALLOCATE(fp)
 ABI_DEALLOCATE(fpp)
 ABI_DEALLOCATE(fppp)
 ABI_DEALLOCATE(fct1)
 ABI_DEALLOCATE(ppsig)
 ABI_DEALLOCATE(fctii)
 ABI_DEALLOCATE(abso)
 ABI_DEALLOCATE(nomega)

end subroutine msig
!!***
