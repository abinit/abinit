!!****m* ABINIT/m_jellium
!! NAME
!!  m_jellium
!!
!! FUNCTION
!!  Routines related to jellium
!!
!! COPYRIGHT
!! Copyright (C) 2007-2019 ABINIT group (SC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_jellium

 use defs_basis
 use m_errors
 use m_abicore

 use defs_abitypes, only : MPI_type
 use m_fft,         only : fourdp

 implicit none

 private
!!***

 public :: jellium
!!***

contains
!!***

!!****f* m_jellium/jellium
!! NAME
!! jellium
!!
!! FUNCTION
!! Optionally compute
!!  option=1 : local ionic (external) potential due to jellium background
!!  option=2 : contribution to the initial density taking into account
!!                the jellium slab
!!
!! INPUTS
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2) (sphere for density and potential)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  option=
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  slabwsrad=Wigner-Seitz radius of jellium background
!!  slabzstart,slabzend=edges of jellium slab
!!
!! OUTPUT
!!  (if option==1) vjell(nfft)=external potential due to jellium background
!!  (if option==1) rhog(2,nfft), rhor(nfft,nspden)=density of positive charge
!!   in reciprocal, real space (only used in setvtr!)
!!
!! SIDE EFFECTS
!!  (if option==2) rhog(2,nfft), rhor(nfft,nspden)=reciprocal, real space
!!   updated initial electronic density
!!
!! PARENTS
!!      extraprho,gstate,setvtr
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

subroutine jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,nspden,&
&  option,slabwsrad,rhog,rhor,rprimd,vjell,slabzstart,slabzend)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden,option
 real(dp),intent(in) :: gsqcut,slabwsrad,slabzend,slabzstart
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3)
 real(dp),intent(inout) :: rhog(2,nfft),rhor(nfft,min(option,nspden))
 real(dp),intent(out) :: vjell(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,id1,id2,id3,ig1,ig2,ig3,ii,ispden,n1,n2,n3,nfftot
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: areaxy,cgz1,cgz2,cutoff,deltac,deltas,gsquar,gz,rhoave,sfi,sfr
 real(dp) :: sgz1,sgz2,zcellength
 character(len=500) :: message
!arrays
 real(dp),allocatable :: rhjg(:,:),rhjr(:),vjelg(:,:)

! *********************************************************************

!Enforce that nspden<=2
 if(nspden>2) then
   MSG_ERROR('Jellium possible only with nspden <= 2.')
 end if

!Make sure option is acceptable
 if (option/=1 .and. option/=2) then
   write(message, '(a,i0,3a)' )&
&   'option=',option,' is not allowed.',ch10,&
&   'Must be 1 or 2.'
   MSG_BUG(message)
 end if

 zcellength=rprimd(3,3)
 areaxy=abs(rprimd(1,1)*rprimd(2,2)-rprimd(1,2)*rprimd(2,1))
 rhoave=-half*three/(four_pi*slabwsrad**3)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2
 nfftot=n1*n2*n3
 cutoff=gsqcut*tolfix

 ABI_ALLOCATE(rhjg,(2,nfft))
 ABI_ALLOCATE(rhjr,(nfft))
 rhjg(:,:)=zero
 if(option==1) then
   ABI_ALLOCATE(vjelg,(2,nfft))
   vjelg(:,:)=zero
 end if

!Produce the potential due to the jellium background
 ii=0
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     do i1=1,n1

       ig1=i1-(i1/id1)*n1-1
       ii=ii+1
       gsquar=gsq_jel(ig1,ig2,ig3)

!      Skip G**2 outside cutoff and use \delta_{G_\|,0}:
       if (gsquar<=cutoff.and.ig1==0.and.ig2==0) then

!        N o t e   t h a t   gz=two_pi*sqrt(gsq(0,0,ig3))
         gz=dble(ig3)*two_pi/zcellength

!        G_z == 0
         if (ig3==0) then
           sfr=two*rhoave*(slabzend-slabzstart)
           sfi=zero
!          G_z /= 0
         else ! of ig3==0
           sgz2=sin(gz*slabzend) ; sgz1=sin(gz*slabzstart)
           cgz2=cos(gz*slabzend) ; cgz1=cos(gz*slabzstart)
           deltas=sgz2-sgz1
           deltac=cgz2-cgz1
           sfr=two*rhoave*deltas/gz
           sfi=two*rhoave*deltac/gz
           if(option==1) then
!            Assemble vjell_G
             vjelg(re,ii)=four_pi*sfr/gz**2
             vjelg(im,ii)=four_pi*sfi/gz**2
           end if
         end if ! of ig3==0
!        Assemble \rho_G
         rhjg(re,ii)=sfr
         rhjg(im,ii)=sfi

       end if ! of gsquar ...

!      End loop on i1
     end do
!    End loop on i2
   end do
!  End loop on i3
 end do

 rhjg(:,:)=rhjg(:,:)/zcellength
 if(option==1) vjelg(:,:)=vjelg(:,:)/zcellength

 call fourdp(1,rhjg,rhjr,1,mpi_enreg,nfft,1,ngfft,0)
 if(option==1) then
   call fourdp(1,vjelg,vjell,1,mpi_enreg,nfft,1,ngfft,0)
   rhog(:,:)=rhjg(:,:)
   rhor(:,1)=rhjr(:)
 else
!  Update the initial electronic density adding -rhjr
   rhog(:,:)=rhog(:,:)-rhjg(:,:)
   do ispden=1,nspden
     rhor(:,ispden)=rhor(:,ispden)-rhjr(:)/dble(ispden)
   end do
 end if

 ABI_DEALLOCATE(rhjg)
 ABI_DEALLOCATE(rhjr)
 if(option==1) then
   ABI_DEALLOCATE(vjelg)
 end if

!DEBUG
!write(std_out,*)' jellium : exit '
!stop
!ENDDEBUG

 contains

   function gsq_jel(i1,i2,i3)

   real(dp) :: gsq_jel
   integer,intent(in) :: i1,i2,i3
   gsq_jel=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_jel

end subroutine jellium
!!***

end module m_jellium
!!***
