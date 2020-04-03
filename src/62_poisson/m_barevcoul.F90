!!****m* ABINIT/m_barevcoul
!! NAME
!!  m_barevcoul
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

module m_barevcoul

 use defs_basis
 use m_errors
 use m_fstrings,        only : sjoin, itoa
 use m_fft,             only : zerosym

 implicit none

 private
!!***

 public :: barevcoul
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/barevcoul
!! NAME
!! barevcoul
!!
!! FUNCTION
!! Compute bare coulomb term in G-space on the FFT mesh i.e. 4pi/(G+q)**2
!!
!! INPUTS
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  gsqcut=cutoff value on G**2 for sphere inside fft box. (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  icutcoul=Option for the Coulomb potential cutoff technique
!!  divgq0= value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq. Used if q = Gamma
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  izero=if 1, unbalanced components of V(q,g) are set to zero # Used by the PAW library
!!  nfft=Total number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  barev(nfft)=4pi/(G+q)**2, G=0 component is set to divgq0/pi if q = Gamma.
!!
!! NOTES
!!  This routine operates on the full FFT mesh. DO NOT PASS MPI_TYPE
!!  One can easily implemente MPI-FFT by just calling this routine and then
!!  extracting the G-vectors treated by the node.
!!
!! PARENTS
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine barevcoul(qphon,gsqcut,gmet,izero,nfft,nkpt_bz,ngfft,ucvol,barev)

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: izero,nfft,nkpt_bz,icutcoul
 real(dp),intent(in)  :: gsqcut,ucvol
!arrays
 integer,intent(in)   :: ngfft(18)
 real(dp),intent(in)  :: gmet(3,3),qphon(3)
 real(dp),intent(out) :: barev(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter    :: cplex1=1
 integer              :: i1,i2,i23,i3,id1,id2,id3
 integer              :: ig,ig1min,ig1,ig1max,ig2,ig2min,ig2max,ig3,ig3min,ig3max
 integer              :: ii,ii1,ing,n1,n2,n3,qeq0,qeq05
 real(dp),parameter   :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp)             :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3,rcut,divgq0
 character(len=100)   :: cutoff_method
!arrays
 integer :: id(3)
 real(dp),allocatable :: gq(:,:)

! Re-use variable defined initially in m_vcoul
 if (icutcoul==0) cutoff_method='SPHERE'
 if (icutcoul==1) cutoff_method='CYLINDER'
 if (icutcoul==2) cutoff_method='SURFACE'
 if (icutcoul==3) cutoff_method='CRYSTAL'
 if (icutcoul==4) cutoff_method='ERF'
 if (icutcoul==5) cutoff_method='ERFC'
 if (icutcoul==6) cutoff_method='AUXILIARY_FUNCTION'
 if (icutcoul==7) cutoff_method='AUX_GB' 

! Treatment of the divergence at q+g=zero

 cutoff_method="SPHERE"

 SELECT CASE (TRIM(cutoff_method))

 CASE ('SPHERE')
    ! Spence&Alavi scheme 
    rcut  = (three*nkpt_bz*ucvol/four_pi)**(one/three)
    divgq0= two_pi*rcut**two
   
 CASE ('CYLINDER')

 CASE ('SURFACE')
    ! Rozzi et al. scheme 

    ! Ismail-Beigi scheme   
    
 CASE ('AUX_FUNCTION')

 CASE ('AUX_GB')

 CASE ('CRYSTAL')
 
 CASE ('ERF')
 
 CASE ('ERFC')

 CASE DEFAULT
   MSG_BUG(sjoin('Unknown cutoff mode: ', cutoff_method))
   MSG_BUG('Please refer to our documentation related to icutcoul variable.')
 END SELECT

!Initialize a few quantities
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 cutoff=gsqcut*tolfix
 barev=zero

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12

     i23=n1*(i2-1 +(n2)*(i3-1))
     ! Do the test that eliminates the Gamma point outside of the inner loop
     ii1=1
     if (i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0) then
       ii1=2
       ! value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq
       barev(1+i23)=divgq0

     end if

     ! Final inner loop on the first dimension (note the lower limit)
     do i1=ii1,n1
       gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
       ii=i1+i23

       if(gs<=cutoff)then
         ! Identify min/max indexes (to cancel unbalanced contributions later)
         ! Count (q+g)-vectors with similar norm
         if ((qeq05==1).and.(izero==1)) then
           ig1=i1-(i1/id1)*n1-1
           ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
           ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
           ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
         end if

         den=piinv/gs

         barev(ii)=barev(ii)+den*(one-cos(rcut*sqrt(four_pi/den)))

       end if ! Cut-off
     end do ! End loop on i1
   end do ! End loop on i2
 end do ! End loop on i3

 ABI_DEALLOCATE(gq)

end subroutine barevcoul
!!***

end module m_barevcoul
!!***
