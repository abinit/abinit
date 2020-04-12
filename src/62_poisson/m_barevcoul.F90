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
 use m_profiling_abi,   only : abimem_record

 use m_gsphere,         only : gsphere_t
 use m_qplusg,          only : cmod_qpg
 use m_cutoff_sphere

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

subroutine barevcoul(rcut,shortrange,qphon,gsqcut,gmet,nfft,nkpt_bz,ngfft,ucvol,barev)

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: nfft,nkpt_bz
 real(dp),intent(in)  :: rcut,gsqcut,ucvol
!arrays
 integer,intent(in)      :: ngfft(18)
 real(dp),intent(in)     :: qphon(3)
 real(dp),intent(inout)  :: gmet(3,3)
 real(dp),intent(inout)  :: barev(nfft)
!Local variables-------------------------------
!scalars
 integer,parameter    :: icutcoul=0,empty(3,3)=0
 integer              :: i1,i2,i23,i3,id1,id2,id3
 integer              :: ig,ig1min,ig1max,ig2min,ig2max,ig3min,ig3max
 integer              :: ii,ing,n1,n2,n3
 real(dp),parameter   :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp)             :: cutoff,gqg2p3,gqgm12,gqgm13,gqgm23,gs2,gs3,divgq0,rcut0,testv(nfft)
 logical              :: shortrange
!arrays
 integer :: id(3)
 real(dp),allocatable :: gq(:,:),gpq(:),gpq2(:)

! Treatment of the divergence at q+g=zero
 rcut0= (three*nkpt_bz*ucvol/four_pi)**(one/three)
 divgq0= two_pi*rcut0**two

!Initialize a few quantities
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 cutoff=gsqcut*tolfix
 barev=zero

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 ABI_ALLOCATE(gpq,(nfft))
 ABI_ALLOCATE(gpq2,(nfft))
 
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

 do i3=1,n3
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2
   do i2=1,n2
     i23=n1*(i2-1 +(n2)*(i3-1))
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12
     do i1=1,n1
        ii=i1+i23
        gpq(ii)=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
        gpq2(ii) = piinv/four/gpq(ii)
     end do
   end do
 end do

 do ig=1,nfft 
     if(abs(gpq(ig))<tol8) then 
        barev(ig)=barev(ig)+divgq0
     else if(gpq(ig)<=cutoff) then
       if(shortrange) then
         barev(ig)=barev(ig)+gpq2(ig)*(one-exp(-pi/(gpq2(ig)*rcut**2)))
       else
         barev(ig)=barev(ig)+gpq(ig)*(one-cos(rcut0*sqrt(four_pi/gpq(ig))))
       end if
    end if
 end do

! if(shortrange) then
!    continue
! else
!    call cutoff_sphere(nfft,gq,1,empty,gmet,rcut0,testv)
! end if

 ABI_DEALLOCATE(gq)
 ABI_DEALLOCATE(gpq)
 ABI_DEALLOCATE(gpq2)

end subroutine barevcoul
!!***

end module m_barevcoul
!!***
