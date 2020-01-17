!!****m* ABINIT/m_sgfft
!! NAME
!!  m_sgfft
!!
!! FUNCTION
!!  This module provides low-level interfaces to Goedecker's FFT library.
!!
!! COPYRIGHT
!! Copyright by Stefan Goedecker, Ithaca, NY USA, July 14, 1993
!! Copyright (C) 1998-2019 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_sgfft

 use defs_basis
 use m_abicore
 use m_errors
 use m_fftcore

 use defs_fftdata,  only : mg

 implicit none

 private

! Public API.
 public :: sg_fft_cc      ! Complex-Complex version (full box)
 public :: sg_fft_rc      ! Real-Complex version (full box)
 public :: sg_fftpad      ! Zero-padding version of "fft".
 public :: sg_fftrisc     ! Fourier transforms of wavefunctions
 public :: sg_fftrisc_2
 public :: sg_poisson     ! Solve the poisson equation in G-space starting from n(r).

CONTAINS  !====================================================================
!!***


!!****f* m_sgfft/sg_fft_cc
!! NAME
!! sg_fft_cc
!!
!! FUNCTION
!! Calculates the discrete Fourier transform:
!!
!!   ftarr(i1,i2,i3)=exp(ris*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) arr(j1,j2,j3)
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  n1,n2,n3=physical dimension of the transform
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  ndat=Number of FFT transforms
!!  isign=+1 for G-->R, -1 for R-->G
!!  arr(2,nd1*nd2*nd3*ndat)=input complex array with alternating real and imaginary
!!  elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  (see SIDE FFECTS).
!!
!! OUTPUT
!!  ftarr(2,nd1*nd2*nd3*ndat)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!  arr(2,nd1*nd2*nd3*ndat) is modified by sg_fftx,sg_ffty,sg_fftz.
!!
!! NOTES
!!  ndi must always be greater or equal to ni.  Recommended choice for nd1
!!  and nd2 is: ni for ni=odd or ni+1 for ni=even (hence 2*(ni/2)+1);
!!  nd3 should always be n3.  Note that choosing nd1 or nd2 larger than
!!  the recommended value can severely degrade efficiency of this routine.
!!  Avoiding even ndi for nd1 and nd2 avoids cache conflicts on cache machines.
!!  Each of n1,n2,n3 must be a
!!  product of the prime factors 2,3,5. If two ni s are equal
!!  it is recommended to place them behind each other.
!!  The largest any of these may be is set by parameter "mg" below.
!!  This fft is particularly efficient for cache architectures.
!!  Note that the meaning of fftcache has changed from the original
!!  ncache of SG (that was the maximum number of COMPLEX*16 in the cache)
!!
!! PARENTS
!!      ccfft,m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fft_cc(fftcache,n1,n2,n3,nd1,nd2,nd3,ndat,isign,arr,ftarr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,n1,n2,n3,nd1,nd2,nd3,ndat,isign
!arrays
 real(dp),intent(inout) :: arr(2,nd1*nd2*nd3*ndat)
 real(dp),intent(inout) :: ftarr(2,nd1*nd2*nd3*ndat)

!Local variables-------------------------------
!scalars
 integer :: idat,start

! *************************************************************************

 do idat=1,ndat
   start = 1 + (idat-1)*nd1*nd2*nd3
   call fft_cc_one_nothreadsafe(fftcache,nd1,nd2,nd3,n1,n2,n3,arr(1,start),ftarr(1,start),real(isign,kind=dp))
 end do

end subroutine sg_fft_cc
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/fft_cc_one_nothreadsafe
!! NAME
!! fft_cc_one_nothreadsafe
!!
!! FUNCTION
!! Calculates the discrete Fourier transform:
!!
!!   ftarr(i1,i2,i3)=exp(ris*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) arr(j1,j2,j3)
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  n1,n2,n3=physical dimension of the transform
!!  arr(2,nd1,nd2,nd3)=input complex array with alternating real and imaginary
!!  elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  (see SIDE FFECTS).
!!  ris=(real(dp)) sign of exponential in transform
!!
!! OUTPUT
!!  ftarr(2,nd1,nd2,nd3)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!  arr(2,nd1,nd2,nd3) is modified by sg_fftx,sg_ffty,sg_fftz.
!!
!! NOTES
!!  ndi must always be greater or equal to ni.  Recommended choice for nd1
!!  and nd2 is: ni for ni=odd or ni+1 for ni=even (hence 2*(ni/2)+1);
!!  nd3 should always be n3.  Note that choosing nd1 or nd2 larger than
!!  the recommended value can severely degrade efficiency of this routine.
!!  Avoiding even ndi for nd1 and nd2 avoids cache conflicts on cache machines.
!!  Each of n1,n2,n3 must be a
!!  product of the prime factors 2,3,5. If two ni s are equal
!!  it is recommended to place them behind each other.
!!  The largest any of these may be is set by parameter "mg" below.
!!  This fft is particularly efficient for cache architectures.
!!  Note that the meaning of fftcache has changed from the original
!!  ncache of SG (that was the maximum number of COMPLEX*16 in the cache)
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine fft_cc_one_nothreadsafe(fftcache,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,n1,n2,n3,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 real(dp),intent(inout) :: arr(2,nd1,nd2,nd3)
 real(dp),intent(inout) :: ftarr(2,nd1,nd2,nd3) !vz_i

!Local variables-------------------------------
!mfac sets maximum number of factors (5, 4, 3, or 2) which may be
!contained within any n1, n2, or n3
!mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
!scalars
 integer,parameter :: mfac=11
 integer :: i2,ic,n1i,n3i
 character(len=500) :: message
!arrays
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: trig(2,mg)

! *************************************************************************

!Check that dimension is not exceeded
 if (n1>mg.or.n2>mg.or.n3>mg) then
   write(message, '(a,3i10,a,i10,a)' )&
&   'one of the dimensions n1,n2,n3=',n1,n2,n3,&
&   'exceeds allowed dimension mg=',mg,ch10
   MSG_BUG(message)
 end if

!transform along x direction
 call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 call sg_fftx(fftcache,mfac,mg,nd1,nd2,nd3,n2,n3,&
& arr,ftarr,trig,aft,now,bef,ris,ind,ic)

!transform along y direction
 if (n2/=n1)then
   call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 end if
 n1i=1 ; n3i=1
 call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,&
& ftarr,arr,trig,aft,now,bef,ris,ind,ic)

!transform along z direction
 if (n3/=n2)then
   call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 end if

!$OMP PARALLEL DO SHARED(aft,arr,bef,ftarr,ind,ic)&
!$OMP SHARED(nd1,nd2,nd3,now,n1,n2,ris,trig)&
!$OMP PRIVATE(i2)
 do i2=1,n2
   call sg_fftz(mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&   trig,aft,now,bef,ris,ind,ic)
 end do
!$OMP END PARALLEL DO

end subroutine fft_cc_one_nothreadsafe
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fft_rc
!! NAME
!! sg_fft_rc
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2. Usually used for density and potentials.
!!
!! There are two different possibilities :
!!  fftalgb=0 means using the complex-to-complex FFT routine,
!!   irrespective of the value of cplex
!!  fftalgb=1 means using a real-to-complex FFT or a complex-to-complex FFT,
!!   depending on the value of cplex.
!!  The only real-to-complex FFT available is from SGoedecker library.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! isign=sign of Fourier transform exponent: current convention uses
!!  +1 for transforming from G to r
!!  -1 for transforming from r to G.
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! PARENTS
!!      fourdp
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fft_rc(cplex,fofg,fofr,isign,nfft,ngfft)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: fofg(2,nfft),fofr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer :: fftalg,fftalga,fftalgb,fftcache,i1,i2,i3,ic1,ic2,ic3,index
 integer :: n1,n1half1,n1halfm,n2,n2half1,n3,n4,n4half1,n5,n5half1,n6
 real(dp) :: ris,xnorm
 character(len=500) :: msg
!arrays
 integer :: aft1(mfac),aft2(mfac),aft3(mfac),bef1(mfac),bef2(mfac),bef3(mfac)
 integer :: ind1(mg),ind2(mg),ind3(mg),now1(mfac),now2(mfac),now3(mfac)
 real(dp) :: trig1(2,mg),trig2(2,mg),trig3(3,mg)
 real(dp),allocatable :: wk2d_a(:,:,:,:),wk2d_b(:,:,:,:),wk2d_c(:,:,:,:)
 real(dp),allocatable :: wk2d_d(:,:,:,:),work1(:,:,:,:),work2(:,:,:,:)

! *************************************************************************

 !DBG_ENTER("COLL")

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 fftcache=ngfft(8)
 fftalg  =ngfft(7)
 fftalga =fftalg/100
 fftalgb =mod(fftalg,100)/10

 ris=dble(isign)
 xnorm=1.0d0/dble(n1*n2*n3)

 if (fftalgb/=0 .and. fftalgb/=1) then
   write(msg, '(a,i4,a,a,a,a,a)' )&
&   'The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   'The second digit (fftalg(B)) must be 0 or 1.',ch10,&
&   'Action: change fftalg in your input file.'
   MSG_BUG(msg)
 end if

 if (fftalgb==1 .and. ALL(fftalga/=(/1,3,4/)) )then
   write(msg,'(a,i4,5a)')&
&   'The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   'When fftalg(B) is 1, the allowed values for fftalg(A) are 1 and 4.',ch10,&
&   'Action: change fftalg in your input file.'
   MSG_BUG(msg)
 end if

 if (n4<n1.or.n5<n2.or.n6<n3) then
   write(msg,'(a,3i8,a,3i8)')'  Each of n4,n5,n6=',n4,n5,n6,'must be >= n1, n2, n3 =',n1,n2,n3
   MSG_BUG(msg)
 end if

!---------------------------------------------------------
!Here sophisticated algorithm based on S. Goedecker routines, only for the REAL case.
!Take advantage of the fact that fofr is real, and that fofg has corresponding symmetry properties.

#ifdef DEBUG_MODE
 if (n1>mg .or. n2>mg .or. n3>mg) then
   write(msg, '(a,3i10,a,a,a,i10,a)' )&
&   'One of the dimensions n1,n2,n3=',n1,n2,n3,',',ch10,&
&   'exceeds allowed dimension mg=',mg,'.'
   MSG_BUG(msg)
 end if
#endif

 n1half1=n1/2+1 ; n1halfm=(n1+1)/2
 n2half1=n2/2+1
!n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
 n4half1=(n1half1/2)*2+1
 n5half1=(n2half1/2)*2+1

!This sophisticated algorithm allows to decrease the memory needs.
 ABI_ALLOCATE(work1,(2,n4,n5half1,n6))
 ABI_ALLOCATE(work2,(2,n4,n5half1,n6))

 if(isign==1)then

!  Compute auxiliary arrays needed for FFTs, here forward FFT
   call sg_ctrig(n1,trig1,aft1,bef1,now1,one,ic1,ind1,mfac,mg)
   call sg_ctrig(n2,trig2,aft2,bef2,now2,one,ic2,ind2,mfac,mg)
   call sg_ctrig(n3,trig3,aft3,bef3,now3,one,ic3,ind3,mfac,mg)

!  Transfer fofg to the expanded fft box (only half of it)

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work1)
   do i3=1,n3
     do i2=1,n2half1
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
         work1(1,i1,i2,i3)=fofg(1,i1+index)
         work1(2,i1,i2,i3)=fofg(2,i1+index)
       end do
     end do
   end do

!$OMP PARALLEL DO SHARED(aft3,bef3,ind3,ic3,now3,n1,n2half1,n4,n5half1,n6,ris,trig3,work1,work2) PRIVATE(i2)
   do i2=1,n2half1
     call sg_fftz(mfac,mg,n4,n5half1,n6,n1,i2,i2,work1,work2,&
&     trig3,aft3,now3,bef3,ris,ind3,ic3)
   end do

!  Loop over x-y planes

!$OMP PARALLEL PRIVATE(i1,i2,i3,index,wk2d_a,wk2d_b,wk2d_c,wk2d_d) &
!$OMP&SHARED(aft1,aft2,bef1,bef2,fftcache,fofg,fofr,ic1,ic2,ind1,ind2) &
!$OMP&SHARED(n1,n1half1,n1halfm,n2,n2half1,n3) &
!$OMP&SHARED(n4,n5,now1,now2,ris,trig1,trig2,work2)

   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_c,(2,2*n1halfm+1,n5,1))
   ABI_ALLOCATE(wk2d_d,(2,2*n1halfm+1,n5,1))

!$OMP DO
   do i3=1,n3

     do i2=1,n2half1
       do i1=1,n1
         wk2d_c(1,i1,i2,1)=work2(1,i1,i2,i3)
         wk2d_c(2,i1,i2,1)=work2(2,i1,i2,i3)
       end do
     end do

     call sg_fftx(fftcache,mfac,mg,2*n1halfm+1,n5,1,n2half1,1,wk2d_c,wk2d_d,&
&     trig1,aft1,now1,bef1,ris,ind1,ic1)

     do i1=1,n1half1-1 ! Compute symmetric and antisymmetric combinations
       wk2d_a(1,i1,1,1)=wk2d_d(1,2*i1-1,1,1)
       wk2d_a(2,i1,1,1)=wk2d_d(1,2*i1  ,1,1)
     end do

     if((2*n1half1-2)/=n1)then  ! If n1 odd, must add last data
       wk2d_a(1,n1half1,1,1)=wk2d_d(1,n1,1,1)
       wk2d_a(2,n1half1,1,1)=0.0d0
     end if

     do i2=2,n2half1
       do i1=1,n1half1-1
         wk2d_a(1,i1,i2,1)     = wk2d_d(1,2*i1-1,i2,1)-wk2d_d(2,2*i1,i2,1)
         wk2d_a(2,i1,i2,1)     = wk2d_d(2,2*i1-1,i2,1)+wk2d_d(1,2*i1,i2,1)
         wk2d_a(1,i1,n2+2-i2,1)= wk2d_d(1,2*i1-1,i2,1)+wk2d_d(2,2*i1,i2,1)
         wk2d_a(2,i1,n2+2-i2,1)=-wk2d_d(2,2*i1-1,i2,1)+wk2d_d(1,2*i1,i2,1)
       end do
       if((2*n1half1-2)/=n1)then
         wk2d_a(1,n1half1,i2,1)     = wk2d_d(1,n1,i2,1)
         wk2d_a(2,n1half1,i2,1)     = wk2d_d(2,n1,i2,1)
         wk2d_a(1,n1half1,n2+2-i2,1)= wk2d_d(1,n1,i2,1)
         wk2d_a(2,n1half1,n2+2-i2,1)=-wk2d_d(2,n1,i2,1)
       end if
     end do

     call sg_ffty(fftcache,mfac,mg,n4,n5,1,1,n1halfm,1,1,wk2d_a,wk2d_b,&
&     trig2,aft2,now2,bef2,ris,ind2,ic2)

     do i2=1,n2  ! Take real part data from expanded box and put it in the original box.
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1half1-1 ! copy data
         fofr(2*i1-1+index)=wk2d_b(1,i1,i2,1)
         fofr(2*i1  +index)=wk2d_b(2,i1,i2,1)
       end do
       if((2*n1half1-2)/=n1)then ! If n1 odd, must add last data
         fofr(n1+index)=wk2d_b(1,n1half1,i2,1)
       end if
     end do

   end do ! loop over x-y planes
!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_c)
   ABI_DEALLOCATE(wk2d_d)
!$OMP END PARALLEL

 else if(isign==-1)then

!  Compute auxiliary arrays needed for FFTs, here backward FFT
   call sg_ctrig(n1,trig1,aft1,bef1,now1,-one,ic1,ind1,mfac,mg)
   call sg_ctrig(n2,trig2,aft2,bef2,now2,-one,ic2,ind2,mfac,mg)
   call sg_ctrig(n3,trig3,aft3,bef3,now3,-one,ic3,ind3,mfac,mg)

!  Treat first x-transform in x-y plane, and multiply
!  by overall normalization factor 1/nfftot

!  Loop over x-y planes

!$OMP PARALLEL PRIVATE(i1,i2,i3,index,wk2d_a,wk2d_b,wk2d_c,wk2d_d) &
!$OMP&SHARED(aft1,aft2,bef1,bef2,fftcache,fofr,ic1,ic2,ind1,ind2) &
!$OMP&SHARED(n1,n1half1,n1halfm,n2,n2half1,n3) &
!$OMP&SHARED(n4,n5,now1,now2,ris,trig1,trig2,work1,xnorm)

   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_c,(2,2*n1halfm+1,n5,1))
   ABI_ALLOCATE(wk2d_d,(2,2*n1halfm+1,n5,1))

!$OMP DO
   do i3=1,n3
     do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1half1-1 ! copy and normalize data
         wk2d_a(1,i1,i2,1)=fofr(2*i1-1+index)*xnorm
         wk2d_a(2,i1,i2,1)=fofr(2*i1  +index)*xnorm
       end do

       if((2*n1half1-2)/=n1)then ! If n1 odd, must add last data
         wk2d_a(1,n1half1,i2,1)=fofr(n1+index)*xnorm
         wk2d_a(2,n1half1,i2,1)=zero
       end if
     end do

     call sg_ffty(fftcache,mfac,mg,n4,n5,1,1,n1halfm,1,1,wk2d_a,wk2d_b,&
&     trig2,aft2,now2,bef2,ris,ind2,ic2)

     do i1=1,n1halfm ! Decompose symmetric and antisymmetric parts
       wk2d_c(1,2*i1-1,1,1)=wk2d_b(1,i1,1,1)
       wk2d_c(2,2*i1-1,1,1)=0.0d0
       wk2d_c(1,2*i1,1,1)=wk2d_b(2,i1,1,1)
       wk2d_c(2,2*i1,1,1)=0.0d0
     end do

     do i2=2,n2half1
       do i1=1,n1halfm
         wk2d_c(1,2*i1-1,i2,1)= (wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+2-i2,1))*0.5d0
         wk2d_c(2,2*i1-1,i2,1)= (wk2d_b(2,i1,i2,1)-wk2d_b(2,i1,n2+2-i2,1))*0.5d0
         wk2d_c(1,2*i1,i2,1)  = (wk2d_b(2,i1,i2,1)+wk2d_b(2,i1,n2+2-i2,1))*0.5d0
         wk2d_c(2,2*i1,i2,1)  =-(wk2d_b(1,i1,i2,1)-wk2d_b(1,i1,n2+2-i2,1))*0.5d0
       end do
     end do

     call sg_fftx(fftcache,mfac,mg,2*n1halfm+1,n5,1,n2half1,1,wk2d_c,wk2d_d,&
&     trig1,aft1,now1,bef1,ris,ind1,ic1)

     do i2=1,n2half1
       do i1=1,n1
         work1(1,i1,i2,i3)=wk2d_d(1,i1,i2,1)
         work1(2,i1,i2,i3)=wk2d_d(2,i1,i2,1)
       end do
     end do

   end do
!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_c)
   ABI_DEALLOCATE(wk2d_d)
!$OMP END PARALLEL

!$OMP PARALLEL DO SHARED(aft3,bef3,ind3,ic3,now3,n1,n2half1,n4,n5half1,n6,ris,trig3,work1,work2) PRIVATE(i2)
   do i2=1,n2half1
     call sg_fftz(mfac,mg,n4,n5half1,n6,n1,i2,i2,work1,work2,&
&     trig3,aft3,now3,bef3,ris,ind3,ic3)
   end do

!  Transfer fft output to the original fft box

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n2half1,n3,work2)
   do i3=1,n3
     do i2=1,n2half1
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
         fofg(1,i1+index)=work2(1,i1,i2,i3)
         fofg(2,i1+index)=work2(2,i1,i2,i3)
       end do
     end do
!    Complete missing values with complex conjugate
!    Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
     if(n2half1>2)then
       do i2=2,n2+1-n2half1
         index=n1*((n2+2-i2)-1)
         if(i3/=1)index=index+n1*n2*((n3+2-i3)-1)
         fofg(1,1+index)= work2(1,1,i2,i3)
         fofg(2,1+index)=-work2(2,1,i2,i3)
         do i1=2,n1
           fofg(1,n1+2-i1+index)= work2(1,i1,i2,i3)
           fofg(2,n1+2-i1+index)=-work2(2,i1,i2,i3)
         end do
       end do
     end if
   end do

 end if ! choice of isign

 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)

 !DBG_EXIT("COLL")

end subroutine sg_fft_rc
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fftpad
!! NAME
!! sg_fftpad
!!
!! FUNCTION
!! Fast Fourier transform. This is the zero-padding version of "fft".
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  mgfft=maximum size of 1D FFTs
!!  n1,n2,n3=physical dimension of the transform
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  ndat=Number of FFT transforms.
!!  isign= sign of exponential in transform
!!  gbound(2*mgfft+8,2)=sphere boundary info
!!
!! OUTPUT
!!  ftarr(2,nd1,nd2,nd3*ndat)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!  arr(2,nd1,nd2,nd3*ndat)=input complex array with alternating real and imaginary
!!    elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  arr(2,nd1,nd2,nd3*ndat) is modified by sg_fftpx,sg_ffty,sg_fftz.
!!
!! NOTES
!!  mfac sets maximum number of factors (5, 4, 3, or 2) which may be
!!  contained within any n1, n2, or n3
!!  mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
!!  XG: the signification of mg is changed with respect to fft3dp !!!
!!
!! PARENTS
!!      fourwf,m_fft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fftpad(fftcache,mgfft,n1,n2,n3,nd1,nd2,nd3,ndat,gbound,isign,arr,ftarr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,mgfft,n1,n2,n3,nd1,nd2,nd3,ndat,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),intent(inout) :: arr(2,nd1,nd2,nd3*ndat)
 real(dp),intent(out) :: ftarr(2,nd1,nd2,nd3*ndat)

!Local variables-------------------------------
!scalars
 integer :: idat,start

! *************************************************************************

 do idat=1,ndat
   start = 1 + (idat-1)*nd3
   call fftpad_one_nothreadsafe(fftcache,mgfft,nd1,nd2,nd3,n1,n2,n3,&
&    arr(1,1,1,start),ftarr(1,1,1,start),real(isign, kind=dp),gbound)
 end do

end subroutine sg_fftpad
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/fftpad_one_nothreadsafe
!! NAME
!! fftpad_one_nothreadsafe
!!
!! FUNCTION
!! Fast Fourier transform. This is the zero-padding version of "fft" for a single array.
!! This version is not thread-safe.
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  mgfft=maximum size of 1D FFTs
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  n1,n2,n3=physical dimension of the transform
!!  arr(2,nd1,nd2,nd3)=input complex array with alternating real and imaginary
!!    elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  ris=(real(dp)) sign of exponential in transform
!!  gbound(2*mgfft+8,2)=sphere boundary info
!!
!! OUTPUT
!!  ftarr(2,nd1,nd2,nd3)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!  arr(2,nd1,nd2,nd3) is modified by sg_fftpx,sg_ffty,sg_fftz.
!!
!! NOTES
!!  mfac sets maximum number of factors (5, 4, 3, or 2) which may be
!!  contained within any n1, n2, or n3
!!  mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
!!  XG: the signification of mg is changed with respect to fft3dp !!!
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine fftpad_one_nothreadsafe(fftcache,mgfft,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris,gbound)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,mgfft,n1,n2,n3,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),intent(inout) :: arr(2,nd1,nd2,nd3)
 real(dp),intent(out) :: ftarr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer :: g3max,g3min,i2,ic,n1i,n3i,n3p
#ifdef DEBUG_MODE
 character(len=500) :: message
#endif
!arrays
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: trig(2,mg)

! *************************************************************************

#ifdef DEBUG_MODE
!Check that dimension is not exceeded
 if (n1>mg.or.n2>mg.or.n3>mg) then
   write(message, '(a,3i10,a,i10)')&
&   'one of the dimensions n1,n2,n3=',n1,n2,n3,' exceeds the allowed dimension mg=',mg
   MSG_BUG(message)
 end if
#endif

 g3min=gbound(3,2)
 g3max=gbound(4,2)

!--------------------------------------------------------------------------

 if (abs(ris-one)<tol12) then

!  Handle G -> r  transform (G sphere to fft box)

!  Transform along x direction
   call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)

!  Zero out the untransformed (0) data part of the work array
!  -- at every (y,z) there are 0 s to be added to the ends of
!  the x data so have to zero whole thing.
   ftarr(:,:,:,:)=0.0d0

!  Note the passing of the relevant part of gbound
   call sg_fftpx(fftcache,mfac,mg,mgfft,nd1,nd2,nd3,n2,n3,&
&   arr,ftarr,trig,aft,now,bef,ris,ind,ic,gbound(3,2))

!  Transform along y direction in two regions of z
   if (n2/=n1)then
     call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  First y transform: z=1..g3max+1
   n3p=g3max+1
   n1i=1 ; n3i=1
   call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3p,ftarr,arr,&
&   trig,aft,now,bef,ris,ind,ic)

!  Zero out the untransformed (0) data part of the work array
!  -- only need to zero specified ranges of z
   arr(:,:,:,n3p+1:g3min+n3)=0.0d0

!  Second y transform: z=g3min+1..0 (wrapped around)
   n3p=-g3min
   if (n3p>0) then
     n3i=1+g3min+n3 ; n1i=1
     call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,ftarr,arr,&
&     trig,aft,now,bef,ris,ind,ic)
   end if

!  Transform along z direction
   if (n3/=n2) then
     call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!$OMP PARALLEL DO
   do i2=1,n2
     call sg_fftz(mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&     trig,aft,now,bef,ris,ind,ic)
   end do

 else

!  *************************************************
!  Handle r -> G transform (from fft box to G sphere)

!  Transform along z direction
   call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)

!$OMP PARALLEL DO
   do i2=1,n2
     call sg_fftz(mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&     trig,aft,now,bef,ris,ind,ic)
   end do

!  Transform along y direction in two regions of z
   if (n2/=n3) then
     call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  First y transform: z=1..g3max+1
   n3p=g3max+1
   n1i=1 ; n3i=1
   call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3p,ftarr,arr,&
&   trig,aft,now,bef,ris,ind,ic)

!  Second y transform: z=g3min+1..0 (wrapped around)
   n3p=-g3min
   if (n3p>0) then
     n1i=1 ; n3i=1+g3min+n3
     call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,ftarr,arr,&
&     trig,aft,now,bef,ris,ind,ic)
   end if

!  Transform along x direction
   if (n1/=n2) then
     call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  Zero out the untransformed (0) data part of the work array
!  -- at every (y,z) there are 0 s to be added to the ends of
!  the x data so have to zero whole thing.
   ftarr(:,:,:,:)=0.0d0

!  Note the passing of the relevant part of gbound
   call sg_fftpx(fftcache,mfac,mg,mgfft,nd1,nd2,nd3,n2,n3,&
&   arr,ftarr,trig,aft,now,bef,ris,ind,ic,gbound(3,2))

!  Data is now ready to be extracted from fft box to sphere
 end if

end subroutine fftpad_one_nothreadsafe
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fftpx
!! NAME
!! sg_fftpx
!!
!! FUNCTION
!! This subroutine is called by the 3-dimensional fft to conduct the
!! "x" transforms for all y and z.
!! Accomodate more optimal treatment of
!! zero padding following the method of fft3dp.
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  mfac = maximum number of factors in 1D FFTs
!!  mg = maximum length of 1D FFTs
!!  mgfft = effective maximum length of 1D FFTs, for dimensioning gbound
!!  nd1=first dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd2=second dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd3=third dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  n2,n3=actual length of y and z transforms
!!  z(2,nd1,nd2,nd3)=INPUT array; destroyed by transformation
!!  trig, aft, now, bef, ind=provided by previous call to ctrig
!!   Note that in this routine (and in ctrig) the values in array trig are
!!   actually cos and tan, not cos and sin.  Use of tan allows advantageous
!!   use of FMA on the ibm rs6000.
!!  ris=sign of exponential in transform (should be 1 or -1; real)
!!  ic=number of (radix) factors of x transform length (from ctrig)
!!  gbound(2*mgfft+4)=sphere boundary info
!!
!! OUTPUT
!!  zbr(2,nd1,nd2,nd3)=OUTPUT transformed array; no scaling applied
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine blocks the x transforms
!! so that all transforms under consideration at one step fit within
!! the cache memory, which is crucial for optimal performance.
!! The blocking factor is set by parameter "fftcache" below, which should
!! be adjusted to be somewhat smaller (say 3/4) than the actual cache size
!! of the machine.
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fftpx(fftcache,mfac,mg,mgfft,nd1,nd2,nd3,n2,n3,&
&    z,zbr,trig,aft,now,bef,ris,ind,ic,gbound)

 implicit none

!Arguments ------------------------------------
!Dimensions of aft, now, bef, ind, and trig should agree with
!those in subroutine ctrig.
!scalars
 integer,intent(in) :: fftcache,ic,mfac,mg,mgfft,n2,n3,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 integer,intent(in) :: aft(mfac),bef(mfac),gbound(2*mgfft+4),ind(mg),now(mfac)
 real(dp),intent(in) :: trig(2,mg)
 real(dp),intent(inout) :: z(2,nd1,nd2,nd3)
 real(dp),intent(inout) :: zbr(2,nd1,nd2,nd3) !vz_i

!Local variables-------------------------------
!scalars
 integer :: g2,g2max,g2min,g3,g3max,g3min,gg3,i,ia,ib,igb,ihalfy,indx,j
 integer :: len3,lot,lowlim,ma,mb,ntb,upplim
!no_abirules
 real(dp),parameter :: &
& cos2=0.3090169943749474d0,&   !cos(2.d0*pi/5.d0)
& cos4=-0.8090169943749474d0,&  !cos(4.d0*pi/5.d0)
& sin42=0.6180339887498948d0    !sin(4.d0*pi/5.d0)/sin(2.d0*pi/5.d0)
 real(dp) :: bb,cr2,cr2s,cr3,cr3p,cr4,cr5,ct2,ct3,ct4,ct5,&
& factor,r,r1,r2,r25,r3,r34,r4,r5,s,sin2,s1,s2,s25,s3,s34,s4,s5

! *************************************************************************

 g3min=gbound(1)
 g3max=gbound(2)
 igb=3
 len3=g3max-g3min+1


!Do x transforms in blocks of size "lot" which is set by how
!many x transform arrays (of size nd1 each) fit into the nominal
!cache size "fftcache".
!Loop over blocks in the loop below.

 factor=0.75d0
 lot=(fftcache*factor*1000d0)/(nd1*8*2)
 if(lot.lt.1) lot=1
!Express loop over y, z in terms of separate z and y loops

!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP SHARED(aft,bef,gbound,g3max,ic,ind,len3,lot)&
!$OMP SHARED(n2,n3,nd2,now,ris,trig,z,zbr)
 do gg3=1,len3

   if (gg3<=g3max+1) then
     g3=gg3
   else
!    wrap around for negative gg3
     g3=gg3-len3+n3
   end if

   igb=gg3*2+1
   g2min=gbound(igb)
   g2max=gbound(igb+1)

!  Split the y loop into positive and wrapped-around negative parts

   do ihalfy=1,2

!    Start at 1 for ihalfy=1; g2min+1+n2 for ihalfy=2
     lowlim=1+(ihalfy-1)*(g2min+n2)
!    End at g2max+1 for ihalfy=1; n2 for ihalfy=2
     upplim=g2max+1+(ihalfy-1)*(n2-g2max-1)

     do g2=lowlim,upplim,lot

!      Find array starting address ma and ending address mb
!      modified xg 980107
!      ma=g2+(g3-1)*nd2
       ma=g2
!      Perform "lot" transforms at a time (until out of data)
!      mb=min(g2+(lot-1),upplim)+(g3-1)*nd2
       mb=min(g2+(lot-1),upplim)

!      -------------------------------------------------------------------------
!
!      Direct transformation

!      Run over all factors except the last (to ic-1), performing
!      x transform

!      Note: fortran should skip this loop if ic=1; beware "onetrip"
!      compiler option which forces each loop at least once

       do i=1,ic-1
         ntb=now(i)*bef(i)

!        Treat radix 4
         if (now(i)==4) then
           ia=0

!          First step of factor 4
           do ib=1,bef(i)
             do j=ma,mb
               r4=z(1,ia*ntb+3*bef(i)+ib,j,g3)
               s4=z(2,ia*ntb+3*bef(i)+ib,j,g3)
               r3=z(1,ia*ntb+2*bef(i)+ib,j,g3)
               s3=z(2,ia*ntb+2*bef(i)+ib,j,g3)
               r2=z(1,ia*ntb+bef(i)+ib,j,g3)
               s2=z(2,ia*ntb+bef(i)+ib,j,g3)
               r1=z(1,ia*ntb+ib,j,g3)
               s1=z(2,ia*ntb+ib,j,g3)

               r=r1 + r3
               s=r2 + r4
               z(1,ia*ntb+ib,j,g3) = r + s
               z(1,ia*ntb+2*bef(i)+ib,j,g3) = r - s
               r=r1 - r3
               s=s2 - s4
               z(1,ia*ntb+bef(i)+ib,j,g3) = r - s*ris
               z(1,ia*ntb+3*bef(i)+ib,j,g3) = r + s*ris
               r=s1 + s3
               s=s2 + s4
               z(2,ia*ntb+ib,j,g3) = r + s
               z(2,ia*ntb+2*bef(i)+ib,j,g3) = r - s
               r=s1 - s3
               s=r2 - r4
               z(2,ia*ntb+bef(i)+ib,j,g3) = r + s*ris
               z(2,ia*ntb+3*bef(i)+ib,j,g3) = r - s*ris
             end do
           end do

!          Second step of factor 4
           do ia=1,aft(i)-1
             indx=ind(ia*4*bef(i)+1)-1
             indx=indx*bef(i)
             cr2=trig(1,indx)
             ct2=trig(2,indx)
             cr3=trig(1,2*indx)
             ct3=trig(2,2*indx)
             cr4=trig(1,3*indx)
             ct4=trig(2,3*indx)
             cr4=cr4/cr2
             cr2s=cr2*ris
             do ib=1,bef(i)
               do j=ma,mb
                 r4=z(1,ia*ntb+3*bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+3*bef(i)+ib,j,g3)*ct4
                 s4=z(1,ia*ntb+3*bef(i)+ib,j,g3)*ct4 + &
&                 z(2,ia*ntb+3*bef(i)+ib,j,g3)
                 r3=z(1,ia*ntb+2*bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+2*bef(i)+ib,j,g3)*ct3
                 s3=z(1,ia*ntb+2*bef(i)+ib,j,g3)*ct3 + &
&                 z(2,ia*ntb+2*bef(i)+ib,j,g3)
                 r2=z(1,ia*ntb+bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)*ct2
                 s2=z(1,ia*ntb+bef(i)+ib,j,g3)*ct2 + &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)
                 r1=z(1,ia*ntb+ib,j,g3)
                 s1=z(2,ia*ntb+ib,j,g3)

                 r=r1 + r3*cr3
                 s=r2 + r4*cr4
                 z(1,ia*ntb+ib,j,g3) = r + s*cr2
                 z(1,ia*ntb+2*bef(i)+ib,j,g3) = r - s*cr2
                 r=r1 - r3*cr3
                 s=s2 - s4*cr4
                 z(1,ia*ntb+bef(i)+ib,j,g3) = r - s*cr2s
                 z(1,ia*ntb+3*bef(i)+ib,j,g3) = r + s*cr2s
                 r=s1 + s3*cr3
                 s=s2 + s4*cr4
                 z(2,ia*ntb+ib,j,g3) = r + s*cr2
                 z(2,ia*ntb+2*bef(i)+ib,j,g3) = r - s*cr2
                 r=s1 - s3*cr3
                 s=r2 - r4*cr4
                 z(2,ia*ntb+bef(i)+ib,j,g3) = r + s*cr2s
                 z(2,ia*ntb+3*bef(i)+ib,j,g3) = r - s*cr2s
               end do
             end do
           end do

!          Treat radix 2
         else if (now(i)==2) then
           ia=0

!          First step of factor 2
           do ib=1,bef(i)
             do j=ma,mb
               r1=z(1,ia*ntb+ib,j,g3)
               s1=z(2,ia*ntb+ib,j,g3)
               r2=z(1,ia*ntb+bef(i)+ib,j,g3)
               s2=z(2,ia*ntb+bef(i)+ib,j,g3)
               z(1,ia*ntb+ib,j,g3) =  r2 + r1
               z(2,ia*ntb+ib,j,g3) =  s2 + s1
               z(1,ia*ntb+bef(i)+ib,j,g3) = -r2 + r1
               z(2,ia*ntb+bef(i)+ib,j,g3) = -s2 + s1
             end do
           end do

!          Second step of radix 2
           do ia=1,aft(i)-1
             indx=ind(ia*2*bef(i)+1)-1
             indx=indx*bef(i)
             cr2=trig(1,indx)
             ct2=trig(2,indx)
             do ib=1,bef(i)
               do j=ma,mb
                 r1=z(1,ia*ntb+ib,j,g3)
                 s1=z(2,ia*ntb+ib,j,g3)
                 r2=z(1,ia*ntb+bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)*ct2
                 s2=z(1,ia*ntb+bef(i)+ib,j,g3)*ct2 + &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)
                 z(1,ia*ntb+ib,j,g3) =  r2*cr2 + r1
                 z(2,ia*ntb+ib,j,g3) =  s2*cr2 + s1
                 z(1,ia*ntb+bef(i)+ib,j,g3) = -r2*cr2 + r1
                 z(2,ia*ntb+bef(i)+ib,j,g3) = -s2*cr2 + s1
               end do
             end do
           end do

!          Treat radix 3
         else if (now(i)==3) then
!          .5d0*sqrt(3.d0)=0.8660254037844387d0
           ia=0
           bb=ris*0.8660254037844387d0

!          First step of radix 3
           do ib=1,bef(i)
             do j=ma,mb
               r1=z(1,ia*ntb+ib,j,g3)
               s1=z(2,ia*ntb+ib,j,g3)
               r2=z(1,ia*ntb+bef(i)+ib,j,g3)
               s2=z(2,ia*ntb+bef(i)+ib,j,g3)
               r3=z(1,ia*ntb+2*bef(i)+ib,j,g3)
               s3=z(2,ia*ntb+2*bef(i)+ib,j,g3)
               r=r2 + r3
               s=s2 + s3
               z(1,ia*ntb+ib,j,g3) = r + r1
               z(2,ia*ntb+ib,j,g3) = s + s1
               r1=r1 - r*.5d0
               s1=s1 - s*.5d0
               r2=r2-r3
               s2=s2-s3
               z(1,ia*ntb+bef(i)+ib,j,g3) = r1 - s2*bb
               z(2,ia*ntb+bef(i)+ib,j,g3) = s1 + r2*bb
               z(1,ia*ntb+2*bef(i)+ib,j,g3) = r1 + s2*bb
               z(2,ia*ntb+2*bef(i)+ib,j,g3) = s1 - r2*bb
             end do
           end do

!          Second step of radix 3
           do ia=1,aft(i)-1
             indx=ind(ia*3*bef(i)+1)-1
             indx=indx*bef(i)
             cr2=trig(1,indx)
             ct2=trig(2,indx)
             cr3=trig(1,2*indx)
             ct3=trig(2,2*indx)
             cr2=cr2/cr3
             cr3p=.5d0*cr3
             bb=ris*cr3*0.8660254037844387d0
             do ib=1,bef(i)
               do j=ma,mb
                 r1=z(1,ia*ntb+ib,j,g3)
                 s1=z(2,ia*ntb+ib,j,g3)
                 r2=z(1,ia*ntb+bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)*ct2
                 s2=z(1,ia*ntb+bef(i)+ib,j,g3)*ct2 + &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)
                 r3=z(1,ia*ntb+2*bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+2*bef(i)+ib,j,g3)*ct3
                 s3=z(1,ia*ntb+2*bef(i)+ib,j,g3)*ct3 + &
&                 z(2,ia*ntb+2*bef(i)+ib,j,g3)
                 r=cr2*r2 + r3
                 s=cr2*s2 + s3
                 z(1,ia*ntb+ib,j,g3) = r*cr3 + r1
                 z(2,ia*ntb+ib,j,g3) = s*cr3 + s1
                 r1=r1 - r*cr3p
                 s1=s1 - s*cr3p
                 r2=cr2*r2-r3
                 s2=cr2*s2-s3
                 z(1,ia*ntb+bef(i)+ib,j,g3) = r1 - s2*bb
                 z(2,ia*ntb+bef(i)+ib,j,g3) = s1 + r2*bb
                 z(1,ia*ntb+2*bef(i)+ib,j,g3) = r1 + s2*bb
                 z(2,ia*ntb+2*bef(i)+ib,j,g3) = s1 - r2*bb
               end do
             end do
           end do

!          Treat radix 5
         else if (now(i)==5) then
!          sin(2.d0*pi/5.d0)
           sin2=ris*0.9510565162951536d0
           ia=0

!          First step of radix 5
           do ib=1,bef(i)
             do j=ma,mb
               r1=z(1,ia*ntb+ib,j,g3)
               s1=z(2,ia*ntb+ib,j,g3)
               r2=z(1,ia*ntb+bef(i)+ib,j,g3)
               s2=z(2,ia*ntb+bef(i)+ib,j,g3)
               r3=z(1,ia*ntb+2*bef(i)+ib,j,g3)
               s3=z(2,ia*ntb+2*bef(i)+ib,j,g3)
               r4=z(1,ia*ntb+3*bef(i)+ib,j,g3)
               s4=z(2,ia*ntb+3*bef(i)+ib,j,g3)
               r5=z(1,ia*ntb+4*bef(i)+ib,j,g3)
               s5=z(2,ia*ntb+4*bef(i)+ib,j,g3)
               r25 = r2 + r5
               r34 = r3 + r4
               s25 = s2 - s5
               s34 = s3 - s4
               z(1,ia*ntb+ib,j,g3) = r1 + r25 + r34
               r = r1 + cos2*r25 + cos4*r34
               s = s25 + sin42*s34
               z(1,ia*ntb+bef(i)+ib,j,g3) = r - sin2*s
               z(1,ia*ntb+4*bef(i)+ib,j,g3) = r + sin2*s
               r = r1 + cos4*r25 + cos2*r34
               s = sin42*s25 - s34
               z(1,ia*ntb+2*bef(i)+ib,j,g3) = r - sin2*s
               z(1,ia*ntb+3*bef(i)+ib,j,g3) = r + sin2*s
               r25 = r2 - r5
               r34 = r3 - r4
               s25 = s2 + s5
               s34 = s3 + s4
               z(2,ia*ntb+ib,j,g3) = s1 + s25 + s34
               r = s1 + cos2*s25 + cos4*s34
               s = r25 + sin42*r34
               z(2,ia*ntb+bef(i)+ib,j,g3) = r + sin2*s
               z(2,ia*ntb+4*bef(i)+ib,j,g3) = r - sin2*s
               r = s1 + cos4*s25 + cos2*s34
               s = sin42*r25 - r34
               z(2,ia*ntb+2*bef(i)+ib,j,g3) = r + sin2*s
               z(2,ia*ntb+3*bef(i)+ib,j,g3) = r - sin2*s
             end do
           end do

!          Second step of radix 5
           do ia=1,aft(i)-1
             indx=ind(ia*5*bef(i)+1)-1
             indx=indx*bef(i)
             cr2=trig(1,indx)
             ct2=trig(2,indx)
             cr3=trig(1,2*indx)
             ct3=trig(2,2*indx)
             cr4=trig(1,3*indx)
             ct4=trig(2,3*indx)
             cr5=trig(1,4*indx)
             ct5=trig(2,4*indx)
             do ib=1,bef(i)
               do j=ma,mb
                 r1=z(1,ia*ntb+ib,j,g3)
                 s1=z(2,ia*ntb+ib,j,g3)
                 r2=cr2*(z(1,ia*ntb+bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+bef(i)+ib,j,g3)*ct2)
                 s2=cr2*(z(1,ia*ntb+bef(i)+ib,j,g3)*ct2 + &
&                 z(2,ia*ntb+bef(i)+ib,j,g3))
                 r3=cr3*(z(1,ia*ntb+2*bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+2*bef(i)+ib,j,g3)*ct3)
                 s3=cr3*(z(1,ia*ntb+2*bef(i)+ib,j,g3)*ct3 + &
&                 z(2,ia*ntb+2*bef(i)+ib,j,g3))
                 r4=z(1,ia*ntb+3*bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+3*bef(i)+ib,j,g3)*ct4
                 s4=z(1,ia*ntb+3*bef(i)+ib,j,g3)*ct4 + &
&                 z(2,ia*ntb+3*bef(i)+ib,j,g3)
                 r5=z(1,ia*ntb+4*bef(i)+ib,j,g3) - &
&                 z(2,ia*ntb+4*bef(i)+ib,j,g3)*ct5
                 s5=z(1,ia*ntb+4*bef(i)+ib,j,g3)*ct5 + &
&                 z(2,ia*ntb+4*bef(i)+ib,j,g3)
                 r25 = r2 + r5*cr5
                 r34 = r3 + r4*cr4
                 s25 = s2 - s5*cr5
                 s34 = s3 - s4*cr4
                 z(1,ia*ntb+ib,j,g3) = r1 + r25 + r34
                 r = r1 + cos2*r25 + cos4*r34
                 s = s25 + sin42*s34
                 z(1,ia*ntb+bef(i)+ib,j,g3) = r - sin2*s
                 z(1,ia*ntb+4*bef(i)+ib,j,g3) = r + sin2*s
                 r = r1 + cos4*r25 + cos2*r34
                 s = sin42*s25 - s34
                 z(1,ia*ntb+2*bef(i)+ib,j,g3) = r - sin2*s
                 z(1,ia*ntb+3*bef(i)+ib,j,g3) = r + sin2*s
                 r25 = r2 - r5*cr5
                 r34 = r3 - r4*cr4
                 s25 = s2 + s5*cr5
                 s34 = s3 + s4*cr4
                 z(2,ia*ntb+ib,j,g3) = s1 + s25 + s34
                 r = s1 + cos2*s25 + cos4*s34
                 s = r25 + sin42*r34
                 z(2,ia*ntb+bef(i)+ib,j,g3) = r + sin2*s
                 z(2,ia*ntb+4*bef(i)+ib,j,g3) = r - sin2*s
                 r = s1 + cos4*s25 + cos2*s34
                 s = sin42*r25 - r34
                 z(2,ia*ntb+2*bef(i)+ib,j,g3) = r + sin2*s
                 z(2,ia*ntb+3*bef(i)+ib,j,g3) = r - sin2*s
               end do
             end do
           end do

         else
!          All radices treated
           MSG_BUG('called with factors other than 2, 3, and 5')
         end if

       end do  ! End of direct transformation (loop over ic)

!      -----------------------------------------------------------------

!      Bitreversal
!      Perform bit reversal on last factor of transformation

!      Treat radix 4
       if (now(ic)==4) then
         ia=0

!        First step of radix 4
         do j=ma,mb
           r4=z(1,ia*4+4,j,g3)
           s4=z(2,ia*4+4,j,g3)
           r3=z(1,ia*4+3,j,g3)
           s3=z(2,ia*4+3,j,g3)
           r2=z(1,ia*4+2,j,g3)
           s2=z(2,ia*4+2,j,g3)
           r1=z(1,ia*4+1,j,g3)
           s1=z(2,ia*4+1,j,g3)

           r=r1 + r3
           s=r2 + r4
           zbr(1,ind(ia*4+1),j,g3) = r + s
           zbr(1,ind(ia*4+3),j,g3) = r - s
           r=r1 - r3
           s=s2 - s4
           zbr(1,ind(ia*4+2),j,g3) = r - s*ris
           zbr(1,ind(ia*4+4),j,g3) = r + s*ris
           r=s1 + s3
           s=s2 + s4
           zbr(2,ind(ia*4+1),j,g3) = r + s
           zbr(2,ind(ia*4+3),j,g3) = r - s
           r=s1 - s3
           s=r2 - r4
           zbr(2,ind(ia*4+2),j,g3) = r + s*ris
           zbr(2,ind(ia*4+4),j,g3) = r - s*ris
         end do

!        Second step of radix 4
         do ia=1,aft(ic)-1
           indx=ind(ia*4+1)-1
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           cr3=trig(1,2*indx)
           ct3=trig(2,2*indx)
           cr4=trig(1,3*indx)
           ct4=trig(2,3*indx)
           cr4=cr4/cr2
           cr2s=cr2*ris
           do j=ma,mb
             r4=z(1,ia*4+4,j,g3) - z(2,ia*4+4,j,g3)*ct4
             s4=z(1,ia*4+4,j,g3)*ct4 + z(2,ia*4+4,j,g3)
             r3=z(1,ia*4+3,j,g3) - z(2,ia*4+3,j,g3)*ct3
             s3=z(1,ia*4+3,j,g3)*ct3 + z(2,ia*4+3,j,g3)
             r2=z(1,ia*4+2,j,g3) - z(2,ia*4+2,j,g3)*ct2
             s2=z(1,ia*4+2,j,g3)*ct2 + z(2,ia*4+2,j,g3)
             r1=z(1,ia*4+1,j,g3)
             s1=z(2,ia*4+1,j,g3)

             r=r1 + r3*cr3
             s=r2 + r4*cr4
             zbr(1,ind(ia*4+1),j,g3) = r + s*cr2
             zbr(1,ind(ia*4+3),j,g3) = r - s*cr2
             r=r1 - r3*cr3
             s=s2 - s4*cr4
             zbr(1,ind(ia*4+2),j,g3) = r - s*cr2s
             zbr(1,ind(ia*4+4),j,g3) = r + s*cr2s
             r=s1 + s3*cr3
             s=s2 + s4*cr4
             zbr(2,ind(ia*4+1),j,g3) = r + s*cr2
             zbr(2,ind(ia*4+3),j,g3) = r - s*cr2
             r=s1 - s3*cr3
             s=r2 - r4*cr4
             zbr(2,ind(ia*4+2),j,g3) = r + s*cr2s
             zbr(2,ind(ia*4+4),j,g3) = r - s*cr2s
           end do
         end do

!        Treat radix 2
       else if (now(ic)==2) then

         ia=0

!        First step of radix 2
         do j=ma,mb
           r1=z(1,ia*2+1,j,g3)
           s1=z(2,ia*2+1,j,g3)
           r2=z(1,ia*2+2,j,g3)
           s2=z(2,ia*2+2,j,g3)
           zbr(1,ind(ia*2+1),j,g3) =  r2 + r1
           zbr(2,ind(ia*2+1),j,g3) =  s2 + s1
           zbr(1,ind(ia*2+2),j,g3) = -r2 + r1
           zbr(2,ind(ia*2+2),j,g3) = -s2 + s1
         end do

!        Second step of radix 2
         do ia=1,aft(ic)-1
           indx=ind(ia*2+1)-1
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           do j=ma,mb
             r1=z(1,ia*2+1,j,g3)
             s1=z(2,ia*2+1,j,g3)
             r2=z(1,ia*2+2,j,g3) - z(2,ia*2+2,j,g3)*ct2
             s2=z(1,ia*2+2,j,g3)*ct2 + z(2,ia*2+2,j,g3)
             zbr(1,ind(ia*2+1),j,g3) =  r2*cr2 + r1
             zbr(2,ind(ia*2+1),j,g3) =  s2*cr2 + s1
             zbr(1,ind(ia*2+2),j,g3) = -r2*cr2 + r1
             zbr(2,ind(ia*2+2),j,g3) = -s2*cr2 + s1
           end do
         end do

!        Treat radix 3
       else if (now(ic)==3) then
!        radix 3
!        .5d0*sqrt(3.d0)=0.8660254037844387d0
         ia=0
         bb=ris*0.8660254037844387d0

!        First step of radix 3
         do j=ma,mb
           r1=z(1,ia*3+1,j,g3)
           s1=z(2,ia*3+1,j,g3)
           r2=z(1,ia*3+2,j,g3)
           s2=z(2,ia*3+2,j,g3)
           r3=z(1,ia*3+3,j,g3)
           s3=z(2,ia*3+3,j,g3)
           r=r2 + r3
           s=s2 + s3
           zbr(1,ind(ia*3+1),j,g3) = r + r1
           zbr(2,ind(ia*3+1),j,g3) = s + s1
           r1=r1 - r*.5d0
           s1=s1 - s*.5d0
           r2=r2-r3
           s2=s2-s3
           zbr(1,ind(ia*3+2),j,g3) = r1 - s2*bb
           zbr(2,ind(ia*3+2),j,g3) = s1 + r2*bb
           zbr(1,ind(ia*3+3),j,g3) = r1 + s2*bb
           zbr(2,ind(ia*3+3),j,g3) = s1 - r2*bb
         end do

         do ia=1,aft(ic)-1
           indx=ind(ia*3+1)-1
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           cr3=trig(1,2*indx)
           ct3=trig(2,2*indx)
           cr2=cr2/cr3
           cr3p=.5d0*cr3
           bb=ris*cr3*0.8660254037844387d0
           do j=ma,mb
             r1=z(1,ia*3+1,j,g3)
             s1=z(2,ia*3+1,j,g3)
             r2=z(1,ia*3+2,j,g3) - z(2,ia*3+2,j,g3)*ct2
             s2=z(1,ia*3+2,j,g3)*ct2 + z(2,ia*3+2,j,g3)
             r3=z(1,ia*3+3,j,g3) - z(2,ia*3+3,j,g3)*ct3
             s3=z(1,ia*3+3,j,g3)*ct3 + z(2,ia*3+3,j,g3)
             r=cr2*r2 + r3
             s=cr2*s2 + s3
             zbr(1,ind(ia*3+1),j,g3) = r*cr3 + r1
             zbr(2,ind(ia*3+1),j,g3) = s*cr3 + s1
             r1=r1 - r*cr3p
             s1=s1 - s*cr3p
             r2=cr2*r2-r3
             s2=cr2*s2-s3
             zbr(1,ind(ia*3+2),j,g3) = r1 - s2*bb
             zbr(2,ind(ia*3+2),j,g3) = s1 + r2*bb
             zbr(1,ind(ia*3+3),j,g3) = r1 + s2*bb
             zbr(2,ind(ia*3+3),j,g3) = s1 - r2*bb
           end do
         end do

!        Treat radix 5
       else if (now(ic)==5) then
!        radix 5
!        sin(2.d0*pi/5.d0)
         sin2=ris*0.9510565162951536d0
         ia=0

!        First step of radix 5
         do j=ma,mb
           r1=z(1,ia*5+1,j,g3)
           s1=z(2,ia*5+1,j,g3)
           r2=z(1,ia*5+2,j,g3)
           s2=z(2,ia*5+2,j,g3)
           r3=z(1,ia*5+3,j,g3)
           s3=z(2,ia*5+3,j,g3)
           r4=z(1,ia*5+4,j,g3)
           s4=z(2,ia*5+4,j,g3)
           r5=z(1,ia*5+5,j,g3)
           s5=z(2,ia*5+5,j,g3)
           r25 = r2 + r5
           r34 = r3 + r4
           s25 = s2 - s5
           s34 = s3 - s4
           zbr(1,ind(ia*5+1),j,g3) = r1 + r25 + r34
           r = r1 + cos2*r25 + cos4*r34
           s = s25 + sin42*s34
           zbr(1,ind(ia*5+2),j,g3) = r - sin2*s
           zbr(1,ind(ia*5+5),j,g3) = r + sin2*s
           r = r1 + cos4*r25 + cos2*r34
           s = sin42*s25 - s34
           zbr(1,ind(ia*5+3),j,g3) = r - sin2*s
           zbr(1,ind(ia*5+4),j,g3) = r + sin2*s
           r25 = r2 - r5
           r34 = r3 - r4
           s25 = s2 + s5
           s34 = s3 + s4
           zbr(2,ind(ia*5+1),j,g3) = s1 + s25 + s34
           r = s1 + cos2*s25 + cos4*s34
           s = r25 + sin42*r34
           zbr(2,ind(ia*5+2),j,g3) = r + sin2*s
           zbr(2,ind(ia*5+5),j,g3) = r - sin2*s
           r = s1 + cos4*s25 + cos2*s34
           s = sin42*r25 - r34
           zbr(2,ind(ia*5+3),j,g3) = r + sin2*s
           zbr(2,ind(ia*5+4),j,g3) = r - sin2*s
         end do

!        Second step of radix 5
         do ia=1,aft(ic)-1
           indx=ind(ia*5+1)-1
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           cr3=trig(1,2*indx)
           ct3=trig(2,2*indx)
           cr4=trig(1,3*indx)
           ct4=trig(2,3*indx)
           cr5=trig(1,4*indx)
           ct5=trig(2,4*indx)
           do j=ma,mb
             r1=z(1,ia*5+1,j,g3)
             s1=z(2,ia*5+1,j,g3)
             r2=cr2*(z(1,ia*5+2,j,g3) - z(2,ia*5+2,j,g3)*ct2)
             s2=cr2*(z(1,ia*5+2,j,g3)*ct2 + z(2,ia*5+2,j,g3))
             r3=cr3*(z(1,ia*5+3,j,g3) - z(2,ia*5+3,j,g3)*ct3)
             s3=cr3*(z(1,ia*5+3,j,g3)*ct3 + z(2,ia*5+3,j,g3))
             r4=z(1,ia*5+4,j,g3) - z(2,ia*5+4,j,g3)*ct4
             s4=z(1,ia*5+4,j,g3)*ct4 + z(2,ia*5+4,j,g3)
             r5=z(1,ia*5+5,j,g3) - z(2,ia*5+5,j,g3)*ct5
             s5=z(1,ia*5+5,j,g3)*ct5 + z(2,ia*5+5,j,g3)
             r25 = r2 + r5*cr5
             r34 = r3 + r4*cr4
             s25 = s2 - s5*cr5
             s34 = s3 - s4*cr4
             zbr(1,ind(ia*5+1),j,g3) = r1 + r25 + r34
             r = r1 + cos2*r25 + cos4*r34
             s = s25 + sin42*s34
             zbr(1,ind(ia*5+2),j,g3) = r - sin2*s
             zbr(1,ind(ia*5+5),j,g3) = r + sin2*s
             r = r1 + cos4*r25 + cos2*r34
             s = sin42*s25 - s34
             zbr(1,ind(ia*5+3),j,g3) = r - sin2*s
             zbr(1,ind(ia*5+4),j,g3) = r + sin2*s
             r25 = r2 - r5*cr5
             r34 = r3 - r4*cr4
             s25 = s2 + s5*cr5
             s34 = s3 + s4*cr4
             zbr(2,ind(ia*5+1),j,g3) = s1 + s25 + s34
             r = s1 + cos2*s25 + cos4*s34
             s = r25 + sin42*r34
             zbr(2,ind(ia*5+2),j,g3) = r + sin2*s
             zbr(2,ind(ia*5+5),j,g3) = r - sin2*s
             r = s1 + cos4*s25 + cos2*s34
             s = sin42*r25 - r34
             zbr(2,ind(ia*5+3),j,g3) = r + sin2*s
             zbr(2,ind(ia*5+4),j,g3) = r - sin2*s
           end do
         end do

       else
!        All radices are treated
         MSG_BUG('called with factors other than 2, 3, and 5')
       end if

!      End of bit reversal

!      -------------------------------------------------------------------
     end do
   end do
 end do
!$OMP END PARALLEL DO

end subroutine sg_fftpx
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fftx
!! NAME
!! sg_fftx
!!
!! FUNCTION
!! This subroutine is called by the 3-dimensional fft to conduct the
!! "x" transforms for all y and z.
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  mfac = maximum number of factors in 1D FFTs
!!  mg = maximum length of 1D FFTs
!!  nd1=first dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd2=second dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd3=third dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  n2,n3=actual length of y and z transforms
!!  z(2,nd1,nd2,nd3)=INPUT array; destroyed by transformation
!!  trig, aft, now, bef, ind=provided by previous call to ctrig
!!   Note that in this routine (and in ctrig) the values in array trig are
!!   actually cos and tan, not cos and sin.  Use of tan allows advantageous
!!   use of FMA on the ibm rs6000.
!!  ris=sign of exponential in transform (should be 1 or -1; real)
!!  ic=number of (radix) factors of x transform length (from ctrig)
!!
!! OUTPUT
!!  zbr(2,nd1,nd2,nd3)=OUTPUT transformed array; no scaling applied
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine blocks the x transforms
!! so that all transforms under consideration at one step fit within
!! the cache memory, which is crucial for optimal performance.
!! The blocking factor is set by parameter "fftcache" below, which should
!! be adjusted to be somewhat smaller (say 3/4) than the actual cache size
!! of the machine.
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fftx(fftcache,mfac,mg,nd1,nd2,nd3,n2,n3,z,zbr,&
& trig,aft,now,bef,ris,ind,ic)

 implicit none

!Arguments ------------------------------------
!Dimensions of aft, now, bef, ind, and trig should agree with
!those in subroutine ctrig.
!scalars
 integer,intent(in) :: fftcache,ic,mfac,mg,n2,n3,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 integer,intent(in) :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp),intent(in) :: trig(2,mg)
 real(dp),intent(inout) :: z(2,nd1,nd2,nd3),zbr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer :: i,i3,ia,ib,indx,j,jj,lot,ma,mb,ntb
 real(dp),parameter :: cos2=0.3090169943749474d0   !cos(2.d0*pi/5.d0)
 real(dp),parameter :: cos4=-0.8090169943749474d0  !cos(4.d0*pi/5.d0)
 real(dp),parameter :: sin42=0.6180339887498948d0  !sin(4.d0*pi/5.d0)/sin(2.d0*pi/5.d0)
 real(dp) :: bb,cr2,cr2s,cr3,cr3p,cr4,cr5,ct2,ct3,ct4,ct5
 real(dp) :: factor,r,r1,r2,r25,r3,r34,r4,r5,s,sin2,s1,s2,s25,s3,s34,s4,s5

! *************************************************************************

!Do x transforms in blocks of size "lot" which is set by how
!many x transform arrays (of size nd1 each) fit into the nominal
!cache size "fftcache".
 factor=0.75d0
 lot=(fftcache*factor*1000d0)/(nd1*8*2)

!XG : due to the dimension problems on the P6, I have slightly
!modified this part of the code, with an external loop
!on n3 ...
!Modifications are indicated explicitely, or
!are related to the increase of the number of dimensions of z and
!zbr ...

 factor=0.75d0
 lot=(fftcache*factor*1000d0)/(nd1*8*2)
 if(lot.lt.1) lot=1 ! this may happen for very large cells
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(aft,bef,ic,ind,lot,n2,n3,now,ris,trig,z,zbr)
 do i3=1,n3
   do jj=1,n2,lot
!    end of modification

!    For each jj, ma and mb give starting and ending addresses for fft
!    ma starts where we left off after last block
     ma=jj
!    mb runs to the end of the block or else to the end of the data
!    modified XG 980107
!    mb=min(jj+(lot-1),n23)
     mb=min(jj+(lot-1),n2)

!    Run over all factors except the last (to ic-1), performing
!    x transform

!    Note: fortran should skip this loop if ic=1; beware "onetrip"
!    compiler option which forces each loop at least once

!    ------------------------------------------------------------------------

!    Direct transformation (to be followed by bit reversal)

     do i=1,ic-1
       ntb=now(i)*bef(i)
!      radix 4

!      Treat radix 4
       if (now(i)==4) then
         ia=0

!        First step of factor 4
         do ib=1,bef(i)
           do j=ma,mb
             r4=z(1,ia*ntb+3*bef(i)+ib,j,i3)
             s4=z(2,ia*ntb+3*bef(i)+ib,j,i3)
             r3=z(1,ia*ntb+2*bef(i)+ib,j,i3)
             s3=z(2,ia*ntb+2*bef(i)+ib,j,i3)
             r2=z(1,ia*ntb+bef(i)+ib,j,i3)
             s2=z(2,ia*ntb+bef(i)+ib,j,i3)
             r1=z(1,ia*ntb+ib,j,i3)
             s1=z(2,ia*ntb+ib,j,i3)

             r=r1 + r3
             s=r2 + r4
             z(1,ia*ntb+ib,j,i3) = r + s
             z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - s
             r=r1 - r3
             s=s2 - s4
             z(1,ia*ntb+bef(i)+ib,j,i3) = r - s*ris
             z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + s*ris
             r=s1 + s3
             s=s2 + s4
             z(2,ia*ntb+ib,j,i3) = r + s
             z(2,ia*ntb+2*bef(i)+ib,j,i3) = r - s
             r=s1 - s3
             s=r2 - r4
             z(2,ia*ntb+bef(i)+ib,j,i3) = r + s*ris
             z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - s*ris
           end do
         end do

!        Second step of factor 4
         do ia=1,aft(i)-1
           indx=ind(ia*4*bef(i)+1)-1
           indx=indx*bef(i)
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           cr3=trig(1,2*indx)
           ct3=trig(2,2*indx)
           cr4=trig(1,3*indx)
           ct4=trig(2,3*indx)
           cr4=cr4/cr2
           cr2s=cr2*ris
           do ib=1,bef(i)
             do j=ma,mb
               r4=z(1,ia*ntb+3*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+3*bef(i)+ib,j,i3)*ct4
               s4=z(1,ia*ntb+3*bef(i)+ib,j,i3)*ct4 + &
&               z(2,ia*ntb+3*bef(i)+ib,j,i3)
               r3=z(1,ia*ntb+2*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3)*ct3
               s3=z(1,ia*ntb+2*bef(i)+ib,j,i3)*ct3 + &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3)
               r2=z(1,ia*ntb+bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+bef(i)+ib,j,i3)*ct2
               s2=z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&               z(2,ia*ntb+bef(i)+ib,j,i3)
               r1=z(1,ia*ntb+ib,j,i3)
               s1=z(2,ia*ntb+ib,j,i3)

               r=r1 + r3*cr3
               s=r2 + r4*cr4
               z(1,ia*ntb+ib,j,i3) = r + s*cr2
               z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - s*cr2
               r=r1 - r3*cr3
               s=s2 - s4*cr4
               z(1,ia*ntb+bef(i)+ib,j,i3) = r - s*cr2s
               z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + s*cr2s
               r=s1 + s3*cr3
               s=s2 + s4*cr4
               z(2,ia*ntb+ib,j,i3) = r + s*cr2
               z(2,ia*ntb+2*bef(i)+ib,j,i3) = r - s*cr2
               r=s1 - s3*cr3
               s=r2 - r4*cr4
               z(2,ia*ntb+bef(i)+ib,j,i3) = r + s*cr2s
               z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - s*cr2s
             end do
           end do
         end do

!        Treat radix 2
       else if (now(i)==2) then
         ia=0

!        First step of factor 2
         do ib=1,bef(i)
           do j=ma,mb
             r1=z(1,ia*ntb+ib,j,i3)
             s1=z(2,ia*ntb+ib,j,i3)
             r2=z(1,ia*ntb+bef(i)+ib,j,i3)
             s2=z(2,ia*ntb+bef(i)+ib,j,i3)
             z(1,ia*ntb+ib,j,i3) =  r2 + r1
             z(2,ia*ntb+ib,j,i3) =  s2 + s1
             z(1,ia*ntb+bef(i)+ib,j,i3) = -r2 + r1
             z(2,ia*ntb+bef(i)+ib,j,i3) = -s2 + s1
           end do
         end do

!        Second step of factor 2
         do ia=1,aft(i)-1
           indx=ind(ia*2*bef(i)+1)-1
           indx=indx*bef(i)
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           do ib=1,bef(i)
             do j=ma,mb
               r1=z(1,ia*ntb+ib,j,i3)
               s1=z(2,ia*ntb+ib,j,i3)
               r2=z(1,ia*ntb+bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+bef(i)+ib,j,i3)*ct2
               s2=z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&               z(2,ia*ntb+bef(i)+ib,j,i3)
               z(1,ia*ntb+ib,j,i3) =  r2*cr2 + r1
               z(2,ia*ntb+ib,j,i3) =  s2*cr2 + s1
               z(1,ia*ntb+bef(i)+ib,j,i3) = -r2*cr2 + r1
               z(2,ia*ntb+bef(i)+ib,j,i3) = -s2*cr2 + s1
             end do
           end do
         end do

!        Treat radix 3
       else if (now(i)==3) then
!        .5d0*sqrt(3.d0)=0.8660254037844387d0
         ia=0
         bb=ris*0.8660254037844387d0

!        First step of factor 3
         do ib=1,bef(i)
           do j=ma,mb
             r1=z(1,ia*ntb+ib,j,i3)
             s1=z(2,ia*ntb+ib,j,i3)
             r2=z(1,ia*ntb+bef(i)+ib,j,i3)
             s2=z(2,ia*ntb+bef(i)+ib,j,i3)
             r3=z(1,ia*ntb+2*bef(i)+ib,j,i3)
             s3=z(2,ia*ntb+2*bef(i)+ib,j,i3)
             r=r2 + r3
             s=s2 + s3
             z(1,ia*ntb+ib,j,i3) = r + r1
             z(2,ia*ntb+ib,j,i3) = s + s1
             r1=r1 - r*.5d0
             s1=s1 - s*.5d0
             r2=r2-r3
             s2=s2-s3
             z(1,ia*ntb+bef(i)+ib,j,i3) = r1 - s2*bb
             z(2,ia*ntb+bef(i)+ib,j,i3) = s1 + r2*bb
             z(1,ia*ntb+2*bef(i)+ib,j,i3) = r1 + s2*bb
             z(2,ia*ntb+2*bef(i)+ib,j,i3) = s1 - r2*bb
           end do
         end do

!        Second step of factor 3
         do ia=1,aft(i)-1
           indx=ind(ia*3*bef(i)+1)-1
           indx=indx*bef(i)
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           cr3=trig(1,2*indx)
           ct3=trig(2,2*indx)
           cr2=cr2/cr3
           cr3p=.5d0*cr3
           bb=ris*cr3*0.8660254037844387d0
           do ib=1,bef(i)
             do j=ma,mb
               r1=z(1,ia*ntb+ib,j,i3)
               s1=z(2,ia*ntb+ib,j,i3)
               r2=z(1,ia*ntb+bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+bef(i)+ib,j,i3)*ct2
               s2=z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&               z(2,ia*ntb+bef(i)+ib,j,i3)
               r3=z(1,ia*ntb+2*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3)*ct3
               s3=z(1,ia*ntb+2*bef(i)+ib,j,i3)*ct3 + &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3)
               r=cr2*r2 + r3
               s=cr2*s2 + s3
               z(1,ia*ntb+ib,j,i3) = r*cr3 + r1
               z(2,ia*ntb+ib,j,i3) = s*cr3 + s1
               r1=r1 - r*cr3p
               s1=s1 - s*cr3p
               r2=cr2*r2-r3
               s2=cr2*s2-s3
               z(1,ia*ntb+bef(i)+ib,j,i3) = r1 - s2*bb
               z(2,ia*ntb+bef(i)+ib,j,i3) = s1 + r2*bb
               z(1,ia*ntb+2*bef(i)+ib,j,i3) = r1 + s2*bb
               z(2,ia*ntb+2*bef(i)+ib,j,i3) = s1 - r2*bb
             end do
           end do
         end do

!        Treat radix 5
       else if (now(i)==5) then
!        sin(2.d0*pi/5.d0)
         sin2=ris*0.9510565162951536d0
         ia=0

!        First step of factor 5
         do ib=1,bef(i)
           do j=ma,mb
             r1=z(1,ia*ntb+ib,j,i3)
             s1=z(2,ia*ntb+ib,j,i3)
             r2=z(1,ia*ntb+bef(i)+ib,j,i3)
             s2=z(2,ia*ntb+bef(i)+ib,j,i3)
             r3=z(1,ia*ntb+2*bef(i)+ib,j,i3)
             s3=z(2,ia*ntb+2*bef(i)+ib,j,i3)
             r4=z(1,ia*ntb+3*bef(i)+ib,j,i3)
             s4=z(2,ia*ntb+3*bef(i)+ib,j,i3)
             r5=z(1,ia*ntb+4*bef(i)+ib,j,i3)
             s5=z(2,ia*ntb+4*bef(i)+ib,j,i3)
             r25 = r2 + r5
             r34 = r3 + r4
             s25 = s2 - s5
             s34 = s3 - s4
             z(1,ia*ntb+ib,j,i3) = r1 + r25 + r34
             r = r1 + cos2*r25 + cos4*r34
             s = s25 + sin42*s34
             z(1,ia*ntb+bef(i)+ib,j,i3) = r - sin2*s
             z(1,ia*ntb+4*bef(i)+ib,j,i3) = r + sin2*s
             r = r1 + cos4*r25 + cos2*r34
             s = sin42*s25 - s34
             z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - sin2*s
             z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + sin2*s
             r25 = r2 - r5
             r34 = r3 - r4
             s25 = s2 + s5
             s34 = s3 + s4
             z(2,ia*ntb+ib,j,i3) = s1 + s25 + s34
             r = s1 + cos2*s25 + cos4*s34
             s = r25 + sin42*r34
             z(2,ia*ntb+bef(i)+ib,j,i3) = r + sin2*s
             z(2,ia*ntb+4*bef(i)+ib,j,i3) = r - sin2*s
             r = s1 + cos4*s25 + cos2*s34
             s = sin42*r25 - r34
             z(2,ia*ntb+2*bef(i)+ib,j,i3) = r + sin2*s
             z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - sin2*s
           end do
         end do

!        Second step of factor 5
         do ia=1,aft(i)-1
           indx=ind(ia*5*bef(i)+1)-1
           indx=indx*bef(i)
           cr2=trig(1,indx)
           ct2=trig(2,indx)
           cr3=trig(1,2*indx)
           ct3=trig(2,2*indx)
           cr4=trig(1,3*indx)
           ct4=trig(2,3*indx)
           cr5=trig(1,4*indx)
           ct5=trig(2,4*indx)
           do ib=1,bef(i)
             do j=ma,mb
               r1=z(1,ia*ntb+ib,j,i3)
               s1=z(2,ia*ntb+ib,j,i3)
               r2=cr2*(z(1,ia*ntb+bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+bef(i)+ib,j,i3)*ct2)
               s2=cr2*(z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&               z(2,ia*ntb+bef(i)+ib,j,i3))
               r3=cr3*(z(1,ia*ntb+2*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3)*ct3)
               s3=cr3*(z(1,ia*ntb+2*bef(i)+ib,j,i3)*ct3 + &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3))
               r4=z(1,ia*ntb+3*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+3*bef(i)+ib,j,i3)*ct4
               s4=z(1,ia*ntb+3*bef(i)+ib,j,i3)*ct4 + &
&               z(2,ia*ntb+3*bef(i)+ib,j,i3)
               r5=z(1,ia*ntb+4*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+4*bef(i)+ib,j,i3)*ct5
               s5=z(1,ia*ntb+4*bef(i)+ib,j,i3)*ct5 + &
&               z(2,ia*ntb+4*bef(i)+ib,j,i3)
               r25 = r2 + r5*cr5
               r34 = r3 + r4*cr4
               s25 = s2 - s5*cr5
               s34 = s3 - s4*cr4
               z(1,ia*ntb+ib,j,i3) = r1 + r25 + r34
               r = r1 + cos2*r25 + cos4*r34
               s = s25 + sin42*s34
               z(1,ia*ntb+bef(i)+ib,j,i3) = r - sin2*s
               z(1,ia*ntb+4*bef(i)+ib,j,i3) = r + sin2*s
               r = r1 + cos4*r25 + cos2*r34
               s = sin42*s25 - s34
               z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - sin2*s
               z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + sin2*s
               r25 = r2 - r5*cr5
               r34 = r3 - r4*cr4
               s25 = s2 + s5*cr5
               s34 = s3 + s4*cr4
               z(2,ia*ntb+ib,j,i3) = s1 + s25 + s34
               r = s1 + cos2*s25 + cos4*s34
               s = r25 + sin42*r34
               z(2,ia*ntb+bef(i)+ib,j,i3) = r + sin2*s
               z(2,ia*ntb+4*bef(i)+ib,j,i3) = r - sin2*s
               r = s1 + cos4*s25 + cos2*s34
               s = sin42*r25 - r34
               z(2,ia*ntb+2*bef(i)+ib,j,i3) = r + sin2*s
               z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - sin2*s
             end do
           end do
         end do

       else
!        All factors have been treated
         MSG_BUG('called with factors other than 2, 3, and 5')
       end if

     end do

!    ---------------------------------------------------------------

!    bitreversal

!    Perform bit reversal on last factor of transformation

!    Treat factor 4
     if (now(ic)==4) then
!      radix 4
       ia=0

!      First step of factor 4
       do j=ma,mb
         r4=z(1,ia*4+4,j,i3)
         s4=z(2,ia*4+4,j,i3)
         r3=z(1,ia*4+3,j,i3)
         s3=z(2,ia*4+3,j,i3)
         r2=z(1,ia*4+2,j,i3)
         s2=z(2,ia*4+2,j,i3)
         r1=z(1,ia*4+1,j,i3)
         s1=z(2,ia*4+1,j,i3)

         r=r1 + r3
         s=r2 + r4
         zbr(1,ind(ia*4+1),j,i3) = r + s
         zbr(1,ind(ia*4+3),j,i3) = r - s
         r=r1 - r3
         s=s2 - s4
         zbr(1,ind(ia*4+2),j,i3) = r - s*ris
         zbr(1,ind(ia*4+4),j,i3) = r + s*ris
         r=s1 + s3
         s=s2 + s4
         zbr(2,ind(ia*4+1),j,i3) = r + s
         zbr(2,ind(ia*4+3),j,i3) = r - s
         r=s1 - s3
         s=r2 - r4
         zbr(2,ind(ia*4+2),j,i3) = r + s*ris
         zbr(2,ind(ia*4+4),j,i3) = r - s*ris
       end do

!      Second step of factor 4
       do ia=1,aft(ic)-1
         indx=ind(ia*4+1)-1
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         cr3=trig(1,2*indx)
         ct3=trig(2,2*indx)
         cr4=trig(1,3*indx)
         ct4=trig(2,3*indx)
         cr4=cr4/cr2
         cr2s=cr2*ris
         do j=ma,mb
           r4=z(1,ia*4+4,j,i3) - z(2,ia*4+4,j,i3)*ct4
           s4=z(1,ia*4+4,j,i3)*ct4 + z(2,ia*4+4,j,i3)
           r3=z(1,ia*4+3,j,i3) - z(2,ia*4+3,j,i3)*ct3
           s3=z(1,ia*4+3,j,i3)*ct3 + z(2,ia*4+3,j,i3)
           r2=z(1,ia*4+2,j,i3) - z(2,ia*4+2,j,i3)*ct2
           s2=z(1,ia*4+2,j,i3)*ct2 + z(2,ia*4+2,j,i3)
           r1=z(1,ia*4+1,j,i3)
           s1=z(2,ia*4+1,j,i3)

           r=r1 + r3*cr3
           s=r2 + r4*cr4
           zbr(1,ind(ia*4+1),j,i3) = r + s*cr2
           zbr(1,ind(ia*4+3),j,i3) = r - s*cr2
           r=r1 - r3*cr3
           s=s2 - s4*cr4
           zbr(1,ind(ia*4+2),j,i3) = r - s*cr2s
           zbr(1,ind(ia*4+4),j,i3) = r + s*cr2s
           r=s1 + s3*cr3
           s=s2 + s4*cr4
           zbr(2,ind(ia*4+1),j,i3) = r + s*cr2
           zbr(2,ind(ia*4+3),j,i3) = r - s*cr2
           r=s1 - s3*cr3
           s=r2 - r4*cr4
           zbr(2,ind(ia*4+2),j,i3) = r + s*cr2s
           zbr(2,ind(ia*4+4),j,i3) = r - s*cr2s
         end do
       end do

!      Treat factor 2
     else if (now(ic)==2) then
!      radix 2
       ia=0

!      First step of factor 2
       do j=ma,mb
         r1=z(1,ia*2+1,j,i3)
         s1=z(2,ia*2+1,j,i3)
         r2=z(1,ia*2+2,j,i3)
         s2=z(2,ia*2+2,j,i3)
         zbr(1,ind(ia*2+1),j,i3) =  r2 + r1
         zbr(2,ind(ia*2+1),j,i3) =  s2 + s1
         zbr(1,ind(ia*2+2),j,i3) = -r2 + r1
         zbr(2,ind(ia*2+2),j,i3) = -s2 + s1
       end do

!      Second step of factor 2
       do ia=1,aft(ic)-1
         indx=ind(ia*2+1)-1
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         do j=ma,mb
           r1=z(1,ia*2+1,j,i3)
           s1=z(2,ia*2+1,j,i3)
           r2=z(1,ia*2+2,j,i3) - z(2,ia*2+2,j,i3)*ct2
           s2=z(1,ia*2+2,j,i3)*ct2 + z(2,ia*2+2,j,i3)
           zbr(1,ind(ia*2+1),j,i3) =  r2*cr2 + r1
           zbr(2,ind(ia*2+1),j,i3) =  s2*cr2 + s1
           zbr(1,ind(ia*2+2),j,i3) = -r2*cr2 + r1
           zbr(2,ind(ia*2+2),j,i3) = -s2*cr2 + s1
         end do
       end do

!      Treat factor 3
     else if (now(ic)==3) then
!      radix 3
!      .5d0*sqrt(3.d0)=0.8660254037844387d0
       ia=0
       bb=ris*0.8660254037844387d0

!      First step of factor 3
       do j=ma,mb
         r1=z(1,ia*3+1,j,i3)
         s1=z(2,ia*3+1,j,i3)
         r2=z(1,ia*3+2,j,i3)
         s2=z(2,ia*3+2,j,i3)
         r3=z(1,ia*3+3,j,i3)
         s3=z(2,ia*3+3,j,i3)
         r=r2 + r3
         s=s2 + s3
         zbr(1,ind(ia*3+1),j,i3) = r + r1
         zbr(2,ind(ia*3+1),j,i3) = s + s1
         r1=r1 - r*.5d0
         s1=s1 - s*.5d0
         r2=r2-r3
         s2=s2-s3
         zbr(1,ind(ia*3+2),j,i3) = r1 - s2*bb
         zbr(2,ind(ia*3+2),j,i3) = s1 + r2*bb
         zbr(1,ind(ia*3+3),j,i3) = r1 + s2*bb
         zbr(2,ind(ia*3+3),j,i3) = s1 - r2*bb
       end do

!      Second step of factor 3
       do ia=1,aft(ic)-1
         indx=ind(ia*3+1)-1
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         cr3=trig(1,2*indx)
         ct3=trig(2,2*indx)
         cr2=cr2/cr3
         cr3p=.5d0*cr3
         bb=ris*cr3*0.8660254037844387d0
         do j=ma,mb
           r1=z(1,ia*3+1,j,i3)
           s1=z(2,ia*3+1,j,i3)
           r2=z(1,ia*3+2,j,i3) - z(2,ia*3+2,j,i3)*ct2
           s2=z(1,ia*3+2,j,i3)*ct2 + z(2,ia*3+2,j,i3)
           r3=z(1,ia*3+3,j,i3) - z(2,ia*3+3,j,i3)*ct3
           s3=z(1,ia*3+3,j,i3)*ct3 + z(2,ia*3+3,j,i3)
           r=cr2*r2 + r3
           s=cr2*s2 + s3
           zbr(1,ind(ia*3+1),j,i3) = r*cr3 + r1
           zbr(2,ind(ia*3+1),j,i3) = s*cr3 + s1
           r1=r1 - r*cr3p
           s1=s1 - s*cr3p
           r2=cr2*r2-r3
           s2=cr2*s2-s3
           zbr(1,ind(ia*3+2),j,i3) = r1 - s2*bb
           zbr(2,ind(ia*3+2),j,i3) = s1 + r2*bb
           zbr(1,ind(ia*3+3),j,i3) = r1 + s2*bb
           zbr(2,ind(ia*3+3),j,i3) = s1 - r2*bb
         end do
       end do

!      Treat factor 5
     else if (now(ic)==5) then
!      radix 5
!      sin(2.d0*pi/5.d0)
       sin2=ris*0.9510565162951536d0
       ia=0

!      First step of factor 5
       do j=ma,mb
         r1=z(1,ia*5+1,j,i3)
         s1=z(2,ia*5+1,j,i3)
         r2=z(1,ia*5+2,j,i3)
         s2=z(2,ia*5+2,j,i3)
         r3=z(1,ia*5+3,j,i3)
         s3=z(2,ia*5+3,j,i3)
         r4=z(1,ia*5+4,j,i3)
         s4=z(2,ia*5+4,j,i3)
         r5=z(1,ia*5+5,j,i3)
         s5=z(2,ia*5+5,j,i3)
         r25 = r2 + r5
         r34 = r3 + r4
         s25 = s2 - s5
         s34 = s3 - s4
         zbr(1,ind(ia*5+1),j,i3) = r1 + r25 + r34
         r = r1 + cos2*r25 + cos4*r34
         s = s25 + sin42*s34
         zbr(1,ind(ia*5+2),j,i3) = r - sin2*s
         zbr(1,ind(ia*5+5),j,i3) = r + sin2*s
         r = r1 + cos4*r25 + cos2*r34
         s = sin42*s25 - s34
         zbr(1,ind(ia*5+3),j,i3) = r - sin2*s
         zbr(1,ind(ia*5+4),j,i3) = r + sin2*s
         r25 = r2 - r5
         r34 = r3 - r4
         s25 = s2 + s5
         s34 = s3 + s4
         zbr(2,ind(ia*5+1),j,i3) = s1 + s25 + s34
         r = s1 + cos2*s25 + cos4*s34
         s = r25 + sin42*r34
         zbr(2,ind(ia*5+2),j,i3) = r + sin2*s
         zbr(2,ind(ia*5+5),j,i3) = r - sin2*s
         r = s1 + cos4*s25 + cos2*s34
         s = sin42*r25 - r34
         zbr(2,ind(ia*5+3),j,i3) = r + sin2*s
         zbr(2,ind(ia*5+4),j,i3) = r - sin2*s
       end do

!      Second step of factor 5
       do ia=1,aft(ic)-1
         indx=ind(ia*5+1)-1
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         cr3=trig(1,2*indx)
         ct3=trig(2,2*indx)
         cr4=trig(1,3*indx)
         ct4=trig(2,3*indx)
         cr5=trig(1,4*indx)
         ct5=trig(2,4*indx)
         do j=ma,mb
           r1=z(1,ia*5+1,j,i3)
           s1=z(2,ia*5+1,j,i3)
           r2=cr2*(z(1,ia*5+2,j,i3) - z(2,ia*5+2,j,i3)*ct2)
           s2=cr2*(z(1,ia*5+2,j,i3)*ct2 + z(2,ia*5+2,j,i3))
           r3=cr3*(z(1,ia*5+3,j,i3) - z(2,ia*5+3,j,i3)*ct3)
           s3=cr3*(z(1,ia*5+3,j,i3)*ct3 + z(2,ia*5+3,j,i3))
           r4=z(1,ia*5+4,j,i3) - z(2,ia*5+4,j,i3)*ct4
           s4=z(1,ia*5+4,j,i3)*ct4 + z(2,ia*5+4,j,i3)
           r5=z(1,ia*5+5,j,i3) - z(2,ia*5+5,j,i3)*ct5
           s5=z(1,ia*5+5,j,i3)*ct5 + z(2,ia*5+5,j,i3)
           r25 = r2 + r5*cr5
           r34 = r3 + r4*cr4
           s25 = s2 - s5*cr5
           s34 = s3 - s4*cr4
           zbr(1,ind(ia*5+1),j,i3) = r1 + r25 + r34
           r = r1 + cos2*r25 + cos4*r34
           s = s25 + sin42*s34
           zbr(1,ind(ia*5+2),j,i3) = r - sin2*s
           zbr(1,ind(ia*5+5),j,i3) = r + sin2*s
           r = r1 + cos4*r25 + cos2*r34
           s = sin42*s25 - s34
           zbr(1,ind(ia*5+3),j,i3) = r - sin2*s
           zbr(1,ind(ia*5+4),j,i3) = r + sin2*s
           r25 = r2 - r5*cr5
           r34 = r3 - r4*cr4
           s25 = s2 + s5*cr5
           s34 = s3 + s4*cr4
           zbr(2,ind(ia*5+1),j,i3) = s1 + s25 + s34
           r = s1 + cos2*s25 + cos4*s34
           s = r25 + sin42*r34
           zbr(2,ind(ia*5+2),j,i3) = r + sin2*s
           zbr(2,ind(ia*5+5),j,i3) = r - sin2*s
           r = s1 + cos4*s25 + cos2*s34
           s = sin42*r25 - r34
           zbr(2,ind(ia*5+3),j,i3) = r + sin2*s
           zbr(2,ind(ia*5+4),j,i3) = r - sin2*s
         end do
       end do

     else
!      All factors treated
       MSG_BUG('called with factors other than 2, 3, and 5')
     end if

!    ---------------------------------------------------------------

   end do ! do i3=1,n3
 end do  ! do jj=1,n2,lot
!$OMP END PARALLEL DO

end subroutine sg_fftx
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_ffty
!! NAME
!! sg_ffty
!!
!! FUNCTION
!! This subroutine is called by the 3-dimensional fft to conduct the
!! "y" transforms for all x and z.
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  mfac = maximum number of factors in 1D FFTs
!!  mg = maximum length of 1D FFTs
!!  nd1=first dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd2=second dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd3=third dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  n1i=lower i1 index, used for blocking : the do-loop will be i1=n1i,n1
!!   put to 1 for usual ffty
!!  n1=upper i1 index, used for blocking, put usual n1 for usual ffty
!!  n3i=lower i3 index, used for blocking : the do-loop will be i3=n3i,n3
!!   put to 1 for usual ffty
!!  n3=upper i3 index, used for blocking, put usual n3 for usual ffty
!!  z(2,nd1,nd2,nd3)=INPUT array; destroyed by transformation
!!  trig, aft, now, bef, ind=provided by previous call to ctrig
!!   Note that in this routine (and in ctrig) the values in array trig are
!!   actually cos and tan, not cos and sin.  Use of tan allows advantageous
!!   use of FMA on the ibm rs6000.
!!  ris=sign of exponential in transform (should be 1 or -1; real)
!!  ic=number of (radix) factors of x transform length (from ctrig)
!!
!! OUTPUT
!!  zbr(2,nd1,nd2,nd3)=OUTPUT transformed array; no scaling applied
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,&
&          z,zbr,trig,aft,now,bef,ris,ind,ic)

 implicit none

!Arguments ------------------------------------
!Dimensions of aft, now, bef, ind, and trig should agree with
!those in subroutine ctrig.
!scalars
 integer,intent(in) :: fftcache,ic,mfac,mg,n1,n1i,n3,n3i,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 integer,intent(in) :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp),intent(in) :: trig(2,mg)
 real(dp),intent(inout) :: z(2,nd1,nd2,nd3),zbr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer :: i,ia,ib,indx,j1,j2,ntb
 real(dp),parameter :: cos2=0.3090169943749474d0   !cos(2.d0*pi/5.d0)
 real(dp),parameter :: cos4=-0.8090169943749474d0  !cos(4.d0*pi/5.d0)
 real(dp),parameter :: sin42=0.6180339887498948d0  !sin(4.d0*pi/5.d0)/sin(2.d0*pi/5.d0)
 real(dp) :: bb,cr2,cr2s,cr3,cr3p,cr4,cr5,ct2,ct3,ct4,ct5
 real(dp) :: r,r1,r2,r25,r3,r34,r4,r5,s,sin2,s1,s2,s25,s3,s34,s4,s5

! *************************************************************************

 if (fftcache<0) then
   MSG_ERROR('fftcache must be positive')
 end if

!Outer loop over z planes (j2)--note range from n3i to n3

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(aft,bef,ic,ind,n1,n1i,n3,n3i,now,ris,trig,z,zbr)
 do j2=n3i,n3

!  Direct transformation
   do i=1,ic-1
     ntb=now(i)*bef(i)

!    Treat radix 4
     if (now(i)==4) then
       ia=0

!      First step of radix 4
       do ib=1,bef(i)
!        Inner loop over all x values (j1) -- note range from n1i to n1
!        y transform is performed for this range of x values repeatedly
!        below

         do j1=n1i,n1
           r4=z(1,j1,ia*ntb+3*bef(i)+ib,j2)
           s4=z(2,j1,ia*ntb+3*bef(i)+ib,j2)
           r3=z(1,j1,ia*ntb+2*bef(i)+ib,j2)
           s3=z(2,j1,ia*ntb+2*bef(i)+ib,j2)
           r2=z(1,j1,ia*ntb+bef(i)+ib,j2)
           s2=z(2,j1,ia*ntb+bef(i)+ib,j2)
           r1=z(1,j1,ia*ntb+ib,j2)
           s1=z(2,j1,ia*ntb+ib,j2)

           r=r1 + r3
           s=r2 + r4
           z(1,j1,ia*ntb+ib,j2) = r + s
           z(1,j1,ia*ntb+2*bef(i)+ib,j2) = r - s
           r=r1 - r3
           s=s2 - s4
           z(1,j1,ia*ntb+bef(i)+ib,j2) = r - s*ris
           z(1,j1,ia*ntb+3*bef(i)+ib,j2) = r + s*ris
           r=s1 + s3
           s=s2 + s4
           z(2,j1,ia*ntb+ib,j2) = r + s
           z(2,j1,ia*ntb+2*bef(i)+ib,j2) = r - s
           r=s1 - s3
           s=r2 - r4
           z(2,j1,ia*ntb+bef(i)+ib,j2) = r + s*ris
           z(2,j1,ia*ntb+3*bef(i)+ib,j2) = r - s*ris
         end do ! j1
       end do ! ib

!      Second step of radix 4
       do ia=1,aft(i)-1
         indx=ind(ia*4*bef(i)+1)-1
         indx=indx*bef(i)
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         cr3=trig(1,2*indx)
         ct3=trig(2,2*indx)
         cr4=trig(1,3*indx)
         ct4=trig(2,3*indx)
         cr4=cr4/cr2
         cr2s=cr2*ris
         do ib=1,bef(i)
!          Range of x array again (also appears many times below)
           do j1=n1i,n1
             r4=z(1,j1,ia*ntb+3*bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+3*bef(i)+ib,j2)*ct4
             s4=z(1,j1,ia*ntb+3*bef(i)+ib,j2)*ct4 + &
&             z(2,j1,ia*ntb+3*bef(i)+ib,j2)
             r3=z(1,j1,ia*ntb+2*bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+2*bef(i)+ib,j2)*ct3
             s3=z(1,j1,ia*ntb+2*bef(i)+ib,j2)*ct3 + &
&             z(2,j1,ia*ntb+2*bef(i)+ib,j2)
             r2=z(1,j1,ia*ntb+bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)*ct2
             s2=z(1,j1,ia*ntb+bef(i)+ib,j2)*ct2 + &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)
             r1=z(1,j1,ia*ntb+ib,j2)
             s1=z(2,j1,ia*ntb+ib,j2)

             r=r1 + r3*cr3
             s=r2 + r4*cr4
             z(1,j1,ia*ntb+ib,j2) = r + s*cr2
             z(1,j1,ia*ntb+2*bef(i)+ib,j2) = r - s*cr2
             r=r1 - r3*cr3
             s=s2 - s4*cr4
             z(1,j1,ia*ntb+bef(i)+ib,j2) = r - s*cr2s
             z(1,j1,ia*ntb+3*bef(i)+ib,j2) = r + s*cr2s
             r=s1 + s3*cr3
             s=s2 + s4*cr4
             z(2,j1,ia*ntb+ib,j2) = r + s*cr2
             z(2,j1,ia*ntb+2*bef(i)+ib,j2) = r - s*cr2
             r=s1 - s3*cr3
             s=r2 - r4*cr4
             z(2,j1,ia*ntb+bef(i)+ib,j2) = r + s*cr2s
             z(2,j1,ia*ntb+3*bef(i)+ib,j2) = r - s*cr2s
           end do ! j1
         end do ! ib
       end do ! ia

!      Treat radix 2
     else if (now(i)==2) then
       ia=0

!      First step of radix 2
       do ib=1,bef(i)
         do j1=n1i,n1
           r1=z(1,j1,ia*ntb+ib,j2)
           s1=z(2,j1,ia*ntb+ib,j2)
           r2=z(1,j1,ia*ntb+bef(i)+ib,j2)
           s2=z(2,j1,ia*ntb+bef(i)+ib,j2)
           z(1,j1,ia*ntb+ib,j2) =  r2 + r1
           z(2,j1,ia*ntb+ib,j2) =  s2 + s1
           z(1,j1,ia*ntb+bef(i)+ib,j2) = -r2 + r1
           z(2,j1,ia*ntb+bef(i)+ib,j2) = -s2 + s1
         end do
       end do

!      Second step of radix 2
       do ia=1,aft(i)-1
         indx=ind(ia*2*bef(i)+1)-1
         indx=indx*bef(i)
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         do ib=1,bef(i)
           do j1=n1i,n1
             r1=z(1,j1,ia*ntb+ib,j2)
             s1=z(2,j1,ia*ntb+ib,j2)
             r2=z(1,j1,ia*ntb+bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)*ct2
             s2=z(1,j1,ia*ntb+bef(i)+ib,j2)*ct2 + &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)
             z(1,j1,ia*ntb+ib,j2) =  r2*cr2 + r1
             z(2,j1,ia*ntb+ib,j2) =  s2*cr2 + s1
             z(1,j1,ia*ntb+bef(i)+ib,j2) = -r2*cr2 + r1
             z(2,j1,ia*ntb+bef(i)+ib,j2) = -s2*cr2 + s1
           end do
         end do
       end do

!      Treat radix 3
     else if (now(i)==3) then
!      .5d0*sqrt(3.d0)=0.8660254037844387d0
       ia=0
       bb=ris*0.8660254037844387d0

!      First step of radix 3
       do ib=1,bef(i)
         do j1=n1i,n1
           r1=z(1,j1,ia*ntb+ib,j2)
           s1=z(2,j1,ia*ntb+ib,j2)
           r2=z(1,j1,ia*ntb+bef(i)+ib,j2)
           s2=z(2,j1,ia*ntb+bef(i)+ib,j2)
           r3=z(1,j1,ia*ntb+2*bef(i)+ib,j2)
           s3=z(2,j1,ia*ntb+2*bef(i)+ib,j2)
           r=r2 + r3
           s=s2 + s3
           z(1,j1,ia*ntb+ib,j2) = r + r1
           z(2,j1,ia*ntb+ib,j2) = s + s1
           r1=r1 - r*.5d0
           s1=s1 - s*.5d0
           r2=r2-r3
           s2=s2-s3
           z(1,j1,ia*ntb+bef(i)+ib,j2) = r1 - s2*bb
           z(2,j1,ia*ntb+bef(i)+ib,j2) = s1 + r2*bb
           z(1,j1,ia*ntb+2*bef(i)+ib,j2) = r1 + s2*bb
           z(2,j1,ia*ntb+2*bef(i)+ib,j2) = s1 - r2*bb
         end do
       end do

!      Second step of radix 3
       do ia=1,aft(i)-1
         indx=ind(ia*3*bef(i)+1)-1
         indx=indx*bef(i)
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         cr3=trig(1,2*indx)
         ct3=trig(2,2*indx)
         cr2=cr2/cr3
         cr3p=.5d0*cr3
         bb=ris*cr3*0.8660254037844387d0
         do ib=1,bef(i)
           do j1=n1i,n1
             r1=z(1,j1,ia*ntb+ib,j2)
             s1=z(2,j1,ia*ntb+ib,j2)
             r2=z(1,j1,ia*ntb+bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)*ct2
             s2=z(1,j1,ia*ntb+bef(i)+ib,j2)*ct2 + &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)
             r3=z(1,j1,ia*ntb+2*bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+2*bef(i)+ib,j2)*ct3
             s3=z(1,j1,ia*ntb+2*bef(i)+ib,j2)*ct3 + &
&             z(2,j1,ia*ntb+2*bef(i)+ib,j2)
             r=cr2*r2 + r3
             s=cr2*s2 + s3
             z(1,j1,ia*ntb+ib,j2) = r*cr3 + r1
             z(2,j1,ia*ntb+ib,j2) = s*cr3 + s1
             r1=r1 - r*cr3p
             s1=s1 - s*cr3p
             r2=cr2*r2-r3
             s2=cr2*s2-s3
             z(1,j1,ia*ntb+bef(i)+ib,j2) = r1 - s2*bb
             z(2,j1,ia*ntb+bef(i)+ib,j2) = s1 + r2*bb
             z(1,j1,ia*ntb+2*bef(i)+ib,j2) = r1 + s2*bb
             z(2,j1,ia*ntb+2*bef(i)+ib,j2) = s1 - r2*bb
           end do
         end do
       end do

!      Treat radix 5
     else if (now(i)==5) then
!      sin(2.d0*pi/5.d0)
       sin2=ris*0.9510565162951536d0
       ia=0

!      First step of radix 5
       do ib=1,bef(i)
         do j1=n1i,n1
           r1=z(1,j1,ia*ntb+ib,j2)
           s1=z(2,j1,ia*ntb+ib,j2)
           r2=z(1,j1,ia*ntb+bef(i)+ib,j2)
           s2=z(2,j1,ia*ntb+bef(i)+ib,j2)
           r3=z(1,j1,ia*ntb+2*bef(i)+ib,j2)
           s3=z(2,j1,ia*ntb+2*bef(i)+ib,j2)
           r4=z(1,j1,ia*ntb+3*bef(i)+ib,j2)
           s4=z(2,j1,ia*ntb+3*bef(i)+ib,j2)
           r5=z(1,j1,ia*ntb+4*bef(i)+ib,j2)
           s5=z(2,j1,ia*ntb+4*bef(i)+ib,j2)
           r25 = r2 + r5
           r34 = r3 + r4
           s25 = s2 - s5
           s34 = s3 - s4
           z(1,j1,ia*ntb+ib,j2) = r1 + r25 + r34
           r = r1 + cos2*r25 + cos4*r34
           s = s25 + sin42*s34
           z(1,j1,ia*ntb+bef(i)+ib,j2) = r - sin2*s
           z(1,j1,ia*ntb+4*bef(i)+ib,j2) = r + sin2*s
           r = r1 + cos4*r25 + cos2*r34
           s = sin42*s25 - s34
           z(1,j1,ia*ntb+2*bef(i)+ib,j2) = r - sin2*s
           z(1,j1,ia*ntb+3*bef(i)+ib,j2) = r + sin2*s
           r25 = r2 - r5
           r34 = r3 - r4
           s25 = s2 + s5
           s34 = s3 + s4
           z(2,j1,ia*ntb+ib,j2) = s1 + s25 + s34
           r = s1 + cos2*s25 + cos4*s34
           s = r25 + sin42*r34
           z(2,j1,ia*ntb+bef(i)+ib,j2) = r + sin2*s
           z(2,j1,ia*ntb+4*bef(i)+ib,j2) = r - sin2*s
           r = s1 + cos4*s25 + cos2*s34
           s = sin42*r25 - r34
           z(2,j1,ia*ntb+2*bef(i)+ib,j2) = r + sin2*s
           z(2,j1,ia*ntb+3*bef(i)+ib,j2) = r - sin2*s
         end do
       end do

!      Second step of radix 5
       do ia=1,aft(i)-1
         indx=ind(ia*5*bef(i)+1)-1
         indx=indx*bef(i)
         cr2=trig(1,indx)
         ct2=trig(2,indx)
         cr3=trig(1,2*indx)
         ct3=trig(2,2*indx)
         cr4=trig(1,3*indx)
         ct4=trig(2,3*indx)
         cr5=trig(1,4*indx)
         ct5=trig(2,4*indx)
         do ib=1,bef(i)
           do j1=n1i,n1
             r1=z(1,j1,ia*ntb+ib,j2)
             s1=z(2,j1,ia*ntb+ib,j2)
             r2=cr2*(z(1,j1,ia*ntb+bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+bef(i)+ib,j2)*ct2)
             s2=cr2*(z(1,j1,ia*ntb+bef(i)+ib,j2)*ct2 + &
&             z(2,j1,ia*ntb+bef(i)+ib,j2))
             r3=cr3*(z(1,j1,ia*ntb+2*bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+2*bef(i)+ib,j2)*ct3)
             s3=cr3*(z(1,j1,ia*ntb+2*bef(i)+ib,j2)*ct3 + &
&             z(2,j1,ia*ntb+2*bef(i)+ib,j2))
             r4=z(1,j1,ia*ntb+3*bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+3*bef(i)+ib,j2)*ct4
             s4=z(1,j1,ia*ntb+3*bef(i)+ib,j2)*ct4 + &
&             z(2,j1,ia*ntb+3*bef(i)+ib,j2)
             r5=z(1,j1,ia*ntb+4*bef(i)+ib,j2) - &
&             z(2,j1,ia*ntb+4*bef(i)+ib,j2)*ct5
             s5=z(1,j1,ia*ntb+4*bef(i)+ib,j2)*ct5 + &
&             z(2,j1,ia*ntb+4*bef(i)+ib,j2)
             r25 = r2 + r5*cr5
             r34 = r3 + r4*cr4
             s25 = s2 - s5*cr5
             s34 = s3 - s4*cr4
             z(1,j1,ia*ntb+ib,j2) = r1 + r25 + r34
             r = r1 + cos2*r25 + cos4*r34
             s = s25 + sin42*s34
             z(1,j1,ia*ntb+bef(i)+ib,j2) = r - sin2*s
             z(1,j1,ia*ntb+4*bef(i)+ib,j2) = r + sin2*s
             r = r1 + cos4*r25 + cos2*r34
             s = sin42*s25 - s34
             z(1,j1,ia*ntb+2*bef(i)+ib,j2) = r - sin2*s
             z(1,j1,ia*ntb+3*bef(i)+ib,j2) = r + sin2*s
             r25 = r2 - r5*cr5
             r34 = r3 - r4*cr4
             s25 = s2 + s5*cr5
             s34 = s3 + s4*cr4
             z(2,j1,ia*ntb+ib,j2) = s1 + s25 + s34
             r = s1 + cos2*s25 + cos4*s34
             s = r25 + sin42*r34
             z(2,j1,ia*ntb+bef(i)+ib,j2) = r + sin2*s
             z(2,j1,ia*ntb+4*bef(i)+ib,j2) = r - sin2*s
             r = s1 + cos4*s25 + cos2*s34
             s = sin42*r25 - r34
             z(2,j1,ia*ntb+2*bef(i)+ib,j2) = r + sin2*s
             z(2,j1,ia*ntb+3*bef(i)+ib,j2) = r - sin2*s
           end do
         end do
       end do

     else
!      All radices treated
       MSG_BUG('called with factors other than 2, 3, and 5')
     end if

   end do

!  ---------------------------------------------------------------

!  bitreversal

!  Treat radix 4
   if (now(ic)==4) then
     ia=0

!    First step of radix 4
     do j1=n1i,n1
       r4=z(1,j1,ia*4+4,j2)
       s4=z(2,j1,ia*4+4,j2)
       r3=z(1,j1,ia*4+3,j2)
       s3=z(2,j1,ia*4+3,j2)
       r2=z(1,j1,ia*4+2,j2)
       s2=z(2,j1,ia*4+2,j2)
       r1=z(1,j1,ia*4+1,j2)
       s1=z(2,j1,ia*4+1,j2)

       r=r1 + r3
       s=r2 + r4
       zbr(1,j1,ind(ia*4+1),j2) = r + s
       zbr(1,j1,ind(ia*4+3),j2) = r - s
       r=r1 - r3
       s=s2 - s4
       zbr(1,j1,ind(ia*4+2),j2) = r - s*ris
       zbr(1,j1,ind(ia*4+4),j2) = r + s*ris
       r=s1 + s3
       s=s2 + s4
       zbr(2,j1,ind(ia*4+1),j2) = r + s
       zbr(2,j1,ind(ia*4+3),j2) = r - s
       r=s1 - s3
       s=r2 - r4
       zbr(2,j1,ind(ia*4+2),j2) = r + s*ris
       zbr(2,j1,ind(ia*4+4),j2) = r - s*ris
     end do

!    Second step of radix 4
     do ia=1,aft(ic)-1
       indx=ind(ia*4+1)-1
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       cr3=trig(1,2*indx)
       ct3=trig(2,2*indx)
       cr4=trig(1,3*indx)
       ct4=trig(2,3*indx)
       cr4=cr4/cr2
       cr2s=cr2*ris
       do j1=n1i,n1
         r4=z(1,j1,ia*4+4,j2) - z(2,j1,ia*4+4,j2)*ct4
         s4=z(1,j1,ia*4+4,j2)*ct4 + z(2,j1,ia*4+4,j2)
         r3=z(1,j1,ia*4+3,j2) - z(2,j1,ia*4+3,j2)*ct3
         s3=z(1,j1,ia*4+3,j2)*ct3 + z(2,j1,ia*4+3,j2)
         r2=z(1,j1,ia*4+2,j2) - z(2,j1,ia*4+2,j2)*ct2
         s2=z(1,j1,ia*4+2,j2)*ct2 + z(2,j1,ia*4+2,j2)
         r1=z(1,j1,ia*4+1,j2)
         s1=z(2,j1,ia*4+1,j2)

         r=r1 + r3*cr3
         s=r2 + r4*cr4
         zbr(1,j1,ind(ia*4+1),j2) = r + s*cr2
         zbr(1,j1,ind(ia*4+3),j2) = r - s*cr2
         r=r1 - r3*cr3
         s=s2 - s4*cr4
         zbr(1,j1,ind(ia*4+2),j2) = r - s*cr2s
         zbr(1,j1,ind(ia*4+4),j2) = r + s*cr2s
         r=s1 + s3*cr3
         s=s2 + s4*cr4
         zbr(2,j1,ind(ia*4+1),j2) = r + s*cr2
         zbr(2,j1,ind(ia*4+3),j2) = r - s*cr2
         r=s1 - s3*cr3
         s=r2 - r4*cr4
         zbr(2,j1,ind(ia*4+2),j2) = r + s*cr2s
         zbr(2,j1,ind(ia*4+4),j2) = r - s*cr2s
       end do
     end do

!    Treat radix 2
   else if (now(ic)==2) then
     ia=0

!    First step of radix 2
     do j1=n1i,n1
       r1=z(1,j1,ia*2+1,j2)
       s1=z(2,j1,ia*2+1,j2)
       r2=z(1,j1,ia*2+2,j2)
       s2=z(2,j1,ia*2+2,j2)
       zbr(1,j1,ind(ia*2+1),j2) =  r2 + r1
       zbr(2,j1,ind(ia*2+1),j2) =  s2 + s1
       zbr(1,j1,ind(ia*2+2),j2) = -r2 + r1
       zbr(2,j1,ind(ia*2+2),j2) = -s2 + s1
     end do

!    Second step of radix 2
     do ia=1,aft(ic)-1
       indx=ind(ia*2+1)-1
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       do j1=n1i,n1
         r1=z(1,j1,ia*2+1,j2)
         s1=z(2,j1,ia*2+1,j2)
         r2=z(1,j1,ia*2+2,j2) - z(2,j1,ia*2+2,j2)*ct2
         s2=z(1,j1,ia*2+2,j2)*ct2 + z(2,j1,ia*2+2,j2)
         zbr(1,j1,ind(ia*2+1),j2) =  r2*cr2 + r1
         zbr(2,j1,ind(ia*2+1),j2) =  s2*cr2 + s1
         zbr(1,j1,ind(ia*2+2),j2) = -r2*cr2 + r1
         zbr(2,j1,ind(ia*2+2),j2) = -s2*cr2 + s1
       end do
     end do

!    Treat radix 3
   else if (now(ic)==3) then
!    .5d0*sqrt(3.d0)=0.8660254037844387d0
     ia=0
     bb=ris*0.8660254037844387d0

!    First step of radix 3
     do j1=n1i,n1
       r1=z(1,j1,ia*3+1,j2)
       s1=z(2,j1,ia*3+1,j2)
       r2=z(1,j1,ia*3+2,j2)
       s2=z(2,j1,ia*3+2,j2)
       r3=z(1,j1,ia*3+3,j2)
       s3=z(2,j1,ia*3+3,j2)
       r=r2 + r3
       s=s2 + s3
       zbr(1,j1,ind(ia*3+1),j2) = r + r1
       zbr(2,j1,ind(ia*3+1),j2) = s + s1
       r1=r1 - r*.5d0
       s1=s1 - s*.5d0
       r2=r2-r3
       s2=s2-s3
       zbr(1,j1,ind(ia*3+2),j2) = r1 - s2*bb
       zbr(2,j1,ind(ia*3+2),j2) = s1 + r2*bb
       zbr(1,j1,ind(ia*3+3),j2) = r1 + s2*bb
       zbr(2,j1,ind(ia*3+3),j2) = s1 - r2*bb
     end do

!    Second step of radix 3
     do ia=1,aft(ic)-1
       indx=ind(ia*3+1)-1
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       cr3=trig(1,2*indx)
       ct3=trig(2,2*indx)
       cr2=cr2/cr3
       cr3p=.5d0*cr3
       bb=ris*cr3*0.8660254037844387d0
       do j1=n1i,n1
         r1=z(1,j1,ia*3+1,j2)
         s1=z(2,j1,ia*3+1,j2)
         r2=z(1,j1,ia*3+2,j2) - z(2,j1,ia*3+2,j2)*ct2
         s2=z(1,j1,ia*3+2,j2)*ct2 + z(2,j1,ia*3+2,j2)
         r3=z(1,j1,ia*3+3,j2) - z(2,j1,ia*3+3,j2)*ct3
         s3=z(1,j1,ia*3+3,j2)*ct3 + z(2,j1,ia*3+3,j2)
         r=cr2*r2 + r3
         s=cr2*s2 + s3
         zbr(1,j1,ind(ia*3+1),j2) = r*cr3 + r1
         zbr(2,j1,ind(ia*3+1),j2) = s*cr3 + s1
         r1=r1 - r*cr3p
         s1=s1 - s*cr3p
         r2=cr2*r2-r3
         s2=cr2*s2-s3
         zbr(1,j1,ind(ia*3+2),j2) = r1 - s2*bb
         zbr(2,j1,ind(ia*3+2),j2) = s1 + r2*bb
         zbr(1,j1,ind(ia*3+3),j2) = r1 + s2*bb
         zbr(2,j1,ind(ia*3+3),j2) = s1 - r2*bb
       end do
     end do

!    Treat radix 5
   else if (now(ic)==5) then
!    sin(2.d0*pi/5.d0)
     sin2=ris*0.9510565162951536d0
     ia=0

!    First step of radix 5
     do j1=n1i,n1
       r1=z(1,j1,ia*5+1,j2)
       s1=z(2,j1,ia*5+1,j2)
       r2=z(1,j1,ia*5+2,j2)
       s2=z(2,j1,ia*5+2,j2)
       r3=z(1,j1,ia*5+3,j2)
       s3=z(2,j1,ia*5+3,j2)
       r4=z(1,j1,ia*5+4,j2)
       s4=z(2,j1,ia*5+4,j2)
       r5=z(1,j1,ia*5+5,j2)
       s5=z(2,j1,ia*5+5,j2)
       r25 = r2 + r5
       r34 = r3 + r4
       s25 = s2 - s5
       s34 = s3 - s4
       zbr(1,j1,ind(ia*5+1),j2) = r1 + r25 + r34
       r = r1 + cos2*r25 + cos4*r34
       s = s25 + sin42*s34
       zbr(1,j1,ind(ia*5+2),j2) = r - sin2*s
       zbr(1,j1,ind(ia*5+5),j2) = r + sin2*s
       r = r1 + cos4*r25 + cos2*r34
       s = sin42*s25 - s34
       zbr(1,j1,ind(ia*5+3),j2) = r - sin2*s
       zbr(1,j1,ind(ia*5+4),j2) = r + sin2*s
       r25 = r2 - r5
       r34 = r3 - r4
       s25 = s2 + s5
       s34 = s3 + s4
       zbr(2,j1,ind(ia*5+1),j2) = s1 + s25 + s34
       r = s1 + cos2*s25 + cos4*s34
       s = r25 + sin42*r34
       zbr(2,j1,ind(ia*5+2),j2) = r + sin2*s
       zbr(2,j1,ind(ia*5+5),j2) = r - sin2*s
       r = s1 + cos4*s25 + cos2*s34
       s = sin42*r25 - r34
       zbr(2,j1,ind(ia*5+3),j2) = r + sin2*s
       zbr(2,j1,ind(ia*5+4),j2) = r - sin2*s
     end do

!    Second step of radix 5
     do ia=1,aft(ic)-1
       indx=ind(ia*5+1)-1
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       cr3=trig(1,2*indx)
       ct3=trig(2,2*indx)
       cr4=trig(1,3*indx)
       ct4=trig(2,3*indx)
       cr5=trig(1,4*indx)
       ct5=trig(2,4*indx)
       do j1=n1i,n1
         r1=z(1,j1,ia*5+1,j2)
         s1=z(2,j1,ia*5+1,j2)
         r2=cr2*(z(1,j1,ia*5+2,j2) - z(2,j1,ia*5+2,j2)*ct2)
         s2=cr2*(z(1,j1,ia*5+2,j2)*ct2 + z(2,j1,ia*5+2,j2))
         r3=cr3*(z(1,j1,ia*5+3,j2) - z(2,j1,ia*5+3,j2)*ct3)
         s3=cr3*(z(1,j1,ia*5+3,j2)*ct3 + z(2,j1,ia*5+3,j2))
         r4=z(1,j1,ia*5+4,j2) - z(2,j1,ia*5+4,j2)*ct4
         s4=z(1,j1,ia*5+4,j2)*ct4 + z(2,j1,ia*5+4,j2)
         r5=z(1,j1,ia*5+5,j2) - z(2,j1,ia*5+5,j2)*ct5
         s5=z(1,j1,ia*5+5,j2)*ct5 + z(2,j1,ia*5+5,j2)
         r25 = r2 + r5*cr5
         r34 = r3 + r4*cr4
         s25 = s2 - s5*cr5
         s34 = s3 - s4*cr4
         zbr(1,j1,ind(ia*5+1),j2) = r1 + r25 + r34
         r = r1 + cos2*r25 + cos4*r34
         s = s25 + sin42*s34
         zbr(1,j1,ind(ia*5+2),j2) = r - sin2*s
         zbr(1,j1,ind(ia*5+5),j2) = r + sin2*s
         r = r1 + cos4*r25 + cos2*r34
         s = sin42*s25 - s34
         zbr(1,j1,ind(ia*5+3),j2) = r - sin2*s
         zbr(1,j1,ind(ia*5+4),j2) = r + sin2*s
         r25 = r2 - r5*cr5
         r34 = r3 - r4*cr4
         s25 = s2 + s5*cr5
         s34 = s3 + s4*cr4
         zbr(2,j1,ind(ia*5+1),j2) = s1 + s25 + s34
         r = s1 + cos2*s25 + cos4*s34
         s = r25 + sin42*r34
         zbr(2,j1,ind(ia*5+2),j2) = r + sin2*s
         zbr(2,j1,ind(ia*5+5),j2) = r - sin2*s
         r = s1 + cos4*s25 + cos2*s34
         s = sin42*r25 - r34
         zbr(2,j1,ind(ia*5+3),j2) = r + sin2*s
         zbr(2,j1,ind(ia*5+4),j2) = r - sin2*s
       end do
     end do

   else
!    All radices done
     MSG_BUG('Called with factors other than 2, 3, and 5')
   end if
 end do
!$OMP END PARALLEL DO

end subroutine sg_ffty
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fftz
!! NAME
!! sg_fftz
!!
!! FUNCTION
!! This subroutine is called by the 3-dimensional fft to conduct the
!! "z" transforms for all x and y.
!!
!! INPUTS
!!  mfac = maximum number of factors in 1D FFTs
!!  mg = maximum length of 1D FFTs
!!  nd1=first dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd2=second dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd3=third dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  n1=actual length of x and y transforms
!!  n2i=lower i2 index, used for blocking : the do-loop will be i2=n2i,n2
!!   put to 1 for usual ffty
!!  n2=upper i2 index, used for blocking, put usual n2 for usual ffty
!!  z(2,nd1,nd2,nd3)=INPUT array; destroyed by transformation
!!  trig, aft, now, bef, ind=provided by previous call to ctrig
!!   Note that in this routine (and in ctrig) the values in array trig are
!!   actually cos and tan, not cos and sin.  Use of tan allows advantageous
!!   use of FMA on the ibm rs6000.
!!  ris=sign of exponential in transform (should be 1 or -1; real)
!!  ic=number of (radix) factors of x transform length (from ctrig)
!!
!! OUTPUT
!!  zbr(2,nd1,nd2,nd3)=OUTPUT transformed array; no scaling applied
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fftz(mfac,mg,nd1,nd2,nd3,n1,n2i,n2,z,zbr,trig,aft,now,bef,ris,ind,ic)

 implicit none

!Arguments ------------------------------------
!Dimensions of aft, now, bef, ind, and trig should agree with
!those in subroutine ctrig.
!scalars
 integer,intent(in) :: ic,mfac,mg,n1,n2,n2i,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 integer,intent(in) :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp),intent(in) :: trig(2,mg)
 real(dp),intent(inout) :: z(2,nd1,nd2,nd3),zbr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer :: b_i,i,i2,ia,ib,indx,j,ntb
 real(dp),parameter :: cos2=0.3090169943749474d0   !cos(2.d0*pi/5.d0)
 real(dp),parameter :: cos4=-0.8090169943749474d0  !cos(4.d0*pi/5.d0)
 real(dp),parameter :: sin42=0.6180339887498948d0  !sin(4.d0*pi/5.d0)/sin(2.d0*pi/5.d0)
 real(dp) :: bb,cr2,cr2s,cr3,cr3p,cr4,cr5,ct2,ct3,ct4,ct5
 real(dp) :: r,r1,r2,r25,r3,r34,r4,r5,s,sin2,s1,s2,s25,s3,s34,s4,s5

! *************************************************************************

!n12 occurs as a loop index repeated below; do z transform while
!looping over all n12 lines of data

!Direct transformation (to ic-1), bitreversal will be in second part
!of routine

 do i=1,ic-1
   ntb=now(i)*bef(i)
   b_i=bef(i)

!  Treat radix 4
   if (now(i)==4) then
     ia=0

!    First step of radix 4
     do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,ia,ib,n1,n2i,n2,ntb,ris,z)
       do i2=n2i,n2
         do j=1,n1
           r4=z(1,j,i2,ia*ntb+3*b_i+ib)
           s4=z(2,j,i2,ia*ntb+3*b_i+ib)
           r3=z(1,j,i2,ia*ntb+2*b_i+ib)
           s3=z(2,j,i2,ia*ntb+2*b_i+ib)
           r2=z(1,j,i2,ia*ntb+b_i+ib)
           s2=z(2,j,i2,ia*ntb+b_i+ib)
           r1=z(1,j,i2,ia*ntb+ib)
           s1=z(2,j,i2,ia*ntb+ib)

           r=r1 + r3
           s=r2 + r4
           z(1,j,i2,ia*ntb+ib) = r + s
           z(1,j,i2,ia*ntb+2*b_i+ib) = r - s
           r=r1 - r3
           s=s2 - s4
           z(1,j,i2,ia*ntb+b_i+ib) = r - s*ris
           z(1,j,i2,ia*ntb+3*b_i+ib) = r + s*ris
           r=s1 + s3
           s=s2 + s4
           z(2,j,i2,ia*ntb+ib) = r + s
           z(2,j,i2,ia*ntb+2*b_i+ib) = r - s
           r=s1 - s3
           s=r2 - r4
           z(2,j,i2,ia*ntb+b_i+ib) = r + s*ris
           z(2,j,i2,ia*ntb+3*b_i+ib) = r - s*ris
         end do ! j
       end do ! i2
!$OMP END PARALLEL DO
     end do ! ib

!    Second step of radix 4
     do ia=1,aft(i)-1
       indx=ind(ia*4*b_i+1)-1
       indx=indx*b_i
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       cr3=trig(1,2*indx)
       ct3=trig(2,2*indx)
       cr4=trig(1,3*indx)
       ct4=trig(2,3*indx)
       cr4=cr4/cr2
       cr2s=cr2*ris
       do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,cr2,cr3,cr4,ct2,cr2s,ct3,ct4,i,ia,ib,n1,n2i,n2,ntb,ris,z)
         do i2=n2i,n2
           do j=1,n1
             r4=z(1,j,i2,ia*ntb+3*b_i+ib) - &
&             z(2,j,i2,ia*ntb+3*b_i+ib)*ct4
             s4=z(1,j,i2,ia*ntb+3*b_i+ib)*ct4 + &
&             z(2,j,i2,ia*ntb+3*b_i+ib)
             r3=z(1,j,i2,ia*ntb+2*b_i+ib) - &
&             z(2,j,i2,ia*ntb+2*b_i+ib)*ct3
             s3=z(1,j,i2,ia*ntb+2*b_i+ib)*ct3 + &
&             z(2,j,i2,ia*ntb+2*b_i+ib)
             r2=z(1,j,i2,ia*ntb+b_i+ib) - &
&             z(2,j,i2,ia*ntb+b_i+ib)*ct2
             s2=z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&             z(2,j,i2,ia*ntb+b_i+ib)
             r1=z(1,j,i2,ia*ntb+ib)
             s1=z(2,j,i2,ia*ntb+ib)

             r=r1 + r3*cr3
             s=r2 + r4*cr4
             z(1,j,i2,ia*ntb+ib) = r + s*cr2
             z(1,j,i2,ia*ntb+2*b_i+ib) = r - s*cr2
             r=r1 - r3*cr3
             s=s2 - s4*cr4
             z(1,j,i2,ia*ntb+b_i+ib) = r - s*cr2s
             z(1,j,i2,ia*ntb+3*b_i+ib) = r + s*cr2s
             r=s1 + s3*cr3
             s=s2 + s4*cr4
             z(2,j,i2,ia*ntb+ib) = r + s*cr2
             z(2,j,i2,ia*ntb+2*b_i+ib) = r - s*cr2
             r=s1 - s3*cr3
             s=r2 - r4*cr4
             z(2,j,i2,ia*ntb+b_i+ib) = r + s*cr2s
             z(2,j,i2,ia*ntb+3*b_i+ib) = r - s*cr2s
           end do ! j
         end do ! i2
!$OMP END PARALLEL DO
       end do ! ib

     end do ! ia

!    Treat radix 2
   else if (now(i)==2) then
     ia=0

!    First step of radix 2
     do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,ia,ib,n1,n2,n2i,ntb,z)
       do i2=n2i,n2
         do j=1,n1
           r1=z(1,j,i2,ia*ntb+ib)
           s1=z(2,j,i2,ia*ntb+ib)
           r2=z(1,j,i2,ia*ntb+b_i+ib)
           s2=z(2,j,i2,ia*ntb+b_i+ib)
           z(1,j,i2,ia*ntb+ib) =  r2 + r1
           z(2,j,i2,ia*ntb+ib) =  s2 + s1
           z(1,j,i2,ia*ntb+b_i+ib) = -r2 + r1
           z(2,j,i2,ia*ntb+b_i+ib) = -s2 + s1
         end do ! j
       end do ! i2
!$OMP END PARALLEL DO
     end do ! ib

!    Second step of radix 2
     do ia=1,aft(i)-1
       indx=ind(ia*2*b_i+1)-1
       indx=indx*b_i
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,cr2,ct2,ia,ib,n1,n2,n2i,ntb,z)
         do i2=n2i,n2
           do j=1,n1
             r1=z(1,j,i2,ia*ntb+ib)
             s1=z(2,j,i2,ia*ntb+ib)
             r2=z(1,j,i2,ia*ntb+b_i+ib) - &
&             z(2,j,i2,ia*ntb+b_i+ib)*ct2
             s2=z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&             z(2,j,i2,ia*ntb+b_i+ib)
             z(1,j,i2,ia*ntb+ib) =  r2*cr2 + r1
             z(2,j,i2,ia*ntb+ib) =  s2*cr2 + s1
             z(1,j,i2,ia*ntb+b_i+ib) = -r2*cr2 + r1
             z(2,j,i2,ia*ntb+b_i+ib) = -s2*cr2 + s1
           end do ! j
         end do ! i2
!$OMP END PARALLEL DO
       end do ! ib

     end do ! ia

!    Treat radix 3
   else if (now(i)==3) then
!    .5d0*sqrt(3.d0)=0.8660254037844387d0
     ia=0
     bb=ris*0.8660254037844387d0

!    First step of radix 3
     do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,b_i,ia,ib,n1,n2,n2i,ntb,z)
       do i2=n2i,n2
         do j=1,n1
           r1=z(1,j,i2,ia*ntb+ib)
           s1=z(2,j,i2,ia*ntb+ib)
           r2=z(1,j,i2,ia*ntb+b_i+ib)
           s2=z(2,j,i2,ia*ntb+b_i+ib)
           r3=z(1,j,i2,ia*ntb+2*b_i+ib)
           s3=z(2,j,i2,ia*ntb+2*b_i+ib)
           r=r2 + r3
           s=s2 + s3
           z(1,j,i2,ia*ntb+ib) = r + r1
           z(2,j,i2,ia*ntb+ib) = s + s1
           r1=r1 - r*.5d0
           s1=s1 - s*.5d0
           r2=r2-r3
           s2=s2-s3
           z(1,j,i2,ia*ntb+b_i+ib) = r1 - s2*bb
           z(2,j,i2,ia*ntb+b_i+ib) = s1 + r2*bb
           z(1,j,i2,ia*ntb+2*b_i+ib) = r1 + s2*bb
           z(2,j,i2,ia*ntb+2*b_i+ib) = s1 - r2*bb
         end do ! j
       end do ! i2
!$OMP END PARALLEL DO
     end do ! ib

!    Second step of radix 3
     do ia=1,aft(i)-1
       indx=ind(ia*3*b_i+1)-1
       indx=indx*b_i
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       cr3=trig(1,2*indx)
       ct3=trig(2,2*indx)
       cr2=cr2/cr3
       cr3p=.5d0*cr3
       bb=ris*cr3*0.8660254037844387d0
       do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,b_i,cr2,cr3,cr3p,ct2,ct3,ia,ib,n1,n2,n2i,ntb,z)
         do i2=n2i,n2
           do j=1,n1
             r1=z(1,j,i2,ia*ntb+ib)
             s1=z(2,j,i2,ia*ntb+ib)
             r2=z(1,j,i2,ia*ntb+b_i+ib) - &
&             z(2,j,i2,ia*ntb+b_i+ib)*ct2
             s2=z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&             z(2,j,i2,ia*ntb+b_i+ib)
             r3=z(1,j,i2,ia*ntb+2*b_i+ib) - &
&             z(2,j,i2,ia*ntb+2*b_i+ib)*ct3
             s3=z(1,j,i2,ia*ntb+2*b_i+ib)*ct3 + &
&             z(2,j,i2,ia*ntb+2*b_i+ib)
             r=cr2*r2 + r3
             s=cr2*s2 + s3
             z(1,j,i2,ia*ntb+ib) = r*cr3 + r1
             z(2,j,i2,ia*ntb+ib) = s*cr3 + s1
             r1=r1 - r*cr3p
             s1=s1 - s*cr3p
             r2=cr2*r2-r3
             s2=cr2*s2-s3
             z(1,j,i2,ia*ntb+b_i+ib) = r1 - s2*bb
             z(2,j,i2,ia*ntb+b_i+ib) = s1 + r2*bb
             z(1,j,i2,ia*ntb+2*b_i+ib) = r1 + s2*bb
             z(2,j,i2,ia*ntb+2*b_i+ib) = s1 - r2*bb
           end do ! j
         end do ! i2
!$OMP END PARALLEL DO
       end do ! ib

     end do ! ia

!    Treat radix 5
   else if (now(i)==5) then
     sin2=ris*0.9510565162951536d0
     ia=0

!    First step of radix 5
     do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,ia,ib,n1,n2,n2i,ntb,sin2,z)
       do i2=n2i,n2
         do j=1,n1
           r1=z(1,j,i2,ia*ntb+ib)
           s1=z(2,j,i2,ia*ntb+ib)
           r2=z(1,j,i2,ia*ntb+b_i+ib)
           s2=z(2,j,i2,ia*ntb+b_i+ib)
           r3=z(1,j,i2,ia*ntb+2*b_i+ib)
           s3=z(2,j,i2,ia*ntb+2*b_i+ib)
           r4=z(1,j,i2,ia*ntb+3*b_i+ib)
           s4=z(2,j,i2,ia*ntb+3*b_i+ib)
           r5=z(1,j,i2,ia*ntb+4*b_i+ib)
           s5=z(2,j,i2,ia*ntb+4*b_i+ib)
           r25 = r2 + r5
           r34 = r3 + r4
           s25 = s2 - s5
           s34 = s3 - s4
           z(1,j,i2,ia*ntb+ib) = r1 + r25 + r34
           r = r1 + cos2*r25 + cos4*r34
           s = s25 + sin42*s34
           z(1,j,i2,ia*ntb+b_i+ib) = r - sin2*s
           z(1,j,i2,ia*ntb+4*b_i+ib) = r + sin2*s
           r = r1 + cos4*r25 + cos2*r34
           s = sin42*s25 - s34
           z(1,j,i2,ia*ntb+2*b_i+ib) = r - sin2*s
           z(1,j,i2,ia*ntb+3*b_i+ib) = r + sin2*s
           r25 = r2 - r5
           r34 = r3 - r4
           s25 = s2 + s5
           s34 = s3 + s4
           z(2,j,i2,ia*ntb+ib) = s1 + s25 + s34
           r = s1 + cos2*s25 + cos4*s34
           s = r25 + sin42*r34
           z(2,j,i2,ia*ntb+b_i+ib) = r + sin2*s
           z(2,j,i2,ia*ntb+4*b_i+ib) = r - sin2*s
           r = s1 + cos4*s25 + cos2*s34
           s = sin42*r25 - r34
           z(2,j,i2,ia*ntb+2*b_i+ib) = r + sin2*s
           z(2,j,i2,ia*ntb+3*b_i+ib) = r - sin2*s
         end do ! j
       end do ! i2
!$OMP END PARALLEL DO
     end do ! ib

!    Second step of radix 5
     do ia=1,aft(i)-1
       indx=ind(ia*5*b_i+1)-1
       indx=indx*b_i
       cr2=trig(1,indx)
       ct2=trig(2,indx)
       cr3=trig(1,2*indx)
       ct3=trig(2,2*indx)
       cr4=trig(1,3*indx)
       ct4=trig(2,3*indx)
       cr5=trig(1,4*indx)
       ct5=trig(2,4*indx)
       do ib=1,b_i
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,cr2,cr3,cr4,cr5,ct2,ct3,ct4,ct5,ia,ib,n1,n2,n2i,ntb,sin2,z)
         do i2=n2i,n2
           do j=1,n1
             r1=z(1,j,i2,ia*ntb+ib)
             s1=z(2,j,i2,ia*ntb+ib)
             r2=cr2*(z(1,j,i2,ia*ntb+b_i+ib) - &
&             z(2,j,i2,ia*ntb+b_i+ib)*ct2)
             s2=cr2*(z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&             z(2,j,i2,ia*ntb+b_i+ib))
             r3=cr3*(z(1,j,i2,ia*ntb+2*b_i+ib) - &
&             z(2,j,i2,ia*ntb+2*b_i+ib)*ct3)
             s3=cr3*(z(1,j,i2,ia*ntb+2*b_i+ib)*ct3 + &
&             z(2,j,i2,ia*ntb+2*b_i+ib))
             r4=z(1,j,i2,ia*ntb+3*b_i+ib) - &
&             z(2,j,i2,ia*ntb+3*b_i+ib)*ct4
             s4=z(1,j,i2,ia*ntb+3*b_i+ib)*ct4 + &
&             z(2,j,i2,ia*ntb+3*b_i+ib)
             r5=z(1,j,i2,ia*ntb+4*b_i+ib) - &
&             z(2,j,i2,ia*ntb+4*b_i+ib)*ct5
             s5=z(1,j,i2,ia*ntb+4*b_i+ib)*ct5 + &
&             z(2,j,i2,ia*ntb+4*b_i+ib)
             r25 = r2 + r5*cr5
             r34 = r3 + r4*cr4
             s25 = s2 - s5*cr5
             s34 = s3 - s4*cr4
             z(1,j,i2,ia*ntb+ib) = r1 + r25 + r34
             r = r1 + cos2*r25 + cos4*r34
             s = s25 + sin42*s34
             z(1,j,i2,ia*ntb+b_i+ib) = r - sin2*s
             z(1,j,i2,ia*ntb+4*b_i+ib) = r + sin2*s
             r = r1 + cos4*r25 + cos2*r34
             s = sin42*s25 - s34
             z(1,j,i2,ia*ntb+2*b_i+ib) = r - sin2*s
             z(1,j,i2,ia*ntb+3*b_i+ib) = r + sin2*s
             r25 = r2 - r5*cr5
             r34 = r3 - r4*cr4
             s25 = s2 + s5*cr5
             s34 = s3 + s4*cr4
             z(2,j,i2,ia*ntb+ib) = s1 + s25 + s34
             r = s1 + cos2*s25 + cos4*s34
             s = r25 + sin42*r34
             z(2,j,i2,ia*ntb+b_i+ib) = r + sin2*s
             z(2,j,i2,ia*ntb+4*b_i+ib) = r - sin2*s
             r = s1 + cos4*s25 + cos2*s34
             s = sin42*r25 - r34
             z(2,j,i2,ia*ntb+2*b_i+ib) = r + sin2*s
             z(2,j,i2,ia*ntb+3*b_i+ib) = r - sin2*s
           end do ! j
         end do ! i2
!$OMP END PARALLEL DO
       end do ! ib

     end do ! ia

!    All radices treated
   else
     MSG_BUG('called with factors other than 2, 3, and 5')
   end if

!  End of direct transformation
 end do

!------------------------------------------------------------
!bitreversal  (zbr is for z"bit-reversed")

!Treat radix 4
 if (now(ic)==4) then
   ia=0

!  First step of radix 4
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,ind,n1,n2,n2i,ntb,ris,z,zbr)
   do i2=n2i,n2
     do j=1,n1
       r4=z(1,j,i2,ia*4+4)
       s4=z(2,j,i2,ia*4+4)
       r3=z(1,j,i2,ia*4+3)
       s3=z(2,j,i2,ia*4+3)
       r2=z(1,j,i2,ia*4+2)
       s2=z(2,j,i2,ia*4+2)
       r1=z(1,j,i2,ia*4+1)
       s1=z(2,j,i2,ia*4+1)

       r=r1 + r3
       s=r2 + r4
       zbr(1,j,i2,ind(ia*4+1)) = r + s
       zbr(1,j,i2,ind(ia*4+3)) = r - s
       r=r1 - r3
       s=s2 - s4
       zbr(1,j,i2,ind(ia*4+2)) = r - s*ris
       zbr(1,j,i2,ind(ia*4+4)) = r + s*ris
       r=s1 + s3
       s=s2 + s4
       zbr(2,j,i2,ind(ia*4+1)) = r + s
       zbr(2,j,i2,ind(ia*4+3)) = r - s
       r=s1 - s3
       s=r2 - r4
       zbr(2,j,i2,ind(ia*4+2)) = r + s*ris
       zbr(2,j,i2,ind(ia*4+4)) = r - s*ris
     end do ! j
   end do ! i2
!$OMP END PARALLEL DO

!  Second step of radix 4
   do ia=1,aft(ic)-1
     indx=ind(ia*4+1)-1
     cr2=trig(1,indx)
     ct2=trig(2,indx)
     cr3=trig(1,2*indx)
     ct3=trig(2,2*indx)
     cr4=trig(1,3*indx)
     ct4=trig(2,3*indx)
     cr4=cr4/cr2
     cr2s=cr2*ris
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,cr2,cr2s,cr3,cr4,ct2,ct3,ct4,ind,n1,n2,n2i,z,zbr)
     do i2=n2i,n2
       do j=1,n1
         r4=z(1,j,i2,ia*4+4) - z(2,j,i2,ia*4+4)*ct4
         s4=z(1,j,i2,ia*4+4)*ct4 + z(2,j,i2,ia*4+4)
         r3=z(1,j,i2,ia*4+3) - z(2,j,i2,ia*4+3)*ct3
         s3=z(1,j,i2,ia*4+3)*ct3 + z(2,j,i2,ia*4+3)
         r2=z(1,j,i2,ia*4+2) - z(2,j,i2,ia*4+2)*ct2
         s2=z(1,j,i2,ia*4+2)*ct2 + z(2,j,i2,ia*4+2)
         r1=z(1,j,i2,ia*4+1)
         s1=z(2,j,i2,ia*4+1)

         r=r1 + r3*cr3
         s=r2 + r4*cr4
         zbr(1,j,i2,ind(ia*4+1)) = r + s*cr2
         zbr(1,j,i2,ind(ia*4+3)) = r - s*cr2
         r=r1 - r3*cr3
         s=s2 - s4*cr4
         zbr(1,j,i2,ind(ia*4+2)) = r - s*cr2s
         zbr(1,j,i2,ind(ia*4+4)) = r + s*cr2s
         r=s1 + s3*cr3
         s=s2 + s4*cr4
         zbr(2,j,i2,ind(ia*4+1)) = r + s*cr2
         zbr(2,j,i2,ind(ia*4+3)) = r - s*cr2
         r=s1 - s3*cr3
         s=r2 - r4*cr4
         zbr(2,j,i2,ind(ia*4+2)) = r + s*cr2s
         zbr(2,j,i2,ind(ia*4+4)) = r - s*cr2s
       end do ! j
     end do ! i2
!$OMP END PARALLEL DO

   end do ! ia

!  Treat radix 2
 else if (now(ic)==2) then
   ia=0

!  First step of radix 2
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,ind,n1,n2,n2i,z,zbr)
   do i2=n2i,n2
     do j=1,n1
       r1=z(1,j,i2,ia*2+1)
       s1=z(2,j,i2,ia*2+1)
       r2=z(1,j,i2,ia*2+2)
       s2=z(2,j,i2,ia*2+2)
       zbr(1,j,i2,ind(ia*2+1)) =  r2 + r1
       zbr(2,j,i2,ind(ia*2+1)) =  s2 + s1
       zbr(1,j,i2,ind(ia*2+2)) = -r2 + r1
       zbr(2,j,i2,ind(ia*2+2)) = -s2 + s1
     end do ! j
   end do ! i2
!$OMP END PARALLEL DO

!  Second step of radix 2
   do ia=1,aft(ic)-1
     indx=ind(ia*2+1)-1
     cr2=trig(1,indx)
     ct2=trig(2,indx)
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(cr2,ct2,ia,ind,n1,n2,n2i,z,zbr)
     do i2=n2i,n2
       do j=1,n1
         r1=z(1,j,i2,ia*2+1)
         s1=z(2,j,i2,ia*2+1)
         r2=z(1,j,i2,ia*2+2) - z(2,j,i2,ia*2+2)*ct2
         s2=z(1,j,i2,ia*2+2)*ct2 + z(2,j,i2,ia*2+2)
         zbr(1,j,i2,ind(ia*2+1)) =  r2*cr2 + r1
         zbr(2,j,i2,ind(ia*2+1)) =  s2*cr2 + s1
         zbr(1,j,i2,ind(ia*2+2)) = -r2*cr2 + r1
         zbr(2,j,i2,ind(ia*2+2)) = -s2*cr2 + s1
       end do ! j
     end do ! i2
!$OMP END PARALLEL DO
   end do ! ia

!  Treat radix 3
 else if (now(ic)==3) then
!  .5d0*sqrt(3.d0)=0.8660254037844387d0
   ia=0
   bb=ris*0.8660254037844387d0

!  First step of radix 3
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,ia,ind,n1,n2,n2i,z,zbr)
   do i2=n2i,n2
     do j=1,n1
       r1=z(1,j,i2,ia*3+1)
       s1=z(2,j,i2,ia*3+1)
       r2=z(1,j,i2,ia*3+2)
       s2=z(2,j,i2,ia*3+2)
       r3=z(1,j,i2,ia*3+3)
       s3=z(2,j,i2,ia*3+3)
       r=r2 + r3
       s=s2 + s3
       zbr(1,j,i2,ind(ia*3+1)) = r + r1
       zbr(2,j,i2,ind(ia*3+1)) = s + s1
       r1=r1 - r*.5d0
       s1=s1 - s*.5d0
       r2=r2-r3
       s2=s2-s3
       zbr(1,j,i2,ind(ia*3+2)) = r1 - s2*bb
       zbr(2,j,i2,ind(ia*3+2)) = s1 + r2*bb
       zbr(1,j,i2,ind(ia*3+3)) = r1 + s2*bb
       zbr(2,j,i2,ind(ia*3+3)) = s1 - r2*bb
     end do ! j
   end do ! i2
!$OMP END PARALLEL DO

!  Second step of radix 3
   do ia=1,aft(ic)-1
     indx=ind(ia*3+1)-1
     cr2=trig(1,indx)
     ct2=trig(2,indx)
     cr3=trig(1,2*indx)
     ct3=trig(2,2*indx)
     cr2=cr2/cr3
     cr3p=.5d0*cr3
     bb=ris*cr3*0.8660254037844387d0
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,cr2,cr3,cr3p,ct2,ct3,ia,ind,n1,n2,n2i,z,zbr)
     do i2=n2i,n2
       do j=1,n1
         r1=z(1,j,i2,ia*3+1)
         s1=z(2,j,i2,ia*3+1)
         r2=z(1,j,i2,ia*3+2) - z(2,j,i2,ia*3+2)*ct2
         s2=z(1,j,i2,ia*3+2)*ct2 + z(2,j,i2,ia*3+2)
         r3=z(1,j,i2,ia*3+3) - z(2,j,i2,ia*3+3)*ct3
         s3=z(1,j,i2,ia*3+3)*ct3 + z(2,j,i2,ia*3+3)
         r=cr2*r2 + r3
         s=cr2*s2 + s3
         zbr(1,j,i2,ind(ia*3+1)) = r*cr3 + r1
         zbr(2,j,i2,ind(ia*3+1)) = s*cr3 + s1
         r1=r1 - r*cr3p
         s1=s1 - s*cr3p
         r2=cr2*r2-r3
         s2=cr2*s2-s3
         zbr(1,j,i2,ind(ia*3+2)) = r1 - s2*bb
         zbr(2,j,i2,ind(ia*3+2)) = s1 + r2*bb
         zbr(1,j,i2,ind(ia*3+3)) = r1 + s2*bb
         zbr(2,j,i2,ind(ia*3+3)) = s1 - r2*bb
       end do ! j
     end do ! i2
!$OMP END PARALLEL DO
   end do ! ia

!  Treat radix 5
 else if (now(ic)==5) then
!  sin(2.d0*pi/5.d0)
   sin2=ris*0.9510565162951536d0
   ia=0

!  First step of radix 5
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,ind,n1,n2,n2i,sin2,z,zbr)
   do i2=n2i,n2
     do j=1,n1
       r1=z(1,j,i2,ia*5+1)
       s1=z(2,j,i2,ia*5+1)
       r2=z(1,j,i2,ia*5+2)
       s2=z(2,j,i2,ia*5+2)
       r3=z(1,j,i2,ia*5+3)
       s3=z(2,j,i2,ia*5+3)
       r4=z(1,j,i2,ia*5+4)
       s4=z(2,j,i2,ia*5+4)
       r5=z(1,j,i2,ia*5+5)
       s5=z(2,j,i2,ia*5+5)
       r25 = r2 + r5
       r34 = r3 + r4
       s25 = s2 - s5
       s34 = s3 - s4
       zbr(1,j,i2,ind(ia*5+1)) = r1 + r25 + r34
       r = r1 + cos2*r25 + cos4*r34
       s = s25 + sin42*s34
       zbr(1,j,i2,ind(ia*5+2)) = r - sin2*s
       zbr(1,j,i2,ind(ia*5+5)) = r + sin2*s
       r = r1 + cos4*r25 + cos2*r34
       s = sin42*s25 - s34
       zbr(1,j,i2,ind(ia*5+3)) = r - sin2*s
       zbr(1,j,i2,ind(ia*5+4)) = r + sin2*s
       r25 = r2 - r5
       r34 = r3 - r4
       s25 = s2 + s5
       s34 = s3 + s4
       zbr(2,j,i2,ind(ia*5+1)) = s1 + s25 + s34
       r = s1 + cos2*s25 + cos4*s34
       s = r25 + sin42*r34
       zbr(2,j,i2,ind(ia*5+2)) = r + sin2*s
       zbr(2,j,i2,ind(ia*5+5)) = r - sin2*s
       r = s1 + cos4*s25 + cos2*s34
       s = sin42*r25 - r34
       zbr(2,j,i2,ind(ia*5+3)) = r + sin2*s
       zbr(2,j,i2,ind(ia*5+4)) = r - sin2*s
     end do ! j
   end do ! i2
!$OMP END PARALLEL DO

!  Second step of radix 5
   do ia=1,aft(ic)-1
     indx=ind(ia*5+1)-1
     cr2=trig(1,indx)
     ct2=trig(2,indx)
     cr3=trig(1,2*indx)
     ct3=trig(2,2*indx)
     cr4=trig(1,3*indx)
     ct4=trig(2,3*indx)
     cr5=trig(1,4*indx)
     ct5=trig(2,4*indx)
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(cr2,cr3,cr4,cr5,ct2,ct3,ct4,ct5,ia,ind,n1,n2,n2i,sin2,z,zbr)
     do i2=n2i,n2
       do j=1,n1
         r1=z(1,j,i2,ia*5+1)
         s1=z(2,j,i2,ia*5+1)
         r2=cr2*(z(1,j,i2,ia*5+2) - z(2,j,i2,ia*5+2)*ct2)
         s2=cr2*(z(1,j,i2,ia*5+2)*ct2 + z(2,j,i2,ia*5+2))
         r3=cr3*(z(1,j,i2,ia*5+3) - z(2,j,i2,ia*5+3)*ct3)
         s3=cr3*(z(1,j,i2,ia*5+3)*ct3 + z(2,j,i2,ia*5+3))
         r4=z(1,j,i2,ia*5+4) - z(2,j,i2,ia*5+4)*ct4
         s4=z(1,j,i2,ia*5+4)*ct4 + z(2,j,i2,ia*5+4)
         r5=z(1,j,i2,ia*5+5) - z(2,j,i2,ia*5+5)*ct5
         s5=z(1,j,i2,ia*5+5)*ct5 + z(2,j,i2,ia*5+5)
         r25 = r2 + r5*cr5
         r34 = r3 + r4*cr4
         s25 = s2 - s5*cr5
         s34 = s3 - s4*cr4
         zbr(1,j,i2,ind(ia*5+1)) = r1 + r25 + r34
         r = r1 + cos2*r25 + cos4*r34
         s = s25 + sin42*s34
         zbr(1,j,i2,ind(ia*5+2)) = r - sin2*s
         zbr(1,j,i2,ind(ia*5+5)) = r + sin2*s
         r = r1 + cos4*r25 + cos2*r34
         s = sin42*s25 - s34
         zbr(1,j,i2,ind(ia*5+3)) = r - sin2*s
         zbr(1,j,i2,ind(ia*5+4)) = r + sin2*s
         r25 = r2 - r5*cr5
         r34 = r3 - r4*cr4
         s25 = s2 + s5*cr5
         s34 = s3 + s4*cr4
         zbr(2,j,i2,ind(ia*5+1)) = s1 + s25 + s34
         r = s1 + cos2*s25 + cos4*s34
         s = r25 + sin42*r34
         zbr(2,j,i2,ind(ia*5+2)) = r + sin2*s
         zbr(2,j,i2,ind(ia*5+5)) = r - sin2*s
         r = s1 + cos4*s25 + cos2*s34
         s = sin42*r25 - r34
         zbr(2,j,i2,ind(ia*5+3)) = r + sin2*s
         zbr(2,j,i2,ind(ia*5+4)) = r - sin2*s
       end do ! j
     end do ! i2
!$OMP END PARALLEL DO
   end do ! ia

 else !  All radices treated
   MSG_BUG('called with factors other than 2, 3, and 5')
 end if

end subroutine sg_fftz
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_ctrig
!! NAME
!! sg_ctrig
!!
!! FUNCTION
!! Precalculates trigonometric expressions and bitreversal key IND (Stefan Goedecker lib).
!!
!! INPUTS
!! n=Number of FFT points for 1D FFT.
!! ris  = sign of exponential in transform (should be 1 or -1; real)
!! mfac = maximum number of factors in 1D FFTs
!! mg   = maximum length of 1D FFTs
!!
!! OUTPUT
!! trig(2,mg) TO BE DESCRIBED SB 090902
!! aft(mfac) TO BE DESCRIBED SB 090902
!! bef(mfac) TO BE DESCRIBED SB 090902
!! now(mfac) TO BE DESCRIBED SB 090902
!! ic = number of (radix) factors of x transform length (from ctrig)
!! ind(mg) TO BE DESCRIBED SB 090902
!!
!! NOTES
!! * This version of sg_ctrig produces cos and tan instead of sin and cos--
!!   this allows for much greater efficiency on the superscalar architecture
!!   of ibm rs6000 where floating point multiply and add (FMA) is used.
!!
!! * This routine is not thread-safe due to the presence of variables with the save attribute!
!!   DO NOT CALL THIS ROUTINE INSIDE A OPENMP PARALLEL REGION
!!
!! TODO
!! Should describe arguments
!! Should suppress one-letter variables
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_ctrig(n,trig,aft,bef,now,ris,ic,ind,mfac,mg)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mfac,mg,n
 integer,intent(out) :: ic
 real(dp),intent(in) :: ris
!arrays
 integer,intent(out) :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp),intent(out) :: trig(2,mg)

!Local variables-------------------------------
!scalars
 integer,save :: nextmx=4
 integer :: i,ii,inc,irep,j,k,l,next,nh
 integer,save :: prime(4)=(/5,4,3,2/)  !"prime" is the set of radices coded elsewhere for fft
 real(dp) :: angle,trigc,trigs,twopi
 character(len=500) :: message

! *************************************************************************

!**Note**
!2*Pi must not be defined too accurately here or else
!cos(twopi/2) will be exactly 0 and sin/cos below will be
!infinite; if a small error is left in Pi, then sin/cos will
!be about 10**14 and later cos * (sin/cos) will be 1 to within
!about 10**(-14) and the fft routines will work
!The precision on sgi causes the algorithm to fail if
!twopi is defined as 8.d0*atan(1.0d0).

 twopi=6.2831853071795867d0

 angle=ris*twopi/n
!trig(1,0)=1.d0
!trig(2,0)=0.d0
 if (mod(n,2)==0) then
   nh=n/2
   trig(1,nh)=-1.d0
   trig(2,nh)=0.d0
   do i=1,nh-1
     trigc=cos(i*angle)
     trigs=sin(i*angle)
     trig(1,i)=trigc
     trig(2,i)=trigs/trigc
     trig(1,n-i)=trigc
     trig(2,n-i)=-trigs/trigc
   end do
 else
   nh=(n-1)/2
   do i=1,nh
     trigc=cos(i*angle)
     trigs=sin(i*angle)
     trig(1,i)=trigc
     trig(2,i)=trigs/trigc
     trig(1,n-i)=trigc
     trig(2,n-i)=-trigs/trigc
   end do
 end if

 ic=1
 aft(ic)=1
 bef(ic)=n
 next=1

!An infinite loop, with exit or cycle instructions
 do
   if( (bef(ic)/prime(next))*prime(next)<bef(ic) ) then
     next=next+1
     if (next<=nextmx) then
       cycle
     else
       now(ic)=bef(ic)
       bef(ic)=1
     end if
   else
     now(ic)=prime(next)
     bef(ic)=bef(ic)/prime(next)
   end if
   aft(ic+1)=aft(ic)
   now(ic+1)=now(ic)
   bef(ic+1)=bef(ic)
   ic=ic+1
   if (ic>mfac) then
     write(message, '(a,i0,2a,i0)' )&
&     'number of factors ic=',ic,ch10,&
&     'exceeds dimensioned mfac=',mfac
     MSG_BUG(message)
   end if
   if (bef(ic)/=1) then
     aft(ic)=aft(ic)*now(ic)
     cycle
   end if
!  If not cycled, exit
   exit
 end do

 ic=ic-1

!DEBUG
!write(std_out,*) 'now',(now(i),i=1,ic)
!write(std_out,*) 'aft',(aft(i),i=1,ic)
!write(std_out,*) 'bef',(bef(i),i=1,ic)
!ENDDEBUG

 do i=1,n
   ind(i)=1
 end do

 irep=1
 inc=n
 do l=ic,1,-1
   inc=inc/now(l)
   ii=0
   do k=1,1+(n-1)/(now(l)*irep)
     do j=0,now(l)-1
       do i=1,irep
         ii=ii+1
         ind(ii)=ind(ii)+j*inc
       end do
     end do
   end do
   irep=irep*now(l)
 end do

 if (irep/=n) then
   write(message,'(a,i0,a,i0)')'  irep should equal n ; irep=',irep,' n=',n
   MSG_BUG(message)
 end if

 if (inc/=1) then
   write(message, '(a,i0)' )' inc should equal 1 in sg_ctrig; inc=',inc
   MSG_BUG(message)
 end if

end subroutine sg_ctrig
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fftrisc
!! NAME
!! sg_fftrisc
!!
!! FUNCTION
!!  Wrapper around fftrisc_one_nothreadsafe that supports ndat transforms.
!!
!! * This routine is not thread-safe due to the presence of variables with the save attribute!
!!   DO NOT CALL THIS ROUTINE INSIDE A OPENMP PARALLEL REGION
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fftrisc(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
& kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,n4,n5,n6,option,weight_r, weight_i)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,n4,n5,n6,ndat,npwin,npwout,option
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(in) :: fofgin(2,npwin*ndat)
 real(dp),intent(inout) :: denpot(cplex*n4*n5*n6),fofr(2,n4*n5*n6*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer :: idat,fofgin_p,fofr_p,fofgout_p
!arrays
 real(dp) :: dum_fofgin(0,0),dum_fofr(0,0),dum_fofgout(0,0)

! *************************************************************************

 do idat=1,ndat
   fofgin_p = 1 + (idat-1) * npwin
   fofr_p = 1 + (idat - 1) * n4*n5*n6
   fofgout_p = 1 + (idat-1) * npwout

   select case (option)
   case (0)
     call fftrisc_one_nothreadsafe(&
&      cplex,denpot,fofgin(1,fofgin_p),dum_fofgout,fofr(1,fofr_p),&
&      gboundin,gboundout,istwf_k,&
&      kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)

   case (1)
     ! Don't know why but fofr is not touched by this option.
     call fftrisc_one_nothreadsafe(&
&      cplex,denpot,fofgin(1,fofgin_p),dum_fofgout,dum_fofr,&
&      gboundin,gboundout,istwf_k,&
&      kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)

   case (2)
     call fftrisc_one_nothreadsafe(&
&      cplex,denpot,fofgin(1,fofgin_p),fofgout(1,fofgout_p),dum_fofr,&
&      gboundin,gboundout,istwf_k,&
&      kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)

   case (3)
     call fftrisc_one_nothreadsafe(&
&      cplex,denpot,dum_fofgin,fofgout(1,fofgout_p),fofr(1,fofr_p),&
&      gboundin,gboundout,istwf_k,&
&      kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)

   case default
      MSG_ERROR("Wrong option")
   end select
 end do

end subroutine sg_fftrisc
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/fftrisc_one_nothreadsafe
!! NAME
!! fftrisc_one_nothreadsafe
!!
!! FUNCTION
!! Carry out Fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!!
!! NOTES
!! Specifically uses rather sophisticated algorithms, based on S Goedecker
!! routines, specialized for superscalar RISC architecture.
!! Zero padding : saves 7/12 execution time
!! Bi-dimensional data locality in most of the routine : cache reuse
!! For k-point (0 0 0) : takes advantage of symmetry of data.
!! Note however that no blocking is used, in both 1D z-transform
!! or subsequent 2D transform. This should be improved.
!!
!! * This routine is not thread-safe due to the presence of variables with the save attribute!
!!   DO NOT CALL THIS ROUTINE INSIDE A OPENMP PARALLEL REGION
!!
!! INPUTS
!!  cplex= if 1 , denpot is real, if 2 , denpot is complex
!!     (cplex=2 only allowed for option=2 when istwf_k=1)
!!     one can also use cplex=0 if option=0 or option=3
!!  fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!  gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!!  gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_kin(3,npwin)=reduced planewave coordinates, input
!!  kg_kout(3,npwout)=reduced planewave coordinates, output
!!  mgfft=maximum size of 1D FFTs
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwin=number of elements in fofgin array (for option 0, 1 and 2)
!!  npwout=number of elements in fofgout array (for option 2 and 3)
!!  n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!!  option= if 0: do direct FFT
!!          if 1: do direct FFT, then sum the density
!!          if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!          if 3: do reverse FFT only
!!  weight=weight to be used for the accumulation of the density in real space
!!          (needed only when option=1)
!!
!! OUTPUT
!!  (see side effects)
!!
!! OPTIONS
!!  The different options are:
!!  - reciprocal to real space and output the result (when option=0),
!!  - reciprocal to real space and accumulate the density (when option=1) or
!!  - reciprocal to real space, apply the local potential to the wavefunction
!!    in real space and produce the result in reciprocal space (when option=2)
!!  - real space to reciprocal space (when option=3).
!!  option=0 IS NOT ALLOWED when istwf_k>2
!!  option=3 IS NOT ALLOWED when istwf_k>=2
!!
!! SIDE EFFECTS
!!  for option==0, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 fofr(2,n4,n5,n6) contains the Fourier Transform of fofgin;
!!                 no use of denpot, fofgout and npwout.
!!  for option==1, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*n4,n5,n6) contains the input density at input,
!!                 and the updated density at output;
!!                 no use of fofgout and npwout.
!!  for option==2, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*n4,n5,n6) contains the input local potential;
!!                 fofgout(2,npwout) contains the output function;
!!  for option==3, fofr(2,n4,n5,n6) contains the real space wavefunction;
!!                 fofgout(2,npwout) contains its Fourier transform;
!!                 no use of fofgin and npwin.
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine fftrisc_one_nothreadsafe(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
& kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,n4,n5,n6,npwin,npwout,option
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofr(2,n4,n5,n6)
 real(dp),intent(out) :: fofgout(2,npwout)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer,save :: ic1,ic2,ic3,ic4,ic5,ic6,n1_save=0,n2_save=0,n3_save=0
 integer :: fftcache,g2max,g2min,i1,i1max,i2,i3,i3inv,ig,igb
 integer :: igb_inv,igbmax,ii2,lot,lotin,lotout,mgb,n1
 integer :: n1half1,n1halfm,n1i,n2,n2half1,n3,n4half1,n5half1,nfftot,ngbin
 integer :: ngbout,nlot,nproc_omp
 real(dp) :: ai,ar,fraction,norm,phai,phar,wkim,wkre
 character(len=500) :: message
!arrays
 integer,save :: aft1(mfac),aft2(mfac),aft3(mfac),aft4(mfac),aft5(mfac)
 integer,save :: aft6(mfac),bef1(mfac),bef2(mfac),bef3(mfac),bef4(mfac)
 integer,save :: bef5(mfac),bef6(mfac),ind1(mg),ind2(mg),ind3(mg),ind4(mg)
 integer,save :: ind5(mg),ind6(mg),now1(mfac),now2(mfac),now3(mfac),now4(mfac)
 integer,save :: now5(mfac),now6(mfac)
 integer :: gbound_dum(4)
 integer,allocatable :: indpw_kin(:,:),indpw_kout(:,:)
 real(dp),save :: trig1(2,mg),trig2(2,mg),trig3(2,mg),trig4(2,mg),trig5(2,mg)
 real(dp),save :: trig6(2,mg)
 real(dp),allocatable :: pha1(:,:),pha2(:,:),pha3(:,:),wk1d_a(:,:,:,:)
 real(dp),allocatable :: wk1d_b(:,:,:,:),wk2d_a(:,:,:,:),wk2d_b(:,:,:,:)
 real(dp),allocatable :: wk2d_c(:,:,:,:),wk2d_d(:,:,:,:)
#if defined HAVE_OPENMP
 integer,external :: OMP_GET_NUM_THREADS
#endif

! *************************************************************************

 if(istwf_k>2 .and. option==0)then
   write(message,'(a,i0)')' option=0 is not allowed with istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(istwf_k>=2 .and. option==3)then
   write(message,'(a,i0)')' option=3 is not allowed with istwf_k=',istwf_k
   MSG_BUG(message)
 end if

!For all other tests of validity of inputs, assume that they
!have been done in the calling routine

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3) ; nfftot=n1*n2*n3
 fftcache=ngfft(8)

 if(option/=3)then
   ABI_ALLOCATE(indpw_kin,(4,npwin))
   call indfftrisc(gboundin(3:3+2*mgfft+4,1),indpw_kin,kg_kin,mgfft,ngbin,ngfft,npwin)
 end if
 if(option==2 .or. option==3)then
   ABI_ALLOCATE(indpw_kout,(4,npwout))
   call indfftrisc(gboundout(3:3+2*mgfft+4,1),indpw_kout,kg_kout,mgfft,ngbout,ngfft,npwout)
 end if

!Define the dimension of the first work arrays, for 1D transforms along z ,
!taking into account the need to avoid the cache trashing
 if(option==2)then
   mgb=max(ngbin,ngbout)
 else if(option==0 .or. option==1)then
   mgb=ngbin ; ngbout=1
 else if(option==3)then
   mgb=ngbout ; ngbin=1
 end if

 if(mod(mgb,2)/=1)mgb=mgb+1

!Initialise openmp, if needed
!$OMP PARALLEL
!$OMP SINGLE
 nproc_omp=1
#if defined HAVE_OPENMP
 nproc_omp=OMP_GET_NUM_THREADS()
#endif
!$OMP END SINGLE
!$OMP END PARALLEL

!For the treatment of the z transform,
!one tries to use only a fraction of the cache, since the
!treatment of the array wk1d_a will not involve contiguous segments
 fraction=0.25
!First estimation of lot and nlot
 lot=(fftcache*fraction*1000)/(n3*8*2)+1
!Select the smallest integer multiple of nproc_omp, larger
!or equal to nlot. In this way, the cache size is not exhausted,
!and one takes care correctly of the number of processors.
!Treat separately the in and out cases
 nlot=(ngbin-1)/lot+1
 nlot=nproc_omp*((nlot-1)/nproc_omp+1)
 lotin=(ngbin-1)/nlot+1
 nlot=(ngbout-1)/lot+1
 nlot=nproc_omp*((nlot-1)/nproc_omp+1)
 lotout=(ngbout-1)/nlot+1
!The next line impose only one lot. Usually, comment it.
!lotin=mgb ; lotout=mgb

!Compute auxiliary arrays needed for FFTs
 if(n1/=n1_save)then
   call sg_ctrig(n1,trig1,aft1,bef1,now1,one,ic1,ind1,mfac,mg)
   call sg_ctrig(n1,trig4,aft4,bef4,now4,-one,ic4,ind4,mfac,mg)
   n1_save=n1
 end if
 if(n2/=n2_save)then
   call sg_ctrig(n2,trig2,aft2,bef2,now2,one,ic2,ind2,mfac,mg)
   call sg_ctrig(n2,trig5,aft5,bef5,now5,-one,ic5,ind5,mfac,mg)
   n2_save=n2
 end if
 if(n3/=n3_save)then
   call sg_ctrig(n3,trig3,aft3,bef3,now3,one,ic3,ind3,mfac,mg)
   call sg_ctrig(n3,trig6,aft6,bef6,now6,-one,ic6,ind6,mfac,mg)
   n3_save=n3
 end if

!------------------------------------------------------------------
!Here, call general k-point code

 if(istwf_k==1)then

!  Note that the z transform will appear as a y transform
   ABI_ALLOCATE(wk1d_a,(2,mgb,n3,1))
   ABI_ALLOCATE(wk1d_b,(2,mgb,n3,1))

   if(option/=3)then

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(n3,ngbin,wk1d_a)
     do i3=1,n3
       do igb=1,ngbin
         wk1d_a(1,igb,i3,1)=zero
         wk1d_a(2,igb,i3,1)=zero
       end do
     end do
!$OMP END PARALLEL DO

!    Insert fofgin into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(fofgin,indpw_kin,npwin,wk1d_a)
     do ig=1,npwin
       igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
       wk1d_a(1,igb,i3,1)=fofgin(1,ig)
       wk1d_a(2,igb,i3,1)=fofgin(2,ig)
     end do
!$OMP END PARALLEL DO

!    Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbin,lotin
       igbmax=min(igb+lotin-1,ngbin)
!      Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a,wk1d_b, &
&       trig3,aft3,now3,bef3,one,ind3,ic3)
     end do
!$OMP END PARALLEL DO

   end if !  if(option/=3)

!  Do-loop on the planes stacked in the z direction
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP&SHARED(aft1,aft2,aft4,aft5,bef1,bef2,bef4,bef5,cplex,denpot) &
!$OMP&SHARED(fftcache,fofr,gboundin,gboundout)&
!$OMP&SHARED(ic1,ic2,ic4,ic5,ind1,ind2,ind4) &
!$OMP&SHARED(ind5,indpw_kin,indpw_kout,mgb,n1,n2,n3,n4,n5,ngbin) &
!$OMP&SHARED(ngbout,now1,now2,now4,now5,option,trig1,trig2,trig4,trig5) &
!$OMP&SHARED(weight_r,weight_i,wk1d_a,wk1d_b)

!  Allocate two 2-dimensional work arrays
   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
!$OMP DO
   do i3=1,n3

     if(option/=3)then
!      Zero the values on the current plane
!      wk2d_a(1:2,1:n1,1:n2,1)=zero
       do i2=1,n2
         do i1=1,n1
           wk2d_a(1,i1,i2,1)=zero
           wk2d_a(2,i1,i2,1)=zero
         end do
       end do
!      Copy the data in the current plane
       do igb=1,ngbin
         i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
         wk2d_a(1,i1,i2,1)=wk1d_b(1,igb,i3,1)
         wk2d_a(2,i1,i2,1)=wk1d_b(2,igb,i3,1)
       end do
!      Perform x transform, taking into account arrays of zeros
       g2min=gboundin(3,1) ; g2max=gboundin(4,1)
       if ( g2min+n2 >= g2max+2 ) then
         do i2=g2max+2,g2min+n2
           do i1=1,n1
             wk2d_b(1,i1,i2,1)=zero
             wk2d_b(2,i1,i2,1)=zero
           end do
         end do
       end if
       gbound_dum(1)=1 ; gbound_dum(2)=1
       gbound_dum(3)=g2min ; gbound_dum(4)=g2max
       call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_a,wk2d_b,&
&       trig1,aft1,now1,bef1,one,ind1,ic1,gbound_dum)
!      Perform y transform
       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_b,wk2d_a, &
&       trig2,aft2,now2,bef2,one,ind2,ic2)
!      The wave function is now in real space, for the current plane
     end if

     if(option==0)then ! Copy the transformed function at the right place
       do i2=1,n2
         do i1=1,n1
           fofr(1,i1,i2,i3)=wk2d_a(1,i1,i2,1)
           fofr(2,i1,i2,i3)=wk2d_a(2,i1,i2,1)
         end do
       end do
     end if

     if(option==1)then ! Accumulate density
       do i2=1,n2
         do i1=1,n1
           denpot(i1,i2,i3)=denpot(i1,i2,i3)+weight_r*wk2d_a(1,i1,i2,1)**2+weight_i*wk2d_a(2,i1,i2,1)**2
         end do
       end do
     end if

     if(option==2)then ! Apply local potential
       if(cplex==1)then
         do i2=1,n2
           do i1=1,n1
             wk2d_a(1,i1,i2,1)=denpot(i1,i2,i3)*wk2d_a(1,i1,i2,1)
             wk2d_a(2,i1,i2,1)=denpot(i1,i2,i3)*wk2d_a(2,i1,i2,1)
           end do
         end do
       else
         do i2=1,n2
           do i1=1,n1
             wkre=wk2d_a(1,i1,i2,1)
             wkim=wk2d_a(2,i1,i2,1)
             wk2d_a(1,i1,i2,1)=denpot(2*i1-1,i2,i3)*wkre -denpot(2*i1  ,i2,i3)*wkim
             wk2d_a(2,i1,i2,1)=denpot(2*i1-1,i2,i3)*wkim +denpot(2*i1  ,i2,i3)*wkre
           end do
         end do
       end if
     end if

     if(option==3)then ! Copy the function to be tranformed at the right place
       do i2=1,n2
         do i1=1,n1
           wk2d_a(1,i1,i2,1)=fofr(1,i1,i2,i3)
           wk2d_a(2,i1,i2,1)=fofr(2,i1,i2,i3)
         end do
       end do
     end if

     if(option==2 .or. option==3)then  ! Perform y transform
       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_a,wk2d_b, &
&       trig5,aft5,now5,bef5,-one,ind5,ic5)
!      Perform x transform, taking into account arrays of zeros
       gbound_dum(1)=1 ; gbound_dum(2)=1
       gbound_dum(3)=gboundout(3,1) ; gbound_dum(4)=gboundout(4,1)
       call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_b,wk2d_a,&
&       trig4,aft4,now4,bef4,-one,ind4,ic4,gbound_dum)
!      Copy the data from the current plane to wk1d_b
       do igb=1,ngbout
         i1=indpw_kout(1,igb) ; i2=indpw_kout(2,igb)
         wk1d_b(1,igb,i3,1)=wk2d_a(1,i1,i2,1)
         wk1d_b(2,igb,i3,1)=wk2d_a(2,i1,i2,1)
       end do
     end if

!    End loop on planes
   end do
!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
!$OMP END PARALLEL

   if(option==2 .or. option==3)then

!    Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft6,bef6,fftcache,ind6,ic6,lotout,mgb)&
!$OMP&SHARED(ngbout,now6,n3,trig6,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbout,lotout
       igbmax=min(igb+lotout-1,ngbout)
!      Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_b,wk1d_a, &
&       trig6,aft6,now6,bef6,-one,ind6,ic6)

     end do
!$OMP END PARALLEL DO

!    Transfer the data in the output array, after normalization
     norm=1.d0/dble(nfftot)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(fofgout,indpw_kout,norm,npwout,wk1d_a)
     do ig=1,npwout
       igb=indpw_kout(4,ig) ; i3=indpw_kout(3,ig)
       fofgout(1,ig)=wk1d_a(1,igb,i3,1)*norm
       fofgout(2,ig)=wk1d_a(2,igb,i3,1)*norm
     end do
!$OMP END PARALLEL DO
   end if

   ABI_DEALLOCATE(wk1d_a)
   ABI_DEALLOCATE(wk1d_b)

!  End general k-point part
 end if

!------------------------------------------------------------------
!Here, use of time-reversal symmetry

 if(istwf_k>=2)then

   n1half1=n1/2+1 ; n1halfm=(n1+1)/2
   n2half1=n2/2+1
!  n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
   n4half1=(n1half1/2)*2+1
   n5half1=(n2half1/2)*2+1
!  Note that the z transform will appear as a y transform
   ABI_ALLOCATE(wk1d_a,(2,mgb,n3,1))
   ABI_ALLOCATE(wk1d_b,(2,mgb,n3,1))

   if(istwf_k/=2)then
     ABI_ALLOCATE(pha1,(2,n1))
     ABI_ALLOCATE(pha2,(2,n2))
     ABI_ALLOCATE(pha3,(3,n3))
     do i1=1,n1
       pha1(1,i1)=cos(dble(i1-1)*pi/dble(n1))
       pha1(2,i1)=sin(dble(i1-1)*pi/dble(n1))
     end do
     do i2=1,n2
       pha2(1,i2)=cos(dble(i2-1)*pi/dble(n2))
       pha2(2,i2)=sin(dble(i2-1)*pi/dble(n2))
     end do
     do i3=1,n3
       pha3(1,i3)=cos(dble(i3-1)*pi/dble(n3))
       pha3(2,i3)=sin(dble(i3-1)*pi/dble(n3))
     end do
   end if

   if(option/=3)then

!    Zero the components of wk1d_a
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(n3,ngbin,wk1d_a)
     do i3=1,n3
       do igb=1,ngbin
         wk1d_a(1,igb,i3,1)=zero
         wk1d_a(2,igb,i3,1)=zero
       end do
     end do
!$OMP END PARALLEL DO

!    Insert fofgin into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(fofgin,indpw_kin,npwin,wk1d_a)
     do ig=1,npwin
       igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
       wk1d_a(1,igb,i3,1)=fofgin(1,ig)
       wk1d_a(2,igb,i3,1)=fofgin(2,ig)
     end do
!$OMP END PARALLEL DO

!    Must complete the i2=1 plane when $k_y \equiv 0$

!    Take care of i1=1 when $k_x \equiv 0$
     if(istwf_k==2)then
!      Take care of i1=1
       do i3=n3/2+1,n3
         i3inv=n3+2-i3
         wk1d_a(1,1,i3,1)= wk1d_a(1,1,i3inv,1)
         wk1d_a(2,1,i3,1)=-wk1d_a(2,1,i3inv,1)
       end do
     else if(istwf_k==4)then
!      Take care of i1=1
       do i3=n3/2+1,n3
         i3inv=n3+1-i3
         wk1d_a(1,1,i3,1)= wk1d_a(1,1,i3inv,1)
         wk1d_a(2,1,i3,1)=-wk1d_a(2,1,i3inv,1)
       end do
     end if

!    Now, take care of other i1 values, except i3==1 when $k_z \equiv 0$
     i1max=gboundin(6,1)+1
     if(istwf_k==2)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(i1max,n3,wk1d_a)
       do igb=2,2*i1max-1
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+2-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     else if(istwf_k==3)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(i1max,n3,wk1d_a)
       do igb=1,2*i1max
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+2-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     else if(istwf_k==4)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(i1max,n3,wk1d_a)
       do igb=2,2*i1max-1
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+1-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     else if(istwf_k==5)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(i1max,n3,wk1d_a)
       do igb=1,2*i1max
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+1-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     end if

!    Now, i3==1
     if(istwf_k==2)then
       do igb=2,i1max
         igb_inv=2*i1max+1-igb
         wk1d_a(1,igb_inv,1,1)= wk1d_a(1,igb,1,1)
         wk1d_a(2,igb_inv,1,1)=-wk1d_a(2,igb,1,1)
       end do
     else if(istwf_k==3)then
       do igb=1,i1max
         igb_inv=2*i1max+1-igb
         wk1d_a(1,igb_inv,1,1)= wk1d_a(1,igb,1,1)
         wk1d_a(2,igb_inv,1,1)=-wk1d_a(2,igb,1,1)
       end do
     end if

!    Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbin,lotin
       igbmax=min(igb+lotin-1,ngbin)
!      Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a,wk1d_b, &
&       trig3,aft3,now3,bef3,one,ind3,ic3)
     end do
!$OMP END PARALLEL DO

!    Change the phase if $k_z \neq 0$
     if(istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ngbin,n3,pha3,wk1d_b)
       do i3=1,n3
         phar=pha3(1,i3)
         phai=pha3(2,i3)
         do igb=1,ngbin
           ar=wk1d_b(1,igb,i3,1)
           ai=wk1d_b(2,igb,i3,1)
           wk1d_b(1,igb,i3,1)=phar*ar-phai*ai
           wk1d_b(2,igb,i3,1)=phai*ar+phar*ai
         end do
       end do
!$OMP END PARALLEL DO
     end if

   end if !  if(option/=3)

!  Do-loop on the planes stacked in the z direction

!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP&SHARED(aft1,aft2,aft4,aft5,bef1,bef2,bef4,bef5,denpot) &
!$OMP&SHARED(fftcache,fofr,gboundin,ic1,ic2,ic4,ic5,ind1,ind2,ind4,ind5) &
!$OMP&SHARED(indpw_kin,indpw_kout,istwf_k,mgb,n1,n1half1) &
!$OMP&SHARED(n1halfm,n2,n2half1,n3,n4,n5,ngbin,ngbout) &
!$OMP&SHARED(now1,now2,now4,now5,option,pha1,pha2,trig1) &
!$OMP&SHARED(trig2,trig4,trig5,weight_r,weight_i,wk1d_a,wk1d_b)

!  Allocate two 2-dimensional work arrays
   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_c,(2,2*n1halfm,n5,1))
   ABI_ALLOCATE(wk2d_d,(2,2*n1halfm,n5,1))
!$OMP DO
   do i3=1,n3

     g2max=gboundin(4,1)

     if(option/=3)then
!      Zero the values on the current plane : need only from i2=1 to g2max+1
       do i2=1,g2max+1
         do i1=1,n1
           wk2d_a(1,i1,i2,1)=zero
           wk2d_a(2,i1,i2,1)=zero
         end do
       end do

!      Copy the data in the current plane
       do igb=1,ngbin
         i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
         wk2d_a(1,i1,i2,1)=wk1d_b(1,igb,i3,1)
         wk2d_a(2,i1,i2,1)=wk1d_b(2,igb,i3,1)
       end do

!      Perform x transform, taking into account arrays of zeros
       call sg_fftx(fftcache,mfac,mg,n4,n5,1,g2max+1,1,wk2d_a,wk2d_b,&
&       trig1,aft1,now1,bef1,one,ind1,ic1)

!      Change the phase if $k_x \neq 0$
       if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9)then
         do i1=1,n1
           phar=pha1(1,i1)
           phai=pha1(2,i1)
           do i2=1,g2max+1
             ar=wk2d_b(1,i1,i2,1)
             ai=wk2d_b(2,i1,i2,1)
             wk2d_b(1,i1,i2,1)=phar*ar-phai*ai
             wk2d_b(2,i1,i2,1)=phai*ar+phar*ai
           end do
         end do
       end if

!      Compute symmetric and antisymmetric combinations
       if(istwf_k>=2 .and. istwf_k<=5)then
         do i1=1,n1half1-1
           wk2d_a(1,i1,1,1)=wk2d_b(1,2*i1-1,1,1)
           wk2d_a(2,i1,1,1)=wk2d_b(1,2*i1  ,1,1)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           wk2d_a(1,n1half1,1,1)=wk2d_b(1,n1,1,1)
           wk2d_a(2,n1half1,1,1)=zero
         end if
         ii2=2
       else
         ii2=1
       end if
       if( g2max+1 >= ii2)then
         do i2=ii2,g2max+1
           do i1=1,n1half1-1
             wk2d_a(1,i1,i2,1)=        wk2d_b(1,2*i1-1,i2,1)-wk2d_b(2,2*i1,i2,1)
             wk2d_a(2,i1,i2,1)=        wk2d_b(2,2*i1-1,i2,1)+wk2d_b(1,2*i1,i2,1)
             wk2d_a(1,i1,n2+ii2-i2,1)= wk2d_b(1,2*i1-1,i2,1)+wk2d_b(2,2*i1,i2,1)
             wk2d_a(2,i1,n2+ii2-i2,1)=-wk2d_b(2,2*i1-1,i2,1)+wk2d_b(1,2*i1,i2,1)
           end do
           if((2*n1half1-2)/=n1)then
             wk2d_a(1,n1half1,i2,1)=        wk2d_b(1,n1,i2,1)
             wk2d_a(2,n1half1,i2,1)=        wk2d_b(2,n1,i2,1)
             wk2d_a(1,n1half1,n2+ii2-i2,1)= wk2d_b(1,n1,i2,1)
             wk2d_a(2,n1half1,n2+ii2-i2,1)=-wk2d_b(2,n1,i2,1)
           end if
         end do
       end if
       if ( n2half1 >= g2max+2 ) then
         do i2=g2max+2,n2half1
           do i1=1,n1half1-1
             wk2d_a(1,i1,i2,1)=zero
             wk2d_a(2,i1,i2,1)=zero
             wk2d_a(1,i1,n2+ii2-i2,1)=zero
             wk2d_a(2,i1,n2+ii2-i2,1)=zero
           end do
           if((2*n1half1-2)/=n1)then
             wk2d_a(1,n1half1,i2,1)=zero
             wk2d_a(2,n1half1,i2,1)=zero
             wk2d_a(1,n1half1,n2+ii2-i2,1)=zero
             wk2d_a(2,n1half1,n2+ii2-i2,1)=zero
           end if
         end do
       end if

       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1halfm,1,1,wk2d_a,wk2d_b,&
&       trig2,aft2,now2,bef2,one,ind2,ic2)

!      Change the phase if $k_y \neq 0$
       if(istwf_k>=6 .and. istwf_k<=9)then
         do i2=1,n2
           phar=pha2(1,i2)
           phai=pha2(2,i2)
           do i1=1,n1halfm
             ar=wk2d_b(1,i1,i2,1)
             ai=wk2d_b(2,i1,i2,1)
             wk2d_b(1,i1,i2,1)= phar*ar-phai*ai
             wk2d_b(2,i1,i2,1)= phai*ar+phar*ai
           end do
         end do
       end if

     end if ! option/=3

!    The wave function is now in real space, for the current plane,
!    represented by REAL numbers, although packed in the complex array wk2d_b

     g2max=gboundin(4,1)

     if(option==0)then
!      This option is only permitted for istwf_k==2 (Gamma point)
!      Copy the transformed function at the right place
       do i2=1,n2
         do i1=1,n1half1-1
           fofr(1,2*i1-1,i2,i3)=wk2d_b(1,i1,i2,1)
           fofr(1,2*i1  ,i2,i3)=wk2d_b(2,i1,i2,1)
           fofr(2,2*i1-1,i2,i3)=zero
           fofr(2,2*i1  ,i2,i3)=zero
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           fofr(1,n1,i2,i3)=wk2d_b(1,n1half1,i2,1)
           fofr(2,n1,i2,i3)=zero
         end if
       end do
     end if

     if(option==1)then ! Accumulate density
       do i2=1,n2
         do i1=1,n1half1-1
           denpot(2*i1-1,i2,i3)=denpot(2*i1-1,i2,i3)+weight_r*wk2d_b(1,i1,i2,1)**2
           denpot(2*i1  ,i2,i3)=denpot(2*i1  ,i2,i3)+weight_i*wk2d_b(2,i1,i2,1)**2
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           denpot(n1,i2,i3)=denpot(n1,i2,i3)+weight_r*wk2d_b(1,n1half1,i2,1)**2
         end if
       end do
     end if

     if(option==2)then ! Apply local potential
       do i2=1,n2
         do i1=1,n1half1-1
           wk2d_a(1,i1,i2,1)=denpot(2*i1-1,i2,i3)*wk2d_b(1,i1,i2,1)
           wk2d_a(2,i1,i2,1)=denpot(2*i1  ,i2,i3)*wk2d_b(2,i1,i2,1)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           wk2d_a(1,n1half1,i2,1)=denpot(n1,i2,i3)*wk2d_b(1,n1half1,i2,1)
           wk2d_a(2,n1half1,i2,1)=zero
         end if
       end do
     end if

     if(option==3)then
!      This option is only permitted for istwf_k==2 (Gamma point)
!      Copy the transformed function at the right place
       do i2=1,n2
         do i1=1,n1half1-1
           wk2d_b(1,i1,i2,1)=fofr(1,2*i1-1,i2,i3)
           wk2d_b(2,i1,i2,1)=fofr(1,2*i1  ,i2,i3)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           wk2d_b(1,n1half1,i2,1)=fofr(1,n1,i2,i3)
         end if
       end do
     end if

     if(option==2 .or. option==3)then  ! Change the phase if $k_y \neq 0$
       if(istwf_k>=6 .and. istwf_k<=9)then
         do i2=1,n2
           phar=pha2(1,i2)
           phai=pha2(2,i2)
           do i1=1,n1halfm
             ar=wk2d_a(1,i1,i2,1)
             ai=wk2d_a(2,i1,i2,1)
             wk2d_a(1,i1,i2,1)= phar*ar+phai*ai
             wk2d_a(2,i1,i2,1)=-phai*ar+phar*ai
           end do
         end do
       end if

!      Perform y transform
       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1halfm,1,1,wk2d_a,wk2d_b, &
&       trig5,aft5,now5,bef5,-one,ind5,ic5)

!      Decompose symmetric and antisymmetric parts
       if(istwf_k>=2 .and. istwf_k<=5)then
         do i1=1,n1halfm
           wk2d_c(1,2*i1-1,1,1)=wk2d_b(1,i1,1,1)
           wk2d_c(2,2*i1-1,1,1)=zero
           wk2d_c(1,2*i1,1,1)=wk2d_b(2,i1,1,1)
           wk2d_c(2,2*i1,1,1)=zero
         end do
         ii2=2
       else
         ii2=1
       end if
       do i2=ii2,g2max+1
         do i1=1,n1halfm
           wk2d_c(1,2*i1-1,i2,1)=(wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+ii2-i2,1))*0.5d0
           wk2d_c(2,2*i1-1,i2,1)=(wk2d_b(2,i1,i2,1)-wk2d_b(2,i1,n2+ii2-i2,1))*0.5d0
           wk2d_c(1,2*i1,i2,1)= ( wk2d_b(2,i1,i2,1)+wk2d_b(2,i1,n2+ii2-i2,1))*0.5d0
           wk2d_c(2,2*i1,i2,1)= (-wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+ii2-i2,1))*0.5d0
         end do
       end do

!      Change the phase if $k_x \neq 0$
       if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9 )then
         do i1=1,n1
           phar=pha1(1,i1)
           phai=pha1(2,i1)
           do i2=1,g2max+1
             ar=wk2d_c(1,i1,i2,1)
             ai=wk2d_c(2,i1,i2,1)
             wk2d_c(1,i1,i2,1)= phar*ar+phai*ai
             wk2d_c(2,i1,i2,1)=-phai*ar+phar*ai
           end do
         end do
       end if

!      Perform x transform : for y=1 to g2max+1, to benefit from zeros
       call sg_fftx(fftcache,mfac,mg,2*n1halfm,n5,1,g2max+1,1,wk2d_c,wk2d_d,&
&       trig4,aft4,now4,bef4,-one,ind4,ic4)

!      Copy the data from the current plane to wk1d_b
       do igb=1,ngbout
         i1=indpw_kout(1,igb) ; i2=indpw_kout(2,igb)
         wk1d_b(1,igb,i3,1)=wk2d_d(1,i1,i2,1)
         wk1d_b(2,igb,i3,1)=wk2d_d(2,i1,i2,1)
       end do

     end if ! option==2 or 3

!    End loop on planes
   end do

!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_c)
   ABI_DEALLOCATE(wk2d_d)
!$OMP END PARALLEL

   if(option==2 .or. option==3)then

!    Change the phase if $k_z \neq 0$
     if(istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ngbout,n3,pha3,wk1d_b)
       do i3=1,n3
         phar=pha3(1,i3)
         phai=pha3(2,i3)
         do igb=1,ngbout
           ar=wk1d_b(1,igb,i3,1)
           ai=wk1d_b(2,igb,i3,1)
           wk1d_b(1,igb,i3,1)= phar*ar+phai*ai
           wk1d_b(2,igb,i3,1)=-phai*ar+phar*ai
         end do
       end do
!$OMP END PARALLEL DO
     end if

!    Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft6,bef6,fftcache,ind6,ic6,lotout,mgb)&
!$OMP&SHARED(ngbout,now6,n3,trig6,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbout,lotout
       igbmax=min(igb+lotout-1,ngbout)
!      Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_b,wk1d_a, &
&       trig6,aft6,now6,bef6,-one,ind6,ic6)

     end do
!$OMP END PARALLEL DO

!    Transfer the data in the output array, after normalization
     norm=1.d0/dble(nfftot)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(fofgout,indpw_kout,norm,npwout,wk1d_a)
     do ig=1,npwout
       igb=indpw_kout(4,ig) ; i3=indpw_kout(3,ig)
       fofgout(1,ig)=wk1d_a(1,igb,i3,1)*norm
       fofgout(2,ig)=wk1d_a(2,igb,i3,1)*norm
     end do
!$OMP END PARALLEL DO

   end if

   ABI_DEALLOCATE(wk1d_a)
   ABI_DEALLOCATE(wk1d_b)

   if(istwf_k/=2)then
     ABI_DEALLOCATE(pha1)
     ABI_DEALLOCATE(pha2)
     ABI_DEALLOCATE(pha3)
   end if

 end if !  End time-reversal symmetry

!------------------------------------------------------------------

 if(option/=3) then
   ABI_DEALLOCATE(indpw_kin)
 end if
 if(option==2 .or. option==3) then
   ABI_DEALLOCATE(indpw_kout)
 end if

end subroutine fftrisc_one_nothreadsafe
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_fftrisc_2
!! NAME
!! sg_fftrisc_2
!!
!! FUNCTION
!! Carry out Fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!! if luse_ndo is activated, do two FFT, and compute the density with
!! non-diagonal occupations.
!!
!! NOTES
!! * Specifically uses rather sophisticated algorithms, based on S Goedecker
!!   routines, specialized for superscalar RISC architecture.
!!   Zero padding : saves 7/12 execution time
!!   Bi-dimensional data locality in most of the routine : cache reuse
!!   For k-point (0 0 0) : takes advantage of symmetry of data.
!!   Note however that no blocking is used, in both 1D z-transform
!!   or subsequent 2D transform. This should be improved.
!!
!! * This routine is not thread-safe due to the presence of variables with the save attribute!
!!   DO NOT CALL THIS ROUTINE INSIDE A OPENMP PARALLEL REGION
!!
!! INPUTS
!!  cplex= if 1 , denpot is real, if 2 , denpot is complex
!!     (cplex=2 only allowed for option=2 when istwf_k=1)
!!     one can also use cplex=0 if option=0 or option=3
!!  fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!  fofgin_p(2,npwin) (optional) =holds second input wavefunction in G vector basis sphere.
!!  gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!!  gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_kin(3,npwin)=reduced planewave coordinates, input
!!  kg_kout(3,npwout)=reduced planewave coordinates, output
!!  luse_ndo (optional) = use non diagonal occup (in this case, exists fofgin_p)
!!  npwin=number of elements in fofgin array (for option 0, 1 and 2)
!!  npwout=number of elements in fofgout array (for option 2 and 3)
!!  mgfft=maximum size of 1D FFTs
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!!  option= if 0: do direct FFT
!!          if 1: do direct FFT, then sum the density
!!          if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!          if 3: do reverse FFT only
!!  weight=weight to be used for the accumulation of the density in real space
!!          (needed only when option=1)
!!
!! OUTPUT
!!  (see side effects)
!!
!! OPTIONS
!!  The different options are:
!!  - reciprocal to real space and output the result (when option=0),
!!  - reciprocal to real space and accumulate the density (when option=1) or
!!  - reciprocal to real space, apply the local potential to the wavefunction
!!    in real space and produce the result in reciprocal space (when option=2)
!!  - real space to reciprocal space (when option=3).
!!  option=0 IS NOT ALLOWED when istwf_k>2
!!  option=3 IS NOT ALLOWED when istwf_k>=2
!!  (this version can be used to compute fft of two wavefunction and
!!     compute the product in denpot)
!!
!! SIDE EFFECTS
!!  for option==0, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 fofr(2,n4,n5,n6) contains the Fourier Transform of fofgin;
!!                 no use of denpot, fofgout and npwout.
!!  for option==1, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*n4,n5,n6) contains the input density at input,
!!                 and the updated density at output;
!!                 fofr(2,n4,n5,n6) contains the Fourier transform of fofgin,
!!                 except in the case of the hp library subroutine;
!!                 no use of fofgout and npwout.
!!  for option==2, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*n4,n5,n6) contains the input local potential;
!!                 fofgout(2,npwout) contains the output function;
!!                 fofr(2,n4,n5,n6) contains the Fourier transform of fofgin,
!!                 except in the case of the hp library subroutine.
!!  for option==3, fofr(2,n4,n5,n6) contains the real space wavefunction;
!!                 fofgout(2,npwout) contains its Fourier transform;
!!                 no use of fofgin and npwin.
!!
!! TODO
!! Complete input and output list.
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_fftrisc_2(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
& istwf_k,kg_kin,kg_kout,&
& mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_2,&
& luse_ndo,fofgin_p) ! optional

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,n4,n5,n6,npwin,npwout,option
 real(dp),intent(in) :: weight_r
 real(dp),intent(in),optional :: weight_2
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 logical,intent(in),optional :: luse_ndo
 real(dp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(in),optional :: fofgin_p(:,:)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofr(2,n4,n5,n6)
 real(dp),intent(out) :: fofgout(2,npwout)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer,save :: ic1,ic2,ic3,ic4,ic5,ic6,n1_save=0,n2_save=0,n3_save=0
 integer :: fftcache,g2max,g2min,i1,i1max,i2,i3,i3inv,ig,igb
 integer :: igb_inv,igbmax,ii2,lot,lotin,lotout,mgb,n1
 integer :: n1half1,n1halfm,n1i,n2,n2half1,n3,n4half1,n5half1,nfftot,ngbin
 integer :: ngbout,nlot,nproc_omp
 integer :: weight_i
 real(dp) :: ai,ar,fraction,norm,phai,phar,wkim,wkre
 character(len=500) :: message
!arrays
 integer,save :: aft1(mfac),aft2(mfac),aft3(mfac),aft4(mfac),aft5(mfac)
 integer,save :: aft6(mfac),bef1(mfac),bef2(mfac),bef3(mfac),bef4(mfac)
 integer,save :: bef5(mfac),bef6(mfac),ind1(mg),ind2(mg),ind3(mg),ind4(mg)
 integer,save :: ind5(mg),ind6(mg),now1(mfac),now2(mfac),now3(mfac),now4(mfac)
 integer,save :: now5(mfac),now6(mfac)
 integer :: gbound_dum(4)
 integer,allocatable :: indpw_kin(:,:),indpw_kout(:,:)
 logical :: lluse_ndo
 real(dp),save :: trig1(2,mg),trig2(2,mg),trig3(2,mg),trig4(2,mg),trig5(2,mg)
 real(dp),save :: trig6(2,mg)
 real(dp),allocatable :: pha1(:,:),pha2(:,:),pha3(:,:),wk1d_a(:,:,:,:)
 real(dp),allocatable :: wk1d_b(:,:,:,:),wk2d_a(:,:,:,:),wk2d_b(:,:,:,:)
 real(dp),allocatable :: wk2d_c(:,:,:,:),wk2d_d(:,:,:,:)
 real(dp),allocatable :: wk1d_a_p(:,:,:,:),wk1d_b_p(:,:,:,:)
 real(dp),allocatable :: wk2d_a_p(:,:,:,:),wk2d_b_p(:,:,:,:)
#if defined HAVE_OPENMP
 integer,external :: OMP_GET_NUM_THREADS
#endif

! *************************************************************************

 !DBG_ENTER("COLL")

!DEBUG
!write(std_out,*)' sg_fftrisc_2 : enter, istwf_k= ',istwf_k
!write(std_out,*)' sg_fftrisc_2 : option,mgfft=',option,mgfft
!write(std_out,*)' sg_fftrisc_2 : gboundin(3:2*mgfft+6,1)='
!do ii=1,mgfft+2
!write(std_out,*)gboundin(2*ii+1,1),gboundin(2*ii+2,1)
!end do
!stop
!ENDDEBUG
!
 lluse_ndo=.true.
 if(istwf_k/=1)then
   write(message,'(a,i0)' )' It is not yet allowed to use dmft with istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(istwf_k>2 .and. option==0)then
   write(message, '(a,i0)' )' It is not allowed to use option=0 with istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(istwf_k>=2 .and. option==3)then
   write(message, '(a,i0)' )'  It is not allowed to use option=3 with istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 lluse_ndo=.false.
 if(present(luse_ndo).and.present(fofgin_p)) then
   if(luse_ndo) lluse_ndo=.true.
 end if
 if(lluse_ndo) then
   if((size(fofgin_p,2)==0).and.(luse_ndo)) then
     write(message, '(a,a,a,i4,i5)' )&
&     'fofgin_p has a dimension equal to zero and luse_ndo true',ch10,&
&     'Action: check dimension of fofgin_p',size(fofgin_p,2),luse_ndo
     MSG_BUG(message)
   end if
 end if

 weight_i= weight_r
 if ( present (weight_2 )) then
     weight_i= weight_2
     if ( present(luse_ndo) .and. (luse_ndo) )weight_i=weight_r
 end if

!For all other tests of validity of inputs, assume that they
!have been done in the calling routine

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3) ; nfftot=n1*n2*n3
 fftcache=ngfft(8)

 if(option/=3)then
   ABI_ALLOCATE(indpw_kin,(4,npwin))
   call indfftrisc(gboundin(3:3+2*mgfft+4,1),indpw_kin,kg_kin,mgfft,ngbin,ngfft,npwin)
 end if
 if(option==2 .or. option==3)then
   ABI_ALLOCATE(indpw_kout,(4,npwout))
   call indfftrisc(gboundout(3:3+2*mgfft+4,1),indpw_kout,kg_kout,mgfft,ngbout,ngfft,npwout)
 end if

!Define the dimension of the first work arrays, for 1D transforms along z ,
!taking into account the need to avoid the cache trashing
 if(option==2)then
   mgb=max(ngbin,ngbout)
 else if(option==0 .or. option==1)then
   mgb=ngbin ; ngbout=1
 else if(option==3)then
   mgb=ngbout ; ngbin=1
 end if

 if(mod(mgb,2)/=1)mgb=mgb+1

!Initialise openmp, if needed
!$OMP PARALLEL
!$OMP SINGLE
 nproc_omp=1
#if defined HAVE_OPENMP
 nproc_omp=OMP_GET_NUM_THREADS()
#endif
!$OMP END SINGLE
!$OMP END PARALLEL

!For the treatment of the z transform,
!one tries to use only a fraction of the cache, since the
!treatment of the array wk1d_a will not involve contiguous segments
 fraction=0.25
!First estimation of lot and nlot
 lot=(fftcache*fraction*1000)/(n3*8*2)+1
!Select the smallest integer multiple of nproc_omp, larger
!or equal to nlot. In this way, the cache size is not exhausted,
!and one takes care correctly of the number of processors.
!Treat separately the in and out cases
 nlot=(ngbin-1)/lot+1
 nlot=nproc_omp*((nlot-1)/nproc_omp+1)
 lotin=(ngbin-1)/nlot+1
 nlot=(ngbout-1)/lot+1
 nlot=nproc_omp*((nlot-1)/nproc_omp+1)
 lotout=(ngbout-1)/nlot+1
!The next line impose only one lot. Usually, comment it.
!lotin=mgb ; lotout=mgb

!Compute auxiliary arrays needed for FFTs
 if(n1/=n1_save)then
   call sg_ctrig(n1,trig1,aft1,bef1,now1,one,ic1,ind1,mfac,mg)
   call sg_ctrig(n1,trig4,aft4,bef4,now4,-one,ic4,ind4,mfac,mg)
   n1_save=n1
 end if
 if(n2/=n2_save)then
   call sg_ctrig(n2,trig2,aft2,bef2,now2,one,ic2,ind2,mfac,mg)
   call sg_ctrig(n2,trig5,aft5,bef5,now5,-one,ic5,ind5,mfac,mg)
   n2_save=n2
 end if
 if(n3/=n3_save)then
   call sg_ctrig(n3,trig3,aft3,bef3,now3,one,ic3,ind3,mfac,mg)
   call sg_ctrig(n3,trig6,aft6,bef6,now6,-one,ic6,ind6,mfac,mg)
   n3_save=n3
 end if

!------------------------------------------------------------------
!Here, call general k-point code

 if(istwf_k==1)then

!  Note that the z transform will appear as a y transform
   ABI_ALLOCATE(wk1d_a,(2,mgb,n3,1))
   ABI_ALLOCATE(wk1d_b,(2,mgb,n3,1))
   ABI_ALLOCATE(wk1d_a_p,(2,mgb,n3,1))
   ABI_ALLOCATE(wk1d_b_p,(2,mgb,n3,1))

   if(option/=3)then

     if(lluse_ndo)  then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(n3,ngbin,wk1d_a_p)
       do i3=1,n3
         do igb=1,ngbin
           wk1d_a_p(1,igb,i3,1)=zero
           wk1d_a_p(2,igb,i3,1)=zero
         end do
       end do
!$OMP END PARALLEL DO

!      Insert fofgin_p into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(fofgin_p,indpw_kin,npwin,wk1d_a_p)
       do ig=1,npwin
         igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
         wk1d_a_p(1,igb,i3,1)=fofgin_p(1,ig)
         wk1d_a_p(2,igb,i3,1)=fofgin_p(2,ig)
       end do
!$OMP END PARALLEL DO

!      Go from wk1d_a_p to wk1d_b_p, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a_p,wk1d_b_p)&
!$OMP&PRIVATE(igb,igbmax)
       do igb=1,ngbin,lotin
         igbmax=min(igb+lotin-1,ngbin)
!        Go from wk1d_a_p to wk1d_b_p, using 1D FFTs on the z direction
!        However, due to special packing of data, use routine ffty
         call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a_p,wk1d_b_p, &
&         trig3,aft3,now3,bef3,one,ind3,ic3)
       end do
!$OMP END PARALLEL DO

     end if ! lluse_ndo

!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(n3,ngbin,wk1d_a)
     do i3=1,n3
       do igb=1,ngbin
         wk1d_a(1,igb,i3,1)=zero
         wk1d_a(2,igb,i3,1)=zero
       end do
     end do
!$OMP END PARALLEL DO

!    Insert fofgin into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(fofgin,indpw_kin,npwin,wk1d_a)
     do ig=1,npwin
       igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
       wk1d_a(1,igb,i3,1)=fofgin(1,ig)
       wk1d_a(2,igb,i3,1)=fofgin(2,ig)
     end do
!$OMP END PARALLEL DO

!    Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbin,lotin
       igbmax=min(igb+lotin-1,ngbin)
!      Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a,wk1d_b, &
&       trig3,aft3,now3,bef3,one,ind3,ic3)
     end do
!$OMP END PARALLEL DO

   end if !  if(option/=3)

!  Do-loop on the planes stacked in the z direction
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP&SHARED(aft1,aft2,aft4,aft5,bef1,bef2,bef4,bef5,cplex,denpot) &
!$OMP&SHARED(fftcache,fofr,gboundin,gboundout)&
!$OMP&SHARED(ic1,ic2,ic4,ic5,ind1,ind2,ind4) &
!$OMP&SHARED(ind5,indpw_kin,indpw_kout,lluse_ndo,mgb,n1,n2,n3,n4,n5,ngbin) &
!$OMP&SHARED(ngbout,now1,now2,now4,now5,option,trig1,trig2,trig4,trig5) &
!$OMP&SHARED(weight_r,weight_i,weight_2,wk1d_a,wk1d_b,wk1d_b_p)

!  Allocate two 2-dimensional work arrays
   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_a_p,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b_p,(2,n4,n5,1))
!$OMP DO
   do i3=1,n3

     if(option/=3)then
       if(lluse_ndo)  then
!        Zero the values on the current plane
!        wk2d_a_p(1:2,1:n1,1:n2,1)=zero
         do i2=1,n2
           do i1=1,n1
             wk2d_a_p(1,i1,i2,1)=zero
             wk2d_a_p(2,i1,i2,1)=zero
           end do
         end do
!        Copy the data in the current plane
         do igb=1,ngbin
           i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
           wk2d_a_p(1,i1,i2,1)=wk1d_b_p(1,igb,i3,1)
           wk2d_a_p(2,i1,i2,1)=wk1d_b_p(2,igb,i3,1)
         end do
!        Perform x transform, taking into account arrays of zeros
         g2min=gboundin(3,1) ; g2max=gboundin(4,1)
         if ( g2min+n2 >= g2max+2 ) then
           do i2=g2max+2,g2min+n2
             do i1=1,n1
               wk2d_b_p(1,i1,i2,1)=zero
               wk2d_b_p(2,i1,i2,1)=zero
             end do
           end do
         end if
         gbound_dum(1)=1 ; gbound_dum(2)=1
         gbound_dum(3)=g2min ; gbound_dum(4)=g2max
         call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_a_p,wk2d_b_p,&
&         trig1,aft1,now1,bef1,one,ind1,ic1,gbound_dum)
!        Perform y transform
         n1i=1
         call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_b_p,wk2d_a_p, &
&         trig2,aft2,now2,bef2,one,ind2,ic2)
!        The wave function is now in real space, for the current plane
       end if  ! lluse_ndo

!      Zero the values on the current plane
!      wk2d_a(1:2,1:n1,1:n2,1)=zero
       do i2=1,n2
         do i1=1,n1
           wk2d_a(1,i1,i2,1)=zero
           wk2d_a(2,i1,i2,1)=zero
         end do
       end do
!      Copy the data in the current plane
       do igb=1,ngbin
         i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
         wk2d_a(1,i1,i2,1)=wk1d_b(1,igb,i3,1)
         wk2d_a(2,i1,i2,1)=wk1d_b(2,igb,i3,1)
       end do
!      Perform x transform, taking into account arrays of zeros
       g2min=gboundin(3,1) ; g2max=gboundin(4,1)
       if ( g2min+n2 >= g2max+2 ) then
         do i2=g2max+2,g2min+n2
           do i1=1,n1
             wk2d_b(1,i1,i2,1)=zero
             wk2d_b(2,i1,i2,1)=zero
           end do
         end do
       end if
       gbound_dum(1)=1 ; gbound_dum(2)=1
       gbound_dum(3)=g2min ; gbound_dum(4)=g2max
       call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_a,wk2d_b,&
&       trig1,aft1,now1,bef1,one,ind1,ic1,gbound_dum)
!      Perform y transform
       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_b,wk2d_a, &
&       trig2,aft2,now2,bef2,one,ind2,ic2)
!      The wave function is now in real space, for the current plane
     end if

     if(option==0)then
!      Copy the transformed function at the right place
       do i2=1,n2
         do i1=1,n1
           fofr(1,i1,i2,i3)=wk2d_a(1,i1,i2,1)
           fofr(2,i1,i2,i3)=wk2d_a(2,i1,i2,1)
         end do
       end do
     end if

     if(option==1)then
!      Accumulate density
       do i2=1,n2
         do i1=1,n1
           if(lluse_ndo)  then
             denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&             weight_r*(wk2d_a(1,i1,i2,1)*wk2d_a_p(1,i1,i2,1)&
&             +wk2d_a(2,i1,i2,1)*wk2d_a_p(2,i1,i2,1))
             if(present(weight_2)) then
               denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&               weight_2*(wk2d_a_p(2,i1,i2,1)*wk2d_a(1,i1,i2,1)&
&               -wk2d_a_p(1,i1,i2,1)*wk2d_a(2,i1,i2,1))
             end if
           else
             denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&             weight_r*wk2d_a(1,i1,i2,1)**2+ weight_i*wk2d_a(2,i1,i2,1)**2
           end if
         end do
       end do
     end if

     if(option==2)then
!      Apply local potential
       if(cplex==1)then
         do i2=1,n2
           do i1=1,n1
             wk2d_a(1,i1,i2,1)=denpot(i1,i2,i3)*wk2d_a(1,i1,i2,1)
             wk2d_a(2,i1,i2,1)=denpot(i1,i2,i3)*wk2d_a(2,i1,i2,1)
           end do
         end do
       else
         do i2=1,n2
           do i1=1,n1
             wkre=wk2d_a(1,i1,i2,1)
             wkim=wk2d_a(2,i1,i2,1)
             wk2d_a(1,i1,i2,1)=denpot(2*i1-1,i2,i3)*wkre &
&             -denpot(2*i1  ,i2,i3)*wkim
             wk2d_a(2,i1,i2,1)=denpot(2*i1-1,i2,i3)*wkim &
&             +denpot(2*i1  ,i2,i3)*wkre
           end do
         end do
       end if
     end if

     if(option==3)then
!      Copy the function to be tranformed at the right place
       do i2=1,n2
         do i1=1,n1
           wk2d_a(1,i1,i2,1)=fofr(1,i1,i2,i3)
           wk2d_a(2,i1,i2,1)=fofr(2,i1,i2,i3)
         end do
       end do
     end if

     if(option==2 .or. option==3)then
!      Perform y transform
       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_a,wk2d_b, &
&       trig5,aft5,now5,bef5,-one,ind5,ic5)
!      Perform x transform, taking into account arrays of zeros
       gbound_dum(1)=1 ; gbound_dum(2)=1
       gbound_dum(3)=gboundout(3,1) ; gbound_dum(4)=gboundout(4,1)
       call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_b,wk2d_a,&
&       trig4,aft4,now4,bef4,-one,ind4,ic4,gbound_dum)
!      Copy the data from the current plane to wk1d_b
       do igb=1,ngbout
         i1=indpw_kout(1,igb) ; i2=indpw_kout(2,igb)
         wk1d_b(1,igb,i3,1)=wk2d_a(1,i1,i2,1)
         wk1d_b(2,igb,i3,1)=wk2d_a(2,i1,i2,1)
       end do
     end if

!    End loop on planes
   end do
!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_a_p)
   ABI_DEALLOCATE(wk2d_b_p)
!$OMP END PARALLEL

   if(option==2 .or. option==3)then

!    Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft6,bef6,fftcache,ind6,ic6,lotout,mgb)&
!$OMP&SHARED(ngbout,now6,n3,trig6,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbout,lotout
       igbmax=min(igb+lotout-1,ngbout)
!      Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_b,wk1d_a, &
&       trig6,aft6,now6,bef6,-one,ind6,ic6)

     end do
!$OMP END PARALLEL DO

!    Transfer the data in the output array, after normalization
     norm=1.d0/dble(nfftot)
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(fofgout,indpw_kout,norm,npwout,wk1d_a)
     do ig=1,npwout
       igb=indpw_kout(4,ig) ; i3=indpw_kout(3,ig)
       fofgout(1,ig)=wk1d_a(1,igb,i3,1)*norm
       fofgout(2,ig)=wk1d_a(2,igb,i3,1)*norm
     end do
!$OMP END PARALLEL DO
   end if

   ABI_DEALLOCATE(wk1d_a)
   ABI_DEALLOCATE(wk1d_b)
   ABI_DEALLOCATE(wk1d_a_p)
   ABI_DEALLOCATE(wk1d_b_p)

!  End general k-point part
 end if

!------------------------------------------------------------------
!Here, use of time-reversal symmetry

 if(istwf_k>=2)then

   n1half1=n1/2+1 ; n1halfm=(n1+1)/2
   n2half1=n2/2+1
!  n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
   n4half1=(n1half1/2)*2+1
   n5half1=(n2half1/2)*2+1
!  Note that the z transform will appear as a y transform
   ABI_ALLOCATE(wk1d_a,(2,mgb,n3,1))
   ABI_ALLOCATE(wk1d_b,(2,mgb,n3,1))

   if(istwf_k/=2)then
     ABI_ALLOCATE(pha1,(2,n1))
     ABI_ALLOCATE(pha2,(2,n2))
     ABI_ALLOCATE(pha3,(3,n3))
     do i1=1,n1
       pha1(1,i1)=cos(dble(i1-1)*pi/dble(n1))
       pha1(2,i1)=sin(dble(i1-1)*pi/dble(n1))
     end do
     do i2=1,n2
       pha2(1,i2)=cos(dble(i2-1)*pi/dble(n2))
       pha2(2,i2)=sin(dble(i2-1)*pi/dble(n2))
     end do
     do i3=1,n3
       pha3(1,i3)=cos(dble(i3-1)*pi/dble(n3))
       pha3(2,i3)=sin(dble(i3-1)*pi/dble(n3))
     end do
   end if

   if(option/=3)then

!    Zero the components of wk1d_a
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(n3,ngbin,wk1d_a)
     do i3=1,n3
       do igb=1,ngbin
         wk1d_a(1,igb,i3,1)=zero
         wk1d_a(2,igb,i3,1)=zero
       end do
     end do
!$OMP END PARALLEL DO

!    Insert fofgin into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(fofgin,indpw_kin,npwin,wk1d_a)
     do ig=1,npwin
       igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
       wk1d_a(1,igb,i3,1)=fofgin(1,ig)
       wk1d_a(2,igb,i3,1)=fofgin(2,ig)
     end do
!$OMP END PARALLEL DO

!    Must complete the i2=1 plane when $k_y \equiv 0$

!    Take care of i1=1 when $k_x \equiv 0$
     if(istwf_k==2)then
!      Take care of i1=1
       do i3=n3/2+1,n3
         i3inv=n3+2-i3
         wk1d_a(1,1,i3,1)= wk1d_a(1,1,i3inv,1)
         wk1d_a(2,1,i3,1)=-wk1d_a(2,1,i3inv,1)
       end do
     else if(istwf_k==4)then
!      Take care of i1=1
       do i3=n3/2+1,n3
         i3inv=n3+1-i3
         wk1d_a(1,1,i3,1)= wk1d_a(1,1,i3inv,1)
         wk1d_a(2,1,i3,1)=-wk1d_a(2,1,i3inv,1)
       end do
     end if

!    Now, take care of other i1 values, except i3==1 when $k_z \equiv 0$
     i1max=gboundin(6,1)+1
     if(istwf_k==2)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(i1max,n3,wk1d_a)
       do igb=2,2*i1max-1
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+2-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     else if(istwf_k==3)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(i1max,n3,wk1d_a)
       do igb=1,2*i1max
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+2-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     else if(istwf_k==4)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(i1max,n3,wk1d_a)
       do igb=2,2*i1max-1
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+1-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     else if(istwf_k==5)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(i1max,n3,wk1d_a)
       do igb=1,2*i1max
         igb_inv=2*i1max+1-igb
         do i3=n3/2+1,n3
           i3inv=n3+1-i3
           wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
           wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
         end do
       end do
!$OMP END PARALLEL DO

     end if

!    Now, i3==1
     if(istwf_k==2)then
       do igb=2,i1max
         igb_inv=2*i1max+1-igb
         wk1d_a(1,igb_inv,1,1)= wk1d_a(1,igb,1,1)
         wk1d_a(2,igb_inv,1,1)=-wk1d_a(2,igb,1,1)
       end do
     else if(istwf_k==3)then
       do igb=1,i1max
         igb_inv=2*i1max+1-igb
         wk1d_a(1,igb_inv,1,1)= wk1d_a(1,igb,1,1)
         wk1d_a(2,igb_inv,1,1)=-wk1d_a(2,igb,1,1)
       end do
     end if

!    Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbin,lotin
       igbmax=min(igb+lotin-1,ngbin)
!      Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a,wk1d_b, &
&       trig3,aft3,now3,bef3,one,ind3,ic3)
     end do
!$OMP END PARALLEL DO

!    Change the phase if $k_z \neq 0$
     if(istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(ngbin,n3,pha3,wk1d_b)
       do i3=1,n3
         phar=pha3(1,i3)
         phai=pha3(2,i3)
         do igb=1,ngbin
           ar=wk1d_b(1,igb,i3,1)
           ai=wk1d_b(2,igb,i3,1)
           wk1d_b(1,igb,i3,1)=phar*ar-phai*ai
           wk1d_b(2,igb,i3,1)=phai*ar+phar*ai
         end do
       end do
!$OMP END PARALLEL DO
     end if

   end if !  if(option/=3)

!  Do-loop on the planes stacked in the z direction

!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP&SHARED(aft1,aft2,aft4,aft5,bef1,bef2,bef4,bef5,denpot) &
!$OMP&SHARED(fftcache,fofr,gboundin,ic1,ic2,ic4,ic5,ind1,ind2,ind4,ind5) &
!$OMP&SHARED(indpw_kin,indpw_kout,istwf_k,mgb,n1,n1half1) &
!$OMP&SHARED(n1halfm,n2,n2half1,n3,n4,n5,ngbin,ngbout) &
!$OMP&SHARED(now1,now2,now4,now5,option,pha1,pha2,trig1) &
!$OMP&SHARED(trig2,trig4,trig5,weight_r,weight_i,wk1d_a,wk1d_b)

!  Allocate two 2-dimensional work arrays
   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_c,(2,2*n1halfm,n5,1))
   ABI_ALLOCATE(wk2d_d,(2,2*n1halfm,n5,1))
!$OMP DO
   do i3=1,n3

     g2max=gboundin(4,1)

     if(option/=3)then
!      Zero the values on the current plane : need only from i2=1 to g2max+1
       do i2=1,g2max+1
         do i1=1,n1
           wk2d_a(1,i1,i2,1)=zero
           wk2d_a(2,i1,i2,1)=zero
         end do
       end do

!      Copy the data in the current plane
       do igb=1,ngbin
         i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
         wk2d_a(1,i1,i2,1)=wk1d_b(1,igb,i3,1)
         wk2d_a(2,i1,i2,1)=wk1d_b(2,igb,i3,1)
       end do

!      Perform x transform, taking into account arrays of zeros
       call sg_fftx(fftcache,mfac,mg,n4,n5,1,g2max+1,1,wk2d_a,wk2d_b,&
&       trig1,aft1,now1,bef1,one,ind1,ic1)

!      Change the phase if $k_x \neq 0$
       if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9)then
         do i1=1,n1
           phar=pha1(1,i1)
           phai=pha1(2,i1)
           do i2=1,g2max+1
             ar=wk2d_b(1,i1,i2,1)
             ai=wk2d_b(2,i1,i2,1)
             wk2d_b(1,i1,i2,1)=phar*ar-phai*ai
             wk2d_b(2,i1,i2,1)=phai*ar+phar*ai
           end do
         end do
       end if

!      Compute symmetric and antisymmetric combinations
       if(istwf_k>=2 .and. istwf_k<=5)then
         do i1=1,n1half1-1
           wk2d_a(1,i1,1,1)=wk2d_b(1,2*i1-1,1,1)
           wk2d_a(2,i1,1,1)=wk2d_b(1,2*i1  ,1,1)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           wk2d_a(1,n1half1,1,1)=wk2d_b(1,n1,1,1)
           wk2d_a(2,n1half1,1,1)=zero
         end if
         ii2=2
       else
         ii2=1
       end if
       if( g2max+1 >= ii2)then
         do i2=ii2,g2max+1
           do i1=1,n1half1-1
             wk2d_a(1,i1,i2,1)=        wk2d_b(1,2*i1-1,i2,1)-wk2d_b(2,2*i1,i2,1)
             wk2d_a(2,i1,i2,1)=        wk2d_b(2,2*i1-1,i2,1)+wk2d_b(1,2*i1,i2,1)
             wk2d_a(1,i1,n2+ii2-i2,1)= wk2d_b(1,2*i1-1,i2,1)+wk2d_b(2,2*i1,i2,1)
             wk2d_a(2,i1,n2+ii2-i2,1)=-wk2d_b(2,2*i1-1,i2,1)+wk2d_b(1,2*i1,i2,1)
           end do
           if((2*n1half1-2)/=n1)then
             wk2d_a(1,n1half1,i2,1)=        wk2d_b(1,n1,i2,1)
             wk2d_a(2,n1half1,i2,1)=        wk2d_b(2,n1,i2,1)
             wk2d_a(1,n1half1,n2+ii2-i2,1)= wk2d_b(1,n1,i2,1)
             wk2d_a(2,n1half1,n2+ii2-i2,1)=-wk2d_b(2,n1,i2,1)
           end if
         end do
       end if
       if ( n2half1 >= g2max+2 ) then
         do i2=g2max+2,n2half1
           do i1=1,n1half1-1
             wk2d_a(1,i1,i2,1)=zero
             wk2d_a(2,i1,i2,1)=zero
             wk2d_a(1,i1,n2+ii2-i2,1)=zero
             wk2d_a(2,i1,n2+ii2-i2,1)=zero
           end do
           if((2*n1half1-2)/=n1)then
             wk2d_a(1,n1half1,i2,1)=zero
             wk2d_a(2,n1half1,i2,1)=zero
             wk2d_a(1,n1half1,n2+ii2-i2,1)=zero
             wk2d_a(2,n1half1,n2+ii2-i2,1)=zero
           end if
         end do
       end if

       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1halfm,1,1,wk2d_a,wk2d_b,&
&       trig2,aft2,now2,bef2,one,ind2,ic2)

!      Change the phase if $k_y \neq 0$
       if(istwf_k>=6 .and. istwf_k<=9)then
         do i2=1,n2
           phar=pha2(1,i2)
           phai=pha2(2,i2)
           do i1=1,n1halfm
             ar=wk2d_b(1,i1,i2,1)
             ai=wk2d_b(2,i1,i2,1)
             wk2d_b(1,i1,i2,1)= phar*ar-phai*ai
             wk2d_b(2,i1,i2,1)= phai*ar+phar*ai
           end do
         end do
       end if

     end if ! option/=3

!    The wave function is now in real space, for the current plane,
!    represented by REAL numbers, although packed in the complex array wk2d_b

     g2max=gboundin(4,1)

     if(option==0)then
!      This option is only permitted for istwf_k==2 (Gamma point)
!      Copy the transformed function at the right place
       do i2=1,n2
         do i1=1,n1half1-1
           fofr(1,2*i1-1,i2,i3)=wk2d_b(1,i1,i2,1)
           fofr(1,2*i1  ,i2,i3)=wk2d_b(2,i1,i2,1)
           fofr(2,2*i1-1,i2,i3)=zero
           fofr(2,2*i1  ,i2,i3)=zero
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           fofr(1,n1,i2,i3)=wk2d_b(1,n1half1,i2,1)
           fofr(2,n1,i2,i3)=zero
         end if
       end do
     end if

     if(option==1)then
!      Accumulate density
       do i2=1,n2
         do i1=1,n1half1-1
           denpot(2*i1-1,i2,i3)=denpot(2*i1-1,i2,i3)+weight_r*wk2d_b(1,i1,i2,1)**2
           denpot(2*i1  ,i2,i3)=denpot(2*i1  ,i2,i3)+weight_i*wk2d_b(2,i1,i2,1)**2
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           denpot(n1,i2,i3)=denpot(n1,i2,i3)+weight_r*wk2d_b(1,n1half1,i2,1)**2
!          not use in DMFT because istwfk required to be one.
         end if
       end do
     end if

     if(option==2)then
!      Apply local potential
       do i2=1,n2
         do i1=1,n1half1-1
           wk2d_a(1,i1,i2,1)=denpot(2*i1-1,i2,i3)*wk2d_b(1,i1,i2,1)
           wk2d_a(2,i1,i2,1)=denpot(2*i1  ,i2,i3)*wk2d_b(2,i1,i2,1)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           wk2d_a(1,n1half1,i2,1)=denpot(n1,i2,i3)*wk2d_b(1,n1half1,i2,1)
           wk2d_a(2,n1half1,i2,1)=zero
         end if
       end do
     end if

     if(option==3)then
!      This option is only permitted for istwf_k==2 (Gamma point)
!      Copy the transformed function at the right place
       do i2=1,n2
         do i1=1,n1half1-1
           wk2d_b(1,i1,i2,1)=fofr(1,2*i1-1,i2,i3)
           wk2d_b(2,i1,i2,1)=fofr(1,2*i1  ,i2,i3)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           wk2d_b(1,n1half1,i2,1)=fofr(1,n1,i2,i3)
         end if
       end do
     end if

     if(option==2 .or. option==3)then
!      Change the phase if $k_y \neq 0$
       if(istwf_k>=6 .and. istwf_k<=9)then
         do i2=1,n2
           phar=pha2(1,i2)
           phai=pha2(2,i2)
           do i1=1,n1halfm
             ar=wk2d_a(1,i1,i2,1)
             ai=wk2d_a(2,i1,i2,1)
             wk2d_a(1,i1,i2,1)= phar*ar+phai*ai
             wk2d_a(2,i1,i2,1)=-phai*ar+phar*ai
           end do
         end do
       end if

!      Perform y transform
       n1i=1
       call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1halfm,1,1,wk2d_a,wk2d_b, &
&       trig5,aft5,now5,bef5,-one,ind5,ic5)

!      Decompose symmetric and antisymmetric parts
       if(istwf_k>=2 .and. istwf_k<=5)then
         do i1=1,n1halfm
           wk2d_c(1,2*i1-1,1,1)=wk2d_b(1,i1,1,1)
           wk2d_c(2,2*i1-1,1,1)=zero
           wk2d_c(1,2*i1,1,1)=wk2d_b(2,i1,1,1)
           wk2d_c(2,2*i1,1,1)=zero
         end do
         ii2=2
       else
         ii2=1
       end if
       do i2=ii2,g2max+1
         do i1=1,n1halfm
           wk2d_c(1,2*i1-1,i2,1)=(wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+ii2-i2,1))*0.5d0
           wk2d_c(2,2*i1-1,i2,1)=(wk2d_b(2,i1,i2,1)-wk2d_b(2,i1,n2+ii2-i2,1))*0.5d0
           wk2d_c(1,2*i1,i2,1)= ( wk2d_b(2,i1,i2,1)+wk2d_b(2,i1,n2+ii2-i2,1))*0.5d0
           wk2d_c(2,2*i1,i2,1)= (-wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+ii2-i2,1))*0.5d0
         end do
       end do

!      Change the phase if $k_x \neq 0$
       if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9 )then
         do i1=1,n1
           phar=pha1(1,i1)
           phai=pha1(2,i1)
           do i2=1,g2max+1
             ar=wk2d_c(1,i1,i2,1)
             ai=wk2d_c(2,i1,i2,1)
             wk2d_c(1,i1,i2,1)= phar*ar+phai*ai
             wk2d_c(2,i1,i2,1)=-phai*ar+phar*ai
           end do
         end do
       end if

!      Perform x transform : for y=1 to g2max+1, to benefit from zeros
       call sg_fftx(fftcache,mfac,mg,2*n1halfm,n5,1,g2max+1,1,wk2d_c,wk2d_d,&
&       trig4,aft4,now4,bef4,-one,ind4,ic4)

!      Copy the data from the current plane to wk1d_b
       do igb=1,ngbout
         i1=indpw_kout(1,igb) ; i2=indpw_kout(2,igb)
         wk1d_b(1,igb,i3,1)=wk2d_d(1,i1,i2,1)
         wk1d_b(2,igb,i3,1)=wk2d_d(2,i1,i2,1)
       end do

     end if ! option==2 or 3

!    End loop on planes
   end do

!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_c)
   ABI_DEALLOCATE(wk2d_d)
!$OMP END PARALLEL

   if(option==2 .or. option==3)then

!    Change the phase if $k_z \neq 0$
     if(istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(ngbout,n3,pha3,wk1d_b)
       do i3=1,n3
         phar=pha3(1,i3)
         phai=pha3(2,i3)
         do igb=1,ngbout
           ar=wk1d_b(1,igb,i3,1)
           ai=wk1d_b(2,igb,i3,1)
           wk1d_b(1,igb,i3,1)= phar*ar+phai*ai
           wk1d_b(2,igb,i3,1)=-phai*ar+phar*ai
         end do
       end do
!$OMP END PARALLEL DO
     end if

!    Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!    However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft6,bef6,fftcache,ind6,ic6,lotout,mgb)&
!$OMP&SHARED(ngbout,now6,n3,trig6,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
     do igb=1,ngbout,lotout
       igbmax=min(igb+lotout-1,ngbout)
!      Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!      However, due to special packing of data, use routine ffty
       call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_b,wk1d_a, &
&       trig6,aft6,now6,bef6,-one,ind6,ic6)

     end do
!$OMP END PARALLEL DO

!    Transfer the data in the output array, after normalization
     norm=1.d0/dble(nfftot)
!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP&SHARED(fofgout,indpw_kout,norm,npwout,wk1d_a)
     do ig=1,npwout
       igb=indpw_kout(4,ig) ; i3=indpw_kout(3,ig)
       fofgout(1,ig)=wk1d_a(1,igb,i3,1)*norm
       fofgout(2,ig)=wk1d_a(2,igb,i3,1)*norm
     end do
!$OMP END PARALLEL DO

   end if

   ABI_DEALLOCATE(wk1d_a)
   ABI_DEALLOCATE(wk1d_b)

   if(istwf_k/=2)then
     ABI_DEALLOCATE(pha1)
     ABI_DEALLOCATE(pha2)
     ABI_DEALLOCATE(pha3)
   end if

!  End time-reversal symmetry
 end if

 if(option/=3) then
   ABI_DEALLOCATE(indpw_kin)
 end if
 if(option==2 .or. option==3) then
   ABI_DEALLOCATE(indpw_kout)
 end if

 !DBG_EXIT("COLL")

end subroutine sg_fftrisc_2
!!***

!----------------------------------------------------------------------

!!****f* m_sgfft/sg_poisson
!! NAME
!! sg_poisson
!!
!! FUNCTION
!!  Solve the Poisson equation in G-space given the density, n(r),
!!  in real space of the FFT box.
!!
!! INPUTS
!! fftcache=size of the cache (kB)
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nx,ny,nz=Number of FFT points along the three directions.
!! ldx,ldy,ldz=Leading dimension of the array nr and vg.
!! ndat = Number of densities
!! vg(nx*ny*nz)=Potential in reciprocal space.
!!
!! SIDE EFFECTS
!! nr(cplex*ldx*ldy*ldz*ndat)
!!    input: n(r) (real or complex)
!!    output: the hartree potential in real space
!!
!! NOTES
!!   vg is given on the FFT mesh instead of the augmented mesh [ldx,ldy,ldz]
!!   in order to simplify the interface with the other routines operating of vg
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      sg_fft_cc
!!
!! SOURCE

subroutine sg_poisson(fftcache,cplex,nx,ny,nz,ldx,ldy,ldz,ndat,vg,nr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,cplex,nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 real(dp),intent(inout) :: nr(cplex*ldx*ldy*ldz*ndat)
 real(dp),intent(in) :: vg(nx*ny*nz)

!Local variables-------------------------------
 integer,parameter :: ndat1=1
 integer :: ii,jj,kk,ifft,dat,ptr,ig
 real(dp) :: fft_fact
!arrays
 real(dp),allocatable :: work(:,:)

! *************************************************************************

 fft_fact = one/(nx*ny*nz)

 ABI_CHECK(cplex==2,"cplex!=2 not coded")

 ABI_MALLOC(work, (2,ldx*ldy*ldz))

 do dat=1,ndat
   ! n(r) --> n(G)
   ptr = 1 + (dat-1)*cplex*ldx*ldy*ldz
   call sg_fft_cc(fftcache,nx,ny,nz,ldx,ldy,ldz,ndat1,-1,nr(ptr),work)

   ! Multiply by v(G)
   ig = 0
   do kk=1,nz
     do jj=1,ny
       do ii=1,nx
         ig = ig + 1
         ifft = ii + (jj-1)*ldx + (kk-1)*ldx*ldy
         work(1:2,ifft) = work(1:2,ifft) * vg(ig) * fft_fact
      end do
     end do
   end do

   ! compute vh(r)
   call sg_fft_cc(fftcache,nx,ny,nz,ldx,ldy,ldz,ndat1,+1,work,nr(ptr))
 end do

 ABI_FREE(work)

end subroutine sg_poisson
!!***

END MODULE m_sgfft
