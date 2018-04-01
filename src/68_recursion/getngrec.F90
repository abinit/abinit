!{\src2tex{textfont=tt}}
!!****f* ABINIT/getngrec
!! NAME
!! getngrec
!!
!! FUNCTION
!! This routine computes the fft box for the recursion method, accordingly to the troncation radius.
!! It is quite similar to getng, but : 
!!     - there is no xboxcut and ecut consistency
!!     - ngfft (the initial fft box) is the maximum fft box
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ngfft(18)=non truncated fft box
!!  mgfft=maximum of ngfft(1:3)
!!  inf_rmet=define the infinitesimal metric : rprimd*(transpose(rprimd)), divided by the number of discretisation point
!!  recrcut=truncating
!!  delta=to obtain radius of truncation
!!
!! OUTPUT
!!  ngfftrec=truncated fft box
!!  nfftrec= truncated nfft 
!!  tronc=True if truncation is made 
!!
!! SIDE EFFECTS
!! 
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      sort_int,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getngrec(ngfft,rmet,ngfftrec,nfftrec,recrcut,delta,tronc)

 use defs_basis
 use m_profiling_abi
 use m_sort

 use m_time,        only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getngrec'
 use interfaces_18_timing
!End of the abilint section

implicit none

!Arguments -------------------------------
!scalars
real(dp),intent(in) :: recrcut,delta
integer,intent(out) :: nfftrec
logical,intent(out) :: tronc
!arrays
integer,intent(in)  :: ngfft(18)
real(dp),intent(in) :: rmet(3,3)
integer,intent(out) :: ngfftrec(18)

!Local variables-------------------------------
!scalars
integer :: ii,iimin,index,jj,jjmin,kk,kkmin,largest_ngfftrec,maxpow11,maxpow2
integer :: maxpow3,maxpow5,maxpow7,mmsrch,plane
real(dp) :: dsm,dsp,dsqmin,rtroncat
!arrays
integer :: get_ngfftrec(3),maxsrch(3),minsrch(3)
integer,allocatable :: iperm(:),srch(:)
real(dp) :: tsec(2)
real(dp) :: inf_rmet(3,3)
! *************************************************************************

 call timab(602,1,tsec)

 ngfftrec(:) = ngfft(:) 

 if (recrcut>tol14) then  !default value dtset%recrcut = zero means no troncation 
   rtroncat = recrcut+delta
   get_ngfftrec(:)=1
   plane = 1 

   do ii=1,3
     inf_rmet(ii,:) = rmet(ii,:)/(real(ngfft(1:3)*ngfft(ii),dp))
   end do


!  minimum value of ngfftrec
   do ii = 1,3
     ngfftrec(ii)=floor(2*rtroncat/sqrt(inf_rmet(ii,ii)))+1  !minimum value
     if(ngfftrec(ii)>=ngfft(ii))then
       ngfftrec(ii)=ngfft(ii)
       get_ngfftrec(ii)=0
     end if
   end do


   if(sum(get_ngfftrec)/=0)then
     largest_ngfftrec=maxval(ngfft(1:3))
     maxpow2=int(log(largest_ngfftrec+0.5d0)/log(two))
     maxpow3=int(log(largest_ngfftrec+0.5d0)/log(three))
     maxpow5=int(log(largest_ngfftrec+0.5d0)/log(five))
     maxpow7=0
     maxpow11=0
     mmsrch=(maxpow2+1)*(maxpow3+1)*(maxpow5+1)*(maxpow7+1)*(maxpow11+1)
     ABI_ALLOCATE(srch,(mmsrch))
     ABI_ALLOCATE(iperm,(mmsrch))
!    Factors of 2
     srch(1)=1
     do ii=1,maxpow2
       srch(ii+1)=srch(ii)*2
     end do
!    Factors of 3
     index=maxpow2+1
     if(maxpow3>0)then
       do ii=1,maxpow3
         srch(1+ii*index:(ii+1)*index)=3*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    Factors of 5
     index=(maxpow3+1)*index
     if(maxpow5>0)then
       do ii=1,maxpow5
         srch(1+ii*index:(ii+1)*index)=5*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    Factors of 7
     index=(maxpow5+1)*index
     if(maxpow7>0)then
       do ii=1,maxpow7
         srch(1+ii*index:(ii+1)*index)=7*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    Factors of 11
     index=(maxpow7+1)*index
     if(maxpow11>0)then
       do ii=1,maxpow11
         srch(1+ii*index:(ii+1)*index)=11*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    srch is the set of allowed ngfftrec values

     call sort_int(mmsrch,srch,iperm)
     ABI_DEALLOCATE(iperm)

     do ii=1,3
       if(get_ngfftrec(ii)==1)then
         do jj=1,mmsrch
           if(srch(jj)>=ngfftrec(ii))then
             minsrch(ii)=jj
             ngfftrec(ii)=srch(jj)
             exit
           end if
         end do
         do jj=minsrch(ii),mmsrch
           if(srch(jj)>ngfft(ii))then 
!            since ngfftrec(ii)<ngfft(ii) for get_ngfftrec(ii)==1, 
!            and srch(mmsrch)maxval(ngfft(1:3)), 
!            that will appens in the range minsrch(ii),mmsrch
             maxsrch(ii)=jj-1
             exit
           end if
         end do
       end if
!      since ngfft(ii) is in srch, we have here srch(maxsrch(ii))=ngfft(ii)
!      minsrch(ii), maxsrch(ii) is the range of index of srch in which we can 
!      search ngfftrec(ii)

       if(ngfftrec(ii)>=ngfft(ii))then
         ngfftrec(ii)=ngfft(ii)
         get_ngfftrec(ii)=0
       end if
     end do
   end if

!  verify that the entiere truncation sphere is in the fft box ; 
!  but only in the dimension in which we do not consider the entiere fft box
   do while(sum(get_ngfftrec)/=0)  !again...

!    determining the minimum distance between 0 and the boundary 
!    of the fft box
!    quite similar to the subroutine "bound", but only over the plane which 
!    are not the whole fft box
     dsqmin=dsq_rec(ngfftrec(1)/2,-ngfftrec(2)/2,-ngfftrec(3)/2,inf_rmet)+0.01d0

     if(get_ngfftrec(1)/=0)then
!      look at +/- g1 planes:
       do jj=-ngfftrec(2)/2,ngfftrec(2)/2
         do kk=-ngfftrec(3)/2,ngfftrec(3)/2
           dsp = dsq_rec(ngfftrec(1)/2, jj, kk,inf_rmet)
           dsm = dsq_rec( - ngfftrec(1)/2, jj, kk,inf_rmet)
           if (dsp<dsqmin) then
             dsqmin = dsp
             iimin = ngfftrec(1)/2
             jjmin = jj
             kkmin = kk
             plane=1
           end if
           if (dsm<dsqmin) then
             dsqmin = dsm
             iimin =  - ngfftrec(1)/2
             jjmin = jj
             kkmin = kk
             plane=1
           end if
         end do
       end do
     end if

     if(get_ngfftrec(2)/=0)then
!      +/- g2 planes:
       do ii=-ngfftrec(1)/2,ngfftrec(1)/2
         do kk=-ngfftrec(3)/2,ngfftrec(3)/2
           dsp = dsq_rec(ii,ngfftrec(2)/2,kk,inf_rmet)
           dsm = dsq_rec(ii,-ngfftrec(2)/2,kk,inf_rmet)
           if (dsp<dsqmin) then
             dsqmin = dsp
             iimin = ii
             jjmin = ngfftrec(2)/2
             kkmin = kk
             plane=2
           end if
           if (dsm<dsqmin) then
             dsqmin = dsm
             iimin = ii
             jjmin =  - ngfftrec(2)/2
             kkmin = kk
             plane=2
           end if
         end do
       end do
     end if

     if(get_ngfftrec(3)/=0)then
!      +/- g3 planes:
       do ii=-ngfftrec(1)/2,ngfftrec(1)/2
         do jj=-ngfftrec(2)/2,ngfftrec(2)/2
           dsp = dsq_rec(ii,jj,ngfftrec(3)/2,inf_rmet)
           dsm = dsq_rec(ii,jj,-ngfftrec(3)/2,inf_rmet)
           if (dsp<dsqmin) then
             dsqmin = dsp
             iimin = ii
             jjmin = jj
             kkmin = ngfftrec(3)/2
             plane=3
           end if
           if (dsm<dsqmin) then
             dsqmin = dsm
             iimin = ii
             jjmin = jj
             kkmin =  - ngfftrec(3)/2
             plane=3
           end if
         end do
       end do
     end if

     if(dsqmin>=rtroncat)then
       get_ngfftrec=0
       exit
     end if

!    Fix nearest boundary
     do ii=minsrch(plane),maxsrch(plane)
       if (srch(ii)>=ngfftrec(plane)) then
!        redefine ngfft(plane) to next higher choice
         ngfftrec(plane)=srch(ii+1)
!        verify if we cover the whole box
         if(ngfftrec(plane)>=ngfft(plane))then
           ngfftrec(plane)=ngfft(plane)
           get_ngfftrec(plane)=0
         end if
!        Exit the loop over ii
         exit
       end if
     end do

   end do

   if (allocated(srch)) then
     ABI_DEALLOCATE(srch)
   end if

!  if(mod(ngfftrec(1),16)/=0) then
!  ngfftrec(1) = ngfftrec(1)+(16-mod(ngfftrec(1),16))
!  ngfftrec(2:3) = ngfftrec(1)
!  endif

   ngfftrec(4)=2*(ngfftrec(1)/2)+1
   ngfftrec(5)=2*(ngfftrec(2)/2)+1
   ngfftrec(6)=ngfftrec(3)

!  --algorithm
   ngfftrec(7)=ngfft(7)   ! to be improved for a better non-parallel algorithm - here it is automatically 401
   ngfftrec(8)=ngfft(8)

 end if

!--For now, recursion method doesn't use paralelism on FFT - which would require a great number of processors 
 nfftrec = product(ngfftrec(1:3))
 ngfftrec(9:11) = (/0,1,0/)   !--(/ paral, nproc, %me \)
 ngfftrec(12:13) = ngfftrec(2:3)   ! n2proc ! n3proc

 tronc  = all(ngfftrec(:3)/=ngfft(:3))
 call timab(602,2,tsec)

 contains

   function dsq_rec(ii,jj,kk,inf_rmet)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dsq_rec'
!End of the abilint section

   real(dp) :: dsq_rec
   integer,intent(in) :: ii,jj,kk
   real(dp),intent(in) :: inf_rmet(3,3)
   dsq_rec=sqrt(&
&   inf_rmet(1,1)*dble(ii**2)&
&   +inf_rmet(2,2)*dble(jj**2)&
&   +inf_rmet(3,3)*dble(kk**2)&
&   +two*(inf_rmet(1,2)*dble(ii*jj)&
&   +inf_rmet(2,3)*dble(jj*kk)&
&   +inf_rmet(3,1)*dble(kk*ii)))
 end function dsq_rec


end subroutine getngrec
!!***
