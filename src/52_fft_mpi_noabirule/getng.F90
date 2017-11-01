!{\src2tex{textfont=tt}}
!!****f* ABINIT/getng
!! NAME
!! getng
!!
!! FUNCTION
!! From ecut and metric tensor in reciprocal space, computes recommended ngfft(1:3)
!! Also computes the recommended value of nfft and mgfft
!! Pay attention that the FFT grid must be compatible with the symmetry operations (see irrzg.f).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! boxcutmin=minimum value of boxcut admitted (boxcut is the ratio
!!  between the radius of the sphere contained in the FFT box, and the
!!  radius of the planewave sphere) : usually 2.0 .
!! ecut=energy cutoff in Hartrees
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! kpt(3)=input k vector in terms of reciprocal lattice primitive translations
!! me_fft=index of the processor in the FFT set (use 0 if sequential)
!! nproc_fft=number of processors in the FFT set (use 1 if sequential)
!! nsym=number of symmetry elements in group
!! paral_fft=0 if no FFT parallelisation ; 1 if FFT parallelisation
!! symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!
!! OUTPUT
!! mgfft= max(ngfft(1),ngfft(2),ngfft(3))
!! nfft=number of points in the FFT box=ngfft(1)*ngfft(2)*ngfft(3)/nproc_fft
!!
!! SIDE EFFECTS
!! Input/Output
!! ngfft(1:18)=integer array with FFT box dimensions and other
!!   information on FFTs. On input ngfft(1:3) contains
!!   optional trial values. If ngfft(1:3)/minbox is greater than value
!!   calculated to avoid wrap-around error and ngfft obeys constraint
!!   placed by the FFT routine that is used
!!   then ngfft(1:3) is left unchanged. Otherwise set to value computed in now.
!!
!! Note that there is the possibility of an undetected error if we
!! are dealing with a cubic unit cell and ngfft(1), ngfft(2) and ngfft(3)
!! are different. In the future we should handle this case.
!!
!! ngfft(4),ngfft(5),ngfft(6)= modified values to avoid cache trashing,
!!        presently : ngfft(4)=ngfft(1)+1 if ngfft(1) is even ;
!!                    ngfft(5)=ngfft(2)+1 if ngfft(2) is even.
!!           in the other cases, ngfft(4:6)=ngfft(1:3).
!!   Other choices may better, but this is left for the future.
!! ngfft(7)=choice for FFT algorithm, see the input variable fftalg
!! ngfft(8)=size of the cache, in bytes (not used here presently).!!
!!   other ngfft slots are used for parallelism see ~abinit/doc/input_variables/vargs.htm#ngfft
!! [unit] = Output Unit number (DEFAULT std_out)
!!
!! PARENTS
!!      fftprof,m_fft,m_fft_prof,memory_eval,mpi_setup,mrgscr,scfcv
!!
!! CHILDREN
!!      bound,print_ngfft,sort_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getng(boxcutmin,ecut,gmet,kpt,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym,paral_fft,symrel,&
&                ngfftc,use_gpu_cuda,unit) ! optional

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_sort

 use defs_fftdata,  only : mg
 use m_fftcore,     only : print_ngfft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getng'
 use interfaces_14_hidewrite
 use interfaces_52_fft_mpi_noabirule, except_this_one => getng
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: me_fft,nproc_fft,nsym,paral_fft
 integer,intent(out) :: mgfft,nfft
 integer,optional,intent(in) :: unit,use_gpu_cuda
 real(dp),intent(in) :: boxcutmin,ecut
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(in),optional :: ngfftc(3)
 integer,intent(inout) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer,save :: first=1,msrch(3),previous_paral_mode=0
 integer :: element,ii,index,isrch,isrch1,isrch2,isrch3,isym,jj,mu,paral_fft_
 integer :: plane,testok,tobechecked,ount,fftalga
 real(dp),parameter :: minbox=0.75_dp
 real(dp) :: dsqmax,dsqmin,ecutmx,prodcurrent,prodtrial,xx,yy
 logical :: testdiv
 character(len=500) :: message
 integer,parameter :: largest_ngfft=mg ! Goedecker FFT: any powers of 2, 3, and 5 - must be coherent with defs_fftdata.F90
 integer,parameter :: maxpow2 =16      ! int(log(largest_ngfft+half)/log(two))
 integer,parameter :: maxpow3 =6       ! int(log(largest_ngfft+half)/log(three))
 integer,parameter :: maxpow5 =6       ! int(log(largest_ngfft+half)/log(five))
 integer,parameter :: maxpow7 =0
 integer,parameter :: maxpow11=0
 integer,parameter :: mmsrch=(maxpow2+1)*(maxpow3+1)*(maxpow5+1)*(maxpow7+1)*(maxpow11+1)
 integer,save :: iperm(mmsrch),srch(mmsrch,3)
 integer(i8b) :: li_srch(mmsrch)
 integer :: divisor(3,3),gbound(3),imax(3),imin(3),ngcurrent(3)
 integer :: ngmax(3),ngsav(3),ngtrial(3)

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

 fftalga = ngfft(7)/100

!If not yet done, compute recommended (boxcut>=2) fft grid dimensions
!In case we switch for paral to sequential mode, recompute srch.
!This is the case e.g. when computing ngfftdiel in sequential mode
!after an initial computation of ngfft in parallel

 paral_fft_=paral_fft;if (nproc_fft==0) paral_fft_=0

 if(first==1.or.paral_fft_ /= previous_paral_mode) then
   first=0; previous_paral_mode=paral_fft_
   srch(:,:)=0

   ! Factors of 2
   srch(1,1)=1
   do ii=1,maxpow2
     srch(ii+1,1)=srch(ii,1)*2
   end do

   ! Factors of 3
   index=maxpow2+1
   if(maxpow3>0)then
     do ii=1,max(1,maxpow3)
       srch(1+ii*index:(ii+1)*index,1)=3*srch(1+(ii-1)*index:ii*index,1)
     end do
   end if

   ! Factors of 5
   index=(maxpow3+1)*index
   if(maxpow5>0)then
     do ii=1,max(1,maxpow5)
       li_srch = 0
       li_srch(1+ii*index:(ii+1)*index)=5_i8b*srch(1+(ii-1)*index:ii*index,1)
       where (li_srch > huge(maxpow3)) li_srch = huge(maxpow3)
       srch(1+ii*index:(ii+1)*index,1)=li_srch(1+ii*index:(ii+1)*index)
     end do
   end if

   ! Factors of 7
   index=(maxpow5+1)*index
   if(maxpow7>0)then
     do ii=1,max(1,maxpow7)
       srch(1+ii*index:(ii+1)*index,1)=7*srch(1+(ii-1)*index:ii*index,1)
     end do
   end if

   ! Factors of 11
   index=(maxpow7+1)*index
   if(maxpow11>0)then
     do ii=1,max(1,maxpow11)
       srch(1+ii*index:(ii+1)*index,1)=11*srch(1+(ii-1)*index:ii*index,1)
     end do
   end if

   call sort_int(mmsrch,srch(:,1),iperm)

   do ii=1,mmsrch
     if(srch(ii,1)>largest_ngfft)exit
   end do
   msrch(1)=ii-1

   ! In case of FFT parallelism, one need ngfft(2) and ngfft(3) to be multiple of nproc_fft
   if(paral_fft_==1)then
     msrch(2)=0
     do ii=1,msrch(1)
       if(modulo(srch(ii,1),nproc_fft)==0) then
         msrch(2)=msrch(2)+1
         srch(msrch(2),2)=srch(ii,1)
       end if
     end do
     !write(message,'(a,i0,a,i0,2a,i0)')&
     ! 'The second and third dimension of the FFT grid: ',ngfft(2),", ",ngfft(3),ch10,&
     ! 'were imposed to be multiple of the number of processors for the FFT: ', nproc_fft
     !if (ount /= dev_null) MSG_COMMENT(message)
   else
     msrch(2)=msrch(1)
     srch(:,2)=srch(:,1)
   end if

   ! The second and third search list have the same constraint
   msrch(3)=msrch(2)
   srch(:,3)=srch(:,2)

!  The set of allowed ngfft values has been found
 end if ! first==1

!Save input values of ngfft
 ngsav(1:3) = ngfft(1:3)

!As an initial guess for ngfft, use the provided coarse mesh grid
 if (PRESENT(ngfftc)) then
   ngfft(1:3)=ngfftc(1:3)
   call wrtout(ount,' Using supplied coarse mesh as initial guess.','COLL')
 else
   ngfft(1:3)=2
 end if

!Enlarge the initial guess until the set of ngfft entirely comprises the sphere
 do

   call bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)

   ! Exit the infinite do-loop if the sphere is inside the FFT box
   if (dsqmin>=(half*boxcutmin**2*ecut/pi**2)) exit

   ! Fix nearest boundary
   do ii=1,msrch(plane)-1
     if (srch(ii,plane)>=ngfft(plane)) then
!      redefine ngfft(plane) to next higher choice
       ngfft(plane)=srch(ii+1,plane)
       exit ! Exit the loop over ii
     end if

     if (ii==msrch(plane)-1)then
       ! Here, we are in trouble
       write(message, '(a,i12,5a)' ) &
&       'ngfft is bigger than allowed value =',ngfft(plane),'.',ch10,&
&       'This indicates that desired ngfft is larger than getng',ch10,&
&       'can handle. The code has to be changed and compiled.'
       MSG_BUG(message)
     end if
   end do

 end do ! End of the infinite do-loop : will either "exit", or abort

!ecutmx=maximum ecut consistent with chosen ngfft
 ecutmx=0.5_dp*pi**2*dsqmin

!Print results
 write(message, '(a,1p,e14.6,a,3i8,a,a,e14.6)' ) &
& ' For input ecut=',ecut,' best grid ngfft=',ngfft(1:3),ch10,&
& '       max ecut=',ecutmx
 call wrtout(ount,message,'COLL')

! The FFT grid is compatible with the symmetries if for each
! symmetry isym, each ii and each jj, the quantity
! (ngfft(jj)*symrel(jj,ii,isym))/ngfft(ii) is an integer.
! This relation is immediately verified for diagonal elements, since
! symrel is an integer. It is also verified if symrel(ii,jj,isym) is zero.

!Compute the biggest (positive) common divisor of each off-diagonal element of the symmetry matrices
 divisor(:,:)=0; tobechecked=0

 do ii=1,3
   do jj=1,3
     if(jj==ii)cycle
     do isym=1,nsym
       if(symrel(jj,ii,isym)==0 .or. divisor(jj,ii)==1 )cycle
       tobechecked=1
       element=abs(symrel(jj,ii,isym))
       testdiv= ( divisor(jj,ii)==0 .or. divisor(jj,ii)==element .or. element==1)
       if(testdiv)then
         divisor(jj,ii)=element
       else
!        Must evaluate common divisor between non-trivial numbers
         do
           if(divisor(jj,ii)<element)element=element-divisor(jj,ii)
           if(divisor(jj,ii)>element)divisor(jj,ii)=divisor(jj,ii)-element
           if(divisor(jj,ii)==element)exit
         end do
       end if
     end do
   end do
 end do

!Check whether there is a problem
 testok=1
 if(tobechecked==1)then
   do ii=1,3
     do jj=1,3
       xx=divisor(jj,ii)*ngfft(jj)
       yy=xx/ngfft(ii)
       if(abs(yy-nint(yy))>tol8)testok=0
     end do
   end do
 end if

!There is definitely a problem
 if(testok==0)then
!  Use a brute force algorithm
!  1) Because one knows that three identical numbers will satisfy
!  the constraint, use the maximal ngfft value to define current triplet
!  and associate total number of grid points
   ngcurrent(1:3)=maxval(ngfft(1:3))
!  Takes into account the fact that ngfft(2) and ngfft(3) must
!  be multiple of nproc_fft
   if(mod(ngcurrent(1),nproc_fft)/=0)ngcurrent(1:3)=ngcurrent(1:3)*max(1,nproc_fft)
   prodcurrent=ngcurrent(1)**3+1.0d-3
!  2) Define maximal values for each component, limited
!  by the maximal value of the list
   ngmax(1)=min(int(prodcurrent/(ngfft(2)*ngfft(3))),srch(msrch(1),1))
   ngmax(2)=min(int(prodcurrent/(ngfft(1)*ngfft(3))),srch(msrch(2),2))
   ngmax(3)=min(int(prodcurrent/(ngfft(1)*ngfft(2))),srch(msrch(3),3))
!  3) Get minimal and maximal search indices
   do ii=1,3
     do isrch=1,msrch(ii)
       index=srch(isrch,ii)
       if(index==ngfft(ii))imin(ii)=isrch
!      One cannot suppose that imax belongs to the allowed list,
!      so must use <= instead of == , to determine largest index
       if(index<=ngmax(ii))imax(ii)=isrch
     end do
   end do
!  4) Compute product of trial ngffts
!  DEBUG
!  write(ount,*)' getng : enter triple loop '
!  write(ount,*)'imin',imin(1:3)
!  write(ount,*)'imax',imax(1:3)
!  write(ount,*)'ngcurrent',ngcurrent(1:3)
!  ENDDEBUG
   do isrch1=imin(1),imax(1)
     ngtrial(1)=srch(isrch1,1)
     do isrch2=imin(2),imax(2)
       ngtrial(2)=srch(isrch2,2)
       do isrch3=imin(3),imax(3)
         ngtrial(3)=srch(isrch3,3)
         prodtrial=real(ngtrial(1))*real(ngtrial(2))*real(ngtrial(3))+1.0d-3
         if(prodtrial>prodcurrent-1.0d-4)exit
!        The trial product is lower or equal to the current product,
!        so now, checks whether the symmetry constraints are OK
         testok=1
         do ii=1,3
           do jj=1,3
             xx=divisor(jj,ii)*ngtrial(jj)
             yy=xx/ngtrial(ii)
             if(abs(yy-nint(yy))>tol8)testok=0
           end do
         end do
!        DEBUG
!        write(ount,'(a,3i6,a,i3,a,es16.6)' )' getng : current trial triplet',ngtrial(1:3),&
!        &     ' testok=',testok,' prodtrial=',prodtrial
!        ENDDEBUG
         if(testok==0)cycle
!        The symmetry constraints are fulfilled, so update current values
         ngcurrent(1:3)=ngtrial(1:3)
         prodcurrent=prodtrial
       end do
     end do
   end do

   ngfft(1:3)=ngcurrent(1:3)
   call bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)
!  ecutmx=maximum ecut consistent with chosen ngfft
   ecutmx=0.5_dp*pi**2*dsqmin
!  Give results
   write(message, '(a,3i8,a,a,e14.6)' ) &
&   ' However, must be changed due to symmetry =>',ngfft(1:3),ch10,&
&   '       with max ecut=',ecutmx
   call wrtout(ount,message,'COLL')

   if (prodcurrent>huge(ii)) then
     write(message, '(5a)' )&
&     'The best FFT grid will lead to indices larger',ch10,&
&     'than the largest representable integer on this machine.',ch10,&
&     'Action: try to deal with smaller problems. Also contact ABINIT group.'
     MSG_ERROR(message)
   end if

 end if ! testok==0

!Possibly use the input values of ngfft
 if ( int( dble(ngsav(1)) / minbox ) >= ngfft(1) .and.&
&     int( dble(ngsav(2)) / minbox ) >= ngfft(2) .and.&
&     int( dble(ngsav(3)) / minbox ) >= ngfft(3) ) then

!  Must check whether the values are in the allowed list
   testok=0
   do mu=1,3
     do ii=1,msrch(mu)
       if(srch(ii,mu)==ngsav(mu))then
         testok=testok+1
         exit
       end if
     end do
   end do
   if(testok==3)then
     write(message,'(a,3(a,i1,a,i3),a)') ' input values of',&
&     (' ngfft(',mu,') =',ngsav(mu),mu=1,3),' are alright and will be used'
     call wrtout(ount,message,'COLL')
     do mu = 1,3
       ngfft(mu) = ngsav(mu)
     end do
   end if

 end if

!mgfft needs to be set to the maximum of ngfft(1),ngfft(2),ngfft(3)
 mgfft = maxval(ngfft(1:3))

 if (paral_fft_==1) then
   ! For the time being, one need ngfft(2) and ngfft(3) to be multiple of nproc_fft
   if(modulo(ngfft(2),nproc_fft)/=0)then
     write(message,'(3a,i5,a,i5)')&
&     'The second dimension of the FFT grid, ngfft(2), should be',&
&     'a multiple of the number of processors for the FFT, nproc_fft.',&
&     'However, ngfft(2)=',ngfft(2),' and nproc_fft=',nproc_fft
     MSG_BUG(message)
   end if
   if(modulo(ngfft(3),nproc_fft)/=0)then
     write(message,'(3a,i5,a,i5)')&
&     'The third dimension of the FFT grid, ngfft(3), should be',&
&     'a multiple of the number of processors for the FFT, nproc_fft.',&
&     'However, ngfft(3)=',ngfft(3),' and nproc_fft=',nproc_fft
     MSG_BUG(message)
   end if

 else if (paral_fft_/=0) then
   write(message,'(a,i0)')'paral_fft_ should be 0 or 1, but its value is ',paral_fft_
   MSG_BUG(message)
 end if

! Compute effective number of FFT points (for this MPI node if parall FFT)
 nfft=product(ngfft(1:3))/max(1,nproc_fft)

!Set up fft array dimensions ngfft(4,5,6) to avoid cache conflicts
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=ngfft(3)
 if (any(fftalga == [FFT_FFTW3, FFT_DFTI])) then
   ! FFTW3 supports leading dimensions but at the price of a larger number of FFTs
   ! to be executed along z when the zero-padded version is used.
   ! One should check whether the augmentation is beneficial for FFTW3.
   ngfft(4)=2*(ngfft(1)/2)+1
   ngfft(5)=2*(ngfft(2)/2)+1
   !ngfft(4)=ngfft(1)
   !ngfft(5)=ngfft(2)
   ngfft(6)=ngfft(3)
 end if

 if (present(use_gpu_cuda)) then
   if (use_gpu_cuda==1) then
     ngfft(4)=ngfft(1)
     ngfft(5)=ngfft(2)
     ngfft(6)=ngfft(3)
   end if
 end if

 ngfft(14:18)=0 ! ngfft(14) to be filled outside of getng

 if (paral_fft_==0) then
   ngfft(9)=0     ! paral_fft_
   ngfft(10)=1    ! nproc_fft
   ngfft(11)=0    ! me_fft
   ngfft(12)=0    ! n2proc
   ngfft(13)=0    ! n3proc
 else
   ngfft(9)=1     ! paral_fft_
   ngfft(10)=nproc_fft
   ngfft(11)=me_fft
   ngfft(12)=ngfft(2)/nproc_fft
   ngfft(13)=ngfft(3)/nproc_fft
 end if


 call print_ngfft(ngfft,header="FFT mesh",unit=ount,mode_paral="COLL")

end subroutine getng
!!***
