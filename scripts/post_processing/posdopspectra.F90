program posdopspectra

!! Script to analyze the electron-positron momentum distributions
!! included in the _DOPPLER file. Calculated 1D distributions in
!! three different directions, normal to the (001), (011) and (111)
!! planes in the reciprocal space.
!! This corresponds to [001], [011] and [111] directions
!! for cubic systems. The script convolutes the data with
!! a Gaussian and interpolates them on an uniform grid with
!! 0.1 mrad spacing. Calculates the S and W parameters 
!! of the Doppler broadening of annihilation radiation.

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!scalars
 integer, parameter:: dp=kind(0.d0) 
 integer :: i1,i2,i3,ifft,ii,ikpt,jj
 integer :: nfft,nkpt,npoints1,npoints2,npoints3
 integer :: nd1,nd2,nd3,nd4,nd5,nd6
 integer :: nspline
 real(dp) :: der1,der2,der3,der4,der5,der6,fwhm,lambda,norm,S_dop
 real(dp) :: zero,one,two,three,two_pi,InvFineStruct
 real(dp) :: vec,vec0,ucvol,W_dop
 logical :: file_exists
 character(len=264) :: filnam
!arrays
 integer,allocatable :: ipoint1(:),ipoint2(:),ipoint3(:)
 integer,allocatable :: weight1(:),weight2(:),weight3(:)
 real(dp) :: a(3),b(3),c(3),dir(3),rprim(3,3),smax,wmin,wmax
 real(dp),allocatable :: der001(:),der011(:),der111(:)
 real(dp),allocatable :: gauss(:),mesh_spline(:)
 real(dp),allocatable :: pcart(:,:,:),rho_moment(:,:),rho_001(:),rho001(:),rho001conv(:)
 real(dp),allocatable :: rho_011(:),rho011(:),rho011conv(:),rho_111(:),rho111(:),rho111conv(:)
 real(dp),allocatable :: spectrum1(:,:),spectrum3(:,:),spectrum2(:,:)
 real(dp),allocatable :: vec001(:),vec011(:),vec111(:)
 character(len=500) :: msg

!******************************************************************
 zero=0._dp;one=1._dp;two=2._dp;three=3._dp
 two_pi=8*atan(1.d0)
 InvFineStruct=137.035999679_dp  ! Inverse of fine structure constant
! Get name of the momentum distribution file
 write(*,*)
 write(*,*) ' What is the name of the 3D electron-positron momentum distribution file?'
 read(*,'(a)')filnam
 write(*,*)
!  Checking the existence of data file
 INQUIRE(FILE=filnam, EXIST=file_exists)   ! file_exists will be TRUE if the file
 if (.not.file_exists) then
   write(*,*) 'Missing data file: '//TRIM(filnam)
   stop
 end if

 write(*,'(a,a,a,i4)')'  posdopspectra : read file ',trim(filnam)
 write(*,*)
 write(*,'(a,a)') 'Opening file ', filnam
 open(unit=19,file=filnam,form='unformatted',status='old')
 read(19) nfft,nkpt,ucvol,rprim
 ALLOCATE(pcart(3,nfft,nkpt))
 ALLOCATE(rho_moment(nfft,nkpt))
 do ikpt=1,nkpt
   read(19) pcart(1:3,1:nfft,ikpt),rho_moment(1:nfft,ikpt)
 end do
 write(*,*) ' => Your momentum distribution file is : ',trim(filnam)
 write(*,*) ' Choose FWHM (in mrad) for the convolution.'
 read(*,*) fwhm

! Calculate vectors normal to three different planes in the reciprocal space - (001), (011) and (111)
! They are equivalent to [001], [011] and [111] directions in the real space
! and to [001], [011] and [111] directions in the reciprocal space for cubic systems

! Vector normal to the (001) plane
 dir=(/zero ,zero ,one/)
 a(:)=matmul(rprim(:,:),dir)
 a(:)=a(:)/(sqrt(dot_product(a,a)))

! Vector normal to the (011) plane
 dir=(/zero ,one ,one/)
 b(:)=matmul(rprim(:,:),dir)
 b(:)=b(:)/(sqrt(dot_product(b,b)))

! Vector normal to the (111) plane
 dir=(/one ,one ,one/)
 c(:)=matmul(rprim(:,:),dir)
 c(:)=c(:)/(sqrt(dot_product(c,c)))


! Calculate the number of different points in the three directions (or equivalent)
! and attribute the lenght of the corresponding vector

 ALLOCATE(vec001 (nfft*nkpt))
 ALLOCATE(vec011 (nfft*nkpt))
 ALLOCATE(vec111 (nfft*nkpt))
 ALLOCATE(ipoint1(nfft*nkpt))
 ALLOCATE(ipoint2(nfft*nkpt))
 ALLOCATE(ipoint3(nfft*nkpt))

 vec001=zero;vec011=zero;vec111=zero
 npoints1=1;npoints2=1;npoints3=1
 vec001(1)=a(1)*pcart(1,1,1)+a(2)*pcart(2,1,1)+a(3)*pcart(3,1,1)
 vec011(1)=b(1)*pcart(1,1,1)+b(2)*pcart(2,1,1)+b(3)*pcart(3,1,1)
 vec111(1)=c(1)*pcart(1,1,1)+c(2)*pcart(2,1,1)+c(3)*pcart(3,1,1)

 ipoint1(1)=1;ipoint2(1)=1;ipoint3(1)=1

 outer1: do jj=2,nkpt*nfft
   ikpt=int((jj-1)/nfft)+1;ifft=mod(jj-1,nfft)+1
   do ii=1,npoints1
     vec0=vec001(ii)-a(1)*pcart(1,ifft,ikpt)-a(2)*pcart(2,ifft,ikpt)-a(3)*pcart(3,ifft,ikpt)
     if (abs(vec0)<0.00001) then
       ipoint1(jj)=ii
       cycle outer1
     end if
   end do
   npoints1 = npoints1 + 1
   ipoint1(jj) = maxval(ipoint1)+1
   vec001 (npoints1) = a(1)*pcart(1,ifft,ikpt)+a(2)*pcart(2,ifft,ikpt)+a(3)*pcart(3,ifft,ikpt)
 end do outer1

 outer2: do jj=2,nkpt*nfft
   ikpt=int((jj-1)/nfft)+1;ifft=mod(jj-1,nfft)+1
   do ii=1,npoints2
     vec0=vec011(ii)-b(1)*pcart(1,ifft,ikpt)-b(2)*pcart(2,ifft,ikpt)-b(3)*pcart(3,ifft,ikpt)
     if (abs(vec0)<0.00001) then
       ipoint2(jj)=ii
       cycle outer2
     end if
   end do
   npoints2 = npoints2 + 1
   ipoint2(jj) = maxval(ipoint2)+1
   vec011 (npoints2) = b(1)*pcart(1,ifft,ikpt)+b(2)*pcart(2,ifft,ikpt)+b(3)*pcart(3,ifft,ikpt)
 end do outer2

 outer3: do jj=2,nkpt*nfft
   ikpt=int((jj-1)/nfft)+1;ifft=mod(jj-1,nfft)+1
   do ii=1,npoints3
     vec0=vec111(ii)-c(1)*pcart(1,ifft,ikpt)-c(2)*pcart(2,ifft,ikpt)-c(3)*pcart(3,ifft,ikpt)
     if (abs(vec0)<0.00001) then
       ipoint3(jj)=ii
       cycle outer3
     end if
   end do
   npoints3 = npoints3 + 1
   ipoint3(jj) = maxval(ipoint3)+1
   vec111 (npoints3) = c(1)*pcart(1,ifft,ikpt)+c(2)*pcart(2,ifft,ikpt)+c(3)*pcart(3,ifft,ikpt)
 end do outer3
! Calculate the weight of each point in the three directions

 ALLOCATE(weight1(npoints1))
 ALLOCATE(weight2(npoints2))
 ALLOCATE(weight3(npoints3))
 do ii=1,npoints1
   weight1(ii)=count((ipoint1(1:nfft*nkpt)==ii))
 end do
 do ii=1,npoints2
   weight2(ii)=count((ipoint2(1:nfft*nkpt)==ii))
 end do
 do ii=1,npoints3
   weight3(ii)=count((ipoint3(1:nfft*nkpt)==ii))
 end do

! Transform vector lenght to mrad

 vec001(1:npoints1)=vec001(1:npoints1)*two_pi*1000_dp/InvFineStruct
 vec011(1:npoints2)=vec011(1:npoints2)*two_pi*1000_dp/InvFineStruct
 vec111(1:npoints3)=vec111(1:npoints3)*two_pi*1000_dp/InvFineStruct
 ALLOCATE(rho001(npoints1))
 ALLOCATE(rho011(npoints2))
 ALLOCATE(rho111(npoints3))

 do jj=1,nkpt*nfft
   ikpt=int((jj-1)/nfft)+1;ifft=mod(jj-1,nfft)+1
   i1=ipoint1(jj);i2=ipoint2(jj);i3=ipoint3(jj)
   rho001(i1)=rho001(i1)+rho_moment(ifft,ikpt)/weight1(i1)
   rho011(i2)=rho011(i2)+rho_moment(ifft,ikpt)/weight2(i2)
   rho111(i3)=rho111(i3)+rho_moment(ifft,ikpt)/weight3(i3)
 end do

 DEALLOCATE(weight1)
 DEALLOCATE(weight2)
 DEALLOCATE(weight3)
 DEALLOCATE(ipoint1)
 DEALLOCATE(ipoint2)
 DEALLOCATE(ipoint3)


 DEALLOCATE(pcart)
 DEALLOCATE(rho_moment)

! Sort the vectors in ascending order

 nd1=-count(vec001(1:npoints1)<0);nd2=npoints1+nd1-1
 nd3=-count(vec011(1:npoints2)<0);nd4=npoints2+nd3-1
 nd5=-count(vec111(1:npoints3)<0);nd6=npoints3+nd5-1
 ALLOCATE(spectrum1(2,nd1:nd2))
 ALLOCATE(spectrum2(2,nd3:nd4))
 ALLOCATE(spectrum3(2,nd5:nd6))

 do ii=1,npoints1
   i1=count(vec001(ii)>vec001(1:npoints1))+nd1
   spectrum1(1,i1)=vec001(ii)
   spectrum1(2,i1)=rho001(ii)
 end do
 do ii=1,npoints2
   i1=count(vec011(ii)>vec011(1:npoints2))+nd3
   spectrum2(1,i1)=vec011(ii)
   spectrum2(2,i1)=rho011(ii)
 end do
 do ii=1,npoints3
   i1=count(vec111(ii)>vec111(1:npoints3))+nd5
   spectrum3(1,i1)=vec111(ii)
   spectrum3(2,i1)=rho111(ii)
 end do

 DEALLOCATE(vec001)
 DEALLOCATE(rho001)
 DEALLOCATE(vec011)
 DEALLOCATE(vec111)
 DEALLOCATE(rho011)
 DEALLOCATE(rho111)
! Convolve the 001 spectrum with a Gaussian 

 ALLOCATE(rho001conv(nd1:nd2))
 ALLOCATE(rho011conv(nd3:nd4))
 ALLOCATE(rho111conv(nd5:nd6))
 rho001conv=zero; rho011conv=zero;rho111conv=zero
 call conv( nd1, nd2, spectrum1, fwhm, rho001conv)
 call conv( nd3, nd4, spectrum2, fwhm, rho011conv)
 call conv( nd5, nd6, spectrum3, fwhm, rho111conv)

! Write the raw and convoluted spectra to the doppler_out file

 open(unit=114, file='doppler_out', status='replace')
 write(114,*)'Raw momentum distribution in [001] direction:'
 write(114,*)
 write(114,*)'  p001 (mrad):             Probalility:'
 do ii=0,nd2
   write(114,*) spectrum1(1,ii),spectrum1(2,ii)
 end do
 write(114,*)
 write(114,*)'Convoluted momentum distribution in [001] direction:'
 write(114,*)'FWHM (in mrad) =',fwhm
 write(114,*)
 write(114,*)'  p001 (mrad):             Probalility:'
 do ii=0,nd2
   write(114,*) spectrum1(1,ii),rho001conv(ii)
 end do
 write(114,*)'Raw momentum distribution in [011] direction:'
 write(114,*)
 write(114,*)'  p011 (mrad):             Probalility:'
 do ii=0,nd4
   write(114,*) spectrum2(1,ii),spectrum2(2,ii)
 end do
 write(114,*)
 write(114,*)'Convoluted momentum distribution in [011] direction:'
 write(114,*)'FWHM (in mrad) =',fwhm
 write(114,*)
 write(114,*)'  p011 (mrad):             Probalility:'
 do ii=0,nd4
   write(114,*) spectrum2(1,ii),rho011conv(ii)
 end do
 write(114,*)'Raw momentum distribution in [111] direction:'
 write(114,*)
 write(114,*)'  p111 (mrad):             Probalility:'
 do ii=0,nd6
   write(114,*) spectrum3(1,ii),spectrum3(2,ii)
 end do
 write(114,*)
 write(114,*)'Convoluted momentum distribution in [111] direction:'
 write(114,*)'FWHM (in mrad) =',fwhm
 write(114,*)
 write(114,*)'  p111 (mrad):             Probalility:'
 do ii=0,nd6
   write(114,*) spectrum3(1,ii),rho111conv(ii)
 end do

! Interpolate the spectra

 der1=rho001conv(nd1);der2=rho001conv(nd2)
 der3=rho011conv(nd3);der4=rho011conv(nd4)
 der5=rho111conv(nd5);der6=rho111conv(nd6)
 ALLOCATE(der001(nd1:nd2))
 ALLOCATE(der011(nd3:nd4))
 ALLOCATE(der111(nd5:nd6))
 call spline(spectrum1(1,:),rho001conv,nd2-nd1,der1,der2,der001)
 call spline(spectrum2(1,:),rho011conv,nd4-nd3,der3,der4,der011)
 call spline(spectrum3(1,:),rho111conv,nd6-nd5,der5,der6,der111)
 nspline=1001
 ALLOCATE(mesh_spline(nspline))
 ALLOCATE(rho_001(nspline))
 ALLOCATE(rho_011(nspline))
 ALLOCATE(rho_111(nspline))

 do ii=1,nspline
   mesh_spline(ii)=0.1_dp*(ii-1)
 end do
 call splint(nd2-nd1,spectrum1(1,:),rho001conv,der001,nspline,mesh_spline,rho_001)
 call splint(nd4-nd3,spectrum2(1,:),rho011conv,der011,nspline,mesh_spline,rho_011)
 call splint(nd6-nd5,spectrum3(1,:),rho111conv,der111,nspline,mesh_spline,rho_111)

 DEALLOCATE(rho001conv)
 DEALLOCATE(rho011conv)
 DEALLOCATE(rho111conv)

 DEALLOCATE(der001)
 DEALLOCATE(der011)
 DEALLOCATE(der111)

 DEALLOCATE(spectrum1)
 DEALLOCATE(spectrum2)
 DEALLOCATE(spectrum3)


! Write the interpolated and normalized 001 spectrum to the rho_001 file

 lambda=sum(rho_001(:))*(mesh_spline(2)-mesh_spline(1))
 rho_001(:)=rho_001(:)/lambda
 lambda=sum(rho_011(:))*(mesh_spline(2)-mesh_spline(1))
 rho_011(:)=rho_011(:)/lambda
 lambda=sum(rho_111(:))*(mesh_spline(2)-mesh_spline(1))
 rho_111(:)=rho_111(:)/lambda

 open(unit=111, file='rho_001', status='replace')
 do ii=1,nspline
   write(111,'(f10.2,es24.15)') mesh_spline(ii),rho_001(ii)
 end do
 close(111)

 open(unit=112, file='rho_111', status='replace')
 do ii=1,nspline
   write(112,'(f10.2,es24.15)') mesh_spline(ii),rho_111(ii)
 end do
 close(112)

 open(unit=113, file='rho_011', status='replace')
 do ii=1,nspline
   write(113,'(f10.2,es24.15)') mesh_spline(ii),rho_011(ii)
 end do
 close(113)

 write(*,*)
 write(*,*) 'Convoluted and normalized spectra on a uniform grid&
& have been written to rho_001, rho_011 and rho_111 files.' 
 write(*,*)

! Calculate the S and W parameters using the 001 projection
 write(*,*) ' Give the upper limit  for S calculations (in mrad)'
 read(*,*) smax
 write(*,*) ' Give ranges for W calculations (in mrad)'
 read(*,*) wmin,wmax
 S_dop=zero;W_dop=zero
 do ii=2,nint(smax*10+1)
   S_dop=S_dop+(rho_001(ii)+rho_011(ii)+rho_111(ii))/three
 end do
 S_dop=(S_dop+rho_001(1)/two)*0.1_dp
 do ii=nint(wmin*10+1),nint(wmax*10+1)
   W_dop=W_dop+(rho_001(ii)+rho_011(ii)+rho_111(ii))/three
 end do
 W_dop=W_dop*0.1_dp
 write(*,*) 'S parameter:', S_dop
 write(*,*) 'W parameter:', W_dop

 DEALLOCATE(mesh_spline)
 DEALLOCATE(rho_001)
 DEALLOCATE(rho_011)
 DEALLOCATE(rho_111)

!Write results
 write(114,*)'S parameter:',S_dop
 write(114,*)'W parameter:',W_dop
 close(114)

 write(*,*)
 write(*,*) 'S and W parameters and momentum distribution spectra have been written to the doppler_out file.' 
 write(*,*)

end program posdopspectra

!!***
!!***

!! NAME
!!  conv
!!
!! FUNCTION
!!  Computes convolution with a gaussian with fwhm.


subroutine conv( ng1, ng2, spectrum, fwhm, rhoconv)


  implicit none

  integer, parameter:: dp=kind(0.d0) 
  integer, intent(in) :: ng1
  integer, intent(in) :: ng2
  real(dp), intent(in) :: fwhm
  real(dp), intent(in) :: spectrum(2,ng1:ng2)
  real(dp), intent(out) :: rhoconv(ng1:ng2)

  integer :: i1,i2,iz
  real(dp) :: two,zero
  real(dp), allocatable :: gauss(:)

  zero=0._dp;two=2._dp

 ALLOCATE(gauss(ng1:ng2))
 gauss=zero

 do iz=ng1,ng2
   gauss(iz)=exp(-dble(spectrum(1,iz)*spectrum(1,iz))/(two*((fwhm/2.35482_dp)**2)))
 end do

 do iz=ng1,0
   rhoconv(iz) = zero
   i1=ng1
   i2=iz-ng1
   do while (i2 >= ng1)
     if (i2<=ng2) then
       rhoconv(iz) = rhoconv(iz) + spectrum(2,i1)*gauss(i2)
     end if
     i1 = i1+1
     i2 = i2-1
   end do
 end do
 do iz=1,ng2
   rhoconv(iz) = zero
   i1=iz-ng2
   i2=ng2
   do while (i1 <= ng2)
     if (i1>=ng1) then
       rhoconv(iz) = rhoconv(iz) + spectrum(2,i1)*gauss(i2)
     end if
     i1 = i1+1
     i2 = i2-1
   end do
 end do
 DEALLOCATE(gauss)

end subroutine conv
!!***
!!***


!   subroutines spline and splint copied from the ABINIT source
!   from the m_splines module
!   (/src/28_numeric_noabirule/m_splines.F90)

subroutine spline( t, y, n, ybcbeg, ybcend, ypp )
!
!  Author:
!
!    John Burkardt
!    (XGonze got it from http://www.psc.edu/~burkardt/src/spline/spline.html)
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2. 
!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!    spline will actually be linear. 
!
!    Input, double precision T(N), the knot values, that is, the points where data
!    is specified.  The knot values should be distinct, and increasing.
!
!    Input, double precision Y(N), the data values to be interpolated.
!
!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!
!    Output, double precision YPP(N), the second derivatives of the cubic spline.


  implicit none

  integer, parameter:: dp=kind(0.d0) 
  integer, intent(in) :: n
  real(dp), intent(in) :: t(n)
  real(dp), intent(in) :: y(n)
  real(dp), intent(in) :: ybcbeg
  real(dp), intent(in) :: ybcend

  real(dp), intent(out) :: ypp(n)

  integer :: ibcbeg
  integer :: ibcend
  integer :: i,k
  real(dp) :: ratio,pinv
  real(dp), allocatable :: tmp(:)

  ALLOCATE(tmp(n))

!
!  XG041127
  ibcbeg=1 ; ibcend=1
  if(ybcbeg>1.0d+30)ibcbeg=0
  if(ybcend>1.0d+30)ibcend=0
!
!  Set the first and last equations.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.d0
    tmp(1) = 0.d0
  else if ( ibcbeg == 1 ) then
    ypp(1) = -0.5d0
    tmp(1) = (3.d0/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1))-ybcbeg)
  end if
  if ( ibcend == 0 ) then
    ypp(n) = 0.d0
    tmp(n) = 0.d0
  else if ( ibcend == 1 ) then
    ypp(n) = 0.5d0
    tmp(n) = (3.d0/(t(n)-t(n-1)))*(ybcend-(y(n)-y(n-1))/(t(n)-t(n-1)))
  end if

!
!  Set the intermediate equations.
!
  do i=2,n-1
   ratio=(t(i)-t(i-1))/(t(i+1)-t(i-1))
   pinv = 1.0d0/(ratio*ypp(i-1) + 2.0d0)
   ypp(i) = (ratio-1.0d0)*pinv
   tmp(i)=(6.0d0*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1)) &
&    /(t(i)-t(i-1)))/(t(i+1)-t(i-1))-ratio*tmp(i-1))*pinv
   if (abs(tmp(i))<1.d5*tiny(0.d0)) tmp(i)=0.d0   !MT20050927
  enddo

! Solve the equations
  ypp(n) = (tmp(n)-ypp(n)*tmp(n-1))/(ypp(n)*ypp(n-1)+1.0d0)
  do k=n-1,1,-1
   ypp(k)=ypp(k)*ypp(k+1)+tmp(k)
  enddo

  DEALLOCATE(tmp)

  return
end subroutine spline
!!***
!!***

!----------------------------------------------------------------------

!! NAME
!!  splint
!!
!!  Compute spline interpolation. There is no hypothesis
!!  about the spacing of the input grid points.
!!
!! INPUTS
!!  nspline: number of grid points of input mesh
!!  xspline(nspline): input mesh
!!  yspline(nspline): function on input mesh
!!  ysplin2(nspline): second derivative of yspline on input mesh
!!  nfit: number of points of output mesh
!!  xfit(nfit): output mesh
!!
!! OUTPUT
!!  yfit(nfit): function on output mesh
!!  [ierr]=A non-zero value is used to signal that some points in xfit exceed xspline(nspline).
!!    The input value is incremented by the number of such points.

subroutine splint(nsplin,xspline,yspline,ysplin2,nfit,xfit,yfit)

 implicit none

 integer, parameter:: dp=kind(0.d0) 
 integer, intent(in) :: nfit, nsplin
 real(dp), intent(in) :: xspline(nsplin)
 real(dp), intent(in) :: yspline(nsplin)
 real(dp), intent(in) :: ysplin2(nsplin)
 real(dp), intent(in) :: xfit(nfit)

 real(dp), intent(out) :: yfit(nfit)

!local
 integer :: left,i,k,right,my_err
 real(dp) :: delarg,invdelarg,aa,bb

!source
 my_err=0

 left = 1
 do i=1, nfit
   yfit(i)=0.d0  ! Initialize for the unlikely event that rmax exceed r(mesh)
   !
   do k=left+1, nsplin

     if(xspline(k) >= xfit(i)) then
       if(xspline(k-1) <= xfit(i)) then
         right = k
         left = k-1
       end if
       delarg= xspline(right) - xspline(left)
       invdelarg= 1.0d0/delarg
       aa= (xspline(right)-xfit(i))*invdelarg
       bb= (xfit(i)-xspline(left))*invdelarg

       yfit(i) = aa*yspline(left) + bb*yspline(right)    &
&               +( (aa*aa*aa-aa)*ysplin2(left) +         &
&                  (bb*bb*bb-bb)*ysplin2(right) ) *delarg*delarg/6.0d0
       exit
     end if
   end do ! k
   !

   if (k==nsplin+1) my_err=my_err+1 ! xfit not found 
 end do ! i

end subroutine splint
