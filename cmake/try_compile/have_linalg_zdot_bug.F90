program have_linalg_zdot_bug

  implicit none

  real :: x(2)=[3.,4.], y(2)=[1.,1.]
  complex :: w(2)=[(4.,3.),(3.,4.)], z(2)=[(5.,6.),(7.,8.)]
  double complex :: zw(2)=[(4.d0,3.d0),(3.d0,4.d0)]
  double complex :: zz(2)=[(5.d0,6.d0),(7.d0,8.d0)]
  real, external :: sdot, sdsdot, snrm2, scasum
  complex, external :: cdotu, cdotc
  double complex, external :: zdotu, zdotc
  logical :: blas_ok =.true.

  if (abs(sdot(2,x,1,y,1)-real(7.0))>0.00001) blas_ok=.false.
  if (abs(sdsdot(2,0.0,x,1,y,1)-real(7.0))>0.00001) blas_ok=.false.
  if (abs(snrm2(2,x,1)-real(5.0))>0.00001) blas_ok=.false.
  if (abs(scasum(2,w,1)-real(14.0))>0.00001) blas_ok=.false.
  if (abs(real(cdotu(2,w,1,z,1))-real(-9.0))>0.00001) blas_ok=.false.
  if (abs(aimag(cdotu(2,w,1,z,1))-real(91.0))>0.00001) blas_ok=.false.
  if (abs(real(cdotc(2,w,1,z,1))-real(91.0))>0.00001) blas_ok=.false.
  if (abs(aimag(cdotc(2,w,1,z,1))-real(5.0))>0.00001) blas_ok=.false.
  if (abs(dble(zdotu(2,zw,1,zz,1))-dble(-9.0))>0.00001d0) blas_ok=.false.
  if (abs(aimag(zdotu(2,zw,1,zz,1))-dble(91.0))>0.00001d0) blas_ok=.false.
  if (abs(dble(zdotc(2,zw,1,zz,1))-dble(91.0))>0.00001d0) blas_ok=.false.
  if (abs(aimag(zdotc(2,zw,1,zz,1))-dble(5.0))>0.00001d0) blas_ok=.false.

  if (blas_ok) then
    print *,"OK"
  else
    print *,"FAILED"
  end if

end program have_linalg_zdot_bug
