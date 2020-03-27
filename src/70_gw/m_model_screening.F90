!!****m* ABINIT/m_model_screening
!! NAME
!! m_model_screening
!!
!! FUNCTION
!!  Module containing functions for calculating and fitting model dielectric functions
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MS)
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

MODULE m_model_screening

 use defs_basis
 use m_errors
 use m_abicore

 use m_io_tools,       only : open_file

 implicit none

 private

 public :: im_screening         ! Calc. Drude-Lorentz model function from parameters.
 public :: re_screening         ! Calc. Drude-Lorentz model function from parameters.
 public :: re_and_im_screening  ! Calc. Drude-Lorentz model function from parameters.
 public :: re_and_im_screening_with_phase  ! Calc. Drude-Lorentz model function from parameters.
                                           !  with the addition of a phase
 public :: sequential_fitting   ! Fit poles one by one
 public :: init_peaks_from_grid ! find approximate expression for parameters from
                                !  chi0 or eps^-1 on a grid in the complex plane.
 public :: init_peaks_even_dist ! Initial guess from even distributuin of peaks
 public :: init_single_peak     ! Initialise a single peak from the maximum, the
                                !  origin, and the second value along the
                                !  imaginary axis
! public :: int_screening        ! Find the integral along real or complex axis
                                !  from parameters.
 public :: remove_phase

CONTAINS  !==============================================================================
!!***

!!****f* m_model_screening/im_screening
!! NAME
!!  im_screening
!!
!! FUNCTION
!!  Return the imaginary part of model dielectric / inverse dielectric
!!  function as given by a Drude-Lorentx model
!!
!!  The function is a sum in the complex plane:
!!      f(z) = Sum_n  f_n * Im[ ((w_n^2-z^2) - i*gamma*z)^-1 ], z=a-i*b
!!
!!  Here, f_n is the oscillator strength, w_n the location of the peak for
!!   the imaginary function, and gamma is related to the width
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  coeff      = The coefficients in order: f_1,w_1,gamma_1,f_2,w_2,gamma_2,
!!                                          ...,f_n,w_n,gamma_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine im_screening(omega,fval,nomega,coeff,ncoeff)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,ncoeff
!arrays
  complex(dpc),intent(in)  :: omega(nomega)
  real(gwp)  ,intent(in)  :: coeff(ncoeff)
  real(gwp)  ,intent(out) :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: io,ip
  real(gwp) :: rez,imz,realp,imagp
  real(gwp) :: fn,wn,gamman
! *********************************************************************

! The expression is: -f_n*(2*rez*imz-rez*gamma_n)
!    /( (-imz*gamma_n+w_n^2+imz^2-rez^2)^2 + (2*rez*imz-rez*gamma_n)^2 )

  do io=1,nomega
    fval(io) = 0.0
    rez =  REAL(omega(io))
    imz = AIMAG(omega(io))
    do ip=1,ncoeff,3
      fn     = coeff(ip)
      wn     = coeff(ip+1)
      gamman = coeff(ip+2)
      realp  = -imz*gamman+wn*wn+imz*imz-rez*rez
      imagp  = rez*(two*imz-gamman)

        fval(io) = fval(io)-fn*imagp/((realp*realp)+(imagp*imagp))

    end do
  end do

end subroutine im_screening
!!***

!!****f* m_model_screening/re_screening
!! NAME
!!  re_screening
!!
!! FUNCTION
!!  Return the real part of model dielectric / inverse dielectric
!!  function as evaluated from pole coefficients.
!!
!!  The function is a sum of poles in the complex plane:
!!      f(z) = Sum_n[ A/(z-(B-iC)) - A/(z-(-B+iC)) ],
!!   where each pole occurs twice in a time-ordered fashion.
!!
!!  Here, the A are the oscillator strengths, B the real component of the position
!!  of the pole, and C the imaginary component.
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine re_screening(omega,fval,nomega,coeff,ncoeff)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,ncoeff
!arrays
  complex(dpc),intent(in)  :: omega(nomega)
  real(gwp)  ,intent(in)  :: coeff(ncoeff)
  real(gwp)  ,intent(out) :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: io,ip
  real(gwp) :: rez,imz,realp,imagp
  real(gwp) :: fn,wn,gamman
! *********************************************************************

! The expression is: fn*(-imz*gamma+w_n^2+imz^2-rez^2)
!    /( (-imz*gamma+w_n^2+imz^2-rez^2)^2 + (2*rez*imz-rez*gamma)^2 )

  do io=1,nomega
    fval(io) = 0.0
    rez =  REAL(omega(io))
    imz = AIMAG(omega(io))
    do ip=1,ncoeff,3
      fn     = coeff(ip)
      wn     = coeff(ip+1)
      gamman = coeff(ip+2)
      realp  = -imz*gamman+wn*wn+imz*imz-rez*rez
      imagp  = rez*(two*imz-gamman)

        fval(io) = fval(io)-fn*realp/((realp*realp)+(imagp*imagp))

    end do
  end do

end subroutine re_screening
!!***

!!****f* m_model_screening/re_and_im_screening
!! NAME
!!  re_and_im_screening
!!
!! FUNCTION
!!  Return the real and imaginary part of model dielectric / inverse dielectric
!!  function as evaluated from pole coefficients.
!!
!!  The function is a sum of poles in the complex plane:
!!      f(z) = Sum_n[ A/(z-(B-iC)) - A/(z-(-B+iC)) ],
!!   where each pole occurs twice in a time-ordered fashion.
!!
!!  Here, the A are the oscillator strengths, B the real component of the position
!!  of the pole, and C the imaginary component.
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_model_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine re_and_im_screening(omega,fval,nomega,coeff,ncoeff)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,ncoeff
!arrays
  complex(dpc) ,intent(in)  :: omega(nomega)
  real(gwp)   ,intent(in)  :: coeff(ncoeff)
  complex(gwp),intent(out) :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: io,ip
  real(gwp) :: rez,imz,realp,imagp,refval,imfval
  real(gwp) :: fn,wn,gamman
! *********************************************************************

! The expression is: fn*(-imz*gamma+w_n^2+imz^2-rez^2)
!    /( (-imz*gamma+w_n^2+imz^2-rez^2)^2 + (2*rez*imz-rez*gamma)^2 )

  do io=1,nomega
    fval(io) = 0.0
    rez =  REAL(omega(io))
    imz = AIMAG(omega(io))
    do ip=1,ncoeff,3
      fn     = coeff(ip)
      wn     = coeff(ip+1)
      gamman = coeff(ip+2)
      realp  = -imz*gamman+wn*wn+imz*imz-rez*rez
      imagp  = rez*(two*imz-gamman)

        refval   = fn*realp/((realp*realp)+(imagp*imagp))
        imfval   = fn*imagp/((realp*realp)+(imagp*imagp))

        fval(io) = fval(io)-CMPLX(refval,imfval)

    end do
  end do

end subroutine re_and_im_screening
!!***

!!****f* m_model_screening/re_and_im_screening_with_phase
!! NAME
!!  re_and_im_screening_with_phase
!!
!! FUNCTION
!!  Return the real and imaginary part of model dielectric / inverse dielectric
!!  function as evaluated from pole coefficients.
!!
!!  The function is a sum of poles in the complex plane:
!!      f(z) = Sum_n[ A/(z-(B-iC)) - A/(z-(-B+iC)) ],
!!   where each pole occurs twice in a time-ordered fashion.
!!
!!  Here, the A are the oscillator strengths, B the real component of the position
!!  of the pole, and C the imaginary component.
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      calc_sigc_pole_cd
!!
!! CHILDREN
!!
!! SOURCE

subroutine re_and_im_screening_with_phase(omega,fval,nomega,coeff,ncoeff)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,ncoeff
!arrays
  complex(dpc) ,intent(in)  :: omega(nomega)
  real(gwp)   ,intent(in)  :: coeff(ncoeff)
  complex(gwpc),intent(out) :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: io,ip,npoles
  real(gwp) :: rez,imz,realp,imagp,refval,imfval,retemp,imtemp
  real(gwp) :: fn,wn,gamman,imrot,rerot
! *********************************************************************

! The expression is: fn*(-imz*gamma+w_n^2+imz^2-rez^2)
!    /( (-imz*gamma+w_n^2+imz^2-rez^2)^2 + (2*rez*imz-rez*gamma)^2 )
  npoles = (ncoeff-1)/3

  do io=1,nomega
    fval(io) = 0.0
    rez =  REAL(omega(io))
    imz = AIMAG(omega(io))
    do ip=1,(ncoeff-1),3
      fn     = coeff(ip)
      wn     = coeff(ip+1)
      gamman = coeff(ip+2)
      realp  = -imz*gamman+wn*wn+imz*imz-rez*rez
      imagp  = rez*(two*imz-gamman)

        refval   = fn*realp/((realp*realp)+(imagp*imagp))
        imfval   = fn*imagp/((realp*realp)+(imagp*imagp))

        fval(io) = fval(io)-CMPLX(refval,imfval)

    end do
    ! Restore phase
    rerot = COS(coeff(npoles*3+1))
    imrot = SIN(coeff(npoles*3+1))
    retemp = REAL(fval(io))
    imtemp = AIMAG(fval(io))
    fval(io) = CMPLX(rerot*retemp-imrot*imtemp,rerot*imtemp + imrot*retemp)
  end do

end subroutine re_and_im_screening_with_phase
!!***


!!****f* m_model_screening/sequential_fitting
!! NAME
!!  sequential_fitting
!!
!! FUNCTION
!!  Fit a function in the complex plane pole-by-pole in such
!!  a way as to increasingly minimise the error
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  refval     = Real part of function to be fitted
!!  imfval     = imaginary part of function to be fitted
!!  nomega     = Total number of points in the complex plane
!!  nfreqre    = Number of points along real axis
!!  nfreqim    = Number of points along imaginary axis
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  ncoeff     = number of coefficients
!!  prtvol     = Diagnostics verbose level
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine sequential_fitting(omega,refval,imfval,nomega,nfreqre,coeff,&
& ncoeff,prtvol,startcoeff)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,nfreqre,ncoeff,prtvol
!arrays
  complex(dpc),intent(in)     :: omega(nomega)
  real(gwp)  ,intent(out)    :: coeff(ncoeff)
  real(gwp)  ,intent(inout)  :: refval(nomega),imfval(nomega)
  real(gwp),optional,intent(out) :: startcoeff(ncoeff)

!Local variables-------------------------------
!scalars
  integer :: ip,npoles,idx
  real(gwp) :: thiscoeff(3),norm,invnorm
  real(dpc)  :: re_zvals(nomega),im_zvals(nomega)
!  real(dpc)  :: orig_refval(nomega),orig_imfval(nomega)
  complex(gwp) :: pole_func(nomega)
! *********************************************************************

  npoles = ncoeff/3
  re_zvals(:) = REAL(omega(:))
  im_zvals(:) = AIMAG(omega(:))

  ! Normalise
  norm    = MAXVAL(ABS(imfval))
  invnorm = 1.0_gwp/norm
  refval  = invnorm*refval
  imfval  = invnorm*imfval

  ! Loop over poles to fit
  do ip=1,npoles
    idx = 3*(ip-1)+1
    ! Initialise pole
    call init_single_peak(omega,refval,imfval,nomega,nfreqre,thiscoeff,prtvol)
    if (present(startcoeff)) then
      startcoeff(idx:idx+2) = thiscoeff(1:3)
    end if
    ! Make fit
#ifdef HAVE_LEVMAR
    call dfit_re_and_im_screening(re_zvals,im_zvals,imfval,refval,&
&    nomega,3,thiscoeff,prtvol)
#else
    MSG_ERROR(' ABINIT was not compiled with the levmar library!')
#endif
    ! Remove current fit
    call re_and_im_screening(omega,pole_func,nomega,thiscoeff,3)
    refval(:) = refval(:) -  REAL(pole_func(:))
    imfval(:) = imfval(:) - AIMAG(pole_func(:))
    coeff(idx:idx+2) = thiscoeff(1:3)
    coeff(idx) = norm*coeff(idx)
  end do

end subroutine sequential_fitting
!!***

!!****f* m_model_screening/init_peaks_from_grid
!! NAME
!!  init_peaks_from_grid
!!
!! FUNCTION
!!
!!  Find an initial guess of coefficents from the "valleys" and "hills" in
!!   the complex plane.
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  yvals      = The function to be fitted in the complex plane
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!  prtvol   = Verbosity of diagnostics
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_peaks_from_grid(omega,fval,nomega,nfreqre,nfreqim,coeff,ncoeff,prtvol)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,nfreqre,nfreqim,ncoeff,prtvol
!arrays
  complex(dpc) ,intent(in)  :: omega(nomega)
  real(gwp)   ,intent(out) :: coeff(ncoeff)
  complex(gwpc),intent(in)  :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: npoles,iline,idx,ip
  real(gwp) :: pol,gam,maxv,df,dk,val2,b2,osc,val1
  real(gwp) :: temp1,temp2,temp3

!arrays
  integer :: ploc(ncoeff/3)

! *********************************************************************

  npoles = ncoeff/3

! Map the evolution of peaks if prtvol>10
  if (prtvol>10) then
    call print_peaks(omega,fval,nomega,nfreqre,nfreqim)
  end if ! prtvol>10

! Count the number of peaks per line and find the location of the
! constant-imaginary frequency line wich has at least a number of
! peaks commensurate with the requested number of poles
  call find_peaks(fval,nomega,nfreqre,nfreqim,ploc,npoles,iline)
  write(std_out,*) ' Optimum peak locations:',ploc
  write(std_out,*) '               on iline:',iline
! Now fit the peaks. A linear interpolation along the imaginary
!  direction is used to get a rough estimate of the width of
!  the peak.
  do ip=1,npoles
    pol  = REAL(omega(ploc(ip)))
    maxv = AIMAG(fval(ploc(ip)))
    write(std_out,*) ' maxv:',maxv
    if (ploc(ip)<nfreqre+1) then ! We are right on the real axis
      if (ploc(ip)==1) then ! Peak is at origin (Drude peak, i.e. metal)
        b2   = AIMAG(omega(nfreqre+1))
        val2 = AIMAG(fval(nfreqre+1))
        write(std_out,*) '1: ploc:',ploc(ip),' b2:',b2,' val2:',val2
      else ! Second value will be in c-plane
        idx  = nfreqre+nfreqim+ploc(ip)-1
        b2   = AIMAG(omega(idx))
        val2 = AIMAG(fval(idx))
        write(std_out,*) '2: ploc:',ploc(ip),' b2:',b2,' val2:',val2
      end if
    else if (ploc(ip)<nfreqre+nfreqim+1) then ! We are right on the imaginary axis
      if (ploc(ip)==nfreqre+nfreqim) then
        MSG_ERROR(' Peak in upper left corner. This should never happen')
      end if
      b2   = AIMAG(omega(ploc(ip)+1))
      val2 = AIMAG(fval(ploc(ip)+1))
        write(std_out,*) '3: ploc:',ploc(ip),' b2:',b2,' val2:',val2
    else ! We are in the complex plane
      idx  = ploc(ip)+nfreqre-1
      b2   = AIMAG(omega(idx))
      val2 = AIMAG(fval(idx))
      write(std_out,*) '4: ploc:',ploc(ip),' idx:',idx,' b2:',b2,' val2:',val2
    end if
    df = ABS(val2 - maxv)
    dk = df/b2
    gam = -ABS(val2/dk)
    !temp1 = SQRT(-b2*b2*val2*(val2+maxv)+(maxv*maxv)**2)
    !temp2 = b2*(two*pol*pol+b2*b2)*val2-b2*maxv*pol*pol
    !temp3 = (pol*pol+b2*b2)*val2+maxv*pol*pol
    !gam = -((temp1-temp2)/temp3)
    if (gam>zero) gam = ((temp1+temp2)/temp3)
    osc = maxv*gam*pol
    idx = 3*(ip-1)+1
    coeff(idx  ) = osc ! Oscillator strength
    coeff(idx+1) = pol ! Position of maximum
    coeff(idx+2) = gam ! Spread of function
    if (prtvol>9) then
      write(std_out,'(a,a,i0)') ch10,' Pole no,: ',ip
      write(std_out,'(a,ES16.8)')    '  Osc. strength:',osc
      write(std_out,'(a,ES16.8,a)')  '  Peak location:',pol*Ha_eV,' eV'
      write(std_out,'(a,ES16.8,a)')  '     Peak width:',gam*Ha_eV,' eV'
      val2 = gam*half
      val1 = SIGN(1.0_gwp,pol*pol-val2*val2)*SQRT(ABS(pol*pol-val2*val2))
      write(std_out,'(a,ES16.8,a)')  ' Re[z] for pole:',val1*Ha_eV,' eV'
      write(std_out,'(a,ES16.8,a)')  ' Im[z] for pole:',val2*Ha_eV,' eV'
      write(std_out,'(a,ES16.8)')    '      Amplitude:',osc*half/ABS(val1)
    end if
  end do

end subroutine init_peaks_from_grid
!!***

!!****f* m_model_screening/init_single_peak
!! NAME
!!  init_single_peak
!!
!! FUNCTION
!!  Initialise a single peak by using the behaviour along the imaginary axis
!!   and the main peak
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  yvals      = The function to be fitted in the complex plane
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!  prtvol     = Verbosity of diagnostics
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_model_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_single_peak(omega,refval,imfval,nomega,nfreqre,coeff,prtvol)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,nfreqre,prtvol
!arrays
  complex(dpc),intent(in)  :: omega(nomega)
  real(gwp)  ,intent(out) :: coeff(3)
  real(gwp)  ,intent(in) :: refval(nomega),imfval(nomega)

!Local variables-------------------------------
!scalars
  integer :: maxpos,idx
  real(gwp) :: pol,osc,gam,val1,val2,pol_sq

! *********************************************************************

  maxpos = MAXLOC(ABS(imfval(1:nfreqre)),1)
  if (maxpos==1) maxpos = MAXLOC(ABS(imfval(2:nfreqre)),1)
  pol = REAL(omega(maxpos))
  pol_sq = pol*pol
  osc = -refval(1)*pol_sq
  idx = nfreqre+1
  val2 = refval(idx)
  val1 = osc+val2*pol_sq+val2*AIMAG(omega(idx))*AIMAG(omega(idx))
  gam  = -ABS(val1/(AIMAG(omega(idx))*val2))

  coeff(1) = osc
  coeff(2) = pol
  coeff(3) = gam

  if (prtvol>9) then
    write(std_out,'(a,ES16.8)')    '  Osc. strength:',osc
    write(std_out,'(a,ES16.8,a)')  '  Peak location:',pol*Ha_eV,' eV'
    write(std_out,'(a,ES16.8,a)')  '     Peak width:',gam*Ha_eV,' eV'
    val2 = gam*half
    val1 = SIGN(1.0_gwp,pol*pol-val2*val2)*SQRT(ABS(pol*pol-val2*val2))
    write(std_out,'(a,ES16.8,a)')  ' Re[z] for pole:',val1*Ha_eV,' eV'
    write(std_out,'(a,ES16.8,a)')  ' Im[z] for pole:',val2*Ha_eV,' eV'
    write(std_out,'(a,ES16.8)')    '      Amplitude:',osc*half/ABS(val1)
  end if

end subroutine init_single_peak
!!***

!!****f* m_model_screening/init_peaks_even_dist
!! NAME
!!  init_peaks_even_dist
!!
!! FUNCTION
!!
!!  Distribute the peaks evenly along a line in the complex plane and
!!   normalise.
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  yvals      = The function to be fitted in the complex plane
!!  coeff      = The coefficients in order: A_1,B_1,C_1,A_2,B_2,C_2,...,A_n,B_n,C_n
!!  nomega     = number of fit points
!!  ncoeff     = number of coefficients
!!  prtvol     = Verbosity of diagnostics
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_peaks_even_dist(omega,fval,nomega,nfreqre,coeff,ncoeff,prtvol)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,nfreqre,ncoeff,prtvol
!arrays
  complex(dpc) ,intent(in)  :: omega(nomega)
  real(gwp)   ,intent(out) :: coeff(ncoeff)
  complex(gwpc),intent(in)  :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: npoles,ip,idx,iw
  real(gwp) :: delta,norm,div,val1,val2,osc,pol,gam

! *********************************************************************

  npoles = ncoeff/3
  div = real(npoles,gwp)

  delta = (omega(nfreqre)-omega(1))/(div+1.0_gwp)
  ! Integrate function along real axis (trapezoid rule) and have normalised
  ! oscillator strengths
  norm = fval(1)*half
  do iw=2,nfreqre-1
    norm = norm + fval(iw)
  end do
  norm = norm + fval(nfreqre)*half
  norm = norm*(omega(nfreqre)-omega(1))/real(nfreqre,gwp)
  norm = norm/div

  do ip=1,npoles
    idx = 3*(ip-1)+1
    pol = delta*ip   ! Position of maximum
    gam = 0.1_gwp ! Spread of function
    val2 = gam*half
    val1 = SQRT(pol*pol-val2*val2)
    osc = norm*val1*two
    coeff(idx  ) = osc!*(-1.0_gwp)**(ip-1)
    coeff(idx+1) = pol   ! Position of maximum
    coeff(idx+2) = -gam ! Spread of function
    if (prtvol>9) then
      write(std_out,'(a,a,i0)') ch10,' Pole no,: ',ip
      write(std_out,'(a,ES16.8)')    '  Osc. strength:',osc
      write(std_out,'(a,ES16.8,a)')  '  Peak location:',pol*Ha_eV,' eV'
      write(std_out,'(a,ES16.8,a)')  '     Peak width:',gam*Ha_eV,' eV'
      val2 = gam*half
      val1 = SIGN(1.0_gwp,pol*pol-val2*val2)*SQRT(ABS(pol*pol-val2*val2))
      write(std_out,'(a,ES16.8,a)')  ' Re[z] for pole:',val1*Ha_eV,' eV'
      write(std_out,'(a,ES16.8,a)')  ' Im[z] for pole:',val2*Ha_eV,' eV'
      write(std_out,'(a,ES16.8)')    '      Amplitude:',osc*half/ABS(val1)
    end if
  end do

end subroutine init_peaks_even_dist
!!***

!!****f* m_model_screening/print_peaks
!! NAME
!!  print_peaks
!!
!! FUNCTION
!!
!!  Find and output the location of peaks on the grid in a file
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  fvals      = The function to be fitted in the complex plane
!!  nomega     = number of fit points
!!  nfreqre    = number or imaginary gridlines
!!  nfreqim    = number of real gridlines
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_model_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_peaks(omega,fval,nomega,nfreqre,nfreqim)

 implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)   :: nomega,nfreqre,nfreqim
!arrays
  complex(dpc) ,intent(in)  :: omega(nomega)
  complex(gwpc),intent(in)  :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer :: ire,iim,unt_tmp
  integer :: idx1,idx2,idx3
  real(gwp) :: rez,imz,val1,val2,val3
  character(len=500) :: msg

! *********************************************************************

  if (open_file("grid_peak_tree.dat", msg, newunit=unt_tmp) /= 0) then
    MSG_ERROR(msg)
  end if

  do iim=nfreqim,1,-1
!    write(std_out,*) ' iim:',iim
!   Check first points
    idx1 = nfreqre+iim
    idx2 = nfreqre+nfreqim+1+(nfreqre-1)*(iim-1)
!    write(std_out,*) ' idx1:',idx1,' idx2:',idx2
    val1 = AIMAG(fval(idx1))
    val2 = AIMAG(fval(idx2))
    if (ABS(val1)>ABS(val2)) then
      rez = REAL(omega(idx1))
      imz = AIMAG(omega(idx1))
      write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val1
    end if
!   Do all but the last
    do ire=1,nfreqre-4
       idx1 = nfreqre+nfreqim+1+(nfreqre-1)*(iim-1)+ire
       idx2 = idx1+1
       idx3 = idx1+2
!       write(std_out,*) ' idx1:',idx1,' idx2:',idx2,' idx3:',idx3
       rez = REAL(omega(idx2))
       imz = AIMAG(omega(idx2))
       val1 = AIMAG(fval(idx1))
       val2 = AIMAG(fval(idx2))
       val3 = AIMAG(fval(idx3))
       if (((val1<val2).AND.(val2>val3))) then
         if (sign(1.0_gwp,val2)<zero) CYCLE
         write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val2
       else if (((val1>val2).AND.(val2<val3))) then
         if (sign(1.0_gwp,val2)>zero) CYCLE
         write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val2
       end if
    end do
!   Check last point
    idx1 = nfreqre+nfreqim+1+(nfreqre-1)*(iim-1)+nfreqre-3
    idx2 = idx1 + 1
!   write(std_out,*) ' idx1:',idx1,' idx2:',idx2
    rez = REAL(omega(idx2))
    imz = AIMAG(omega(idx2))
    val1 = AIMAG(fval(idx1))
    val2 = AIMAG(fval(idx2))
    if (ABS(val1)<ABS(val2)) then
      write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val2
    end if
  end do
! finally, do the purely real axis
! Check first points
  idx1 = 1; idx2 = 2
  val1 = AIMAG(fval(idx1)); val2 = AIMAG(fval(idx2))
  if (ABS(val1)>ABS(val2)) then
    rez = REAL(omega(idx1))
    imz = AIMAG(omega(idx1))
    write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val1
  end if
  do ire=2,nfreqre-3
    idx1 = ire; idx2 = idx1+1; idx3 = idx1+2
    rez = REAL(omega(idx2))
    imz = AIMAG(omega(idx2))
    val1 = AIMAG(fval(idx1))
    val2 = AIMAG(fval(idx2))
    val3 = AIMAG(fval(idx3))
    if (((val1<val2).AND.(val2>val3))) then
      if (sign(1.0_gwp,val2)<zero) CYCLE
      write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val2
    else if (((val1>val2).AND.(val2<val3))) then
      if (sign(1.0_gwp,val2)>zero) CYCLE
      write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val2
    end if
  end do
! Check last point
  idx1 = nfreqre-2
  idx2 = idx1 + 1
!  write(std_out,*) ' idx1:',idx1,' idx2:',idx2
  rez = REAL(omega(idx2))
  imz = AIMAG(omega(idx2))
  val1 = AIMAG(fval(idx1))
  val2 = AIMAG(fval(idx2))
  if (ABS(val1)<ABS(val2)) then
    write(unt_tmp,'(2f8.2,4x,ES16.8)') rez*Ha_eV,imz*Ha_eV,val2
  end if

  close(unt_tmp)

end subroutine print_peaks
!!***

!!****f* m_model_screening/find_peaks
!! NAME
!!  find_peaks
!!
!! FUNCTION
!!
!!  Find the location of the highest peaks along gridlines starting at the real axis
!!  and then moving towards higher imaginary frequencies. Stop when enough
!!  peaks to satisfy the number of poles needed has been found
!!
!! INPUTS
!!  omega      = (complex) Real and imaginary part of the frequency points
!!  fvals      = The function to be fitted in the complex plane
!!  nomega     = number of fit points
!!  nfreqre    = number or imaginary gridlines
!!  nfreqim    = number of real gridlines
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_model_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine find_peaks(fval,nomega,nfreqre,nfreqim,ploc,npoles,iline)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in)     :: nomega,nfreqre,nfreqim,npoles
  integer, intent(inout) :: iline
!arrays
  integer    ,intent(inout) :: ploc(npoles)
  complex(gwpc), intent(in) :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer    :: ire,iim,ipoles
  integer    :: idx1,idx2,idx3,ipol
  real(gwp) :: val1,val2,val3

!arrays
  integer :: ploc_prev(npoles)
  real    :: pval(npoles),pval_prev(npoles)

! *********************************************************************

  ploc=-1; ploc_prev=-1
  pval=zero; pval_prev=zero; ipol=1; ipoles=0

! First do a line along the real axis
  idx1 = 1; idx2 = 2
  val1 = AIMAG(fval(idx1)); val2 = AIMAG(fval(idx2))
  if (ABS(val1)>ABS(val2)) then
    ipoles = ipoles + 1
    ploc(1)=idx1; pval(1)=val1
    write(std_out,*) ' pval:',pval
    write(std_out,*) ' ploc:',ploc
  end if
  do ire=2,nfreqre-3
    idx1 = ire; idx2 = idx1+1; idx3 = idx1+2
    val1 = AIMAG(fval(idx1))
    val2 = AIMAG(fval(idx2))
    val3 = AIMAG(fval(idx3))
    if (((val1<val2).AND.(val2>val3))) then
      if (sign(1.0_gwp,val2)<zero) CYCLE
      ipoles = ipoles + 1
      if (ANY( ABS(pval(:))<ABS(val2) )) then
        ipol = MINLOC(ABS(pval(:)),1)
        ploc(ipol)=idx2; pval(ipol)=val2
        write(std_out,*) ' pval:',pval
        write(std_out,*) ' ploc:',ploc
      end if
    else if (((val1>val2).AND.(val2<val3))) then
      if (sign(1.0_gwp,val2)>zero) CYCLE
      ipoles = ipoles + 1
      if (ANY( ABS(pval(:))<ABS(val2) )) then
        ipol = MINLOC(ABS(pval(:)),1)
        ploc(ipol)=idx2; pval(ipol)=val2
        write(std_out,*) ' pval:',pval
        write(std_out,*) ' ploc:',ploc
      end if
    end if
  end do
! Check last point
  idx1 = nfreqre-2
  idx2 = idx1 + 1
!  write(std_out,*) ' idx1:',idx1,' idx2:',idx2
  val1 = AIMAG(fval(idx1))
  val2 = AIMAG(fval(idx2))
  if (ABS(val1)<ABS(val2)) then
    ipoles = ipoles + 1
    if (ANY( ABS(pval(:))<ABS(val2) )) then
      ipol = MINLOC(ABS(pval(:)),1)
      ploc(ipol)=idx2; pval(ipol)=val2
       write(std_out,*) ' pval:',pval
       write(std_out,*) ' ploc:',ploc
    end if
  end if
  write(std_out,'(a,i0)') ' Number of poles real axis:',ipoles


  ploc_prev = ploc; pval_prev = pval

! Do the rest of the imaginary grid until total
!  number of peaks found equals npoles or less
  do iim=1,nfreqim-1
    ploc=-1; pval=zero; ipol=1; ipoles=0
    idx1 = nfreqre+iim
    idx2 = nfreqre+nfreqim+1+(nfreqre-1)*(iim-1)
!    write(std_out,*) ' idx1:',idx1,' idx2:',idx2
    val1 = AIMAG(fval(idx1))
    val2 = AIMAG(fval(idx2))
    if (ABS(val1)>ABS(val2)) then
      ipoles = ipoles + 1
      ploc(1)=idx1; pval(1)=val1
      write(std_out,*) ' pval:',pval
      write(std_out,*) ' ploc:',ploc
    end if
!   Do all but the last
    do ire=1,nfreqre-4
       idx1 = nfreqre+nfreqim+1+(nfreqre-1)*(iim-1)+ire
       idx2 = idx1+1
       idx3 = idx1+2
!       write(std_out,*) ' idx1:',idx1,' idx2:',idx2,' idx3:',idx3
       val1 = AIMAG(fval(idx1))
       val2 = AIMAG(fval(idx2))
       val3 = AIMAG(fval(idx3))
       if (((val1<val2).AND.(val2>val3))) then
         if (sign(1.0_gwp,val2)<zero) CYCLE
         ipoles = ipoles + 1
         if (ANY( ABS(pval(:))<ABS(val2) )) then
           ipol = MINLOC(ABS(pval(:)),1)
           ploc(ipol)=idx2; pval(ipol)=val2
           write(std_out,*) ' pval:',pval
           write(std_out,*) ' ploc:',ploc
         end if
       else if (((val1>val2).AND.(val2<val3))) then
         if (sign(1.0_gwp,val2)>zero) CYCLE
         ipoles = ipoles + 1
         if (ANY( ABS(pval(:))<ABS(val2) )) then
           ipol = MINLOC(ABS(pval(:)),1)
           ploc(ipol)=idx2; pval(ipol)=val2
           write(std_out,*) ' pval:',pval
           write(std_out,*) ' ploc:',ploc
         end if
       end if
    end do
!   Check last point
    idx1 = nfreqre+nfreqim+1+(nfreqre-1)*(iim-1)+nfreqre-3
    idx2 = idx1 + 1
!   write(std_out,*) ' idx1:',idx1,' idx2:',idx2
    val1 = AIMAG(fval(idx1))
    val2 = AIMAG(fval(idx2))
    if (ABS(val1)<ABS(val2)) then
      ipoles = ipoles + 1
      if (ANY( ABS(pval(:))<ABS(val2) )) then
        ipol = MINLOC(ABS(pval(:)),1)
        ploc(ipol)=idx2; pval(ipol)=val2
         write(std_out,*) ' pval:',pval
         write(std_out,*) ' ploc:',ploc
      end if
    end if
    write(std_out,'(2(a,i0))') ' Line,:',iim,' ipoles:',ipoles
    if (ipoles<=npoles) then
      iline = iim - 1
      ploc = ploc_prev
      EXIT
    end if

    ploc_prev = ploc; pval_prev = pval

 end do

end subroutine find_peaks
!!***

!!****f* m_model_screening/remove_phase
!! NAME
!!  remove_phase
!!
!! FUNCTION
!!  Find out what the complex phase factor is for off-diagonal elements
!!   and unmix the components.
!!
!! INPUTS
!!  fvals  = The function to be fitted in the complex plane.
!!  nomega = Number of fit points.
!!  nfreqre    = number or imaginary gridlines
!!  phase  = The phase angle.
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine remove_phase(fval,nomega,phase)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,    intent(in)  :: nomega
  real(gwp), intent(out) :: phase
!arrays
  complex(gwpc), intent(inout) :: fval(nomega)

!Local variables-------------------------------
!scalars
  integer       :: io
  real(gwp)    :: a,b,retemp,imtemp

! *********************************************************************

! The phase can be found by checking when the function is
!  identically zero along the imaginary axis
  if (ABS(AIMAG(fval(1)))<tol14) then ! Phase is zero
    phase = zero
    RETURN
  else if (ABS(REAL(fval(1)))<tol14) then ! Phase is exactly pi/2
    phase = pi*half
    a = zero
    b = -1.0_gwp
  else
    phase = ATAN(AIMAG(fval(1))/REAL(fval(1)))
    a = COS(phase)
    b = -SIN(phase)
  end if
! Rotate values
  do io=1,nomega
    retemp = REAL(fval(io))
    imtemp = AIMAG(fval(io))
    fval(io) = CMPLX(a*retemp-b*imtemp,a*imtemp+b*retemp)
  end do

end subroutine remove_phase
!!***

END MODULE m_model_screening
!!***
