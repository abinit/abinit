!!****f* ABINIT/m_polynomial_filter
!! NAME
!! m_polynomial_filter
!!
!! FUNCTION
!! This module contains routines to implement various polynomial filters (scalar) and 
!! damping coefficients. 
!!
!! COPYRIGHT
!! Copyright (C) 2018-2026 ABINIT group (IML)
!! This file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

module m_polynomial_filter

    use defs_basis
    use defs_abitypes
    use m_abicore
    use m_errors
    use m_time, only : timab
    use m_sort, only: sort_dp

    use m_xmpi
    use m_xomp
#ifdef HAVE_OPENMP
    use omp_lib
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize
#endif

#if defined(HAVE_GPU_MARKERS)
    use m_nvtx_data
#endif

    implicit none

    private

    ! Timers
    !---------------------------------------------------
    integer, parameter :: tim_swap        = 1761
    integer, parameter :: tim_RR_q        = 1759
    integer, parameter :: tim_barrier     = 1764
    integer, parameter :: tim_copy        = 1765
    integer, parameter :: tim_Bortho_X    = 1641
    integer, parameter :: tim_getAX_BX    = 1754
    integer, parameter :: tim_invovl      = 1755
    integer, parameter :: tim_lanczos     = 2167
    integer, parameter :: tim_trace       = 2168

    ! Public methods
    !-------------------------------------------------
    public :: jackson_step_coeffs
    public :: erf_step_coeffs
    public :: smooth_step_coeffs
    public :: lanczos_step_coeffs
    public :: buildChebyshevJacksonCoeffs ! left from previous implementation
    public :: cheb_poly
    public :: bandpass_sca
    public :: print_scalar_filter

    CONTAINS  
!=====================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/jackson_step_coeffs
!! NAME
!! jackson_step_coeffs
!! 
!! FUNCTION
!! Jackson damped coefficients
!! 
!! SOURCE

  subroutine jackson_step_coeffs(a,b,lmin,lmax,deg,ctilde)
      
      implicit none

      integer, intent(in) :: deg
      real(dp), intent(in) :: a,b,lmin,lmax
      real(dp), intent(out) :: ctilde(1:deg+1)
      integer :: k
      real(dp) :: c, r, a_, b_, ck, gk, alpha, cotv

      c  = (lmin + lmax)/2.d0
      r  = (lmax - lmin)/2.d0
      a_ = (a - c)/r
      b_ = (b - c)/r
      alpha = pi/(deg+2)
      cotv = cos(alpha)/sin(alpha)

      do k = 0, deg
         ! Chebyshev coefficient for step function
         if (k == 0) then
            ck = (acos(a_) - acos(b_))/pi
         else
            ck = 2.0d0/pi*(sin(k*acos(a_)) - sin(k*acos(b_)))/k
         end if

         ! Jackson damping factor
         gk = ((deg - k + 1)*cos(pi*k/(deg+1)) + sin(pi*k/(deg+1))*cotv) / (deg+1)

         ctilde(k+1) = ck * gk
      end do

  end subroutine jackson_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/erf_step_coeffs
!! NAME
!! erf_step_coeffs
!! 
!! FUNCTION
!! Erf damped coefficients
!! 
!! SOURCE

  function erf_step_coeffs(b, ndeg, sigma, Ngrid) result(coeffs)

      implicit none
      real(dp), intent(in) :: b, sigma
      integer, intent(in) :: ndeg, Ngrid

      real(dp) :: coeffs(ndeg+1)
      integer :: i, k
      real(dp) :: x(Ngrid+1)
      real(dp) :: f(Ngrid+1)

      x = (/ ( cos(Pi*(i-2)/(Ngrid-1) ) , i=1,Ngrid+1) /)
      f = (/ ( 0.5 * (1.0 - fast_erf((x(i) - b)/sigma)), i=1,Ngrid+1) /)
      
      do k=0, ndeg
        coeffs(k+1) = (2.0/Ngrid) * sum( (/ (f(i)*cos(k*acos(x(i))), i=1,Ngrid+1 ) /) )
      end do
      coeffs(1) = coeffs(1) / 2.d0 

  end function erf_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/fast_erf
!! NAME
!! fast_erf
!! 
!! FUNCTION
!! Fast approximation of erf(x) using Abramowitz & Stegun 7.1.26

  function fast_erf(x) result(erf_val)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: erf_val
    real(dp) :: t, tau, ax
    real(dp), parameter :: p  = 0.3275911_dp
    real(dp), parameter :: a1 = 0.254829592_dp
    real(dp), parameter :: a2 = -0.284496736_dp
    real(dp), parameter :: a3 = 1.421413741_dp
    real(dp), parameter :: a4 = -1.453152027_dp
    real(dp), parameter :: a5 = 1.061405429_dp

    ax = abs(x)
    t = 1.0_dp / (1.0_dp + p * ax)
    tau = (((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t) * exp(-ax*ax)
    erf_val = 1.0_dp - tau
    if (x < 0.0_dp) erf_val = -erf_val
  end function fast_erf
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/smooth_step_coeffs
!! NAME
!! smooth_step_coeffs
!! 
!! FUNCTION
!! Smooth damped coefficients
!! 
!! SOURCE

  function smooth_step_coeffs(b, ndeg, alpha, Ngrid) result(coeffs)

      implicit none
      real(dp), intent(in) :: b, alpha
      integer, intent(in) :: ndeg, Ngrid

      real(dp) :: coeffs(ndeg+1)
      integer :: i, k
      real(dp) :: x(Ngrid+1)
      real(dp) :: f(Ngrid+1)

      x = (/ ( cos(Pi*(i-2)/(Ngrid-1) ) , i=1,Ngrid+1) /)
      f = (/ ( 0.5 * (1.0 - tanh( alpha*(x(i) - b) )) , i=1,Ngrid+1) /)
      
      do k=0, ndeg
        coeffs(k+1) = (2.0/Ngrid) * sum( (/ (f(i)*cos(k*acos(x(i))), i=1,Ngrid+1 ) /) )
      end do
      coeffs(1) = coeffs(1) / 2.d0 

  end function smooth_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/lanczos_step_coeffs
!! NAME
!! lanczos_step_coeffs
!! 
!! FUNCTION
!! Chebyshev coefficients with Lanczos damping
!! 
!! SOURCE

function lanczos_step_coeffs(b, ndeg) result(coeffs)

    implicit none
    real(dp), intent(in) :: b
    integer, intent(in) :: ndeg

    real(dp) :: coeffs(ndeg+1)
    integer :: j
    real(dp) :: x(ndeg+1)
    real(dp) :: f(ndeg+1)
    real(dp) :: k_array(ndeg+1)

    x = (/ ( cos(Pi*(2*j-1)/(2*(ndeg+1)) ) , j=1,ndeg+1) /)
    f = merge(1.0, 0.0, x < b) ! f=1 if x<b else 0
    k_array = (/ (j, j=1,ndeg+1) /)
    do j = 1, ndeg+1
        coeffs(j) = 2.0d0 / (ndeg+1) * &
            sum( f(:) * cos( pi*(j-1)*(2.0d0*k_array(:)-1.0d0) / (2.0d0*(ndeg+1)) ) )
    end do
    coeffs(1) = coeffs(1) / 2.d0 
    do j = 2, ndeg+1
        coeffs(j) = coeffs(j) * sin(pi*(j-1)/(ndeg+1)) / (pi*(j-1)/(ndeg+1))
    end do 

end function lanczos_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/buildChebyshevJacksonCoeffs
!! NAME
!! buildChebyshevJacksonCoeffs
!!
!! SOURCE
subroutine buildChebyshevJacksonCoeffs(ls, us, ndeg_filter, cja)

    implicit none

    integer, intent(in) :: ndeg_filter
    real(dp), intent(in) :: ls, us ! scaled to [-1,1)
    real(dp), intent(inout) :: cja(ndeg_filter+1)

    integer :: ideg
    real(dp) :: cdeg
    real(dp) :: mu, damp
    
    ! *********************************************************************

    cdeg = Pi/(ndeg_filter+2)
    mu = 1.d0/Pi*(ACOS(ls)-ACOS(us))
    damp = 1.d0 ! Jackson damping

    cja(1) = mu*damp

    do ideg = 0, ndeg_filter - 1
        
        ! Accumulate X with weight in Xsum for bandpass filters
        mu = 2/Pi * (SIN((ideg+1)*ACOS(ls)) - SIN((ideg+1)*ACOS(us)))/(ideg+1)
        damp = ((1 - (ideg+1)/(ndeg_filter+2))*SIN(cdeg)*COS((ideg+1)*cdeg) + &
                1/(ndeg_filter+2)*COS(cdeg)*SIN((ideg+1)*cdeg))/SIN(cdeg)
        
        cja(ideg+2) = mu*damp
        
    end do

end subroutine buildChebyshevJacksonCoeffs
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/cheb_poly
!! NAME
!! cheb_poly
!!
!! FUNCTION
!! Compute Chebyshev polynomial
!!
!! INPUTS
!!  xx= input variable
!!  aa= left bound of the interval
!!  bb= right bound of the interval
!!  nn=
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

function cheb_poly(xx,nn,aa,bb) result(yy)

    implicit none

    ! Arguments ------------------------------------
    integer,  intent(in) :: nn
    real(dp), intent(in) :: xx, aa, bb
    real(dp)             :: yy

    ! Local variables-------------------------------
    integer  :: ii
    real(dp) :: xred,yim1,temp

    ! *************************************************************************

    xred = (xx-(aa+bb)/2)/(bb-aa)*2
    yy = xred
    yim1 = 1
    do ii= 2, nn
        temp = yy
        yy = 2*xred*yy - yim1
        yim1 = temp
    end do

end function cheb_poly
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/bandpass_sca
!! NAME
!! bandpass_sca
!!
!! FUNCTION
!! Computes delta-Dirac polynomial filter f(x) approximated by a Chebyshev
!! expansion of order deg centered at gamma, evaluated at point x=t
!!
!! INPUTS
!!  t=      scalar to evaluate filter on
!!  deg=    order of Chebyshev expansion
!!  gam=    center of Chebyshev expansion
!!
!! OUTPUT
!!  res
!!
!! SOURCE

function bandpass_sca(t, deg, gam) result(f_t)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in ) :: t, gam
    integer , intent(in ) :: deg

    real(dp) :: f_t
    
    !Local variables-------------------------------
    real(dp) :: yt0, yt, yg0, yg, yt_swap, yg_swap
    real(dp) :: mu, damp, rho, rhog, theta
    integer  :: i
    
    ! *********************************************************************

    ! init cheby of deg=0,1 eval at t,gamma
    yt0 = 1.d0
    yt = t
    
    yg0 = 1.d0
    yg = gam

    ! init delta-Dirac filters of deg=0
    theta = Pi/(deg + 1)
    damp = SIN(theta) / theta
    rho = 0.5d0 + gam * damp * yt
    rhog = 0.5d0 + gam * damp * yg

    do i=2,deg 

        ! Update Chebyshev polynomials
        yt_swap = yt
        yt = 2 * t * yt - yt0
        yt0 = yt_swap
        
        yg_swap = yg
        yg = 2 * gam * yg - yg0
        yg0 = yg_swap

        ! Update delta-Dirac filters
        mu = COS(i * ACOS(gam))
        damp = SIN(i * theta) / (i * theta)
        rho = rho + mu * damp * yt
        rhog = rhog + mu * damp * yg
        
    end do

    f_t = rho / rhog

end function bandpass_sca
!!***

!----------------------------------------------------------------------

!!****f* m_polynomial_filter/print_scalar_filter
!! NAME
!! print_scalar_filter
!! 
!! FUNCTION
!! Print x,f(x) for every x in (a,b).
!! 
!! INPUTS
!! lb,ub= interval to amplify/vanish
!! glb,gub= guaranteed spectral bounds used to scale to [-1,1]
!! ndeg= polynomial degree of filter
!! is_lowpass= flag. If true then Chebyshev if false Chebyshev-Jackson

subroutine print_scalar_filter(a, b, lb, ub, glb, gub, ndeg, is_lowpass)

    implicit none

    real(dp), intent(in) :: a,b,lb,ub,glb,gub
    integer, intent(in) :: ndeg
    logical, intent(in) :: is_lowpass

    integer :: npt, ipt
    real(dp) :: c, r, pt, fun_pt

    npt = 100
    c = (gub + glb)/2.d0
    r = (gub - glb)/2.d0
    write(std_out,*) ' '
    write(std_out,*) 'Plot filter ==== x | f(x)'
    if (is_lowpass) then
        do ipt=1,npt
            pt = a + (ipt-1)*(b-a)/npt ! unscaled!
            fun_pt = cheb_poly(pt,ndeg,ub,gub)
            write(std_out,*) pt, fun_pt
        end do
    else
        do ipt=1,npt
            pt = (a + (ipt-1)*(b-a)/npt - c)/r ! scaled!
            fun_pt = bandpassIndicator_sca(pt,(a-c)/r,(b-c)/r,ndeg)
            write(std_out,*) pt, fun_pt
        end do
    end if
    write(std_out,*) ' '

end subroutine print_scalar_filter
!!***

end module m_polynomial_filter
!!***
