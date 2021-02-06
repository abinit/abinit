!!****m* ABINIT/m_eliashberg_1d
!! NAME
!!  m_eliashberg_1d
!!
!! FUNCTION
!!  Solve the Eliashberg equations in the isotropic case
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MVer)
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

module m_eliashberg_1d

 use defs_basis
 use defs_elphon
 use m_errors
 use m_abicore
 use m_io_tools
 use m_abicore

 use m_numeric_tools,   only : simpson_int

 implicit none

 private
!!***

 public :: eliashberg_1d
!!***

contains
!!***

!!****f* m_eliashberg_1d/eliashberg_1d
!!
!! NAME
!! eliashberg_1d
!!
!! FUNCTION
!!  Solve the Eliashberg equations in the isotropic case
!!  First the linearized case, which allows the estimation of Tc
!!  then the full case which gives the gap as a function of temperature.
!!
!! INPUTS
!!  a2f_1d = alpha^2F function averaged over the FS (only energy dependence)
!!  elph_ds = datastructure with phonon matrix elements
!!  gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!  natom = number of atoms
!!
!! OUTPUT
!!
!! PARENTS
!!      m_elphon
!!
!! CHILDREN
!!
!! NOTES
!!  na2f = number of frequency points for alpha^2F function
!!
!! SOURCE

subroutine eliashberg_1d(a2f_1d,elph_ds,mustar)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: mustar
 type(elph_type),intent(in) :: elph_ds
!arrays
 real(dp),intent(in) :: a2f_1d(elph_ds%na2f)

!Local variables-------------------------------
  ! for diagonalization of gammma matrix
  ! output variables for gtdyn9+dfpt_phfrq
!scalars
 integer :: iiter,imatsu
 integer :: maxiter,nmatsu,unit_del,unit_lam,unit_z
 real(dp) :: maxeigval,omega_cutoff
 real(dp) :: tc
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: delta_1d(:),lambda_1d(:),mm_1d(:,:),z_1d(:)

! *********************************************************************

 call wrtout(std_out,'Solving the 1-D Eliashberg equation (isotropic case)',"COLL")

 if (elph_ds%nsppol /= 1) then
   write(std_out,*) 'eliashberg_1d is not coded for nsppol > 1 yet'
   return
 end if

!maximum number of iterations to find T_c
 maxiter=30

!Fix nmatsu. Should add test in iiter loop to check if
!omega_cutoff is respected
 nmatsu = 50
!write(std_out,*) ' eliashberg_1d : nmatsu = ', nmatsu

 ABI_MALLOC(lambda_1d,(-nmatsu:nmatsu))
 ABI_MALLOC(z_1d,(-nmatsu:nmatsu))
 ABI_MALLOC(delta_1d,(-nmatsu:nmatsu))
 ABI_MALLOC(mm_1d,(-nmatsu:nmatsu,-nmatsu:nmatsu))

 unit_lam=get_unit()
 fname=trim(elph_ds%elph_base_name) // "_LAM"
 open (UNIT=unit_lam,FILE=fname,STATUS='REPLACE')
 unit_z=get_unit()
 fname=trim(elph_ds%elph_base_name) // "_Z"
 open (UNIT=unit_z,FILE=fname,STATUS='REPLACE')
 unit_del=get_unit()
 fname=trim(elph_ds%elph_base_name) // "_DEL"
 open (UNIT=unit_del,FILE=fname,STATUS='REPLACE')

!
!1) use linearized Eliashberg equation to find Tc
!! \sum_j \mathbf{M}_{ij} \Delta_j = \zeta \cdot \Delta_i $ $i,j = 1 .. n_{\mathrm{Matsubara}}$
!$\zeta = 1$ gives T$_c$ $\beta = \frac{1}{\mathrm{T}}$ $\omega_i = (2 i + 1) \pi \mathrm{T}$
!! \mathbf{M}_{ij} = \frac{\pi}{\beta} \frac{\lambda (\omega_i - \omega_j)}{Z (\omega_i)}$
!! Z (\omega_i) = 1 + \frac{\pi}{\beta \omega_i} \sum_j \lambda(\omega_i - \omega_j) \mathrm{sgn}(\omega_j)$
!

!initial guess for T$_c$ in Hartree (1Ha =3.067e5 K)
 tc = 0.0001
!
!big iterative loop
!
 do iiter=1,maxiter

   omega_cutoff = (two*nmatsu+one) * pi * tc

!
!  calculate array of lambda values
!
   call eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)
   write (unit_lam,'(a)') '#'
   write (unit_lam,'(a)') '# ABINIT package : lambda file'
   write (unit_lam,'(a)') '#'
   write (unit_lam,'(a,I10,a)') '# lambda_1d array containing 2*', nmatsu, '+1 Matsubara frequency points'
   write (unit_lam,'(a,E16.6,a,E16.6)') '#  from ', -omega_cutoff, ' to ', omega_cutoff
   write (unit_lam,'(a)') '#  lambda_1d is the frequency dependent coupling constant '
   write (unit_lam,'(a)') '#  in the Eliashberg equations '
   write (unit_lam,'(a)') '#'
   do imatsu=-nmatsu,nmatsu
     write (unit_lam,*) imatsu,lambda_1d(imatsu)
   end do
   write (unit_lam,*)

!
!  calculate array of z values
!
   call eli_z_1d (lambda_1d,nmatsu,z_1d)
   write (unit_z,'(a)') '#'
   write (unit_z,'(a)') '# ABINIT package : Z file'
   write (unit_z,'(a)') '#'
   write (unit_z,'(a,I10,a)') '# z_1d array containing 2*', nmatsu, '+1 Matsubara frequency points'
   write (unit_z,'(a,E16.6,a,E16.6)') '# from ', -omega_cutoff, ' to ', omega_cutoff
   write (unit_z,'(a)') '# z_1d is the renormalization factor in the Eliashberg equations'
   write (unit_z,'(a)') '#'
   do imatsu=-nmatsu,nmatsu
     write (unit_z,*) imatsu,z_1d(imatsu)
   end do

!  !
!  ! apply M matrix until a maximal eigenvalue is found.
!  !
!  call eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,z_1d)

!
!  diagonalize M brute forcefully
!
   call eli_diag_m_1d(delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)

   write (unit_del,'(a)') '#'
   write (unit_del,'(a)') '# eliashberg_1d : delta_1d = '
   write (unit_del,'(a)') '#'
   write (unit_del,'(a,i6,a)') '# delta_1d array containing 2*', nmatsu, '+1 Matsubara frequency points'
   write (unit_z,'(a,E16.6,a,E16.6)') '# from ', -omega_cutoff, ' to ', omega_cutoff
   write (unit_z,'(a)') '# delta_1d is the gap function in the Eliashberg equations'
   write (unit_z,'(a)') '#'
   do imatsu=-nmatsu,nmatsu
     write (unit_del,*) imatsu,delta_1d(imatsu)
   end do
   write (unit_del,*)

!  if eigenvalue is < 1 increase T
!  else if eigenvalue is > 1 decrease T
!  if eigenvalue ~= 1 stop
!
   if (abs(maxeigval-one) < tol8) then
     write(std_out,*) 'Eliashberg Tc found = ', tc, ' (Ha) = ', tc/kb_HaK, ' (K)'
     exit
   else if (maxeigval > 0.001_dp) then
     tc = tc * maxeigval
   else
     write(std_out,*) 'maxeigval is very small'
     tc = tc * 1000.0_dp
   end if


 end do
!end iiter do

 if (abs(maxeigval-one) > tol8) then
   write(std_out,*) 'eliashberg_1d : Tc not converged. ', maxeigval, ' /= 1'
   write(std_out,*) 'Eliashberg Tc nonetheless = ', tc, ' (Ha) = ', tc/kb_HaK, ' (K)'
 end if

 ABI_FREE(lambda_1d)
 ABI_FREE(z_1d)
 ABI_FREE(delta_1d)
 ABI_FREE(mm_1d)

 close (UNIT=unit_z)
 close (UNIT=unit_lam)
 close (UNIT=unit_del)

 write(std_out,*) ' eliashberg_1d : end '

end subroutine eliashberg_1d
!!***

!!****f* m_eliashberg_1d/eli_app_m_1d
!!
!! NAME
!! eli_app_m_1d
!!
!! FUNCTION
!!   Apply the linearized Eliashberg matrix once to the input vector.
!!
!! INPUTS
!!   lambda_1d = coupling constant as a function of frequency
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!   z_1d = renormalization Z as a function of frequency
!!
!! SIDE EFFECTS
!!   delta_1d = imaginary gap function as a function of frequency changed
!!
!! PARENTS
!!      m_eliashberg_1d
!!
!! CHILDREN
!!
!! SOURCE


subroutine eli_app_m_1d (delta_1d,lambda_1d,nmatsu,z_1d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu),z_1d(-nmatsu:nmatsu)
 real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,jmatsu,miguelflag
 real(dp) :: zfact
!arrays
 real(dp) :: delta_tmp(-nmatsu:nmatsu),freqfact(-nmatsu:nmatsu)

! *********************************************************************

 miguelflag = 0


 do imatsu=-nmatsu,nmatsu
   freqfact(imatsu) = one / abs(two*imatsu+one)
 end do

 delta_tmp(:) = delta_1d(:)

 if (miguelflag == 1) then
   do imatsu=-nmatsu,nmatsu
!    zfact = pi*tc / z_1d(imatsu)
     zfact = one / z_1d(imatsu)

     do jmatsu=max(-nmatsu,-nmatsu+imatsu),min(nmatsu,nmatsu+imatsu)
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + delta_1d(jmatsu) &
&       * lambda_1d(imatsu-jmatsu) &
&       * freqfact(jmatsu)
     end do
     delta_tmp(imatsu) = delta_tmp(imatsu)*zfact
   end do

 else

!  i < 0
   do imatsu=-nmatsu,-1

!    j < 0
     do jmatsu=max(-nmatsu,-nmatsu+imatsu),-1
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       - lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do
!    j > 0
     do jmatsu=0,min(nmatsu,nmatsu+imatsu)
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do

   end do

!  i > 0
   do imatsu=0,nmatsu

!    j < 0
     do jmatsu=max(-nmatsu,-nmatsu+imatsu),-1
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do
!    j > 0
     do jmatsu=0,min(nmatsu,nmatsu+imatsu)
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       - lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do

   end do

 end if

 delta_1d(:) = delta_tmp(:)

end subroutine eli_app_m_1d
!!***

!!****f* m_eliashberg_1d/eli_diag_m_1d
!!
!! NAME
!! eli_diag_m_1d
!!
!! FUNCTION
!!  diagonalize M matrix. Heavy and should be avoided for production.
!!  Actually, since M is not symmetrical, diagonalize M^{t} M and
!!  get right-eigenvalues and vectors
!!
!! INPUTS
!!   lambda_1d = coupling constant as a function of frequency
!!   mustar = Coulomb potential parameter in Eliashberg equation
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!   z_1d = renormalization Z as a function of frequency
!!
!! OUTPUT
!!   maxeigval = estimation for maximum eigenvalue of M
!!
!! SIDE EFFECTS
!!   delta_1d = imaginary gap function as a function of frequency
!!
!! PARENTS
!!      m_eliashberg_1d
!!
!! CHILDREN
!!
!! SOURCE

subroutine eli_diag_m_1d (delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)

 use m_linalg_interfaces

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: mustar,tc
 real(dp),intent(out) :: maxeigval
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu),z_1d(-nmatsu:nmatsu)
 real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,info,jmatsu,kmatsu,lwork,tmiguel
 real(dp) :: si,sj,sqtimat,sqtjmat
!arrays
 real(dp) :: mtm_eig(2*nmatsu+1),symm_mtm(-nmatsu:nmatsu,-nmatsu:nmatsu)
 real(dp) :: work(3*(2*nmatsu+1))

! *********************************************************************

 tmiguel = 0

 if (tmiguel == 1) then
   do imatsu=-nmatsu,nmatsu
     do jmatsu=-nmatsu,nmatsu

       symm_mtm(imatsu,jmatsu) = zero
       do kmatsu=max(-nmatsu+imatsu,-nmatsu+jmatsu,-nmatsu),min(nmatsu+imatsu,nmatsu+jmatsu,nmatsu)
         symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) &
&         +  lambda_1d(kmatsu-imatsu)*lambda_1d(kmatsu-jmatsu) &
&         /  ( z_1d(kmatsu)*z_1d(kmatsu) )
       end do
!      symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) / ((two*imatsu+one)*(two*jmatsu+one))
       symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) * pi * tc * pi * tc
     end do
   end do

 else

   symm_mtm(:,:) = -mustar

   si = -one
   do imatsu=-nmatsu,nmatsu
     sqtimat = one / sqrt(two*abs(imatsu)+one)
     if (imatsu == 0) si = one
     sj = -one
     do jmatsu=max(-nmatsu,-nmatsu+imatsu),min(nmatsu,nmatsu+imatsu)
       sqtjmat = one / sqrt(two*abs(jmatsu)+one)
       if (jmatsu == 0) sj = one

       symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) &
&       + lambda_1d(imatsu-jmatsu)*sqtimat*sqtjmat

       symm_mtm(imatsu,imatsu) = symm_mtm(imatsu,imatsu) &
&       - lambda_1d(imatsu-jmatsu)*si*sj*sqtimat*sqtimat
     end do
   end do

 end if

 lwork = 3*(2*nmatsu+1)
 call DSYEV('V', 'U', 2*nmatsu+1, symm_mtm, 2*nmatsu+1, mtm_eig, work, lwork, info )

 write(std_out,*) 'last eigenvalues = '
 write(std_out,*) mtm_eig(2*nmatsu-9:2*nmatsu+1)

 do imatsu=-nmatsu,nmatsu
   delta_1d(imatsu) = symm_mtm(imatsu,nmatsu)*sqrt(two*abs(imatsu)+one)
 end do

 maxeigval = mtm_eig(2*nmatsu+1)

end subroutine eli_diag_m_1d
!!***

!!****f* m_eliashberg_1d/eli_lambda_1d
!!
!! NAME
!! eli_lambda_1d
!!
!! FUNCTION
!!  In the solving of the 1D (energy only) Eliashberg equations, calculate
!!  the lambda, which is the e-p coupling strength. See Allen and Mitrovic
!!  Solid State Physics vol 37 ed Ehrenreich Seitz and Turnbull, p.45 [[cite:Allen1983c]]
!!
!! INPUTS
!!   a2f_1d = 1D alpha2F function
!!   elph_ds = elphon dataset
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!
!! OUTPUT
!!   lambda_1d = coupling constant as a function of frequency
!!
!! PARENTS
!!      m_eliashberg_1d
!!
!! CHILDREN
!!
!! NOTES
!!  lambda is used at points which are differences of Matsubara freqs,
!!  and hence is tabulated on points going through 0.
!!
!! SOURCE

subroutine eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: tc
 type(elph_type),intent(in) :: elph_ds
!arrays
 real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
 real(dp),intent(out) :: lambda_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,iomega
 real(dp) :: nu_matsu,nu_matsu2,omega,domega
!arrays
 real(dp) :: lambda_int(elph_ds%na2f),tmplambda(elph_ds%na2f)

! *********************************************************************
!
!MG: the step should be calculated locally using nomega and the extrema of the spectrum.
!One should not rely on previous calls for the setup of elph_ds%domega
!I will remove elph_ds%domega since mka2f.F90 will become a method of gamma_t
 domega =elph_ds%domega

 do imatsu=-nmatsu,nmatsu
   nu_matsu = (two*imatsu)*pi*tc
   nu_matsu2 = nu_matsu*nu_matsu

   tmplambda(:) = zero
   omega=domega
   do iomega=2,elph_ds%na2f
     tmplambda(iomega) = a2f_1d(iomega) * two * omega / (nu_matsu2 + omega*omega)
     omega=omega+domega
   end do
   call simpson_int(elph_ds%na2f,domega,tmplambda,lambda_int)

   lambda_1d(imatsu) = lambda_int(elph_ds%na2f)
 end do

end subroutine eli_lambda_1d
!!***

!!****f* m_eliashberg_1d/eli_m_iter_1d
!! NAME
!! eli_m_iter_1d
!!
!! FUNCTION
!!  Find largest eigenvalue of M matrix, to deduce superconducting Tc
!!
!! INPUTS
!!   lambda_1d = coupling constant as a function of frequency
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!   z_1d = renormalization Z as a function of frequency
!!
!! OUTPUT
!!   maxeigval = estimation for maximum eigenvalue of M
!!
!! SIDE EFFECTS
!!   delta_1d = imaginary gap function as a function of frequency
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,z_1d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(out) :: maxeigval
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu),z_1d(-nmatsu:nmatsu)
 real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: iiterm,imatsu,nfilter,ngeteig
 real(dp) :: dnewnorm,dnorm,fact
!arrays
 real(dp) :: delta_old(-nmatsu:nmatsu)

! *********************************************************************

 nfilter = 10
 ngeteig = 10


!
!1) apply M matrix enough times to filter out largest eigenvalue
!
 do iiterm=1,nfilter
   call eli_app_m_1d (delta_1d,lambda_1d,nmatsu,z_1d)

!  DEBUG
!  dnorm=zero
!  do imatsu=-nmatsu,nmatsu
!  dnorm = dnorm + delta_1d(imatsu)*delta_1d(imatsu)/(two*imatsu+one)
!  end do
!  dnorm = sqrt(dnorm)
!  write(std_out,*) 'eli_m_iter_1d : dnorm ', dnorm
!  ENDDEBUG
 end do

!
!2) calculate norm
!
 dnorm=zero
 do imatsu=-nmatsu,nmatsu
   dnorm = dnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
 end do
 dnorm = sqrt(dnorm)

!normalize delta_1d
 delta_1d(:) = delta_1d(:) / dnorm

 delta_old = delta_1d

!DEBUG
!dnewnorm=zero
!do imatsu=-nmatsu,nmatsu
!dnewnorm = dnewnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
!end do
!dnewnorm = sqrt(dnewnorm)
!write(std_out,*) 'eli_m_iter_1d : dnewnorm1 ', dnewnorm
!ENDDEBUG

!
!3) re-apply M matrix ngeteig times
!
 do iiterm=1,ngeteig
   call eli_app_m_1d (delta_1d,lambda_1d,nmatsu,z_1d)
!  DEBUG
!  dnewnorm=zero
!  do imatsu=-nmatsu,nmatsu
!  dnewnorm = dnewnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
!  end do
!  dnewnorm = sqrt(dnewnorm)
!  write(std_out,*) 'eli_m_iter_1d : dnewnorm ', dnewnorm
!
!  do imatsu=-nmatsu,nmatsu
!  write (112,*) imatsu,delta_1d(imatsu)/delta_old(imatsu)
!  end do
!  write (112,*)
!  delta_old = delta_1d
!  ENDDEBUG
 end do

!
!4) calculate new norm and estimate eigenvalue
!
 dnewnorm=zero
 do imatsu=-nmatsu,nmatsu
   dnewnorm = dnewnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
 end do
 dnewnorm = sqrt(dnewnorm)

!maxeigval = exp ( log(dnewnorm/dnorm) / ngeteig )
 maxeigval = exp ( log(dnewnorm) / ngeteig )

 write(std_out,*) 'eli_m_iter_1d : maxeigval =' , maxeigval
!fact = exp(-log(maxeigval) * (ngeteig+nfilter))
 fact = exp(-log(maxeigval) * (ngeteig))
 do imatsu=-nmatsu,nmatsu
   delta_1d(imatsu) = delta_1d(imatsu) * fact
 end do

end subroutine eli_m_iter_1d
!!***

!!****f* m_eliashberg_1d/eli_z_1d
!!
!! NAME
!! eli_z_1d
!!
!! FUNCTION
!!  In the solving of the 1D (energy only) Eliashberg equations, calculate
!!  the Z function, which is the renormalization factor. See Allen and Mitrovic
!!  Solid State Physics vol 37 ed Ehrenreich Seitz and Turnbull [[cite:Allen1983c]]
!!
!! INPUTS
!!   lambda_1d = coupling constant as a function of frequency
!!   nmatsu = number of Matsubara frequencies
!!
!! OUTPUT
!!   z_1d = renormalizing Z as a function of frequency
!!
!! PARENTS
!!      m_eliashberg_1d
!!
!! CHILDREN
!!
!! NOTES
!!  Because Z only depends on lambda(n-n'), and lambda(omega)
!!   is an even function, Z is symmetrical in n and -n
!!   hence only calculate for n>0 and complete the rest
!!
!! SOURCE

subroutine eli_z_1d (lambda_1d,nmatsu,z_1d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
 real(dp),intent(out) :: z_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,jmatsu

! *********************************************************************


 do imatsu=0,nmatsu

   z_1d(imatsu) = zero
!  count $\mathrm{sign}(omega_{Matsubara})$
   do jmatsu=-nmatsu+imatsu,-1
     z_1d(imatsu) = z_1d(imatsu) - lambda_1d(imatsu-jmatsu)
   end do
   do jmatsu=0,nmatsu
     z_1d(imatsu) = z_1d(imatsu) + lambda_1d(imatsu-jmatsu)
   end do

!  NOTE: the pi*Tc factor in Z cancels the one in the Matsubara frequency.
   z_1d(imatsu) = one + z_1d(imatsu) / (two*imatsu+one)
   z_1d(-imatsu) = z_1d(imatsu)
 end do

end subroutine eli_z_1d
!!***

end module m_eliashberg_1d
!!***
