!!****m* ABINIT/m_gwls_lineqsolver
!! NAME
!! m_gwls_lineqsolver
!!
!! FUNCTION
!!  .
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (JLJ, BR, MC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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



!---------------------------------------------------------------------
!  Module to solve A.x = b efficiently, where A will involve
!  the Hamiltonian.
!---------------------------------------------------------------------


module m_gwls_lineqsolver
!----------------------------------------------------------------------------------------------------
! This module contains routines to solve A x = b interatively to solve the Sternheimer equation
! in various contexts...
!----------------------------------------------------------------------------------------------------


! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_hamiltonian

! abinit modules
use defs_basis
use m_abicore
use m_xmpi
use m_bandfft_kpt
use m_cgtools

use m_time,      only : timab
use m_io_tools,  only : get_unit

implicit none
save
private
!!***

logical :: activate_inf_shift_poles = .false.
real(dp) :: inf_shift_poles = 1.0d-4
!!***

public :: sqmr, qmr, activate_inf_shift_poles, inf_shift_poles
!!***

contains 

!!****f* m_hamiltonian/sqmr
!! NAME
!!  sqmr
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_polarisability
!!
!! CHILDREN
!!      hpsikc,precondition_cplx,unset_precondition,xmpi_sum
!!
!! SOURCE

subroutine sqmr(b,x,lambda,project_on_what,omega,omega_imaginary,kill_Pc_x)
!--------------------------------------------------------------------------------
! This subroutine solves the linear algebra problem
!
!                                 A x = b
! 
! Where:
!                INPUT
!                -----
!        real(dp) b                          right-hand-side of the equation to be solved
!     real(dp) omega                         *OPTIONAL* frequency used in building A
!     logical  omega_imaginary               *OPTIONAL* is the frequency imaginary?
!     real(dp) lambda                        value to be subtracted from the Hamiltonian
!     integer  project_on_what               flag which determines the projection scheme.
!
!                OUTPUT
!                -----
!        real(dp) x                          solution
!
!    Note that blocksize corresponds to the number of band processors; it is a global
!    variable defined in gwls_hamiltonian. The name is inspired from lobpcgwf.F90.
!
! with:
!        omega     omega_imaginary      Operator
!        ------------------------------------------------------
!     absent          -            A =   (H - lambda)  
!     present         -            A =   (H - lambda)^2 - omega^2
!     present    present, true     A =   (H - lambda)^2 + omega^2
!
!        project_on_what                        action
!        ------------------------------------------------------
!                0                no projection
!                1                projection on conduction states
!                2                projection out of subspace degenerate with lambda
!                3                projection on states beyond all the states explicitly stored
! 
! NOTE: It is the developper's responsibility to apply (H-ev) on the input 
!       if the frequency is not zero.
!--------------------------------------------------------------------------------
implicit none

!External variables
real(dp), intent(in)  :: b(2,npw_g) 
real(dp), intent(in)  :: lambda
real(dp), intent(out) :: x(2,npw_g)
integer, intent(in)   :: project_on_what
real(dp), intent(in), optional :: omega
logical, optional     :: omega_imaginary, kill_Pc_x

!Local variables
real(dp) :: norm, tmp(2), residual
real(dp), allocatable :: g(:), theta(:), rho(:), sigma(:), c(:)
real(dp), allocatable :: t(:,:), delta(:,:), r(:,:), d(:,:), w(:,:), wmb(:,:)
integer :: k, l
real(dp):: signe
real(dp):: norm_Axb


real(dp):: norm_b, tol14

integer :: min_index
logical :: singular
logical :: precondition_on


integer :: pow

logical :: imaginary

integer,save :: counter = 0
integer      :: io_unit
character(128) :: filename
logical        :: file_exists
logical        :: head_node

integer      :: ierr

integer      :: mpi_communicator, mpi_rank, mpi_group



! timing
real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB

! *************************************************************************

! The processors communicate over FFT!
mpi_communicator = mpi_enreg%comm_fft

! what is the rank of this processor, within its group?
mpi_rank  = mpi_enreg%me_fft

! Which group does this processor belong to, given the communicator?
mpi_group = mpi_enreg%me_band


! Test if the input has finite norm
tol14 = 1.0D-14
tmp  = cg_zdotc(npw_g,b,b)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
norm_b = tmp(1)

if (norm_b < tol14) then
  ! Because of band parallelism, it is possible that sqmr gets a zero norm argument.
  ! A | x>  = 0 implies |x > = 0.
  x(:,:) = zero
  return
end if

GWLS_TIMAB   = 1523
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

! only the head node should write to the log file
head_node = ( mpi_rank == 0 )

!Memory allocation for local variables
ABI_ALLOCATE(g,    (nline))
ABI_ALLOCATE(theta,(nline))
ABI_ALLOCATE(rho,  (nline))
ABI_ALLOCATE(sigma,(nline))
ABI_ALLOCATE(c,    (nline))

ABI_ALLOCATE(t,    (2,npw_g))
ABI_ALLOCATE(delta,(2,npw_g))
ABI_ALLOCATE(r,    (2,npw_g))
ABI_ALLOCATE(d,    (2,npw_g))
ABI_ALLOCATE(w,    (2,npw_g))
ABI_ALLOCATE(wmb,  (2,npw_g))



!Some vectors won't be filled (first iteration missing) so it's useful to initialise them.
g        = zero
theta    = zero
rho      = zero
sigma    = zero
c        = zero
x        = zero
delta    = zero
r        = zero
d        = zero
w        = zero


! Determine if the frequency is imaginary
if (present(omega_imaginary) .and. present(omega)) then
  imaginary = omega_imaginary
else
  imaginary = .false.
end if


! Define the sign in front of (H-eig(v))**2. 
! If omega_imaginary is not given, we assume that omega is real (and sign=-1). 
if (present(omega)) then
  if ( imaginary ) then
    signe = one
  else
    signe =-one
  end if
end if


!Check for singularity problems
if (present(omega)) then
  norm      = minval(abs((eig(1:nbandv)-lambda)**2 + signe*(omega)**2))
  min_index = minloc(abs((eig(1:nbandv)-lambda)**2 + signe*(omega)**2),1)
else
  norm      = minval(abs(eig(1:nbandv)-lambda))
  min_index = minloc(abs(eig(1:nbandv)-lambda),1)
end if
singular = norm < 1.0d-12

!--------------------------------------------------------------------------------
! If the linear operator has a kernel, then the intermediate vectors obtained in 
! SQMR must be projected out of this subspace several time at each iterations, 
! otherwise SQMR is unstable. 
!
! This is true even if the seed vector has been initially projected out of this 
! subspace, since the preconditionning will re-introduce a non-zero component in 
! the subspace of the kernel of the linear operator. 
!                 ===>   Use project_on_what==2 in such cases.
!
! Here, test if the operator is singular and if we are NOT projecting out of 
! the kernel. 
!                ===>   If true, stop the code. 
!
!--------------------------------------------------------------------------------

! Quit if the operator has an uncontrolled kernel, a sign that the routine is being 
! misused by a developper...
if (singular .and. ( (project_on_what==1 .and. (min_index > nbandv)) .or. project_on_what==0 ))  then
  write(std_out,*) "ERROR - SQMR: Quasi-singuar problem treated, min. eigenvalue of A is ", norm," < 1d-12."
  write(std_out,*) "              Yet, there is no projection out of the kernel of A.                      "

  if (project_on_what==1 .and. (min_index > nbandv)) then
    write(std_out,*) " "
    write(std_out,*) "              There is a projection on the conduction states, but A is singular in this "
    write(std_out,*) "              subspace (the kernel contains state i=",min_index," > ",nbandv,"=# of valence states)."
  end if

  write(std_out,*) " "
  write(std_out,*) "              In this situation, SQMR will be unstable. Use project_on_what==2 as an   "
  write(std_out,*) "              input argument of SQMR."
  write(std_out,*) " "
  write(std_out,*) "                                      Decision taken to exit..."
  stop
end if

!--------------------------------------------------------------------------------
! Open a log file for the output of SQMR; only write if head of group!
!--------------------------------------------------------------------------------
if (head_node) then

  io_unit  = get_unit()

  write(filename,'(A,I0.4,A)') "SQMR_GROUP=",mpi_group,".log"

  inquire(file=filename,exist=file_exists)

  if (file_exists) then
    open( io_unit,file=filename,position='append',status=files_status_old)
  else
    open( io_unit,file=filename,status=files_status_new)
    write(io_unit,10) "#======================================================================================="
    write(io_unit,10) "#                                                                                       "
    write(io_unit,10) "#   This file contains information regarding the application of the SQMR scheme,        "
    write(io_unit,10) "#   for this MPI group.                                                                 "
    write(io_unit,10) "#======================================================================================="
    flush(io_unit)
  end if        

  counter = counter + 1
  write(io_unit,10) "#                                                                                       "
  write(io_unit,11) "#   Call # ", counter
  write(io_unit,12) "#                 lambda = ",lambda," Ha                                                "
  if (present(omega)) then
    write(io_unit,12) "#                 omega  = ",omega," Ha                                          "
    if (imaginary) then
      write(io_unit,10) "#                 omega is imaginary                                             "
    else
      write(io_unit,10) "#                 omega is real                                                  "
    end if
  else
    write(io_unit,10) "#                 omega is absent                                                 "
  end if

  write(io_unit,13) "#        project_on_what = ",project_on_what,"                                                 "
  write(io_unit,13) "#                                                                                               "
  if (present(omega) ) then
    if (imaginary) then
      write(io_unit,10) "#                SOLVE ((H-lambda)^2 + omega^2) x = b"
    else
      write(io_unit,10) "#                SOLVE ((H-lambda)^2 - omega^2) x = b"
    end if
  else
    write(io_unit,10) "#                SOLVE  (H-lambda) x = b"
  end if

  flush(io_unit)
end if ! head_node
!--------------------------------------------------------------------------------
! Precondition to accelerate convergence
!--------------------------------------------------------------------------------
precondition_on = .true. 
if(imaginary) then
  if(omega > 10.0_dp) then
    precondition_on = .false.
  end if
end if

! DEBUG
pow             = project_on_what
!precondition_on = .false. 

!Prepare to precondition
if (precondition_on) then
  if ( imaginary ) then
    call set_precondition(lambda,omega)
  else
    call set_precondition()
  end if
else
  call unset_precondition()
end if

!--------------------------------------------------------------------------------
!Initialisation
!--------------------------------------------------------------------------------

k = 1
l = 1
r = b

if (head_node) then
  write(io_unit,10) "# "
  write(io_unit,10) "# iteration          approximate residual"
  write(io_unit,10) "#----------------------------------------"
  flush(io_unit)
end if

do ! outer loop
call precondition(d,r)

! g(k)   = norm_k(r)
tmp    = cg_zdotc(npw_g,r,r)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
g(k)   = dsqrt(tmp(1))


!tmp    = scprod_k(r,d)
tmp    = cg_zdotc(npw_g,r,d)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
rho(k) = tmp(1)

if (head_node)  then
  write(io_unit,16) k, g(k)**2
  flush(io_unit)
end if

do ! inner loop
k=k+1
l=l+1

! Apply the A operator
if (present(omega)) then 
  call Hpsik(w,d,lambda)
  call Hpsik(w,cte=lambda)
  w = w + d*signe*omega**2
else
  call Hpsik(w,d,lambda)
end if 

! Apply projections, if requested
!if(dtset%gwcalctyp /= 2) then !This is a test to obtain the time taken by the orthos.
if(pow == 1) call pc_k_valence_kernel(w)
!if(pow == 2) call pc_k(w,eig_e=lambda)
!if(pow == 3) call pc_k(w,n=nband,above=.true.)
!end if

! Apply SQMR scheme
!tmp        = scprod_k(d,w)
tmp    = cg_zdotc(npw_g,d,w)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!

sigma(k-1) = tmp(1)
r          = r-(rho(k-1)/sigma(k-1))*w

! The following two lines must have a bug! We cannot distribute the norm this way!
! theta(k)   = norm_k(r)/g(k-1)
! call xmpi_sum(theta(k), mpi_communicator,ierr) ! sum on all processors working on FFT!
tmp    = cg_zdotc(npw_g,r,r)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
theta(k)   = dsqrt(tmp(1))/g(k-1)

c(k)       = one/dsqrt(one+theta(k)**2)
g(k)       = g(k-1)*theta(k)*c(k)
delta      = delta*(c(k)*theta(k-1))**2+d*rho(k-1)/sigma(k-1)*c(k)**2

x          = x+delta


if (head_node) then
  write(io_unit,16) k, g(k)**2
  flush(io_unit)
end if

! Test exit condition
if(g(k)**2<tolwfr .or. k>= nline) exit
!if(k>=nline) exit

! Safety test every 100 iterations, check that estimated residual is of the right order of magnitude. 
! If not, restart SQMR.
if(mod(l,100)==0) then
  if(present(omega)) then
    call Hpsik(w,x,lambda)
    call Hpsik(w,cte=lambda)
    w = w + x*signe*omega**2
  else
    call Hpsik(w,x,lambda)
  end if

  if(pow == 1) call pc_k_valence_kernel(w)
  !if(pow == 2) call pc_k(w,eig_e=lambda)
  !if(pow == 3) call pc_k(w,n=nband,above=.true.)

  !if(norm_k(w-b)**2 > 10*g(k)**2) exit
  wmb   = w-b
  tmp   = cg_zdotc(npw_g,wmb,wmb)
  call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
  if(tmp(1) > 10*g(k)**2) exit

end if

! Get ready for next cycle
call precondition(w,r)
!if(dtset%gwcalctyp /= 2) then
if(pow == 1) call pc_k_valence_kernel(w)
!if(pow == 2) call pc_k(w,eig_e=lambda)
!if(pow == 3) call pc_k(w,n=nband,above=.true.)
!end if

!tmp    = scprod_k(r,w)
tmp   = cg_zdotc(npw_g,r,w)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
rho(k) = tmp(1)
d      = w+d*rho(k)/rho(k-1)

end do ! end inner loop

! Exit condition
if(g(k)**2<tolwfr .or. k>=nline) exit
!if(k>=nline) exit

if (head_node) write(io_unit,10) "  ----     RESTART of SQMR -----"

!norm_Axb = norm_k(w-b)**2
wmb  = w-b     
tmp  = cg_zdotc(npw_g,wmb,wmb)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
norm_Axb = tmp(1)

if (head_node) then
  write(io_unit,*) "|Ax-b|^2 :",norm_Axb 
  write(io_unit,*) "g(k)^2 :",g(k)**2
  flush(io_unit)
end if

k = k+1
l = 1

! Apply the operator
if(present(omega)) then
  call Hpsik(r,x,lambda)
  call Hpsik(r,cte=lambda)
  r = r + x*signe*omega**2
else
  call Hpsik(r,x,lambda)
end if

if(pow == 1) call pc_k_valence_kernel(r)
!if(pow == 2) call pc_k(r,eig_e=lambda)
!if(pow == 3) call pc_k(r,n=nband,above=.true.)

r=b-r
end do


ktot = ktot+k
if(k >= nline .and. head_node ) then 
  write(io_unit,10) " **** Iterations were not enough to converge!  ****"
end if

if( present(kill_Pc_x) ) then
  if (.not. kill_Pc_x .and. pow == 1) call pc_k_valence_kernel(x)
end if

if( .not. present(kill_Pc_x) .and. pow == 1 ) call pc_k_valence_kernel(x)




!if(pow == 2) call pc_k(x,eig_e=lambda)
!if(pow == 3) call pc_k(x,n=nband,above=.true.)

if(present(omega)) then
  call Hpsik(w,x,lambda)
  call Hpsik(w,cte=lambda)
  w = w + x*signe*omega**2
else
  call Hpsik(w,x,lambda)
end if
if(pow == 1) call pc_k_valence_kernel(w)
!if(pow == 2) call pc_k(w,eig_e=lambda)
!if(pow == 3) call pc_k(w,n=nband,above=.true.)

!residual = norm_k(w-b)**2
wmb      = w-b
tmp      = cg_zdotc(npw_g,wmb,wmb)
call xmpi_sum(tmp, mpi_communicator,ierr) ! sum on all processors working on FFT!
residual = tmp(1)

if (head_node) then
  write(io_unit,15) "iterations            :", k
  write(io_unit,14) "tolwfr                :", tolwfr
  write(io_unit,14) "residuals (estimated) :", g(k)**2
  write(io_unit,14) "residuals : |Ax-b|^2  :", residual
  close(io_unit)
end if



ABI_DEALLOCATE(g)
ABI_DEALLOCATE(theta)
ABI_DEALLOCATE(rho)
ABI_DEALLOCATE(sigma)
ABI_DEALLOCATE(c)
ABI_DEALLOCATE(t)
ABI_DEALLOCATE(delta)
ABI_DEALLOCATE(r)
ABI_DEALLOCATE(d)
ABI_DEALLOCATE(w)
ABI_DEALLOCATE(wmb)


OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



10 format(A)
11 format(A,I8)
12 format(A,E24.16,A)
13 format(A,I2,A)
14 format(20X,A,E24.16)
15 format(20X,A,I8)
16 format(5X,I5,15X,E12.3)

end subroutine sqmr
!!***

!!****f* m_hamiltonian/qmr
!! NAME
!!  qmr
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_polarisability
!!
!! CHILDREN
!!      hpsikc,precondition_cplx,unset_precondition,xmpi_sum
!!
!! SOURCE

subroutine qmr(b,x,lambda) !,project_on_what)
!--------------------------------------------------------------------------------
! This subroutine solves the linear algebra problem
!
!                                 A x = b
!
! where A :=  (H - lambda)  can be non-hermitian. 
! Thus, complex values of lambda are allowed and 
! non-hermitian operators could be handled instead of H. 
! 
! Arguments :
!                INPUT
!                -----
!        real(dp) b(2,npw_k)                right-hand-side of the equation to be solved
!     real(dp) lambda(2)                value to be subtracted from the Hamiltonian (complex).
!     integer  project_on_what        flag which determines the projection scheme.
!
!                OUTPUT
!                -----
!        real(dp) x(2,npw_k)                solution
!
!        project_on_what                        action
!        ------------------------------------------------------
!                0                no projection
!                1                projection on conduction states
!                2                projection out of subspace degenerate with lambda
!                3                projection on states beyond all the states explicitly stored
!--------------------------------------------------------------------------------
implicit none

!External variables
real(dp), intent(in)  :: b(2,npw_k)
real(dp), intent(in)  :: lambda(2)
real(dp), intent(out) :: x(2,npw_k)
!integer, intent(in)   :: project_on_what !Unused yet, no projections done.

!Local variables
complex(dpc), allocatable :: xc(:), r(:), v(:), w(:), z(:), p(:), q(:), y(:), t(:), d(:), s(:)
complex(dpc), allocatable :: beta(:), eta(:), delta(:), epsilonn(:)
complex(dpc) :: lambdac
real(dp), allocatable :: rho(:), zeta(:), gama(:), theta(:), resid(:)
integer :: i
integer :: ierr

integer :: mpi_communicator
!logical :: precondition_on

! *************************************************************************

if(sum(b**2) < tol12) then
  x=zero
else
  
  !Allocation
  ABI_ALLOCATE(xc,(npw_k))
  ABI_ALLOCATE(r ,(npw_k))
  ABI_ALLOCATE(v ,(npw_k))
  ABI_ALLOCATE(w ,(npw_k))
  ABI_ALLOCATE(z ,(npw_k))
  ABI_ALLOCATE(p ,(npw_k))
  ABI_ALLOCATE(q ,(npw_k))
  ABI_ALLOCATE(y ,(npw_k))
  ABI_ALLOCATE(t ,(npw_k))
  ABI_ALLOCATE(d ,(npw_k))
  ABI_ALLOCATE(s ,(npw_k))
  
  ABI_ALLOCATE(beta    ,(nline))
  ABI_ALLOCATE(rho     ,(nline+1))
  ABI_ALLOCATE(zeta    ,(nline+1))
  ABI_ALLOCATE(gama    ,(nline+1))
  ABI_ALLOCATE(eta     ,(nline+1))
  ABI_ALLOCATE(theta   ,(nline+1))
  ABI_ALLOCATE(delta   ,(nline))
  ABI_ALLOCATE(epsilonn,(nline))
  ABI_ALLOCATE(resid   ,(nline+1))
  
  !Initialization
  x  = zero
  xc = zero
  r  = zero
  v  = zero
  w  = zero
  z  = zero
  p  = zero
  q  = zero
  y  = zero
  t  = zero
  d  = zero
  s  = zero
  
  beta     = zero
  rho      = zero
  zeta     = zero
  gama     = zero
  eta      = zero
  theta    = zero
  delta    = zero
  epsilonn = zero
  resid    = zero
  
  call unset_precondition()
  
  
  !mpi_communicator = mpi_enreg%comm_fft
  mpi_communicator = mpi_enreg%comm_bandfft
  
  lambdac  = dcmplx(lambda(1),lambda(2))
  
  i = 1
  r = dcmplx(b(1,:),b(2,:))
  v = r
  
  rho(i) = norm_kc(v)
  call xmpi_sum(rho(i),mpi_communicator ,ierr) ! sum on all processors working on FFT!
  
  w = r
  call precondition_cplx(z,w)
  zeta(i) = norm_kc(z)
  call xmpi_sum(zeta(i),mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  gama(i) = one
  eta(i) = -one
  !theta(i) = zero
  !p = zero
  !q = zero
  
  do i=1,nline
  v = v/rho(i)
  w = w/zeta(i)
  z = z/zeta(i)
  delta(i) = scprod_kc(z,v)
  call xmpi_sum(delta(i),mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  call precondition_cplx(y,v)
  if(i/=1) then
    p = y - (zeta(i)*delta(i)/epsilonn(i-1))*p
    q = z - ( rho(i)*delta(i)/epsilonn(i-1))*q
  else
    p = y
    q = z
  end if
  call Hpsikc(t,p,lambdac)
  epsilonn(i) = scprod_kc(q,t)
  call xmpi_sum(epsilonn(i),mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  beta(i) = epsilonn(i)/delta(i)
  v = t - beta(i)*v
  rho(i+1) = norm_kc(v)
  call xmpi_sum(rho(i+1),mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  call Hpsikc(z,q,conjg(lambdac)) ; w = z - beta(i)*w
  call precondition_cplx(z,w)
  zeta(i+1) = norm_kc(z)
  call xmpi_sum(zeta(i+1), mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  theta(i+1) = rho(i+1)/(gama(i)*abs(beta(i)))
  gama(i+1) = 1./sqrt(1+theta(i+1)**2)
  eta(i+1) = -eta(i)*rho(i)*gama(i+1)**2/(beta(i)*gama(i)**2)
  d = eta(i+1)*p + ((theta(i)*gama(i+1))**2)*d
  s = eta(i+1)*t + ((theta(i)*gama(i+1))**2)*s
  xc = xc + d
  r = r - s
  resid(i) = norm_kc(r)**2
  call xmpi_sum(resid(i), mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  !write(std_out,*) "QMR residual**2 = ",resid(i),"; i = ",i-1
  if(resid(i) < tolwfr) exit
  end do
  
  if(i>=nline) then
    write(std_out,*) " **** Iterations were not enough to converge!  ****"
  end if
  
  call Hpsikc(r,xc,lambdac)
  v = dcmplx(b(1,:),b(2,:))
  resid(nline+1) = norm_kc(r - v)**2
  call xmpi_sum(resid(nline+1),  mpi_communicator,ierr) ! sum on all processors working on FFT!
  
  write(std_out,*) "QMR residual**2 (at end) = ",resid(nline+1),"; # iterations = ",i-1
  
  x(1,:) = dble(xc)
  x(2,:) = dimag(xc)
  
  !Deallocate
  ABI_DEALLOCATE(xc) 
  ABI_DEALLOCATE(r) 
  ABI_DEALLOCATE(v)
  ABI_DEALLOCATE(w)
  ABI_DEALLOCATE(z)
  ABI_DEALLOCATE(p)
  ABI_DEALLOCATE(q)
  ABI_DEALLOCATE(y)
  ABI_DEALLOCATE(t)
  ABI_DEALLOCATE(d)
  ABI_DEALLOCATE(s)
  
  ABI_DEALLOCATE(beta    )  
  ABI_DEALLOCATE(rho     )
  ABI_DEALLOCATE(zeta    )
  ABI_DEALLOCATE(gama    )
  ABI_DEALLOCATE(eta     )
  ABI_DEALLOCATE(theta   )
  ABI_DEALLOCATE(delta   )
  ABI_DEALLOCATE(epsilonn)
  ABI_DEALLOCATE(resid   ) 

end if

end subroutine qmr
!!***

end module m_gwls_lineqsolver
!!***
