!!****m* ABINIT/m_gwls_wf
!! NAME
!! m_gwls_wf
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
!  Modules to handle low-level Abinit entities, like the Hamiltonian
!  and wavefunctions.
!---------------------------------------------------------------------

module m_gwls_wf

! local modules
use m_gwls_utility

! abinit modules
use defs_basis
use m_abicore
use m_cgtools
use m_xmpi

implicit none
save
private
!!***

real(dp) :: dv=0.
integer :: nx=0, ny=0, nz=0, nk=0, nb=0, ng=0, i=0, cbf=1, cf=1, cb=1
!integer :: j=0, k=0
!!***

public :: set_wf, norm_k, norm_kc, scprod_kc, scprod_k
! public :: contribution, contribution_bloc
!!***
contains

!!****f* m_hamiltonian/set_wf
!! NAME
!!  set_wf
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_hamiltonian
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine set_wf(dv2,nx2,ny2,nz2,nk2,nb2,ng2,cbf2,cf2,cb2)

implicit none
real(dp), intent(in) :: dv2
integer, intent(in) :: nx2, ny2, nz2, nk2, nb2, ng2, cbf2, cf2, cb2
! *************************************************************************
dv=dv2
nx=nx2
ny=ny2
nz=nz2
nk=nk2
nb=nb2
ng=ng2
cbf=cbf2
cf=cf2
cb=cb2
end subroutine set_wf
!!***

!!****f* m_hamiltonian/norm_k
!! NAME
!!  norm_k
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

real(dp) function norm_k(v)

implicit none
real(dp), intent(in) :: v(2,nk)
! *************************************************************************
norm_k = 0.0_dp

do i=1,nk
norm_k = norm_k + v(1,i)**2 + v(2,i)**2
end do
call xmpi_sum(norm_k,cbf,i) ! sum on all processors

norm_k = dsqrt(norm_k)

end function norm_k
!!***

!!****f* m_hamiltonian/norm_kc
!! NAME
!!  norm_kc
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

real(dp) function norm_kc(v)

implicit none
complex(dpc), intent(in) :: v(nk)
! *************************************************************************
norm_kc = zero

do i=1,nk
norm_kc = norm_kc + dble(v(i))**2+dimag(v(i))**2
end do
call xmpi_sum(norm_kc,cbf,i) ! sum on all processors
norm_kc = dsqrt(norm_kc)

end function norm_kc
!!***

!!****f* m_hamiltonian/scprod_kc
!! NAME
!!  scprod_kc
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

complex(dpc) function scprod_kc(v1,v2)

implicit none
complex(dpc), intent(in) :: v1(nk), v2(nk)
! *************************************************************************
scprod_kc = zero
do i=1,nk
scprod_kc = scprod_kc + conjg(v1(i))*v2(i)
end do
call xmpi_sum(scprod_kc,cbf,i) ! sum on all processors

!scprod_kc = sum(conjg(v1)*v2)  !Eliminated to avoid functions that return vectors embeded in another function/subroutine call.

end function scprod_kc
!!***

!!****f* m_hamiltonian/contribution
!! NAME
!!  contribution
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

! real(dp) function contribution(alpha,beta,k,norm_svne)
!
!  implicit none
!  real(dp), intent(in) :: alpha(k), beta(k-1), norm_svne
!  integer, intent(in) :: k
!
!  integer :: i
!  real(dp), allocatable :: eq_lin(:,:)
!! *************************************************************************
!  ABI_ALLOCATE(eq_lin,(4,k))
!
!  !Copy the input data into 4 vectors that will be overwritten by dgtsv()
!  eq_lin = zero
!  eq_lin(1,1:k-1) = beta       !sub-diagonal (only the 1:kmax-1 elements are used)
!  eq_lin(2,:) = alpha      !diagonal
!  eq_lin(3,1:k-1) = beta       !supra-diagonal (only the 1:kmax-1 elements are used)
!  eq_lin(4,1) = 1.0        !the RHS vector to the linear equation (here, |1,0,0,...>)
!
!  !DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!  !dgtsv(matrix size, # of column of RHS, sub-diagonal elements, diagonal elements, super-diagonal elements,
!  !      RHS of equation (solution of equation at end of routine), size of RHS vector,
!  !      error message integer (0: success, -i: illegal ith argument, i: ith factorization failed))
!  call dgtsv(k,1,eq_lin(1,1:k-1),eq_lin(2,:),eq_lin(3,1:k-1),eq_lin(4,:),k,i)
!  !Do the scalar product between the <qr_k|sv|ne>=(1,0,0,...) vector and the solution x to T*x=(1,0,0,...) to obtain contribution
!  !to SEX, at order k, from orbital n.
!  !We now replace the screened coulomb interaction by the coulomb hole...
!  contribution = (eq_lin(4,1)-1.0)*norm_svne**2
!  ABI_DEALLOCATE(eq_lin)
! end function contribution
!!***

!!****f* m_hamiltonian/contribution_bloc
!! NAME
!!  contribution_bloc
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

! function contribution_bloc(alpha,beta,kmax,norm_svne,nseeds)
!
!  implicit none
!  real(dp) :: contribution_bloc(2)
!  integer, intent(in) :: kmax, nseeds
!  real(dp), intent(in) :: alpha(2,nseeds,nseeds,kmax), beta(2,nseeds,nseeds,kmax-1), norm_svne
!
!  integer :: i=0, j=0, k=0
!  integer, allocatable :: ipiv(:)
!  complex(dpc), allocatable :: a(:,:),b(:,:)  !For the (non-banded) solver of AX=B
!! *************************************************************************
!  !write(std_out,*)  "Allocating..."
!  ABI_ALLOCATE(a,(kmax*nseeds,kmax*nseeds))
!  ABI_ALLOCATE(b,(kmax*nseeds,1))
!  ABI_ALLOCATE(ipiv,(kmax*nseeds))
!  !write(std_out,*)  "Zeroing..."
!  a=zero
!  b=zero
!  ipiv=zero
!
!  !Copy the input data into 4 vectors that will be overwritten by dgtsv()
!   !do j=1,k*nseeds
!   ! do i=max(1,j-nseeds),min(k*nseeds,j+nseeds)
!   !  ab(2*nseeds+1+i-j,j) = cmplx(alpha
!   ! end do
!   !end do
!  !write(std_out,*)  "Copying alpha..."
!  do k=1,kmax
!   do j=1,nseeds
!    do i=1,nseeds
!     a((k-1)*nseeds+i,(k-1)*nseeds+j) = cmplx(alpha(1,i,j,k),alpha(2,i,j,k),dpc)
!    end do
!   end do
!  end do
!  !write(std_out,*)  "Copying beta..."
!  do k=1,kmax-1
!   do j=1,nseeds
!    do i=1,j
!      a(k*nseeds+i,(k-1)*nseeds+j) = cmplx(beta(1,i,j,k),beta(2,i,j,k),dpc)
!      a((k-1)*nseeds+i,k*nseeds+j) = cmplx(beta(1,j,i,k),-beta(2,j,i,k),dpc)
!    end do
!   end do
!  end do
!  !write(std_out,*)  "Setting RHS..."
!  b(1,1) = (1.0,0.0)
!
!  !write(std_out,*)  "Solving..."
!  call zgesv(kmax*nseeds,1,a,kmax*nseeds,ipiv,b,kmax*nseeds,i)
!  !Do the scalar product between the <qr_k|sv|ne>=(1,0,0,...) vector and the solution x to T*x=(1,0,0,...) to obtain contribution
!  !to SEX, at order k, from orbital n
!  !write(std_out,*)  "Obtaining contribution..."
!  contribution_bloc(1) = real(b(1,1))*norm_svne**2
!  contribution_bloc(2) = aimag(b(1,1))*norm_svne**2
!  !write(std_out,*)  "Deallocating..."
!  ABI_DEALLOCATE(a)
!  ABI_DEALLOCATE(b)
!  ABI_DEALLOCATE(ipiv)
! end function contribution_bloc
!!***

!!****f* m_hamiltonian/scprod_k
!! NAME
!!  scprod_k
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

function scprod_k(v1,v2)
!--------------------------------------------------------------------------------
! This function computes the inner product of two "vectors" (typically
! wavefunctions), < v1 | v2 >.
!--------------------------------------------------------------------------------
implicit none
real(dp) :: scprod_k(2)
real(dp), intent(in) :: v1(2,nk), v2(2,nk)
! *************************************************************************
scprod_k = zero

! The mitaine way
!  do i=1,nk
!   scprod_k(1) = scprod_k(1) + v1(1,i)*v2(1,i) + v1(2,i)*v2(2,i)
!   scprod_k(2) = scprod_k(2) + v1(1,i)*v2(2,i) - v1(2,i)*v2(1,i)
!  end do

! The ABINIT way
scprod_k = cg_zdotc(nk,v1,v2)

! Collect from every processor
call xmpi_sum(scprod_k,cbf,i) ! sum on all processors

end function scprod_k
!!***

end module m_gwls_wf
!!***
