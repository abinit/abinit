!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_exp_mat
!! NAME
!!  m_exp_mat
!!
!! FUNCTION
!!  This subroutine calculate the exponential of a  matrix
!!
!! COPYRIGHT
!! Copyright (C) 2002-2019 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_exp_mat

 use defs_basis
 use m_abicore
 use m_errors
 use m_linalg_interfaces

 implicit none

 private

 public :: exp_mat  ! exponential of a complex matrix


 interface exp_mat
  module procedure exp_mat_cx
 end interface exp_mat


CONTAINS  !===========================================================
 !!***

 !!****f* m_exp_mat/exp_mat_cx
 !! NAME
 !! exp_mat_cx
 !!
 !! FUNCTION
 !! Returns the exponential of a complex matrix multiplied a scalar
 !!
 !! INPUTS
 !! mat_a = the complex matrix
 !! mat_a = the size of mat_a
 !! factor =  a real factor multiplying at the exponent
 !!
 !! OUTPUT
 !!  exp_mat_cx = its exponential is returned in the same matrix
 !! PARENTS
!!
 !! CHILDREN
!!      zgeev,zgetrf,zgetri
!!
 !! SOURCE

 subroutine exp_mat_cx(mat_a,mat_a_size,factor)

  !Arguments ------------------------------------
  ! scalars
  real(dp),intent(in) ::  factor
  ! arrays
  complex(dpc),intent(inout) :: mat_a(:,:)
  !complex(dpc) :: exp_mat_cx(mat_a_size,mat_a_size)

  !Local ------------------------------------------
  ! scalars
  integer :: info,mat_a_size,ii
  integer,parameter :: maxsize=3
  integer,parameter :: lwork=(1+32)*maxsize

  ! arrays
  integer :: ipvt(mat_a_size)
  character(len=500) :: msg
  real(dp) :: rwork(2*maxsize)
  complex(dpc),allocatable :: ww(:),uu(:,:)
  complex(dpc) :: work(lwork),vl(1,1)
  ! *********************************************************************

  mat_a_size = max(1,size(mat_a,dim=1))
  ABI_ALLOCATE(ww,(mat_a_size))
  ABI_ALLOCATE(uu,(mat_a_size,mat_a_size))

  !Now it calculates the eigenvalues and eigenvectors of the matrix
  call ZGEEV('No left vectors','Vectors (right)',mat_a_size, mat_a, mat_a_size,ww,&
    vl,1,uu, mat_a_size, work, lwork, rwork, info)
  if (info/=0) then
   write(msg,'(a,i4)')'Wrong value for rwork ',info
   MSG_BUG(msg)
  end if

  !!debbug
  !    write(std_out,*)'mat_a',mat_a
  !    write(std_out,*)'mat_a_size',mat_a_size

  !    write(std_out,*)'eigenvalues'
  !    write(std_out,*)ww
  !    write(std_out,*)'eigenvectors'
  !    do ii=1,g_mat_size
  !     write(std_out,*)'autov',ii
  !     write(std_out,*)uu(:,ii)
  !    end do
  !    write(std_out,*)'optimal workl=',work(1)

  !    write(std_out,*)'info', info
  !    write(std_out,*) '------------------'
  !!enddebug



  !-----------------------------------------------------------
  !Now it calculates the exponential of the eigenvectors (times the factor)

  !--exponential of the diagonal
  ww(:) = exp( ww(:)*factor )

  !--construction exponential matrix
  mat_a = zero
  mat_a(:,1) = ww(:)
  mat_a(:,:) = cshift(array=mat_a,shift=(/ (-ii,ii=0,mat_a_size) /), dim=2 )


  !uu.exp(ww*factor)
  mat_a(:,:) = matmul(uu,mat_a)

  !the inverse of the eigenvectors matrix
  call ZGETRF( mat_a_size, mat_a_size, uu,mat_a_size, ipvt, info )
  if (info/=0) then
   write(msg,'(a,i4)')'Wrong value for rwork ',info
   MSG_BUG(msg)
  end if

  call ZGETRI( mat_a_size, uu, mat_a_size, ipvt, work, lwork, info )
  if (info/=0) then
   write(msg,'(a,i4)')'Wrong value for rwork ',info
   MSG_BUG(msg)
  end if

  !(uu.exp(ww*factor)).uu-1
  mat_a = matmul(mat_a,uu)

  ABI_DEALLOCATE(ww)
  ABI_DEALLOCATE(uu)

 end subroutine exp_mat_cx
 !!***




END MODULE m_exp_mat
!!***
