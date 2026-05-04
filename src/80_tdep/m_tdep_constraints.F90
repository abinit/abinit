
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_constraints

  use defs_basis
  use m_errors
  use m_xmpi
 
  implicit none
 
  type S_product
    double precision, allocatable :: SS  (:,:,:)
    double precision, allocatable :: SSS (:,:,:,:)
    double precision, allocatable :: SSSS(:,:,:,:,:)
  end type S_product
 
  type Asr_Rot
    double precision, allocatable :: ABG (:,:,:)
    double precision, allocatable :: ABGD(:,:,:,:)
    double precision, allocatable :: ABGDE(:,:,:,:,:)
  end type Asr_Rot
 
  type,public :: Constraints_type
    type(S_product),allocatable :: Sprod(:,:)
    type(Asr_Rot),allocatable :: AsrRot3(:,:,:)
    type(Asr_Rot),allocatable :: AsrRot4(:,:,:,:)
  end type Constraints_type
 
  public :: tdep_calc_orthonorm

contains

!====================================================================================================

subroutine tdep_calc_orthonorm(dim1,dim2,nindep,vect)

  integer, intent(in) :: dim1,dim2
  integer, intent(out) :: nindep
  double precision, intent(inout) :: vect(dim1,dim2)

  integer :: ii,jj,kk
  double precision :: prod_scal

! Filter non-zero vectors
  ii=0
  do kk=1,dim2
    if (sum(abs(vect(:,kk))).gt.tol8) then
      ii=ii+1
      vect(:,ii)=vect(:,kk)
    end if
  end do
  nindep=ii

! Gram-Schmidt orthogonalization
  do kk=2,nindep
    do jj=1,kk-1
      prod_scal=sum(vect(:,jj)*vect(:,jj))
      if (abs(prod_scal).gt.tol8) then
        vect(:,kk)=vect(:,kk)-sum(vect(:,kk)*vect(:,jj))/prod_scal*vect(:,jj)
      end if
    end do
  end do

! Store the non-zero vectors and normalize
  ii=0
  do kk=1,nindep
    prod_scal=sum(vect(:,kk)*vect(:,kk))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      vect(:,ii)=vect(:,kk)/dsqrt(prod_scal)
    end if
  end do
  do kk=nindep+1,dim2
    vect(:,kk)=zero
  end do  
  nindep=ii

end subroutine tdep_calc_orthonorm

!====================================================================================================

end module m_tdep_constraints
