!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_matrix
!! NAME
!! m_matrix
!!
!! FUNCTION
!! Module containing some function acting on a matrix 
!!  (sqrt root)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

MODULE m_matrix

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_abilasi,  only : xginv

 implicit none

 private

! public :: init_matrix         ! Main creation method
 public :: invsqrt_matrix         ! inv of Sqrt of Matrix
! public :: inverse_matrix      ! Inverse matrix
! public :: nullify_matrix      ! Nullify the object
! public :: destroy_matrix      ! Frees the allocated memory
! public :: print_matrix        ! Printout of the basic info


CONTAINS  !===========================================================

!! FUNCTION
!!  Initialize matrix
!!
!! INPUTS
!!  ndim = dimension of matrix
!!  matrix= matrix
!!
!! OUTPUT
!!  matrix= square root of the matrix
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine invsqrt_matrix(matrix,tndim)


 use defs_basis
! use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invsqrt_matrix'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tndim 
 complex(dpc),intent(inout) :: matrix(tndim,tndim)
!arrays

!Local variables-------------------------------
!scalars
 integer :: im,im1,im2,info,lwork
 character(len=500) :: message
 real(dp) :: pawprtvol
!arrays
 real(dp),allocatable :: eig(:),rwork(:)
 complex(dpc),allocatable :: zwork(:),diag(:,:)
 complex(dpc),allocatable :: sqrtmat(:,:),zhdp2(:,:),sqrtmatinv(:,:)
 complex(dpc),allocatable :: initialmatrix(:,:)
 
! *************************************************************************

!Do not remove this silly print instruction. Seems needed to avoid floating
!point exception on vm1_gcc51 ...
#if __GFORTRAN__ == 1 && __GNUC__ == 5 && (__GNUC_MINOR__ == 1 || __GNUC_MINOR__ == 2)
 write(std_out,'(a)')' invsqrt_matrix at m_matrix.F90 : enter ( needed to avoid FPE with GCC5[1,2] )'
#endif

 DBG_ENTER("COLL")
 pawprtvol=2

 ABI_ALLOCATE(initialmatrix,(tndim,tndim))
 initialmatrix=matrix
!  == First diagonalize matrix and keep the matrix for the change of basis
 lwork=2*tndim-1
 ABI_ALLOCATE(rwork,(3*tndim-2))
 ABI_ALLOCATE(zwork,(lwork))
 ABI_ALLOCATE(eig,(tndim))
 
 call zheev('v','u',tndim,matrix,tndim,eig,zwork,lwork,rwork,info)
 
 
 ABI_DEALLOCATE(zwork)
 ABI_DEALLOCATE(rwork)
 if(info/=0) then
  message = 'Error in diagonalization of zmat (zheev) ! - '
  MSG_ERROR(message)
 end if

!  == Secondly Compute sqrt(diagonalized matrix)
 ABI_ALLOCATE(diag,(tndim,tndim))
 diag=czero
 do im=1,tndim

   if(eig(im)<tol16) then
     message = "  - Eigenvalues from zheev are negative or zero ! - "
     MSG_ERROR(message)
   else
     diag(im,im)=cmplx(sqrt(eig(im)),zero,kind=dp)
   endif
 enddo
 ABI_DEALLOCATE(eig)
! write(6,*) "sqrt(eig)                , diag(1,1)",sqrt(eig(1)),diag(1,1)
! write(6,*) "cmplx(sqrt(eig(1)),zero,dp) , diag(1,1)",cmplx(sqrt(eig(1)),zero,dp),diag(1,1)
! write(6,*) "sqrt(cmplx(eig(1),zero,dp)) , diag(1,1)",sqrt(cmplx(eig(1),zero,dp)),diag(1,1)

!  == Thirdly Multiply by  matrix for the change of basis
 ABI_ALLOCATE(sqrtmat,(tndim,tndim))
 ABI_ALLOCATE(zhdp2,(tndim,tndim))
 if(pawprtvol>3) then
   write(message,'(2a)') ch10,'  - sqrt(Eigenmatrix) - '
   call wrtout(std_out,message,'COLL')
   do im1=1,tndim
     write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))')&
!     write(message,'(12(1x,18(1x,"(",f20.16,",",f20.16,")")))')&
&     (diag(im1,im2),im2=1,tndim)
    call wrtout(std_out,message,'COLL')
   end do
 endif
!zgemm(A,B,C) : C = op(A) op(B)
 call zgemm('n','t',tndim,tndim,tndim,cone,diag,tndim,conjg(matrix),tndim,czero,zhdp2,tndim)
 call zgemm('n','n',tndim,tndim,tndim,cone,matrix,tndim,zhdp2,tndim,czero,sqrtmat,tndim)
! if(abs(pawprtvol)>=3) then
 if(pawprtvol>3) then
   write(message,'(3a)') ch10,"  - Sqrt root of matrix is - "
   call wrtout(std_out,message,'COLL')
   do im1=1,tndim
     write(message,'(12(1x,18(1x,"(",f20.16,",",f20.16,")")))')&
&     (sqrtmat(im1,im2),im2=1,tndim)
     call wrtout(std_out,message,'COLL')
   end do
 endif
! endif
 ABI_DEALLOCATE(diag)

!  == Forthly Compute the inverse of the square root
! call matcginv_dpc(sqrtmat,tndim,tndim)
 call xginv(sqrtmat,tndim)
 ABI_ALLOCATE(sqrtmatinv,(tndim,tndim))
 sqrtmatinv=sqrtmat
 if(pawprtvol>3) then
   write(message,'(2a)') ch10,"  - inverse Sqrt root of matrix is - "
   call wrtout(std_out,message,'COLL')
   do im1=1,tndim
     write(message,'(12(1x,18(1x,"(",f20.16,",",f20.16,")")))')&
&     (sqrtmatinv(im1,im2),im2=1,tndim)
     call wrtout(std_out,message,'COLL')
   end do
 endif
 ABI_DEALLOCATE(sqrtmat)

!  == Fifthly Check that O^{-0/5} O O{-0/5}=I
!  zgemm(A,B,C) : C = op(A) op(B)
 call zgemm('n','n',tndim,tndim,tndim,cone,initialmatrix,tndim,sqrtmatinv,tndim,czero,zhdp2,tndim)
 call zgemm('n','n',tndim,tndim,tndim,cone,sqrtmatinv,tndim,zhdp2,tndim,czero,initialmatrix,tndim)
 if(pawprtvol>3) then
   write(message,'(3a)') ch10,"  - O^{-0/5} O O^{-0/5}=I - "
   call wrtout(std_out,message,'COLL')
   do im1=1,tndim
     write(message,'(12(1x,18(1x,"(",f10.6,",",f4.1,")")))')&
!     write(message,'(12(1x,18(1x,"(",f20.16,",",f20.16,")")))')&
&     (initialmatrix(im1,im2),im2=1,tndim)
     call wrtout(std_out,message,'COLL')
   end do
 endif
 ABI_DEALLOCATE(zhdp2)
 matrix=sqrtmatinv
 ABI_DEALLOCATE(sqrtmatinv)
 ABI_DEALLOCATE(initialmatrix)

 DBG_EXIT("COLL")

end subroutine invsqrt_matrix
!!***

END MODULE m_matrix
