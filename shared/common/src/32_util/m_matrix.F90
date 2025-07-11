!!****m* ABINIT/m_matrix
!! NAME
!! m_matrix
!!
!! FUNCTION
!! Module containing some function acting on a matrix 
!!  (sqrt root)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2025 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_matrix

 use defs_basis
 use m_errors
 use m_abicore

 use m_hide_lapack,  only : xginv

 implicit none

 private

 public :: invsqrt_matrix       ! inv of Sqrt of Matrix
 public :: blockdiago_fordsyev  ! inv of Sqrt of Matrix
 public :: blockdiago_forzheev  ! inv of Sqrt of Matrix
 public :: mat33det             ! Determinant of a 3x3 matrix
 public :: mati3inv             ! Invert and transpose orthogonal 3x3 matrix of INTEGER elements.
 public :: mati3det             ! Compute the determinant of a 3x3 matrix of INTEGER elements.
 public :: matr3inv             ! Invert and TRANSPOSE general 3x3 matrix of real*8 elements.


 ! the determinant of a 3*3 matrix
 interface mat33det
    procedure  real_mat33det
    procedure  int_mat33det
 end interface mat33det


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
!!  force_diag = 0 if it no 0 on diagonal
!!             = nb of zeros found otherwise
!!
!! SOURCE

subroutine invsqrt_matrix(matrix,tndim,force_diag)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tndim 
 complex(dpc),intent(inout) :: matrix(tndim,tndim)
 integer, intent(out) :: force_diag
!arrays

!Local variables-------------------------------
!scalars
 integer :: im,im1,im2,info,lwork,nb_of_zero
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

 ABI_MALLOC(initialmatrix,(tndim,tndim))
 initialmatrix=matrix
!  == First diagonalize matrix and keep the matrix for the change of basis
 lwork=2*tndim-1
 ABI_MALLOC(rwork,(3*tndim-2))
 ABI_MALLOC(zwork,(lwork))
 ABI_MALLOC(eig,(tndim))
 
 call zheev('v','u',tndim,matrix,tndim,eig,zwork,lwork,rwork,info)
 if(pawprtvol>3) then
   write(message,'(2a)') ch10,'  - rotation matrix - '
   call wrtout(std_out,message,'COLL')
   do im1=1,tndim
     write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))')&
!     write(message,'(12(1x,18(1x,"(",f20.16,",",f20.16,")")))')&
&     (matrix(im1,im2),im2=1,tndim)
    call wrtout(std_out,message,'COLL')
   end do
 endif
 
 
 ABI_FREE(zwork)
 ABI_FREE(rwork)
 if(info/=0) then
  message = 'Error in diagonalization of zmat (zheev) ! - '
  ABI_ERROR(message)
 end if

!  == Secondly Compute sqrt(diagonalized matrix)
 ABI_MALLOC(diag,(tndim,tndim))
 diag=czero
 nb_of_zero=0
 do im=1,tndim

   if(eig(im)< -tol8) then
     message = "  - Eigenvalues from zheev are negative or zero ! - "
     write(std_out,*)
     write(std_out,*) "    Eigenvalue=",eig(im)
     write(std_out,*) "    Matrix is"
     do im1=1,tndim
       write(std_out,'(100f7.3)') (initialmatrix(im1,im2),im2=1,tndim)
     enddo
     ABI_ERROR(message)
   else if(abs(eig(im))<tol8) then
     nb_of_zero=nb_of_zero+1
   else
     diag(im,im)=cmplx(one/sqrt(eig(im)),zero,kind=dp)
   endif
 enddo
 force_diag=nb_of_zero
 ABI_FREE(eig)
! write(std_out,*) "sqrt(eig)                , diag(1,1)",sqrt(eig(1)),diag(1,1)
! write(std_out,*) "cmplx(sqrt(eig(1)),zero,dp) , diag(1,1)",cmplx(sqrt(eig(1)),zero,dp),diag(1,1)
! write(std_out,*) "sqrt(cmplx(eig(1),zero,dp)) , diag(1,1)",sqrt(cmplx(eig(1),zero,dp)),diag(1,1)

!  == Thirdly Multiply by  matrix for the change of basis
 ABI_MALLOC(sqrtmat,(tndim,tndim))
 ABI_MALLOC(zhdp2,(tndim,tndim))
 if(pawprtvol>3) then
   write(message,'(2a)') ch10,'  - 1.0/sqrt(Eigenmatrix) - '
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
   write(message,'(3a)') ch10,"  - inverse Sqrt root of matrix is - "
   call wrtout(std_out,message,'COLL')
   do im1=1,tndim
     write(message,'(12(1x,18(1x,"(",f20.16,",",f20.16,")")))')&
&     (sqrtmat(im1,im2),im2=1,tndim)
     call wrtout(std_out,message,'COLL')
   end do
 endif
! endif
 ABI_FREE(diag)

!  == Forthly Compute the inverse of the square root
! call matcginv_dpc(sqrtmat,tndim,tndim)
 !call xginv(sqrtmat,tndim)
 ABI_MALLOC(sqrtmatinv,(tndim,tndim))
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
 ABI_FREE(sqrtmat)

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
 ABI_FREE(zhdp2)
 matrix=sqrtmatinv
 ABI_FREE(sqrtmatinv)
 ABI_FREE(initialmatrix)

 DBG_EXIT("COLL")

end subroutine invsqrt_matrix
!!***

!! FUNCTION
!!  Transform matrix into block diagonal form before diagonalisation
!!
!! INPUTS
!!  ndim = dimension of matrix
!!  matrix= matrix
!!
!! OUTPUT
!!  matrix= square root of the matrix
!!
!! SOURCE

subroutine blockdiago_fordsyev(matrix,tndim,eig)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tndim
 real(dp),intent(inout) :: matrix(tndim,tndim)
 real(dp),intent(inout) :: eig(tndim)
!arrays

!Local variables-------------------------------
!scalars
 integer :: im1,im2,im3,info,lwork,im4,indice_formax,shift !im5,
 character(len=500) :: message
 real(dp):: tmpx,maxvalue
 integer(dp):: tmpi,newstarting,current_dege,prtopt
!arrays
 real(dp),allocatable :: work(:)
 real(dp),allocatable :: Permutcol(:,:)
 real(dp),allocatable :: Apermutcol(:,:)
 real(dp),allocatable :: Apermutline(:,:)
 real(dp),allocatable :: Apermutlineback(:,:)
 real(dp),allocatable :: Permutline(:,:)
 real(dp),allocatable :: matrix_save(:,:) !,W(:)
 integer,allocatable :: nonnul(:)
 integer,allocatable :: nonnuldege(:)
 logical :: testdege,swap
 
! *************************************************************************

!!!Do not remove this silly print instruction. Seems needed to avoid floating
!!!point exception on vm1_gcc51 ...
!!#if __GFORTRAN__ == 1 && __GNUC__ == 5 && (__GNUC_MINOR__ == 1 || __GNUC_MINOR__ == 2)
!! write(std_out,'(a)')' invsqrt_matrix at m_matrix.F90 : enter ( needed to avoid FPE with GCC5[1,2] )'
!!#endif
 DBG_ENTER("COLL")

 lwork=10*tndim
 ABI_MALLOC(work,(lwork))
 work = zero

 ABI_MALLOC(matrix_save,(tndim,tndim))
 matrix_save=matrix
 
 ABI_MALLOC(Permutcol,(tndim,tndim))

 Permutcol=zero
 do im1=1,tndim
   Permutcol(im1,im1)=1.d0
 end do

 prtopt=0

 ABI_MALLOC(nonnul,(tndim))
 do im1=1,tndim
   if(im1==1) nonnul(im1)=0
   if(im1>1) nonnul(im1)=nonnul(im1-1)
   do im2=1,tndim
     if (abs(matrix(im1,im2))>0.000000000001.and.im2>nonnul(im1)) then
       nonnul(im1)=nonnul(im1)+1
     !  write(std_out,*) "im2,nonnul(im1)",im2,nonnul(im1)
       ! permute
       do im3=1,tndim
         tmpx=matrix(im3,im2)
         matrix(im3,im2)=matrix(im3,nonnul(im1))
         matrix(im3,nonnul(im1))=tmpx
         tmpi=Permutcol(im3,im2)
         Permutcol(im3,im2)=Permutcol(im3,nonnul(im1))
         Permutcol(im3,nonnul(im1))=tmpi
       enddo
     elseif (abs(matrix(im1,im2))<0.000000000001) then
         matrix(im1,im2)=zero
     endif
   enddo
 enddo
 if(prtopt==1) then
   write(std_out,*) "MATRIX AFTER COLUMN PERMUT"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
   end do
   write(std_out,*) "Permutcol MATRIX AFTER"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (Permutcol(im1,im2),im2=1,tndim)
   end do
 endif

 ABI_MALLOC(Apermutcol,(tndim,tndim))
 if(prtopt==1) then
   write(std_out,*) "Check product of original matrix by permutation matrix "
 endif 
 Apermutcol=zero
 do im1=1,tndim
  do im2=1,tndim
   Apermutcol(im1,im2)=zero
   do im3=1,tndim
    Apermutcol(im1,im2)=matrix_save(im1,im3)*Permutcol(im3,im2)+Apermutcol(im1,im2)
   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
   end do
  end do
 end do
 if(prtopt==1) then
   write(std_out,*) "Asave*Permutcol"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (Apermutcol(im1,im2),im2=1,tndim)
   end do
 endif
 


 ABI_MALLOC(Permutline,(tndim,tndim))
 Permutline=zero
 do im1=1,tndim
   Permutline(im1,im1)=1.d0
 end do

 do im1=1,tndim
   if(im1==1) nonnul(im1)=0
   if(im1>1) nonnul(im1)=nonnul(im1-1)
   do im2=1,tndim
     ! write(std_out,*) "im1,im2, abs matrix(im2,im1),nonnul(im1)",im1,im2,abs(B(im1,im2)),nonnul(im1)
     if (abs(matrix(im2,im1))>0.000000000001.and.im2>nonnul(im1)) then
     !  write(std_out,*) "im2,nonnul(im1)",im2,nonnul(im1)
       nonnul(im1)=nonnul(im1)+1
     !  write(std_out,*) "im2,nonnul(im1)",im2,nonnul(im1)
       ! permute
       do im3=1,tndim
         tmpx=matrix(im2,im3)
         matrix(im2,im3)=matrix(nonnul(im1),im3)
         matrix(nonnul(im1),im3)=tmpx
         tmpi=Permutline(im2,im3)
         Permutline(im2,im3)=Permutline(nonnul(im1),im3)
         Permutline(nonnul(im1),im3)=tmpi
       enddo
     elseif (abs(matrix(im2,im1))<0.000000000001) then
         matrix(im2,im1)=zero
     endif
   enddo
 enddo
 if(prtopt==1) then
   write(std_out,*) "matrix AFTER"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
   end do
   write(std_out,*) "Permutline MATRIX AFTER"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (Permutline(im1,im2),im2=1,tndim)
   end do
 endif

 if(prtopt==1) then
   write(std_out,*) "Check product of Apermutcol matrix by permutation matrix of the line "
 endif
 ABI_MALLOC(Apermutline,(tndim,tndim))
 Apermutline=zero
 do im1=1,tndim
  do im2=1,tndim
   Apermutline(im1,im2)=zero
   do im3=1,tndim
    Apermutline(im1,im2)=Apermutcol(im3,im2)*Permutline(im1,im3)+Apermutline(im1,im2)
   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
   end do
  end do
 end do
 if(prtopt==1) then
   write(std_out,*) "Permutline*Apermutcol"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (Apermutline(im1,im2),im2=1,tndim)
   end do
 endif
 work=0.d0
 call dsyev('v','u',tndim,matrix_save,tndim,eig,work,lwork,info)
 if(info/=0) then
  message = 'Error in diagonalization of matrix (dsyev) ! - '
  ABI_ERROR(message)
 end if
 if(prtopt==1) then
   write(std_out,*) 'output',INFO
   write(std_out,*) "Eigenvalues"
   write(std_out,'(2x,20f20.15) ') (eig(im1),im1=1,tndim)
   write(std_out,*) "Eigenvectors"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f20.15,f20.15)))') (matrix_save(im1,im2),im2=1,tndim)
   end do
 endif


! call dsyev('v','u',tndim,A,LDA,W,WORKTMP,LWORK,INFO)
! write(std_out,*) "optimal lwork",worktmp(1)

 work=0.d0
 call dsyev('v','u',tndim,matrix,tndim,eig,work,lwork,info)
 if(info/=0) then
  message = 'Error in diagonalization of matrix (dsyev) ! - '
  ABI_ERROR(message)
 end if
 if(prtopt==1) then
  write(std_out,*) 'output',INFO
  write(std_out,*) "Eigenvalues"
  write(std_out,'(2x,20f20.15) ') (eig(im1),im1=1,tndim)
  write(std_out,*) "Eigenvectors"
  do im1=1,tndim
     write(std_out,'(2(1x,18(1x,f20.15,f20.15)))') (matrix(im1,im2),im2=1,tndim)
  end do
 endif


!! REORDER EIGENVECTORS
 ABI_MALLOC(nonnuldege,(tndim))
 newstarting=1
 current_dege=1
 do im4=2,tndim
  if(im4<tndim) testdege=((eig(im4)-eig(im4-1))<tol12)
  if(im4==tndim) then
   testdege=.false.
   current_dege=current_dege+1
  endif
  if(testdege) then
   current_dege=current_dege+1
  else 
   !new set of degenerate state: reorder it: put it into block diagonal
   !form for column
     if(prtopt==1) write(std_out,*) "newstarting, current_dege",newstarting, current_dege
     shift=0
     do im1=1,tndim ! balaye les premiers coefficients puis les autres

     !  if(im1==1) nonnuldege(im1)=0
     !  if(im1>1) nonnuldege(im1)=nonnuldege(im1-1)
       maxvalue=0.00000001
       swap=.false.
       do im2=newstarting+shift,newstarting+current_dege-1
         if(abs(matrix(im1,im2))>maxvalue) then
           maxvalue=abs(matrix(im1,im2))
           indice_formax=im2
           swap=.true.
         endif
       enddo
      ! found max value: permute
       if(swap) then
        do im3=1,tndim
          tmpx=matrix(im3,indice_formax)
          matrix(im3,indice_formax)=matrix(im3,newstarting+shift)
          matrix(im3,newstarting+shift)=tmpx
        enddo
        shift=shift+1
       endif
       !write(std_out,*) "Eigenvectors after m1"
       !do im3=1,tndim
       !   write(std_out,'(2(1x,18(1x,f20.15,f20.15)))') (matrix(im3,im5),im5=1,tndim)
       !end do

     enddo
     if(prtopt==1) then
       write(std_out,*) "Eigenvectors after set of dege"
       do im2=1,tndim
          write(std_out,'(2(1x,18(1x,f20.15,f20.15)))') (matrix(im2,im3),im3=1,tndim)
       end do
     endif
     newstarting=im4
     current_dege=1
  endif
 enddo
 ABI_FREE(nonnuldege)
 if(prtopt==1) then
   write(std_out,*) "Ordered Eigenvectors"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f20.15,f20.15)))') (matrix(im1,im2),im2=1,tndim)
   end do
 endif

 if(prtopt==1) then
   write(std_out,*) "inverse operation: reconstitute original matrix: only the line here"
 endif
 ABI_MALLOC(Apermutlineback,(tndim,tndim))
 Apermutlineback=zero
 do im1=1,tndim
  do im2=1,tndim
   Apermutlineback(im1,im2)=zero
   do im3=1,tndim
    Apermutlineback(im1,im2)=matrix(im3,im2)*Permutline(im3,im1)+Apermutlineback(im1,im2)
   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
   end do
  end do
 end do
 matrix=Apermutlineback
 if(prtopt==1) then
   write(std_out,*) "t(Permutline)*Apermutcol"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
   end do
 endif

! Now, set the first coefficient of eigenvectors positive.

 do im2=1,tndim ! loop over eigenvectors
   do im1=1,tndim ! loop over components
     if(abs(matrix(im1,im2))>tol8) then
       if(matrix(im1,im2)<0) then
         do im3=1,tndim
           if(abs(matrix(im3,im2))>tol8) then
             matrix(im3,im2)=-matrix(im3,im2)
           endif
         enddo
       endif
       exit
     endif
   enddo
 enddo
 if(prtopt==1) then
   write(std_out,*) "Impose first component of eigenvectors is positive"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
   end do
 endif
! write(std_out,*) "inverse operation: reconstitute original matrix: then the column"
! Apermutcolback=zero
! do im1=1,tndim
!  do im2=1,tndim
!   Apermutcolback(im1,im2)=zero
!   do im3=1,tndim
!    Apermutcolback(im1,im2)=Apermutlineback(im1,im3)*Permutcol(im2,im3)+Apermutcolback(im1,im2)
!   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
!   end do
!  end do
! end do
! write(std_out,*) "Apermutlineback*t(Permutcol)"
! do im1=1,10
!    write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (Apermutcolback(im1,im2),im2=1,10)
! end do


 ABI_FREE(Apermutlineback)
 ABI_FREE(Apermutline)
 ABI_FREE(matrix_save)
 ABI_FREE(Apermutcol)
 ABI_FREE(work)
 ABI_FREE(Permutcol)
 ABI_FREE(nonnul)
 ABI_FREE(Permutline)

 DBG_EXIT("COLL")

end subroutine blockdiago_fordsyev
!!***

!! FUNCTION
!!  Transform matrix into block diagonal form before diagonalisation
!!
!! INPUTS
!!  ndim = dimension of matrix
!!  matrix= matrix
!!
!! OUTPUT
!!  matrix= square root of the matrix
!!
!! SOURCE

subroutine blockdiago_forzheev(matrix,tndim,eig)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tndim
 complex(dpc),intent(inout) :: matrix(tndim,tndim)
 real(dp),intent(inout) :: eig(tndim)
!arrays

!Local variables-------------------------------
!scalars
 integer :: im1,im2,im3,info,lwork
 character(len=500) :: message
 complex(dpc):: tmpx
 integer(dp):: tmpi
!arrays
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)
 real(dp),allocatable :: Permutcol(:,:)
 complex(dpc),allocatable :: Apermutcol(:,:)
 complex(dpc),allocatable :: Apermutline(:,:)
 complex(dpc),allocatable :: Apermutlineback(:,:)
 real(dp),allocatable :: Permutline(:,:)
 complex(dpc),allocatable :: matrix_save(:,:) !,W(:)
 integer,allocatable :: nonnul(:)
 
! *************************************************************************

!!!Do not remove this silly print instruction. Seems needed to avoid floating
!!!point exception on vm1_gcc51 ...
!!#if __GFORTRAN__ == 1 && __GNUC__ == 5 && (__GNUC_MINOR__ == 1 || __GNUC_MINOR__ == 2)
!! write(std_out,'(a)')' invsqrt_matrix at m_matrix.F90 : enter ( needed to avoid FPE with GCC5[1,2] )'
!!#endif
 DBG_ENTER("COLL")

 lwork=10*tndim
 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(rwork,(3*tndim-2))

 ABI_MALLOC(matrix_save,(tndim,tndim))
 matrix_save=matrix
 
 ABI_MALLOC(Permutcol,(tndim,tndim))

 Permutcol=zero
 do im1=1,tndim
   Permutcol(im1,im1)=1.d0
 end do
 write(std_out,*) "MATRIX"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (matrix_save(im1,im2),im2=1,tndim)
 end do

 ABI_MALLOC(nonnul,(tndim))
 do im1=1,tndim
   if(im1==1) nonnul(im1)=0
   if(im1>1) nonnul(im1)=nonnul(im1-1)
   do im2=1,tndim
     if (abs(matrix(im1,im2))>0.000000000001.and.im2>nonnul(im1)) then
       nonnul(im1)=nonnul(im1)+1
     !  write(std_out,*) "im2,nonnul(im1)",im2,nonnul(im1)
       ! permute
       do im3=1,tndim
         tmpx=matrix(im3,im2)
         matrix(im3,im2)=matrix(im3,nonnul(im1))
         matrix(im3,nonnul(im1))=tmpx
         tmpi=Permutcol(im3,im2)
         Permutcol(im3,im2)=Permutcol(im3,nonnul(im1))
         Permutcol(im3,nonnul(im1))=tmpi
       enddo
     elseif (abs(matrix(im1,im2))<0.000000000001) then
         matrix(im1,im2)=czero
     endif
   enddo
 enddo
 write(std_out,*) "MATRIX AFTER COLUMN PERMUT"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
 end do
 write(std_out,*) "Permutcol MATRIX AFTER"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (Permutcol(im1,im2),im2=1,tndim)
 end do

 ABI_MALLOC(Apermutcol,(tndim,tndim))
 write(std_out,*) "Check product of original matrix by permutation matrix "
 Apermutcol=czero
 do im1=1,tndim
  do im2=1,tndim
   Apermutcol(im1,im2)=czero
   do im3=1,tndim
    Apermutcol(im1,im2)=matrix_save(im1,im3)*Permutcol(im3,im2)+Apermutcol(im1,im2)
   ! write(std_out,*) im1,im2,im3,Apermutcol(im1,im2)
   end do
  end do
 end do
 write(std_out,*) "Asave*Permutcol"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (Apermutcol(im1,im2),im2=1,tndim)
 end do
 


 ABI_MALLOC(Permutline,(tndim,tndim))
 Permutline=zero
 do im1=1,tndim
   Permutline(im1,im1)=1.d0
 end do

 do im1=1,tndim
   if(im1==1) nonnul(im1)=0
   if(im1>1) nonnul(im1)=nonnul(im1-1)
   do im2=1,tndim
     ! write(std_out,*) "im1,im2, abs matrix(im2,im1),nonnul(im1)",im1,im2,abs(B(im1,im2)),nonnul(im1)
     if (abs(matrix(im2,im1))>0.000000000001.and.im2>nonnul(im1)) then
     !  write(std_out,*) "im2,nonnul(im1)",im2,nonnul(im1)
       nonnul(im1)=nonnul(im1)+1
     !  write(std_out,*) "im2,nonnul(im1)",im2,nonnul(im1)
       ! permute
       do im3=1,tndim
         tmpx=matrix(im2,im3)
         matrix(im2,im3)=matrix(nonnul(im1),im3)
         matrix(nonnul(im1),im3)=tmpx
         tmpi=Permutline(im2,im3)
         Permutline(im2,im3)=Permutline(nonnul(im1),im3)
         Permutline(nonnul(im1),im3)=tmpi
       enddo
     elseif (abs(matrix(im2,im1))<0.000000000001) then
         matrix(im2,im1)=czero
     endif
   enddo
 enddo
 write(std_out,*) "matrix AFTER"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
 end do
 write(std_out,*) "Permutline MATRIX AFTER"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (Permutline(im1,im2),im2=1,tndim)
 end do

 write(std_out,*) "Check product of Apermutcol matrix by permutation matrix of the line "
 ABI_MALLOC(Apermutline,(tndim,tndim))
 Apermutline=czero
 do im1=1,tndim
  do im2=1,tndim
   Apermutline(im1,im2)=czero
   do im3=1,tndim
    Apermutline(im1,im2)=Apermutcol(im3,im2)*Permutline(im1,im3)+Apermutline(im1,im2)
   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
   end do
  end do
 end do
 write(std_out,*) "Permutline*Apermutcol"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (Apermutline(im1,im2),im2=1,tndim)
 end do
 work=czero
 call zheev('v','u',tndim,matrix_save,tndim,eig,work,lwork,rwork,info)
 if(info/=0) then
  message = 'Error in diagonalization of matrix (zheev) ! - '
  ABI_ERROR(message)
 end if
 write(std_out,*) 'output',INFO
 write(std_out,*) "Eigenvalues"
 write(std_out,'(2x,20f20.15) ') (eig(im1),im1=1,tndim)
 write(std_out,*) "Eigenvectors"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f20.15,f20.15)))') (matrix_save(im1,im2),im2=1,tndim)
 end do


! call dsyev('v','u',tndim,A,LDA,W,WORKTMP,LWORK,INFO)
! write(std_out,*) "optimal lwork",worktmp(1)

 work=czero
 call zheev('v','u',tndim,matrix,tndim,eig,work,lwork,rwork,info)
 if(info/=0) then
  message = 'Error in diagonalization of matrix (zheev) ! - '
  ABI_ERROR(message)
 end if
 write(std_out,*) 'output',INFO
 write(std_out,*) "Eigenvalues"
 write(std_out,'(2x,20f20.15) ') (eig(im1),im1=1,tndim)
 write(std_out,*) "Eigenvectors"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f20.15,f20.15)))') (matrix(im1,im2),im2=1,tndim)
 end do

 write(std_out,*) "inverse operation: reconstitute original matrix: first the line"
 ABI_MALLOC(Apermutlineback,(tndim,tndim))
 Apermutlineback=czero
 do im1=1,tndim
  do im2=1,tndim
   Apermutlineback(im1,im2)=czero
   do im3=1,tndim
    Apermutlineback(im1,im2)=matrix(im3,im2)*Permutline(im3,im1)+Apermutlineback(im1,im2)
   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
   end do
  end do
 end do
 matrix=Apermutlineback
 write(std_out,*) "t(Permutline)*Apermutcol"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (matrix(im1,im2),im2=1,tndim)
 end do

! write(std_out,*) "inverse operation: reconstitute original matrix: then the column"
! Apermutcolback=zero
! do im1=1,tndim
!  do im2=1,tndim
!   Apermutcolback(im1,im2)=zero
!   do im3=1,tndim
!    Apermutcolback(im1,im2)=Apermutlineback(im1,im3)*Permutcol(im2,im3)+Apermutcolback(im1,im2)
!   ! write(std_out,*) PROD(im1,im2),A(im3,im1),A(im3,im2)
!   end do
!  end do
! end do
! write(std_out,*) "Apermutlineback*t(Permutcol)"
! do im1=1,10
!    write(std_out,'(2(1x,18(1x,f22.18,f22.18)))') (Apermutcolback(im1,im2),im2=1,10)
! end do


 ABI_FREE(Apermutlineback)
 ABI_FREE(Apermutline)
 ABI_FREE(matrix_save)
 ABI_FREE(Apermutcol)
 ABI_FREE(work)
 ABI_FREE(Permutcol)
 ABI_FREE(nonnul)
 ABI_FREE(Permutline)
 ABI_FREE(rwork)

 DBG_EXIT("COLL")

end subroutine blockdiago_forzheev
!!***

!! FUNCTION
!!  Compute the determinant of a 3x3 real matrix
!!
!! INPUTS
!!  A = 3x3 matrix
!!
!! OUTPUT
!!  det = The determinant
!!
!! SOURCE

function real_mat33det(A) result(det)
  real(dp), intent(in) :: A(3,3)
  real(dp) :: det
  DET =  A(1,1)*A(2,2)*A(3,3)  &
       - A(1,1)*A(2,3)*A(3,2)  &
       - A(1,2)*A(2,1)*A(3,3)  &
       + A(1,2)*A(2,3)*A(3,1)  &
       + A(1,3)*A(2,1)*A(3,2)  &
       - A(1,3)*A(2,2)*A(3,1)
end function real_mat33det
!!***

!! FUNCTION
!!  Compute the determinant of a 3x3 integer matrix
!!
!! INPUTS
!!  A = 3x3 matrix
!!
!! OUTPUT
!!  det = The determinant
!!
!! SOURCE

function int_mat33det(A) result(det)
  integer, intent(in) :: A(3,3)
  integer :: det
  DET =  A(1,1)*A(2,2)*A(3,3)  &
       - A(1,1)*A(2,3)*A(3,2)  &
       - A(1,2)*A(2,1)*A(3,3)  &
       + A(1,2)*A(2,3)*A(3,1)  &
       + A(1,3)*A(2,1)*A(3,2)  &
       - A(1,3)*A(2,2)*A(3,1)
end function int_mat33det
!!***

!!****f* m_matrix/mati3inv
!! NAME
!! mati3inv
!!
!! FUNCTION
!! Invert and transpose orthogonal 3x3 matrix of INTEGER elements.
!!
!! INPUTS
!! mm = integer matrix to be inverted
!!
!! OUTPUT
!! mit = inverse of mm input matrix
!!
!! NOTES
!! Used for symmetry operations.
!! This routine applies to ORTHOGONAL matrices only.
!! Since these form a group, inverses are also integer arrays.
!! Returned array is TRANSPOSE of inverse, as needed.
!! Note use of integer arithmetic.
!!
!! SOURCE

subroutine mati3inv(mm, mit)

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: mit(3,3)

!Local variables-------------------------------
!scalars
 integer :: dd
 character(len=500) :: msg
!arrays
 integer :: tt(3,3)

! *************************************************************************

 tt(1,1) = mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)
 tt(2,1) = mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)
 tt(3,1) = mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3)
 tt(1,2) = mm(3,1) * mm(2,3) - mm(2,1) * mm(3,3)
 tt(2,2) = mm(1,1) * mm(3,3) - mm(3,1) * mm(1,3)
 tt(3,2) = mm(2,1) * mm(1,3) - mm(1,1) * mm(2,3)
 tt(1,3) = mm(2,1) * mm(3,2) - mm(3,1) * mm(2,2)
 tt(2,3) = mm(3,1) * mm(1,2) - mm(1,1) * mm(3,2)
 tt(3,3) = mm(1,1) * mm(2,2) - mm(2,1) * mm(1,2)
 dd = mm(1,1) * tt(1,1) + mm(2,1) * tt(2,1) + mm(3,1) * tt(3,1)

 ! Make sure matrix is not singular
 if (dd /= 0) then
   mit(:,:)=tt(:,:)/dd
 else
   write(msg, '(2a,2x,9(i0,1x),a)' )'Attempting to invert integer array',ch10,mm,' ==> determinant is zero.'
   ABI_ERROR(msg)
 end if

 ! If matrix is orthogonal, determinant must be 1 or -1
 if (abs(dd) /= 1) then
   write(msg, '(3a,i0)' )'Absolute value of determinant should be one',ch10,'but determinant= ',dd
   ABI_ERROR(msg)
 end if

end subroutine mati3inv
!!***

!!****f* m_matrix/mati3det
!! NAME
!! mati3det
!!
!! FUNCTION
!! Compute the determinant of a 3x3 matrix of INTEGER elements.
!!
!! INPUTS
!! mm = integer matrix
!!
!! OUTPUT
!! det = determinant of the matrix
!!
!! SOURCE

subroutine mati3det(mm, det)

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: det

! *************************************************************************
 det=mm(1,1)*(mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)) &
   + mm(2,1)*(mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)) &
   + mm(3,1)*(mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3))

end subroutine mati3det
!!***

!!****f* m_matrix/matr3inv
!! NAME
!! matr3inv
!!
!! FUNCTION
!! Invert and transpose general 3x3 matrix of real*8 elements.
!!
!! INPUTS
!! aa = 3x3 matrix to be inverted
!!
!! OUTPUT
!! ait = inverse of aa input matrix
!!
!! NOTES
!! Returned array is TRANSPOSE of inverse, as needed to get g from r.
!!
!! SOURCE

subroutine matr3inv(aa, ait)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3)
 real(dp),intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd,det,t1,t2,t3
 character(len=500) :: msg

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 det  = aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3

!Make sure matrix is not singular
 if (abs(det)>tol16) then
   dd=one/det
 else
   write(msg, '(2a,2x,9es16.8,a,a,es16.8,a)' )&
     'Attempting to invert real(8) 3x3 array',ch10,aa(:,:),ch10,'   ==> determinant=',det,' is zero.'
   ABI_BUG(msg)
 end if

 ait(1,1) = t1 * dd
 ait(2,1) = t2 * dd
 ait(3,1) = t3 * dd
 ait(1,2) = (aa(3,1)*aa(2,3)-aa(2,1)*aa(3,3)) * dd
 ait(2,2) = (aa(1,1)*aa(3,3)-aa(3,1)*aa(1,3)) * dd
 ait(3,2) = (aa(2,1)*aa(1,3)-aa(1,1)*aa(2,3)) * dd
 ait(1,3) = (aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)) * dd
 ait(2,3) = (aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)) * dd
 ait(3,3) = (aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)) * dd

end subroutine matr3inv
!!***


END MODULE m_matrix
