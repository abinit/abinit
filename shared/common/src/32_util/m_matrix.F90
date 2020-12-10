!!****m* ABINIT/m_matrix
!! NAME
!! m_matrix
!!
!! FUNCTION
!! Module containing some function acting on a matrix 
!!  (sqrt root)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (BA)
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
 use m_abicore

 use m_hide_lapack,  only : xginv

 implicit none

 private

! public :: init_matrix         ! Main creation method
 public :: invsqrt_matrix         ! inv of Sqrt of Matrix
 public :: blockdiago_fordsyev         ! inv of Sqrt of Matrix
 public :: blockdiago_forzheev         ! inv of Sqrt of Matrix
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
!!  force_diag = 0 if it no 0 on diagonal
!!             = nb of zeros found otherwise
!!
!! PARENTS
!!
!! CHILDREN
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

 ABI_ALLOCATE(initialmatrix,(tndim,tndim))
 initialmatrix=matrix
!  == First diagonalize matrix and keep the matrix for the change of basis
 lwork=2*tndim-1
 ABI_ALLOCATE(rwork,(3*tndim-2))
 ABI_ALLOCATE(zwork,(lwork))
 ABI_ALLOCATE(eig,(tndim))
 
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
 
 
 ABI_DEALLOCATE(zwork)
 ABI_DEALLOCATE(rwork)
 if(info/=0) then
  message = 'Error in diagonalization of zmat (zheev) ! - '
  ABI_ERROR(message)
 end if

!  == Secondly Compute sqrt(diagonalized matrix)
 ABI_ALLOCATE(diag,(tndim,tndim))
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
 ABI_DEALLOCATE(eig)
! write(std_out,*) "sqrt(eig)                , diag(1,1)",sqrt(eig(1)),diag(1,1)
! write(std_out,*) "cmplx(sqrt(eig(1)),zero,dp) , diag(1,1)",cmplx(sqrt(eig(1)),zero,dp),diag(1,1)
! write(std_out,*) "sqrt(cmplx(eig(1),zero,dp)) , diag(1,1)",sqrt(cmplx(eig(1),zero,dp)),diag(1,1)

!  == Thirdly Multiply by  matrix for the change of basis
 ABI_ALLOCATE(sqrtmat,(tndim,tndim))
 ABI_ALLOCATE(zhdp2,(tndim,tndim))
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
 ABI_DEALLOCATE(diag)

!  == Forthly Compute the inverse of the square root
! call matcginv_dpc(sqrtmat,tndim,tndim)
 !call xginv(sqrtmat,tndim)
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
!! PARENTS
!!
!! CHILDREN
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
 ABI_ALLOCATE(work,(lwork))
 work = zero

 ABI_ALLOCATE(matrix_save,(tndim,tndim))
 matrix_save=matrix
 
 ABI_ALLOCATE(Permutcol,(tndim,tndim))

 Permutcol=zero
 do im1=1,tndim
   Permutcol(im1,im1)=1.d0
 end do

 prtopt=0

 ABI_ALLOCATE(nonnul,(tndim))
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

 ABI_ALLOCATE(Apermutcol,(tndim,tndim))
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
 


 ABI_ALLOCATE(Permutline,(tndim,tndim))
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
 ABI_ALLOCATE(Apermutline,(tndim,tndim))
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
 ABI_ALLOCATE(nonnuldege,(tndim))
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
 ABI_DEALLOCATE(nonnuldege)
 if(prtopt==1) then
   write(std_out,*) "Ordered Eigenvectors"
   do im1=1,tndim
      write(std_out,'(2(1x,18(1x,f20.15,f20.15)))') (matrix(im1,im2),im2=1,tndim)
   end do
 endif

 if(prtopt==1) then
   write(std_out,*) "inverse operation: reconstitute original matrix: only the line here"
 endif
 ABI_ALLOCATE(Apermutlineback,(tndim,tndim))
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


 ABI_DEALLOCATE(Apermutlineback)
 ABI_DEALLOCATE(Apermutline)
 ABI_DEALLOCATE(matrix_save)
 ABI_DEALLOCATE(Apermutcol)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(Permutcol)
 ABI_DEALLOCATE(nonnul)
 ABI_DEALLOCATE(Permutline)

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
!! PARENTS
!!
!! CHILDREN
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
 ABI_ALLOCATE(work,(lwork))
 ABI_ALLOCATE(rwork,(3*tndim-2))

 ABI_ALLOCATE(matrix_save,(tndim,tndim))
 matrix_save=matrix
 
 ABI_ALLOCATE(Permutcol,(tndim,tndim))

 Permutcol=zero
 do im1=1,tndim
   Permutcol(im1,im1)=1.d0
 end do
 write(std_out,*) "MATRIX"
 do im1=1,tndim
    write(std_out,'(2(1x,30(1x,f22.18,f22.18)))') (matrix_save(im1,im2),im2=1,tndim)
 end do

 ABI_ALLOCATE(nonnul,(tndim))
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

 ABI_ALLOCATE(Apermutcol,(tndim,tndim))
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
 


 ABI_ALLOCATE(Permutline,(tndim,tndim))
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
 ABI_ALLOCATE(Apermutline,(tndim,tndim))
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
 ABI_ALLOCATE(Apermutlineback,(tndim,tndim))
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


 ABI_DEALLOCATE(Apermutlineback)
 ABI_DEALLOCATE(Apermutline)
 ABI_DEALLOCATE(matrix_save)
 ABI_DEALLOCATE(Apermutcol)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(Permutcol)
 ABI_DEALLOCATE(nonnul)
 ABI_DEALLOCATE(Permutline)
 ABI_DEALLOCATE(rwork)

 DBG_EXIT("COLL")

end subroutine blockdiago_forzheev
!!***

END MODULE m_matrix
