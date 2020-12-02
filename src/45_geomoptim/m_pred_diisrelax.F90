!!****m* ABINIT/m_pred_diisrelax
!! NAME
!!  m_pred_diisrelax
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, JCC, SE)
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

module m_pred_diisrelax

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist
 use m_linalg_interfaces

 use m_geometry,  only : xcart2xred, xred2xcart
 use m_bfgs, only : hessinit, hessupdt

 implicit none

 private
!!***

 public :: pred_diisrelax
!!***

contains
!!***

!!****f* ABINIT/pred_diisrelax
!! NAME
!! pred_diisrelax
!!
!! FUNCTION
!! Ionmov predictor (20) Direct inversion of the iterative subspace
!!
!! IONMOV 20:
!! Given a starting point xred that is a vector of length 3*natom
!! (reduced nuclei coordinates), and unit cell parameters (rprimd)
!! this routine uses the DIIS (direct inversion of the iterative
!! subspace) to minize the gradient (forces) on atoms. The preconditioning
!! used to compute errors from gradients is using an inversed hessian
!! matrix obtained by a BFGS algorithm.
!! This method is known to converge to the nearest point where gradients
!! vanish. This is efficient to refine positions around a saddle point
!! for instance.
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! itime  : Index of the present iteration
!! ntime  : Maximal number of iterations
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!!
!! PARENTS
!!      m_precpred_1geo
!!
!! CHILDREN
!!      dcopy,dgemv,dsysv,hessinit,hessupdt,hist2var,var2hist,xcart2xred
!!      xred2xcart
!!
!! SOURCE

subroutine pred_diisrelax(ab_mover,hist,itime,ntime,zDEBUG,iexit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(inout),target :: hist

!Local variables-------------------------------
!scalars
 integer  :: ihist_prev,ndim,nhist,shift,diisSize
 integer  :: ii,jj,kk,info
 real(dp) :: etotal
 real(dp) :: suma
!arrays
 real(dp) :: acell(3)
 real(dp) :: rprimd(3,3)
!real(dp) :: fred_corrected(3,ab_mover%natom)
 real(dp) :: xred(3,ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom)
 real(dp) :: strten(6)
 real(dp) ::  ident(3,3)
 real(dp),allocatable,save :: hessin(:,:)
! DIISRELAX SPECIFIC
! error:          Store the supposed error
!                 steps. it is required to compute the DIIS matrix.
! diisMatrix:     Store the matrix used to compute the coefficients.
! diisCoeff:      Store the coefficients computed from diisMatrix.
! workMatrix:     Lapack work array.
! workArray:      Lapack work array.
! ipiv:           Lapack work array.
 integer,  allocatable :: ipiv(:)
 real(dp) :: fcart_tmp(3*ab_mover%natom)
 real(dp) :: error_tmp(3*ab_mover%natom)
 real(dp) :: xcart_tmp(3*ab_mover%natom)
 real(dp), allocatable,save :: error(:, :, :)
 real(dp), allocatable :: xcart_hist(:,:,:)
 real(dp), allocatable :: diisMatrix(:, :)
 real(dp), allocatable :: diisCoeff(:)
 real(dp), allocatable :: workArray(:)
 real(dp), allocatable :: workMatrix(:, :)
 real(dp),pointer :: fcart_hist(:,:,:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   if (allocated(ipiv))         then
     ABI_DEALLOCATE(ipiv)
   end if
   if (allocated(error))        then
     ABI_DEALLOCATE(error)
   end if
   if (allocated(diisMatrix))   then
     ABI_DEALLOCATE(diisMatrix)
   end if
   if (allocated(diisCoeff))    then
     ABI_DEALLOCATE(diisCoeff)
   end if
   if (allocated(workArray))    then
     ABI_DEALLOCATE(workArray)
   end if
   if (allocated(workMatrix))   then
     ABI_DEALLOCATE(workMatrix)
   end if
   if (allocated(hessin))       then
     ABI_DEALLOCATE(hessin)
   end if
   return
 end if

!write(std_out,*) 'diisrelax 01'
!##########################################################
!### 01. Debugging and Verbose

 if(zDEBUG)then
   write(std_out,'(a,3a,40a,37a)') ch10,('-',kk=1,3),&
&   'Debugging and Verbose for pred_diisrelax',('-',kk=1,37)
   write(std_out,*) 'ionmov: ',20
   write(std_out,*) 'itime:  ',itime
 end if

!write(std_out,*) 'diisrelax 02'
!##########################################################
!### 02. Compute the dimension of vectors (ndim=3*natom)

 ndim=3*ab_mover%natom
 nhist=hist%mxhist

 if(zDEBUG) write(std_out,*) 'Dimension of vin, vout and hessian (ndim): ',ndim

!write(std_out,*) 'diisrelax 03'
!##########################################################
!### 03. Allocate the arrays

!Notice that the arrays could be allocated
!From a previous dataset with a different ndim
 if(itime==1)then
   if (allocated(error))  then
     ABI_DEALLOCATE(error)
   end if
   if (allocated(hessin))  then
     ABI_DEALLOCATE(hessin)
   end if

   ABI_ALLOCATE(error,(3, ab_mover%natom, nhist))
   ABI_ALLOCATE(hessin,(ndim,ndim))

   ident(:, :) = zero
   ident(1, 1) = -one
   ident(2, 2) = -one
   ident(3, 3) = -one
   call hessinit(ab_mover,hessin,ident,ndim,zero)

 end if

!write(std_out,*) 'diisrelax 04'
!##########################################################
!### 04. Compute the shift in the history and the size
!###     of the diisMatrix

!When itime > diismemory we need to shift the records
!in the history to obtain the right values, the variable
!'shift' contains the actual shift to be applied on each
!iteration

 shift=max(0,itime-ab_mover%diismemory)

!Initially the diisMatrix grows with the iteration itime
!(itime+1) but when it arrives diismemory, the value of the
!matrix will be fixed on (diismemory+1)

 if (ab_mover%diismemory>itime)then
   diisSize=itime
 else
   diisSize=ab_mover%diismemory
 end if


!write(std_out,*) 'diisrelax 05'
!##########################################################
!### 05. Obtain the present values from the history

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 strten(:)=hist%strten(:,hist%ihist)
 etotal   =hist%etot(hist%ihist)

 if(zDEBUG)then
   write (std_out,*) 'fcart:'
   do kk=1,ab_mover%natom
     write (std_out,*) hist%fcart(:,kk,hist%ihist)
   end do
   write (std_out,*) 'strten:'
   write (std_out,*) strten(1:3),ch10,strten(4:6)
   write (std_out,*) 'etotal:'
   write (std_out,*) etotal
 end if

!Need also history in cartesian coordinates
 ABI_ALLOCATE(xcart_hist,(3,ab_mover%natom,diisSize))
 do ii=1,diisSize
   call xred2xcart(ab_mover%natom,rprimd,xcart_hist(:,:,ii),hist%xred(:,:,ii+shift))
 end do
 fcart_hist => hist%fcart(:,:,1+shift:diisSize+shift)

!Get rid of mean force on whole unit cell, but only if no
!generalized constraints are in effect
!  call fcart2fred(hist%fcart(:,:,hist%ihist),fred_corrected,rprimd,ab_mover%natom)
!  if(ab_mover%nconeq==0)then
!    do ii=1,3
!      if (ii/=3.or.ab_mover%jellslab==0) then
!        favg=sum(fred_corrected(ii,:))/dble(ab_mover%natom)
!        fred_corrected(ii,:)=fred_corrected(ii,:)-favg
!      end if
!    end do
!  end if

!write(std_out,*) 'diisrelax 06'
!##########################################################
!### 06. Precondition the error using the hessian matrix.

!Precondition the error using the hessian matrix.
!In the quadratic approximation, we have:
!e = H^-1.g, where e is the error vectors, H the hessian
!and g the gradient.

 if(zDEBUG)then
   write (std_out,*) 'Stored xcart:'
   do ii = 1, diisSize, 1
     write (std_out,*) 'ii,diisSize,shift',ii,diisSize,shift
     do kk=1,ab_mover%natom
       write (std_out,*) xcart_hist(:,kk,ii)
     end do
   end do
   write (std_out,*) 'Stored fcart'
   do ii = 1, diisSize, 1
     write (std_out,*) 'ii,diisSize,shift',ii,diisSize,shift
     do kk=1,ab_mover%natom
       write (std_out,*) fcart_hist(:,kk,ii)
     end do
   end do
 end if

 do ii=1,diisSize,1
   fcart_tmp(:)=RESHAPE( fcart_hist(:,:,ii), (/ ndim /) )
!  *  BLAS ROUTINE LEVEL 2
!  *  DGEMV  performs one of the matrix-vector operations
!  *
!  *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!  *
!  *  where alpha and beta are scalars, x and y are vectors and A is an
!  *  m by n matrix.

!  Here we are computing:
!  error(ndim) := 1*hessin(ndim x ndim)*fcart(ndim) + 0*error(ndim)
!
   call DGEMV('N',ndim,ndim,one,hessin,&
&   ndim,fcart_tmp,1,zero,error_tmp,1)
   error(:,:,ii)=RESHAPE( error_tmp, (/ 3, ab_mover%natom /) )

   if(zDEBUG)then
     write (std_out,*) 'Precondition ',ii
     write (std_out,*) 'fcart_tmp:'
     do kk=1,3*ab_mover%natom
       write (std_out,*) fcart_tmp(kk)
     end do
     write (std_out,*) 'error:'
     do kk=1,ab_mover%natom
       write (std_out,*) error(:,kk,ii)
     end do
     write(std_out,*) 'Hessian matrix (hessin):',ndim,'x',ndim
     do kk=1,ndim
       do jj=1,ndim,3
         if (jj+2<=ndim)then
           write(std_out,'(I3,1p,3e22.14)') jj,hessin(jj:jj+2,kk)
         else
           write(std_out,'(I3,1p,3e22.14)') jj,hessin(jj:ndim,kk)
         end if
       end do
     end do
   end if

 end do

 if(zDEBUG)then
   write (std_out,*) 'Computed error'
   do ii = 1, diisSize, 1
     write (std_out,*) 'ii,diisSize,shift',ii,diisSize,shift
     do kk=1,ab_mover%natom
       write (std_out,*) error(:,kk,ii+shift)
     end do
   end do
 end if

!write(std_out,*) 'diisrelax 07'
!##########################################################
!### 07. Create the DIIS Matrix

 ABI_ALLOCATE(diisMatrix,(diisSize + 1, diisSize + 1))
 diisMatrix(:,:) = zero
 if(zDEBUG) write(std_out,*) "DIIS matrix", diisSize+1,'x',diisSize+1
 do ii = 1, diisSize, 1
   do jj = ii, diisSize, 1
     diisMatrix(jj, ii) = ddot(ndim, error(1,1,ii),&
&     1, error(1,1,jj),1)
     diisMatrix(ii, jj) = diisMatrix(jj, ii)
   end do
   diisMatrix(ii, diisSize + 1) = -one
   diisMatrix(diisSize + 1, ii) = -one
   if(zDEBUG) write(std_out,*) diisMatrix(1:diisSize + 1, ii)
 end do

!write(std_out,*) 'diisrelax 08'
!##########################################################
!### 08. Solve the system using Lapack

 ABI_ALLOCATE(diisCoeff,(diisSize + 1))
 diisCoeff(:) = zero
 diisCoeff(diisSize + 1) = -one
 if(zDEBUG) write(std_out,*) "B vector:", diisCoeff(1:diisSize + 1)
 ABI_ALLOCATE(workMatrix,(diisSize + 1, diisSize + 1))
 ABI_ALLOCATE(workArray,((diisSize + 1) ** 2))
 ABI_ALLOCATE(ipiv,(diisSize + 1))
!*     DCOPY(N,DX,INCX,DY,INCY)
!*     copies a vector, x, to a vector, y.
!*     uses unrolled loops for increments equal to one.
 call DCOPY((diisSize + 1) ** 2, diisMatrix(1:diisSize + 1, 1:diisSize + 1), 1, workMatrix, 1)

!*     DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
!*
!*  Purpose
!*  =======
!*
!*  DSYSV computes the solution to a real system of linear equations
!*     A * X = B,
!*  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!*  matrices.
!*
!*  The diagonal pivoting method is used to factor A as
!*     A = U * D * U**T,  if UPLO = 'U', or
!*     A = L * D * L**T,  if UPLO = 'L',
!*  where U (or L) is a product of permutation and unit upper (lower)
!*  triangular matrices, and D is symmetric and block diagonal with
!*  1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!*  used to solve the system of equations A * X = B.
!*
!*  Arguments
!*  =========
!*
!*  UPLO    (input) CHARACTER*1
!*          = 'U':  Upper triangle of A is stored;
!*          = 'L':  Lower triangle of A is stored.
!*
!*  N       (input) INTEGER
!*          The number of linear equations, i.e., the order of the
!*          matrix A.  N >= 0.
!*
!*  NRHS    (input) INTEGER
!*          The number of right hand sides, i.e., the number of columns
!*          of the matrix B.  NRHS >= 0.
!*
!*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!*          N-by-N upper triangular part of A contains the upper
!*          triangular part of the matrix A, and the strictly lower
!*          triangular part of A is not referenced.  If UPLO = 'L', the
!*          leading N-by-N lower triangular part of A contains the lower
!*          triangular part of the matrix A, and the strictly upper
!*          triangular part of A is not referenced.
!*
!*          On exit, if INFO = 0, the block diagonal matrix D and the
!*          multipliers used to obtain the factor U or L from the
!*          factorization A = U*D*U**T or A = L*D*L**T as computed by
!*          DSYTRF.
!*
!*  LDA     (input) INTEGER
!*          The leading dimension of the array A.  LDA >= max(1,N).
!*
!*  IPIV    (output) INTEGER array, dimension (N)
!*          Details of the interchanges and the block structure of D, as
!*          determined by DSYTRF.  If IPIV(k) > 0, then rows and columns
!*          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
!*          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
!*          then rows and columns k-1 and -IPIV(k) were interchanged and
!*          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
!*          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
!*          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
!*          diagonal block.
!*
!*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!*          On entry, the N-by-NRHS right hand side matrix B.
!*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!*
!*  LDB     (input) INTEGER
!*          The leading dimension of the array B.  LDB >= max(1,N).
!*
!*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!*
!*  LWORK   (input) INTEGER
!*          The length of WORK.  LWORK >= 1, and for best performance
!*          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!*          DSYTRF.
!*
!*          If LWORK = -1, then a workspace query is assumed; the routine
!*          only calculates the optimal size of the WORK array, returns
!*          this value as the first entry of the WORK array, and no error
!*          message related to LWORK is issued by XERBLA.
!*
!*  INFO    (output) INTEGER
!*          = 0: successful exit
!*          < 0: if INFO = -i, the i-th argument had an illegal value
!*          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!*               has been completed, but the block diagonal matrix D is
!*               exactly singular, so the solution could not be computed.
!*
!*  =====================================================================

 if(zDEBUG)then
   write(std_out,*) "(A*X=B) A Matrix:", diisSize + 1,'x',diisSize + 1
   do ii=1,diisSize + 1
     write(std_out,*) workMatrix(1:diisSize + 1,ii)
   end do

   write(std_out,*) "(A*X=B) B Vector:", diisSize + 1
   do ii=1,diisSize + 1
     write(std_out,*) diisCoeff(ii)
   end do
 end if

 call DSYSV('L', diisSize + 1, 1, workMatrix, &
& diisSize + 1, ipiv, diisCoeff, diisSize + 1, &
& workArray, (diisSize + 1) ** 2, info)

 if (info /= 0) then
   write(std_out,*) "error solving DIIS matrix", info
   do ii=1,diisSize + 1
     write(std_out,*) workMatrix(1:diisSize + 1,ii)
   end do
 end if

 if(zDEBUG)then
   write(std_out,*) "(A*X=B) X Vector:",diisSize+1
   suma=0.0
   do ii=1,diisSize+1
     write(std_out,*) ii,diisCoeff(ii)
     suma=suma+diisCoeff(ii)
   end do

   suma=suma-diisCoeff(diisSize+1)
   write(std_out,*) 'Sum of coefficients=',suma
 end if

 ABI_DEALLOCATE(ipiv)
 ABI_DEALLOCATE(workArray)
 ABI_DEALLOCATE(workMatrix)
 ABI_DEALLOCATE(diisMatrix)

!write(std_out,*) 'diisrelax 09'
!##########################################################
!### 09. Build the new coordinates

!Build the new coordinates, to do it, we compute a new error e,
!using the linear coefficients (temporary store it in error) applied
!on previous gradient: e=H^-1(sum_i c_i.g_i)
 xcart(:, :) = zero
 error(:, :, diisSize) = zero
 do ii = 1, diisSize, 1
   xcart(:, :) = xcart(:, :) +&
&   xcart_hist(:, :, ii) * diisCoeff(ii)
   error(:, :, diisSize) = error(:, :, diisSize)+&
&   fcart_hist(:, :, ii) * diisCoeff(ii)

   if(zDEBUG)then
     write (std_out,*) 'Building new coordinates (ii):',ii
     write (std_out,*) 'diisCoeff(ii)',diisCoeff(ii)
     write (std_out,*) 'xcart_hist(:, :, ii)'
     do kk=1,ab_mover%natom
       write (std_out,*) xcart_hist(:,kk,ii)
     end do
     write (std_out,*) 'fcart_hist(:, :, ii)'
     do kk=1,ab_mover%natom
       write (std_out,*) fcart_hist(:,kk,ii)
     end do
     write (std_out,*) 'xcart:'
     do kk=1,ab_mover%natom
       write (std_out,*) xcart(:,kk)
     end do
     write (std_out,*) 'error:'
     do kk=1,ab_mover%natom
       write (std_out,*) error(:,kk,diisSize)
     end do
   end if

 end do
 ABI_DEALLOCATE(diisCoeff)

 error_tmp(:)=RESHAPE( error(:,:,diisSize), (/ 3*ab_mover%natom /) )
 xcart_tmp(:)=RESHAPE( xcart(:,:), (/ 3*ab_mover%natom /) )

!*  BLAS ROUTINE LEVEL 2
!*  DGEMV  performs one of the matrix-vector operations
!*
!*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!*
!*  where alpha and beta are scalars, x and y are vectors and A is an
!*  m by n matrix.

!Here we are computing:
!xcart_tmp(ndim) := -1*hessin(ndim x ndim)*error_tmp(ndim) + 1*xcart_tmp(ndim)
!
 call DGEMV('N', ndim, ndim, -one , hessin, &
& ndim, error_tmp, 1, one, xcart_tmp, 1)
 xcart(:,:)=RESHAPE( xcart_tmp(:), (/ 3, ab_mover%natom /) )

!write(std_out,*) 'diisrelax 10'
!##########################################################
!### 10. Update the hessian matrix using a BFGS algorithm.
 if (itime > 1) then
   call hessupdt(hessin, ab_mover%iatfix, ab_mover%natom, ndim, &
&   reshape(xcart_hist(:,:,diisSize)  , (/ ndim /)), &
&   reshape(xcart_hist(:,:,diisSize-1), (/ ndim /)), &
&   reshape(fcart_hist(:,:,diisSize)  , (/ ndim /)), &
&   reshape(fcart_hist(:,:,diisSize-1), (/ ndim /)))
 end if

 ABI_DEALLOCATE(xcart_hist)

!write(std_out,*) 'diisrelax 11'
!##########################################################
!### 11. Update the history with the prediction

!Increase indexes
 hist%ihist = abihist_findIndex(hist,+1)

!Compute xred from xcart and rprimd
 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)

!Fill the history with the variables
!xred, acell, rprimd, vel
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 ihist_prev = abihist_findIndex(hist,-1)
 hist%vel(:,:,hist%ihist)=hist%vel(:,:,ihist_prev)

 if (.false.) write(std_out,*) ntime

end subroutine pred_diisrelax
!!***

end module m_pred_diisrelax
!!***
