!!****m* ABINIT/m_iterative_solvers
!! NAME
!!  m_iterative_solvers
!!
!! FUNCTION
!!  This module contains basic (matrix-free) iterative solvers (GMRES, CG).
!!
!! COPYRIGHT
!! TODO
!!
!! SOURCE

#include "abi_common.h"

module m_iterative_solvers

    use defs_basis
#if defined HAVE_LINALG_MKL_OMATCOPY
    use mkl_rci, only : dfgmres, dfgmres_check, dfgmres_get, dfgmres_init
#endif

    implicit none
    private
    public :: gmres_linear_solver, cg_eigen_solver_treshold
 
    contains

!-------------------------------------------------------------------------------------------
! Eigensolvers
!-------------------------------------------------------------------------------------------

    !****f* m_iterative_solvers/cg_eigen_solver_treshold
    !! NAME
    !!  cg_eigen_solver_treshold
    !!
    !! FUNCTION
    !!  Compute the smallest eigenvalues and corresponding eigenvectors of a matrix using
    !!  the Conjugate Gradient method to minimize the Rayleigh quotient. 
    !!  The stopping criterion is : at least one of the computed eigenvalues is above 'treshold'.
    !!
    !! INPUTS
    !!  n              = Size of the matrix.
    !!  matvec         = Subroutine that performs matrix-vector multiplication.
    !!  x0             = Initial guess for the eigenvector.
    !!  tol            = Convergence tolerance for the residual norm.
    !!  max_iter       = Maximum number of iterations for the Conjugate Gradient method.
    !!  max_neig       = Maximum number of eigenvalues to compute.
    !!  eigenvalue_threshold = Threshold below which to stop computing further eigenvalues.
    !!
    !! OUTPUTS
    !!  eigenvalues    = Array of computed smallest eigenvalues.
    !!  eigenvectors   = Matrix of corresponding eigenvectors.
    !!  n_eig          = Number of computed eigenvalues.
    !!
    !! SOURCE
    subroutine cg_eigen_solver_treshold(n, matvec, x0, tol, max_iter, max_neig, eigenvalue_threshold, eigenvalues, eigenvectors, n_eig)
        ! Input parameters
        integer, intent(in) :: n, max_iter, max_neig
        real(dp), intent(in) :: tol, eigenvalue_threshold
        real(dp), intent(in) :: x0(n)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface

        ! Output parameters
        real(dp), intent(out) :: eigenvalues(max_neig)
        real(dp), intent(out) :: eigenvectors(n, max_neig)
        integer, intent(out) :: n_eig

        ! Local variables
        real(dp) :: r(n), p(n), Ap(n), x(n)
        real(dp) :: residual_norm, rayleigh_quotient
        integer :: i, j, iter
        
        ! *************************************************************************
        
        ! Initialize variables
        x = x0 / sqrt(sum(x0**2))   ! Normalize the initial guess
        r = x                       ! Residual vector
        p = r                       ! Search direction
        residual_norm = sqrt(sum(r**2))

        ! Conjugate Gradient iterations to minimize Rayleigh quotient
        do iter = 1, max_iter
            call cg_update(n, matvec, x, r, p, Ap, rayleigh_quotient, residual_norm)
            if (residual_norm < tol) exit
        end do

        ! Store the computed eigenvalue and eigenvector
        eigenvalues(1) = rayleigh_quotient
        eigenvectors(:, 1) = x
        n_eig = 1

        ! If more eigenvalues are required, use deflation to compute subsequent eigenvalues
        do j = 2, max_neig
            ! Check if the eigenvalue is below the threshold
            if (eigenvalues(j-1) < eigenvalue_threshold) exit

            ! If it is not, we need to compute more eigenvalues
            ! Orthogonalize the initial guess against previously computed eigenvectors
            x = x0
            do i = 1, j - 1
                x = x - dot_product(x, eigenvectors(:, i)) * eigenvectors(:, i)
            end do
            x = x / sqrt(sum(x**2))

            ! Reset residual and search direction for the next eigenvalue
            r = x
            p = r
            residual_norm = sqrt(sum(r**2))

            do iter = 1, max_iter
                call cg_update(n, matvec, x, r, p, Ap, rayleigh_quotient, residual_norm)
                if (residual_norm < tol) exit
            end do

            eigenvalues(j) = rayleigh_quotient
            eigenvectors(:, j) = x
            n_eig = n_eig + 1
        end do

    end subroutine cg_eigen_solver_treshold

    !****f* m_iterative_solvers/cg_update
    !! NAME
    !!  cg_update
    !!
    !! FUNCTION
    !!  Perform a single Conjugate Gradient update step to minimize the Rayleigh quotient.
    !!
    !! INPUTS
    !!  n              = Size of the matrix.
    !!  matvec         = Subroutine that performs matrix-vector multiplication.
    !!
    !! INPUT/OUTPUTS
    !!  x              = Current solution vector (eigenvector approximation).
    !!  r              = Residual vector.
    !!  p              = Search direction vector.
    !!  Ap             = Result of matrix-vector multiplication (A * p).
    !!
    !! OUTPUTS
    !!  rayleigh_quotient = Current approximation of the eigenvalue.
    !!  residual_norm     = Norm of the residual vector.
    !!
    !! SOURCE
    subroutine cg_update(n, matvec, x, r, p, Ap, rayleigh_quotient, residual_norm)
        ! Arguments
        integer, intent(in) :: n
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
        real(dp), intent(inout) :: x(n), r(n), p(n), Ap(n)
        real(dp), intent(out) :: rayleigh_quotient, residual_norm

        ! Local variables
        real(dp) :: alpha, beta

        ! *************************************************************************

        ! Apply the matrix-vector multiplication
        call matvec(n, p, Ap)

        ! Compute Rayleigh quotient (approximation of eigenvalue)
        rayleigh_quotient = dot_product(x, Ap) / dot_product(x, x)

        ! Compute alpha (step size)
        alpha = dot_product(r, r) / dot_product(p, Ap)

        ! Update the solution vector
        x = x + alpha * p

        ! Update the residual vector
        r = r - alpha * Ap

        ! Compute residual norm
        residual_norm = sqrt(sum(r**2))

        ! Compute beta (update factor for search direction)
        beta = dot_product(r, r) / dot_product(r - alpha * Ap, r - alpha * Ap)

        ! Update the search direction
        p = r + beta * p

    end subroutine cg_update

!-------------------------------------------------------------------------------------------
! Linear solvers (GMRES)
!-------------------------------------------------------------------------------------------

    !****f* m_iterative_solvers/call_FGMRES
    !! NAME
    !!  call_FGMRES
    !!
    !! FUNCTION
    !!  Call the MKL FGMRES routine to solve a linear system.
    !!
    !! INPUTS
    !!  n              = Size of the matrix.
    !!  matvec         = Subroutine that performs matrix-vector multiplication.
    !!  rhs            = Right-hand side vector of the linear system.
    !!  gmres_maxiter  = Maximum number of iterations for the FGMRES algorithm.
    !!  gmres_rtol     = Relative tolerance for convergence.
    !!
    !! INPUT/OUTPUTS
    !!  est            = Initial guess for the solution vector, updated with the computed solution.
    !!
    !! SOURCE
    subroutine call_FGMRES(n, matvec, rhs, est, gmres_maxiter, gmres_rtol)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp),intent(in) :: rhs(:)
        real(dp),intent(inout) :: est(:)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
        !Local variables-------------------------------
        !MKL FGMRES
        integer :: RCI_request, itercount, size_vres
        integer :: ipar(128)
        real(dp) :: dpar(128)
        real(dp), allocatable :: tmp(:)

        ! *************************************************************************
        
        !FGMRES initialization

        ABI_MALLOC(tmp, ((2*gmres_maxiter+1)*n + gmres_maxiter*(gmres_maxiter+9)/2 + 1))
        call dfgmres_init(n, est, rhs, RCI_request, ipar, dpar, tmp)
        !setting FGMRES parameters
        ipar(7) = 0              ! control verbosity : no warning message
        ipar(5) = gmres_maxiter  ! maximum number of iterations
        ipar(8) = 1              ! dfgmres routine performs the stopping test for the maximum number of iterations ipar(4)≤ipar(5)
        ipar(9) = 1              ! dfgmres routine performs the residual stopping test dpar(5)≤dpar(4)=dpar(1)*dpar(3)+dpar(2)
        ipar(10) = 0             ! no user defined stopping tests
        ipar(11) = 0             ! non-preconditioned GMRES
        ipar(12) = 1             ! dfgmres routine performs the automatic test dpar(7)≤dpar(8)
        ipar(15) = gmres_maxiter ! number of the non-restarted FGMRES iterations (no restart here)
        dpar(1) = gmres_rtol     ! relative tolerance
        !dpar(2) = 0.01          ! absolute tolerance
        
        !FGMRES iterations
        
        call dfgmres_check(n, est, rhs, RCI_request, ipar, dpar, tmp)
        call dfgmres(n, est, rhs, RCI_request, ipar, dpar, tmp)
        
        do
            if (RCI_request==-1) then
            !    maximum number of iterations is reached
                call dfgmres_get(n, est, rhs, RCI_request, ipar, dpar, tmp, itercount)
                exit
            else if (RCI_request==0) then
            !    successful completion of the task
                call dfgmres_get(n, est, rhs, RCI_request, ipar, dpar, tmp, itercount)
                exit
            else  if (RCI_request==1) then
            !    multiply the matrix P by tmp(ipar(22)) and put the result in tmp(ipar(23))
                call matvec(n, tmp(ipar(22):ipar(22)+2*size_vres-1), tmp(ipar(23):ipar(23)+2*size_vres-1))
            !    proceed with FGMRES iterations
                call dfgmres(2*size_vres, est, rhs, RCI_request, ipar, dpar, tmp)
        !---------------------------------------------------------------------
        !  FGMRES Errors
            else if (RCI_request==-10) then
                ABI_BUG('FGMRES : attempt to divide by zero')
                exit
            else if (RCI_request==-11) then
                ABI_BUG('FGMRES : infinite cycle')
                exit
            else if (RCI_request==-12) then
                ABI_BUG('FGMRES : errors were found in the method parameters')
                exit
            ! RCI_request = 2, 3, 4 should not happen with this choice of parameters
            else
                ABI_BUG('FGMRES : RCI_request has unexpected value')
            end if
        !---------------------------------------------------------------------
        end do
        ABI_FREE(tmp)
    end subroutine call_FGMRES

    !****f* m_iterative_solvers/call_gmresm
    !! NAME
    !!  call_gmresm
    !!
    !! FUNCTION
    !!  Call the gmresm (Willis, A. (2017) SoftwareX 6, 124-127, code at the end of this file) routine to solve a linear system.
    !!
    !! INPUTS
    !!  n              = Size of the matrix.
    !!  matvec         = Subroutine that performs matrix-vector multiplication.
    !!  rhs            = Right-hand side vector of the linear system.
    !!  gmres_maxiter  = Maximum number of iterations for the GMRES algorithm.
    !!  gmres_rtol     = Relative tolerance for convergence.
    !!
    !! INPUT/OUTPUTS
    !!  est            = Initial guess for the solution vector, updated with the computed solution.
    !!
    !! SOURCE
    subroutine call_gmresm(n, matvec, est, rhs, gmres_maxiter, gmres_rtol)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp),intent(in) :: rhs(n)
        real(dp),intent(inout) :: est(n)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
        !Local variables-------------------------------
        integer :: its, info, m
        real(dp) :: res, del
        real(dp), allocatable :: h(:, :), v(:, :)

        ! *************************************************************************

        m = gmres_maxiter
        ABI_MALLOC(h, (m+1, m))
        ABI_MALLOC(v, (n, m+1))
        res = gmres_rtol
        del = 0
        its = gmres_maxiter  ! No restart
        info = 1
        call gmresm(m, n, est, rhs, matvec, psolve, dotprd, h, v, res, del, its, info)
        ABI_FREE(h)
        ABI_FREE(v)

        contains
        ! Dummy :  No preconditioning
        subroutine psolve(n_, x)
            integer, intent(in) :: n_
            real(dp), intent(inout) :: x
        end subroutine psolve
        ! Dot product
        function dotprd(n_, a, b) result(c)
            integer, intent(in) :: n_
            real(dp), intent(inout) :: a(n_), b(n_)
            real(dp) :: c
            ! ***********************
            c = dot_product(a, b)
        end function dotprd

    end subroutine call_gmresm

    !****f* m_iterative_solvers/gmres_linear_solver
    !! NAME
    !!  gmres_linear_solver
    !!
    !! FUNCTION
    !!  Solve a linear system using GMRES. Depending on the availability of MKL, 
    !!  it either calls the MKL FGMRES routine or the gmresm routine.
    !!
    !! INPUTS
    !!  n              = Size of the matrix.
    !!  matvec         = Subroutine that performs matrix-vector multiplication.
    !!  rhs            = Right-hand side vector of the linear system.
    !!  gmres_maxiter  = Maximum number of iterations for the GMRES algorithm.
    !!  gmres_rtol     = Relative tolerance for convergence.
    !!
    !! INPUT/OUTPUTS
    !!  est            = Initial guess for the solution vector, updated with the computed solution.
    !!
    !! SOURCE
    subroutine gmres_linear_solver(n, matvec, rhs, est, gmres_maxiter, gmres_rtol)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp), intent(in) :: rhs(n)
        real(dp), intent(inout) :: est(n)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
      
        ! *************************************************************************
        
        !TODO : dirty check of MKL availability
#if defined HAVE_LINALG_MKL_OMATCOPY
        call call_FGMRES(n, matvec, rhs, est, gmres_maxiter, gmres_rtol)
#else
        call call_gmresm(n, matvec, est, rhs, gmres_maxiter, gmres_rtol)
#endif
      
    end subroutine gmres_linear_solver

!-------------------------------------------------------------------------------------------

!----------------------------------------------------------------------
! Openpipeflow.org.  If used in your work, please cite
! Willis, A. (2017) SoftwareX 6, 124-127.
! https://doi.org/10.1016/j.softx.2017.05.003 (open access)
!                                      Thanks in advance! Ashley 2019.
!----------------------------------------------------------------------
! solve A x = b for x ;  
! minimise |Ax-b| subject to constraint |x| < delta .
! requires lapack routines dgelsy, dgesvd.
!----------------------------------------------------------------------
! m	  gmres dimension
! n 	  dimension of x
! x	  on input:  guess for x, can be 0
!         on exit:  solution x, subject to constraint if del>0
! b	  input b
! matvec  performs y := A x, call matvec(N,x, y)
! psolve  preconditioner, solve M x_out = x_in, call psolve(N,x)
! dotprd  dot product, d = dotprd(n,a,b)
! h       Hessian matrix,  size (m+1)*m
! v       Krylov subspace, size n*(m+1)
! res	  on input: |Ax-b|/|b|<res; 
!         on exit:  residual reached
! del     on input: if(del>0) then the x returned is the hookstep
!         on exit:  norm of next b predicted by hook
! its	  on input: max num its; 
!         on exit:  number of its taken
! info	  on input: if(info==1) print* residuals
!                   if(info==2) recalc hookstep with new del
! 	  on exit:  0 sucessful, 1 method breakdown, 2 max its
!							A.P.Willis 2008
!----------------------------------------------------------------------

 subroutine gmresm(m,n,x,b,matvec,psolve,dotprd,h,v,res,del,its,info)
   implicit none
   integer,          intent(in)    :: m
   integer,          intent(in)    :: n
   double precision, intent(inout) :: x(n)
   double precision, intent(in)    :: b(n)
   external                        :: matvec,psolve
   double precision, external      :: dotprd
   double precision, intent(inout) :: h(m+1,m)
   double precision, intent(inout) :: v(n,m+1)
   double precision, intent(inout) :: res
   double precision, intent(inout) :: del
   integer,          intent(inout) :: its
   integer,          intent(inout) :: info
   double precision :: tol,res_,stgn, w(n), z(n)
   double precision :: h_(m+1,m), y(m+1), p(m+1), work(4*m+1)
   integer :: imx, piv(m), rank, i
   double precision, save :: beta
   integer, save :: j
   logical :: done   

   if(info==2) then
      call hookstep(j,h,m,beta,del, y)
      z = matmul(v(:,1:j),y(1:j))
      call psolve(n, z)
      x = z
      info = 0
      return
   end if	 

   tol = res
   imx = its
   its = 0
   v   = 0d0

 1 continue
   res_ = 1d99
   stgn = 1d0 - 1d-14
 
   beta = dsqrt(dotprd(n,x,x)) 
   if(beta==0d0)  w = 0d0
   if(beta/=0d0)  call matvec(n,x, w)
   w = b - w
   beta = dsqrt(dotprd(n,w,w)) 
   v(:,1) = w / beta
     
   h = 0d0
   do j = 1, m
      its = its + 1
      z = v(:,j)      
      call psolve(n, z)
      call matvec(n, z, w)
      do i = 1, j
         h(i,j) = dotprd(n,w,v(1,i))
         w = w - h(i,j)*v(:,i)
      end do
      h(j+1,j) = dsqrt(dotprd(n,w,w))
      v(:,j+1) = w / h(j+1,j)
          
      p(1) = beta
      p(2:j+1) = 0d0
      h_(1:j+1,1:j) = h(1:j+1,1:j)
      !call dgelsy(j+1,j,1,h_(1:m+1, 1:j),m+1,p,m+1,piv,m,rank,work,4*m+1,i)
      call dgelsy(j+1,j,1,h_,m+1,p,m+1,piv,m,rank,work,4*m+1,i)
      if(i/=0) stop 'gmresm: dgelsy'
      y = p

      p(1:j+1) = - matmul(h(1:j+1,1:j),y(1:j))
      p(1) = p(1) + beta
      res = dsqrt(dot_product(p(1:j+1),p(1:j+1)))
      if(info==1) print*, 'gmresm: it=', its,' res=', real(res)
      
      done = (res<=tol .or. its==imx .or. res>res_)
      if(done .or. j==m) then
        if(del>0d0)  call hookstep(j,h,m,beta,del, y)
         z = matmul(v(:,1:j),y(1:j))
         call psolve(n, z)
         x = x + z
        if(its==imx) info = 2
        if(res>res_) info = 1
        if(res<=tol) info = 0
         if(done)     return
        if(del>0d0)  print*, 'gmres: warning! restart affects hookstep'
         goto 1       ! (j==m) restart
      end if
      res_ = res*stgn

   end do   
 
 end subroutine gmresm
 
 
!-----------------------------------------------------------------
! replace y with a vector that generates a hookstep
! c.f. Viswanath (2008) arXiv:0809.1498
!-----------------------------------------------------------------
 subroutine hookstep(j,h,m,beta,del, y)
   implicit none
   integer,          intent(in)    :: j, m
   double precision, intent(in)    :: h(m+1,j), beta
   double precision, intent(inout) :: del
   double precision, intent(out)   :: y(j)
   double precision :: a(j+1,j), s(j), u(j+1,j+1), vt(j,j), work(5*(j+1))
   double precision :: p(j+1), q(j), mu, qn
   integer :: info
   
   a = h(1:j+1,1:j)
   
   call dgesvd('A','A',j+1,j,a,j+1,s,u,j+1,vt,j,work,5*(j+1),info)
   if(info/=0) stop 'hookstep: dgesvd'
   
   p(1:j) = beta * u(1,1:j)   

   mu = max(s(j)*s(j)*1d-6,1d-99)
   qn = 1d99
   do while(qn>del)
      mu = mu * 1.1d0
      q = p(1:j)*s/(mu+s*s)
      qn = dsqrt(dot_product(q,q))
   end do

   y = matmul(q,vt)

   p = - matmul(h(1:j+1,1:j),y(1:j))
   p(1) = p(1) + beta
   del = dsqrt(dot_product(p,p))
 
 end subroutine hookstep

end module m_iterative_solvers