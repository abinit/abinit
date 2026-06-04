!!****m* ABINIT/m_iterative_solvers
!! NAME
!!  m_iterative_solvers
!!
!! FUNCTION
!!  This module contains basic (matrix-free) MPI-aware iterative solvers (GMRES, CG).
!!
!! COPYRIGHT
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_iterative_solvers
    
    use m_errors
    use defs_basis
    use m_xmpi
    use m_specialmsg

    implicit none
    private
    public :: cg_linear_solver, gmres_linear_solver, cg_eigen_solver_treshold
 
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
        real(dp) :: x(n), Ax(n), p(n), g(n), beta
        real(dp) :: residual_norm, rayleigh_quotient
        integer :: i, j, iter
        logical :: do_exit
        integer :: ierr
        
        ! *************************************************************************
        
        ! Initialize variables
        x = x0
        call matvec(n, x, Ax)                                   ! Compute initial matrix-vector product
        rayleigh_quotient = dot_product(x, Ax) / dot_product(x, x)
        g = 2/dot_product(x, x) * (Ax - rayleigh_quotient * x)  ! Gradient of Rayleigh quotient
        p = 0
        beta = 0

        ! Conjugate Gradient iterations to minimize Rayleigh quotient
        do iter = 1, max_iter
            call cg_update(n, matvec, x, Ax, p, g, beta, rayleigh_quotient, residual_norm)
            do_exit = (residual_norm < tol)
            call xmpi_bcast(do_exit, 0, xmpi_world, ierr)
            if (do_exit) exit
        end do

        ! Store the computed eigenvalue and eigenvector
        eigenvalues(1) = rayleigh_quotient
        eigenvectors(:, 1) = x/norm2(x)
        n_eig = 1

        ! If more eigenvalues are required, use deflation to compute subsequent eigenvalues
        do j = 2, max_neig

            ! Check if the eigenvalue is below the threshold
            do_exit = (eigenvalues(j-1) < eigenvalue_threshold)
            call xmpi_bcast(do_exit, 0, xmpi_world, ierr)
            if (do_exit) exit

            ! If it is not, we need to compute more eigenvalues :
            ! Orthogonalize the initial guess against previously computed eigenvectors
            x = x0
            do i = 1, j - 1
                x = x - dot_product(p, eigenvectors(:, i))/dot_product(eigenvectors(:, i), eigenvectors(:, i)) * eigenvectors(:, i)
            end do

            ! Initialization
            call matvec(n, x, Ax)                                   ! Compute initial matrix-vector product
            rayleigh_quotient = dot_product(x, Ax) / dot_product(x, x)
            g = 2/dot_product(x, x) * (Ax - rayleigh_quotient * x)  ! Gradient of Rayleigh quotient
            p = 0
            beta = 0

            ! Perform conjugate gradient iterations
            do iter = 1, max_iter
                call cg_update(n, matvec, x, Ax, p, g, beta, rayleigh_quotient, residual_norm, eigenvectors(:, 1:j-1))
                do_exit = (residual_norm < tol)
                call xmpi_bcast(do_exit, 0, xmpi_world, ierr)   ! MPI aware: avoid desynchronization.
                if (do_exit) exit
            end do

            ! Store the computed eigenvalue and eigenvector
            eigenvalues(j) = rayleigh_quotient
            eigenvectors(:, j) = x/norm2(x)
            n_eig = n_eig + 1

        end do

    end subroutine cg_eigen_solver_treshold

    !****f* m_iterative_solvers/cg_update
    !! NAME
    !!  cg_update
    !!
    !! FUNCTION
    !!  Perform a single Conjugate Gradient update step to minimize the Rayleigh quotient 
    !!  of a self-adjoint operator given by matvec. This method is used to approximate 
    !!  the eigenvalue and eigenvector of the operator.
    !!
    !! INPUTS
    !!  n              = Integer, size of the matrix (number of rows/columns).
    !!  matvec         = Subroutine, performs matrix-vector multiplication (A * v).
    !!
    !! INPUT/OUTPUTS
    !!  x              = Current solution vector (eigenvector approximation).
    !!  Ax             = A * x, where A is the operator defined by matvec.
    !!  p              = Previous search direction vector.
    !!  g              = Gradient of the Rayleigh-quotient at x.
    !!  beta           = Scalar, used to update the search direction.
    !!
    !! OUTPUTS
    !!  rayleigh_quotient = Current approximation of the eigenvalue (lambda).
    !!  residual_norm     = Norm of the residual.
    !!
    !! NOTES
    !!  - The input vectors (x, r, p, Ap) must be properly initialized before calling 
    !!    this subroutine.
    !!  - The subroutine assumes that the operator is self-adjoint (Hermitian).
    !!
    !! SOURCE
    subroutine cg_update(n, matvec, x, Ax, p, g, beta, rayleigh_quotient, residual_norm, eigenvectors)
        ! Arguments
        integer, intent(in) :: n
        real(dp), optional, intent(in) :: eigenvectors(:,:)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
        real(dp), intent(inout) :: x(n), Ax(n), g(n), p(n), beta
        real(dp), intent(inout) :: rayleigh_quotient, residual_norm

        ! Local variables
        real(dp) :: Ap(n), a, b, c, d, e, f, alpha, alpha_(2)
        integer :: n_eig, i

        ! *************************************************************************
        if (present(eigenvectors)) then
            n_eig = size(eigenvectors, 2)   ! Number of eigenvectors already computed
        else
            n_eig = 0
        end if

        ! Update search direction
        p = -g + beta*p
        ! Orthogonalize the search direction against previously computed eigenvectors
        do i = 1, n_eig
            p = p - dot_product(p, eigenvectors(:, i))/dot_product(eigenvectors(:, i), eigenvectors(:, i)) * eigenvectors(:, i)
        end do

        ! Compute A*p
        call matvec(n, p, Ap)

        ! Compute alpha (step size = minimizer of R(x+alpha*p) that is the solution (+) of a quadratic problem)
        a = dot_product(p, Ap)
        b = 2*dot_product(x, Ap)
        c = dot_product(x, Ax)
        d = dot_product(p, p)
        e = 2*dot_product(x, p)
        f = dot_product(x, x)
        alpha_ = quadratic_roots(a*e-b*d, 2*(f*a-d*c), b*f-c*e)
        alpha = alpha_(1)

        ! Update solution vector
        x = x + alpha * p
        Ax = Ax + alpha * Ap
        
        ! Update the Rayleigh-quotient
        rayleigh_quotient = dot_product(x, Ax)/dot_product(x, x)

        ! Update the (x-normalized) gradient and beta = dot(g, g)/dot(g_old, g_old)
        beta = 1/dot_product(g, g)
        g = 2/norm2(x) * (Ax - rayleigh_quotient * x)
        beta = beta * dot_product(g, g)

        ! Update the residual norm (= norm of the gradient)
        residual_norm = norm2(g)

    end subroutine cg_update

    function quadratic_roots(a, b, c) result(r)
        ! Arguments
        real(dp) :: a, b, c
        real(dp) :: r(2)

        ! Local variables
        real(dp) :: discriminant, sqrt_discriminant

        ! *************************************************************************

        if (a == 0.0_dp) then
            r(1) = 0.0_dp
            r(2) = 0.0_dp
            return
        end if

        discriminant = b**2 - 4.0_dp * a * c
        if (discriminant < 0.0_dp) then
            r(1) = 0.0_dp
            r(2) = 0.0_dp
            return
        end if

        sqrt_discriminant = sqrt(discriminant)
        r(1) = (-b + sqrt_discriminant) / (2.0_dp * a)
        r(2) = (-b - sqrt_discriminant) / (2.0_dp * a)

    end function quadratic_roots

!-------------------------------------------------------------------------------------------
! Linear solvers (CG and GMRES)
!-------------------------------------------------------------------------------------------

    subroutine cg_linear_solver(n, matvec, rhs, est, cg_maxiter, cg_rtol, verbose)
        
        !Arguments ------------------------------------
        integer, intent(in) :: n, cg_maxiter
        real(dp), intent(in) :: cg_rtol
        real(dp), intent(in) :: rhs(:)
        logical, intent(in) :: verbose
        real(dp),intent(inout) :: est(:)
        character(len=500) :: msg
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface

        !Local variables-------------------------------
        real(dp) :: r(n), p(n), Ap(n)
        real(dp) :: alpha, beta, rsold, rsnew
        integer :: iter
        logical :: do_exit
        integer :: ierr

        ! *************************************************************************

        ! Initialize
        call matvec(n, est, Ap)
        r = rhs - Ap
        p = r
        rsold = dot_product(r, r)

        ! Conjugate Gradient iterations
        do iter = 1, cg_maxiter
            call matvec(n, p, Ap)
            alpha = rsold / dot_product(p, Ap)
            est = est + alpha * p
            r = r - alpha * Ap
            rsnew = dot_product(r, r)
            
            ! Check for convergence
            do_exit = (sqrt(rsnew) < cg_rtol)
            call xmpi_bcast(do_exit, 0, xmpi_world, ierr)   ! MPI aware: avoid desynchronization.
            if (do_exit) exit

            beta = rsnew / rsold
            p = r + beta * p
            write(msg, *)'cg: it=', iter,' res=', sqrt(rsnew)
            if (verbose) call wrtout(std_out, msg)
            !if (verbose) write(std_out,*) 'cg: it=', iter,' res=', sqrt(rsnew)
            rsold = rsnew
        end do

    end subroutine cg_linear_solver

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
    !!  verbose        = Logical, if true, print residuals at each iteration.
    !!
    !! INPUT/OUTPUTS
    !!  est            = Initial guess for the solution vector, updated with the computed solution.
    !!
    !! SOURCE
    subroutine call_gmresm(n, matvec, est, rhs, gmres_maxiter, gmres_rtol, verbose)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp), intent(in) :: rhs(n)
        logical, intent(in) :: verbose
        real(dp), intent(inout) :: est(n)
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
        res = gmres_rtol * norm2(rhs)
        del = 0
        its = gmres_maxiter  ! No restart
        info = 0
        if (verbose) then
            info = 1
        end if
        call gmresm(m, n, est, rhs, matvec, psolve, dotprd, h, v, res, del, its, info)
        ABI_FREE(h)
        ABI_FREE(v)

        contains
        ! Dummy :  No preconditioning
        subroutine psolve(n_, x)
            integer, intent(in) :: n_
            real(dp), intent(inout) :: x(n_)
            ! ***********************
            ! We do nothing here but don't wan't to be flashed by abirule.
            if (.false.) then
                x = zero
            end if
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
    subroutine gmres_linear_solver(n, matvec, rhs, est, gmres_maxiter, gmres_rtol, verbose)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp), intent(in) :: rhs(n)
        logical :: verbose
        real(dp), intent(inout) :: est(n)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
      
        ! *************************************************************************
        
        call call_gmresm(n, matvec, est, rhs, gmres_maxiter, gmres_rtol, verbose)
      
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
   integer :: ierr  
   character(len=500) :: msg


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
      ! MPI aware: broadcast the 'res' value of master to avoid desynchronization.
      call xmpi_bcast(res, 0, xmpi_world, ierr)
      !if(info==1) print*, 'gmresm: it=', its,' res=', real(res)
      write(msg, *)'gmresm: it=', its,' res=', real(res)
      if(info==1) call wrtout(std_out, msg)

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
        !if(del>0d0)  print*, 'gmres: warning! restart affects hookstep'
        if(del>0d0) call wrtout(std_out, 'gmres: warning! restart affects hookstep')
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