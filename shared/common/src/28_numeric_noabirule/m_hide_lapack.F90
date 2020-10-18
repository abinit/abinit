!!****m* ABINIT/m_hide_lapack
!! NAME
!!  m_hide_lapack
!!
!! FUNCTION
!!  ABINIT Linear Algebra Subroutine Interfaces.
!!
!!  This modules provides interfaces performing the overloading of commonly used Lapack routines.
!!  The main purpose of this module is to create a layer between abinit routines and Lapack procedures.
!!  This layer can be used to hide the parallel Scalapack version. In this case, only the MPI commutator
!!  has to be provided in input as the wrapper will take care of the initialization of the Scalapack grid as
!!  well as of the distribution of the matrix. Note that this allows one to reduce
!!  the CPU time per processor but not the memory allocated since the entire matrix has to be provided in input.
!!  The interfaces are very similar to the Lapack F77 version (neither F90 constructs nor
!!  F90 assumed size arrays are used). The main simplification with respect to the F77 version
!!  of Lapack is that the work arrays are allocated inside the wrapper with optimal size
!!  thus reducing the number of input argcomm_scalapackuments that has to be passed.
!!  Leading dimensions have been removed from the interface whenever possible.
!!  In F90 one can pass the array descriptor if the routines should operate on a slice
!!  of the local array (seldom done in abinit).
!!  Using array descriptor is OK but will it likely slow-down the calculation as
!!  some compilers perform a copy of the input-output data.
!!  If efficiency is a concern, then the F77 call should be used
!!
!! COPYRIGHT
!! Copyright (C) 1992-2020 ABINIT group (MG, GMR, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! TODO
!!  1) Use a function to define the size of the Scalapack block according to some heuristic method.
!!  2) Define a threshold below which Scalapack is not used although the MPI communicator is passed.
!!  3) Split MPI communicator for Scalapack (.i.e. use a max size for the Scalapack comm; see abi_linalg_init).
!!  4) On certain networks, xmpi_sum might crash due to the size of the MPI packet.
!!     This problem should be solved in hide_mpi (Module containing a private global variable
!!     defining a threshold above which the input array is split into smaller chunks.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_hide_lapack

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_slk
 use m_linalg_interfaces

 use m_time,       only : cwtime
 use m_fstrings,   only : firstchar

 implicit none

 private

! ==========================
! Complex Hermitian matrices
! ==========================
! The _cplex version receives matrices declared as arr(cplex,N,M)
! and calls the complex/real version depending on cplex (used e.g. for istwf_k = 2)

 public :: xheev   ! Computes all the eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix.
 public :: xheev_cplex

 public :: xhpev   ! Computes all the eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix
                   ! in packed storage (Scalapack version not available)

 public :: xhegv   ! Compute all the eigenvalues, and optionally, the eigenvectors of a complex generalized
                   ! Hermitian-definite eigenproblem, of the form:
                   ! A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x

 public :: xheevx  ! Computes selected eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.
                   ! Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
                   ! indices for the desired eigenvalues.

 public :: xheevx_cplex

 public :: xhegvx  ! Computes selected eigenvalues, and optionally, eigenvectors
                   ! of a complex generalized Hermitian-definite eigenproblem,
                   ! of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x.
                   ! Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
                   ! indices for the desired eigenvalues.
 public :: xhegvx_cplex


 public :: xhesv_cplex  ! Solve A * X = B, where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS matrices.


! ==============================
! Complex non-symmetric matrices
! ==============================

 public :: xgeev   ! Computes for a complex nonsymmetric matrix A, the eigenvalues and, optionally,
                   ! the left and/or right eigenvectors.

 public :: xginv   ! Invert a general matrix of complex elements by means of LU factorization.


 public :: xhdp_invert   ! Invert a Hermitian positive definite matrix.

 interface xheev
   module procedure wrap_CHEEV
   module procedure wrap_ZHEEV
 end interface xheev

 interface xhpev
   module procedure wrap_CHPEV
   module procedure wrap_ZHPEV
 end interface xhpev

 interface xhegv
   module procedure wrap_ZHEGV
 end interface xhegv

 public :: xhegv_cplex

 interface xheevx
   module procedure wrap_ZHEEVX
 end interface xheevx

 interface xhegvx
   module procedure wrap_ZHEGVX
 end interface xhegvx

 interface xgeev
   module procedure wrap_CGEEV
   module procedure wrap_ZGEEV
 end interface xgeev

 interface xginv
   module procedure cginv
   module procedure zginv
 end interface xginv

 interface xhdp_invert
   module procedure zhpd_invert
 end interface xhdp_invert

 public :: matrginv      ! Invert a general matrix of real*8 elements.
 public :: matr3eigval   ! Find the eigenvalues of a real symmetric 3x3 matrix, entered in full storage mode.

 !FIXME This procedures are deprecated, use lapack API
 public :: jacobi        ! Computes all eigenvalues and eigenvectors of a real symmetric matrix a,
 public :: ludcmp
 public :: lubksb
 public :: dzgedi
 public :: dzgefa

!----------------------------------------------------------------------
! support for unitary tests and profiling.

 type,public :: latime_t
   character(len=500) :: testname
   integer :: msize
   real(dp) :: ctime
   real(dp) :: wtime
   real(dp) :: max_abserr=-one
   real(dp) :: gflops
 end type latime_t

 public :: test_xginv

!----------------------------------------------------------------------
! private variables

 integer,private,parameter :: SLK_BLOCK_SIZE = 24
 ! Default block size for Scalapack distribution.
 ! As recommended by Intel MKL, a more sensible default than the previous value of 40

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_hide_lapack/wrap_CHEEV
!! NAME
!!  wrap_CHEEV
!!
!! FUNCTION
!!  wrap_CHEEV computes the eigenvalues and, optionally, the eigenvectors of a
!!  complex Hermitian matrix in single precision. [PRIVATE]
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!! OUTPUT
!!  W       (output) REAL(SP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) COMPLEX(SPC) array, dimension (N, N)
!!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = 'L',
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!!          orthonormal eigenvectors of the matrix A.
!!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!!          or the upper triangle (if UPLO='U') of A, including the
!!          diagonal, is destroyed.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_CHEEV(jobz, uplo, n, a, w)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 character(len=*),intent(in) :: jobz,uplo
!scalars
 real(sp),intent(out) :: w(n)
 complex(spc),intent(inout) :: a(n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info
 character(len=500) :: msg
!arrays
 real(sp),allocatable :: rwork(:)
 complex(spc),allocatable :: work(:)

!************************************************************************

 lwork = MAX(1,2*n-1)

 ABI_MALLOC(work, (lwork))
 ABI_MALLOC(rwork, (MAX(1,3*n-2)))

 call CHEEV(jobz,uplo,n,a,n,w,work,lwork,rwork,info)

 if (info < 0) then
   write(msg,'(a,i0,a)')"The ",-info,"-th argument of CHEEV had an illegal value."
   MSG_ERROR(msg)
 end if

 if (info > 0) then
   write(msg,'(2a,i0,a)')&
    "CHEEV: the algorithm failed to converge; ",ch10,&
    info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
   MSG_ERROR(msg)
 end if

 ABI_FREE(rwork)
 ABI_FREE(work)

 !TODO scaLAPACK version (complex single precision buffer is needed in matrix_scalapack)

end subroutine wrap_CHEEV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_ZHEEV
!! NAME
!!  wrap_ZHEEV
!!
!! FUNCTION
!!  wrap_ZHEEV computes the eigenvalues and, optionally, the eigenvectors of a
!!  complex Hermitian matrix in double precision. [PRIVATE]
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!! OUTPUT
!!  W       (output) REAL(DP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) COMPLEX(DPC) array, dimension (N, N)
!!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = 'L',
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!!          orthonormal eigenvectors of the matrix A.
!!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!!          or the upper triangle (if UPLO='U') of A, including the
!!          diagonal, is destroyed.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_ZHEEV(jobz, uplo, n, a, w, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,optional,intent(in) :: comm
 character(len=*),intent(in) :: jobz,uplo
!arrays
 complex(dpc),intent(inout) :: a(n,n)
 real(dp),intent(out) :: w(n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_mat,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif
!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)

   lwork = MAX(1,2*n-1)
   ABI_MALLOC(work, (lwork))
   ABI_MALLOC(rwork, (MAX(1,3*n-2)))

   call ZHEEV(jobz,uplo,n,a,n,w,work,lwork,rwork,info)

   if (info < 0) then
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEEV had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(2a,i0,a)')&
     "ZHEEV: the algorithm failed to converge; ",ch10,&
     info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    MSG_ERROR(msg)
   end if

   ABI_FREE(rwork)
   ABI_FREE(work)
   RETURN

 CASE (.TRUE.)
#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE
   call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call slk_matrix_from_global_dpc_2D(Slk_mat,uplo,a)

   want_eigenvectors = firstchar(jobz,(/"V","v"/))
   if (want_eigenvectors) then ! Initialize the distributed vectors.
    call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   end if

   ! Solve the problem with scaLAPACK.
   call slk_pzheev(jobz,uplo,Slk_mat,Slk_vec,w)

   call destruction_matrix_scalapack(Slk_mat)

   if (want_eigenvectors) then ! A is overwritten with the eigenvectors
    a = czero
    call slk_matrix_to_global_dpc_2D(Slk_vec,"All",a) ! Fill the entries calculated by this node.
    call destruction_matrix_scalapack(Slk_vec)
    call xmpi_sum(a,comm,ierr)                        ! Fill the remaing entries of the global matrix
   end if

   call end_scalapack(Slk_processor)
   RETURN
#endif

  MSG_BUG("You should not be here!")
 END SELECT

end subroutine wrap_ZHEEV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/xheev_cplex
!! NAME
!!  xheev_cplex
!!
!! FUNCTION
!!  xheev_cplex computes the eigenvalues and, optionally, the eigenvectors of a
!!  (complex Hermitian| real symmetric) matrix in double precision.
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  CPLEX= Size of the first dimension of the A matrix.
!!          1 for a real symmetric matrix.
!!          2 for complex Hermitian matrix stored in a real array with real and imaginary part.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  W       (output) REAL(DP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) REAL(DP) array, dimension (CPLEX, N, N)
!!          On entry, the (complex Hermitian|Real symmetric) matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = 'L',
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!!          orthonormal eigenvectors of the matrix A.
!!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!!          or the upper triangle (if UPLO='U') of A, including the
!!          diagonal, is destroyed.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xheev_cplex(jobz, uplo, cplex, n, a, w, msg, ierr, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,cplex
 integer,optional,intent(in) :: comm
 character(len=*),intent(in) :: jobz,uplo
 integer,intent(out) :: ierr
 character(len=*),intent(out) :: msg
!arrays
 real(dp),intent(inout) :: a(cplex,n,n)
 real(dp),intent(out) :: w(n)

!Local variables ------------------------------
!scalars
 integer :: lwork,nprocs
 logical :: use_scalapack
!arrays
 real(dp),allocatable :: rwork(:)
 real(dp),allocatable :: work_real(:)
 complex(dpc),allocatable :: work_cplx(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_mat,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif
!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 if (ALL(cplex/= [1, 2])) then
   write(msg,'(a,i0)')" Wrong value for cplex: ",cplex
   ierr = 1; return
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
  if (cplex==1) then
    ! Real symmetric case.
    lwork = MAX(1,3*n-1)
    ABI_MALLOC(work_real,(lwork))

    call DSYEV(jobz,uplo,n,a,n,w,work_real,lwork,ierr)

    if (ierr < 0) then
     write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of DSYEV had an illegal value."
    end if

    if (ierr > 0) then
     write(msg,'(2a,i0,a)')&
       "DSYEV: the algorithm failed to converge; ",ch10,&
       ierr," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    end if

    ABI_FREE(work_real)
    RETURN

  else
    ! Hermitian case.
    lwork = MAX(1,2*n-1)

    ABI_MALLOC(work_cplx, (lwork))
    ABI_MALLOC(rwork, (MAX(1,3*n-2)))

    call ZHEEV(jobz,uplo,n,a,n,w,work_cplx,lwork,rwork,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of ZHEEV had an illegal value."
    end if

    if (ierr > 0) then
      write(msg,'(2a,i0,a)')&
       "ZHEEV: the algorithm failed to converge; ",ch10,&
       ierr," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    end if

    ABI_FREE(rwork)
    ABI_FREE(work_cplx)
    RETURN
  end if ! cplex

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   MSG_ERROR("Not coded yet")

   !call init_scalapack(Slk_processor,comm)
   !istwf_k=1
   !
   !! Initialize and fill Scalapack matrix from the global one.
   !tbloc=SLK_BLOCK_SIZE
   !call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   !
   !write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   !call wrtout(std_out,msg,"PERS")
   !
   !call slk_matrix_from_global_dpc_2D(Slk_mat,uplo,a)
   !
   !want_eigenvectors = firstchar(jobz,(/"V","v"/))
   !if (want_eigenvectors) then ! Initialize the distributed vectors.
   ! call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   !end if
   !
   !! Solve the problem with scaLAPACK.
   !call slk_pzheev(jobz,uplo,Slk_mat,Slk_vec,w)
   !
   !call destruction_matrix_scalapack(Slk_mat)
   !
   !if (want_eigenvectors) then ! A is overwritten with the eigenvectors
   ! a = czero
   ! call slk_matrix_to_global_dpc_2D(Slk_vec,"All",a) ! Fill the entries calculated by this node.
   ! call destruction_matrix_scalapack(Slk_vec)
   ! call xmpi_sum(a,comm,ierr)                        ! Fill the remaing entries of the global matrix
   !end if
   !
   !call end_scalapack(Slk_processor)

   RETURN
#endif

   MSG_BUG("You should not be here!")
 END SELECT

end subroutine xheev_cplex
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_CHPEV
!! NAME
!!  wrap_CHPEV
!!
!! FUNCTION
!!  wrap_CHPEV computes all the eigenvalues and, optionally, eigenvectors of a
!!  complex Hermitian matrix in packed storage. Scalapack version is not available. [PRIVATE].
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!!  LDZ     (input) INTEGER
!!          The leading dimension of the array Z.  LDZ >= 1, and if
!!          JOBZ = 'V', LDZ >= max(1,N).
!!
!! OUTPUT
!!  W       (output) REAL(SP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!!  Z       (output) COMPLEX(SPC) array, dimension (LDZ, N)
!!          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!!          eigenvectors of the matrix A, with the i-th column of Z
!!          holding the eigenvector associated with W(i).
!!          If JOBZ = 'N', then Z is not referenced.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!
!!  AP      (input/output) COMPLEX(SPC) array, dimension (N*(N+1)/2)
!!          On entry, the upper or lower triangle of the Hermitian matrix
!!          A, packed columnwise in a linear array.  The j-th column of A
!!          is stored in the array AP as follows:
!!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!!
!!          On exit, AP is overwritten by values generated during the
!!          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!!          and first superdiagonal of the tridiagonal matrix T overwrite
!!          the corresponding elements of A, and if UPLO = 'L', the
!!          diagonal and first subdiagonal of T overwrite the
!!          corresponding elements of A.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_CHPEV(jobz, uplo, n, ap, w, z, ldz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,ldz
 character(len=*),intent(in) :: jobz,uplo
!arrays
 real(sp),intent(out) :: w(n)
 complex(spc),intent(inout) :: ap(n*(n+1)/2)
 complex(spc),intent(out) :: z(ldz,n)

!Local variables ------------------------------
!scalars
 integer :: info
 character(len=500) :: msg
!arrays
 real(sp),allocatable :: rwork(:)
 complex(spc),allocatable :: work(:)

!************************************************************************

 ABI_MALLOC(work, (MAX(1,2*n-1)))
 ABI_MALLOC(rwork, (MAX(1,3*n-2)))

 call CHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, INFO )

 if (info < 0) then
   write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEEV had an illegal value."
   MSG_ERROR(msg)
 end if

 if (info > 0) then
   write(msg,'(2a,i0,a)')&
    "ZHPEV: the algorithm failed to converge; ",ch10,&
    info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
   MSG_ERROR(msg)
 end if

 ABI_FREE(rwork)
 ABI_FREE(work)

end subroutine wrap_CHPEV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_ZHPEV
!! NAME
!!  wrap_ZHPEV
!!
!! FUNCTION
!!  wrap_ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a
!!  complex Hermitian matrix in packed storage. Scalapack version is not available. [PRIVATE].
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!!  LDZ     (input) INTEGER
!!          The leading dimension of the array Z.  LDZ >= 1, and if
!!          JOBZ = 'V', LDZ >= max(1,N).
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called. Note that scalapack does not provide native
!!        support for packed symmetric matrices. Threfore we have to distribute the full matrix among the nodes.
!!        in order to perform the calculation in parallel.
!!
!! OUTPUT
!!  W       (output) REAL(DP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!!  Z       (output) COMPLEX(DPC) array, dimension (LDZ, N)
!!          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!!          eigenvectors of the matrix A, with the i-th column of Z
!!          holding the eigenvector associated with W(i).
!!          If JOBZ = 'N', then Z is not referenced.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!
!!  AP      (input/output) COMPLEX(DPC) array, dimension (N*(N+1)/2)
!!          On entry, the upper or lower triangle of the Hermitian matrix
!!          A, packed columnwise in a linear array.  The j-th column of A
!!          is stored in the array AP as follows:
!!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!!
!!          On exit, AP is overwritten by values generated during the
!!          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!!          and first superdiagonal of the tridiagonal matrix T overwrite
!!          the corresponding elements of A, and if UPLO = 'L', the
!!          diagonal and first subdiagonal of T overwrite the
!!          corresponding elements of A. Unchanged if ScaLAPACK is used.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_ZHPEV(jobz, uplo, n, ap, w, z, ldz, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,ldz
 integer,optional,intent(in) :: comm
 character(len=*),intent(in) :: jobz,uplo
!arrays
 real(dp),intent(out) :: w(n)
 complex(dpc),intent(inout) :: ap(n*(n+1)/2)
 complex(dpc),intent(out) :: z(ldz,n)

!Local variables ------------------------------
!scalars
 integer :: info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_mat,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
   ABI_MALLOC(work, (MAX(1,2*n-1)))
   ABI_MALLOC(rwork, (MAX(1,3*n-2)))

   call ZHPEV(jobz,uplo,n,ap,w,z,ldz,work,rwork,info)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHPEV had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
     write(msg,'(2a,i0,a)')&
      "ZHPEV: the algorithm failed to converge; ",ch10,&
      info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
     MSG_ERROR(msg)
   end if

   ABI_FREE(rwork)
   ABI_FREE(work)
   RETURN

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE
   call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call slk_matrix_from_global_dpc_1Dp(Slk_mat,uplo,ap)

   want_eigenvectors = firstchar(jobz,(/"V","v"/))
   if (want_eigenvectors) then ! Initialize the distributed vectors.
    call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   end if

   ! Solve the problem with scaLAPACK.
   call slk_pzheev(jobz,uplo,Slk_mat,Slk_vec,w)

   call destruction_matrix_scalapack(Slk_mat)

   if (want_eigenvectors) then ! Collect the eigenvectors.
    z = zero
    call slk_matrix_to_global_dpc_2D(Slk_vec,"All",z) ! Fill the entries calculated by this node.
    call destruction_matrix_scalapack(Slk_vec)
    call xmpi_sum(z,comm,ierr)                        ! Fill the remaing entries of the global matrix
   end if

   call end_scalapack(Slk_processor)

   RETURN
#endif

   MSG_BUG("You should not be here!")
 END SELECT

end subroutine wrap_ZHPEV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_ZHEGV
!! NAME
!!  wrap_ZHEGV
!!
!! FUNCTION
!!  wrap_ZHEGV computes all the  eigenvalues, and  optionally, the eigenvectors of a  complex generalized
!!  Hermitian-definite eigenproblem, of  the form
!!        A*x=(lambda)*B*x  (1),
!!       A*Bx=(lambda)*x,   (2), or
!!      B*A*x=(lambda)*x    (3).
!!  Here A and B are assumed to be Hermitian and B is also positive definite.
!!
!! INPUTS
!!
!!  ITYPE   (input) INTEGER Specifies the problem type to be solved:
!!          = 1:  A*x = (lambda)*B*x
!!          = 2:  A*B*x = (lambda)*x
!!          = 3:  B*A*x = (lambda)*x
!!
!!  JOBZ    (input) CHARACTER*1
!!          = "N":  Compute eigenvalues only;
!!          = "V":  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = "U":  Upper triangle of A and B are stored;
!!          = "L":  Lower triangle of A and B are stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrices A and B.  N >= 0.
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  W       (output) REAL(DP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) COMPLEX(DPC) array, dimension (N, N)
!!          On  entry, the Hermitian matrix A.  If UPLO = "U", the leading N-by-N upper triangular part of A
!!          <S-F1>contains the upper triangular part of the matrix A.
!!          If UPLO = "L", the leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A.
!!
!!          On exit, if JOBZ = "V", then A contains the matrix Z of eigenvectors.
!!          The eigenvectors are normalized as follows: if ITYPE = 1  or  2,
!!          Z**H*B*Z  = I; if ITYPE = 3, Z**H*inv(B)*Z = I.
!!          If JOBZ = "N", then on exit the upper triangle (if UPLO="U") or the lower triangle
!!          (if UPLO="L") of A, including the diagonal, is destroyed.
!!
!!
!!  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
!!          On entry, the Hermitian positive definite matrix B.
!!          If UPLO = "U", the leading N-by-N upper triangular part of B contains the upper triangular part of the matrix B.
!!          If UPLO = "L", the leading N-by-N lower triangular part of B contains the lower triangular part of the matrix B.
!!
!!          On exit, if INFO <= N, the part of B containing the matrix is overwritten by the triangular
!!          factor U or L from the Cholesky factorization B = U**H*U or B = L*L**H.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_ZHEGV(itype, jobz, uplo, n, a, b, w, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,itype
 integer,optional,intent(in) :: comm
 character(len=*),intent(in) :: jobz,uplo
!arrays
 complex(dpc),intent(inout) :: a(n,n),b(n,n)
 real(dp),intent(out) :: w(n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs,ii
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 type(matrix_scalapack)    :: Slk_matA,Slk_matB
 type(processor_scalapack) :: Slk_processor
#endif
!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
   lwork = MAX(1,2*n-1)

   ABI_MALLOC(work, (lwork))
   ABI_MALLOC(rwork, (MAX(1,3*n-2)))

   call ZHEGV(itype,jobz,uplo,n,a,n,b,n,w,work,lwork,rwork,info)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEGV had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
     if (info<= n) then
       write(msg,'(2a,i0,a)')&
        "ZHEGV failed to converge: ",ch10,&
        info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
     else
       ii = info -n
       write(msg,'(3a,i0,3a)')&
        "ZHEGV failed to converge: ",ch10,&
        "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
        "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
     end if
     MSG_ERROR(msg)
   end if

   ABI_FREE(rwork)
   ABI_FREE(work)
   RETURN

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call init_matrix_scalapack(Slk_matA,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   call slk_matrix_from_global_dpc_2D(Slk_matA,uplo,a)

   call init_matrix_scalapack(Slk_matB,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   call slk_matrix_from_global_dpc_2D(Slk_matB,uplo,b)

   ! Solve the problem with scaLAPACK.
   MSG_ERROR("slk_pZHEGV not yet coded")
   ! TODO
   !% call slk_pzhegv(itype,jobz,uplo,Slk_matA,Slk_matB,w)

   call destruction_matrix_scalapack(Slk_matB)

   if (firstchar(jobz,(/"V","v"/))) then ! A is overwritten with the eigenvectors
     a = czero
     call slk_matrix_to_global_dpc_2D(Slk_matA,"All",a) ! Fill the entries calculated by this node.
     call xmpi_sum(a,comm,ierr)                         ! Fill the remaing entries of the global matrix
   end if

   call destruction_matrix_scalapack(Slk_matA)
   call end_scalapack(Slk_processor)
   RETURN
#endif

   MSG_BUG("You should not be here!")
 END SELECT

end subroutine wrap_ZHEGV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/xhegv_cplex
!! NAME
!!  xhegv_cplex
!!
!! FUNCTION
!!  xhegv_cplex computes all the  eigenvalues, and  optionally, the eigenvectors of a
!!  (real generalized symmetric-definite| complex generalized  Hermitian-definite)
!!  eigenproblem, of  the form
!!        A*x=(lambda)*B*x  (1),
!!       A*Bx=(lambda)*x,   (2), or
!!      B*A*x=(lambda)*x    (3).
!!  Here A and B are assumed to be (symmetric|Hermitian) and B is also positive definite.
!!
!! INPUTS
!!
!!  ITYPE   (input) INTEGER Specifies the problem type to be solved:
!!          = 1:  A*x = (lambda)*B*x
!!          = 2:  A*B*x = (lambda)*x
!!          = 3:  B*A*x = (lambda)*x
!!
!!  JOBZ    (input) CHARACTER*1
!!          = "N":  Compute eigenvalues only;
!!          = "V":  Compute eigenvalues and eigenvectors.
!!
!!  UPLO    (input) CHARACTER*1
!!          = "U":  Upper triangle of A and B are stored;
!!          = "L":  Lower triangle of A and B are stored.
!!
!!  CPLEX   Size of the first dimension of the A and B matrices.
!!          1 for a real symmetric matrix.
!!          2 for a complex Hermitian matrix.
!!
!!  N       (input) INTEGER
!!          The order of the matrices A and B.  N >= 0.
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  W       (output) REAL(DP) array, dimension (N)
!!          If INFO = 0, the eigenvalues in ascending order.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) REAL(DP) array, dimension (CPLEX,N, N)
!!          On  entry, the (real symmetric|Hermitian) matrix A.  If UPLO = "U", the leading N-by-N upper triangular part of A
!!          <S-F1>contains the upper triangular part of the matrix A.
!!          If UPLO = "L", the leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A.
!!
!!          On exit, if JOBZ = "V", then A contains the matrix Z of eigenvectors.
!!          The eigenvectors are normalized as follows:
!!          if ITYPE =1 or 2:
!!            Z**T*B*Z = I if CPLEX=1
!!            Z**H*B*Z = I if CPLEX=2.
!!          if ITYPE = 3,
!!             Z**T*inv(B)*Z = I if CPLEX=1
!!             Z**H*inv(B)*Z = I if CPLEX=2
!!
!!          If JOBZ = "N", then on exit the upper triangle (if UPLO="U") or the lower triangle
!!          (if UPLO="L") of A, including the diagonal, is destroyed.
!!
!!
!!  B       (input/output) REAL(DP) array, dimension (CPLEX,N, N)
!!          On entry, the (real symmetric|Hermitian) positive definite matrix B.
!!          If UPLO = "U", the leading N-by-N upper triangular part of B contains the upper triangular part of the matrix B.
!!          If UPLO = "L", the leading N-by-N lower triangular part of B contains the lower triangular part of the matrix B.
!!
!!          On exit, if INFO <= N, the part of B containing the matrix is overwritten by the triangular
!!          factor U or L from the Cholesky factorization
!!          B = U**T*U or B = L*L**T if CPLEX=1
!!          B = U**H*U or B = L*L**H if CPLEX=2
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xhegv_cplex(itype, jobz, uplo, cplex, n, a, b, w, msg, ierr, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,itype,cplex
 character(len=*),intent(in) :: jobz, uplo
 character(len=*),intent(out) :: msg
 integer,intent(out) :: ierr
 integer,optional,intent(in) :: comm
!arrays
 real(dp),intent(inout) :: a(cplex,n,n), b(cplex,n,n)
 real(dp),intent(out) :: w(n)

!Local variables ------------------------------
!scalars
 integer :: lwork,nprocs,ii
 logical :: use_scalapack
!arrays
 real(dp),allocatable :: rwork(:), work_real(:)
 complex(dpc),allocatable :: work_cplx(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: istwf_k, tbloc
 type(matrix_scalapack)    :: Slk_matA,Slk_matB
 type(processor_scalapack) :: Slk_processor
#endif
!************************************************************************

 use_scalapack = .FALSE.
 if (present(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = nprocs > 1
#endif
 end if

 if (all(cplex /= [1, 2])) then
   write(msg,'(a,i0)')"Wrong value for cplex: ",cplex
   ierr = 1; return
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)

  if (cplex==1) then
    ! Real symmetric case.
    lwork = MAX(1,3*n-1)

    ABI_MALLOC(work_real,(lwork))
    call DSYGV(itype,jobz,uplo,n,a,n,b,n,w,work_real,lwork,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of DSYGV had an illegal value."
    end if

    if (ierr > 0) then
      if (ierr <= n) then
        write(msg,'(2a,i0,a)')&
          " DSYGV failed to converge: ",ch10,&
        ierr," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
      else
        ii = ierr - n
        write(msg,'(3a,i0,3a)')&
         "DSYGV failed to converge: ",ch10,&
         "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
         "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
      end if
    end if

    ABI_FREE(work_real)
    return

  else
    ! complex Hermitian case
    lwork = MAX(1,2*n-1)

    ABI_MALLOC(work_cplx,(lwork))
    ABI_MALLOC(rwork,(MAX(1,3*n-2)))

    call ZHEGV(itype,jobz,uplo,n,a,n,b,n,w,work_cplx,lwork,rwork,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of ZHEGV had an illegal value."
    end if

    if (ierr > 0) then
      if (ierr <= n) then
        write(msg,'(2a,i0,a)')&
         "ZHEGV failed to converge: ",ch10,&
         ierr," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
      else
        ii = ierr -n
        write(msg,'(3a,i0,3a)')&
        "ZHEGV failed to converge: ",ch10,&
        "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
        "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
      end if
    end if

    ABI_FREE(rwork)
    ABI_FREE(work_cplx)
    return
  end if ! cplex

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK

  MSG_ERROR("Not coded yet")
  ! call init_scalapack(Slk_processor,comm)
  ! istwf_k=1

  ! ! Initialize and fill Scalapack matrix from the global one.
  ! tbloc=SLK_BLOCK_SIZE

  ! write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
  ! call wrtout(std_out,msg,"COLL")

  ! call init_matrix_scalapack(Slk_matA,n,n,Slk_processor,istwf_k,tbloc=tbloc)
  ! call slk_matrix_from_global_dpc_2D(Slk_matA,uplo,a)

  ! call init_matrix_scalapack(Slk_matB,n,n,Slk_processor,istwf_k,tbloc=tbloc)
  ! call slk_matrix_from_global_dpc_2D(Slk_matB,uplo,b)

  ! ! Solve the problem with scaLAPACK.
  ! MSG_ERROR("slk_pZHEGV not yet coded")
  ! ! TODO
  ! call slk_pzhegv(itype,jobz,uplo,Slk_matA,Slk_matB,w)

  ! call destruction_matrix_scalapack(Slk_matB)
  !
  ! if (firstchar(jobz,(/"V","v"/))) then ! A is overwritten with the eigenvectors
  !  a = czero
  !  call slk_matrix_to_global_dpc_2D(Slk_matA,"All",a) ! Fill the entries calculated by this node.
  !  call xmpi_sum(a,comm,ierr)                         ! Fill the remaing entries of the global matrix
  ! end if

  ! call destruction_matrix_scalapack(Slk_matA)

  ! call end_scalapack(Slk_processor)

  RETURN
#endif

  MSG_BUG("You should not be here!")
 END SELECT

end subroutine xhegv_cplex
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_ZHEEVX
!! NAME
!!  wrap_ZHEEVX
!!
!! FUNCTION
!!  wrap_ZHEEVX computes selected eigenvalues and, optionally, eigenvectors
!!  of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can
!!  be selected by specifying either a range of values or a range of
!!  indices for the desired eigenvalues.
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  RANGE   (input) CHARACTER*1
!!          = 'A': all eigenvalues will be found.
!!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!!                 will be found.
!!          = 'I': the IL-th through IU-th eigenvalues will be found.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!!  LDA     (input) INTEGER
!!          The leading dimension of the array A.  LDA >= max(1,N).
!!
!!  VL      (input) REAL(DP)
!!  VU      (input) REAL(DP)
!!          If RANGE='V', the lower and upper bounds of the interval to
!!          be searched for eigenvalues. VL < VU.
!!          Not referenced if RANGE = 'A' or 'I'.
!!
!!  IL      (input) INTEGER
!!  IU      (input) INTEGER
!!          If RANGE='I', the indices (in ascending order) of the
!!          smallest and largest eigenvalues to be returned.
!!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!!          Not referenced if RANGE = 'A' or 'V'.
!!
!!  ABSTOL  (input) REAL(DP)
!!          The absolute error tolerance for the eigenvalues.
!!          An approximate eigenvalue is accepted as converged
!!          when it is determined to lie in an interval [a,b]
!!          of width less than or equal to
!!
!!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!!
!!          where EPS is the machine precision.  If ABSTOL is less than
!!          or equal to zero, then  EPS*|T|  will be used in its place,
!!          where |T| is the 1-norm of the tridiagonal matrix obtained
!!          by reducing A to tridiagonal form.
!!
!!          Eigenvalues will be computed most accurately when ABSTOL is
!!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!!          If this routine returns with INFO>0, indicating that some
!!          eigenvectors did not converge, try setting ABSTOL to
!!          2*DLAMCH('S').
!!
!!          See "Computing Small Singular Values of Bidiagonal Matrices
!!          with Guaranteed High Relative Accuracy," by Demmel and
!!          Kahan, LAPACK Working Note #3.
!!
!!  LDZ     (input) INTEGER
!!          The leading dimension of the array Z.  LDZ >= 1, and if
!!          JOBZ = 'V', LDZ >= max(1,N).
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  M       (output) INTEGER
!!          The total number of eigenvalues found.  0 <= M <= N.
!!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!!
!!  W       (output) REAL(DP) array, dimension (N)
!!          On normal exit, the first M elements contain the selected
!!          eigenvalues in ascending order.
!!
!!  Z       (output) COMPLEX(DPC) array, dimension (LDZ, max(1,M))
!!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!!          contain the orthonormal eigenvectors of the matrix A
!!          corresponding to the selected eigenvalues, with the i-th
!!          column of Z holding the eigenvector associated with W(i).
!!          If an eigenvector fails to converge, then that column of Z
!!          contains the latest approximation to the eigenvector, and the
!!          index of the eigenvector is returned in IFAIL.
!!          If JOBZ = 'N', then Z is not referenced.
!!          Note: the user must ensure that at least max(1,M) columns are
!!          supplied in the array Z; if RANGE = 'V', the exact value of M
!!          is not known in advance and an upper bound must be used.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) COMPLEX(DPC) array, dimension (N, N)
!!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = 'L',
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!          On exit, the lower triangle (if UPLO='L') or the upper
!!          triangle (if UPLO='U') of A, including the diagonal, is
!!          destroyed.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_ZHEEVX(jobz,range,uplo,n,a,vl,vu,il,iu,abstol,m,w,z,ldz,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,iu,ldz,n
 integer,optional,intent(in) :: comm
 integer,intent(inout) :: m
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo
!arrays
 real(dp),intent(out) :: w(n)
 complex(dpc),intent(out) :: z(ldz,m)
 complex(dpc),intent(inout) :: a(n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 integer,allocatable :: ifail(:),iwork(:)
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_mat,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.) ! Standard LAPACK call.

  lwork = MAX(1,2*n)
  ABI_MALLOC(work,(lwork))
  ABI_MALLOC(rwork,(7*n))
  ABI_MALLOC(iwork,(5*n))
  ABI_MALLOC(ifail,(n))

  call ZHEEVX(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,ldz,work,lwork,rwork,iwork,ifail,info)

  if (info < 0) then
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEEVX had an illegal value."
    MSG_ERROR(msg)
  end if

  if (info > 0) then
    write(msg,'(2a,i0,a)')"ZHEEVX: the algorithm failed to converge; ",ch10,&
     info,"eigenvectors failed to converge. "
    MSG_ERROR(msg)
  end if

  ABI_FREE(iwork)
  ABI_FREE(ifail)
  ABI_FREE(rwork)
  ABI_FREE(work)
  RETURN

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE
   call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call slk_matrix_from_global_dpc_2D(Slk_mat,uplo,a)

   want_eigenvectors = firstchar(jobz,(/"V","v"/))
   if (want_eigenvectors) then ! Initialize the distributed vectors.
     call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   end if

   ! Solve the problem.
   call slk_pzheevx(jobz,range,uplo,Slk_mat,vl,vu,il,iu,abstol,Slk_vec,m,w)

   call destruction_matrix_scalapack(Slk_mat)

   if (want_eigenvectors) then ! A is overwritten with the eigenvectors
    z = czero
    call slk_matrix_to_global_dpc_2D(Slk_vec,"All",z) ! Fill the entries calculated by this node.
    call destruction_matrix_scalapack(Slk_vec)
    call xmpi_sum(z,comm,ierr)                        ! Fill the remaing entries of the global matrix
   end if

   call end_scalapack(Slk_processor)
   RETURN
#endif

   MSG_BUG("You should not be here!")
 END SELECT

end subroutine wrap_ZHEEVX
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/xheevx_cplex
!! NAME
!!  xheevx_cplex
!!
!! FUNCTION
!!  xheevx_cplex computes selected eigenvalues and, optionally, eigenvectors
!!  of a (real symmetric|complex Hermitian) matrix A.  Eigenvalues and eigenvectors can
!!  be selected by specifying either a range of values or a range of
!!  indices for the desired eigenvalues.
!!
!! INPUTS
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  RANGE   (input) CHARACTER*1
!!          = 'A': all eigenvalues will be found.
!!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!!                 will be found.
!!          = 'I': the IL-th through IU-th eigenvalues will be found.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  CPLEX   Size of the first dimension of the matrix A.
!!          1 for real symmetric matrix
!!          2 for complex Hermitian matrix.
!!
!!  N       (input) INTEGER
!!          The order of the matrix A.  N >= 0.
!!
!!  LDA     (input) INTEGER
!!          The leading dimension of the array A.  LDA >= max(1,N).
!!
!!  VL      (input) REAL(DP)
!!  VU      (input) REAL(DP)
!!          If RANGE='V', the lower and upper bounds of the interval to
!!          be searched for eigenvalues. VL < VU.
!!          Not referenced if RANGE = 'A' or 'I'.
!!
!!  IL      (input) INTEGER
!!  IU      (input) INTEGER
!!          If RANGE='I', the indices (in ascending order) of the
!!          smallest and largest eigenvalues to be returned.
!!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!!          Not referenced if RANGE = 'A' or 'V'.
!!
!!  ABSTOL  (input) REAL(DP)
!!          The absolute error tolerance for the eigenvalues.
!!          An approximate eigenvalue is accepted as converged
!!          when it is determined to lie in an interval [a,b]
!!          of width less than or equal to
!!
!!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!!
!!          where EPS is the machine precision.  If ABSTOL is less than
!!          or equal to zero, then  EPS*|T|  will be used in its place,
!!          where |T| is the 1-norm of the tridiagonal matrix obtained
!!          by reducing A to tridiagonal form.
!!
!!          Eigenvalues will be computed most accurately when ABSTOL is
!!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!!          If this routine returns with INFO>0, indicating that some
!!          eigenvectors did not converge, try setting ABSTOL to
!!          2*DLAMCH('S').
!!
!!          See "Computing Small Singular Values of Bidiagonal Matrices
!!          with Guaranteed High Relative Accuracy," by Demmel and
!!          Kahan, LAPACK Working Note #3.
!!
!!  LDZ     (input) INTEGER
!!          The leading dimension of the array Z.  LDZ >= 1, and if
!!          JOBZ = 'V', LDZ >= max(1,N).
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  M       (output) INTEGER
!!          The total number of eigenvalues found.  0 <= M <= N.
!!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!!
!!  W       (output) REAL(DP) array, dimension (N)
!!          On normal exit, the first M elements contain the selected
!!          eigenvalues in ascending order.
!!
!!  Z       (output) REAL(DP) array, dimension (CPLEX, LDZ, max(1,M))
!!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!!          contain the orthonormal eigenvectors of the matrix A
!!          corresponding to the selected eigenvalues, with the i-th
!!          column of Z holding the eigenvector associated with W(i).
!!          If an eigenvector fails to converge, then that column of Z
!!          contains the latest approximation to the eigenvector, and the
!!          index of the eigenvector is returned in IFAIL.
!!          If JOBZ = 'N', then Z is not referenced.
!!          Note: the user must ensure that at least max(1,M) columns are
!!          supplied in the array Z; if RANGE = 'V', the exact value of M
!!          is not known in advance and an upper bound must be used.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) REAL(DP) array, dimension (CPLEX, N, N)
!!          On entry, the (real symmetric|complex Hermitian) matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = 'L',
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!          On exit, the lower triangle (if UPLO='L') or the upper
!!          triangle (if UPLO='U') of A, including the diagonal, is
!!          destroyed.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xheevx_cplex(jobz, range, uplo, cplex, n, a, vl, vu, il, iu, &
                        abstol, m, w, z, ldz, msg, ierr, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,iu,ldz,n,cplex
 integer,optional,intent(in) :: comm
 integer,intent(inout) :: m
 integer,intent(out) :: ierr
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo
 character(len=*),intent(out) :: msg
!arrays
 real(dp),intent(out) :: w(n)
 !real(dp),intent(out) :: z(cplex,ldz,n)
 real(dp),intent(out) :: z(cplex,ldz,m)
 real(dp),intent(inout) :: a(cplex,n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,nprocs
 logical :: use_scalapack
!arrays
 integer,allocatable :: ifail(:),iwork(:)
 real(dp),allocatable :: rwork(:)
 real(dp),allocatable :: work_real(:)
 complex(dpc),allocatable :: work_cplx(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_mat,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 if (ALL(cplex/=(/1,2/))) then
   write(msg,'(a,i0)')" Wrong value for cplex: ",cplex
   ierr = 1; return
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
  ! Standard LAPACK call.

  if (cplex==1) then
    ! Real symmetric case
    lwork = MAX(1,8*n)
    ABI_MALLOC(work_real,(lwork))
    ABI_MALLOC(iwork,(5*n))
    ABI_MALLOC(ifail,(n))

    call DSYEVX(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_real,lwork,iwork,ifail,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of DSYEVX had an illegal value."
    end if

    if (ierr > 0) then
      write(msg,'(2a,i0,a)')&
       "DSYEVX: the algorithm failed to converge; ",ch10,ierr,"eigenvectors failed to converge. "
    end if

    ABI_FREE(work_real)
    ABI_FREE(iwork)
    ABI_FREE(ifail)
    RETURN

  else
    ! Complex Hermitian case.
    lwork = MAX(1,2*n)
    ABI_MALLOC(work_cplx,(lwork))
    ABI_MALLOC(rwork,(7*n))
    ABI_MALLOC(iwork,(5*n))
    ABI_MALLOC(ifail,(n))

    call ZHEEVX(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_cplx,lwork,rwork,iwork,ifail,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of ZHEEVX had an illegal value."
    end if

    if (ierr > 0) then
      write(msg,'(2a,i0,a)')&
      "ZHEEVX: the algorithm failed to converge; ",ch10,ierr,"eigenvectors failed to converge. "
    end if

    ABI_FREE(iwork)
    ABI_FREE(ifail)
    ABI_FREE(rwork)
    ABI_FREE(work_cplx)
    RETURN
  end if

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
  MSG_ERROR("Not coded yet")
  ! call init_scalapack(Slk_processor,comm)
  ! istwf_k=1

  ! ! Initialize and fill Scalapack matrix from the global one.
  ! tbloc=SLK_BLOCK_SIZE
  ! call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

  ! write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
  ! call wrtout(std_out,msg,"COLL")

  ! call slk_matrix_from_global_dpc_2D(Slk_mat,uplo,a)

  ! want_eigenvectors = firstchar(jobz,(/"V","v"/))
  ! if (want_eigenvectors) then ! Initialize the distributed vectors.
  !  call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
  ! end if

  ! ! Solve the problem.
  ! call slk_pzheevx(jobz,range,uplo,Slk_mat,vl,vu,il,iu,abstol,Slk_vec,m,w)

  ! call destruction_matrix_scalapack(Slk_mat)
  !
  ! if (want_eigenvectors) then ! A is overwritten with the eigenvectors
  !  z = czero
  !  call slk_matrix_to_global_dpc_2D(Slk_vec,"All",z) ! Fill the entries calculated by this node.
  !  call destruction_matrix_scalapack(Slk_vec)
  !  call xmpi_sum(z,comm,ierr)                        ! Fill the remaing entries of the global matrix
  ! end if

  ! call end_scalapack(Slk_processor)

  RETURN
#endif

  MSG_BUG("You should not be here!")
 END SELECT

end subroutine xheevx_cplex
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_ZHEGVX
!! NAME
!!  wrap_ZHEGVX
!!
!! FUNCTION
!!  wrap_ZHEGVX  - compute selected eigenvalues, and optionally, eigenvectors of a
!!  complex generalized Hermitian-definite eigenproblem, of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x.
!!  Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
!!  indices for the desired eigenvalues.
!!
!! INPUTS
!!
!!  ITYPE   (input) INTEGER Specifies the problem type to be solved:
!!          = 1:  A*x = (lambda)*B*x
!!          = 2:  A*B*x = (lambda)*x
!!          = 3:  B*A*x = (lambda)*x
!!
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  RANGE   (input) CHARACTER*1
!!          = 'A': all eigenvalues will be found.
!!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!!                 will be found.
!!          = 'I': the IL-th through IU-th eigenvalues will be found.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  N       (input) INTEGER
!!          The order of the matrices A and B.  N >= 0.
!!
!!  LDA     (input) INTEGER
!!          The leading dimension of the array A.  LDA >= max(1,N).
!!
!!  VL      (input) REAL(DP)
!!  VU      (input) REAL(DP)
!!          If RANGE='V', the lower and upper bounds of the interval to
!!          be searched for eigenvalues. VL < VU.
!!          Not referenced if RANGE = 'A' or 'I'.
!!
!!  IL      (input) INTEGER
!!  IU      (input) INTEGER
!!          If RANGE='I', the indices (in ascending order) of the
!!          smallest and largest eigenvalues to be returned.
!!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!!          Not referenced if RANGE = 'A' or 'V'.
!!
!!  ABSTOL  (input) REAL(DP)
!!          The absolute error tolerance for the eigenvalues.
!!          An approximate eigenvalue is accepted as converged
!!          when it is determined to lie in an interval [a,b]
!!          of width less than or equal to
!!
!!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!!
!!          where EPS is the machine precision.  If ABSTOL is less than
!!          or equal to zero, then  EPS*|T|  will be used in its place,
!!          where |T| is the 1-norm of the tridiagonal matrix obtained
!!          by reducing A to tridiagonal form.
!!
!!          Eigenvalues will be computed most accurately when ABSTOL is
!!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!!          If this routine returns with INFO>0, indicating that some
!!          eigenvectors did not converge, try setting ABSTOL to
!!          2*DLAMCH('S').
!!
!!  LDZ     (input) INTEGER
!!          The leading dimension of the array Z.  LDZ >= 1, and if
!!          JOBZ = 'V', LDZ >= max(1,N).
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  M       (output) INTEGER
!!          The total number of eigenvalues found.  0 <= M <= N.
!!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!!
!!  W       (output) REAL(DP) array, dimension (N)
!!          On normal exit, the first M elements contain the selected
!!          eigenvalues in ascending order.
!!
!!  Z       (output) COMPLEX(DPC) array, dimension (LDZ, max(1,M))
!!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!!          contain the orthonormal eigenvectors of the matrix A
!!          corresponding to the selected eigenvalues, with the i-th
!!          column of Z holding the eigenvector associated with W(i).
!!          The eigenvectors are normalized as follows: if ITYPE = 1 or 2, Z**T*B*Z = I; if ITYPE = 3, Z**T*inv(B)*Z = I.
!!          If an eigenvector fails to converge, then that column of Z
!!          contains the latest approximation to the eigenvector, and the
!!          index of the eigenvector is returned in IFAIL.
!!          If JOBZ = 'N', then Z is not referenced.
!!          Note: the user must ensure that at least max(1,M) columns are
!!          supplied in the array Z; if RANGE = 'V', the exact value of M
!!          is not known in advance and an upper bound must be used.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) COMPLEX(DPC) array, dimension (N, N)
!!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = "L",
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!
!!          On exit, the lower triangle (if UPLO="L") or the upper
!!          triangle (if UPLO="U") of A, including the diagonal, is
!!          destroyed.
!!
!!   B      (input/output) COMPLEX(DPC) array, dimension (LDB, N)
!!          On entry, the Hermitian matrix B.  If UPLO = "U", the leading N-by-N upper triangular part
!!          of B contains the upper triangular part  of the matrix B.
!!          If UPLO = "L", the leading N-by-N lower triangular part of B contains the lower triangular part of the matrix B.
!!
!!          On exit, if INFO <= N, the part of B containing the matrix is overwritten by the triangular factor
!!          U or L from the Cholesky factorization B = U**H*U or B = L*L**H.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_ZHEGVX(itype,jobz,range,uplo,n,a,b,vl,vu,il,iu,abstol,m,w,z,ldz,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,iu,ldz,n,itype
 integer,optional,intent(in) :: comm
 integer,intent(inout) :: m
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo
!arrays
 real(dp),intent(out) :: w(n)
 !complex(dpc),intent(out) :: z(ldz,n)
 complex(dpc),intent(out) :: z(ldz,m)
 complex(dpc),intent(inout) :: a(n,n),b(n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs,ii
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 integer,allocatable :: ifail(:),iwork(:)
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_matA,Slk_matB,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
   ! Standard LAPACK call.
   lwork = MAX(1,2*n)
   ABI_MALLOC(work,(lwork))
   ABI_MALLOC(rwork,(7*n))
   ABI_MALLOC(iwork,(5*n))
   ABI_MALLOC(ifail,(n))

   call ZHEGVX(itype,jobz,range,uplo,n,a,n,b,n,vl,vu,il,iu,abstol,m,w,z,ldz,work,lwork,rwork,iwork,ifail,info)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEGVX had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
     if (info<= n) then
       write(msg,'(a,i0,a)')"ZHEGVX failed to converge: ",info," eigenvectors failed to converge. "
     else
       ii = info -n
       write(msg,'(3a,i0,3a)')&
        "ZHEEVX failed to converge: ",ch10,&
        "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
        "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
     end if
     MSG_ERROR(msg)
   end if

   ABI_FREE(iwork)
   ABI_FREE(ifail)
   ABI_FREE(rwork)
   ABI_FREE(work)
   RETURN

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call init_matrix_scalapack(Slk_matA,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   call slk_matrix_from_global_dpc_2D(Slk_matA,uplo,a)

   call init_matrix_scalapack(Slk_matB,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   call slk_matrix_from_global_dpc_2D(Slk_matB,uplo,b)

   want_eigenvectors = firstchar(jobz,(/"V","v"/))
   if (want_eigenvectors) then ! Initialize the distributed vectors.
     call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
   end if

   ! Solve the problem.
   MSG_ERROR("slk_pZHEGVX not coded yet")
   ! TODO write the scaLAPACK wrapper.
   !call slk_pZHEGVX(itype,jobz,range,uplo,Slk_matA,Slk_matB,vl,vu,il,iu,abstol,Slk_vec,m,w)

   call destruction_matrix_scalapack(Slk_matA)
   call destruction_matrix_scalapack(Slk_matB)

   if (want_eigenvectors) then ! A is overwritten with the eigenvectors
     z = czero
     call slk_matrix_to_global_dpc_2D(Slk_vec,"All",z) ! Fill the entries calculated by this node.
     call destruction_matrix_scalapack(Slk_vec)
     call xmpi_sum(z,comm,ierr)                        ! Fill the remaing entries of the global matrix
   end if

   call end_scalapack(Slk_processor)

   RETURN
#endif

  MSG_BUG("You should not be here!")
 END SELECT

end subroutine wrap_ZHEGVX
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/xhegvx_cplex
!! NAME
!!  xhegvx_cplex
!!
!! FUNCTION
!!  xhegvx_cplex  - compute selected eigenvalues, and optionally, eigenvectors of a
!!  (real symmetric-definite|complex generalized Hermitian-definite) eigenproblem, of the form
!!  A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x.
!!  Here A and B are assumed to be (real symmetric|complex Hermitian) and B is also positive definite.
!!  Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
!!  indices for the desired eigenvalues.
!!
!! INPUTS
!!
!!  ITYPE   (input) INTEGER Specifies the problem type to be solved:
!!          = 1:  A*x = (lambda)*B*x
!!          = 2:  A*B*x = (lambda)*x
!!          = 3:  B*A*x = (lambda)*x
!!
!!  JOBZ    (input) CHARACTER*1
!!          = 'N':  Compute eigenvalues only;
!!          = 'V':  Compute eigenvalues and eigenvectors.
!!
!!  RANGE   (input) CHARACTER*1
!!          = 'A': all eigenvalues will be found.
!!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!!                 will be found.
!!          = 'I': the IL-th through IU-th eigenvalues will be found.
!!
!!  UPLO    (input) CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!!
!!  CPLEX   Size of the first dimension of the matrices A and B
!!          1 for Real symmetric matrices
!!          2 for complex Hermitianmatrices
!!
!!  N       (input) INTEGER
!!          The order of the matrices A and B.  N >= 0.
!!
!!  LDA     (input) INTEGER
!!          The leading dimension of the array A.  LDA >= max(1,N).
!!
!!  VL      (input) REAL(DP)
!!  VU      (input) REAL(DP)
!!          If RANGE='V', the lower and upper bounds of the interval to
!!          be searched for eigenvalues. VL < VU.
!!          Not referenced if RANGE = 'A' or 'I'.
!!
!!  IL      (input) INTEGER
!!  IU      (input) INTEGER
!!          If RANGE='I', the indices (in ascending order) of the
!!          smallest and largest eigenvalues to be returned.
!!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!!          Not referenced if RANGE = 'A' or 'V'.
!!
!!  ABSTOL  (input) REAL(DP)
!!          The absolute error tolerance for the eigenvalues.
!!          An approximate eigenvalue is accepted as converged
!!          when it is determined to lie in an interval [a,b]
!!          of width less than or equal to
!!
!!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!!
!!          where EPS is the machine precision.  If ABSTOL is less than
!!          or equal to zero, then  EPS*|T|  will be used in its place,
!!          where |T| is the 1-norm of the tridiagonal matrix obtained
!!          by reducing A to tridiagonal form.
!!
!!          Eigenvalues will be computed most accurately when ABSTOL is
!!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!!          If this routine returns with INFO>0, indicating that some
!!          eigenvectors did not converge, try setting ABSTOL to
!!          2*DLAMCH('S').
!!
!!  LDZ     (input) INTEGER
!!          The leading dimension of the array Z.  LDZ >= 1, and if
!!          JOBZ = 'V', LDZ >= max(1,N).
!!
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! OUTPUT
!!  M       (output) INTEGER
!!          The total number of eigenvalues found.  0 <= M <= N.
!!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!!
!!  W       (output) REAL(DP) array, dimension (N)
!!          On normal exit, the first M elements contain the selected
!!          eigenvalues in ascending order.
!!
!!  Z       (output) REAL(DP) array, dimension (CPLEX ,LDZ, max(1,M))
!!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!!          contain the orthonormal eigenvectors of the matrix A
!!          corresponding to the selected eigenvalues, with the i-th
!!          column of Z holding the eigenvector associated with W(i).
!!          The eigenvectors are normalized as follows:
!!           if ITYPE = 1 or 2, Z**T*B*Z = I; if ITYPE = 3, Z**T*inv(B)*Z = I.
!!
!!          If an eigenvector fails to converge, then that column of Z
!!          contains the latest approximation to the eigenvector, and the
!!          index of the eigenvector is returned in IFAIL.
!!          If JOBZ = 'N', then Z is not referenced.
!!          Note: the user must ensure that at least max(1,M) columns are
!!          supplied in the array Z; if RANGE = 'V', the exact value of M
!!          is not known in advance and an upper bound must be used.
!!
!! See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  A       (input/output) REAL(DP) array, dimension (CPLEX, N, N)
!!          On entry, the (real symmetric| complex Hermitian) matrix A.  If UPLO = 'U', the
!!          leading N-by-N upper triangular part of A contains the
!!          upper triangular part of the matrix A.  If UPLO = "L",
!!          the leading N-by-N lower triangular part of A contains
!!          the lower triangular part of the matrix A.
!!
!!          On exit, the lower triangle (if UPLO="L") or the upper
!!          triangle (if UPLO="U") of A, including the diagonal, is
!!          destroyed.
!!
!!   B      (input/output) REAL(DP) array, dimension (CPLEX, LDB, N)
!!          On entry, the (real symmetric| complex Hermitian) matrix B.  If UPLO = "U", the leading N-by-N upper triangular part
!!          of B contains the upper triangular part  of the matrix B.
!!          If UPLO = "L", the leading N-by-N lower triangular part of B contains the lower triangular part of the matrix B.
!!
!!          On exit, if INFO <= N, the part of B containing the matrix is overwritten by the triangular factor
!!          U or L from the Cholesky factorization B = U**H*U or B = L*L**H.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xhegvx_cplex(itype, jobz, range, uplo, cplex, n, a, b, &
                        vl, vu, il, iu, abstol, m, w, z, ldz, msg, ierr, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,iu,ldz,n,itype,cplex
 integer,optional,intent(in) :: comm
 integer,intent(inout) :: m
 integer,intent(out) :: ierr
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo
 character(len=*),intent(out) :: msg
!arrays
 real(dp),intent(out) :: w(n)
 !real(dp),intent(out) :: z(cplex,ldz,n)
 real(dp),intent(out) :: z(cplex,ldz,m)
 real(dp),intent(inout) :: a(cplex,n,n),b(cplex,n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,nprocs,ii
 logical :: use_scalapack
!arrays
 integer,allocatable :: ifail(:),iwork(:)
 real(dp),allocatable :: rwork(:)
 real(dp),allocatable :: work_real(:)
 complex(dpc),allocatable :: work_cplx(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,tbloc
 logical :: want_eigenvectors
 type(matrix_scalapack)    :: Slk_matA,Slk_matB,Slk_vec
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 if (ALL(cplex/=(/1,2/))) then
   write(msg,'(a,i0)')" Wrong value for cplex: ",cplex
   ierr = 1; return
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.)
  ! Standard LAPACK call.
  if (cplex==1) then
    ! Real symmetric case
    lwork = MAX(1,8*n)

    ABI_MALLOC(work_real,(lwork))
    ABI_MALLOC(iwork,(5*n))
    ABI_MALLOC(ifail,(n))

    call DSYGVX(itype,jobz,range,uplo,n,a,n,b,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_real,lwork,iwork,ifail,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')" The ",-ierr,"-th argument of DSYGVX had an illegal value."
    end if

    if (ierr > 0) then
      if (ierr<= n) then
       write(msg,'(a,i0,a)')" DSYGVX failed to converge: ",ierr," eigenvectors failed to converge. "
      else
       ii = ierr - n
       write(msg,'(3a,i0,3a)')&
        " DSYGVX failed to converge: ",ch10,&
        " The leading minor of order ",ii," of B is not positive definite. ",ch10,&
        " The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
      end if
    end if

    ABI_FREE(iwork)
    ABI_FREE(ifail)
    ABI_FREE(work_real)
    RETURN

  else
    ! Complex Hermitian case.
    lwork = MAX(1,2*n)

    ABI_MALLOC(work_cplx,(lwork))
    ABI_MALLOC(rwork,(7*n))
    ABI_MALLOC(iwork,(5*n))
    ABI_MALLOC(ifail,(n))

    !write(std_out,*)"Calling ZHEGVX"
    call ZHEGVX(itype,jobz,range,uplo,n,a,n,b,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_cplx,lwork,rwork,iwork,ifail,ierr)

    if (ierr < 0) then
      write(msg,'(a,i0,a)')"The ",-ierr,"-th argument of ZHEGVX had an illegal value."
    end if

    if (ierr > 0) then
      if (ierr<= n) then
        write(msg,'(a,i0,a)')"ZHEGVX failed to converge: ",ierr," eigenvectors failed to converge. "
      else
        ii = ierr -n
        write(msg,'(3a,i0,3a)')&
         "ZHEEVX failed to converge: ",ch10,&
         "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
         "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
      end if
    end if

    ABI_FREE(iwork)
    ABI_FREE(ifail)
    ABI_FREE(rwork)
    ABI_FREE(work_cplx)
    RETURN
  end if ! cplex

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
  MSG_ERROR("not coded yet")
  ! call init_scalapack(Slk_processor,comm)
  ! istwf_k=1

  ! ! Initialize and fill Scalapack matrix from the global one.
  ! tbloc=SLK_BLOCK_SIZE

  ! write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
  ! call wrtout(std_out,msg,"COLL")

  ! call init_matrix_scalapack(Slk_matA,n,n,Slk_processor,istwf_k,tbloc=tbloc)
  ! call slk_matrix_from_global_dpc_2D(Slk_matA,uplo,a)

  ! call init_matrix_scalapack(Slk_matB,n,n,Slk_processor,istwf_k,tbloc=tbloc)
  ! call slk_matrix_from_global_dpc_2D(Slk_matB,uplo,b)

  ! want_eigenvectors = firstchar(jobz,(/"V","v"/))
  ! if (want_eigenvectors) then ! Initialize the distributed vectors.
  !  call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc)
  ! end if

  ! ! Solve the problem.
  ! MSG_ERROR("slk_pZHEGVX not coded yet")
  ! ! TODO write the scaLAPACK wrapper.
  ! call slk_pZHEGVX(itype,jobz,range,uplo,Slk_matA,Slk_matB,vl,vu,il,iu,abstol,Slk_vec,m,w)

  ! call destruction_matrix_scalapack(Slk_matA)
  ! call destruction_matrix_scalapack(Slk_matB)
  !
  ! if (want_eigenvectors) then ! A is overwritten with the eigenvectors
  !  z = czero
  !  call slk_matrix_to_global_dpc_2D(Slk_vec,"All",z) ! Fill the entries calculated by this node.
  !  call destruction_matrix_scalapack(Slk_vec)
  !  call xmpi_sum(z,comm,ierr)                        ! Fill the remaing entries of the global matrix
  ! end if

  ! call end_scalapack(Slk_processor)

  ! RETURN
#endif

  MSG_BUG("You should not be here!")
 END SELECT

end subroutine xhegvx_cplex
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_CGEEV
!! NAME
!!  wrap_CGEEV
!!
!! FUNCTION
!!  wrap_CGEEV computes for an N-by-N complex nonsymmetric matrix A, the
!!  eigenvalues and, optionally, the left and/or right eigenvectors using single precision arithmetic. [PRIVATE]
!!
!!  The right eigenvector v(j) of A satisfies: A * v(j) = lambda(j) * v(j)
!!  where lambda(j) is its eigenvalue.
!!  The left eigenvector u(j) of A satisfies u(j)**H * A = lambda(j) * u(j)**H
!!  where u(j)**H denotes the conjugate transpose of u(j).
!!
!!  The computed eigenvectors are normalized to have Euclidean norm
!!  equal to 1 and largest component real.
!!
!! INPUTS
!!   JOBVL   (input) CHARACTER*1
!!           = 'N': left eigenvectors of A are not computed;
!!           = 'V': left eigenvectors of are computed.
!!
!!   JOBVR   (input) CHARACTER*1
!!           = 'N': right eigenvectors of A are not computed;
!!           = 'V': right eigenvectors of A are computed.
!!
!!   N       (input) INTEGER
!!           The order of the matrix A. N >= 0.
!!
!!   LDA     (input) INTEGER
!!           The leading dimension of the array A.  LDA >= max(1,N).
!!
!!   LDVL    (input) INTEGER
!!           The leading dimension of the array VL.  LDVL >= 1; if
!!           JOBVL = 'V', LDVL >= N.
!!
!!   LDVR    (input) INTEGER
!!           The leading dimension of the array VR.  LDVR >= 1; if
!!           JOBVR = 'V', LDVR >= N.
!!
!! OUTPUT
!!   W       (output) COMPLEX(SPC) array, dimension (N)
!!           W contains the computed eigenvalues.
!!   VL      (output) COMPLEX(SCP) array, dimension (LDVL,N)
!!           If JOBVL = 'V', the left eigenvectors u(j) are stored one
!!           after another in the columns of VL, in the same order
!!           as their eigenvalues.
!!           If JOBVL = 'N', VL is not referenced.
!!           u(j) = VL(:,j), the j-th column of VL.
!!   VR      (output) COMPLEX(SPC) array, dimension (LDVR,N)
!!           If JOBVR = 'V', the right eigenvectors v(j) are stored one
!!           after another in the columns of VR, in the same order
!!           as their eigenvalues.
!!           If JOBVR = 'N', VR is not referenced.
!!           v(j) = VR(:,j), the j-th column of VR.
!!
!!  See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!   A       (input/output) COMPLEX(SPC) array, dimension (LDA,N)
!!           On entry, the N-by-N matrix A.
!!           On exit, A has been overwritten.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_CGEEV(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,lda,ldvl,ldvr
 character(len=*),intent(in) ::  jobvl,jobvr
!arrays
 complex(spc),intent(inout) :: a(lda,n)
 complex(spc),intent(out) :: w(n)
 complex(spc),intent(out) :: vl(ldvl,n)
 complex(spc),intent(out) :: vr(ldvr,n)

!Local variables ------------------------------
!scalars
 integer :: info,lwork
 character(len=500) :: msg
!arrays
 real(sp),allocatable :: rwork(:)
 complex(spc),allocatable :: work(:)

!************************************************************************

 lwork = MAX(1,2*n)

 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(rwork,(2*n))

 call CGEEV(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

 if (info < 0) then
   write(msg,'(a,i0,a)')" The ",-info,"-th argument of CGEEV had an illegal value."
   MSG_ERROR(msg)
 end if

 if (info > 0) then
   write(msg,'(3a,i0,a,i0,a)')&
     "CGEEV: The QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed;",ch10,&
     "Elements ",info+1,":",n," of W contain eigenvalues which have converged. "
   MSG_ERROR(msg)
 end if

 ABI_FREE(work)
 ABI_FREE(rwork)

end subroutine wrap_CGEEV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/wrap_ZGEEV
!! NAME
!!  wrap_ZGEEV
!!
!! FUNCTION
!!  wrap_ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
!!  eigenvalues and, optionally, the left and/or right eigenvectors using double precision arithmetic. [PRIVATE]
!!
!!  The right eigenvector v(j) of A satisfies: A * v(j) = lambda(j) * v(j)
!!  where lambda(j) is its eigenvalue.
!!  The left eigenvector u(j) of A satisfies u(j)**H * A = lambda(j) * u(j)**H
!!  where u(j)**H denotes the conjugate transpose of u(j).
!!
!!  The computed eigenvectors are normalized to have Euclidean norm
!!  equal to 1 and largest component real.
!!  No scalapack version is available (PZGEEV is not provided by the Scalapack team)
!!
!! INPUTS
!!   JOBVL   (input) CHARACTER*1
!!           = 'N': left eigenvectors of A are not computed;
!!           = 'V': left eigenvectors of are computed.
!!
!!   JOBVR   (input) CHARACTER*1
!!           = 'N': right eigenvectors of A are not computed;
!!           = 'V': right eigenvectors of A are computed.
!!
!!   N       (input) INTEGER
!!           The order of the matrix A. N >= 0.
!!
!!   LDA     (input) INTEGER
!!           The leading dimension of the array A.  LDA >= max(1,N).
!!
!!   LDVL    (input) INTEGER
!!           The leading dimension of the array VL.  LDVL >= 1; if
!!           JOBVL = 'V', LDVL >= N.
!!
!!   LDVR    (input) INTEGER
!!           The leading dimension of the array VR.  LDVR >= 1; if
!!           JOBVR = 'V', LDVR >= N.
!!
!! OUTPUT
!!   W       (output) COMPLEX(DPC) array, dimension (N)
!!           W contains the computed eigenvalues.
!!   VL      (output) COMPLEX(DPC) array, dimension (LDVL,N)
!!           If JOBVL = 'V', the left eigenvectors u(j) are stored one
!!           after another in the columns of VL, in the same order
!!           as their eigenvalues.
!!           If JOBVL = 'N', VL is not referenced.
!!           u(j) = VL(:,j), the j-th column of VL.
!!   VR      (output) COMPLEX(DPC) array, dimension (LDVR,N)
!!           If JOBVR = 'V', the right eigenvectors v(j) are stored one
!!           after another in the columns of VR, in the same order
!!           as their eigenvalues.
!!           If JOBVR = 'N', VR is not referenced.
!!           v(j) = VR(:,j), the j-th column of VR.
!!
!!  See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!   A       (input/output) COMPLEX(DPC) array, dimension (LDA,N)
!!           On entry, the N-by-N matrix A.
!!           On exit, A has been overwritten.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrap_ZGEEV(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,lda,ldvl,ldvr
 character(len=*),intent(in) ::  jobvl,jobvr
!arrays
 complex(dpc),intent(inout) :: a(lda,n)
 complex(dpc),intent(out) :: w(n)
 complex(dpc),intent(out) :: vl(ldvl,n)
 complex(dpc),intent(out) :: vr(ldvr,n)

!Local variables ------------------------------
!scalars
 integer :: info,lwork
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rwork(:)
 complex(dpc),allocatable :: work(:)

!************************************************************************

 use_scalapack=.FALSE.

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)

   lwork = MAX(1,2*n)
   ABI_MALLOC(work,(lwork))
   ABI_MALLOC(rwork,(2*n))

   call ZGEEV(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

   if (info < 0) then
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZGEEV had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(3a,i0,a,i0,a)')&
     "ZGEEV: The QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed;",ch10,&
     "Elements ",info+1,":",n," of W contain eigenvalues which have converged. "
    MSG_ERROR(msg)
   end if

   ABI_FREE(work)
   ABI_FREE(rwork)
   RETURN

 CASE (.TRUE.)
   MSG_BUG("You should not be here!")
 END SELECT

end subroutine wrap_ZGEEV
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/cginv
!! NAME
!! cginv
!!
!! FUNCTION
!! Invert a general matrix of complex elements in single precision.
!!  CGETRF computes an LU factorization of a general N-by-N matrix A using partial pivoting with row interchanges.
!!  CGETRI computes the inverse of a matrix using the LU factorization computed by CGETRF.
!!
!! INPUTS
!! n=size of complex matrix a
!! a=matrix of complex elements
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1,
!!        in this case the sequential LAPACK routine is called.
!!
!! SIDE EFFECTS
!! a(n,n)= array of complex elements, input, inverted at output
!!
!! TODO
!!  Add Scalapack version, matrix_scalapack has to be modified by adding a single precision complex buffer.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cginv(a, n, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,optional,intent(in) :: comm
!arrays
 complex(spc),intent(inout) :: a(n,n)

!Local variables-------------------------------
!scalars
 integer :: lwork,info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 integer,allocatable :: ipiv(:)
 complex(spc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: ierr,istwf_k,ipiv_size,liwork,tbloc
 integer,allocatable :: iwork(:)
 type(matrix_scalapack)    :: Slk_mat
 type(processor_scalapack) :: Slk_processor
#endif

! *************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
  nprocs = xmpi_comm_size(comm)
  ! TODO
!#ifdef HAVE_LINALG_SCALAPACK
!  use_scalapack = (nprocs>1)
!#endif
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.)
   ABI_MALLOC(ipiv, (n))

   call CGETRF(n,n,a,n,ipiv,info) ! P* L* U  Factorization.

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of CGETRF had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
     write(msg,'(3a,i0,4a)')&
      "The matrix that has been passed in argument is probably either singular or nearly singular.",ch10,&
      "U(i,i) in the P*L*U factorization is exactly zero for i = ",info,ch10,&
      "The factorization has been completed but the factor U is exactly singular.",ch10,&
      "Division by zero will occur if it is used to solve a system of equations."
     MSG_ERROR(msg)
   end if

   lwork=MAX(1,n)
   ABI_MALLOC(work,(lwork))

   call CGETRI(n,a,n,ipiv,work,lwork,info) ! Inverts U and the computes inv(A)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of CGETRI had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
     write(msg,'(3a,i0,a)')&
      "The matrix that has been passed to this subroutine is probably either singular or nearly singular.",ch10,&
      "U(i,i) for i= ",info," is exactly zero; the matrix is singular and its inverse could not be computed."
     MSG_ERROR(msg)
   end if

   ABI_FREE(ipiv)
   ABI_FREE(work)
   RETURN

 CASE (.TRUE.)

#if 0
! FIXME matrix_scalapack does not have a single precision complex buffer

#ifdef HAVE_LINALG_SCALAPACK
  call init_scalapack(Slk_processor,comm)
  istwf_k=1

  ! Initialize and fill Scalapack matrix from the global one.
  tbloc=SLK_BLOCK_SIZE
  call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

  write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
  call wrtout(std_out,msg,"COLL")

  ! IMPORTANT NOTE: PZGETRF requires square block decomposition i.e.,  MB_A = NB_A.
  if ( Slk_mat%descript%tab(MB_)/=Slk_mat%descript%tab(NB_) ) then
   msg ="PZGETRF requires square block decomposition i.e.,  MB_A = NB_A."
   MSG_ERROR(msg)
  end if

  !!call slk_matrix_from_global_dpc_2D(Slk_mat,"All",a)

  ipiv_size = my_locr(Slk_mat) + Slk_mat%descript%tab(MB_)
  ABI_MALLOC(ipiv,(ipiv_size))

  call PCGETRF(Slk_mat%sizeb_global(1),Slk_mat%sizeb_global(2),Slk_mat%buffer_cplx_sp,&
&   1,1,Slk_mat%descript%tab,ipiv,info) ! P * L * U  Factorization.

  if (info/=0) then
   write(msg,'(a,i0)')"PCGETRF returned info= ",info
   MSG_ERROR(msg)
  end if

  ! Get optimal size of workspace for PCGETRI.
  lwork=-1; liwork=-1
  ABI_MALLOC(work,(1))
  ABI_MALLOC(iwork,(1))

  call PCGETRI(Slk_mat%sizeb_global(1),Slk_mat%buffer_cplx_sp,1,1,Slk_mat%descript%tab,ipiv,&
&  work,lwork,iwork,liwork,info)

  ABI_CHECK(info==0,"PZGETRI: Error during compuation of workspace size")

  lwork = NINT(DBLE(work(1))); liwork=iwork(1)
  ABI_FREE(work)
  ABI_FREE(iwork)

  ! Solve the problem.
  ABI_MALLOC(work,(lwork))
  ABI_MALLOC(iwork,(liwork))

  call PCGETRI(Slk_mat%sizeb_global(1),Slk_mat%buffer_cplx_sp,1,1,Slk_mat%descript%tab,ipiv,&
&  work,lwork,iwork,liwork,info)

  if (info/=0) then
   write(msg,'(a,i0)')"PZGETRI returned info= ",info
   MSG_ERROR(msg)
  end if

  ABI_FREE(work)
  ABI_FREE(iwork)
  ABI_FREE(ipiv)

  ! Reconstruct the global matrix from the distributed one.
  a = czero
  !! call slk_matrix_to_global_dpc_2D(Slk_mat,"All",a)  ! Fill the entries calculated by this node.
  call destruction_matrix_scalapack(Slk_mat)

  call xmpi_sum(a,comm,ierr)                         ! Fill the remaing entries of the global matrix
  call end_scalapack(Slk_processor)

  RETURN
#endif

#endif

  MSG_BUG("You should not be here!")

 END SELECT

end subroutine cginv
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/zginv
!! NAME
!! zginv
!!
!! FUNCTION
!! Invert a general matrix of complex elements in double precision.
!!  ZGETRF computes an LU factorization of a general N-by-N matrix A using partial pivoting with row interchanges.
!!  ZGETRI computes the inverse of a matrix using the LU factorization computed by ZGETRF.
!!
!! INPUTS
!! n=size of complex matrix a
!! a=matrix of complex elements
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1.
!!        In this case the sequential LAPACK routine is called.
!!
!! SIDE EFFECTS
!! a(n,n)= array of complex elements, input, inverted at output
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine zginv(a, n, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,optional,intent(in) :: comm
!arrays
 complex(dpc),intent(inout) :: a(n,n)

!Local variables-------------------------------
!scalars
 integer :: lwork,info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 integer,allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:)
#ifdef HAVE_LINALG_SCALAPACK
 integer :: istwf_k,tbloc,ierr
 type(matrix_scalapack)    :: Slk_mat
 type(processor_scalapack) :: Slk_processor
#endif

! *************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
   ABI_MALLOC(ipiv, (n))
   call ZGETRF(n,n,a,n,ipiv,info) ! P* L* U  Factorization.

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZGETRF had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(3a,i0,4a)')&
      "The matrix that has been passed in argument is probably either singular or nearly singular.",ch10,&
      "U(i,i) in the P*L*U factorization is exactly zero for i = ",info,ch10,&
      "The factorization has been completed but the factor U is exactly singular.",ch10,&
      "Division by zero will occur if it is used to solve a system of equations."
    MSG_ERROR(msg)
   end if

   lwork=MAX(1,n)
   ABI_MALLOC(work,(lwork))

   call ZGETRI(n,a,n,ipiv,work,lwork,info) ! Invert U and then compute inv(A)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZGETRI had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(3a,i0,a)')&
      "The matrix that has been passed to this subroutine is probably either singular or nearly singular.",ch10,&
      "U(i,i) for i= ",info," is exactly zero; the matrix is singular and its inverse could not be computed."
    MSG_ERROR(msg)
   end if

   ABI_FREE(ipiv)
   ABI_FREE(work)
   RETURN

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE
   call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call slk_matrix_from_global_dpc_2D(Slk_mat,"All",a)

   ! Perform the calculation with scaLAPACK.
   call slk_zinvert(Slk_mat)

   ! Reconstruct the global matrix from the distributed one.
   a = czero
   call slk_matrix_to_global_dpc_2D(Slk_mat,"All",a)  ! Fill the entries calculated by this node.
   call destruction_matrix_scalapack(Slk_mat)

   call xmpi_sum(a,comm,ierr)                         ! Fill the remaing entries of the global matrix
   call end_scalapack(Slk_processor)

   RETURN
#endif

  MSG_BUG("You should not be here!")
 END SELECT

end subroutine zginv
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/zhpd_invert
!! NAME
!! zhpd_invert
!!
!! FUNCTION
!! Invert a Hermitian positive definite matrix of complex elements in double precision.
!!
!! INPUTS
!! uplo= 'U':  Upper triangle of A is stored;
!!     = 'L':  Lower triangle of A is stored.
!! n=size of complex matrix a
!! a=matrix of complex elements
!! [comm]=MPI communicator for ScaLAPACK inversion. Only available if the code has been compiled with Scalapack support.
!!        To avoid wasting CPU time the scalapack initialization is avoided if the number of processors in 1.
!!        In this case the sequential LAPACK routine is called.
!!
!! SIDE EFFECTS
!! a(n,n)=
!!    On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!!    N-by-N upper triangular part of A contains the upper
!!    triangular part of the matrix A, and the strictly lower
!!    triangular part of A is not referenced.  If UPLO = 'L', the
!!    leading N-by-N lower triangular part of A contains the lower
!!    triangular part of the matrix A, and the strictly upper
!!    triangular part of A is not referenced.
!!    On exit, the upper or lower triangle of the (Hermitian) inverse of A
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine zhpd_invert(uplo, a, n, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: uplo
 integer,intent(in) :: n
 integer,optional,intent(in) :: comm
!arrays
 complex(dpc),intent(inout) :: a(n,n)

!Local variables-------------------------------
!scalars
 integer :: info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
#ifdef HAVE_LINALG_SCALAPACK
 integer :: istwf_k,tbloc,ierr
 type(matrix_scalapack)    :: Slk_mat
 type(processor_scalapack) :: Slk_processor
#endif

! *************************************************************************

 use_scalapack=.FALSE.
 if (PRESENT(comm)) then
   nprocs = xmpi_comm_size(comm)
#ifdef HAVE_LINALG_SCALAPACK
   use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)
 CASE (.FALSE.)
   ! *  ZPOTRF computes the Cholesky factorization of a complex Hermitian positive definite.
   ! *     A = U**H * U,  if UPLO = 'U', or
   ! *     A = L  * L**H,  if UPLO = 'L',
   call ZPOTRF(uplo,n,a,n,info)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZPOTRF had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(a,i0,3a)')&
      "The leading minor of order ",info," is not positive definite, ",ch10,&
      "and the factorization could not be completed."
    MSG_ERROR(msg)
   end if
   !
   ! *  ZPOTRI computes the inverse of a complex Hermitian positive definite
   ! *  matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
   ! *  computed by ZPOTRF.
   ! *  On exit, the upper or lower triangle of the (Hermitian)
   ! *  inverse of A, overwriting the input factor U or L.
   call ZPOTRI(uplo,n,a,n,info)

   if (info < 0) then
     write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZPOTRI had an illegal value."
     MSG_ERROR(msg)
   end if

   if (info > 0) then
     write(msg,'(a,2(1x,i0),a)')&
       "The ( ",info,info,")element of the factor U or L is zero, and the inverse could not be computed."
     MSG_ERROR(msg)
   end if

   RETURN

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK
   call init_scalapack(Slk_processor,comm)
   istwf_k=1

   ! Initialize and fill Scalapack matrix from the global one.
   tbloc=SLK_BLOCK_SIZE
   call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)

   write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
   call wrtout(std_out,msg,"COLL")

   call slk_matrix_from_global_dpc_2D(Slk_mat,uplo,a)

   ! Perform the calculation with scaLAPACK.
   call slk_zdhp_invert(Slk_mat,uplo)

   ! Reconstruct the global matrix from the distributed one.
   a = czero
   call slk_matrix_to_global_dpc_2D(Slk_mat,uplo,a)  ! Fill the entries calculated by this node.
   call destruction_matrix_scalapack(Slk_mat)

   call xmpi_sum(a,comm,ierr)                         ! Fill the remaing entries of the global matrix
   call end_scalapack(Slk_processor)

   RETURN
#endif

   MSG_BUG("You should not be here!")
 END SELECT

end subroutine zhpd_invert
!!***

!----------------------------------------------------------------------

!!****f* m_hide_lapack/matrginv
!! NAME
!! matrginv
!!
!! FUNCTION
!! Invert a general matrix of real*8 elements.
!!
!! INPUTS
!! lda=leading dimension of complex matrix a
!! n=size of complex matrix a
!! a=matrix of real elements
!! OUTPUT
!! a=inverse of a input matrix
!!
!! SIDE EFFECTS
!! a(lda,n)= array of real elements, input, inverted at output
!!
!! PARENTS
!!      m_a2ftr,m_bethe_salpeter,m_ddb_elast,m_ddb_piezo,m_geometry,m_haydock
!!      m_mlwfovlp,m_paw_optics,m_symtk,m_vcoul,m_wfd_optic
!!
!! CHILDREN
!!
!! SOURCE

subroutine matrginv(a,lda,n)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lda,n
!arrays
 real(dp),intent(inout) :: a(lda,n)

!Local variables-------------------------------
!scalars
 integer :: ierr,nwork
#if defined HAVE_LINALG_ESSL
 real(dp) :: rcond
#endif
 character(len=500) :: message
!arrays
 integer,allocatable :: ipvt(:)
#if defined HAVE_LINALG_ESSL
 real(dp) :: det(2)
#elif defined HAVE_LINALG_ASL
 real(dp) :: det(2)
#endif
 real(dp),allocatable :: work(:)

! *************************************************************************

#if defined HAVE_LINALG_ESSL
 nwork=200*n
#else
 nwork=n
#endif

 ABI_ALLOCATE(work,(nwork))
 ABI_ALLOCATE(ipvt,(n))

#if defined HAVE_LINALG_ESSL

 call dgeicd(a,lda,n,0,rcond,det,work,nwork)
 if(abs(rcond)==zero) then
   write(message, '(7a)' )&
   '  The matrix that has been passed in argument of this subroutine',ch10,&
   '  is probably either singular or nearly singular.',ch10,&
   '  The ESSL routine dgeicd failed.',ch10,&
   '  Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

#elif defined HAVE_LINALG_ASL

 call dbgmlu(a,lda,n,ipvt,ierr)
 if(ierr /= 0) then
   write(message, '(7a)' ) ch10,&
   '  The matrix that has been passed in argument of this subroutine',ch10,&
   '  is probably either singular or nearly singular.',ch10,&
   '  The ASL routine dbgmlu failed.',ch10,&
   '  Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

 call dbgmdi(a,lda,n,ipvt,det,-1,work,ierr)

 if(ierr /= 0) then
   write(message, '(7a)' ) &
   '  The matrix that has been passed in argument of this subroutine',ch10,&
   '  is probably either singular or nearly singular.',ch10,&
   '  The ASL routine dbgmdi failed.',ch10,&
   '  Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

#else

 call dgetrf(n,n,a,lda,ipvt,ierr)
 if(ierr /= 0) then
   write(message, '(7a)' ) &
   '  The matrix that has been passed in argument of this subroutine',ch10,&
   '  is probably either singular or nearly singular.',ch10,&
   '  The LAPACK routine dgetrf failed.',ch10,&
   '  Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

 call dgetri(n,a,lda,ipvt,work,n,ierr)

 if(ierr /= 0) then
   write(message, '(7a)' ) &
   '  The matrix that has been passed in argument of this subroutine',ch10,&
   '  is probably either singular or nearly singular.',ch10,&
   '  The LAPACK routine dgetri failed.',ch10,&
   '  Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

#endif

 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(ipvt)

end subroutine matrginv
!!***

!!****f* m_hide_lapack/matr3eigval
!! NAME
!! matr3eigval
!!
!! FUNCTION
!! Find the eigenvalues of a real symmetric 3x3 matrix, entered in full storage mode.
!!
!! INPUTS
!!  matr(3,3)=real symmetric 3x3 matrix
!!
!! OUTPUT
!!  eigval(3)=three eigenvalues
!!
!! PARENTS
!!      m_geometry
!!
!! CHILDREN
!!
!! SOURCE

subroutine matr3eigval(eigval,matr)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: matr(3,3)
 real(dp),intent(out) :: eigval(3)

!Local variables-------------------------------
!scalars
 integer :: ier
!arrays
 real(dp) :: eigvec(2,3,3),matrx(2,6),zhpev1(2,2*3-1),zhpev2(3*3-2)

! *************************************************************************

 matrx(1,1)=matr(1,1)
 matrx(1,2)=matr(1,2)
 matrx(1,3)=matr(2,2)
 matrx(1,4)=matr(1,3)
 matrx(1,5)=matr(2,3)
 matrx(1,6)=matr(3,3)
 matrx(2,:)=zero

 call ZHPEV ('V','U',3,matrx,eigval,eigvec,3,zhpev1,zhpev2,ier)
!write(std_out,*)' eigval=',eigval

end subroutine matr3eigval
!!***


!!****f* ABINIT/jacobi
!! NAME
!!  jacobi
!!
!! FUNCTION
!!  Computes all eigenvalues and eigenvectors of a real symmetric matrix a,
!!  which is of size n by n, stored in a physical np by np array. On output,
!!  elements of a above the diagonal are destroyed. d returns the
!!  eigenvalues of a in its first n elements. v is a matrix with the same
!!  logical and physical dimensions as a, whose columns contain, on output,
!!  the normalized eigenvectors of a. nrot returns the number of Jacobi
!!  rotations that were required.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!  This routine is deprecated, use Lapack API
!!
!! PARENTS
!!      m_bader,m_conducti
!!
!! CHILDREN
!!
!! SOURCE

subroutine jacobi(a,n,np,d,v,nrot)

!Arguments
 integer :: n,np,nrot
 real*8 :: a(np,np),d(np),v(np,np)
!Local variables
 integer, parameter :: NMAX=500
 integer i,ip,iq,j
 real*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
 do ip=1,n
   do iq=1,n
     v(ip,iq)=0.
   enddo
   v(ip,ip)=1.
 enddo
 do ip=1,n
   b(ip)=a(ip,ip)
   d(ip)=b(ip)
   z(ip)=0.
 enddo
 nrot=0
 do i=1,50
   sm=0.
   do ip=1,n-1
     do iq=ip+1,n
       sm=sm+abs(a(ip,iq))
     enddo
   enddo
   if(sm.eq.0.)return
   if(i.lt.4)then
     tresh=0.2*sm/n**2
   else
     tresh=0.
   endif
   do ip=1,n-1
     do iq=ip+1,n
       g=100.*abs(a(ip,iq))
       if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
&          .and.(abs(d(iq))+g.eq.abs(d(iq))))then
            a(ip,iq)=0.
       else if(abs(a(ip,iq)).gt.tresh)then
         h=d(iq)-d(ip)
         if(abs(h)+g.eq.abs(h))then
           t=a(ip,iq)/h
         else
           theta=0.5*h/a(ip,iq)
           t=1./(abs(theta)+sqrt(1.+theta**2))
           if(theta.lt.0.)t=-t
         endif
         c=1./sqrt(1+t**2)
         s=t*c
         tau=s/(1.+c)
         h=t*a(ip,iq)
         z(ip)=z(ip)-h
         z(iq)=z(iq)+h
         d(ip)=d(ip)-h
         d(iq)=d(iq)+h
         a(ip,iq)=0.
         do j=1,ip-1
           g=a(j,ip)
           h=a(j,iq)
           a(j,ip)=g-s*(h+g*tau)
           a(j,iq)=h+s*(g-h*tau)
         enddo
         do j=ip+1,iq-1
           g=a(ip,j)
           h=a(j,iq)
           a(ip,j)=g-s*(h+g*tau)
           a(j,iq)=h+s*(g-h*tau)
         enddo
         do j=iq+1,n
           g=a(ip,j)
           h=a(iq,j)
           a(ip,j)=g-s*(h+g*tau)
           a(iq,j)=h+s*(g-h*tau)
         enddo
         do j=1,n
           g=v(j,ip)
           h=v(j,iq)
           v(j,ip)=g-s*(h+g*tau)
           v(j,iq)=h+s*(g-h*tau)
         enddo
         nrot=nrot+1
       endif
     enddo
   enddo
   do ip=1,n
     b(ip)=b(ip)+z(ip)
     d(ip)=b(ip)
     z(ip)=0.
   enddo
 enddo
 write(std_out,*) 'too many iterations in jacobi'

end subroutine jacobi
!!***

!!****f* m_hide_lapack/ludcmp
!! NAME
!!  ludcmp
!!
!! FUNCTION
!!  Given a matrix a(1:n,1:n), with physical dimension np by np, this
!!  routine replaces it by the LU decomposition of a rowwise permutation of
!!  itself. a and n are input. a is output, arranged as in equation (2.3.14)
!!  above; indx(1:n) is an output vector that records the row permutation
!!  effected by the partial pivoting; id is output as +- 1 depending on
!!  whether the number of row interchanges was even or odd,
!!  respectively. This routine is used in combination with lubksb to solve
!!  linear equations or invert a matrix.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!   This routine is depreacted, use lapack API
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

SUBROUTINE ludcmp(a,n,np,indx,id,info)

      INTEGER n,np,indx(n),NMAX,id,info
      REAL*8 a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)

      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)

!      write(std_out,*) 'ENTERING LUDCMP...'
!      write(std_out,*) 'in ludcmp n=',n,' np=',np
!      write(std_out,201) ((a(i,j),j=1,n),i=1,n)
! 201  FORMAT('A in ludcmp ',/,3F16.8,/,3F16.8,/,3F16.8)
      id=1
      info=0
      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) then
          write(std_out,*) 'LUDCMP: singular matrix !!!'
          do j=1,3
            write(std_out,*) (a(j,k),k=1,3)
          enddo
          info=1
          return
!          stop 'singular matrix in ludcmp'
        endif
        vv(i)=1./aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax)then
          do  k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          id=-id
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
!      write(std_out,*) 'LEAVING LUDCMP...'
      return
END SUBROUTINE ludcmp
!!***

!!****f* m_hide_lapack/lubksb
!! NAME
!!  lubksb
!!
!! FUNCTION
!!  Solves the set of n linear equations A . X = B. Here a is input, not as
!!  the matrix A but rather as its LU decomposition, determined by the
!!  routine ludcmp. indx is input as the permutation vector returned by
!!  ludcmp. b(1:n) is input as the right-hand side vector B, and returns
!!  with the solution vector X. a, n, np, and indx are not modified by this
!!  routine and can be left in place for successive calls with different
!!  right-hand sides b. This routine takes into account the possibility that
!!  b will begin with many zero elements, so it is efficient for use in
!!  matrix inversion.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  This routine is deprecated, use lapack API
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

SUBROUTINE lubksb(a,n,np,indx,b)

      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)

      INTEGER i,ii,j,ll
      REAL*8 sum
!      write(std_out,*) 'ENTERING LUBKSB...'
!      write(std_out,201) ((a(i,j),j=1,n),i=1,n)
! 201  FORMAT('A in lubksb ',/,3F16.8,/,3F16.8,/,3F16.8)

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
!      write(std_out,*) 'LEAVING LUBKSB...'
      return

END SUBROUTINE LUBKSB
!!***

!!****f* m_hide_lapack/dzegdi
!! NAME
!!  dzgedi
!!
!! FUNCTION
!!  This routine is the clone of zgefa.F90 using real*8 a(2) instead of complex*16
!!  for the purpose of ABINIT
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_berryphase,m_berrytk,m_dfpt_fef,m_elpolariz,m_pead_nl_loop,m_relaxpol
!!
!! CHILDREN
!!
!! SOURCE

subroutine dzgedi(a,lda,n,ipvt,det,work,job)

      integer :: lda,n,ipvt(n),job
      real*8 :: a(2,lda,n),det(2,2),work(2,n)
!
!     zgedi computes the determinant and inverse of a matrix
!     using the factors computed by zgeco or zgefa.
!
!     on entry
!
!        a       complex*16(lda, n)
!                the output from zgeco or zgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from zgeco or zgefa.
!
!        work    complex*16(n)
!                work vector.  contents destroyed.
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!
!        det     complex*16(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. cabs1(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if zgeco has set rcond .gt. 0.0 or zgefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     internal variables
!
      double precision :: r(2),rk(2),rkj(2)
      double precision :: ten,rinv2,rabs
      integer :: i,j,k,kb,kp1,l,nm1
!
!     compute determinant
!
      if (job/10 .eq. 0) go to 70
         det(1,1) = 1.0d0; det(2,1) = 0.0d0
         det(1,2) = 0.0d0; det(2,2) = 0.0d0
         ten = 10.0d0
         do i = 1, n
            if (ipvt(i) .ne. i) then
                det(1,1) = -det(1,1)
                det(2,1) = -det(2,1)
            end if
            r(1)=det(1,1); r(2)=det(2,1)
            det(1,1) = r(1)*a(1,i,i)-r(2)*a(2,i,i)
            det(2,1) = r(2)*a(1,i,i)+r(1)*a(2,i,i)
!        ...exit
            rabs = abs(det(1,1))+abs(det(2,1))
            if (rabs .eq. 0.0d0) go to 60
   10       continue
            rabs = abs(det(1,1))+abs(det(2,1))
            if (rabs .ge. 1.0d0) go to 20
               det(1,1) = ten*det(1,1); det(2,1) = ten*det(2,1)
               det(1,2) = det(1,2) - 1.0d0
            go to 10
   20       continue
   30       continue
            rabs = abs(det(1,1))+abs(det(2,1))
            if (rabs .lt. ten) go to 40
               det(1,1) = det(1,1)/ten; det(2,1) = det(2,1)/ten
               det(1,2) = det(1,2) + 1.0d0
            go to 30
   40       continue
         end do
   60    continue
   70 continue
!
!     compute inverse(u)
!
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            !a(k,k) = (1.0d0,0.0d0)/a(k,k)
            !t = -a(k,k)
            !call zscal(k-1,t,a(1,k),1)
            rinv2 = 1.d0/(a(1,k,k)**2+a(2,k,k)**2)
            a(1,k,k) =  rinv2*a(1,k,k)
            a(2,k,k) = -rinv2*a(2,k,k)
            rk(1) = -a(1,k,k); rk(2) = -a(2,k,k)
            do i=1,k-1
               r(1)=a(1,i,k)
               r(2)=a(2,i,k)
               a(1,i,k)=rk(1)*r(1)-rk(2)*r(2)
               a(2,i,k)=rk(1)*r(2)+rk(2)*r(1)
            end do
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               !t = a(k,j)
               !a(k,j) = (0.0d0,0.0d0)
               !call zaxpy(k,t,a(1,k),1,a(1,j),1)
               rkj(1) = a(1,k,j); rkj(2) = a(2,k,j)
               a(1,k,j) = 0.d0; a(2,k,j) = 0.d0
               do i=1,k
                  a(1,i,j)=rkj(1)*a(1,i,k)-rkj(2)*a(2,i,k)+a(1,i,j)
                  a(2,i,j)=rkj(2)*a(1,i,k)+rkj(1)*a(2,i,k)+a(2,i,j)
               end do
   80       continue
   90       continue
  100    continue
  do i=1,n
  end do
!
!        form inverse(u)*inverse(l)
!
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(1,i) = a(1,i,k); work(2,i) = a(2,i,k)
               a(1,i,k) = 0.0d0; a(2,i,k) = 0.d0
  110       continue
            do 120 j = kp1, n
               r(1) = work(1,j); r(2) = work(2,j)
               !call zaxpy(n,t,a(1,j),1,a(1,k),1)
               do i=1,n
                  a(1,i,k)=r(1)*a(1,i,j)-r(2)*a(2,i,j)+a(1,i,k)
                  a(2,i,k)=r(2)*a(1,i,j)+r(1)*a(2,i,j)+a(2,i,k)
               end do
  120       continue
            l = ipvt(k)
            if (l .ne. k) then
               !call zswap(n,a(1,k),1,a(1,l),1)
               do i=1,n
                  r(1) = a(1,i,k); r(2) = a(2,i,k)
                  a(1,i,k) = a(1,i,l); a(2,i,k) = a(2,i,l)
                  a(1,i,l) = r(1); a(2,i,l) = r(2)
               end do
            end if
  130    continue
  140    continue
  150 continue

end subroutine dzgedi
!!***

!!****f* m_hide_lapack/dzgefa
!! NAME
!!  dzgefa
!!
!! FUNCTION
!!   This routine is the clone of zgefa.F90 using real*8 a(2) instead of complex*16
!!   for the purpose of ABINIT (2008,TD)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_berryphase,m_berrytk,m_dfpt_fef,m_elpolariz,m_pead_nl_loop,m_relaxpol
!!
!! CHILDREN
!!
!! SOURCE

subroutine dzgefa(a,lda,n,ipvt,info)

 use m_linalg_interfaces

!Arguments
 integer :: lda,n,ipvt(n),info
 real*8  :: a(2,lda,n)
!
!     zgefa factors a complex*16 matrix by gaussian elimination.
!
!     dzgefa is usually called by zgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
!
!     on entry
!
!        a       complex*16(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that zgesl or zgedi will divide by zero
!                     if called.  use  rcond  in zgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     internal variables
!
!Local variables
 real*8 :: r(2),rk(2),rlj(2)
 real*8 :: rinv2,rmax,rabs
 integer :: i,j,k,kp1,l,nm1

!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         !l = izamax(n-k+1,a(k,k),1) + k - 1
         rmax=0.d0
         l=0
         do i=k,n
            rabs=abs(a(1,i,k))+abs(a(2,i,k))
            if(rmax<=rabs) then
              rmax=rabs
              l=i
            end if
         end do
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (abs(a(1,l,k))+abs(a(2,l,k)) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               r(1) = a(1,l,k); r(2) = a(2,l,k)
               a(1,l,k) = a(1,k,k); a(2,l,k) = a(2,k,k)
               a(1,k,k) = r(1); a(2,k,k) = r(2)
   10       continue
!
!           compute multipliers
!
            rinv2 = 1.d0/(a(1,k,k)**2+a(2,k,k)**2)
            rk(1) = -rinv2*a(1,k,k)
            rk(2) =  rinv2*a(2,k,k)
            !call zscal(n-k,t,a(k+1,k),1)
            do i=k+1,n
               r(1)=a(1,i,k)
               r(2)=a(2,i,k)
               a(1,i,k)=rk(1)*r(1)-rk(2)*r(2)
               a(2,i,k)=rk(1)*r(2)+rk(2)*r(1)
            end do
!
!           row elimination with column indexing
!
            do j = kp1, n
               rlj(1) = a(1,l,j); rlj(2) = a(2,l,j)
               if (l .eq. k) go to 20
                  a(1,l,j) = a(1,k,j); a(2,l,j) = a(2,k,j)
                  a(1,k,j) = rlj(1); a(2,k,j) = rlj(2)
   20          continue
               !call zaxpy(n-k,t,a(1,k+1,k),1,a(1,k+1,j),1)
               do i=k+1,n
                  a(1,i,j)=rlj(1)*a(1,i,k)-rlj(2)*a(2,i,k)+a(1,i,j)
                  a(2,i,j)=rlj(2)*a(1,i,k)+rlj(1)*a(2,i,k)+a(2,i,j)
               end do
            end do
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (abs(a(1,n,n))+abs(a(2,n,n)) .eq. 0.0d0) info = n

end subroutine dzgefa
!!***

!!****f* m_hide_lapack/test_xginv
!! NAME
!!  test_xginv
!!
!! FUNCTION

subroutine test_xginv(msize,skinds,do_check,Tres,comm)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: msize,comm
 logical,intent(in) :: do_check
 character(len=*),intent(in) :: skinds
 type(latime_t),intent(out) :: Tres
!arrays
 !complex(spc),allocatable :: cmat_spc(:,:)
 !complex(spc),allocatable :: cmat_spc_check(:,:)
 complex(dpc),allocatable :: cmat_dpc(:,:)
 complex(dpc),allocatable :: cmat_dpc_check(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: max_abserr

! *************************************************************************

 if (.FALSE.) write(std_out,*)skinds

 if (do_check) then
   ABI_MALLOC(cmat_dpc_check,(msize,msize))
   cmat_dpc_check = czero
   do ii=1,msize
    cmat_dpc_check(ii,ii) = cone
   end do
   !call xginv(cmat_dpc_check,msize,comm=xmpi_comm_self)
 end if

 ABI_MALLOC(cmat_dpc,(msize,msize))
 do ii=1,msize
  cmat_dpc(ii,ii) = cone
 end do

 call cwtime(Tres%ctime,Tres%wtime,Tres%gflops,"start")

 call xginv(cmat_dpc,msize,comm)

 call cwtime(Tres%ctime,Tres%wtime,Tres%gflops,"stop")
 Tres%testname  = 'test_xginv'
 Tres%msize     = msize

 max_abserr = -one
 if (do_check) then
   max_abserr = MAXVAL( ABS(cmat_dpc - cmat_dpc_check) )
 end if
 Tres%max_abserr = max_abserr

 ABI_FREE(cmat_dpc)

 ABI_SFREE(cmat_dpc_check)

end subroutine test_xginv
!!***

!!****f* m_cgtools/xhesv_cplex
!! NAME
!!   xhesv_cplex
!!
!! FUNCTION
!! ZHESV computes the solution to a complex (real) system of linear equations
!!    A * X = B,
!!
!! where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS matrices.
!! The value of cplex (1 or 2) defines whether we have a complex Hermitian or real symmetric matrix
!!
!! The diagonal pivoting method is used to factor A as
!!    A = U * D * U**H,  if UPLO = 'U', or
!!    A = L * D * L**H,  if UPLO = 'L',
!!
!! where U (or L) is a product of permutation and unit upper (lower)
!! triangular matrices, and D is Hermitian and block diagonal with
!! 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!! used to solve the system of equations A * X = B.
!!
!! INPUTS
!!
!![in]    UPLO
!!          UPLO is CHARACTER*1
!!          = 'U':  Upper triangle of A is stored;
!!          = 'L':  Lower triangle of A is stored.
!![in]    N
!!          N is INTEGER
!!          The number of linear equations, i.e., the order of the
!!          matrix A.  N >= 0.
!![in]    NRHS
!!          NRHS is INTEGER
!!          The number of right hand sides, i.e., the number of columns
!!          of the matrix B.  NRHS >= 0.
!![in,out]    A
!!          A is COMPLEX*16 array, dimension (LDA,N)
!!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!!          N-by-N upper triangular part of A contains the upper
!!          triangular part of the matrix A, and the strictly lower
!!          triangular part of A is not referenced.  If UPLO = 'L', the
!!          leading N-by-N lower triangular part of A contains the lower
!!          triangular part of the matrix A, and the strictly upper
!!          triangular part of A is not referenced.
!!
!!          On exit, if INFO = 0, the block diagonal matrix D and the
!!          multipliers used to obtain the factor U or L from the
!!          factorization A = U*D*U**H or A = L*D*L**H as computed by
!!          ZHETRF.
!![in,out]    B
!!          B is COMPLEX*16 array, dimension (LDB,NRHS)
!!          On entry, the N-by-NRHS right hand side matrix B.
!!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!![out]   INFO
!!          INFO is INTEGER
!!          = 0: successful exit
!!          < 0: if INFO = -i, the i-th argument had an illegal value
!!          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!!               has been completed, but the block diagonal matrix D is
!!               exactly singular, so the solution could not be computed.
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xhesv_cplex(UPLO, cplex, N, NRHS, A, B, msg, info)

!Arguments ------------------------------------
 character(len=1),intent(in) :: UPLO
 integer,intent(in) :: cplex, N, NRHS
 real(dp),intent(inout) :: A(cplex, N, N)
 real(dp),intent(inout) :: B(cplex, N, NRHS)
 character(len=*),intent(out) :: msg
 integer,intent(out) :: info

!Local variables ------------------------------
!scalars
 integer :: lwork, lda, ldb
 integer,allocatable :: ipiv(:)
 real(dp),allocatable :: work(:,:)

!************************************************************************

 if (all(cplex /= [1, 2])) then
   write(msg,'(a,i0)')" Wrong value for cplex: ",cplex
   info = 1; return
 end if

 lda = N; ldb = N

 ABI_MALLOC(ipiv, (N))
 ABI_MALLOC(work, (cplex, 1))
 lwork = -1

 if (cplex == 2) then
   ! Complex version
   call zhesv(uplo, N, NRHS, A, lda, ipiv, B, ldb, work, lwork, info)
   lwork = int(work(1, 1))
   ABI_REMALLOC(work, (cplex, lwork))

   call zhesv(uplo, N, NRHS, A, lda, ipiv, B, ldb, work, lwork, info)
 else
   ! Read version
   call dsysv(uplo, N, NRHS, A, lda, ipiv, B, ldb, work, lwork, info)
   lwork = int(work(1, 1))
   ABI_REMALLOC(work, (cplex, lwork))

   call dsysv(uplo, N, NRHS, A, lda, ipiv, B, ldb, work, lwork, info)
 end if

 ABI_FREE(ipiv)
 ABI_FREE(work)

 if (info < 0) then
   write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZPOTRI had an illegal value."
 else if (info > 0) then
   write(msg, "(a,i0,4a)") &
     " D(i,i) is exactly zero for i= ", info, ch10, &
     "The factorization has been completed, but the block diagonal matrix D is ", ch10, &
     "exactly singular, so the solution could not be computed."
 end if

end subroutine xhesv_cplex
!!***

end module m_hide_lapack
!!***
