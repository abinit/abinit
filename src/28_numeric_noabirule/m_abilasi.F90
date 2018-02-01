!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abilasi
!! NAME
!!  m_abilasi
!!
!! FUNCTION
!!  ABInit Linear Algebra Subroutine Interfaces.
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
!!  In F90 one can pass the array descriptor if the routines should operate on a slice of the local array (seldom done in abinit).
!!  Using array descriptor is OK but will it likely slow-down the calculation as some compilers perform a copy of the input-output data.
!!  If efficiency is a concern, then the F77 call should be used 
!!
!! COPYRIGHT
!! Copyright (C) 1992-2018 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN

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

MODULE m_abilasi

 use defs_basis
 use m_profiling_abi
 use m_xmpi
 use m_errors 
 use m_slk
 use m_linalg_interfaces

 use m_time,       only : cwtime
 use m_fstrings,   only : firstchar

 implicit none

 private

! Procedures for complex Hermitian matrices

 public :: xheev   ! Computes all the eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix.

 public :: xhpev   ! Computes all the eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix 
                   !   in packed storage (Scalapack version not available)

 public :: xhegv   ! Compute all the eigenvalues, and optionally, the eigenvectors of a complex generalized 
                   !   Hermitian-definite eigenproblem, of the form:
                   !   A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x

 public :: xheevx  ! Computes selected eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  
                   !   Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
                   !   indices for the desired eigenvalues.

 public :: xhegvx  ! Computes selected eigenvalues, and optionally, eigenvectors of a complex generalized Hermitian-definite eigenproblem, 
                   !   of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. 
                   !   Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
                   !   indices for the desired eigenvalues.

! Procedures for complex non-symmetric matrices

 public :: xgeev   ! Computes for a complex nonsymmetric matrix A, the eigenvalues and, optionally, 
                   ! the left and/or right eigenvectors.

 public :: xginv   ! Invert a general matrix of complex elements by means of LU factorization.


 public :: xhdp_invert   ! Invert a Hermitian positive definite matrix.

 interface xheev
  module procedure wrap_CHEEV
  module procedure wrap_ZHEEV
  module procedure wrap_DSYEV_ZHEEV
 end interface xheev

 interface xhpev
  module procedure wrap_CHPEV
  module procedure wrap_ZHPEV
 end interface xhpev

 interface xhegv
  module procedure wrap_ZHEGV
  module procedure wrap_DSYGV_ZHEGV
 end interface xhegv

 interface xheevx
  module procedure wrap_ZHEEVX 
  module procedure wrap_DSYEVX_ZHEEVX
 end interface xheevx

 interface xhegvx
  module procedure wrap_ZHEGVX 
  module procedure wrap_DSYGVX_ZHEGVX
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

!!****f* m_abilasi/wrap_CHEEV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_CHEEV(jobz,uplo,n,a,w)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_CHEEV'
!End of the abilint section

 implicit none

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
&  "CHEEV: the algorithm failed to converge; ",ch10,&
&  info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
  MSG_ERROR(msg)
 end if

 ABI_FREE(rwork)
 ABI_FREE(work)

 !TODO scaLAPACK version (complex single precision buffer is needed in matrix_scalapack)

end subroutine wrap_CHEEV
!!***

!----------------------------------------------------------------------

!!****f* m_abilasi/wrap_ZHEEV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_ZHEEV(jobz,uplo,n,a,w,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_ZHEEV'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
&   "ZHEEV: the algorithm failed to converge; ",ch10,&
&   info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
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

!!****f* m_abilasi/wrap_DSYEV_ZHEEV
!! NAME
!!  wrap_DSYEV_ZHEEV
!!
!! FUNCTION
!!  wrap_DSYEV_ZHEEV computes the eigenvalues and, optionally, the eigenvectors of a 
!!  (complex Hermitian| real symmetric) matrix in double precision. [PRIVATE]
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_DSYEV_ZHEEV(jobz,uplo,cplex,n,a,w,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_DSYEV_ZHEEV'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,cplex
 integer,optional,intent(in) :: comm
 character(len=*),intent(in) :: jobz,uplo 
!arrays
 real(dp),intent(inout) :: a(cplex,n,n)
 real(dp),intent(out) :: w(n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
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

 if (ALL(cplex/=(/1,2/))) then 
  write(msg,'(a,i0)')" Wrong value for cplex: ",cplex
  MSG_BUG(msg)
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.)

  if (cplex==1) then  ! Real symmetric case.

   lwork = MAX(1,3*n-1)

   ABI_MALLOC(work_real,(lwork))

   call DSYEV(jobz,uplo,n,a,n,w,work_real,lwork,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of DSYEV had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(2a,i0,a)')&
&    "DSYEV: the algorithm failed to converge; ",ch10,&
&    info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    MSG_ERROR(msg)
   end if

   ABI_FREE(work_real)

   RETURN

  else                ! Hermitian case.

   lwork = MAX(1,2*n-1)

   ABI_MALLOC(work_cplx, (lwork))
   ABI_MALLOC(rwork, (MAX(1,3*n-2)))

   call ZHEEV(jobz,uplo,n,a,n,w,work_cplx,lwork,rwork,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEEV had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(2a,i0,a)')&
&    "ZHEEV: the algorithm failed to converge; ",ch10,&
&    info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    MSG_ERROR(msg)
   end if

   ABI_FREE(rwork)
   ABI_FREE(work_cplx)

   RETURN
  end if ! cplex

 CASE (.TRUE.)

#ifdef HAVE_LINALG_SCALAPACK

  MSG_ERROR("Not coded yet")

!  call init_scalapack(Slk_processor,comm)
!  istwf_k=1
!
!  ! Initialize and fill Scalapack matrix from the global one.
!  tbloc=SLK_BLOCK_SIZE
!  call init_matrix_scalapack(Slk_mat,n,n,Slk_processor,istwf_k,tbloc=tbloc)
!
!  write(msg,'(2(a,i0))')" Using scaLAPACK version with nprocs = ",nprocs,"; block size = ",tbloc
!  call wrtout(std_out,msg,"PERS")
!
!  call slk_matrix_from_global_dpc_2D(Slk_mat,uplo,a)
!
!  want_eigenvectors = firstchar(jobz,(/"V","v"/))
!  if (want_eigenvectors) then ! Initialize the distributed vectors.
!   call init_matrix_scalapack(Slk_vec,n,n,Slk_processor,istwf_k,tbloc=tbloc) 
!  end if
!
!  ! Solve the problem with scaLAPACK.
!  call slk_pzheev(jobz,uplo,Slk_mat,Slk_vec,w)
!
!  call destruction_matrix_scalapack(Slk_mat)
!  
!  if (want_eigenvectors) then ! A is overwritten with the eigenvectors 
!   a = czero  
!   call slk_matrix_to_global_dpc_2D(Slk_vec,"All",a) ! Fill the entries calculated by this node.
!   call destruction_matrix_scalapack(Slk_vec)
!   call xmpi_sum(a,comm,ierr)                        ! Fill the remaing entries of the global matrix
!  end if
!
!  call end_scalapack(Slk_processor)

  RETURN 
#endif

  MSG_BUG("You should not be here!")

 END SELECT

end subroutine wrap_DSYEV_ZHEEV
!!***

!----------------------------------------------------------------------

!!****f* m_abilasi/wrap_CHPEV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_CHPEV(jobz,uplo,n,ap,w,z,ldz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_CHPEV'
!End of the abilint section

 implicit none

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
&  "ZHPEV: the algorithm failed to converge; ",ch10,&
&  info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
  MSG_ERROR(msg)
 end if

 ABI_FREE(rwork)
 ABI_FREE(work)

end subroutine wrap_CHPEV
!!***

!----------------------------------------------------------------------

!!****f* m_abilasi/wrap_ZHPEV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_ZHPEV(jobz,uplo,n,ap,w,z,ldz,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_ZHPEV'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
&   "ZHPEV: the algorithm failed to converge; ",ch10,&
&   info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
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


!!****f* m_abilasi/wrap_ZHEGV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_ZHEGV(itype,jobz,uplo,n,a,b,w,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_ZHEGV'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
&    "ZHEGV failed to converge: ",ch10,&
&    info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
   else
    ii = info -n 
    write(msg,'(3a,i0,3a)')&
&   "ZHEGV failed to converge: ",ch10,&
&   "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
&   "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
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

!!****f* m_abilasi/wrap_DSYGV_ZHEGV
!! NAME
!!  wrap_DSYGV_ZHEGV
!!
!! FUNCTION
!!  wrap_DSYGV_ZHEGV computes all the  eigenvalues, and  optionally, the eigenvectors of a  
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_DSYGV_ZHEGV(itype,jobz,uplo,cplex,n,a,b,w,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_DSYGV_ZHEGV'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n,itype,cplex
 integer,optional,intent(in) :: comm
 character(len=*),intent(in) :: jobz,uplo 
!arrays
 real(dp),intent(inout) :: a(cplex,n,n),b(cplex,n,n)
 real(dp),intent(out) :: w(n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs,ii
 logical :: use_scalapack
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rwork(:)
 real(dp),allocatable :: work_real(:)
 complex(dpc),allocatable :: work_cplx(:)
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

 if (ALL(cplex/=(/1,2/))) then 
  write(msg,'(a,i0)')"Wrong value for cplex: ",cplex
  MSG_ERROR(msg)
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.)

  if (cplex==1) then ! Real symmetric case.

   lwork = MAX(1,3*n-1)

   ABI_MALLOC(work_real,(lwork))

   call DSYGV(itype,jobz,uplo,n,a,n,b,n,w,work_real,lwork,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of DSYGV had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    if (info<= n) then 
     write(msg,'(2a,i0,a)')&
&     " DSYGV failed to converge: ",ch10,&
&     info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    else
     ii = info -n 
     write(msg,'(3a,i0,3a)')&
&    "DSYGV failed to converge: ",ch10,&
&    "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
&    "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
    end if
    MSG_ERROR(msg)
   end if

   ABI_FREE(work_real)

   RETURN

  else               ! complex Hermitian case

   lwork = MAX(1,2*n-1)

   ABI_MALLOC(work_cplx,(lwork))
   ABI_MALLOC(rwork,(MAX(1,3*n-2)))
 
   call ZHEGV(itype,jobz,uplo,n,a,n,b,n,w,work_cplx,lwork,rwork,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEGV had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    if (info<= n) then 
     write(msg,'(2a,i0,a)')&
&     "ZHEEV failed to converge: ",ch10,&
&     info," off-diagonal elements of an intermediate tridiagonal form did not converge to zero. "
    else
     ii = info -n 
     write(msg,'(3a,i0,3a)')&
&    "ZHEEV failed to converge: ",ch10,&
&    "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
&    "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
    end if
    MSG_ERROR(msg)
   end if

   ABI_FREE(rwork)
   ABI_FREE(work_cplx)

   RETURN

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

end subroutine wrap_DSYGV_ZHEGV
!!***

!----------------------------------------------------------------------


!!****f* m_abilasi/wrap_ZHEEVX
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_ZHEEVX(jobz,range,uplo,n,a,vl,vu,il,iu,abstol,m,w,z,ldz,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_ZHEEVX'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
   write(msg,'(2a,i0,a)')&
&   "ZHEEVX: the algorithm failed to converge; ",ch10,&
&   info,"eigenvectors failed to converge. "
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


!!****f* m_abilasi/wrap_DSYEVX_ZHEEVX
!! NAME
!!  wrap_DSYEVX_ZHEEVX
!!
!! FUNCTION
!!  wrap_DSYEVX_ZHEEVX computes selected eigenvalues and, optionally, eigenvectors
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_DSYEVX_ZHEEVX(jobz,range,uplo,cplex,n,a,vl,vu,il,iu,abstol,m,w,z,ldz,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_DSYEVX_ZHEEVX'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,iu,ldz,n,cplex
 integer,optional,intent(in) :: comm
 integer,intent(inout) :: m
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo 
!arrays
 real(dp),intent(out) :: w(n)
 !real(dp),intent(out) :: z(cplex,ldz,n)
 real(dp),intent(out) :: z(cplex,ldz,m)
 real(dp),intent(inout) :: a(cplex,n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs
 logical :: use_scalapack
 character(len=500) :: msg
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
  MSG_ERROR(msg)
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.) ! Standard LAPACK call.

  if (cplex==1) then      ! Real symmetric case

   lwork = MAX(1,8*n)

   ABI_MALLOC(work_real,(lwork))
   ABI_MALLOC(iwork,(5*n))
   ABI_MALLOC(ifail,(n))

   call DSYEVX(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_real,lwork,iwork,ifail,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of DSYEVX had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(2a,i0,a)')&
&    "DSYEVX: the algorithm failed to converge; ",ch10,&
&    info,"eigenvectors failed to converge. "
    MSG_ERROR(msg)
   end if

   ABI_FREE(work_real)
   ABI_FREE(iwork)
   ABI_FREE(ifail)

   RETURN

  else                    ! Complex Hermitian case.

   lwork = MAX(1,2*n)

   ABI_MALLOC(work_cplx,(lwork))
   ABI_MALLOC(rwork,(7*n))
   ABI_MALLOC(iwork,(5*n))
   ABI_MALLOC(ifail,(n))

   call ZHEEVX(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_cplx,lwork,rwork,iwork,ifail,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZHEEVX had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    write(msg,'(2a,i0,a)')&
&    "ZHEEVX: the algorithm failed to converge; ",ch10,&
&    info,"eigenvectors failed to converge. "
    MSG_ERROR(msg)
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

end subroutine wrap_DSYEVX_ZHEEVX
!!***

!----------------------------------------------------------------------

!!****f* m_abilasi/wrap_ZHEGVX
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_ZHEGVX(itype,jobz,range,uplo,n,a,b,vl,vu,il,iu,abstol,m,w,z,ldz,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_ZHEGVX'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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

 CASE (.FALSE.) ! Standard LAPACK call.

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
    write(msg,'(a,i0,a)')&
&    "ZHEGVX failed to converge: ",info," eigenvectors failed to converge. "
   else
    ii = info -n 
    write(msg,'(3a,i0,3a)')&
&   "ZHEEVX failed to converge: ",ch10,&
&   "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
&   "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
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

!!****f* m_abilasi/wrap_DSYGVX_ZHEGVX
!! NAME
!!  wrap_DSYGVX_ZHEGVX
!!
!! FUNCTION
!!  wrap_DSYGVX_ZHEGVX  - compute selected eigenvalues, and optionally, eigenvectors of a 
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_DSYGVX_ZHEGVX(itype,jobz,range,uplo,cplex,n,a,b,vl,vu,il,iu,abstol,m,w,z,ldz,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_DSYGVX_ZHEGVX'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,iu,ldz,n,itype,cplex
 integer,optional,intent(in) :: comm
 integer,intent(inout) :: m
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo 
!arrays
 real(dp),intent(out) :: w(n)
 !real(dp),intent(out) :: z(cplex,ldz,n)
 real(dp),intent(out) :: z(cplex,ldz,m)
 real(dp),intent(inout) :: a(cplex,n,n),b(cplex,n,n)

!Local variables ------------------------------
!scalars
 integer :: lwork,info,nprocs,ii
 logical :: use_scalapack
 character(len=500) :: msg
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
  MSG_ERROR(msg)
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.) ! Standard LAPACK call.

  if (cplex==1) then  ! Real symmetric case 

   lwork = MAX(1,8*n)

   ABI_MALLOC(work_real,(lwork))
   ABI_MALLOC(iwork,(5*n))
   ABI_MALLOC(ifail,(n))

   call DSYGVX(itype,jobz,range,uplo,n,a,n,b,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_real,lwork,iwork,ifail,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')" The ",-info,"-th argument of DSYGVX had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    if (info<= n) then 
     write(msg,'(a,i0,a)')&
&     " DSYGVX failed to converge: ",info," eigenvectors failed to converge. "
    else
     ii = info -n 
     write(msg,'(3a,i0,3a)')&
&    " DSYGVX failed to converge: ",ch10,&
&    " The leading minor of order ",ii," of B is not positive definite. ",ch10,&
&    " The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
    end if
    MSG_ERROR(msg)
   end if

   ABI_FREE(iwork)
   ABI_FREE(ifail)
   ABI_FREE(work_real)

   RETURN

  else                ! Complex Hermitian case.

   lwork = MAX(1,2*n)

   ABI_MALLOC(work_cplx,(lwork))
   ABI_MALLOC(rwork,(7*n))
   ABI_MALLOC(iwork,(5*n))
   ABI_MALLOC(ifail,(n))

   !write(std_out,*)"Calling ZHEGVX"

   call ZHEGVX(itype,jobz,range,uplo,n,a,n,b,n,vl,vu,il,iu,abstol,m,w,z,ldz,work_cplx,lwork,rwork,iwork,ifail,info)

   if (info < 0) then 
    write(msg,'(a,i0,a)')"The ",-info,"-th argument of ZHEGVX had an illegal value."
    MSG_ERROR(msg)
   end if

   if (info > 0) then
    if (info<= n) then 
     write(msg,'(a,i0,a)')&
&     "ZHEGVX failed to converge: ",info," eigenvectors failed to converge. "
    else
     ii = info -n 
     write(msg,'(3a,i0,3a)')&
&    "ZHEEVX failed to converge: ",ch10,&
&    "The leading minor of order ",ii," of B is not positive definite. ",ch10,&
&    "The factorization of B could not be completed and no eigenvalues or eigenvectors were computed."
    end if
    MSG_ERROR(msg)
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

end subroutine wrap_DSYGVX_ZHEGVX
!!***

!----------------------------------------------------------------------

!!****f* m_abilasi/wrap_CGEEV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_CGEEV(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_CGEEV'
!End of the abilint section

 implicit none

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
&  "CGEEV: The QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed;",ch10,&
&  "Elements ",info+1,":",n," of W contain eigenvalues which have converged. "
  MSG_ERROR(msg)
 end if
      
 ABI_FREE(work)
 ABI_FREE(rwork)

end subroutine wrap_CGEEV
!!***

!----------------------------------------------------------------------

!!****f* m_abilasi/wrap_ZGEEV
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine wrap_ZGEEV(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap_ZGEEV'
!End of the abilint section

 implicit none

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
&   "ZGEEV: The QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed;",ch10,&
&   "Elements ",info+1,":",n," of W contain eigenvalues which have converged. "
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

!!****f* m_abilasi/cginv
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine cginv(a,n,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cginv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
#ifdef HAVE_LINALG_SCALAPACK
  use_scalapack = (nprocs>1)
#endif
 end if

 SELECT CASE(use_scalapack)

 CASE (.FALSE.)

  ABI_MALLOC(ipiv,(n))

  call CGETRF(n,n,a,n,ipiv,info) ! P* L* U  Factorization.

  if (info < 0) then 
   write(msg,'(a,i0,a)')" The ",-info,"-th argument of CGETRF had an illegal value."
   MSG_ERROR(msg)
  end if

  if (info > 0) then
   write(msg,'(3a,i0,4a)')&
&   "The matrix that has been passed in argument is probably either singular or nearly singular.",ch10,&
&   "U(i,i) in the P*L*U factorization is exactly zero for i = ",info,ch10,&
&   "The factorization has been completed but the factor U is exactly singular.",ch10,&
&   "Division by zero will occur if it is used to solve a system of equations."
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
&   "The matrix that has been passed to this subroutine is probably either singular or nearly singular.",ch10,&
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

!!****f* m_abilasi/zginv
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine zginv(a,n,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zginv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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

  ABI_MALLOC(ipiv,(n))
  call ZGETRF(n,n,a,n,ipiv,info) ! P* L* U  Factorization.

  if (info < 0) then 
   write(msg,'(a,i0,a)')" The ",-info,"-th argument of ZGETRF had an illegal value."
   MSG_ERROR(msg)
  end if

  if (info > 0) then
   write(msg,'(3a,i0,4a)')&
&   "The matrix that has been passed in argument is probably either singular or nearly singular.",ch10,&
&   "U(i,i) in the P*L*U factorization is exactly zero for i = ",info,ch10,&
&   "The factorization has been completed but the factor U is exactly singular.",ch10,&
&   "Division by zero will occur if it is used to solve a system of equations."
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
&   "The matrix that has been passed to this subroutine is probably either singular or nearly singular.",ch10,&
&   "U(i,i) for i= ",info," is exactly zero; the matrix is singular and its inverse could not be computed."
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

!!****f* m_abilasi/zhpd_invert
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
!!      cwtime,xginv
!!
!! SOURCE

subroutine zhpd_invert(uplo,a,n,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zhpd_invert'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
&     "The leading minor of order ",info," is not positive definite, ",ch10,&
&     "and the factorization could not be completed."
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
&     "The ( ",info,info,")element of the factor U or L is zero, and the inverse could not be computed."
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


!!****f* m_abilasi/test_xginv
!! NAME
!!  test_xginv
!!
!! FUNCTION

subroutine test_xginv(msize,skinds,do_check,Tres,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'test_xginv'
!End of the abilint section

 implicit none

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
 Tres%testname  = ABI_FUNC
 Tres%msize     = msize

 max_abserr = -one
 if (do_check) then
   max_abserr = MAXVAL( ABS(cmat_dpc - cmat_dpc_check) )
 end if
 Tres%max_abserr = max_abserr

 ABI_FREE(cmat_dpc)

 if (allocated(cmat_dpc_check)) then
   ABI_FREE(cmat_dpc_check)
 end if

end subroutine test_xginv
!!***

END MODULE m_abilasi
!!***
