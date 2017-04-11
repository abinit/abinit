!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_elpa
!! NAME
!! m_elpa
!!
!! FUNCTION
!! This module contains interfaces to ELPA library methods
!! See http://elpa.mpcdf.mpg.de
!!
!! COPYRIGHT
!! Copyright (C) 2016-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_elpa

 use defs_basis
 use m_errors
#ifdef HAVE_LINALG_ELPA
 use elpa1
#endif

 implicit none

 private

!Public procedures
#ifdef HAVE_LINALG_ELPA
 public :: elpa_get_elpa_communicators
 public :: elpa_solve_evp_real_1stage
 public :: elpa_solve_evp_complex_1stage
 public :: elpa_cholesky_real
 public :: elpa_cholesky_complex
 public :: elpa_invert_trm_real
 public :: elpa_invert_trm_complex
 public :: elpa_mult_at_b_real
 public :: elpa_mult_ah_b_complex
#endif

CONTAINS  !==============================================================================
!!***

#ifdef HAVE_LINALG_ELPA

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_get_elpa_communicators
!! NAME
!!  elpa_get_elpa_communicators
!!
!! FUNCTION
!!  Wrapper to get_elpa_communicators ELPA function
!!
!! INPUTS
!!  mpi_comm_global=Global communicator for the calculations (in)
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! OUTPUT
!!  mpi_comm_rows=Communicator for communicating within rows of processes (out)
!!  mpi_comm_cols=Communicator for communicating within columns of processes (out)
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_get_elpa_communicators(mpi_comm_global,my_prow,my_pcol,mpi_comm_rows,mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_get_elpa_communicators'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: mpi_comm_global, my_prow, my_pcol
 integer,intent(out) :: mpi_comm_rows, mpi_comm_cols

!Local variables-------------------------------
 integer  :: mpierr

! *********************************************************************

 mpierr=0
 
#if (defined HAVE_LINALG_ELPA_2013) || (defined HAVE_LINALG_ELPA_2014)
 call get_elpa_row_col_comms(mpi_comm_global,my_prow,my_pcol,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015)
 mpierr=get_elpa_row_col_comms(mpi_comm_global,my_prow,my_pcol,mpi_comm_rows,mpi_comm_cols)
#else
 mpierr=get_elpa_communicators(mpi_comm_global,my_prow,my_pcol,mpi_comm_rows,mpi_comm_cols)
#endif

 ABI_CHECK_MPI(mpierr,"elpa_get_elpa_communicators")

end subroutine elpa_get_elpa_communicators
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_solve_evp_real_1stage
!! NAME
!!  elpa_solve_evp_real_1stage
!!
!! FUNCTION
!!  Wrapper to solve_evp_real_1stage ELPA function
!!
!! INPUTS
!!  na=Order of matrix aa
!!  nev=Number of eigenvalues needed.
!!  lda=Leading dimension of aa
!!  ldq=Leading dimension of qq
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  matrixCols=Distributed number of matrix columns
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(ldq,matrixCols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na,na))
!!                      even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!! aa(lda,matrixCols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_solve_evp_real_1stage(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols, &
&                                     mpi_comm_rows, mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_solve_evp_real_1stage'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,nev,lda,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!arrays
 real(dp),intent(inout) :: aa(lda,matrixCols)
 real(dp),intent(out) :: ev(na),qq(ldq,matrixCols)

!Local variables-------------------------------
 logical  :: success

! *********************************************************************

 success=.true.

#if (defined HAVE_LINALG_ELPA_2013)
  call solve_evp_real(na,nev,aa,lda,ev,qq,ldq,nblk,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2015)
  success=solve_evp_real(na,nev,aa,lda,ev,qq,ldq,nblk,mpi_comm_rows,mpi_comm_cols)
#else
  success=solve_evp_real_1stage(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (solve_evp_real)!')
 end if

end subroutine elpa_solve_evp_real_1stage
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_solve_evp_complex_1stage
!! NAME
!!  elpa_solve_evp_complex_1stage
!!
!! FUNCTION
!!  Wrapper to solve_evp_complex_1stage ELPA function
!!
!! INPUTS
!!  na=Order of matrix aa
!!  nev=Number of eigenvalues needed.
!!  lda=Leading dimension of aa
!!  ldq=Leading dimension of qq
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  matrixCols=Distributed number of matrix columns
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(ldq,matrixCols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na,na))
!!                      even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!! aa(lda,matrixCols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_solve_evp_complex_1stage(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols, &
&                                     mpi_comm_rows, mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_solve_evp_complex_1stage'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,nev,lda,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!arrays
 complex(dpc),intent(inout) :: aa(lda,matrixCols)
 real(dp),intent(out) :: ev(na)
 complex(dpc),intent(out) :: qq(ldq,matrixCols)

!Local variables-------------------------------
 logical  :: success

! *********************************************************************

 success=.true.

#if (defined HAVE_LINALG_ELPA_2013)
  call solve_evp_complex(na,nev,aa,lda,ev,qq,ldq,nblk,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2015)
  success=solve_evp_complex(na,nev,aa,lda,ev,qq,ldq,nblk,mpi_comm_rows,mpi_comm_cols)
#else
  success=solve_evp_complex_1stage(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (solve_evp_complex)!')
 end if

end subroutine elpa_solve_evp_complex_1stage
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_cholesky_real
!! NAME
!!  elpa_cholesky_real
!!
!! FUNCTION
!!  Wrapper to cholesky_real ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     On return, the upper triangle contains the Cholesky factor
!!                     and the lower triangle is set to 0.
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_cholesky_real(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_cholesky_real'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!arrays
 real(dp),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 logical :: success

! *********************************************************************

 success=.true.

#if (defined HAVE_LINALG_ELPA_2013)
 call cholesky_real(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014)
 call cholesky_real(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2015)
 call cholesky_real(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,.false.,success)
#else
 call cholesky_real(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols,.false.,success)
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (cholesky_real)!')
 end if

end subroutine elpa_cholesky_real
!!***

!----------------------------------------------------------------------
!!****f* m_elpa/elpa_cholesky_complex
!! NAME
!!  elpa_cholesky_complex
!!
!! FUNCTION
!!  Wrapper to cholesky_complex ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     On return, the upper triangle contains the Cholesky factor
!!                     and the lower triangle is set to 0.
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_cholesky_complex(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_cholesky_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!arrays
 complex(dpc),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 logical :: success

! *********************************************************************

 success=.true.

#if (defined HAVE_LINALG_ELPA_2013)
 call cholesky_complex(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014)
 call cholesky_complex(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2015)
 call cholesky_complex(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,.false.,success)
#else
 call cholesky_complex(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols,.false.,success)
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (cholesky_complex)!')
 end if

end subroutine elpa_cholesky_complex
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_invert_trm_real
!! NAME
!!  elpa_invert_trm_real
!!
!! FUNCTION
!!  Wrapper to invert_trm_real ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     The lower triangle is not referenced.
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_invert_trm_real(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_invert_trm_real'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!arrays
 real(dp),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 logical :: success

! *********************************************************************

 success=.true.

#if (defined HAVE_LINALG_ELPA_2013)
 call invert_trm_real(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014)
 call invert_trm_real(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2015)
 call invert_trm_real(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,.false.,success)
#else
 call invert_trm_real(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols,.false.,success)
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (invert_trm_real)!')
 end if

end subroutine elpa_invert_trm_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_invert_trm_complex
!! NAME
!!  elpa_invert_trm_complex
!!
!! FUNCTION
!!  Wrapper to invert_trm_complex ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     The lower triangle is not referenced.
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_invert_trm_complex(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_invert_trm_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!arrays
 complex(dpc),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 logical :: success

! *********************************************************************

 success=.true.

#if (defined HAVE_LINALG_ELPA_2013)
 call invert_trm_complex(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014)
 call invert_trm_complex(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2015)
 call invert_trm_complex(na,aa,lda,nblk,mpi_comm_rows,mpi_comm_cols,.false.,success)
#else
 call invert_trm_complex(na,aa,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols,.false.,success)
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (invert_trm_complex)!')
 end if

end subroutine elpa_invert_trm_complex
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_mult_at_b_real
!! NAME
!!  elpa_mult_at_b_real
!!
!! FUNCTION
!!  Wrapper to mult_at_b_real ELPA function
!!  Performs C := A**T * B
!!
!! INPUTS
!! uplo_a='U' if A is upper triangular
!!        'L' if A is lower triangular
!!        anything else if A is a full matrix
!!        Please note: This pertains to the original A (as set in the calling program)
!!          whereas the transpose of A is used for calculations
!!          If uplo_a is 'U' or 'L', the other triangle is not used at all,
!!          i.e. it may contain arbitrary numbers
!!  uplo_c='U' if only the upper diagonal part of C is needed
!!         'L' if only the upper diagonal part of C is needed
!!         anything else if the full matrix C is needed
!!         Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!!         written to a certain extent, i.e. one shouldn't rely on the content there!
!!  aa=Matrix A
!!  bb=Matrix B
!!  lda=Leading dimension of aa
!!  ldb=Leading dimension of bb
!!  ldc=Leading dimension of cc
!!  na=Number of rows/columns of A, number of rows of B and C
!!  ncb=Number of columns  of B and C
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  matrixCols=local columns of matrix a
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! OUTPUT
!!  cc=Matrix C
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_mult_at_b_real(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,nblk,matrixCols, &
&                              mpi_comm_rows,mpi_comm_cols,cc,ldc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_mult_at_b_real'
!End of the abilint section

 implicit none
 
!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,ncb,lda,ldb,ldc,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
 character*1 :: uplo_a, uplo_c
!arrays
 real(dp),intent(in) :: aa(lda,matrixCols),bb(ldb,matrixCols)
 real(dp),intent(out) :: cc(ldc,matrixCols)

!Local variables-------------------------------
 logical :: success

! *********************************************************************

 call mult_at_b_real(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,nblk, &
&                    mpi_comm_rows,mpi_comm_cols,cc,ldc)

end subroutine elpa_mult_at_b_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_mult_ah_b_complex
!! NAME
!!  elpa_mult_ah_b_complex
!!
!! FUNCTION
!!  Wrapper to mult_ah_b_complex ELPA function
!!  Performs C := A**H * B
!!
!! INPUTS
!! uplo_a='U' if A is upper triangular
!!        'L' if A is lower triangular
!!        anything else if A is a full matrix
!!        Please note: This pertains to the original A (as set in the calling program)
!!          whereas the transpose of A is used for calculations
!!          If uplo_a is 'U' or 'L', the other triangle is not used at all,
!!          i.e. it may contain arbitrary numbers
!!  uplo_c='U' if only the upper diagonal part of C is needed
!!         'L' if only the upper diagonal part of C is needed
!!         anything else if the full matrix C is needed
!!         Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!!         written to a certain extent, i.e. one shouldn't rely on the content there!
!!  aa=Matrix A
!!  bb=Matrix B
!!  lda=Leading dimension of aa
!!  ldb=Leading dimension of bb
!!  ldc=Leading dimension of cc
!!  na=Number of rows/columns of A, number of rows of B and C
!!  ncb=Number of columns  of B and C
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  matrixCols=local columns of matrix a
!!  mpi_comm_rows=MPI-Communicator for rows
!!  mpi_comm_cols=MPI-Communicator for columns
!!
!! OUTPUT
!!  cc=Matrix C
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_mult_ah_b_complex(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,nblk,matrixCols, &
&                              mpi_comm_rows,mpi_comm_cols,cc,ldc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_mult_ah_b_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,ncb,lda,ldb,ldc,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
 character*1 :: uplo_a, uplo_c
!arrays
 complex(dpc),intent(in) :: aa(lda,matrixCols),bb(ldb,matrixCols)
 complex(dpc),intent(out) :: cc(ldc,matrixCols)

!Local variables-------------------------------
 logical :: success

! *********************************************************************

 call mult_ah_b_complex(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,nblk, &
&                    mpi_comm_rows,mpi_comm_cols,cc,ldc)

end subroutine elpa_mult_ah_b_complex
!!***

!----------------------------------------------------------------------
#endif

END MODULE m_elpa
!!***
