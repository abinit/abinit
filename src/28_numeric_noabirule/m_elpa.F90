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

#if (defined HAVE_LINALG_ELPA) && (defined HAVE_LINALG_ELPA_2017)
#define HAVE_ELPA_FORTRAN2008
#else
#undef HAVE_ELPA_FORTRAN2008
#endif

module m_elpa

 use defs_basis
 use m_errors

#ifdef HAVE_LINALG_ELPA
#ifdef HAVE_ELPA_FORTRAN2008
 use elpa
#else
 use elpa1
#endif
#endif

 implicit none

 private

#ifdef HAVE_LINALG_ELPA

!Public procedures
!Had to choose names different from those provided by elpa
 public :: elpa_func_init
 public :: elpa_func_uninit
 public :: elpa_func_get_communicators
 public :: elpa_func_solve_evp_1stage
 public :: elpa_func_cholesky
 public :: elpa_func_invert_triangular
 public :: elpa_func_hermitian_multiply

 interface elpa_func_solve_evp_1stage
   module procedure elpa_func_solve_evp_1stage_real
   module procedure elpa_func_solve_evp_1stage_complex
 end interface elpa_func_solve_evp_1stage

 interface elpa_func_cholesky
   module procedure elpa_func_cholesky_real
   module procedure elpa_func_cholesky_complex
 end interface elpa_func_cholesky

 interface elpa_func_invert_triangular
   module procedure elpa_func_invert_triangular_real
   module procedure elpa_func_invert_triangular_complex
 end interface elpa_func_invert_triangular

 interface elpa_func_hermitian_multiply
   module procedure elpa_func_hermitian_multiply_real
   module procedure elpa_func_hermitian_multiply_complex
 end interface elpa_func_hermitian_multiply

#ifdef HAVE_ELPA_FORTRAN2008
!Handle for ELPA type
 class(elpa_t),pointer,private :: elpa_hdl
#else
!MPI-Communicator for rows
 integer,private,save :: elpa_comm_rows
!MPI-Communicator for columns
 integer,private,save :: elpa_comm_cols
#endif

#endif

CONTAINS  !==============================================================================
!!***

#ifdef HAVE_LINALG_ELPA

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_init
!! NAME
!!  elpa_func_init
!!
!! FUNCTION
!!  Wrapper to elpa_init ELPA function
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abi_linalg
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_init()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
#if (defined HAVE_LINALG_ELPA_2017)
 integer,parameter :: elpa_min_version=20170403
#endif
 logical :: success
!
! *********************************************************************

 success=.true.

#if   (defined HAVE_LINALG_ELPA_2017)
 success=(elpa_init(elpa_min_version)==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016) || (defined HAVE_LINALG_ELPA_2015)
!No init function
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
!No init function
#else
 success=.true.
#endif

 if (.not.success) then
   MSG_ERROR('Problem with ELPA (elpa_init)!')
 end if

end subroutine elpa_func_init
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_uninit
!! NAME
!!  elpa_func_uninit
!!
!! FUNCTION
!!  Wrapper to elpa_uninit ELPA function
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abi_linalg
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_uninit()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_uninit'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

! *********************************************************************

#if   (defined HAVE_LINALG_ELPA_2017)
 call elpa_uninit()
#elif (defined HAVE_LINALG_ELPA_2016) || (defined HAVE_LINALG_ELPA_2015)
!No uninit function
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
!No uninit function
#else
!No uninit function
#endif

end subroutine elpa_func_uninit
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_get_communicators
!! NAME
!!  elpa_func_get_communicators
!!
!! FUNCTION
!!  Wrapper to elpa_get_communicators ELPA function
!!
!! INPUTS
!!  mpi_comm_elpa=Global communicator for the calculations (in)
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! OUTPUT
!!  No output; store private variables elpa_comm_rows and elpa_comm_cols
!!
!! PARENTS
!!      m_slk
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_get_communicators(mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_get_communicators'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol

!Local variables-------------------------------
 integer  :: mpierr

! *********************************************************************

 mpierr=0

#if (defined HAVE_LINALG_ELPA_2017)
!This function doesnt exist anymore in the F2008 interface
#elif (defined HAVE_LINALG_ELPA_2016)
 mpierr=elpa_get_communicators(mpi_comm_elpa,my_prow,my_pcol,elpa_comm_rows,elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015)
 mpierr=get_elpa_row_col_comms(mpi_comm_elpa,my_prow,my_pcol,elpa_comm_rows,elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
 call get_elpa_row_col_comms(mpi_comm_elpa,my_prow,my_pcol,elpa_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 mpierr=get_elpa_communicators(mpi_comm_elpa,my_prow,my_pcol,elpa_comm_rows,elpa_comm_cols)
#endif

 if (mpierr/=0) then
   MSG_ERROR('Problem with ELPA (elpa_get_communicators)!')
 end if

end subroutine elpa_func_get_communicators
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_solve_evp_1stage_real
!! NAME
!!  elpa_func_solve_evp_1stage_real
!!
!! FUNCTION
!!  Wrapper to elpa_solve_evp_real_1stage ELPA function
!!
!! INPUTS
!!  na=Order of matrix aa
!!  nev=Number of eigenvalues needed.
!!  lda=Leading dimension of aa
!!  ldq=Leading dimension of qq
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  matrixCols=Distributed number of matrix columns
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
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
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_solve_evp_1stage_real(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,&
&                                          mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_solve_evp_1stage_real'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,nev,lda,ldq,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
!arrays
 real(dp),intent(inout) :: aa(lda,matrixCols)
 real(dp),intent(out) :: ev(na),qq(ldq,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical  :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol,&
&                                    nev=nev,ldq=ldq)
 if (success) then
   call elpa_hdl%set("solver",ELPA_SOLVER_1STAGE,error)
   if (error==0) call elpa_hdl%eigenvectors(aa,ev,qq,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif  (defined HAVE_LINALG_ELPA_2016)
 success=elpa_solve_evp_real_1stage(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,&
&                                   elpa_comm_rows,elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015) || (defined HAVE_LINALG_ELPA_2014)
 success=solve_evp_real(na,nev,aa,lda,ev,qq,ldq,nblk,elpa_comm_rows,elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2013)
 call solve_evp_real(na,nev,aa,lda,ev,qq,ldq,nblk,elpneva_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success=elpa_solve_evp_real_1stage_double(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,&
&                                          elpa_comm_rows,elpa_comm_cols,mpi_comm_elpa,.false.)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (solve_evp_1stage_real)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_solve_evp_1stage_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_solve_evp_1stage_complex
!! NAME
!!  elpa_func_solve_evp_1stage_complex
!!
!! FUNCTION
!!  Wrapper to elpa_solve_evp_complex_1stage ELPA function
!!
!! INPUTS
!!  na=Order of matrix aa
!!  nev=Number of eigenvalues needed.
!!  lda=Leading dimension of aa
!!  ldq=Leading dimension of qq
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  matrixCols=Distributed number of matrix columns
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
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
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_solve_evp_1stage_complex(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,&
&                                             mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_solve_evp_1stage_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,nev,lda,ldq,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
!arrays
 complex(dpc),intent(inout) :: aa(lda,matrixCols)
 real(dp),intent(out) :: ev(na)
 complex(dpc),intent(out) :: qq(ldq,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical  :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol,&
&                                    nev=nev,ldq=ldq)
 if (success) then
   call elpa_hdl%set("solver",ELPA_SOLVER_1STAGE,error)
   if (error==0) call elpa_hdl%eigenvectors(aa,ev,qq,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
 success=elpa_solve_evp_complex_1stage(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,&
&                                      elpa_comm_rows,elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015) || (defined HAVE_LINALG_ELPA_2014)
 success=solve_evp_complex(na,nev,aa,lda,ev,qq,ldq,nblk,elpa_comm_rows,elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2013)
 call solve_evp_complex(na,nev,aa,lda,ev,qq,ldq,nblk,elpa_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success=elpa_solve_evp_complex_1stage_double(na,nev,aa,lda,ev,qq,ldq,nblk,matrixCols,&
&                                             elpa_comm_rows,elpa_comm_cols,mpi_comm_elpa,.false.)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (solve_evp_1stage_complex)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_solve_evp_1stage_complex
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_cholesky_real
!! NAME
!!  elpa_func_cholesky_real
!!
!! FUNCTION
!!  Wrapper to elpa_cholesky_real ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     On return, the upper triangle contains the Cholesky factor
!!                     and the lower triangle is set to 0.
!!
!! PARENTS
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_cholesky_real(na,aa,lda,nblk,matrixCols,&
&                                  mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_cholesky_real'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
!arrays
 real(dp),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol)
 if (success) then
   call elpa_hdl%cholesky(aa,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_cholesky_real(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015)
 call cholesky_real(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call cholesky_real(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2013)
 call cholesky_real(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_cholesky_real_double(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (cholesky_real)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_cholesky_real
!!***

!----------------------------------------------------------------------
!!****f* m_elpa/elpa_func_cholesky_complex
!! NAME
!!  elpa_func_cholesky_complex
!!
!! FUNCTION
!!  Wrapper to elpa_cholesky_complex ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     On return, the upper triangle contains the Cholesky factor
!!                     and the lower triangle is set to 0.
!!
!! PARENTS
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_cholesky_complex(na,aa,lda,nblk,matrixCols,&
&                                     mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_cholesky_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
!arrays
 complex(dpc),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol)
 if (success) then
   call elpa_hdl%cholesky(aa,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_cholesky_complex(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015)
 call cholesky_complex(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call cholesky_complex(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2013)
 call cholesky_complex(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_cholesky_complex_double(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (cholesky_complex)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_cholesky_complex
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_invert_triangular_real
!! NAME
!!  elpa_func_invert_triangular_real
!!
!! FUNCTION
!!  Wrapper to elpa_invert_triangular_real ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     The lower triangle is not referenced.
!!
!! PARENTS
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_invert_triangular_real(na,aa,lda,nblk,matrixCols,&
&                                           mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_invert_triangular_real'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
!arrays
 real(dp),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol)
 if (success) then
   call elpa_hdl%invert_triangular(aa,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_invert_trm_real(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015)
 call invert_trm_real(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call invert_trm_real(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,success)
 success=.true. ! Sometimes get unexpected success=false
#elif (defined HAVE_LINALG_ELPA_2013)
 call invert_trm_real(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_invert_trm_real_double(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (invert_triangular_real)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_invert_triangular_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_invert_triangular_complex
!! NAME
!!  elpa_func_invert_triangular_complex
!!
!! FUNCTION
!!  Wrapper to elpa_invert_triangular_complex ELPA function
!!
!! INPUTS
!!  na=Order of matrix
!!  lda=Leading dimension of aa
!!  matrixCols=local columns of matrix a
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! SIDE EFFECTS
!!  aa(lda,matrixCols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     The lower triangle is not referenced.
!!
!! PARENTS
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_invert_triangular_complex(na,aa,lda,nblk,matrixCols,&
&                                              mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_invert_triangular_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,lda,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
!arrays
 complex(dpc),intent(inout) :: aa(lda,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol)
 if (success) then
   call elpa_hdl%invert_triangular(aa,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_invert_trm_complex(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015)
 call invert_trm_complex(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call invert_trm_complex(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols,success)
 success=.true. ! Sometimes get unexpected success=false
#elif (defined HAVE_LINALG_ELPA_2013)
 call invert_trm_complex(na,aa,lda,nblk,elpa_comm_rows,elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_invert_trm_complex_double(na,aa,lda,nblk,matrixCols,elpa_comm_rows,elpa_comm_cols,.false.)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (invert_triangular_complex)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_invert_triangular_complex
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_hermitian_multiply_real
!! NAME
!!  elpa_func_hermitian_multiply_real
!!
!! FUNCTION
!!  Wrapper to elpa_hermitian_multiply_real ELPA function
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
!!  mpi_comm_elpa=Global communicator for the calculations
!!  my_prow=Row coordinate of the calling process in the process grid (in)
!!  my_pcol=Column coordinate of the calling process in the process grid (in)
!!
!! OUTPUT
!!  cc=Matrix C
!!
!! PARENTS
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_hermitian_multiply_real(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,cc,ldc,&
&                                            nblk,matrixCols,mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_hermitian_multiply_real'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,ncb,lda,ldb,ldc,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
 character*1 :: uplo_a, uplo_c
!arrays
 real(dp),intent(in) :: aa(lda,matrixCols),bb(ldb,matrixCols)
 real(dp),intent(out) :: cc(ldc,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol)
 if (success) then
   call elpa_hdl%hermitian_multiply(uplo_a,uplo_c,ncb,aa,bb,ldb,matrixCols,cc,ldc,matrixCols,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_mult_at_b_real(uplo_a,uplo_c,na,ncb,aa,lda,matrixCols,bb,ldb,matrixCols,nblk, &
&                              elpa_comm_rows,elpa_comm_cols,cc,ldc,matrixCols)
#elif (defined HAVE_LINALG_ELPA_2015) || (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
 call mult_at_b_real(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,nblk, &
&                    elpa_comm_rows,elpa_comm_cols,cc,ldc)
#else
!ELPA-LEGACY-2017
 success =  elpa_mult_at_b_real_double(uplo_a,uplo_c,na,ncb,aa,lda,matrixCols,bb,ldb,matrixCols,nblk, &
&                                      elpa_comm_rows,elpa_comm_cols,cc,ldc,matrixCols)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (hermitian_multiply_real)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_hermitian_multiply_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_hermitian_multiply_complex
!! NAME
!!  elpa_func_hermitian_multiply_complex
!!
!! FUNCTION
!!  Wrapper to elpa_hermitian_multiply_complex ELPA function
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
!!  elpa_comm_rows=MPI-Communicator for rows
!!  elpa_comm_cols=MPI-Communicator for columns
!!
!! OUTPUT
!!  cc=Matrix C
!!
!! PARENTS
!!
!! CHILDREN
!!      elpa_deallocate,elpa_hdl%hermitian_multiply,elpa_hdl%set
!!      mult_ah_b_complex
!!
!! SOURCE

subroutine elpa_func_hermitian_multiply_complex(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,cc,ldc, &
&                                               nblk,matrixCols,mpi_comm_elpa,my_prow,my_pcol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_hermitian_multiply_complex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: na,ncb,lda,ldb,ldc,nblk,matrixCols
 integer,intent(in)  :: mpi_comm_elpa,my_prow,my_pcol
 character*1 :: uplo_a, uplo_c
!arrays
 complex(dpc),intent(in) :: aa(lda,matrixCols),bb(ldb,matrixCols)
 complex(dpc),intent(out) :: cc(ldc,matrixCols)

!Local variables-------------------------------
 integer :: error
 logical :: success
 character(len=100) :: errmsg

! *********************************************************************

 success=.true. ; error=0 ; errmsg=""

#if  (defined HAVE_LINALG_ELPA_2017)
 elpa_hdl => elpa_allocate()
 success=elpa_func_set_matrix_params(errmsg,na,nblk,lda,matrixCols,mpi_comm_elpa,my_prow,my_pcol)
 if (success) then
   call elpa_hdl%hermitian_multiply(uplo_a,uplo_c,ncb,aa,bb,ldb,matrixCols,cc,ldc,matrixCols,error)
   success=(error==0)
 end if
 call elpa_deallocate(elpa_hdl)
#elif (defined HAVE_LINALG_ELPA_2016)
  success = elpa_mult_ah_b_complex(uplo_a,uplo_c,na,ncb,aa,lda,matrixCols,bb,ldb,matrixCols,nblk, &
&                                  elpa_comm_rows,elpa_comm_cols,cc,ldc,matrixCols)
#elif (defined HAVE_LINALG_ELPA_2015) || (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
 call mult_ah_b_complex(uplo_a,uplo_c,na,ncb,aa,lda,bb,ldb,nblk, &
&                       elpa_comm_rows,elpa_comm_cols,cc,ldc)
#else
!ELPA-LEGACY-2017
  success = elpa_mult_ah_b_complex_double(uplo_a,uplo_c,na,ncb,aa,lda,matrixCols,bb,ldb,matrixCols,nblk, &
&                                         elpa_comm_rows,elpa_comm_cols,cc,ldc,matrixCols)
#endif

 if (.not.success) then
   if (trim(errmsg)=="") errmsg="Problem with ELPA (hermitian_multiply_complex)!"
   MSG_ERROR(errmsg)
 end if

end subroutine elpa_func_hermitian_multiply_complex
!!***

!----------------------------------------------------------------------

#ifdef HAVE_ELPA_FORTRAN2008

!!****f* m_elpa/elpa_func_set_matrix_params
!! NAME
!!  elpa_func_set_matrix_params
!!
!! FUNCTION
!!  Set parameters decribing a matrix and it's MPI distribution
!!
!! INPUTS
!!  na=Order of matrix A
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  local_nrows=Leading dimension of A
!!  local_ncols=Local columns of matrixes A and Q (eigenvectors)
!!  mpi_comm_parent=Global communicator for the calculations
!!  process_row=Row coordinate of the calling process in the process grid
!!  process_col=Column coordinate of the calling process in the process grid
!!  [nev]= -- optional -- Number of (smallest) eigenvalues/eigenvectors to be computed
!!  [ldq]= -- optional -- For testing purpose: leading dimension of Q
!!  [gpu]= -- optional -- Flag (0 or 1): use GPU version
!!
!! OUTPUT
!!  errmsg=error message if any
!!  success=.TRUE. if successful
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function elpa_func_set_matrix_params(errmsg,na,nblk,local_nrows,local_ncols,&
&             mpi_comm_parent,process_row,process_col,nev,ldq,gpu) &
& result(success)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpa_func_set_matrix_params'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: na,nblk,local_nrows,local_ncols
 integer,intent(in) :: mpi_comm_parent,process_row,process_col
 integer,intent(in),optional :: nev,ldq,gpu
 logical :: success
 character(len=*) :: errmsg
!arrays

!Local variables-------------------------------
 integer :: err

! *********************************************************************

 success=.true. ; errmsg=""

 if (present(ldq)) then
   if (local_nrows/=ldq) then
     MSG_BUG('Problem with ELPA (lda/=ldq)!')
   end if
 end if

 call elpa_hdl%set("na",na,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (na)!";success=.false.;return
 endif

 call elpa_hdl%set("nblk",nblk,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (nblk)!";success=.false.;return
 endif

 call elpa_hdl%set("local_nrows",local_nrows,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (local_nrows)!";success=.false.;return
 endif

 call elpa_hdl%set("local_ncols",local_ncols,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (local_ncols)!";success=.false.;return
 endif

 call elpa_hdl%set("process_row",process_row,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (process_row)!";success=.false.;return
 endif

 call elpa_hdl%set("process_col",process_col,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (process_col)!";success=.false.;return
 endif

 call elpa_hdl%set("mpi_comm_parent",mpi_comm_parent,err)
 if (err/=0) then
   errmsg="Cannot setup ELPA instance (mpi_comm_parent)!";success=.false.;return
 endif

 if (present(nev)) then
   call elpa_hdl%set("nev",nev,err)
   if (err/=0) then
     errmsg="Cannot setup ELPA instance (nev)!";success=.false.;return
   endif
 end if

 if (present(gpu)) then
   call elpa_hdl%set("gpu",gpu,err)
   if (err/=0) then
     errmsg="Cannot setup ELPA instance (gpu)!";success=.false.;return
   endif
 end if

 if (elpa_hdl%setup()/=ELPA_OK) then
   errmsg="Cannot setup ELPA instance!"
   success=.false.;return
 endif

end function elpa_func_set_matrix_params
!!***

!----------------------------------------------------------------------
#endif

#endif


END MODULE m_elpa
!!***
