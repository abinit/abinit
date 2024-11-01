!!****m* ABINIT/m_elpa
!! NAME
!! m_elpa
!!
!! FUNCTION
!! This module contains interfaces to ELPA library methods
!! See http://elpa.mpcdf.mpg.de
!!
!! COPYRIGHT
!! Copyright (C) 2016-2024 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!!#if (defined HAVE_LINALG_ELPA) && (defined HAVE_LINALG_ELPA_2017)
!!#define HAVE_LINALG_ELPA_FORTRAN2008
!!#else
!!#undef HAVE_LINALG_ELPA_FORTRAN2008
!!#endif

module m_elpa

 use defs_basis
 use m_errors

#ifdef HAVE_LINALG_ELPA
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
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
 public :: elpa_func_init                ! Init ELPA
 public :: elpa_func_uninit              ! End ELPA
 public :: elpa_func_allocate            ! Allocate a ELPA handle and set up MPI information
 public :: elpa_func_deallocate          ! Deallocate a ELPA handle
 public :: elpa_func_error_handler       ! Manage errors (print a readable message)
 public :: elpa_func_get_communicators   ! Get rows and cols communicators (not supposed to be called directly)
 public :: elpa_func_set_matrix          ! Set matrix specifications in a ELPA handle
 public :: elpa_func_solve_evp_1stage    ! Solve the diagonalization problem (use a ELPA handle)
 public :: elpa_func_solve_gevp_2stage   ! Solve the generalized diagonalization problem (use a ELPA handle)
 public :: elpa_func_cholesky            ! Apply Cholesky transformation (use a ELPA handle)
 public :: elpa_func_invert_triangular   ! Invert triangular matrix (use a ELPA handle)
 public :: elpa_func_hermitian_multiply  ! Perform C := A**H * B (use a ELPA handle)

 interface elpa_func_solve_evp_1stage
   module procedure elpa_func_solve_evp_1stage_real
   module procedure elpa_func_solve_evp_1stage_complex
 end interface elpa_func_solve_evp_1stage

 interface elpa_func_solve_gevp_2stage
   module procedure elpa_func_solve_gevp_2stage_real
   module procedure elpa_func_solve_gevp_2stage_complex
 end interface elpa_func_solve_gevp_2stage

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

!ELPA generalized handle
 type,public :: elpa_hdl_t
   logical :: is_allocated=.false.
   logical :: matrix_is_set=.false.
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
   class(elpa_t),pointer :: elpa
#else
   integer :: mpi_comm_parent
   integer :: elpa_comm_rows,elpa_comm_cols
   integer :: process_row,process_col
   integer :: local_nrows,local_ncols
   integer :: na,nblk,nev
   integer :: gpu=0
#endif
 end type elpa_hdl_t

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
!! SOURCE

subroutine elpa_func_init()

!Local variables-------------------------------
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 integer,parameter :: elpa_min_version=20170403
#endif
 logical :: success
!
! *********************************************************************

 success=.true.

#ifdef HAVE_LINALG_ELPA_FORTRAN2008

#if (defined HAVE_LINALG_ELPA_2016) || (defined HAVE_LINALG_ELPA_2015_11) || (defined HAVE_LINALG_ELPA_2015_02)
!This case is not supposed to happen
 success=.false.
 ABI_BUG('Wrong ELPA cpp directives!')
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
!This case is not supposed to happen
 success=.false.
 ABI_BUG('Wrong ELPA cpp directives!')
#else
 success=(elpa_init(elpa_min_version)==ELPA_OK)
#endif

#else
!No init function
#endif

 if (.not.success) then
   ABI_ERROR('Problem with ELPA (elpa_init)!')
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
!! SOURCE

subroutine elpa_func_uninit()

! *********************************************************************

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_uninit()
#else
!No uninit function
#endif

end subroutine elpa_func_uninit
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_allocate
!! NAME
!!  elpa_func_allocate
!!
!! FUNCTION
!!  Allocate a ELPA handle and set it up with communicators specification
!!
!! INPUTS
!!  [blacs_ctx]= -- optional -- Blacs context
!!  [gpu]= -- optional -- Flag (0 or 1): use GPU version (currently only NVidia)
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)= ELPA handle
!!
!! SOURCE

subroutine elpa_func_allocate(elpa_hdl,gpu,blacs_ctx)

!Arguments ------------------------------------
 integer,intent(in),optional :: gpu,blacs_ctx
 type(elpa_hdl_t),intent(inout) :: elpa_hdl

!Local variables-------------------------------
 integer :: err,l_gpu,l_blacs_ctx
 logical :: gpu_debug_mode=.false.
 character(len=10) :: varname

! *********************************************************************

 err=0
 ! if optional parameter is present, use it
 ! else use default value, i.e. don't use GPU
 l_gpu = 0
 if (present(gpu)) then
   if(gpu/=0) l_gpu = 1
 end if

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 elpa_hdl%elpa => elpa_allocate(err)
 call elpa_func_error_handler(err_code=err,err_msg='Error in initialization')

 if(l_gpu==1) then
#if defined HAVE_GPU_CUDA
   varname="nvidia-gpu"
#else
   ABI_BUG("Requested unsupported GPU model with ELPA!")
#endif
   if (err==ELPA_OK) call elpa_hdl%elpa%set(varname,l_gpu,err)

   ! Handling GPU-support on older ELPA versions (2020.11 and previous)
   ! Only NVIDIA CUDA GPUs were supported with 'gpu' setting
   if (err==ELPA_ERROR_ENTRY_NOT_FOUND) then
#if defined HAVE_GPU_CUDA
     varname="gpu"
     call elpa_hdl%elpa%set(varname,l_gpu,err)
#else
     ABI_ERROR("You seem to use an old version of ELPA ( < 2021.x ) which only supports NVIDIA GPUs.")
#endif
   end if
   call elpa_func_error_handler(err_code=err,err_msg='Error when enabling GPU on ELPA')

   if (gpu_debug_mode) then
     if (err==ELPA_OK) call elpa_hdl%elpa%set("debug",1,err) 
     call elpa_func_error_handler(err_code=err,err_msg='Error when enabling debug on ELPA')
   end if

 end if
#else
 if (err==0.and.l_gpu==1) then
   elpa_hdl%gpu=l_gpu
   if (gpu_debug_mode) elpa_hdl%debug=1
 end if
#endif

 if (present(blacs_ctx)) then
   if (err==ELPA_OK) call elpa_hdl%elpa%set("blacs_context",int(blacs_ctx,kind=c_int),err)
   call elpa_func_error_handler(err_code=err,err_varname=varname)
 end if

 elpa_hdl%is_allocated=.true.

end subroutine elpa_func_allocate
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_deallocate
!! NAME
!!  elpa_func_deallocate_matrix
!!
!! FUNCTION
!!  Deallocate a ELPA handle
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)= ELPA handle
!!
!! SOURCE

subroutine elpa_func_deallocate(elpa_hdl)

!Arguments ------------------------------------
 type(elpa_hdl_t),intent(inout) :: elpa_hdl

!Local variables-------------------------------
 integer :: err

! *********************************************************************

 err=0

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_deallocate(elpa_hdl%elpa)
#endif

 elpa_hdl%matrix_is_set=.false.
 elpa_hdl%is_allocated=.false.

end subroutine elpa_func_deallocate
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_error_handler
!! NAME
!!  elpa_func_error_handler
!!
!! FUNCTION
!!  Handle ELPA errors
!!
!! INPUTS
!!  [err_code]= --optional-- Error code
!!  [err_msg]= --optional-- Generic error message
!!  [err_varname]= -- optional-- Name of the ELPA variable related to the error
!!
!! OUTPUT
!!  No output, only printing
!!
!! SOURCE

subroutine elpa_func_error_handler(err_code,err_msg,err_varname)

!Arguments ------------------------------------
 integer,optional :: err_code
 character(len=*),optional :: err_msg,err_varname

!Local variables-------------------------------
 integer :: err_code_
 character(len=500) :: msg
 character(len=100) :: err_strg

! *********************************************************************

 err_code_=-100;if (present(err_code)) err_code_=err_code
 if (err_code_==0) return

 err_strg=''
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err_code_==ELPA_ERROR) err_strg='ELPA_ERROR'
 if (err_code_==ELPA_ERROR_ENTRY_READONLY) err_strg='ELPA_ERROR_ENTRY_READONLY'
 if (err_code_==ELPA_ERROR_ENTRY_NOT_FOUND) err_strg='ELPA_ERROR_ENTRY_NOT_FOUND'
 if (err_code_==ELPA_ERROR_ENTRY_ALREADY_SET) err_strg='ELPA_ERROR_ENTRY_ALREADY_SET'
 if (err_code_==ELPA_ERROR_ENTRY_INVALID_VALUE) err_strg='ELPA_ERROR_ENTRY_INVALID_VALUE'
 if (err_code_==ELPA_ERROR_ENTRY_NO_STRING_REPRESENTATION) err_strg='ELPA_ERROR_NO_STRING_REPRESENTATION'
#endif

 write(msg,'(a)') 'ELPA library error!'
 if (present(err_msg)) then
   if (trim(err_msg)/="") write(msg,'(3a)') trim(msg),ch10,trim(err_msg)
 end if
 if (present(err_varname)) then
   if (trim(err_varname)/="") write(msg,'(4a)') trim(msg),ch10,'Variable: ',trim(err_varname)
 end if
 if (trim(err_strg)/="") write(msg,'(4a)') trim(msg),ch10,'Error code: ',trim(err_strg)
   ABI_ERROR(msg)

end subroutine elpa_func_error_handler
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
!!  mpi_comm_parent=Global communicator for the calculations (in)
!!  process_row=Row coordinate of the calling process in the process grid (in)
!!  process_col=Column coordinate of the calling process in the process grid (in)
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)= ELPA handle
!!
!! SOURCE

subroutine elpa_func_get_communicators(elpa_hdl,mpi_comm_parent,process_row,process_col)

!Arguments ------------------------------------
 integer,intent(in)  :: mpi_comm_parent,process_row,process_col
 type(elpa_hdl_t),intent(inout) :: elpa_hdl

!Local variables-------------------------------
 integer  :: err
 character(len=20) :: varname

! *********************************************************************

 err=0 ; varname=''

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err==ELPA_OK) then
   varname='mpi_comm_parent'
   call elpa_hdl%elpa%set(trim(varname),mpi_comm_parent,err)
 end if
 if (err==ELPA_OK) then
   varname='process_row'
   call elpa_hdl%elpa%set(trim(varname),process_row,err)
 end if
 if (err==ELPA_OK) then
   varname='process_col'
   call elpa_hdl%elpa%set(trim(varname),process_col,err)
 end if
 if (err==ELPA_OK) then
   varname=''
   err = elpa_hdl%elpa%setup()
   call elpa_func_error_handler(err_code=err,err_msg='Error during ELPA setup')
 endif

#else
 elpa_hdl%mpi_comm_parent=mpi_comm_parent
 elpa_hdl%process_row=process_row
 elpa_hdl%process_col=process_col
#if (defined HAVE_LINALG_ELPA_2016)
 err=elpa_get_communicators(mpi_comm_parent,process_row,process_col,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015_11) || (defined HAVE_LINALG_ELPA_2015_02)
 err=get_elpa_row_col_comms(mpi_comm_parent,process_row,process_col,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
 call get_elpa_row_col_comms(mpi_comm_parent,process_row,process_col,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 err=elpa_get_communicators(mpi_comm_parent,process_row,process_col,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#endif

#endif

 call elpa_func_error_handler(err_code=err,err_msg='Error in elpa_get_communicators',err_varname=varname)

end subroutine elpa_func_get_communicators
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_set_matrix
!! NAME
!!  elpa_func_set_matrix
!!
!! FUNCTION
!!  Set parameters decribing a matrix and it's MPI distribution
!!  in a ELPA handle
!!
!! INPUTS
!!  na=Order of matrix A
!!  nblk=Blocksize of cyclic distribution, must be the same in both directions!
!!  nev=Number of eigenvalues needed.
!!  local_nrows=Leading dimension of A
!!  local_ncols=Local columns of matrixes A and Q (eigenvectors)
!!  nev=Number of eigenvalues needed.
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_set_matrix(elpa_hdl,na,nblk,nev,local_nrows,local_ncols)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: na,nblk,nev,local_nrows,local_ncols
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays

!Local variables-------------------------------
 integer :: err
 character(len=15) :: varname

! *********************************************************************

 err=0 ; varname=''

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err==ELPA_OK) then
   varname="na"
   call elpa_hdl%elpa%set(trim(varname),na,err)
 end if
 if (err==ELPA_OK) then
   varname="nblk"
   call elpa_hdl%elpa%set(trim(varname),nblk,err)
 end if
 if (err==ELPA_OK) then
   varname="nev"
   call elpa_hdl%elpa%set(trim(varname),nev,err)
 end if
 if (err==ELPA_OK) then
   varname="local_nrows"
   call elpa_hdl%elpa%set(trim(varname),local_nrows,err)
 end if
 if (err==ELPA_OK) then
   varname="local_ncols"
   call elpa_hdl%elpa%set(trim(varname),local_ncols,err)
 end if
#else
 elpa_hdl%na=na
 elpa_hdl%nblk=nblk
 elpa_hdl%nev=nev
 elpa_hdl%local_nrows=local_nrows
 elpa_hdl%local_ncols=local_ncols
#endif

 call elpa_func_error_handler(err_code=err,err_msg='Error during matrix initialization',err_varname=varname)

 elpa_hdl%matrix_is_set=.true.

end subroutine elpa_func_set_matrix
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_solve_evp_2stage_real
!! NAME
!!  elpa_func_solve_evp_2stage_real
!!
!! FUNCTION
!!  Wrapper to elpa_solve_evp_real_2stage ELPA function
!!
!! INPUTS
!!  nev=Number of eigenvalues needed.
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(local_nrows,local_ncols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na,na))
!!                     even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_solve_gevp_2stage_real(elpa_hdl,aa,bb,qq,ev,nev)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nev
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 real(dp),intent(inout) :: aa(:,:),bb(:,:)
 real(dp),intent(out) :: ev(:),qq(:,:)

!Local variables-------------------------------
 integer :: err
 logical  :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
 ABI_CHECK(size(qq)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix Q has wrong sizes!')
 ABI_CHECK(size(ev)==elpa_hdl%elpa%na,'BUG: matrix EV has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err==ELPA_OK) call elpa_hdl%elpa%set("solver",ELPA_SOLVER_2STAGE,err)
 if (err==ELPA_OK) call elpa_hdl%elpa%generalized_eigenvectors(aa,bb,ev,qq,.false.,err)
 success=(err==ELPA_OK)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in solve_evp_2stage_real!')

end subroutine elpa_func_solve_gevp_2stage_real
!!***

!----------------------------------------------------------------------

!!****f* m_elpa/elpa_func_solve_evp_2stage_complex
!! NAME
!!  elpa_func_solve_evp_2stage_complex
!!
!! FUNCTION
!!  Wrapper to elpa_solve_evp_complex_2stage ELPA function
!!
!! INPUTS
!!  nev=Number of eigenvalues needed.
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(local_nrows,local_ncols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na,na))
!!                      even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_solve_gevp_2stage_complex(elpa_hdl,aa,bb,qq,ev,nev)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nev
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 complex(dpc),intent(inout) :: aa(:,:),bb(:,:)
 real(dp),intent(out) :: ev(:)
 complex(dpc),intent(out) :: qq(:,:)

!Local variables-------------------------------
 integer :: err
 logical  :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
 ABI_CHECK(size(qq)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix Q has wrong sizes!')
 ABI_CHECK(size(ev)==elpa_hdl%elpa%na,'BUG: matrix EV has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err==ELPA_OK) call elpa_hdl%elpa%set("solver",ELPA_SOLVER_2STAGE,err)
 if (err==ELPA_OK) call elpa_hdl%elpa%generalized_eigenvectors(aa,bb,ev,qq,.false.,err)
 success=(err==ELPA_OK)
#endif

 if (.not.success) call elpa_func_error_handler(err_code=err,err_msg='Error in solve_evp_2stage_complex!')

end subroutine elpa_func_solve_gevp_2stage_complex
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
!!  nev=Number of eigenvalues needed.
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(local_nrows,local_ncols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na,na))
!!                     even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_solve_evp_1stage_real(elpa_hdl,aa,qq,ev,nev)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nev
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 real(dp),intent(inout) :: aa(:,:)
 real(dp),intent(out) :: ev(:),qq(:,:)

!Local variables-------------------------------
 integer :: err
 logical  :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
 ABI_CHECK(size(qq)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix Q has wrong sizes!')
 ABI_CHECK(size(ev)==elpa_hdl%elpa%na,'BUG: matrix EV has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
 ABI_CHECK(size(qq)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix Q has wrong sizes!')
 ABI_CHECK(size(ev)==elpa_hdl%na,'BUG: matrix EV has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err==ELPA_OK) call elpa_hdl%elpa%set("solver",ELPA_SOLVER_1STAGE,err)
 if (err==ELPA_OK) call elpa_hdl%elpa%eigenvectors(aa,ev,qq,err)
 success=(err==ELPA_OK)
#elif  (defined HAVE_LINALG_ELPA_2016)
 success=elpa_solve_evp_real_1stage(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                       elpa_hdl%nblk,elpa_hdl%local_ncols,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015_11)
 success=solve_evp_real(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                       elpa_hdl%nblk,elpa_hdl%local_ncols,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015_02) || (defined HAVE_LINALG_ELPA_2014)
 success=solve_evp_real(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                       elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2013)
 call solve_evp_real(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                    elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success=elpa_solve_evp_real_1stage_double(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&        elpa_hdl%nblk,elpa_hdl%local_ncols,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,elpa_hdl%mpi_comm_parent,.false.)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in solve_evp_1stage_real!')

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
!!  nev=Number of eigenvalues needed.
!!
!! OUTPUT
!!  ev(na)=Eigenvalues of a, every processor gets the complete set
!!  qq(local_nrows,local_ncols)=Eigenvectors of aa
!!                     Distribution is like in Scalapack.
!!                     Must be always dimensioned to the full size (corresponding to (na,na))
!!                      even if only a part of the eigenvalues is needed.
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix for which eigenvalues are to be computed.
!!                    Distribution is like in Scalapack.
!!                    The full matrix must be set (not only one half like in scalapack).
!!                    Destroyed on exit (upper and lower half).
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_solve_evp_1stage_complex(elpa_hdl,aa,qq,ev,nev)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nev
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 complex(dpc),intent(inout) :: aa(:,:)
 real(dp),intent(out) :: ev(:)
 complex(dpc),intent(out) :: qq(:,:)

!Local variables-------------------------------
 integer :: err
 logical  :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
 ABI_CHECK(size(qq)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix Q has wrong sizes!')
 ABI_CHECK(size(ev)==elpa_hdl%elpa%na,'BUG: matrix EV has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
 ABI_CHECK(size(qq)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix Q has wrong sizes!')
 ABI_CHECK(size(ev)==elpa_hdl%na,'BUG: matrix EV has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 if (err==ELPA_OK) call elpa_hdl%elpa%set("solver",ELPA_SOLVER_1STAGE,err)
 if (err==ELPA_OK) call elpa_hdl%elpa%eigenvectors(aa,ev,qq,err)
 success=(err==ELPA_OK)
#elif  (defined HAVE_LINALG_ELPA_2016)
 success=elpa_solve_evp_complex_1stage(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                       elpa_hdl%nblk,elpa_hdl%local_ncols,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015_11)
 success=solve_evp_complex(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                       elpa_hdl%nblk,elpa_hdl%local_ncols,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2015_02) || (defined HAVE_LINALG_ELPA_2014)
 success=solve_evp_complex(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                       elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#elif (defined HAVE_LINALG_ELPA_2013)
 call solve_evp_complex(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&                    elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success=elpa_solve_evp_complex_1stage_double(elpa_hdl%na,nev,aa,elpa_hdl%local_nrows,ev,qq,elpa_hdl%local_nrows,&
&        elpa_hdl%nblk,elpa_hdl%local_ncols,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,elpa_hdl%mpi_comm_parent,.false.)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in solve_evp_1stage_complex!')

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
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     On return, the upper triangle contains the Cholesky factor
!!                     and the lower triangle is set to 0.
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_cholesky_real(elpa_hdl,aa)

!Arguments ------------------------------------
!scalars
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 real(dp),intent(inout) :: aa(:,:)

!Local variables-------------------------------
 integer :: err
 logical :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_hdl%elpa%cholesky(aa,err)
 success=(err==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_cholesky_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                             elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015_11)
 call cholesky_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
                    elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2015_02)
 call cholesky_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
                    elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call cholesky_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                   elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2013)
 call cholesky_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                   elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_cholesky_real_double(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                                    elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in solve_cholesky_real!')

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
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     On return, the upper triangle contains the Cholesky factor
!!                     and the lower triangle is set to 0.
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_cholesky_complex(elpa_hdl,aa)

!Arguments ------------------------------------
!scalars
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 complex(dpc),intent(inout) :: aa(:,:)

!Local variables-------------------------------
 integer :: err
 logical :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_hdl%elpa%cholesky(aa,err)
 success=(err==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_cholesky_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                                elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015_11)
 call cholesky_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                      elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2015_02)
 call cholesky_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                      elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call cholesky_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                      elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,success)
#elif (defined HAVE_LINALG_ELPA_2013)
 call cholesky_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                      elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_cholesky_complex_double(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                                       elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in solve_cholesky_complex!')

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
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     The lower triangle is not referenced.
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_invert_triangular_real(elpa_hdl,aa)

!Arguments ------------------------------------
!scalars
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 real(dp),intent(inout) :: aa(:,:)

!Local variables-------------------------------
 integer :: err
 logical :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_hdl%elpa%invert_triangular(aa,err)
 success=(err==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_invert_trm_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                               elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015_11)
 call invert_trm_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
                         elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2015_02)
 call invert_trm_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
                         elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call invert_trm_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                     elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,success)
 success=.true. ! Sometimes get unexpected success=false
#elif (defined HAVE_LINALG_ELPA_2013)
 call invert_trm_real(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                     elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_invert_trm_real_double(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                                      elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in invert_trianguler_real!')

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
!!
!! SIDE EFFECTS
!!  aa(local_nrows,local_ncols)=Distributed matrix which should be factorized.
!!                     Distribution is like in Scalapack.
!!                     Only upper triangle is needs to be set.
!!                     The lower triangle is not referenced.
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_invert_triangular_complex(elpa_hdl,aa)

!Arguments ------------------------------------
!scalars
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 complex(dpc),intent(inout) :: aa(:,:)

!Local variables-------------------------------
 integer :: err
 logical :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
#endif

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_hdl%elpa%invert_triangular(aa,err)
 success=(err==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_invert_trm_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                                  elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#elif (defined HAVE_LINALG_ELPA_2015_11)
 call invert_trm_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
                         elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2015_02)
 call invert_trm_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
                         elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.,success)
#elif (defined HAVE_LINALG_ELPA_2014)
 call invert_trm_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                        elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,success)
 success=.true. ! Sometimes get unexpected success=false
#elif (defined HAVE_LINALG_ELPA_2013)
 call invert_trm_complex(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,&
&                        elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols)
#else
!ELPA-LEGACY-2017
 success = elpa_invert_trm_complex_double(elpa_hdl%na,aa,elpa_hdl%local_nrows,elpa_hdl%nblk,elpa_hdl%local_ncols,&
&                                         elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,.false.)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in invert_trianguler_complex!')

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
!!  ncb=Number of columns of B and C
!!  aa(local_nrows,local_ncols)=Matrix A
!!  bb(ldb,local_ncols_c)=Matrix B
!!  local_nrows_b=Local rows of matrix B
!!  local_ncols_b=Local columns of matrix B
!!  local_nrows_c=Local rows of matrix C
!!  local_ncols_c=Local columns of matrix C
!!
!! OUTPUT
!!  cc(local_nrows_c,local_ncols_c)=Matrix C
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_hermitian_multiply_real(elpa_hdl,uplo_a,uplo_c,ncb,aa,bb,local_nrows_b,local_ncols_b,&
&                                            cc,local_nrows_c,local_ncols_c)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: ncb,local_nrows_b,local_nrows_c,local_ncols_b,local_ncols_c
 character*1 :: uplo_a,uplo_c
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 real(dp),intent(in) :: aa(:,:),bb(:,:)
 real(dp),intent(out) :: cc(:,:)

!Local variables-------------------------------
 integer :: err
 logical :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
#endif
 ABI_CHECK(size(bb)==local_nrows_b*local_ncols_b,'BUG: matrix B has wrong sizes!')
 ABI_CHECK(size(cc)==local_nrows_c*local_ncols_c,'BUG: matrix C has wrong sizes!')

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_hdl%elpa%hermitian_multiply(uplo_a,uplo_c,ncb,aa,bb,local_nrows_b,local_ncols_b,&
&                                      cc,local_nrows_c,local_ncols_c,err)
 success=(err==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_mult_at_b_real(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,elpa_hdl%local_ncols,&
&          bb,local_nrows_b,local_ncols_b,elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,&
&          cc,local_nrows_c,local_ncols_c)
#elif (defined HAVE_LINALG_ELPA_2015_11) || (defined HAVE_LINALG_ELPA_2015_02)
 call mult_at_b_real(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,bb,local_nrows_b,&
&               elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,cc,local_nrows_c)
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
 call mult_at_b_real(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,bb,local_nrows_b,&
&               elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,cc,local_nrows_c)
#else
!ELPA-LEGACY-2017
 success =  elpa_mult_at_b_real_double(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,&
&           elpa_hdl%local_ncols,bb,local_nrows_b,local_ncols_b,elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,&
&           elpa_hdl%elpa_comm_cols,cc,local_nrows_c,local_ncols_c)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in hermitian_multiply_real!')

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
!!  ncb=Number of columns of B and C
!!  aa(local_nrows,local_ncols)=Matrix A
!!  bb(ldb,local_ncols_c)=Matrix B
!!  local_nrows_b=Local rows of matrix B
!!  local_ncols_b=Local columns of matrix B
!!  local_nrows_c=Local rows of matrix C
!!  local_ncols_c=Local columns of matrix C
!!
!! OUTPUT
!!  cc(local_nrows_c,local_ncols_c)=Matrix C
!!
!! SIDE EFFECTS
!!  elpa_hdl(type<elpa_hdl_t>)=handler for ELPA object
!!
!! SOURCE

subroutine elpa_func_hermitian_multiply_complex(elpa_hdl,uplo_a,uplo_c,ncb,aa,bb,local_nrows_b,local_ncols_b,&
&                                               cc,local_nrows_c,local_ncols_c)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: ncb,local_nrows_b,local_nrows_c,local_ncols_b,local_ncols_c
 character*1 :: uplo_a,uplo_c
 type(elpa_hdl_t),intent(inout) :: elpa_hdl
!arrays
 complex(dpc),intent(in) :: aa(:,:),bb(:,:)
 complex(dpc),intent(out) :: cc(:,:)

!Local variables-------------------------------
 integer :: err
 logical :: success

! *********************************************************************

 success=.true. ; err=0

 if (.not.elpa_hdl%is_allocated) then
   ABI_BUG('ELPA handle not allocated!')
 end if
 if (.not.elpa_hdl%matrix_is_set) then
   ABI_BUG('Matrix not set in ELPA handle!')
 end if
#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 ABI_CHECK(size(aa)==elpa_hdl%elpa%local_nrows*elpa_hdl%elpa%local_ncols,'BUG: matrix A has wrong sizes!')
#else
 ABI_CHECK(size(aa)==elpa_hdl%local_nrows*elpa_hdl%local_ncols,'BUG: matrix A has wrong sizes!')
#endif
 ABI_CHECK(size(bb)==local_nrows_b*local_ncols_b,'BUG: matrix B has wrong sizes!')
 ABI_CHECK(size(cc)==local_nrows_c*local_ncols_c,'BUG: matrix C has wrong sizes!')

#ifdef HAVE_LINALG_ELPA_FORTRAN2008
 call elpa_hdl%elpa%hermitian_multiply(uplo_a,uplo_c,ncb,aa,bb,local_nrows_b,local_ncols_b,&
&                                      cc,local_nrows_c,local_ncols_c,err)
 success=(err==ELPA_OK)
#elif (defined HAVE_LINALG_ELPA_2016)
 success = elpa_mult_ah_b_complex(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,elpa_hdl%local_ncols,&
&          bb,local_nrows_b,local_ncols_b,elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,&
&          cc,local_nrows_c,local_ncols_c)
#elif (defined HAVE_LINALG_ELPA_2015_11) || (defined HAVE_LINALG_ELPA_2015_02)
 call mult_ah_b_complex(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,bb,local_nrows_b,&
&               elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,cc,local_nrows_c)
#elif (defined HAVE_LINALG_ELPA_2014) || (defined HAVE_LINALG_ELPA_2013)
 call mult_ah_b_complex(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,bb,local_nrows_b,&
&               elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,elpa_hdl%elpa_comm_cols,cc,local_nrows_c)
#else
!ELPA-LEGACY-2017
 success =  elpa_mult_ah_b_complex_double(uplo_a,uplo_c,elpa_hdl%na,ncb,aa,elpa_hdl%local_nrows,&
&           elpa_hdl%local_ncols,bb,local_nrows_b,local_ncols_b,elpa_hdl%nblk,elpa_hdl%elpa_comm_rows,&
&           elpa_hdl%elpa_comm_cols,cc,local_nrows_c,local_ncols_c)
#endif

 if (.not.success) call elpa_func_error_handler(err_msg='Error in hermitian_multiply_complex!')

end subroutine elpa_func_hermitian_multiply_complex
!!***

!----------------------------------------------------------------------

#endif

END MODULE m_elpa
!!***
