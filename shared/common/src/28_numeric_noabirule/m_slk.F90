!!****m* ABINIT/m_slk
!! NAME
!! m_slk
!!
!! FUNCTION
!!  High-level objects and wrappers around the ScaLAPACK and ELPA API.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2024 ABINIT group (CS,GZ,FB,MG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_slk

 use defs_basis
 use m_xmpi
 use m_errors
 use m_abicore
#ifdef HAVE_LINALG_ELPA
 use m_elpa
#endif
#ifdef HAVE_MPI2
 use mpi
#endif

 use m_fstrings,      only : firstchar, toupper, itoa, sjoin, ltoa
 use m_time,          only : cwtime, cwtime_report
 use m_numeric_tools, only : blocked_loop !, print_arr

 implicit none

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

 private

 ! scaLAPACK array descriptor.
 integer,private,parameter :: DLEN_ = 9    ! length
 integer,private,parameter :: Dtype_ = 1   ! type
 integer,private,parameter :: CTXT_ = 2    ! BLACS context
 integer,private,parameter :: M_ = 3       ! nb global lines
 integer,private,parameter :: N_ = 4       ! nb global columns
 integer,private,parameter :: MB_ = 5      ! nb lines of a block
 integer,private,parameter :: NB_ = 6      ! nb columns of a block
 integer,private,parameter :: RSRC_ = 7    ! line of processors at the beginning
 integer,private,parameter :: CSRC_ = 8    ! column of processors at the beginning
 integer,private,parameter :: LLD_ = 9     ! local number of lines
!!***

!----------------------------------------------------------------------

!!****t* m_slk/grid_scalapack
!! NAME
!!  grid_scalapack
!!
!! FUNCTION
!!  Grid of ScaLAPACK processors.
!!
!! SOURCE

 type,public :: grid_scalapack

   integer :: nbprocs = -1
   ! Total number of processors

   integer :: dims(2) = -1
   ! Number of procs for rows/columns

   integer :: ictxt = xmpi_comm_null
   ! BLACS context i.e. MPI communicator.

   logical :: use_gpu = .false.
   ! Wether GPU is used, relevant for determining matrix block size.

 contains
   procedure :: init =>  grid_init  ! Set up the processor grid for ScaLAPACK.
 end type grid_scalapack
!!***

!----------------------------------------------------------------------

!!****t* m_slk/processor_scalapack
!! NAME
!!  processor_scalapack
!!
!! FUNCTION
!!  One processor in the grid.
!!
!! SOURCE

 type,public :: processor_scalapack

   integer :: myproc = -1
   ! rank the processor

   integer :: comm = xmpi_comm_null
   ! MPI communicator underlying the BLACS grid.

   integer :: coords(2) = -1
   ! Coordinates of the processor in the grid.

   type(grid_scalapack) :: grid
   ! the grid to which the processor is associated to.

 contains
   procedure :: init => processor_init        ! Initializes an instance of processor ScaLAPACK from an MPI communicator.
   procedure :: free => processor_free        ! Free the object
 end type processor_scalapack

 private :: build_processor_scalapack     ! Build a ScaLAPACK processor descriptor
!!***

!----------------------------------------------------------------------

!!****t* m_slk/descript_scalapack
!! NAME
!!  descript_scalapack
!!
!! FUNCTION
!! ScaLAPACK matrix descriptor.
!!
!! SOURCE

 type,public :: descript_scalapack
   integer :: tab(DLEN_)
 end type descript_scalapack
!!***

!----------------------------------------------------------------------

!!****t* m_slk/basemat_t
!! NAME
!!  basemat_t
!!
!! FUNCTION
!!  Base class for scalapack matrices
!!
!! SOURCE

 type,private :: basemat_t

   integer :: sizeb_local(2) = -1
    ! dimensions of the local buffer

   integer :: sizeb_global(2) = -1
    ! dimensions of the global matrix

   integer :: sizeb_blocs(2) = -1
    ! size of the block of consecutive data

   integer :: istwf_k = -1

   type(processor_scalapack),pointer :: processor => null()

   type(descript_scalapack) :: descript

 contains

   procedure :: init => init_matrix_scalapack
    ! Constructor

   procedure :: glob2loc => slk_glob2loc
    ! Determine the local indices of an element from its global indices and return haveit bool flag.

   procedure :: loc2glob => slk_matrix_loc2glob
    ! Return global indices of a matrix element from the local indices.

   procedure :: loc2grow => slk_matrix_loc2grow
    ! Determine the global row index from the local index

   procedure :: loc2gcol => slk_matrix_loc2gcol
    ! Determine the global column index from the local index

   procedure :: locmem_mb => locmem_mb
    ! Return memory allocated for the local buffer in Mb.

   procedure :: print => slkmat_print
    ! Print info on the object.

   procedure :: check_local_shape => slkmat_check_local_shape
   !  Debugging tool to test the local shape `lshape` of the local buffer.

   procedure :: free => matrix_scalapack_free
    ! Free memory

   procedure :: change_size_blocs => slk_change_size_blocs
    ! Change the block sizes, return new object.

   procedure :: get_trace => slk_get_trace
    ! Compute the trace of an N-by-N distributed matrix.

   procedure :: set_imag_diago_to_zero => slk_set_imag_diago_to_zero
    ! Set the imaginary part of the diagonal to zero.

   procedure :: invert => slk_invert
    ! Inverse of a complex matrix.

 end type basemat_t
!!***

!----------------------------------------------------------------------

!!****t* m_slk/matrix_scalapack
!! NAME
!!  matrix_scalapack
!!
!! FUNCTION
!!  high-level interface to ScaLAPACK matrix (double precision version)
!!
!! SOURCE

 type, public, extends(basemat_t) :: matrix_scalapack

   real(dp),allocatable :: buffer_real(:,:)
    ! local part of the (real) matrix.
    ! The istwf_k option passed to the constructor defines whether we have a real or complex matrix

   complex(dpc),allocatable :: buffer_cplx(:,:)
    ! local part of the (complex) matrix

 contains

   procedure :: get_head_and_wings => slk_get_head_and_wings
    ! Return global arrays with the head and the wings of the matrix.

   procedure :: set_head_and_wings => slk_set_head_and_wings
    ! Set head and the wings of the matrix starting from global arrays.

   procedure :: copy => matrix_scalapack_copy
    ! Copy object

   procedure :: hpd_invert => slk_hpd_invert
    ! Inverse of a Hermitian positive definite matrix.

   procedure :: ptrans => slk_ptrans
    ! Transpose matrix

   procedure :: cut => slk_cut
    ! Extract submatrix and create new matrix with `size_blocs` and `processor`

   procedure :: take_from => slk_take_from
    ! Take values from source

   procedure :: collect_cplx => slk_collect_cplx
    ! Return on all processors the submatrix of shape (mm, nn) starting at position ija.

   procedure :: heev => slk_heev
    ! Compute eigenvalues and, optionally, eigenvectors of an Hermitian matrix A. A * X = lambda * X

   procedure :: pzheevx => slk_pzheevx
    ! Compute Eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A. ! A * X = lambda *  X

   procedure :: pzhegvx => slk_pzhegvx
    ! Eigenvalues and, optionally, eigenvectors of a complex
    ! generalized Hermitian-definite eigenproblem, of the form
    ! sub( A )*x=(lambda)*sub( B )*x,  sub( A )*sub( B )x=(lambda)*x,
    ! or sub( B )*sub( A )*x=(lambda)*x.

   procedure :: symmetrize => slk_symmetrize
    ! Symmetrizes a square scaLAPACK matrix.

   procedure :: bsize_and_type => slk_bsize_and_type
    ! Returns the byte size and the MPI datatype

 end type matrix_scalapack
!!***

!!****t* m_slk/slkmat_sp_t
!! NAME
!!  slkmat_sp_t
!!
!! FUNCTION
!!  high-level interface to ScaLAPACK matrix (single precision version)
!!
!! SOURCE

 type, public, extends(basemat_t) :: slkmat_sp_t

   real(sp),allocatable :: buffer_real(:,:)
    ! local part of the (real) matrix.
    ! The istwf_k option passed to the constructor defines whether we have a real or complex matrix

   complex(sp),allocatable :: buffer_cplx(:,:)
    ! local part of the (complex) matrix

 contains

   procedure :: copy => slkmat_sp_copy
    ! Copy object

   procedure :: take_from => slkmat_sp_take_from
    ! Take values from source

   procedure :: ptrans => slkmat_sp_ptrans
    ! Transpose matrix

   procedure :: set_head_and_wings => slkmat_sp_set_head_and_wings
    ! Set head and the wings of the matrix starting from global arrays.

   !procedure :: cut => slk_cut
    ! Extract submatrix and create new matrix with `size_blocs` and `processor`

   procedure :: collect_cplx => slkmat_sp_collect_cplx
    ! Return on all processors the submatrix of shape (mm, nn) starting at position ija.

   procedure :: heev => slkmat_sp_heev
    ! Compute eigenvalues and, optionally, eigenvectors of an Hermitian matrix A. A * X = lambda * X

   procedure :: hpd_invert => slkmat_sp_hpd_invert
    ! Inverse of a Hermitian positive definite matrix.

 end type slkmat_sp_t
!!***

 public :: block_dist_1d                   ! Return block size for one-dimensional block column/row distribution
 public :: slk_has_elpa                    ! Return True if ELPA support is activated

 public :: matrix_get_local_cplx           ! Return a local matrix coefficient of complex type.
 public :: matrix_get_local_real           ! Return a local matrix coefficient of double precision type.
 public :: matrix_set_local_cplx           ! Sets a local matrix coefficient of complex type.
 public :: matrix_set_local_real           ! Sets a local matrix coefficient of double precision type.
 public :: idx_loc                         ! Local indices of an entry
                                           ! from its global indices, independently of the processor.
 !public :: glob_loc__                     ! Return global location of a matrix coefficient.
 !public :: loc_glob__                     ! Return global index from a local index (row or column)
                                           ! as a function of a given processor
 public :: matrix_from_global              ! Fills SCALAPACK matrix from full matrix.
 public :: matrix_from_global_sym          ! Fills SCALAPACK matrix from a full matrix.
 public :: matrix_from_realmatrix          ! Fills SCALAPACK matrix from full matrix.
 public :: matrix_from_complexmatrix       ! Fills SCALAPACK matrix from a full matrix.
 public :: matrix_to_global                ! Inserts a ScaLAPACK matrix into a global one.
 public :: matrix_to_realmatrix            ! Inserts a ScaLAPACK matrix into a real matrix.
 public :: matrix_to_complexmatrix         ! Inserts a ScaLAPACK matrix into a complex matrix.
 public :: matrix_to_reference             ! Fill a full matrix with respect to a SCALAPACK matrix.
 public :: slk_matrix_from_global_dpc_2D   ! Fill a complex SCALAPACK matrix with respect to a global matrix.
 public :: slk_matrix_from_global_dpc_1Dp  ! Fill a complex SCALAPACK matrix with respect to a global matrix.
                                           ! target: double precision complex matrix in packed form.
 public :: slk_matrix_to_global_dpc_2D     ! Fill a global matrix with respect to a SCALAPACK matrix.
                                           ! target: Two-dimensional Double precision complex matrix.

 public :: slk_pgemm                         ! Compute: C := alpha*A*B + beta*C
 interface slk_pgemm
   module procedure slk_pgemm_dp
   module procedure slk_pgemm_sp
 end interface slk_pgemm


 public :: compute_eigen_problem             ! Compute eigenvalues and eigenvectors of: A * X = lambda * X.
                                             ! complex and real cases.
 public :: compute_generalized_eigen_problem ! Compute generalized eigenvalue problem
 public :: compute_eigen1                    ! Compute eigenvalues and eigenvectors.  complex and real cases.
 public :: compute_eigen2                    ! Compute eigenvalues and eigenvectors: A * X = lambda * B * X
                                             ! complex and real cases.

 public :: slk_write                         ! Writes a square scaLAPACK distributed matrix on an external file using MPI-IO.
 public :: slk_read                          ! Read a square scaLAPACK distributed matrix from an external file using MPI-IO.
 public :: slk_single_fview_read_mask        ! Returns an MPI datatype that can be used to read a scaLAPACK matrix from
                                             ! a binary file using MPI-IO.
                                             ! The view is created using the user-defined mask function
 public :: slk_single_fview_read             ! Returns an MPI datatype to read a scaLAPACK distributed matrix
                                             ! from a binary file using MPI-IO.
 public :: slk_single_fview_write            ! Returns an MPI datatype to write a scaLAPACK distributed matrix
                                             ! to a binary file using MPI-IO.


 public :: slk_array_free                    !  Deallocate array of matrix_scalapack elements
 interface slk_array_free
   module procedure slk_array1_free
   module procedure slk_array2_free
   module procedure slk_array3_free
   module procedure slk_array4_free
 end interface slk_array_free

 public :: slk_array_set                       ! Elemental routine to set the value of the buffer to a costant value `cvalue`.
 public :: slk_array_locmem_mb                 ! Compute memory allocated for an array of matrix_scalapack elements

CONTAINS  !==============================================================================
!!***

!!****f* m_slk/grid_init
!! NAME
!!  grid_init
!!
!! FUNCTION
!!  Set up the ScaLAPACKgrid given the total number of processors.
!!
!! INPUTS
!!  nbprocs= total number of processors
!!  comm= MPI communicator
!!  [grid_dims]=Number of procs for each dimension.
!!
!! OUTPUT
!!  grid= the grid of processors used by Scalapack
!!
!! SOURCE

subroutine grid_init(grid, nbprocs, comm, use_gpu, grid_dims)

!Arguments ------------------------------------
 class(grid_scalapack),intent(out) :: grid
 integer,intent(in) :: nbprocs,comm
 logical,intent(in) :: use_gpu
 integer,optional,intent(in) :: grid_dims(2)

!Local variables-------------------------------
 integer :: i

! *********************************************************************

 grid%nbprocs = nbprocs

 if (.not. present(grid_dims)) then
   ! Search for a rectangular grid of processors
   i=INT(SQRT(float(nbprocs)))
   do while (MOD(nbprocs,i) /= 0)
     i = i-1
   end do
   i=max(i,1)

   grid%dims(1) = i
   grid%dims(2) = INT(nbprocs/i)

 else
   grid%dims = grid_dims
 end if

 ABI_CHECK(product(grid%dims) == nbprocs, sjoin("grid%dims:", ltoa(grid%dims), "does not agree with nprocs:", itoa(nbprocs)))

 grid%ictxt = comm
 grid%use_gpu = use_gpu

#ifdef HAVE_LINALG_SCALAPACK
 ! 'R': Use row-major natural ordering
 call BLACS_GRIDINIT(grid%ictxt, 'R', grid%dims(1), grid%dims(2))
#endif

end subroutine grid_init
!!***

!----------------------------------------------------------------------

!!****f* m_slk/build_processor_scalapack
!! NAME
!!  build_processor_scalapack
!!
!! FUNCTION
!!  Builds a ScaLAPACK processor descriptor.
!!  Build of the data related to one processor in a grid
!!
!! INPUTS
!!  grid= array representing the grid of processors.
!!  myproc= selected processor
!!  comm= MPI communicator
!!
!! OUTPUT
!!  processor= descriptor of a processor
!!
!! SOURCE

subroutine build_processor_scalapack(processor, grid, myproc, comm)

!Arguments ------------------------------------
 class(processor_scalapack),intent(inout) :: processor
 integer,intent(in) :: myproc,comm
 class(grid_scalapack),intent(in) :: grid

! *********************************************************************

 processor%grid = grid
 processor%myproc = myproc
 processor%comm = comm

#ifdef HAVE_LINALG_SCALAPACK
 call BLACS_GRIDINFO(grid%ictxt, processor%grid%dims(1), processor%grid%dims(2), &
                     processor%coords(1), processor%coords(2))
#endif

 ! These values are the same as those computed by BLACS_GRIDINFO
 ! except in the case where the myproc argument is not the local proc
 processor%coords(1) = INT((myproc) / grid%dims(2))
 processor%coords(2) = MOD((myproc), grid%dims(2))

end subroutine build_processor_scalapack
!!***

!----------------------------------------------------------------------

!!****f* m_slk/processor_init
!! NAME
!!  processor_init
!!
!! FUNCTION
!!  Initializes an instance of processor ScaLAPACK from an MPI communicator.
!!
!! INPUTS
!!  comm= MPI communicator
!!  [grid_dims]=Number of procs for each dimension.
!!
!! OUTPUT
!!  processor= descriptor of a processor
!!
!! SOURCE

subroutine processor_init(processor, comm, grid_dims)

!Arguments ------------------------------------
 class(processor_scalapack),intent(out) :: processor
 integer, intent(in) :: comm
 integer,optional,intent(in) :: grid_dims(2)

!Local variables-------------------------------
 type(grid_scalapack) :: grid
 integer :: nbproc,myproc

! *********************************************************************

 nbproc = xmpi_comm_size(comm)
 myproc = xmpi_comm_rank(comm)

 if (present(grid_dims)) then
   call grid%init(nbproc, comm, .false., grid_dims=grid_dims)
 else
   call grid%init(nbproc, comm, .false.)
 end if

 call build_processor_scalapack(processor, grid, myproc, comm)

end subroutine processor_init
!!***

!----------------------------------------------------------------------

!!****f* m_slk/processor_free
!! NAME
!!  processor_free
!!
!! FUNCTION
!!  Removes a processor from the ScaLAPACK grid.
!!
!! SOURCE

subroutine processor_free(processor)

!Arguments ------------------------------------
 class(processor_scalapack),intent(inout) :: processor

! *********************************************************************

#ifdef HAVE_LINALG_SCALAPACK
 if (processor%grid%ictxt /= xmpi_comm_null) then
   call BLACS_GRIDEXIT(processor%grid%ictxt)
   !call BLACS_EXIT(0)
 end if
#endif

end subroutine processor_free
!!***

!----------------------------------------------------------------------

!!****f* m_slk/init_matrix_scalapack
!! NAME
!!  init_matrix_scalapack
!!
!! FUNCTION
!!  Initialisation of a SCALAPACK matrix (each proc initializes its own part of the matrix)
!!
!! INPUTS
!!  processor= descriptor of a processor
!!  nbli_global= total number of lines
!!  nbco_global= total number of columns
!!  istwf_k= 2 if we have a real matrix else complex.
!!  [size_blocs]= custom block sizes. Use -1 to use global size along that direction.
!!    Useful to distribute only rows or columns. Obviously, [-1, -1] is not allowed.
!!
!! OUTPUT
!!  matrix= the matrix to process
!!
!! SOURCE

subroutine init_matrix_scalapack(matrix, nbli_global, nbco_global, processor, istwf_k, size_blocs)

!Arguments ------------------------------------
 class(basemat_t),intent(inout) :: matrix
 integer,intent(in) :: nbli_global, nbco_global, istwf_k
 type(processor_scalapack),target,intent(in) :: processor
 integer,intent(in),optional :: size_blocs(2)

#ifdef HAVE_LINALG_SCALAPACK
!Local variables-------------------------------
#ifdef HAVE_LINALG_ELPA
 integer, parameter :: DEFAULT_SIZE_BLOCS = 1
#else
 ! As recommended by Intel MKL, a more sensible default than the previous value of 40
 integer, parameter :: DEFAULT_SIZE_BLOCS = 24
#endif
 ! As recommanded in ELPA, which advises distributions as squared as possible using powers of 2
 integer, parameter :: DEFAULT_SIZE_BLOCS_GPU = 16
 integer :: info,sizeb
 integer,external :: NUMROC
 !character(len=500) :: msg

! *********************************************************************

 !call matrix%free()

 sizeb = DEFAULT_SIZE_BLOCS
 if(processor%grid%use_gpu) sizeb = DEFAULT_SIZE_BLOCS_GPU

 !Records of the matrix type
 matrix%processor => processor
 matrix%sizeb_blocs(1) = MIN(sizeb, nbli_global)
 matrix%sizeb_blocs(2) = MIN(sizeb, nbco_global)

#ifdef HAVE_LINALG_ELPA
 if(matrix%sizeb_blocs(1) .ne. matrix%sizeb_blocs(2)) then
    matrix%sizeb_blocs(1) = MIN(matrix%sizeb_blocs(1), matrix%sizeb_blocs(2))
    matrix%sizeb_blocs(2) = matrix%sizeb_blocs(1)
 end if
#endif

 ! Use custom block sizes.
 if (present(size_blocs)) then
   ABI_CHECK(.not. all(size_blocs == -1), "size_blocs [-1, -1]  is not allowed")
   if (size_blocs(1) == -1) then
     matrix%sizeb_blocs(1) = nbli_global
   else
     matrix%sizeb_blocs(1) = MIN(size_blocs(1), nbli_global)
   end if
   if (size_blocs(2) == -1) then
     matrix%sizeb_blocs(2) = nbco_global
   else
     matrix%sizeb_blocs(2) = MIN(size_blocs(2), nbco_global)
   end if
 end if

 matrix%sizeb_global(1) = nbli_global
 matrix%sizeb_global(2) = nbco_global

 ! Size of the local buffer
 ! NUMROC computes the NUMber of Rows Or Columns of a distributed matrix owned by the process indicated by IPROC.
 ! NUMROC (n, nb, iproc, isrcproc, nprocs)
 matrix%sizeb_local(1) = NUMROC(nbli_global, matrix%sizeb_blocs(1), &
                                processor%coords(1), 0, processor%grid%dims(1))

 matrix%sizeb_local(2) = NUMROC(nbco_global,matrix%sizeb_blocs(2), &
                                processor%coords(2), 0, processor%grid%dims(2))

 call idx_loc(matrix, matrix%sizeb_global(1), matrix%sizeb_global(2), &
              matrix%sizeb_local(1), matrix%sizeb_local(2))

 ! Initialisation of the SCALAPACK description of the matrix
 ! (desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info)
 call DESCINIT(matrix%descript%tab, nbli_global, nbco_global, &
               matrix%sizeb_blocs(1), matrix%sizeb_blocs(2), 0, 0, &
               processor%grid%ictxt, MAX(1, matrix%sizeb_local(1)), info)

 if (info /= 0) then
   ABI_ERROR(sjoin("Error while initializing scalapack matrix. info:", itoa(info)))
 end if

 ! Allocate local buffer.
 matrix%istwf_k = istwf_k
 select type (matrix)
 class is (matrix_scalapack)
   if (istwf_k /= 2) then
     ABI_MALLOC(matrix%buffer_cplx, (matrix%sizeb_local(1), matrix%sizeb_local(2)))
     matrix%buffer_cplx = czero
   else
     ABI_MALLOC(matrix%buffer_real, (matrix%sizeb_local(1), matrix%sizeb_local(2)))
     matrix%buffer_real = zero
   end if
 class is (slkmat_sp_t)
   if (istwf_k /= 2) then
     ABI_MALLOC(matrix%buffer_cplx, (matrix%sizeb_local(1), matrix%sizeb_local(2)))
     matrix%buffer_cplx = czero_sp
   else
     ABI_MALLOC(matrix%buffer_real, (matrix%sizeb_local(1), matrix%sizeb_local(2)))
     matrix%buffer_real = zero_sp
   end if
 class default
   ABI_ERROR("Wrong class")
 end select
#endif

end subroutine init_matrix_scalapack
!!***

!----------------------------------------------------------------------

!!****f* m_slk/locmem_mb
!! NAME
!!  locmem_mb
!!
!! FUNCTION
!!  Returns memory allocated for the local buffer in Mb.
!!
!! SOURCE

pure real(dp) function locmem_mb(mat)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: mat

! *********************************************************************

 locmem_mb = zero
 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_real)) locmem_mb = product(int(shape(mat%buffer_real))) * dp
   if (allocated(mat%buffer_cplx)) locmem_mb = product(int(shape(mat%buffer_cplx))) * two * dp
 class is (slkmat_sp_t)
   if (allocated(mat%buffer_real)) locmem_mb = product(int(shape(mat%buffer_real))) * sp
   if (allocated(mat%buffer_cplx)) locmem_mb = product(int(shape(mat%buffer_cplx))) * two * sp
 end select
 locmem_mb = locmem_mb * b2Mb

end function locmem_mb
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_print
!! NAME
!!  slkmat_print
!!
!! FUNCTION
!!  Print info on scalapack matrix.
!!
!! INPUTS
!!  [unit]=Unit number (default: std_out)
!!  [header]=title for info
!!  [prtvol]=Verbosity level (default: 0)
!!
!! SOURCE

subroutine slkmat_print(mat, header, unit, prtvol)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: mat
 character(len=*),optional,intent(in) :: header
 integer,optional,intent(in) :: prtvol, unit

!Local variables-------------------------------
 integer :: unt, my_prtvol, grid_dims(2)
 character(len=50) :: matrix_dtype
 character(len=5000) :: msg

! *********************************************************************

 unt = std_out; if (present(unit)) unt =unit
 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 msg = ' ==== Info on matrix_scalapack ==== '
 if (present(header)) msg=' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(unt, msg)

 matrix_dtype = "undefined"
 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_real)) matrix_dtype = "real"
   if (allocated(mat%buffer_cplx)) matrix_dtype = "complex"
 class is (slkmat_sp_t)
   if (allocated(mat%buffer_real)) matrix_dtype = "real"
   if (allocated(mat%buffer_cplx)) matrix_dtype = "complex"
 class default
   ABI_ERROR("Wrong class")
 end select

 grid_dims = [-1, -1]
 if (associated(mat%processor)) grid_dims = mat%processor%grid%dims

 write(msg,'(5(3a),a,f8.1,a)') &
   "  matrix_dtype ..... ", trim(matrix_dtype), ch10, &
   "  sizeb_global ..... ", trim(ltoa(mat%sizeb_global)), ch10, &
   "  sizeb_local ...... ", trim(ltoa(mat%sizeb_local)), ch10, &
   "  sizeb_blocs ...... ", trim(ltoa(mat%sizeb_blocs)), ch10, &
   "  processor grid  .. ", trim(ltoa(grid_dims)), ch10, &
   "  memory (Mb) ...... ", mat%locmem_mb(), ch10
 call wrtout(unt, msg)

 !if (prtvol > 10) call mat%write(unit)

end subroutine slkmat_print
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_check_shape
!! NAME
!!  slkmat_check_shape
!!
!! FUNCTION
!!  Debugging tool to test the local shape `lshape` of the local buffer.
!!  Return exit status in `ok` and error message in `msg`.
!!
!! SOURCE

logical function slkmat_check_local_shape(mat, lshape, msg) result (ok)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: mat
 integer,intent(in) :: lshape(2)
 character(len=*),intent(out) :: msg

! *********************************************************************

 msg = ""
 ok = all(mat%sizeb_local == lshape)
 if (.not. ok) then
   msg = sjoin("mat%sizeb_local:", ltoa(mat%sizeb_local), " not equal to input local lshape ", ltoa(lshape))
   return
 end if

 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_cplx)) then
     ok = all(shape(mat%buffer_cplx) == lshape)
     if (.not. ok) then
       msg = sjoin("shape(buffer_cplx):", ltoa(shape(mat%buffer_cplx)), " != input local lshape ", ltoa(lshape)); return
     end if
   else if (allocated(mat%buffer_real)) then
     ok = all(shape(mat%buffer_real) == lshape)
     if (.not. ok) then
       msg = sjoin("shape(buffer_real):", ltoa(shape(mat%buffer_real)), " != input local lshape ", ltoa(lshape)); return
     end if
   end if

 class is (slkmat_sp_t)
   ! Same piece of code as above. May use include file!
   if (allocated(mat%buffer_cplx)) then
     ok = all(shape(mat%buffer_cplx) == lshape)
     if (.not. ok) then
       msg = sjoin("shape(buffer_cplx):", ltoa(shape(mat%buffer_cplx)), " != input local lshape ", ltoa(lshape)); return
     end if
   else if (allocated(mat%buffer_real)) then
     ok = all(shape(mat%buffer_real) == lshape)
     if (.not. ok) then
       msg = sjoin("shape(buffer_real):", ltoa(shape(mat%buffer_real)), " != input local lshape ", ltoa(lshape)); return
     end if
   end if

 class default
   ABI_ERROR("Wrong class")
 end select

end function slkmat_check_local_shape
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_get_head_and_wings
!! NAME
!!  slk_get_head_and_wings
!!
!! FUNCTION
!!  Return global arrays with head and wings of the matrix.
!!  If call_mpi if False, global MPI sum is postponed.
!!  Useful to reduce the number of MPI calls if one has to operate on multiple matrices.
!!
!! SOURCE

subroutine slk_get_head_and_wings(mat, head, low_wing, up_wing, call_mpi)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: mat
 complex(dp),intent(out) :: head, low_wing(mat%sizeb_global(1)), up_wing(mat%sizeb_global(2))
 logical,intent(in) :: call_mpi

!Local variables-------------------------------
 integer :: ierr, il_g1, il_g2, iglob1, iglob2
 logical :: is_cplx

! *********************************************************************

 head = zero; low_wing = zero; up_wing = zero

 is_cplx = allocated(mat%buffer_cplx)

 do il_g2=1,mat%sizeb_local(2)
   iglob2 = mat%loc2gcol(il_g2)
   do il_g1=1,mat%sizeb_local(1)
     iglob1 = mat%loc2grow(il_g1)

     if (iglob1 == 1 .or. iglob2 == 1) then
       if (iglob1 == 1 .and. iglob2 == 1) then
         if (is_cplx) then
           head = mat%buffer_cplx(il_g1, il_g2)
         else
           head = mat%buffer_real(il_g1, il_g2)
         end if
       else if (iglob1 == 1) then
         if (is_cplx) then
           up_wing(iglob2) = mat%buffer_cplx(il_g1, il_g2)
         else
           up_wing(iglob2) = mat%buffer_real(il_g1, il_g2)
         end if
       else if (iglob2 == 1) then
         if (is_cplx) then
           low_wing(iglob1) = mat%buffer_cplx(il_g1, il_g2)
         else
           low_wing(iglob1) = mat%buffer_real(il_g1, il_g2)
         end if
       end if
     end if

   end do
 end do

 if (call_mpi) then
   call xmpi_sum(head, mat%processor%grid%ictxt, ierr)
   call xmpi_sum(low_wing, mat%processor%grid%ictxt, ierr)
   call xmpi_sum(up_wing, mat%processor%grid%ictxt, ierr)
 end if

end subroutine slk_get_head_and_wings
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_set_head_and_wings
!! NAME
!!  slk_set_head_and_wings
!!
!! FUNCTION
!!  Set head and the wings of the matrix starting from global arrays.
!!
!! SOURCE

subroutine slk_set_head_and_wings(mat, head, low_wing, up_wing)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: mat
 complex(dp),intent(in) :: head, low_wing(mat%sizeb_global(1)), up_wing(mat%sizeb_global(2))

!Local variables-------------------------------
 integer :: il_g1, il_g2, iglob1, iglob2
 logical :: is_cplx

! *********************************************************************

 is_cplx = allocated(mat%buffer_cplx)

 do il_g2=1,mat%sizeb_local(2)
   iglob2 = mat%loc2gcol(il_g2)
   do il_g1=1,mat%sizeb_local(1)
     iglob1 = mat%loc2grow(il_g1)

     if (iglob1 == 1 .or. iglob2 == 1) then
       if (iglob1 == 1 .and. iglob2 == 1) then
         if (is_cplx) then
           mat%buffer_cplx(il_g1, il_g2) = head
         else
           mat%buffer_real(il_g1, il_g2) = real(head)
         end if
       else if (iglob1 == 1) then
         if (is_cplx) then
           mat%buffer_cplx(il_g1, il_g2) = up_wing(iglob2)
         else
           mat%buffer_real(il_g1, il_g2) = real(up_wing(iglob2))
         end if
       else if (iglob2 == 1) then
         if (is_cplx) then
           mat%buffer_cplx(il_g1, il_g2) = low_wing(iglob1)
         else
           mat%buffer_real(il_g1, il_g2) = real(low_wing(iglob1))
         end if
       end if
     end if

   end do
 end do

end subroutine slk_set_head_and_wings
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_sp_set_head_and_wings
!! NAME
!!  slkmat_sp_set_head_and_wings
!!
!! FUNCTION
!!  Set head and the wings of the matrix starting from global arrays.
!!
!! SOURCE

subroutine slkmat_sp_set_head_and_wings(mat, head, low_wing, up_wing)

!Arguments ------------------------------------
 class(slkmat_sp_t),intent(inout) :: mat
 complex(sp),intent(in) :: head, low_wing(mat%sizeb_global(1)), up_wing(mat%sizeb_global(2))

!Local variables-------------------------------
 integer :: il_g1, il_g2, iglob1, iglob2
 logical :: is_cplx

! *********************************************************************

 is_cplx = allocated(mat%buffer_cplx)

 do il_g2=1,mat%sizeb_local(2)
   iglob2 = mat%loc2gcol(il_g2)
   do il_g1=1,mat%sizeb_local(1)
     iglob1 = mat%loc2grow(il_g1)

     if (iglob1 == 1 .or. iglob2 == 1) then
       if (iglob1 == 1 .and. iglob2 == 1) then
         if (is_cplx) then
           mat%buffer_cplx(il_g1, il_g2) = head
         else
           mat%buffer_real(il_g1, il_g2) = real(head)
         end if
       else if (iglob1 == 1) then
         if (is_cplx) then
           mat%buffer_cplx(il_g1, il_g2) = up_wing(iglob2)
         else
           mat%buffer_real(il_g1, il_g2) = real(up_wing(iglob2))
         end if
       else if (iglob2 == 1) then
         if (is_cplx) then
           mat%buffer_cplx(il_g1, il_g2) = low_wing(iglob1)
         else
           mat%buffer_real(il_g1, il_g2) = real(low_wing(iglob1))
         end if
       end if
     end if

   end do
 end do

end subroutine slkmat_sp_set_head_and_wings
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_scalapack_copy
!! NAME
!!  matrix_scalapack_copy
!!
!! FUNCTION
!!  Copy in_mat to out_mat. If empty is True, the values in the local buffer are not copied. Default: False
!!
!! SOURCE

subroutine matrix_scalapack_copy(in_mat, out_mat, empty)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: in_mat
 class(matrix_scalapack),intent(out) :: out_mat
 logical,optional,intent(in) :: empty

!Local variables-------------------------------
 logical :: empty__

! *********************************************************************

 call out_mat%init(in_mat%sizeb_global(1), in_mat%sizeb_global(2), in_mat%processor, in_mat%istwf_k, &
                   size_blocs=in_mat%sizeb_blocs)

 empty__ = .False.; if (present(empty)) empty__ = empty
 if (.not. empty__) then
   if (in_mat%istwf_k == 1) then
     out_mat%buffer_cplx = in_mat%buffer_cplx
   else
     out_mat%buffer_real = in_mat%buffer_real
   end if
 end if

end subroutine matrix_scalapack_copy
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_sp_copy
!! NAME
!!  slkmat_sp_copy
!!
!! FUNCTION
!!  Copy in_mat to out_mat. If empty is True, the values in the local buffer are not copied. Default: False
!!
!! SOURCE

subroutine slkmat_sp_copy(in_mat, out_mat, empty)

!Arguments ------------------------------------
 class(slkmat_sp_t),intent(in) :: in_mat
 class(slkmat_sp_t),intent(out) :: out_mat
 logical,optional,intent(in) :: empty

!Local variables-------------------------------
 logical :: empty__

! *********************************************************************

 call out_mat%init(in_mat%sizeb_global(1), in_mat%sizeb_global(2), in_mat%processor, in_mat%istwf_k, &
                   size_blocs=in_mat%sizeb_blocs)

 empty__ = .False.; if (present(empty)) empty__ = empty
 if (.not. empty__) then
   if (in_mat%istwf_k == 1) then
     out_mat%buffer_cplx = in_mat%buffer_cplx
   else
     out_mat%buffer_real = in_mat%buffer_real
   end if
 end if

end subroutine slkmat_sp_copy
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_scalapack_free
!! NAME
!!  matrix_scalapack_free
!!
!! FUNCTION
!!  Free dynamic memory
!!
!! SOURCE

subroutine matrix_scalapack_free(mat)

!Arguments ------------------------------------
 class(basemat_t),intent(inout) :: mat

! *********************************************************************

 ! Don't free the grid. Just nullify the pointer as there might be other objects keeping a ref to processor.
 mat%processor => null()

 mat%sizeb_global = 0
 mat%sizeb_blocs = 0
 mat%sizeb_local = 0
 mat%descript%tab = 0

 select type (mat)
 class is (matrix_scalapack)
   ABI_SFREE(mat%buffer_cplx)
   ABI_SFREE(mat%buffer_real)
 class is (slkmat_sp_t)
   ABI_SFREE(mat%buffer_cplx)
   ABI_SFREE(mat%buffer_real)
 class default
   ABI_ERROR("Wrong class")
 end select

end subroutine matrix_scalapack_free
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_array1_free
!! NAME
!!  slk_array1_free
!!
!! FUNCTION
!!  Deallocate 1d array of matrix_scalapack elements
!!
!! SOURCE

subroutine slk_array1_free(slk_arr1)
  class(basemat_t),intent(inout) :: slk_arr1(:)
  integer :: i1
  do i1=1,size(slk_arr1, dim=1)
    call slk_arr1(i1)%free()
  end do
end subroutine slk_array1_free
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_array2_free
!! NAME
!!  slk_array2_free
!!
!! FUNCTION
!!  Deallocate 2d array of matrix_scalapack elements
!!
!! SOURCE

subroutine slk_array2_free(slk_arr2)
  class(basemat_t),intent(inout) :: slk_arr2(:,:)
  integer :: i1, i2
  do i2=1,size(slk_arr2, dim=2)
    do i1=1,size(slk_arr2, dim=1)
      call slk_arr2(i1, i2)%free()
    end do
  end do
end subroutine slk_array2_free
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_array3_free
!! NAME
!!  slk_array3_free
!!
!! FUNCTION
!!  Deallocate 3d array of matrix_scalapack elements
!!
!! SOURCE

subroutine slk_array3_free(slk_arr3)
  class(basemat_t),intent(inout) :: slk_arr3(:,:,:)
  integer :: i1, i2, i3
  do i3=1,size(slk_arr3, dim=3)
    do i2=1,size(slk_arr3, dim=2)
      do i1=1,size(slk_arr3, dim=1)
        call slk_arr3(i1, i2, i3)%free()
      end do
    end do
  end do
end subroutine slk_array3_free
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_array4_free
!! NAME
!!  slk_array4_free
!!
!! FUNCTION
!!  Deallocate 4d array of matrix_scalapack elements
!!
!! SOURCE

subroutine slk_array4_free(slk_arr4)
  class(basemat_t),intent(inout) :: slk_arr4(:,:,:,:)
  integer :: i1, i2, i3, i4
  do i4=1,size(slk_arr4, dim=4)
    do i3=1,size(slk_arr4, dim=3)
      do i2=1,size(slk_arr4, dim=2)
        do i1=1,size(slk_arr4, dim=1)
          call slk_arr4(i1, i2, i3, i4)%free()
        end do
      end do
    end do
  end do
end subroutine slk_array4_free
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_array_set
!! NAME
!!  slk_array_set
!!
!! FUNCTION
!!  Elemental routine to set the value of the PBLAS buffer to a costant value `cvalue`.
!!  Usually used to zero all the buffers in an array of matrix_scalapack objects.
!!
!! SOURCE

elemental subroutine slk_array_set(mat, cvalue)

!Arguments ------------------------------------
 class(basemat_t),intent(inout) :: mat
 complex(dp),intent(in) :: cvalue

 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_cplx)) mat%buffer_cplx = cvalue
   if (allocated(mat%buffer_real)) mat%buffer_real = real(cvalue, kind=dp)
 class is (slkmat_sp_t)
   if (allocated(mat%buffer_cplx)) mat%buffer_cplx = cmplx(cvalue, kind=sp)
   if (allocated(mat%buffer_real)) mat%buffer_real = real(cvalue, kind=sp)
 end select

end subroutine slk_array_set
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_array_locmem_mb
!! NAME
!!  slk_array_locmem_mb
!!
!! FUNCTION
!!  Elemental function to compute the memory allocated for an array of matrix_scalapack elements
!!  Usage: mem_mb = sum(mat_array)
!!
!! SOURCE

elemental real(dp) function slk_array_locmem_mb(mat) result(mem_mb)
  class(basemat_t),intent(in) :: mat
  mem_mb = mat%locmem_mb()
end function slk_array_locmem_mb
!!***

!----------------------------------------------------------------------

!!****f* m_slk/block_dist_1d
!! NAME
!!  block_dist_1d
!!
!! FUNCTION
!!  Return block size for one-dimensional block column (row) distribution.
!!  Mainly used to assign blocks of contiguous columns (rows) of a matrix to successive processes
!!  when a 1d grid is employed.
!!
!!  It is usually interfaced with CPP macros, e.g:
!!
!!    ABI_CHECK(block_dist_1d(mat_size, nproc, block_size, msg), msg)
!!
!! INPUTS
!!  mat_size=Size of the matrix (either number of rows or number of colums)
!!  nproc=Number of processoes in the 1D scalapack grid
!!
!! OUTPUT
!!  ok= Boolean flag with exit status (idle processes are not allowed).
!!  block_size=Size of the block along this axis needed for one-dimensional block distribution
!!  msg=Error message (if not ok)
!!
!! SOURCE

logical function block_dist_1d(mat_size, nproc, block_size, msg) result (ok)

!Arguments ------------------------------------
 integer, intent(in) :: mat_size, nproc
 integer,intent(out) :: block_size
 character(len=*),intent(out) :: msg

! *********************************************************************

 ok = .True.; msg = ""
 !block_size = 1; return

 block_size = mat_size / nproc
 if (block_size == 0) then
   ok = .False.
   write(msg, "(2(a,i0), 2a)") &
     "The number of MPI processors: ", nproc, " exceeeds the number of rows (columms) of the matrix: ", mat_size, ch10, &
     "Decrease the number of MPI processes for the scalapack level."
   return
 end if

 if (mod(mat_size, nproc) /= 0) block_size = block_size + 1

end function block_dist_1d
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_has_elpa
!! NAME
!!  slk_has_elpa
!!
!! FUNCTION
!!  Return True if ELPA support is activated
!!
!! SOURCE

pure logical function slk_has_elpa() result (ans)

! *********************************************************************

 ans = .False.
#ifdef HAVE_LINALG_ELPA
 ans = .True.
#endif

end function slk_has_elpa
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_get_local_cplx
!! NAME
!!  matrix_get_local_cplx
!!
!! FUNCTION
!!  Returns a local matrix coefficient of complex type.
!!  Access to a component thanks to its local indices
!!
!! INPUTS
!!  matrix= the matrix to process
!!  i= row in the matrix
!!  j= column in the matrix
!!
!! OUTPUT
!!  The value of the local matrix.
!!
!! SOURCE

pure complex(dpc) function matrix_get_local_cplx(matrix, i, j)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: matrix
 integer, intent(in) :: i,j

! *********************************************************************

 matrix_get_local_cplx = matrix%buffer_cplx(i,j)

end function matrix_get_local_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_get_local_real
!! NAME
!!  matrix_get_local_real
!!
!! FUNCTION
!!  Returns a local matrix coefficient of double precision type.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  i= row in the matrix
!!  j= column in the matrix
!!
!! SOURCE

pure real(dp) function matrix_get_local_real(matrix,i,j)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: matrix
 integer, intent(in) :: i,j

! *********************************************************************

 matrix_get_local_real = matrix%buffer_real(i,j)

end function matrix_get_local_real
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_set_local_cplx
!! NAME
!!  matrix_set_local_cplx
!!
!! FUNCTION
!!  Sets a local matrix coefficient of complex type.
!! -------------------------------------------------------
!!  Positioning of a component of a matrix thanks to its local indices
!! -------------------------------------------------------
!!
!! INPUTS
!!  i= row in the matrix
!!  j= column in the matrix
!!  value= the value to set
!!
!! SIDE EFFECTS
!!  matrix%buffer_cplx(i,j) filled with value
!!
!! SOURCE

pure subroutine matrix_set_local_cplx(matrix,i,j,value)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: matrix
 integer, intent(in) :: i,j
 complex(dp), intent(in) :: value

! *********************************************************************

 matrix%buffer_cplx(i,j) = value

end subroutine matrix_set_local_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_set_local_real
!! NAME
!!  matrix_set_local_real
!!
!! FUNCTION
!!  Sets a local matrix coefficient of double precision type.
!!
!! INPUTS
!!  i= row in the matrix
!!  j= column in the matrix
!!  value= the value to set
!!
!! SIDE EFFECTS
!!  matrix%buffer_real(i,j) set to value
!!
!! SOURCE

pure subroutine matrix_set_local_real(matrix, i, j, value)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: matrix
 integer, intent(in) :: i,j
 real(dp), intent(in) :: value

! *********************************************************************

 matrix%buffer_real(i,j) = value

end subroutine matrix_set_local_real
!!***

!----------------------------------------------------------------------

!!****f* m_slk/idx_loc
!! NAME
!!  idx_loc
!!
!! FUNCTION
!!  Return local indices from global indices, **independently** of the processor.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  i= row in the matrix
!!  j= column in the matrix
!!
!! OUTPUT
!!  iloc= local row of the coefficient
!!  jloc= local column of the coefficient
!!
!! SOURCE

subroutine idx_loc(matrix, i, j, iloc, jloc)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: matrix
 integer, intent(in) :: i,j
 integer, intent(out) :: iloc,jloc

! *********************************************************************

 iloc = glob_loc__(matrix, i, 1)
 jloc = glob_loc__(matrix, j, 2)

end subroutine idx_loc
!!***

!----------------------------------------------------------------------

!!****f* m_slk/glob_loc__
!! NAME
!!  glob_loc__
!!
!! FUNCTION
!!  Returns the global location of a matrix coefficient.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  idx= number of rows in the distributed matrix
!!  lico= block size index
!!
!! SOURCE

integer function glob_loc__(matrix, idx, lico)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: matrix
 integer, intent(in) :: idx, lico

#ifdef HAVE_LINALG_SCALAPACK
!Local variables-------------------------------
 integer,external :: NUMROC

! *********************************************************************

 glob_loc__ = NUMROC(idx, matrix%sizeb_blocs(lico), &
                   matrix%processor%coords(lico), 0, matrix%processor%grid%dims(lico))
#endif

end function glob_loc__
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_glob2loc
!! NAME
!!  slk_glob2loc
!!
!! FUNCTION
!!  Determine the local indices of an element from its global indices and return haveit bool flag.
!!
!! INPUTS
!!  iloc= local row index.
!!  jloc= local column index.
!!
!! OUTPUT
!!  iloc= row in the matrix
!!  jloc= column in the matrix
!!  haveit= True if (iglob, jglob) is stored on this proc
!!
!! SOURCE

subroutine slk_glob2loc(mat, iglob, jglob, iloc, jloc, haveit)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: mat
 integer, intent(in) :: iglob, jglob
 integer, intent(out) :: iloc, jloc
 logical,intent(out) :: haveit

!Local variables-------------------------------
 integer :: row_src, col_src

! *********************************************************************

#ifdef HAVE_LINALG_SCALAPACK
 ! SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW, MYCOL, LRINDX, LCINDX, RSRC, CSRC)

 call INFOG2L(iglob, jglob, mat%descript%tab, mat%processor%grid%dims(1), mat%processor%grid%dims(2), &
   mat%processor%coords(1), mat%processor%coords(2), iloc, jloc, row_src, col_src)

 haveit = all(mat%processor%coords == [row_src, col_src])
#endif

end subroutine slk_glob2loc
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_matrix_loc2glob
!! NAME
!!  slk_matrix_loc2glob
!!
!! FUNCTION
!!  Determine the global indices of an element from its local indices.
!!
!! INPUTS
!!  matrix= the matrix to process.
!!  iloc= local row index.
!!  jloc= local column index.
!!
!! OUTPUT
!!  i= row in the matrix
!!  j= column in the matrix
!!
!! SOURCE

pure subroutine slk_matrix_loc2glob(matrix, iloc, jloc, i, j)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: matrix
 integer, intent(in) :: iloc,jloc
 integer, intent(out) :: i,j

! *********************************************************************

 i = loc_glob__(matrix, matrix%processor, iloc, 1)
 j = loc_glob__(matrix, matrix%processor, jloc, 2)

end subroutine slk_matrix_loc2glob
!!***

!!****f* m_slk/slk_matrix_loc2grow
!! NAME
!!  slk_matrix_loc2grow
!!
!! FUNCTION
!!  Determine the global row index from the local index
!!
!! INPUTS
!!  matrix= the matrix to process.
!!  iloc= local row index.
!!
!! SOURCE

integer pure function slk_matrix_loc2grow(matrix, iloc) result(iglob)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: matrix
 integer, intent(in) :: iloc

! *********************************************************************

 iglob = loc_glob__(matrix, matrix%processor, iloc, 1)

end function slk_matrix_loc2grow
!!***

!!****f* m_slk/slk_matrix_loc2gcol
!! NAME
!!  slk_matrix_loc2gcol
!!
!! FUNCTION
!!  Determine the global column index of an element from the local index
!!
!! INPUTS
!!  matrix= the matrix to process.
!!  jloc= local column index.
!!
!! SOURCE

integer pure function slk_matrix_loc2gcol(matrix, jloc) result(jglob)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: matrix
 integer, intent(in) :: jloc

! *********************************************************************

 jglob = loc_glob__(matrix, matrix%processor, jloc, 2)

end function slk_matrix_loc2gcol
!!***

!----------------------------------------------------------------------

!!****f* m_slk/loc_glob__
!! NAME
!!  loc_glob__
!!
!! FUNCTION
!!  Determine the global index from a local index (row or column) as a function of a given processor
!!
!! INPUTS
!!  matrix= the matrix to process
!!  proc= descriptor of a processor
!!  idx= number of rows in the distributed matrix
!!  lico= block size index. 1 for rows. 2 for columns
!!
!! SOURCE

integer pure function loc_glob__(matrix, proc, idx, lico)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: matrix
 class(processor_scalapack),intent(in) :: proc
 integer, intent(in) :: idx,lico

!Local variables-------------------------------
 integer :: nbcyc, rest, nblocs

! *********************************************************************

 nbcyc = INT((idx-1) / matrix%sizeb_blocs(lico))
 rest = MOD(idx-1, matrix%sizeb_blocs(lico))
 nblocs = nbcyc * proc%grid%dims(lico) + proc%coords(lico)

 loc_glob__ = nblocs * matrix%sizeb_blocs(lico) + rest + 1

end function loc_glob__
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_from_global
!! NAME
!!  matrix_from_global
!!
!! FUNCTION
!!  Routine to fill a SCALAPACK matrix from a global PACKED matrix.
!!
!! INPUTS
!!  istwf_k= 2 if we have a real matrix else complex.
!!  reference= one-dimensional array with packed matrix.
!!
!! SIDE EFFECTS
!!  matrix= the matrix to process
!!
!! SOURCE

subroutine matrix_from_global(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: matrix
 integer,intent(in) :: istwf_k
 real(dp),intent(in) :: reference(*)

!Local variables-------------------------------
 integer :: i,j,iglob,jglob,ind
 real(dp) :: val_real
 complex(dp) :: val_cplx

! *********************************************************************

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)

     if (istwf_k/=2) then
        ind = jglob*(jglob-1)+2*iglob-1
        val_cplx = dcmplx(reference(ind),reference(ind+1))
        call matrix_set_local_cplx(matrix,i,j,val_cplx)
    else
       ind = (jglob*(jglob-1))/2 + iglob
       val_real = reference(ind)
       call matrix_set_local_real(matrix,i,j,val_real)
     end if

   end do
 end do

end subroutine matrix_from_global
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_from_global_sym
!! NAME
!!  matrix_from_global_sym
!!
!! FUNCTION
!!
!! INPUTS
!!  istwf_k= 2 if we have a real matrix else complex.
!!  reference= one-dimensional array
!!
!! SIDE EFFECTS
!!  matrix= the matrix to process
!!
!! SOURCE

subroutine matrix_from_global_sym(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout)  :: matrix
 real(dp),intent(in) :: reference(:)
 integer,intent(in) :: istwf_k

!Local variables-------------------------------
 integer :: i,j,iglob,jglob,ind
 complex(dp):: val_cplx
 real(dp) ::val_real

! *********************************************************************

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i,j,iglob,jglob)
     if(jglob < iglob) then
        ind = iglob*(iglob-1)+2*jglob-1
     else
        ind = jglob*(jglob-1)+2*iglob-1
     end if
     if (istwf_k/=2) then
        val_cplx = dcmplx(reference(ind),reference(ind+1))
        if(jglob < iglob) then
           call matrix_set_local_cplx(matrix,i,j,conjg(val_cplx))
        else
           call matrix_set_local_cplx(matrix,i,j,val_cplx)
        end if
     else
        ind = (ind + 1) / 2
        val_real = reference(ind)
        call matrix_set_local_real(matrix,i,j,val_real)
     end if
   end do
 end do

end subroutine matrix_from_global_sym
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_from_realmatrix
!! NAME
!!  matrix_from_realmatrix
!!
!! FUNCTION
!!  Routine to fill a SCALAPACK matrix from a real global matrix (FULL STORAGE MODE)
!!
!! INPUTS
!!  istwf_k= 2 if we have a real matrix else complex.
!!  reference= a real matrix
!!
!! SIDE EFFECTS
!!  matrix= the matrix to process
!!
!! SOURCE

subroutine matrix_from_realmatrix(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: matrix
 integer,intent(in) :: istwf_k
!arrays
 real(dp),intent(in) :: reference(:,:)

!Local variables-------------------------------
 integer :: i,j,iglob,jglob
 real(dp) :: val

! *********************************************************************

 ABI_UNUSED(istwf_k)

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)
     val = reference(iglob, jglob)
     call matrix_set_local_real(matrix,i,j,val)
   end do
 end do

end subroutine matrix_from_realmatrix
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_from_complexmatrix
!! NAME
!!  matrix_from_complexmatrix
!!
!! FUNCTION
!!  Routine to fill a SCALAPACK matrix from a global matrix (FULL STORAGE MODE)
!!
!! INPUTS
!!  istwf_k= 2 if we have a real matrix else complex.
!!  reference= a complex matrix
!!
!! SIDE EFFECTS
!!  matrix= the matrix to process
!!
!! SOURCE

subroutine matrix_from_complexmatrix(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: matrix
 integer,intent(in) :: istwf_k
!arrays
 real(dp),intent(in) :: reference(:,:)

!Local variables-------------------------------
 integer :: i,j,iglob,jglob
 complex(dpc) :: val

! *********************************************************************

 ABI_UNUSED(istwf_k)

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)
     val = dcmplx(reference(2*iglob-1, jglob),reference(2*iglob, jglob))
     call matrix_set_local_cplx(matrix,i,j,val)
   end do
 end do

end subroutine matrix_from_complexmatrix
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_to_global
!! NAME
!!  matrix_to_global
!!
!! FUNCTION
!!  Inserts a ScaLAPACK matrix into a global one in PACKED storage mode.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  istwf_k= 2 if we have a real matrix else complex.
!!  nband_k= number of bands at this k point for that spin polarization
!!
!! SIDE EFFECTS
!!  reference= one-dimensional array
!!
!! SOURCE

subroutine matrix_to_global(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: matrix
 integer,intent(in) :: istwf_k          !,nband_k
 real(dp),intent(inout) :: reference(*) !(nband_k*(nband_k+1))

!Local variables-------------------------------
 integer  :: i,j,iglob,jglob,ind

! *********************************************************************

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)

     ind = jglob*(jglob-1)+2*iglob-1
     if (ind <= matrix%sizeb_global(2)*(matrix%sizeb_global(2)+1)) then
        if (istwf_k/=2) then
         reference(ind)   = real(matrix_get_local_cplx(matrix,i,j))
         reference(ind+1) = aimag(matrix_get_local_cplx(matrix,i,j))
       else
          ind=(ind+1)/2 !real packed storage
          reference(ind) = matrix_get_local_real(matrix,i,j)
       end if
     end if
   end do
 end do

end subroutine matrix_to_global
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_to_realmatrix
!! NAME
!!  matrix_to_realmatrix
!!
!! FUNCTION
!!  Inserts a ScaLAPACK matrix into a real matrix in FULL STORAGE MODE.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  istwf_k= 2 if we have a real matrix else complex.
!!
!! SIDE EFFECTS
!!  reference= the matrix to fill
!!
!! SOURCE

subroutine matrix_to_realmatrix(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: matrix
 integer,intent(in) :: istwf_k
!arrays
 real(dp),intent(inout) :: reference(:,:)

!Local variables-------------------------------
 integer :: i,j,iglob,jglob
 !complex(dpc) :: zvar

! *********************************************************************

 ABI_UNUSED(istwf_k)

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)
     reference(iglob,jglob) = matrix_get_local_real(matrix,i,j)
   end do
 end do

end subroutine matrix_to_realmatrix
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_to_complexmatrix
!! NAME
!!  matrix_to_complexmatrix
!!
!! FUNCTION
!!  Inserts a ScaLAPACK matrix into a complex matrix in FULL STORAGE MODE.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  istwf_k= 2 if we have a real matrix else complex.
!!
!! SIDE EFFECTS
!!  reference= the matrix to fill
!!
!! SOURCE

subroutine matrix_to_complexmatrix(matrix, reference, istwf_k)

!Arguments ------------------------------------
 integer,intent(in) :: istwf_k
 class(matrix_scalapack),intent(in) :: matrix
!arrays
 complex(dpc),intent(inout) :: reference(:,:)

!Local variables-------------------------------
 integer  :: i,j,iglob,jglob

! *********************************************************************

 ABI_UNUSED(istwf_k)

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)
     reference(iglob,jglob) = matrix_get_local_cplx(matrix,i,j)
   end do
 end do

end subroutine matrix_to_complexmatrix
!!***

!----------------------------------------------------------------------

!!****f* m_slk/matrix_to_reference
!! NAME
!!  matrix_to_reference
!!
!! FUNCTION
!!  Routine to fill a full matrix with respect to a SCALAPACK matrix.
!!
!! INPUTS
!!  matrix= the matrix to process
!!  istwf_k= 2 if we have a real matrix else complex.
!!
!! SIDE EFFECTS
!!  reference= one-dimensional array
!!
!! SOURCE

subroutine matrix_to_reference(matrix, reference, istwf_k)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: matrix
 integer,intent(in) :: istwf_k
!arrays
 real(dp),intent(inout) :: reference(:,:)

!Local variables-------------------------------
 integer  :: i,j,iglob,jglob,ind

! *********************************************************************

 do i=1,matrix%sizeb_local(1)
   do j=1,matrix%sizeb_local(2)
     call matrix%loc2glob(i, j, iglob, jglob)

     if (istwf_k/=2) then
       ind=(iglob-1)*2+1
       reference(ind,jglob)   = real(matrix_get_local_cplx(matrix,i,j))
       reference(ind+1,jglob) = aimag(matrix_get_local_cplx(matrix,i,j))
     else
        ind=iglob
        reference(ind,jglob)   = matrix_get_local_real(matrix,i,j)
        !reference(ind+1,jglob) = 0._dp
     end if

   end do
 end do

end subroutine matrix_to_reference
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_matrix_from_global_dpc_2D
!! NAME
!!  slk_matrix_from_global_dpc_2D
!!
!! FUNCTION
!!  Routine to fill a complex SCALAPACK matrix with respect to a global matrix.
!!  target: Two-dimensional double precision complex matrix
!!
!! INPUTS
!!  glob_mat=Two-dimensional array containing the global matrix.
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix is used:
!!    = "U":  Upper triangular
!!    = "L":  Lower triangular
!!    = "A":  Full matrix (used for general complex matrices)
!!
!! SIDE EFFECTS
!!  mat<matrix_scalapack>=The distributed matrix.
!!    %buffer_cplx=Local buffer containg the value this node is dealing with.
!!
!! SOURCE

subroutine slk_matrix_from_global_dpc_2D(mat, uplo, glob_mat)

!Arguments ------------------------------------
!scalars
 class(matrix_scalapack),intent(inout)  :: mat
 character(len=*),intent(in) :: uplo
!array
 complex(dpc),intent(in) :: glob_mat(:,:)

!Local variables-------------------------------
 integer :: ii, jj, iglob, jglob

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "%buffer_cplx not allocated")

 select case (uplo(1:1))

 case ("A", "a")
   ! Full global matrix is used.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)
       mat%buffer_cplx(ii,jj) = glob_mat(iglob,jglob)
     end do
   end do

 case ("U", "u")
   ! Only the upper triangle of the global matrix is used.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)
       if (jglob>=iglob) then
         mat%buffer_cplx(ii,jj) =        glob_mat(iglob,jglob)
       else
         mat%buffer_cplx(ii,jj) = DCONJG(glob_mat(jglob,iglob))
       end if
     end do
   end do

 case ("L", "l")
   ! Only the lower triangle of the global matrix is used.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)
       if (jglob<=iglob) then
         mat%buffer_cplx(ii,jj) =        glob_mat(iglob,jglob)
       else
         mat%buffer_cplx(ii,jj) = DCONJG(glob_mat(jglob,iglob))
       end if
     end do
   end do

 case default
   ABI_BUG(" Wrong uplo: "//TRIM(uplo))
 end select

end subroutine slk_matrix_from_global_dpc_2D
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_matrix_from_global_dpc_1Dp
!! NAME
!!  slk_matrix_from_global_dpc_1Dp
!!
!! FUNCTION
!!  Routine to fill a complex SCALAPACK matrix with respect to a global matrix.
!!  target: double precision complex matrix in packed form.
!!
!! INPUTS
!!  glob_pmat(n*(n+1)/2)=One-dimensional array containing the global matrix A packed columnwise in a linear array.
!!    The j-th column of A is stored in the array glob_pmat as follows:
!!      if uplo = "U", glob_pmat(i + (j-1)*j/2)       = A(i,j) for 1<=i<=j;
!!      if uplo = "L", glob_pmat(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!!      where n is the number of rows or columns in the global matrix.
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix is used:
!!    = "U":  Upper triangular
!!    = "L":  Lower triangular
!!
!! SIDE EFFECTS
!!  mat<matrix_scalapack>=The distributed matrix.
!!    %buffer_cplx=Local buffer containg the value this node is dealing with.
!!
!! SOURCE

subroutine slk_matrix_from_global_dpc_1Dp(mat,uplo,glob_pmat)

!Arguments ------------------------------------
!scalars
 class(matrix_scalapack),intent(inout)  :: mat
 character(len=*),intent(in) :: uplo
!array
 complex(dpc),intent(in) :: glob_pmat(:)

!Local variables-------------------------------
 integer :: ii,jj,iglob,jglob,ind,n
 real(dp) :: szm

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "%buffer_cplx not allocated")

 szm = SIZE(glob_pmat)
 n = NINT( (-1 + SQRT(one+8*szm) )*half )
 if (n*(n+1)/2 /= SIZE(glob_pmat)) then
   ABI_ERROR("Buggy compiler")
 end if

 select case (uplo(1:1))

 case ("U", "u")
   ! Only the upper triangle of the global matrix is used.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)

       if (jglob>=iglob) then
         ind = iglob + jglob*(jglob-1)/2
         mat%buffer_cplx(ii,jj) =        glob_pmat(ind)
       else
         ind = jglob + iglob*(iglob-1)/2
         mat%buffer_cplx(ii,jj) = DCONJG( glob_pmat(ind) )
       end if

     end do
   end do

 case ("L", "l")
   ! Only the lower triangle of the global matrix is used.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)

       if (jglob<=iglob) then
         ind = iglob + (jglob-1)*(2*n-jglob)/2
         mat%buffer_cplx(ii,jj) =        glob_pmat(ind)
       else
         ind = jglob + (iglob-1)*(2*n-iglob)/2
         mat%buffer_cplx(ii,jj) = DCONJG( glob_pmat(ind) )
       end if

     end do
   end do

 case default
   ABI_BUG(" Wrong uplo: "//TRIM(uplo))
 end select

end subroutine slk_matrix_from_global_dpc_1Dp
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_matrix_to_global_dpc_2D
!! NAME
!!  slk_matrix_to_global_dpc_2D
!!
!! FUNCTION
!!  Fill a global matrix with respect to a SCALAPACK matrix.
!!  target: Two-dimensional double precision complex matrix.
!!
!! INPUTS
!!  Slk_mat<matrix_scalapack>=The distributed matrix.
!!  uplo=String specifying whether the upper or lower triangular part of the global matrix has to be filled:
!!    = "U":  Upper triangular
!!    = "L":  Lower triangular
!!    = "A":  Full matrix is filled (used for general complex matrices)
!!
!! SIDE EFFECTS
!!  glob_mat=The global matrix where the entries owned by this processors have been overwritten.
!!  Note that the remaing entries not treated by this node are not changed.
!!
!! SOURCE

subroutine slk_matrix_to_global_dpc_2D(mat, uplo, glob_mat)

!Arguments ------------------------------------
!scalaras
 class(matrix_scalapack),intent(in) :: mat
 character(len=*),intent(in) :: uplo
!arrays
 complex(dpc),intent(inout) :: glob_mat(:,:)

!Local variables-------------------------------
 integer :: ii,jj,iglob,jglob

!************************************************************************

 select case (uplo(1:1))

 case ("A", "a")
   ! Full global matrix has to be filled.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)
       glob_mat(iglob,jglob) = mat%buffer_cplx(ii,jj)
     end do
   end do

 case ("U", "u")
   ! Only the upper triangle of the global matrix is filled.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)
       if (jglob>=iglob) glob_mat(iglob,jglob) = mat%buffer_cplx(ii,jj)
     end do
   end do

 case ("L", "l")
   ! Only the lower triangle of the global matrix is filled.
   do jj=1,mat%sizeb_local(2)
     do ii=1,mat%sizeb_local(1)
       call mat%loc2glob(ii, jj, iglob, jglob)
       if (jglob<=iglob) glob_mat(iglob,jglob) = mat%buffer_cplx(ii,jj)
     end do
   end do

 case default
   ABI_BUG(" Wrong uplo: "//TRIM(uplo))
 end select

end subroutine slk_matrix_to_global_dpc_2D
!!***

!----------------------------------------------------------------------

!!****f* m_slk/my_locr
!! NAME
!! my_locr
!!
!! FUNCTION
!!  Method of matrix_scalapack wrapping the scaLAPACK tool LOCr.
!!
!! OUTPUT
!!  my_locr= For the meaning see NOTES below.
!!
!! NOTES
!!  Let K be the number of rows or columns of a distributed matrix, and assume that its process grid has dimension p x q.
!!  LOCr( K ) denotes the number of elements of K that a process would receive if K were distributed over the p
!!  processes of its process column.
!!  Similarly,  LOCc(  K  ) denotes the number of elements of K that a process would receive if K were distributed over
!!  the q processes of its process row.
!!  The values of LOCr() and LOCc() may be determined via a call to the ScaLAPACK tool function, NUMROC:
!!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).  An upper bound for these quantities may  be  computed
!!  by:
!!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!!
!! SOURCE

integer function my_locr(mat)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: mat

!Local variables-------------------------------
#ifdef HAVE_LINALG_SCALAPACK
 integer :: M, MB_A, MYROW, RSRC_A, NPROW
 integer,external :: NUMROC

! *************************************************************************

 M      = mat%descript%tab(M_ )      ! The number of rows in the global matrix.
 MB_A   = mat%descript%tab(MB_)      ! The number of rows in a block.
 MYROW  = mat%processor%coords(1)    ! The row index of my processor
 RSRC_A = mat%descript%tab(RSRC_)    ! The row of the processors at the beginning.
 NPROW  = mat%processor%grid%dims(1) ! The number of processors per row in the Scalapack grid.

 my_locr = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW )
#endif

end function my_locr
!!***

!----------------------------------------------------------------------

!!****f* m_slk/my_locc
!! NAME
!! my_locc
!!
!! FUNCTION
!!  Method of matrix_scalapack wrapping the scaLAPACK tool LOCc.
!!
!! OUTPUT
!!  my_locc= For the meaning see NOTES below.
!!
!! NOTES
!!  Let K be the number of rows or columns of a distributed matrix, and assume that its process grid has dimension p x q.
!!  LOCr( K ) denotes the number of elements of K that a process would receive if K were distributed over the p
!!  processes of its process column.
!!  Similarly,  LOCc(  K  ) denotes the number of elements of K that a process would receive if K were distributed over
!!  the q processes of its process row.
!!  The values of LOCr() and LOCc() may be determined via a call to the ScaLAPACK tool function, NUMROC:
!!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).  An upper bound for these quantities may  be  computed
!!  by:
!!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!!
!! SOURCE

integer function my_locc(mat)

!Arguments ------------------------------------
 class(basemat_t),intent(in) :: mat

#ifdef HAVE_LINALG_SCALAPACK
!Local variables-------------------------------
 integer :: N, NB_A, MYCOL, CSRC_A, NPCOL
 integer,external :: NUMROC

! *************************************************************************

 N      = mat%descript%tab(N_ )      ! The number of columns in the global matrix.
 NB_A   = mat%descript%tab(NB_)      ! The number of columns in a block.
 MYCOL  = mat%processor%coords(2)    ! The column index of my processor
 CSRC_A = mat%descript%tab(CSRC_)    ! The column of the processors at the beginning.
 NPCOL  = mat%processor%grid%dims(2) ! The number of processors per column in the Scalapack grid.

 my_locc = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL )
#endif

end function my_locc
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_pgemm_dp
!! NAME
!!  slk_pgemm_dp
!!
!! FUNCTION
!!  Extended matrix * matrix product: C := alpha*A*B + beta*C
!!  For a simple matrix vector product, one can simply pass alpha = cone and beta = czero
!!
!! INPUTS
!!  matrix1= first ScaLAPACK matrix (matrix A)
!!  matrix2= second ScaLAPACK matrix (matrix B)
!!  alpha= scalar multiplicator for the A*B product
!!  beta= scalar multiplicator for the C matrix
!!  [ija, ijb, ijc]:  (global). The row and column indices in the distributed matrix A/B/C
!!    indicating the first row and the first column of the submatrix, respectively. Default [1, 1]
!!
!! OUTPUT
!!  results= ScaLAPACK matrix coming out of the operation
!!
!! NOTES
!! The ESLL manual says that
!!      "matrices matrix1 and matrix2 must have no common elements otherwise, results are unpredictable."
!! However the official scaLAPACK documentation does not report this (severe) limitation.
!!
!! SOURCE

subroutine slk_pgemm_dp(transa, transb, matrix1, alpha, matrix2, beta, results, &
                        ija, ijb, ijc) ! optional

!Arguments ------------------------------------
 character(len=1),intent(in) :: transa, transb
 class(matrix_scalapack),intent(in) :: matrix1, matrix2
 class(matrix_scalapack),intent(inout) :: results
 complex(dpc),intent(in) :: alpha, beta
 integer,optional,intent(in) :: ija(2), ijb(2), ijc(2)

!Local variables-------------------------------
 integer :: mm, nn, kk, ija__(2), ijb__(2), ijc__(2)

!************************************************************************

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijb__ = [1, 1]; if (present(ijb)) ijb__ = ijb
 ijc__ = [1, 1]; if (present(ijc)) ijc__ = ijc

 mm  = matrix1%sizeb_global(1)
 nn  = matrix2%sizeb_global(2)
 kk  = matrix1%sizeb_global(2)

 if (toupper(transa) /= 'N') then
   mm = matrix1%sizeb_global(2)
   kk = matrix1%sizeb_global(1)
 end if
 if (toupper(transb) /= 'N') nn = matrix2%sizeb_global(1)

#ifdef HAVE_LINALG_SCALAPACK
 ! pzgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
 if (matrix1%istwf_k /= 2) then
   call PZGEMM(transa, transb, mm, nn, kk, alpha, &
               matrix1%buffer_cplx, ija__(1), ija__(2), matrix1%descript%tab, &
               matrix2%buffer_cplx, ijb__(1), ijb__(2), matrix2%descript%tab, &
               beta, results%buffer_cplx, ijc__(1), ijc__(2), results%descript%tab)
 else
   call PDGEMM(transa, transb, mm, nn, kk, real(alpha, kind=dp), &
               matrix1%buffer_cplx, ija__(1), ija__(2), matrix1%descript%tab, &
               matrix2%buffer_cplx, ijb__(1), ijb__(2), matrix2%descript%tab, &
               real(beta, kind=dp), results%buffer_cplx, ijc__(1), ijc__(2), results%descript%tab)
 end if
#endif

end subroutine slk_pgemm_dp
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_pgemm_sp
!! NAME
!!  slk_pgemm_sp
!!
!! FUNCTION
!!  Extended matrix * matrix product: C := alpha*A*B + beta*C
!!  For a simple matrix vector product, one can simply pass alpha = cone and beta = czero
!!
!! INPUTS
!!  matrix1= first ScaLAPACK matrix (matrix A)
!!  matrix2= second ScaLAPACK matrix (matrix B)
!!  alpha= scalar multiplicator for the A*B product
!!  beta= scalar multiplicator for the C matrix
!!  [ija, ijb, ijc]:  (global). The row and column indices in the distributed matrix A/B/C
!!    indicating the first row and the first column of the submatrix, respectively. Default [1, 1]
!!
!! OUTPUT
!!  results= ScaLAPACK matrix coming out of the operation
!!
!! NOTES
!! The ESLL manual says that
!!      "matrices matrix1 and matrix2 must have no common elements otherwise, results are unpredictable."
!! However the official scaLAPACK documentation does not report this (severe) limitation.
!!
!! SOURCE

subroutine slk_pgemm_sp(transa, transb, matrix1, alpha, matrix2, beta, results, &
                      ija, ijb, ijc) ! optional

!Arguments ------------------------------------
 character(len=1),intent(in) :: transa, transb
 class(slkmat_sp_t),intent(in) :: matrix1, matrix2
 class(slkmat_sp_t),intent(inout) :: results
 complex(sp),intent(in) :: alpha, beta
 integer,optional,intent(in) :: ija(2), ijb(2), ijc(2)

!Local variables-------------------------------
 integer :: mm, nn, kk, ija__(2), ijb__(2), ijc__(2)

!************************************************************************

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijb__ = [1, 1]; if (present(ijb)) ijb__ = ijb
 ijc__ = [1, 1]; if (present(ijc)) ijc__ = ijc

 mm  = matrix1%sizeb_global(1)
 nn  = matrix2%sizeb_global(2)
 kk  = matrix1%sizeb_global(2)

 if (toupper(transa) /= 'N') then
   mm = matrix1%sizeb_global(2)
   kk = matrix1%sizeb_global(1)
 end if
 if (toupper(transb) /= 'N') nn = matrix2%sizeb_global(1)

#ifdef HAVE_LINALG_SCALAPACK
 ! pzgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
 if (matrix1%istwf_k /= 2) then
   call PCGEMM(transa, transb, mm, nn, kk, alpha, &
               matrix1%buffer_cplx, ija__(1), ija__(2), matrix1%descript%tab, &
               matrix2%buffer_cplx, ijb__(1), ijb__(2), matrix2%descript%tab, &
               beta, results%buffer_cplx, ijc__(1), ijc__(2), results%descript%tab)
 else
   call PSGEMM(transa, transb, mm, nn, kk, real(alpha, kind=sp), &
               matrix1%buffer_cplx, ija__(1), ija__(2), matrix1%descript%tab, &
               matrix2%buffer_cplx, ijb__(1), ijb__(2), matrix2%descript%tab, &
               real(beta, kind=sp), results%buffer_cplx, ijc__(1), ijc__(2), results%descript%tab)
 end if
#endif

end subroutine slk_pgemm_sp
!!***

!----------------------------------------------------------------------

!!****f* m_slk/compute_eigen_problem
!! NAME
!!  compute_eigen_problem
!!
!! FUNCTION
!!  Calculation of eigenvalues and eigenvectors: A * X = lambda * X, complex and real cases.
!!
!! INPUTS
!!  processor= descriptor of a processor
!!  matrix= the matrix to process
!!  comm= MPI communicator
!!  istwf_k= 2 if we have a real matrix else complex.
!!  [nev]= Number of eigenvalues needed. Default: full set
!!  [use_gpu_elpa]= Flag to activate the use of GPU (ELPA only)
!!
!! OUTPUT
!!  results= ScaLAPACK matrix coming out of the operation (global dimensions must be equal to matrix
!!         even if only a part of the eigenvectors is needed.
!!  eigen= eigenvalues of the matrix dimensioned as the global size of the square matrix
!!         even if only a part of the eigenvalues is needed.
!!
!! SOURCE

subroutine compute_eigen_problem(processor, matrix, results, eigen, comm, istwf_k, &
&                                nev, use_gpu_elpa) ! Optional arguments

#ifdef HAVE_LINALG_ELPA
  !Arguments ------------------------------------
  class(processor_scalapack),intent(in) :: processor
  class(matrix_scalapack),intent(inout) :: matrix
  class(matrix_scalapack),intent(inout) :: results
  DOUBLE PRECISION,intent(inout) :: eigen(:)
  integer,intent(in)  :: comm,istwf_k
  integer,optional,intent(in) :: nev
  integer,optional,intent(in) :: use_gpu_elpa

  !Local variables ------------------------------
  type(elpa_hdl_t) :: elpa_hdl
  integer :: nev__,use_gpu_elpa_

!************************************************************************

  nev__ = matrix%sizeb_global(1); if (present(nev)) nev__ = nev
  use_gpu_elpa_=0
#ifdef HAVE_LINALG_ELPA
  if (present(use_gpu_elpa)) use_gpu_elpa_=use_gpu_elpa
#endif

  call elpa_func_allocate(elpa_hdl,gpu=use_gpu_elpa_)
  call elpa_func_set_matrix(elpa_hdl,matrix%sizeb_global(1),matrix%sizeb_blocs(1),nev__,&
&                           matrix%sizeb_local(1),matrix%sizeb_local(2))
  call elpa_func_get_communicators(elpa_hdl,processor%comm,processor%coords(1),processor%coords(2))

  if (istwf_k/=2) then
    call elpa_func_solve_evp_1stage(elpa_hdl,matrix%buffer_cplx,results%buffer_cplx,eigen,nev__)
  else
    call elpa_func_solve_evp_1stage(elpa_hdl,matrix%buffer_real,results%buffer_real,eigen,nev__)
  end if

  call elpa_func_deallocate(elpa_hdl)

#else
  !Arguments ------------------------------------
  class(processor_scalapack),intent(in)       :: processor
  class(matrix_scalapack),intent(in)          :: matrix
  class(matrix_scalapack),intent(inout)       :: results
  DOUBLE PRECISION,intent(inout) :: eigen(:)
  integer,intent(in)  :: comm,istwf_k
  integer,optional,intent(in) :: nev
  integer,optional,intent(in) :: use_gpu_elpa

#ifdef HAVE_LINALG_SCALAPACK
  !Local variables-------------------------------
  integer            :: LRWORK,LIWORK,LCWORK,INFO
  character(len=500) :: msg

  integer         , dimension(1) :: IWORK_tmp
  DOUBLE PRECISION, dimension(1) :: RWORK_tmp
  complex(dpc)     , dimension(1) :: CWORK_tmp

  integer         , allocatable  :: IWORK(:)
  DOUBLE PRECISION, allocatable  :: RWORK(:)
  complex(dpc)     , allocatable  :: CWORK(:)

  integer,          allocatable :: ICLUSTR(:)
  integer,          allocatable :: IFAIL(:)
  DOUBLE PRECISION, allocatable :: GAP(:)

  DOUBLE PRECISION            :: ABSTOL,ORFAC
  integer,          parameter :: IZERO=0

  integer ::  M,NZ,ierr,TWORK_tmp(3),TWORK(3) ! IA,JA,IZ,JZ,
  integer :: nev__, il, iu
  character(len=1) :: range

  DOUBLE PRECISION, external :: PDLAMCH

! *************************************************************************

  ABI_UNUSED(use_gpu_elpa) ! No GPU implementation is using scaLAPACK
  nev__ = matrix%sizeb_global(1); range = "A"; il = 0; iu = 0
  if (present(nev)) then
    nev__ = nev; range = "I"; il = 1; iu = nev
  end if

  ! Initialisation
  INFO   = 0
  ABSTOL = zero
  ORFAC  = -1.D+0

  ! Allocation of the variables for the results of the calculations
  ABI_MALLOC(IFAIL,(matrix%sizeb_global(2)))
  ABI_MALLOC(ICLUSTR,(2*processor%grid%dims(1)*processor%grid%dims(2)))
  ABI_MALLOC(GAP,(processor%grid%dims(1)*processor%grid%dims(2)))

  ! Get the size of the work arrays
  if (istwf_k/=2) then
     call PZHEEVX('V', range, 'U',&
&      matrix%sizeb_global(2),&
&      matrix%buffer_cplx,1,1,matrix%descript%tab, &
&      ZERO,ZERO,il,iu,ABSTOL,&
&      m,nz,eigen,ORFAC, &
&      results%buffer_cplx,1,1,results%descript%tab, &
&      CWORK_tmp,-1,RWORK_tmp,-1,IWORK_tmp,-1,&
&      IFAIL,ICLUSTR,GAP,INFO)
  else
     call PDSYEVX('V', range, 'U',&
&      matrix%sizeb_global(2),&
&      matrix%buffer_real,1,1,matrix%descript%tab, &
&      ZERO,ZERO,il,iu,ABSTOL,&
&      m,nz,eigen,ORFAC, &
&      results%buffer_real,1,1,results%descript%tab, &
&      RWORK_tmp,-1,IWORK_tmp,-1,&
&      IFAIL,ICLUSTR,GAP,INFO)
  end if

  if (INFO/=0) then
     write(msg,'(A,I6)') "Problem to compute workspace to use ScaLAPACK, INFO=",INFO
     ABI_ERROR(msg)
  endif

  TWORK_tmp(1) = IWORK_tmp(1)
  TWORK_tmp(2) = INT(RWORK_tmp(1))
  TWORK_tmp(3) = INT(real(CWORK_tmp(1)))

  ! Get the maximum of the size of the work arrays processor%comm
  call MPI_ALLREDUCE(TWORK_tmp,TWORK,3,MPI_integer,MPI_MAX,comm,ierr)

  LIWORK = TWORK(1)
  LRWORK = TWORK(2) + matrix%sizeb_global(2) *(matrix%sizeb_global(2)-1)
  LCWORK = TWORK(3)

  ! Allocation of the work arrays
  if (LIWORK>0) then
    ABI_MALLOC(IWORK,(LIWORK))
    IWORK(:) = 0
  else
    ABI_MALLOC(IWORK,(1))
  end if
  if (LRWORK>0) then
    ABI_MALLOC(RWORK,(LRWORK))
    RWORK(:) = 0._dp
  else
    ABI_MALLOC(RWORK,(1))
  end if
  if (LCWORK>0) then
    ABI_MALLOC(CWORK,(LCWORK))
    CWORK(:) = (0._dp,0._dp)
  else
    ABI_MALLOC(CWORK,(1))
  end if

  ! prototype
  !call pzheevx(jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w,
  !             orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)

  ! Call the calculation routine
  if (istwf_k/=2) then
    ! write(std_out,*) 'I am using PZHEEVX'
    call PZHEEVX('V', range, 'U',&
     matrix%sizeb_global(2),&
     matrix%buffer_cplx,1,1,matrix%descript%tab, &
     ZERO,ZERO,il,iu,ABSTOL,&
     m,nz,eigen,ORFAC, &
     results%buffer_cplx,1,1,results%descript%tab, &
     CWORK,LCWORK,RWORK,LRWORK,IWORK,LIWORK,&
     IFAIL,ICLUSTR,GAP,INFO)
  else
    ! write(std_out,*) ' I am using PDSYEVX'
    call PDSYEVX('V', range, 'U',&
     matrix%sizeb_global(2),&
     matrix%buffer_real,1,1,matrix%descript%tab, &
     ZERO,ZERO,il,iu,ABSTOL,&
     m,nz,eigen,ORFAC, &
     results%buffer_real,1,1,results%descript%tab, &
     RWORK,LRWORK,IWORK,LIWORK,&
     IFAIL,ICLUSTR,GAP,INFO)
  endif

  ! MG: TODO: Recheck the computation of the workspace as I got INFO 2 with a 5x5x5 si supercell.
  if (INFO/=0) then
    write(msg,'(A,I0)') "Problem to compute eigenvalues and eigenvectors with ScaLAPACK, INFO=",INFO
    ABI_ERROR(msg)
  endif

  ABI_FREE(IFAIl)
  ABI_FREE(ICLUSTR)
  ABI_FREE(GAP)

  ABI_SFREE(IWORK)
  ABI_SFREE(RWORK)
  ABI_SFREE(CWORK)
#endif
#endif
  return

end subroutine compute_eigen_problem
!!***

!----------------------------------------------------------------------

!!****f* m_slk/solve_gevp_complex
!! NAME
!!  solve_gevp_complex
!!
!! FUNCTION
!!  Calculation of eigenvalues and eigenvectors: A * X = lambda * B * X
!!  complex and real cases.
!!
!! INPUTS
!!  processor= descriptor of a processor
!!  matrix1= first ScaLAPACK matrix (matrix A)
!!  matrix2= second ScaLAPACK matrix (matrix B)
!!  comm= MPI communicator
!!  istwf_k= 2 if we have a real matrix else complex.
!!  [use_gpu_elpa]= Flag to activate the use of GPU (ELPA only)
!!
!! SIDE EFFECTS
!!  results= ScaLAPACK matrix coming out of the operation
!!  eigen= eigenvalues of the matrix
!!
!! SOURCE

#ifdef HAVE_LINALG_ELPA

subroutine solve_gevp_complex(na,nev,na_rows,na_cols,nblk,a,b,ev,z,tmp1,tmp2, &
&                             my_prow,my_pcol,np_rows,np_cols,sc_desc,comm,&
&                             use_gpu_elpa) ! Optional parameter

  !-Arguments
  integer,intent(in) :: na
  integer,intent(in) :: nev
  integer,intent(in) :: na_rows,na_cols
  integer,intent(in) :: nblk
  integer,intent(in) :: my_pcol,my_prow
  integer,intent(in) :: np_cols,np_rows
  integer,intent(in) :: sc_desc(9)
  integer,intent(in) :: comm
  integer,optional,intent(in) :: use_gpu_elpa
  real*8 :: ev(na)
  complex*16 :: a(na_rows,na_cols),b(na_rows,na_cols),z(na_rows,na_cols)
  complex*16 :: tmp1(na_rows,na_cols),tmp2(na_rows,na_cols)
  !-Local variables
  integer :: i, n_col, n_row, use_gpu_elpa_
  integer,external :: indxl2g,numroc
  complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)
  type(elpa_hdl_t) :: elpa_hdl

! *************************************************************************
  
  use_gpu_elpa_=0
#ifdef HAVE_LINALG_ELPA
  if (present(use_gpu_elpa)) use_gpu_elpa_=use_gpu_elpa
#endif

! Allocate ELPA handle
  call elpa_func_allocate(elpa_hdl,blacs_ctx=sc_desc(CTXT_),gpu=use_gpu_elpa_)
  call elpa_func_set_matrix(elpa_hdl,na,nblk,nev,na_rows,na_cols)
  call elpa_func_get_communicators(elpa_hdl,comm,my_prow,my_pcol)

  call elpa_func_solve_gevp_2stage(elpa_hdl,a,b,z,ev,nev)

  call elpa_func_deallocate(elpa_hdl)

end subroutine solve_gevp_complex

!----------------------------------------------------------------------

subroutine solve_gevp_real(na,nev,na_rows,na_cols,nblk,a,b,ev,z,tmp1,tmp2, &
&                          my_prow,my_pcol,np_rows,np_cols,sc_desc,comm, &
&                          use_gpu_elpa) ! Optional argument

  !-Arguments
  integer,intent(in) :: na
  integer,intent(in) :: nev
  integer,intent(in) :: na_rows,na_cols
  integer,intent(in) :: nblk
  integer,intent(in) :: my_pcol,my_prow
  integer,intent(in) :: np_cols,np_rows
  integer,intent(in) :: sc_desc(9)
  integer,intent(in) :: comm
  integer,optional,intent(in) :: use_gpu_elpa
  real*8 :: ev(na)
  real*8 :: a(na_rows,na_cols),b(na_rows,na_cols),z(na_rows,na_cols)
  real*8::tmp1(na_rows,na_cols),tmp2(na_rows,na_cols)
  !-Local variables
  integer :: i, n_col, n_row, use_gpu_elpa_
  integer,external :: indxl2g,numroc
  type(elpa_hdl_t) :: elpa_hdl

! *************************************************************************

  use_gpu_elpa_=0
#ifdef HAVE_LINALG_ELPA
  if (present(use_gpu_elpa)) use_gpu_elpa_=use_gpu_elpa
#endif

! Allocate ELPA handle
  call elpa_func_allocate(elpa_hdl,blacs_ctx=sc_desc(CTXT_),gpu=use_gpu_elpa_)
  call elpa_func_set_matrix(elpa_hdl,na,nblk,nev,na_rows,na_cols)
  call elpa_func_get_communicators(elpa_hdl,comm,my_prow,my_pcol)

  !FIXME Need to figure out why generalized_eigenvectors doesn't work in this real case
  !      while it is fine with complex case
  if(.false.) then
    call elpa_func_solve_gevp_2stage(elpa_hdl,a,b,tmp1,ev,nev)
  else
    ! 1. Calculate Cholesky factorization of Matrix B = U**T * U
    !    and invert triangular matrix U
    call elpa_func_cholesky(elpa_hdl,b)
    call elpa_func_invert_triangular(elpa_hdl,b)
    ! 2. Calculate U**-T * A * U**-1
    ! 2a. tmp1 = U**-T * A
    call elpa_func_hermitian_multiply(elpa_hdl,'U','L',na,b,a,na_rows,na_cols,tmp1,na_rows,na_cols)
    ! 2b. tmp2 = tmp1**T
    call pdtran(na,na,1.d0,tmp1,1,1,sc_desc,0.d0,tmp2,1,1,sc_desc)
    ! 2c. A =  U**-T * tmp2 ( = U**-T * Aorig * U**-1 )
    call elpa_func_hermitian_multiply(elpa_hdl,'U','U',na,b,tmp2,na_rows,na_cols,a,na_rows,na_cols)
    ! A is only set in the upper half, solve_evp_real needs a full matrix
    ! Set lower half from upper half
    call pdtran(na,na,1.d0,a,1,1,sc_desc,0.d0,tmp1,1,1,sc_desc)
    do i=1,na_cols
       ! Get global column corresponding to i and number of local rows up to
       ! and including the diagonal, these are unchanged in A
       n_col = indxl2g(i,     nblk, my_pcol, 0, np_cols)
       n_row = numroc (n_col, nblk, my_prow, 0, np_rows)
       a(n_row+1:na_rows,i) = tmp1(n_row+1:na_rows,i)
    enddo
    ! 3. Calculate eigenvalues/eigenvectors of U**-T * A * U**-1
    !    Eigenvectors go to tmp1
    call elpa_func_solve_evp_1stage(elpa_hdl,a,tmp1,ev,nev)
    ! 4. Backtransform eigenvectors: Z = U**-1 * tmp1
    !    hermitian_multiply needs the transpose of U**-1, thus tmp2 = (U**-1)**T
    call pdtran(na,na,1.d0,b,1,1,sc_desc,0.d0,tmp2,1,1,sc_desc)
    call elpa_func_hermitian_multiply(elpa_hdl,'L','N',nev,tmp2,tmp1,na_rows,na_cols,z,na_rows,na_cols)
  end if

  call elpa_func_deallocate(elpa_hdl)

 end subroutine solve_gevp_real
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_slk/compute_generalized_eigen_problem
!! NAME
!!  compute_generalized_eigen_problem
!!
!! FUNCTION
!!  Calculation of eigenvalues and eigenvectors of the generalized eigenvalue problem: A * X = lambda * B X
!!  complex and real cases.
!!
!! INPUTS
!!  processor= descriptor of a processor
!!  matrix1=  A matrix
!!  matrix2=  B matrix
!!  comm= MPI communicator
!!  istwf_k= 2 if we have a real matrix else complex.
!!  [nev]= Number of eigenvalues needed. Default: full set
!!  [use_gpu_elpa]= Flag to activate the use of GPU (ELPA only)
!!
!! OUTPUT
!!  results= ScaLAPACK matrix coming out of the operation (global dimensions must be equal to matrix
!!         even if only a part of the eigenvectors is needed.
!!  eigen= eigenvalues of the matrix dimensioned as the global size of the square matrix
!!         even if only a part of the eigenvalues is needed.
!!
!! SOURCE

subroutine compute_generalized_eigen_problem(processor,matrix1,matrix2,results,eigen,comm,istwf_k,&
&                                            nev,use_gpu_elpa) ! Optional arguments

#ifdef HAVE_LINALG_ELPA
!Arguments ------------------------------------
  class(processor_scalapack),intent(in)       :: processor
  class(matrix_scalapack),intent(in)          :: matrix1,matrix2
  class(matrix_scalapack),intent(inout)       :: results
  DOUBLE PRECISION,intent(inout) :: eigen(:)
  integer,intent(in)  :: comm,istwf_k
  integer,optional,intent(in) :: nev
  integer,optional,intent(in) :: use_gpu_elpa
!Local
  type(matrix_scalapack) :: tmp1, tmp2
  integer :: i,n_col, n_row, nev__,use_gpu_elpa__
  integer,external :: indxl2g,numroc

  nev__ = matrix1%sizeb_global(2); if (present(nev)) nev__ = nev
  use_gpu_elpa__ = 0
#ifdef HAVE_LINALG_ELPA
  if (present(use_gpu_elpa)) use_gpu_elpa__ = use_gpu_elpa
#endif

  call tmp1%init(matrix1%sizeb_global(1),matrix1%sizeb_global(2),processor,istwf_k)
  call tmp2%init(matrix1%sizeb_global(1),matrix1%sizeb_global(2),processor,istwf_k)

  if (istwf_k/=2) then
     call solve_gevp_complex(matrix1%sizeb_global(1), nev__, &
&          matrix1%sizeb_local(1),matrix1%sizeb_local(2),matrix1%sizeb_blocs(1), &
&          matrix1%buffer_cplx,matrix2%buffer_cplx,eigen,results%buffer_cplx, &
&          tmp1%buffer_cplx,tmp2%buffer_cplx, &
&          processor%coords(1),processor%coords(2), &
&          processor%grid%dims(1),processor%grid%dims(2), &
&          matrix1%descript%tab,processor%comm,use_gpu_elpa=use_gpu_elpa__)
  else
     call solve_gevp_real(matrix1%sizeb_global(1), nev__, &
&          matrix1%sizeb_local(1),matrix1%sizeb_local(2),matrix1%sizeb_blocs(1), &
&          matrix1%buffer_real,matrix2%buffer_real,eigen,results%buffer_real, &
&          tmp1%buffer_real,tmp2%buffer_real, &
&          processor%coords(1),processor%coords(2), &
&          processor%grid%dims(1),processor%grid%dims(2), &
&          matrix1%descript%tab,processor%comm,use_gpu_elpa=use_gpu_elpa__)
  end if
  call tmp1%free()
  call tmp2%free()

#else
!Arguments ------------------------------------
  class(processor_scalapack),intent(in)       :: processor
  class(matrix_scalapack),intent(in)          :: matrix1,matrix2
  class(matrix_scalapack),intent(inout)       :: results
  DOUBLE PRECISION,intent(inout) :: eigen(:)
  integer,intent(in)  :: comm,istwf_k
  integer,optional,intent(in) :: nev
  integer,optional,intent(in) :: use_gpu_elpa

#ifdef HAVE_LINALG_SCALAPACK
!Local variables-------------------------------
  integer            :: LRWORK,LIWORK,LCWORK,INFO
  character(len=500) :: msg
  integer         , dimension(1) :: IWORK_tmp
  DOUBLE PRECISION, dimension(1) :: RWORK_tmp
  complex(dpc)     , dimension(1) :: CWORK_tmp

  integer         , allocatable  :: IWORK(:)
  DOUBLE PRECISION, allocatable  :: RWORK(:)
  complex(dpc)     , allocatable  :: CWORK(:)
  integer,          allocatable :: ICLUSTR(:)
  integer,          allocatable :: IFAIL(:)
  DOUBLE PRECISION, allocatable :: GAP(:)
  DOUBLE PRECISION            :: ABSTOL,ORFAC
  integer         , parameter :: IZERO=0
  integer ::  M,NZ,ierr,TWORK_tmp(3),TWORK(3) ! IA,JA,IZ,JZ,
  character(len=1) :: range
  integer :: nev__, il, iu
  DOUBLE PRECISION, external :: PDLAMCH

! *************************************************************************

  ABI_UNUSED(use_gpu_elpa) ! No GPU implementation is using scaLAPACK
  nev__ = matrix1%sizeb_global(2); range = "A"; il = 0; iu = 0
  if (present(nev)) then
    nev__ = nev; range = "I"; il = 1; iu = nev
  end if

  ! Initialisation
  INFO   = 0
  ABSTOL = zero
  ORFAC  = -1.D+0

  ! Allocate the arrays for the results of the calculation
  ABI_MALLOC(IFAIL  ,(matrix1%sizeb_global(2)))
  ABI_MALLOC(ICLUSTR,(2*processor%grid%dims(1)*processor%grid%dims(2)))
  ABI_MALLOC(GAP    ,(  processor%grid%dims(1)*processor%grid%dims(2)))

  ! Get the size of the work arrays
  if (istwf_k /= 2) then
     call PZHEGVX(1, 'V', range, 'U',&
&      matrix1%sizeb_global(2),&
&      matrix1%buffer_cplx,1,1,matrix1%descript%tab, &
&      matrix2%buffer_cplx,1,1,matrix2%descript%tab, &
&      ZERO,ZERO,il,iu,ABSTOL,&
&      m,nz,eigen,ORFAC, &
&      results%buffer_cplx,1,1,results%descript%tab, &
&      CWORK_tmp,-1,RWORK_tmp,-1,IWORK_tmp,-1,&
&      IFAIL,ICLUSTR,GAP,INFO)
  else
     call PDSYGVX(1,'V',range,'U',&
&      matrix1%sizeb_global(2),&
&      matrix1%buffer_real,1,1,matrix1%descript%tab, &
&      matrix2%buffer_real,1,1,matrix2%descript%tab, &
&      ZERO,ZERO,il,iu,ABSTOL,&
&      m,nz,eigen,ORFAC, &
&      results%buffer_real,1,1,results%descript%tab, &
&      RWORK_tmp,-1,IWORK_tmp,-1,&
&      IFAIL,ICLUSTR,GAP,INFO)
  endif

  if (INFO/=0) then
     write(msg,'(A,I6)') "Problem to compute workspace to use ScaLAPACK, INFO=",INFO
     ABI_ERROR(msg)
  endif

  TWORK_tmp(1) = IWORK_tmp(1)
  TWORK_tmp(2) = INT(RWORK_tmp(1)) + matrix1%sizeb_global(2) *(matrix1%sizeb_global(2)-1)
  TWORK_tmp(3) = INT(real(CWORK_tmp(1)))

 ! Get the maximum of sizes of the work arrays processor%comm
  call MPI_ALLREDUCE(TWORK_tmp,TWORK,3,MPI_integer,MPI_MAX,comm,ierr)

  LIWORK = TWORK(1)
  LRWORK = TWORK(2)
  LCWORK = TWORK(3)

 ! Allocate the work arrays
  if (LIWORK>0) then
    ABI_MALLOC(IWORK,(LIWORK))
    IWORK(:) = 0
  else
    ABI_MALLOC(IWORK,(1))
  end if
  if (LRWORK>0) then
    ABI_MALLOC(RWORK,(LRWORK))
    RWORK(:) = 0._dp
  else
    ABI_MALLOC(RWORK,(1))
  end if
  if (LCWORK>0) then
    ABI_MALLOC(CWORK,(LCWORK))
    CWORK(:) = (0._dp,0._dp)
  else
    ABI_MALLOC(CWORK,(1))
  end if

  ! Call the calculation routine
  if (istwf_k/=2) then
     ! write(std_out,*) 'I am using PZHEGVX'
     call PZHEGVX(1,'V',range,'U',&
&      matrix1%sizeb_global(2),&
&      matrix1%buffer_cplx,1,1,matrix1%descript%tab, &
&      matrix2%buffer_cplx,1,1,matrix2%descript%tab, &
&      ZERO,ZERO,il,iu,ABSTOL,&
&      m,nz,eigen,ORFAC, &
&      results%buffer_cplx,1,1,results%descript%tab, &
&      CWORK,LCWORK,RWORK,LRWORK,IWORK,LIWORK,&
&      IFAIL,ICLUSTR,GAP,INFO)
  else
     ! write(std_out,*) 'I am using PDSYGVX'
     call PDSYGVX(1,'V',range,'U',&
&      matrix1%sizeb_global(2),&
&      matrix1%buffer_real,1,1,matrix1%descript%tab, &
&      matrix2%buffer_real,1,1,matrix2%descript%tab, &
&      ZERO,ZERO,il,iu,ABSTOL,&
&      m,nz,eigen,ORFAC, &
&      results%buffer_real,1,1,results%descript%tab, &
&      RWORK,LRWORK,IWORK,LIWORK,&
&      IFAIL,ICLUSTR,GAP,INFO)
  endif

  if (INFO/=0) then
     write(msg,'(A,I6)') "Problem to compute eigen problem with ScaLAPACK, INFO=",INFO
     ABI_ERROR(msg)
  endif

  ABI_FREE(IFAIl)
  ABI_FREE(ICLUSTR)
  ABI_FREE(GAP)
  ABI_SFREE(IWORK)
  ABI_SFREE(RWORK)
  ABI_SFREE(CWORK)
#endif
#endif
  return

end subroutine compute_generalized_eigen_problem
!!***

!----------------------------------------------------------------------

!!****f* m_slk/compute_eigen1
!! NAME
!!  compute_eigen1
!!
!! FUNCTION
!!  Calculation of eigenvalues and eigenvectors. complex and real cases.
!!
!! INPUTS
!!  comm= MPI communicator
!!  cplex=1 if matrix is real, 2 if complex
!!  nbli_global number of lines
!!  nbco_global number of columns
!!  matrix= the matrix to process
!!  vector= eigenvalues of the matrix
!!  istwf_k= 2 if we have a real matrix else complex.
!!  [use_gpu_elpa]= Flag to activate the use of GPU (ELPA only)
!!
!! OUTPUT
!!  vector
!!
!! SIDE EFFECTS
!!  results= ScaLAPACK matrix coming out of the operation
!!  eigen= eigenvalues of the matrix
!!
!! SOURCE

subroutine compute_eigen1(comm,processor,cplex,nbli_global,nbco_global,matrix,vector,istwf_k,&
&                         use_gpu_elpa) ! Optional argument

!Arguments ------------------------------------
!scalaras
 integer,intent(in) :: comm
 integer,intent(in) :: cplex,nbli_global,nbco_global
 integer,intent(in) :: istwf_k
 class(processor_scalapack),intent(in) :: processor
 integer,intent(in),optional :: use_gpu_elpa
!arrays
 real(dp),intent(inout) :: matrix(cplex*nbli_global,nbco_global)
 real(dp),intent(inout) :: vector(:)

!Local variables-------------------------------
#ifdef HAVE_LINALG_ELPA
 integer :: i,j
#endif
 integer :: ierr,use_gpu_elpa_
 type(matrix_scalapack) :: sca_matrix1
 type(matrix_scalapack) :: sca_matrix2
 real(dp),allocatable :: r_tmp_evec(:,:)
 complex(dpc),allocatable :: z_tmp_evec(:,:)

! *************************************************************************

 use_gpu_elpa_=0
#ifdef HAVE_LINALG_ELPA
 if (present(use_gpu_elpa)) use_gpu_elpa_=use_gpu_elpa
#endif

 ! ================================
 ! INITIALISATION SCALAPACK MATRIX
 ! ================================
 call sca_matrix1%init(nbli_global,nbco_global,processor,istwf_k)
 call sca_matrix2%init(nbli_global,nbco_global,processor,istwf_k)

 ! ==============================
 ! FILLING SCALAPACK MATRIX
 ! ==============================
 if ( istwf_k /= 2 ) then
   ABI_CHECK_IEQ(cplex, 2, "cplex != 2")
   ABI_MALLOC(z_tmp_evec,(nbli_global,nbco_global))
   z_tmp_evec=cmplx(0._DP,0._DP)
#ifdef HAVE_LINALG_ELPA
   ! The full matrix must be set (not only one half like in scalapack).
   do j=1,nbco_global
      do i=j+1,nbli_global
         matrix(2*(i-1)+1,j) = matrix(2*(j-1)+1,i)
         matrix(2*(i-1)+2,j) = -matrix(2*(j-1)+2,i)
      end do
   end do
#endif
   call matrix_from_complexmatrix(sca_matrix1,matrix,istwf_k)
 else
   ABI_CHECK_IEQ(cplex, 1, "cplex != 2")
   ABI_MALLOC(r_tmp_evec,(nbli_global,nbco_global))
   r_tmp_evec(:,:)=0._DP
#ifdef HAVE_LINALG_ELPA
   ! The full matrix must be set (not only one half like in scalapack).
   do j=1,nbco_global
      do i=j+1,nbli_global
         matrix(i,j) = matrix(j,i)
      end do
   end do
#endif
   call matrix_from_realmatrix(sca_matrix1,matrix,istwf_k)
 endif

 ! ================================
 ! COMPUTE EIGEN VALUES AND VECTORS : A * X = lambda  * X
 ! ================================
 call compute_eigen_problem(processor,sca_matrix1, sca_matrix2,vector, comm,istwf_k, &
&                           use_gpu_elpa=use_gpu_elpa_)

 ! ==============================
 ! CONCATENATE EIGEN VECTORS
 ! ==============================
#ifdef HAVE_MPI
 if (istwf_k /= 2) then
   call matrix_to_complexmatrix(sca_matrix2,z_tmp_evec,istwf_k)
   call MPI_ALLREDUCE(z_tmp_evec, matrix, nbli_global*nbco_global, MPI_DOUBLE_complex, MPI_SUM,comm,ierr)
 else
   call matrix_to_realmatrix(sca_matrix2,r_tmp_evec,istwf_k)
   call MPI_ALLREDUCE(r_tmp_evec, matrix, nbli_global*nbco_global, MPI_DOUBLE_PRECISION, MPI_SUM,comm,ierr)
 endif
#endif

 ! ====================================
 ! DESTRUCTION SCALAPACK AND TMP MATRICES
 ! ====================================
 call sca_matrix1%free()
 call sca_matrix2%free()

 ABI_SFREE(z_tmp_evec)
 ABI_SFREE(r_tmp_evec)

#ifndef HAVE_LINALG_ELPA
 ABI_UNUSED(use_gpu_elpa)
#endif

end subroutine compute_eigen1
!!***

!----------------------------------------------------------------------

!!****f* m_slk/compute_eigen2
!! NAME
!!  compute_eigen2
!!
!! FUNCTION
!!  Calculation of eigenvalues and eigenvectors: A * X = lambda * B * X
!!  complex and real cases.
!!
!! INPUTS
!!  comm= MPI communicator
!!  cplex=1 if matrix is real, 2 if complex
!!  nbli_global number of lines
!!  nbco_global number of columns
!!  matrix1= first ScaLAPACK matrix (matrix A)
!!  matrix2= second ScaLAPACK matrix (matrix B)
!!  vector=
!!  istwf_k= 2 if we have a real matrix else complex.
!!  [use_gpu_elpa]= Flag to activate the use of GPU (ELPA only)
!!
!! SIDE EFFECTS
!!  results= ScaLAPACK matrix coming out of the operation
!!  eigen= eigenvalues of the matrix
!!
!! SOURCE

subroutine compute_eigen2(comm,processor,cplex,nbli_global,nbco_global,matrix1,matrix2,vector,istwf_k, &
&                         use_gpu_elpa) ! Optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nbli_global,nbco_global
 integer,intent(in) :: comm
 integer,intent(in) :: istwf_k
 class(processor_scalapack),intent(in) :: processor
 integer,optional,intent(in) :: use_gpu_elpa
!arrays
 real(dp),intent(inout) :: matrix1(cplex*nbli_global,nbco_global)
 real(dp),intent(inout) :: matrix2(cplex*nbli_global,nbco_global)
 real(dp),intent(inout) :: vector(:)

!Local variables-------------------------------
#ifdef HAVE_LINALG_ELPA
 integer :: i,j
#endif
 integer :: ierr,use_gpu_elpa_
 type(matrix_scalapack) :: sca_matrix1, sca_matrix2, sca_matrix3
 real(dp),allocatable :: r_tmp_evec(:,:)
 complex(dpc),allocatable :: z_tmp_evec(:,:)

! *************************************************************************

 use_gpu_elpa_=0
#if defined HAVE_LINALG_ELPA
 if (present(use_gpu_elpa)) use_gpu_elpa_=use_gpu_elpa
#endif

 ! ================================
 ! INITIALISATION SCALAPACK MATRIX
 ! ================================
 call sca_matrix1%init(nbli_global,nbco_global,processor,istwf_k)
 call sca_matrix2%init(nbli_global,nbco_global,processor,istwf_k)
 call sca_matrix3%init(nbli_global,nbco_global,processor,istwf_k)

 ! ==============================
 ! FILLING SCALAPACK MATRIX
 ! ==============================
 if ( istwf_k /= 2 ) then
   ABI_CHECK_IEQ(cplex, 2, "cplex != 2")
   ABI_MALLOC(z_tmp_evec,(nbli_global,nbco_global))
   z_tmp_evec=cmplx(0._DP,0._DP)
#ifdef HAVE_LINALG_ELPA
   ! The full matrix must be set (not only one half like in scalapack).
   do j=1,nbco_global
      do i=j+1,nbli_global
         matrix1(2*(i-1)+1,j) = matrix1(2*(j-1)+1,i)
         matrix1(2*(i-1)+2,j) = -matrix1(2*(j-1)+2,i)
         matrix2(2*(i-1)+1,j) = matrix2(2*(j-1)+1,i)
         matrix2(2*(i-1)+2,j) = -matrix2(2*(j-1)+2,i)
      end do
   end do
#endif
   call matrix_from_complexmatrix(sca_matrix1,matrix1,istwf_k)
   call matrix_from_complexmatrix(sca_matrix2,matrix2,istwf_k)
 else
   ABI_CHECK_IEQ(cplex, 1, "cplex != 1")
   ABI_MALLOC(r_tmp_evec,(nbli_global,nbco_global))
   r_tmp_evec(:,:)=0._DP
#ifdef HAVE_LINALG_ELPA
   ! The full matrix must be set (not only one half like in scalapack).
   do j=1,nbco_global
      do i=j+1,nbli_global
         matrix1(i,j) = matrix1(j,i)
         matrix2(i,j) = matrix2(j,i)
      end do
   end do
#endif
   call matrix_from_realmatrix(sca_matrix1,matrix1,istwf_k)
   call matrix_from_realmatrix(sca_matrix2,matrix2,istwf_k)
 endif

 ! ================================
 ! COMPUTE EIGEN VALUES AND VECTORS : A * X = lambda * B * X
 ! ================================
 call compute_generalized_eigen_problem(processor,sca_matrix1,sca_matrix2,&
&             sca_matrix3,vector,comm,istwf_k,use_gpu_elpa=use_gpu_elpa_)

 ! ==============================
 ! CONCATENATE EIGEN VECTORS
 ! ==============================
#ifdef HAVE_MPI
 if ( istwf_k /= 2 ) then
   call matrix_to_complexmatrix(sca_matrix3,z_tmp_evec,istwf_k)
   call MPI_ALLREDUCE(z_tmp_evec, matrix1, nbli_global*nbco_global, MPI_DOUBLE_complex,&
&    MPI_SUM,comm,ierr)
 else
   call matrix_to_realmatrix(sca_matrix3,r_tmp_evec,istwf_k)
   call MPI_ALLREDUCE(r_tmp_evec, matrix1, nbli_global*nbco_global, MPI_DOUBLE_PRECISION,&
&    MPI_SUM,comm,ierr)
 endif
#endif

 ! ====================================
 ! DESTRUCTION SCALAPACK AND TMP MATRICES
 ! ====================================
 call sca_matrix1%free()
 call sca_matrix2%free()
 call sca_matrix3%free()


#ifndef HAVE_LINALG_ELPA
 ABI_UNUSED(use_gpu_elpa)
#endif

end subroutine compute_eigen2
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_heev
!! NAME
!! slk_heev
!!
!! FUNCTION
!!  slk_heev computes selected eigenvalues and, optionally, eigenvectors of an Hermitian matrix A.
!!   A * X = lambda * X
!!
!! INPUTS
!!  JOBZ    (global input) CHARACTER*1
!!          Specifies whether or not to compute the eigenvectors:
!!          = "N":  Compute eigenvalues only.
!!          = "V":  Compute eigenvalues and eigenvectors.
!!  UPLO    (global input) CHARACTER*1
!!          Specifies whether the upper or lower triangular part of the symmetric matrix A is stored:
!!          = "U":  Upper triangular
!!          = "L":  Lower triangular
!!
!!  mat=The object storing the local buffer in DOUBLE PRECISION, the array descriptor, the PBLAS context.
!!  vec=The distributed eigenvectors. Not referenced if JOBZ="N"
!!
!! OUTPUT
!!  W       (global output) array, dimension (N) where N is the rank of the global matrix.
!!          On normal exit, the first M entries contain the selected eigenvalues in ascending order.
!!
!! SIDE EFFECTS
!!  If JOBZ="V", the local buffer vec%buffer_cplx will contain part of the distributed eigenvectors.
!!  On exit, the lower triangle (if UPLO='L') or the upper triangle (if UPLO='U') of A, including the diagonal, is destroyed.
!!
!! SOURCE

subroutine slk_heev(mat, jobz, uplo, vec, w, &
                      mat_size, ija, ijz) ! Optional

!Arguments ------------------------------------
!scalars
 class(matrix_scalapack),intent(inout) :: mat
 character(len=*),intent(in) :: jobz, uplo
 class(matrix_scalapack),intent(inout) :: vec
!arrays
 real(dp),intent(out) :: w(:)
 integer,optional,intent(in) :: mat_size, ija(2), ijz(2)

#ifdef HAVE_LINALG_SCALAPACK
!Local variables ------------------------------
!scalars
 integer :: lwork, lrwork, info, nn
 !character(len=500) :: msg
!arrays
 integer :: ija__(2), ijz__(2)
 real(dp),allocatable :: rwork_dp(:)
 complex(dp),allocatable :: work_dp(:)

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "buffer_cplx not allocated")

 nn = mat%sizeb_global(2); if (present(mat_size)) nn = mat_size
 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijz__ = [1, 1]; if (present(ijz)) ijz__ = ijz

 ! Get optimal size of workspace.
 lwork = - 1; lrwork = -1
 ABI_MALLOC(work_dp, (1))
 ABI_MALLOC(rwork_dp, (1))

 !call pzheev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, info)

 call PZHEEV(jobz, uplo, nn, mat%buffer_cplx, ija__(1), ija__(2), mat%descript%tab, &
             w, vec%buffer_cplx, ijz__(1), ijz__(2), vec%descript%tab, work_dp, lwork, rwork_dp, lrwork, info)
 ABI_CHECK(info == 0, sjoin("Error in the calculation of the workspace size, info:", itoa(info)))

 lwork = NINT(real(work_dp(1))); lrwork= NINT(rwork_dp(1)) !*2
 ABI_FREE(work_dp)
 ABI_FREE(rwork_dp)

 ! MG: Nov 23 2011. On my mac with the official scalapack package, rwork(1) is not large enough and causes a SIGFAULT.
 if (firstchar(jobz, ['V'])) then
   if (lrwork < 2*nn + 2*nn-2) lrwork = 2*nn + 2*nn-2
 else if (firstchar(jobz, ['N'])) then
   if (lrwork < 2*nn) lrwork = 2*nn
 end if
 !write(std_out,*)lwork,lrwork

 ! Solve the problem.
 ABI_MALLOC(work_dp, (lwork))
 ABI_MALLOC(rwork_dp, (lrwork))

 call PZHEEV(jobz, uplo, nn, mat%buffer_cplx, ija__(1), ija__(2), mat%descript%tab, &
             w, vec%buffer_cplx, ijz__(1), ijz__(2), vec%descript%tab, work_dp, lwork, rwork_dp, lrwork, info)
 ABI_CHECK(info == 0, sjoin("PZHEEV returned info:", itoa(info)))
 ABI_FREE(work_dp)
 ABI_FREE(rwork_dp)
#endif

end subroutine slk_heev
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_sp_heev
!! NAME
!! slkmat_sp_heev
!!
!! FUNCTION
!!  slkmat_sp_heev computes selected eigenvalues and, optionally, eigenvectors of an Hermitian matrix A.
!!   A * X = lambda * X
!!
!! INPUTS
!!  JOBZ    (global input) CHARACTER*1
!!          Specifies whether or not to compute the eigenvectors:
!!          = "N":  Compute eigenvalues only.
!!          = "V":  Compute eigenvalues and eigenvectors.
!!  UPLO    (global input) CHARACTER*1
!!          Specifies whether the upper or lower triangular part of the symmetric matrix A is stored:
!!          = "U":  Upper triangular
!!          = "L":  Lower triangular
!!
!!  mat=The object storing the local buffer in SINGLE PRECISION, the array descriptor, the PBLAS context.
!!  vec=The distributed eigenvectors. Not referenced if JOBZ="N"
!!
!! OUTPUT
!!  W       (global output) array, dimension (N) where N is the rank of the global matrix.
!!          On normal exit, the first M entries contain the selected eigenvalues in ascending order.
!!
!! SIDE EFFECTS
!!  If JOBZ="V", the local buffer vec%buffer_cplx will contain part of the distributed eigenvectors.
!!  On exit, the lower triangle (if UPLO='L') or the upper triangle (if UPLO='U') of A, including the diagonal, is destroyed.
!!
!! SOURCE

subroutine slkmat_sp_heev(mat, jobz, uplo, vec, w, &
                          mat_size, ija, ijz) ! Optional

!Arguments ------------------------------------
!scalars
 class(slkmat_sp_t),intent(inout) :: mat
 character(len=*),intent(in) :: jobz, uplo
 class(slkmat_sp_t),intent(inout) :: vec
!arrays
 real(sp),intent(out) :: w(:)
 integer,optional,intent(in) :: mat_size, ija(2), ijz(2)

#ifdef HAVE_LINALG_SCALAPACK
!Local variables ------------------------------
!scalars
 integer :: lwork, lrwork, info, nn
!arrays
 integer :: ija__(2), ijz__(2)
 real(sp),allocatable :: rwork_sp(:)
 complex(sp),allocatable :: work_sp(:)

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "buffer_cplx not allocated")

 nn = mat%sizeb_global(2); if (present(mat_size)) nn = mat_size
 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijz__ = [1, 1]; if (present(ijz)) ijz__ = ijz

 ! Get optimal size of workspace.
 lwork = - 1; lrwork = -1
 ABI_MALLOC(work_sp, (1))
 ABI_MALLOC(rwork_sp, (1))

 !call pzheev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, info)

 call PCHEEV(jobz, uplo, nn, mat%buffer_cplx, ija__(1), ija__(2), mat%descript%tab, &
             w, vec%buffer_cplx, ijz__(1), ijz__(2), vec%descript%tab, work_sp, lwork, rwork_sp, lrwork, info)
 ABI_CHECK(info == 0, sjoin("Error in the calculation of the workspace size, info:", itoa(info)))

 lwork = NINT(real(work_sp(1))); lrwork= NINT(rwork_sp(1)) !*2
 ABI_FREE(work_sp)
 ABI_FREE(rwork_sp)

 ! MG: Nov 23 2011. On my mac with the official scalapack package, rwork(1) is not large enough and causes a SIGFAULT.
 if (firstchar(jobz, ['V'])) then
   if (lrwork < 2*nn + 2*nn-2) lrwork = 2*nn + 2*nn-2
 else if (firstchar(jobz, ['N'])) then
   if (lrwork < 2*nn) lrwork = 2*nn
 end if
 !write(std_out,*)lwork,lrwork

 ! Solve the problem.
 ABI_MALLOC(work_sp, (lwork))
 ABI_MALLOC(rwork_sp, (lrwork))

 call PCHEEV(jobz, uplo, nn, mat%buffer_cplx, ija__(1), ija__(2), mat%descript%tab, &
             w, vec%buffer_cplx, ijz__(1), ijz__(2), vec%descript%tab, work_sp, lwork, rwork_sp, lrwork, info)
 ABI_CHECK(info == 0, sjoin("PCHEEV returned info:", itoa(info)))
 ABI_FREE(work_sp)
 ABI_FREE(rwork_sp)
#endif

end subroutine slkmat_sp_heev
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_pzheevx
!! NAME
!!  slk_pzheevx
!!
!! FUNCTION
!!  slk_pzheevx computes selected eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.
!!  A * X = lambda * X
!!
!! INPUTS
!!  mat=ScaLAPACK matrix (matrix A)
!!
!!  vec=The distributed eigenvectors X. Not referenced if JOBZ="N"
!!
!!  JOBZ    (global input) CHARACTER*1
!!          Specifies whether or not to compute the eigenvectors:
!!          = "N":  Compute eigenvalues only.
!!          = "V":  Compute eigenvalues and eigenvectors.
!!
!!  RANGE   (global input) CHARACTER*1
!!          = "A": all eigenvalues will be found.
!!          = "V": all eigenvalues in the interval [VL,VU] will be found.
!!          = "I": the IL-th through IU-th eigenvalues will be found.
!!
!!  UPLO    (global input) CHARACTER*1
!!          Specifies whether the upper or lower triangular part of the Hermitian matrix A is stored:
!!          = "U":  Upper triangular
!!          = "L":  Lower triangular
!!
!!  VL      (global input) DOUBLE PRECISION
!!          If  RANGE="V",the lower bound of the interval to be searched for eigenvalues. Not referenced if RANGE =
!!          "A" or "I"
!!
!!  VU      (global input) DOUBLE PRECISION
!!          If RANGE="V", the upper bound of the interval to be searched for eigenvalues.  Not referenced  if  RANGE  =
!!          "A" or "I".
!!
!!  IL     (global input) integer
!!         If  RANGE="I",  the  index  (from smallest to largest) of the smallest eigenvalue to be returned.  IL >= 1.
!!         Not referenced if RANGE = "A" or "V".
!!
!!  IU     (global input) integer
!!         If RANGE="I", the index (from smallest to largest) of the largest eigenvalue to be returned.  min(IL,N)  <=
!!         IU <= N.  Not referenced if RANGE = "A" or "V"
!!
!!  ABSTOL  (global input) DOUBLE PRECISION
!!          If JOBZ="V", setting ABSTOL to PDLAMCH( CONTEXT, "U") yields the most orthogonal eigenvectors.
!!          The  absolute error tolerance for the eigenvalues.  An approximate eigenvalue is accepted as converged when
!!          it is determined to lie in an interval [a,b] of width less than or equal to
!!
!!           ABSTOL + EPS *   max( |a|,|b| ) ,
!!
!!          where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then EPS*norm(T) will be used
!!          in  its  place, where norm(T) is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.
!!          Eigenvalues will be computed  most  accurately  when  ABSTOL  is  set  to  twice  the  underflow  threshold
!!          2*PDLAMCH("S")  not  zero.   If  this routine returns with ((MOD(INFO,2).NE.0) .OR.  (MOD(INFO/8,2).NE.0)),
!!          indicating that some eigenvalues or eigenvectors did not converge, try setting ABSTOL to 2*PDLAMCH("S").
!!
!! OUTPUT
!!  mene_found= (global output) Total number of eigenvalues found.  0 <= mene_found <= N.
!!  eigen(N)= (global output) Eigenvalues of A where N is the dimension of M
!!            On normal exit, the first mene_found entries contain the selected eigenvalues in ascending order.
!!
!! SIDE EFFECTS
!!  If JOBZ="V", the local buffer vec%buffer_cplx will contain part of the distributed eigenvectors.
!!  Slk%mat%buffer_cplx is destroyed when the routine returns
!!
!! SOURCE

subroutine slk_pzheevx(mat, jobz, range, uplo, vl, vu, il, iu, abstol, vec, mene_found, eigen)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: mat
 integer,intent(in) :: il, iu
 integer,intent(out) :: mene_found
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo
 class(matrix_scalapack),intent(inout) :: vec
!arrays
 real(dp),intent(out) :: eigen(*)

#ifdef HAVE_LINALG_SCALAPACK
!Local variables-------------------------------
!scalars
 integer  :: lwork,lrwork,liwork,info,nvec_calc !,ierr
 real(dp) :: orfac
 character(len=500) :: msg
!arrays
 !integer :: ibuff(3),max_ibuff(3)
 integer,allocatable  :: iwork(:),iclustr(:),ifail(:)
 real(dp),allocatable  :: rwork(:),gap(:)
 complex(dpc),allocatable :: work(:)

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "buffer_cplx is not allocated!")

 ! abstol = PDLAMCH(vec%processor%grid%ictxt,'U')

 orfac  = -one ! Only for eigenvectors: use default value 10d-3.
 ! Vectors within orfac*norm(A) will be reorthogonalized.

 ! Allocate the arrays for the results of the calculation
 ABI_MALLOC(gap, (mat%processor%grid%dims(1) * mat%processor%grid%dims(2)))

 if (firstchar(jobz, ["V","v"])) then
   ABI_MALLOC(ifail, (mat%sizeb_global(2)))
   ABI_MALLOC(iclustr, (2*mat%processor%grid%dims(1) * mat%processor%grid%dims(2)))
 end if

 ! Get the optimal size of the work arrays.
 lwork=-1; lrwork=-1; liwork=-1
 ABI_MALLOC(work, (1))
 ABI_MALLOC(iwork, (1))
 ABI_MALLOC(rwork, (3))
 ! This is clearly seen in the source in which rwork(1:3) is accessed
 ! in the calculation of the workspace size.

  ! prototype
  !call pzheevx(jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w,
  !             orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)

  call PZHEEVX(jobz,range,uplo, mat%sizeb_global(2),mat%buffer_cplx,1,1,mat%descript%tab,&
    vl,vu,il,iu,abstol,mene_found,nvec_calc,eigen,orfac,&
    vec%buffer_cplx,1,1,vec%descript%tab,&
    work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)

  ABI_CHECK(info == 0, sjoin("Problem to compute workspace, info:", itoa(info)))

  lwork  = NINT(real(work(1)),kind=dp)
  lrwork = NINT(rwork(1))
  liwork = iwork(1)

  ABI_FREE(work)
  ABI_FREE(rwork)
  ABI_FREE(iwork)
  !
  ! FROM THE SCALAPACK MAN PAGE:
  ! The computed eigenvectors may not be orthogonal if the minimal workspace is supplied and ORFAC is too
  ! small. If you  want to guarantee orthogonality (at the cost of potentially poor performance) you should
  ! add the following to LRWORK: (CLUSTERSIZE-1)*N where CLUSTERSIZE is  the  number  of  eigenvalues  in  the
  ! largest cluster, where a cluster is defined as a set of close eigenvalues: { W(K),...,W(K+CLUSTERSIZE-1) |
  ! W(J+1) <= W(J) + ORFAC*2*norm(A) }.

  if (firstchar(jobz, ["V","v"])) then
    lrwork = INT( lrwork + mat%sizeb_global(2) *(mat%sizeb_global(2)-1) )
  end if

  ! ibuff(1) = lwork
  ! ibuff(2) = lrwork !INT(lrwork + mat%sizeb_global(2) *(mat%sizeb_global(2)-1)
  ! ibuff(3) = liwork

  ! Get the maximum of sizes of the work arrays processor%comm
  ! call MPI_ALLREDUCE(ibuff,max_ibuff,3,MPI_integer,MPI_MAX,comm,ierr)

  ! lwork  = max_ibuff(1)
  ! lrwork = max_ibuff(2)
  ! liwork = max_ibuff(3)

  ABI_MALLOC(work , (lwork ))
  ABI_MALLOC(rwork, (lrwork))
  ABI_MALLOC(iwork, (liwork))

  ! prototype
  !call pzheevx(jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w,
  !             orfac, z, iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)

  ! Call the scaLAPACK routine.
  ! write(std_out,*) 'I am using PZHEEVX'
  call PZHEEVX(jobz,range,uplo, mat%sizeb_global(2),mat%buffer_cplx,1,1,mat%descript%tab,&
    vl,vu,il,iu,abstol,mene_found,nvec_calc, eigen,orfac,&
    vec%buffer_cplx,1,1,vec%descript%tab,&
    work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)

 ! TODO
 !call pxheevx_info_to_msg(info, jobz, nvec_calc, mene_found, msg)
 !ABI_CHECK(info == 0, msg)

 ! Handle possible error.
 if (info < 0) then
   write(msg,'(a,i0,a)')" The ", -info, "-th argument of P?HEEVX had an illegal value."
   if (info == -25) msg = " LRWORK is too small to compute all the eigenvectors requested, no computation is performed"
   ABI_ERROR(msg)
 end if

 if (info > 0) then
   write(msg,'(a,i0)') " P?HEEVX returned info: ",info
   call wrtout(std_out, msg)
   if (MOD(info, 2) /= 0)then
     write(msg,'(3a)')&
     " One or more eigenvectors failed to converge. ",ch10,&
     " Their indices are stored in IFAIL. Ensure ABSTOL=2.0*PDLAMCH('U')"
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 2, 2) /= 0) then
     write(msg,'(5a)')&
     " Eigenvectors corresponding to one or more clusters of eigenvalues ",ch10,&
     " could not be reorthogonalized because of insufficient workspace. ",ch10,&
     " The indices of the clusters are stored in the array ICLUSTR."
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 4, 2) /= 0) then
     write(msg,'(3a)')" Space limit prevented PZHEEVX from computing all of the eigenvectors between VL and VU. ",ch10,&
      " The number of eigenvectors computed is returned in NZ."
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 8, 2) /= 0) then
     call wrtout(std_out, "PZSTEBZ failed to compute eigenvalues. Ensure ABSTOL=2.0*PDLAMCH('U')")
   end if
   ABI_ERROR("Cannot continue")
 end if

 ! Check the number of eigenvalues found wrt to the number of vectors calculated.
 if ( firstchar(jobz, ['V','v']) .and. mene_found /= nvec_calc) then
   write(msg,'(5a)')&
   " The user supplied insufficient space and PZHEEVX is not able to detect this before beginning computation. ",ch10,&
   " To get all the eigenvectors requested, the user must supply both sufficient space to hold the ",ch10,&
   " eigenvectors in Z (M .LE. DESCZ(N_)) and sufficient workspace to compute them. "
   !ierr = huge(1)
   ABI_ERROR(msg)
 end if

 ABI_FREE(work)
 ABI_FREE(rwork)
 ABI_FREE(iwork)
 ABI_FREE(gap)

 ABI_SFREE(ifail)
 ABI_SFREE(iclustr)
#endif

end subroutine slk_pzheevx
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_pzhegvx
!! NAME
!!  slk_pzhegvx
!!
!! FUNCTION
!!  slk_pzhegvx provides an object-oriented interface to the ScaLAPACK routine PZHEGVX that
!!  computes selected eigenvalues and, optionally, eigenvectors of a complex generalized
!!  Hermitian-definite eigenproblem, of the form
!!  sub( A )*x=(lambda)*sub( B )*x,  sub( A )*sub( B )x=(lambda)*x,  or sub( B )*sub( A )*x=(lambda)*x.
!!  Here sub( A ) denoting A( IA:IA+N-1, JA:JA+N-1 ) is assumed to be
!!  Hermitian, and sub( B ) denoting B( IB:IB+N-1, JB:JB+N-1 ) is assumed
!!  to be Hermitian positive definite.
!!
!! INPUTS
!!  Slk_matA<matrix_scalapack>=ScaLAPACK matrix (matrix A)
!!  Slk_matB<matrix_scalapack>=ScaLAPACK matrix (matrix B)
!!  Slk_vec<matrix_scalapack>=The distributed eigenvectors X. Not referenced if JOBZ="N"
!!
!!  IBtype   (global input) integer
!!          Specifies the problem type to be solved:
!!          = 1:  sub( A )*x = (lambda)*sub( B )*x
!!          = 2:  sub( A )*sub( B )*x = (lambda)*x
!!          = 3:  sub( B )*sub( A )*x = (lambda)*x
!!
!!  JOBZ    (global input) CHARACTER*1
!!          Specifies whether or not to compute the eigenvectors:
!!          = "N":  Compute eigenvalues only.
!!          = "V":  Compute eigenvalues and eigenvectors.
!!
!!  RANGE   (global input) CHARACTER*1
!!          = "A": all eigenvalues will be found.
!!          = "V": all eigenvalues in the interval [VL,VU] will be found.
!!          = "I": the IL-th through IU-th eigenvalues will be found.
!!
!!  UPLO    (global input) CHARACTER*1
!!          Specifies whether the upper or lower triangular part of the Hermitian matrix sub(A) and sub(B) is stored:
!!          = "U":  Upper triangular
!!          = "L":  Lower triangular
!!
!!  VL      (global input) DOUBLE PRECISION
!!          If  RANGE="V",the lower bound of the interval to be searched for eigenvalues. Not referenced if RANGE =
!!          "A" or "I"
!!
!!  VU      (global input) DOUBLE PRECISION
!!          If RANGE="V", the upper bound of the interval to be searched for eigenvalues.  Not referenced  if  RANGE  =
!!          "A" or "I".
!!
!!  IL     (global input) integer
!!         If  RANGE="I",  the  index  (from smallest to largest) of the smallest eigenvalue to be returned.  IL >= 1.
!!         Not referenced if RANGE = "A" or "V".
!!
!!  IU     (global input) integer
!!         If RANGE="I", the index (from smallest to largest) of the largest eigenvalue to be returned.  min(IL,N)  <=
!!         IU <= N.  Not referenced if RANGE = "A" or "V"
!!
!!  ABSTOL  (global input) DOUBLE PRECISION
!!          If JOBZ="V", setting ABSTOL to PDLAMCH( CONTEXT, "U") yields the most orthogonal eigenvectors.
!!          The  absolute error tolerance for the eigenvalues.  An approximate eigenvalue is accepted as converged when
!!          it is determined to lie in an interval [a,b] of width less than or equal to
!!
!!           ABSTOL + EPS *   max( |a|,|b| ) ,
!!
!!          where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then EPS*norm(T) will be used
!!          in  its  place, where norm(T) is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.
!!          Eigenvalues will be computed  most  accurately  when  ABSTOL  is  set  to  twice  the  underflow  threshold
!!          2*PDLAMCH("S")  not  zero.   If  this routine returns with ((MOD(INFO,2).NE.0) .OR.  (MOD(INFO/8,2).NE.0)),
!!          indicating that some eigenvalues or eigenvectors did not converge, try setting ABSTOL to 2*PDLAMCH("S").
!!
!! OUTPUT
!!  mene_found= (global output) Total number of eigenvalues found.  0 <= mene_found <= N.
!!
!!  eigen(N)= (global output) Eigenvalues of A where N is the dimension of M
!!            On normal exit, the first mene_found entries contain the selected eigenvalues in ascending order.
!!
!! SIDE EFFECTS
!!  Slk_vec<matrix_scalapack>:
!!   %buffer_cplx local output (global dimension (N,N)
!!     If JOBZ = 'V', then on normal exit the first M columns of Z
!!     contain the orthonormal eigenvectors of the matrix
!!     corresponding to the selected eigenvalues.
!!     If JOBZ = 'N', then Z is not referenced.
!!
!!  Slk_matA<matrix_scalapack>:
!!    %buffer_cplx
!!      (local input/local output) complex(DPC) pointer into the
!!      local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
!!      On entry, this array contains the local pieces of the
!!      N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
!!      the leading N-by-N upper triangular part of sub( A ) contains
!!      the upper triangular part of the matrix.  If UPLO = 'L', the
!!      leading N-by-N lower triangular part of sub( A ) contains
!!      the lower triangular part of the matrix.
!!
!!      On exit, if JOBZ = 'V', then if INFO = 0, sub( A ) contains
!!      the distributed matrix Z of eigenvectors.  The eigenvectors
!!      are normalized as follows:
!!      if IBtype = 1 or 2, Z**H*sub( B )*Z = I;
!!      if IBtype = 3, Z**H*inv( sub( B ) )*Z = I.
!!      If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
!!      or the lower triangle (if UPLO='L') of sub( A ), including
!!      the diagonal, is destroyed.
!!
!!  Slk_matB=
!!    %buffer_cplx
!!      (local input/local output) complex*(DPC) pointer into the
!!      local memory to an array of dimension (LLD_B, LOCc(JB+N-1)).
!!      On entry, this array contains the local pieces of the
!!      N-by-N Hermitian distributed matrix sub( B ). If UPLO = 'U',
!!      the leading N-by-N upper triangular part of sub( B ) contains
!!      the upper triangular part of the matrix.  If UPLO = 'L', the
!!      leading N-by-N lower triangular part of sub( B ) contains
!!      the lower triangular part of the matrix.
!!
!!      On exit, if INFO <= N, the part of sub( B ) containing the
!!      matrix is overwritten by the triangular factor U or L from
!!      the Cholesky factorization sub( B ) = U**H*U or
!!      sub( B ) = L*L**H.
!!
!! SOURCE

subroutine slk_pzhegvx(Slk_matA, ibtype, jobz, range, uplo, Slk_matB, vl, vu, il, iu, abstol, Slk_vec, mene_found, eigen)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: Slk_matA
 integer,intent(in) :: il,iu,ibtype
 integer,intent(out) :: mene_found
 real(dp),intent(in) :: abstol,vl,vu
 character(len=*),intent(in) :: jobz,range,uplo
 class(matrix_scalapack),intent(inout) :: Slk_matB
 class(matrix_scalapack),intent(inout) :: Slk_vec
!arrays
 real(dp),intent(out) :: eigen(*)

#ifdef HAVE_LINALG_SCALAPACK
!Local variables-------------------------------
!scalars
 integer  :: lwork,lrwork,liwork,info,nvec_calc !,ierr
 real(dp) :: orfac
 logical :: ltest
 character(len=500) :: msg
!arrays
 !integer :: ibuff(3),max_ibuff(3)
 integer :: desca(DLEN_),descb(DLEN_),descz(DLEN_)
 integer,allocatable  :: iwork(:),iclustr(:),ifail(:)
 real(dp),allocatable  :: rwork(:),gap(:)
 complex(dpc),allocatable :: work(:)

!************************************************************************

 ABI_CHECK(allocated(Slk_matA%buffer_cplx), "buffer_cplx is not allocated!")

 ! abstol = PDLAMCH(Slk_vecprocessor%grid%ictxt,'U')

 orfac  = -one ! Only for eigenvectors: use default value 10d-3.
 ! Vectors within orfac*norm(A) will be reorthogonalized.

 ! ======================
 ! Alignment requirements
 ! ======================
 ! The distributed submatrices A(IA:*, JA:*), C(IC:IC+M-1,JC:JC+N-1),
 ! and B( IB:IB+N-1, JB:JB+N-1 ) must verify some alignment properties,

 desca = Slk_matA%descript%tab
 descb = Slk_matB%descript%tab
 if (firstchar(jobz, ["V", "v"])) then
   descz = Slk_vec%descript%tab
 else
   descz = Slk_matA%descript%tab
 end if

 ltest = .TRUE.
 ltest = ltest .and. (DESCA(MB_) == DESCA(NB_))
 !IA = IB = IZ
 !JA = IB = JZ
 ltest = ltest .and.ALL(DESCA(M_   ) == [DESCB(M_   ), DESCZ(M_   )])
 ltest = ltest .and.ALL(DESCA(N_   ) == [DESCB(N_   ), DESCZ(N_   )])
 ltest = ltest .and.ALL(DESCA(MB_  ) == [DESCB(MB_  ), DESCZ(MB_  )])
 ltest = ltest .and.ALL(DESCA(NB_  ) == [DESCB(NB_  ), DESCZ(NB_  )])
 ltest = ltest .and.ALL(DESCA(RSRC_) == [DESCB(RSRC_), DESCZ(RSRC_)])
 ltest = ltest .and.ALL(DESCA(CSRC_) == [DESCB(CSRC_), DESCZ(CSRC_)])
 !MOD( IA-1, DESCA( MB_ ) ) = 0
 !MOD( JA-1, DESCA( NB_ ) ) = 0
 !MOD( IB-1, DESCB( MB_ ) ) = 0
 !MOD( JB-1, DESCB( NB_ ) ) = 0

 if (.not.ltest) then
   ABI_ERROR("Alignment requirements not satisfied, check the caller")
 end if

 !Allocate the arrays for the results of the calculation
 ABI_MALLOC(gap, (Slk_matA%processor%grid%dims(1) * Slk_matA%processor%grid%dims(2)))

 if (firstchar(jobz, ["V","v"])) then
   ABI_MALLOC(ifail,(Slk_matA%sizeb_global(2)))
   ABI_MALLOC(iclustr,( 2*Slk_matA%processor%grid%dims(1) * Slk_matA%processor%grid%dims(2)))
 else
   ABI_MALLOC(ifail,(1))
 end if

 ! Get the optimal size of the work arrays.
 lwork=-1; lrwork=-1; liwork=-1
 ABI_MALLOC(work,(1))
 ABI_MALLOC(iwork,(1))
 ABI_MALLOC(rwork,(3))
 ! This is clearly seen in the source in which rwork(1:3) is accessed
 ! in the calcuation of the workspace size.

 call pzhegvx(ibtype,jobz,range,uplo, Slk_matA%sizeb_global(2),Slk_matA%buffer_cplx,1,1,Slk_matA%descript%tab,&
   Slk_matB%buffer_cplx,1,1,Slk_matB%descript%tab,&
   vl,vu,il,iu,abstol,mene_found,nvec_calc,eigen,orfac,&
   Slk_vec%buffer_cplx,1,1,Slk_vec%descript%tab,&
   work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)

 ABI_CHECK(info == 0, sjoin("Problem to compute workspace, info:", itoa(info)))

 lwork  = NINT(real(work(1)),kind=dp)
 lrwork = NINT(rwork(1))
 liwork = iwork(1)

 ABI_FREE(work)
 ABI_FREE(rwork)
 ABI_FREE(iwork)

 !FROM THE SCALAPACK MAN PAGE:
 !The computed eigenvectors may not be orthogonal if the minimal workspace is supplied and ORFAC is too
 !small. If you  want to guarantee orthogonality (at the cost of potentially poor performance) you should
 !add the following to LRWORK: (CLUSTERSIZE-1)*N where CLUSTERSIZE is  the  number  of  eigenvalues  in  the
 !largest cluster, where a cluster is defined as a set of close eigenvalues: { W(K),...,W(K+CLUSTERSIZE-1) |
 !W(J+1) <= W(J) + ORFAC*2*norm(A) }.

 if (firstchar(jobz, ["V","v"])) then
   lrwork = INT( lrwork + Slk_matA%sizeb_global(2) *(Slk_matA%sizeb_global(2)-1) )
 end if

 !ibuff(1) = lwork
 !ibuff(2) = lrwork !INT(lrwork + Slk_matA%sizeb_global(2) *(Slk_matA%sizeb_global(2)-1)
 !ibuff(3) = liwork

 !Get the maximum of sizes of the work arrays processor%comm
 !call MPI_ALLREDUCE(ibuff,max_ibuff,3,MPI_integer,MPI_MAX,comm,ierr)

 !lwork  = max_ibuff(1)
 !lrwork = max_ibuff(2)
 !liwork = max_ibuff(3)

 ABI_MALLOC(work , (lwork ))
 ABI_MALLOC(rwork, (lrwork))
 ABI_MALLOC(iwork, (liwork))

 ! Call the scaLAPACK routine.
 ! write(std_out,*) 'I am using PZHEGVX'
 call pzhegvx(ibtype,jobz,range,uplo, Slk_matA%sizeb_global(2),Slk_matA%buffer_cplx,1,1,Slk_matA%descript%tab,&
    Slk_matB%buffer_cplx,1,1,Slk_matB%descript%tab,&
    vl,vu,il,iu,abstol,mene_found,nvec_calc, eigen,orfac,&
    Slk_vec%buffer_cplx,1,1,Slk_vec%descript%tab,&
   work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)

 ! Handle the possible error.
 if (info < 0) then
   write(msg,'(a,i7,a)')" The ",-info,"-th argument of PZHEGVX had an illegal value."
   if (info==-25) msg = " LRWORK is too small to compute all the eigenvectors requested, no computation is performed"
   ABI_ERROR(msg)
 end if

 if (info > 0) then
   write(msg,'(a,i0)') " PZHEGVX returned info: ",info
   call wrtout(std_out, msg)
   if (MOD(info,2)/=0)then
     write(msg,'(3a)')&
     " One or more eigenvectors failed to converge. ",ch10,&
     " Their indices are stored in IFAIL. Ensure ABSTOL=2.0*PDLAMCH('U')"
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 2, 2) /= 0) then
     write(msg,'(5a)')&
     " Eigenvectors corresponding to one or more clusters of eigenvalues ",ch10,&
     " could not be reorthogonalized because of insufficient workspace. ",ch10,&
     " The indices of the clusters are stored in the array ICLUSTR."
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 4, 2) /= 0) then
     write(msg,'(3a)')&
     " Space limit prevented PZHEGVX from computing all of the eigenvectors between VL and VU. ",ch10,&
     " The number of eigenvectors  computed  is returned in NZ."
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 8, 2) /= 0) then
     msg = " PZSTEBZ  failed to compute eigenvalues. Ensure ABSTOL=2.0*PDLAMCH('U')"
     call wrtout(std_out, msg)
   end if
   if (MOD(info / 16, 2) /= 0) then
     write(msg,'(3a)')&
     " B was not positive definite.",ch10,&
     " IFAIL(1) indicates the order of the smallest minor which is not positive definite."
     call wrtout(std_out, msg)
   end if
   ABI_ERROR("Cannot continue")
 end if

 ! Check the number of eigenvalues found wrt to the number of vectors calculated.
 if ( firstchar(jobz, ['V','v']) .and. mene_found/=nvec_calc) then
   write(msg,'(5a)')&
   " The user supplied insufficient space and PZHEGVX is not able to detect this before beginning computation. ",ch10,&
   " To get all the  eigenvectors requested, the user must supply both sufficient space to hold the ",ch10,&
   " eigenvectors in Z (M .LE. DESCZ(N_)) and sufficient workspace to compute them. "
   ABI_ERROR(msg)
 end if

 ABI_FREE(work)
 ABI_FREE(rwork)
 ABI_FREE(iwork)
 ABI_FREE(gap)
 ABI_FREE(ifail)
 ABI_SFREE(iclustr)
#endif

end subroutine slk_pzhegvx
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_invert
!! NAME
!! slk_invert
!!
!! FUNCTION
!!  Compute the inverse of a complex matrix.
!!
!! SIDE EFFECTS
!! mat
!!    In input, the matrix to invert.
!!    In output the matrix inverted and distributed among the nodes.
!!
!! SOURCE

subroutine slk_invert(mat)

!Arguments ------------------------------------
 class(basemat_t),intent(inout) :: mat

#ifdef HAVE_LINALG_SCALAPACK
!Local variables ------------------------------
!scalars
 integer :: lwork,info,ipiv_size,liwork
!array
 integer,allocatable :: ipiv(:), iwork(:)
 complex(dp),allocatable :: work_dp(:)
 complex(sp),allocatable :: work_sp(:)

!************************************************************************

 ! IMPORTANT NOTE: PZGETRF requires square block decomposition i.e., MB_A = NB_A.
 if (mat%descript%tab(MB_) /= mat%descript%tab(NB_)) then
   ABI_ERROR(" PZGETRF requires square block decomposition i.e MB_A = NB_A.")
 end if

 ipiv_size = my_locr(mat) + mat%descript%tab(MB_)
 ABI_MALLOC(ipiv, (ipiv_size))

 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_cplx)) then
     ! P * L * U  Factorization.
     call PZGETRF(mat%sizeb_global(1), mat%sizeb_global(2), mat%buffer_cplx, 1, 1, mat%descript%tab,ipiv, info)
     ABI_CHECK(info == 0, sjoin(" PZGETRF returned info:", itoa(info)))

     ! Get optimal size of workspace for PZGETRI.
     lwork = -1; liwork = -1
     ABI_MALLOC(work_dp,(1))
     ABI_MALLOC(iwork,(1))

     call PZGETRI(mat%sizeb_global(1), mat%buffer_cplx, 1, 1, mat%descript%tab, ipiv, work_dp, lwork, iwork, liwork, info)
     ABI_CHECK(info == 0, "PZGETRI: Error while computing workspace size")

     lwork = nint(real(work_dp(1))); liwork=iwork(1)
     ABI_FREE(work_dp)
     ABI_FREE(iwork)

     ! Solve the problem.
     ABI_MALLOC(work_dp, (lwork))
     ABI_MALLOC(iwork, (liwork))

     call PZGETRI(mat%sizeb_global(1), mat%buffer_cplx, 1, 1, mat%descript%tab, ipiv, work_dp, lwork, iwork, liwork, info)
     ABI_CHECK(info == 0, sjoin("PZGETRI returned info:", itoa(info)))
     ABI_FREE(work_dp)

   else if (allocated(mat%buffer_real)) then
     ABI_ERROR("Inversion for real matrices not coded!")
   end if

 class is (slkmat_sp_t)
   if (allocated(mat%buffer_cplx)) then
     ! P * L * U  Factorization.
     call PCGETRF(mat%sizeb_global(1), mat%sizeb_global(2), mat%buffer_cplx, 1, 1, mat%descript%tab,ipiv, info)
     ABI_CHECK(info == 0, sjoin(" PCGETRF returned info:", itoa(info)))

     ! Get optimal size of workspace for PZGETRI.
     lwork = -1; liwork = -1
     ABI_MALLOC(work_sp,(1))
     ABI_MALLOC(iwork,(1))

     call PCGETRI(mat%sizeb_global(1), mat%buffer_cplx, 1, 1, mat%descript%tab, ipiv, work_sp, lwork, iwork, liwork, info)
     ABI_CHECK(info == 0, "PZGETRI: Error while computing workspace size")

     lwork = nint(real(work_sp(1))); liwork=iwork(1)
     ABI_FREE(work_sp)
     ABI_FREE(iwork)

     ! Solve the problem.
     ABI_MALLOC(work_sp, (lwork))
     ABI_MALLOC(iwork, (liwork))

     call PCGETRI(mat%sizeb_global(1), mat%buffer_cplx, 1, 1, mat%descript%tab, ipiv, work_sp, lwork, iwork, liwork, info)
     ABI_CHECK(info == 0, sjoin("PZGETRI returned info:", itoa(info)))
     ABI_FREE(work_sp)

   else if (allocated(mat%buffer_real)) then
     ABI_ERROR("Inversion for real matrices not coded!")
   end if

 class default
   ABI_ERROR("Wrong class")
 end select

 ABI_FREE(iwork)
 ABI_FREE(ipiv)
#endif

end subroutine slk_invert
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_hpd_invert
!! NAME
!! slk_hpd_invert
!!
!! FUNCTION
!!  Compute the inverse of an Hermitian positive definite matrix.
!!
!! INPUTS
!!  uplo: global input
!!    = 'U':  Upper triangle of sub( A ) is stored;
!!    = 'L':  Lower triangle of sub( A ) is stored.
!!  [full]: If full PBLAS matrix is neeeded. Default: True
!!
!! SIDE EFFECTS
!!  mat= The object storing the local buffer, the array descriptor, the context, etc.
!!    On entry, this array contains the local pieces of the N-by-N Hermitian distributed matrix sub( A ) to be factored.
!!    If UPLO = 'U', the leading N-by-N upper triangular part of sub( A ) contains the upper triangular part of the matrix,
!!    and its strictly lower triangular part is not referenced.
!!    If UPLO = 'L', the leading N-by-N lower triangular part of sub( A ) contains the lower triangular part of the distribu-
!!    ted matrix, and its strictly upper triangular part is not referenced.
!!    On exit, the local pieces of the upper or lower triangle of the (Hermitian) inverse of sub( A )
!!
!! SOURCE

subroutine slk_hpd_invert(mat, uplo, full)

!Arguments ------------------------------------
 character(len=*),intent(in) :: uplo
 class(matrix_scalapack),intent(inout) :: mat
 logical,optional,intent(in) :: full

#ifdef HAVE_LINALG_SCALAPACK
!Local variables ------------------------------
!scalars
 integer :: info, mm, il1, il2, iglob1, iglob2
 type(matrix_scalapack) :: work_mat
 logical :: full__

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "buffer_cplx not allocated")

 ! ZPOTRF computes the Cholesky factorization of a complex Hermitian positive definite.
 !  A = U**H * U,   if UPLO = 'U', or
 !  A = L  * L**H,  if UPLO = 'L',
 mm = mat%sizeb_global(1)
 call PZPOTRF(uplo, mm, mat%buffer_cplx, 1, 1, mat%descript%tab, info)
 ABI_CHECK(info == 0, sjoin("PZPOTRF returned info:", itoa(info)))

 ! PZPOTRI computes the inverse of a complex Hermitian positive definite
 ! distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the
 ! Cholesky factorization sub( A ) = U**H*U or L*L**H computed by PZPOTRF.
 call PZPOTRI(uplo, mm, mat%buffer_cplx, 1, 1, mat%descript%tab, info)
 ABI_CHECK(info == 0, sjoin("PZPOTRI returned info:", itoa(info)))

 full__ = .True.; if (present(full)) full__ = full
 if (full__) then
   ! Only the uplo part contains the inverse so we need to fill the other triangular part.
   !     1)  Fill the missing triangle with zeros and copy results to work_mat
   !     2)  Call pzgeadd to compute: sub(C) := beta*sub(C) + alpha*op(sub(A))
   !     3)  Divide diagonal elements by two.

   do il2=1,mat%sizeb_local(2)
     iglob2 = mat%loc2gcol(il2)
     do il1=1,mat%sizeb_local(1)
       iglob1 = mat%loc2grow(il1)
       if (uplo == "L" .and. iglob2 > iglob1) mat%buffer_cplx(il1, il2) = zero
       if (uplo == "U" .and. iglob2 < iglob1) mat%buffer_cplx(il1, il2) = zero
     end do
   end do

   call mat%copy(work_mat, empty=.False.)

   ! call pzgeadd(trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
   ! sub(C) := beta*sub(C) + alpha*op(sub(A))
   call pzgeadd("C", mm, mm, cone, work_mat%buffer_cplx, 1, 1, work_mat%descript%tab, &
         cone, mat%buffer_cplx, 1, 1, mat%descript%tab)
   call work_mat%free()

   do il2=1,mat%sizeb_local(2)
     iglob2 = mat%loc2gcol(il2)
     do il1=1,mat%sizeb_local(1)
       iglob1 = mat%loc2grow(il1)
       if (iglob2 == iglob1) mat%buffer_cplx(il1, il2) = half * mat%buffer_cplx(il1, il2)
     end do
   end do
 end if ! full__
#endif

end subroutine slk_hpd_invert
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_sp_hpd_invert
!! NAME
!! slkmat_sp_hpd_invert
!!
!! FUNCTION
!!  Compute the inverse of an Hermitian positive definite matrix.
!!
!! INPUTS
!!  uplo: global input
!!    = 'U':  Upper triangle of sub( A ) is stored;
!!    = 'L':  Lower triangle of sub( A ) is stored.
!!  [full]: If full PBLAS matrix is neeeded. Default: True
!!
!! SIDE EFFECTS
!!  mat= The object storing the local buffer, the array descriptor, the context, etc.
!!    On entry, this array contains the local pieces of the N-by-N Hermitian distributed matrix sub( A ) to be factored.
!!    If UPLO = 'U', the leading N-by-N upper triangular part of sub( A ) contains the upper triangular part of the matrix,
!!    and its strictly lower triangular part is not referenced.
!!    If UPLO = 'L', the leading N-by-N lower triangular part of sub( A ) contains the lower triangular part of the distribu-
!!    ted matrix, and its strictly upper triangular part is not referenced.
!!    On exit, the local pieces of the upper or lower triangle of the (Hermitian) inverse of sub( A )
!!
!! SOURCE

subroutine slkmat_sp_hpd_invert(mat, uplo, full)

!Arguments ------------------------------------
 character(len=*),intent(in) :: uplo
 class(slkmat_sp_t),intent(inout) :: mat
 logical,optional,intent(in) :: full

#ifdef HAVE_LINALG_SCALAPACK
!Local variables ------------------------------
!scalars
 integer :: info, mm, il1, il2, iglob1, iglob2
 type(slkmat_sp_t) :: work_mat
 logical :: full__

!************************************************************************

 ABI_CHECK(allocated(mat%buffer_cplx), "buffer_cplx not allocated")

 ! ZPOTRF computes the Cholesky factorization of a complex Hermitian positive definite.
 !  A = U**H * U,   if UPLO = 'U', or
 !  A = L  * L**H,  if UPLO = 'L',
 mm = mat%sizeb_global(1)
 call PCPOTRF(uplo, mm, mat%buffer_cplx, 1, 1, mat%descript%tab, info)
 ABI_CHECK(info == 0, sjoin("PCPOTRF returned info:", itoa(info)))

 ! PZPOTRI computes the inverse of a complex Hermitian positive definite
 ! distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the
 ! Cholesky factorization sub( A ) = U**H*U or L*L**H computed by PZPOTRF.
 call PCPOTRI(uplo, mm, mat%buffer_cplx, 1, 1, mat%descript%tab, info)
 ABI_CHECK(info == 0, sjoin("PCPOTRI returned info:", itoa(info)))

 full__ = .True.; if (present(full)) full__ = full
 if (full__) then
   ! Only the uplo part contains the inverse so we need to fill the other triangular part.
   !     1)  Fill the missing triangle with zeros and copy results to work_mat
   !     2)  Call pzgeadd to compute: sub(C) := beta*sub(C) + alpha*op(sub(A))
   !     3)  Divide diagonal elements by two.

   do il2=1,mat%sizeb_local(2)
     iglob2 = mat%loc2gcol(il2)
     do il1=1,mat%sizeb_local(1)
       iglob1 = mat%loc2grow(il1)
       if (uplo == "L" .and. iglob2 > iglob1) mat%buffer_cplx(il1, il2) = zero_sp
       if (uplo == "U" .and. iglob2 < iglob1) mat%buffer_cplx(il1, il2) = zero_sp
     end do
   end do

   call mat%copy(work_mat, empty=.False.)

   ! call pzgeadd(trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
   ! sub(C) := beta*sub(C) + alpha*op(sub(A))
   call pcgeadd("C", mm, mm, cone_sp, work_mat%buffer_cplx, 1, 1, work_mat%descript%tab, &
         cone_sp, mat%buffer_cplx, 1, 1, mat%descript%tab)
   call work_mat%free()

   do il2=1,mat%sizeb_local(2)
     iglob2 = mat%loc2gcol(il2)
     do il1=1,mat%sizeb_local(1)
       iglob1 = mat%loc2grow(il1)
       if (iglob2 == iglob1) mat%buffer_cplx(il1, il2) = 0.5_sp * mat%buffer_cplx(il1, il2)
     end do
   end do
 end if ! full__
#endif

end subroutine slkmat_sp_hpd_invert
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_ptrans
!! NAME
!!  slk_ptrans
!!
!! FUNCTION
!! Transposes a matrix
!!
!!    sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
!!
!! where
!!
!!    sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
!!    sub( A ) denotes A(IA:IA+N-1,JA:JA+M-1), and, op( X ) = X'.
!!
!! Thus, op( sub( A ) ) denotes A(IA:IA+N-1,JA:JA+M-1)'.
!! Beta is a scalar, sub( C ) is an m by n submatrix, and sub( A ) is an n by m submatrix.
!!
!! INPUTS
!!  [ija(2)]: (global) The row and column indices in the distributed matrix in_mat indicating
!!       the first row and the first column of the submatrix sub(A), respectively.
!!  [ijc(2)]: (global) The row and column indices in the distributed matrix out_mat
!!   indicating the first row and the first column of the submatrix sub(C), respectively.
!!  [free]: True to deallocate in_mat. Default: False
!!
!! SOURCE

subroutine slk_ptrans(in_mat, trans, out_mat, &
                      out_gshape, ija, ijc, size_blocs, alpha, beta, free) ! optional

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: in_mat
 character(len=1),intent(in) :: trans
 class(matrix_scalapack),intent(inout) :: out_mat
 integer,optional,intent(in) :: out_gshape(2), size_blocs(2), ija(2), ijc(2)
 complex(dp),optional,intent(in) :: alpha, beta
 logical,optional,intent(in) :: free

!Local variables-------------------------------
 integer :: sb, mm, nn, size_blocs__(2)
 real(dp) :: ralpha__, rbeta__
 integer :: ija__(2), ijc__(2)
 complex(dp) :: calpha__, cbeta__

! *************************************************************************

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijc__ = [1, 1]; if (present(ijc)) ijc__ = ijc

 ! transposed output (sub)matrix has shape (nn, mm)
 if (present(out_gshape)) then
   nn = out_gshape(1)
   mm = out_gshape(2)
 else
   nn = in_mat%sizeb_global(2)
   mm = in_mat%sizeb_global(1)
 end if

 if (present(size_blocs)) then
   size_blocs__ = size_blocs
 else
   ! FIXME: This can cause problems if I start to use round-robin distribution in GWR!!!!!
   size_blocs__(1) = in_mat%sizeb_global(2)
   sb = in_mat%sizeb_global(1) / in_mat%processor%grid%dims(2)
   if (mod(in_mat%sizeb_global(1), in_mat%processor%grid%dims(2)) /= 0) sb = sb + 1
   size_blocs__(2) = sb
   !size_blocs__(2) = in_mat%sizeb_blocs(1); size_blocs__(1) = in_mat%sizeb_blocs(2)
 end if

 call out_mat%init(nn, mm, in_mat%processor, in_mat%istwf_k, size_blocs=size_blocs__)

 ! prototype: call pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

 if (allocated(in_mat%buffer_cplx)) then
#ifdef HAVE_LINALG_SCALAPACK
   select case (trans)
   case ("N")
     ! sub(C) := beta*sub(C) + alpha*sub(A)',
     calpha__ = cone; if (present(alpha)) calpha__ = alpha
     cbeta__ = czero; if (present(beta)) cbeta__ = beta
     call pztranu(nn, mm, calpha__, in_mat%buffer_cplx, ija__(1), ija__(2), &
                  in_mat%descript%tab, cbeta__, out_mat%buffer_cplx, ijc__(1), ijc__(2), out_mat%descript%tab)

   case ("C")
     ! sub(C) := beta * sub(C) + alpha * conjg(sub(A)')
     calpha__ = cone; if (present(alpha)) calpha__ = alpha
     cbeta__ = czero; if (present(beta)) cbeta__ = beta
     call pztranc(nn, mm, calpha__, in_mat%buffer_cplx, ija__(1), ija__(2), &
                  in_mat%descript%tab, cbeta__, out_mat%buffer_cplx, ijc__(1), ijc__(2), out_mat%descript%tab)

   case default
     ABI_ERROR(sjoin("Invalid value for trans:", trans))
   end select

 else if (allocated(in_mat%buffer_real)) then
     ralpha__ = one; if (present(alpha)) ralpha__ = real(alpha)
     rbeta__ = zero; if (present(beta)) rbeta__ = real(beta)
     call pdtran(nn, mm, ralpha__, in_mat%buffer_real, ija__(1), ija__(2), &
                 in_mat%descript%tab, rbeta__, out_mat%buffer_real, ijc__(1), ijc__(2), out_mat%descript%tab)
#endif
 else
   ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
 end if

 if (present(free)) then
   if (free) call in_mat%free()
 end if

end subroutine slk_ptrans
!!***

!!****f* m_slk/slkmat_sp_ptrans
!! NAME
!!  slkmat_sp_ptrans
!!
!! FUNCTION
!! Transposes a matrix
!!
!!    sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
!!
!! where
!!
!!    sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
!!    sub( A ) denotes A(IA:IA+N-1,JA:JA+M-1), and, op( X ) = X'.
!!
!! Thus, op( sub( A ) ) denotes A(IA:IA+N-1,JA:JA+M-1)'.
!! Beta is a scalar, sub( C ) is an m by n submatrix, and sub( A ) is an n by m submatrix.
!!
!! INPUTS
!!  [ija(2)]: (global) The row and column indices in the distributed matrix in_mat indicating
!!       the first row and the first column of the submatrix sub(A), respectively.
!!  [ijc(2)]: (global) The row and column indices in the distributed matrix out_mat
!!   indicating the first row and the first column of the submatrix sub(C), respectively.
!!  [free]: True to deallocate in_mat. Default: False
!!
!! SOURCE

subroutine slkmat_sp_ptrans(in_mat, trans, out_mat, &
                      out_gshape, ija, ijc, size_blocs, alpha, beta, free) ! optional

!Arguments ------------------------------------
 class(slkmat_sp_t),intent(inout) :: in_mat
 character(len=1),intent(in) :: trans
 class(slkmat_sp_t),intent(inout) :: out_mat
 integer,optional,intent(in) :: out_gshape(2), size_blocs(2), ija(2), ijc(2)
 complex(sp),optional,intent(in) :: alpha, beta
 logical,optional,intent(in) :: free

!Local variables-------------------------------
 integer :: sb, mm, nn, size_blocs__(2)
 real(sp) :: ralpha__, rbeta__
 integer :: ija__(2), ijc__(2)
 complex(sp) :: calpha__, cbeta__

! *************************************************************************

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijc__ = [1, 1]; if (present(ijc)) ijc__ = ijc

 ! transposed output (sub)matrix has shape (nn, mm)
 if (present(out_gshape)) then
   nn = out_gshape(1)
   mm = out_gshape(2)
 else
   nn = in_mat%sizeb_global(2)
   mm = in_mat%sizeb_global(1)
 end if

 if (present(size_blocs)) then
   size_blocs__ = size_blocs
 else
   ! FIXME: This can cause problems if I start to use round-robin distribution in GWR!!!!!
   size_blocs__(1) = in_mat%sizeb_global(2)
   sb = in_mat%sizeb_global(1) / in_mat%processor%grid%dims(2)
   if (mod(in_mat%sizeb_global(1), in_mat%processor%grid%dims(2)) /= 0) sb = sb + 1
   size_blocs__(2) = sb
   !size_blocs__(2) = in_mat%sizeb_blocs(1); size_blocs__(1) = in_mat%sizeb_blocs(2)
 end if

 call out_mat%init(nn, mm, in_mat%processor, in_mat%istwf_k, size_blocs=size_blocs__)

 ! prototype: call pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

 if (allocated(in_mat%buffer_cplx)) then
#ifdef HAVE_LINALG_SCALAPACK
   select case (trans)
   case ("N")
     ! sub(C) := beta*sub(C) + alpha*sub(A)',
     calpha__ = cone_sp; if (present(alpha)) calpha__ = alpha
     cbeta__ = czero_sp; if (present(beta)) cbeta__ = beta
     call pctranu(nn, mm, calpha__, in_mat%buffer_cplx, ija__(1), ija__(2), &
                  in_mat%descript%tab, cbeta__, out_mat%buffer_cplx, ijc__(1), ijc__(2), out_mat%descript%tab)

   case ("C")
     ! sub(C) := beta * sub(C) + alpha * conjg(sub(A)')
     calpha__ = cone_sp; if (present(alpha)) calpha__ = alpha
     cbeta__ = czero_sp; if (present(beta)) cbeta__ = beta
     call pctranc(nn, mm, calpha__, in_mat%buffer_cplx, ija__(1), ija__(2), &
                  in_mat%descript%tab, cbeta__, out_mat%buffer_cplx, ijc__(1), ijc__(2), out_mat%descript%tab)

   case default
     ABI_ERROR(sjoin("Invalid value for trans:", trans))
   end select

 else if (allocated(in_mat%buffer_real)) then
     ralpha__ = one_sp; if (present(alpha)) ralpha__ = real(alpha)
     rbeta__ = zero_sp; if (present(beta)) rbeta__ = real(beta)
     call pstran(nn, mm, ralpha__, in_mat%buffer_real, ija__(1), ija__(2), &
                 in_mat%descript%tab, rbeta__, out_mat%buffer_real, ijc__(1), ijc__(2), out_mat%descript%tab)
#endif
 else
   ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
 end if

 if (present(free)) then
   if (free) call in_mat%free()
 end if

end subroutine slkmat_sp_ptrans
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_change_size_blocs
!! NAME
!!  slk_change_size_blocs
!!
!! FUNCTION
!!  Change the block sizes, return new matrix in out_mat
!!
!! INPUTS
!!  [free]: True if `in_mat` should be deallocated. Default: False
!!
!! OUTPUT
!!
!! SOURCE

subroutine slk_change_size_blocs(in_mat, out_mat, &
                                 size_blocs, processor, free)  ! Optional

!Arguments ------------------------------------
 class(basemat_t),target,intent(inout) :: in_mat
 class(basemat_t),intent(out) :: out_mat
 integer,optional,intent(in) :: size_blocs(2)
 class(processor_scalapack), target, optional,intent(in) :: processor
 logical,optional,intent(in) :: free

!Local variables-------------------------------
 type(processor_scalapack), pointer :: processor__

! *************************************************************************

 processor__ => in_mat%processor; if (present(processor)) processor__ => processor

 if (present(size_blocs)) then
   call out_mat%init(in_mat%sizeb_global(1), in_mat%sizeb_global(2), processor__, in_mat%istwf_k, size_blocs=size_blocs)
 else
   call out_mat%init(in_mat%sizeb_global(1), in_mat%sizeb_global(2), processor__, in_mat%istwf_k)
 end if
 !call out_mat%print(header="output matrix generated by slk_change_size_blocs")

 ABI_CHECK(same_type_as(in_mat, out_mat), "in_mat and out_mat should have same type!")

 ! p?gemr2d: Copies a submatrix from one general rectangular matrix to another.
 ! prototype
 !call pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)

#ifdef HAVE_LINALG_SCALAPACK
 select type (in_mat)
 class is (matrix_scalapack)
   select type (out_mat)
   class is (matrix_scalapack)
   if (allocated(in_mat%buffer_cplx)) then
     call pzgemr2d(in_mat%sizeb_global(1), in_mat%sizeb_global(2),  &
                   in_mat%buffer_cplx, 1, 1, in_mat%descript%tab,   &
                   out_mat%buffer_cplx, 1, 1, out_mat%descript%tab, &
                   processor__%grid%ictxt)

   else if (allocated(in_mat%buffer_real)) then
     call pdgemr2d(in_mat%sizeb_global(1), in_mat%sizeb_global(2),  &
                   in_mat%buffer_real, 1, 1, in_mat%descript%tab,   &
                   out_mat%buffer_real, 1, 1, out_mat%descript%tab, &
                   processor__%grid%ictxt)
   else
     ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
   end if
   end select

 class is (slkmat_sp_t)
   select type (out_mat)
   class is (slkmat_sp_t)
   if (allocated(in_mat%buffer_cplx)) then
     call pcgemr2d(in_mat%sizeb_global(1), in_mat%sizeb_global(2),  &
                   in_mat%buffer_cplx, 1, 1, in_mat%descript%tab,   &
                   out_mat%buffer_cplx, 1, 1, out_mat%descript%tab, &
                   processor__%grid%ictxt)

   else if (allocated(in_mat%buffer_real)) then
     call psgemr2d(in_mat%sizeb_global(1), in_mat%sizeb_global(2),  &
                   in_mat%buffer_real, 1, 1, in_mat%descript%tab,   &
                   out_mat%buffer_real, 1, 1, out_mat%descript%tab, &
                   processor__%grid%ictxt)
   else
     ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
   end if
   end select

 class default
   ABI_ERROR("Wrong class")
 end select
#endif

 if (present(free)) then
   if (free) call in_mat%free()
 end if

end subroutine slk_change_size_blocs
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_cut
!! NAME
!!  slk_cut
!!
!! FUNCTION
!!  Extract submatrix of shape (glob_nrows, glob_ncols) starting at `ija` from `in_mat`
!!  and create new matrix with `size_blocs` and `processor`
!!
!! INPUTS
!!  [free]: True if `in_mat` should be deallocated. Default: False
!!
!! OUTPUT
!!
!! SOURCE

subroutine slk_cut(in_mat, glob_nrows, glob_ncols, out_mat, &
                   size_blocs, processor, ija, ijb, free)  ! Optional

!Arguments ------------------------------------
 class(matrix_scalapack),target,intent(inout) :: in_mat
 integer,intent(in) :: glob_nrows, glob_ncols
 class(matrix_scalapack),intent(out) :: out_mat
 integer,optional,intent(in) :: size_blocs(2)
 class(processor_scalapack), target, optional,intent(in) :: processor
 integer,optional,intent(in) :: ija(2), ijb(2)
 logical,optional,intent(in) :: free

!Local variables-------------------------------
 type(processor_scalapack), pointer :: processor__
 integer :: ija__(2), ijb__(2)

! *************************************************************************

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijb__ = [1, 1]; if (present(ijb)) ijb__ = ijb

 processor__ => in_mat%processor; if (present(processor)) processor__ => processor

 if (present(size_blocs)) then
   call out_mat%init(glob_nrows, glob_ncols, processor__, in_mat%istwf_k, size_blocs=size_blocs)
 else
   call out_mat%init(glob_nrows, glob_ncols, processor__, in_mat%istwf_k)
 end if
 !call out_mat%print(header="output matrix generated by slk_cut")

 ! p?gemr2d: Copies a submatrix from one general rectangular matrix to another.
 ! prototype
 !call pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)

 if (allocated(in_mat%buffer_cplx)) then
#ifdef HAVE_LINALG_SCALAPACK
   call pzgemr2d(glob_nrows, glob_ncols, &
                 in_mat%buffer_cplx, ija__(1), ija__(2), in_mat%descript%tab,   &
                 out_mat%buffer_cplx, ijb__(1), ijb__(2), out_mat%descript%tab, &
                 processor__%grid%ictxt)

 else if (allocated(in_mat%buffer_real)) then
   call pdgemr2d(glob_nrows, glob_ncols, &
                 in_mat%buffer_real, ija__(1), ija__(2), in_mat%descript%tab,   &
                 out_mat%buffer_real, ijb__(1), ijb__(2), out_mat%descript%tab, &
                 processor__%grid%ictxt)
#endif
 else
   ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
 end if

 if (present(free)) then
   if (free) call in_mat%free()
 end if

end subroutine slk_cut
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_take_from
!! NAME
!!  slk_take_from
!!
!! FUNCTION
!!  Take values from source
!!  NB: This routine should be called by all procs owning mat and source.
!!
!! INPUTS
!!  [free]: True if source should be deallocated. Default: False
!!
!! OUTPUT
!!
!! SOURCE

subroutine slk_take_from(out_mat, source, &
                         ija, ijb, free) ! optional

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: out_mat
 class(matrix_scalapack),intent(inout) :: source
 integer,optional,intent(in) :: ija(2), ijb(2)
 logical,optional,intent(in) :: free

!Local variables-------------------------------
 integer :: mm, nn
 character(len=500) :: msg
 integer :: ija__(2), ijb__(2)

! *************************************************************************

 ! prototype
 !call pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)

 ! Take care when context A is disjoint from context B. The general rules for which parameters need to be set are:
 !
 !   - All calling processes must have the correct m and n.
 !   - Processes in context A must correctly define all parameters describing A.
 !   - Processes in context B must correctly define all parameters describing B.
 !   - Processes which are not members of context A must pass ctxt_a = -1 and need not set other parameters describing A.
 !   - Processes which are not members of contextB must pass ctxt_b = -1 and need not set other parameters describing B.

 mm = source%sizeb_global(1)
 nn = source%sizeb_global(2)

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijb__ = [1, 1]; if (present(ijb)) ijb__ = ijb

 if (all(out_mat%sizeb_global == -1)) then
   out_mat%descript%tab(CTXT_) = -1
 else
   ABI_CHECK_IEQ(out_mat%istwf_k, source%istwf_k, "istwfk_mat /= istwfk_source")
   if (any(out_mat%sizeb_global /= source%sizeb_global)) then
     msg = sjoin("Matrices should have same global shape but out_mat:", ltoa(out_mat%sizeb_global), &
                 "source:", ltoa(source%sizeb_global))
     ABI_ERROR(msg)
   end if
 end if

 if (allocated(source%buffer_cplx)) then
#ifdef HAVE_LINALG_SCALAPACK
   call pzgemr2d(mm, nn,  &
                 source%buffer_cplx, ija__(1), ija__(2), source%descript%tab,   &
                 out_mat%buffer_cplx, ijb__(1), ijb__(2), out_mat%descript%tab, &
                 source%processor%grid%ictxt)

 else if (allocated(source%buffer_real)) then
   call pdgemr2d(mm, nn,  &
                 source%buffer_real, ija__(1), ija__(2), source%descript%tab,   &
                 out_mat%buffer_real,ijb__(1), ijb__(2), out_mat%descript%tab, &
                 source%processor%grid%ictxt)
#endif
 else
   ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
 end if

 if (present(free)) then
   if (free) call source%free()
 end if

end subroutine slk_take_from
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_sp_take_from
!! NAME
!!  slkmat_sp_take_from
!!
!! FUNCTION
!!  Take values from source
!!  NB: This routine should be called by all procs owning mat and source.
!!
!! INPUTS
!!  [free]: True if source should be deallocated. Default: False
!!
!! OUTPUT
!!
!! SOURCE

subroutine slkmat_sp_take_from(out_mat, source, &
                         ija, ijb, free) ! optional

!Arguments ------------------------------------
 class(slkmat_sp_t),intent(inout) :: out_mat
 class(slkmat_sp_t),intent(inout) :: source
 integer,optional,intent(in) :: ija(2), ijb(2)
 logical,optional,intent(in) :: free

!Local variables-------------------------------
 integer :: mm, nn
 character(len=500) :: msg
 integer :: ija__(2), ijb__(2)

! *************************************************************************

 ! prototype
 !call pzgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)

 ! Take care when context A is disjoint from context B. The general rules for which parameters need to be set are:
 !
 !   - All calling processes must have the correct m and n.
 !   - Processes in context A must correctly define all parameters describing A.
 !   - Processes in context B must correctly define all parameters describing B.
 !   - Processes which are not members of context A must pass ctxt_a = -1 and need not set other parameters describing A.
 !   - Processes which are not members of contextB must pass ctxt_b = -1 and need not set other parameters describing B.

 mm = source%sizeb_global(1)
 nn = source%sizeb_global(2)

 ija__ = [1, 1]; if (present(ija)) ija__ = ija
 ijb__ = [1, 1]; if (present(ijb)) ijb__ = ijb

 if (all(out_mat%sizeb_global == -1)) then
   out_mat%descript%tab(CTXT_) = -1
 else
   ABI_CHECK_IEQ(out_mat%istwf_k, source%istwf_k, "istwfk_mat /= istwfk_source")
   if (any(out_mat%sizeb_global /= source%sizeb_global)) then
     msg = sjoin("Matrices should have same global shape but out_mat:", ltoa(out_mat%sizeb_global), &
                 "source:", ltoa(source%sizeb_global))
     ABI_ERROR(msg)
   end if
 end if

 if (allocated(source%buffer_cplx)) then
#ifdef HAVE_LINALG_SCALAPACK
   call pcgemr2d(mm, nn,  &
                 source%buffer_cplx, ija__(1), ija__(2), source%descript%tab,   &
                 out_mat%buffer_cplx, ijb__(1), ijb__(2), out_mat%descript%tab, &
                 source%processor%grid%ictxt)

 else if (allocated(source%buffer_real)) then
   call psgemr2d(mm, nn,  &
                 source%buffer_real, ija__(1), ija__(2), source%descript%tab,   &
                 out_mat%buffer_real,ijb__(1), ijb__(2), out_mat%descript%tab, &
                 source%processor%grid%ictxt)
#endif
 else
   ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
 end if

 if (present(free)) then
   if (free) call source%free()
 end if

end subroutine slkmat_sp_take_from
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_collect_cplx
!! NAME
!!  slk_collect_cplx
!!
!! FUNCTION
!!  Return on all processors the complex submatrix of shape (mm, nn) starting at position ija.
!!  NB: `out_carr` is allocated by the routine.
!!  If optional argument request is use, the routine uses non-blocking BCAST and client code is
!!  supposed to wait before accessing out_carr.
!!
!! SOURCE

subroutine slk_collect_cplx(in_mat, mm, nn, ija, out_carr, request)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(in) :: in_mat
 integer,intent(in) :: mm, nn, ija(2)
 complex(dp) ABI_ASYNC, allocatable,intent(out) :: out_carr(:,:)
 integer ABI_ASYNC, optional,intent(out) :: request

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: ierr
 type(processor_scalapack) :: self_processor
 type(matrix_scalapack) :: out_mat

! *************************************************************************

 ABI_CHECK(allocated(in_mat%buffer_cplx), "buffer_cplx is not allocated")

 if (in_mat%processor%grid%nbprocs == 1) then
   ! Copy buffer and return
   ABI_MALLOC(out_carr, (mm, nn))
   out_carr(:,:) = in_mat%buffer_cplx(ija(1):ija(1)+mm-1, ija(2):ija(2)+nn-1); return
 end if

 ! Two-step algorithm:
 !     1) Use pzgemr2d to collect submatrix on master.
 !     2) Master brodacasts submatrix.

 if (in_mat%processor%myproc == master) then
   call self_processor%init(xmpi_comm_self)
   call out_mat%init(mm, nn, self_processor, in_mat%istwf_k, size_blocs=[mm, nn])
 else
   out_mat%descript%tab(CTXT_) = -1
 end if

#ifdef HAVE_LINALG_SCALAPACK
 call pzgemr2d(mm, nn,  &
               in_mat%buffer_cplx, ija(1), ija(2), in_mat%descript%tab,   &
               out_mat%buffer_cplx, 1, 1, out_mat%descript%tab, &
               in_mat%processor%grid%ictxt)
#endif

 if (in_mat%processor%myproc == master) then
   ABI_MOVE_ALLOC(out_mat%buffer_cplx, out_carr)
   call out_mat%free()
   call self_processor%free()
 else
   ABI_MALLOC(out_carr, (mm, nn))
 end if

 if (present(request)) then
   call xmpi_ibcast(out_carr, master, in_mat%processor%comm, request, ierr)
 else
   call xmpi_bcast(out_carr, master, in_mat%processor%comm, ierr)
 end if

end subroutine slk_collect_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slkmat_sp_collect_cplx
!! NAME
!!  slkmat_sp_collect_cplx
!!
!! FUNCTION
!!  Return on all processors the complex submatrix of shape (mm, nn) starting at position ija.
!!  NB: `out_carr` is allocated by the routine.
!!  If optional argument request is use, the routine uses non-blocking BCAST and client code is
!!  supposed to wait before accessing out_carr.
!!
!! SOURCE

subroutine slkmat_sp_collect_cplx(in_mat, mm, nn, ija, out_carr, request)

!Arguments ------------------------------------
 class(slkmat_sp_t),intent(in) :: in_mat
 integer,intent(in) :: mm, nn, ija(2)
 complex(sp) ABI_ASYNC, allocatable,intent(out) :: out_carr(:,:)
 integer ABI_ASYNC, optional,intent(out) :: request

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: ierr
 type(processor_scalapack) :: self_processor
 type(slkmat_sp_t) :: out_mat

! *************************************************************************

 ABI_CHECK(allocated(in_mat%buffer_cplx), "buffer_cplx is not allocated")

 if (in_mat%processor%grid%nbprocs == 1) then
   ! Copy buffer and return
   ABI_MALLOC(out_carr, (mm, nn))
   out_carr(:,:) = in_mat%buffer_cplx(ija(1):ija(1)+mm-1, ija(2):ija(2)+nn-1); return
 end if

 ! Two-step algorithm:
 !     1) Use pzgemr2d to collect submatrix on master.
 !     2) Master brodacasts submatrix.

 if (in_mat%processor%myproc == master) then
   call self_processor%init(xmpi_comm_self)
   call out_mat%init(mm, nn, self_processor, in_mat%istwf_k, size_blocs=[mm, nn])
 else
   out_mat%descript%tab(CTXT_) = -1
 end if

#ifdef HAVE_LINALG_SCALAPACK
 call pcgemr2d(mm, nn,  &
               in_mat%buffer_cplx, ija(1), ija(2), in_mat%descript%tab,   &
               out_mat%buffer_cplx, 1, 1, out_mat%descript%tab, &
               in_mat%processor%grid%ictxt)
#endif

 if (in_mat%processor%myproc == master) then
   ABI_MOVE_ALLOC(out_mat%buffer_cplx, out_carr)
   call out_mat%free()
   call self_processor%free()
 else
   ABI_MALLOC(out_carr, (mm, nn))
 end if

 if (present(request)) then
   call xmpi_ibcast(out_carr, master, in_mat%processor%comm, request, ierr)
 else
   call xmpi_bcast(out_carr, master, in_mat%processor%comm, ierr)
 end if

end subroutine slkmat_sp_collect_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_get_trace
!! NAME
!!  slk_get_trace
!!
!! FUNCTION
!!  Compute the trace of an N-by-N distributed matrix.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

complex(dp) function slk_get_trace(mat) result(ctrace)

!Arguments ------------------------------------
 class(basemat_t), intent(in) :: mat

!Local variables-------------------------------
#ifdef HAVE_LINALG_SCALAPACK
 real(dp) :: rtrace
 real(sp) :: rtrace_sp
 complex(sp) :: ctrace_sp
 real(dp),external :: PDLATRA
 real(sp),external :: PSLATRA
 complex(sp),external :: PCLATRA
 complex(dp),external :: PZLATRA
#endif

! *************************************************************************

 ABI_CHECK_IEQ(mat%sizeb_global(1), mat%sizeb_global(2), "get_trace assumes square matrix!")

 ! prototype for complex version.
 ! COMPLEX*16 FUNCTION PZLATRA( N, A, IA, JA, DESCA )
#ifdef HAVE_LINALG_SCALAPACK
 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_cplx)) then
     ctrace = PZLATRA(mat%sizeb_global(1), mat%buffer_cplx, 1, 1, mat%descript%tab)
   else if (allocated(mat%buffer_real)) then
     rtrace = PDLATRA(mat%sizeb_global(1), mat%buffer_real, 1, 1, mat%descript%tab)
     ctrace = rtrace
   else
     ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
   end if
 class is (slkmat_sp_t)
    if (allocated(mat%buffer_cplx)) then
      ctrace_sp = PCLATRA(mat%sizeb_global(1), mat%buffer_cplx, 1, 1, mat%descript%tab)
      ctrace = ctrace_sp
    else if (allocated(mat%buffer_real)) then
      rtrace_sp = PSLATRA(mat%sizeb_global(1), mat%buffer_real, 1, 1, mat%descript%tab)
      ctrace = rtrace_sp
    else
      ABI_ERROR("Neither buffer_cplx nor buffer_real are allocated!")
    end if
 class default
   ABI_ERROR("Wrong class")
 end select
#endif

end function slk_get_trace
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_set_imag_diago_to_zero
!! NAME
!!  slk_set_imag_diago_to_zero
!!
!! FUNCTION
!!  Set the imaginary part of the diagonal to zero.
!!  Return in local_max the max of the imaginar part in the local buffer.
!!  No MPI communication is performed inside the routine. Client code can easily reduce
!!  local_max within the PBLAS communicator if needed.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine slk_set_imag_diago_to_zero(mat, local_max)

!Arguments ------------------------------------
 class(basemat_t), intent(inout) :: mat
 real(dp),intent(out) :: local_max

!Local variables-------------------------------
 integer :: il1, iglob1, il2, iglob2

! *************************************************************************

 local_max = -huge(one)

 select type (mat)
 class is (matrix_scalapack)
   if (allocated(mat%buffer_real)) return
   do il2=1,mat%sizeb_local(2)
     iglob2 = mat%loc2gcol(il2)
     do il1=1,mat%sizeb_local(1)
       iglob1 = mat%loc2grow(il1)
       if (iglob1 == iglob2) then
         local_max = max(local_max, aimag(mat%buffer_cplx(il1, il2)))
         mat%buffer_cplx(il1, il2) = real(mat%buffer_cplx(il1, il2))
       end if
     end do
   end do

 class is (slkmat_sp_t)
   if (allocated(mat%buffer_real)) return
   do il2=1,mat%sizeb_local(2)
     iglob2 = mat%loc2gcol(il2)
     do il1=1,mat%sizeb_local(1)
       iglob1 = mat%loc2grow(il1)
       if (iglob1 == iglob2) then
         local_max = max(local_max, aimag(mat%buffer_cplx(il1, il2)))
         mat%buffer_cplx(il1, il2) = real(mat%buffer_cplx(il1, il2))
       end if
     end do
   end do
 class default
   ABI_ERROR("Wrong class")
 end select


end subroutine slk_set_imag_diago_to_zero
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_write
!! NAME
!!  slk_write
!!
!! FUNCTION
!!  Routine to write a square scaLAPACK-distributed matrix to an external file using MPI-IO.
!!
!! INPUTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK distribution with the local buffer
!!    containing the distributed matrix.
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix is used:
!!    = "U":  Upper triangular
!!    = "L":  Lower triangular
!!    = "A":  Full matrix (used for general complex matrices)
!!  is_fortran_file=.FALSE. is C stream is used. .TRUE. for writing Fortran binary files.
!!  [fname]= Mutually exclusive with mpi_fh. The name of the external file on which the matrix will be written.
!!           The file is open and closed inside the routine with MPI flags specified by flags.
!!  [mpi_fh]=File handler associated to the file (already open in the caller). Not compatible with fname.
!!  [flags]=MPI-IO flags used to open the file in MPI_FILE_OPEN.
!!    Default is MPI_MODE_CREATE + MPI_MODE_WRONLY + MPI_MODE_EXCL.
!!  [glob_subarray(2,2)] = Used to select the subarray of the global matrix. Used only when uplo="All"
!!     NOTE that each node should call the routine with the same value.
!!     glob_subarray(:,1)=starting global coordinates of the subarray in each dimension
!!       (array of nonnegative integers >=1, <=array_of_sizes)
!!     glob_subarray(:,2)=Number of elements in each dimension of the subarray (array of positive integers)
!!
!! OUTPUT
!!  Only writing. The global scaLAPACK matrix is written to file fname.
!!  If fname is present then the file is open and closed inside the routine. Any exception is fatal.
!!
!! SIDE EFFECTS
!!  [offset]=
!!    input:  Offset used to access the content of the file. Default is zero.
!!    output: New offset incremented with the byte size of the matrix that has been read (Fortran
!!            markers are included if is_fortran_file=.TRUE.)
!! TODO
!!  - Generalize the implementation adding the writing the real buffer.
!!  - This routine should be removed and replaced by hdf5 + mpi-io
!!
!! SOURCE

subroutine slk_write(Slk_mat, uplo, is_fortran_file, fname,mpi_fh, offset, flags, glob_subarray)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: flags
 integer,optional,intent(inout) :: mpi_fh
 integer(XMPI_OFFSET_KIND),optional,intent(inout) :: offset
 logical,intent(in) :: is_fortran_file
 character(len=*),optional,intent(in) :: fname
 character(len=*),intent(in) :: uplo
 class(matrix_scalapack),intent(in) :: Slk_mat
!array
 integer,optional,intent(in) :: glob_subarray(2,2)

!Local variables ------------------------------
!scalars
#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
 integer :: jloc,iloc,nrows_glob,ncols_glob,elw,nrows_w,ncols_w ! iglob,jglob,
 integer :: slk_type,offset_err,etype,nfrec,bsize_elm,mpi_type_elm
 integer(XMPI_OFFSET_KIND) :: my_offset
 logical :: do_open
 integer :: comm,my_flags,my_fh,buffer_size
 integer :: ierr,nelw,col_glob ! ij_loc,
!arrays
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecord(:)
 integer,pointer :: elw2slk(:,:)
 complex(dpc),allocatable :: buffer1_cplx(:)
 character(len=500) :: msg

!************************************************************************

 ABI_CHECK(allocated(Slk_mat%buffer_cplx), "buffer_cplx not allocated")

 if (firstchar(uplo, ["U","L"]) .and. Slk_mat%sizeb_global(1) /= Slk_mat%sizeb_global(2) ) then
   ABI_ERROR("rectangular matrices are not compatible with the specified uplo")
 end if

 if (PRESENT(glob_subarray).and. .not. firstchar(uplo, ["A"])) then
   ABI_ERROR("glob_subarray should not be used when uplo/=All")
 end if

 do_open = PRESENT(fname)
 if (do_open) then
   ABI_CHECK(.not.PRESENT(fname),"fname should not be present")
 else
   ABI_CHECK(PRESENT(mpi_fh),"mpi_fh should be present")
 end if

 my_offset=0; if (PRESENT(offset)) my_offset=offset

 comm = Slk_mat%processor%comm

 nrows_glob=Slk_mat%sizeb_global(1)
 ncols_glob=Slk_mat%sizeb_global(1)
 buffer_size= PRODUCT(Slk_mat%sizeb_local(1:2))

 call slk_mat%bsize_and_type(bsize_elm, mpi_type_elm)

 if (do_open) then !Open the file.
   my_flags=MPI_MODE_CREATE + MPI_MODE_WRONLY + MPI_MODE_APPEND
   if (PRESENT(flags)) my_flags = flags

   call MPI_FILE_OPEN(comm, fname, my_flags, MPI_INFO_NULL, my_fh, ierr)
   ABI_CHECK_MPI(ierr, "MPI_FILE_OPEN "//TRIM(fname))
 else
   my_fh = mpi_fh
 end if

 if (PRESENT(glob_subarray)) then
   call slk_single_fview_write(Slk_mat,uplo,nelw,elw2slk,etype,slk_type,offset_err,&
     is_fortran_file=is_fortran_file,glob_subarray=glob_subarray)
 else
   call slk_single_fview_write(Slk_mat,uplo,nelw,elw2slk,etype,slk_type,offset_err,&
     is_fortran_file=is_fortran_file)
 end if

 if (offset_err /= 0) then
   write(msg,"(3a)")&
    " Global position index cannot be stored in standard Fortran integer ",ch10,&
    " scaLAPACK matrix cannot be read with a single MPI-IO call."
   ABI_ERROR(msg)
 end if

 call MPI_FILE_SET_VIEW(my_fh, my_offset, etype, slk_type, 'native', MPI_INFO_NULL, ierr)
 ABI_CHECK_MPI(ierr,"SET_VIEW")

 call MPI_TYPE_FREE(slk_type,ierr)
 ABI_CHECK_MPI(ierr,"MPI_type_FREE")

 if (nelw == buffer_size) then
   ! Dump Slk_mat% immediately.
   call MPI_FILE_WRITE_ALL(my_fh, Slk_mat%buffer_cplx, buffer_size, MPI_DOUBLE_complex, MPI_STATUS_IGNORE, ierr)
   ABI_CHECK_MPI(ierr,"WRITE_ALL")
 else
   ! Have to extract the data to be written.
   ABI_MALLOC(buffer1_cplx,(nelw))
   do elw=1,nelw
     iloc = elw2slk(1,elw)
     jloc = elw2slk(2,elw)
     buffer1_cplx(elw) = Slk_mat%buffer_cplx(iloc,jloc)
   end do
   call MPI_FILE_WRITE_ALL(my_fh, buffer1_cplx, nelw, MPI_DOUBLE_complex, MPI_STATUS_IGNORE, ierr)
   ABI_CHECK_MPI(ierr,"WRITE_ALL")
   ABI_FREE(buffer1_cplx)
 end if

 ABI_FREE(elw2slk)
 !
 ! Number of columns and rows that have been written.
 ! Used to write the Fortran markers and to increment the offset.
 nrows_w = nrows_glob
 ncols_w = ncols_glob
 if (PRESENT(glob_subarray)) then
   nrows_w = glob_subarray(1,2) - glob_subarray(1,1) + 1
   ncols_w = glob_subarray(2,2) - glob_subarray(2,1) + 1
   if (.not.firstchar(uplo, ["A"])) then
     ABI_ERROR("glob_subarray should not be used when uplo/=All")
   end if
 end if

 !TODO check whether slk_single_fview_write can report an offset to reduce the extent.
 if (is_fortran_file) then ! Collective writing of the Fortran markers.
   nfrec = ncols_w
   ABI_MALLOC(bsize_frecord,(nfrec))
   if (firstchar(uplo, ["A"])) then
     bsize_frecord = nrows_w * bsize_elm
   else if (firstchar(uplo, ["U"])) then
     bsize_frecord = (/(col_glob * bsize_elm, col_glob=1,nfrec)/)
   else if (firstchar(uplo, ["L"])) then
     bsize_frecord = (/(col_glob * bsize_elm, col_glob=nfrec,1,-1)/)
   else
     ABI_ERROR("Wrong uplo")
   end if
   call xmpio_write_frmarkers(mpi_fh,my_offset,xmpio_collective,nfrec,bsize_frecord,ierr)
   ABI_CHECK(ierr==0,"Error while writing Fortran markers")
   ABI_FREE(bsize_frecord)
 end if

 if (do_open) then
   ! Close the file.
   call MPI_FILE_CLOSE(my_fh, ierr)
   ABI_CHECK_MPI(ierr,"FILE_CLOSE")
 end if

 ! Increment the offset
 if (PRESENT(offset)) then
   if (firstchar(uplo, ["A"])) then
     offset = offset + nrows_w*ncols_w*bsize_elm
     if (is_fortran_file) offset = offset + ncols_w*2*xmpio_bsize_frm
   else if (firstchar(uplo, ["U","L"])) then
     offset = offset + ( (Slk_mat%sizeb_global(2) * (Slk_mat%sizeb_global(2))+1)/2 ) * bsize_elm
     if (is_fortran_file) offset = offset + Slk_mat%sizeb_global(2)*2*xmpio_bsize_frm
   else
     ABI_ERROR("Wrong uplo")
   end if
 end if

 call xmpi_barrier(comm)
 RETURN

#else
  ABI_ERROR("MPI-IO support not activated")
#endif

end subroutine slk_write
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_read
!! NAME
!!  slk_read
!!
!! FUNCTION
!!  Routine to read a square scaLAPACK distributed matrix from an external file using MPI-IO.
!!
!! INPUTS
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix is stored on disk:
!!    = "U":  Upper triangular is stored
!!    = "L":  Lower triangular is stored
!!    = "A":  Full matrix (used for general complex matrices)
!!  symtype=Symmetry type of the matrix stored on disk (used only if uplo = "L" or "A").
!!    = "H" for Hermitian matrix
!!    = "S" for symmetric matrix.
!!    = "N" if matrix has no symmetry (not compatible with uplo="L" or uplo="U".
!!  is_fortran_file=.FALSE. is C stream is used. .TRUE. for writing Fortran binary files.
!!  [fname]= Mutually exclusive with mpi_fh. The name of the external file from which the matrix will be read.
!!           The file is open and closed inside the routine with MPI flags specified by flags.
!!  [mpi_fh]=File handler associated to the file (already open in the caller). Not compatible with fname.
!!  [flags]=MPI-IO flags used to open the file in MPI_FILE_OPEN. Default is MPI_MODE_RDONLY. Referenced only when fname is used.
!!
!! SIDE EFFECTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK distribution with the local buffer
!!    supposed to be allocated.
!!    %buffer_cplx=Local buffer containg the distributed matrix stored on the external file.
!!  If fname is present then the file is opened and closed inside the routine. Any exception is fatal.
!!  [offset]=
!!    input:  Offset used to access the content of the file. Default is zero.
!!    output: New offset incremented with the byte size of the matrix that has been read (Fortran
!!            markers are included if is_fortran_file=.TRUE.)
!!
!! TODO
!!  - Generalize the implementation adding the reading of the real buffer.
!!
!!  - This routine is not portable as this kind of access pattern is not supported by all MPI implementations
!!    E.g. with MPICH we have
!!
!! --- !ERROR
!! src_file: m_slk.F90
!! src_line: 3780
!! mpi_rank: 1
!! message: |
!!     SET_VIEW
!!     Other I/O error , error stack:
!!     ADIO_Set_view(48):  **iobadoverlap displacements of filetype must be in a monotonically nondecreasing order
!! ...
!!
!! - This routine should be removed and replaced by hdf5 + mpi-io
!!
!! SOURCE

subroutine slk_read(Slk_mat,uplo,symtype,is_fortran_file,fname,mpi_fh,offset,flags)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: flags,mpi_fh
 integer(XMPI_OFFSET_KIND),optional,intent(inout) :: offset
 character(len=*),optional,intent(in) :: fname
 character(len=*),intent(in) :: uplo,symtype
 logical,intent(in) :: is_fortran_file
 class(matrix_scalapack),intent(inout) :: Slk_mat

!Local variables ------------------------------
#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
!scalars
 integer :: nrows_glob,offset_err,slk_type,etype
 integer(XMPI_OFFSET_KIND) :: my_offset
 logical :: do_open
 integer :: comm,my_flags,my_fh,buffer_size,ierr,col_glob
 integer :: nfrec,bsize_elm,mpi_type_elm
 !complex(dpc) :: ctest
 logical,parameter :: check_frm=.TRUE.
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecord(:)
!arrays
 character(len=500) :: msg

!************************************************************************

 do_open = PRESENT(fname)
 if (do_open) then
   ABI_CHECK(.not.PRESENT(fname), "fname should not be present")
 else
   ABI_CHECK(PRESENT(mpi_fh), "mpi_fh should be present")
 end if

 my_offset=0; if (PRESENT(offset)) my_offset=offset

 ABI_CHECK(allocated(Slk_mat%buffer_cplx), "%buffer_cplx not allocated")
 if (firstchar(uplo, ["U","L"]) .and. Slk_mat%sizeb_global(1) /= Slk_mat%sizeb_global(2) ) then
   ABI_ERROR("rectangular matrices are not compatible with the specified uplo")
 end if

 nrows_glob = Slk_mat%sizeb_global(1)

 buffer_size= PRODUCT(Slk_mat%sizeb_local(1:2))

 call wrtout(std_out, "slk_read: Using MPI-IO")

 comm = Slk_mat%processor%comm

 if (do_open) then ! Open the file.
   my_flags=MPI_MODE_RDONLY; if (PRESENT(flags)) my_flags = flags
   call MPI_FILE_OPEN(comm, fname, my_flags, MPI_INFO_NULL, my_fh, ierr)
   ABI_CHECK_MPI(ierr,"FILE_OPEN "//TRIM(fname))
 else
   my_fh = mpi_fh
 end if

 call slk_single_fview_read(Slk_mat,uplo,etype,slk_type,offset_err,is_fortran_file=is_fortran_file)

 if (offset_err/=0) then
   write(msg,"(3a)")&
    "Global position index cannot be stored in standard Fortran integer ",ch10,&
    "scaLAPACK matrix cannot be read with a single MPI-IO call."
   ABI_ERROR(msg)
 end if

 call MPI_FILE_SET_VIEW(my_fh, my_offset, etype, slk_type, 'native', MPI_INFO_NULL, ierr)
 ABI_CHECK_MPI(ierr,"SET_VIEW")

 call MPI_FILE_READ_ALL(my_fh, Slk_mat%buffer_cplx, buffer_size, MPI_DOUBLE_complex, MPI_STATUS_IGNORE, ierr)
 ABI_CHECK_MPI(ierr,"READ_ALL")

 ! Symmetrize local buffer if uplo /= "All"
 call slk_symmetrize(Slk_mat, uplo, symtype)

!BEGINDEBUG
!call MPI_FILE_READ_AT(mpi_fh,my_offset+xmpio_bsize_frm,ctest,1,MPI_DOUBLE_complex,MPI_STATUS_IGNORE,ierr)
!write(std_out,*)"ctest",ctest
!call MPI_FILE_READ_AT(mpi_fh,my_offset+2*xmpio_bsize_frm,ctest,1,MPI_DOUBLE_complex,MPI_STATUS_IGNORE,ierr)
!write(std_out,*)"ctest",ctest
!ENDDEBUG

 !call print_arr(Slk_mat%buffer_cplx,max_r=10,max_c=10,unit=std_out)
 !
 ! Close the file and release the MPI filetype.
 call MPI_type_FREE(slk_type,ierr)
 ABI_CHECK_MPI(ierr,"MPI_type_FREE")

 call slk_mat%bsize_and_type(bsize_elm, mpi_type_elm)

!It seems that personal call makes the code stuck
!if (is_fortran_file .and. check_frm .and. Slk_mat%Processor%myproc==0) then ! Master checks the Fortran markers.
 if (is_fortran_file .and. check_frm) then ! Master checks the Fortran markers.
   call wrtout(std_out,"Checking Fortran record markers...", do_flush=.True.)
   nfrec = Slk_mat%sizeb_global(2)
   ABI_MALLOC(bsize_frecord,(nfrec))
   if (firstchar(uplo, ["A"])) then
     bsize_frecord = Slk_mat%sizeb_global(1) * bsize_elm
   else if (firstchar(uplo, ["U"])) then
     bsize_frecord = (/(col_glob * bsize_elm, col_glob=1,nfrec)/)
   else if (firstchar(uplo, ["L"])) then
     bsize_frecord = (/(col_glob * bsize_elm, col_glob=nfrec,1,-1)/)
   else
     ABI_ERROR("Wrong uplo")
   end if
   call xmpio_check_frmarkers(my_fh,my_offset,xmpio_collective,nfrec,bsize_frecord,ierr)
   ABI_CHECK(ierr==0,"Wrong Fortran record markers")
   ABI_FREE(bsize_frecord)
 end if

 if (do_open) then ! Close the file.
   call MPI_FILE_CLOSE(my_fh, ierr)
   ABI_CHECK_MPI(ierr,"FILE_CLOSE")
 end if

!Increment the offset
 if (PRESENT(offset)) then
   if (firstchar(uplo, ["A"])) then
     offset = offset + PRODUCT(Slk_mat%sizeb_global(1:2)) * bsize_elm
     if (is_fortran_file) offset = offset + Slk_mat%sizeb_global(2)*2*xmpio_bsize_frm
   else if (firstchar(uplo, ["U","L"])) then
     offset = offset + ( (Slk_mat%sizeb_global(2) * (Slk_mat%sizeb_global(2))+1)/2 ) * bsize_elm
     if (is_fortran_file) offset = offset + Slk_mat%sizeb_global(2)*2*xmpio_bsize_frm
   else
     ABI_ERROR("Wrong uplo")
   end if
 end if

 call xmpi_barrier(comm)
 RETURN

#else
 ABI_ERROR("MPI-IO support not enabled")
#endif

end subroutine slk_read
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_single_fview_read_mask
!! NAME
!!  slk_single_fview_read_mask
!!
!! FUNCTION
!!  Return an MPI datatype that can be used to read a scaLAPACK distributed matrix from
!!  a binary file using MPI-IO. The view is created using the user-defined mask function
!!  mask_of_glob. The storage of the data on file is described via the user-defined function offset_of_glob.
!!
!! INPUTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK matrix.
!!  mask_of_glob(row_glob,col_glob,size_glob) is an integer function that accepts in input
!!     the global indices of the matrix size_glob(1:2) are the global dimensions.
!!     Return 0 if (row_glob,col_glob) should not be read.
!!  offset_of_glob(row_glob,col_glob,size_glob,nsblocks,sub_block,bsize_elm,bsize_frm)
!!  nsblocks=Number of sub-blocks (will be passed to offset_of_glob)
!!  sub_block(2,2,nsblocks)=Global coordinates of the extremal points delimiting the sub-blocs
!!   e.g. sub_block(:,1,1) gives the coordinates of the left upper corner of the first block.
!!        sub_block(:,2,1) gives the coordinates of the right lower corner of the first block.
!!  [is_fortran_file]=.FALSE. is C stream is used. Set to .TRUE. for writing Fortran binary files
!!    with record marker.
!!
!! OUTPUT
!!  my_nel=Number of elements that will be read by this node.
!!  etype=Elementary data type (handle) defining the elementary unit used to access the file.
!!    This is the elementary type that must be used to creae the view (MPI_BYTE is used).
!!  slk_type=New MPI type that can be used to instantiate the MPI-IO view for the Fortran file.
!!     Note that the view assumes that the file pointer points to the FIRST Fortran record marker.
!!  offset_err=Error code. A returned non-zero value signals that the global matrix is too large
!!    for a single MPI-IO access. See notes in other slk_single_fview_* routines.
!!
!! SIDE EFFECTS
!!  myel2loc(:,:)
!!    input: pointer to NULL
!!    output: myel2loc(2,my_nel):  myel2loc(:,el) gives (iloc,jloc) for el=1,my_nel.
!!
!! SOURCE

subroutine slk_single_fview_read_mask(Slk_mat,mask_of_glob,offset_of_glob,nsblocks,sub_block,&
                                      my_nel,myel2loc,etype,slk_type,offset_err,is_fortran_file)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsblocks
 integer,intent(out) :: my_nel,offset_err,slk_type,etype
 logical,optional,intent(in) :: is_fortran_file
 class(matrix_scalapack),intent(in) :: Slk_mat
!arrays
 integer,intent(in) :: sub_block(2,2,nsblocks)
 integer,pointer :: myel2loc(:,:)

 interface
   function mask_of_glob(row_glob,col_glob,size_glob)
     use defs_basis
     integer :: mask_of_glob
     integer,intent(in) :: row_glob,col_glob
     integer,intent(in) :: size_glob(2)
   end function mask_of_glob
 end interface

 interface
   function offset_of_glob(row_glob,col_glob,size_glob,nsblocks,sub_block,bsize_elm,bsize_frm)
     use defs_basis
     use m_xmpi
     integer(XMPI_OFFSET_KIND) :: offset_of_glob
     integer,intent(in) :: row_glob,col_glob,bsize_elm,bsize_frm,nsblocks
     integer,intent(in) :: size_glob(2),sub_block(2,2,nsblocks)
   end function offset_of_glob
 end interface

!Local variables ------------------------------
!scalars
 integer :: el,jloc,iloc,iglob,jglob,mpi_err,sweep
 integer :: bsize_frm,mpi_type_elm,bsize_elm
 integer(XMPI_OFFSET_KIND) :: tmp_off,max_displ
!arrays
 character(len=500) :: msg
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

#ifdef HAVE_MPI_IO
 bsize_frm = xmpio_bsize_frm  ! Byte size of the Fortran record marker.
 if (PRESENT(is_fortran_file)) then
   if (.not.is_fortran_file) bsize_frm = 0
 end if

 ! Byte size of the matrix element.
 call slk_mat%bsize_and_type(bsize_elm, mpi_type_elm)

 ! Find the number of local matrix elements to be read, then create the table myel2loc.
 do sweep=1,2
   if (sweep==2) then
      ABI_MALLOC(myel2loc,(2,my_nel))
   end if
   my_nel=0

   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       if ( mask_of_glob(iglob,jglob,Slk_mat%sizeb_global)/= 0) then ! Will fill this entry.
         my_nel  = my_nel+1
         if (sweep==2) myel2loc(:,my_nel) = (/iloc,jloc/)
       end if
     end do
   end do
 end do

 etype = MPI_BYTE

 ! Define the mapping between scaLAPACK buffer and the storage on file.
 ! Note that the view assumes that the file pointer points to the first Fortran record marker.
 ABI_MALLOC(block_length,(my_nel+2))
 ABI_MALLOC(block_displ,(my_nel+2))
 ABI_MALLOC(block_type,(my_nel+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 offset_err=0; max_displ=0
 do el=1,my_nel
   iloc = myel2loc(1,el)
   jloc = myel2loc(2,el)
   call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
   tmp_off = offset_of_glob(iglob,jglob,Slk_mat%sizeb_global,nsblocks,sub_block,bsize_elm,bsize_frm)
   if (xmpio_max_address(tmp_off)) offset_err=1   ! Test for possible wraparounds.
   max_displ = MAX(max_displ,tmp_off)
   block_displ (el+1) = tmp_off
   block_type  (el+1) = mpi_type_elm
   block_length(el+1) = 1
   !write(std_out,*)" iglob, jglob, tmp_off ",iglob, jglob, tmp_off
 end do
 !write(std_out,*)" MAX displ is ",MAXVAL(block_displ)

 if (offset_err/=0) then  ! just warn, let the caller handle the exception.
   write(msg,"(3a)")&
    " Global position index cannot be stored in standard Fortran integer ",ch10,&
    " scaLAPACK matrix cannot be read with a single MPI-IO call ."
   ABI_WARNING(msg)
 end if

 block_length(my_nel+2) = 1
 block_displ (my_nel+2) = max_displ
 block_type  (my_nel+2) = MPI_UB

 call xmpio_type_struct(my_nel+2,block_length,block_displ,block_type,slk_type,mpi_err)
 ABI_CHECK_MPI(mpi_err,"MPI_type_STRUCT")

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

 call MPI_type_COMMIT(slk_type,mpi_err)
 ABI_CHECK_MPI(mpi_err,"MPI_type_COMMIT")

#else
 ABI_ERROR("MPI-IO is mandatatory in slk_single_fview_read_mask")
#endif

end subroutine slk_single_fview_read_mask
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_symmetrize
!! NAME
!!  slk_symmetrize
!!
!! FUNCTION
!!  Symmetrize a square scaLAPACK matrix.
!!
!! INPUTS
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix has been read
!!    = "U":  Upper triangular has been read.
!!    = "L":  Lower triangular has been read.
!!    = "A":  Full matrix (used for general complex matrices)
!!  symtype=Symmetry type of the matrix (used only if uplo = "L" or "A").
!!    = "H" for Hermitian matrix
!!    = "S" for symmetric matrix.
!!    = "N" if matrix has no symmetry (not compatible with uplo="L" or uplo="U".
!!
!! SIDE EFFECTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK distribution with the local buffer
!!    supposed to be allocated.
!!    %buffer_cplx=Local buffer containg the distributed matrix stored on the external file.
!!
!! SOURCE

subroutine slk_symmetrize(Slk_mat, uplo, symtype)

!Arguments ------------------------------------
!scalars
 class(matrix_scalapack),intent(inout) :: Slk_mat
 character(len=*),intent(in) :: symtype
 character(len=*),intent(in) :: uplo

!Local variables ------------------------------
!scalars
 integer :: jloc,iloc,iglob,jglob,ij_loc
 logical :: is_hermitian,is_real,is_cplx,is_symmetric
 character(len=500) :: msg

!************************************************************************

 is_cplx = (allocated(Slk_mat%buffer_cplx))
 is_real = (allocated(Slk_mat%buffer_real))

 ! One and only one buffer should be allocated.
 if (is_real .and. is_cplx) then
   write(msg,'(a,2l1)')" ScaLAPACK buffers are not allocated correctly, is_real=, is_cplx ",is_real,is_cplx
   ABI_ERROR(msg)
 end if

 if (is_real) RETURN

 is_hermitian=.FALSE.; is_symmetric=.FALSE.
 select case (symtype(1:1))
 case ("H", "h")
   is_hermitian = .TRUE.
 case ("S","s")
   is_symmetric = .TRUE.
 case("N","n")
   if (ALL(uplo(1:1) /= ["A","a"])) then
     msg = " Found symtype= "//TRIM(symtype)//", but uplo= "//TRIM(uplo)
     ABI_ERROR(msg)
   end if
   RETURN  ! Nothing to do.
 case default
   ABI_ERROR("Wrong symtype "//TRIM(symtype))
 end select

 !write(std_out,*)"is_cplx",is_cplx
 !write(std_out,*)"is_hermitian",is_hermitian

 select case (uplo(1:1))
 case ("A","a")
   ! Full global matrix has been read, nothing to do.
   return

 case ("U", "u")
   ! Only the upper triangle of the global matrix was read.
   if (is_cplx .and. is_hermitian) then
     ij_loc=0
     do jloc=1,Slk_mat%sizeb_local(2)
       do iloc=1,Slk_mat%sizeb_local(1)
         call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
         ij_loc = ij_loc+1
         if (jglob < iglob) then
          ! Diagonal elements are not forced to be real.
           Slk_mat%buffer_cplx(iloc,jloc) = DCONJG(Slk_mat%buffer_cplx(iloc,jloc))
         end if
         !if (iglob==jglob) Slk_mat%buffer_cplx(iloc,jloc) =  real(Slk_mat%buffer_cplx(iloc,jloc))
       end do
     end do
   end if

 case ("L", "l")
   ! Only the lower triangle of the global matrix was read.
   if (is_cplx .and. is_hermitian) then
     ij_loc=0
     do jloc=1,Slk_mat%sizeb_local(2)
       do iloc=1,Slk_mat%sizeb_local(1)
         call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
         ij_loc = ij_loc+1
         if (jglob>iglob) then ! diagonal elements are not forced to be real.
           Slk_mat%buffer_cplx(iloc,jloc) =  DCONJG(Slk_mat%buffer_cplx(iloc,jloc))
         end if
         !if (iglob==jglob) Slk_mat%buffer_cplx(iloc,jloc) =  real(Slk_mat%buffer_cplx(iloc,jloc))
       end do
     end do
   end if

 case default
   ABI_BUG(" Wrong uplo: "//TRIM(uplo))
 end select

end subroutine slk_symmetrize
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_single_fview_read
!! NAME
!!  slk_single_fview_read
!!
!! FUNCTION
!!  Return an MPI datatype that can be used to read a scaLAPACK distributed matrix from
!!  a binary file using MPI-IO.
!!
!! INPUTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK distribution with the local buffer.
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix is stored on disk:
!!    = "U":  Upper triangular is stored
!!    = "L":  Lower triangular is stored
!!    = "A":  Full matrix (used for general complex matrices)
!!  [is_fortran_file]=.FALSE. is C stream is used. .TRUE. for writing Fortran binary files
!!    with record markers. In this case etype is set xmpio_mpi_type_frm provided that
!!    the mpi_type of the matrix element is commensurate with xmpio_mpi_type_frm. Defaults to .TRUE.
!!
!! OUTPUT
!!  etype=Elementary data type (handle) defining the elementary unit used to access the file.
!!  slk_type=New MPI type that can be used to instantiate the MPI-IO view for the Fortran file.
!!     Note that the view assumes that the file pointer points to the FIRST Fortran record marker.
!!  offset_err=Error code. A non-zero value signals that the global matrix is too large
!!    for a single MPI-IO access (see notes below).
!!
!! NOTES
!!  With (signed) Fortran integers, the maximum size of the file that
!!  that can be read in one-shot is around 2Gb when etype is set to byte.
!!  Using a larger etype might create portability problems (real data on machines using
!!  integer*16 for the marker) since etype must be a multiple of the Fortran record marker
!!  Due to the above reason, block_displ is given in bytes and must be stored in a integer
!!  of kind XMPI_ADDRESS_KIND. If the displacement is too large, the routine returns
!!  offset_err=1 so that the caller will know that several MPI-IO reads are needed to
!!  read the local buffer.
!!
!! SOURCE

subroutine slk_single_fview_read(Slk_mat,uplo,etype,slk_type,offset_err,is_fortran_file)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: offset_err,slk_type,etype
 character(len=*),intent(in) :: uplo
 logical,optional,intent(in) :: is_fortran_file
 class(matrix_scalapack),intent(in) :: Slk_mat

!Local variables ------------------------------
!scalars
 integer :: jloc,iloc,iglob,jglob,nrows_glob,ncols_glob,mpi_err,nel
 integer :: bsize_frm,mpi_type_elm,ij_loc,bsize_etype,bsize_elm
 integer(XMPI_OFFSET_KIND) :: ijp_glob,my_offset,cpad_frm
!arrays
 character(len=500) :: msg
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

#ifdef HAVE_MPI_IO
!@matrix_scalapack
 bsize_frm = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 if (PRESENT(is_fortran_file)) then
   if (.not.is_fortran_file) bsize_frm = 0
 end if

 call slk_mat%bsize_and_type(bsize_elm, mpi_type_elm)

 ! Global dimensions.
 nrows_glob=Slk_mat%sizeb_global(1)
 ncols_glob=Slk_mat%sizeb_global(2)

 ! Number of matrix elements treated by this node.
 nel = PRODUCT(Slk_mat%sizeb_local(1:2))

 !Cannot use MPI_type_CREATE_INDEXED_BLOCK since it is not correctly implemented in several MPI libraries.
 !etype has to be set to MPI_BYTE, since the displacement in MPI structures is always in byte.
 ! ABI_WARNING("Using MPI_type_STRUCT for the MPI-IO file view")

 etype = MPI_BYTE
 call MPI_type_SIZE(etype,bsize_etype,mpi_err)

 ! Define the mapping between scaLAPACK buffer and the storage on file.
 ABI_MALLOC(block_length, (nel+2))
 ABI_MALLOC(block_displ, (nel+2))
 ABI_MALLOC(block_type, (nel+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 ! Note that the view assumes that the file pointer points to the first Fortran record marker.
 offset_err=0
 select case (uplo(1:1))
 case ("A","a") ! The entire global matrix is stored on disk. TODO can use contigous vectors for better access.
   ij_loc=0
   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       ij_loc  = ij_loc+1
       my_offset = 2*(jglob-1)*bsize_frm + bsize_frm + (jglob-1)*nrows_glob*bsize_elm + (iglob-1) * bsize_elm
       my_offset = my_offset / bsize_etype
       if (xmpio_max_address(my_offset)) offset_err=1   ! Test for possible wraparounds
       block_displ (ij_loc+1) = my_offset
       block_type  (ij_loc+1) = mpi_type_elm
       block_length(ij_loc+1) = 1
     end do
   end do

 case ("U","u") ! Only the upper triangle of the global matrix is stored on disk.
   ij_loc=0
   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       if (jglob>=iglob) then
         ijp_glob = iglob + jglob*(jglob-1)/2  ! Index for packed form
         cpad_frm = 2*(jglob-1)*bsize_frm
       else
         ijp_glob = jglob + iglob*(iglob-1)/2  ! Index for packed form
         cpad_frm = 2*(iglob-1)*bsize_frm
       end if
       ij_loc = ij_loc+1
       my_offset = cpad_frm + bsize_frm + (ijp_glob-1) * bsize_elm
       my_offset = my_offset / bsize_etype
       if (xmpio_max_address(my_offset)) offset_err=1  ! Test for possible wraparounds
       block_displ (ij_loc+1) = my_offset
       block_type  (ij_loc+1) = mpi_type_elm
       block_length(ij_loc+1) = 1
     end do
   end do

 case ("L","l") ! Only the lower triangle of the global matrix is stored on disk.
   ij_loc=0
   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       if (jglob<=iglob) then
         ijp_glob = iglob + (jglob-1)*(2*nrows_glob-jglob)/2 ! Index for packed form
         cpad_frm = 2*(jglob-1)*bsize_frm
       else
         ijp_glob = jglob + (iglob-1)*(2*nrows_glob-iglob)/2 ! Index for packed form
         cpad_frm = 2*(iglob-1)*bsize_frm
       end if
       ij_loc = ij_loc+1
       my_offset = cpad_frm + bsize_frm + (ijp_glob-1) * bsize_elm
       my_offset = my_offset / bsize_etype
       if (xmpio_max_address(my_offset)) offset_err=1   ! block_displ is usually integer*4. Test for possible wraparounds
       block_displ  (ij_loc+1) = my_offset
       block_type  (ij_loc+1) = mpi_type_elm
       block_length(ij_loc+1) = 1
     end do
   end do

   if (offset_err/=0) then  ! just warn, let the caller handle the exception.
     write(msg,"(3a)")&
      " Global position index cannot be stored in standard Fortran integer ",ch10,&
      " scaLAPACK matrix cannot be read with a single MPI-IO call ."
     ABI_WARNING(msg)
   end if

 case default
   ABI_BUG(" Wrong uplo: "//TRIM(uplo))
 end select

 block_length(nel+2)= 1
 block_displ (nel+2)= ncols_glob * (nrows_glob*bsize_elm + 2*bsize_frm) / bsize_etype
 block_type  (nel+2)= MPI_UB

 call xmpio_type_struct(nel+2,block_length,block_displ,block_type,slk_type,mpi_err)
 ABI_CHECK_MPI(mpi_err,"MPI_type_STRUCT")

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

 call MPI_type_COMMIT(slk_type,mpi_err)
 ABI_CHECK_MPI(mpi_err,"MPI_type_COMMIT")

#else
 ABI_ERROR("MPI-IO is mandatatory in slk_single_fview_read")
#endif

end subroutine slk_single_fview_read
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_single_fview_write
!! NAME
!!  slk_single_fview_write
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to write a scaLAPACK distributed matrix to
!!  a binary file using MPI-IO.
!!
!! INPUTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK distribution with the local buffer.
!!  uplo=String specifying whether only the upper or lower triangular part of the global matrix is stored on disk:
!!    = "U":  Upper triangular is stored
!!    = "L":  Lower triangular is stored
!!    = "A":  Full matrix (used for general complex matrices)
!!  [is_fortran_file]=.FALSE. is C stream is used. .TRUE. for writing Fortran binary files
!!    with record marker. In this case etype is set xmpio_mpi_type_frm provided that
!!    the mpi_type of the matrix element is commensurate with xmpio_mpi_type_frm. Defaults to .TRUE.
!!  glob_subarray(2,2) = Used to select the subarray of the global matrix. Used only when uplo="All"
!!     glob_subarray(:,1)=starting global coordinates of the subarray in each dimension
!!     (array of nonnegative integers >=1, <=array_of_sizes)
!!     glob_subarray(:,2)=Number of elements in each dimension of the subarray (array of positive integers)
!!
!! OUTPUT
!!  nelw=Number of elements to be written.
!!  etype=Elementary data type (handle) defining the elementary unit used to access the file.
!!  slk_type=New MPI type that can be used to instantiate the MPI-IO view for the Fortran file.
!!     Note that the view assumes that the file pointer points to the FIRST Fortran record marker.
!!  offset_err=Error code. A non-zero value signals that the global matrix is too large
!!    for a single MPI-IO access (see notes below).
!!
!! SIDE EFFECTS
!!  elw2slk(:,:) =
!!    input:  pointer to null().
!!    output: elw2slk(2,nelw) contains the local coordinates of the matrix elements to be written.
!!     (useful only if the upper or lower triangle of the global matrix has to be written or when
!!      uplo="all" but a global subarray is written.
!!
!! NOTES
!!  With (signed) Fortran integers, the maximum size of the file that
!!  that can be read in one-shot is around 2Gb when etype is set to byte.
!!  Using a larger etype might create portability problems (real data on machines using
!!  integer*16 for the marker) since etype must be a multiple of the Fortran record marker
!!  Due to the above reason, block_displ is given in bytes and must be stored in a Fortran
!!  integer of kind XMPI_ADDRESS_KIND. If the displacement is too large, the routine returns
!!  offset_err=1 so that the caller will know that several MPI-IO reads are needed to
!!  write the local buffer.
!!
!! SOURCE

subroutine slk_single_fview_write(Slk_mat,uplo,nelw,elw2slk,etype,slk_type,offset_err,is_fortran_file,glob_subarray)

!Arguments ------------------------------------
!scalars
 class(matrix_scalapack),intent(in) :: Slk_mat
 integer,intent(out) :: offset_err,slk_type,etype,nelw
 character(len=*),intent(in) :: uplo
 logical,optional,intent(in) :: is_fortran_file
!arrays
 integer,pointer :: elw2slk(:,:)
 integer,optional,intent(in) :: glob_subarray(2,2)

!Local variables ------------------------------
!scalars
 integer :: jloc,iloc,iglob,jglob,nrows_glob,ncols_glob,mpi_err,nel_max
 integer :: grow_min,grow_max,gcol_min,gcol_max
 integer :: bsize_frm,mpi_type_elm,ij_loc,bsize_elm
 integer(XMPI_OFFSET_KIND) :: ijp_glob,my_offset,cpad_frm
!arrays
 character(len=500) :: msg
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

#ifdef HAVE_MPI_IO
!@matrix_scalapack
 bsize_frm = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 if (PRESENT(is_fortran_file)) then
   if (.not.is_fortran_file) bsize_frm = 0
 end if

 if (PRESENT(glob_subarray).and..not.firstchar(uplo, ["A"])) then
   ABI_ERROR("glob_subarray should not be used when uplo/=All")
 end if

 call slk_mat%bsize_and_type(bsize_elm, mpi_type_elm)

 ! Global dimensions.
 nrows_glob=Slk_mat%sizeb_global(1)
 ncols_glob=Slk_mat%sizeb_global(2)

 ! Number of matrix elements treated by this node.
 nel_max = PRODUCT(Slk_mat%sizeb_local(1:2))

 ABI_MALLOC(elw2slk,(2,nel_max))
 elw2slk=0

 ! Cannot use MPI_type_CREATE_INDEXED_BLOCK since it is not correctly implemented in several MPI libraries.
 ! etype has to be set to MPI_BYTE, since the displacement in MPI structures is always in byte.
 etype = MPI_BYTE

 ! Define the mapping between scaLAPACK buffer and the storage on file.
 ABI_MALLOC(block_length, (nel_max+2))
 ABI_MALLOC(block_displ, (nel_max+2))
 ABI_MALLOC(block_type, (nel_max+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 ! Note that the view assumes that the file pointer points to the first Fortran record marker.
 offset_err=0

 select case (uplo(1:1))
 case ("A","a")
   ! The entire global matrix is written on disk. TODO can use contigous vectors for better access.
   grow_min=1; grow_max=nrows_glob
   gcol_min=1; gcol_max=ncols_glob
   if (PRESENT(glob_subarray)) then ! subarray access.
     grow_min = glob_subarray(1,1)
     gcol_min = glob_subarray(2,1)
     grow_max = grow_min + glob_subarray(1,2) -1
     gcol_max = gcol_min + glob_subarray(2,2) -1
   end if

   ij_loc=0
   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       if (iglob>=grow_min.and.iglob<=grow_max .and. &  ! glob_subarray element.
           jglob>=gcol_min.and.jglob<=gcol_max) then
         ij_loc  = ij_loc+1
         my_offset = 2*(jglob-1)*bsize_frm + bsize_frm + (jglob-1)*nrows_glob*bsize_elm + (iglob-1) * bsize_elm
         if (xmpio_max_address(my_offset)) offset_err=1   ! Test for possible wraparounds
         block_displ (ij_loc+1) = my_offset
         block_type  (ij_loc+1) = mpi_type_elm
         block_length(ij_loc+1) = 1
         elw2slk(:,ij_loc) = (/iloc,jloc/) ! useless when subarray are not used but oh well!
       end if
     end do
   end do

 case ("U","u")
   ! Only the upper triangle of the global matrix is stored on disk.
   ij_loc=0
   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       if (jglob>=iglob) then
         ijp_glob = iglob + jglob*(jglob-1)/2  ! Index for packed form
         cpad_frm = 2*(jglob-1)*bsize_frm
         ij_loc = ij_loc+1
         my_offset = cpad_frm + bsize_frm + (ijp_glob-1) * bsize_elm
         if (xmpio_max_address(my_offset)) offset_err=1   ! Test for possible wraparounds
         block_displ (ij_loc+1) = my_offset
         block_type  (ij_loc+1) = mpi_type_elm
         block_length(ij_loc+1) = 1
         elw2slk(:,ij_loc) = (/iloc,jloc/)
       end if
     end do
   end do

 case ("L","l")
   ! Only the lower triangle of the global matrix is stored on disk.
   ij_loc=0
   do jloc=1,Slk_mat%sizeb_local(2)
     do iloc=1,Slk_mat%sizeb_local(1)
       call slk_mat%loc2glob(iloc, jloc, iglob, jglob)
       if (jglob<=iglob) then
         ijp_glob = iglob + (jglob-1)*(2*nrows_glob-jglob)/2 ! Index for packed form
         cpad_frm = 2*(jglob-1)*bsize_frm
         ij_loc = ij_loc+1
         my_offset = cpad_frm + bsize_frm + (ijp_glob-1) * bsize_elm
         if (xmpio_max_address(my_offset)) offset_err=1   ! block_displ is usually integer*4. Test for possible wraparounds
         block_displ (ij_loc+1) = my_offset
         block_type  (ij_loc+1) = mpi_type_elm
         block_length(ij_loc+1) = 1
         elw2slk(:,ij_loc) = (/iloc,jloc/)
       end if
     end do
   end do

 case default
   ABI_BUG(" Wrong uplo: "//TRIM(uplo))
 end select

 if (offset_err/=0) then  ! just warn, let the caller handle the exception.
   write(msg,"(3a)")&
    "Global position index cannot be stored in standard Fortran integer ",ch10,&
    "scaLAPACK matrix cannot be read with a single MPI-IO call ."
   ABI_WARNING(msg)
 end if

 ! Final number of matrix elements that will be written by this node.
 nelw = ij_loc

 block_length(nelw+2)= 1
 block_displ (nelw+2)= ncols_glob * (nrows_glob*bsize_elm + 2*bsize_frm)
 block_type  (nelw+2)= MPI_UB

 call xmpio_type_struct(nelw+2,block_length,block_displ,block_type,slk_type,mpi_err)
 ABI_CHECK_MPI(mpi_err, "MPI_type_STRUCT")

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

 call MPI_type_COMMIT(slk_type,mpi_err)
 ABI_CHECK_MPI(mpi_err, "MPI_type_COMMIT")

#else
 ABI_ERROR("MPI-IO is mandatatory in slk_single_fview_read")
#endif

end subroutine slk_single_fview_write
!!***

!----------------------------------------------------------------------

!!****f* m_slk/slk_bsize_and_type
!! NAME
!!  slk_bsize_and_type
!!
!! FUNCTION
!!  Returns the byte size and the MPI datatype associated to the matrix elements
!!  that are stored in the ScaLAPACK_matrix
!!
!! INPUTS
!!  Slk_mat<matrix_scalapack>=Structured datatype defining the scaLAPACK distribution with the local buffer
!!
!! OUTPUT
!!  bsize_elm=Byte size of the matrix element.
!!  mpi_type_elm=MPI datatype of the matrix element.
!!
!! SOURCE

subroutine slk_bsize_and_type(Slk_mat, bsize_elm, mpi_type_elm)

!Arguments ------------------------------------
!scalars
 class(matrix_scalapack),intent(in) :: Slk_mat
 integer,intent(out) :: bsize_elm,mpi_type_elm

!Local variables ------------------------------
!scalars
 integer :: ierr
 character(len=500) :: msg

! ************************************************************************

 ! @matrix_scalapack
 ierr=0
#ifdef HAVE_MPI
 if (allocated(Slk_mat%buffer_cplx)) then
   ierr = ierr + 1
   mpi_type_elm = MPI_DOUBLE_COMPLEX
   bsize_elm    = xmpi_bsize_dpc
 end if

 if (allocated(Slk_mat%buffer_real)) then
   ierr = ierr + 1
   mpi_type_elm = MPI_DOUBLE_PRECISION
   bsize_elm    = xmpi_bsize_dp
 end if
#endif

 ! One and only one buffer should be allocated.
 if (ierr /= 1) then
   write(msg,'(a,i0)')" ScaLAPACK buffers are not allocated correctly, ierr= ",ierr
   ABI_ERROR(msg)
 end if

end subroutine slk_bsize_and_type
!!***

end module m_slk
!!***
