!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi_linalg
!! NAME
!!  m_abi_linalg
!!
!! FUNCTION
!!  management of Linear Algebra wrappers routines
!! with support of different external library (scalapack, elpa, plasma, magma, ... )
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (LNguyen,FDahm)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abi_linalg

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_xomp
  use m_slk
  use iso_c_binding
#ifdef HAVE_LINALG_ELPA
 use m_elpa
#endif
#ifdef HAVE_LINALG_PLASMA
 use plasma, except_dp => dp, except_sp => sp
#endif
 use m_time,  only : timab
 

 implicit none

 private
!!***

 ! Working arrays for eigen problem
 integer,save :: eigen_s_maxsize=0
 integer,save :: eigen_d_maxsize=0
 integer,save :: eigen_c_maxsize=0
 integer,save :: eigen_z_maxsize=0

 integer,save :: eigen_s_lwork=0
 integer,save :: eigen_d_lwork=0
 integer,save :: eigen_c_lwork=0
 integer,save :: eigen_z_lwork=0

 real(sp),save,allocatable :: eigen_c_rwork(:)
 real(dp),save,allocatable :: eigen_z_rwork(:)
 real(sp),save,allocatable :: eigen_s_work(:)
 real(dp),save,allocatable :: eigen_d_work(:)

 complex(spc),save,allocatable :: eigen_c_work(:)
 complex(dpc),save,allocatable :: eigen_z_work(:)

#ifdef HAVE_LINALG_SCALAPACK
 type(processor_scalapack),save,public :: abi_processor
 integer,save,public :: abi_communicator
 integer,save,public :: abi_complement_communicator
#endif

! This flag enables the PLASMA version in the BLAS/LAPACK wrappers.
! Set to True in linalg_allow_plasma
 public :: linalg_allow_plasma

 logical,private,save :: XPLASMA_ISON = .False.

!----------------------------------------------------------------------
!!***


!!****t* m_abi_linalg/laparams_t
!! NAME
!! laparams_t
!!
!! FUNCTION
!!  Gather the parameters and the options to be passed to the
!!  linear algebra routines.
!!
!! SOURCE

 integer,public,parameter :: LIB_LINALG=1
 integer,public,parameter :: LIB_MAGMA=2
 integer,public,parameter :: LIB_PLASMA=3
 integer,public,parameter :: LIB_SLK=4

 type,public :: laparams_t

   integer :: lib = LIB_LINALG
     ! Flag defining the library to use.

   integer :: comm = xmpi_comm_null
     ! MPI communicator (used if lib == LIB_SLK)

   !integer :: timopt, tim
   ! Options for timing.
 end type laparams_t
!!***


 !Procedures ------------------------------------
 public :: abi_linalg_init             ! Initialization and allocation routine
 public :: abi_linalg_finalize         ! CleanuUp routine
 public :: abi_linalg_eigen_setmaxsize ! Allocate work array for eigen problem
 !----------------------------------------------------------------------

 !BLAS INTERFACE
 public :: abi_xgemm
 interface abi_xgemm
    module procedure abi_zgemm_2d
    module procedure abi_d2zgemm
 end interface
 !----------------------------------------------------------------------
 public :: abi_xcopy
 interface abi_xcopy
    module procedure abi_zcopy
    module procedure abi_zcopy_1d
    module procedure abi_dcopy
    module procedure abi_dcopy_1d
    module procedure abi_dcopy_2d     ! FIXME To be removed. One can pass the base adress of the array!
    module procedure abi_dcopy_0d_1d
    module procedure abi_dcopy_1d_0d
    module procedure abi_d2zcopy_2d  ! FIXME To be removed. One can pass the base adress of the array!
    module procedure abi_z2dcopy_2d  ! FIXME To be removed. One can pass the base adress of the array!
 end interface
 !----------------------------------------------------------------------
 public :: abi_xtrsm
 interface abi_xtrsm
    module procedure abi_ztrsm
    module procedure abi_dtrsm
    module procedure abi_d2ztrsm
    !module procedure abi_d2ztrsm_3d
 end interface

 public :: abi_d2ztrsm_3d ! Used in bestwfk TODO to be Removed
 !----------------------------------------------------------------------

 !LAPACK INTERFACE
 public :: abi_xheev
 interface abi_xheev
    !module procedure abi_dheev
    !module procedure abi_cheev
    !module procedure abi_zheev
    module procedure abi_dheev_alloc
    module procedure abi_cheev_alloc
    module procedure abi_zheev_alloc
 end interface
 !----------------------------------------------------------------------
 public :: abi_xhegv
 interface abi_xhegv
    !module procedure abi_dhegv
    !module procedure abi_chegv
    !module procedure abi_zhegv
    module procedure abi_dhegv_alloc
    module procedure abi_chegv_alloc
    module procedure abi_zhegv_alloc
 end interface
 !----------------------------------------------------------------------
 public :: abi_xhpev
 interface abi_xhpev
    !module procedure abi_dhpev
    !module procedure abi_chpev
    !module procedure abi_zhpev
    module procedure abi_dhpev_alloc_1d
    module procedure abi_dhpev_alloc_2d
    module procedure abi_chpev_alloc
    module procedure abi_zhpev_alloc
 end interface
 !----------------------------------------------------------------------
 public :: abi_xhpgv
 interface abi_xhpgv
    !module procedure abi_dhpgv
    !module procedure abi_chpgv
    !module procedure abi_zhpgv
    module procedure abi_dhpgv_alloc_1d
    module procedure abi_dhpgv_alloc_2d
    module procedure abi_chpgv_alloc
    module procedure abi_zhpgv_alloc
 end interface
 !----------------------------------------------------------------------
 interface abi_xpotrf
    !module procedure abi_dpotrf
    module procedure abi_d2zpotrf
    module procedure abi_zpotrf_2d
 end interface
 !----------------------------------------------------------------------
 public :: abi_xorthonormalize
 interface abi_xorthonormalize
    module procedure xorthonormalize
    module procedure zorthonormalize
 end interface

 ! This version operates on arrays in which the real and the imaginary part
 ! are packed together (real parts first, them imaginary parts), used when gamma-point and istwfk=2
 public :: ortho_reim
 !----------------------------------------------------------------------

#ifndef HAVE_GPU_CUDA
 !dummy routines replace gpu helper routines
 public :: alloc_on_gpu
 public :: copy_from_gpu
 public :: copy_on_gpu
 public :: dealloc_on_gpu
 public :: gpu_linalg_init
 public :: gpu_linalg_shutdown
 public :: gpu_xgemm
 public :: gpu_xtrsm
#endif

 public :: gpu_xorthonormalize

 logical,external :: LSAME

 ! Timab slots, used if we want to profile BLAS calls
 ! Fine-grained profiling, must be enabled with the CPP option DEV_LINALG_TIMING
 ! For the time being, I use the same slots employed in lobpcgwf although
 ! one should define specialized entries.If the index of the slot 0, no profiling is done.

 integer,parameter,private :: TIMAB_XCOPY=584
 integer,parameter,private :: TIMAB_XGEMM=532
 integer,parameter,private :: TIMAB_XORTHO=535
 integer,parameter,private :: TIMAB_XEIGEN=587
 integer,parameter,private :: TIMAB_XPRECO=536
 integer,parameter,private :: TIMAB_WFCOPY=584
 integer,parameter,private :: TIMAB_XTRSM=535

#define DEV_LINALG_TIMING 1

! Support for [Z,C]GEMM3M routines
 logical,save,private :: XGEMM3M_ISON=.False.
 ! True if [Z,C]GEMM3M can be used (can be set with linalg_allow_gemm3m)

 public :: linalg_allow_gemm3m

 ! Thresholds for the activation of [Z,C]GEMM3M
 integer,parameter,private :: ZGEMM3M_LIMIT=325000
 integer,parameter,private :: CGEMM3M_LIMIT=200000

! Handy macros
#ifdef HAVE_LINALG_GEMM3M
#define _ZGEMM3M ZGEMM3M
#define _CGEMM3M CGEMM3M
#else
#define _ZGEMM3M ZGEMM
#define _CGEMM3M CGEMM
#endif


CONTAINS  !===========================================================
!!***

!!****f* m_abi_linalg/abi_linalg_init
!! NAME
!! abi_linalg_init
!!
!! FUNCTION
!! Initalization of linear algebra environnement
!!
!! INPUTS
!!
!! PARENTS
!!      compute_kgb_indicator,driver
!!
!! CHILDREN
!!
!! SOURCE
!!
 subroutine abi_linalg_init(comm_scalapack,eigen_group_size,max_eigen_pb_size,my_rank,&
&                           only_scalapack) ! optional parameter

#if defined HAVE_MPI2
  use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_linalg_init'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: comm_scalapack,eigen_group_size,max_eigen_pb_size,my_rank
 logical,intent(in),optional :: only_scalapack

!Local variables ------------------------------
#ifdef HAVE_LINALG_SCALAPACK
 integer :: abi_info1,rank,commsize
 integer :: sizecart(2)
 logical :: periodic(2), reorder, keepdim(2)
 integer :: commcart
#endif
#ifdef HAVE_LINALG_PLASMA
 integer :: abi_info2,core_id,num_cores,num_cores_node
 integer,allocatable :: affinity(:)
#endif

!******************************************************************

#ifdef HAVE_LINALG_SCALAPACK
!Scalapack initalization
 if (eigen_group_size>0) then
   abi_communicator = comm_scalapack
   rank=xmpi_comm_rank(comm_scalapack)

   ! We create abi_communicator using a cartesian grid, and store its complement
   commsize = MIN(eigen_group_size, xmpi_comm_size(comm_scalapack))
   sizecart = (/commsize, xmpi_comm_size(comm_scalapack)/commsize/)
   periodic = (/.true.,.true./)
   reorder = .false.
   call MPI_CART_CREATE(comm_scalapack,2,sizecart,periodic,reorder,commcart,abi_info1)
   keepdim = (/.true., .false./)
   call MPI_CART_SUB(commcart, keepdim, abi_communicator,abi_info1)
   keepdim = (/.false., .true./)
   call MPI_CART_SUB(commcart, keepdim, abi_complement_communicator,abi_info1)

   call init_scalapack(abi_processor,abi_communicator)
 else
   abi_communicator=xmpi_comm_null
   abi_complement_communicator = xmpi_comm_null
 end if
#endif
 if (present(only_scalapack)) then
   if (only_scalapack) return
 end if

#ifdef HAVE_LINALG_ELPA
 call elpa_func_init()
#endif

#ifdef HAVE_LINALG_PLASMA
!Plasma Initialization
!Because use of hybrid use of mpi+openmp+plasma,
!  we need to set manually the thread bindings policy
!  to avoid conflicts between mpi process due to plasma
 num_cores=xomp_get_max_threads()
 num_cores_node=xomp_get_num_cores_node()

 if (num_cores_node == 0) then ! This means that OMP is not enabled.
   num_cores_node = 1
   MSG_WARNING("You are using PLASMA but OpenMP is not enabled in Abinit!")
 end if

 ABI_ALLOCATE(affinity,(num_cores))
 do core_id =1,num_cores
   affinity(core_id) = MOD(my_rank*num_cores + (core_id-1), num_cores_node)
 end do
 call PLASMA_Init_Affinity(num_cores,affinity(1),abi_info2)
 ABI_DEALLOCATE(affinity)
 XPLASMA_ISON = .True.
#endif

#ifdef HAVE_LINALG_MAGMA
#ifdef HAVE_LINALG_MAGMA_15
 call magmaf_init()
#endif
#endif

#ifdef HAVE_GPU_CUDA
!Cublas initialization
 call gpu_linalg_init()
#endif

!Allocation of work sizes
 eigen_s_maxsize = 0
 eigen_d_maxsize = 0
 eigen_c_maxsize = 0
 eigen_z_maxsize = 0
 eigen_s_lwork   = 0
 eigen_d_lwork   = 0
 eigen_c_lwork   = 0
 eigen_z_lwork   = 0

 call abi_linalg_eigen_setmaxsize(max_eigen_pb_size)

 return

 ABI_UNUSED(comm_scalapack)
 ABI_UNUSED(eigen_group_size)
 ABI_UNUSED(my_rank)

 end subroutine abi_linalg_init
!!***

!!****f* m_abi_linalg/abi_linalg_eigen_setmaxsize
!! NAME
!! abi_linalg_eigen_setmaxsize
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_abi_linalg
!!
!! CHILDREN
!!
!! SOURCE
!!
 subroutine abi_linalg_eigen_setmaxsize(max_eigen_pb_size)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_linalg_eigen_setmaxsize'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: max_eigen_pb_size

!******************************************************************

!Single precision work arrays (real and complex)
 eigen_s_maxsize = max_eigen_pb_size
 eigen_s_lwork = max(1,3*eigen_s_maxsize)
 if(allocated(eigen_s_work)) then
   ABI_DEALLOCATE(eigen_s_work)
 end if
 ABI_ALLOCATE(eigen_s_work,(max(1,eigen_s_lwork)))

 eigen_c_maxsize = max_eigen_pb_size
 eigen_c_lwork = max(1,2*eigen_c_maxsize-1)
 if(allocated(eigen_c_work)) then
   ABI_DEALLOCATE(eigen_c_work)
 end if

 ABI_ALLOCATE(eigen_c_work,(max(1,eigen_c_lwork)))
 if(allocated(eigen_c_rwork)) then
   ABI_DEALLOCATE(eigen_c_rwork)
 end if
 ABI_ALLOCATE(eigen_c_rwork,(max(1, 3*eigen_c_maxsize-2)))

!Double precision work arrays (real and complex)
 eigen_d_maxsize = max_eigen_pb_size
 eigen_d_lwork = max(eigen_d_lwork,2*(2*eigen_d_maxsize - 1)) !for zh[ep][eg]v
 eigen_d_lwork = max(eigen_d_lwork,(3*eigen_d_maxsize))       !for ds[ep][eg]v
#ifdef HAVE_LINALG_MAGMA
 eigen_d_lwork = max(eigen_d_lwork,2*(2*(eigen_d_maxsize**2) + 6*eigen_d_maxsize +1))!for magma_zhe[eg]v
 eigen_d_lwork = max(eigen_d_lwork, eigen_d_maxsize**2 + 33*eigen_d_maxsize)         !for magma_dsy[eg]v
#endif
 if(allocated(eigen_d_work)) then
   ABI_DEALLOCATE(eigen_d_work)
 end if
 ABI_ALLOCATE(eigen_d_work,(max(1,eigen_d_lwork)))

 eigen_z_maxsize = max_eigen_pb_size
 eigen_z_lwork = max(1,2*eigen_z_maxsize-1)
 if(allocated(eigen_z_work)) then
   ABI_DEALLOCATE(eigen_z_work)
 end if
 ABI_ALLOCATE(eigen_z_work,(max(1,eigen_z_lwork)))
 if(allocated(eigen_z_rwork)) then
   ABI_DEALLOCATE(eigen_z_rwork)
 end if

#ifdef HAVE_LINALG_MAGMA
 ABI_ALLOCATE(eigen_z_rwork,(max(1, 2*(eigen_z_maxsize**2) + 5*eigen_z_maxsize + 1)))
#else
 ABI_ALLOCATE(eigen_z_rwork,(max(1, 3*eigen_z_maxsize-2)))
#endif

 end subroutine abi_linalg_eigen_setmaxsize
!!***

!!****f* m_abi_linalg/abi_linalg_finalize
!! NAME
!! abi_linalg_finalize
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      compute_kgb_indicator,driver
!!
!! CHILDREN
!!
!! SOURCE
!!
 subroutine abi_linalg_finalize(only_scalapack) ! optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_linalg_finalize'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical,intent(in),optional :: only_scalapack

!Local variables ------------------------------
#ifdef HAVE_LINALG_PLASMA
 integer :: info
#endif

!******************************************************************

#ifdef HAVE_LINALG_SCALAPACK
 if (abi_communicator/=xmpi_comm_null) then
   call xmpi_comm_free(abi_communicator)
   call end_scalapack(abi_processor)
 end if
#endif
 if (present(only_scalapack)) then
   if (only_scalapack) return
 end if

#ifdef HAVE_LINALG_ELPA
 call elpa_func_uninit()
#endif

#ifdef HAVE_LINALG_PLASMA
 call PLASMA_Finalize(info)
#endif

#ifdef HAVE_LINALG_MAGMA
#ifdef HAVE_LINALG_MAGMA_15
 call magmaf_finalize()
#endif
#endif

#ifdef HAVE_GPU_CUDA
 call gpu_linalg_shutdown()
#endif

!Memory freeing
 if(allocated(eigen_s_work)) then
   ABI_DEALLOCATE(eigen_s_work)
 end if
 if(allocated(eigen_d_work)) then
   ABI_DEALLOCATE(eigen_d_work)
 end if
 if(allocated(eigen_c_work)) then
   ABI_DEALLOCATE(eigen_c_work)
 end if
 if(allocated(eigen_c_rwork)) then
   ABI_DEALLOCATE(eigen_c_rwork)
 end if
 if(allocated(eigen_z_work)) then
   ABI_DEALLOCATE(eigen_z_work)
 end if
 if(allocated(eigen_z_rwork)) then
   ABI_DEALLOCATE(eigen_z_rwork)
 end if

 eigen_s_maxsize = 0
 eigen_d_maxsize = 0
 eigen_c_maxsize = 0
 eigen_z_maxsize = 0
 eigen_s_lwork   = 0
 eigen_d_lwork   = 0
 eigen_c_lwork   = 0
 eigen_z_lwork   = 0

 end subroutine abi_linalg_finalize
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/linalg_allow_gemm3m
!! NAME
!!
!! FUNCTION
!!  Programmatic interface to enable the use of [Z,C]GEMM3M calls
!!
!! SOURCE

subroutine linalg_allow_gemm3m(bool)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linalg_allow_gemm3m'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: bool

! *************************************************************************

 XGEMM3M_ISON = bool

end subroutine linalg_allow_gemm3m
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/use_zgemm3m
!! NAME
!!  use_zgemm3m
!!
!! FUNCTION
!!  Enable the use of ZGEMM3M
!!
!! PARENTS
!!
!! NOTES
!!  The CGEMM3M and ZGEMM3M routines use an algorithm requiring 3 real matrix
!!  multiplications and 5 real matrix additions to compute the complex matrix
!!  product; CGEMM(3S) and ZGEMM(3S) use 4 real matrix multiplications and 2
!!  real matrix additions. Because the matrix multiplication time is usually
!!  the limiting performance factor in these routines, CGEMM3M and ZGEMM3M
!!  may run up to 33 percent faster than CGEMM and ZGEMM.  Because of other
!!  overhead associated with the 3M routines, however, these performance
!!  improvements may not always be realized.  For example, on one processor
!!  the 3M routines will generally run more slowly than the standard complex
!!  matrix multiplication routines when m * n * k < FACTOR, where m, n, and k
!!  are the input matrix dimensions and FACTOR is approximately 200000 for
!!  CGEMM3M and 325000 for ZGEMM3M.
!!  from: http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?coll=0650&db=man&raw=1&fname=/usr/share/catman/p_man/cat3/SCSL/ZGEMM3M.z
!!
!! SOURCE

pure logical function use_zgemm3m(m,n,k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'use_zgemm3m'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m,n,k

! *************************************************************************

 use_zgemm3m = .False.
 if (XGEMM3M_ISON) use_zgemm3m = ((m * n * k) > ZGEMM3M_LIMIT)

#ifndef HAVE_LINALG_GEMM3M
 use_zgemm3m = .False.
#endif

end function use_zgemm3m
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/use_cgemm3m
!! NAME
!!  use_cgemm3m
!!
!! FUNCTION
!!  Enable the use of CGEMM3M
!!
!! PARENTS
!!
!! NOTES
!!  See use_zgemm3m
!!
!! SOURCE

pure logical function use_cgemm3m(m,n,k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'use_cgemm3m'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m,n,k

! *************************************************************************

 use_cgemm3m = .False.
 if (XGEMM3M_ISON) use_cgemm3m = ((m * n * k) > CGEMM3M_LIMIT)
#ifndef HAVE_LINALG_GEMM3M
 use_cgemm3m = .False.
#endif

end function use_cgemm3m
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/linalg_allow_plasma
!! NAME
!!
!! FUNCTION
!!  Programmatic interface to enable the use of PLASMA
!!  False to disable PLASMA version.
!!
!! SOURCE

subroutine linalg_allow_plasma(bool)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linalg_allow_plasma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: bool

! *************************************************************************

 XPLASMA_ISON = bool
#ifndef HAVE_LINALG_PLASMA
 ! Just to be on the safe-side.
 ! I have to use a weird set of branches to make abirules happy in the BLAS/LAPACK
 ! wrappers, and one cannot set XPLASMA_MODE to .True. if PLASMA is not available.
 XPLASMA_ISON = .False.
#endif

end subroutine linalg_allow_plasma
!!***

#ifdef HAVE_LINALG_PLASMA

!!****f* m_abi_linalg/uplo_plasma
!! NAME
!!
!! FUNCTION
!!  Convert uplo character to PLASMA integer
!!
!! SOURCE

integer function uplo_plasma(uplo)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'uplo_plasma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: uplo

! *************************************************************************

 if (LSAME(uplo,'U')) then
    uplo_plasma = PlasmaUpper
 else
    uplo_plasma = PlasmaLower
 end if

end function uplo_plasma
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/trans_plasma
!! NAME
!!
!! FUNCTION
!!  Convert trans character to PLASMA integer
!!
!! PARENTS
!!
!! SOURCE

integer function trans_plasma(trans)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'trans_plasma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: trans

! *************************************************************************

 if (LSAME(trans,'C')) then
   trans_plasma = PlasmaConjTrans
 else if (LSAME(trans,'T')) then
   trans_plasma = PlasmaTrans
 else
   trans_plasma = PlasmaNoTrans
 end if

end function trans_plasma
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/side_plasma
!! NAME
!!
!! FUNCTION
!!  Convert side character to PLASMA integer
!!
!! PARENTS
!!
!! SOURCE

integer function side_plasma(side)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'side_plasma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: side

! *************************************************************************

 if(LSAME(side,'L')) then
    side_plasma = PlasmaLeft
 else
    side_plasma = PlasmaRight
 end if

end function side_plasma
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/diag_plasma
!! NAME
!!
!! FUNCTION
!!  Convert diag character to PLASMA integer
!!
!! PARENTS
!!
!! SOURCE

integer function diag_plasma(diag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diag_plasma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: diag

! *************************************************************************

 if (LSAME(diag,'U')) then
   diag_plasma = PlasmaUnit
 else
   diag_plasma = PlasmaNonUnit
 end if

end function diag_plasma
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/jobz_plasma
!! NAME
!!
!! FUNCTION
!!  Convert jobz character to PLASMA integer
!!
!! PARENTS
!!
!! SOURCE

integer function jobz_plasma(jobz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jobz_plasma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: jobz

! *************************************************************************

 if (LSAME(jobz,'N')) then
   jobz_plasma = PlasmaNoVec
 else
   jobz_plasma = PlasmaVec
 end if

end function jobz_plasma
!!***

#endif

! Include files providing wrappers for some of the most commonly used BLAS & LAPACK routines

! *********************************************************************!
! ******************* BLAS_XGEMM interface ****************************!
! *********************************************************************!
#include "abi_xgemm.f90"

! *********************************************************************!
! ******************** BLAS_XCOPY interface ***************************!
! *********************************************************************!
#include "abi_xcopy.f90"

! *********************************************************************!
! ******************** BLAS_XTRSM interface ***************************!
! *********************************************************************!
#include "abi_xtrsm.f90"

! *********************************************************************!
! ********************* LAPACK XHEEV  interface ***********************!
! *********************************************************************!
#include "abi_xheev.f90"

! *********************************************************************!
! ********************* LAPACK XHEGV  interface ***********************!
! *********************************************************************!
#include "abi_xhegv.f90"

! *********************************************************************!
! ********************* LAPACK XHPEV  interface ***********************!
! *********************************************************************!
#include "abi_xhpev.f90"

! *********************************************************************!
! ********************* LAPACK XHPGV  interface ***********************!
! *********************************************************************!
#include "abi_xhpgv.f90"

! *********************************************************************!
! ****************** LAPACK XPOTRF  interface *************************!
! *********************************************************************!
#include "abi_xpotrf.f90"

! *********************************************************************!
! ************  ABINIT Orthonormalization interface *******************!
! *********************************************************************!
#include "abi_xorthonormalize.f90"

! *********************************************************************!
! ******************    GPU LINALG interface  *************************!
! *********************************************************************!
#include "abi_gpu_linalg.f90"

!----------------------------------------------------------------------

end module m_abi_linalg
!!***
