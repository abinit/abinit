!!****m* ABINIT/m_abi_linalg
!! NAME
!!  m_abi_linalg
!!
!! FUNCTION
!!  management of Linear Algebra wrappers routines
!!  with support of different external library (scalapack, elpa, plasma, magma, ... )
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2020 ABINIT group (LNguyen,FDahm,MT)
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

#if defined HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

#if defined HAVE_MPI2
 use mpi
#endif

 use m_time,  only : timab

 implicit none

 private
!!***

!This flag is ON if abi_linalg functions are in use (if this module is in use)
 logical :: abi_linalg_in_use=.true.

!These flags enable the different versions in the BLAS/LAPACK wrappers
 logical,private,save :: ABI_LINALG_SCALAPACK_ISON=.False.
 logical,private,save :: ABI_LINALG_MAGMA_ISON=.False.
 logical,private,save :: ABI_LINALG_PLASMA_ISON=.False.

!Working arrays for eigen problem

 logical,save :: lapack_single_precision=.false.
 logical,save :: lapack_double_precision=.false.
 logical,save :: lapack_full_storage    =.false.
 logical,save :: lapack_packed_storage  =.false.
 logical,save :: lapack_divide_conquer  =.false.

 integer,save :: eigen_s_maxsize=0
 integer,save :: eigen_d_maxsize=0
 integer,save :: eigen_c_maxsize=0
 integer,save :: eigen_z_maxsize=0
 integer,save :: eigen_s_lwork=0
 integer,save :: eigen_d_lwork=0
 integer,save :: eigen_c_lwork=0
 integer,save :: eigen_z_lwork=0
 integer,save :: eigen_c_lrwork=0
 integer,save :: eigen_z_lrwork=0
 integer,save :: eigen_liwork=0

 integer,save,target,allocatable :: eigen_iwork(:)
 real(sp),save,target,allocatable :: eigen_c_rwork(:)
 real(dp),save,target,allocatable :: eigen_z_rwork(:)
 real(sp),save,target,allocatable :: eigen_s_work(:)
 real(dp),save,target,allocatable :: eigen_d_work(:)
 complex(spc),save,target,allocatable :: eigen_c_work(:)
 complex(dpc),save,target,allocatable :: eigen_z_work(:)

 integer,save,private :: slk_minsize=1
 integer,save,public :: slk_communicator=xmpi_comm_null
 integer,save,public :: slk_complement_communicator=xmpi_comm_null
#ifdef HAVE_LINALG_SCALAPACK
 type(processor_scalapack),save,public :: slk_processor
#endif

!Plasma can be activated via command line
!Use XPLASMA_ISON flag and modifiy it with linalg_allow_plasma
 logical,private,save :: XPLASMA_ISON=.false.
 public :: linalg_allow_plasma
#ifdef HAVE_LINALG_PLASMA
 type(c_ptr) :: plasma_work
#endif

!----------------------------------------------------------------------
!!***

 !Procedures ------------------------------------
 public :: abi_linalg_init          ! Initialization routine
 public :: abi_linalg_finalize      ! CleanuUp routine
 public :: abi_linalg_work_allocate ! Allocate work arrays
 !----------------------------------------------------------------------

!BLAS INTERFACE
 !public :: abi_zgemm
 public :: abi_xgemm

 interface abi_xgemm
    module procedure abi_zgemm_2d
    module procedure abi_d2zgemm
 end interface abi_xgemm

 public :: abi_zgemm
 public :: abi_zgemm_2d
 public :: abi_zgemm_2r
 interface abi_zgemm  ! No x_cplx stuff here!
    module procedure abi_zgemm_2d
    module procedure abi_zgemm_2r
 end interface abi_zgemm

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
 end interface abi_xcopy

 !----------------------------------------------------------------------
 public :: abi_xtrsm
 interface abi_xtrsm
    module procedure abi_ztrsm
    module procedure abi_dtrsm
    module procedure abi_d2ztrsm
    !module procedure abi_d2ztrsm_3d
 end interface abi_xtrsm

 public :: abi_d2ztrsm_3d ! Used in bestwfk TODO to be Removed
 !----------------------------------------------------------------------

!LAPACK INTERFACE
 public :: abi_xheev
 interface abi_xheev
    module procedure abi_dheev
    module procedure abi_cheev
    module procedure abi_zheev
 end interface
 !----------------------------------------------------------------------
 public :: abi_xhegv
 interface abi_xhegv
    module procedure abi_dhegv
    module procedure abi_chegv
    module procedure abi_zhegv
 end interface
 !----------------------------------------------------------------------
 public :: abi_xhpev
 interface abi_xhpev
    module procedure abi_dhpev
    module procedure abi_chpev
    module procedure abi_zhpev
 end interface
 !----------------------------------------------------------------------
 public :: abi_xhpgv
 interface abi_xhpgv
    module procedure abi_dhpgv
    module procedure abi_chpgv
    module procedure abi_zhpgv
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

! Define this variable to activate timing routines
!#define DEV_LINALG_TIMING 1

! Support for [Z,C]GEMM3M routines
 logical,save,private :: XGEMM3M_ISON = .False.
 !logical,save,private :: XGEMM3M_ISON = .True.
 ! True if [Z,C]GEMM3M can be used (can be set with linalg_allow_gemm3m)

 public :: linalg_allow_gemm3m

 ! Thresholds for the activation of [Z,C]GEMM3M
 integer,parameter,private :: ZGEMM3M_LIMIT = 325000
 integer,parameter,private :: CGEMM3M_LIMIT = 200000

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
!! max_eigen_pb_size= max. size of eigenproblem during calculation
!! optdriver= type of calculation (ground-state, response function, GW, ...)
!! wfoptalg= wave functions optimization algorithm (CG, LOBPCG, CHEBFI, ...)
!! paral_kgb= 1 if (k,g,b) parallelism is on
!! use_gpu_cuda= 1 if Cuda (GPU) is on
!! use_slk= 1 if use of Scalapack is on
!! np_slk= max. number of processes to be used in Scalapack calls
!! comm_scalapack= global communicator to be used in case of Scalapack
!!
!! PARENTS
!!      m_driver,m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE
!!

 subroutine abi_linalg_init(max_eigen_pb_size,optdriver,wfoptalg,paral_kgb,&
&                           use_gpu_cuda,use_slk,np_slk,comm_scalapack)

!Arguments ------------------------------------
 integer,intent(in) :: max_eigen_pb_size
 integer,intent(in) :: optdriver,wfoptalg,paral_kgb
 integer,intent(in) :: comm_scalapack,np_slk
 integer,intent(in) :: use_gpu_cuda,use_slk

!Local variables ------------------------------
 integer :: max_eigen_pb_size_eff=0
 logical :: need_work_space=.true.
#ifdef HAVE_LINALG_SCALAPACK
 integer :: abi_info1,rank,commsize,commcart,sizecart(2)
 logical :: reorder,periodic(2),keepdim(2)
#endif
#ifdef HAVE_LINALG_PLASMA
 integer :: abi_info2,core_id,rank
 integer :: num_cores=0,num_cores_node=0
 integer,allocatable :: affinity(:)
#endif

!******************************************************************

!Use only abi_linalg in case of GS calculations
 abi_linalg_in_use=(optdriver==RUNL_GSTATE.or.optdriver==RUNL_GWLS)

 max_eigen_pb_size_eff=0
 lapack_single_precision=.false.
 lapack_double_precision=.false.
 lapack_full_storage    =.false.
 lapack_packed_storage  =.false.
 lapack_divide_conquer  =.false.
 eigen_s_maxsize=0 ; eigen_d_maxsize=0
 eigen_c_maxsize=0 ; eigen_z_maxsize=0
 eigen_s_lwork=0   ; eigen_d_lwork=0
 eigen_c_lwork=0   ; eigen_z_lwork=0
 eigen_c_lrwork=0  ; eigen_z_lrwork=0
 eigen_liwork=0
 ABI_LINALG_SCALAPACK_ISON=.False.
 ABI_LINALG_MAGMA_ISON=.False.
 ABI_LINALG_PLASMA_ISON=.False.
 slk_communicator=xmpi_comm_null
 slk_complement_communicator=xmpi_comm_null
 slk_minsize=1

!Exit here if we don't use this abi_linalg module
 if (.not.abi_linalg_in_use) return

!Set Lapack parameters
 max_eigen_pb_size_eff=max_eigen_pb_size
 if (wfoptalg==4.or.wfoptalg==14.or.use_gpu_cuda/=0) max_eigen_pb_size_eff=3*max_eigen_pb_size_eff
 lapack_full_storage=(wfoptalg==4.or.wfoptalg==14.or.use_gpu_cuda/=0)
 lapack_packed_storage=.true.
 lapack_single_precision=.false.
 lapack_double_precision=.true.
 lapack_divide_conquer=.false.

!Set maximum sizes
 eigen_s_maxsize = max_eigen_pb_size_eff
 eigen_d_maxsize = max_eigen_pb_size_eff
 eigen_c_maxsize = max_eigen_pb_size_eff
 eigen_z_maxsize = max_eigen_pb_size_eff
 need_work_space=.true.

#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI
 if ((paral_kgb==1.or.use_slk==1).and.np_slk>0) then
   rank=xmpi_comm_rank(comm_scalapack)
   ! We create slk_communicator using a cartesian grid, and store its complement
   commsize = MIN(np_slk, xmpi_comm_size(comm_scalapack))
   sizecart = (/commsize, xmpi_comm_size(comm_scalapack)/commsize/)
   periodic = (/.true.,.true./) ; reorder = .false.
   call MPI_CART_CREATE(comm_scalapack,2,sizecart,periodic,reorder,commcart,abi_info1)
   keepdim = (/.true., .false./)
   call MPI_CART_SUB(commcart, keepdim, slk_communicator,abi_info1)
   keepdim = (/.false., .true./)
   call MPI_CART_SUB(commcart, keepdim, slk_complement_communicator,abi_info1)
   call init_scalapack(slk_processor,slk_communicator)
   slk_minsize=maxval(slk_processor%grid%dims(1:2))
   need_work_space=(use_slk/=1) ! In this case we never use the work arrays
   ABI_LINALG_SCALAPACK_ISON = .true.
 end if
#else
 ABI_UNUSED(comm_scalapack)
 ABI_UNUSED(paral_kgb)
 ABI_UNUSED(use_slk)
 ABI_UNUSED(np_slk)
#endif

#ifdef HAVE_LINALG_ELPA
 call elpa_func_init()
#endif

#ifdef HAVE_LINALG_PLASMA
!Plasma Initialization
!Because use of hybrid use of mpi+openmp+plasma,
!  we need to set manually the thread bindings policy
!  to avoid conflicts between mpi process due to plasma
 if (XPLASMA_ISON) then
   num_cores=xomp_get_max_threads()
   num_cores_node=xomp_get_num_cores_node()
   rank=xmpi_comm_rank(xmpi_world)
   if (num_cores_node == 0) then ! This means that OMP is not enabled.
     num_cores_node = 1
     MSG_WARNING("You are using PLASMA but OpenMP is not enabled in Abinit!")
   end if
   ABI_ALLOCATE(affinity,(num_cores))
   do core_id =1,num_cores
     affinity(core_id) = MOD(rank*num_cores + (core_id-1), num_cores_node)
   end do
   call PLASMA_Init_Affinity(num_cores,affinity(1),abi_info2)
   ABI_DEALLOCATE(affinity)
   ABI_LINALG_PLASMA_ISON = .True.
   lapack_divide_conquer=.true.
 end if
#endif

#ifdef HAVE_LINALG_MAGMA
#ifdef HAVE_LINALG_MAGMA_15
 call magmaf_init()
 ABI_LINALG_MAGMA_ISON = .true.
 lapack_divide_conquer=.true.
#endif
#endif

#ifdef HAVE_GPU_CUDA
!Cublas initialization
 call gpu_linalg_init()
#endif

 if (need_work_space) call abi_linalg_work_allocate()

 end subroutine abi_linalg_init
!!***

!!****f* m_abi_linalg/abi_linalg_work_allocate
!! NAME
!! abi_linalg_work_allocate
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
 subroutine abi_linalg_work_allocate()

!Arguments ------------------------------------

!Local variables ------------------------------
#ifdef HAVE_LINALG_MAGMA
 integer :: nb
 integer :: magmaf_get_ssytrd_nb
 integer :: magmaf_get_dsytrd_nb
 integer :: magmaf_get_chetrd_nb
 integer :: magmaf_get_zhetrd_nb
#endif

!******************************************************************

!Single precision WORK
 eigen_s_lwork = 0
 if (eigen_s_maxsize>0) then
   if (lapack_single_precision) then
     if (lapack_full_storage) then
       eigen_s_lwork = max(eigen_s_lwork,3*eigen_s_maxsize-1) ! SSYEV, SSYGV
     end if
     if (lapack_packed_storage) then
       eigen_s_lwork = max(eigen_s_lwork,3*eigen_s_maxsize) ! SSPEV[D], SSPGV[D]
     end if
     if (lapack_divide_conquer) then
       eigen_s_lwork = max(eigen_s_lwork,1+6*eigen_s_maxsize+2*eigen_s_maxsize**2) ! SSYEVD, SSYGVD
     end if
     if (ABI_LINALG_MAGMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
#if defined HAVE_LINALG_MAGMA
         nb=magmaf_get_ssytrd_nb(eigen_s_maxsize)
         eigen_s_lwork = max(eigen_s_lwork,eigen_s_maxsize*(nb+2)) ! MAGMAF_SSYEVD, MAGMAF_SSYGVD
#endif
       end if
     end if
     if (ABI_LINALG_PLASMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
         eigen_s_lwork = max(eigen_s_lwork,eigen_s_maxsize**2) ! PLASMA_SSYEV
       end if
     end if
   end if
 end if
 ABI_SFREE(eigen_s_work)
 ABI_ALLOCATE(eigen_s_work,(eigen_s_lwork))

!Double precision WORK
 eigen_d_lwork = 0
 if (eigen_d_maxsize>0) then
   if (lapack_double_precision) then
     if (lapack_full_storage) then
       eigen_d_lwork = max(eigen_d_lwork,3*eigen_d_maxsize-1) ! DSYEV, DSYGV
     end if
     if (lapack_packed_storage) then
       eigen_d_lwork = max(eigen_d_lwork,3*eigen_d_maxsize) ! DSPEV[D], DSPGV[D]
     end if
     if (lapack_divide_conquer) then
       eigen_d_lwork = max(eigen_d_lwork,1+6*eigen_d_maxsize+2*eigen_d_maxsize**2) ! DSYEVD, DSYGVD
     end if
     if (ABI_LINALG_MAGMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
#if defined HAVE_LINALG_MAGMA
         nb=magmaf_get_dsytrd_nb(eigen_d_maxsize)
         eigen_d_lwork = max(eigen_d_lwork,eigen_d_maxsize*(nb+2)) ! MAGMAF_DSYEVD, MAGMAF_DSYGVD
#endif
       end if
     end if
     if (ABI_LINALG_PLASMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
         eigen_d_lwork = max(eigen_d_lwork,eigen_d_maxsize**2) ! PLASMA_DSYEV
       end if
     end if
   end if
 end if
 ABI_SFREE(eigen_d_work)
 ABI_ALLOCATE(eigen_d_work,(eigen_d_lwork))

!Single complex WORK
 eigen_c_lwork = 0
 if (eigen_c_maxsize>0) then
   if (lapack_single_precision) then
     if (lapack_full_storage) then
       eigen_c_lwork = max(eigen_c_lwork,2*eigen_c_maxsize-1) ! CHEEV, CHEGV
     end if
     if (lapack_packed_storage) then
       eigen_c_lwork = max(eigen_c_lwork,2*eigen_c_maxsize) ! CHPEV[D], CHPGV[D]
     end if
     if (lapack_divide_conquer) then
       eigen_c_lwork = max(eigen_c_lwork,2*eigen_c_maxsize+eigen_c_maxsize**2) ! CHEEVD, CHEGVD
     end if
     if (ABI_LINALG_MAGMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
#if defined HAVE_LINALG_MAGMA
         nb=magmaf_get_chetrd_nb(eigen_c_maxsize)
         eigen_c_lwork = max(eigen_c_lwork,eigen_c_maxsize*(nb+1)) ! MAGMAF_CHEEVD, MAGMAF_CHEGVD
#endif
       end if
     end if
     if (ABI_LINALG_PLASMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
         eigen_c_lwork = max(eigen_c_lwork,eigen_c_maxsize**2) ! PLASMA_CHEEV
       end if
     end if
   end if
 end if
 ABI_SFREE(eigen_c_work)
 ABI_ALLOCATE(eigen_c_work,(eigen_c_lwork))

!Double complex WORK
 eigen_z_lwork = 0
 if (eigen_z_maxsize>0) then
   if (lapack_double_precision) then
     if (lapack_full_storage) then
       eigen_z_lwork = max(eigen_z_lwork,2*eigen_z_maxsize-1) ! ZHEEV, ZHEGV
     end if
     if (lapack_packed_storage) then
       eigen_z_lwork = max(eigen_z_lwork,2*eigen_z_maxsize) ! ZHPEV[D], ZHPGV[D]
     end if
     if (lapack_divide_conquer) then
       eigen_z_lwork = max(eigen_z_lwork,2*eigen_z_maxsize+eigen_z_maxsize**2) ! ZHEEVD, ZHEGVD
     end if
     if (ABI_LINALG_MAGMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
#if defined HAVE_LINALG_MAGMA
         nb=magmaf_get_zhetrd_nb(eigen_z_maxsize)
         eigen_z_lwork = max(eigen_z_lwork,eigen_z_maxsize*(nb+1)) ! MAGMAF_ZHEEVD, MAGMAF_ZHEGVD
#endif
       end if
     end if
     if (ABI_LINALG_PLASMA_ISON) then
       if (lapack_full_storage.and.lapack_divide_conquer) then
         eigen_z_lwork = max(eigen_z_lwork,eigen_z_maxsize**2) ! PLASMA_ZHEEV
       end if
     end if
   end if
 end if
 ABI_SFREE(eigen_z_work)
 ABI_ALLOCATE(eigen_z_work,(eigen_z_lwork))

!Single precision RWORK
 eigen_c_lrwork = 0
 if (eigen_c_maxsize>0) then
   if (lapack_single_precision) then
     if (lapack_full_storage.or.lapack_packed_storage) then
       eigen_c_lrwork = max(eigen_c_lrwork,3*eigen_c_maxsize-2) ! CHEEV, CHEGV, CHPEV, CHPGV
     end if
     if (lapack_divide_conquer) then
       eigen_c_lrwork = max(eigen_c_lrwork,1+5*eigen_c_maxsize+2*eigen_c_maxsize**2) ! CHEEVD, CHEGVD, CHPEVD, CHPGVD
     end if
   end if
 end if
 ABI_SFREE(eigen_c_rwork)
 ABI_ALLOCATE(eigen_c_rwork,(eigen_c_lrwork))

!Double precision RWORK
 eigen_z_lrwork = 0
 if (eigen_z_maxsize>0) then
   if (lapack_double_precision) then
     if (lapack_full_storage.or.lapack_packed_storage) then
       eigen_z_lrwork = max(eigen_z_lrwork,3*eigen_z_maxsize-2) ! ZHEEV, ZHEGV, ZHPEV, ZHPGV
     end if
     if (lapack_divide_conquer) then
       eigen_z_lrwork = max(eigen_z_lrwork,1+5*eigen_z_maxsize+2*eigen_z_maxsize**2) ! ZHEEVD, ZHEGVD, ZHPEVD, ZHPGVD
     end if
   end if
 end if
 ABI_SFREE(eigen_z_rwork)
 ABI_ALLOCATE(eigen_z_rwork,(eigen_z_lrwork))

!Integer IWORK
 eigen_liwork = 0
 if (lapack_divide_conquer) then
   if (lapack_single_precision) then
     if (eigen_s_maxsize>0) eigen_liwork = max(eigen_liwork,3+5*eigen_s_maxsize)
     if (eigen_c_maxsize>0) eigen_liwork = max(eigen_liwork,3+5*eigen_c_maxsize)
   end if
   if (lapack_double_precision) then
     if (eigen_d_maxsize>0) eigen_liwork = max(eigen_liwork,3+5*eigen_d_maxsize)
     if (eigen_z_maxsize>0) eigen_liwork = max(eigen_liwork,3+5*eigen_z_maxsize)
   end if
 end if
 ABI_SFREE(eigen_iwork)
 ABI_ALLOCATE(eigen_iwork,(eigen_liwork))

 end subroutine abi_linalg_work_allocate
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
!!      m_driver,m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE
!!
 subroutine abi_linalg_finalize()

!Arguments ------------------------------------

!Local variables ------------------------------
#ifdef HAVE_LINALG_PLASMA
 integer :: info
#endif

!******************************************************************

 if (.not.abi_linalg_in_use) return

 eigen_s_maxsize = 0
 eigen_d_maxsize = 0
 eigen_c_maxsize = 0
 eigen_z_maxsize = 0
 eigen_s_lwork   = 0
 eigen_d_lwork   = 0
 eigen_c_lwork   = 0
 eigen_z_lwork   = 0
 eigen_c_lrwork  = 0
 eigen_z_lrwork  = 0
 eigen_liwork    = 0

 lapack_full_storage=.False.
 lapack_packed_storage=.False.
 lapack_single_precision=.False.
 lapack_double_precision=.False.
 lapack_divide_conquer=.false.

#ifdef HAVE_LINALG_SCALAPACK
 if (ABI_LINALG_SCALAPACK_ISON) then
   call end_scalapack(slk_processor)
   call xmpi_comm_free(slk_communicator)
   call xmpi_comm_free(slk_complement_communicator)
   slk_communicator=xmpi_comm_null
   slk_complement_communicator=xmpi_comm_null
   slk_minsize=1
 end if
#endif

#ifdef HAVE_LINALG_ELPA
 call elpa_func_uninit()
#endif

#ifdef HAVE_LINALG_PLASMA
   call PLASMA_Finalize(info)
   ABI_LINALG_PLASMA_ISON=.False.
 end if
#endif

#ifdef HAVE_LINALG_MAGMA
#ifdef HAVE_LINALG_MAGMA_15
 if (ABI_LINALG_PLASMA_ISON) then
   call magmaf_finalize()
 end if
#endif
#endif

#ifdef HAVE_GPU_CUDA
 call gpu_linalg_shutdown()
#endif

!Memory freeing
 ABI_SFREE(eigen_s_work)
 ABI_SFREE(eigen_d_work)
 ABI_SFREE(eigen_c_work)
 ABI_SFREE(eigen_z_work)
 ABI_SFREE(eigen_c_rwork)
 ABI_SFREE(eigen_z_rwork)
 ABI_SFREE(eigen_iwork)

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

subroutine linalg_allow_gemm3m(bool, write_msg)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: bool, write_msg

! *************************************************************************

 XGEMM3M_ISON = bool
 if (write_msg) then
#ifdef HAVE_LINALG_GEMM3M
   if (bool) then
     MSG_COMMENT("Activating ZGEMM3M version instead of ZGEMM")
   else
     MSG_COMMENT("Using ZGEMM instead of ZGEMM3M")
   end if
#else
   if (bool) then
     MSG_WARNING("Cannot activate ZGEMM3M as HAVE_LINALG_GEMM3M is not defined!")
   end if
#endif
 endif

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

pure logical function use_zgemm3m(m, n, k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m,n,k

! *************************************************************************

 use_zgemm3m = .False.
 if (XGEMM3M_ISON) use_zgemm3m = ((m * n * k) > ZGEMM3M_LIMIT)
 !if (XGEMM3M_ISON) use_zgemm3m = .True.

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

pure logical function use_cgemm3m(m, n, k)

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
