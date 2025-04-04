!!****m* ABINIT/m_xg_ortho_RR
!! NAME
!!  m_xg_ortho_RR
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2024-2025 ABINIT group (J. Bieder, L. Baguet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

#include "nvtx_macros.h"

module m_xg_ortho_RR

  use m_errors
  use m_abicore
  use defs_basis
  use m_time, only : timab,abi_wtime
  use m_xmpi

  use m_xg
  use m_xg_nonlop
  use m_xgScalapack

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

#if defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

  implicit none

  private

  integer, parameter :: VAR_X   = 1000
  integer, parameter :: VAR_XW  = 1010
  integer, parameter :: VAR_XWP = 1100

  integer, parameter :: EIGENVX = 1
  integer, parameter :: EIGENVD = 2
  integer, parameter :: EIGENV = 3
  integer, parameter :: EIGENPVX = 4
  integer, parameter :: EIGENPVD = 5
  integer, parameter :: EIGENPV = 6
  integer, parameter :: EIGENEVD = 7
  integer, parameter :: EIGENEV = 8
  integer, parameter :: EIGENPEVD = 9
  integer, parameter :: EIGENPEV = 10
  integer, parameter :: EIGENSLK = 11
  logical, parameter :: EIGPACK(11) = &
    (/ .false.,.false.,.false., &
       .true. ,.true. ,.true. ,&
       .false.,.false.,&
       .true. ,.true., .false.  /)

  integer, parameter :: tim_RR_diago  = 1795
  integer, parameter :: tim_RR_gemm_1 = 1796
  integer, parameter :: tim_RR_gemm_2 = 1797

  public :: xg_Borthonormalize
  public :: xg_Borthonormalize_cprj
  public :: xg_RayleighRitz
  public :: xg_RayleighRitz_cprj

  contains
!!***

!!****f* m_xg_ortho_RR/xg_Borthonormalize
!!
!! NAME
!! xg_Borthonormalize
  subroutine xg_Borthonormalize(X,BX,info,timer,gpu_option,AX)

    integer        , intent(in   ) :: timer
    integer        , intent(in   ) :: gpu_option
    type(xgBlock_t), intent(inout) :: X
    type(xgBlock_t), intent(inout) :: BX
    type(xgBlock_t), intent(inout),optional :: AX
    integer       , intent(  out) :: info
    integer :: space_buf
    type(xg_t) :: buffer
    double precision :: tsec(2)

    ABI_NVTX_START_RANGE(NVTX_B_ORTHO)

    call timab(timer,1,tsec)

    call xgBlock_check(X,BX)
    if (present(AX)) then
      call xgBlock_check(X,AX)
    end if

    if (space(X)/=SPACE_CR) then
      space_buf = SPACE(X)
    else
      space_buf = SPACE_R
    end if
    call xg_init(buffer,space_buf,cols(X),cols(X),comm(X),gpu_option=gpu_option)

    ! If space(X)==SPACE_CR : set imaginary part of G=0 component to zero to improve numerical stability
    call xgBlock_zero_im_g0(X)
    call xgBlock_zero_im_g0(BX)
    if (present(AX)) then
      call xgBlock_zero_im_g0(AX)
    end if

    ! Compute X^TBX
    call xgBlock_gemm('t','n',1.d0,X,BX,0.d0,buffer%self,comm=comm(X))

    ! Compute Cholesky decomposition (Upper part)
    call xgBlock_potrf(buffer%self,'u',info)

    if ( info /= 0 ) then
      ABI_COMMENT("Cholesky decomposition did not work. Orthonormalization not done")
      call xg_free(buffer)
      ABI_NVTX_END_RANGE()
      return
    end if

    ! Solve YU=X
    call xgBlock_trsm('r','u',buffer%normal,'n',1.d0,buffer%self,X)

    ! Solve BYU=BX
    call xgBlock_trsm('r','u',buffer%normal,'n',1.d0,buffer%self,BX)

    if (present(AX)) then
      ! Solve AYU=AX
      call xgBlock_trsm('r','u',buffer%normal,'n',1.d0,buffer%self,AX)
    end if

    call xg_free(buffer)

    ABI_NVTX_END_RANGE()
    call timab(timer,2,tsec)

  end subroutine xg_Borthonormalize
!!***

!!****f* m_xg_ortho_RR/xg_Borthonormalize_cprj
!!
!! NAME
!! xg_Borthonormalize_cprj
subroutine xg_Borthonormalize_cprj(xg_nonlop,blockdim,X,cprjX,info,timer,gpu_option,AX)

    integer          , intent(in   ) :: timer,gpu_option
    integer          , intent(in   ) :: blockdim
    integer          , intent(  out) :: info
    type(xg_nonlop_t), intent(in   ) :: xg_nonlop
    type(xgBlock_t)  , intent(inout) :: X
    type(xgBlock_t)  , intent(inout) :: cprjX
    type(xgBlock_t)  , intent(inout),optional :: AX

    type(xg_t) :: buffer,cprj_work
    type(xgBlock_t) :: cprjX_spinor,cprj_work_spinor
    integer :: space_buf
    integer :: spacecom,ncols_cprj,nn,nspinor
    double precision :: tsec(2)

    call timab(timer,1,tsec)

    if (gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(X,cprjX)
    if (present(AX)) then
      call xgBlock_check_gpu_option(X,AX)
    end if

    nn = cols(X)
    ncols_cprj = cols(cprjX)
    nspinor = xg_nonlop%nspinor

    spacecom = comm(X)

    if (space(X)/=SPACE_CR) then
      space_buf = SPACE(X)
    else
      space_buf = SPACE_R
    end if
    call xg_init(buffer,space_buf,nn,nn,spacecom)

    ! Compute X^TX
    call xgBlock_gemm('t','n',1.d0,X,X,0.d0,buffer%self,comm=spacecom)
    call xg_init(cprj_work,space(cprjX),rows(cprjX),ncols_cprj,spacecom)
    if (xg_nonlop%paw) then
      call xg_nonlop_getXSX(xg_nonlop,cprjX,cprjX,cprj_work%self,buffer%self,blockdim)
    end if

    ! Compute Cholesky decomposition (Upper part)
    call xgBlock_potrf(buffer%self,'u',info)

    if ( info /= 0 ) then
      ABI_COMMENT("Cholesky decomposition did not work. Orthonormalization not done")
      call xg_free(buffer)
      return
    end if

    ! Solve YU=X
    call xgBlock_trsm('r','u','n','n',1.d0,buffer%self,X)

    if (present(AX)) then
      ! Solve AYU=AX
      call xgBlock_trsm('r','u','n','n',1.d0,buffer%self,AX)
    end if

    ! Solve (cprjY)U=(cprjX)
    call xgBlock_reshape_spinor(cprjX,cprjX_spinor,nspinor,COLS2ROWS)
    if (ncols_cprj==nspinor*nn) then
      call xgBlock_trsm('r','u','n','n',1.d0,buffer%self,cprjX_spinor)
    else
      call xgBlock_invert_tri('u','n',buffer%self)
      call xgBlock_zerotri(buffer%self,'u')
      call xgBlock_zero(cprj_work%self)
      call xgBlock_reshape_spinor(cprj_work%self,cprj_work_spinor,nspinor,COLS2ROWS)
      call xgBlock_gemm_mpi_cyclic_permutation(cprjX_spinor,buffer%self,cprj_work_spinor,&
        & xg_nonlop%me_band,blockdim/nspinor,comm=xg_nonlop%comm_band)
      call xgBlock_copy(cprj_work%self,cprjX)
    end if

    call xg_free(cprj_work)

    call xg_free(buffer)

    call timab(timer,2,tsec)

  end subroutine xg_Borthonormalize_cprj
!!***

!!****f* m_xg_ortho_RR/xg_RayleighRitz
!!
!! NAME
!! xg_RayleighRitz
  subroutine xg_RayleighRitz(X,AX,BX,eigenvalues,info,prtvol,timer,gpu_option,&
    & tolerance,XW,AW,BW,P,AP,BP,WP,AWP,BWP,XWP,solve_ax_bx)

    integer        , intent(in   ) :: timer
    integer        , intent(in   ) :: gpu_option
    integer        , intent(in   ) :: prtvol
    type(xgBlock_t), intent(inout) :: eigenvalues
    type(xgBlock_t), intent(inout) :: X
    type(xgBlock_t), intent(inout) :: AX
    type(xgBlock_t), intent(inout) :: BX
    integer        , intent(  out) :: info
    double precision, optional, intent(in) :: tolerance
    ! LOBPCG only :
    type(xgBlock_t), intent(inout),optional :: XW
    type(xgBlock_t), intent(inout),optional :: AW
    type(xgBlock_t), intent(inout),optional :: BW
    type(xgBlock_t), intent(inout),optional :: P
    type(xgBlock_t), intent(inout),optional :: AP
    type(xgBlock_t), intent(inout),optional :: BP
    type(xgBlock_t), intent(inout),optional :: WP
    type(xgBlock_t), intent(inout),optional :: AWP
    type(xgBlock_t), intent(inout),optional :: BWP
    type(xgBlock_t), intent(inout),optional :: XWP
    ! End LOBPCG only
    logical, intent(in),optional :: solve_ax_bx
    integer :: var
    integer :: spacedim
    integer :: blockdim
    integer :: space_buf
    integer :: subdim
    integer :: spacecom
    integer :: eigenSolver
    double precision :: abstol
#ifdef HAVE_LINALG_SCALAPACK
    logical :: use_slk
#endif
    type(xg_t) :: vec
    type(xg_t) :: subA
    type(xg_t) :: subB
    type(xgBlock_t) :: subsub
    type(xgBlock_t) :: Cwp
    type(xgScalapack_t) :: scalapack
    double precision :: tsec(2)
    logical :: solve_ax_bx_

    call timab(timer , 1, tsec)

    solve_ax_bx_ = .false.
    if (present(solve_ax_bx)) then
      solve_ax_bx_ = solve_ax_bx
    end if

    var = VAR_X
    call xgBlock_check(X,AX)
    call xgBlock_check(X,BX)
    if (.not.solve_ax_bx_) then
      eigenSolver = EIGENEVD
    else
      eigenSolver = EIGENVD
    end if
    spacedim = rows(X)
    blockdim = cols(X)
    space_buf = space(X)
    if (space(X) == space_CR) then
      space_buf = space_R
    end if
    spacecom = comm(X)
    subdim   = blockdim

    if (present(XW)) then
      if (solve_ax_bx_) then
        ABI_ERROR('solve_ax_bx is not implemented in that case')
      end if
      var = VAR_XW
      eigenSolver = EIGENVD
      call xgBlock_check(X,AW)
      call xgBlock_check(X,BW)
      call xgBlock_check(X,P)
      call xgBlock_check(X,AP)
      call xgBlock_check(X,BP)
      call xgBlock_check(X,XW,fact_col=2)
      call xgBlock_check(X,WP,fact_col=2)
      call xgBlock_check(X,AWP,fact_col=2)
      call xgBlock_check(X,BWP,fact_col=2)
      subdim = 2*blockdim
      if (present(XWP)) then
        var = VAR_XWP
        call xgBlock_check(X,XWP,fact_col=3)
        subdim = 3*blockdim
      end if
    end if

#ifdef HAVE_LINALG_SCALAPACK
    call xgScalapack_init(scalapack,spacecom,subdim,prtvol-2,(gpu_option/=ABI_GPU_DISABLED),use_slk)
    if ( use_slk) then
      eigenSolver = EIGENSLK
    end if
#endif

    ! Select diago algorithm

    abstol = 0d0 ; if ( present(tolerance) ) abstol = tolerance

    call xg_init(subA,space_buf,subdim,subdim,spacecom,gpu_option=gpu_option)
    if ( solve_ax_bx_ .or. var /= VAR_X ) then
      call xg_init(subB,space_buf,subdim,subdim,spacecom,gpu_option=gpu_option)
    end if

    if ( eigenSolver == EIGENVX .or. eigenSolver == EIGENPVX ) then
      call xg_init(vec,space_buf,subdim,blockdim,gpu_option=gpu_option)
    else if ( EIGPACK(eigenSolver) ) then
      call xg_init(vec,space_buf,subdim,subdim,gpu_option=gpu_option)
    else
      call xg_setBlock(subA,vec%self,subdim,blockdim)
    endif

    ! Compute subA and subB by part
    !--- begin
    ! |  E  |  XAW  | XAP |  |  I  |  XBW  | XBP |
    ! |  *  |  WAW  | WAP |  |  *  |   I   | WBP |
    ! |  *  |   *   | PAP |  |  *  |   *   |  I  |

    call timab(tim_RR_gemm_1,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_RR_GEMM_1)

    ! If space(X)==SPACE_CR : set imaginary part of G=0 component to zero to improve numerical stability
    if (var == VAR_X)   call xgBlock_zero_im_g0(X)
    if (var == VAR_XW)  call xgBlock_zero_im_g0(XW)
    if (var == VAR_XWP) call xgBlock_zero_im_g0(XWP)
    call xgBlock_zero_im_g0(AX)
    call xgBlock_zero_im_g0(BX)

    call xg_setBlock(subA,subsub,blockdim,blockdim)
    call xgBlock_gemm('t','n',1.0d0,X,AX,0.d0,subsub,comm=spacecom)

    if ( solve_ax_bx_ .or. var /= VAR_X ) then
      call xg_setBlock(subB,subsub,blockdim,blockdim)
      call xgBlock_gemm('t','n',1.0d0,X,BX,0.d0,subsub,comm=spacecom)
    endif

    if ( var == VAR_XW .or. var == VAR_XWP ) then

      ! If space(X)==SPACE_CR : set imaginary part of G=0 component to zero to improve numerical stability
      call xgBlock_zero_im_g0(AW)
      call xgBlock_zero_im_g0(BW)

      ! subA
      call xg_setBlock(subA,subsub,2*blockdim,blockdim,fcol=blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XW,AW,0.d0,subsub,comm=spacecom)

      ! subB
      call xg_setBlock(subB,subsub,2*blockdim,blockdim,fcol=blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XW,BW,0.d0,subsub,comm=spacecom)

    end if

    if ( var == VAR_XWP ) then

      ! If space(X)==SPACE_CR : set imaginary part of G=0 component to zero to improve numerical stability
      call xgBlock_zero_im_g0(AP)
      call xgBlock_zero_im_g0(BP)

      ! subA
      call xg_setBlock(subA,subsub,3*blockdim,blockdim,fcol=2*blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XWP,AP,0.d0,subsub,comm=spacecom)

      ! subB
      call xg_setBlock(subB,subsub,3*blockdim,blockdim,fcol=2*blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XWP,BP,0.d0,subsub,comm=spacecom)

    end if

    call timab(tim_RR_gemm_1,2,tsec)
    ABI_NVTX_END_RANGE()

    if ( EIGPACK(eigenSolver) ) then
      call xgBlock_pack(subA%self,subA%self,'u')
      if ( solve_ax_bx_ .or. var /= VAR_X ) then
        call xgBlock_pack(subB%self,subB%self,'u')
      end if
    end if

    call timab(tim_RR_diago,1,tsec)
    tsec(2) = abi_wtime()
    if ( .not.solve_ax_bx_ .and. var == VAR_X ) then
      ABI_NVTX_START_RANGE(NVTX_RR_HEEV)
    ! Solve Hermitian eigen problem
      select case (eigenSolver)
      case (EIGENEVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using heevd"
        call xgBlock_heevd('v','u',subA%self,eigenvalues,info) ! work with GPU
      case (EIGENEV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using heev"
        call xgBlock_heev('v','u',subA%self,eigenvalues,info)
      case (EIGENPEVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpevd"
        call xgBlock_hpevd('v','u',subA%self,eigenvalues,vec%self,info)
      case (EIGENPEV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpev"
        call xgBlock_hpev('v','u',subA%self,eigenvalues,vec%self,info)
      case (EIGENSLK)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using pheev"
        call xgScalapack_heev(scalapack,subA%self,eigenvalues,gpu_option=gpu_option) ! work with GPU
        info = 0 ! No error code returned for the moment
      case default
        ABI_ERROR("Error for Eigen Solver HEEV")
      end select
    else
      ABI_NVTX_START_RANGE(NVTX_RR_HEGV)
      ! Solve Hermitian general eigen problem only for first blockdim eigenvalues
      select case (eigenSolver)
      case (EIGENVX)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegvx"
        call xgBlock_hegvx(1,'v','i','u',subA%self,subB%self,0.d0,0.d0,1,blockdim,abstol,&
          eigenvalues,vec%self,info)
      case (EIGENVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegvd"
        call xgBlock_hegvd(1,'v','u',subA%self,subB%self,eigenvalues,info) ! work with GPU
      case (EIGENV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegv"
        call xgBlock_hegv(1,'v','u',subA%self,subB%self,eigenvalues,info)
      case (EIGENPVX)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgvx"
        call xgBlock_hpgvx(1,'v','i','u',subA%self,subB%self,0.d0,0.d0,1,blockdim,abstol,&
          eigenvalues,vec%self,info)
      case (EIGENPVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgvd"
        call xgBlock_hpgvd(1,'v','u',subA%self,subB%self,eigenvalues,vec%self,info)
      case (EIGENPV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgv"
        call xgBlock_hpgv(1,'v','u',subA%self,subB%self,eigenvalues,vec%self,info)
      case (EIGENSLK)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using phegv"
        call xgScalapack_hegv(scalapack,subA%self,subB%self,eigenvalues,gpu_option=gpu_option) ! work with GPU
        info = 0 ! No error code returned for the moment
      case default
        ABI_ERROR("Error for Eigen Solver HEGV")
      end select
    end if
    if ( eigenSolver == EIGENSLK ) then
      call xgScalapack_free(scalapack)
    end if
    tsec(2) = abi_wtime() - tsec(2)
    if ( prtvol == 4 ) write(std_out,*) tsec(2)

    call timab(tim_RR_diago,2,tsec)
    ABI_NVTX_END_RANGE()

    if ( eigenSolver == EIGENVX .or. EIGPACK(eigenSolver)) then
      call xg_free(subA)
    end if
    call xg_free(subB)

    call timab(tim_RR_gemm_2,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_RR_GEMM_2)

    !FIXME Avoid those transfers
    if ( info == 0 ) then
      call xg_init(subB,space(X),spacedim,blockdim,comm=comm(X),me_g0=me_g0(X),gpu_option=gpu_option)

      !/* Easy basic solution */
      !/* Compute first part of X here */
      ! Use subB as buffer
      !lobpcg%XWP (:,X+1:X+blockdim) = matmul(lobpcg%XWP (:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_setBlock(vec%self,Cwp,blockdim,blockdim)
      call xgBlock_gemm('n','n',1.0d0,X,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,X)

      !lobpcg%AXWP(:,X+1:X+blockdim) = matmul(lobpcg%AXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm('n','n',1.0d0,AX,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,AX)

      !lobpcg%BXWP(:,X+1:X+blockdim) = matmul(lobpcg%BXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm('n','n',1.0d0,BX,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,BX)

      if ( var /= VAR_X ) then
        ! Cost to pay to avoid temporary array in xgemm
        if(gpu_option==ABI_GPU_OPENMP) call xgBlock_copy_from_gpu(vec%self) !FIXME Avoid that transfer
        call xgBlock_cshift(vec%self,blockdim,1) ! Bottom 2*blockdim lines are now at the top
        if(gpu_option==ABI_GPU_OPENMP) call xgBlock_copy_to_gpu(vec%self) !FIXME Avoid that transfer
        call xgBlock_setBlock(vec%self,Cwp,subdim-blockdim,blockdim)

        !lobpcg%XWP (:,P+1:P+blockdim) = matmul(lobpcg%XWP (:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm('n','n',1.0d0,WP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,P)

        !lobpcg%AXWP(:,P+1:P+blockdim) = matmul(lobpcg%AXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm('n','n',1.0d0,AWP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,AP)

        !lobpcg%BXWP(:,P+1:P+blockdim) = matmul(lobpcg%BXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm('n','n',1.0d0,BWP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,BP)

        !/* Maybe faster solution
        ! * Sum previous contribution plus P direction
        ! */
        call xgBlock_add(X,P)
        call xgBlock_add(AX,AP)
        call xgBlock_add(BX,BP)
      end if
    end if

    call timab(tim_RR_gemm_2,2,tsec)

    ! Doing free on an already free object does not do anything
    call xg_free(vec)
    call xg_free(subA)
    call xg_free(subB)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    if (gpu_option==ABI_GPU_KOKKOS) then
      call gpu_device_synchronize()
    end if
#endif

    ABI_NVTX_END_RANGE()
    call timab(timer , 2, tsec)

  end subroutine xg_RayleighRitz
!!***

!!****f* m_xg_ortho_RR/xg_RayleighRitz_cprj
!!
!! NAME
!! xg_RayleighRitz_cprj
subroutine xg_RayleighRitz_cprj(xg_nonlop,X,cprjX,AX,eigenvalues,blockdim_cprj,info,prtvol,timer,gpu_option,&
     tolerance,XW,W,cprjXW,cprjW,AW,P,cprjP,AP,WP,cprjWP,AWP,XWP,cprjXWP,solve_ax_bx,add_Anl)

    integer          , intent(in   ) :: timer,gpu_option
    integer          , intent(in   ) :: prtvol
    integer          , intent(in   ) :: blockdim_cprj
    type(xgBlock_t)  , intent(inout) :: eigenvalues
    type(xgBlock_t)  , intent(inout) :: X
    type(xgBlock_t)  , intent(inout) :: cprjX
    type(xgBlock_t)  , intent(inout) :: AX
    type(xg_nonlop_t), intent(in   ) :: xg_nonlop
    integer        , intent(  out) :: info
    double precision, optional, intent(in) :: tolerance
    ! LOBPCG only :
    type(xgBlock_t), intent(inout),optional :: XW
    type(xgBlock_t), intent(inout),optional :: W
    type(xgBlock_t), intent(inout),optional :: cprjXW
    type(xgBlock_t), intent(inout),optional :: cprjW
    type(xgBlock_t), intent(inout),optional :: AW
    type(xgBlock_t), intent(inout),optional :: P
    type(xgBlock_t), intent(inout),optional :: cprjP
    type(xgBlock_t), intent(inout),optional :: AP
    type(xgBlock_t), intent(inout),optional :: WP
    type(xgBlock_t), intent(inout),optional :: cprjWP
    type(xgBlock_t), intent(inout),optional :: AWP
    type(xgBlock_t), intent(inout),optional :: XWP
    type(xgBlock_t), intent(inout),optional :: cprjXWP
    ! End LOBPCG only
    logical, intent(in),optional :: solve_ax_bx
    logical, intent(in),optional :: add_Anl
    integer :: var
    integer :: spacedim
    integer :: blockdim
    integer :: space_buf
    integer :: subdim
    integer :: spacecom
    integer :: eigenSolver
    integer :: nspinor,cprjdim,ncols_cprj
    !integer :: neigen
    double precision :: abstol
#ifdef HAVE_LINALG_SCALAPACK
    logical :: use_slk
#endif
    type(xg_t) :: vec
    type(xg_t) :: subA
    type(xg_t) :: subB
    type(xg_t) :: Xwork
    type(xg_t) :: cprjXwork
    type(xg_t) :: cprj_work
    type(xgBlock_t) :: cprj_workX,cprj_workX_spinor,cprjX_spinor,cprjW_spinor,cprjWP_spinor
    type(xgBlock_t) :: subsub
    type(xgBlock_t) :: Cwp
    type(xgScalapack_t) :: scalapack
    double precision :: tsec(2)
    logical :: solve_ax_bx_,add_Anl_

    call timab(timer , 1, tsec)

    if (gpu_option/=ABI_GPU_DISABLED) then
      ABI_ERROR('Not implemented for GPU')
    end if
    call xgBlock_check_gpu_option(X,cprjX)
    call xgBlock_check_gpu_option(X,AX)

    solve_ax_bx_ = .false.
    if (present(solve_ax_bx)) then
      solve_ax_bx_ = solve_ax_bx
    end if
    add_Anl_ = .false.
    if (present(add_Anl)) then
      add_Anl_ = add_Anl
    end if

    call timab(timer , 1, tsec)

    var = VAR_X
    call xgBlock_check(X,AX)
    if (.not.solve_ax_bx_) then
      eigenSolver = EIGENEVD
    else
      eigenSolver = EIGENVD
    end if
    spacedim  = rows(X)
    blockdim  = cols(X)
    space_buf = space(X)
    if (space(X)==SPACE_CR) then
      space_buf = SPACE_R
    end if
    spacecom = comm(X)
    subdim   = blockdim
    nspinor  = xg_nonlop%nspinor
    cprjdim  = xg_nonlop%cprjdim

    if (present(XW)) then
      if (solve_ax_bx_) then
        ABI_ERROR('solve_ax_bx is not implemented in that case')
      end if
      var = VAR_XW
      eigenSolver = EIGENVD
      !TODO Do checks on other optional arguments
      call xgBlock_check(X,AW)
      call xgBlock_check_gpu_option(X,AW)
      call xgBlock_check(X,P)
      call xgBlock_check_gpu_option(X,P)
      call xgBlock_check(X,AP)
      call xgBlock_check_gpu_option(X,AP)
      call xgBlock_check(X,XW,fact_col=2)
      call xgBlock_check_gpu_option(X,XW)
      call xgBlock_check(X,WP,fact_col=2)
      call xgBlock_check_gpu_option(X,WP)
      call xgBlock_check(X,AWP,fact_col=2)
      call xgBlock_check_gpu_option(X,AWP)
      subdim = 2*blockdim
      if (present(XWP)) then
        var = VAR_XWP
        call xgBlock_check(X,XWP,fact_col=3)
        subdim = 3*blockdim
      end if
    end if

#ifdef HAVE_LINALG_SCALAPACK
    call xgScalapack_init(scalapack,spacecom,subdim,prtvol-2,(gpu_option/=ABI_GPU_DISABLED),use_slk)
    if (use_slk) then
      eigenSolver = EIGENSLK
    end if
#endif

    abstol = 0d0 ; if ( present(tolerance) ) abstol = tolerance

    call xg_init(subA,space_buf,subdim,subdim,spacecom)
    if ( solve_ax_bx_ .or. var /= VAR_X ) then
      call xg_init(subB,space_buf,subdim,subdim,spacecom)
    end if

    if ( eigenSolver == EIGENVX .or. eigenSolver == EIGENPVX ) then
      call xg_init(vec,space_buf,subdim,blockdim)
    else if ( EIGPACK(eigenSolver) ) then
      call xg_init(vec,space_buf,subdim,subdim)
    else
      call xg_setBlock(subA,vec%self,subdim,blockdim)
    endif

    ! Compute subA and subB by part
    !--- begin
    ! |  E  |  XAW  | XAP |  |  I  |  XBW  | XBP |
    ! |  *  |  WAW  | WAP |  |  *  |   I   | WBP |
    ! |  *  |   *   | PAP |  |  *  |   *   |  I  |

    call timab(tim_RR_gemm_1,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_RR_GEMM_1)

    ! Compute XAX
    call xg_setBlock(subA,subsub,blockdim,blockdim)

    call xgBlock_gemm('t','n',1.0d0,X,AX,0.d0,subsub,comm=spacecom)

    ! Add the nonlocal part (H)
    call xg_init(cprj_work,space(cprjX),rows(cprjX),cols(cprjX),spacecom)
    if (add_Anl_) then
      call xg_nonlop_getXHX(xg_nonlop,cprjX,cprjX,cprj_work%self,subsub,blockdim_cprj)
    end if

    if ( solve_ax_bx_ .or. var /= VAR_X ) then
      ! Compute XBX
      call xg_setBlock(subB,subsub,blockdim,blockdim)
      call xgBlock_gemm('t','n',1.0d0,X,X,0.d0,subsub,comm=spacecom)
      ! Add the nonlocal part (S)
      if (xg_nonlop%paw) then
        call xg_nonlop_getXSX(xg_nonlop,cprjX,cprjX,cprj_work%self,subsub,blockdim_cprj)
      end if
    end if

    if ( var == VAR_XW .or. var == VAR_XWP ) then

      ! Compute XAW and WAW
      call xg_setBlock(subA,subsub,2*blockdim,blockdim,fcol=blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XW,AW,0.d0,subsub,comm=spacecom)

      ! Add the nonlocal part (H)
      if (add_Anl_) then
        call xg_nonlop_getXHX(xg_nonlop,cprjXW,cprjW,cprj_work%self,subsub,blockdim_cprj)
      end if

      ! Compute XBW and WBW
      call xg_setBlock(subB,subsub,2*blockdim,blockdim,fcol=blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XW,W,0.d0,subsub,comm=spacecom)
      ! Add the nonlocal part (S)
      if (xg_nonlop%paw) then
        call xg_nonlop_getXSX(xg_nonlop,cprjXW,cprjW,cprj_work%self,subsub,blockdim_cprj)
      end if

    end if

    if ( var == VAR_XWP ) then
      ! Compute XAP, WAP and PAP
      call xg_setBlock(subA,subsub,3*blockdim,blockdim,fcol=2*blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XWP,AP,0.d0,subsub,comm=spacecom)

      ! Add the nonlocal part (H)
      if (add_Anl_) then
        call xg_nonlop_getXHX(xg_nonlop,cprjXWP,cprjP,cprj_work%self,subsub,blockdim_cprj)
      end if

      call xg_setBlock(subB,subsub,3*blockdim,blockdim,fcol=2*blockdim+1)
      call xgBlock_gemm('t','n',1.0d0,XWP,P,0.d0,subsub,comm=spacecom)
      ! Add the nonlocal part (S)
      if (xg_nonlop%paw) then
        call xg_nonlop_getXSX(xg_nonlop,cprjXWP,cprjP,cprj_work%self,subsub,blockdim_cprj)
      end if
    end if

    call xg_free(cprj_work)

    call timab(tim_RR_gemm_1,2,tsec)
    ABI_NVTX_END_RANGE()

    if ( EIGPACK(eigenSolver) ) then
      call xgBlock_pack(subA%self,subA%self,'u')
    end if

    call timab(tim_RR_diago,1,tsec)
    tsec(2) = abi_wtime()
    if ( .not.solve_ax_bx_ .and. var == VAR_X ) then
    ! Solve Hermitian eigen problem
      select case (eigenSolver)
      case (EIGENEVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using heevd"
        call xgBlock_heevd('v','u',subA%self,eigenvalues,info)
      case (EIGENEV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using heev"
        call xgBlock_heev('v','u',subA%self,eigenvalues,info)
      case (EIGENPEVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpevd"
        call xgBlock_hpevd('v','u',subA%self,eigenvalues,vec%self,info)
      case (EIGENPEV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpev"
        call xgBlock_hpev('v','u',subA%self,eigenvalues,vec%self,info)
      case (EIGENSLK)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using pheev"
        call xgScalapack_heev(scalapack,subA%self,eigenvalues)
        info = 0 ! No error code returned for the moment
      case default
        ABI_ERROR("Error for Eigen Solver HEEV")
      end select
    else
      ! Solve Hermitian general eigen problem only for first blockdim eigenvalues
      select case (eigenSolver)
      case (EIGENVX)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegvx"
        call xgBlock_hegvx(1,'v','i','u',subA%self,subB%self,0.d0,0.d0,1,blockdim,abstol,&
          eigenvalues,vec%self,info)
      case (EIGENVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegvd"
        call xgBlock_hegvd(1,'v','u',subA%self,subB%self,eigenvalues,info)
      case (EIGENV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegv"
        call xgBlock_hegv(1,'v','u',subA%self,subB%self,eigenvalues,info)
      case (EIGENPVX)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgvx"
        call xgBlock_hpgvx(1,'v','i','u',subA%self,subB%self,0.d0,0.d0,1,blockdim,abstol,&
          eigenvalues,vec%self,info)
      case (EIGENPVD)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgvd"
        call xgBlock_hpgvd(1,'v','u',subA%self,subB%self,eigenvalues,vec%self,info)
      case (EIGENPV)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgv"
        call xgBlock_hpgv(1,'v','u',subA%self,subB%self,eigenvalues,vec%self,info)
      case (EIGENSLK)
        if ( prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using phegv"
        call xgScalapack_hegv(scalapack,subA%self,subB%self,eigenvalues)
        info = 0 ! No error code returned for the moment
      case default
        ABI_ERROR("Error for Eigen Solver HEGV")
      end select
    end if
    if ( eigenSolver == EIGENSLK ) then
      call xgScalapack_free(scalapack)
    end if
    tsec(2) = abi_wtime() - tsec(2)
    if ( prtvol == 4 ) write(std_out,*) tsec(2)
    call timab(tim_RR_diago,2,tsec)

    if ( eigenSolver == EIGENVX .or. EIGPACK(eigenSolver)) then
      call xg_free(subA)
    end if
    call xg_free(subB)

    call timab(tim_RR_gemm_2,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_RR_GEMM_2)

    if ( info == 0 ) then
      call xg_init(Xwork,space(X),spacedim,blockdim,me_g0=me_g0(X))
      ncols_cprj = cols(cprjX)
      call xg_init(cprjXwork,space(cprjX),cprjdim,ncols_cprj,comm(cprjX))
      call xg_setBlock(cprjXwork,cprj_workX,cprjdim,ncols_cprj)
      call xgBlock_reshape_spinor(cprj_workX,cprj_workX_spinor,nspinor,COLS2ROWS)

      call xgBlock_setBlock(vec%self,Cwp,blockdim,blockdim)

      !/* Easy basic solution */
      !/* Compute first part of X here */
      !XWP (:,X+1:X+blockdim) = matmul(XWP (:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm('n','n',1.0d0,X,Cwp,0.d0,Xwork%self)
      call xgBlock_copy(Xwork%self,X)

      !AXWP(:,X+1:X+blockdim) = matmul(AXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm('n','n',1.0d0,AX,Cwp,0.d0,Xwork%self)
      call xgBlock_copy(Xwork%self,AX)

      !cprjXWP(:,X+1:X+blockdim) = matmul(cprjXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_reshape_spinor(cprjX,cprjX_spinor,nspinor,COLS2ROWS)
      if (ncols_cprj==nspinor*cols(X)) then
        call xgBlock_gemm('n','n',1.0d0,cprjX_spinor,Cwp,0.d0,cprj_workX_spinor)
      else
        call xgBlock_zero(cprj_workX)
        call xgBlock_gemm_mpi_cyclic_permutation(cprjX_spinor,Cwp,cprj_workX_spinor,&
          & xg_nonlop%me_band,blockdim_cprj/nspinor,comm=xg_nonlop%comm_band)
      end if
      call xgBlock_copy(cprj_workX,cprjX)

      if ( var /= VAR_X ) then
        ! Cost to pay to avoid temporary array in xgemm
        call xgBlock_cshift(vec%self,blockdim,1) ! Bottom 2*blockdim lines are now at the top
        call xgBlock_setBlock(vec%self,Cwp,subdim-blockdim,blockdim)

        !XWP (:,P+1:P+blockdim) = matmul(XWP (:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm('n','n',1.0d0,WP,Cwp,0.d0,Xwork%self)
        call xgBlock_copy(Xwork%self,P)

        !AXWP(:,P+1:P+blockdim) = matmul(AXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm('n','n',1.0d0,AWP,Cwp,0.d0,Xwork%self)
        call xgBlock_copy(Xwork%self,AP)

        !cprjXWP(:,P+1:P+blockdim) = matmul(cprjXWP(:,W+1:X+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_reshape_spinor(cprjWP,cprjWP_spinor,nspinor,COLS2ROWS)
        if (ncols_cprj==cols(WP)) then
          call xgBlock_gemm('n','n',1.0d0,cprjWP_spinor,Cwp,0.d0,cprj_workX_spinor)
        else
          call xgBlock_zero(cprj_workX)
          if ( var==VAR_XW ) then
            call xgBlock_reshape_spinor(cprjW,cprjW_spinor,nspinor,COLS2ROWS)
            call xgBlock_gemm_mpi_cyclic_permutation(cprjW_spinor,Cwp,cprj_workX_spinor,&
              & xg_nonlop%me_band,blockdim_cprj/nspinor,comm=xg_nonlop%comm_band)
          else if ( var==VAR_XWP ) then
            call xgBlock_gemm_mpi_cyclic_permutation(cprjWP_spinor,Cwp,cprj_workX_spinor,&
              & xg_nonlop%me_band,blockdim_cprj/nspinor,comm=xg_nonlop%comm_band)
          else
            ABI_ERROR('not implemented')
          end if
        end if
        call xgBlock_copy(cprj_workX,cprjP)

        !/* Maybe faster solution
        ! * Sum previous contribution plus P direction
        ! */
        call xgBlock_add(X,P)
        call xgBlock_add(AX,AP)
        call xgBlock_add(cprjX,cprjP)
      end if

    end if

    call timab(tim_RR_gemm_2,2,tsec)
    ABI_NVTX_END_RANGE()

    ! Doing free on an already free object does not do anything
    call xg_free(vec)
    call xg_free(subA)
    call xg_free(Xwork)
    call xg_free(cprjXwork)

    call timab(timer , 2, tsec)

  end subroutine xg_RayleighRitz_cprj
!!***

end module m_xg_ortho_RR
!!***
