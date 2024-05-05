!!****m* ABINIT/m_xg_ortho_RR
!! NAME
!!  m_xg_ortho_RR
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2024-2024 ABINIT group (L. Baguet)
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
  use m_time, only : timab
  use m_xmpi

  use m_xg
  use m_xgScalapack

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
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
  public :: xg_RayleighRitz

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
    type(xg_t) :: buffer
    double precision :: tsec(2)

    ABI_NVTX_START_RANGE(NVTX_B_ORTHO)

    call timab(timer,1,tsec)

    call xgBlock_check(X,BX)
    if (present(AX)) then
      call xgBlock_check(X,AX)
    end if

    call xg_init(buffer,space(X),cols(X),cols(X),comm(X),gpu_option=gpu_option)

    ! Compute X^TBX
    call xgBlock_gemm(X%trans,BX%normal,1.d0,X,BX,0.d0,buffer%self)

    ! Compute Cholesky decomposition (Upper part)
    call xgBlock_potrf(buffer%self,'u',info)

    if ( info /= 0 ) then
      ABI_COMMENT("Cholesky decomposition did not work. Orthonormalization not done")
      call xg_free(buffer)
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

!!****f* m_xg_ortho_RR/xg_RayleighRitz
!!
!! NAME
!! xg_RayleighRitz
  subroutine xg_RayleighRitz(X,AX,BX,eigenvalues,info,prtvol,timer,gpu_option,&
    & tolerance,XW,AW,BW,P,AP,BP,WP,AWP,BWP,XWP,solve_ax_bx,istwf_k,usepaw,me_g0)

    use m_time
    integer        , intent(in   ) :: timer
    integer        , intent(in   ) :: gpu_option
    integer        , intent(in   ) :: prtvol
    type(xgBlock_t), intent(inout) :: eigenvalues
    type(xgBlock_t), intent(inout) :: X
    type(xgBlock_t), intent(inout) :: AX
    type(xgBlock_t), intent(inout) :: BX
    integer        , intent(  out) :: info
    double precision, optional, intent(in) :: tolerance
    ! CHEBFI only
    integer        , optional, intent(in)  :: istwf_k,me_g0
    logical        , optional, intent(in)  :: usepaw
    ! End CHEBFI only
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
    integer :: Xspace
    integer :: subdim
    integer :: spacecom
    integer :: eigenSolver
    integer :: istwf_k_,me_g0_
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
    logical :: usepaw_
    type(xgBlock_t) :: X_first_row
    type(xgBlock_t) :: AX_first_row
    type(xgBlock_t) :: BX_first_row


    call timab(timer , 1, tsec)

    solve_ax_bx_ = .false.
    if (present(solve_ax_bx)) then
      solve_ax_bx_ = solve_ax_bx
    end if

    istwf_k_=1; if(present(istwf_k)) istwf_k_=istwf_k
    usepaw_=.false.; if(present(usepaw)) usepaw_=usepaw
    me_g0_=1; if(present(me_g0)) me_g0_=me_g0

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
    Xspace   = space(X)
    spacecom = comm(X)
    subdim   = blockdim

    if (present(XW)) then
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

    call xg_init(subA,Xspace,subdim,subdim,spacecom,gpu_option=gpu_option)
    if ( solve_ax_bx_ .or. var /= VAR_X ) then
      call xg_init(subB,Xspace,subdim,subdim,spacecom,gpu_option=gpu_option)
    end if

    if ( eigenSolver == EIGENVX .or. eigenSolver == EIGENPVX ) then
      call xg_init(vec,Xspace,subdim,blockdim,gpu_option=gpu_option)
    else if ( EIGPACK(eigenSolver) ) then
      call xg_init(vec,Xspace,subdim,subdim,gpu_option=gpu_option)
    else
      call xg_setBlock(subA,vec%self,1,subdim,blockdim)
    endif

     ! Compute subA and subB by part
    !--- begin
    ! |  E  |  XAW  | XAP |  |  I  |  XBW  | XBP |
    ! |  *  |  WAW  | WAP |  |  *  |   I   | WBP |
    ! |  *  |   *   | PAP |  |  *  |   *   |  I  |

    call timab(tim_RR_gemm_1,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_RR_GEMM_1)

    call xg_setBlock(subA,subsub,1,blockdim,blockdim)
    call xgBlock_gemm(X%trans,AX%normal,1.0d0,X,AX,0.d0,subsub)

    if ( solve_ax_bx_ .or. var /= VAR_X ) then
      call xg_setBlock(subB,subsub,1,blockdim,blockdim)
      call xgBlock_gemm(X%trans,BX%normal,1.0d0,X,BX,0.d0,subsub)
    endif

    if ( var == VAR_XW .or. var == VAR_XWP ) then
      ! subA
      call xg_setBlock(subA,subsub,blockdim+1,2*blockdim,blockdim)
      call xgBlock_gemm(XW%trans,AW%normal,1.0d0,XW,AW,0.d0,subsub)

      ! subB
      call xg_setBlock(subB,subsub,blockdim+1,2*blockdim,blockdim)
      call xgBlock_gemm(XW%trans,BW%normal,1.0d0,XW,BW,0.d0,subsub)

    end if

    if ( var == VAR_XWP ) then
      ! subA
      call xg_setBlock(subA,subsub,2*blockdim+1,3*blockdim,blockdim)
      call xgBlock_gemm(XWP%trans,AP%normal,1.0d0,XWP,AP,0.d0,subsub)

      ! subB
      call xg_setBlock(subB,subsub,2*blockdim+1,3*blockdim,blockdim)
      call xgBlock_gemm(XWP%trans,BP%normal,1.0d0,XWP,BP,0.d0,subsub)

    end if

    call timab(tim_RR_gemm_1,2,tsec)
    if(istwf_k_ == 2) then
      call xgBlock_scale(X, 1/sqrt2, 1)
      if (me_g0_ == 1)  then
        call xgBlock_setBlock(X, X_first_row, 1, 2, subdim) !has to be 2 rows in SPACE_CR
        call xgBlock_scale(X_first_row, sqrt2, 1)
      end if
      call xgBlock_scale(AX, 1/sqrt2, 1)
      if (me_g0_ == 1)  then
        call xgBlock_setBlock(AX, AX_first_row, 1, 2, subdim)
        call xgBlock_scale(AX_first_row, sqrt2, 1)
      end if
      if (usepaw_)  then
        call xgBlock_scale(BX, 1/sqrt2, 1)
        if (me_g0_ == 1)  then
          call xgBlock_setBlock(BX, BX_first_row, 1, 2, subdim)
          call xgBlock_scale(BX_first_row, sqrt2, 1)
        end if
      end if
    end if

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
        call xgScalapack_heev(scalapack,subA%self,eigenvalues) ! work with GPU
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
        call xgScalapack_hegv(scalapack,subA%self,subB%self,eigenvalues) ! work with GPU
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
      call xg_init(subB,Xspace,spacedim,blockdim,comm=comm(X),gpu_option=gpu_option)

      !/* Easy basic solution */
      !/* Compute first part of X here */
      ! Use subB as buffer
      !lobpcg%XWP (:,X+1:X+blockdim) = matmul(lobpcg%XWP (:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_setBlock(vec%self,Cwp,1,blockdim,blockdim)
      call xgBlock_gemm(X%normal,Cwp%normal,1.0d0,X,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,X)

      !lobpcg%AXWP(:,X+1:X+blockdim) = matmul(lobpcg%AXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm(AX%normal,Cwp%normal,1.0d0,AX,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,AX)

      !lobpcg%BXWP(:,X+1:X+blockdim) = matmul(lobpcg%BXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm(BX%normal,Cwp%normal,1.0d0,BX,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,BX)

      if ( var /= VAR_X ) then
        ! Cost to pay to avoid temporary array in xgemm
        if(gpu_option==ABI_GPU_OPENMP) call xgBlock_copy_from_gpu(vec%self) !FIXME Avoid that transfer
        call xgBlock_cshift(vec%self,blockdim,1) ! Bottom 2*blockdim lines are now at the top
        if(gpu_option==ABI_GPU_OPENMP) call xgBlock_copy_to_gpu(vec%self) !FIXME Avoid that transfer
        call xgBlock_setBlock(vec%self,Cwp,1,subdim-blockdim,blockdim)

        !lobpcg%XWP (:,P+1:P+blockdim) = matmul(lobpcg%XWP (:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm(WP%normal,Cwp%normal,1.0d0,WP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,P)

        !lobpcg%AXWP(:,P+1:P+blockdim) = matmul(lobpcg%AXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm(AWP%normal,Cwp%normal,1.0d0,AWP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,AP)

        !lobpcg%BXWP(:,P+1:P+blockdim) = matmul(lobpcg%BXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm(BWP%normal,Cwp%normal,1.0d0,BWP,Cwp,0.d0,subB%self)
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

    ! Doing free on an already free object does not doe anything
    call xg_free(vec)
    call xg_free(subA)
    call xg_free(subB)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    if (chebfi%gpu_option==ABI_GPU_KOKKOS) then
      call gpu_device_synchronize()
    end if
#endif

    ABI_NVTX_END_RANGE()
    call timab(timer , 2, tsec)

  end subroutine xg_RayleighRitz
!!***

end module m_xg_ortho_RR
!!***
