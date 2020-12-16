!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xtrsm
!! NAME
!!  abi_xtrsm
!!
!! FUNCTION
!!  abi_xtrsm is the generic function that solve :
!! *     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!! *
!! *  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!! *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!! *
!! *     op( A ) = A   or   op( A ) = A**T.
!! *
!! *  The matrix X is overwritten on B.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (LNguyen,FDahm (CS))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/abi_ztrsm
!! NAME
!! abi_ztrsm
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

!Arguments-------------------------------------
 character(len=1), intent(in) :: side
 character(len=1), intent(in) :: uplo
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: diag
 integer, intent(in) :: m,n,ldb,lda
 complex(dpc), intent(in) :: alpha
 complex(dpc),target,intent(in) :: a(lda,*)
 complex(dpc),target,intent(inout) :: b(ldb,*)

!Local variables-------------------------------
#ifdef HAVE_LINALG_PLASMA
 integer :: info
#endif

#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XTRSM,1,tsec)
#endif

 if (ABI_LINALG_PLASMA_ISON) then
#ifdef HAVE_LINALG_PLASMA
   info = PLASMA_ztrsm_c(side_plasma(side),uplo_plasma(uplo),trans_plasma(transa),diag_plasma(diag),&
&     m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
#endif
 else
   call ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XTRSM,2,tsec)
#endif

end subroutine abi_ztrsm
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_dtrsm
!! NAME
!! abi_dtrsm
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb,&
&       x_cplx)

 !Arguments-------------------------------------
 character(len=1), intent(in) :: side,uplo,transa,diag
 integer, intent(in) :: m,n,lda,ldb
 real(dp), intent(in) :: alpha
 real(dp),target, intent(in) :: a(lda,*)       ! FIXME should be lda * x_cplx
 real(dp),target, intent(inout) :: b(ldb,*)
 !Only for lobpcgwf
 integer, intent(in), optional :: x_cplx

 !Local variables-------------------------------
 integer  :: cplx_
#ifdef HAVE_LINALG_PLASMA
 integer :: info
#endif

#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XTRSM,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 if (ABI_LINALG_PLASMA_ISON) then
#ifdef HAVE_LINALG_PLASMA
   if(cplx_ == 2) then
      info = PLASMA_ztrsm_c(side_plasma(side),uplo_plasma(uplo),trans_plasma(TRANSA),diag_plasma(diag),&
&       m,n,cmplx(alpha,0.d0,dpc),c_loc(a),lda,c_loc(b),ldb)
   else
      info = PLASMA_dtrsm_c(side_plasma(side),uplo_plasma(uplo),trans_plasma(TRANSA),diag_plasma(diag),&
&       m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
   end if
#endif
 else
   if(cplx_ == 2) then
      call ztrsm(side,uplo,transa,diag,m,n,cmplx(alpha,0.d0,dpc),a,lda,b,ldb)
   else
      call dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
   end if
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XTRSM,2,tsec)
#endif

end subroutine abi_dtrsm
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2ztrsm
!! NAME
!! abi_d2ztrsm
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

 subroutine abi_d2ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb,&
&  x_cplx)

!Arguments-------------------------------------
 character(len=1), intent(in) :: side,uplo,transa,diag
 integer, intent(in) :: m,n,lda,ldb
 complex(dpc), intent(in) :: alpha
 real(dp),target, intent(in) :: a(lda,*)           ! FIXME should be lda * x_cplx
 real(dp),target, intent(inout) :: b(ldb,*)
 !Only for lobpcgwf
 integer, intent(in), optional :: x_cplx

!Local variables-------------------------------
 integer  :: cplx_
#ifdef HAVE_LINALG_PLASMA
 integer :: info
#endif

#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XTRSM,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 if (ABI_LINALG_PLASMA_ISON) then
#ifdef HAVE_LINALG_PLASMA
   if(cplx_ == 2) then
      info = PLASMA_ztrsm_c(side_plasma(side),uplo_plasma(uplo),trans_plasma(TRANSA),diag_plasma(diag),&
&       m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
   else
      info = PLASMA_dtrsm_c(side_plasma(side),uplo_plasma(uplo),trans_plasma(TRANSA),diag_plasma(diag),&
&       m,n,real(alpha,dp),c_loc(a),lda,c_loc(b),ldb)
   end if
#endif
 else
   if(cplx_ == 2) then
      call ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
   else
      call dtrsm(side,uplo,transa,diag,m,n,real(alpha,dp),a,lda,b,ldb)
   end if
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XTRSM,2,tsec)
#endif

end subroutine abi_d2ztrsm
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2ztrsm_3d
!! NAME
!! abi_d2ztrsm_3d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!

  subroutine abi_d2ztrsm_3d(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

!Arguments-------------------------------------
 character(len=1), intent(in) :: side,uplo,transa,diag
 integer, intent(in) :: m,n,lda,ldb
 complex(dpc), intent(in) :: alpha
 real(dp), target,intent(in) :: a(2,lda,*)
 real(dp), target,intent(inout) :: b(2,ldb,*)

!Local variables-------------------------------
#ifdef HAVE_LINALG_PLASMA
 integer :: info
#endif

#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XTRSM,1,tsec)
#endif

 if (ABI_LINALG_PLASMA_ISON) then
#ifdef HAVE_LINALG_PLASMA
   info = PLASMA_ztrsm_c(side_plasma(side),uplo_plasma(uplo),trans_plasma(TRANSA),diag_plasma(diag),&
&    m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
#endif
 else
   call ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XTRSM,2,tsec)
#endif

end subroutine abi_d2ztrsm_3d
!!***

!----------------------------------------------------------------------
