!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xgemm
!! NAME
!!  abi_xgemm
!!
!! FUNCTION
!!  abi_xgemm is the generic function that solve :
!! *     C := alpha*op( A )*op( B ) + beta*C,
!! *
!! *  where  op( X ) is one of
!! *
!! *     op( X ) = X   or   op( X ) = X**T,
!! *
!! *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!! *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (LNguyen,FDahm (CS))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zgemm_2d
!! NAME
!! abi_zgemm_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!!
!! PARENTS
!!
!! SOURCE

 subroutine abi_zgemm_2d(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

 !Arguments------------------------------------
 character(len=1),intent(in) :: TRANSA
 character(len=1),intent(in) :: TRANSB
 integer,intent(in) :: K
 integer,intent(in) :: LDA
 integer,intent(in) :: LDB
 integer,intent(in) :: LDC
 integer,intent(in) :: M
 integer,intent(in) :: N
 complex(dpc),intent(in) :: ALPHA
 complex(dpc),intent(in) :: BETA
 complex(dpc),target,intent(in) :: A(lda,*)
 complex(dpc),target,intent(in) :: B(ldb,*)
 complex(dpc),target,intent(inout) :: C(ldc,*)

 integer :: info
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XGEMM,1,tsec)
#endif

 if (ABI_LINALG_PLASMA_ISON) then
   info = -1
#ifdef HAVE_LINALG_PLASMA
   !write(std_out,*)"Will call PLASMA_zgemm_c"
   info = PLASMA_zgemm_c(trans_plasma(TRANSA),trans_plasma(TRANSB),M,N,K,ALPHA,&
&    c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)
#endif
   ABI_CHECK(info==0,"PLASMA_zgemm_c returned info !=0")
 else
   if (use_zgemm3m(m,n,k)) then
     call _ZGEMM3M(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   else
     call zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   end if
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XGEMM,2,tsec)
#endif

end subroutine abi_zgemm_2d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm
!! NAME
!! abi_d2zgemm
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine abi_d2zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dpc),intent(in) :: alpha
 complex(dpc),intent(in) :: beta
 real(dp),target, intent(in)    :: a(lda,*)    ! FIXME should be lda * x_cplx
 real(dp),target, intent(in)    :: b(ldb,*)
 real(dp),target, intent(inout) :: c(ldc,*)
!Only for lobpcgwf
 integer, intent(in), optional :: x_cplx

!Local variables-------------------------------
 integer  :: cplx_,info
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XGEMM,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 if (ABI_LINALG_PLASMA_ISON) then
   info = -1
#ifdef HAVE_LINALG_PLASMA
   !write(std_out,*) "Will call plasma_[zd]gemm"
   if (cplx_ == 2) then
      info = PLASMA_zgemm_c(trans_plasma(TRANSA),trans_plasma(TRANSB),M,N,K,&
&      ALPHA,c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)
   else
      info = PLASMA_dgemm_c(trans_plasma(TRANSA),trans_plasma(TRANSB),M,N,K,&
&      real(ALPHA,dp),c_loc(A),LDA,c_loc(B),LDB,real(BETA,dp),c_loc(C),LDC)
   end if
#endif
   ABI_CHECK(info==0,"PLASMA_[z,d]gemm_c returned info !=0")
 else
   if (cplx_ == 2) then
     if (use_zgemm3m(m,n,k)) then
       call _ZGEMM3M(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     else
       call zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     end if
   else
     call dgemm(TRANSA,TRANSB,M,N,K,real(ALPHA,dp),A,LDA,B,LDB,real(BETA,dp),C,LDC)
   end if
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XGEMM,2,tsec)
#endif

end subroutine abi_d2zgemm
!!***

!----------------------------------------------------------------------

