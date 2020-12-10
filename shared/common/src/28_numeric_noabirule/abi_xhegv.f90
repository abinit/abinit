!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xhegv
!! NAME
!!  abi_xhegv
!!
!! FUNCTION
!!  abi_xhegv is the generic function that compute
!!  all eigenvalues and, optionally, eigenvectors of a
!!  generalized symmetric-definite eigenproblem, of the form
!!  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
!!  Here A and B are assumed to be symmetric (or hermitian) and
!!  B is also positive definite.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (LNguyen,FDahm,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/abi_dhegv
!! NAME
!! abi_dhegv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!

  subroutine abi_dhegv(itype,jobz,uplo,n,a,lda,b,ldb,w, &
&            x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

!Arguments ------------------------------------
 integer, intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda,ldb
 real(dp), intent(inout) :: a(n,*) ! FIXME should be cplex * lda
 real(dp), intent(inout) :: b(n,*)
 real(dp),target,intent(out) :: w(n)
 integer, optional, intent(in) :: x_cplx
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: timopt,tim_xeigen
 integer, optional, intent(in) :: use_gpu,use_slk

!Local variables ------------------------------------
 integer :: cplx_,istwf_k_,usegpu_,use_slk_
 integer :: info
 real(dp) :: tsec(2)

! *********************************************************************

 ABI_CHECK(lapack_full_storage,"BUG(1) in abi_dhegv (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_dhegv (precision)!")
 ABI_CHECK(n<=eigen_d_maxsize,"BUG(3) in abi_dhegv (maxsize)!")

 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,1,tsec)
 end if

 cplx_=1 ; if(present(x_cplx)) cplx_ = x_cplx
 usegpu_=0; if (present(use_gpu)) usegpu_=use_gpu
 use_slk_=0; if (present(use_slk)) use_slk_=use_slk
 istwf_k_=1; if (present(istwf_k)) istwf_k_=istwf_k

!===== MAGMA
 if (ABI_LINALG_MAGMA_ISON.and.usegpu_==1) then
#if defined HAVE_LINALG_MAGMA
   ABI_CHECK((lapack_divide_conquer),"BUG(4) in abi_dhegv (d&c)!")
   if (cplx_ == 2) then
     call magmaf_zhegvd(itype,jobz,uplo,n,a,lda,b,ldb,w,eigen_z_work,eigen_z_lwork, &
&           eigen_z_rwork,eigen_z_lrwork,eigen_iwork,eigen_liwork,info)
   else
     call magmaf_dsygvd(jobz,uplo,n,a,lda,w,eigen_d_work,eigen_d_lwork,&
&                        eigen_iwork,eigen_liwork,info)
   endif
#endif

!===== SCALAPACK
 else if (ABI_LINALG_SCALAPACK_ISON.and.use_slk_==1.and.n>slk_minsize) then
#if defined HAVE_LINALG_SCALAPACK
   ABI_CHECK(present(x_cplx),"BUG(5) in abi_dhegv (x_cplx)!")
   call compute_eigen2(slk_communicator,slk_processor,cplx_,n,n,a,b,w,istwf_k_)
   info = 0 ! This is to avoid unwanted warning but it's not clean
#endif

!===== PLASMA
 !FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !  the eigenvalues computation (jobz=N) and not the
 ! full eigenvectors bases determination (jobz=V)
 else if (ABI_LINALG_PLASMA_ISON.and.LSAME(jobz,'N')) then
#if defined HAVE_LINALG_PLASMA
   if (cplx_ == 2) then
     call PLASMA_Alloc_Workspace_zhegv(n,n,plasma_work,info)
     info = PLASMA_zhegv_c(itype,jobz_plasma(jobz),uplo_plasma(uplo),n,c_loc(a),lda,&
&                          c_loc(b),ldb,c_loc(w),plasma_work,c_loc(eigen_z_work),n)
   else
     call PLASMA_Alloc_Workspace_dsygv(n,n,plasma_work,info)
     info = PLASMA_dsygv_c(itype,jobz_plasma(jobz),uplo_plasma(uplo),n,c_loc(a),lda,&
&                          c_loc(b),ldb,c_loc(w),plasma_work,c_loc(eigen_d_work),n)
   endif
   call PLASMA_Dealloc_handle(plasma_work,info)
#endif

!===== LAPACK
 else
   if (cplx_ == 2) then
     call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,eigen_z_work,eigen_z_lwork,eigen_z_rwork,info)
   else
     call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,eigen_d_work,eigen_d_lwork,info)
   endif
 end if

 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,2,tsec)
 end if

 ABI_CHECK(info==0,"abi_dhegv returned info!=0!")

end subroutine abi_dhegv
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_chegv
!! NAME
!! abi_chegv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
subroutine abi_chegv(itype,jobz,uplo,n,a,lda,b,ldb,w)

 implicit none

 !Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda,ldb
 complex(spc), intent(inout) :: a(lda,*)
 complex(spc), intent(inout) :: b(ldb,*)
 real(sp), intent(out) :: w(n)

!Local variables-------------------------------
 integer :: info,lwork
 real(sp),pointer :: rwork(:)
 complex(spc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_full_storage,"BUG(1) in abi_chegv (storage)!")
 ABI_CHECK(lapack_single_precision,"BUG(2) in abi_chegv (precision)!")
 ABI_CHECK(n<=eigen_c_maxsize,"BUG(3) in abi_chegv (maxsize)!")

 work => eigen_c_work ; rwork => eigen_c_rwork
 lwork=eigen_c_lwork

!===== PLASMA
 !FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !  the eigenvalues computation (jobz=N) and not the
 !  full eigenvectors bases determination (jobz=V)
 if (ABI_LINALG_PLASMA_ISON.and.LSAME(jobz,'N')) then
#if defined HAVE_LINALG_PLASMA
   if (eigen_c_lwork==0) then
     ABI_MALLOC(work,(n**2))
   end if
   call PLASMA_Alloc_Workspace_chegv(n,n,plasma_work,info)
   info = call PLASMA_chegv(itype,jobz_plasma(jobz),uplo_plasma(uplo),n,a,lda,b,ldb,w,&
&                           plasma_work,c_loc(work),n)
   call PLASMA_Dealloc_handle(plasma_work,info)
   if (eigen_c_lwork==0) then
     ABI_FREE(work)
   end if
#endif

!===== LAPACK
 else
   if (eigen_c_lwork==0) then
     lwork=2*n-1
     ABI_MALLOC(work,(lwork))
   end if
   if (eigen_c_lrwork==0) then
     ABI_MALLOC(rwork,(3*n-2))
   end if
   call chegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
   if (eigen_c_lwork==0) then
     ABI_FREE(work)
   end if
   if (eigen_c_lrwork==0) then
     ABI_FREE(rwork)
   end if
 end if

 ABI_CHECK(info==0,"abi_chegv returned info!=0!")

end subroutine abi_chegv
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zhegv
!! NAME
!! abi_zhegv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!!
!! SOURCE

subroutine abi_zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda,ldb
 complex(dpc), intent(inout) :: a(lda,*)
 complex(dpc), intent(inout) :: b(ldb,*)
 real(dp), intent(out) :: w(n)

!Local variables-------------------------------
 integer :: info,lwork
 real(dp),pointer :: rwork(:)
 complex(dpc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_full_storage,"BUG(1) in abi_zhegv (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_zhegv (precision)!")
 ABI_CHECK(n<=eigen_z_maxsize,"BUG(3) in abi_zhegv (maxsize)!")

 work => eigen_z_work ; rwork => eigen_z_rwork
 lwork=eigen_z_lwork

!===== PLASMA
 !FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !  the eigenvalues computation (jobz=N) and not the
 ! full eigenvectors bases determination (jobz=V)
 if (ABI_LINALG_PLASMA_ISON.and.LSAME(jobz,'N')) then
#if defined HAVE_LINALG_PLASMA
   if (eigen_z_lwork==0) then
     ABI_MALLOC(work,(n**2))
   end if
   call PLASMA_Alloc_Workspace_zhegv(n,n,plasma_work,info)
   call PLASMA_zhegv(itype,jobz_plasma(jobz),uplo_plasma(uplo),n,a,lda,b,ldb,w,&
&                    plasma_work,c_loc(work),n)
   call PLASMA_Dealloc_handle(plasma_work,info)
   if (eigen_z_lwork==0) then
     ABI_FREE(work)
   end if
#endif

!===== LAPACK
 else
   if (eigen_z_lwork==0) then
     lwork=2*n-1
     ABI_MALLOC(work,(lwork))
   end if
   if (eigen_z_lrwork==0) then
     ABI_MALLOC(rwork,(3*n-2))
   end if
   call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
   if (eigen_z_lwork==0) then
     ABI_FREE(work)
   end if
   if (eigen_z_lrwork==0) then
     ABI_FREE(rwork)
   end if
 end if

 ABI_CHECK(info==0,"abi_zhegv returned info!=0!")

end subroutine abi_zhegv
!!***
