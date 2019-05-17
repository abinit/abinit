!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xheev
!! NAME
!!  abi_xheev
!!
!! FUNCTION
!!  abi_xheev is the generic function that compute
!!  all eigenvalues and, optionally, eigenvectors of a
!!  symmetric or hermitian matrix A.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2019 ABINIT group (LNguyen,FDahm,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/abi_dheev
!! NAME
!! abi_dheev
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
  subroutine abi_dheev(jobz,uplo,n,a,lda,w,&
&            x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dheev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda
 real(dp), intent(inout) :: a(n,*)  ! FIXME should be cplex * lda
 real(dp), intent(out) :: w(n)
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: x_cplx
 integer, optional, intent(in) :: timopt,tim_xeigen
 integer, optional, intent(in) :: use_gpu,use_slk

!Local variables-------------------------------
 integer :: cplx_,istwf_k_,usegpu_,use_slk_
 integer :: info
 real(dp) :: tsec(2)

! *********************************************************************

 ABI_CHECK(lapack_full_storage,"BUG(1) in abi_dheev (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_dheev (precision)!")
 ABI_CHECK(n<=eigen_d_maxsize,"BUG(3) in abi_dheev (maxsize)!")

 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,1,tsec)
 end if

 cplx_=1 ; if(present(x_cplx)) cplx_ = x_cplx
 usegpu_=0;if (present(use_gpu)) usegpu_=use_gpu
 istwf_k_=1;if (present(istwf_k)) istwf_k_=istwf_k
 use_slk_ = 0; if(present(use_slk)) use_slk_ = use_slk

!===== MAGMA
 if (ABI_LINALG_MAGMA_ISON.and.usegpu_==1) then
#if defined HAVE_LINALG_MAGMA
   ABI_CHECK((lapack_divide_conquer),"BUG(4) in abi_dheev (d&c)!")
   if (cplx_ == 2) then
     call magmaf_zheevd(jobz,uplo,n,a,lda,w,eigen_z_work,eigen_z_lwork, &
&           eigen_z_rwork,eigen_z_lrwork,eigen_iwork,eigen_liwork,info)
   else
     call magmaf_dsyevd(jobz,uplo,n,a,lda,w,eigen_d_work,eigen_d_lwork,&
&                       eigen_iwork,eigen_liwork,info)
   endif
#endif

!===== SCALAPACK
 else if (ABI_LINALG_SCALAPACK_ISON.and.use_slk_==1.and.n>slk_minsize) then
#if defined HAVE_LINALG_SCALAPACK
   ABI_CHECK(present(x_cplx),"BUG(5) in abi_dheev (x_cplx)!")
   call compute_eigen1(slk_communicator,slk_processor,cplx_,n,n,a,w,istwf_k_)
   info = 0 ! This is to avoid unwanted warning but it's not clean
#endif

!===== PLASMA
 !FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !  the eigenvalues computation (jobz=N) and not the
 !  full eigenvectors bases determination (jobz=V)
 else if (ABI_LINALG_PLASMA_ISON.and.LSAME(jobz,'N')) then
#if defined HAVE_LINALG_PLASMA
   jobz_plasma_a = jobz_plasma(jobz)
   if (cplx_ == 2) then
     call PLASMA_Alloc_Workspace_zheev(n,n,plasma_work,info)
     info = PLASMA_zheev_c(jobz_plasma(jobz),uplo_plasma(uplo),n,c_loc(a),lda,c_loc(w),&
&                          plasma_work,c_loc(eigen_z_work),n)
   else
     call PLASMA_Alloc_Workspace_dsyev(n,n,plasma_work,info)
     info = PLASMA_dsyev_c(jobz_plasma(jobz),uplo_plasma(uplo),n,c_loc(a),lda,c_loc(w),&
&                          plasma_work,c_loc(eigen_d_work),n)
   endif
   call PLASMA_Dealloc_handle(plasma_work,info)
#endif

!===== LAPACK
 else
   if (cplx_ == 2) then
     call zheev(jobz,uplo,n,a,lda,w,eigen_z_work,eigen_z_lwork,eigen_z_rwork,info)
   else
     call dsyev(jobz,uplo,n,a,lda,w,eigen_d_work,eigen_d_lwork,info)
   endif
 end if

 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,2,tsec)
 end if

 ABI_CHECK(info==0,"abi_dheev returned info!=0!")

end subroutine abi_dheev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_cheev
!! NAME
!! abi_cheev
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_cheev(jobz,uplo,n,a,lda,w)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_cheev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda
 complex(spc), intent(inout) :: a(lda,*)
 real(sp), intent(out) :: w(n)

!Local variables-------------------------------
 integer :: info,lwork
 real(sp),pointer :: rwork(:)
 complex(spc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_full_storage,"BUG(1) in abi_cheev (storage)!")
 ABI_CHECK(lapack_single_precision,"BUG(2) in abi_cheev (precision)!")
 ABI_CHECK(n<=eigen_c_maxsize,"BUG(3) in abi_cheev (maxsize)!")

 work => eigen_c_work ; rwork => eigen_c_rwork
 lwork=eigen_c_lwork

!===== PLASMA
 !FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !  the eigenvalues computation (jobz=N) and not the
 ! full eigenvectors bases determination (jobz=V)
 if (ABI_LINALG_PLASMA_ISON.and.LSAME(jobz,'N')) then
#if defined HAVE_LINALG_PLASMA
   if (eigen_c_lwork==0) then
     ABI_ALLOCATE(work,(n**2))
   end if
   call PLASMA_Alloc_Workspace_cheev(n,n,plasma_work,info)
   info = PLASMA_cheev_c(jobz_plasma(jobz),uplo_plasma(uplo),n,c_loc(a),lda,c_loc(w),&
&                        plasma_work,c_loc(work),n)
   call PLASMA_Dealloc_handle(plasma_work,info)
   if (eigen_c_lwork==0) then
     ABI_DEALLOCATE(work)
   end if
#endif

!===== LAPACK
 else
   if (eigen_c_lwork==0) then
     lwork=2*n-1
     ABI_ALLOCATE(work,(lwork))
   end if
   if (eigen_c_lrwork==0) then
     ABI_ALLOCATE(rwork,(3*n-2))
   end if
   call cheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
   if (eigen_c_lwork==0) then
     ABI_DEALLOCATE(work)
   end if
   if (eigen_c_lrwork==0) then
     ABI_DEALLOCATE(rwork)
   end if
 end if

 ABI_CHECK(info==0,"abi_cheev returned info!=!0")

end subroutine abi_cheev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zheev
!! NAME
!! abi_zheev
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_zheev(jobz,uplo,n,a,lda,w)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zheev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda
 complex(dpc), intent(inout) :: a(lda,*)
 real(dp), intent(out) :: w(n)

!Local variables-------------------------------
 integer :: info,lwork
 real(dp),pointer :: rwork(:)
 complex(dpc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_full_storage,"BUG(1) in abi_zheev (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_zheev (precision)!")
 ABI_CHECK(n<=eigen_z_maxsize,"BUG(3) in abi_zheev (maxsize)!")

 work => eigen_z_work ; rwork => eigen_z_rwork
 lwork=eigen_z_lwork

!===== PLASMA
 !FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 ! the eigenvalues computation (jobz=N) and not the
 ! full eigenvectors bases determination (jobz=V)
 if (ABI_LINALG_PLASMA_ISON.and.LSAME(jobz,'N')) then
#if defined HAVE_LINALG_PLASMA
   if (eigen_z_lwork==0) then
     ABI_ALLOCATE(work,(n**2))
   end if
   call PLASMA_Alloc_Workspace_zheev(n,n,plasma_work,info)
   info = PLASMA_zheev_c(jobz_plasma(jobz),uplo_plasma(uplo),&
&                        plasma_work,c_loc(work),n)
   call PLASMA_Dealloc_handle(plasma_work,info)
   if (eigen_z_lwork==0) then
     ABI_DEALLOCATE(work)
   end if
#endif

!===== LAPACK
 else
   if (eigen_z_lwork==0) then
     lwork=2*n-1
     ABI_ALLOCATE(work,(lwork))
   end if
   if (eigen_z_lrwork==0) then
     ABI_ALLOCATE(rwork,(3*n-2))
   end if
  call zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
   if (eigen_z_lwork==0) then
     ABI_DEALLOCATE(work)
   end if
   if (eigen_z_lrwork==0) then
     ABI_DEALLOCATE(rwork)
   end if
 end if

 ABI_CHECK(info==0,"abi_zheev returned info !=0!")

end subroutine abi_zheev
!!***
