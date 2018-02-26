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
!!  Copyright (C) 2001-2018 ABINIT group (LNguyen,FDahm (CS))
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

  subroutine abi_dhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info, &
&       x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhegv'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
 integer :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: lwork
 real(dp), target, intent(inout) :: a(lda,*) ! FIXME should be cplex * lda
 real(dp), target, intent(inout) :: b(ldb,*)
 real(dp), target, intent(inout) :: work(*)
 real(dp),target,optional, intent(inout) :: rwork(lwork)
 real(dp),target,intent(out) :: w(n)
 integer, intent(out) :: info

!Optional Arguments ------------------------------------
 integer, optional, intent(in) :: x_cplx
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: use_slk
 integer, optional, intent(in) :: timopt,tim_xeigen,use_gpu

!Local variables ------------------------------------
 character(len=500) :: msg
 integer :: istwf_k_ , usegpu_, use_slk_, cplx_
 real(dp) :: tsec(2)

#ifdef HAVE_LINALG_MAGMA
 integer :: liwork,lwork_, lrwork_
 integer, dimension(:), allocatable :: iwork
#endif
#ifdef HAVE_LINALG_PLASMA
 integer :: jobz_plasma_a
 type(c_ptr) :: plasma_work
#endif

! *********************************************************************

 if (present(tim_xeigen).and.present(timopt)) then
    if(abs(timopt)==3) then
      call timab(tim_xeigen,1,tsec)
    end if
 end if

 usegpu_=0; if (present(use_gpu)) usegpu_=use_gpu
 use_slk_=0; if (present(use_slk)) use_slk_=use_slk
 istwf_k_=1; if (present(istwf_k)) istwf_k_=istwf_k
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 if( n > eigen_d_maxsize ) then
    write(msg,'(a,2i3)')' Eigen size higher than max size set!!',n,eigen_d_maxsize
    MSG_ERROR(msg)
 endif

#ifdef HAVE_LINALG_MAGMA
 if (usegpu_==1) then
    lrwork_= 1 + 5*n + 2*(n**2)
    liwork= 3 + 5*n
    ABI_ALLOCATE(iwork,(liwork))
    if (cplx_ == 2 .and. present(rwork)) then
       lwork_=33*n + (n**2)
       call magmaf_zhegvd(itype,jobz,uplo,n,a,lda,b,ldb,w, &
&            work(1:2*lwork_),lwork_,rwork(1:lrwork_),lrwork_,iwork(1:liwork),liwork,info)
    else
       lwork_= 1+n*(6 + 2*n)
       call magmaf_dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,&
&            work(1:lwork_),lwork_,iwork(1:liwork),liwork,info)
    endif
    ABI_DEALLOCATE(iwork)
 else
#endif

#ifdef HAVE_LINALG_SCALAPACK
 if(use_slk_ == 1) then
   ABI_CHECK(PRESENT(x_cplx),"x_cplx must be present")
   call compute_eigen2(abi_communicator,abi_processor,cplx_,n,n,a,b,w,istwf_k_)
   info = 0 ! This is to avoid unwanted warning but it's not clean
 else
#endif

#ifdef HAVE_LINALG_PLASMA
 ! FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !the eigenvalues computation (jobz=N) and not the
 !full eigenvectors bases determination (jobz=V)
 if (LSAME(jobz,'N')) then
    jobz_plasma_a = jobz_plasma(jobz)

    MSG_ERROR("This code is broken")
    if (cplx_ == 2 .and. present(rwork)) then
       call PLASMA_Alloc_Workspace_zhegv(n,n,plasma_work,info)

!       info = PLASMA_zhegv_c(itype,jobz_plasma_a,uplo_plasma(uplo),&
!&        n,c_loc(a),lda,c_loc(b),ldb,c_loc(w),plasma_work,c_loc(rwork),lwork)

       call PLASMA_Dealloc_handle(plasma_work,info)
    else
       call PLASMA_Alloc_Workspace_dsygv(n,n,plasma_work,info)

!       info = PLASMA_dsygv_c(itype,jobz_plasma_a,uplo_plasma(uplo),&
!&        n,c_loc(a),lda,c_loc(b),ldb,c_loc(w),plasma_work,c_loc(lwork),lwork)

       call PLASMA_Dealloc_handle(plasma_work,info)
    endif
 else
#endif
 if ( cplx_ == 2 .and. present(rwork)) then
    call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork/2,rwork,info)
 else
    call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
 endif
#ifdef HAVE_LINALG_PLASMA
 end if
#endif
#ifdef HAVE_LINALG_SCALAPACK
 end if
#endif
#ifdef HAVE_LINALG_MAGMA
 end if
#endif

 if(info/=0) then
    write(msg,'(a,i3)')' Problem in abi_xhegv, info= ',info
    MSG_ERROR(msg)
 endif

 if (present(tim_xeigen).and.present(timopt)) then
    if(abs(timopt)==3) then
      call timab(tim_xeigen,2,tsec)
    end if
 end if

end subroutine abi_dhegv
!!***

!!****f* m_abi_linalg/abi_dhegv_alloc
!! NAME
!! abi_dhegv_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
 subroutine abi_dhegv_alloc(itype,jobz,uplo,n,a,b,w, &
&       x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhegv_alloc'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), DEV_CONTARRD intent(inout) :: a(:,:)  ! here I need lda to be consistent
 real(dp), DEV_CONTARRD intent(inout) :: b(:,:)
 real(dp), intent(out) :: w(n)
 integer :: info

!Optional Arguments ------------------------------------
 integer, optional, intent(in) :: x_cplx
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: timopt,tim_xeigen,use_slk, use_gpu

 call  abi_dhegv(itype,jobz,uplo,n,a,n,b,n,w,eigen_d_work,eigen_d_lwork,eigen_z_rwork,info, &
&         x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

end subroutine abi_dhegv_alloc
!!***

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
subroutine abi_chegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_chegv'
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: lwork
 complex(spc),target, intent(inout) :: a(lda,*)
 complex(spc),target, intent(inout) :: b(ldb,*)
 complex(spc),target, intent(inout) :: work(*)
 real(sp), intent(inout) :: rwork(lwork)
 real(sp), intent(out) :: w(n)
 integer, intent(out) :: info

#ifdef HAVE_LINALG_PLASMA
 integer :: jobz_plasma_a
 type(c_ptr) :: plasma_work
#endif
    ! *********************************************************************

#ifdef HAVE_LINALG_PLASMA
 ! FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !the eigenvalues computation (jobz=N) and not the
 !full eigenvectors bases determination (jobz=V)
 if (LSAME(jobz,'N')) then
    jobz_plasma_a = jobz_plasma(jobz)

    call PLASMA_Alloc_Workspace_chegv(n,n,plasma_work,info)

    MSG_ERROR("This code is broken")
    !call PLASMA_chegv(itype,jobz_plasma_a,uplo_plasma(uplo),n,a,lda,b,ldb,w,plasma_work,rwork,lwork,info)

    call PLASMA_Dealloc_handle(plasma_work,info)
 else
#endif
   call chegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
#ifdef HAVE_LINALG_PLASMA
 end if
#endif

end subroutine abi_chegv
!!***

!!****f* m_abi_linalg/abi_chegv_alloc
!! NAME
!! abi_chegv_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!

 subroutine abi_chegv_alloc(itype,jobz,uplo,n,a,b,w)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_chegv_alloc'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(spc), DEV_CONTARRD intent(inout) :: a(:,:) ! here I need lda to be consistent
 complex(spc), DEV_CONTARRD intent(inout) :: b(:,:)
 real(sp), intent(out) :: w(n)
 integer :: info

! *********************************************************************

 call abi_chegv(itype,jobz,uplo,n,a,n,b,n,w,eigen_c_work,eigen_c_lwork,eigen_c_rwork,info)

end subroutine abi_chegv_alloc
!!***

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

subroutine abi_zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zhegv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: lwork
 complex(dpc), target,intent(inout) :: a(lda,*)
 complex(dpc), target,intent(inout) :: b(ldb,*)
 complex(dpc), target,intent(inout) :: work(*)
 real(dp), intent(inout) :: rwork(lwork)
 real(dp), intent(out) :: w(n)
 integer, intent(out) :: info

#ifdef HAVE_LINALG_PLASMA
 !Optional Arguments ------------------------------------
 integer :: jobz_plasma_a
 type(c_ptr) :: plasma_work
#endif

 ! *********************************************************************

#ifdef HAVE_LINALG_PLASMA
 ! FDahm & LNGuyen  (November 2012) :
 !  In Plasma v 2.4.6, eigen routines support only
 !the eigenvalues computation (jobz=N) and not the
 !full eigenvectors bases determination (jobz=V)
 if (LSAME(jobz,'N')) then
    jobz_plasma_a = jobz_plasma(jobz)

    call PLASMA_Alloc_Workspace_zhegv(n,n,plasma_work,info)

    MSG_ERROR("This code is broken")
    !call PLASMA_zhegv(itype,jobz_plasma_a,uplo_plasma(uplo),n,a,lda,b,ldb,w,plasma_work,rwork,lwork,info)

    call PLASMA_Dealloc_handle(plasma_work,info)
 else
#endif
   call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
#ifdef HAVE_LINALG_PLASMA
 end if
#endif

end subroutine abi_zhegv
!!***

!!****f* m_abi_linalg/abi_zhegv_alloc
!! NAME
!! abi_zhegv_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_zhegv_alloc(itype,jobz,uplo,n,a,b,w)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zhegv_alloc'

!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(dpc), DEV_CONTARRD intent(inout) :: a(:,:) ! here I need lda to be consistent
 complex(dpc), DEV_CONTARRD intent(inout) :: b(:,:)
 real(dp), intent(out) :: w(n)
 integer :: info

! *********************************************************************
 call abi_zhegv(itype,jobz,uplo,n,a,n,b,n,w,eigen_z_work,eigen_z_lwork,eigen_z_rwork,info)

end subroutine abi_zhegv_alloc
!!***
