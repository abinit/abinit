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
!!  Copyright (C) 2001-2017 ABINIT group (LNguyen,FDahm (CS))
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
  subroutine abi_dheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info, &
&       x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dheev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: lda
 integer, intent(in) :: lwork
 integer, intent(out) :: info
 real(dp),target, intent(inout) :: a(lda,*)  ! FIXME should be cplex * lda
 real(dp),target,intent(inout) :: work(*)
 real(dp),optional,target,intent(inout) :: rwork(lwork)
 real(dp), target,intent(out) :: w(n)
 integer, optional, intent(in) :: x_cplx
 integer, optional, intent(in) :: timopt,tim_xeigen,use_gpu
 integer, optional, intent(in) :: use_slk
 integer, optional, intent(in) :: istwf_k

!Local variables-------------------------------
 character(len=500) :: msg
 integer :: istwf_k_ , usegpu_,use_slk_
 real(dp) :: tsec(2)
 integer  :: cplx_
#ifdef HAVE_LINALG_MAGMA
 integer :: liwork_,lwork_, lrwork_,iwk_(1)
 integer, dimension(:), allocatable :: iwork
 real(dp) :: rwk_(1)
 complex(dpc) :: cwk_(1)
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

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 usegpu_=0;if (present(use_gpu)) usegpu_=use_gpu
 istwf_k_=1;if (present(istwf_k)) istwf_k_=istwf_k
 use_slk_ = 0; if(present(use_slk)) use_slk_ = 1
 
 if( n > eigen_d_maxsize ) then
    write(msg,'(a,2i3)')' Eigen size higher than max size set!!',n,eigen_d_maxsize
    MSG_ERROR(msg)
 endif
    
#ifdef HAVE_LINALG_MAGMA
 if (usegpu_==1) then
    !work only if  lwork=n**2 + 33*n?
    if (cplx_ == 2 .and. present(rwork)) then
       !old
       !lwork_=33*n + (n**2)
       !lrwork_= 1 + 5*n + 2*(n**2)
       !liwork_ = 3 + 5*n
       !new
       call magmaf_zheevd(jobz,uplo,n,a,lda,w,cwk_(1),-1,rwk_(1),-1,iwk_(1),-1,info)
       lwork_=int(real(cwk_(1))) ; lrwork_=int(rwk_(1)) ; liwork_=iwk_(1)
       ABI_ALLOCATE(iwork,(liwork_))
       call magmaf_zheevd(jobz,uplo,n,a,lda,w,work(1:2*lwork_),lwork_,rwork(1:lrwork_),lrwork_,iwork,liwork_,info)
       ABI_DEALLOCATE(iwork)
    else
       !old
       !lwork_= 1+n*(6 + 2*n)
       !lrwork_= 1 + 5*n + 2*(n**2)
       !liwork_ = 3 + 5*n
       !new
       call magmaf_dsyevd(jobz,uplo,n,a,lda,w,rwk_(1),-1,iwk_(1),-1,info)
       lwork_=int(rwk_(1)) ; liwork_=iwk_(1)
       ABI_ALLOCATE(iwork,(liwork_))
       call magmaf_dsyevd(jobz,uplo,n,a,lda,w,work(1:lwork_),lwork_,iwork(1:liwork_),liwork_,info)
       ABI_DEALLOCATE(iwork)
    endif
 else
#endif

#ifdef HAVE_LINALG_SCALAPACK
 if( use_slk_ == 1.and.( n > maxval(abi_processor%grid%dims(1:2))) )  then
    ABI_CHECK(present(x_cplx),"x_cplx must be present")
    call compute_eigen1(abi_communicator,abi_processor,cplx_,n,n,a,w,istwf_k_)
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

    if ( cplx_ == 2 .and. present(rwork)) then
       call PLASMA_Alloc_Workspace_zheev(n,n,plasma_work,info)
       info = PLASMA_zheev_c(jobz_plasma_a,uplo_plasma(uplo),n,c_loc(a),lda,c_loc(w),plasma_work,c_loc(rwork),lwork)
       ABI_CHECK(info==0,"PLASMA_zheev_c returned info !=0")
       call PLASMA_Dealloc_handle(plasma_work,info)
    else
       call PLASMA_Alloc_Workspace_dsyev(n,n,plasma_work,info)
       info = PLASMA_dsyev_c(jobz_plasma_a,uplo_plasma(uplo),n,c_loc(a),lda,c_loc(w),plasma_work,c_loc(rwork),lwork)
       ABI_CHECK(info==0,"PLASMA_dsyev_c returned info !=0")
       call PLASMA_Dealloc_handle(plasma_work,info)
    endif
 else
#endif
   if ( cplx_ == 2 .and. present(rwork)) then
      call zheev(jobz,uplo,n,a,lda,w,work,lwork/2,rwork,info)
   else
#ifdef FC_NAG
      ! MG: This hack needed to pass paral[25] paral[29] and mpiio on nag@petrus with np=4
      if (n < 0) write(std_out, *)"work: ",work(1:3)
#endif
      call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
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
    write(msg,'(a,i0)')' Problem in abi_xheev, info= ',info
    !MSG_WARNING(msg)
    MSG_ERROR(msg)
 endif
 
 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_xeigen,2,tsec)
   end if
 end if

end subroutine abi_dheev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_dheev_alloc
!! NAME
!! abi_dheev_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_dheev_alloc(jobz,uplo,n,a,w, &
&       x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dheev_alloc'

!End of the abilint section

 implicit none
 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), DEV_CONTARRD intent(inout) :: a(:,:)  ! Here I need lda to be consistent
 real(dp), intent(out) :: w(n) 
 integer, optional, intent(in) :: x_cplx
 integer, optional, intent(in) :: timopt,tim_xeigen,use_slk,use_gpu
 integer, optional, intent(in) :: istwf_k

 integer :: info
    
! *********************************************************************
    
 call abi_dheev(jobz,uplo,n,a,n,w,eigen_d_work,eigen_d_lwork,eigen_z_rwork,info, &
&         x_cplx,istwf_k,timopt,tim_xeigen,use_slk,use_gpu)

end subroutine abi_dheev_alloc
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

subroutine abi_cheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_cheev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: lda
 integer, intent(in) :: lwork
 complex(spc), target, intent(inout) :: a(lda,*)
 complex(spc), target, intent(inout) :: work(*)
 real(sp),target, intent(inout) :: rwork(lwork)
 real(sp),target, intent(out) :: w(n)
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

    call PLASMA_Alloc_Workspace_cheev(n,n,plasma_work,info)

    info = PLASMA_cheev_c(jobz_plasma_a,uplo_plasma(uplo),&
&     n,c_loc(a),lda,c_loc(w),plasma_work,c_loc(rwork),lwork)

    call PLASMA_Dealloc_handle(plasma_work,info)
 else
#endif
   call cheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
#ifdef HAVE_LINALG_PLASMA
 end if
#endif

end subroutine abi_cheev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_cheev_alloc
!! NAME
!! abi_cheev_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
  subroutine abi_cheev_alloc(jobz,uplo,n,a,w)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_cheev_alloc'

!End of the abilint section

 implicit none
 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(spc), DEV_CONTARRD intent(inout) :: a(:,:) ! Here I need lda to be consistent
 real(sp), intent(out) :: w(n)

 integer :: info

! *********************************************************************
    
 call abi_cheev(jobz,uplo,n,a,n,w,eigen_c_work,eigen_c_lwork,eigen_c_rwork,info)

end subroutine abi_cheev_alloc
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

subroutine abi_zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zheev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: lda
 integer, intent(in) :: lwork
 complex(dpc),target,intent(inout) :: a(lda,*)
 complex(dpc), target, intent(inout) :: work(*)
 real(dp),target, intent(inout) :: rwork(lwork)        ! TODO Check this!
 real(dp),target,intent(out) :: w(n)
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

    call PLASMA_Alloc_Workspace_zheev(n,n,plasma_work,info)

    info = PLASMA_zheev_c(jobz_plasma_a,uplo_plasma(uplo),&
&     n,c_loc(a),lda,c_loc(w),plasma_work,c_loc(rwork),lwork)

    call PLASMA_Dealloc_handle(plasma_work,info)
 else
#endif
   call zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
#ifdef HAVE_LINALG_PLASMA
 end if
#endif

end subroutine abi_zheev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zheev_alloc
!! NAME
!! abi_zheev_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_zheev_alloc(jobz,uplo,n,a,w)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zheev_alloc'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(dpc), DEV_CONTARRD intent(inout) :: a(:,:)   ! Here I need lda to be consistent
 real(dp), intent(out) :: w(n)

 integer :: info

! *********************************************************************
    
 call abi_zheev(jobz,uplo,n,a,n,w,eigen_z_work,eigen_z_lwork,eigen_z_rwork,info)

end subroutine abi_zheev_alloc
!!***

!----------------------------------------------------------------------
