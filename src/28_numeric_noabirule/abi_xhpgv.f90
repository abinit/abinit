!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xhpgv
!! NAME
!!  abi_xhpgv
!!
!! FUNCTION
!!  abi_xhpgv is the generic function that compute
!!  all eigenvalues and, optionally, eigenvectors of a
!!  generalized symmetric-definite eigenproblem, of the form
!!  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
!!  Here A and B are assumed to be symmetric (or hermitian),
!!  stored in packed format  and B is also positive definite.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2018 ABINIT group (LNguyen,FDahm (CS))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
!!***


!!****f* m_abi_linalg/abi_dhpgv
!! NAME
!! abi_dhpgv
!!
!! PARENTS
!!
!! SOURCE
!!
  subroutine abi_dhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,info,rwork,istwf_k)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhpgv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: ldz
 real(dp) :: a(*)    !vz_d
 real(dp) :: b(*)    !vz_d
 real(dp) :: z(*)    !vz_d
 real(dp), intent(inout) :: work(*)
 real(dp), optional,intent(inout) :: rwork(*)
 real(dp) :: w(*)    !vz_d
 integer, intent(out) :: info
 integer, optional, intent(in) :: istwf_k

!Local Arguments ------------------------------------
 integer :: istwf_k_

! *********************************************************************

 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

 !MG: FIXME This is clearly wrong but tests are OK!
 !if ( present(istwf_k) .and. istwf_k == 2 .and. present(rwork)) then
 if (istwf_k_ /= 2) then
    ABI_CHECK(present(rwork),"rwork must be present")
    call zhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)
 else
    call dspgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,info)
 endif

 ABI_CHECK(info==0,"[z,d]hpgv returned info !=0")

end subroutine abi_dhpgv
!!***

!!****f* m_abi_linalg/abi_dhpgv_alloc_1d
!! NAME
!! abi_dhpgv_alloc_1d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
  subroutine abi_dhpgv_alloc_1d(itype,jobz,uplo,n,a,b,w,z,istwf_k,use_slk)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhpgv_alloc_1d'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(:)
 real(dp), intent(inout) :: b(:)
 real(dp), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(:)
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: use_slk

!Local Arguments ------------------------------------
 character(len=500) :: msg
 integer :: info,use_slk_,istwf_k_
#ifdef HAVE_LINALG_SCALAPACK
 type(matrix_scalapack)    :: sca_a,sca_b,sca_ev
 integer :: ierr
#endif

 use_slk_ = 0; if (present(use_slk)) use_slk_ = use_slk
 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

 if( n > eigen_d_maxsize ) then
   write(msg,'(a,2i3)')' Eigen size higher than max size set!!',n,eigen_d_maxsize
   MSG_ERROR(msg)
 endif
 info = 0 !to avoid unwanted warning when info is not set by scalapack

#ifdef HAVE_LINALG_SCALAPACK
 if (use_slk_ == 1) then
   z = 0._dp
   call init_matrix_scalapack(sca_a,n,n,abi_processor,istwf_k_,10)
   call init_matrix_scalapack(sca_b,n,n,abi_processor,istwf_k_,10)
   call init_matrix_scalapack(sca_ev,n,n,abi_processor,istwf_k_,10)
#ifdef HAVE_LINALG_ELPA
   call matrix_from_global_sym(sca_a,a,istwf_k_)
   call matrix_from_global_sym(sca_b,b,istwf_k_)
#else
   call matrix_from_global(sca_a,a,istwf_k_)
   call matrix_from_global(sca_b,b,istwf_k_)
#endif
   call compute_generalized_eigen_problem(abi_processor,sca_a,sca_b,&
&       sca_ev,w,abi_communicator,istwf_k_)

   call matrix_to_global(sca_a,a,istwf_k_)
   call matrix_to_global(sca_b,b,istwf_k_)
   call matrix_to_reference(sca_ev,z,istwf_k_)

   call xmpi_sum(z,abi_communicator,ierr)

   CALL destruction_matrix_scalapack(sca_a)
   CALL destruction_matrix_scalapack(sca_ev)
 else
#endif
   call  abi_dhpgv(itype,jobz,uplo,n,a,b,w,z,n, &    !vz_d
&            eigen_d_work,info, rwork=eigen_z_rwork,istwf_k=istwf_k_)
#ifdef HAVE_LINALG_SCALAPACK
 end if
#endif

 if(info/=0) then
   write(msg,'(a,i0)')' Problem in abi_xhpgv, info= ',info
   MSG_ERROR(msg)
 endif

end subroutine abi_dhpgv_alloc_1d
!!***

!!****f* m_abi_linalg/abi_dhpgv_alloc_2d
!! NAME
!! abi_dhpgv_alloc_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_dhpgv_alloc_2d(itype,jobz,uplo,n,a,b,w,z,istwf_k)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhpgv_alloc_2d'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(:,:)
 real(dp), intent(inout) :: b(:,:)
 real(dp), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(:)
 integer :: info
 integer, optional, intent(in) :: istwf_k

!Local Arguments ------------------------------------
 integer :: istwf_k_

! *********************************************************************

 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

 call  abi_dhpgv(itype,jobz,uplo,n,a,b,w,z,n,eigen_d_work,info,&
&         rwork=eigen_z_rwork,istwf_k=istwf_k_)

end subroutine abi_dhpgv_alloc_2d
!!***

!!****f* m_abi_linalg/abi_chpgv
!! NAME
!! abi_chpgv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
  subroutine abi_chpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_chpgv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: ldz
 complex(spc), intent(inout) :: a(:,:)
 complex(spc), intent(inout) :: b(:,:)
 complex(spc), intent(out) :: z(:,:)
 complex(spc), intent(inout) :: work(:)
 real(sp), intent(inout) :: rwork(:)
 real(sp), intent(out) :: w(n)
 integer, intent(out) :: info

! *********************************************************************

 call chpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)

end subroutine abi_chpgv
!!***

!!****f* m_abi_linalg/abi_chpgv_alloc
!! NAME
!! abi_chpgv_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_chpgv_alloc(itype,jobz,uplo,n,a,b,w,z)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_chpgv_alloc'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer,intent(in) :: itype
 integer, intent(in) :: n
 complex(spc), intent(inout) :: a(:,:)
 complex(spc), intent(inout) :: b(:,:)
 complex(spc), intent(out) :: z(:,:)
 real(sp), intent(out) :: w(n)

 integer :: info

! *********************************************************************

 call abi_chpgv(itype,jobz,uplo,n,a,b,w,z,n,eigen_c_work,eigen_c_rwork,info)

end subroutine abi_chpgv_alloc
!!***

!!****f* m_abi_linalg/abi_zhpgv
!! NAME
!! abi_zhpgv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_zhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zhpgv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 integer, intent(in) :: n
 integer, intent(in) :: ldz
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 complex(dpc), intent(inout) :: a(:,:)
 complex(dpc), intent(inout) :: b(:,:)
 complex(dpc), intent(out) :: z(:,:)
 complex(dpc), intent(inout) :: work(:)
 real(dp), intent(inout) :: rwork(:)
 real(dp), intent(out) :: w(n)
 integer, intent(out) :: info

 ! *********************************************************************

 call zhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)

end subroutine abi_zhpgv
!!***

!!****f* m_abi_linalg/abi_zhpgv_alloc
!! NAME
!! abi_zhpgv_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_zhpgv_alloc(itype,jobz,uplo,n,a,b,w,z)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zhpgv_alloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(dpc), intent(inout) :: a(:,:)
 complex(dpc), intent(inout) :: b(:,:)
 complex(dpc), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(n)

 integer :: info

! *********************************************************************
 call abi_zhpgv(itype,jobz,uplo,n,a,b,w,z,n,eigen_z_work,eigen_z_rwork,info)

end subroutine abi_zhpgv_alloc
!!***
