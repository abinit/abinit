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
!!  Copyright (C) 2001-2024 ABINIT group (LNguyen,FDahm,MT)
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
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
  subroutine abi_dhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,istwf_k,use_slk,use_gpu_elpa)

    use m_fstrings,     only : sjoin, itoa

!Arguments ------------------------------------
 integer :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,ldz
 real(dp), intent(inout) :: a(:)
 real(dp), intent(inout) :: b(:)
 real(dp), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(:)
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: use_slk,use_gpu_elpa

!Local variables-------------------------------
 integer :: info,use_slk_,use_gpu_elpa_,istwf_k_
#ifdef HAVE_LINALG_SCALAPACK
 type(matrix_scalapack) :: sca_a,sca_b,sca_ev
 integer :: ierr
#endif

! *********************************************************************

 ABI_CHECK(lapack_packed_storage,"BUG(1) in abi_dhpgv (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_dhpgv (precision)!")
 ABI_CHECK(n<=eigen_d_maxsize,"BUG(3) in abi_dhpgv (maxsize)!")

 info = 0 !to avoid unwanted warning when info is not set by scalapack

 use_slk_ = 0; if (present(use_slk)) use_slk_ = use_slk
 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k
 use_gpu_elpa_=0
#ifdef HAVE_LINALG_ELPA
 if (present(use_gpu_elpa)) use_gpu_elpa_=use_gpu_elpa
#endif

!===== SCALAPACK
 if (ABI_LINALG_SCALAPACK_ISON.and.use_slk_==1.and.n>slk_minsize)  then
#if defined HAVE_LINALG_SCALAPACK
   z = zero
   call sca_a%init(n,n,slk_processor,istwf_k_)
   call sca_b%init(n,n,slk_processor,istwf_k_)
   call sca_ev%init(n,n,slk_processor,istwf_k_)
#ifdef HAVE_LINALG_ELPA
   call matrix_from_global_sym(sca_a,a,istwf_k_)
   call matrix_from_global_sym(sca_b,b,istwf_k_)
#else
   call matrix_from_global(sca_a,a,istwf_k_)
   call matrix_from_global(sca_b,b,istwf_k_)
#endif
   call compute_generalized_eigen_problem(slk_processor,sca_a,sca_b,&
&       sca_ev,w,slk_communicator,istwf_k_,use_gpu_elpa=use_gpu_elpa_)
   call matrix_to_global(sca_a,a,istwf_k_)
   call matrix_to_global(sca_b,b,istwf_k_)
   call matrix_to_reference(sca_ev,z,istwf_k_)
   call xmpi_sum(z,slk_communicator,ierr)
   call sca_a%free()
   call sca_ev%free()
#endif

!===== LAPACK
 else
   if (istwf_k_/=2) then
     call zhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,eigen_z_work,eigen_z_rwork,info)
   else
     call dspgv(itype,jobz,uplo,n,a,b,w,z,ldz,eigen_d_work,info)
   endif
 end if

  if (info < 0) then
    ABI_COMMENT(sjoin("argument #", itoa(-info), "had an illegal value"))
 end if

 if (info > 0) then
    ABI_COMMENT("DSPEV failed to converge")
    if (info <= n) then
       ABI_COMMENT(sjoin("DSPEV failed to converge;", itoa(info), " off-diagonal elements of"))
       ABI_COMMENT(" an intermediate tridiagonal form did not converge to zero.")
    else
       ABI_COMMENT("The factorization of B could not be completed and no eigenvalues or eigenvectors were computed.")
    endif
 end if

 ABI_CHECK(info==0,"abi_dhpgv returned info!=0!")

#ifndef HAVE_LINALG_ELPA
 ABI_UNUSED(use_gpu_elpa)
#endif

end subroutine abi_dhpgv
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_chpgv
!! NAME
!! abi_chpgv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
  subroutine abi_chpgv(itype,jobz,uplo,n,a,b,w,z,ldz)

!Arguments ------------------------------------
 integer,intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,ldz
 complex(spc), intent(inout) :: a(:,:)
 complex(spc), intent(inout) :: b(:,:)
 complex(spc), intent(out) :: z(:,:)
 real(sp), intent(out) :: w(:)

!Local variables-------------------------------
 integer :: info
 real(sp),pointer :: rwork(:)
 complex(spc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_packed_storage,"BUG(1) in abi_chpgv (storage)!")
 ABI_CHECK(lapack_single_precision,"BUG(2) in abi_chpgv (precision)!")
 ABI_CHECK(n<=eigen_c_maxsize,"BUG(3) in abi_chpgv (maxsize)!")

 work => eigen_c_work ; rwork => eigen_c_rwork

!===== LAPACK
 if (eigen_c_lwork==0) then
   ABI_MALLOC(work,(2*n-1))
 end if
 if (eigen_c_lrwork==0) then
   ABI_MALLOC(rwork,(3*n-2))
 end if
 call chpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)
 if (eigen_c_lwork==0) then
   ABI_FREE(work)
 end if
 if (eigen_c_lrwork==0) then
   ABI_FREE(rwork)
 end if

 ABI_CHECK(info==0,"abi_chpgv returned info!=0!")

end subroutine abi_chpgv
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zhpgv
!! NAME
!! abi_zhpgv
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zhpgv(itype,jobz,uplo,n,a,b,w,z,ldz)

!Arguments ------------------------------------
 integer,intent(in) :: itype
 integer, intent(in) :: n,ldz
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 complex(dpc), intent(inout) :: a(:,:)
 complex(dpc), intent(inout) :: b(:,:)
 complex(dpc), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(:)

!Local variables-------------------------------
 integer :: info
 real(dp),pointer :: rwork(:)
 complex(dpc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_packed_storage,"BUG(1) in abi_zhpgv (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_zhpgv (precision)!")
 ABI_CHECK(n<=eigen_z_maxsize,"BUG(3) in abi_zhpgv (maxsize)!")

 work => eigen_z_work ; rwork => eigen_z_rwork

!===== LAPACK
 if (eigen_z_lwork==0) then
   ABI_MALLOC(work,(2*n-1))
 end if
 if (eigen_z_lrwork==0) then
   ABI_MALLOC(rwork,(3*n-2))
 end if
 call zhpgv(itype,jobz,uplo,n,a,b,w,z,ldz,work,rwork,info)
 if (eigen_z_lwork==0) then
   ABI_FREE(work)
 end if
 if (eigen_z_lrwork==0) then
   ABI_FREE(rwork)
 end if

 ABI_CHECK(info==0,"abi_zhpgv returned info!=0!")

end subroutine abi_zhpgv
!!***
