!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xhpev
!! NAME
!!  abi_xhpev
!!
!! FUNCTION
!!  abi_xhpev is the generic function that compute
!!  all eigenvalues and, optionally, eigenvectors of a
!!  symmetric or hermitian matrix A in packed storage
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (LNguyen,FDahm,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/abi_dhpev
!! NAME
!! abi_dhpev
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_dhpev(jobz,uplo,n,a,w,z,ldz,istwf_k,use_slk)

 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,ldz
 real(dp), intent(inout) :: a(:)
 real(dp), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(:)
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: use_slk

!Local variables-------------------------------
 integer :: info,use_slk_,istwf_k_
#ifdef HAVE_LINALG_SCALAPACK
 type(matrix_scalapack) :: sca_a,sca_ev
 real(dp),allocatable :: tmp_evec(:,:)
 integer :: dim_evec1,ierr
#endif

! *********************************************************************

 ABI_CHECK(lapack_packed_storage,"BUG(1) in abi_dhpev (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_dhpev (precision)!")
 ABI_CHECK(n<=eigen_d_maxsize,"BUG(3) in abi_dhpev (maxsize)!")

 info = 0 ! to avoid unwanted warnings but scalapack doesn't check return

 use_slk_ = 0; if (present(use_slk)) use_slk_ = use_slk
 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

!===== SCALAPACK
 if (ABI_LINALG_SCALAPACK_ISON.and.use_slk_==1.and.n>slk_minsize)  then
#if defined HAVE_LINALG_SCALAPACK
   ! if istwfk=1, then dim_evec1=2*n and if istwfk=2, dim_evec1=n
   dim_evec1= 2*n/istwf_k_
   ABI_ALLOCATE(tmp_evec,(dim_evec1,n))
   tmp_evec = zero
   call init_matrix_scalapack(sca_a,n,n,slk_processor,istwf_k_, tbloc=10)
   call init_matrix_scalapack(sca_ev,n,n,slk_processor,istwf_k_, tbloc=10)
#ifdef HAVE_LINALG_ELPA
   call matrix_from_global_sym(sca_a,a,istwf_k_)
#else
   call matrix_from_global(sca_a,a,istwf_k_)
#endif
   call compute_eigen_problem(slk_processor,sca_a,sca_ev,w,slk_communicator,istwf_k_)
   call matrix_to_global(sca_a,a,istwf_k_)
   call matrix_to_reference(sca_ev,tmp_evec,istwf_k_)
   call xmpi_sum(tmp_evec,z,dim_evec1*n,slk_communicator,ierr)
   call sca_a%free()
   call sca_ev%free()
   ABI_DEALLOCATE(tmp_evec)
#endif

!===== LAPACK
 else
   if (istwf_k_/=2) then
     call zhpev(jobz,uplo,n,a,w,z,ldz,eigen_z_work,eigen_z_rwork,info)
   else
     call dspev(jobz,uplo,n,a,w,z,ldz,eigen_d_work,info)
   end if
 end if

 ABI_CHECK(info==0,"dhpev returned info!=0")

end subroutine abi_dhpev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_chpev
!! NAME
!! abi_chpev
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_chpev(jobz,uplo,n,a,w,z,ldz)

 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,ldz
 complex(spc), intent(inout) :: a(:,:)
 complex(spc), intent(out) :: z(:,:)
 real(sp), intent(out) :: w(:)

!Local variables-------------------------------
 integer :: info
 real(sp),pointer :: rwork(:)
 complex(spc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_packed_storage,"BUG(1) in abi_chpev (storage)!")
 ABI_CHECK(lapack_single_precision,"BUG(2) in abi_chpev (precision)!")
 ABI_CHECK(n<=eigen_c_maxsize,"BUG(3) in abi_chpev (maxsize)!")

 work => eigen_c_work ; rwork => eigen_c_rwork

!===== LAPACK
 if (eigen_c_lwork==0) then
   ABI_ALLOCATE(work,(2*n-1))
 end if
 if (eigen_c_lrwork==0) then
   ABI_ALLOCATE(rwork,(3*n-2))
 end if
 call chpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)
 if (eigen_c_lwork==0) then
   ABI_DEALLOCATE(work)
 end if
 if (eigen_c_lrwork==0) then
   ABI_DEALLOCATE(rwork)
 end if

 ABI_CHECK(info==0,"abi_chpev returned info!=0!")

end subroutine abi_chpev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zhpev
!! NAME
!! abi_zhpev
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_zhpev(jobz,uplo,n,a,w,z,ldz)

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,ldz
 complex(dpc), intent(inout) :: a(:,:)
 complex(dpc), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(:)

!Local variables-------------------------------
 integer :: info
 real(dp),pointer :: rwork(:)
 complex(dpc),pointer :: work(:)

! *********************************************************************

 ABI_CHECK(lapack_packed_storage,"BUG(1) in abi_zhpev (storage)!")
 ABI_CHECK(lapack_double_precision,"BUG(2) in abi_zhpev (precision)!")
 ABI_CHECK(n<=eigen_z_maxsize,"BUG(3) in abi_zhpev (maxsize)!")

 work => eigen_z_work ; rwork => eigen_z_rwork

!===== LAPACK
 if (eigen_z_lwork==0) then
   ABI_ALLOCATE(work,(2*n-1))
 end if
 if (eigen_z_lrwork==0) then
   ABI_ALLOCATE(rwork,(3*n-2))
 end if
 call zhpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)
 if (eigen_z_lwork==0) then
   ABI_DEALLOCATE(work)
 end if
 if (eigen_z_lrwork==0) then
   ABI_DEALLOCATE(rwork)
 end if

 ABI_CHECK(info==0,"abi_zhpev returned info!=0!")

end subroutine abi_zhpev
!!***
