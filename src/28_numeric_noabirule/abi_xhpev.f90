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
!!  Copyright (C) 2001-2018 ABINIT group (LNguyen,FDahm (CS))
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

  subroutine abi_dhpev(jobz,uplo,n,a,w,z,ldz,work,info,rwork,istwf_k)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhpev'
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: ldz
 real(dp) :: a(*)    !vz_d
 real(dp) :: z(*)    !vz_d
 real(dp), intent(inout) :: work(*)
 real(dp),optional, intent(inout) :: rwork(*)
 real(dp) :: w(*)    !vz_d
 integer, intent(out) :: info
 !Optional Arguments ------------------------------------
 integer, optional, intent(in) :: istwf_k

!Local Arguments ------------------------------------
 integer :: istwf_k_

 ! *********************************************************************

 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

 ! MG: FIXME: THIS IS BROKEN: rwork and istwfk_k!
 !if ( present(istwf_k) .and. (istwf_k .ne. 2) .and. present(rwork)) then
 if (istwf_k_ /= 2) then
    ABI_CHECK(present(rwork),"rwork must be present")
    call zhpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)
 else
    call dspev(jobz,uplo,n,a,w,z,ldz,work(1:3*n),info)
 endif

end subroutine abi_dhpev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_dhpev_alloc_1d
!! NAME
!! abi_dhpev_alloc_1d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_dhpev_alloc_1d(jobz,uplo,n,a,w,z,istwf_k,use_slk)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhpev_alloc_1d'

!End of the abilint section

 implicit none
 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(:)
 real(dp), intent(out) :: w(:)
 real(dp), intent(out) :: z(:,:)
 integer, optional, intent(in) :: istwf_k
 integer, optional, intent(in) :: use_slk

 !Local Arguments ------------------------------------
 character(len=500) :: msg
 integer :: info,use_slk_,istwf_k_
#ifdef HAVE_LINALG_SCALAPACK
 type(matrix_scalapack)    :: sca_a,sca_ev
 real(dp),allocatable :: tmp_evec(:,:)
 integer :: dim_evec1,ierr
#endif

! *********************************************************************
 ! write(*,*) "ENTERING abi_dhpev_alloc_1d n=",n
 if( n > eigen_d_maxsize ) then
    write(msg,'(a,2i3)')' Eigen size higher than max size set!!',n,eigen_d_maxsize
    MSG_ERROR(msg)
 endif
 info = 0 ! to avoid unwanted warnings but scalapack doesn't check return

 use_slk_ = 0; if (present(use_slk)) use_slk_ = use_slk
 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

#ifdef HAVE_LINALG_SCALAPACK
 if (use_slk_ == 1) then
   ! if istwfk=1, then dim_evec1=2*n and if istwfk=2, dim_evec1=n
   dim_evec1= 2*n/istwf_k_
   ABI_ALLOCATE(tmp_evec,(dim_evec1,n))
   tmp_evec = 0._dp
   call init_matrix_scalapack(sca_a,n,n,abi_processor,istwf_k_,10)
   call init_matrix_scalapack(sca_ev,n,n,abi_processor,istwf_k_,10)

#ifdef HAVE_LINALG_ELPA
   call matrix_from_global_sym(sca_a,a,istwf_k_)
#else
   call matrix_from_global(sca_a,a,istwf_k_)
#endif
   call compute_eigen_problem(abi_processor,sca_a,&
&        sca_ev,w,abi_communicator,istwf_k_)

   call matrix_to_global(sca_a,a,istwf_k_)
   call matrix_to_reference(sca_ev,tmp_evec,istwf_k_)

   call xmpi_sum(tmp_evec,z,dim_evec1*n,abi_communicator,ierr)

   CALL destruction_matrix_scalapack(sca_a)
   CALL destruction_matrix_scalapack(sca_ev)
   ABI_DEALLOCATE(tmp_evec)
  else
#endif
    call abi_dhpev(jobz,uplo,n,a,w,z,n,eigen_d_work,info,&
&     rwork=eigen_z_rwork,istwf_k=istwf_k_)
#ifdef HAVE_LINALG_SCALAPACK
  end if
#endif
  if(info/=0) then
     write(msg,'(a,i3)')' Problem in abi_xhpev, info= ',info
     MSG_WARNING(msg)
  endif

end subroutine abi_dhpev_alloc_1d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_dhpev_alloc_2d
!! NAME
!! abi_dhpev_alloc_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

 subroutine abi_dhpev_alloc_2d(jobz,uplo,n,a,w,z,istwf_k)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_dhpev_alloc_2d'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(:,:)
 real(dp), intent(out) :: w(:)
 real(dp), intent(out) :: z(:,:)
 integer :: info
 integer, optional, intent(in) :: istwf_k

!Local Arguments ------------------------------------
 integer :: istwf_k_

 ! *********************************************************************

 istwf_k_ = 1; if (present(istwf_k)) istwf_k_ = istwf_k

 call abi_dhpev(jobz,uplo,n,a,w,z,n,eigen_d_work,info,&
&  rwork=eigen_z_rwork,istwf_k=istwf_k_)

end subroutine abi_dhpev_alloc_2d
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

  subroutine abi_chpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_chpev'
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: ldz
 complex(spc), intent(inout) :: a(:,:)
 complex(spc), intent(out) :: z(:,:)
 complex(spc), intent(inout) :: work(:)
 real(sp), intent(inout) :: rwork(:)
 real(sp), intent(out) :: w(:)
 integer, intent(out) :: info

! *********************************************************************

 call chpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)

end subroutine abi_chpev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_chpev_alloc
!! NAME
!! abi_chpev_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

 subroutine abi_chpev_alloc(jobz,uplo,n,a,w,z)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_chpev_alloc'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(spc), intent(inout) :: a(:,:)
 complex(spc), intent(out) :: z(:,:)
 real(sp), intent(out) :: w(n)
 integer :: info

! *********************************************************************

 call abi_chpev(jobz,uplo,n,a,w,z,n,eigen_c_work,eigen_c_rwork,info)

end subroutine abi_chpev_alloc
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

  subroutine abi_zhpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zhpev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 integer, intent(in) :: ldz
 complex(dpc), intent(inout) :: a(:,:)
 complex(dpc), intent(out) :: z(:,:)
 complex(dpc), intent(inout) :: work(:)
 real(dp), intent(inout) :: rwork(:)
 real(dp), intent(out) :: w(n)
 integer, intent(out) :: info

! *********************************************************************

 call zhpev(jobz,uplo,n,a,w,z,ldz,work,rwork,info)

end subroutine abi_zhpev
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zhpev_alloc
!! NAME
!! abi_zhpev_alloc
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

 subroutine abi_zhpev_alloc(jobz,uplo,n,a,w,z)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_zhpev_alloc'

!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 complex(dpc), intent(inout) :: a(:,:)
 complex(dpc), intent(out) :: z(:,:)
 real(dp), intent(out) :: w(n)
 integer :: info

! *********************************************************************

 call abi_zhpev(jobz,uplo,n,a,w,z,n,eigen_z_work,eigen_z_rwork,info)

end subroutine abi_zhpev_alloc
!!***

!----------------------------------------------------------------------
