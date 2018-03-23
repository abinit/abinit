!{\src2tex{textfont=tt}}
!!****f* ABINIT/matrginv
!! NAME
!! matrginv
!!
!! FUNCTION
!! Invert a general matrix of real*8 elements.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! lda=leading dimension of complex matrix a
!! n=size of complex matrix a
!! a=matrix of real elements
!! OUTPUT
!! a=inverse of a input matrix
!!
!! SIDE EFFECTS
!! a(lda,n)= array of real elements, input, inverted at output
!!
!!
!! PARENTS
!!      calc_optical_mels,ddb_elast,ddb_piezo,get_tau_k,linear_optics_paw
!!      m_haydock,m_vcoul,matpointsym,mka2f_tr,mlwfovlp_ylmfar,setup_bse
!!      strainsym
!!
!! CHILDREN
!!      dbgmdi,dbgmlu,dgeicd,dgetrf,dgetri
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine matrginv(a,lda,n)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matrginv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lda,n
!arrays
 real(dp),intent(inout) :: a(lda,n)

!Local variables-------------------------------
!scalars
 integer :: ierr,nwork
#if defined HAVE_LINALG_ESSL
 real(dp) :: rcond
#endif
 character(len=500) :: message
!arrays
 integer,allocatable :: ipvt(:)
#if defined HAVE_LINALG_ESSL
 real(dp) :: det(2)
#elif defined HAVE_LINALG_ASL
 real(dp) :: det(2)
#endif
 real(dp),allocatable :: work(:)

! *************************************************************************

#if defined HAVE_LINALG_ESSL
 nwork=200*n
#else
 nwork=n
#endif

 ABI_ALLOCATE(work,(nwork))
 ABI_ALLOCATE(ipvt,(n))


#if defined HAVE_LINALG_ESSL

 call dgeicd(a,lda,n,0,rcond,det,work,nwork)
 if(abs(rcond)==zero) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ESSL routine dgeicd failed.',ch10,&
&   '  Action : Contact ABINIT group '
   MSG_ERROR(message)
 end if

#elif defined HAVE_LINALG_ASL

 call dbgmlu(a,lda,n,ipvt,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ASL routine dbgmlu failed.',ch10,&
&   '  Action : Contact ABINIT group '
   MSG_ERROR(message)
 end if
 call dbgmdi(a,lda,n,ipvt,det,-1,work,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ASL routine dbgmdi failed.',ch10,&
&   '  Action : Contact ABINIT group '
   MSG_ERROR(message)
 end if

#else

 call dgetrf(n,n,a,lda,ipvt,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine dgetrf failed.',ch10,&
&   '  Action : Contact ABINIT group '
   MSG_ERROR(message)
 end if
 call dgetri(n,a,lda,ipvt,work,n,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine dgetri failed.',ch10,&
&   '  Action : Contact ABINIT group '
   MSG_ERROR(message)
 end if

#endif

 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(ipvt)

end subroutine matrginv
!!***
