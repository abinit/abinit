!{\src2tex{textfont=tt}}
!!****f* ABINIT/hermit
!! NAME
!! hermit
!!
!! FUNCTION
!! Take a matrix in hermitian storage mode (lower triangle stored)
!! and redefine diagonal elements to impose Hermiticity
!! (diagonal terms have to be real).
!! If abs(Im(H(i,i)))>4096*machine precision, print error warning.
!! (Typical 64 bit machine precision is 2^-52 or 2.22e-16)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  chmin(n*n+n)=complex hermitian matrix with numerical noise possibly
!!   rendering Im(diagonal elements) approximately 1e-15 or so
!!  ndim=size of complex hermitian matrix
!!
!! OUTPUT
!!  chmout(n*n+n)=redefined matrix with strictly real diagonal elements.
!!   May be same storage location as chmin.
!!  ierr=0 if no problem, 1 if the imaginary part of some element
!!   too large (at present, stop in this case).
!!
!! TODO
!!  Name is misleading, perhaps hermit_force_diago?
!!  Interface allows aliasing
!!
!! PARENTS
!!      extrapwf,subdiago
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine hermit(chmin,chmout,ierr,ndim)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hermit'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim
 integer,intent(out) :: ierr
!arrays
 real(dp),intent(inout) :: chmin(ndim*ndim+ndim)
 real(dp),intent(inout) :: chmout(ndim*ndim+ndim)

!Local variables-------------------------------
!scalars
 integer,save :: mmesgs=20,nmesgs=0
 integer :: idim,merrors,nerrors
 real(dp),parameter :: eps=epsilon(0.0d0)
 real(dp) :: ch_im,ch_re,moduls,tol
 character(len=500) :: message

! *************************************************************************

 tol=4096.0d0*eps

 ierr=0
 merrors=0

!Copy matrix into possibly new location
 chmout(:)=chmin(:)

!Loop over diagonal elements of matrix (off-diag not altered)
 do idim=1,ndim

   ch_im=chmout(idim*idim+idim  )
   ch_re=chmout(idim*idim+idim-1)

!  check for large absolute Im part and print warning when
!  larger than (some factor)*(machine precision)
   nerrors=0
   if( abs(ch_im) > tol .and. abs(ch_im) > tol8*abs(ch_re)) nerrors=2
   if( abs(ch_im) > tol .or. abs(ch_im) > tol8*abs(ch_re)) nerrors=1

   if( (abs(ch_im) > tol .and. nmesgs<mmesgs) .or. nerrors==2)then
     write(message, '(3a,i0,a,es20.12,a,es20.12,a)' )&
&     'Input Hermitian matrix has nonzero relative Im part on diagonal:',ch10,&
&     'for component',idim,' Im part is',ch_im,', Re part is',ch_re,'.'
     call wrtout(std_out,message,'PERS')
     nmesgs=nmesgs+1
   end if

   if( ( abs(ch_im) > tol8*abs(ch_re) .and. nmesgs<mmesgs) .or. nerrors==2)then
     write(message, '(3a,i0,a,es20.12,a,es20.12,a)' )&
&     'Input Hermitian matrix has nonzero relative Im part on diagonal:',ch10,&
&     'for component',idim,' Im part is',ch_im,', Re part is',ch_re,'.'
     call wrtout(std_out,message,'PERS')
     nmesgs=nmesgs+1
   end if

!  compute modulus $= (\Re^2+\Im^2)^{1/2}$
   moduls=sqrt(ch_re**2+ch_im**2)

!  set Re part to modulus with sign of original Re part
   chmout(idim*idim+idim-1)=sign(moduls,ch_re)

!  set Im part to 0
   chmout(idim*idim+idim)=zero

   merrors=max(merrors,nerrors)

 end do

 if(merrors==2)then
   ierr=1
   write(message, '(3a)' )&
&   'Imaginary part(s) of diagonal Hermitian matrix element(s) is too large.',ch10,&
&   'See previous messages.'
   MSG_BUG(message)
 end if

end subroutine hermit
!!***
