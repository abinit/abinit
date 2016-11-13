!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_skw
!! NAME
!!  m_skw
!!
!! FUNCTION
!!  Shankland-Koelling-Wood Fourier interpolation scheme.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_skw

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_crystal
 use m_sort

 use m_fstrings,       only : itoa, sjoin
 use m_numeric_tools,  only : vdiff_t, vdiff_eval, vdiff_print
 !use m_geometry,      only : normv


 implicit none

 private

!!***

!----------------------------------------------------------------------

!!****t* m_skw/skw_t
!! NAME
!! skw_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: skw_t

  integer :: nr
   ! Number of real-space lattice points.

  integer :: nk
   ! Number of ab-initio k-points

  integer :: nband
   ! Number of bands

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: cplex
  ! 1 if time-reversal symmetry can be used, 2 otherwise.

  integer(dp),allocatable :: rpts(:,:)
  ! rpts(3, nr)
  ! Real-space lattice points ordered with non-decreasing length.

  complex(dpc),allocatable :: coef(:,:,:)
  ! coef(nr,nband,nsppol)

 end type skw_t

 public :: skw_init     ! Initialize the object.
 public :: skw_evalk    ! Evaluate the interpolated eigenvalues.
 public :: skw_free     ! Free memory.
!!***

CONTAINS  !=====================================================================================
!!***

!!****f* m_skw/skw_init
!! NAME
!!  skw_init
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  cplex=1 if time reversal can be used, 2 otherwise.
!!  nband=Number of bands
!!  nkpt=Number of ab-initio k-points.
!!  nsppol=Number of independent spin polarizations.
!!  kpts(3,nkpt)=ab-initio k-points in reduced coordinates.
!!  eig(nband,nkpt,nsppol)=ab-initio eigenvalues.
!!
!! OUTPUT
!!  skw<skw_t>=New object
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine skw_init(skw, cryst, cplex, nband, nkpt, nsppol, kpts, eig)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkpt,nsppol,cplex
 type(crystal_t),intent(in) :: cryst
 type(skw_t),intent(out) :: skw
!arrays
 real(dp),intent(in) :: kpts(3,nkpt)
 real(dp),intent(in) :: eig(nband,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ir,ik,ii,jj,nr,nk,band,spin,isym,nsh!,info
 integer :: ratio
 real(dp) :: arg,ecut,r2
 real(dp),parameter :: c1=0.25_dp,c2=0.25_dp
!arrays
 integer,allocatable :: ipiv(:)
 real(dp) :: oeig(nband)
 real(dp),allocatable :: delta_eig(:,:,:),rhor(:)
 complex(dpc),allocatable :: srk(:,:),hmat(:,:),lambda(:,:,:)

! *********************************************************************

 skw%nk = nkpt; skw%nband = nband; skw%nsppol = nsppol; skw%cplex = cplex

 ! Select R-points then ordered them according to |R|
 ratio = 5
 nr = nkpt * ratio

 !call skw_genrpts(cryst, rmax, nr, rpts)

 ecut = zero; nsh = 0; nr = 1
 !call setshells(ecut,nr,nsh,cryst%nsym,gmet,gprimd,symrel,tag,ucvol)
 !call kpgcount(ecut,exchn2n3d,gmet,istwfk,kpt,ngmax,ngmin,nkpt)
 ABI_MALLOC(skw%rpts, (3, nr))
 skw%nr = nr

 ! Construct star functions.
 ABI_CALLOC(srk, (nr,nk))

 do ir=1,nr
   do ik=1,nkpt
     do isym=1,cryst%nsym
       !arg = two_pi * dot_product(kpts(:,ik), matmul(cryst%symrel(:,:,isym), skw%rpts(:,ir)))
       srk(ir, ik) = srk(ir, ik) + exp(j_dpc * arg)
     end do
     srk(ir, ik) = srk(ir, ik) / cryst%nsym
   end do
 end do

 !do ik=1,nkpt
 !  call skw_mkstar(skw, cryst, kpt, srk(:,ik))
 !end do

 ABI_MALLOC(rhor, (nr))
 do ir=1,nr
   r2 = dot_product(skw%rpts(:,ir), matmul(cryst%rmet, skw%rpts(:,ir)))
   rhor(ir) = c1 * r2 + c2 * r2**2
 end do

 ABI_CALLOC(hmat, (nk-1, nk-1))
 do jj=1,nk-1
   do ii=1,jj
     do ir=2,nr
       hmat(ii, jj) = hmat(ii, jj) + &
          (srk(ir, ii) - srk(ir, nk)) * conjg(srk(ir, jj) - srk(ir, nk)) / rhor(ir)
     end do
   end do
 end do

 ABI_MALLOC(delta_eig, (nk-1,nband,nsppol))
 do spin=1,nsppol
   do band=1,nband
     delta_eig(:,band,spin) = eig(band,1:nk-1,spin) - eig(band,nk,spin)
   end do
 end do

 ! Solve system of linear equations to get lambda coeffients.
 ABI_MALLOC(lambda, (nk-1, nband, nsppol))
 ABI_MALLOC(ipiv, (nk-1))

 ! Solve all the equations for all the bands at once (eq. 10 of PRB 38 2721)

 ! LU factorization of hmat matrix, kept in 2 triangular blocks of the matrix
 !call DGETRF(nk-1,nk-1,hmat,nk-1,ipiv,info)
 !CALL DGETRS('N',bs%nkpt-1,bs%icut2-bs%icut1+1,hmat,bs%nkpt-1,ipiv,lambda(1,bs%icut1),bs%nkpt-1,inf)

 !call dgesv(nk-1, nband*nsppol, hmat, nk-1, ipiv, lambda, nk-1, info)
 !call zgesv(nk-1, nband*nsppol, hmat, nk-1, ipiv, lambda, nk-1, info)
 !ABI_CHECK(info /= 0, sjoin("DGESV returned info", itoa(info)))
 ABI_FREE(ipiv)

 ! Compute coefficients
 ABI_MALLOC(skw%coef, (nr,nband,nsppol))

 do spin=1,nsppol
   do band=1,nband
     !do ir=2,nr
     !  skw%coef(ir,band,spin) = (one/rhor(ir)) * dot_product(lambda(:nk-1), conjg(srk(ir,:) - srk(ir,nk))
     !end do
     skw%coef(1,band,spin) = eig(band,nk,spin) - dot_product(conjg(skw%coef(2:, band,spin)), srk(2:, nk))
   end do
 end do

 ABI_FREE(rhor)
 ABI_FREE(hmat)
 ABI_FREE(delta_eig)
 ABI_FREE(lambda)

 ! Compare ab-initio data with interpolated results.
 do spin=1,nsppol
   do ik=1,nkpt
     call skw_evalk(skw, cryst, kpts(:,ik), spin, oeig)
     write(std_out,*)eig(:,ik,spin) - oeig
     !call vdiff_print(vdiff_eval(1,nband,eig(:,ik,spin),oeig,one))
   end do
 end do

end subroutine skw_init
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_mkstar
!! NAME
!!  skw_mkstar
!!
!! FUNCTION
!!  Compute the star function at the given k-point kpt
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine skw_mkstar(skw, cryst, kpt, srk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_mkstar'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(out) :: srk(skw%nr)

!Local variables-------------------------------
!scalars
 integer :: ir,isym
 real(dp) :: arg

! *********************************************************************

 srk = zero
 do ir=1,skw%nr
   do isym=1,cryst%nsym
     arg = two_pi * dot_product(kpt, matmul(cryst%symrel(:,:,isym), skw%rpts(:,ir)))
     srk(ir) = srk(ir) + exp(j_dpc * arg)
   end do
   srk(ir) = srk(ir) / cryst%nsym
 end do

end subroutine skw_mkstar
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_evalk
!! NAME
!!  skw_evalk
!!
!! FUNCTION
!!  Interpolate the energies for a given k-point and spin.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine skw_evalk(skw, cryst, kpt, spin, oeig)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_evalk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 type(skw_t),intent(in) :: skw
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: oeig(skw%nband)

!Local variables-------------------------------
!scalars
 integer :: band
!arrays
 complex(dpc) :: srk(skw%nr)

! *********************************************************************

 ! Compute star function for this k-point
 call skw_mkstar(skw, cryst, kpt, srk)

 do band=1,skw%nband
   oeig(band) = dot_product(conjg(skw%coef(:,band,spin)), srk)
 end do

end subroutine skw_evalk
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_free
!! NAME
!!  skw_free
!!
!! FUNCTION
!!  Free memory
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine skw_free(skw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(inout) :: skw

! *********************************************************************

 if (allocated(skw%rpts)) then
   ABI_FREE(skw%rpts)
 end if
 if (allocated(skw%coef)) then
   ABI_FREE(skw%coef)
 end if

end subroutine skw_free
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_genrpts
!! NAME
!!  skw_genrpts
!!
!! FUNCTION
!!  Generate the list of real-space points.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine skw_genrpts(cryst, rmax, nr, rpts)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_genrpts'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: nr
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: rmax(3)
 real(dp),allocatable :: rpts(:,:)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: tolr=tol12
 integer :: i1,i2,i3,cnt,msize
!arrays
 !integer :: ngmax(3),ngmin(3)
 integer,allocatable :: iperm(:)
 real(dp),allocatable :: rtmp(:,:),r2(:)

! *********************************************************************

 msize = product(2*rmax+1)
 ABI_MALLOC(rtmp, (3, msize))
 ABI_MALLOC(r2, (msize))

 !ecut = ten * (2*pi**2)
 !do
 !  call kpgcount(ecut,0,gmet,[1],[zero,zero,zero],ngmax,ngmin,1)
 !  if ((2*ngmax + 1) > nr) exit
 !  ecut = two * ecut
 !end do

 cnt = 0
 do i3=-rmax(3),rmax(3)
   do i2=-rmax(2),rmax(2)
     do i1=-rmax(1),rmax(1)
       cnt = cnt + 1
       rtmp(:, cnt) = [i1,i2,i3]
       r2(cnt) = dot_product(rtmp(:,cnt), matmul(cryst%rmet, rtmp(:,cnt)))
     end do
   end do
 end do

 ! Sort norm
 ABI_MALLOC(iperm, (msize))
 iperm = [(i1, i1=1,msize)]
 call sort_dp(msize,r2,iperm,tolr)

 ! Find nr closing the shell
 do i1=nr+1,msize
   if (abs(r2(i1) - r2(nr)) > tolr) exit
 end do
 nr = (i1 - 1)

 ! Copy lattice points sorted by norm.
 ABI_MALLOC(rpts, (3, nr))
 do i1=1,nr
   rpts(:,i1) = rtmp(:,iperm(i1))
 end do

 ABI_FREE(iperm)
 ABI_FREE(rtmp)
 ABI_FREE(r2)

end subroutine skw_genrpts
!!***

!----------------------------------------------------------------------

END MODULE m_skw
!!***
