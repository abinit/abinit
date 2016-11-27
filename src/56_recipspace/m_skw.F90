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

 use m_fstrings,       only : itoa, sjoin, ktoa
 use m_numeric_tools,  only : vdiff_t, vdiff_eval, vdiff_print
 use m_bz_mesh,        only : isamek

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

  integer,allocatable :: rpts(:,:)
  ! rpts(3, nr)
  ! Real-space lattice points (reduced coordinates) ordered with non-decreasing length.

  complex(dpc),allocatable :: coef(:,:,:)
  ! coef(nr,nband,nsppol)

 end type skw_t

 public :: skw_new      ! Create new object.
 public :: skw_print    ! Print info about object.
 public :: skw_evalk    ! Evaluate the interpolated eigenvalues.
 public :: skw_free     ! Free memory.
!!***

CONTAINS  !=====================================================================================
!!***

!!****f* m_skw/skw_new
!! NAME
!!  skw_new
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
!!  comm=MPI communicator
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

type(skw_t) function skw_new(cryst, cplex, nband, nkpt, nsppol, kpts, eig, comm) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkpt,nsppol,cplex,comm
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpts(3,nkpt)
 real(dp),intent(in) :: eig(nband,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: my_rank,nprocs,cnt
 integer :: nk,ir,ik,ii,jj,nr,band,spin,isym,nsh,ierr
 real(dp) :: arg,ecut,r2,rmin2
 real(dp),parameter :: c1=0.25_dp,c2=0.25_dp
 character(len=500) :: fmt
!arrays
 integer,allocatable :: ipiv(:)
 real(dp) :: oeig(nband)
 real(dp),allocatable :: delta_eig(:,:,:),rhor(:)
 complex(dpc),allocatable :: srk(:,:),hmat(:,:),lambda(:,:,:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 new%nk = nkpt; new%nband = nband; new%nsppol = nsppol; new%cplex = cplex

 ! Select R-points and order them according to |R|
 call genrpts_(cryst, new%nk, new%nr, new%rpts)
 rmin2 = dot_product(new%rpts(:,2), matmul(cryst%rmet, new%rpts(:,2)))

 ! Construct star functions for the ab-initio k-points.
 nk = new%nk; nr = new%nr
 ABI_MALLOC(srk, (nr, nk))
 do ik=1,nkpt
   call mkstar_(new, cryst, kpts(:,ik), srk(:,ik))
 end do

 ! Compute roughness function.
 ABI_MALLOC(rhor, (nr))
 do ir=1,nr
   r2 = dot_product(new%rpts(:,ir), matmul(cryst%rmet, new%rpts(:,ir)))
   !rhor(ir) = (one - c1 * r2/rmin2)**2 + c2 * (r2 / rmin2)**3
   rhor(ir) = c1 * r2 + c2 * r2**2
 end do

 ! Build H matrix.
 ABI_CALLOC(hmat, (nk-1, nk-1))
 cnt = 0
 do jj=1,nk-1
   !do ii=1,jj
   do ii=1,nk-1
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi parallelism.
     do ir=2,nr
       hmat(ii, jj) = hmat(ii, jj) + &
         (srk(ir, ii) - srk(ir, nk)) * conjg(srk(ir, jj) - srk(ir, nk)) / rhor(ir)
     end do
   end do
 end do
 call xmpi_sum(hmat, comm, ierr)

 ABI_MALLOC(delta_eig, (nk-1,nband,nsppol))
 do spin=1,nsppol
   do band=1,nband
     delta_eig(:,band,spin) = eig(band,1:nk-1,spin) - eig(band,nk,spin)
   end do
 end do

 ! Solve system of linear equations to get lambda coeffients.
 ABI_MALLOC(lambda, (nk-1, nband, nsppol))
 ABI_MALLOC(ipiv, (nk-1))
 lambda = delta_eig

 ! Solve all the equations for all the bands at once (eq. 10 of PRB 38 2721)

 ! LU factorization of hmat matrix, kept in 2 triangular blocks of the matrix
 !call DGETRF(nk-1,nk-1,hmat,nk-1,ipiv,ierr)
 !CALL DGETRS('N',bs%nkpt-1,bs%icut2-bs%icut1+1,hmat,bs%nkpt-1,ipiv,lambda(1,bs%icut1),bs%nkpt-1,inf)

 call zgesv(nk-1, nband*nsppol, hmat, nk-1, ipiv, lambda, nk-1, ierr)
 !call zhesv(nk-1, nband*nsppol, hmat, nk-1, ipiv, lambda, nk-1, ierr)
 ABI_CHECK(ierr == 0, sjoin("ZGESV returned info:", itoa(ierr)))
 !lambda = conjg(lambda)

 ! Compute coefficients
 ABI_MALLOC(new%coef, (nr,nband,nsppol))

 do spin=1,nsppol
   do band=1,nband
     do ir=2,nr
       new%coef(ir,band,spin) = (one/rhor(ir)) * dot_product(srk(ir,:nk-1) - srk(ir,nk), lambda(:nk-1, band, spin))
       !new%coef(ir,band,spin) = (one/rhor(ir)) * dot_product(lambda(:nk-1, band, spin), conjg(srk(ir,:) - srk(ir,nk)))
       !new%coef(ir,band,spin) = (one/rhor(ir)) * dot_product(lambda(:nk-1, band, spin), conjg(srk(ir,:) - srk(ir,1)))
     end do
     new%coef(1,band,spin) = eig(band,nk,spin) - dot_product(conjg(new%coef(2:nr, band,spin)), srk(2:nr, nk))
   end do
 end do

 ABI_FREE(srk)
 ABI_FREE(rhor)
 ABI_FREE(hmat)
 ABI_FREE(delta_eig)
 ABI_FREE(lambda)
 ABI_FREE(ipiv)

 ! Compare ab-initio data with interpolated results.
 fmt = sjoin("(a,",itoa(nband),"(f8.3))")
 do spin=1,nsppol
   do ik=1,nkpt
     call skw_evalk(new, cryst, kpts(:,ik), spin, oeig)
     write(std_out,fmt)"-- ref ", eig(:,ik,spin) * Ha_eV
     write(std_out,fmt)"-- int ", oeig * Ha_eV
     write(std_out,*)"maxerr:", maxval(eig(:,ik,spin) - oeig) * Ha_eV, "[eV], ", trim(ktoa(kpts(:,ik)))
     call vdiff_print(vdiff_eval(1, nband, eig(:,ik,spin), oeig, one))
   end do
 end do

 call skw_print(new, std_out)

end function skw_new
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_print
!! NAME
!!  skw_print
!!
!! FUNCTION
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

subroutine skw_print(skw, unt)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
 integer,intent(in) :: unt

! *********************************************************************

 write(unt,"(a)")sjoin("Number of real-space lattice points:", itoa(skw%nr))
 write(unt,"(a)")sjoin("Number of ab-initio k-points:", itoa(skw%nk))
 write(unt,"(a)")sjoin("nband:", itoa(skw%nband), "nsppol", itoa(skw%nsppol), "cplex:", itoa(skw%cplex))

end subroutine skw_print
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_evalk
!! NAME
!!  skw_evalk
!!
!! FUNCTION
!!  Interpolate the energies for a given k-point and spin with slow FT.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  kpt(3)=K-point in reduced coordinates.
!!  spin=Spin index.
!!
!! OUTPUT
!!  oeig(skw%nband)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine skw_evalk(skw, cryst, kpt, spin, oeig)


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
 call mkstar_(skw, cryst, kpt, srk)

 do band=1,skw%nband
   oeig(band) = dot_product(conjg(skw%coef(:,band,spin)), srk)
 end do

 ! TODO: Derivatives

end subroutine skw_evalk
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_eval_kgrid
!! NAME
!!  skw_eval_kgrid
!!
!! FUNCTION
!!  Interpolate the energies on a grid with FFT.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  spin=Spin index.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine skw_eval_kgrid(skw, cryst, ngkpt, spin)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_kgrid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 type(skw_t),intent(in) :: skw
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ngkpt(3)

!Local variables-------------------------------
!scalars
 !integer :: band
!arrays
 !complex(dpc) :: srk(skw%nr)

! *********************************************************************

 ! Compute star function for this k-point
 !call mkstar_(skw, cryst, kpt, srk)

 !do band=1,skw%nband
 !  oeig(band) = dot_product(conjg(skw%coef(:,band,spin)), srk)
 !end do

 ! TODO: Derivatives

end subroutine skw_eval_kgrid
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

!!****f* m_skw/genrpts_
!! NAME
!!  genrpts_
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

subroutine genrpts_(cryst, nk, nr, rpts)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'genrpts_'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nk
 integer,intent(out) :: nr
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,allocatable,intent(out) :: rpts(:,:)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: tolr=tol12
 integer :: i1,i2,i3,cnt,msize
 real(dp) :: ratio
!arrays
 !integer :: ngmax(3),ngmin(3)
 integer :: rmax(3)
 integer,allocatable :: iperm(:)
 real(dp),allocatable :: rtmp(:,:),r2(:)

! *********************************************************************

 ratio = 120 !; nr = int(nkpt * ratio)
 !ecut = zero; nsh = 0; nr = 1
 !call setshells(ecut,nr,nsh,cryst%nsym,gmet,gprimd,symrel,tag,ucvol)
 !call kpgcount(ecut,exchn2n3d,gmet,istwfk,kpt,ngmax,ngmin,nkpt)
 !ecut = ten * (2*pi**2)
 !do
 !  call kpgcount(ecut,0,gmet,[1],[zero,zero,zero],ngmax,ngmin,1)
 !  if ((2*ngmax + 1) > nr) exit
 !  ecut = two * ecut
 !end do

 rmax = [20, 20, 20]
 msize = product(2*rmax + 1)
 ABI_MALLOC(rtmp, (3, msize))
 ABI_MALLOC(r2, (msize))

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

 ! Initial guess for nr
 do i1=1,msize
   if (i1 > nk * ratio) exit
 end do
 nr = i1 - 1

 ! Find nr closing the shell
 do i1=nr+1,msize
   if (abs(r2(i1) - r2(nr)) > tolr) exit
 end do
 nr = i1 - 1

 ! Copy lattice points sorted by norm.
 ABI_MALLOC(rpts, (3, nr))
 do i1=1,nr
   rpts(:,i1) = rtmp(:,iperm(i1))
 end do

 ABI_FREE(iperm)
 ABI_FREE(rtmp)
 ABI_FREE(r2)

end subroutine genrpts_
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar_
!! NAME
!!  mkstar_
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

subroutine mkstar_(skw, cryst, kpt, srk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkstar_'
 use interfaces_32_util
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
 integer :: ik,ir,isym,my_nsym,nkstar
 real(dp) :: arg
 logical :: found
!arrays
 integer :: g0(3),mit(3,3)
 real(dp) :: kstar(3,cryst%nsym),new_sk(3)

! *********************************************************************

 srk = zero

#if 1
 nkstar = 1; kstar(:, 1) = kpt
 do isym=1,cryst%nsym
   if (cryst%symafm(isym) == -1) cycle
   !call mati3inv(cryst%symrel(:,:,isym), mit)
   !new_sk = matmul(transpose(mit), kpt)
   new_sk = matmul(transpose(cryst%symrel(:,:,isym)), kpt)
   do ik=1,nkstar
     found = isamek(new_sk, kstar(:,ik), g0)
     if (any(g0 /= 0)) found = .False.
     if (found) exit
   end do
   if (.not. found) then
      nkstar = nkstar + 1
      kstar(:, nkstar) = new_sk
   end if
 end do

 do ir=1,skw%nr
   do ik=1,nkstar
     arg = two_pi * dot_product(kstar(:,ik), skw%rpts(:,ir))
     srk(ir) = srk(ir) + exp(j_dpc * arg)
   end do
   srk(ir) = srk(ir) / nkstar
   !srk(ir) = srk(ir) / cryst%nsym
 end do

#else
 my_nsym = count(cryst%symafm == 1)
 do ir=1,skw%nr
   do isym=1,cryst%nsym
     if (cryst%symafm(isym) == -1) cycle
     arg = two_pi * dot_product(kpt, matmul(cryst%symrel(:,:,isym), skw%rpts(:,ir)))
     srk(ir) = srk(ir) + exp(j_dpc * arg)
   end do
   srk(ir) = srk(ir) / my_nsym
 end do
#endif

end subroutine mkstar_
!!***

!----------------------------------------------------------------------

end module m_skw
!!***
