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
!!  Object implementing the Shankland-Koelling-Wood Fourier interpolation scheme.
!!  It can be used to interpolate functions in k-space with the periodicity of the
!!  reciprocal lattice and satisfying F(k) = F(Sk) for each rotation S
!!  belonging to the point group of the crystal. For readability reason,
!!  the names of the variables are chosen assuming we are interpolating electronic eigenvalues
!!  but the same object can be use to interpolate phonons as well. Just use nsppol=1 and nband = 3 * natom
!!
!! SOURCE

 !type :: skwcoefs_t
 !  complex(dpc),allocatable :: vals(:,:,:)
 !end type skwcoefs_t

 type,public :: skw_t

  integer :: cplex
   ! 1 if time-reversal symmetry can be used, 2 otherwise.

  integer :: nr
   ! Number of real-space lattice points.

  integer :: nkpt
   ! Number of ab-initio k-points

  integer :: band_block(2)
   ! Initial and final band index treated by this processor

  integer :: spin_block(2)
   ! Initial and final spin index treated by this processor

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer,allocatable :: rpts(:,:)
   ! rpts(3, nr)
   ! Real-space lattice points (reduced coordinates) ordered with non-decreasing length.

  complex(dpc),allocatable :: coefs(:,:,:)
   ! coefs(nr, mband, nsppol)

   !type(skcoefs_t),allocatable :: coefs(:,:)
   ! coefs(mband, nsppol)

  complex(dpc),allocatable :: cached_srk(:)
   ! cached_srk(skw%nr)
   ! The star function for cached_kpt (used in skw_eval_bks)
  real(dp) :: cached_kpt(3)

  complex(dpc),allocatable :: cached_srk_dk1(:,:)
   ! cached_srk_dk1(skw%nr, 3)
   ! The 1d derivative wrt k of the star function for cached_kpt_dk1 (used in skw_eval_bks)
  real(dp),private :: cached_kpt_dk1(3)

  complex(dpc),allocatable :: cached_srk_dk2(:,:,:)
   ! cached_srk_dk2(skw%nr,3,3)
   ! The 2d derivatives wrt k of the star function for cached_kpt_dk2 (used in skw_eval_bks)
  real(dp) :: cached_kpt_dk2(3)

 end type skw_t

 public :: skw_new          ! Create new object.
 public :: skw_print        ! Print info about object.
 public :: skw_eval_bks     ! Interpolate eigenvalues, 1st, 2nd derivates wrt k, at an arbitrary k-point.
 public :: skw_free         ! Free memory.
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
!!  nband=Total Number of bands in eig array.
!!  nkpt=Number of ab-initio k-points.
!!  nsppol=Number of independent spin polarizations.
!!  kpts(3,nkpt)=ab-initio k-points in reduced coordinates.
!!  eig(nband,nkpt,nsppol)=ab-initio eigenvalues.
!!  band_block(2)=Initial and final band index. If [0,0], all bands are used
!!  spin_block(2)=Initial and final spin index. If [0,0], all bands are used
!!  comm=MPI communicator
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

type(skw_t) &
function skw_new(cryst, cplex, nband, nkpt, nsppol, kpts, eig, band_block, spin_block, comm) result(new)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nband,nkpt,nsppol,comm
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: band_block(2),spin_block(2)
 real(dp),intent(in) :: kpts(3,nkpt)
 real(dp),intent(in) :: eig(nband,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: my_rank,nprocs,cnt,bstop,bstart,bcount
 integer :: ir,ik,ib,ii,jj,nr,band,spin,isym,nsh,ierr,i1,i2,i3,msize
 real(dp),parameter :: tolr=tol12
 real(dp) :: arg,ecut,r2,r2min,mare,mae,adiff_meV,rel_err,ratio
 real(dp),parameter :: c1=0.25_dp,c2=0.25_dp
 character(len=500) :: fmt
!arrays
 integer :: rmax(3) !,ngmax(3),ngmin(3)
 integer,allocatable :: ipiv(:),iperm(:)
 real(dp),allocatable :: rtmp(:,:),r2vals(:)
 real(dp),allocatable :: delta_eig(:,:,:),rhor(:),oeig(:)
 complex(dpc),allocatable :: srk(:,:),hmat(:,:),lambda(:,:,:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Get slice of bands to be treated.
 new%band_block = band_block; if (all(band_block == 0)) new%band_block = [1, nband]
 new%spin_block = spin_block; if (all(spin_block == 0)) new%spin_block = [1, nsppol]
 call xmpi_min(new%band_block(1), bstart, comm, ierr)
 call xmpi_max(new%band_block(2), bstop, comm, ierr)
 bcount = bstop - bstart + 1

 new%cplex = cplex; new%nkpt = nkpt; new%nsppol = nsppol

 !call crystal_point_group(cryst, new%ptg_nsym, new%ptg_symrel, new%ptg_symrec, has_inversion)

 ! -----------------------------------------------
 ! Select R-points and order them according to |R|
 ! -----------------------------------------------

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

 rmax = [20, 20, 20]; msize = product(2*rmax + 1)
 ABI_MALLOC(rtmp, (3, msize))
 ABI_MALLOC(r2vals, (msize))

 cnt = 0
 do i3=-rmax(3),rmax(3)
   do i2=-rmax(2),rmax(2)
     do i1=-rmax(1),rmax(1)
       cnt = cnt + 1
       rtmp(:, cnt) = [i1,i2,i3]
       r2vals(cnt) = dot_product(rtmp(:,cnt), matmul(cryst%rmet, rtmp(:,cnt)))
     end do
   end do
 end do

 ! Sort norm
 ABI_MALLOC(iperm, (msize))
 iperm = [(i1, i1=1,msize)]
 call sort_dp(msize, r2vals, iperm, tolr)

 ! Initial guess for nr
 do i1=1,msize
   if (i1 > nkpt * ratio) exit
 end do
 nr = i1 - 1

 ! Find nr closing the shell
 do i1=nr+1,msize
   if (abs(r2vals(i1) - r2vals(nr)) > tolr) exit
 end do
 nr = i1 - 1

 ! Copy lattice points sorted by norm.
 new%nr = nr
 ABI_MALLOC(new%rpts, (3, nr))
 do i1=1,nr
   new%rpts(:,i1) = rtmp(:,iperm(i1))
 end do
 r2min = dot_product(new%rpts(:,2), matmul(cryst%rmet, new%rpts(:,2)))

 ABI_FREE(iperm)
 ABI_FREE(rtmp)
 ABI_FREE(r2vals)

 ! Construct star functions for the ab-initio k-points.
 ABI_MALLOC(srk, (nr, nkpt))
 do ik=1,nkpt
   call mkstar(new, cryst, kpts(:,ik), srk(:,ik))
 end do

 ! Compute roughness function.
 ABI_MALLOC(rhor, (nr))
 do ir=1,nr
   r2 = dot_product(new%rpts(:,ir), matmul(cryst%rmet, new%rpts(:,ir)))
   !rhor(ir) = (one - c1 * r2/r2min)**2 + c2 * (r2 / r2min)**3
   rhor(ir) = c1 * r2 + c2 * r2**2
 end do

 ! Build H(k,k') matrix.
 ABI_CALLOC(hmat, (nkpt-1, nkpt-1))
 cnt = 0
 do jj=1,nkpt-1
   !do ii=1,jj
   do ii=1,nkpt-1
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi parallelism.
     do ir=2,nr
       hmat(ii, jj) = hmat(ii, jj) + &
         (srk(ir, ii) - srk(ir, nkpt)) * conjg(srk(ir, jj) - srk(ir, nkpt)) / rhor(ir)
     end do
   end do
 end do
 call xmpi_sum(hmat, comm, ierr)

 ABI_MALLOC(delta_eig, (nkpt-1, bcount, nsppol))
 do spin=1,nsppol
   do ib=1,bcount
     band = ib + bstart - 1
     delta_eig(:,ib,spin) = eig(band,1:nkpt-1,spin) - eig(band,nkpt,spin)
   end do
 end do

 ! Solve system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)
 ! Solve all bands and spins at once
 ABI_MALLOC(lambda, (nkpt-1, bcount, nsppol))
 ABI_MALLOC(ipiv, (nkpt-1))
 lambda = delta_eig

 !call DGETRF(nkpt-1,nkpt-1,hmat,nkpt-1,ipiv,ierr)
 !CALL DGETRS('N',bs%nkpt-1,bs%icut2-bs%icut1+1,hmat,bs%nkpt-1,ipiv,lambda(1,bs%icut1),bs%nkpt-1,inf)

 call zgesv(nkpt-1, bcount*nsppol, hmat, nkpt-1, ipiv, lambda, nkpt-1, ierr)
 !call zhesv(nkpt-1, bcount*nsppol, hmat, nkpt-1, ipiv, lambda, nkpt-1, ierr)
 ABI_CHECK(ierr == 0, sjoin("ZGESV returned info:", itoa(ierr)))

 ! Compute coefficients
 ABI_MALLOC(new%coefs, (nr,bcount,nsppol))

 do spin=1,nsppol
   do ib=1,bcount
     band = ib + bstart - 1
     do ir=2,nr
       new%coefs(ir,ib,spin) = (one/rhor(ir)) * dot_product(srk(ir,:nkpt-1) - srk(ir,nkpt), lambda(:nkpt-1, ib, spin))
       !new%coefs(ir,ib,spin) = (one/rhor(ir)) * dot_product(lambda(:nkpt-1, ib, spin), conjg(srk(ir,:) - srk(ir,nkpt)))
       !new%coefs(ir,ib,spin) = (one/rhor(ir)) * dot_product(lambda(:nkpt-1, ib, spin), conjg(srk(ir,:) - srk(ir,1)))
     end do
     new%coefs(1,ib,spin) = eig(band,nkpt,spin) - dot_product(conjg(new%coefs(2:nr, ib,spin)), srk(2:nr, nkpt))
   end do
 end do

 ! Prepare workspace arrays for star functions.
 new%cached_kpt = huge(one)
 ABI_MALLOC(new%cached_srk, (new%nr))
 new%cached_kpt_dk1 = huge(one)
 ABI_MALLOC(new%cached_srk_dk1, (new%nr, 3))
 new%cached_kpt_dk2 = huge(one)
 ABI_MALLOC(new%cached_srk_dk2, (new%nr, 3, 3))

 ABI_FREE(srk)
 ABI_FREE(rhor)
 ABI_FREE(hmat)
 ABI_FREE(delta_eig)
 ABI_FREE(lambda)
 ABI_FREE(ipiv)

 ! Compare ab-initio data with interpolated results.
 ABI_MALLOC(oeig, (bcount))
 fmt = sjoin("(a,",itoa(bcount),"(f8.3))")
 bstop = bstart + bcount - 1
 mare = zero; mae = zero
 do spin=1,nsppol
   do ik=1,nkpt
     do ib=1,bcount
       band = ib + new%band_block(1) - 1
       call skw_eval_bks(new, cryst, band, kpts(:,ik), spin, oeig(ib))

       adiff_meV = abs(eig(band,ik,spin) - oeig(ib)); rel_err = zero
       if (abs(eig(band,ik,spin)) > tol16) rel_err = adiff_meV / abs(eig(band,ik,spin))
       rel_err = 100 * rel_err; adiff_meV = adiff_meV * Ha_meV
       mae = mae + adiff_meV; mare = mare + rel_err
     end do

     write(std_out,fmt)"-- ref ", eig(bstart:bstop,ik,spin) * Ha_meV
     write(std_out,fmt)"-- int ", oeig * Ha_meV
     write(std_out,*)"maxerr:", maxval(eig(bstart:bstop,ik,spin) - oeig) * Ha_meV, " [meV], ", trim(ktoa(kpts(:,ik)))
     !call vdiff_print(vdiff_eval(1, bcount, eig(bstart:bstop,ik,spin), oeig, one))
   end do
 end do

 !list2 = [mare/cnt, mae/cnt]
 !call xmpi_sum(list2, comm, ierr)

 call skw_print(new, std_out)

 ABI_FREE(oeig)

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
 write(unt,"(a)")sjoin("Number of ab-initio k-points:", itoa(skw%nkpt))
 write(unt,"(a)")sjoin("nsppol", itoa(skw%nsppol), "cplex:", itoa(skw%cplex))

end subroutine skw_print
!!***

!----------------------------------------------------------------------

!!****f* m_skw/skw_eval_bks
!! NAME
!!  skw_eval_bks
!!
!! FUNCTION
!!  Interpolate the energies for an arbitrary k-point and spin with slow FT.
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  band=Band index.
!!  kpt(3)=K-point in reduced coordinates.
!!  spin=Spin index.
!!
!! OUTPUT
!!  oeig=interpolated eigenvalues
!!    Note that oeig is not necessarily sorted in ascending order.
!!    The routine does not reorder the interpolated eigenvalues
!!    to be consistent with the interpolation of the derivatives.
!!  [oder1(3)]=First-order derivatives wrt k in reduced coordinates.
!!  [oder2(3,3)]=Second-order derivatives wrt k in reduced coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine skw_eval_bks(skw, cryst, band, kpt, spin, oeig, oder1, oder2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skw_eval_bks'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,spin
 type(skw_t),intent(inout) :: skw
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: oeig
 real(dp),optional,intent(out) :: oder1(3)
 real(dp),optional,intent(out) :: oder2(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
! *********************************************************************

#ifdef DEBUG_MODE
 !ABI_CHECK(allocated(ebspl%coefs(band, spin)%vals), sjoin("Unallocated (band, spin):", ltoa([band, spin])))
 !ABI_CHECK(ib >= 1 .and. ib <= skw%bcount, sjoin("out of range band:", itoa(band)))
#endif

 ! Compute star function for this k-point (if not already in memory)
 if (any(kpt /= skw%cached_kpt)) then
   call mkstar(skw, cryst, kpt, skw%cached_srk)
   skw%cached_kpt = kpt
 end if

 oeig = dot_product(conjg(skw%coefs(:,band,spin)), skw%cached_srk)

 ! TODO: Finalize Derivatives
 if (present(oder1)) then
   ! Compute first-order derivatives.
   if (any(kpt /= skw%cached_kpt_dk1)) then
     call mkstar_dk1(skw, cryst, kpt, skw%cached_srk_dk1)
     skw%cached_kpt_dk1 = kpt
   end if
   do ii=1,3
     oder1(ii) = dot_product(conjg(skw%coefs(:,band,spin)), skw%cached_srk_dk1(:,ii))
   end do
 end if

 if (present(oder2)) then
   ! Compute second-order derivatives.
   if (any(kpt /= skw%cached_kpt_dk2)) then
     call mkstar_dk2(skw, cryst, kpt, skw%cached_srk_dk2)
     skw%cached_kpt_dk2 = kpt
   end if

   oder2 = zero
   do jj=1,3
     do ii=1,jj
       oder2(ii, jj) = dot_product(conjg(skw%coefs(:,band,spin)), skw%cached_srk_dk2(:,ii,jj))
       if (ii /= jj) oder2(jj, ii) = oder2(ii, jj)
     end do
   end do
 end if

end subroutine skw_eval_bks
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
 if (allocated(skw%coefs)) then
   ABI_FREE(skw%coefs)
 end if

 if (allocated(skw%cached_srk)) then
   ABI_FREE(skw%cached_srk)
 end if
 skw%cached_kpt = huge(one)

 if (allocated(skw%cached_srk_dk1)) then
   ABI_FREE(skw%cached_srk_dk1)
 end if
 skw%cached_kpt_dk1 = huge(one)

 if (allocated(skw%cached_srk_dk2)) then
   ABI_FREE(skw%cached_srk_dk2)
 end if
 skw%cached_kpt_dk2 = huge(one)

end subroutine skw_free
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar
!! NAME
!!  mkstar
!!
!! FUNCTION
!!  Compute the star function for k-point kpt
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  kpt(3)=K-point in reduced coordinates.
!!
!! OUTPUT
!!  srk(%nr)=Star function for this k-point.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkstar(skw, cryst, kpt, srk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkstar'
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

end subroutine mkstar
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar_dk1
!! NAME
!!  mkstar_dk1
!!
!! FUNCTION
!!  Compute the 1st derivative of the star function wrt k
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  kpt(3)=K-point in reduced coordinates.
!!
!! OUTPUT
!!  srk_dk1(%nr,3)=Derivative of the star function wrt k
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkstar_dk1(skw, cryst, kpt, srk_dk1)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkstar_dk1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(out) :: srk_dk1(skw%nr,3)

!Local variables-------------------------------
!scalars
 !integer :: ik,ir,isym,my_nsym,nkstar
!arrays
 !integer :: g0(3),mit(3,3)
 !real(dp) :: kstar(3,cryst%nsym),new_sk(3)

! *********************************************************************

 srk_dk1 = zero

end subroutine mkstar_dk1
!!***

!----------------------------------------------------------------------

!!****f* m_skw/mkstar_dk2
!! NAME
!!  mkstar_dk2
!!
!! FUNCTION
!!  Compute the 2st derivatives of the star function wrt k
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure.
!!  kpt(3)=K-point in reduced coordinates.
!!
!! OUTPUT
!!  srk_dk2(%nr,3,3)=2nd derivatives of the star function wrt k
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkstar_dk2(skw, cryst, kpt, srk_dk2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkstar_dk2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(skw_t),intent(in) :: skw
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(out) :: srk_dk2(skw%nr,3,3)

!Local variables-------------------------------
!scalars
 !integer :: ik,ir,isym,my_nsym,nkstar
!arrays
 !integer :: g0(3),mit(3,3)
 !real(dp) :: kstar(3,cryst%nsym),new_sk(3)

! *********************************************************************

 srk_dk2 = zero

end subroutine mkstar_dk2
!!***

!----------------------------------------------------------------------

end module m_skw
!!***
