!!****m* ABINIT/m_krank
!! NAME
!! m_krank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!!
!! COPYRIGHT
!! Copyright (C) 2010-2025 ABINIT group (MVer, HM, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! TODO: Remove file
!#include "libtetra.h"

module m_krank

 use defs_basis
 use m_abicore
 use m_errors

 use m_fstrings, only : itoa, sjoin

 implicit none

 private
!!***

!!****t* m_krank/krank_t
!! NAME
!! krank_t
!!
!! FUNCTION
!!  structure to contain a rank/inverse rank pair of arrays, with dimensions
!!
!! SOURCE

 type,public :: krank_t

   integer :: max_linear_density

   integer :: min_rank

   integer :: max_rank

   integer :: npoints

   logical :: time_reversal

   logical :: kpts_owns_memory = .False.

   integer,allocatable :: invrank(:)

   real(dp),ABI_CONTIGUOUS pointer :: kpts(:,:)
    ! Reference to input k-points or copy of the array depending on kpts_owns_memory

    ! Internal tables used by krank_get_mapping
   integer,allocatable :: rank2ikpt_(:), rank2symtime_(:)

 contains

   procedure :: get_rank
    ! Calculates the rank for one kpt

   procedure :: get_index => krank_get_index
    ! Return the index of the k-point `kpt` in the initial set. -1 if not found.

   procedure :: copy => krank_copy
    ! Deep copy of the object

   procedure :: free => krank_free
    ! Free memory

   procedure :: print => krank_print
    ! Prints the arrays and dimensions of a krank_t structure

   procedure :: get_mapping => krank_get_mapping
    ! Use symmetries to map input kptn2 to the list of k-points used to generate krank_t.
    ! Similar to listkk but, unlike listkk, this algo does not try to minimize the distance.
    ! Must faster than listkk for dense meshes

 end type krank_t

 public :: krank_from_kptrlatt  ! Initialize object from kptrlatt
 public :: krank_new            ! Sets up the kpt ranks for comparing kpts

 public :: get_ibz2bz           ! Return array with the index of the IBZ wave vectors in the BZ.
 public :: star_from_ibz_idx    ! Return array with the indices of the star of the ik_ibz wavevector in the IBZ.
!!***

contains


!!****f* m_krank/krank_from_kptrlatt
!! NAME
!! krank_from_kptrlatt
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! INPUTS
!!  npt = number of kpoints (eventually irreducible)
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!  krank = object containing ranking and inverse ranking
!!
!! NOTES
!!  By default, the object holds a reference to kpts so do not change/deallocate this array
!!  while using krank.
!!
!! SOURCE

type(krank_t) function krank_from_kptrlatt(nkpt, kpts, kptrlatt, compute_invrank) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 logical,optional,intent(in) :: compute_invrank
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),target,intent(in) :: kpts(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: ii, jj, ikpt, max_linear_density, opt=0
 logical :: compute_invrank_
 real(dp) :: min_kpt

! *********************************************************************

 opt=0
 do jj=1,3
   do ii=1,3
     if (ii == jj .and. kptrlatt(ii, ii) == 0) then
       ABI_ERROR("kptrlatt with zero matrix element on the diagonal!")
     end if
     if (ii /= jj .and. kptrlatt(ii, jj) /= 0) then
       ! ABI_WARNING("kptrlatt with non-zero off-diagonal matrix elements is not supported")
       opt=1
     end if
   end do
 end do

 compute_invrank_ = .True.; if (present(compute_invrank)) compute_invrank_ = compute_invrank

 min_kpt = 1
 if (opt == 1) then
   do ikpt=1,nkpt
     do ii=1,3
       if (abs(kpts(ii,ikpt)) < min_kpt .and. abs(kpts(ii,ikpt)) /= 0) then
         min_kpt = abs(kpts(ii,ikpt)) ! used as tmp variable
       end if
     end do
   end do
   max_linear_density = ceiling(2/min_kpt)
 else
   max_linear_density = maxval([kptrlatt(1,1), kptrlatt(2,2), kptrlatt(3,3)])
 end if

 new = krank_new(nkpt, kpts, max_linear_density=max_linear_density, compute_invrank=compute_invrank_)

end function krank_from_kptrlatt
!!***

!!****f* m_krank/krank_new
!! NAME
!! krank_new
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! INPUTS
!!  npt = number of kpoints (eventually irreducible)
!!  kpt = coordinates of kpoints
!!  time_reversal = true or false to use time reversal symmetry.
!!     Default is true, but only important if nsym and symrec are present
!!
!! SOURCE

type(krank_t) function krank_new(nkpt, kpts, &
                                 nsym, symrec, time_reversal, max_linear_density, compute_invrank) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 integer,intent(in), optional :: nsym
 logical,intent(in), optional :: time_reversal
 integer,optional,intent(in) :: max_linear_density
 logical,optional,intent(in) :: compute_invrank
!arrays
 real(dp),target,intent(in) :: kpts(3,nkpt)
 integer,intent(in), optional :: symrec(3,3, *)

!Local variables -------------------------
!scalars
 integer :: ikpt, isym, symkptrank, irank, timrev, itim
 logical :: compute_invrank_
 real(dp) :: smallestlen
 character(len=500) :: msg
!arrays
 real(dp) :: symkpt(3)

! *********************************************************************

 compute_invrank_ = .True.; if (present(compute_invrank)) compute_invrank_ = compute_invrank

 new%kpts => kpts
 new%kpts_owns_memory = .False.

 if (.not. present(max_linear_density)) then
   ! Find smallest linear length from input kpts
   smallestlen = one
   do ikpt=1, nkpt
     if (abs(kpts(1,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpts(1,ikpt)))
     if (abs(kpts(2,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpts(2,ikpt)))
     if (abs(kpts(3,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpts(3,ikpt)))
   end do
   new%max_linear_density = nint(one/smallestlen)
 else
   ! Get it from input
   new%max_linear_density = max_linear_density
 end if

 new%npoints = nkpt
 new%min_rank = nint(real(new%max_linear_density)*(half+tol8 +&
                     real(new%max_linear_density)*(half+tol8 +&
                     real(new%max_linear_density)*(half+tol8))))

 new%max_rank = nint(real(new%max_linear_density)*(1+half+tol8 +&
                     real(new%max_linear_density)*(1+half+tol8 +&
                     real(new%max_linear_density)*(1+half+tol8))))

 timrev = 2
 new%time_reversal = .true.
 if (present(time_reversal)) then
   if (.not. time_reversal) timrev = 1
   new%time_reversal = .false.
 end if

 ! Ensure kpt(i)+one is positive, and the smallest difference between kpts should be larger than 1/100 ie ngkpt < 100.
 ! the following fills invrank for the k-points in the list provided (may be only the irred kpts)
 if (compute_invrank_) then
   ABI_MALLOC(new%invrank, (new%min_rank:new%max_rank))
   new%invrank(:) = -1

   do ikpt=1,nkpt
     irank = new%get_rank(kpts(:,ikpt))
     if (irank > new%max_rank .or. irank < new%min_rank) then
       write(msg,'(a,2i0)')" rank above max_rank or below min_rank, ikpt, rank ", ikpt, irank
       ABI_ERROR(msg)
     end if
     new%invrank(irank) = ikpt
   end do
 end if

 ! if symrec is provided, fill invrank with appropriate irred kpt indices
 ! for symmetry completion: kptrank_t%invrank points to the irred k-point
 ! equivalent to the k-point whose rank is provided
 if (present(symrec)) then
   if(.not. present(nsym)) then
     ABI_ERROR("need both symrec and nsym arguments together")
   end if
   do ikpt=1,nkpt
     ! itim == 1 for positive, and itim==2 gives Kramers opposite of k-point
     ! favor the former by looping it last
     do itim = timrev, 1, -1
       do isym = 1, nsym
         symkpt = (-1)**(itim+1) * matmul(symrec(:,:,isym), kpts(:, ikpt))
         symkptrank = new%get_rank(symkpt(:))
         new%invrank(symkptrank) = ikpt
       end do
     end do
   end do
 end if

end function krank_new
!!***

!----------------------------------------------------------------------

!!****f* m_krank/get_rank
!! NAME
!! get_rank
!!
!! FUNCTION
!! This routine calculates the rank for one kpt
!!
!! INPUTS
!!  kpt = coordinates of kpoints
!!  krank = rank object for the k-grid we are using
!!
!! OUTPUT
!!  rank = rank of the kpoint
!!
!! SOURCE

integer function get_rank(krank, kpt) result(rank)

!Arguments ------------------------------------
!scalars
 class(krank_t), intent(in) :: krank
!arrays
 real(dp),intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays
 real(dp) :: redkpt(3)
! *************************************************************************

 ! wrap to [0, 1[ -> replaced call to wrap2_zero2one inline, to encapsulate this module
 if (kpt(1)>zero) then
   redkpt(1)=mod((kpt(1)+tol12),one)-tol12
 else
   redkpt(1)=-mod(-(kpt(1)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(1))<tol12)redkpt(1)=zero

 if (kpt(2)>zero) then
   redkpt(2)=mod((kpt(2)+tol12),one)-tol12
 else
   redkpt(2)=-mod(-(kpt(2)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(2))<tol12)redkpt(2)=zero

 if (kpt(3)>zero) then
   redkpt(3)=mod((kpt(3)+tol12),one)-tol12
 else
   redkpt(3)=-mod(-(kpt(3)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(3))<tol12)redkpt(3)=zero

! rank = int(real(krank%max_linear_density)*(redkpt(3)+half+tol8 + &
!            real(krank%max_linear_density)*(redkpt(2)+half+tol8 + &
!            real(krank%max_linear_density)*(redkpt(1)+half+tol8))))
 rank = nint(real(krank%max_linear_density)*(redkpt(1)+half+tol8 + &
             real(krank%max_linear_density)*(redkpt(2)+half+tol8 + &
             real(krank%max_linear_density)*(redkpt(3)+half+tol8))))

 if (rank > krank%max_rank) then
   write(msg,'(2(a,i0))') ' Rank should be <= max_rank: ', krank%max_rank, ' but got: ', rank
   ABI_ERROR(msg)
 end if
 if (rank < krank%min_rank) then
   !print *, "redkpt", redkpt
   write(msg,'(2(a,i0))') ' Rank should be >= min_rank ', krank%min_rank, ' but got: ', rank
   ABI_ERROR(msg)
 end if

end function get_rank
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_get_index
!! NAME
!! krank_get_index
!!
!! FUNCTION
!!  Return the index of the k-point `kpt` in the initial set. -1 if not found.
!!
!! INPUTS
!!  krank = rank object for the k-grid we are using
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!
!! SOURCE

integer function krank_get_index(krank, kpt) result(ikpt)

!Arguments ------------------------------------
!scalars
 class(krank_t), intent(in) :: krank
!arrays
 real(dp),intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 integer :: kpt_rank

! *************************************************************************

 kpt_rank = krank%get_rank(kpt)
 ikpt = -1
 if (kpt_rank < krank%max_rank) ikpt = krank%invrank(kpt_rank)

end function krank_get_index
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_copy
!!
!! NAME
!! krank_copy
!!
!! FUNCTION
!! Deep copy of the object
!!
!! INPUTS
!!
!! OUTPUT
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! SOURCE

type(krank_t) function krank_copy(krank_in) result(krank_out)

!Arguments ------------------------------------
!scalars
 class(krank_t), intent(in) :: krank_in

! *********************************************************************

 krank_out%max_linear_density = krank_in%max_linear_density
 krank_out%min_rank = krank_in%min_rank
 krank_out%max_rank = krank_in%max_rank
 krank_out%npoints = krank_in%npoints

 ABI_MALLOC(krank_out%invrank, (krank_out%min_rank:krank_out%max_rank))
 krank_out%invrank = krank_in%invrank

 ! This is why I call it deep copy!
 ABI_MALLOC(krank_out%kpts, (3, size(krank_in%kpts, dim=2)))
 krank_out%kpts = krank_in%kpts
 krank_out%kpts_owns_memory = .True.

end function krank_copy
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_free
!! NAME
!! krank_free
!!
!! FUNCTION
!! This routine deallocates the arrays in a krank_t structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! SOURCE

subroutine krank_free(krank)

!Arguments ------------------------------------
 class(krank_t), intent(inout) :: krank

! *********************************************************************

 ABI_SFREE(krank%invrank)
 ABI_SFREE(krank%rank2ikpt_)
 ABI_SFREE(krank%rank2symtime_)

 if (krank%kpts_owns_memory) then
   ABI_SFREE_PTR(krank%kpts)
   krank%kpts_owns_memory = .False.
   krank%kpts => null()
 else
   krank%kpts => null()
 end if

end subroutine krank_free
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_print
!!
!! NAME
!! krank_print
!!
!! FUNCTION
!! This routine prints the arrays and dimensions of a krank_t structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking
!!  unout = unit for open file to print to
!!
!! SOURCE

subroutine krank_print(krank, unout)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: unout
!arrays
 class(krank_t), intent(in) :: krank

! *********************************************************************

  write(unout, *)
  write(unout, '(a)') ' Dump of the contents of a krank_t structure with k-point rank information'
  write(unout, '(a,i0)') ' max linear density of points in 3 directions: max_linear_density = ',  krank%max_linear_density
  write(unout, '(a,i0)') ' maximum rank for any point in grid: max_rank = ',  krank%max_rank
  write(unout, '(a,i0)') ' number of points in input grid: npoints = ',  krank%npoints
  write(unout, *)
  write(unout, '(a)') ' invrank array = '
  write(unout, '(i0)') krank%invrank(:)
  write(unout, *)

end subroutine krank_print
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_get_mapping
!! NAME
!! krank_get_mapping
!!
!! FUNCTION
!! Use symmetries to map input kptn2 to the list of k-points used to generate krank_t.
!! Similar to listkk but, unlike listkk, this algo does not try to minimize the distance
!! Mainly used to map two set of k-points associated to the same grid (e.g. BZ --> IBZ, IBZ(q) --> IBZ etc.
!! Must be faster than listkk for dense meshes
!! although this routine requires the allocation of temporary array of shape (2, self%min_rank:self%max_rank)
!! Returns indirect indexing list indkk.
!!
!! INPUTS
!!  kptns2(3,nkpt2)=list of final k points
!!  nkpt2=number of final k points
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  nsym=number of symmetry elements
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symmat(3,3,nsym)=symmetry operations (symrel or symrec, depending on value of use_symrec)
!!  timrev=1 if the use of time-reversal is allowed; 0 otherwise
!!  [use_symrec]: if present and true, symmat assumed to be symrec, otherwise assumed to be symrel (default)
!!  [qpt]: q-point to be added to kptns2. 0 if not specified.
!!
!! OUTPUT
!!  dksqmax=maximal value of the norm**2 of the difference between
!!    a kpt2 vector and the closest k-point found from the kptns1 set, using symmetries.
!!  indkk(6, nkpt2)=
!!    indkk(:,1)=k point index of kptns1
!!    indkk(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!    indkk(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!      to give kpt1b, that is the closest to kpt2.
!!    indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!
!! SOURCE

subroutine krank_get_mapping(self, nkpt2, kptns2, dksqmax, gmet, indkk, nsym, symafm, symmat, timrev, &
                             use_symrec, qpt) ! optional

!Arguments ------------------------------------
!scalars
 class(krank_t),intent(inout) :: self
 integer,intent(in) :: nkpt2, nsym, timrev
 real(dp),intent(out) :: dksqmax
 logical,optional,intent(in) :: use_symrec
!arrays
 integer,intent(in) :: symafm(nsym), symmat(3,3,nsym)
 integer,intent(out) :: indkk(6, nkpt2)
 real(dp),intent(in) :: gmet(3,3), kptns2(3,nkpt2)
 real(dp),optional,intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer :: irank, ikpt1, ikpt2, itimrev, isym, ii
 logical :: my_use_symrec
!arrays
 integer :: dkint(3), my_symmat(3, 3, nsym)
 !integer,allocatable :: rank2ikpt(:), rank2symtime(:)
 real(dp) :: kpt1a(3), dk(3), my_qpt(3)

! *************************************************************************

 my_use_symrec = .False.; if (present(use_symrec)) my_use_symrec = use_symrec
 my_qpt = zero; if (present(qpt)) my_qpt = qpt

 if (my_use_symrec) then
   ! Symrec k
   my_symmat = symmat
 else
   ! Symrel^T k
   do isym=1,nsym
     my_symmat(:,:,isym) = transpose(symmat(:,:,isym))
   end do
 end if

 ABI_MALLOC_IFNOT(self%rank2symtime_, (self%min_rank:self%max_rank))
 ABI_MALLOC_IFNOT(self%rank2ikpt_, (self%min_rank:self%max_rank))
 self%rank2ikpt_ = -1 !; self%rank2symtime_ = -1

 do ikpt1=1,self%npoints

   do itimrev=0,timrev
     do isym=1,nsym
       ! Do not use magnetic symmetries.
       if (symafm(isym) == -1) cycle

       kpt1a = (1 - 2*itimrev) * matmul(my_symmat(:, :, isym), self%kpts(:, ikpt1))
       irank = self%get_rank(kpt1a)
       if (self%rank2ikpt_(irank) == -1) then
         self%rank2ikpt_(irank) = ikpt1
         self%rank2symtime_(irank) = isym + itimrev * nsym
       end if
     end do
   end do

 end do

 dksqmax = zero
 do ikpt2=1,nkpt2
   irank = self%get_rank(kptns2(:, ikpt2) + my_qpt)
   ikpt1 = self%rank2ikpt_(irank)
   ii = self%rank2symtime_(irank)
   isym = 1 + mod(ii - 1, nsym)
   itimrev = (ii - 1) / nsym
   kpt1a = (1 - 2 * itimrev) * matmul(my_symmat(:, :, isym), self%kpts(:, ikpt1))
   dk(:) = kptns2(:,ikpt2) + my_qpt - kpt1a(:)
   dkint(:) = nint(dk(:) + tol12)

   indkk(1, ikpt2) = ikpt1
   indkk(2, ikpt2) = isym
   indkk(3:5, ikpt2) = dkint(:)
   indkk(6, ikpt2) = itimrev

   ! Compute norm of the difference vector.
   dk(:) = dk(:) - dkint(:)

   dksqmax = max(dksqmax, &
                 gmet(1,1)*dk(1)**2 + gmet(2,2)*dk(2)**2 + gmet(3,3)*dk(3)**2 + &
                 two * (gmet(2,1)*dk(2)*dk(1) + gmet(3,2)*dk(3)*dk(2)+gmet(3,1)*dk(3)*dk(1)))

   !if (dksqmax > tol8) then
   !  print *, "kbase:", self%kpts(:, ikpt1)
   !  print *, "k", kptns2(:, ikpt2)
   !  print *, "k + q", kptns2(:, ikpt2) + my_qpt
   !  print *, "dk", dk, "dksqmax", dksqmax
   !end if
 end do

 !ABI_FREE(self%rank2ikpt_)
 !ABI_FREE(self%rank2symtime_)

end subroutine krank_get_mapping
!!***

!----------------------------------------------------------------------

!!****f* m_krank/get_ibz2bz
!! NAME
!! get_ibz2bz
!!
!! FUNCTION
!!  Return array with the index of the IBZ wave vectors in the BZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine get_ibz2bz(nibz, nbz, bz2ibz, ibz2bz, err_msg, ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nibz, nbz
 integer,intent(in) :: bz2ibz(6, nbz)
 integer,intent(out) :: ierr
 character(len=*),intent(out) :: err_msg
!arrays
 integer,allocatable,intent(out) :: ibz2bz(:)

!Local variables-------------------------------
!scalars
 integer :: iq_bz, iq_ibz, isym_q, trev_q, cnt, g0_q(3)
 logical :: isirr_q

!----------------------------------------------------------------------

 ABI_MALLOC(ibz2bz, (nibz))

 cnt = 0
 do iq_bz=1,nbz
   iq_ibz = bz2ibz(1, iq_bz); isym_q = bz2ibz(2, iq_bz)
   trev_q = bz2ibz(6, iq_bz); g0_q = bz2ibz(3:5,iq_bz)
   isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
   if (isirr_q) then
     cnt = cnt + 1
     ibz2bz(iq_ibz) = iq_bz
   end if
 end do

 ierr = merge(0, 1, cnt == nibz)
 err_msg = ""
 if (ierr /= 0) then
   err_msg = sjoin("The number of points in the IBZ computed from symmetry table is: ", itoa(cnt), " while it should be: ", itoa(nibz))
 end if

end subroutine get_ibz2bz
!!***

!----------------------------------------------------------------------

!!****f* m_krank/star_from_ibz_idx
!! NAME
!! star_from_ibz_idx
!!
!! FUNCTION
!!  Return array with the indices of the star of the ik_ibz wavevector in the IBZ.
!!  Return number of points in the star, allocate array with the indices of the
!!  star points in the BZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine star_from_ibz_idx(ik_ibz, nkbz, bz2ibz, nk_in_star, kstar_bz_inds)

!Arguments ------------------------------------
 integer,intent(in) :: ik_ibz, nkbz
 integer,intent(out) :: nk_in_star
 integer,intent(in) :: bz2ibz(6, nkbz)
 integer,allocatable,intent(out) :: kstar_bz_inds(:)

!Local variables-------------------------------
!scalars
 integer :: iq_bz

!----------------------------------------------------------------------

 nk_in_star = count(bz2ibz(1, :) == ik_ibz)
 ABI_MALLOC(kstar_bz_inds, (nk_in_star))

 nk_in_star = 0
 do iq_bz=1,nkbz
   if (bz2ibz(1, iq_bz) /= ik_ibz) cycle
   nk_in_star = nk_in_star + 1
   kstar_bz_inds(nk_in_star) = iq_bz
 end do

end subroutine star_from_ibz_idx
!!***

!----------------------------------------------------------------------

end module m_krank
!!***
