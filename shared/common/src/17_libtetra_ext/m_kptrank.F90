!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kptrank
!! NAME
!! m_kptrank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (MVer,HM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "libtetra.h"

module m_kptrank

 USE_MEMORY_PROFILING
 USE_MSG_HANDLING

 implicit none

 private

 double precision :: zero = 0.0d0, half = 0.5d0, one = 1.0d0, tol8 = 1.d-8,  tol10 = 1.d-10, tol12 = 1.d-12
!!***

!!****t* m_kptrank/kptrank_type
!! NAME
!! kptrank_type
!!
!! FUNCTION
!!  structure to contain a rank/inverse rank pair of arrays, with dimensions
!!
!! SOURCE

 type,public :: kptrank_type
   integer :: max_linear_density
   integer :: min_rank
   integer :: max_rank
   integer :: npoints
   logical :: time_reversal
   integer,allocatable :: invrank(:)
 end type kptrank_type

 public :: mkkptrank       ! Sets up the kpt ranks for comparing kpts
 public :: get_rank_1kpt   ! Calculates the rank for one kpt
 public :: kptrank_index   ! Return the index of the k-point `kpt` in the initial set. -1 if not found.
 public :: copy_kptrank    ! Copy the object
 public :: destroy_kptrank ! Free memory
 public :: dump_kptrank    ! Prints the arrays and dimensions of a kptrank_type structure
!!***

contains
!!***

!!****f* m_kptrank/mkkptrank
!!
!! NAME
!! mkkptrank
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
!! OUTPUT
!!  krank = object containing ranking and inverse ranking
!!
!! PARENTS
!!      get_full_kgrid,m_ddk,m_ebands,m_fstab,m_nesting,m_phgamma,m_pptools
!!      m_tetrahedron,mkfskgrid,mkqptequiv,order_fs_kpts,outelph,read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkkptrank (kpt,nkpt,krank,nsym,symrec, time_reversal)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 integer,intent(in), optional :: nsym
 logical,intent(in), optional :: time_reversal
!arrays
 type(kptrank_type), intent(out) :: krank
 double precision,intent(in) :: kpt(3,nkpt)
 integer,intent(in), optional :: symrec(3,3, *)

!Local variables -------------------------
!scalars
 integer :: ikpt, isym, symkptrank, irank
 integer :: timrev, itim
 double precision :: smallestlen
 character(len=500) :: msg
!arrays
 double precision :: symkpt(3)

! *********************************************************************

! find smallest linear length
 smallestlen = one
 do ikpt=1, nkpt
   if (abs(kpt(1,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpt(1,ikpt)))
   if (abs(kpt(2,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpt(2,ikpt)))
   if (abs(kpt(3,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpt(3,ikpt)))
 end do

 krank%max_linear_density = nint(one/smallestlen)
 krank%npoints = nkpt
 krank%min_rank = nint(real(krank%max_linear_density)*(half+tol8 +&
&                      real(krank%max_linear_density)*(half+tol8 +&
&                      real(krank%max_linear_density)*(half+tol8))))

 krank%max_rank = nint(real(krank%max_linear_density)*(1+half+tol8 +&
&                      real(krank%max_linear_density)*(1+half+tol8 +&
&                      real(krank%max_linear_density)*(1+half+tol8))))

 TETRA_ALLOCATE(krank%invrank, (krank%min_rank:krank%max_rank))
 krank%invrank(:) = -1

 timrev = 2
 krank%time_reversal = .true.
 if (present(time_reversal)) then
   if (.not. time_reversal) timrev = 1
 end if

!Ensure kpt(i)+one is positive, and the smallest
!difference between kpts should be larger than 1/100
!ie ngkpt < 100.
! the following fills invrank for the k-points in the list provided (may be only the irred kpts)
 do ikpt=1,nkpt
   call get_rank_1kpt (kpt(:,ikpt), irank, krank)

   if (irank > krank%max_rank .or. irank < krank%min_rank) then
     write(msg,'(a,2i0)')" rank above max_rank or bellow min_rank, ikpt, rank ", ikpt, irank
     TETRA_ERROR(msg)
   end if
   krank%invrank(irank) = ikpt
 end do

! if symrec is provided, fill invrank with appropriate irred kpt indices
! for symmetry completion: kptrank_t%invrank points to the irred k-point
! equivalent to the k-point whose rank is provided
 if (present(symrec)) then
   if(.not. present(nsym)) then
     TETRA_ERROR("need both symrec and nsym arguments together")
   end if
   do ikpt=1,nkpt
     ! itim == 1 for positive, and itim==2 gives Kramers opposite of k-point
     ! favor the former by looping it last
     do itim = timrev, 1, -1
       do isym = 1, nsym
         symkpt = (-1)**(timrev+1) * matmul(symrec(:,:,isym), kpt(:, ikpt))
         call get_rank_1kpt (symkpt(:), symkptrank, krank)
         krank%invrank(symkptrank) = ikpt
       end do
     end do
   end do
 end if

end subroutine mkkptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/get_rank_1kpt
!!
!! NAME
!! get_rank_1kpt
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
!! PARENTS
!!      elphon,get_full_kgrid,integrate_gamma,integrate_gamma_alt,k_neighbors
!!      m_ddk,m_kptrank,m_nesting,m_pptools,m_tetrahedron,mkfskgrid,mkqptequiv
!!      read_el_veloc,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_rank_1kpt(kpt,rank,krank)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: rank
 type(kptrank_type), intent(in) :: krank
!arrays
 double precision,intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays
 double precision :: redkpt(3)

! *************************************************************************

! wrap to [0, 1[ -> replaced call to wrap2_zeroone inline, to encapsulate this module
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

! rank = int(real(krank%max_linear_density)*(redkpt(3)+half+tol8 +&
!&           real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
!&           real(krank%max_linear_density)*(redkpt(1)+half+tol8))))
 rank = nint(real(krank%max_linear_density)*(redkpt(1)+half+tol8 +&
&           real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
&           real(krank%max_linear_density)*(redkpt(3)+half+tol8))))

 if (rank > krank%max_rank) then
   write(msg,'(a,i0,a,i0)') ' rank should be inferior to ', krank%max_rank, ' got ', rank
   TETRA_ERROR(msg)
 end if

end subroutine get_rank_1kpt
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/kptrank_index
!!
!! NAME
!! kptrank_index
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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function kptrank_index(krank, kpt) result(ikpt)

!Arguments ------------------------------------
!scalars
 type(kptrank_type), intent(in) :: krank
!arrays
 double precision,intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 integer :: kpt_rank

! *************************************************************************

 call get_rank_1kpt(kpt, kpt_rank, krank)
 ikpt = -1
 if (kpt_rank < krank%max_rank) ikpt = krank%invrank(kpt_rank)

end function kptrank_index
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/copy_kptrank
!!
!! NAME
!! copy_kptrank
!!
!! FUNCTION
!! Copy the object
!!
!! INPUTS
!!
!! OUTPUT
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! PARENTS
!!      defs_elphon,elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_kptrank (krank_in, krank_out)

!Arguments ------------------------------------
!scalars
 type(kptrank_type), intent(in) :: krank_in
 type(kptrank_type), intent(out) :: krank_out

! *********************************************************************
 krank_out%max_linear_density = krank_in%max_linear_density
 krank_out%min_rank = krank_in%min_rank
 krank_out%max_rank = krank_in%max_rank
 krank_out%npoints = krank_in%npoints

 TETRA_ALLOCATE(krank_out%invrank, (krank_out%min_rank:krank_out%max_rank))
 krank_out%invrank = krank_in%invrank

end subroutine copy_kptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/destroy_kptrank
!!
!! NAME
!! destroy_kptrank
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! PARENTS
!!      defs_elphon,get_full_kgrid,m_ddk,m_ebands,m_fstab,m_nesting,m_phgamma
!!      m_pptools,m_tetrahedron,mkfskgrid,mkqptequiv,order_fs_kpts,outelph
!!      read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_kptrank (krank)

!Arguments ------------------------------------
!scalars
 type(kptrank_type), intent(inout) :: krank

! *********************************************************************

 if (allocated(krank%invrank))  then
   TETRA_DEALLOCATE(krank%invrank)
 end if

end subroutine destroy_kptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/dump_kptrank
!!
!! NAME
!! dump_kptrank
!!
!! FUNCTION
!! This routine prints the arrays and dimensions of a kptrank_type structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking
!!  unout = unit for open file to print to
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dump_kptrank (krank, unout)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: unout
!arrays
 type(kptrank_type), intent(in) :: krank

! *********************************************************************

  write(unout, *)
  write(unout, '(a)') ' Dump of the contents of a kptrank_type structure with k-point rank information'
  write(unout, '(a,I8)') ' max linear density of points in 3 directions: max_linear_density = ',  krank%max_linear_density
  write(unout, '(a,I8)') ' maximum rank for any point in grid: max_rank = ',  krank%max_rank
  write(unout, '(a,I8)') ' number of points in input grid: npoints = ',  krank%npoints
  write(unout, *)
  write(unout, '(a)') ' invrank array = '
  write(unout, '(I4)') krank%invrank(:)
  write(unout, *)

end subroutine dump_kptrank
!!***

!----------------------------------------------------------------------

end module m_kptrank
!!***
