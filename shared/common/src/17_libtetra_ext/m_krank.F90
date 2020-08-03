!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_krank
!! NAME
!! m_krank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!!
!! COPYRIGHT
!! Copyright (C) 2010-2020 ABINIT group (MVer,HM)
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

module m_krank

 USE_MEMORY_PROFILING
 USE_MSG_HANDLING

 implicit none

 private

 integer, parameter :: dp=kind(1.0d0)

 real(dp) :: zero = 0.0d0, half = 0.5d0, one = 1.0d0, tol8 = 1.d-8,  tol10 = 1.d-10, tol12 = 1.d-12
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
   integer,allocatable :: invrank(:)
   !real(dp),pointer kpoints(:,:)

 contains

   procedure :: get_rank
    ! Calculates the rank for one kpt

   procedure :: get_index => krank_get_index
    ! Return the index of the k-point `kpt` in the initial set. -1 if not found.

   procedure :: copy => krank_copy
    ! Copy the object

   procedure :: free => krank_free
    ! Free memory

   procedure :: dump => krank_dump
    ! Prints the arrays and dimensions of a krank_t structure

 end type krank_t

 public :: krank_new       ! Sets up the kpt ranks for comparing kpts
!!***

contains
!!***

!!****f* m_krank/krank_new
!!
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

type(krank_t) function krank_new(nkpt, kpt, nsym, symrec, time_reversal) result(krank)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 integer,intent(in), optional :: nsym
 logical,intent(in), optional :: time_reversal
!arrays
 real(dp),intent(in) :: kpt(3,nkpt)
 integer,intent(in), optional :: symrec(3,3, *)

!Local variables -------------------------
!scalars
 integer :: ikpt, isym, symkptrank, irank, timrev, itim
 real(dp) :: smallestlen
 character(len=500) :: msg
!arrays
 real(dp) :: symkpt(3)

! *********************************************************************

 !krank%kpts => kpt

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
                       real(krank%max_linear_density)*(half+tol8 +&
                       real(krank%max_linear_density)*(half+tol8))))

 krank%max_rank = nint(real(krank%max_linear_density)*(1+half+tol8 +&
                       real(krank%max_linear_density)*(1+half+tol8 +&
                       real(krank%max_linear_density)*(1+half+tol8))))

 TETRA_ALLOCATE(krank%invrank, (krank%min_rank:krank%max_rank))
 krank%invrank(:) = -1

 timrev = 2
 krank%time_reversal = .true.
 if (present(time_reversal)) then
   if (.not. time_reversal) timrev = 1
   krank%time_reversal = .false.
 end if

 ! Ensure kpt(i)+one is positive, and the smallest difference between kpts should be larger than 1/100 ie ngkpt < 100.
 ! the following fills invrank for the k-points in the list provided (may be only the irred kpts)
 do ikpt=1,nkpt
   irank = krank%get_rank(kpt(:,ikpt))

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
         symkpt = (-1)**(itim+1) * matmul(symrec(:,:,isym), kpt(:, ikpt))
         symkptrank = krank%get_rank(symkpt(:))
         krank%invrank(symkptrank) = ikpt
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
!! PARENTS
!!      elphon,get_full_kgrid,integrate_gamma,integrate_gamma_alt,k_neighbors
!!      m_ddk,m_krank,m_nesting,m_pptools,m_tetrahedron,mkfskgrid,mkqptequiv
!!      read_el_veloc,read_gkk
!!
!! CHILDREN
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
!            real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
!            real(krank%max_linear_density)*(redkpt(1)+half+tol8))))
 rank = nint(real(krank%max_linear_density)*(redkpt(1)+half+tol8 +&
             real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
             real(krank%max_linear_density)*(redkpt(3)+half+tol8))))

 if (rank > krank%max_rank) then
   write(msg,'(a,i0,a,i0)') ' rank should be inferior to ', krank%max_rank, ' got ', rank
   TETRA_ERROR(msg)
 end if

end function get_rank
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_get_index
!!
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
!! PARENTS
!!
!! CHILDREN
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

type(krank_t) function krank_copy(krank_in) result(krank_out)

!Arguments ------------------------------------
!scalars
 class(krank_t), intent(in) :: krank_in

! *********************************************************************

 krank_out%max_linear_density = krank_in%max_linear_density
 krank_out%min_rank = krank_in%min_rank
 krank_out%max_rank = krank_in%max_rank
 krank_out%npoints = krank_in%npoints

 TETRA_ALLOCATE(krank_out%invrank, (krank_out%min_rank:krank_out%max_rank))
 krank_out%invrank = krank_in%invrank
 !krank_out%kpts => krank_in

end function krank_copy
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_free
!!
!! NAME
!! krank_free
!!
!! FUNCTION
!! This routine deallocates the arrays in a krank_t structure
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

subroutine krank_free(krank)

!Arguments ------------------------------------
 class(krank_t), intent(inout) :: krank

! *********************************************************************

 if (allocated(krank%invrank))  then
   TETRA_DEALLOCATE(krank%invrank)
 end if

end subroutine krank_free
!!***

!----------------------------------------------------------------------

!!****f* m_krank/krank_dump
!!
!! NAME
!! krank_dump
!!
!! FUNCTION
!! This routine prints the arrays and dimensions of a krank_t structure
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

subroutine krank_dump (krank, unout)

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

end subroutine krank_dump
!!***

!----------------------------------------------------------------------

subroutine krank_map(self, nkpt1, kptns1, dksqmax, gmet, indkk, nsym, symafm, symmat, timrev, comm, &
                     use_symrec) ! optional

!Arguments ------------------------------------
!scalars
 class(krank_t),intent(in) :: self
 integer,intent(in) :: nkpt1, nsym, timrev, comm !nkpt2,
 real(dp),intent(out) :: dksqmax
 logical,optional,intent(in) :: use_symrec
!arrays
 integer,intent(in) :: symafm(nsym),symmat(3,3,nsym)
 integer,intent(out) :: indkk(6, self%npoints)
 real(dp),intent(in) :: gmet(3,3),kptns1(3,nkpt1) !,kptns2(3,nkpt2)

!Local variables-------------------------------
!scalars
 integer :: ikpt1, itimrev, isym, ik2_rank, isk, ierr, my_rank, nprocs
 real(dp) :: kpt1a(3)

! *************************************************************************

 !my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 indkk = 0

 do ikpt1=1,nkpt1
   !if (mod(ikpt1, nprocs) /= my_rank) cycle ! MPI parallelism.

   itimrev_loop: &
   do itimrev=0,timrev
     do isym=1,nsym
       ! Select magnetic characteristic of symmetries
       if (symafm(isym) == -1) cycle
       !if (isppol == 1 .and. symafm(isym) == -1) cycle
       !if (isppol == 2 .and. symafm(isym) == 1) cycle

       ! original code only used transpose(symrel)
       if (present(use_symrec)) then
         if (use_symrec) then
           kpt1a(:) = MATMUL(symmat(:,:,isym), kptns1(:,ikpt1))
         else
           kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)), kptns1(:,ikpt1))
         end if
       else
         kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)), kptns1(:,ikpt1))
       end if
       kpt1a(:) = (1-2*itimrev)*kpt1a(:)

       !ik2_rank = self%get_rank(kpt1a)
       !if ik2_rank
       isk = self%get_index(kpt1a)
       if (isk /= -1) then
         !indkk(1, isk) = ikpt1
         !indkk(2, isk) = isym
         !indkk(3:5, isk) = jdkint(:)
         !indkk(6, isk) = itimrev

         ! Compute norm of the difference vector, and update kpt1 if better.
         !dksq=gmet(1,1)*dk(1)**2+gmet(2,2)*dk(2)**2+ &
         !     gmet(3,3)*dk(3)**2+two*(gmet(2,1)*dk(2)*dk(1)+ &
         !     gmet(3,2)*dk(3)*dk(2)+gmet(3,1)*dk(3)*dk(1))

         !dksqmax = max(dksqmax, dksqmn)
         exit itimrev_loop
       end if

     end do
   end do itimrev_loop
 end do

 !if (nprocs > 1) call xmpi_sum(indkk, comm, ierr)

end subroutine krank_map
!!***

end module m_krank
!!***
