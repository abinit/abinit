!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_polynomial_term
!!
!! NAME
!! m_polynomial_term
!!
!! FUNCTION
!! COPYRIGHT
!! Copyright (C) 2010-2015 ABINIT group (AM)
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

module m_polynomial_term

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 public :: polynomial_term_init
 public :: polynomial_term_free
 public :: polynomial_term_dot
!!***

!!****t* m_polynomial_term/polynomial_term_type
!! NAME
!! polynomial_term_type
!!
!! FUNCTION
!! structure for specific displacements of term
!!
!! SOURCE

 type, public :: polynomial_term_type

   integer,allocatable :: atindx(:,:)
!     atindx(2,ndisp)
!     Indexes of the atoms a and b in the unit cell

   integer,allocatable :: cell(:,:,:)
!     cell(3,2,ndisp)
!     indexes of the cell of the atom a and b

   integer,allocatable :: direction(:)
!     direction(ndisp)
!     direction of the displacement

   integer :: ndisp = zero
!     Number of displacement for this terms
!     1 for (X_y-O_y)^3, 2 for (X_y-O_y)(X_x-O_y)^2...

   integer,allocatable :: power(:)
!     power(ndisp)
!     power of the displacement 2 (X_z-O_z)^2 or 1 for (X_y-O_y)^1

   real(dp) :: weight = zero
!     weight of the term


 end type polynomial_term_type
!!***

interface operator (==)
  module procedure terms_compare
end interface

CONTAINS  !===========================================================================================


!!****f* m_polynomial_term/polynomial_term_init
!!
!! NAME
!! polynomial_term_init
!!
!! FUNCTION
!! Initialize polynomial_term_init for given set of displacements
!!
!! INPUTS
!! atindx(2) = Indexes of the atoms a and b in the unit cell
!! cell(3,2) = Indexes of the cell of the atom a and b
!! ndisp     = Number of displacement for this terms
!! power     = Power of the displacement 2 (X_z-O_z)^2 or 1 for (X_y-O_y)^1
!! weight    = Weight of the term
!!
!! OUTPUT
!! polynomial_term = polynomial_term structure to be initialized
!!
!! PARENTS
!!      m_effective_potential_file,m_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_term_init(atindx,cell,direction,ndisp,polynomial_term,power,weight)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_term_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ndisp
 real(dp),intent(in) :: weight
!arrays
 integer, intent(in) :: atindx(2,ndisp)
 integer, intent(in) :: cell(3,2,ndisp)
 integer, intent(in) :: direction(ndisp),power(ndisp)
 type(polynomial_term_type), intent(out) :: polynomial_term
!Local variables-------------------------------
!scalar
 integer :: idisp1,idisp2,ndisp_tmp
!arrays
 integer :: power_tmp(ndisp)
 character(500) :: msg

! *************************************************************************

!Do some checks
 if (size(atindx,2) /= ndisp) then
   write(msg,'(a)')' atindx and ndisp have not the same size'
   MSG_ERROR(msg)
 end if

 if (size(cell,3) /= ndisp) then
   write(msg,'(a)')' cell and ndisp have not the same size'
   MSG_ERROR(msg)
 end if

 if (size(direction) /= ndisp) then
   write(msg,'(a)')' direction and ndisp have not the same size'
   MSG_ERROR(msg)
 end if

 if (size(power) /= ndisp) then
   write(msg,'(a)')' power and ndisp have not the same size'
   MSG_ERROR(msg)
 end if

!First free structure before init
 call polynomial_term_free(polynomial_term)

!Check if displacement are identical, in this case
!increase the power 
 power_tmp(:) = power(:)

 do idisp1=1,ndisp
   do idisp2=idisp1,ndisp
     if (idisp1/=idisp2.and.&
&        atindx(1,idisp1)   == atindx(1,idisp2).and.&
&        atindx(2,idisp1)   == atindx(2,idisp2).and.&
&        direction(idisp1)  == direction(idisp2).and.&
&         all(cell(:,1,idisp1)==cell(:,1,idisp2)).and.&
&         all(cell(:,2,idisp1)==cell(:,2,idisp2)).and.&
&        power_tmp(idisp2) > 0 )then
       power_tmp(idisp1) = power_tmp(idisp1) + 1
       power_tmp(idisp2) = 0
     end if
   end do
 end do

!Count the number of power avec the previous check
!or just remove the power equal to zero
 ndisp_tmp=zero
 do idisp1=1,ndisp
   if(power_tmp(idisp1) > zero) then
     ndisp_tmp = ndisp_tmp + 1
   end if
 end do

!init the values
 polynomial_term%ndisp  = ndisp_tmp
 polynomial_term%weight = weight

 ABI_ALLOCATE(polynomial_term%atindx,(2,polynomial_term%ndisp))
 ABI_ALLOCATE(polynomial_term%direction,(polynomial_term%ndisp))
 ABI_ALLOCATE(polynomial_term%cell,(3,2,polynomial_term%ndisp))
 ABI_ALLOCATE(polynomial_term%power,(polynomial_term%ndisp))

 idisp2 = zero
 do idisp1=1,ndisp
   if(power_tmp(idisp1) > zero)then
     idisp2 =  idisp2 + 1
     polynomial_term%atindx(:,idisp2) = atindx(:,idisp1) 
     polynomial_term%direction(idisp2) = direction(idisp1)
     polynomial_term%cell(:,:,idisp2) = cell(:,:,idisp1)
     polynomial_term%power(idisp2) = power_tmp(idisp2)
   end if
 end do

end subroutine polynomial_term_init
!!***


!!****f* m_polynomial_term/polynomial_term_free
!!
!! NAME
!! polynomial_term_free
!!
!! FUNCTION
!! Free polynomial_term
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_term = polynomial_term structure to be free
!!
!! PARENTS
!!      m_effective_potential_file,m_polynomial_coeff,m_polynomial_term
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_term_free(polynomial_term)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_term_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(polynomial_term_type), intent(inout) :: polynomial_term
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************

 polynomial_term%ndisp     = zero
 polynomial_term%weight    = zero

 if(allocated(polynomial_term%atindx))then
   polynomial_term%atindx(:,:) = zero
   ABI_DEALLOCATE(polynomial_term%atindx)
 end if

 if(allocated(polynomial_term%cell))then
   polynomial_term%cell(:,:,:) = zero
   ABI_DEALLOCATE(polynomial_term%cell)
 end if

 if(allocated(polynomial_term%direction))then
   polynomial_term%direction(:) = zero
   ABI_DEALLOCATE(polynomial_term%direction)
 end if

 if(allocated(polynomial_term%power))then
   polynomial_term%power(:) = zero
   ABI_DEALLOCATE(polynomial_term%power)
 end if


end subroutine polynomial_term_free
!!***

!!****f* m_polynomial_term/polynomial_term_dot
!!
!! NAME
!! polynomial_term_dot
!!
!! FUNCTION
!! Return the multiplication of two terms
!!
!! INPUTS
!! term1_in = Firts term
!! term2_in = Second term
!!
!! OUTPUT
!! term_out = multiplication of the two input terms
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_term_dot(term_out,term1_in,term2_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_term_dot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(polynomial_term_type), intent(in) :: term1_in
 type(polynomial_term_type), intent(in) :: term2_in
 type(polynomial_term_type), intent(out):: term_out
!Local variables-------------------------------
!scalar
 integer :: idisp1,idisp2,ndisp,new_ndisp
!arrays

! *************************************************************************

!Get the number of displacement for this new term
 new_ndisp = term1_in%ndisp + term2_in%ndisp

 ndisp = 0

! do while (ndisp < new_ndisp)
   do idisp1=1,term1_in%ndisp
     do idisp2=1,term2_in%ndisp
       if (term1_in%atindx(1,idisp1)/=&
&          term2_in%atindx(1,idisp2).or.&
&          term1_in%atindx(2,idisp1)/=&
&          term2_in%atindx(2,idisp2).or.&
&          term1_in%direction(idisp1)/=&
&          term2_in%direction(idisp2))  then
         ndisp = ndisp + 1
       end if
     end do
   end do
! end do

end subroutine polynomial_term_dot
!!***

!!****f* m_polynomial_term/terms_compare
!! NAME
!!  equal
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function terms_compare(t1,t2) result (res)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'terms_compare'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  type(polynomial_term_type), intent(in) :: t1,t2
  logical :: res
!local
!variable
  integer :: ndisp
  integer :: ia,idisp1,idisp2,mu
  logical :: found
!array
  integer :: blkval(2,t1%ndisp)
! *************************************************************************
  res = .true.
  blkval = zero
  ndisp = 0
  if(t1%ndisp==t2%ndisp)then
    do idisp1=1,t1%ndisp
      if(blkval(1,idisp1)==one)cycle!already found
      do idisp2=1,t2%ndisp
        if(blkval(2,idisp2)==one)cycle!already found
        found = .false.
        if(t1%atindx(1,idisp1)  ==  t2%atindx(1,idisp2).and.&
&          t1%atindx(2,idisp1)  ==  t2%atindx(2,idisp2).and.&
&          t1%direction(idisp1) ==  t2%direction(idisp2).and.&
&          t1%power(idisp1) == t2%power(idisp2))then
          found=.true.
          do ia=1,2
            do mu=1,3
              if(t1%cell(mu,ia,idisp1) /= t2%cell(mu,ia,idisp2))then
                found = .false.
                cycle
              end if
            end do
          end do
          if(found)then
            blkval(1,idisp1)=one
            blkval(2,idisp2)=one
          end if
        end if
      end do
    end do
   if(any(blkval(:,:)==zero))res = .false.
  else
    res = .false.
  end if

!     do idisp1=1,t1%ndisp
!       do idisp2=1,t2%ndisp
!         if(t1%atindx(1,idisp1) ==  t2%atindx(1,idisp2).and.&
! &          t1%atindx(2,idisp1) ==  t2%atindx(2,idisp2).and.&
! &          t1%direction(idisp1) == t2%direction(idisp2).and.&
! &          t1%power(idisp1) == t2%power(idisp2))then!.and.&
! !&          t1%weight == t2%weight)then
!           found = .true.
!           do ia=1,2
!             do mu=1,3
!               if(t1%cell(mu,ia,idisp1) /= t2%cell(mu,ia,idisp2))then
!                 found = .false.
!                 cycle
!               end if
!             end do
!           end do
!           if (found)then
!             ndisp = ndisp +1 
!           end if
!         end if
!       end do
!     end do

!     if(ndisp == t1%ndisp)then
!       res = .true.
!     end if
!  end if

end function terms_compare
!!***


end module m_polynomial_term
!!***
