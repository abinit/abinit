!!****f* ABINIT/m_polynomial_term
!!
!! NAME
!! m_polynomial_term
!!
!! FUNCTION
!! Module with the datatype polynomial terms
!!
!! COPYRIGHT
!! Copyright (C) 2010-2025 ABINIT group (AM)
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
 use m_abicore

 implicit none

 public :: polynomial_term_init
 public :: polynomial_term_free
 public :: polynomial_term_copy
 public :: terms_compare
 public :: terms_compare_inverse
!!***

!!****t* m_polynomial_term/polynomial_term_type
!! NAME
!! polynomial_term_type
!!
!! FUNCTION
!! Datatype for a terms  (displacements or strain)
!! related to a polynomial coefficient
!!
!! SOURCE

 type, public :: polynomial_term_type

   integer :: ndisp = 0
!     Number of displacement for this terms
!     1 for (X_y-O_y)^3, 2 for (X_y-O_y)(X_x-O_y)^2...

   integer :: nstrain = 0
!     Number of strain for this terms

   integer :: nindex = -1
!     Number of index

   integer,allocatable :: atindx(:,:)
!     atindx(2,ndisp)
!     Indexes of the atoms a and b in the unit cell

   integer,allocatable :: cell(:,:,:)
!     cell(3,2,ndisp)
!     indexes of the cell of the atom a and b

   integer,allocatable :: direction(:)
!     direction(ndisp)
!     direction of the displacement

   integer,allocatable :: strain(:)
!     strain(nstrain)
!     strain

   integer,allocatable :: power_disp(:)
!     power_disp(ndisp)
!     power of the displacement 2 (X_z-O_z)^2 or 1 for (X_y-O_y)^1

   integer,allocatable :: power_strain(:)
!     power_strain(nstrain)
!     power of the strain 2 (\eta_1)^2 or 1 for (\eta_1)^1

   real(dp) :: weight = zero
!     weight of the term

   integer,allocatable  ::  index_coeff(:)

   ! a string to identify the term
    character(len=500) :: debug_str = ''

   contains
           procedure :: get_nbody
           procedure :: get_total_power
           !final :: polynomial_term_finalizer
 end type polynomial_term_type
!!***


 interface  operator (==)
  module procedure terms_compare
end interface

CONTAINS  !===========================================================================================


function get_nbody(self) result(nbody)
  class(polynomial_term_type), intent(in) :: self
  integer :: nbody
  nbody = self%ndisp + self%nstrain
end function get_nbody


function get_total_power(self) result(total_power)
  class(polynomial_term_type), intent(in) :: self
  integer :: total_power
  total_power = sum(self%power_disp) + sum(self%power_strain)
end function get_total_power


!!****f* m_polynomial_term/polynomial_term_init
!!
!! NAME
!! polynomial_term_init
!!
!! FUNCTION
!! Initialize a polynomial_term_init for given set of displacements/strain
!!
!! INPUTS
!! atindx(2) = Indexes of the atoms a and b in the unit cell
!! cell(3,2) = Indexes of the cell of the atom a and b
!! direction = direction of the perturbation => 1,2,3 for atomic displacement
!! ndisp     = Number of displacement for this terms
!! nstrain   = Number of strain for this terms
!! power_disp   = power_disp of the displacement 2 (X_z-O_z)^2 or 1 for (X_y-O_y)^1
!! power_strain = power_strain of the strain 2 (\eta)^2 or 1 for (\eta)^1
!! strain(nstrain) = index of strain, 1 2 3 4 5 or 6
!! weight     = Weight of the term
!! check      = optional,logical => if TRUE, the term will be check, same displacement/strain
!!                                          are gathered in the same displacement but with an higher
!!                                          power_disp. For example:
!!                                          ((Sr_y-O1_y)^1(Sr_y-O1_y)^1 => (Sr_y-O1_y)^2)
!!                                 if FALSE, default, do nothing
!!
!! OUTPUT
!! polynomial_term<type(polynomial_term)> = polynomial_term datatype is now initialized
!!
!! SOURCE

subroutine polynomial_term_init(atindx,cell,direction,ndisp,nstrain,polynomial_term,power_disp,&
&                               power_strain,strain,weight,check, index_coeff, debug_str)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ndisp,nstrain
 real(dp),intent(in) :: weight
 logical,optional,intent(in)  :: check
!arrays
 integer, intent(in) :: atindx(2,ndisp)
 integer, intent(in) :: cell(3,2,ndisp)
 integer, intent(in) :: direction(ndisp),power_disp(ndisp)
 integer, intent(in) :: strain(nstrain),power_strain(nstrain)
 integer, optional,  intent(in) :: index_coeff(:)
 type(polynomial_term_type), intent(out) :: polynomial_term
!Local variables-------------------------------
!scalar
 integer :: idisp1,idisp2,ndisp_tmp,nstrain_tmp
 logical :: check_in
!arrays
 integer :: power_disp_tmp(ndisp),power_strain_tmp(nstrain)
 character(500) :: msg

character(*), optional :: debug_str

! *************************************************************************

if (present(debug_str)) then
  polynomial_term%debug_str = debug_str
else
  polynomial_term%debug_str = ' '
end if

!Do some checks
 if (size(atindx,2) /= ndisp) then
   write(msg,'(a)')' atindx and ndisp have not the same size'
   ABI_ERROR(msg)
 end if

 if (size(cell,3) /= ndisp) then
   write(msg,'(a)')' cell and ndisp have not the same size'
   ABI_ERROR(msg)
 end if

 if (size(direction) /= ndisp) then
   write(msg,'(a)')' direction and ndisp have not the same size'
   ABI_ERROR(msg)
 end if

 if (size(power_disp) /= ndisp) then
   write(msg,'(a)')' power_disp and ndisp have not the same size'
   ABI_ERROR(msg)
 end if

 if (size(power_strain) /= nstrain) then
   write(msg,'(a)')' power_strain and nstrain have not the same size'
   ABI_ERROR(msg)
 end if

 if (size(strain) /= nstrain) then
   write(msg,'(a)')' strain and nstrain have not the same size'
   ABI_ERROR(msg)
 end if
!First free datatype before init
 call polynomial_term_free(polynomial_term)
 check_in = .false.

!Copy the power array before check
 power_disp_tmp(:) = power_disp(:)
 power_strain_tmp(:) = power_strain(:)

 if(present(check)) check_in = check

 if(check_in)then
!Check if displacement are identical, in this case
!increase the power_disp
   do idisp1=1,ndisp
     do idisp2=idisp1,ndisp
       if (idisp1/=idisp2.and.&
&        atindx(1,idisp1)   == atindx(1,idisp2).and.&
&        atindx(2,idisp1)   == atindx(2,idisp2).and.&
&        direction(idisp1)  == direction(idisp2).and.&
&         all(cell(:,1,idisp1)==cell(:,1,idisp2)).and.&
&         all(cell(:,2,idisp1)==cell(:,2,idisp2)).and.&
&        power_disp_tmp(idisp2) > 0 )then
         power_disp_tmp(idisp1) = power_disp_tmp(idisp1) + 1
         power_disp_tmp(idisp2) = 0
       end if
     end do
   end do

! Count the number of power_disp avec the previous check
! or just remove the power_disp equal to zero
   ndisp_tmp = 0
   do idisp1=1,ndisp
     if(power_disp_tmp(idisp1) > zero) then
       ndisp_tmp = ndisp_tmp + 1
     end if
   end do

!Check if strain are identical, in this case
!increase the power_strain
   do idisp1=1,nstrain
     do idisp2=idisp1,nstrain
       if (idisp1/=idisp2.and.&
&        strain(idisp1)  == strain(idisp2).and.&
&        power_strain_tmp(idisp2) > 0 )then
         power_strain_tmp(idisp1) = power_strain_tmp(idisp1) + 1
         power_strain_tmp(idisp2) = 0
       end if
     end do
   end do

! Count the number of power_strain avec the previous check
! or just remove the power_strain equal to zero
   nstrain_tmp = 0
   do idisp1=1,nstrain
     if(power_strain_tmp(idisp1) > zero) then
       nstrain_tmp = nstrain_tmp + 1
     end if
   end do

 else
   ndisp_tmp   = ndisp
   nstrain_tmp = nstrain
 end if!end check

!init the values
 polynomial_term%ndisp    = ndisp_tmp
 polynomial_term%nstrain  = nstrain_tmp
 polynomial_term%weight   = weight

 ABI_MALLOC(polynomial_term%atindx,(2,polynomial_term%ndisp))
 ABI_MALLOC(polynomial_term%direction,(polynomial_term%ndisp))
 ABI_MALLOC(polynomial_term%cell,(3,2,polynomial_term%ndisp))
 ABI_MALLOC(polynomial_term%power_disp,(polynomial_term%ndisp))
 ABI_MALLOC(polynomial_term%power_strain,(polynomial_term%nstrain))
 ABI_MALLOC(polynomial_term%strain,(polynomial_term%nstrain))

!Transfert displacement
 idisp2 = 0
 do idisp1=1,ndisp
   if(power_disp_tmp(idisp1) > zero)then
     idisp2 =  idisp2 + 1
     polynomial_term%direction(idisp2)  =  direction(idisp1)
     polynomial_term%power_disp(idisp2) =  power_disp_tmp(idisp1)
     polynomial_term%atindx(:,idisp2)   =  atindx(:,idisp1)
     polynomial_term%cell(:,:,idisp2)   =  cell(:,:,idisp1)
     polynomial_term%power_disp(idisp2) =  power_disp_tmp(idisp1)
   end if
 end do

!Transfert strain
 idisp2 = 0
 do idisp1=1,nstrain
   if(power_strain_tmp(idisp1) > zero)then
     idisp2 =  idisp2 + 1
     polynomial_term%power_strain(idisp2) = power_strain_tmp(idisp1)
     polynomial_term%strain(idisp2) = strain(idisp1)
   end if
 end do

 if (present(index_coeff)) then
   polynomial_term%nindex = size(index_coeff)
   ABI_MALLOC(polynomial_term%index_coeff, (polynomial_term%nindex))
   polynomial_term%index_coeff(:)=index_coeff(:)
 else
   polynomial_term%nindex = -1
   ABI_MALLOC(polynomial_term%index_coeff, (0))
 end if




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
!! polynomial_term<type(polynomial_term)> =  datatype to free
!!
!! OUTPUT
!! polynomial_term<type(polynomial_term)> =  datatype to free
!!
!! SOURCE

subroutine polynomial_term_free(polynomial_term)

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(polynomial_term_type), intent(inout) :: polynomial_term
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************

 polynomial_term%ndisp     = 0
 polynomial_term%nstrain   = 0
 polynomial_term%weight    = zero

 ABI_SFREE(polynomial_term%atindx)
 ABI_SFREE(polynomial_term%cell)
 ABI_SFREE(polynomial_term%direction)
 ABI_SFREE(polynomial_term%power_disp)
 ABI_SFREE(polynomial_term%power_strain)
 ABI_SFREE(polynomial_term%strain)
 ABI_SFREE(polynomial_term%index_coeff)
 
end subroutine polynomial_term_free
!!***

!!****f* m_polynomial_term/polynomial_term_finalizer
!!
!! NAME
!! polynomial_term_finalizer
!!
!! FUNCTION
!! Finalizer procedure for polynomial_term_type to automatically free allocated memory
!!
!! SOURCE

subroutine polynomial_term_finalizer(self)
  type(polynomial_term_type), intent(inout) :: self
  !print *, "Warning: polynomial_term_finalizer called. Debug str: ", self%debug_str
  call polynomial_term_free(self)
end subroutine polynomial_term_finalizer
!!***

! function polynomial_term_type_get_index_coeff(term, ndisp, nstrain) result(list)
!   type(polynomial_term_type),  intent(inout) :: term
!   integer, intent(in) :: size
!   integer :: list()
!   integer :: i, counter, ip
!   counter=1
!   do i=1, term%ndisp
!     do ip=1, term%power_disp
!       list(counter) = term%
!     end do
!   end do
! end function polynomial_term_type_get_index_coeff

subroutine polynomial_term_list_free(terms)
  type(polynomial_term_type), allocatable, intent(inout) :: terms(:)
  integer :: iterm
  do iterm=1, size(terms)
    call polynomial_term_free(terms(iterm))
  end do
  ABI_SFREE(terms)
end subroutine polynomial_term_list_free

!!****f* m_polynomial_term/terms_compare
!! NAME
!!  equal
!!
!! FUNCTION
!!  Compare two polynomial_term_dot
!!
!! INPUTS
!! t1<type(polynomial_term)> =  datatype of the first term
!! t2<type(polynomial_term)> =  datatype of the second term
!!
!! OUTPUT
!! res = logical
!!
!! SOURCE

pure function terms_compare(t1,t2) result (res)
!Arguments ------------------------------------
 implicit none

!Arguments ------------------------------------
  type(polynomial_term_type), intent(in) :: t1,t2
  logical :: res
!local
!variable
  integer :: ia,idisp1,idisp2,mu
  logical :: found
!array
  integer :: blkval(2,t1%ndisp+t1%nstrain)
! *************************************************************************
  res = .true.
  blkval(:,:) = 0
  if(t1%ndisp==t2%ndisp.and.t1%nstrain==t2%nstrain)then
!   Check strain
    blkval(:,:) = 0
    do idisp1=1,t1%nstrain
      if(blkval(1,t1%ndisp+idisp1)==1)cycle!already found
      do idisp2=1,t2%nstrain
        if(blkval(2,t1%ndisp+idisp2)==1)cycle!already found
        found = .false.
        if(t1%strain(idisp1) ==  t2%strain(idisp2).and.&
&          t1%power_strain(idisp1) == t2%power_strain(idisp2))then
          found=.true.
        end if
        if(found)then
          blkval(1,t1%ndisp+idisp1) = 1
          blkval(2,t1%ndisp+idisp2) = 1
        end if
      end do
    end do
    if(any(blkval(:,t1%ndisp+1:t1%ndisp+t1%nstrain) == 0))then
      res = .false.
      return
    end if
!   Check displacement
    do idisp1=1,t1%ndisp
      if(blkval(1,idisp1)==1)cycle!already found
      do idisp2=1,t2%ndisp
        if(blkval(2,idisp2)==1)cycle!already found
        found = .false.
        if(t1%atindx(1,idisp1)  ==  t2%atindx(1,idisp2).and.&
&          t1%atindx(2,idisp1)  ==  t2%atindx(2,idisp2).and.&
&          t1%direction(idisp1) ==  t2%direction(idisp2).and.&
&          t1%power_disp(idisp1) == t2%power_disp(idisp2))then
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
            blkval(1,idisp1) = 1
            blkval(2,idisp2) = 1
          end if
        end if
      end do
    end do
    if(any(blkval(:,:)==0))res = .false.
  else
    res = .false.
  end if
end function terms_compare
!!***


function terms_compare_inverse(t1, t2) result(res)
  type(polynomial_term_type), intent(in) :: t1,t2
  logical :: res
  type(polynomial_term_type) :: t3
  call polynomial_term_copy(t2, t3)
  t3%atindx(1,:) =t2%atindx(2, :)
  t3%atindx(2,:) =t2%atindx(1, :)
  t3%cell(1,:, :) = 0
  t3%cell(2, :, :) = -t2%cell(2,:, : )
  res=terms_compare(t1, t3)
  call polynomial_term_free(t3)
end function terms_compare_inverse


subroutine polynomial_term_copy(in, out)
  type(polynomial_term_type), intent(in) :: in
  type(polynomial_term_type), intent(out) ::  out
  if (in%nindex>-1) then
    call polynomial_term_init(in%atindx,in%cell,in%direction,in%ndisp,&
      &                              in%nstrain,out,in%power_disp,&
      &                              in%power_strain,in%strain,in%weight, &
      &                             check=.True., index_coeff=in%index_coeff)
  else
    call polynomial_term_init(in%atindx,in%cell,in%direction,in%ndisp,&
      &                              in%nstrain,out,in%power_disp,&
      &                              in%power_strain,in%strain,in%weight, &
      &                             check=.True.)
  endif
end subroutine polynomial_term_copy

end module m_polynomial_term
!!***
