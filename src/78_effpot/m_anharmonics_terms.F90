!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_anharmonics_terms
!!
!! NAME
!! m_anharmonics_term
!!
!! FUNCTION
!! Module with structure and tools for the anharmonics terms
!! Compute strain phonon coupling by finite differences
!! Return the effective_potential with the third order
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

module m_anharmonics_terms

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_polynomial_coeff
 use m_ifc

 implicit none

 public :: anharmonics_terms_init
 public :: anharmonics_terms_free
 public :: anharmonics_terms_setCoeffs

!!***

!!****t* defs_abitypes/anharmonics_terms_type
!! NAME
!! anharmonics_terms_type
!!
!! FUNCTION
!! structure for a effective potential constructed.
!!
!! SOURCE

 type, public :: anharmonics_terms_type

   integer :: ncoeff
!  nterm store the number of coefficients

   type(polynomial_coeff_type),dimension(:),allocatable :: coefficients
!  array with all the coefficients from fited polynome

   real(dp) :: elastics(6,6,6)
!     elastic_constant(6,6,6)
!     Elastic tensor Hartree

   real(dp), allocatable :: internal_strain(:,:,:,:)    
!    internal_strain(6,6,natom,3)
!    internal strain tensor 

   type(ifc_type),dimension(:),allocatable :: phonon_strain
!   Array of ifc with phonon_strain coupling for each strain 

 end type anharmonics_terms_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_anharmonics_terms/anharmonics_terms_init
!!
!! NAME
!! anharmonics_terms_init
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, qpoint chosen, and atoms
!!
!! INPUTS
!! natom  = number of atoms in primitive cell
!! ncoeff = number of coefficient for the fited polynome 
!! nrpt   = number of cell into ifc
!!
!! OUTPUT
!! anharmonics_terms = anharmonics_terms structure to be initialized
!!
!! PARENTS
!!    multibinit
!!
!! CHILDREN
!!    effective_potential_free
!!
!! SOURCE

subroutine anharmonics_terms_init(anharmonics_terms,natom,ncoeff,nrpt,&
&                                 elastics,internal_strain,phonon_strain,coeffs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,ncoeff,nrpt
 type(anharmonics_terms_type), intent(out) :: anharmonics_terms
 real(dp),optional,intent(in) :: internal_strain(6,6,natom,3)
 real(dp),optional,intent(in) :: elastics(6,6,6)
 type(polynomial_coeff_type),optional :: coeffs(ncoeff)
 type(ifc_type),intent(in) :: phonon_strain(6)
!arrays
!Local variables-------------------------------
!scalar
 integer :: ii
!arrays
 character(len=500) :: msg

! *************************************************************************

 call anharmonics_terms_free(anharmonics_terms)

! Check the number of atoms
 if (natom < 1) then
   write(msg, '(a,a,a,i10,a)' )&
&   'The cell must have at least one atom.',ch10,&
&   'The number of atom is  ',natom,'.'
   MSG_BUG(msg)
 end if

 if (nrpt < 1) then
   write(msg, '(a,a,a,i10,a)' )&
&   'The cell must have at least one rpt point.',ch10,&
&   'The number of rpt points is  ',nrpt,'.'
   MSG_BUG(msg)
 end if

 anharmonics_terms%elastics = zero
 if(present(elastics))then
   anharmonics_terms%elastics = elastics
 end if
 
!Allocation of phonon strain coupling array (3rd order)
 ABI_DATATYPE_ALLOCATE(anharmonics_terms%phonon_strain,(6))
 do ii = 1,6
   ABI_ALLOCATE(anharmonics_terms%phonon_strain(ii)%atmfrc,(2,3,natom,3,natom,nrpt))
   ABI_ALLOCATE(anharmonics_terms%phonon_strain(ii)%cell,(3,nrpt))
   anharmonics_terms%phonon_strain(ii)%nrpt   = nrpt
   anharmonics_terms%phonon_strain(ii)%atmfrc = zero
   anharmonics_terms%phonon_strain(ii)%cell   = zero
 end do

!Allocation of elastic 3rd order
 ABI_ALLOCATE(anharmonics_terms%internal_strain,(6,6,natom,3))
 anharmonics_terms%internal_strain = zero

!Allocation of the coefficient
  if(present(coeffs))then
   if(ncoeff /= size(coeffs))then
     write(msg, '(a)' )&
&        ' ncoeff has not the same size than coeffs array, '
     MSG_BUG(msg)
   end if
   call anharmonics_terms_setCoeffs(coeffs,anharmonics_terms,ncoeff)
 end if

end subroutine anharmonics_terms_init 
!!***


!****f* m_anharmonics_terms/anharmonics_terms_free
!!
!! NAME
!! anharmonics_terms_free
!!
!! FUNCTION
!! deallocate all dynamic memory for this anharmonics_terms structure
!!
!! INPUTS
!!
!! OUTPUT
!! anharmonics_terms = structure with anharmonics terms
!!
!! PARENTS
!!   anharmonics_terms_init
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine anharmonics_terms_free(anharmonics_terms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_free'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  type(anharmonics_terms_type), intent(inout) :: anharmonics_terms

!Local variables-------------------------------
!scalars
  integer :: ii
!array

! *************************************************************************

  if(allocated(anharmonics_terms%internal_strain)) then
    anharmonics_terms%internal_strain=zero
    ABI_DEALLOCATE(anharmonics_terms%internal_strain)
  end if

  if(allocated(anharmonics_terms%phonon_strain))then
    do ii = 1,6
       call ifc_free(anharmonics_terms%phonon_strain(ii))
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%phonon_strain)
  end if

  if(allocated(anharmonics_terms%coefficients))then
    do ii=1,anharmonics_terms%ncoeff
      call polynomial_coeff_free(anharmonics_terms%coefficients(ii))
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%coefficients)
  end if

  anharmonics_terms%elastics = zero
  anharmonics_terms%ncoeff = zero

end subroutine anharmonics_terms_free
!!***

!****f* m_anharmonics_terms/anharmonics_terms_setCoeffs
!!
!! NAME
!! anharmonics_terms_setCoeffs
!!
!! FUNCTION
!! Set the coefficients
!!
!! INPUTS
!! coeffs = polynomial_coeff_type
!! anharmonics = anharmonics type
!! ncoeff = number of coefficient
!!
!! OUTPUT
!! anharmonics = set the coefficient from the fited polynome 
!!
!!
!! PARENTS
!!   multibinit
!!
!! CHILDREN
!!   wrtout
!!
!! SOURCE
 
subroutine anharmonics_terms_setCoeffs(coeffs,anharmonics,ncoeff)

 use m_polynomial_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_setCoeffs'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: ncoeff
!array
  type(anharmonics_terms_type),intent(inout) :: anharmonics
  type(polynomial_coeff_type),intent(in) :: coeffs(ncoeff)
!Local variables-------------------------------
!scalar
  integer :: ii
  character(len=500) :: msg
!array
! *************************************************************************

  if(ncoeff /= size(coeffs))then
    write(msg, '(a)' )&
&        ' ncoeff has not the same size than coeffs array, '
    MSG_BUG(msg)
  end if

! 1-deallocation of the previous value
  if(allocated(anharmonics%coefficients))then
    do ii=1,anharmonics%ncoeff
      call polynomial_coeff_free(anharmonics%coefficients(ii))
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics%coefficients)
  end if

! Allocation of the new array
  anharmonics%ncoeff = ncoeff  
  ABI_DATATYPE_ALLOCATE(anharmonics%coefficients,(ncoeff))
  do ii=1,anharmonics%ncoeff
    call polynomial_coeff_init(coeffs(ii)%coefficient,coeffs(ii)%name,coeffs(ii)%nterm,&
&                              anharmonics%coefficients(ii),coeffs(ii)%terms)
  end do

end subroutine anharmonics_terms_setCoeffs
!!***

end module m_anharmonics_terms
!!***
