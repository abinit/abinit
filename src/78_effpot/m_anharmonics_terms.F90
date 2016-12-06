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
 public :: anharmonics_terms_setElastic3rd
 public :: anharmonics_terms_setElasticDispCoupling
 public :: anharmonics_terms_setStrainPhononCoupling

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

   integer :: ncoeff = 0
!    nterm store the number of coefficients

   logical ::  has_elastic3rd
!   Flag to know if the 3rd derivatives with respect to strain is present

   logical ::  has_strain_coupling
!   Flag to know if the 3rd derivatives with respect to strain and 2 atom disp is present

   logical ::  has_elastic_displ
!   Flag to know if the 3rd derivatives with respect to 2 strain and 3 atom disp is present

   type(polynomial_coeff_type),dimension(:),allocatable :: coefficients
!    array with all the coefficients from fited polynome

   real(dp) :: elastic3rd(6,6,6)
!    elastic_constant(6,6,6)
!    Elastic tensor Hartree

   real(dp), allocatable :: elastic_displacement(:,:,:,:)    
!    elastic_displacement(6,6,3,natom)
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

subroutine anharmonics_terms_init(anharmonics_terms,natom,ncoeff,&
&                                 elastic3rd,elastic_displacement,phonon_strain,coeffs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,ncoeff
 type(anharmonics_terms_type), intent(out) :: anharmonics_terms
 real(dp),optional,intent(in) :: elastic_displacement(6,6,3,natom)
 real(dp),optional,intent(in) :: elastic3rd(6,6,6)
 type(polynomial_coeff_type),optional :: coeffs(ncoeff)
 type(ifc_type),optional,intent(in) :: phonon_strain(6)
!arrays
!Local variables-------------------------------
!scalar
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

!Allocation of phonon strain coupling array (3rd order)
 if(present(phonon_strain)) then
   call anharmonics_terms_setStrainPhononCoupling(anharmonics_terms,natom,phonon_strain)
 end if

!Set the 3rd order elastic tensor
 anharmonics_terms%elastic3rd = zero
 if(present(elastic3rd))then
   call anharmonics_terms_setElastic3rd(anharmonics_terms,elastic3rd)
 end if

!Allocation of 3rd order with respecto to 2 strain and 1 atomic displacement
 if(present(elastic_displacement))then
   call anharmonics_terms_setElasticDispCoupling(anharmonics_terms,natom,elastic_displacement)
 end if

 anharmonics_terms%ncoeff = zero
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

   anharmonics_terms%has_elastic3rd = .FALSE.
   anharmonics_terms%has_strain_coupling  = .FALSE.
   anharmonics_terms%has_elastic_displ = .FALSE.

  if(allocated(anharmonics_terms%elastic_displacement)) then
    anharmonics_terms%elastic_displacement=zero
    ABI_DEALLOCATE(anharmonics_terms%elastic_displacement)
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

  anharmonics_terms%ncoeff = zero

  anharmonics_terms%elastic3rd = zero


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
 
subroutine anharmonics_terms_setCoeffs(coeffs,anharmonics_terms,ncoeff)

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
  type(anharmonics_terms_type),intent(inout) :: anharmonics_terms
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
  if(allocated(anharmonics_terms%coefficients))then
    do ii=1,anharmonics_terms%ncoeff
      call polynomial_coeff_free(anharmonics_terms%coefficients(ii))
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%coefficients)
  end if

! Allocation of the new array
  anharmonics_terms%ncoeff = ncoeff  
  ABI_DATATYPE_ALLOCATE(anharmonics_terms%coefficients,(ncoeff))
  do ii=1,anharmonics_terms%ncoeff
    call polynomial_coeff_init(coeffs(ii)%coefficient,coeffs(ii)%name,coeffs(ii)%nterm,&
&                              anharmonics_terms%coefficients(ii),coeffs(ii)%terms)
  end do

end subroutine anharmonics_terms_setCoeffs
!!***

!****f* m_anharmonics_terms/anharmonics_terms_setElastic3rd
!!
!! NAME
!! anharmonics_terms_setElastic3rd
!!
!! FUNCTION
!! Set the 3rd order derivative of with respect to 3 strain
!!
!! INPUTS
!! elastics = 3d order of elastics constant
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
 
subroutine anharmonics_terms_setElastic3rd(anharmonics_terms,elastics)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_setElastic3rd'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  type(anharmonics_terms_type),intent(inout) :: anharmonics_terms
  real(dp),intent(in) :: elastics(6,6,6)
!Local variables-------------------------------
!scalar
!array
! *************************************************************************

! 1-reinitialise the previous value
  anharmonics_terms%elastic3rd(:,:,:) = zero

! 2-Allocation of the new array
  anharmonics_terms%elastic3rd(:,:,:) = elastics(:,:,:)

! 3-Set the flag
  if(any(abs(anharmonics_terms%elastic3rd)> tol15)) then
    anharmonics_terms%has_elastic3rd = .TRUE. 
  end if

end subroutine anharmonics_terms_setElastic3rd
!!***

!****f* m_anharmonics_terms/anharmonics_terms_setStrainPhononCoupling
!!
!! NAME
!! anharmonics_terms_setStrainPhononCoupling
!!
!! FUNCTION
!! Set the coefficients
!!
!! INPUTS
!! strain_phonon = (size 6) array of type ifc
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
 
subroutine anharmonics_terms_setStrainPhononCoupling(anharmonics_terms,natom,phonon_strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_setStrainPhononCoupling'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(anharmonics_terms_type),intent(inout) :: anharmonics_terms
  type(ifc_type),intent(in) :: phonon_strain(6)
!Local variables-------------------------------
!scalar
  integer :: ii,nrpt
  character(500) :: msg
!array
! *************************************************************************

! 1-Do some check
  do ii=1,6
    if(natom /= size(phonon_strain(ii)%atmfrc,3).or.&
&      phonon_strain(ii)%nrpt < 0)then
      write(msg, '(a)' )&
&        ' natom or/and nrpt have not the same size than phonon_strain array. '
      MSG_BUG(msg)
    end if
  end do

! 1-reinitialise the previous value
  if(allocated(anharmonics_terms%phonon_strain))then
    do ii = 1,6
       call ifc_free(anharmonics_terms%phonon_strain(ii))
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%phonon_strain)
  end if

! 2-Allocation of the new array and filling
 ABI_DATATYPE_ALLOCATE(anharmonics_terms%phonon_strain,(6))
 do ii = 1,6
   nrpt = phonon_strain(ii)%nrpt
   ABI_ALLOCATE(anharmonics_terms%phonon_strain(ii)%atmfrc,(2,3,natom,3,natom,nrpt))
   ABI_ALLOCATE(anharmonics_terms%phonon_strain(ii)%cell,(3,nrpt))
   anharmonics_terms%phonon_strain(ii)%nrpt   = phonon_strain(ii)%nrpt
   anharmonics_terms%phonon_strain(ii)%atmfrc(:,:,:,:,:,:) = phonon_strain(ii)%atmfrc(:,:,:,:,:,:)
   anharmonics_terms%phonon_strain(ii)%cell(:,:)   = phonon_strain(ii)%cell(:,:)

!3-Set the flag
   if(any(abs(anharmonics_terms%phonon_strain(ii)%atmfrc)> tol15)) then
     anharmonics_terms%has_strain_coupling  = .TRUE.
!  If there is no value inside the array,
!  We don't need to store it
   else
     anharmonics_terms%has_strain_coupling  = .FALSE.
     ABI_DEALLOCATE(anharmonics_terms%phonon_strain(ii)%atmfrc)
     ABI_DEALLOCATE(anharmonics_terms%phonon_strain(ii)%cell)
     anharmonics_terms%phonon_strain(ii)%nrpt = zero
   end if
 end do


end subroutine anharmonics_terms_setStrainPhononCoupling
!!***

!****f* m_anharmonics_terms/anharmonics_terms_setElasticDispCoupling
!!
!! NAME
!! anharmonics_terms_setElasticDispCoupling
!!
!! FUNCTION
!! Set the coefficients
!!
!! INPUTS
!! strain_phonon = (size 6) array of type ifc
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
 
subroutine anharmonics_terms_setElasticDispCoupling(anharmonics_terms,natom,elastic_displacement)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_setElasticDispCoupling'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(anharmonics_terms_type),intent(inout) :: anharmonics_terms
  real(dp),intent(in) :: elastic_displacement(6,6,3,natom)
!Local variables-------------------------------
!scalar
  character(500) :: msg
!array
! *************************************************************************

! 1-Do some check
  if(natom /= size(elastic_displacement,4)) then
    write(msg, '(a)' )&
&        ' natom has not the same size elastic_displacement array. '
    MSG_BUG(msg)
  end if

! 1-reinitialise the previous value
  if(allocated(anharmonics_terms%elastic_displacement))then
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%elastic_displacement)
  end if

! 2-Allocation of the new array and filling
  ABI_ALLOCATE(anharmonics_terms%elastic_displacement,(6,6,3,natom))
  anharmonics_terms%elastic_displacement(:,:,:,:) = elastic_displacement(:,:,:,:)

! 3-Set the flag
  if(any(abs(anharmonics_terms%elastic_displacement)> tol15)) then
    anharmonics_terms%has_elastic_displ = .TRUE.
  else
!   If there is no value inside the array,
!   We don't need to store it
    anharmonics_terms%has_elastic_displ = .FALSE.
    ABI_DEALLOCATE(anharmonics_terms%elastic_displacement)
  end if

end subroutine anharmonics_terms_setElasticDispCoupling
!!***

end module m_anharmonics_terms
!!***
