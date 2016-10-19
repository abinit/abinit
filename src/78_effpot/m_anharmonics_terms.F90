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

 use m_ifc

 implicit none

 public :: anharmonics_terms_init
 public :: anharmonics_terms_free

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
!! nrpt   = number of cell into ifc
!!
!! OUTPUT
!! eff_pot = effective_potential structure to be initialized
!!
!! PARENTS
!!    multibinit
!!
!! CHILDREN
!!    effective_potential_free
!!
!! SOURCE

subroutine anharmonics_terms_init(anharmonics_terms,natom,nrpt,&
&                                 elastics,internal_strain,phonon_strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anharmonics_terms_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nrpt
 type(anharmonics_terms_type), intent(out) :: anharmonics_terms
 real(dp),optional,intent(in) :: internal_strain(6,6,natom,3)
 real(dp),optional,intent(in) :: elastics(6,6,6)
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

  anharmonics_terms%elastics = zero

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

end subroutine anharmonics_terms_free
!!***

end module m_anharmonics_terms
!!***
