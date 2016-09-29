!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_harmonics_terms
!!
!! NAME
!! m_harmonics_term
!!
!! FUNCTION
!! Module with structure and tools for the harmonics terms
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

module m_harmonics_terms

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_ifc

 implicit none

 public :: harmonics_terms_init
 public :: harmonics_terms_free

!!***

!!****t* defs_abitypes/harmonics_terms_type
!! NAME
!! harmonics_terms_type
!!
!! FUNCTION
!! structure for harmonic part of effective potential.
!!
!! SOURCE

 type, public :: harmonics_terms_type

   real(dp) :: epsilon_inf(3,3)                     
!     epsilon_inf(3,3)
!     Dielectric tensor

   real(dp) :: elastic_constants(6,6)
!     elastic_constant(6,6)
!     Elastic tensor Hartree

   real(dp), allocatable :: internal_strain(:,:,:)    
!    internal_strain(6,natom,3)
!    internal strain tensor 

   real(dp), allocatable :: zeff(:,:,:)             
!     zeff(3,3,natom) Effective charges

   type(ifc_type) :: ifcs
!   type with ifcs constants (short + ewald) 
!   also contains the number of cell and the indexes

 end type harmonics_terms_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_harmonics_terms/harmonics_terms_init
!!
!! NAME
!! harmonics_terms_init
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
!!    epigene
!!
!! CHILDREN
!!    effective_potential_free
!!
!! SOURCE

subroutine harmonics_terms_init(harmonics_terms,natom,nrpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonics_terms_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nrpt
 type(harmonics_terms_type), intent(out) :: harmonics_terms
!arrays
!Local variables-------------------------------
!scalar
!arrays
 character(len=500) :: msg

! *************************************************************************

 call harmonics_terms_free(harmonics_terms)

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

 harmonics_terms%elastic_constants = zero
 harmonics_terms%epsilon_inf = zero

!Allocation of Effective charges array 
 ABI_ALLOCATE(harmonics_terms%zeff,(3,3,natom))
 harmonics_terms%zeff = zero 

!Allocation of internal strain tensor 
 ABI_ALLOCATE(harmonics_terms%internal_strain,(6,natom,3))
 harmonics_terms%internal_strain = zero

!Allocation of total ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%atmfrc,(2,3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%atmfrc = zero 

!Allocation of ewald part of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%ewald_atmfrc,(2,3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%ewald_atmfrc = zero 

!Allocation of short range part of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%short_atmfrc,(2,3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%short_atmfrc = zero 

!Allocation of cell of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%cell,(nrpt,3))
 harmonics_terms%ifcs%short_atmfrc = zero 

end subroutine harmonics_terms_init 
!!***


!****f* m_harmonics_terms/harmonics_terms_free
!!
!! NAME
!! harmonics_terms_free
!!
!! FUNCTION
!! deallocate all dynamic memory for this harmonic structure
!!
!! INPUTS
!!
!! OUTPUT
!! harmonics_terms = supercell structure with data to be output
!!
!! PARENTS
!!   harmonics_terms_init
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine harmonics_terms_free(harmonics_terms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonics_terms_free'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  type(harmonics_terms_type), intent(inout) :: harmonics_terms
!Local variables-------------------------------
!scalars
!array

! *************************************************************************

  harmonics_terms%elastic_constants = zero
  harmonics_terms%epsilon_inf       = zero

  if(allocated(harmonics_terms%zeff))then
    harmonics_terms%zeff=zero
    ABI_DEALLOCATE(harmonics_terms%zeff)
  end if

  if(allocated(harmonics_terms%internal_strain)) then
    harmonics_terms%internal_strain=zero
    ABI_DEALLOCATE(harmonics_terms%internal_strain)
  end if
  
  call ifc_free(harmonics_terms%ifcs)
  
end subroutine harmonics_terms_free
!!***

end module m_harmonics_terms
!!***
