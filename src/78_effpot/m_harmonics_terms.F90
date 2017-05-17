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
!! Copyright (C) 2010-2017 ABINIT group (AM)
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
 public :: harmonics_terms_setEffectiveCharges
 public :: harmonics_terms_setDynmat
 public :: harmonics_terms_setInternalStrain
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

   integer :: nqpt
!  Number of qpoints

   real(dp) :: epsilon_inf(3,3)                     
!     epsilon_inf(3,3)
!     Dielectric tensor

   real(dp) :: elastic_constants(6,6)
!     elastic_constant(6,6)
!     Elastic tensor Hartree

   real(dp), allocatable :: internal_strain(:,:,:)    
!    internal_strain(6,3,natom)
!    internal strain tensor 

   real(dp), allocatable :: zeff(:,:,:)             
!     zeff(3,3,natom) Effective charges

   type(ifc_type) :: ifcs
!     type with ifcs constants (short + ewald) 
!     also contains the number of cell and the indexes
 
   real(dp), allocatable :: qpoints(:,:) 
!     qph1l(3,nqpt)
!     List of qpoints wavevectors

   real(dp), allocatable :: dynmat(:,:,:,:,:,:)
!     dynmat(2,3,natom,3,natom,nqpt) 
!     dynamical matrix for each q points

   real(dp), allocatable :: phfrq(:,:)
!     phfrq(3*natom,nqpt)
!     array with all phonons frequencies for each q points in Hartree/cm

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
!!      m_effective_potential
!!
!! CHILDREN
!!
!! SOURCE

subroutine harmonics_terms_init(harmonics_terms,ifcs,natom,nrpt,&
&                               dynmat,epsilon_inf,elastic_constants,internal_strain,&
&                               nqpt,phfrq,qpoints,zeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonics_terms_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nrpt
!arrays
 type(ifc_type),intent(in) :: ifcs
 type(harmonics_terms_type), intent(out) :: harmonics_terms
 integer, optional,intent(in) :: nqpt
 real(dp),optional,intent(in) :: epsilon_inf(3,3),dynmat(:,:,:,:,:,:)
 real(dp),optional,intent(in) :: elastic_constants(6,6)
 real(dp),optional,intent(in) :: internal_strain(6,3,natom),zeff(3,3,natom)
 real(dp),optional,intent(in) :: phfrq(:,:),qpoints(:,:)
!Local variables-------------------------------
!scalar
!arrays
 character(len=500) :: msg

! *************************************************************************

 call harmonics_terms_free(harmonics_terms)

! Do some Checks
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

 if (nrpt /= ifcs%nrpt) then
   write(msg, '(3a,i5,a,i5,a)' )&
&   'nrpt must have the same dimension as ifcs.',ch10,&
&   'The number of cell is  ',nrpt,' instead of ',ifcs%nrpt,'.'
   MSG_BUG(msg)
 end if

 if(present(nqpt).and.(.not.present(dynmat).or.&
&                      .not.present(qpoints)   .or.&
&                      .not.present(phfrq)))then
   write(msg, '(a)' )&
&   'nqpt is specified but dynamt,qpoints or phfrq are not.'
   MSG_BUG(msg)
 end if

 if(.not.present(nqpt).and.(present(dynmat).or.&
&                      present(qpoints)   .or.&
&                      present(phfrq)))then
   write(msg, '(a)' )&
&   ' dynamt,qpoints or phfrq are specified but nqpt is not.'
   MSG_BUG(msg)
 end if

!Set number of cell
 harmonics_terms%ifcs%nrpt = nrpt

!Allocation of total ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%atmfrc,(2,3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%atmfrc(:,:,:,:,:,:) = ifcs%atmfrc(:,:,:,:,:,:) 

!Allocation of ewald part of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%ewald_atmfrc,(2,3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%ewald_atmfrc(:,:,:,:,:,:) = ifcs%ewald_atmfrc(:,:,:,:,:,:)

!Allocation of short range part of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%short_atmfrc,(2,3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:,:) = ifcs%short_atmfrc(:,:,:,:,:,:)

!Allocation of cell of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%cell,(3,nrpt))
 harmonics_terms%ifcs%cell(:,:) = ifcs%cell(:,:)

!Allocation of the dynamical matrix
 harmonics_terms%nqpt = zero 
 if(present(nqpt).and.present(dynmat).and.present(qpoints).and.present(phfrq))then
   call harmonics_terms_setDynmat(dynmat,harmonics_terms,natom,nqpt,&
&                                 harmonics_terms%phfrq,harmonics_terms%qpoints)
 end if

!Allocation of the elastic constants
 harmonics_terms%elastic_constants = zero
 if (present(elastic_constants)) then
   harmonics_terms%elastic_constants = elastic_constants
 end if

!Allication of the dielectric tensor
 harmonics_terms%epsilon_inf = zero
 if (present(epsilon_inf)) then
   harmonics_terms%epsilon_inf = epsilon_inf
 end if

!Allocation of Effective charges array 
 ABI_ALLOCATE(harmonics_terms%zeff,(3,3,natom))
 harmonics_terms%zeff = zero 
 if (present(zeff)) then
   call harmonics_terms_setEffectiveCharges(harmonics_terms,natom,zeff)
   harmonics_terms%zeff = zeff
 end if

!Allocation of internal strain tensor 
 ABI_ALLOCATE(harmonics_terms%internal_strain,(6,3,natom))
 harmonics_terms%internal_strain = zero
 if (present(internal_strain)) then
   call harmonics_terms_setInternalStrain(internal_strain,harmonics_terms,natom)
 end if

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
!!      m_effective_potential,m_harmonics_terms
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

  harmonics_terms%nqpt = zero
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

  if(allocated(harmonics_terms%dynmat))then
    harmonics_terms%dynmat=zero
    ABI_DEALLOCATE(harmonics_terms%dynmat)
  end if

  if(allocated(harmonics_terms%phfrq))then
    harmonics_terms%phfrq=zero
    ABI_DEALLOCATE(harmonics_terms%phfrq)
  end if

  if(allocated(harmonics_terms%qpoints))then
    harmonics_terms%qpoints=zero
    ABI_DEALLOCATE(harmonics_terms%qpoints)
  end if
  
  call ifc_free(harmonics_terms%ifcs)
  
end subroutine harmonics_terms_free
!!***

!****f* m_harmonics_terms/harmonics_terms_setInternalStrain
!!
!! NAME
!! harmonics_terms_setInternalStrain
!!
!! FUNCTION
!! Set the internal strain to the harmonics_terms
!!
!! INPUTS
!! natom = number of atoms
!! internal_strain(6,3,natom) = internal strain coupling parameters
!!
!! OUTPUT
!! harmonics_terms = supercell structure with data to be output
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine harmonics_terms_setInternalStrain(internal_strain,harmonics_terms,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonics_terms_setInternalStrain'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: natom
!array
  real(dp),intent(in) :: internal_strain(:,:,:)
  type(harmonics_terms_type), intent(inout) :: harmonics_terms
!Local variables-------------------------------
!scalars
!array
  character(len=500) :: msg

! *************************************************************************

! 0-Checks inputs
  if(natom /= size(internal_strain,3)) then
    write(msg, '(a)' )&
&        ' natom has not the same size internal_strain array. '
    MSG_BUG(msg)
  end if

! 1-deallocate old array
  if(allocated(harmonics_terms%internal_strain))then
    ABI_DEALLOCATE(harmonics_terms%internal_strain)
  end if
   
! 2-allocate and copy the new array
  ABI_ALLOCATE(harmonics_terms%internal_strain,(6,3,natom))
  harmonics_terms%internal_strain(:,:,:) = internal_strain(:,:,:)

end subroutine harmonics_terms_setInternalStrain
!!***


!****f* m_harmonics_terms/harmonics_terms_setEffectiveCharges
!!
!! NAME
!! harmonics_terms_setEffectiveCharges
!!
!! FUNCTION
!! Set the effectives charges to the harmonics_terms
!!
!! INPUTS
!! natom = number of atoms
!! zeff(3,natom) = effective charges
!!
!! OUTPUT
!! harmonics_terms = supercell structure with data to be output
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine harmonics_terms_setEffectiveCharges(harmonics_terms,natom,zeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonics_terms_setEffectiveCharges'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: natom
!array
  real(dp),intent(in) :: zeff(:,:,:)
  type(harmonics_terms_type), intent(inout) :: harmonics_terms
!Local variables-------------------------------
!scalars
!array
  character(len=500) :: msg
! *************************************************************************

! 0-Checks inputs
    if(natom /= size(zeff,3)) then
    write(msg, '(a)' )&
&        ' natom has not the same size zeff array. '
    MSG_BUG(msg)
  end if

! 1-deallocate old array
  if(allocated(harmonics_terms%zeff))then
    ABI_DEALLOCATE(harmonics_terms%zeff)
  end if
   
! 2-allocate and copy the new array
  ABI_ALLOCATE(harmonics_terms%zeff,(3,3,natom))
  harmonics_terms%zeff(:,:,:) = zeff(:,:,:)


end subroutine harmonics_terms_setEffectiveCharges
!!***

!****f* m_harmonics_terms/harmonics_terms_setDynmat
!!
!! NAME
!! harmonics_terms_setDynmat
!!
!! FUNCTION
!! Set the effectives charges to the harmonics_terms
!!
!! INPUTS
!! natom = number of atoms
!! nqpt  = number of qpoints
!! dynmat(2,3,natom,3,natom,nqpt) = dynamical matrix in cartesian coordinates
!! phfrq(3*natom,nqpt) = frequency in hartree 
!! qpoints(3,nqpt) = list of qpoints
!!
!! OUTPUT
!! harmonics_terms = supercell structure with data to be output
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine harmonics_terms_setDynmat(dynmat,harmonics_terms,natom,nqpt,phfrq,qpoints)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonics_terms_setDynmat'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: natom,nqpt
!array
  real(dp),intent(in) :: dynmat(:,:,:,:,:,:)
  real(dp),intent(in) :: qpoints(:,:)
  real(dp),intent(in) :: phfrq(:,:)
  type(harmonics_terms_type), intent(inout) :: harmonics_terms
!Local variables-------------------------------
!scalars
!array
  character(len=500) :: msg
! *************************************************************************

! 0-Checks inputs
    if((natom /= size(dynmat,3)).or.(natom /= size(dynmat,5))) then
    write(msg, '(a)' )&
&        ' natom has not the same size dynmat array. '
    MSG_BUG(msg)
  end if

  if (nqpt /= size(dynmat,6))then
    write(msg, '(a)' )&
&        ' nqpt has not the same size dynmat array. '
    MSG_BUG(msg)
  end if

  if (nqpt /= size(qpoints,2))then
    write(msg, '(a)' )&
&        ' nqpt has not the same size qpoints array. '
    MSG_BUG(msg)
  end if

  if (nqpt /= size(phfrq,2))then
    write(msg, '(a)' )&
&        ' nqpt has not the same size phfrq array. '
    MSG_BUG(msg)
  end if

! 1-deallocate old array
  if(allocated(harmonics_terms%dynmat))then
    ABI_DEALLOCATE(harmonics_terms%dynmat)
  end if

  if(allocated(harmonics_terms%phfrq))then
    ABI_DEALLOCATE(harmonics_terms%phfrq)
  end if

  if(allocated(harmonics_terms%qpoints))then
    ABI_DEALLOCATE(harmonics_terms%qpoints)
  end if
   
! 2-allocate and copy the new array
  harmonics_terms%nqpt = nqpt

  ABI_ALLOCATE(harmonics_terms%dynmat,(2,3,natom,3,natom,nqpt))
  harmonics_terms%dynmat(:,:,:,:,:,:) = dynmat(:,:,:,:,:,:)

  ABI_ALLOCATE(harmonics_terms%phfrq,(3*natom,nqpt))
  harmonics_terms%phfrq(:,:) = phfrq(:,:)

  ABI_ALLOCATE(harmonics_terms%qpoints,(3,nqpt))
  harmonics_terms%qpoints(:,:) = qpoints(:,:)

end subroutine harmonics_terms_setDynmat
!!***

end module m_harmonics_terms
!!***
