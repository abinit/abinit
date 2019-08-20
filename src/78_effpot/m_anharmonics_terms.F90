!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_anharmonics_terms
!!
!! NAME
!! m_anharmonics_term
!!
!! FUNCTION
!! Module with datatype and tools for the anharmonics terms
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (AM)
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
 use m_abicore
 use m_polynomial_coeff
 use m_ifc, only : ifc_type

 implicit none

 public :: anharmonics_terms_init
 public :: anharmonics_terms_free
 public :: anharmonics_terms_freeCoeffs
 public :: anharmonics_terms_evaluateElastic
 public :: anharmonics_terms_evaluateIFCStrainCoupling
 public :: anharmonics_terms_setCoeffs
 public :: anharmonics_terms_setElastic3rd
 public :: anharmonics_terms_setElastic4th
 public :: anharmonics_terms_setElasticDispCoupling
 public :: anharmonics_terms_setStrainPhononCoupling

!!***

!!****t* defs_abitypes/anharmonics_terms_type
!! NAME
!! anharmonics_terms_type
!!
!! FUNCTION
!! datatype for a effective potential constructed.
!!
!! SOURCE

 type, public :: anharmonics_terms_type

   integer :: ncoeff = 0
!    nterm store the number of coefficients

   logical ::  has_elastic3rd
!   Flag to know if the 3rd derivatives with respect to strain is present

   logical ::  has_elastic4th
!   Flag to know if the 3rd derivatives with respect to strain is present

   logical ::  has_strain_coupling
!   Flag to know if the 3rd derivatives with respect to strain and 2 atom disp is present

   logical ::  has_elastic_displ
!   Flag to know if the 3rd derivatives with respect to 2 strain and 3 atom disp is present

   logical :: bounded
!   True : the model is bounded

   type(polynomial_coeff_type),dimension(:),allocatable :: coefficients
!    array with all the coefficients from  polynomial coefficients

   real(dp) :: elastic3rd(6,6,6)
!    elastic_constant(6,6,6)
!    Elastic tensor Hartree

   real(dp) :: elastic4th(6,6,6,6)
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
!! Initialize anharmonics_terms datatype
!!
!! INPUTS
!! natom  = number of atoms in primitive cell
!! ncoeff = number of coefficient for the fited polynome
!! bounded = optional, flag to now if the model in bounded
!! elastic3rd(6,6,6) = optional,3rd order of the elastic constants
!! elastic4th(6,6,6,6) = optional,4st order of the elastic constants
!! elastic_displacement(6,6,3,natom) = optional,elastic constant - force coupling
!! phonon_strain<type(ifc_type)>(6) = optional,phonon strain couling
!! coeffs<type(polynomial_coeff_type)>(ncoeff) = optional,datatype with polynomial coefficients
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype to be initialized
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_init(anharmonics_terms,natom,ncoeff,&
&                                 bounded,elastic3rd,elastic4th,elastic_displacement,&
&                                 phonon_strain,coeffs)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,ncoeff
 type(anharmonics_terms_type), intent(out) :: anharmonics_terms
 real(dp),optional,intent(in) :: elastic_displacement(6,6,3,natom)
 real(dp),optional,intent(in) :: elastic3rd(6,6,6),elastic4th(6,6,6,6)
 type(polynomial_coeff_type),optional :: coeffs(ncoeff)
 type(ifc_type),optional,intent(in) :: phonon_strain(6)
 logical,optional,intent(in) :: bounded
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

!Set the 3rd order elastic tensor
 anharmonics_terms%elastic4th = zero
 if(present(elastic4th))then
   call anharmonics_terms_setElastic4th(anharmonics_terms,elastic4th)
 end if

!Allocation of 3rd order with respecto to 2 strain and 1 atomic displacement
 if(present(elastic_displacement))then
   call anharmonics_terms_setElasticDispCoupling(anharmonics_terms,natom,elastic_displacement)
 end if

 anharmonics_terms%ncoeff = 0

!Allocation of the coefficient
  if(present(coeffs))then
   if(ncoeff /= size(coeffs))then
     write(msg, '(a)' )&
&        ' ncoeff has not the same size than coeffs array, '
     MSG_BUG(msg)
   end if
   call anharmonics_terms_setCoeffs(coeffs,anharmonics_terms,ncoeff)
 end if

!Set the flag bounded
 if(present(bounded))then
   anharmonics_terms%bounded = bounded
 else
   anharmonics_terms%bounded = .FALSE.
 end if

end subroutine anharmonics_terms_init
!!***

!****f* m_anharmonics_terms/anharmonics_terms_free
!!
!! NAME
!! anharmonics_terms_free
!!
!! FUNCTION
!! deallocate all dynamic memory for this anharmonics_terms datatype
!!
!! INPUTS
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype to be free
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype to be free
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_free(anharmonics_terms)

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
   anharmonics_terms%bounded = .FALSE.

  if(allocated(anharmonics_terms%elastic_displacement)) then
    anharmonics_terms%elastic_displacement=zero
    ABI_DEALLOCATE(anharmonics_terms%elastic_displacement)
  end if

  if(allocated(anharmonics_terms%phonon_strain))then
    do ii = 1,6
       call anharmonics_terms%phonon_strain(ii)%free()
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%phonon_strain)
  end if

  call anharmonics_terms_freeCoeffs(anharmonics_terms)

  anharmonics_terms%elastic3rd = zero


end subroutine anharmonics_terms_free
!!***

!****f* m_anharmonics_terms/anharmonics_terms_freeCoeffs
!!
!! NAME
!! anharmonics_terms_freeCoeffs
!!
!! FUNCTION
!! deallocate all dynamic memory for the coefficients
!! of this  anharmonics_terms datatype
!!
!! INPUTS
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype to be free
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype to be free
!!
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_freeCoeffs(anharmonics_terms)

!Arguments ------------------------------------
!scalars
!array
  type(anharmonics_terms_type), intent(inout) :: anharmonics_terms
!Local variables-------------------------------
!scalars
  integer :: ii
!array

! *************************************************************************

  if(allocated(anharmonics_terms%coefficients))then
    do ii=1,anharmonics_terms%ncoeff
      call polynomial_coeff_free(anharmonics_terms%coefficients(ii))
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%coefficients)
  end if

  anharmonics_terms%ncoeff = 0

end subroutine anharmonics_terms_freeCoeffs
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
!! coeffs(ncoeff)<type(polynomial_coeff_type)> = array with datatype polynomial_coeff_type
!! ncoeff = number of coefficient
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype
!!
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_setCoeffs(coeffs,anharmonics_terms,ncoeff)

 use m_polynomial_coeff

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
    call polynomial_coeff_init(coeffs(ii)%coefficient,coeffs(ii)%nterm,&
&                              anharmonics_terms%coefficients(ii),&
&                              coeffs(ii)%terms,name=coeffs(ii)%name)
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
!! elastics(6,6,6) = 3d order of elastics constant
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_setElastic3rd(anharmonics_terms,elastics)

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
  anharmonics_terms%has_elastic3rd = .FALSE.

! 2-Allocation of the new array
  anharmonics_terms%elastic3rd(:,:,:) = elastics(:,:,:)

! 3-Set the flag
  if(any(abs(anharmonics_terms%elastic3rd)> tol15)) then
    anharmonics_terms%has_elastic3rd = .TRUE.
  end if

end subroutine anharmonics_terms_setElastic3rd
!!***

!****f* m_anharmonics_terms/anharmonics_terms_setElastic4th
!!
!! NAME
!! anharmonics_terms_setElastic4th
!!
!! FUNCTION
!! Set the 4th order derivative of with respect to 4 strain
!!
!! INPUTS
!! elastics = 4th order of elastics constant
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype
!!
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_setElastic4th(anharmonics_terms,elastics)

!Arguments ------------------------------------
!scalars
!array
  type(anharmonics_terms_type),intent(inout) :: anharmonics_terms
  real(dp),intent(in) :: elastics(6,6,6,6)
!Local variables-------------------------------
!scalar
!array
! *************************************************************************

! 1-reinitialise the previous value
  anharmonics_terms%elastic4th(:,:,:,:) = zero
  anharmonics_terms%has_elastic4th = .FALSE.

! 2-Allocation of the new array
  anharmonics_terms%elastic4th(:,:,:,:) = elastics(:,:,:,:)

! 3-Set the flag
  if(any(abs(anharmonics_terms%elastic4th)> tol15)) then
    anharmonics_terms%has_elastic4th = .TRUE.
  end if

end subroutine anharmonics_terms_setElastic4th
!!***


!****f* m_anharmonics_terms/anharmonics_terms_setStrainPhononCoupling
!!
!! NAME
!! anharmonics_terms_setStrainPhononCoupling
!!
!! FUNCTION
!! Set the strain-phonon coupling
!!
!! INPUTS
!! strain_phonon(6)<type(ifc_type) = strain-phonon coupling
!! natom = number of atoms
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_setStrainPhononCoupling(anharmonics_terms,natom,phonon_strain)

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
    if(natom /= size(phonon_strain(ii)%atmfrc,2).or.&
&      phonon_strain(ii)%nrpt < 0)then
      write(msg, '(a)' )&
&        ' natom or/and nrpt have not the same size than phonon_strain array. '
      MSG_BUG(msg)
    end if
  end do

! 1-reinitialise the previous value
  anharmonics_terms%has_strain_coupling  = .FALSE.
  if(allocated(anharmonics_terms%phonon_strain))then
    do ii = 1,6
       call anharmonics_terms%phonon_strain(ii)%free()
    end do
    ABI_DATATYPE_DEALLOCATE(anharmonics_terms%phonon_strain)
  end if

! 2-Allocation of the new array and filling
 ABI_DATATYPE_ALLOCATE(anharmonics_terms%phonon_strain,(6))
 do ii = 1,6
   nrpt = phonon_strain(ii)%nrpt
   ABI_ALLOCATE(anharmonics_terms%phonon_strain(ii)%atmfrc,(3,natom,3,natom,nrpt))
   ABI_ALLOCATE(anharmonics_terms%phonon_strain(ii)%cell,(3,nrpt))
   anharmonics_terms%phonon_strain(ii)%nrpt   = phonon_strain(ii)%nrpt
   anharmonics_terms%phonon_strain(ii)%atmfrc(:,:,:,:,:) = phonon_strain(ii)%atmfrc(:,:,:,:,:)
   anharmonics_terms%phonon_strain(ii)%cell(:,:)   = phonon_strain(ii)%cell(:,:)

!3-Set the flag
   if(any(abs(anharmonics_terms%phonon_strain(ii)%atmfrc)> tol15)) then
     anharmonics_terms%has_strain_coupling  = .TRUE.
!  If there is no value inside the array,
!  We don't need to store it
   else
     ABI_DEALLOCATE(anharmonics_terms%phonon_strain(ii)%atmfrc)
     ABI_DEALLOCATE(anharmonics_terms%phonon_strain(ii)%cell)
     anharmonics_terms%phonon_strain(ii)%nrpt = 0
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
!! Set the Elastic displacement coupling
!!
!! INPUTS
!! elastic_displacement(6,6,3,natom) = elastic displacement coupling
!! natom = number of atom
!!
!! OUTPUT
!! anharmonics_terms<type(anharmonics_terms_type)> = anharmonics_terms datatype
!!
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_setElasticDispCoupling(anharmonics_terms,natom,elastic_displacement)

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
  anharmonics_terms%has_elastic_displ = .FALSE.
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
    ABI_DEALLOCATE(anharmonics_terms%elastic_displacement)
  end if

end subroutine anharmonics_terms_setElasticDispCoupling
!!***

!!****f* m_effective_potential/anharmonics_terms_evaluateElastic
!! NAME
!!  anharmonics_terms_evaluateElastic
!!
!! FUNCTION
!! Compute the energy, stresses and forces related to the application of strain
!!
!! INPUTS
!! disp(3,natom_sc) = atomics displacement between configuration and the reference
!! natom = number of atom in the supercell
!! natom_uc = number of atom in the unit cell
!! ncell  = number of cell
!! strain(6) =  strain to apply
!! elastic3rd(6,6,6) = 3 order derivatives with respect to to 3 strain
!! elastic4th(6,6,66,) = 4 order derivatives with respect to to 4 strain
!! elastic_displacement(6,6,3,natom) = 3 order derivatives with respect to 2 strain and 1 Atom disp
!!
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!!   strten(6) = contribution to the stress tensor
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE
!!
subroutine anharmonics_terms_evaluateElastic(disp,energy,fcart,natom,natom_uc,ncell,strten,strain,&
&                                            elastic3rd,elastic4th,elastic_displacement)

 real(dp),intent(out):: energy
 integer, intent(in) :: natom,natom_uc,ncell
! array
 real(dp),optional,intent(in) :: elastic3rd(6,6,6),elastic4th(6,6,6,6)
 real(dp),optional,intent(in) :: elastic_displacement(6,6,3,natom)
 real(dp),intent(out):: strten(6)
 real(dp),intent(out):: fcart(3,natom)
 real(dp),intent(in) :: disp(3,natom)
 real(dp),intent(in) :: strain(6)

 !Local variables-------------------------------
! scalar
 integer :: ia,ii,mu,alpha,beta,gamma,delta,d1,d2
 real(dp):: cijk
 logical :: has_elastic3rd,has_elastic4th,has_elastic_displ
! array
! *************************************************************************

!Reset output and flags
 energy = zero
 fcart = zero
 strten = zero
 has_elastic3rd    = .FALSE.
 has_elastic4th    = .FALSE.
 has_elastic_displ = .FALSE.
 d1=0;d2=0

!Set the flags
 if(present(elastic3rd)) has_elastic3rd = .TRUE.
 if(present(elastic4th)) then
   has_elastic4th = .TRUE.
   d1=1;d2=6
 end if
 if(present(elastic_displacement)) has_elastic_displ = .TRUE.

!1-Treat 3rd order elastic constants
 if (has_elastic3rd.or.has_elastic4th) then
   do alpha=1,6
     do beta=1,6
       do gamma=1,6
         cijk = ncell*elastic3rd(alpha,beta,gamma)
!        Accumulate energy
         energy = energy + sixth*cijk*strain(alpha)*strain(beta)*strain(gamma)
!        Accumulate stresses contributions
         strten(alpha)=strten(alpha)+ half*cijk*strain(beta)*strain(gamma)
         do delta=d1,d2
           cijk = ncell*elastic4th(alpha,beta,gamma,delta)
!          Accumulate energy
           energy = energy + (1/24.)*cijk*strain(alpha)*strain(beta)*&
&                                                  strain(gamma)*strain(delta)
!          Accumulate stresses contributions
           strten(alpha)=strten(alpha)+ sixth*cijk*strain(beta)*strain(gamma)*&
&                                                           strain(delta)
         end do
       end do
     end do
   end do
 end if

!2-Part due to the internat strain
 if(has_elastic_displ)then
   ii = 1
   do ia = 1,natom
     do mu = 1,3
       do beta=1,6
         do alpha=1,6
           cijk = elastic_displacement(alpha,beta,mu,ii)
!          Accumulte for this atom
           energy = energy + sixth*cijk*strain(alpha)*strain(beta)*disp(mu,ia)
           fcart(mu,ia) = fcart(mu,ia)   +  half*cijk*strain(alpha)*strain(beta)
           strten(alpha) = strten(alpha) +  half*cijk*strain(beta)*disp(mu,ia)
         end do
       end do
     end do
     ii = ii +1
!    Reset to 1 if the number of atoms is superior than in the initial cell
     if(ii==natom_uc+1) ii = 1
   end do
 end if

end subroutine anharmonics_terms_evaluateElastic
!!***

!!****f* m_anharmonics_terms/anharmonics_terms_evaluateIFCStrainCoupling
!! NAME
!!  anharmonics_terms_evaluateIFCStrainCoupling
!!
!! FUNCTION
!!  This fonction compute the harmonic part of the energy
!!  of the supercell in the eff_pot
!! INPUTS
!!  strain_phonon(6)<type(ifc_type) = strain-phonon coupling
!!  disp(3,natom_sc) = atomics displacement between configuration and the reference
!!  natom = number of atoms in the supercell
!!  natom_uc = number of atoms in the unit cell
!!  sc_size(3) = size of the supercell
!!  cells(ncell) = number of the cells into the supercell (1,2,3,4,5)
!!  ncell  = total number of cell to treat
!!  index_cells(3,ncell) = indexes of the cells into  supercell (-1 -1 -1 ,...,1 1 1)
!!  comm=MPI communicator
!!
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!!   strten(6) = contribution to the stress tensor
!!
!! PARENT
!!   effective_potential_evaluate
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      getpbcindexes_supercell,xmpi_sum
!!
!! SOURCE

subroutine anharmonics_terms_evaluateIFCStrainCoupling(phonon_strain,disp,energy,fcart,natom,natom_uc,&
&                                                      sc_size,strain,strten,cells,ncell,&
&                                                      index_cells,comm)

!Arguments -------------------------------
! scalars
  real(dp),intent(out) :: energy
  integer,intent(in) :: natom,natom_uc,ncell
  integer,intent(in) :: comm
! array
  integer,intent(in) ::   cells(ncell),index_cells(ncell,3)
  integer,intent(in) :: sc_size(3)
  type(ifc_type),intent(in) :: phonon_strain(6)
  real(dp),intent(in) :: disp(3,natom)
  real(dp),intent(out) :: fcart(3,natom)
  real(dp),intent(out) :: strten(6)
  real(dp),intent(in) :: strain(6)
!Local variables-------------------------------
! scalar
  integer :: alpha
  integer :: i1,i2,i3,ia,ib,icell,ii
  integer :: irpt,jj,kk,ll,mu,nu
  integer :: ierr
  real(dp):: ifc
! array
  integer :: cell_atom2(3)
  character(500) :: msg

! *************************************************************************

  if (any(sc_size <= 0)) then
    write(msg,'(a,a)')' sc_size can not be inferior or equal to zero'
    MSG_ERROR(msg)
  end if

! Initialisation of variables
  energy   = zero
  fcart(:,:) = zero
  strten(:) = zero

  do icell = 1,ncell
    ii = (cells(icell)-1)*natom_uc
    i1=index_cells(icell,1); i2=index_cells(icell,2); i3=index_cells(icell,3)
    do alpha=1,6
      do irpt = 1,phonon_strain(alpha)%nrpt
!       get the cell of atom2  (0 0 0, 0 0 1...)
        cell_atom2(1) =  i1 + phonon_strain(alpha)%cell(1,irpt)
        cell_atom2(2) =  i2 + phonon_strain(alpha)%cell(2,irpt)
        cell_atom2(3) =  i3 + phonon_strain(alpha)%cell(3,irpt)
        call getPBCIndexes_supercell(cell_atom2(1:3),sc_size(1:3))
!       index of the second atom in the displacement array
        jj = (cell_atom2(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&            (cell_atom2(2)-1)*sc_size(3)*natom_uc+&
&            (cell_atom2(3)-1)*natom_uc
        do ib = 1, natom_uc
          ll = jj + ib
          do nu=1,3
            do ia = 1, natom_uc
              kk = ii + ia
              do mu=1,3
                ifc = phonon_strain(alpha)%atmfrc(mu,ia,nu,ib,irpt)
!               accumule energy
                energy =  energy + sixth*strain(alpha)*disp(mu,kk)*disp(nu,ll)*ifc
!               accumule forces
                fcart(mu,kk) = fcart(mu,kk) + half*strain(alpha)*disp(nu,ll)*ifc
!               accumule stresses
                strten(alpha) = strten(alpha) + half*disp(mu,kk)*disp(nu,ll)*ifc
              end do
            end do
          end do
        end do
      end do
    end do
  end do

! MPI_SUM
  call xmpi_sum(energy, comm, ierr)
  call xmpi_sum(fcart , comm, ierr)
  call xmpi_sum(strten, comm, ierr)

end subroutine anharmonics_terms_evaluateIFCStrainCoupling
!!***
end module m_anharmonics_terms
!!***
