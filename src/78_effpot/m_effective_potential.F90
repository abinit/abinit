!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_effective_potential
!!
!! NAME
!! m_effective_potential
!!
!! FUNCTION
!! Module for using a effective potential
!! Container type is defined, and destruction, print subroutines
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (AM)
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

module m_effective_potential

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_strain
 use m_ifc
 use m_io_tools, only : open_file
 use m_phonon_supercell
 use m_phonons
 use m_ddb
 use m_anharmonics_terms
 use m_harmonics_terms
 use m_special_funcs,  only : factorial
 use m_xmpi
 use m_copy,           only : alloc_copy
 use m_crystal,        only : crystal_t, crystal_init, crystal_free, crystal_t,crystal_print
 use m_anaddb_dataset, only : anaddb_dataset_type, anaddb_dtset_free, outvars_anaddb, invars9
 use m_dynmat,         only : make_bigbox,q0dy3_apply, q0dy3_calc, dfpt_phfrq

 implicit none

 public :: effective_potential_applySumRule
 public :: effective_potential_distributeResidualForces
 public :: effective_potential_effpot2ddb
 public :: effective_potential_effpot2dynmat
 public :: effective_potential_free
 public :: effective_potential_freempi_supercell
 public :: effective_potential_generateDipDip
 public :: effective_potential_getDeltaEnergy
 public :: effective_potential_evaluate
 public :: effective_potential_getForces
 public :: effective_potential_init
 public :: effective_potential_initmpi_supercell
 public :: effective_potential_print
 public :: effective_potential_printPDOS
 public :: effective_potential_printSupercell
 public :: effective_potential_setCoeffs
 public :: effective_potential_setElastic3rd
 public :: effective_potential_setElastic4rd
 public :: effective_potential_setElasticDispCoupling
 public :: effective_potential_setStrainPhononCoupling
 public :: effective_potential_writeAbiInput
 public :: effective_potential_writeXML
 public :: effective_potential_writeNETCDF
 public :: OPERATOR(==)


 private :: index_periodic
 private :: find_bound
!!***

!!****t* defs_abitypes/effective_potential_type
!! NAME
!! effective_potential_type
!!
!! FUNCTION
!! structure for a effective potential constructed.
!!
!! SOURCE

 type, public :: effective_potential_type

   character(len=fnlen) :: name
!     Name of the molecule (CaTiO3,...)

   real(dp):: energy
!     Energy of the system (Hatree)

   real(dp) :: internal_stress(6)
!     internal_stress(6)
!     stress tensor of the structure

   real(dp) :: external_stress(6)
!     internal_stress(6)
!     stress tensor of the structure

   logical :: has_strainCoupling
!     True : the 3rd order derivative is computed

   logical :: has_strain
!     True : strain is apply

   real(dp), allocatable :: forces(:,:)
!    forces(3,natom)
!    initial forces in the system

   type(crystal_t) :: crystal
!    crystal type
!    contains all information of the crystal

   type(strain_type) :: strain
!     strain type
!     Type wich containt strain information (related to reference structure)

   type(harmonics_terms_type) :: harmonics_terms
!     type with all information for anharmonics terms

   type(anharmonics_terms_type) :: anharmonics_terms
!     type with all information for anharmonics terms

   type(supercell_type) :: supercell
!     super cell type
!     Store all the information of the suppercell

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over the supercell
   integer :: comm_supercell
   ! local Communicator over all processors treating the same cell

   integer :: me_supercell
   ! Index of my processor in the comm. over one cell

   integer my_ncell
   ! Number of cell treated by current proc

   integer,allocatable :: my_cells(:)
   ! my_cells(my_ncell)
   ! Number of each cells in the supercell treat by this CPU

   integer,allocatable :: my_index_cells(:,:)
   ! my_cells(my_ncell,3)
   ! indexes of the cells in the supercell treat by this CPU

 end type effective_potential_type
!!***

interface operator (==)
  module procedure effective_potential_compare
end interface

CONTAINS  !===========================================================================================

!!****f* m_phonon_effective_potential/effective_potential_init
!!
!! NAME
!! effective_potential_init
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, qpoint chosen, and atoms
!!
!! INPUTS
!! crytal = structure with all the information for the crystal
!! epsilon_inf(3,3) = dielectric tensor
!! elastic_constants(6,6) = elastic constants tensor
!! energy = energy of the reference structure
!! dynmat(2,3,natom,3,natom,nqpt) = dynamical matrix for each qpoints
!! has_strainCoupling = true if the anharmonic par is present
!! ifcs = ifc type with cell,ewald short and total range of the ifcs
!! internal_strain(6,natom,3) = internal strain tensor
!! nqpt = number of qpoints
!! zeff(3,natom) = effective charges
!! comm  = mpi comunicator
!! elastic3rd(6,6,6) = 3 order derivatives with respect to to 3 strain
!! elastic_displacement(6,6,3,natom) =3 order derivatives with respect to 2 strain and 1 Atom disp
!! external_stress(6) = optional,extrenal strain
!! forces(3,natom) = optional,initial forces in the structure
!! internal_stress(6) = optional,initial strain in the structure
!! phonon_strain(6) = optional,ifc type for the phonon-strain coupling (should be in anharmonics_terms)
!! strain = optional, strain type
!! supercell = optional, supercell type to define
!! name = optional, name of the structure
!!
!! OUTPUT
!! eff_pot = effective_potential structure to be initialized
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_init(crystal,eff_pot,energy,ifcs,ncoeff,nqpt,comm,&
&                                   coeffs,dynmat,elastic_constants,elastic3rd,&
&                                   elastic_displacement,epsilon_inf,external_stress,&
&                                   forces,internal_strain,internal_stress,name,phonon_strain,&
&                                   phfrq,qpoints,strain,has_strainCoupling,supercell,zeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpt
 integer, intent(in) :: comm
 real(dp),intent(in):: energy
 character(len=fnlen), optional,intent(in) :: name
 integer,intent(in) :: ncoeff
 logical,intent(in) :: has_strainCoupling
!arrays
 real(dp),optional,intent(in) :: epsilon_inf(3,3),elastic_constants(6,6)
 real(dp),optional,intent(in) :: dynmat(:,:,:,:,:,:),qpoints(:,:),phfrq(:,:)
 real(dp),optional,intent(in) :: internal_strain(:,:,:),zeff(:,:,:)
 type(crystal_t),intent(in) :: crystal
 type(effective_potential_type), intent(out) :: eff_pot
 type(ifc_type),intent(in) :: ifcs
 type(supercell_type),optional,intent(in) :: supercell
 type(strain_type),optional,intent(in) :: strain
 type(ifc_type),optional,intent(in) :: phonon_strain(6)
 real(dp),optional,intent(in) :: elastic3rd(6,6,6)
 real(dp),optional,intent(in) :: elastic_displacement(6,6,3,crystal%natom)
 type(polynomial_coeff_type),optional :: coeffs(ncoeff)
 real(dp),optional,intent(in) :: external_stress(6),internal_stress(6)
 real(dp),optional,intent(in) :: forces(3,crystal%natom)
!Local variables-------------------------------
!scalar
 character(len=500) :: msg
!arrays

! *************************************************************************

!1-Free the effective potential before filling it
 call effective_potential_free(eff_pot)

 if (present(name)) then
   eff_pot%name = name
 else
   eff_pot%name = ''
 end if

!2-Perform some checks
 if (crystal%natom < 1) then
   write(msg, '(a,a,a,i10,a)' )&
&   'The cell must have at least one atom.',ch10,&
&   'The number of atom is  ',crystal%natom,'.'
   MSG_BUG(msg)
 end if

 if (crystal%ntypat < 1) then
   write(msg, '(a,a,a,i10,a)' )&
&   'The cell must have at least one type of atom.',ch10,&
&   'The number of type of atom is  ',crystal%ntypat,'.'
   MSG_BUG(msg)
 end if

!3-Fill energy of the crystal (hartree)
 eff_pot%energy = energy

!1-Fill the crystal
 eff_pot%crystal%natom  = crystal%natom
 eff_pot%crystal%ntypat = crystal%ntypat
 eff_pot%crystal%ucvol  = crystal%ucvol
 eff_pot%crystal%rprimd = crystal%rprimd


 ABI_ALLOCATE(eff_pot%crystal%xcart,(3,eff_pot%crystal%natom))
 eff_pot%crystal%xcart = crystal%xcart

 ABI_ALLOCATE(eff_pot%crystal%amu,(crystal%ntypat))
 eff_pot%crystal%amu = crystal%amu

 ABI_ALLOCATE(eff_pot%crystal%typat,(crystal%natom))
 eff_pot%crystal%typat = crystal%typat

 ABI_ALLOCATE(eff_pot%crystal%znucl,(crystal%ntypat))
 eff_pot%crystal%znucl = crystal%znucl

!4-Fill harmonic part
 call harmonics_terms_init(eff_pot%harmonics_terms,ifcs,crystal%natom,ifcs%nrpt)

!5-Init the anharmonics_terms to set the flag to false
 call anharmonics_terms_init(eff_pot%anharmonics_terms,crystal%natom,ncoeff)

!5-Fill optional inputs
 ABI_ALLOCATE(eff_pot%forces,(3,eff_pot%crystal%natom))
 eff_pot%forces = zero
 if(present(forces))then
   eff_pot%forces = forces
 end if

 if(present(elastic_constants))then
   eff_pot%harmonics_terms%elastic_constants(:,:) = elastic_constants
 end if

 if(present(epsilon_inf))then
   eff_pot%harmonics_terms%epsilon_inf(:,:) = epsilon_inf(:,:)
 end if

 if(present(dynmat).and.present(qpoints).and.present(phfrq))then
   call harmonics_terms_setDynmat(dynmat,eff_pot%harmonics_terms,crystal%natom,nqpt,phfrq,qpoints)
 end if

 if(present(internal_strain))then
   call harmonics_terms_setInternalStrain(internal_strain,eff_pot%harmonics_terms,crystal%natom)
 end if

 if(present(zeff))then
   call harmonics_terms_setEffectiveCharges(eff_pot%harmonics_terms,crystal%natom,zeff)
 end if

 eff_pot%has_strain = .FALSE.
 if(present(strain)) then
   eff_pot%has_strain = .TRUE.
   eff_pot%strain%delta = strain%delta
   eff_pot%strain%direction = strain%direction
   eff_pot%strain%strain = strain%strain
  end if

 eff_pot%internal_stress = zero
 if(present(internal_stress))then
   eff_pot%internal_stress = internal_stress
 end if

 eff_pot%external_stress = zero
 if(present(external_stress))then
   eff_pot%external_stress = external_stress
 end if

!Set the flag for the strain coupling
 eff_pot%has_strainCoupling = has_strainCoupling

!Allocation of phonon strain coupling array (3rd order)
 if(present(phonon_strain).and.has_strainCoupling) then
   call anharmonics_terms_setStrainPhononCoupling(eff_pot%anharmonics_terms,crystal%natom,phonon_strain)
 end if

!Set the 3rd order elastic tensor
 if(present(elastic3rd).and.has_strainCoupling)then
   call anharmonics_terms_setElastic3rd(eff_pot%anharmonics_terms,elastic3rd)
 end if

!Allocation of 3rd order with respecto to 2 strain and 1 atomic displacement
 if(present(elastic_displacement).and.has_strainCoupling)then
   call anharmonics_terms_setElasticDispCoupling(eff_pot%anharmonics_terms,crystal%natom,&
&                                                elastic_displacement)
 end if

!Allocation of the coefficients
 if(present(coeffs))then
   if(ncoeff /= size(coeffs))then
     write(msg, '(a)' )&
&        ' ncoeff has not the same size than coeffs array, '
     MSG_BUG(msg)
   end if
   call effective_potential_setCoeffs(coeffs,eff_pot,ncoeff)
 end if

 if(present(supercell))then
   call copy_supercell(supercell,eff_pot%supercell)
 else
   call init_supercell(eff_pot%crystal%natom, 0, real((/1,1,1/),dp), eff_pot%crystal%rprimd,&
&                      eff_pot%crystal%typat,eff_pot%crystal%xcart, eff_pot%supercell)
 end if

 call effective_potential_initmpi_supercell(eff_pot,comm)

end subroutine effective_potential_init
!!***


!!****f* m_phonon_effective_potential/effective_potential_initmpi_supercell
!! NAME
!!  effective_potential_initmpi_supercell
!!
!! FUNCTION
!!  Initializes the mpi informations for parallelism over supercell.
!!
!! INPUTS
!!   eff_pot= effective_potential_type
!!   comm = communicator
!!
!! OUTPUT
!! This is for the parallelisation over the supercell
!!
!!  eff_pot%me_supercell =  Index of my processor in the comm. over one cell
!!  eff_pot%my_ncell     =  Number of cell treated by current proc
!!  eff_pot%my_cells(:)  = Number of the cells in the supercell treat by this CPU
!!  eff_pot%my_index_cells(:,:) = indexes of the cells in the supercell treat by this CPU
!!
!! PARENTS
!!      m_effective_potential,mover_effpot
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine effective_potential_initmpi_supercell(eff_pot,comm)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_initmpi_supercell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(effective_potential_type),intent(inout)  :: eff_pot
 integer,intent(in) :: comm
!arrays

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,icell,ii
 integer :: my_rank,ncell_alone,ncell,nproc
 integer :: master = zero
 logical :: iam_master = .false.
 character(len=500) :: msg
!array
 integer :: cell_number(3)
! ***********************************************************************

!Set the number of cell in the supercell
 cell_number(:) = int(eff_pot%supercell%qphon(:))
 ncell = product(cell_number(:))

!Do some checks
 if (any(cell_number <= 0).or.ncell<=0) then
   write(msg,'(a,a)')' No supercell found for setting'
   MSG_ERROR(msg)
 end if

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

 call effective_potential_freempi_supercell(eff_pot)

 eff_pot%comm_supercell = comm
 eff_pot%me_supercell   = my_rank



!Determine the number of cell for each CPU
 ncell_alone = mod(ncell,nproc)
 eff_pot%my_ncell = aint(real(ncell,sp)/nproc)
 if(my_rank >= (nproc-ncell_alone)) then
   eff_pot%my_ncell = eff_pot%my_ncell  + 1
 end if

!Allocation of array
 ABI_ALLOCATE(eff_pot%my_cells,(eff_pot%my_ncell))
 ABI_ALLOCATE(eff_pot%my_index_cells,(eff_pot%my_ncell,3))
 eff_pot%my_cells = zero
 eff_pot%my_index_cells = zero

 do icell = 1,eff_pot%my_ncell
   if(my_rank >= (nproc-ncell_alone))then
     eff_pot%my_cells(icell)=(aint(real(ncell,sp)/nproc))*(my_rank)+&
&                              (my_rank - (nproc-ncell_alone)) + icell
   else
     eff_pot%my_cells(icell)=(eff_pot%my_ncell)*(my_rank)  + icell
   end if
 end do

 icell = 0
 ii = 0
 do i1 = 1,cell_number(1)
   do i2 = 1,cell_number(2)
     do i3 = 1,cell_number(3)
       ii = ii +1
       if(any(eff_pot%my_cells==ii))then
         icell=icell+1
         eff_pot%my_index_cells(icell,1) = i1;
         eff_pot%my_index_cells(icell,2) = i2;
         eff_pot%my_index_cells(icell,3) = i3;
       end if
     end do
   end do
 end do

end subroutine effective_potential_initmpi_supercell
!!***


!****f* m_effective_potential/effective_potential_free
!!
!! NAME
!! effective_potential_free
!!
!! FUNCTION
!! deallocate all dynamic memory for this effective potential structure
!!
!! INPUTS
!! eff_pot = supercell structure with data to be output
!!
!! OUTPUT
!! eff_pot = supercell structure with data to be output
!!
!! PARENTS
!!      compute_anharmonics,m_effective_potential,m_effective_potential_file
!!      multibinit
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_free(eff_pot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_free'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  type(effective_potential_type), intent(inout) :: eff_pot

!Local variables-------------------------------
!scalars
!array

! *************************************************************************

  eff_pot%name   = ''
  eff_pot%energy = zero
  eff_pot%internal_stress = zero
  eff_pot%external_stress = zero
  eff_pot%has_strainCoupling = .false.
  eff_pot%has_strain = .false.

   if(allocated(eff_pot%forces)) then
     eff_pot%forces=zero
     ABI_DEALLOCATE(eff_pot%forces)
   end if


! Free others structures
   call anharmonics_terms_free(eff_pot%anharmonics_terms)
   call strain_free(eff_pot%strain)
   call harmonics_terms_free(eff_pot%harmonics_terms)
   call destroy_supercell(eff_pot%supercell)
   call crystal_free(eff_pot%crystal)
   call effective_potential_freempi_supercell(eff_pot)

end subroutine effective_potential_free
!!***

!****f* m_effective_potential/effective_potential_freempi_supercell
!!
!! NAME
!! effective_potential_freempi_supercell
!!
!! FUNCTION
!! deallocate all dynamic memory for mpi of supercell
!!
!! INPUTS
!! eff_pot = supercell structure with data to be output
!!
!! OUTPUT
!! eff_pot = supercell structure with data to be output
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_freempi_supercell(eff_pot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_freempi_supercell'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
 type(effective_potential_type), intent(inout) :: eff_pot

!Local variables-------------------------------
!scalars
!array

! *************************************************************************

 eff_pot%me_supercell = 0
 eff_pot%my_ncell     = 0

 if (allocated(eff_pot%my_cells)) then
   eff_pot%my_cells(:) = zero
   ABI_DEALLOCATE(eff_pot%my_cells)
 end if

 if (allocated(eff_pot%my_index_cells)) then
   eff_pot%my_index_cells(:,:) = zero
   ABI_DEALLOCATE(eff_pot%my_index_cells)
 end if

end subroutine effective_potential_freempi_supercell
!!***


!****f* m_effective_potential/effective_potential_generateDipDip
!!
!! NAME
!! effective_potential_generateDipDip
!!
!! FUNCTION
!! Generate the supercell of the structure inside of the effective
!! potential and fill the supercell structure. Also adapt the harmonic
!! part for the supercell (compute dipole-dipole interation)
!!
!! INPUTS
!! eff_pot = effective potential structure
!! option  =  0 Just generate supercell and fill effective potential
!!            1 generate supercell and harmonic part if need
!!            2 generate supercell and force to generate new harmonic part
!! n_cell  = (3) number of cell in the direction x, y and z
!! comm=MPI communicator
!!
!! OUTPUT
!! eff_pot
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_generateDipDip(eff_pot,n_cell,option,asr,comm)

 use m_phonon_supercell
 use m_ewald

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_generateDipDip'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,asr
 integer,intent(in) :: comm
!array
 integer,intent(in) :: n_cell(3)
 type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
 integer,parameter :: master=0
 integer :: first_coordinate
 integer :: ia,i1,i2,i3,ii,ierr,irpt,irpt2,irpt_ref,min1,min2,min3
 integer :: min1_cell,min2_cell,min3_cell,max1_cell,max2_cell,max3_cell
 integer :: max1,max2,max3,my_rank,ncell,nproc,second_coordinate,size,sumg0
 integer :: my_nrpt,nrpt_alone
 real(dp) :: ucvol
 character(len=500) :: message
 logical :: iam_master=.FALSE.
!array
 integer,allocatable :: my_index_rpt(:,:)
 integer,allocatable :: bufsize(:),bufdisp(:)
 integer,allocatable :: my_irpt(:)
 real(dp) :: acell(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: buff_ewald(:,:,:,:,:,:),dyew(:,:,:,:,:), dyewq0(:,:,:,:,:)
 real(dp),allocatable :: xred(:,:),xred_tmp(:,:),zeff_tmp(:,:,:)
 type(supercell_type) :: supercell
 type(ifc_type) :: ifc_tmp

! *************************************************************************

!0 MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)
 ierr=0

!0 Check the size of the cell
 do ia=1,3
   if(n_cell(ia)<0.or.n_cell(ia)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'n_cell(',ia,') is ',n_cell(ia),', which is lower than 0 of superior than 20.',&
&     ch10,'Action: correct n_cell(',ia,').'
     MSG_ERROR(message)
   end if
 end do

 call init_supercell(eff_pot%crystal%natom, 0,&
&  real(n_cell,dp),&
&   eff_pot%crystal%rprimd,&
&   eff_pot%crystal%typat,&
&   eff_pot%crystal%xcart,&
&   supercell)

 ncell = product(n_cell(:))

!1 Store the information of the supercell of the reference structure into effective potential
 call copy_supercell(supercell,eff_pot%supercell)
!2 Initialisation of new mpi over supercell
 call effective_potential_initmpi_supercell(eff_pot,comm)

!3 Check if the bound of new cell correspond to the effective potential
!only for option=zero
!set min and max
 min1 = minval(eff_pot%harmonics_terms%ifcs%cell(1,:))
 min2 = minval(eff_pot%harmonics_terms%ifcs%cell(2,:))
 min3 = minval(eff_pot%harmonics_terms%ifcs%cell(3,:))
 max1 = maxval(eff_pot%harmonics_terms%ifcs%cell(1,:))
 max2 = maxval(eff_pot%harmonics_terms%ifcs%cell(2,:))
 max3 = maxval(eff_pot%harmonics_terms%ifcs%cell(3,:))

 if(option==0) then
   if(((max1-min1+1)/=n_cell(1).and.&
&    (max2-min2+1)/=n_cell(2).and.(max3-min3+1)/=n_cell(3))) then
     write(message, '(88a,3I3,5a,3I3,a)' )ch10,('-',i1=1,80),ch10,ch10,&
&      ' WARNING: dipdip is set to zero, the longe range interation might be wrong',ch10,&
&      '          because it is not recompute.',ch10,&
&      '          The previous harmonic part is build for ',(max1-min1+1),(max2-min2+1),(max3-min3+1)&
&,     ' cell.',ch10,'          Be sure than the dipole-dipole interation ',ch10,&
&      '          is correct for the supercell: ',n_cell(:),' or set dipdip to 1'
     call wrtout(std_out,message,"COLL")
   else
     write(message,'(84a)')ch10,('-',i1=1,80),ch10,ch10,&
&    ' WARNING: dipdip is set to zero, the longe range interation is not recompute'
     call wrtout(std_out,message,"COLL")
   end if

   write(message,'(a,(80a))') ch10,('=',i1=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!4-Adapt harmonic part
 else if (option>=1.and.all(n_cell(:)>0)) then

   write(message,'(a,(80a),3a)') ch10,('=',i1=1,80),ch10,' Generation of new ifc',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   irpt_ref = 0
   irpt     = 0
   min1_cell = zero; max1_cell = zero
   min2_cell = zero; max2_cell = zero
   min3_cell = zero; max3_cell = zero

   call find_bound(min1_cell,max1_cell,n_cell(1))
   call find_bound(min2_cell,max2_cell,n_cell(2))
   call find_bound(min3_cell,max3_cell,n_cell(3))

   write(message, '(2a)' )&
&        ' dipdip is set to one, the dipole-dipole interation is recompute.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  Generate new bound
   if(option==1)then
     if ((abs(min1) > abs(min1_cell)).or.(abs(max1) > abs(max1_cell)).or.&
&        (abs(min2) > abs(min2_cell)).or.(abs(max2) > abs(max2_cell)).or.&
&        (abs(min3) > abs(min3_cell)).or.(abs(max3) > abs(max3_cell))) then
       write(message, '(6a,3I4,5a)' )ch10,&
&        ' --- !WARNING',ch10,&
&        '     The previous harmonic part was build for bigger cell,',ch10,&
&        '     ifc is adjust on ',int((/(max1-min1+1),(max2-min2+1),(max3-min3+1)/),dp),' cell',ch10,&
&        '     The simulation might be completely wong',ch10,&
&        ' ---'
       call wrtout(std_out,message,"COLL")
       if (abs(min1) < abs(min1_cell)) min1 = min1_cell
       if (abs(min2) < abs(min2_cell)) min2 = min2_cell
       if (abs(min3) < abs(min3_cell)) min3 = min3_cell
       if (abs(max1) < abs(max1_cell)) max1 = max1_cell
       if (abs(max2) < abs(max2_cell)) max2 = max2_cell
       if (abs(max3) < abs(max3_cell)) max3 = max3_cell
!      If the cell is smaller, we redifine new cell to take into acount all atoms
       call destroy_supercell(supercell)
       call init_supercell(eff_pot%crystal%natom, 0,real((/(max1-min1+1),&
&                          (max2-min2+1),(max3-min3+1)/),dp),&
&                          eff_pot%crystal%rprimd,eff_pot%crystal%typat,&
&                          eff_pot%crystal%xcart,supercell)

!      Store the information of the supercell of the reference structure into effective potential
       call copy_supercell(supercell,eff_pot%supercell)
!      Initialisation of new mpi over supercell
       call effective_potential_initmpi_supercell(eff_pot,comm)
     else
       min1 = min1_cell ; min2 = min2_cell ; min3 = min3_cell
       max1 = max1_cell ; max2 = max2_cell ; max3 = max3_cell
     end if
   end if

!  Print the new boundary
   write(message,'(5a,2I3,a,2I3,a,2I3,4a)') ch10,' New bound for ifc (long+short):',&
&    ch10,ch10, " x=[",min1,max1,"], y=[",min2,max2,"] and z=[",min3,max3,"]",ch10,ch10,&
&    " Computation of new dipole-dipole interaction."
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  Count the new number of ifc
   do i1=min1,max1
     do i2=min2,max2
       do i3=min3,max3
         irpt = irpt +1
         if(i1==0.and.i2==0.and.i3==0) irpt_ref = irpt
       end do
     end do
   end do

   ifc_tmp%nrpt = irpt
   ABI_ALLOCATE(ifc_tmp%cell,(3,ifc_tmp%nrpt))
   ifc_tmp%cell(:,:) = zero

!  Set MPI here and not at the begining because the number of cell is adjust just before
!  Here we store in my_irpt a list with the number of each cell to be treat by this CPU
!  Determine the number of cell for each CPU
   ABI_ALLOCATE(bufsize,(nproc))
   ABI_ALLOCATE(bufdisp,(nproc))

   nrpt_alone = mod(ifc_tmp%nrpt,nproc)
   my_nrpt = aint(real(ifc_tmp%nrpt,sp)/nproc)
   if(my_rank >= (nproc-nrpt_alone)) then
     my_nrpt = my_nrpt  + 1
   end if

!  Initialisation of ifc temporary
   ABI_ALLOCATE(buff_ewald,(2,3,eff_pot%crystal%natom,3,eff_pot%crystal%natom,my_nrpt))
   ABI_ALLOCATE(ifc_tmp%short_atmfrc,(2,3,eff_pot%crystal%natom,3,eff_pot%crystal%natom,ifc_tmp%nrpt))
   ABI_ALLOCATE(ifc_tmp%ewald_atmfrc,(2,3,eff_pot%crystal%natom,3,eff_pot%crystal%natom,ifc_tmp%nrpt))
   ABI_ALLOCATE(ifc_tmp%atmfrc,(2,3,eff_pot%crystal%natom,3,eff_pot%crystal%natom,ifc_tmp%nrpt))
   ABI_ALLOCATE(my_irpt,(my_nrpt))
   ABI_ALLOCATE(my_index_rpt,(3,my_nrpt))

   my_irpt = zero
   my_index_rpt(:,:) = zero
   ifc_tmp%atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%short_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%ewald_atmfrc(:,:,:,:,:,:) = zero
   buff_ewald(:,:,:,:,:,:) = zero

!  Allocation of array
   do irpt = 1,my_nrpt
     if(my_rank >= (nproc-nrpt_alone))then
       my_irpt(irpt)=(aint(real(ifc_tmp%nrpt,sp)/nproc))*(my_rank)+&
&                       (my_rank - (nproc-nrpt_alone)) + irpt
     else
       my_irpt(irpt)=(my_nrpt)*(my_rank) + irpt
     end if
   end do

   irpt = 0
   ii = 0
   do i1=min1,max1
     do i2=min2,max2
       do i3=min3,max3
         ii = ii +1
         ifc_tmp%cell(1,ii) = i1; ifc_tmp%cell(2,ii) = i2; ifc_tmp%cell(3,ii) = i3;
         if(any(my_irpt==ii))then
           irpt=irpt+1
           my_index_rpt(1,irpt) = i1;
           my_index_rpt(2,irpt) = i2;
           my_index_rpt(3,irpt) = i3;
         end if
       end do
     end do
   end do

!  Allocate and initialize some array
   ABI_ALLOCATE(xred_tmp,(3,2*eff_pot%crystal%natom))
   ABI_ALLOCATE(xred,(3,supercell%natom_supercell))
   ABI_ALLOCATE(zeff_tmp,(3,3,2*eff_pot%crystal%natom))
   ABI_ALLOCATE(dyew,(2,3,2*eff_pot%crystal%natom,3,2*eff_pot%crystal%natom))
   ABI_ALLOCATE(dyewq0,(2,3,eff_pot%crystal%natom,3,eff_pot%crystal%natom))

   dyew            = zero
   dyewq0          = zero
   xred(:,:)       = zero
   xred_tmp(:,:)   = zero
   zeff_tmp(:,:,:) = zero
   sumg0           = zero
   acell           = one

   call matr3inv(supercell%rprimd_supercell,gprimd)
   call xcart2xred(supercell%natom_supercell,supercell%rprimd_supercell,&
&                  supercell%xcart_supercell,xred)
   call metric(gmet,gprimd,-1,rmet,supercell%rprimd_supercell,ucvol)

!  Fill the atom position of the first cell (reference cell)
   first_coordinate  = ((irpt_ref-1)*eff_pot%crystal%natom) + 1
   second_coordinate = first_coordinate + eff_pot%crystal%natom-1
   xred_tmp(:,1:eff_pot%crystal%natom) = xred(:,first_coordinate:second_coordinate)
!  Fill fake zeff array for ewald9
   zeff_tmp(:,:,1:eff_pot%crystal%natom) = eff_pot%harmonics_terms%zeff
   zeff_tmp(:,:,eff_pot%crystal%natom+1:2*eff_pot%crystal%natom) = eff_pot%harmonics_terms%zeff

   do irpt=1,my_nrpt
     i1=my_index_rpt(1,irpt); i2=my_index_rpt(2,irpt); i3=my_index_rpt(3,irpt)
!    Compute new dipole-dipole interaction
     dyew = zero
     if (i1==0.and.i2==0.and.i3==0) then
       call ewald9(acell,eff_pot%harmonics_terms%epsilon_inf,dyewq0,&
&                  gmet,gprimd,eff_pot%crystal%natom,real((/0,0,0/),dp),rmet,&
&                  supercell%rprimd_supercell,sumg0,ucvol,xred_tmp(:,1:eff_pot%crystal%natom),&
&                  eff_pot%harmonics_terms%zeff)
       buff_ewald(:,:,:,:,:,irpt) = dyewq0 + tol10
     else
       first_coordinate  = ((my_irpt(irpt)-1)*eff_pot%crystal%natom) + 1
       second_coordinate = first_coordinate + eff_pot%crystal%natom  - 1
       xred_tmp(:,eff_pot%crystal%natom+1:2*eff_pot%crystal%natom)=&
&              xred(:,first_coordinate:second_coordinate)
       call ewald9(acell,eff_pot%harmonics_terms%epsilon_inf,dyew,gmet,gprimd,&
&                  int(2*eff_pot%crystal%natom),real((/0,0,0/),dp),&
&                  rmet,supercell%rprimd_supercell,&
&                  sumg0,ucvol,xred_tmp,zeff_tmp)
       buff_ewald(:,:,:,:,:,irpt) = &
&           dyew(:,:,1:eff_pot%crystal%natom,:,eff_pot%crystal%natom+1:2*eff_pot%crystal%natom) + tol10
     end if
   end do

!  DEALLOCATION OF ARRAYS
   ABI_DEALLOCATE(my_index_rpt)
   ABI_DEALLOCATE(my_irpt)
   ABI_DEALLOCATE(xred_tmp)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(zeff_tmp)
   ABI_DEALLOCATE(dyew)
   ABI_DEALLOCATE(dyewq0)

!  Set the bufsize for mpi allgather
   do ii = 1,nproc
     bufsize(ii) = aint(real(ifc_tmp%nrpt,sp)/nproc)*2*3*eff_pot%crystal%natom*3*eff_pot%crystal%natom
     if(ii > (nproc-nrpt_alone)) then
       bufsize(ii) = bufsize(ii) + 2*3*eff_pot%crystal%natom*3*eff_pot%crystal%natom
     end if
   end do

   bufdisp(1) = 0
   do ii = 2,nproc
     bufdisp(ii) = bufdisp(ii-1) + bufsize(ii-1)
   end do

   size = 2*3*eff_pot%crystal%natom*3*eff_pot%crystal%natom*my_nrpt
   call xmpi_allgatherv(buff_ewald,size,ifc_tmp%ewald_atmfrc,bufsize,bufdisp, comm, ierr)

   ABI_DEALLOCATE(bufsize)
   ABI_DEALLOCATE(bufdisp)
   ABI_DEALLOCATE(buff_ewald)

!  Fill the short range part (calculated previously) only master
   if(iam_master)then
     do irpt=1,ifc_tmp%nrpt
       do irpt2=1,eff_pot%harmonics_terms%ifcs%nrpt
         if(eff_pot%harmonics_terms%ifcs%cell(1,irpt2)==ifc_tmp%cell(1,irpt).and.&
&           eff_pot%harmonics_terms%ifcs%cell(2,irpt2)==ifc_tmp%cell(2,irpt).and.&
&           eff_pot%harmonics_terms%ifcs%cell(3,irpt2)==ifc_tmp%cell(3,irpt).and.&
&           any(eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:,irpt2) > tol9)) then
           ifc_tmp%short_atmfrc(:,:,:,:,:,irpt) = &
&                               eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:,irpt2)
         end if
       end do
     end do
   end if

   call xmpi_bcast (ifc_tmp%short_atmfrc, master, comm, ierr)

!  Compute total ifc
   ifc_tmp%atmfrc = ifc_tmp%short_atmfrc + ifc_tmp%ewald_atmfrc


!  Copy ifc into effective potential
!  !!!Warning eff_pot%harmonics_terms%ifcs only contains atmfrc,short_atmfrc,ewald_atmfrc,nrpt
!    and cell!!  rcan,ifc%rpt,wghatm and other quantities
!    are not needed for effective potential!!!
!  Free ifc before copy
   call ifc_free(eff_pot%harmonics_terms%ifcs)
!  Fill the effective potential with new atmfr
   eff_pot%harmonics_terms%ifcs%nrpt = ifc_tmp%nrpt
   call alloc_copy(ifc_tmp%atmfrc      ,eff_pot%harmonics_terms%ifcs%atmfrc)
   call alloc_copy(ifc_tmp%short_atmfrc,eff_pot%harmonics_terms%ifcs%short_atmfrc)
   call alloc_copy(ifc_tmp%ewald_atmfrc,eff_pot%harmonics_terms%ifcs%ewald_atmfrc)
   call alloc_copy(ifc_tmp%cell        ,eff_pot%harmonics_terms%ifcs%cell)

!  Free temporary ifc
   call ifc_free(ifc_tmp)
 end if

 if(asr >= 0) then
! Impose sum rule
   call effective_potential_applySumRule(asr,eff_pot%harmonics_terms%ifcs,eff_pot%crystal%natom)
 end if

 write(message, '(a,(80a),a)' ) ch10,&
&   ('=',ii=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

! Free suppercell
 call destroy_supercell(supercell)

end subroutine effective_potential_generateDipDip
!!***


!****f* m_effective_potential/effective_potential_applySumRule
!!
!! NAME
!! effective_potential_applySumRule
!!
!! FUNCTION
!! Apply the acoustic sum rule on the effective potential
!!
!! INPUTS
!! ifc = inter-atominc forces constants
!! asr     = acoustic sum rule option (see anaddb help)
!! option  = optional if |no present asr is done on total ifc
!!                       |present and 1 asr is done on short part
!!                       |present and 2 asr is done on ewald part
!!
!! OUTPUT
!! eff_pot
!!
!! PARENTS
!!      compute_anharmonics,m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_applySumRule(asr,ifc,natom,option)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_applySumRule'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: asr
 integer,intent(in) :: natom
 integer,optional,intent(in) :: option
!array
 type(ifc_type),target,intent(inout) :: ifc
!Local variables-------------------------------
!scalar
 integer :: ia,ib,irpt,irpt_ref
 integer :: mu,nu
 real(dp) :: sum
 character(500) :: message
!array
 real(dp),pointer :: atmfrc(:,:,:,:,:,:)
! *************************************************************************

 irpt_ref = 0
! Found the cell of reference
 do irpt = 1,ifc%nrpt
   if(ifc%cell(1,irpt)==0.and.&
&     ifc%cell(2,irpt)==0.and.&
&     ifc%cell(3,irpt)==0) then
     irpt_ref = irpt
     cycle
   end if
 end do

 if (irpt_ref<=0) then
   write(message,'(a,a)')' Unable to find the cell of reference in IFC'
   MSG_ERROR(message)
 end if

 if (present(option)) then
   if (option == 1) then
     atmfrc => ifc%short_atmfrc
     write(message,'(3a)') ch10," Impose acoustic sum rule on short range"
   else if (option == 2) then
     atmfrc => ifc%ewald_atmfrc
     write(message,'(3a)') ch10," Impose acoustic sum rule on long range"
   end if
 else
   atmfrc => ifc%atmfrc
     write(message,'(3a)') ch10," Impose acoustic sum rule on total ifc"
 end if
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!impose acoustic sum rule:
 do mu=1,3
   do nu=1,3
     do ia=1,natom
       sum=zero
       do ib=1,natom
!      Get the sum of interatomic forces acting on the atom ia,
!      either in a symmetrical manner, or an unsymmetrical one.
         if(asr==1)then
           do irpt=1, ifc%nrpt
             sum=sum+atmfrc(1,mu,ia,nu,ib,irpt)
           end do
         else if(asr==2)then
           do irpt=1, ifc%nrpt
              sum=sum+&
&                 (atmfrc(1,mu,ia,nu,ib,irpt)+&
&                  atmfrc(1,nu,ia,mu,ib,irpt))/2
            end do
          end if
        end do
!      Correct the self-interaction in order to fulfill the ASR
        atmfrc(1,mu,ia,nu,ia,irpt_ref)=&
&       atmfrc(1,mu,ia,nu,ia,irpt_ref)-sum
        if(asr==2)then
          atmfrc(1,nu,ia,mu,ia,irpt_ref)=&
&         atmfrc(1,mu,ia,nu,ia,irpt_ref)
        end if
      end do
    end do
  end do

 if (present(option)) then
   if (option == 1) then
     ifc%short_atmfrc = atmfrc
   else if (option == 2) then
     ifc%ewald_atmfrc = atmfrc
   end if
 else
   ifc%atmfrc = atmfrc
 end if

 end subroutine effective_potential_applySumRule
!!***


!****f* m_effective_potential/effective_potential_effpot2dynmat
!!
!! NAME
!! effective_potential_effpot2dynmat
!!
!! FUNCTION
!! Compute finate differences on forces to compute dynmical matrix
!! at gamma for supercell
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!! dynmat   = ddb with all information
!!
!! PARENTS
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_effpot2dynmat(dynmat,delta,eff_pot,natom,n_cell,option)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_effpot2dynmat'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option
 real(dp),intent(in) :: delta
!array
 integer,intent(in) :: n_cell(3)
 type(effective_potential_type),intent(in) :: eff_pot
 real(dp),intent(out) :: dynmat(2,3,natom,3,natom)
!Local variables-------------------------------
!scalar
 character(len=500) :: msg
 integer :: ia,ib,ii,mu,nu,npt
 real(dp):: delt
!array
 real(dp),allocatable :: disp(:,:),diff(:,:,:)
 real(dp),allocatable :: fred(:,:),fcart(:,:)

! *************************************************************************

!Do Some check
 if (natom /=  eff_pot%supercell%natom_supercell) then
   write(msg, '(a,i10,2a,i10,a)' )&
&   'The number of atom ',natom ,ch10,&
&   'is not the same than in supercell  ', eff_pot%supercell%natom_supercell,'.'
   MSG_ERROR(msg)
 end if

 do ii=1,3
   if(n_cell(ii)<0.or.n_cell(ii)>20)then
     write(msg, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'n_cell(',ii,') is ',n_cell(ii),', which is lower than 0 of superior than 20.',&
&     ch10,'Action: correct n_cell(',ii,').'
     MSG_ERROR(msg)
   end if
 end do

 do ii=1,3
   if(eff_pot%supercell%qphon(ii) /= n_cell(ii))then
     write(msg, '(a,i0,a,i0,a,i0,a)' )&
&     'n_cell(',ii,') is ',n_cell(ii),', which not conrespond to n_cell in supercell ',&
&     eff_pot%supercell%qphon(ii)
     MSG_ERROR(msg)
   end if
 end do

 write(msg,'(a,(80a),3a)') ch10,('-',ii=1,80),ch10,' Generation of the dynmical matrix by ',&
& 'finite differences'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 select case (option)
 case (1)
   write(msg,'(2a)') ch10,' Finite differences on 1 points  '
   npt = 1
 case (2)
   write(msg,'(2a)') ch10,' Finite differences on 3 points  '
   npt = 3
 case (3)
   write(msg,'(2a)') ch10,' Finite differences on 5 points  '
   npt = 5
 end select

 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

! Allocation of forces arrays
 ABI_ALLOCATE(disp,(3,natom))
 ABI_ALLOCATE(diff,(npt,3,natom))
 ABI_ALLOCATE(fred,(3,natom))
 ABI_ALLOCATE(fcart,(3,natom))

 dynmat = zero

 do ia=1,eff_pot%supercell%natom_supercell
   do mu=1,3
     diff = zero
     do ii=1,npt
       delt = (-(npt/2+1)+ii) * delta
       disp = zero
       disp(mu,ia) = delt * eff_pot%supercell%rprimd_supercell(mu,mu)
       call effective_potential_getForces(eff_pot,fcart,fred,&
&                                         eff_pot%supercell%natom_supercell,&
&                                         eff_pot%supercell%rprimd_supercell,&
&                                         eff_pot%supercell%xcart_supercell,displacement=disp)
       diff(ii,:,:) = fred(:,:)
     end do
     select case (option)
     case (1)
       dynmat(1,mu,ia,:,:) = (diff(2,:,:)) / (delta)
     case (2)
       dynmat(1,mu,ia,:,:) = (diff(3,:,:)-diff(1,:,:)) / (2*delta)
     case (3)
       dynmat(1,mu,ia,:,:) = (-diff(5,:,:)+8*diff(4,:,:)-8*diff(2,:,:)+diff(1,:,:)) / (12*delta)
     end select
   end do
 end do

!Write the phonon  into ddb format wavevector
 write(999, '(a,3es16.8,f6.1)' )' qpt',real((/0,0,0/),dp),1.0

!Write the matrix elements
 do ib=1,eff_pot%supercell%natom_supercell
   do nu=1,3
     do ia=1,eff_pot%supercell%natom_supercell
       do mu=1,3
          write(999,'(4i4,2d22.14)')nu,ib,mu,ia,&
 &             dynmat(1,mu,ia,nu,ib),dynmat(2,mu,ia,nu,ib)
       end do
     end do
   end do
 end do

! Deallocation of arrays
 ABI_DEALLOCATE(fred)
 ABI_DEALLOCATE(fcart)


 end subroutine effective_potential_effpot2dynmat
!!***

!****f* m_effective_potential/effective_potential_setCoeffs
!!
!! NAME
!! effective_potential_setCoeffs
!!
!! FUNCTION
!! Set the coefficients of  the effective_potential in ouput
!!
!! INPUTS
!! coeffs = polynomial_coeff_type
!! eff_pot = effective potential structure
!! ncoeff = number of coefficient
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_effective_potential,m_effective_potential_file
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_setCoeffs(coeffs,eff_pot,ncoeff)

 use m_polynomial_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_setCoeffs'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: ncoeff
!array
  type(effective_potential_type),intent(inout) :: eff_pot
  type(polynomial_coeff_type),intent(in) :: coeffs(ncoeff)
!Local variables-------------------------------
!scalar
  character(len=500) :: msg
!array
! *************************************************************************

  if(ncoeff /= size(coeffs))then
    write(msg, '(a)' )&
&        ' ncoeff has not the same size than coeffs array, '
    MSG_BUG(msg)
  end if

  call anharmonics_terms_setCoeffs(coeffs,eff_pot%anharmonics_terms,ncoeff)

end subroutine effective_potential_setCoeffs
!!***

!****f* m_effective_potential/effective_potential_setElastic3rd
!!
!! NAME
!! effective_potential_setElastic3rd
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
!!      compute_anharmonics
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
 
subroutine effective_potential_setElastic3rd(eff_pot,elastics)

 use m_polynomial_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_setElastic3rd'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  real(dp),intent(in) :: elastics(6,6,6)
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
!array
! *************************************************************************
  call anharmonics_terms_setElastic3rd(eff_pot%anharmonics_terms,elastics)

end subroutine effective_potential_setElastic3rd
!!***

!****f* m_effective_potential/effective_potential_setElastic4rd
!!
!! NAME
!! effective_potential_setElastic4rd
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
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
 
subroutine effective_potential_setElastic4rd(eff_pot,elastics)

 use m_polynomial_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_setElastic4rd'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  real(dp),intent(in) :: elastics(6,6,6,6)
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
!array
! *************************************************************************
  call anharmonics_terms_setElastic4rd(eff_pot%anharmonics_terms,elastics)

end subroutine effective_potential_setElastic4rd
!!***

!****f* m_effective_potential/effective_potential_setStrainPhononCoupling
!!
!! NAME
!! effective_potential_setStrainPhononCoupling
!!
!! FUNCTION
!! Set the strain phonon coupling of  the effective_potential 
!!
!! INPUTS
!! natom  = number of atoms
!! ncoeff = number of coefficient
!! nrpt   = number of rpt 
!! strain_phonon = (size 6) array of type ifc   
!!
!! OUTPUT
!! eff_pot = effective potential structure
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
 
subroutine effective_potential_setStrainPhononCoupling(eff_pot,natom,nrpt,phonon_strain)

 use m_polynomial_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_setStrainPhononCoupling'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom,nrpt
!array
  type(ifc_type),intent(in) :: phonon_strain(6)
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
!array
! *************************************************************************

  call anharmonics_terms_setStrainPhononCoupling(eff_pot%anharmonics_terms,natom,phonon_strain)
  eff_pot%has_strainCoupling = .True.

end subroutine effective_potential_setStrainPhononCoupling
!!***

!****f* m_effective_potential/effective_potential_setElasticDispCoupling
!!
!! NAME
!! effective_potential_setElasticDispCoupling
!!
!! FUNCTION
!! Set the elastic constant displacement coupling of  the effective_potential 
!!
!! INPUTS
!! natom  = number of atoms
!! elastic_displacement = (6,6,3,natom) array with elastic constant displacement coupling
!!
!! OUTPUT
!! eff_pot = effective potential structure
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
 
subroutine effective_potential_setElasticDispCoupling(eff_pot,natom,elastic_displacement)

 use m_polynomial_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_setElasticDispCoupling'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  real(dp) :: elastic_displacement(6,6,3,natom)
  type(effective_potential_type),intent(inout) :: eff_pot 
!Local variables-------------------------------
!scalar
!array
! *************************************************************************

  call anharmonics_terms_setElasticDispCoupling(eff_pot%anharmonics_terms,natom,elastic_displacement)
  eff_pot%has_strainCoupling = .True.

end subroutine effective_potential_setElasticDispCoupling
!!***


!****f* m_effective_potential/effective_potential_print
!!
!! NAME
!! effective_potential_print
!!
!! FUNCTION
!! write the effective_potential in ouput
!!
!! INPUTS
!! eff_pot = effective potential structure
!! option  =  0 no output
!! option  = -1 only generate xml file
!! option  =  1 only useful information
!! option  =  2 all informations including ifc
!! option  =  3 only useful information + xml file
!! option  =  4 only all informations including ifc + xml file
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_print(eff_pot,option,filename)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_print'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: option
  character(len=*),optional,intent(in) :: filename
!array
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
  integer :: ia,ii
  real(dp):: fact
  character(len=500) :: message
!array
! *************************************************************************

  if(option >= 1) then
    if(present(filename)) then
      write(message, '(a,a,a,a,a,a)' )ch10,' The file ',trim(filename),&
&     ' contains this effective potential for ',trim(eff_pot%name),':'
    else
      write(message, '(a,a,a,a)' )ch10,' This effective potential contains ',&
&     trim(eff_pot%name),':'
    end if

    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

!**********************************************************************
! Write basics values
!**********************************************************************

    write(message,'(a,F20.10,2a,I3,2a,I4,2a,I4,2a,I3,2a)') &
&     '  - Reference energy:  ',eff_pot%energy ,ch10,&
&     '  - Number of types of atoms:  ',eff_pot%crystal%ntypat ,ch10,&
&     '  - Number of atoms:  ',eff_pot%crystal%natom ,ch10,&
&     '  - Number of cells:  ',eff_pot%harmonics_terms%ifcs%nrpt ,ch10,&
&     '  - Number of qpoints:  ',eff_pot%harmonics_terms%nqpt ,ch10,&
&     '  - Primitive vectors (unit:Bohr):  '
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    do ii = 1,3
      write(message,'(3F12.6)') eff_pot%crystal%rprimd(1,ii),&
&                               eff_pot%crystal%rprimd(2,ii),&
&                               eff_pot%crystal%rprimd(3,ii)
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    end do
    write(message,'(2a,3F12.6)') '  - acell (unit:Bohr):',ch10,one,one,one

    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    write(message,'(a)') '  - Dielectric tensor:  '
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    do ii=1,3
      write(message,'(3F12.6)')eff_pot%harmonics_terms%epsilon_inf(1,ii),&
&                              eff_pot%harmonics_terms%epsilon_inf(2,ii),&
&                              eff_pot%harmonics_terms%epsilon_inf(3,ii)
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    end do
      write(message,'(a)') '  - Elastic tensor (unit:10^2GPa):  '
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
      fact = HaBohr3_GPa / eff_pot%crystal%ucvol
    do ii=1,6
      write(message,'(6F12.6)')&
&      eff_pot%harmonics_terms%elastic_constants(1,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(2,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(3,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(4,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(5,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(6,ii)*fact/100
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    end do
    do ia=1,eff_pot%crystal%natom
      write(message,'(a,I4,2a,F10.4,2a,F10.4,2a,3F12.6,2a)')'  - Atoms',ia,ch10,&
&             "    - atomic number:",eff_pot%crystal%znucl(eff_pot%crystal%typat(ia)),ch10,&
&             "    - atomic mass:",eff_pot%crystal%amu(eff_pot%crystal%typat(ia)),ch10,&
&             "    - cartesian position:",eff_pot%crystal%xcart(:,ia),ch10,&
&             "    - Effective charges:"
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
      do ii = 1,3
        write(message,'(a,3(F12.6))') "  ",eff_pot%harmonics_terms%zeff(:,ii,ia)
        call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,message,'COLL')
      end do
    end do
  end if

  if (option >= 3) then
!   to do
  end if
end subroutine effective_potential_print
!!***

!****f* m_effective_potential/effective_potential_printSupercell
!!
!! NAME
!! effective_potential_printSupercell
!!
!! FUNCTION
!! Print the supercell of the effective_potential,
!! or if present the supercell as input
!! WARNING: need to be consistent with eff_pot
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      mover_effpot
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_printSupercell(eff_pot,supercell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_printSupercell'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(supercell_type),optional,target,intent(in) :: supercell
!Local variables-------------------------------
!scalar
 integer :: iatom,ii
 character(len=500) :: message
!array
 real(dp), allocatable :: xred(:,:)
 type(supercell_type),pointer :: supercell_tmp

! *************************************************************************

 if(present(supercell)) then
   supercell_tmp => supercell
 else
   supercell_tmp => eff_pot%supercell
 end if

 if(supercell_tmp%natom_supercell /= eff_pot%supercell%natom_supercell) then
   write(message, '(3a)' )&
&  ' There is not the same numbers of atoms in the two supercell',ch10,&
&   'Action: modify the code'
   MSG_BUG(message)
 end if

 ABI_ALLOCATE(xred,(3,supercell_tmp%natom_supercell))

!**********************************************************************
! Write basics values
!**********************************************************************

 write (message, '(4a,I8,a)') ' Structure parameters of the supercell :',ch10,ch10,&
                       '  natom ', supercell_tmp%natom_supercell,ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write (message, '(a)') '  znucl '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write(message,*) ''
 do iatom = 1, size(eff_pot%crystal%znucl)
   write (message, '(a,I5)') trim(message),int(eff_pot%crystal%znucl(iatom))
   if (mod(iatom,6) == 0) then
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     write(message,*) ''
   end if
 end do
 write (message, '(2a)') trim(message),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write (message, '(a,I7,3a)') '  ntypat', size(eff_pot%crystal%znucl),ch10,ch10, '  typat '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,*) ''
 do iatom = 1, supercell_tmp%natom_supercell
   write (message, '(a,I5)') trim(message),&
&         supercell_tmp%typat_supercell(supercell_tmp%atom_indexing_supercell(iatom))
   if (mod(iatom,12) == 0)then
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     write(message,*) ''
   end if
 end do
 write (message, '(2a)') trim(message),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write (message, '(3a)') '  acell 1.0 1.0 1.0',ch10,ch10
 write (message, '(2a)') trim(message),'  rprim'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 do ii = 1,3
   write(message,'(3E23.14,3E23.14,3E23.14)') supercell_tmp%rprimd_supercell(1,ii),&
&                                             supercell_tmp%rprimd_supercell(2,ii),&
&                                             supercell_tmp%rprimd_supercell(3,ii)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end do

  write (message, '(2a)') ch10,'  xcart'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  do iatom = 1, supercell_tmp%natom_supercell
    write (message, '(3E23.14)') supercell_tmp%xcart_supercell(1,iatom),&
&                                supercell_tmp%xcart_supercell(2,iatom),&
&                                supercell_tmp%xcart_supercell(3,iatom)
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end do
  call xcart2xred(supercell_tmp%natom_supercell,supercell_tmp%rprimd_supercell,&
 &                supercell_tmp%xcart_supercell,xred)
  write (message, '(2a)') ch10,'  xred'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  do iatom = 1, supercell_tmp%natom_supercell
    write (message, '(3E23.14)') xred(1,iatom),&
&                                xred(2,iatom),&
&                                xred(3,iatom)
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end do

! Deallocation array
  ABI_DEALLOCATE(xred)

end subroutine effective_potential_printSupercell
!!***

!!****f* m_effective_potential/effective_potential_writeXML
!! NAME
!! effective_potential_writeXML
!!
!! FUNCTION
!! This routine print the effective potential into xml format
!! Several options are available
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filename = the name of output file
!! eff_pot  = structure contains the effective potential
!! option   = option for the format of the xml file
!!           -1 print both system.xml and coefficients from fitted polynomial 
!!            1 print the xml for a system
!!            2 print the coefficients from the fitted polynomial
!!
!!
!! OUTPUT
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_writeXML(eff_pot,option,filename)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_multibinit_dataset, only : multibinit_dataset_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_writeXML'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: option
  character(len=fnlen),optional,intent(in) :: filename
!arrays
  type(effective_potential_type), intent(in) :: eff_pot

!Local variables-------------------------------
!scalar
 integer :: ii,ia,ib,jj
 integer :: iqpt,irpt,mu,nu
 integer :: unit_xml=22
 character(len=500) :: message
 character(len=fnlen) :: namefile
 character(len=10) :: natom
!arrays
 real :: strain(9,6)

! *************************************************************************

 strain(:,1) = (/1,0,0,0,0,0,0,0,0/)
 strain(:,2) = (/0,0,0,0,1,0,0,0,0/)
 strain(:,3) = (/0,0,0,0,0,0,0,0,1/)
 strain(:,4) = half*(/0,0,0,0,0,1,0,1,0/)
 strain(:,5) = half*(/0,0,1,0,0,0,1,0,0/)
 strain(:,6) = half*(/0,1,0,1,0,0,0,0,0/)

!Print only the reference system in xml format
 if (option==  -1 .or. option == 1) then

!  convert natom in character
   write (natom,'(I9)') eff_pot%crystal%natom
   if(present(filename)) then
     namefile=trim(filename)
   else
     namefile='system.xml'
   end if

   call isfile(namefile,'new')

   if (open_file(namefile,message,unit=unit_xml,form="formatted",&
&      status="new",action="write") /= 0) then
     MSG_ERROR(message)
   end if

   write(message, '(a,(80a),a)' ) ch10,&
&    ('=',ii=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   write(message,'(a,a,a)')ch10,&
 &       ' Generation of the xml file for the reference structure in ',trim(namefile)

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  Write header
   WRITE(unit_xml,'("<?xml version=""1.0"" ?>")')
   WRITE(unit_xml,'("<System_definition>")')

   WRITE(unit_xml,'("  <energy>")')
   WRITE(unit_xml,'(E23.14)') (eff_pot%energy)
   WRITE(unit_xml,'("  </energy>")')

   WRITE(unit_xml,'("  <unit_cell units=""bohrradius"">")')
   WRITE(unit_xml,'(3(E23.14))') (eff_pot%crystal%rprimd)
   WRITE(unit_xml,'("  </unit_cell>")')

   WRITE(unit_xml,'("  <epsilon_inf units=""epsilon0"">")')
   WRITE(unit_xml,'(3(E23.14))') (eff_pot%harmonics_terms%epsilon_inf)
   WRITE(unit_xml,'("  </epsilon_inf>")')

   WRITE(unit_xml,'("  <elastic units=""hartree"">")')
   WRITE(unit_xml,'(6(E23.14))') (eff_pot%harmonics_terms%elastic_constants)
   WRITE(unit_xml,'("  </elastic>")')

   do ia=1,eff_pot%crystal%natom
     WRITE(unit_xml,'("  <atom mass=""",1F10.5,""" massunits=""atomicmassunit"">")') &
&      eff_pot%crystal%amu(eff_pot%crystal%typat(ia))
     WRITE(unit_xml,'("    <position units=""bohrradius"">")')
     WRITE(unit_xml,'(3(E23.14))') (eff_pot%crystal%xcart(:,ia))
     WRITE(unit_xml,'("    </position>")')
     WRITE(unit_xml,'("    <borncharge units=""abs(e)"">")')
     WRITE(unit_xml,'(3(E23.14))') (eff_pot%harmonics_terms%zeff(:,:,ia))
     WRITE(unit_xml,'("    </borncharge>")')
     WRITE(unit_xml,'("  </atom>")')
   end do
!  Print the Ifc short range for each cell data is array 3*natom*3*natom
!  [ [x1 x2 ....]
!    [y1 y2 ....] for atom 1
!    [z1 z2 ....]
!    [x1 x2 ....]
!    [y1 y2 ....] for atom 2
!    [z1 z2 ....]
!    ....       ]
! Warning : The IFC are print in other order 1,mu,ia,nu,ib which is not fortran way
!           We do like that to because the previous script was in python
!           When you read ifc from XML file with fotran, you have to tranpose the matrix
!
   do irpt=1,eff_pot%harmonics_terms%ifcs%nrpt
     if(any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(1,:,:,:,:,irpt))>tol9)) then
       WRITE(unit_xml,'("  <local_force_constant units=""hartree/bohrradius**2"">")')
       WRITE(unit_xml,'("    <data>")')
       do ia=1,eff_pot%crystal%natom
         do mu=1,3
           do ib=1,eff_pot%crystal%natom
             do  nu=1,3
               WRITE(unit_xml,'(e22.14)', advance="no")&
&                         (eff_pot%harmonics_terms%ifcs%short_atmfrc(1,mu,ia,nu,ib,irpt))
             end do
           end do
           WRITE(unit_xml,'(a)')''
         end do
       end do
       WRITE(unit_xml,'("    </data>")')
       WRITE(unit_xml,'("    <cell>")')
       WRITE(unit_xml,'(3(I4))') (eff_pot%harmonics_terms%ifcs%cell(:,irpt))
       WRITE(unit_xml,'("    </cell>")')
       WRITE(unit_xml,'("  </local_force_constant>")')
     end if
!    Print the IFC total for each cell, data is array 3*natom*3*natom
!    [ [x1 x2 ....]
!      [y1 y2 ....] for atom 1
!      [z1 z2 ....]
!      [x1 x2 ....]
!      [y1 y2 ....] for atom 2
!      [z1 z2 ....]
!      ....       ]
     if(all(abs(eff_pot%harmonics_terms%ifcs%atmfrc(1,:,:,:,:,irpt))<tol9)) then
       if(any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(1,:,:,:,:,irpt))>tol9)) then
         write(message, '(a,a,a,a)' )&
&         ' There is no total range but short range in your effective potential',ch10,&
&         'Action: contact abinit group'
         MSG_BUG(message)
       end if
     else
       WRITE(unit_xml,'("  <total_force_constant units=""hartree/bohrradius**2"">")')
       WRITE(unit_xml,'("    <data>")')
       do ia=1,eff_pot%crystal%natom
         do mu=1,3
           do ib=1,eff_pot%crystal%natom
             do nu=1,3
               WRITE(unit_xml,'(e22.14)', advance="no")&
&                   (eff_pot%harmonics_terms%ifcs%atmfrc(1,mu,ia,nu,ib,irpt))
             end do
           end do
           WRITE(unit_xml,'(a)')''
         end do
       end do
       WRITE(unit_xml,'("    </data>")')
       WRITE(unit_xml,'("    <cell>")')
       WRITE(unit_xml,'(3(I4))') (eff_pot%harmonics_terms%ifcs%cell(:,irpt))
       WRITE(unit_xml,'("    </cell>")')
       WRITE(unit_xml,'("    </total_force_constant>")')
     end if
   end do

   do iqpt=1,eff_pot%harmonics_terms%nqpt
     WRITE(unit_xml,'("  <phonon>")')
     WRITE(unit_xml,'("    <qpoint units=""2pi*G0"">")')
     WRITE(unit_xml,'(3(E23.14))') (eff_pot%harmonics_terms%qpoints(:,iqpt))
     WRITE(unit_xml,'("    </qpoint>")')
     WRITE(unit_xml,'("    <frequencies units=""reciprocal cm"">")')
     WRITE(unit_xml,'(3(e22.14))') (eff_pot%harmonics_terms%phfrq(:,iqpt))
     WRITE(unit_xml,'("    </frequencies>")')
     WRITE(unit_xml,'("    <dynamical_matrix units=""hartree/bohrradius**2"">")')
     do ia=1,eff_pot%crystal%natom
       do mu=1,3
         do ib = 1,eff_pot%crystal%natom
           do nu=1,3
             WRITE(unit_xml,'(e22.14)',advance='no')(eff_pot%harmonics_terms%dynmat(1,nu,ib,mu,ia,iqpt))
           end do
         end do
         WRITE(unit_xml,'(a)')''
       end do
     end do
     WRITE(unit_xml,'("    </dynamical_matrix>")')
     WRITE(unit_xml,'("  </phonon>")')
   end do

! if phonon/forces strain is computed
   jj = 1
   do ii = 1,6
     WRITE(unit_xml,'("  <strain_coupling voigt=""",I2,""">")') ii-1
     WRITE(unit_xml,'("    <strain>")')
     WRITE(unit_xml,'(6(e12.4))') (strain(:,jj))
     WRITE(unit_xml,'("    </strain>")')
     WRITE(unit_xml,'("    <correction_force units=""hartree/bohrradius"">")')
     do ia=1,eff_pot%crystal%natom
       do mu=1,3     
         WRITE(unit_xml,'(e22.14)', advance="no")&
&             (eff_pot%harmonics_terms%internal_strain(ii,mu,ia))
       end do
       WRITE(unit_xml,'(a)')''
     end do
     WRITE(unit_xml,'("    </correction_force>")')
     if (eff_pot%anharmonics_terms%has_elastic3rd) then
       WRITE(unit_xml,'("  <elastic3rd units=""hartree"">")')
       WRITE(unit_xml,'(6(E23.14))') (eff_pot%anharmonics_terms%elastic3rd(ii,:,:))
       WRITE(unit_xml,'("  </elastic3rd>")')
     end if
     if (eff_pot%anharmonics_terms%has_elastic_displ) then
       WRITE(unit_xml,'("    <correction_strain_force units=""hartree/bohrradius"">")')
       do ia=1,eff_pot%crystal%natom
         do mu=1,3
           do nu=1,6
             WRITE(unit_xml,'(e22.14)', advance="no")&
&                 (eff_pot%anharmonics_terms%elastic_displacement(ii,nu,mu,ia))
           end do
         end do
         WRITE(unit_xml,'(a)')''
       end do
     end if
     if (eff_pot%anharmonics_terms%has_strain_coupling) then
       WRITE(unit_xml,'("    </correction_strain_force>")')
       do irpt=1,eff_pot%anharmonics_terms%phonon_strain(ii)%nrpt
         WRITE(unit_xml,'("    <correction_force_constant units=""hartree/bohrradius**2"">")')
         WRITE(unit_xml,'("      <data>")')
         do ia=1,eff_pot%crystal%natom
           do mu=1,3
             do ib=1,eff_pot%crystal%natom
               do  nu=1,3
                 WRITE(unit_xml,'(e22.14)', advance="no")&
&                    (eff_pot%anharmonics_terms%phonon_strain(ii)%atmfrc(1,mu,ia,nu,ib,irpt))
               end do
             end do
             WRITE(unit_xml,'(a)')''
           end do
         end do
         WRITE(unit_xml,'("      </data>")')
         WRITE(unit_xml,'("      <cell>")')
         WRITE(unit_xml,'(3(I4))') (eff_pot%anharmonics_terms%phonon_strain(ii)%cell(:,irpt))
         WRITE(unit_xml,'("      </cell>")')
         WRITE(unit_xml,'("    </correction_force_constant>")')
       end do
     end if
     WRITE(unit_xml,'("    </strain_coupling>")')
     jj = jj + 1
   end do
   WRITE(unit_xml,'("</System_definition>")')

! Close file
   CLOSE(unit_xml)
 end if!end option

!Print only the coefficients into XML file
 if (option==  -1 .or. option == 2) then
   namefile='coefficients.xml'
   if(eff_pot%anharmonics_terms%ncoeff > 0) then
     call polynomial_coeff_writeXML(eff_pot%anharmonics_terms%coefficients,&
&                                 eff_pot%anharmonics_terms%ncoeff,namefile)
   end if
 end if!end option

end subroutine effective_potential_writeXML
!!***

!!****f* m_effective_potential/effective_potential_writeNETCDF
!! NAME
!! effective_potential_writeNETCDF
!!
!! FUNCTION
!! This routine print the effective potential into netcdf format
!! Several options are available
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filename = the name of output file
!! eff_pot  = structure contains the effective potential
!! option   = option for the format of the xml file
!!            1 print the xml for a system
!!
!! OUTPUT
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_writeNETCDF(eff_pot,option,filename)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_multibinit_dataset, only : multibinit_dataset_type
 use m_nctk
#if defined HAVE_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_writeNETCDF'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: option
  character(len=fnlen),optional,intent(in) :: filename
!arrays
  type(effective_potential_type), intent(in) :: eff_pot

!Local variables-------------------------------
!scalar
 integer :: amu_id,bec_id,ifccell_id,epsinf_id,elastic_id
 integer :: ifc_id,ifcs_id,natom_id,ntypat_id,nrpt_id,npsp_id,typat_id
 integer :: six_id,two_id,xyz_id,znucl_id
 integer :: ncerr,ncid,npsp
 integer :: dimCids(2),dimEids(2),dimIids(6),dimPids(1),dimRids(2),dimXids(2)
 integer :: etotal_id,rprimd_id,xcart_id
 character(len=500) :: message
 character(len=fnlen) :: namefile
!arrays
 real :: strain(9,6)

! *************************************************************************

#if defined HAVE_NETCDF

 strain(:,1) = (/1,0,0,0,0,0,0,0,0/)
 strain(:,2) = (/0,0,0,0,1,0,0,0,0/)
 strain(:,3) = (/0,0,0,0,0,0,0,0,1/)
 strain(:,4) = half*(/0,0,0,0,0,1,0,1,0/)
 strain(:,5) = half*(/0,0,1,0,0,0,1,0,0/)
 strain(:,6) = half*(/0,1,0,1,0,0,0,0,0/)

!Print only the reference system in xml format
 if (option == 1) then

   if(present(filename)) then
     namefile=filename
   else
     namefile='ref.nc'
   end if

   call isfile(namefile,'new')

   write(message,'(a,a,a)')ch10,&
&   ' Generation of the xml file for the reference structure in ',namefile

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  1. Create netCDF file
   ncerr = nf90_create(path=trim(namefile),cmode=NF90_CLOBBER, ncid=ncid)
   NCF_CHECK_MSG(ncerr,"create netcdf history file")

!  2. Define dimensions
   ncerr = nf90_def_dim(ncid,"natom",eff_pot%crystal%natom,natom_id)
   NCF_CHECK_MSG(ncerr," define dimension natom")

   ncerr = nf90_def_dim(ncid,"ntypat",eff_pot%crystal%ntypat,ntypat_id)
   NCF_CHECK_MSG(ncerr," define dimension ntypat")

   ncerr = nf90_def_dim(ncid,"nrpt",eff_pot%harmonics_terms%ifcs%nrpt,nrpt_id)
   NCF_CHECK_MSG(ncerr," define dimension ntypat")

   ncerr = nf90_def_var(ncid, "typat", NF90_INT, natom_id, typat_id)
   NCF_CHECK_MSG(ncerr," define variable typat")

   npsp = size(eff_pot%crystal%znucl)
   if (npsp /= eff_pot%crystal%ntypat) then
     MSG_WARNING("HIST file does not support alchemical mixing!")
   end if
   ncerr = nf90_def_dim(ncid,"npsp",npsp,npsp_id)
   NCF_CHECK_MSG(ncerr," define dimension npsp")

   ncerr = nf90_def_var(ncid, "znucl", NF90_DOUBLE, npsp_id, znucl_id)
   NCF_CHECK_MSG(ncerr," define variable znucl")

   ncerr = nf90_def_dim(ncid,"xyz",3,xyz_id)
   NCF_CHECK_MSG(ncerr," define dimension xyz")

   ncerr = nf90_def_dim(ncid,"six",6,six_id)
   NCF_CHECK_MSG(ncerr," define dimension six")

   ncerr = nf90_def_dim(ncid,"two",2,two_id)
   NCF_CHECK_MSG(ncerr," define dimension two")

!  Dimensions for xcart,xred,fcart,fred and vel
   dimXids = (/ xyz_id, natom_id /)
!  Dimensions for rprimd
   dimRids = (/ xyz_id, xyz_id /)
!  Dimensions for ifc
   dimIids = (/ 2, xyz_id, natom_id, xyz_id, natom_id, nrpt_id /)
!  Dimensions for position
   dimPids = (/ xyz_id /)
!  Dimension for elastic constant
   dimEids = (/six_id,six_id/)
!  Dimension for cell
   dimCids = (/nrpt_id,3/)

!  3. Define variables and their attributes (units and mnemonics)
   call ab_define_var(ncid, (/1/), etotal_id, NF90_DOUBLE,&
&   "energy","Energy of the reference structure","Ha" )

   call ab_define_var(ncid, dimRids, rprimd_id, NF90_DOUBLE,&
&   "rprimd","Real space PRIMitive translations, Dimensional","bohr" )

   call ab_define_var(ncid, dimRids, epsinf_id, NF90_DOUBLE,&
&   "epsilon_inf","Dielectric tensor, Dimensional","epsilon_inf" )

   call ab_define_var(ncid, dimEids, elastic_id, NF90_DOUBLE,&
&   "elastic","Elastic Constants, Dimensional","Ha" )

   call ab_define_var(ncid, dimRids, bec_id, NF90_DOUBLE,&
&   "bec","Born Effective Charges, Dimensional","abs(e)" )

   call ab_define_var(ncid, dimXids, xcart_id, NF90_DOUBLE,&
&   "xcart","vectors (X) of atom positions in CARTesian coordinates","bohr" )

   call ab_define_var(ncid, dimIids, ifcs_id, NF90_DOUBLE,&
&   "IFCs","Interatomic Forces Constantes in real spaces (short range), Dimensional","Hatree/bohr**2" )

   call ab_define_var(ncid, dimIids, ifc_id, NF90_DOUBLE,&
&   "IFC","Interatomic Forces Constantes in real spaces (total range), Dimensional","Hatree/bohr**2" )

   call ab_define_var(ncid, dimCids, ifccell_id, NF90_DOUBLE,&
&   "cell","cell for the ifc, Dimensional","Dimensionless" )

   call ab_define_var(ncid, [ntypat_id], amu_id, NF90_DOUBLE,&
&   "amu","Masses of each type of atom in atomic mass units", "" )

!  4. End define mode
   ncerr = nf90_enddef(ncid)
   NCF_CHECK_MSG(ncerr," end define mode")

!  5. Write variables
   ncerr = nf90_put_var(ncid,etotal_id, eff_pot%energy)
   NCF_CHECK_MSG(ncerr," write variable energy")

   ncerr = nf90_put_var(ncid,rprimd_id, eff_pot%crystal%rprimd)
   NCF_CHECK_MSG(ncerr," write variable rprimd")

   ncerr = nf90_put_var(ncid,epsinf_id, eff_pot%harmonics_terms%epsilon_inf)
   NCF_CHECK_MSG(ncerr," write variable epsilon_inf")

   ncerr = nf90_put_var(ncid,elastic_id , eff_pot%harmonics_terms%elastic_constants)
   NCF_CHECK_MSG(ncerr," write variable elastic_constant")

   ncerr = nf90_put_var(ncid, bec_id, eff_pot%harmonics_terms%zeff)
   NCF_CHECK_MSG(ncerr," write variable bec")

   ncerr = nf90_put_var(ncid,xcart_id, eff_pot%crystal%xcart)
   NCF_CHECK_MSG(ncerr," write variable xcart")

   ncerr = nf90_put_var(ncid,ifccell_id, eff_pot%harmonics_terms%ifcs%cell)
   NCF_CHECK_MSG(ncerr," write variable cell")

   ncerr = nf90_put_var(ncid,ifcs_id, eff_pot%harmonics_terms%ifcs%short_atmfrc)
   NCF_CHECK_MSG(ncerr," write variable short ifc")

   ncerr = nf90_put_var(ncid,ifc_id, eff_pot%harmonics_terms%ifcs%atmfrc)
   NCF_CHECK_MSG(ncerr," write variable total ifc")


!  6. Close NetCDF file
   ncerr = nf90_close(ncid)
   NCF_CHECK_MSG(ncerr," close netcdf history file")
 end if

#endif
end subroutine effective_potential_writeNETCDF
!!***


!!****f* m_effective_potential/effective_potential_writeAbiInput
!! NAME
!! effective_potential_writeAbiInput
!!
!! FUNCTION
!! This routine print the effective potential into input of abinit
!! to be able to reproduce the DDB
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eff_pot  = structure contains the effective potential
!! filename = the name of input file
!! strain   = strain structure if need to apply strain into rprim
!!
!! OUTPUT
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_writeAbiInput(eff_pot,filename,strain)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_multibinit_dataset, only : multibinit_dataset_type
 use m_fstrings, only : ftoa,itoa,int2char4

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_writeAbiInput'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  type(strain_type),optional,intent(in) :: strain
  character(len=fnlen),optional,intent(in) :: filename
!arrays
  type(effective_potential_type), intent(in) :: eff_pot

!Local variables-------------------------------
!scalar
 integer :: unit = 20
 character(len=500) :: message
 character(len=fnlen) :: namefile
!arrays
 real(dp) :: xred(3,eff_pot%crystal%natom)
 type(strain_type) :: strain_tmp

! ************************************************************************

 if(present(strain)) then
   strain_tmp = strain
 else
   call strain_init(strain_tmp)
 end if

! try to open the file
 if(present(filename)) then
   namefile=filename
 else
   if(eff_pot%name /='') then
     write(namefile,'(a)') 'structure_'//trim(eff_pot%name)
   else
     write(namefile,'(a)') 'structure'
   end if
   if (strain_tmp%name/='') then
     write(namefile,'(a,a,a,a,a,a,a)') trim(namefile)//"_"//trim(strain_tmp%name)//"_"//&
&        trim(itoa(strain_tmp%direction)),"_"//trim(ftoa(strain_tmp%delta))
   end if

   namefile=trim(namefile)//".in"

 end if

 call isfile(namefile,'new')

 if (open_file(namefile,message,unit=unit,form="formatted",status="new",action="write") /= 0) then
   MSG_ERROR(message)
 end if

  write(message,'(a,a,a,a)')ch10,&
 &   ' Generation of the input file in ',namefile,ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(unit,'("#Abinit Input for DFPT, this file contrains the keyword")')
  write(unit,'("#To run DFPT calculation of")')
  write(unit,'(a)') trim(eff_pot%name)
  if (strain_tmp%direction /= 0) then
    write(unit,'("# With a perturbation  ")',advance="no")
    write(unit,'(a)',advance="no") trim(strain_tmp%name)
    write(unit,'(" in the direction : ")',advance="no")
    write(unit,'(a)',advance='no') trim(itoa(strain_tmp%direction))
    write(unit,'(" with the deformation : ")',advance="no")
    write(unit,'(a)') trim(ftoa(strain_tmp%delta))
  end if

  write(unit,'("")')
  write(unit,'("ndtset 1 jdtset 1 2 3")')
  write(unit,'("")')
  write(unit,'("#DATASET1 GROUND  STATE")')
  write(unit,'("tolwfr1 = 1d-15")')
  write(unit,'(" prtwf1 = 1")')
  write(unit,'(" nline1 = 5")')

  write(unit,'("")')
  write(unit,'("#DATASET2 DDK PERTURBATION")')
  write(unit,'("getwfk2 =  1")')
  write(unit,'("  iscf2 = -3")')
  write(unit,'(" nline2 =  15")')
  write(unit,'("nnsclo2 =  5")')
  write(unit,'("kptopt2 =  2")')
  write(unit,'("  nqpt2 =  1")')
  write(unit,'("   qpt2 =  0 0 0 ")')
  write(unit,'("rfelfd2 =  2")')
  write(unit,'(" rfdir2 =  1 1 1 ")')
  write(unit,'("tolwfr2 =  1.0d-20 ")')
  write(unit,'(" prtwf2 =  1 ")')
  write(unit,'(" ")')

  write(unit,'("#DATASET3 RF")')
  write(unit,'(" getddk3 =  2")')
  write(unit,'(" getwfk3 =  1")')
  write(unit,'("   iscf3 =  7")')
  write(unit,'(" kptopt3 =  2")')
  write(unit,'("   nqpt3 =  1")')
  write(unit,'("    qpt3 =  0 0 0")')
  write(unit,'(" rfphon3 =  1")')
  write(unit,'("rfatpol3 =  1 ")',advance='no')
  write(unit,'(a)') itoa(eff_pot%crystal%natom)
  write(unit,'(" rfelfd3 =  3")')
  write(unit,'(" rfstrs3 =  3")')
  write(unit,'("  rfdir3 =  1 1 1")')
  write(unit,'(" tolvrs3 =  1.0d-8")')
  write(unit,'("")')

  write(unit,'("#STRUCTURE")')
  write(unit,'(" natom = ")',advance='no')
  write(unit,'(a)') itoa(eff_pot%crystal%natom)
  write(unit,'(" znucl =")',advance='no')
  write(unit,'(10(F4.0))') (eff_pot%crystal%znucl)
  write(unit,'("ntypat = ")',advance='no')
  write(unit,'(a)') itoa(eff_pot%crystal%ntypat)
  write(unit,'(" typat = ")',advance='no')
  write(unit,'(10(I2))') (eff_pot%crystal%typat)
  write(unit,'(" acell = 1 1 1")')
  write(unit,'(" rprim  ")')
  write(unit,'(3(F20.10))') (matmul(eff_pot%crystal%rprimd,strain%strain))
  write(unit,'("  xred  ")')
  call xcart2xred(eff_pot%crystal%natom,eff_pot%crystal%rprimd,eff_pot%crystal%xcart,xred)
  write(unit,'(3(F15.10))') (xred)
  write(unit,'(" ")')

  write(unit,'("#SCF")')
  write(unit,'("     ecut = ")')
  write(unit,'("pawecutdg = ")')
  write(unit,'("   ecutsm = ")')
  write(unit,'("   tolvrs = ")')
  write(unit,'("    nband = ")')
  write(unit,'("      ixc = ")')
  write(unit,'("   occopt = ")')
  write(unit,'("    nstep = ")')
  write(unit,'("   kptopt = ")')
  write(unit,'("    ngkpt =  ")')
  write(unit,'("")')
  write(unit,'("    prtwf 0 prtden 0 prtdos 0")')

  close(unit)

end subroutine effective_potential_writeAbiInput
!!***


!****f* m_effective_potential/effective_potential_getForces
!!
!! NAME
!! effective_potential_getForces
!!
!! FUNCTION
!! evaluate the gradient of  of effective potential
!! to calculates the forces in reduced coordinates
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_getForces(eff_pot,fcart,fred,natom,rprimd,xcart,displacement)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_getForces'
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
  real(dp),intent(in),optional :: displacement(3,natom)
  real(dp),intent(out) :: fcart(3,natom),fred(3,natom)
!Local variables-------------------------------
!scalar
  integer,parameter :: master=0
  integer :: ii
!array
  real(dp):: disp_tmp1(3,natom),dummy
  integer :: cell_number(3)
  character(500) :: msg

! *************************************************************************

! Do some checks
  if (natom /= eff_pot%supercell%natom_supercell) then
    write(msg,'(a,I7,a,I7,a)')' The number of atoms is not correct :',natom,&
&   ' in argument istead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(msg)
  end if

! Set the number of cell in the supercell
  cell_number(:) = int(eff_pot%supercell%qphon(:))

  if (any(cell_number <= 0)) then
    write(msg,'(a,a)')' No supercell found for getForces'
    MSG_ERROR(msg)
  end if

! Set to zero
  fred = zero
  fcart = zero

  disp_tmp1(:,:) = zero
! Try to compute the displacement
  if (present(displacement)) then
    disp_tmp1(:,:) = displacement(:,:)
  else
    do ii = 1, natom
      disp_tmp1(:,ii) = xcart(:,ii) - eff_pot%supercell%xcart_supercell(:,ii)
    end do
  end if

! ifc contribution of the forces
  call ifc_contribution(eff_pot,disp_tmp1,dummy,fcart,eff_pot%my_cells,&
&                       eff_pot%my_ncell,eff_pot%my_index_cells,eff_pot%comm_supercell)

! Redistribute the residuale of the forces
  call effective_potential_distributeResidualForces(eff_pot,fcart,eff_pot%supercell%natom_supercell)

! convert forces into reduced coordinates and multiply by -1
  fcart = -1 * fcart
  call fcart2fred(fcart,fred,rprimd,natom)

end subroutine effective_potential_getForces
!!***


!****f* m_effective_potential/effective_potential_evaluate
!!
!! NAME
!! effective_potential_evaluate
!!
!! FUNCTION
!! evaluate the harmonic part of the energy and the forces
!! of a structure with the effective potential
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!! energy =  energy of the structure
!! fcart  =  forces in cartesian coordinates
!! fred   =  forces in reduced coordinates
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                       xcart,displacement,strain_in,external_stress)

  use m_strain

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_evaluate'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(inout) :: rprimd(3,3),xcart(3,natom)
  real(dp),intent(out) :: energy
  real(dp),intent(out) :: fcart(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(out) :: fred(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(out) :: strten(6)
  real(dp),intent(in),optional :: strain_in(6)
  real(dp),intent(in),optional :: displacement(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(in),optional :: external_stress(6)

!Local variables-------------------------------
!scalar
  integer :: alpha,beta,gamma,delta,ii,ia,mu,ncell
  real(dp):: cijk,energy_part
  real(dp):: ucvol
  character(len=500) :: message
  logical :: has_strain = .FALSE.
  logical :: iam_master
  integer, parameter:: master = 0
!array
  type(strain_type) :: strain
  integer  :: supercell(3)
  real(dp) :: disp_tmp(3,eff_pot%supercell%natom_supercell)
!  real(dp) :: external_stress_tmp(6)
  real(dp) :: fcart_part(3,eff_pot%supercell%natom_supercell)
  real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
  real(dp) :: strain_tmp(6),strten_part(6)
  real(dp) :: xred_tmp(3,eff_pot%supercell%natom_supercell)
  real(dp) :: xcart_tmp(3,eff_pot%supercell%natom_supercell)

! *************************************************************************

! Set MPI local varibaless
  iam_master = (eff_pot%me_supercell == master)

! Set variables
  supercell(1:3) = eff_pot%supercell%qphon(1:3)
  ncell          = product(eff_pot%supercell%qphon)

!Check some variables
  if (natom /= eff_pot%supercell%natom_supercell) then
    write(message,'(a,I7,a,I7,a)')' The number of atoms is not correct :',natom,&
&   ' in argument istead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(message)
  end if

  if (present(displacement))then
    if(size(displacement(1,:)) /= eff_pot%supercell%natom_supercell) then
      write(message,'(a,I7,a,I7,a)')' The number of atoms is not correct :',size(displacement(1,:)),&
&      ' in displacement array instead of ',eff_pot%supercell%natom_supercell, ' in supercell'
      MSG_ERROR(message)
    end if
  end if
  do ii=1,3
    if(eff_pot%supercell%qphon(ii)<0.or.eff_pot%supercell%qphon(ii)>20)then
      write(message, '(a,i0,a,i2,a,a,a,i0,a)' )&
&     'eff_pot%supercell%qphon(',ii,') is ',eff_pot%supercell%qphon(ii),&
&     ', which is lower than 0 of superior than 10.',ch10,'Action: correct n_cell(',ii,').'
      MSG_ERROR(message)
    end if
  end do


  write(message, '(a,a,a)' ) ch10,' enter get_energy : Calculation of the energy'
  call wrtout(std_out,message,'COLL')
  write(message, '(a,a,a)' ) ch10,' Calculation of the energy with effective potential'
  call wrtout(ab_out,message,'COLL')

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!--------------------------------------------
! 1 - Set the perturbations and intialisation
!--------------------------------------------

  strain_tmp(:) = zero
! Try to find the strain into the input file
  if (eff_pot%has_strain) then
    do ii=1,3
      strain_tmp(ii) = eff_pot%strain%strain(ii,ii)
    end do
    strain_tmp(4) = strain%strain(2,3) + strain%strain(3,2)
    strain_tmp(5) = strain%strain(3,1) + strain%strain(1,3) 
    strain_tmp(6) = strain%strain(2,1) + strain%strain(1,2)
    has_strain = .TRUE.
! Try to find the strain from argument
  else  if (present(strain_in)) then
    strain_tmp(:) = strain_in(:)
    has_strain = .TRUE.
  else
  ! else => calculation of the strain
    call strain_get(strain,rprim=eff_pot%supercell%rprimd_supercell,rprim_def=rprimd)
     write(message,'(80a)')('-',mu=1,80)
     call wrtout(std_out,message,'COLL')
    call strain_print(strain)
    if (strain%name /= "reference")  then
      has_strain = .TRUE.
      do ii=1,3
        strain_tmp(ii) = strain%strain(ii,ii)
      end do
      strain_tmp(4) = strain%strain(2,3) + strain%strain(3,2)
      strain_tmp(5) = strain%strain(3,1) + strain%strain(1,3) 
      strain_tmp(6) = strain%strain(2,1) + strain%strain(1,2)
    else
      strain_tmp(:) = zero
    end if
  end if
  
  if (present(displacement)) then
    disp_tmp(:,:) = displacement(:,:)
  else
    if(has_strain) then
      call xcart2xred(natom,eff_pot%supercell%rprimd_supercell,&
&                     eff_pot%supercell%xcart_supercell,xred_tmp)
      call xred2xcart(natom, rprimd, xcart_tmp, xred_tmp)
      do ii = 1, natom
        disp_tmp(:,ii) = xcart(:,ii) - xcart_tmp(:,ii)
      end do
    else
      do ii = 1, natom
        disp_tmp(:,ii) = xcart(:,ii) - eff_pot%supercell%xcart_supercell(:,ii)
      end do
    end if
  end if

  energy         = zero
  fcart(:,:)     = zero
  strten(:)      = zero


  write(message, '(2a,80a,2a)' ) ch10,('-',mu=1,80),&
&     ch10,' Components of total energy (in Hartree) :'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

!------------------------------------
! 2 - Transfert the reference energy
!------------------------------------

  energy = eff_pot%energy * ncell
                                             
  write(message, '(a,a,1ES24.16,a)' ) ch10,' Energy of the reference strucure          :',&
&                                          energy,' Hartree'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

!------------------------------------
! 3 - Computation of the IFC part :
!------------------------------------

  energy_part    = zero
  fcart_part(:,:)= zero

  call ifc_contribution(eff_pot,disp_tmp,energy_part,fcart_part,&
&                       eff_pot%my_cells,eff_pot%my_ncell,eff_pot%my_index_cells,&
&                       eff_pot%comm_supercell)
                                      
  write(message, '(a,1ES24.16,a)' ) ' Energy of the ifc part                    :',&
&                                     energy_part,' Hartree'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  if(energy_part < zero.and.eff_pot%anharmonics_terms%ncoeff == zero)then
    write(message, '(8a)' )ch10,&
&        ' --- !WARNING!',ch10,&
&        '        The harmonic part is negative, the simulation will diverge',ch10,&
&        '        if the anharmonic part is not used',ch10,&
&        ' ---'
    call wrtout(std_out,message,"COLL")
  end if

  energy = energy + energy_part
  fcart(:,:)= fcart(:,:) + fcart_part(:,:)

!----------------------------------------------------
! 4 - Computation of the elastic part of the energy :
!----------------------------------------------------

! TO DO ADD external stress
! if external stress is present
!  if (present(external_stress)) then
!    external_stress_tmp(:) = external_stress(:)
!    has_strain = .TRUE.
!  else
!    external_stress_tmp(:) = zero
!  end if

  energy_part    = zero
  fcart_part(:,:)= zero
  strten_part(:) = zero

  call elastic_contribution(eff_pot,disp_tmp,energy_part,fcart_part,&
&                           ncell,strten_part,strain_tmp)

  if(has_strain)then                                      
    write(message, '(a,1ES24.16,a)' ) ' Energy of the elastic part                :',&
&                                       energy_part,' Hartree'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end if

  energy = energy + energy_part
  fcart(:,:) = fcart(:,:)  + fcart_part(:,:)
  strten(:) = strten(:) + strten_part(:)

!------------------------------------
! 5 - Treat 3rd order strain-coupling:
!------------------------------------

  if (eff_pot%has_strainCoupling) then

!   1-Treat 3rd order elastic constants
    if (eff_pot%anharmonics_terms%has_elastic3rd) then
      energy_part    = zero
      strten_part(:) = zero

      do alpha=1,6
        do beta=1,6
          do gamma=1,6
            cijk = ncell*eff_pot%anharmonics_terms%elastic3rd(alpha,beta,gamma)
!           Accumulate energy
            energy_part = energy_part + sixth*cijk*strain_tmp(alpha)*strain_tmp(beta)*strain_tmp(gamma)
!           Accumulate stresses contributions
            strten_part(alpha)=strten_part(alpha)+ half*cijk*strain_tmp(beta)*strain_tmp(gamma)
          end do
        end do
      end do

      write(message, '(a,1ES24.16,a)' ) ' Energy of the 3rd elastics constants      :',&
&                                         energy_part,' Hartree'
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

      energy = energy + energy_part
      strten(:) = strten(:) + strten_part(:)
    end if

!   2-Part due to the internat strain
    if(eff_pot%anharmonics_terms%has_elastic_displ)then
      energy_part    = zero
      strten_part(:) = zero
      fcart_part(:,:)= zero

      ii = 1
      do ia = 1,eff_pot%supercell%natom_supercell
        do mu = 1,3
          do beta=1,6
            do alpha=1,6
              cijk = eff_pot%anharmonics_terms%elastic_displacement(alpha,beta,mu,ii)
!             Accumulte for this atom
              energy_part = energy_part + sixth*cijk*strain_tmp(alpha)*strain_tmp(beta)*disp_tmp(mu,ia)
              fcart_part(mu,ia) = fcart_part(mu,ia)  +  half*cijk*strain_tmp(alpha)*strain_tmp(beta)
              strten_part(alpha) = strten_part(alpha) +  half*cijk*strain_tmp(beta)*disp_tmp(mu,ia)
            end do
          end do
        end do
        ii = ii +1
!       Reset to 1 if the number of atoms is superior than in the initial cell
        if(ii==eff_pot%crystal%natom+1) ii = 1
      end do

      energy = energy +  energy_part
      strten(:) = strten(:) + strten_part(:)
      fcart(:,:)= fcart(:,:)+ fcart_part(:,:) 

      write(message, '(a,1ES24.16,a)' ) ' Energy of the 3rd (elastics-disp coupling):',&
&                                         energy_part,' Hartree'
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    end if

!   3-Part due to the strain-phonon coupling
    if (eff_pot%anharmonics_terms%has_strain_coupling) then
      energy_part    = zero
      strten_part(:) = zero
      fcart_part(:,:)= zero

      call ifcStrainCoupling_contribution(eff_pot,disp_tmp,energy_part,fcart_part,strain_tmp,&
&                                      strten_part,eff_pot%my_cells,eff_pot%my_ncell,&
&                                      eff_pot%my_index_cells,eff_pot%comm_supercell)

      write(message, '(a,1ES24.16,a)' ) ' Energy of the 3rd (strain-phonon coupling):',&
&                                         energy_part,' Hartree'
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

      energy = energy + energy_part
      fcart  = fcart  + fcart_part
      strten = strten + strten_part
    end if


!   4-Treat 4rd order elastic constants
    if (eff_pot%anharmonics_terms%has_elastic4rd) then
      energy_part    = zero
      strten_part(:) = zero

      do alpha=1,6
        do beta=1,6
          do gamma=1,6
            do delta=1,6
            cijk = ncell*eff_pot%anharmonics_terms%elastic4rd(alpha,beta,gamma,delta)
!           Accumulate energy
            energy_part = energy_part + (1/24.)*cijk*strain_tmp(alpha)*strain_tmp(beta)*&
&                                                    strain_tmp(gamma)*strain_tmp(delta)
!           Accumulate stresses contributions
            strten_part(alpha)=strten_part(alpha)+ sixth*cijk*strain_tmp(beta)*strain_tmp(gamma)*&
&                                                             strain_tmp(delta)
            end do
          end do
        end do
      end do


      write(message, '(a,1ES24.16,a)' ) ' Energy of the 4rd elastics constants      :',&
&                                         energy_part,' Hartree'
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

      energy = energy + energy_part
      strten(:) = strten(:) + strten_part(:)
    end if
  end if


!----------------------------------
! 6 - Treat polynomial coefficient:
!----------------------------------

  if(eff_pot%anharmonics_terms%ncoeff > zero)then

    energy_part = zero
    fcart(:,:)  = zero

    call coefficients_contribution(eff_pot,disp_tmp,&
&                                  energy_part,fcart_part,eff_pot%supercell%natom_supercell,&
&                                  eff_pot%anharmonics_terms%ncoeff,&
&                                  eff_pot%my_cells,eff_pot%my_ncell,eff_pot%my_index_cells,&
&                                  eff_pot%comm_supercell)

    write(message, '(a,1ES24.16,a)' ) ' Energy of the fitted coefficient          :',&
&                                      energy_part,' Hartree'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')

    energy = energy + energy_part
    fcart  = fcart  + fcart_part

  end if


!---------------------------------
! 7 - Apply factors
!---------------------------------

! divide stess tensor by ucvol
  strten = strten / ucvol

! convert forces into reduced coordinates and multiply by -1
! Redistribute the residuale of the forces
  call effective_potential_distributeResidualForces(eff_pot,fcart,eff_pot%supercell%natom_supercell)

  fcart = -1 * fcart
  call fcart2fred(fcart,fred,rprimd,natom)

!------------------------------------
! 8 - Final Print:
!------------------------------------
  write(message, '(2a,es21.14)' ) ch10,&
&     '    >>>>>>>>> Etotal= ',energy
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(message, '(a,80a,3a)' ) ch10,('-',mu=1,80),ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')

  write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')
  write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')
  write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')

! Also output the pressure (minus one third the trace of the stress
! tensor.
  write(message, '(a,a,es12.4,a)' ) ch10,&
&   '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&   -(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa]'

  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')
  
  write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(1 1)=',strten(1)*HaBohr3_GPa,&
&   '  sigma(3 2)=',strten(4)*HaBohr3_GPa
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')
  write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(2 2)=',strten(2)*HaBohr3_GPa,&
&   '  sigma(3 1)=',strten(5)*HaBohr3_GPa
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')
  write(message, '(a,1p,e16.8,a,1p,e16.8,a)' ) &
&   '- sigma(3 3)=',strten(3)*HaBohr3_GPa,&
&   '  sigma(2 1)=',strten(6)*HaBohr3_GPa
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,  message,'COLL')
  write(message, '(80a,a)' ) ('-',mu=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

end subroutine effective_potential_evaluate
!!***

!!****f* m_effective_potential/ifc_contribution
!! NAME
!!  ifc_contribution
!!
!! FUNCTION
!!  This fonction compute the harmonic part of the energy
!!  of the supercell in the eff_pot
!! INPUTS
!!  eff_pot = effective potential of the structure
!!            also contain supercell information
!!  disp    = diplacement vector (3, cell1 (atm1 atm2 ...) cell2 (atm1 atm2 ...)...)
!!  ncell   = total number of cell to treat
!!  cells(ncell) = number of the cells into the supercell (1,2,3,4,5)
!!  index_cells(3,ncell) = indexes of the cells into  supercell (-1 -1 -1 ,...,1 1 1)
!!  comm=MPI communicator
!!
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!! PARENT
!!   effective_potential_evaluate
!!
!! CHILDREN
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine ifc_contribution(eff_pot,disp,energy,fcart,cells,ncell,index_cells,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_contribution'
!End of the abilint section

 implicit none

!Arguments -------------------------------
! scalars
  real(dp),intent(out) :: energy
  integer,intent(in) :: ncell
  integer,intent(in) :: comm
! array
  integer,intent(in) ::   cells(ncell),index_cells(ncell,3)
  type(effective_potential_type), intent(in) :: eff_pot
  real(dp),intent(in) :: disp(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(out) :: fcart(3,eff_pot%supercell%natom_supercell)

!Local variables-------------------------------
! scalar
  integer :: i1,i2,i3,ia,ib,icell,ierr,irpt,ii,jj,kk,ll
  integer :: mu,nu
  real(dp):: ifc
! array
  integer :: cell_number(3)
  integer :: cell_atom2(3)
  character(500) :: msg

! *************************************************************************

  cell_number(:) = int(eff_pot%supercell%qphon(:))

  if (any(cell_number <= 0)) then
    write(msg,'(a,a)')' No supercell found for getEnergy'
    MSG_ERROR(msg)
  end if

! Initialisation of variables
  energy   = zero
  fcart(:,:) = zero

  do icell = 1,ncell
    ii = (cells(icell)-1)*eff_pot%crystal%natom
    i1=index_cells(icell,1); i2=index_cells(icell,2); i3=index_cells(icell,3)
    do irpt = 1,eff_pot%harmonics_terms%ifcs%nrpt
!     get the cell of atom2  (0 0 0, 0 0 1...)
      cell_atom2(1) =  (i1-1) + eff_pot%harmonics_terms%ifcs%cell(1,irpt)
      call index_periodic(cell_atom2(1),cell_number(1))
      cell_atom2(2) =  (i2-1) + eff_pot%harmonics_terms%ifcs%cell(2,irpt)
      call index_periodic(cell_atom2(2),cell_number(2))
      cell_atom2(3) =  (i3-1) + eff_pot%harmonics_terms%ifcs%cell(3,irpt)
      call index_periodic(cell_atom2(3),cell_number(3))
!     index of the second atom in the displacement array
      jj = cell_atom2(1)*cell_number(2)*cell_number(3)*eff_pot%crystal%natom+&
&          cell_atom2(2)*cell_number(3)*eff_pot%crystal%natom+&
&          cell_atom2(3)*eff_pot%crystal%natom
      do ib = 1, eff_pot%crystal%natom
        ll = jj + ib
        do nu=1,3
          do ia = 1, eff_pot%crystal%natom
            kk = ii + ia
            do mu=1,3
              ifc = eff_pot%harmonics_terms%ifcs%atmfrc(1,mu,ia,nu,ib,irpt)
!             accumule energy
              energy =  energy + disp(mu,kk)*disp(nu,ll)*ifc
!             accumule forces
              fcart(mu,kk) = fcart(mu,kk)+ disp(nu,ll)*ifc
            end do
          end do
        end do
      end do
    end do
  end do

  energy = half * energy

! MPI_SUM
  call xmpi_sum(energy, comm, ierr)
  call xmpi_sum(fcart , comm, ierr)

end subroutine ifc_contribution
!!***

!!****f* m_effective_potential/ifcStrainCoupling_contribution
!! NAME
!!  ifcStrainCoupling_contribution
!!
!! FUNCTION
!!  This fonction compute the harmonic part of the energy
!!  of the supercell in the eff_pot
!! INPUTS
!!  eff_pot = effective potential of the structure
!!            also contain supercell information
!!  disp    = diplacement vector (3, cell1 (atm1 atm2 ...) cell2 (atm1 atm2 ...)...)
!!  ncell   = total number of cell to treat
!!  cells(ncell) = number of the cells into the supercell (1,2,3,4,5)
!!  index_cells(3,ncell) = indexes of the cells into  supercell (-1 -1 -1 ,...,1 1 1)
!!  comm=MPI communicator 
!!
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!! PARENT
!!   effective_potential_evaluate
!!
!! CHILDREN
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine ifcStrainCoupling_contribution(eff_pot,disp,energy,fcart,strain,strten,&
&                                         cells,ncell,index_cells,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifcStrainCoupling_contribution'
!End of the abilint section

 implicit none

!Arguments -------------------------------
! scalars
  real(dp),intent(out) :: energy
  integer,intent(in) :: ncell
  integer,intent(in) :: comm
! array
  integer,intent(in) ::   cells(ncell),index_cells(ncell,3)
  type(effective_potential_type), intent(in) :: eff_pot
  real(dp),intent(in) :: disp(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(out) :: fcart(3,eff_pot%supercell%natom_supercell)
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
  integer :: cell_number(3)
  integer :: cell_atom2(3)
  character(500) :: msg

! *************************************************************************

  cell_number(:) = int(eff_pot%supercell%qphon(:))
  if (any(cell_number <= 0)) then
    write(msg,'(a,a)')' No supercell found for getEnergy'
    MSG_ERROR(msg)
  end if

! Initialisation of variables
  energy   = zero
  fcart(:,:) = zero
  strten(:) = zero

  do icell = 1,ncell
    ii = (cells(icell)-1)*eff_pot%crystal%natom
    i1=index_cells(icell,1); i2=index_cells(icell,2); i3=index_cells(icell,3)
    do alpha=1,6
      do irpt = 1,eff_pot%anharmonics_terms%phonon_strain(alpha)%nrpt
!       get the cell of atom2  (0 0 0, 0 0 1...)
        cell_atom2(1) =  (i1-1) + eff_pot%anharmonics_terms%phonon_strain(alpha)%cell(1,irpt)
        call index_periodic(cell_atom2(1),cell_number(1))
        cell_atom2(2) =  (i2-1) + eff_pot%anharmonics_terms%phonon_strain(alpha)%cell(2,irpt)
        call index_periodic(cell_atom2(2),cell_number(2))
        cell_atom2(3) =  (i3-1) + eff_pot%anharmonics_terms%phonon_strain(alpha)%cell(3,irpt)
        call index_periodic(cell_atom2(3),cell_number(3))
!       index of the second atom in the displacement array
        jj = cell_atom2(1)*cell_number(2)*cell_number(3)*eff_pot%crystal%natom+&
&            cell_atom2(2)*cell_number(3)*eff_pot%crystal%natom+&
&            cell_atom2(3)*eff_pot%crystal%natom
        do ib = 1, eff_pot%crystal%natom
          ll = jj + ib
          do nu=1,3
            do ia = 1, eff_pot%crystal%natom
              kk = ii + ia
              do mu=1,3
                ifc = eff_pot%anharmonics_terms%phonon_strain(alpha)%atmfrc(1,mu,ia,nu,ib,irpt)
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

end subroutine ifcStrainCoupling_contribution
!!***


!!****f* m_effective_potential/elastic_contribution
!! NAME
!!  elastic_contribution
!!
!! FUNCTION
!! Compute the energy related to the application of strain
!!
!! INPUTS
!! eff_pot = effective potential structure
!! ncell   = number of cell
!! strain(6) =  first strain to apply
!! external_strees(6) =  external stress to apply
!!
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
!!
subroutine elastic_contribution(eff_pot,disp,energy,fcart,ncell,strten,strain,&
&                               external_stress)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elastic_contribution'
!End of the abilint section

 real(dp),intent(out):: energy
 integer, intent(in) :: ncell
! array
 type(effective_potential_type), intent(in) :: eff_pot
 real(dp),intent(out):: strten(6)
 real(dp),intent(out):: fcart(3,eff_pot%supercell%natom_supercell)
 real(dp),intent(in) :: disp(3,eff_pot%supercell%natom_supercell)
 real(dp),intent(in) :: strain(6)
 real(dp),optional,intent(in) :: external_stress(6)

 !Local variables-------------------------------
! scalar
 integer :: ia,ii,mu,alpha,beta
 real(dp):: cij
! array
! *************************************************************************

 energy = zero
 fcart = zero
 strten = zero

!1- Part due to elastic constants
 do alpha=1,6
   do beta=1,6
     cij = ncell*eff_pot%harmonics_terms%elastic_constants(alpha,beta)
     energy = energy + half*cij*strain(alpha)*strain(beta)
     strten(alpha) = strten(alpha) + cij*strain(beta)
   end do
 end do

!2-Part due to the internat strain
 ii = 1
 do ia = 1,eff_pot%supercell%natom_supercell
   do mu = 1,3
     do alpha=1,6
       cij = eff_pot%harmonics_terms%internal_strain(alpha,mu,ii)
!      Accumulte for this atom
       energy = energy + half*cij*strain(alpha)*disp(mu,ia)
       fcart(mu,ia)  = fcart(mu,ia)  + cij*strain(alpha)
       strten(alpha) = strten(alpha) + cij*disp(mu,ia)
     end do
   end do
   ii = ii +1
!  Reset to 1 if the number of atoms is superior than in the initial cell
   if(ii==eff_pot%crystal%natom+1) ii = 1
 end do

!3- Part due to extenal stress
  if(present(external_stress)) then
    energy = energy - dot_product(external_stress, strain)
    strten(:) = strten(:) - external_stress(:)
  end if

end subroutine  elastic_contribution
!!***

!!****f* m_effective_potential/coefficients_contribution
!! NAME
!!  coefficients_contribution
!!
!! FUNCTION
!! Compute the energy related to the coefficients from
!! fitted polynome
!!
!! INPUTS
!!  eff_pot = effective potential of the structure
!!            also contain supercell information
!!  ncoeffs   = number of coefficients
!!  ncell   = total number of cell to treat
!!  cells(ncell) = number of the cells into the supercell (1,2,3,4,5)
!!  index_cells(3,ncell) = indexes of the cells into  supercell (-1 -1 -1 ,...,1 1 1)
!!  comm=MPI communicator
!! 
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
!!
subroutine coefficients_contribution(eff_pot,disp,energy,fcart,natom,ncoeff,cells,ncell,&
&                                    index_cells,comm)

!Arguments ------------------------------------
! scalar

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'coefficients_contribution'
!End of the abilint section

  real(dp),intent(out):: energy
  integer, intent(in) :: natom,ncell,ncoeff
  integer, intent(in) :: comm
! array
  integer,intent(in) ::   cells(ncell),index_cells(ncell,3)
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(out):: fcart(3,natom)
  real(dp),intent(in) :: disp(3,natom)

 !Local variables-------------------------------
! scalar
  integer :: i1,i2,i3,ia1,ib1,ia2,ib2,idir1,idir2,ierr,ii
  integer :: icoeff,iterm,idisp1,idisp2,icell,power,weight
  real(dp):: coeff,disp1,disp2,tmp1,tmp2
! array
!  integer :: cell_atoma1(3),cell_atoma2(3)
  integer :: cell_atomb1(3),cell_atomb2(3),cell_number(3)

! *************************************************************************

! Initialisation of variables
  cell_number(:) = int(eff_pot%supercell%qphon(:))
  energy   = zero
  fcart(:,:) = zero

  do icell = 1,ncell
    ii = (cells(icell)-1)*eff_pot%crystal%natom
    i1=index_cells(icell,1); i2=index_cells(icell,2); i3=index_cells(icell,3)
!   Loop over coefficient
    do icoeff=1,ncoeff
!     Set the coefficient
      coeff = eff_pot%anharmonics_terms%coefficients(icoeff)%coefficient
!     Loop over term of this coefficient
      do iterm=1,eff_pot%anharmonics_terms%coefficients(icoeff)%nterm
!       Set the weight of the term
        weight =eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%weight
        tmp1 = one
!       Loop over displacement
        do idisp1=1,eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%ndisp
          tmp2 = one
          idir1 = eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%direction(idisp1)

!         indexes of the cell of the atom a (with PBC)
!           cell_atoma1 = eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%cell(:,1,idisp1)
!           if(any(cell_atoma1(:) /=0)) then
!             cell_atoma1(1) =  (i1-1) + cell_atoma1(1)
!             call index_periodic(cell_atoma1(1),cell_number(1))
!             cell_atoma1(2) =  (i2-1) + cell_atoma1(2)
!             call index_periodic(cell_atoma1(2),cell_number(2))
!             cell_atoma1(3) =  (i3-1) + cell_atoma1(3)
!             call index_periodic(cell_atoma1(3),cell_number(3))
! !           index of the first atom (position in the supercell and direction)
! !           if the cell of the atom a is not 0 0 0 (may happen)
!             ia1 = cell_atoma1(1)*cell_number(2)*cell_number(3)*eff_pot%crystal%natom+&
! &                 cell_atoma1(2)*cell_number(3)*eff_pot%crystal%natom+&
! &                 cell_atoma1(3)*eff_pot%crystal%natom+&
! &                 eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
!           else
!         index of the first atom (position in the supercell and direction)
            ia1 = ii + eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
!          end if
!         indexes of the cell of the atom b  (with PBC)
          cell_atomb1 = eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%cell(:,2,idisp1)
          if(cell_atomb1(1)/=0.or.cell_atomb1(2)/=0.or.cell_atomb1(3)/=0) then
            cell_atomb1(1) =  (i1-1) + cell_atomb1(1)
            call index_periodic(cell_atomb1(1),cell_number(1))
            cell_atomb1(2) =  (i2-1) + cell_atomb1(2)
            call index_periodic(cell_atomb1(2),cell_number(2))
            cell_atomb1(3) =  (i3-1) + cell_atomb1(3)
            call index_periodic(cell_atomb1(3),cell_number(3))

!           index of the second atom in the (position in the supercell) 
!           if the cell of the atom b is not 0 0 0 (may happen)
            ib1 = cell_atomb1(1)*cell_number(2)*cell_number(3)*eff_pot%crystal%natom+&
&                 cell_atomb1(2)*cell_number(3)*eff_pot%crystal%natom+&
&                 cell_atomb1(3)*eff_pot%crystal%natom+&
&                 eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
          else
            ib1 = ii + eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
          end if

!         Get the displacement for the both atoms
          disp1 = disp(idir1,ia1)
          disp2 = disp(idir1,ib1)

!         Set the power of the displacement:
          power = eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%power(idisp1)

!         Accumulate energy fo each displacement (\sum ((A_x-O_x)^Y(A_y-O_c)^Z))
          tmp1 = tmp1 * (disp1-disp2)**power

!         Accumulate forces for each displacement (\sum (Y(A_x-O_x)^Y-1(A_y-O_c)^Z+...))
          tmp2 = tmp2 *  power*(disp1-disp2)**(power-1)

          do idisp2=1,eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%ndisp
            if(idisp1/=idisp2) then

              idir2 = eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%direction(idisp2)

              ! if(any(cell_atoma2(:) /=0)) then
              !   cell_atoma2(1) =  (i1-1) + cell_atoma2(1)
              !   call index_periodic(cell_atoma2(1),cell_number(1))
              !   cell_atoma2(2) =  (i2-1) + cell_atoma2(2)
              !   call index_periodic(cell_atoma2(2),cell_number(2))
              !   cell_atoma2(3) =  (i3-1) + cell_atoma2(3)
              !   call index_periodic(cell_atoma2(3),cell_number(3))
!               index of the first atom (position in the supercell and direction)
!               if the cell of the atom a is not 0 0 0 (may happen)
!                ia2 = cell_atoma2(1)*cell_number(2)*cell_number(3)*eff_pot%crystal%natom+&
!&                     cell_atoma2(2)*cell_number(3)*eff_pot%crystal%natom+&
!&                     cell_atoma2(3)*eff_pot%crystal%natom+&
!&                     eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
!              else
!               index of the first atom (position in the supercell and direction)
                ia2 = ii + eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
!              end if

              cell_atomb2= eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%cell(:,2,idisp2)

              if(cell_atomb2(1)/=0.or.cell_atomb2(2)/=0.or.cell_atomb2(3)/=0) then
!             indexes of the cell2 (with PBC)
                cell_atomb2(1) =  (i1-1) + cell_atomb2(1)
                call index_periodic(cell_atomb2(1),cell_number(1))
                cell_atomb2(2) =  (i2-1) + cell_atomb2(2)
                call index_periodic(cell_atomb2(2),cell_number(2))
                cell_atomb2(3) =  (i3-1) + cell_atomb2(3)
                call index_periodic(cell_atomb2(3),cell_number(3))

!               index of the second atom in the (position in the supercell) 
                ib2 = cell_atomb2(1)*cell_number(2)*cell_number(3)*eff_pot%crystal%natom+&
&                     cell_atomb2(2)*cell_number(3)*eff_pot%crystal%natom+&
&                     cell_atomb2(3)*eff_pot%crystal%natom+&
&                     eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
              else
                ib2 = ii + eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
              end if

              disp1 = disp(idir2,ia2)
              disp2 = disp(idir2,ib2)
!             Set the power of the displacement:
              power = eff_pot%anharmonics_terms%coefficients(icoeff)%terms(iterm)%power(idisp2)

              tmp2 = tmp2 * (disp1-disp2)**power
            end if
          end do

!         Accumule  forces
          fcart(idir1,ia1) =  fcart(idir1,ia1)  + coeff * weight * tmp2

        end do

!       accumule energy
        energy = energy +  coeff * weight * tmp1

      end do
    end do
  end do


! MPI_SUM
  call xmpi_sum(energy, comm, ierr)
  call xmpi_sum(fcart , comm, ierr)

end subroutine coefficients_contribution
!!***


!****f* m_effective_potential/effective_potential_distributeResidualForces
!!
!! NAME
!! effective_potential_distributeResidualForces
!!
!! FUNCTION
!! Distribute the residual forces in a weighted manner
!!
!! INPUTS
!! natom   = number of atoms
!! eff_pot = effective potential structure
!!
!! OUTPUT
!! fcart   = forces in cartesian coordinates
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_distributeResidualForces(eff_pot,fcart,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_distributeResidualForces'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(inout) :: fcart(3,natom)
!Local variables-------------------------------
!scalar
  real(dp):: mass_ia,sum_mass
  integer :: ia
!array
  real(dp):: sum_f(3)

! *************************************************************************

  sum_f(1) = sum(fcart(1,:))
  sum_f(2) = sum(fcart(2,:))
  sum_f(3) = sum(fcart(3,:))
  sum_mass = zero

  do ia=1,natom
    sum_mass = sum_mass + eff_pot%crystal%amu(eff_pot%supercell%typat_supercell(ia))
  end do

  do ia=1,natom
    mass_ia = eff_pot%crystal%amu(eff_pot%supercell%typat_supercell(ia))
    fcart(:,ia) = fcart(:,ia) - (mass_ia/sum_mass) * sum_f(:)
  end do


end subroutine effective_potential_distributeResidualForces
!!***

!****f* m_effective_potential/effective_potential_getDeltaEnergy
!!
!! NAME
!! effective_potential_getDeltaEnergy
!!
!! FUNCTION
!! evaluate the energy due to 1 atomic displacement
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_monte_carlo
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_getDeltaEnergy(eff_pot,energy,iatom,idir,natom,rprimd,displacement)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_getDeltaEnergy'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: iatom,idir,natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in),optional :: displacement(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(out) :: energy
!Local variables-------------------------------
!scalar
  real(dp) :: ifc_contribution
  character(len=500) :: message
!array
! *************************************************************************

!Check some variables
  if (natom /= eff_pot%supercell%natom_supercell) then
    write(message,'(a,I7,a,I7,a)')' The number of atoms is not correct :',natom,&
&   ' in argument istead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(message)
  end if

!------------------------------------
! 2 - Computation of the harmonic part (IFC) of the energy :
!------------------------------------

  ifc_contribution =  0

  energy = energy + ifc_contribution
!TEST_AM
!------------------------------------
! 2 - To be continued
!------------------------------------
!TEST_AM

end subroutine effective_potential_getDeltaEnergy
!!***

!!****f* m_effective_potential/index_periodic
!! NAME
!!  delta ernegy
!!
!! FUNCTION
!!
!! INPUTS
!! index  = index of the cell into the supercell
!! n_cell = number of total cell
!!
!! OUTPUT
!!
!! SOURCE

subroutine index_periodic(index,n_cell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'index_periodic'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
  integer, intent(inout)  :: index
  integer, intent(inout) :: n_cell
!Local variables ---------------------------------------
! *********************************************************************

  do while (index > n_cell-1)
    index = index - n_cell
  end do
  do while (index < 0)
    index = index + n_cell
  end do

end subroutine index_periodic
!!***

!!****f* m_effective_potential/find_bound
!! NAME
!!  find_bound
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine find_bound(min,max,n_cell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_bound'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
  integer, intent(inout) :: min,max
  integer, intent(in) :: n_cell
!Local variables ---------------------------------------
  if(abs(max)>abs(min)) then
    max=(n_cell)/2; min=-max;  if(mod(n_cell,2)==0) max = max -1
  else
    min=-(n_cell)/2; max=-min; if(mod(n_cell,2)==0)  min= min +1
  end if

! *********************************************************************
end subroutine find_bound
!!***


!!****f* m_effective_potential/equal
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

pure function effective_potential_compare(e1,e2) result (res)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_compare'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  type(effective_potential_type), intent(in) :: e1,e2
  logical :: res
! *************************************************************************
  res = .false.
  if(e1%crystal%natom==e2%crystal%natom.and.&
&     e1%harmonics_terms%ifcs%nrpt==e2%harmonics_terms%ifcs%nrpt.and.&
&     e1%crystal%ntypat==e2%crystal%ntypat.and.&
&     e1%harmonics_terms%nqpt==e2%harmonics_terms%nqpt.and.&
&     e1%energy==e2%energy.and.&
&     e1%crystal%ucvol==e2%crystal%ucvol) then
    res = .true.
  end if
end function effective_potential_compare
!!***


!TEST_AM_EXPERIMENTAL SECTION
!****f* m_effective_potential/effective_potential_effpot2ddb
!!
!! NAME
!! effective_potential_effpot2ddb
!!
!! FUNCTION
!! Convert eff_pot into ddb structure
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!! ddb   = ddb with all information
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_effpot2ddb(ddb,crystal,eff_pot,n_cell,nph1l,option,qph1l)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_effpot2ddb'
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: nph1l,option
!array
  integer,intent(in) :: n_cell(3)
  real(dp),intent(in):: qph1l(3,nph1l)
  type(effective_potential_type),intent(inout) :: eff_pot
  type(ddb_type),intent(out) :: ddb
  type(crystal_t),intent(out) :: crystal
!Local variables-------------------------------
!scalar
  integer :: ii,jj,msym
  real(dp):: ucvol

! type(anaddb_dataset_type) :: inp
!array
  real(dp) :: gmet(3,3),rmet(3,3)
  real(dp) :: gprimd(3,3),rprimd(3,3)
  real(dp),allocatable :: xred(:,:)
  character :: title(eff_pot%crystal%ntypat)
  integer,allocatable :: symrel(:,:,:),symafm(:)
  real(dp),allocatable :: tnons(:,:)

! *************************************************************************

  ! Number of 2dte blocks in present object
!  integer,allocatable :: flg(:,:)
  ! flg(msize,nblok)
  ! flag to indicate presence of a given block
!  integer,allocatable :: typ(:)
  ! typ(nblok)
  ! type of each block - ddk, dde, phonon etc...
!  real(dp),allocatable :: amu(:)
  ! amu(ntypat)
  ! mass of the atoms (atomic mass unit)
!  real(dp),allocatable :: nrm(:,:)
  ! nrm(3,nblok)
  ! norm of the q-points for each block - can be 0 to indicate a direction of approach to gamma
!  real(dp),allocatable :: qpt(:,:)
  ! qpt(9,nblok)
  ! q-point vector in reciprocal space (reduced lattice coordinates) for each block
!  real(dp),allocatable :: val(:,:,:)
  ! val(2,msize,nblok)
  ! values of the second energy derivatives in each block

! Useless value
   ddb%nblok = -1

  !option = 1 just print ddb for 1 1 1 cell
  if(option==1) then
!   Compute different matrices in real and reciprocal space, also
!   checks whether ucvol is positive.
    call metric(gmet,gprimd,-1,rmet,eff_pot%crystal%rprimd,ucvol)

!   Convert to rprimd
    do ii=1,3
      do jj=1,3
        rprimd(ii,jj)=eff_pot%crystal%rprimd(ii,jj)
      end do
    end do

!   Obtain reciprocal space primitive transl g from inverse trans of r
!   (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
    call matr3inv(rprimd,gprimd)

!   transfert basic values
    ddb%natom  = eff_pot%crystal%natom
    ddb%mpert  = ddb%natom+6
    ddb%msize  = 3*ddb%mpert*3*ddb%mpert;
    ddb%ntypat = eff_pot%crystal%ntypat
    ddb%occopt = 3 ! default value
    ddb%prtvol = 0 ! default value
    ddb%rprim  = rprimd ! dimensioless real space primitive vectors
    ddb%gprim  = gprimd ! dimensioless reciprocal space primitive vectors
    ddb%acell  = one
        msym   = 1
!  Setup crystal type
    ABI_ALLOCATE(xred,(3,ddb%natom))
!    call xcar2xred(ddb%natom,eff_pot%crystal%rprimd,eff_pot%crystal%xcart,xred)
!Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
    ABI_ALLOCATE(symrel,(3,3,msym))
    ABI_ALLOCATE(symafm,(msym))
    ABI_ALLOCATE(tnons,(3,msym))

    call crystal_init(ddb%amu,crystal,1,ddb%natom,size(eff_pot%crystal%znucl),eff_pot%crystal%ntypat,1,&
&       eff_pot%crystal%rprimd,eff_pot%crystal%typat,xred,eff_pot%crystal%znucl,&
&       eff_pot%crystal%znucl,0,.FALSE.,.FALSE.,title)!,&
!&       symrel=symrel,tnons=tnons,symafm=symafm)
    call crystal_print(crystal)
    stop
!TEST_AM
    ABI_DEALLOCATE(symrel)
    ABI_DEALLOCATE(symafm)
    ABI_DEALLOCATE(tnons)

    ABI_DEALLOCATE(xred)

   else  if (option==2) then
!   Compute different matrices in real and reciprocal space, also
!   checks whether ucvol is positive.
    call metric(gmet,gprimd,-1,rmet,eff_pot%supercell%rprimd_supercell,ucvol)

!   Convert to rprim (dimensionless)
    do ii=1,3
      do jj=1,3
        rprimd(ii,jj)=eff_pot%supercell%rprimd_supercell(ii,jj)
      end do
    end do

!   Obtain reciprocal space primitive transl g from inverse trans of r
!   (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
    call matr3inv(rprimd,gprimd)

!   transfert basic values
    ddb%natom  = eff_pot%supercell%natom_supercell
    ddb%ntypat = eff_pot%crystal%ntypat
    ddb%mpert  = ddb%natom+6
    ddb%msize  = 3*ddb%mpert*3*ddb%mpert;
    ddb%occopt = 3 ! default value
    ddb%prtvol = 0 ! default value
    ddb%rprim  = rprimd ! dimensioless real space primitive vectors
    ddb%gprim  = gprimd ! dimensioless reciprocal space primitive vectors
    ddb%acell  = one

   end if
!TEST_AM
    !print*,"natom ",ddb%natom
    !print*,"ntypat",ddb%ntypat
    !print*,"mpert",ddb%mpert
    !print*,"msize",ddb%msize
    !print*,"occopt",ddb%occopt
    !print*,"prtvol",ddb%prtvol
    !print*,"rprim",ddb%rprim
    !print*,"gprim",ddb%gprim
    !print*,"acell",ddb%acell
!TEST_AM

 end subroutine effective_potential_effpot2ddb
!!***


!****f* m_effective_potential/effective_potential_printPDOS
!!
!! NAME
!! effective_potential_printPDOS
!!
!! FUNCTION
!! Apply the acoustic sum rule on the effective potential
!!
!! INPUTS
!! eff_pot = effective potential structure
!! option  = 0 (default) do nothing
!!         = 1 print PHFRQ for specific qgrid (need nph1l and qph1l)
!!         = 2 print PHFRQ for supercell (q=gamma) (need nph1l and qph1l and ncell)
!! OUTPUT
!! eff_pot
!!
!! PARENTS
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE

subroutine effective_potential_printPDOS(eff_pot,filename,n_cell,nph1l,option,qph1l)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_printPDOS'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: nph1l,option
!array
  integer,intent(in) :: n_cell(3)
  real(dp),intent(in):: qph1l(3,nph1l)
  type(effective_potential_type),intent(inout) :: eff_pot
  character(len=fnlen),intent(in) :: filename
!Local variables-------------------------------
!scalar
 integer :: lenstr
 real(dp) :: tcpui,twalli
 character(len=strlen) :: string
!array
 type(crystal_t) :: Crystal
 type(anaddb_dataset_type) :: inp
 type(ddb_type) :: ddb
 type(asrq0_t) :: asrq0

! *************************************************************************

 if (option > 0) then

!  First: transfer into ddb structure:
   call effective_potential_effpot2ddb(ddb,Crystal,eff_pot,n_cell,nph1l,option,qph1l)

!  Setup fake anaddb_dataset
   string = ''
   lenstr = 0
   call invars9(inp,lenstr,ddb%natom,string)
!  fill it with multibinit_dataset values
   inp%prt_ifc = 1
   inp%ifcflag = 1
   inp%qph1l   = qph1l
   inp%nph1l   = nph1l

   ! In case the interatomic forces are not calculated, the
   ! ASR-correction (asrq0%d2asr) has to be determined here from the Dynamical matrix at Gamma.
   if (inp%ifcflag == 0) then
     asrq0 = ddb_get_asrq0(ddb, inp%asr, inp%rfmeth, crystal%xcart)
   end if

  !MG: Note that I'm passing xmpi_comm_self here.
  call mkphbs(eff_pot%harmonics_terms%ifcs,Crystal,inp,ddb,asrq0,filename,&
&  tcpui,twalli,eff_pot%harmonics_terms%zeff, xmpi_comm_self)

  call asrq0_free(asrq0)

 end if

 end subroutine effective_potential_printPDOS
!!***
!TEST_AM_END_EXPERIMENTAL SECTION

end module m_effective_potential
!!***
