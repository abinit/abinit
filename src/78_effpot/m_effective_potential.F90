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

module m_effective_potential

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_strain
 use m_ifc
 use m_io_tools, only : open_file
 use m_phonon_supercell
 use m_phonons
 use m_ddb
 use m_crystal,        only : crystal_t, crystal_init, crystal_free, crystal_t,crystal_print
 use m_anaddb_dataset, only : anaddb_dataset_type, anaddb_dtset_free, outvars_anaddb, invars9
 use m_dynmat,         only : make_bigbox,q0dy3_apply, q0dy3_calc

 implicit none

 public :: effective_potential_applySumRule
 public :: effective_potential_effpot2ddb
 public :: effective_potential_effpot2dynmat
 public :: effective_potential_free
 public :: effective_potential_generateSupercell
 public :: effective_potential_getDeltaEnergy
 public :: effective_potential_getEnergy
 public :: effective_potential_getForces
 public :: effective_potential_init
 public :: effective_potential_print
 public :: effective_potential_printPDOS
 public :: effective_potential_printSupercell
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

   integer :: natom                                 
!     Number of atoms in primitive cell

   integer :: ntypat                                
!     Number of type of atoms

   integer :: nph1l
!     Number of wavevectors for phonon 

   real(dp):: energy                                
!     Energy of the system (Hatree)

   real(dp):: ucvol                                 
!     Unit cell volume

   real(dp) :: acell(3)
!     Lattice vectors for supercell (Bohr)

   real(dp) :: rprimd(3,3)           
!     Lattice vectors for supercell (Bohr)

   real(dp) :: epsilon_inf(3,3)                     
!     epsilon_inf(3,3)
!     Dielectric tensor

   type(supercell_type) :: supercell                 
!     super cell type
!     Store all the information of the suppercell

   real(dp) :: elastic_constants(6,6)
!     elastic_constant(6,6)
!     Elastic tensor (Hartree)

   real(dp) :: internal_stress(6)
!     internal_stress(6)
!     stress tensor of the structure

   real(dp) :: external_stress(6)
!     internal_stress(6)
!     stress tensor of the structure

   logical :: has_3rd
!     True : the 3rd order derivative is computed 

   type(strain_type) :: strain
!     strain type
!     Type wich containt strain information (related to reference structure)

   integer, allocatable :: typat(:)
!     typat(ntypat) 
!     Type of each atom (primitive cell)

   real(dp), allocatable :: qph1l(:,:) 
!     qph1l(3,nph1l)
!     List of nph1l wavevectors

   real(dp), allocatable :: znucl(:)                
!     znucl(ntype) 
!     Gives the nuclear number for each type of atom (primitive cell)

   real(dp), allocatable :: amu(:)                  
!     amu(ntypat)
!     Mass of the atoms array (atomic mass unit) (primitive cell)

   real(dp), allocatable :: zeff(:,:,:)             
!     zeff(3,3,natom) Effective charges

   real(dp), allocatable :: internal_strain(:,:,:)    
!    internal_strain(6,natom,3)
!    internal strain tensor 

   real(dp), allocatable :: forces(:,:)    
!    internal_strain(6,natom,3)
!    internal strain tensor 

   real(dp), allocatable :: xcart(:,:)
!    xcart(3, natom) 
!    positions of atoms in cartesian coordinates

   real(dp), allocatable :: dynmat(:,:,:,:,:,:)
!    dynmat(2,3,natom,3,natom,nph1l) 
!    dynamical matrix for each q points

   real(dp), allocatable :: phfrq(:,:)
!    phfrq(3*natom,nph1l)
!    array with all phonons frequencies for each q points in Hartree/cm

   type(ifc_type) :: ifcs
!   type with ifcs constants (short + ewald) 
!   also contains the number of cell and the indexes

   type(ifc_type),dimension(:),allocatable :: phonon_strain_coupling
!   Array of ifc with phonon_strain_coupling coupling for each strain 

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
!! natom  = number of atoms in primitive cell
!! ntypat = number of type of atoms
!! nph1l  = number of wavevectors in phonon
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

subroutine effective_potential_init(natom,nph1l,nrpt,ntypat, eff_pot,name)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nph1l,nrpt,ntypat
 character(len=fnlen), optional,intent(in) :: name
 type(effective_potential_type), intent(out) :: eff_pot
!arrays
!Local variables-------------------------------
!scalar
 integer :: ii
!arrays
 character(len=500) :: msg

! *************************************************************************

! free the effective potential in input
  call effective_potential_free(eff_pot)

  eff_pot%natom     = natom  ! Set the number of atoms
  eff_pot%ntypat    = ntypat ! Set the number of type of atoms 
  eff_pot%ifcs%nrpt = nrpt   ! Set the number of real space points 
                             ! used to integrate IFC 
  eff_pot%nph1l     = nph1l  ! Set the list of nph1l wavevector for phonon

  if (present(name)) then
    eff_pot%name = name
  else
    eff_pot%name = ''
  end if

! Check the number of atoms
 if (natom < 1) then
   write(msg, '(a,a,a,i10,a)' )&
&   'The cell must have at least one atom.',ch10,&
&   'The number of atom is  ',natom,'.'
   MSG_BUG(msg)
 end if

 if (ntypat < 1) then
   write(msg, '(a,a,a,i10,a)' )&
&   'The cell must have at least one type of atom.',ch10,&
&   'The number of type of atom is  ',ntypat,'.'
   MSG_BUG(msg)
 end if

 eff_pot%energy = zero
 eff_pot%ucvol  = zero
 eff_pot%epsilon_inf = zero
 eff_pot%elastic_constants = zero
 eff_pot%rprimd = zero
 eff_pot%acell = zero
 eff_pot%strain%name = ""
 eff_pot%strain%delta = zero
 eff_pot%strain%direction = zero
 eff_pot%strain%strain = zero
 eff_pot%internal_stress = zero
 eff_pot%external_stress = zero
 eff_pot%has_3rd = .false.

!Allocation of mass of the atoms array (atomic mass unit)
 ABI_ALLOCATE(eff_pot%amu,(ntypat))
 eff_pot%amu = zero 

!Allocation of Effective charges array 
 ABI_ALLOCATE(eff_pot%typat,(natom))
 eff_pot%typat = zero 

!Allocation of internal strain tensor 
 ABI_ALLOCATE(eff_pot%znucl,(ntypat))
 eff_pot%znucl = zero

!Allocation of Effective charges array 
 ABI_ALLOCATE(eff_pot%zeff,(3,3,natom))
 eff_pot%zeff = zero 

!Allocation of total ifc
 ABI_ALLOCATE(eff_pot%ifcs%atmfrc,(2,3,natom,3,natom,nrpt))
 eff_pot%ifcs%atmfrc = zero 

!Allocation of ewald part of ifc
 ABI_ALLOCATE(eff_pot%ifcs%ewald_atmfrc,(2,3,natom,3,natom,nrpt))
 eff_pot%ifcs%ewald_atmfrc = zero 

!Allocation of short range part of ifc
 ABI_ALLOCATE(eff_pot%ifcs%short_atmfrc,(2,3,natom,3,natom,nrpt))
 eff_pot%ifcs%short_atmfrc = zero 

!Allocation of cell of ifc
 ABI_ALLOCATE(eff_pot%ifcs%cell,(nrpt,3))
 eff_pot%ifcs%short_atmfrc = zero 
 
!Allocation of phonon strain coupling array (3rd order)
 ABI_DATATYPE_ALLOCATE(eff_pot%phonon_strain_coupling,(12))
 do ii = 1,12
   ABI_ALLOCATE(eff_pot%phonon_strain_coupling(ii)%atmfrc,(2,3,natom,3,natom,nrpt))
   ABI_ALLOCATE(eff_pot%phonon_strain_coupling(ii)%cell,(nrpt,3))
   eff_pot%phonon_strain_coupling(ii)%nrpt   = nrpt
   eff_pot%phonon_strain_coupling(ii)%atmfrc = zero
   eff_pot%phonon_strain_coupling(ii)%cell   = zero
 end do

!Allocation of internal strain tensor 
 ABI_ALLOCATE(eff_pot%internal_strain,(6,natom,3))
 eff_pot%internal_strain = zero

!Allocation of forces
 ABI_ALLOCATE(eff_pot%forces,(3,natom))
 eff_pot%forces = zero


!Allocation of list of nph1l wavevectors
 ABI_ALLOCATE(eff_pot%qph1l,(3,nph1l))
 eff_pot%qph1l  = zero

!Allocation of dynamical matrix 
 ABI_ALLOCATE(eff_pot%dynmat,(2,3,natom,3,natom,nph1l))
 eff_pot%dynmat = zero

!Allocation of frequecies arrays
 ABI_ALLOCATE(eff_pot%phfrq,(3*natom,nph1l))
 eff_pot%phfrq = zero

!Allocation of internal strain tensor 
 ABI_ALLOCATE(eff_pot%xcart,(3,natom))
 eff_pot%xcart = zero

!Initialisation of the supercell
 call init_supercell(eff_pot%natom, 0, real((/0,0,0/),dp), eff_pot%rprimd, eff_pot%typat,&
&                    eff_pot%xcart, eff_pot%supercell)

end subroutine effective_potential_init 
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
!!
!! OUTPUT
!! eff_pot = supercell structure with data to be output
!!
!! PARENTS
!!   epigene
!!
!! CHILDREN
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
  integer :: ii
!array

! *************************************************************************

  eff_pot%name   = ''
  eff_pot%natom  = zero
  eff_pot%energy = zero
  eff_pot%ntypat = zero
  eff_pot%nph1l  = zero
  eff_pot%ucvol  = zero
  eff_pot%rprimd = zero
  eff_pot%acell  = zero
  eff_pot%ifcs%nrpt   = zero
  eff_pot%epsilon_inf = zero
  eff_pot%elastic_constants = zero
  eff_pot%strain%name = ""
  eff_pot%strain%delta = zero
  eff_pot%strain%direction = zero
  eff_pot%strain%strain = zero
  eff_pot%internal_stress = zero
  eff_pot%external_stress = zero
  eff_pot%has_3rd = .false.

  if(allocated(eff_pot%typat))then
    eff_pot%typat=zero
    ABI_DEALLOCATE(eff_pot%typat)
  end if

  if(allocated(eff_pot%qph1l))then
    eff_pot%qph1l=zero
    ABI_DEALLOCATE(eff_pot%qph1l)
  end if

  if(allocated(eff_pot%znucl))then
    eff_pot%znucl=zero
    ABI_DEALLOCATE(eff_pot%znucl)
  end if

  if(allocated(eff_pot%amu))then
    eff_pot%amu=zero
    ABI_DEALLOCATE(eff_pot%amu)
  end if

  if(allocated(eff_pot%zeff))then
    eff_pot%zeff=zero
    ABI_DEALLOCATE(eff_pot%zeff)
  end if

  if(allocated(eff_pot%internal_strain)) then
    eff_pot%internal_strain=zero
    ABI_DEALLOCATE(eff_pot%internal_strain)
  end if

  if(allocated(eff_pot%forces)) then
    eff_pot%forces=zero
    ABI_DEALLOCATE(eff_pot%forces)
  end if

  if(allocated(eff_pot%xcart))then
    eff_pot%xcart=zero
    ABI_DEALLOCATE(eff_pot%xcart)
  end if

  if(allocated(eff_pot%dynmat))then
    eff_pot%dynmat=zero
    ABI_DEALLOCATE(eff_pot%dynmat)
  end if

  if(allocated(eff_pot%phfrq))then
    eff_pot%phfrq=zero
    ABI_DEALLOCATE(eff_pot%phfrq)
  end if

  if(allocated(eff_pot%ifcs%short_atmfrc))then
    eff_pot%ifcs%short_atmfrc=zero
    ABI_DEALLOCATE(eff_pot%ifcs%short_atmfrc)
  end if

  if(allocated(eff_pot%phonon_strain_coupling))then
    do ii = 1,12
      eff_pot%phonon_strain_coupling(ii)%atmfrc = zero
      eff_pot%phonon_strain_coupling(ii)%cell   = zero
      ABI_DEALLOCATE(eff_pot%phonon_strain_coupling(ii)%atmfrc)
      ABI_DEALLOCATE(eff_pot%phonon_strain_coupling(ii)%cell)
    end do
    ABI_DATATYPE_DEALLOCATE(eff_pot%phonon_strain_coupling)
  end if

  if(allocated(eff_pot%ifcs%ewald_atmfrc))then
    eff_pot%ifcs%ewald_atmfrc=zero
    ABI_DEALLOCATE(eff_pot%ifcs%ewald_atmfrc)
  end if

  if(allocated(eff_pot%ifcs%atmfrc))then
    eff_pot%ifcs%atmfrc=zero
    ABI_DEALLOCATE(eff_pot%ifcs%atmfrc)
  end if

  call destroy_supercell(eff_pot%supercell)
  
end subroutine effective_potential_free
!!***

!****f* m_effective_potential/effective_potential_generateSupercell
!!
!! NAME
!! effective_potential_generateSupercell
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
!!   mover_effpot
!!
!! CHILDREN
!!   init_supercell
!!
!! SOURCE
 
subroutine effective_potential_generateSupercell(eff_pot,n_cell,option,asr,comm)

 use m_phonon_supercell
 use m_ewald

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_generateSupercell'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,asr
 integer,intent(inout) :: comm
!array
 integer,intent(in) :: n_cell(3)
 type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
 integer,parameter :: master=0
 integer :: first_coordinate
 integer :: ia,i1,i2,i3,irpt,irpt2,irpt_ref,min1,min2,min3
 integer :: min1_cell,min2_cell,min3_cell,max1_cell,max2_cell,max3_cell
 integer :: max1,max2,max3,mu,my_rank,nu,nproc,second_coordinate,sumg0,nrpt
 real(dp) :: sum,ucvol
 character(len=500) :: message
 logical :: short_range = .false.
 logical :: iam_master=.FALSE.
!array
 real(dp) :: acell(3)
 real(dp) :: gmet(3,3),rmet(3,3)
 real(dp) :: gprimd(3,3)
 real(dp),allocatable :: dyew(:,:,:,:,:), dyewq0(:,:,:)
 real(dp),allocatable :: xred(:,:),xred_tmp(:,:),zeff_tmp(:,:,:)

 type(supercell_type) :: super_cell
 type(ifc_type) :: ifc_tmp
! *************************************************************************

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Check the size of the cell
 do ia=1,3
   if(n_cell(ia)<0.or.n_cell(ia)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'n_cell(',ia,') is ',n_cell(ia),', which is lower than 0 of superior than 20.',&
&     ch10,'Action: correct n_cell(',ia,').'
     MSG_ERROR(message)
   end if
 end do

 call init_supercell(eff_pot%natom, 0,&
&  real(n_cell,dp),&
&   eff_pot%rprimd,&
&   eff_pot%typat,&
&   eff_pot%xcart,&
&   super_cell)

!1-Store the information of the supercell of the reference structure into effective potential
 eff_pot%supercell%natom  = eff_pot%natom
 eff_pot%supercell%natom_supercell  = super_cell%natom_supercell
 eff_pot%supercell%typat_supercell  = super_cell%typat_supercell
 eff_pot%supercell%qphon  = super_cell%qphon
 eff_pot%supercell%rprimd_supercell = super_cell%rprimd_supercell
 eff_pot%supercell%xcart_supercell  = super_cell%xcart_supercell
 eff_pot%supercell%uc_indexing_supercell = super_cell%uc_indexing_supercell
 eff_pot%supercell%atom_indexing_supercell = super_cell%atom_indexing_supercell

!2 Check if the bound of new cell correspond to the effective potential 
!only for option=zero
!set min and max
 min1 = minval(eff_pot%ifcs%cell(:,1)) ; max1 = maxval(eff_pot%ifcs%cell(:,1))
 min2 = minval(eff_pot%ifcs%cell(:,2)) ; max2 = maxval(eff_pot%ifcs%cell(:,2))
 min3 = minval(eff_pot%ifcs%cell(:,3)) ; max3 = maxval(eff_pot%ifcs%cell(:,3))

 if(option==0) then
   if(((max1-min1+1)/=n_cell(1).and.&
&    (max2-min2+1)/=n_cell(2).and.(max3-min3+1)/=n_cell(3))) then
     write(message, '(5a,3I3,5a,3I3,2a)' )&
&      'ifcsupercell is set to zero, the longe range interation might be wrong',ch10,&
&      'because it is not recompute.',ch10,&
&      'The previous harmonic part is build for ',(max1-min1+1),(max2-min2+1),(max3-min3+1)&
&,     ' cell.',ch10,'Be sure than the dipole-dipole interation ',ch10,&
&      'is correct for the supercell: ',n_cell(:),ch10,&
&      'or set ifcsupercell to 1'
     MSG_WARNING(message)
   else
     write(message,'(a)')&
&    'ifcsupercell is set to zero, the longe range interation is not recompute'
     MSG_WARNING(message)
   end if

!3-Adapt harmonic part   
 else if (option>=1.and.all(n_cell(:)>0)) then

   write(message,'(a,(80a),3a)') ch10,('=',i1=1,80),ch10,' Generation of new ifc',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   irpt_ref = 0
   irpt     = 0

   call find_bound(min1_cell,max1_cell,n_cell(1))
   call find_bound(min2_cell,max2_cell,n_cell(2))
   call find_bound(min3_cell,max3_cell,n_cell(3))
   write(message, '(2a)' )&
&        ' ifcsupercell is set to one or two, the dipole-dipole interation is recompute.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  Generate new bound
   if(option==1)then
     if ((abs(min1) > abs(min1_cell)).or.(abs(max1) > abs(max1_cell)).or.&
&        (abs(min2) > abs(min2_cell)).or.(abs(max2) > abs(max2_cell)).or.&
&        (abs(min3) > abs(min3_cell)).or.(abs(max3) > abs(max3_cell))) then
       write(message, '(3a,3I3,a)' )&
&        ' The previous harmonic part was build for bigger cell',ch10,&
&         'ifc is adjust on ',int((/(max1-min1+1),(max2-min2+1),(max3-min3+1)/),dp),' cell'
       MSG_WARNING(message)
!      If the cell is smaller, we redifine new cell to take into acount all atoms
       call init_supercell(eff_pot%natom, 0,real((/(max1-min1+1),(max2-min2+1),(max3-min3+1)/),dp),&
&          eff_pot%rprimd,eff_pot%typat,eff_pot%xcart,super_cell)
     else
       min1 = min1_cell ; min2 = min2_cell ; min3 = min3_cell
       max1 = max1_cell ; max2 = max2_cell ; max3 = max3_cell
     end if
   else  if(option==2)then
     if ((abs(min1) > abs(min1_cell)).or.(abs(max1) > abs(max1_cell)).or.&
&        (abs(min2) > abs(min2_cell)).or.(abs(max2) > abs(max2_cell)).or.&
&        (abs(min3) > abs(min3_cell)).or.(abs(max3) > abs(max3_cell))) then
       write(message, '(5a)' )&
&        'ifcsupercell is set to two or three, the dipole-dipole interation is recompute.',ch10,&
&        'The previous harmonic part was build for bigger cell,',ch10,&
&        'The accuracy might be reduced.'
       MSG_WARNING(message)
     end if
     min1 = min1_cell ; min2 = min2_cell ; min3 = min3_cell
     max1 = max1_cell ; max2 = max2_cell ; max3 = max3_cell
   else if (option==3)then
     if ((abs(min1) > abs(min1_cell)).or.(abs(max1) > abs(max1_cell)).or.&
&        (abs(min2) > abs(min2_cell)).or.(abs(max2) > abs(max2_cell)).or.&
&        (abs(min3) > abs(min3_cell)).or.(abs(max3) > abs(max3_cell))) then
       write(message, '(5a)' )&
&        'ifcsupercell is set to two or three, the dipole-dipole interation is recompute.',ch10,&
&        'The previous harmonic part was build for bigger cell,',ch10,&
&        'The accuracy might be reduced.'
       MSG_WARNING(message)
     end if
   end if   

!  Print the new boundary
   write(message,'(5a,2I3,a,2I3,a,2I3,4a)') ch10,' New bound for ifc (long+short):',&
&    ch10,ch10, " x=[",min1,max1,"], y=[",min2,max2,"] and z=[",min3,max3,"]",ch10,ch10,&
&    " Computation of new dipole-dipole interaction."
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')


!TEST_AM_NEW
if(.false.)then
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


!  Initialisation of ifc temporary
   ABI_MALLOC(ifc_tmp%cell,(ifc_tmp%nrpt,3))
   ABI_MALLOC(ifc_tmp%short_atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ABI_MALLOC(ifc_tmp%ewald_atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ABI_MALLOC(ifc_tmp%atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ifc_tmp%atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%short_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%ewald_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%cell(:,:) = zero

!  Allocate and initialize some array
   !ABI_ALLOCATE(xred_tmp,(3,2*eff_pot%natom))
   ABI_ALLOCATE(xred,(3,super_cell%natom_supercell))
   ABI_ALLOCATE(zeff_tmp,(3,3,super_cell%natom_supercell))
   ABI_ALLOCATE(dyew,(2,3,super_cell%natom_supercell,3,super_cell%natom_supercell))
   ABI_ALLOCATE(dyewq0,(3,3,super_cell%natom_supercell))
 
   dyew            = zero
   dyewq0          = zero
   xred(:,:)       = zero
   zeff_tmp(:,:,:) = zero
   sumg0           = zero
   acell           = one


   call matr3inv(super_cell%rprimd_supercell,gprimd)
   call xcart2xred(super_cell%natom_supercell,super_cell%rprimd_supercell,&
&                    super_cell%xcart_supercell,xred)
   call metric(gmet,gprimd,-1,rmet,super_cell%rprimd_supercell,ucvol)

!  Fill fake zeff array for ewald9
   do irpt=1,ifc_tmp%nrpt
     first_coordinate  = ((irpt-1)*eff_pot%natom) + 1
     second_coordinate = first_coordinate + eff_pot%natom-1
     zeff_tmp(:,:,first_coordinate:second_coordinate) = eff_pot%zeff
   end do
   
   write(std_out,*)"enter in big ewald9"
   call ewald9(acell,eff_pot%epsilon_inf,dyew,&
&              gmet,gprimd,super_cell%natom_supercell,real((/0,0,0/),dp),rmet,&
&              super_cell%rprimd_supercell,sumg0,ucvol,xred,&
&              zeff_tmp)
   write(std_out,*)"enter in q0dy3_calc"
   call q0dy3_calc(super_cell%natom_supercell,dyewq0,dyew,2)
   write(std_out,*)"enter in q0dy3_apply"
   call q0dy3_apply(super_cell%natom_supercell,dyewq0,dyew)
   write(std_out,*)"done"

   first_coordinate  = ((irpt_ref-1)*eff_pot%natom) + 1

   irpt=0
   do i1=min1,max1
     do i2=min2,max2
       do i3=min3,max3
         irpt = irpt+1
!        Fill the index of the cell
         ifc_tmp%cell(irpt,1)=i1; ifc_tmp%cell(irpt,2)=i2; ifc_tmp%cell(irpt,3)=i3
!        Fill the short range part (calculated previously)
         second_coordinate = ((irpt-1)*eff_pot%natom) + 1
         ifc_tmp%ewald_atmfrc(:,:,:,:,:,irpt) =&
&          dyew(:,:,first_coordinate:first_coordinate+eff_pot%natom-1,&
&                 :,second_coordinate:second_coordinate+eff_pot%natom-1) + tol10
           do irpt2=1,eff_pot%ifcs%nrpt
             if(eff_pot%ifcs%cell(irpt2,1)==i1.and.&
&               eff_pot%ifcs%cell(irpt2,2)==i2.and.&
&               eff_pot%ifcs%cell(irpt2,3)==i3.and.&
&               any(eff_pot%ifcs%short_atmfrc(:,:,:,:,:,irpt2) > tol9)) then
               ifc_tmp%short_atmfrc(:,:,:,:,:,irpt) = eff_pot%ifcs%short_atmfrc(:,:,:,:,:,irpt2)
             end if
           end do
         end do
       end do
     end do


   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(zeff_tmp)
   ABI_DEALLOCATE(dyew)
   ABI_DEALLOCATE(dyewq0)

else if(.false.)then   
  
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

!  Initialisation of ifc temporary
   ABI_MALLOC(ifc_tmp%cell,(ifc_tmp%nrpt,3))
   ABI_MALLOC(ifc_tmp%short_atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ABI_MALLOC(ifc_tmp%ewald_atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ABI_MALLOC(ifc_tmp%atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ifc_tmp%atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%short_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%ewald_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%cell(:,:) = zero

!  Allocate and initialize some array
   ABI_ALLOCATE(xred_tmp,(3,2*eff_pot%natom))
   ABI_ALLOCATE(xred,(3,super_cell%natom_supercell))
   ABI_ALLOCATE(zeff_tmp,(3,3,2*eff_pot%natom))
   ABI_ALLOCATE(dyew,(2,3,2*eff_pot%natom,3,2*eff_pot%natom))
 
   dyew            = zero
   xred(:,:)       = zero
   xred_tmp(:,:)   = zero
   zeff_tmp(:,:,:) = zero
   sumg0           = zero
   acell           = one


   call matr3inv(super_cell%rprimd_supercell,gprimd)
   call xcart2xred(super_cell%natom_supercell,super_cell%rprimd_supercell,&
&                    super_cell%xcart_supercell,xred)
   call metric(gmet,gprimd,-1,rmet,super_cell%rprimd_supercell,ucvol)

!  Fill the atom position of the first cell (reference cell)
   first_coordinate  = ((irpt_ref-1)*eff_pot%natom) + 1
   second_coordinate = first_coordinate + eff_pot%natom-1
   xred_tmp(:,1:eff_pot%natom) = xred(:,first_coordinate:second_coordinate)
!  Fill fake zeff array for ewald9
   zeff_tmp = zero
   zeff_tmp(:,:,1:eff_pot%natom) = eff_pot%zeff
   zeff_tmp(:,:,eff_pot%natom+1:2*eff_pot%natom) = eff_pot%zeff
   irpt=0
   do i1=min1,max1
     do i2=min2,max2
       do i3=min3,max3
         irpt = irpt+1
!        Fill the index of the cell
         ifc_tmp%cell(irpt,1)=i1; ifc_tmp%cell(irpt,2)=i2; ifc_tmp%cell(irpt,3)=i3

         if(option == 1 .or. option==2) then
!          Compute new dipole-dipole interaction
           dyew = zero
           if (i1==0.and.i2==0.and.i3==0) then
             call ewald9(acell,eff_pot%epsilon_inf,dyew(:,:,1:eff_pot%natom,:,1:eff_pot%natom),&
&                        gmet,gprimd,eff_pot%natom,real((/0,0,0/),dp),rmet,&
&                        super_cell%rprimd_supercell,sumg0,ucvol,xred_tmp(:,1:eff_pot%natom),&
&                        eff_pot%zeff)
             ifc_tmp%ewald_atmfrc(:,:,:,:,:,irpt) = dyew(:,:,1:eff_pot%natom,:,1:eff_pot%natom) + tol10
            else
              first_coordinate  = ((irpt-1)*eff_pot%natom) + 1
              second_coordinate = first_coordinate + eff_pot%natom-1
              xred_tmp(:,eff_pot%natom+1:2*eff_pot%natom)=&
&                 xred(:,first_coordinate:second_coordinate)
              call ewald9(acell,eff_pot%epsilon_inf,dyew,gmet,gprimd,&
&                      int(2*eff_pot%natom),real((/0,0,0/),dp),&
&                      rmet,super_cell%rprimd_supercell,&
&                      sumg0,ucvol,xred_tmp,zeff_tmp)
              ifc_tmp%ewald_atmfrc(:,:,:,:,:,irpt) = &
&              dyew(:,:,1:eff_pot%natom,:,eff_pot%natom+1:2*eff_pot%natom) + tol10
            end if
         end if
!        Fill the short range part (calculated previously)
           do irpt2=1,eff_pot%ifcs%nrpt
             if(eff_pot%ifcs%cell(irpt2,1)==i1.and.&
&               eff_pot%ifcs%cell(irpt2,2)==i2.and.&
&               eff_pot%ifcs%cell(irpt2,3)==i3.and.&
&               any(eff_pot%ifcs%short_atmfrc(:,:,:,:,:,irpt2) > tol9)) then
               ifc_tmp%short_atmfrc(:,:,:,:,:,irpt) = eff_pot%ifcs%short_atmfrc(:,:,:,:,:,irpt2)
               if (option==3)then
                 ifc_tmp%ewald_atmfrc(:,:,:,:,:,irpt) = eff_pot%ifcs%ewald_atmfrc(:,:,:,:,:,irpt2)
               end if
             end if
           end do
         end do
       end do
     end do
   ABI_DEALLOCATE(xred_tmp)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(zeff_tmp)
   ABI_DEALLOCATE(dyew)

!TEST_AM_OLD 
 else

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

!  Initialisation of ifc temporary
   ABI_MALLOC(ifc_tmp%cell,(ifc_tmp%nrpt,3))
   ABI_MALLOC(ifc_tmp%short_atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ABI_MALLOC(ifc_tmp%ewald_atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ABI_MALLOC(ifc_tmp%atmfrc,(2,3,eff_pot%natom,3,eff_pot%natom,ifc_tmp%nrpt))
   ifc_tmp%atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%short_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%ewald_atmfrc(:,:,:,:,:,:) = zero
   ifc_tmp%cell(:,:) = zero

!  Allocate and initialize some array
   ABI_ALLOCATE(xred_tmp,(3,2*eff_pot%natom))
   ABI_ALLOCATE(xred,(3,super_cell%natom_supercell))
   ABI_ALLOCATE(zeff_tmp,(3,3,2*eff_pot%natom))
   ABI_ALLOCATE(dyew,(2,3,2*eff_pot%natom,3,2*eff_pot%natom))
 
   dyew            = zero
   xred(:,:)       = zero
   xred_tmp(:,:)   = zero
   zeff_tmp(:,:,:) = zero
   sumg0           = zero
   acell           = one


   call matr3inv(super_cell%rprimd_supercell,gprimd)
   call xcart2xred(super_cell%natom_supercell,super_cell%rprimd_supercell,&
&                    super_cell%xcart_supercell,xred)
   call metric(gmet,gprimd,-1,rmet,super_cell%rprimd_supercell,ucvol)

!  Fill the atom position of the first cell (reference cell)
   first_coordinate  = ((irpt_ref-1)*eff_pot%natom) + 1
   second_coordinate = first_coordinate + eff_pot%natom-1
   xred_tmp(:,1:eff_pot%natom) = xred(:,first_coordinate:second_coordinate)
!  Fill fake zeff array for ewald9
   zeff_tmp(:,:,1:eff_pot%natom) = eff_pot%zeff
   zeff_tmp(:,:,eff_pot%natom+1:2*eff_pot%natom) = eff_pot%zeff

   irpt=0
   do i1=min1,max1
     do i2=min2,max2
       do i3=min3,max3
         irpt = irpt+1
!        Fill the index of the cell
         ifc_tmp%cell(irpt,1)=i1; ifc_tmp%cell(irpt,2)=i2; ifc_tmp%cell(irpt,3)=i3

!        Compute new dipole-dipole interaction
         dyew = zero
         if (i1==0.and.i2==0.and.i3==0) then
           call ewald9(acell,eff_pot%epsilon_inf,dyew(:,:,1:eff_pot%natom,:,1:eff_pot%natom),&
&                      gmet,gprimd,eff_pot%natom,real((/0,0,0/),dp),rmet,&
&                      super_cell%rprimd_supercell,sumg0,ucvol,xred_tmp(:,1:eff_pot%natom),&
&                      eff_pot%zeff)
           ifc_tmp%ewald_atmfrc(:,:,:,:,:,irpt) = dyew(:,:,1:eff_pot%natom,:,1:eff_pot%natom) + tol10
         else
           first_coordinate  = ((irpt-1)*eff_pot%natom) + 1
           second_coordinate = first_coordinate + eff_pot%natom-1
           xred_tmp(:,eff_pot%natom+1:2*eff_pot%natom)=&
&              xred(:,first_coordinate:second_coordinate)
           call ewald9(acell,eff_pot%epsilon_inf,dyew,gmet,gprimd,&
&                    int(2*eff_pot%natom),real((/0,0,0/),dp),&
&                    rmet,super_cell%rprimd_supercell,&
&                    sumg0,ucvol,xred_tmp,zeff_tmp)
           ifc_tmp%ewald_atmfrc(:,:,:,:,:,irpt) = &
&            dyew(:,:,1:eff_pot%natom,:,eff_pot%natom+1:2*eff_pot%natom) + tol10
         end if

!        Fill the short range part (calculated previously)
         do irpt2=1,eff_pot%ifcs%nrpt
           if(eff_pot%ifcs%cell(irpt2,1)==i1.and.&
&             eff_pot%ifcs%cell(irpt2,2)==i2.and.&
&             eff_pot%ifcs%cell(irpt2,3)==i3.and.&
&             any(eff_pot%ifcs%short_atmfrc(:,:,:,:,:,irpt2) > tol9)) then
             ifc_tmp%short_atmfrc(:,:,:,:,:,irpt) = eff_pot%ifcs%short_atmfrc(:,:,:,:,:,irpt2)
           end if
         end do
       end do
     end do
   end do
   ABI_DEALLOCATE(xred_tmp)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(zeff_tmp)
   ABI_DEALLOCATE(dyew)
 end if

!  Compute total ifc
   ifc_tmp%atmfrc = ifc_tmp%short_atmfrc + ifc_tmp%ewald_atmfrc

!  set new ifc into effective potential
   eff_pot%ifcs = ifc_tmp
 
   ABI_FREE(ifc_tmp%cell)
   ABI_FREE(ifc_tmp%short_atmfrc)
   ABI_FREE(ifc_tmp%ewald_atmfrc)
   ABI_FREE(ifc_tmp%atmfrc)
   
 end if! end if option   

   if(asr >= 0) then
 !Impose sum rule
   call effective_potential_applySumRule(asr,eff_pot%ifcs,eff_pot%natom,1)
   call effective_potential_applySumRule(asr,eff_pot%ifcs,eff_pot%natom,2)
   call effective_potential_applySumRule(asr,eff_pot%ifcs,eff_pot%natom)
 end if

end subroutine effective_potential_generateSupercell
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
!!   effective_potential_generateSupercell
!!
!! CHILDREN
!!   
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

! Found the cell of reference
 do irpt = 1,ifc%nrpt   
   if(ifc%cell(irpt,1)==0.and.&
&     ifc%cell(irpt,2)==0.and.&
&     ifc%cell(irpt,3)==0) then
     irpt_ref = irpt
     cycle
   end if
 end do

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
!!   epigene
!!
!! CHILDREN
!!   
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
&                                         eff_pot%supercell%xcart_supercell,1,displacement=disp)
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
&             dynmat(1,nu,ib,mu,ia),dynmat(2,nu,ib,mu,ia)
       end do
     end do
   end do
 end do



! Deallocation of arrays
 ABI_DEALLOCATE(fred)
 ABI_DEALLOCATE(fcart)


 end subroutine effective_potential_effpot2dynmat
!!***


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
!!   epigene
!!
!! CHILDREN
!!   
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
  character(len=500) :: message
  integer :: ii,irpt,jj,msym
  real(dp):: ucvol

! type(anaddb_dataset_type) :: inp
!array
  real(dp) :: gmet(3,3),rmet(3,3)
  real(dp) :: gprimd(3,3),rprimd(3,3)
  real(dp) :: acell(3),gprim(3,3),rprim(3,3)
  real(dp),allocatable :: xred(:,:)
  character :: title(eff_pot%ntypat) 
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
    call metric(gmet,gprimd,-1,rmet,eff_pot%rprimd,ucvol)

!   Convert to rprim (dimensionless)
    do ii=1,3
      do jj=1,3
        rprim(ii,jj)=eff_pot%rprimd(ii,jj)/eff_pot%acell(jj)
      end do
    end do
    
!   Obtain reciprocal space primitive transl g from inverse trans of r
!   (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
    call matr3inv(rprim,gprim)

!   transfert basic values
    ddb%natom  = eff_pot%natom 
    ddb%mpert  = ddb%natom+6
    ddb%msize  = 3*ddb%mpert*3*ddb%mpert;
    ddb%ntypat = eff_pot%ntypat
    ddb%occopt = 3 ! default value
    ddb%prtvol = 0 ! default value
    ddb%rprim  = rprim ! dimensioless real space primitive vectors
    ddb%gprim  = gprim ! dimensioless reciprocal space primitive vectors
    ddb%acell  = eff_pot%acell 
        msym   = 1
!  Setup crystal type
    ABI_ALLOCATE(xred,(3,ddb%natom))
!    call xcar2xred(ddb%natom,eff_pot%rprimd,eff_pot%xcart,xred)
!Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
    ABI_MALLOC(symrel,(3,3,msym))
    ABI_MALLOC(symafm,(msym))
    ABI_MALLOC(tnons,(3,msym))

    call crystal_init(crystal,1,ddb%natom,size(eff_pot%znucl),eff_pot%ntypat,1,eff_pot%rprimd,&
&       eff_pot%typat,xred,eff_pot%znucl,eff_pot%znucl,0,.FALSE.,.FALSE.,title)!,&
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
        rprim(ii,jj)=eff_pot%supercell%rprimd_supercell(ii,jj)/eff_pot%acell(jj)
      end do
    end do
    
!   Obtain reciprocal space primitive transl g from inverse trans of r
!   (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
    call matr3inv(rprim,gprim)

!   transfert basic values
    ddb%natom  = eff_pot%supercell%natom_supercell 
    ddb%ntypat = eff_pot%ntypat
    ddb%mpert  = ddb%natom+6
    ddb%msize  = 3*ddb%mpert*3*ddb%mpert;
    ddb%occopt = 3 ! default value
    ddb%prtvol = 0 ! default value
    ddb%rprim  = rprim ! dimensioless real space primitive vectors
    ddb%gprim  = gprim ! dimensioless reciprocal space primitive vectors
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
!!   epigene
!!
!! CHILDREN
!!   
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
 integer :: irpt,lenstr
 real(dp) :: tcpui,twalli
 character(len=strlen) :: string
!array
 type(crystal_t) :: Crystal
 type(anaddb_dataset_type) :: inp
 type(ddb_type) :: ddb
 real(dp),allocatable  :: d2asr(:,:,:,:,:),singular(:),uinvers(:,:),vtinvers(:,:)

! *************************************************************************

 if (option > 0) then

!  First: transfer into ddb structure:
   call effective_potential_effpot2ddb(ddb,Crystal,eff_pot,n_cell,nph1l,option,qph1l)
   
!  Setup fake anaddb_dataset
   string = ''
   lenstr = 0
   call invars9(inp,lenstr,ddb%natom,string)
!  fill it with epigene_dataset values
   inp%prt_ifc = 1
   inp%ifcflag = 1
   inp%qph1l   = qph1l
   inp%nph1l   = nph1l

   ABI_ALLOCATE(d2asr,(2,3,ddb%natom,3,ddb%natom))
 ! Pre allocate array used if asr in [3,4]
!  ABI_ALLOCATE(singular,(1:3*ddb%natom*(3*ddb%natom-1)/2))
!  ABI_CALLOC(uinvers,(1:3*ddb%natom*(3*ddb%natom-1)/2,1:3*ddb%natom*(3*ddb%natom-1)/2))
!  ABI_CALLOC(vtinvers,(1:3*ddb%natom*(3*ddb%natom-1)/2,1:3*ddb%natom*(3*ddb%natom-1)/2))

  call mkphbs(eff_pot%ifcs,Crystal,inp,ddb,d2asr,filename,&
&  singular,tcpui,twalli,uinvers,vtinvers,eff_pot%zeff)

   ABI_DEALLOCATE(d2asr)
!  ABI_DEALLOCATE(singular)
!  ABI_DEALLOCATE(vtinvers)
!  ABI_DEALLOCATE(uinvers)

 end if

 end subroutine effective_potential_printPDOS
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
!!   epigene
!!
!! CHILDREN
!!   wrtout
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

    write(message,'(a,F20.10,a,a,I3,a,a,I4,a,a,I4,a,a,I3,a,a)') &
&     '  - Reference energy:  ',eff_pot%energy ,ch10,&
&     '  - Number of types of atoms:  ',eff_pot%ntypat ,ch10,&
&     '  - Number of atoms:  ',eff_pot%natom ,ch10,&
&     '  - Number of cells:  ',eff_pot%ifcs%nrpt ,ch10,&
&     '  - Number of qpoints:  ',eff_pot%nph1l ,ch10,&
&     '  - Primitive vectors:  '
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    write(std_out,'(3(F12.6))') (eff_pot%rprimd) 
    write(ab_out,'(3(F12.6))')  (eff_pot%rprimd)    
    write(message,'(a,a,a)') '  - acell:  '
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    write(std_out,'(3(F12.6))') (eff_pot%acell)
    write(ab_out,'(3(F12.6))')  (eff_pot%acell)
    write(message,'(a,a,a)') '  - Dielectric tensor:  '
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    write(std_out,'(3(F12.6))') (eff_pot%epsilon_inf)
    write(ab_out,'(3(F12.6))')  (eff_pot%epsilon_inf)
    write(message,'(a,a,a)') '  - Elastic tensor:  '
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    write(std_out,'(6(F12.6))') (eff_pot%elastic_constants)
    write(ab_out,'(6(F12.6))')  (eff_pot%elastic_constants)

    do ia=1,eff_pot%natom
      write(message,'(a,I4)') '  - Atoms',ia
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
      write(std_out,'(a,3(F10.4))') "    - atomic number:",eff_pot%znucl(eff_pot%typat(ia))
      write(ab_out,'(a,3(F10.4))')  "    - atomic number:",eff_pot%znucl(eff_pot%typat(ia))
      write(std_out,'(a,3(F10.4))') "    - atomic mass:",eff_pot%amu(eff_pot%typat(ia))
      write(ab_out,'(a,3(F10.4))')  "    - atomic mass:",eff_pot%amu(eff_pot%typat(ia))
      write(std_out,'(a,3(F12.6))') "    - position:",eff_pot%xcart(:,ia)
      write(ab_out,'(a,3(F12.6))')  "    - position:",eff_pot%xcart(:,ia)
      write(std_out,'(a)') "    - Effective charges:"
      write(ab_out,'(a)')  "    - Effective charges:"
      do ii = 1,3
        write(std_out,'(a,3(F12.6))') "  ",eff_pot%zeff(:,ii,ia)
        write(ab_out,'(a,3(F12.6))')  "  ",eff_pot%zeff(:,ii,ia)
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
!! Print the supercell of the effetive_potential
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!   mover_effpot
!!
!! CHILDREN
!!   wrtout
!!
!! SOURCE
 
subroutine effective_potential_printSupercell(eff_pot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_printSupercell'
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
  integer :: iatom
  character(len=500) :: message
!array
  real(dp), allocatable :: xred(:,:)
! *************************************************************************

  ABI_ALLOCATE(xred,(3,eff_pot%supercell%natom_supercell))
!**********************************************************************
! Write basics values 
!**********************************************************************

 write (ab_out, '(2a)') '  Lattice vectors for supercell :',ch10
 write (ab_out, '(a,I7)') ' natom ', eff_pot%supercell%natom_supercell
 write (ab_out, *)
 write (ab_out, '(a)') ' znucl '
 do iatom = 1, size(eff_pot%znucl)
   write (ab_out, '(I5)', ADVANCE="NO") int(eff_pot%znucl(iatom))
   if (mod(iatom,6) == 0) write (ab_out, *)
 end do
 write (ab_out, *)
 write (ab_out, *)
 write (ab_out, '(a,I7)') ' ntypat', size(eff_pot%znucl)
 write (ab_out, '(a)') ' typat '
 do iatom = 1, eff_pot%supercell%natom_supercell
   write (ab_out, '(I5)', ADVANCE="NO")&
&    eff_pot%supercell%typat_supercell(&
&    eff_pot%supercell%atom_indexing_supercell(iatom))
   if (mod(iatom,12) == 0) write (ab_out, *)
 end do
 write (ab_out, *)
 write (ab_out, '(a)') ' acell 1.0 1.0 1.0'
 write (ab_out, '(a)') ' rprim'
 write (ab_out, '(3E23.14)') eff_pot%supercell%rprimd_supercell(:,1)
 write (ab_out, '(3E23.14)') eff_pot%supercell%rprimd_supercell(:,2)
 write (ab_out, '(3E23.14)') eff_pot%supercell%rprimd_supercell(:,3)
 write (ab_out, *)
 write (ab_out, '(a)') ' xcart'
 do iatom = 1, eff_pot%supercell%natom_supercell
   write (ab_out, '(3E23.14)') eff_pot%supercell%xcart_supercell(:,iatom)
 end do
 call xcart2xred(eff_pot%supercell%natom_supercell,eff_pot%supercell%rprimd_supercell,&
&                eff_pot%supercell%xcart_supercell,xred)
 write (ab_out, '(a)') ' xred'
 do iatom = 1, eff_pot%supercell%natom_supercell
   write (ab_out, '(3E23.14)') xred(:,iatom)
 end do
 
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
!! Copyright (C) 2000-2015 ABINIT group (AM)
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
!!      epigene
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_writeXML(eff_pot,option,filename)

 use defs_basis
 use m_errors
 use m_profiling_abi  
 use m_epigene_dataset, only : epigene_dataset_type

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
 integer :: iph1l,irpt,mu,nu
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
 if (option == 1) then

!  convert natom in character
   write (natom,'(I9)') eff_pot%natom
   if(present(filename)) then
     namefile=filename
   else
     namefile='ref.xml'
   end if

   call isfile(namefile,'new')

 if (open_file(namefile,message,unit=unit_xml,form="formatted",status="new",action="write") /= 0) then
   MSG_ERROR(message)
 end if 

 write(message, '(a,(80a),a)' ) ch10,&
&  ('-',ii=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

  write(message,'(a,a,a)')ch10,&
 &   ' Generation of the xml file for the reference structure in ',namefile

  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

!Write header
  WRITE(unit_xml,'("<?xml version=""1.0"" ?>")')
  WRITE(unit_xml,'("<System_definition>")')
  
  WRITE(unit_xml,'("  <energy>")')
  WRITE(unit_xml,'(E23.14)') (eff_pot%energy)
  WRITE(unit_xml,'("  </energy>")')
  
  WRITE(unit_xml,'("  <unit_cell units=""bohrradius"">")')
  WRITE(unit_xml,'(3(E23.14))') (eff_pot%rprimd)
  WRITE(unit_xml,'("  </unit_cell>")')
  
  WRITE(unit_xml,'("  <epsilon_inf units=""epsilon0"">")')
  WRITE(unit_xml,'(3(E23.14))') (eff_pot%epsilon_inf)
  WRITE(unit_xml,'("  </epsilon_inf>")')
  
  WRITE(unit_xml,'("  <elastic units=""hartree"">")')
  WRITE(unit_xml,'(6(E23.14))') (eff_pot%elastic_constants)
  WRITE(unit_xml,'("  </elastic>")')
  
  do ia=1,eff_pot%natom
    WRITE(unit_xml,'("  <atom mass=""",1F10.5,""" massunits=""atomicmassunit"">")') &
&     eff_pot%amu(eff_pot%typat(ia))
    WRITE(unit_xml,'("    <position units=""bohrradius"">")') 
    WRITE(unit_xml,'(3(E23.14))') (eff_pot%xcart(:,ia))
    WRITE(unit_xml,'("    </position>")')
    WRITE(unit_xml,'("    <borncharge units=""abs(e)"">")')
    WRITE(unit_xml,'(3(E23.14))') (eff_pot%zeff(:,:,ia))
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
  do irpt=1,eff_pot%ifcs%nrpt
    if(any(abs(eff_pot%ifcs%short_atmfrc(1,:,:,:,:,irpt))>tol9)) then 
      WRITE(unit_xml,'("  <local_force_constant units=""hartree/bohrradius**2"">")')
      WRITE(unit_xml,'("    <data>")')
      do ia=1,eff_pot%natom
        do mu=1,3
          do ib=1,eff_pot%natom
            do  nu=1,3
              WRITE(unit_xml,'(e22.14)', advance="no")(eff_pot%ifcs%short_atmfrc(1,mu,ia,nu,ib,irpt))
            end do
          end do
          WRITE(unit_xml,'(a)')''
        end do
      end do
      WRITE(unit_xml,'("    </data>")')
      WRITE(unit_xml,'("    <cell>")')
      WRITE(unit_xml,'(3(I4))') (eff_pot%ifcs%cell(irpt,:))
      WRITE(unit_xml,'("    </cell>")')
      WRITE(unit_xml,'("  </local_force_constant>")')
    end if
!TEST_AM
    if(.false.)then
      if(any(abs(eff_pot%ifcs%ewald_atmfrc(1,:,:,:,:,irpt))>tol9)) then 
        WRITE(unit_xml,'("  <ewald_force_constant units=""hartree/bohrradius**2"">")')
        WRITE(unit_xml,'("    <data>")')
        do ia=1,eff_pot%natom
          do mu=1,3
            do ib=1,eff_pot%natom
              do  nu=1,3
                WRITE(unit_xml,'(e22.14)', advance="no")(eff_pot%ifcs%ewald_atmfrc(1,mu,ia,nu,ib,irpt))
              end do
            end do
            WRITE(unit_xml,'(a)')''
          end do
        end do
        WRITE(unit_xml,'("    </data>")')
        WRITE(unit_xml,'("    <cell>")')
        WRITE(unit_xml,'(3(I4))') (eff_pot%ifcs%cell(irpt,:))
        WRITE(unit_xml,'("    </cell>")')
        WRITE(unit_xml,'("  </ewald_force_constant>")')
      end if
    end if
!TEST_AM
!  Print the IFC total for each cell, data is array 3*natom*3*natom
!  [ [x1 x2 ....]
!    [y1 y2 ....] for atom 1
!    [z1 z2 ....]
!    [x1 x2 ....]
!    [y1 y2 ....] for atom 2
!    [z1 z2 ....]
!    ....       ] 
    if(all(abs(eff_pot%ifcs%atmfrc(1,:,:,:,:,irpt))<tol9)) then 
      if(any(abs(eff_pot%ifcs%short_atmfrc(1,:,:,:,:,irpt))>tol9)) then
        write(message, '(a,a,a,a)' )&
&        ' There is no total range but short range in your effective potential',ch10,&
&        'Action: contact abinit group'
        MSG_BUG(message)
      end if
    else
      WRITE(unit_xml,'("  <total_force_constant units=""hartree/bohrradius**2"">")')
      WRITE(unit_xml,'("    <data>")')
      do ia=1,eff_pot%natom
        do mu=1,3
          do ib=1,eff_pot%natom
            do nu=1,3
              WRITE(unit_xml,'(e22.14)', advance="no")(eff_pot%ifcs%atmfrc(1,mu,ia,nu,ib,irpt))
            end do
          end do
          WRITE(unit_xml,'(a)')''
        end do
      end do
      WRITE(unit_xml,'("    </data>")')
      WRITE(unit_xml,'("    <cell>")')
      WRITE(unit_xml,'(3(I4))') (eff_pot%ifcs%cell(irpt,:))
      WRITE(unit_xml,'("    </cell>")')
      WRITE(unit_xml,'("    </total_force_constant>")')
    end if
  end do

  do iph1l=1,eff_pot%nph1l
    WRITE(unit_xml,'("  <phonon>")')
    WRITE(unit_xml,'("    <qpoint units=""2pi*G0"">")')
    WRITE(unit_xml,'(3(E23.14))') (eff_pot%qph1l(:,iph1l))
    WRITE(unit_xml,'("    </qpoint>")')
    WRITE(unit_xml,'("    <frequencies units=""reciprocal cm"">")')
    WRITE(unit_xml,'(3(e22.14))') (eff_pot%phfrq(:,iph1l))
    WRITE(unit_xml,'("    </frequencies>")')
    WRITE(unit_xml,'("    <dynamical_matrix units=""hartree/bohrradius**2"">")')
    do ia=1,eff_pot%natom
      do mu=1,3
        do ib = 1,eff_pot%natom
          do nu=1,3
            WRITE(unit_xml,'(e22.14)',advance='no') (eff_pot%dynmat(1,nu,ib,mu,ia,iph1l))
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
    do ia=1,eff_pot%natom
      do mu=1,3
        WRITE(unit_xml,'(e22.14)', advance="no")&
&            (eff_pot%internal_strain(ii,ia,mu))
      end do
      WRITE(unit_xml,'(a)')''
    end do
    WRITE(unit_xml,'("    </correction_force>")')
    if (eff_pot%has_3rd) then
      do irpt=1,eff_pot%phonon_strain_coupling(ii)%nrpt
        WRITE(unit_xml,'("    <correction_force_constant units=""hartree/bohrradius**2"">")')
        WRITE(unit_xml,'("      <data>")')
        do ia=1,eff_pot%natom
          do mu=1,3
            do ib=1,eff_pot%natom
              do  nu=1,3
                WRITE(unit_xml,'(e22.14)', advance="no")&
&                    (eff_pot%phonon_strain_coupling(ii)%atmfrc(1,mu,ia,nu,ib,irpt))
              end do
            end do
            WRITE(unit_xml,'(a)')''
          end do
        end do
        WRITE(unit_xml,'("      </data>")')
        WRITE(unit_xml,'("      <cell>")')
        WRITE(unit_xml,'(3(I4))') (eff_pot%phonon_strain_coupling(ii)%cell(irpt,:))
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
  
 end if
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
!! Copyright (C) 2000-2015 ABINIT group (AM)
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
!!      epigene
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_writeNETCDF(eff_pot,option,filename)

 use defs_basis
 use m_errors
 use m_profiling_abi  
 use m_epigene_dataset, only : epigene_dataset_type
 use m_nctk
#if defined HAVE_TRIO_NETCDF
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
 integer :: amu_id,bec_id,ifcscell_id,ifccell_id,epsinf_id,elastic_id
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

#if defined HAVE_TRIO_NETCDF

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
   ncerr = nf90_def_dim(ncid,"natom",eff_pot%natom,natom_id)
   NCF_CHECK_MSG(ncerr," define dimension natom")

   ncerr = nf90_def_dim(ncid,"ntypat",eff_pot%ntypat,ntypat_id)
   NCF_CHECK_MSG(ncerr," define dimension ntypat")

   ncerr = nf90_def_dim(ncid,"nrpt",eff_pot%ifcs%nrpt,nrpt_id)
   NCF_CHECK_MSG(ncerr," define dimension ntypat")

   ncerr = nf90_def_var(ncid, "typat", NF90_INT, natom_id, typat_id)
   NCF_CHECK_MSG(ncerr," define variable typat")

   npsp = size(eff_pot%znucl)
   if (npsp /= eff_pot%ntypat) then
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

   ncerr = nf90_put_var(ncid,rprimd_id, eff_pot%rprimd)
   NCF_CHECK_MSG(ncerr," write variable rprimd")

   ncerr = nf90_put_var(ncid,epsinf_id, eff_pot%epsilon_inf)
   NCF_CHECK_MSG(ncerr," write variable epsilon_inf")

   ncerr = nf90_put_var(ncid,elastic_id , eff_pot%elastic_constants)
   NCF_CHECK_MSG(ncerr," write variable elastic_constant")

   ncerr = nf90_put_var(ncid, bec_id, eff_pot%zeff)
   NCF_CHECK_MSG(ncerr," write variable bec")

   ncerr = nf90_put_var(ncid,xcart_id, eff_pot%xcart)
   NCF_CHECK_MSG(ncerr," write variable xcart")

   ncerr = nf90_put_var(ncid,ifccell_id, eff_pot%ifcs%cell)
   NCF_CHECK_MSG(ncerr," write variable cell")

   ncerr = nf90_put_var(ncid,ifcs_id, eff_pot%ifcs%short_atmfrc)
   NCF_CHECK_MSG(ncerr," write variable short ifc")

   ncerr = nf90_put_var(ncid,ifc_id, eff_pot%ifcs%atmfrc)
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
!! Copyright (C) 2000-2015 ABINIT group (AM)
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
!!      epigene
!!
!! CHILDREN
!!
!! SOURCE

subroutine effective_potential_writeAbiInput(eff_pot,filename,strain)

 use defs_basis
 use m_errors
 use m_profiling_abi  
 use m_epigene_dataset, only : epigene_dataset_type
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
 real(dp) :: xred(3,eff_pot%natom)
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
  write(unit,'(a)') itoa(eff_pot%natom)
  write(unit,'(" rfelfd3 =  3")')
  write(unit,'(" rfstrs3 =  3")')
  write(unit,'("  rfdir3 =  1 1 1")')
  write(unit,'(" tolvrs3 =  1.0d-8")')
  write(unit,'("")')
                                                                         
  write(unit,'("#STRUCTURE")')
  write(unit,'(" natom = ")',advance='no')
  write(unit,'(a)') itoa(eff_pot%natom)
  write(unit,'(" znucl =")',advance='no')
  write(unit,'(10(F4.0))') (eff_pot%znucl)
  write(unit,'("ntypat = ")',advance='no')
  write(unit,'(a)') itoa(eff_pot%ntypat)
  write(unit,'(" typat = ")',advance='no')
  write(unit,'(10(I2))') (eff_pot%typat)
  write(unit,'(" acell = 1 1 1")')
  write(unit,'(" rprim  ")')
  write(unit,'(3(F20.10))') (matmul(eff_pot%rprimd,strain%strain))
  write(unit,'("  xred  ")')
  call xcart2xred(eff_pot%natom,eff_pot%rprimd,eff_pot%xcart,xred)
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
!! Evaluate the gradient of  of effective potential
!! to calculates the forces in reduced coordinates
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!   mover
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine effective_potential_getForces(eff_pot,fcart,fred,natom,rprimd,xcart,comm,&
&                   displacement)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_getForces'
 use interfaces_41_geometry
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: comm
  integer, intent(in) :: natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
  real(dp),intent(in),optional :: displacement(3,natom)
  real(dp),intent(out) :: fcart(3,natom),fred(3,natom)
!Local variables-------------------------------
!scalar
  integer,parameter :: master=0
  real(dp):: temp(3)
  integer :: i1,i2,i3,ia,ii,jj,ib,ll,my_rank,nproc
  character(len=500) :: message
  logical :: iam_master=.FALSE.
!array
  real(dp):: disp_tmp1(3,natom)
  integer :: cell_number(3)
  integer :: cell_atom1(3),cell_atom2(3)
  character(500) :: msg

! *************************************************************************

!MPI variables
  nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
  iam_master = (my_rank == master)

  fred = zero
  fcart = zero

  if (natom /= eff_pot%supercell%natom_supercell) then
    write(msg,'(a,I7,a,I7,a)')' The number of atoms is not correct :',natom,&
&   ' in argument istead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(msg)
  end if
  
  cell_number(:) = int(eff_pot%supercell%qphon(:))

  if (any(cell_number <= 0)) then
    write(msg,'(a,a)')' No supercell found for getForces'
    MSG_ERROR(msg)
  end if

  disp_tmp1(:,:) = zero
  if (present(displacement)) then
    disp_tmp1(:,:) = displacement(:,:)
  else
    do ii = 1, natom
      disp_tmp1(:,ii) = xcart(:,ii) - eff_pot%supercell%xcart_supercell(:,ii)
    end do
  end if

  ii = 1
  do i1 = 1,cell_number(1)
    do i2 = 1,cell_number(2)
      do i3 = 1,cell_number(3)
        do ia = 1, eff_pot%natom
          temp = zero
          do jj = 1,eff_pot%ifcs%nrpt
!           get the cell of atom2  (0 0 0, 0 0 1...)
            cell_atom2(1) =  (i1-1) + eff_pot%ifcs%cell(jj,1)
            call index_periodic(cell_atom2(1),cell_number(1))
            cell_atom2(2) =  (i2-1) + eff_pot%ifcs%cell(jj,2)
            call index_periodic(cell_atom2(2),cell_number(2))
            cell_atom2(3) =  (i3-1) + eff_pot%ifcs%cell(jj,3)
            call index_periodic(cell_atom2(3),cell_number(3))
            do ib = 1, eff_pot%natom
!              index of the second atom in the displacement array
               ll = cell_atom2(1)*cell_number(2)*cell_number(3)*eff_pot%natom+&
&                   cell_atom2(2)*cell_number(3)*eff_pot%natom+&
&                   cell_atom2(3)*eff_pot%natom+&
&                   ib

               temp(:) = temp(:) + matmul(eff_pot%ifcs%atmfrc(1,:,ia,:,ib,jj), disp_tmp1(:,ll))

             end do
           end do
          
          fcart(:,ii) = -1 * temp(:)

          ii = ii + 1
        end do
      end do
    end do
  end do
 
  call fcart2fred(fcart,fred,rprimd,natom)

end subroutine effective_potential_getForces
!!***

!****f* m_effective_potential/effective_potential_getEnergy
!!
!! NAME
!! effective_potential_getEnergy
!!
!! FUNCTION
!! Evaluate the effective_potential_getEnergy of effective potential
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!   mover
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine effective_potential_getEnergy(eff_pot,energy,natom,rprimd,xcart,comm,&
&                            displacement,strain1,strain2,external_stress)

  use m_strain

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effective_potential_getEnergy'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: comm
  integer, intent(in) :: natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
  real(dp),intent(in),optional :: strain1(6),strain2(6)
  real(dp),intent(in),optional :: displacement(3,eff_pot%supercell%natom_supercell)
  real(dp),intent(in),optional :: external_stress(6)
  real(dp),intent(out) :: energy
!Local variables-------------------------------
!scalar
  integer :: ii,ncell
  real(dp):: elastic_part,ifc_part
  character(len=500) :: message
  integer :: nproc,my_rank
  logical :: iam_master
  integer, parameter:: master = 0
!array
  type(strain_type) :: strain
  integer  :: supercell(3)
  real(dp) :: strain_tmp1(6),strain_tmp2(6)
  real(dp) :: external_stress_tmp(6)
  real(dp) :: disp_tmp1(3,eff_pot%supercell%natom_supercell)
! *************************************************************************

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Check some variables
  if (natom /= eff_pot%supercell%natom_supercell) then
    write(message,'(a,I7,a,I7,a)')' The number of atoms is not correct :',natom,&
&   ' in argument istead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(message)
  end if

  if (present(displacement).and.(size(displacement(1,:)) /= eff_pot%supercell%natom_supercell)) then
    write(message,'(a,I7,a,I7,a)')' The number of atoms is not correct :',size(displacement(1,:)),&
&   ' in displacement array instead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(message)
  end if

  do ii=1,3
    if(eff_pot%supercell%qphon(ii)<0.or.eff_pot%supercell%qphon(ii)>10)then
      write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'eff_pot%supercell%qphon(',ii,') is ',eff_pot%supercell%qphon(ii),&
&     ', which is lower than 0 of superior than 10.',ch10,'Action: correct n_cell(',ii,').'
      MSG_ERROR(message)
    end if
  end do

  if(iam_master) then
    write(message, '(a,a,a)' ) ch10,' enter get_energy : Calculation of the energy'
    call wrtout(std_out,message,'COLL')
    write(message, '(a,a,a)' ) ch10,' Calculation of the energy with effective potential'
    call wrtout(ab_out,message,'COLL')
  end if

  supercell(1:3) = eff_pot%supercell%qphon(1:3)
  ncell          = product(eff_pot%supercell%qphon)
!------------------------------------
! 1 - Transfert the reference energ
!------------------------------------
  energy = eff_pot%energy * ncell

  write(message, '(a,a,1ES24.16,a)' ) ch10,' Energy of the reference strucure :',energy,' Hartree'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

!------------------------------------
! 2 - Computation of the harmonic part (IFC) of the energy :
!------------------------------------

  disp_tmp1(:,:) = zero
  if (present(displacement)) then
    disp_tmp1(:,:) = displacement(:,:)
  else
    !try to compute the displacement   
    do ii = 1, eff_pot%supercell%natom_supercell
      disp_tmp1(:,ii) = xcart(:,ii) - eff_pot%supercell%xcart_supercell(:,ii)
    end do
  end if

  ifc_part = harmonic_energy(eff_pot,disp_tmp1)

  write(message, '(a,1ES24.16,a)' ) ' Energy of the ifc part :',ifc_part,' Hartree'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  energy = energy + ifc_part

!------------------------------------
! 3 - Computation of the elastic part of the energy :
!------------------------------------
  strain_tmp1(:) = zero
  strain_tmp2(:) = zero
  ! Try to find the strain from argument
  if (present(strain1).and.(present(strain2))) then
    strain_tmp1(:) = strain1(:)
    strain_tmp2(:) = strain2(:)
  else if (present(strain1).and.(.not.present(strain2))) then
    strain_tmp1(:) = strain1(:)
    strain_tmp2(:) = strain1(:)
  else
    ! else => calculation of the strain 
    call strain_get(rprimd,eff_pot%supercell%rprimd_supercell,strain)
    do ii=1,3
      strain_tmp1(ii) = strain%strain(ii,ii)
    end do
      strain_tmp1(4) = strain%strain(2,3)
      strain_tmp1(5) = strain%strain(3,1)
      strain_tmp1(6) = strain%strain(2,1)

      strain_tmp2(:) = strain_tmp1(:)
  end if

  if (present(external_stress)) then
    external_stress_tmp(:) = external_stress(:)
  else
    external_stress_tmp(:) = zero
  end if

  elastic_part = elastic_energy(eff_pot,ncell,strain_tmp1,strain_tmp2,&
&                               external_stress=external_stress_tmp)

  write(message, '(a,1ES24.16,a)' ) ' Energy of the elastic part :',elastic_part,' Hartree'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  energy = energy + elastic_part

! 4 - Print the total energy
  write(message, '(a,a,1ES24.16,a)' ) ch10,' Total energy :',energy,' Hartree'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(message, '(a,a,a)' ) ch10,' end get_energy.F90',ch10
  call wrtout(std_out,message,'COLL')


end subroutine effective_potential_getEnergy
!!***

!!****f* m_effective_potential/elastic_energy
!! NAME
!!  elastic_energy
!!
!! FUNCTION
!! Compute the energy related to the application of strain
!!
!! INPUTS 
!! eff_pot = effective potential structure
!! ncell   = number of cell
!! strain1(6) =  first strain to apply
!! strain2(6) =  second strain to apply
!! external_strees(6) =  external stress to apply
!! 
!! OUTPUT
!! energy 
!!
!! SOURCE
function elastic_energy(eff_pot,ncell,strain1,strain2,external_stress) result(energy)

!Arguments ------------------------------------
! scalar

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elastic_energy'
!End of the abilint section

  real(dp) energy
  integer, intent(in) :: ncell
! array
  type(effective_potential_type), intent(in) :: eff_pot
  real(dp),intent(in) :: strain1(6),strain2(6)
  real(dp),optional,intent(in) :: external_stress(6)
!Local variables-------------------------------
! scalar
! array
! *************************************************************************

  energy = zero

  energy = half*dot_product(matmul(strain1,eff_pot%elastic_constants),strain2)

  energy = energy + dot_product(eff_pot%internal_stress, strain1)

  if(present(external_stress)) then
    energy = energy - dot_product(external_stress, strain1)
  end if

  energy = energy * ncell

end function elastic_energy
!!***

!!****f* m_effective_potential/harmonic_energy
!! NAME
!!  harmonic_energy
!!
!! FUNCTION
!!  This fonction compute the harmonic part of the energy
!!  of the supercell in the eff_pot
!! INPUTS 
!!  eff_pot = effective potential of the structure 
!!            also contain supercell information
!!  disp    = diplacement vector (3, cell1 (atm1 atm2 ...) cell2 (atm1 atm2 ...)...)
!!             
!! OUTPUT
!!
!! SOURCE

function harmonic_energy(eff_pot,disp) result(energy)

!Arguments ------------------------------------
! scalar

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'harmonic_energy'
!End of the abilint section

  real(dp) energy
! array
  type(effective_potential_type), intent(in) :: eff_pot
  real(dp),intent(in) :: disp(3,eff_pot%supercell%natom_supercell)
!Local variables-------------------------------
! scalar
  integer :: i1,i2,i3,ia,ib,ii,jj,kk,ll
! array
  real(dp) :: temp(3)
  integer :: cell_number(3)
  integer :: cell_atom1(3),cell_atom2(3)
  character(500) :: msg
! *************************************************************************

  if ((size(disp(1,:)) /= eff_pot%supercell%natom_supercell)) then
    write(msg,'(a,I7,a,I7,a)')' The number of atoms is not correct :',size(disp(1,:)),&
&   ' in displacement array instead of ',eff_pot%supercell%natom_supercell, ' in supercell'
    MSG_ERROR(msg)
  end if

  cell_number(:) = int(eff_pot%supercell%qphon(:))

  if (any(cell_number <= 0)) then
    write(msg,'(a,a)')' No supercell found for getEnergy'
    MSG_ERROR(msg)
  end if

  energy = zero

  ii = 1
  do i1 = 1,cell_number(1)
    do i2 = 1,cell_number(2)
      do i3 = 1,cell_number(3)
!       get the cell of atom 1 (0 0 0, 0 0 1...)
        cell_atom1 = eff_pot%supercell%uc_indexing_supercell(:,ii)
        call index_periodic(cell_atom1(1),cell_number(1))
        call index_periodic(cell_atom1(2),cell_number(2))
        call index_periodic(cell_atom1(3),cell_number(3))

        do ia = 1, eff_pot%natom
!         index of the first atom in the displacement array
          kk = cell_atom1(1)*cell_number(2)*cell_number(3)*eff_pot%natom+&
&              cell_atom1(2)*cell_number(3)*eff_pot%natom+&
&              cell_atom1(3)*eff_pot%natom+&
&              ia

          temp = zero
          do jj = 1,eff_pot%ifcs%nrpt
!           get the cell of atom2  (0 0 0, 0 0 1...)
            cell_atom2(1) =  (i1-1) + eff_pot%ifcs%cell(jj,1)
            call index_periodic(cell_atom2(1),cell_number(1))
            cell_atom2(2) =  (i2-1) + eff_pot%ifcs%cell(jj,2)
            call index_periodic(cell_atom2(2),cell_number(2))
            cell_atom2(3) =  (i3-1) + eff_pot%ifcs%cell(jj,3)
            call index_periodic(cell_atom2(3),cell_number(3))
            do ib = 1, eff_pot%natom
!             index of the second atom in the displacement array              
               ll = cell_atom2(1)*cell_number(2)*cell_number(3)*eff_pot%natom+&
&                   cell_atom2(2)*cell_number(3)*eff_pot%natom+&
&                   cell_atom2(3)*eff_pot%natom+&
&                   ib

               temp = temp + matmul(eff_pot%ifcs%atmfrc(1,:,ia,:,ib,jj),disp(:,ll))

            end do
          end do

          energy = energy + half * dot_product(temp,disp(:,kk))         
!          energy = energy + half * dot_product(eff_pot%forces(:,ia),disp(:,kk))

          ii = ii + 1
        end do
      end do
    end do
  end do

end function harmonic_energy
!!***

!****f* m_effective_potential/effective_potential_getDeltaEnergy
!!
!! NAME
!! effective_potential_getDeltaEnergy
!!
!! FUNCTION
!! Evaluate the energy due to 1 atomic displacement
!!
!! INPUTS
!! eff_pot = effective potential structure
!!
!! OUTPUT
!!
!!
!! PARENTS
!!   mover
!!
!! CHILDREN
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
  real(dp) :: ifc_part
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

  ifc_part =  0

  energy = energy + ifc_part


end subroutine effective_potential_getDeltaEnergy
!!***

!!****f* m_effective_potential/index_periodic
!! NAME
!!  harmonic_energy
!!
!! FUNCTION
!!
!! INPUTS 
!!
!! OUTPUT
!!
!! SOURCE

recursive subroutine index_periodic(index,n_cell)


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

  if (index < 0) then
    index = index + n_cell
    call index_periodic(index,n_cell)
  else
    if(index > n_cell-1) then
      index = index - n_cell
      call index_periodic(index,n_cell)
    end if
  end if

! *********************************************************************
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
    max=(n_cell)/2; min=-max; if(mod(n_cell,2)==0) min= min + 1
  else
    min=-(n_cell)/2; max=-min; if(mod(n_cell,2)==0) max= max - 1
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
  if(e1%natom==e2%natom.and.&
&     e1%ifcs%nrpt==e2%ifcs%nrpt.and.&
&     e1%ntypat==e2%ntypat.and.&
&     e1%nph1l==e2%nph1l.and.&
&     e1%energy==e2%energy.and.&
&     e1%ucvol==e2%ucvol) then
    res = .true.
  end if
end function effective_potential_compare
!!***

end module m_effective_potential
