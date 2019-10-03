!!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_effective_potential
!!
!! NAME
!! m_effective_potential
!!
!! FUNCTION
!! Module for the effective potential
!! Container type is defined, and destruction, print subroutines
!! Contain also routine to evaluate the energy,forces and stresses
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

module m_effective_potential

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_strain
 use m_ifc
 use m_supercell
 use m_phonons
 use m_ddb
 use m_polynomial_conf
 use m_polynomial_coeff
 use m_anharmonics_terms
 use m_harmonics_terms
 use m_xmpi
 use m_ewald
 use m_nctk
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,       only : replace, ftoa, itoa
 use m_io_tools,       only : open_file, get_unit
 use m_dtfil,          only : isfile
 use m_symtk,          only : matr3inv
 use m_effpot_mpi,     only : effpot_mpi_init,effpot_mpi_type,effpot_mpi_free
 use m_abihist,        only : abihist
 use m_geometry,       only : fred2fcart,fcart2fred, xcart2xred, xred2xcart, metric
 use m_crystal,        only : crystal_t, crystal_init
 !use m_anaddb_dataset, only : anaddb_dataset_type, anaddb_dtset_free, outvars_anaddb, invars9

 implicit none

 public :: effective_potential_distributeResidualForces
 public :: effective_potential_evaluate
 public :: effective_potential_free
 public :: effective_potential_freeCoeffs
 public :: effective_potential_freempi
 public :: effective_potential_generateDipDip
 public :: effective_potential_getDisp
 public :: effective_potential_init
 public :: effective_potential_initmpi
 public :: effective_potential_print
 public :: effective_potential_printSupercell
 public :: effective_potential_setCoeffs
 public :: effective_potential_setConfinement
 public :: effective_potential_setElastic3rd
 public :: effective_potential_setElastic4th
 public :: effective_potential_setElasticDispCoupling
 public :: effective_potential_setStrainPhononCoupling
 public :: effective_potential_setSupercell
 public :: effective_potential_writeAbiInput
 public :: effective_potential_writeXML
 !AM_EXPERIMENTAL
 public :: effective_potential_computeGradient
! public :: effective_potential_effpot2ddb
 ! public :: effective_potential_printPDOS
 public :: effective_potential_checkDEV
 public :: effective_potential_writeNETCDF
 public :: OPERATOR(==)
 !AM_EXPERIMENTAL
!!***

!!****t* m_effective_potential/effective_potential_type
!! NAME
!! effective_potential_type
!!
!! FUNCTION
!! datatype for a effective potential constructed.
!!
!! SOURCE

 type, public :: effective_potential_type

   character(len=fnlen) :: name
!     Name of the molecule (CaTiO3,...)

   type(crystal_t) :: crystal
!    crystal type
!    contains all information of the crystal

   real(dp):: energy
!     Energy of the system (Hatree)

   real(dp), allocatable :: fcart(:,:)
!    forces(3,natom)
!    initial cartesian forces of the system

   real(dp) :: strten(6)
!    strten(6)
!    initial stresses (Ha/bohr^3) of the system

   type(harmonics_terms_type) :: harmonics_terms
!     type with all information for harmonics terms

   type(anharmonics_terms_type) :: anharmonics_terms
!     type with all information for anharmonics terms

   type(polynomial_conf_type) :: confinement
!     type with all the informations for the confinement

   type(supercell_type) :: supercell
!     super cell type
!     Store all the information of the suppercell

   logical :: has_anharmonicsTerms
!     True : the aharmonic part is present

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over the supercell
   type(effpot_mpi_type) :: mpi_ifc
!  effpot_mpi_type with all the information for the IFC paralellisation

   type(effpot_mpi_type) :: mpi_coeff
!  effpot_mpi_type with all the information for the polynomial coefficients paralellisation

 end type effective_potential_type
!!***

interface operator (==)
  module procedure effective_potential_compare
end interface

CONTAINS  !===========================================================================================

!!****f* m_effective_potential/effective_potential_init
!!
!! NAME
!! effective_potential_init
!!
!! FUNCTION
!! Initialize effective_potential datatype
!!
!! INPUTS
!! crytal<type(crystal_t)> = datatype with all the information for the crystal
!! energy = energy of the reference structure
!! ifcs <type(ifc_type)> = ifc type with cell,ewald short and total range of the ifcs
!! ncoeff = number of coefficients in the polynomial
!! nqpt = number of qpoints
!! comm  = mpi comunicator
!! coeffs(ncoeff)<type(polynomial_coeff_type)> = optional,list of coefficients for the anharmonic part
!! dynmat(2,3,natom,3,natom,nqpt) = optional,dynamical matrix for each qpoints
!! epsilon_inf(3,3) = optional,dielectric tensor
!! elastic_constants(6,6) = optional,elastic constants tensor
!! elastic3rd(6,6,6) = optional,3 order derivatives with respect to to 3 strain
!! elastic_displacement(6,6,3,natom)=optional, 3 order derivatives with respect to 2 strain and
!!                                             1 Atom disp
!! fcart(3,natom) = optional,optional,initial fcart in the structure
!! strain_coupling(6,natom,3) = optional, internal strain coupling parameters
!! strten(6) = optional,optional,initial strain in the structure
!! name = optional, name of the structure
!! phonon_strain(6) = optional,ifc type for the phonon-strain coupling (should be in anharmonics_terms)
!! phfrq(3*natom,nqpt) = optional,phonons frequencies for each q points in Hartree/cm
!! qpoints(3,nqpt) = optional,list of qpoints wavevectors
!! has_anharmonicsTerms = optional, (default false) flag to set the anharmonics terms
!! supercell<type(supercell_type)> = optional, supercell type to define
!! zeff(3,natom) = optional,effective charges
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = effective_potential datatype to be initialized
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_init(crystal,eff_pot,energy,ifcs,ncoeff,nqpt,comm,&
&                                   coeffs,dynmat,elastic_constants,elastic3rd,&
&                                   elastic_displacement,epsilon_inf,&
&                                   fcart,strain_coupling,strten,name,phonon_strain,&
&                                   polynomial_conf,phfrq,qpoints,has_anharmonicsTerms,&
&                                   supercell,zeff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,ncoeff,nqpt
 real(dp),intent(in):: energy
 logical,optional,intent(in) :: has_anharmonicsTerms
!arrays
 type(crystal_t),intent(in) :: crystal
 type(effective_potential_type), intent(out) :: eff_pot
 type(ifc_type),intent(in) :: ifcs
 character(len=fnlen), optional,intent(in) :: name
 real(dp),optional,intent(in) :: epsilon_inf(3,3),elastic_constants(6,6)
 real(dp),optional,intent(in) :: dynmat(:,:,:,:,:,:),qpoints(:,:),phfrq(:,:)
 real(dp),optional,intent(in) :: strain_coupling(:,:,:),zeff(:,:,:)
 type(supercell_type),optional,intent(in) :: supercell
 type(ifc_type),optional,intent(in) :: phonon_strain(6)
 real(dp),optional,intent(in) :: elastic3rd(6,6,6)
 real(dp),optional,intent(in) :: elastic_displacement(6,6,3,crystal%natom)
 type(polynomial_coeff_type),optional :: coeffs(ncoeff)
 real(dp),optional,intent(in) :: strten(6)
 real(dp),optional,intent(in) :: fcart(3,crystal%natom)
 type(polynomial_conf_type),optional,intent(in) :: polynomial_conf

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
!Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
 call crystal_init(crystal%amu,eff_pot%crystal,crystal%space_group,crystal%natom,&
&                  crystal%npsp,crystal%ntypat,crystal%nsym,crystal%rprimd,&
&                  crystal%typat,crystal%xred,crystal%zion,crystal%znucl,&
&                  crystal%timrev,.FALSE.,.FALSE.,crystal%title,&
&                  symrel=crystal%symrel,tnons=crystal%tnons,symafm=crystal%symafm)

!4-Fill harmonic part
 call harmonics_terms_init(eff_pot%harmonics_terms,ifcs,crystal%natom,ifcs%nrpt)

!5-Init the anharmonics_terms to set the flag to false
 call anharmonics_terms_init(eff_pot%anharmonics_terms,crystal%natom,ncoeff)

!5-Fill optional inputs
 ABI_ALLOCATE(eff_pot%fcart,(3,eff_pot%crystal%natom))
 eff_pot%fcart = zero
 if(present(fcart))then
   eff_pot%fcart = fcart
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

 if(present(strain_coupling))then
   call harmonics_terms_setInternalStrain(eff_pot%harmonics_terms,crystal%natom,strain_coupling)
 end if

 if(present(zeff))then
   call harmonics_terms_setEffectiveCharges(eff_pot%harmonics_terms,crystal%natom,zeff)
 end if

 eff_pot%strten = zero
 if(present(strten))then
   eff_pot%strten(:) = strten(:)
 end if

!Set the flag for the strain coupling
 if(present(has_anharmonicsTerms)) then
   eff_pot%has_anharmonicsTerms = has_anharmonicsTerms
 else
   eff_pot%has_anharmonicsTerms = .false.
 end if

!Allocation of phonon strain coupling array (3rd order)
 if(present(phonon_strain).and.has_anharmonicsTerms) then
   call anharmonics_terms_setStrainPhononCoupling(eff_pot%anharmonics_terms,crystal%natom,&
&                                                 phonon_strain)
 end if

!Set the 3rd order elastic tensor
 if(present(elastic3rd).and.has_anharmonicsTerms)then
   call anharmonics_terms_setElastic3rd(eff_pot%anharmonics_terms,elastic3rd)
 end if

!Allocation of 3rd order with respecto to 2 strain and 1 atomic displacement
 if(present(elastic_displacement).and.has_anharmonicsTerms)then
   call anharmonics_terms_setElasticDispCoupling(eff_pot%anharmonics_terms,crystal%natom,&
&                                                elastic_displacement)
 end if

!Allocation of the coefficients
 if(present(coeffs))then
   if(ncoeff /= size(coeffs))then
     MSG_BUG('ncoeff has not the same size than coeffs array')
   end if
   call effective_potential_setCoeffs(coeffs,eff_pot,ncoeff)
 end if

 if(present(supercell))then
   call effective_potential_setSupercell(eff_pot,comm,supercell=supercell)
 else
!   call init_supercell(eff_pot%crystal%natom, (/1,0,0, 0,1,0,  0,0,1/), eff_pot%crystal%rprimd,&
!&                      eff_pot%crystal%typat, eff_pot%crystal%xcart, eff_pot%crystal%znucl, eff_pot%supercell)
   call effective_potential_setSupercell(eff_pot,comm,ncell=(/1,1,1/))
 end if

!Set the confinement potential
 if(present(polynomial_conf)) then
   call effective_potential_setConfinement(polynomial_conf%cutoff_disp,polynomial_conf%cutoff_strain,&
&                                          eff_pot,polynomial_conf%factor_disp,&
&                                          polynomial_conf%factor_strain,polynomial_conf%ndisp,&
&                                          polynomial_conf%power_disp,polynomial_conf%power_strain,&
&                                          polynomial_conf%need_confinement)
 end if

end subroutine effective_potential_init
!!***


!!****f* m_effective_potential/effective_potential_initmpi
!! NAME
!!  effective_potential_initmpi
!!
!! FUNCTION
!!  Initializes the mpi informations for parallelism over supercell.
!!  Only the parallelisation over cell is done here.
!!  The parallelisation over cell and coeff is disable for now (experimental)
!!
!! INPUTS
!!  eff_pot<type(effective_potential_type)> = datatype for the effective potential
!!  comm = MPI communicator
!!
!! OUTPUT
!! This is for the parallelisation over the supercell
!!  eff_pot%mpi_ifc%me_supercell =  Index of my processor in the comm. over one cell
!!  eff_pot%mpi_ifc%my_ncell     =  Number of cell treated by current proc
!!  eff_pot%mpi_ifc%my_cells(:)  = Number of the cells in the supercell treat by this CPU
!!  eff_pot%mpi_ifc%my_index_cells(:,:) = indexes of the cells in the supercell treat by this CPU
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_initmpi(eff_pot,comm)

!Arguments ------------------------------------
!scalars
 type(effective_potential_type),intent(inout)  :: eff_pot
 integer,intent(in) :: comm
!arrays

!Local variables-------------------------------
!scalars
 integer :: ndiv
 integer :: ncell
!array
 integer :: cell_number(3)
 character(len=500) :: msg
! ***********************************************************************

!Set the number of cell in the supercell
 cell_number(1) = eff_pot%supercell%rlatt(1,1)
 cell_number(2) = eff_pot%supercell%rlatt(2,2)
 cell_number(3) = eff_pot%supercell%rlatt(3,3)
 ncell = product(cell_number(:))

!Do some checks
 if (any(cell_number <= 0).or.ncell<=0) then
   write(msg,'(a,a)')' No supercell found for setting'
   MSG_ERROR(msg)
 end if

!First mpi_ifc
 ndiv = 1
 call effpot_mpi_free(eff_pot%mpi_ifc)
 call effpot_mpi_init(eff_pot%harmonics_terms%ifcs%cell,cell_number,eff_pot%mpi_ifc,&
&                     eff_pot%crystal%natom,ndiv,eff_pot%harmonics_terms%ifcs%nrpt,comm)

!Second mpi_coeff
 ndiv = 1
 call effpot_mpi_free(eff_pot%mpi_coeff)
 call effpot_mpi_init(eff_pot%harmonics_terms%ifcs%cell,cell_number,eff_pot%mpi_coeff,&
&                     eff_pot%crystal%natom,ndiv,eff_pot%harmonics_terms%ifcs%nrpt,comm)

end subroutine effective_potential_initmpi
!!***


!****f* m_effective_potential/effective_potential_free
!!
!! NAME
!! effective_potential_free
!!
!! FUNCTION
!! deallocate all dynamic memory for this effective potential datatype
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)>  = effective_potential datatype
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)>  = effective_potential datatype
!!
!! PARENTS
!!      compute_anharmonics,m_effective_potential,m_effective_potential_file
!!      multibinit
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_free(eff_pot)

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
  eff_pot%strten = zero
  eff_pot%has_anharmonicsTerms = .FALSE.
  eff_pot%anharmonics_terms%bounded = .FALSE.

   if(allocated(eff_pot%fcart)) then
     eff_pot%fcart=zero
     ABI_DEALLOCATE(eff_pot%fcart)
   end if

! Free others datatypes
   call anharmonics_terms_free(eff_pot%anharmonics_terms)
   call harmonics_terms_free(eff_pot%harmonics_terms)
   call destroy_supercell(eff_pot%supercell)
   call eff_pot%crystal%free()
   call effective_potential_freempi(eff_pot)
   call polynomial_conf_free(eff_pot%confinement)

end subroutine effective_potential_free
!!***

!****f* m_effective_potential/effective_potential_freeCoeffs
!!
!! NAME
!! effective_potential_freeCoeffs
!!
!! FUNCTION
!! deallocate all dynamic memory for the coefficients of
!! the anharmonic part of the effective potential datatype
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)>  = effective_potential datatype
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)>  = effective_potential datatype
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_freeCoeffs(eff_pot)

!Arguments ------------------------------------
!scalars
!array
  type(effective_potential_type), intent(inout) :: eff_pot
!Local variables-------------------------------
!scalars
!array
! *************************************************************************

  call anharmonics_terms_freeCoeffs(eff_pot%anharmonics_terms)

end subroutine effective_potential_freeCoeffs
!!***

!****f* m_effective_potential/effective_potential_freempi
!!
!! NAME
!! effective_potential_freempi
!!
!! FUNCTION
!! deallocate all dynamic memory for mpi of supercell
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)>  = effective_potential datatype
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)>  = effective_potential datatype
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_freempi(eff_pot)

!Arguments ------------------------------------
!scalars
!array
 type(effective_potential_type), intent(inout) :: eff_pot

!Local variables-------------------------------
!scalars
!array

! *************************************************************************
 call effpot_mpi_free(eff_pot%mpi_ifc)
 call effpot_mpi_free(eff_pot%mpi_coeff)

end subroutine effective_potential_freempi
!!***

!****f* m_effective_potential/effective_potential_generateDipDip
!!
!! NAME
!! effective_potential_generateDipDip
!!
!! FUNCTION
!! Generate the supercell of the structure inside of the effective
!! potential and fill the supercell datatype. Also adapt the harmonic
!! part for the supercell (compute dipole-dipole interation)
!!
!! INPUTS
!! eff_pot = effective potential datatype
!! option  =  0 Just generate supercell and fill effective potential
!!            1 Regenerate the dipole dipo;le interaction
!! ncell(3) = number of cell in the direction x, y and z
!! comm=MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = effective potential datatype
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_generateDipDip(eff_pot,ncell,option,asr,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,asr
 integer,intent(in) :: comm
!array
 integer,intent(in) :: ncell(3)
 type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
 integer,parameter :: master=0
 integer :: first_coordinate
 integer :: ia,i1,i2,i3,ii,ierr,irpt,irpt2,irpt_ref,min1,min2,min3
 integer :: min1_cell,min2_cell,min3_cell,max1_cell,max2_cell,max3_cell
 integer :: max1,max2,max3,my_rank,natom_uc
 integer :: nproc,second_coordinate,size_tmp,sumg0
 integer :: my_nrpt,nrpt_alone
 real(dp) :: ucvol
 character(len=500) :: msg
 logical :: iam_master=.FALSE.
!array
 integer,allocatable :: my_index_rpt(:,:)
 integer,allocatable :: bufsize(:),bufdisp(:)
 integer,allocatable :: my_irpt(:)
 real(dp) :: acell(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: buff_ewald(:,:,:,:,:,:),dyew(:,:,:,:,:), dyewq0(:,:,:,:,:)
 real(dp),allocatable :: xred(:,:),xred_tmp(:,:),zeff_tmp(:,:,:),qdrp_cart(:,:,:,:)
 type(supercell_type) :: supercell
 type(ifc_type) :: ifc_tmp
 integer, allocatable :: full_cell(:,:)
 integer :: full_nrpt
 real(dp), allocatable :: full_cell_atmfrc(:,:,:,:,:), full_cell_short_atmfrc(:,:,:,:,:), full_cell_ewald_atmfrc(:,:,:,:,:)

! *************************************************************************

!0 MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)
 ierr=0

!0 Check the size of the cell
 do ia=1,3
   if(ncell(ia)<0.or.ncell(ia)>100)then
     write(msg, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ncell(',ia,') is ',ncell(ia),', which is lower than 0 of superior than 50.',&
&     ch10,'Action: correct ncell(',ia,').'
     MSG_ERROR(msg)
   end if
 end do

 call init_supercell(eff_pot%crystal%natom, &
&   (/ncell(1),0,0,  0,ncell(2),0,  0,0,ncell(3)/),&
&   eff_pot%crystal%rprimd,&
&   eff_pot%crystal%typat,&
&   eff_pot%crystal%xcart,&
&   eff_pot%crystal%znucl,&
&   supercell)

!set variables
 natom_uc = eff_pot%crystal%natom

!1 Store the information of the supercell of the reference structure into effective potential
 call copy_supercell(supercell,eff_pot%supercell)
!2 Initialisation of new mpi over supercell
 call effective_potential_initmpi(eff_pot,comm)

!3 Check if the bound of new cell correspond to the effective potential
!only for option=zero
!set min and max
 min1 = minval(eff_pot%harmonics_terms%ifcs%cell(1,:))
 min2 = minval(eff_pot%harmonics_terms%ifcs%cell(2,:))
 min3 = minval(eff_pot%harmonics_terms%ifcs%cell(3,:))
 max1 = maxval(eff_pot%harmonics_terms%ifcs%cell(1,:))
 max2 = maxval(eff_pot%harmonics_terms%ifcs%cell(2,:))
 max3 = maxval(eff_pot%harmonics_terms%ifcs%cell(3,:))
   write(msg,'(5a,2I3,a,2I3,a,2I3,a)') ch10,' Bound for ifc SR:',&
&    ch10,ch10, " x=[",min1,max1,"], y=[",min2,max2,"] and z=[",min3,max3,"]"
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

 if(option==0) then
   if(((max1-min1+1)/=ncell(1).and.&
&    (max2-min2+1)/=ncell(2).and.(max3-min3+1)/=ncell(3))) then
     write(msg, '(90a,3I3,5a,3I3,3a)' )ch10,('-',i1=1,80),ch10,ch10,&
&      ' --- !WARNING:',ch10,&
&      '     dipdip is set to zero, the longe range interation might be wrong',ch10,&
&      '     because it is not recompute.',ch10,&
&      '     The previous harmonic part is build for ',(max1-min1+1),(max2-min2+1),(max3-min3+1)&
&,     ' cell.',ch10,&
&      '     Be sure than the dipole-dipole interation is correct ',ch10,&
&      '     for the supercell: ',ncell(:),' or set dipdip to 1',ch10,&
&      ' ---'
     call wrtout(std_out,msg,"COLL")
   else
     write(msg,'(89a)')ch10,('-',i1=1,80),ch10,ch10,&
&          ' --- !WARNING:',ch10,&
&          '     dipdip is set to zero, the longe range interation is not recompute.',ch10,&
&          ' ---',ch10
     call wrtout(std_out,msg,"COLL")
   end if


   write(msg,'(a,(80a))') ch10,('=',i1=1,80)
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

!4-Adapt harmonic part
 else if (option>=1.and.all(ncell(:)>0)) then

   write(msg,'(a,(80a),3a)') ch10,('=',i1=1,80),ch10,' Generation of new ifc',ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

   irpt_ref = 0
   irpt     = 0
   min1_cell = 0; max1_cell = 0
   min2_cell = 0; max2_cell = 0
   min3_cell = 0; max3_cell = 0

   ! ncell
   call findBound_supercell(min1_cell,max1_cell,ncell(1))
   call findBound_supercell(min2_cell,max2_cell,ncell(2))
   call findBound_supercell(min3_cell,max3_cell,ncell(3))

   write(msg, '(2a)' )&
&        ' dipdip is set to one, the dipole-dipole interation is recompute.'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

!!  Generate new bound
!   !if(option==1)then
!     ! If NGQPT > NCELL "=" dipdip range
!     !if ((abs(min1) > abs(min1_cell)).or.(abs(max1) > abs(max1_cell)).or.&
!&    !    (abs(min2) > abs(min2_cell)).or.(abs(max2) > abs(max2_cell)).or.&
!&    !    (abs(min3) > abs(min3_cell)).or.(abs(max3) > abs(max3_cell))) then
!     !  write(msg, '(6a,3I4,3a)' )ch10,&
!&    !    ' --- !WARNING',ch10,&
!&    !    '     The range of the dipole-dipole interaction is the same than the short-range.',ch10,&
!&    !    '     So the range of the total ifc is ',int((/(max1-min1+1),(max2-min2+1),(max3-min3+1)/),dp),' cell',ch10,&
!&    !    ' ---'
!     !  call wrtout(std_out,msg,"COLL")
!     !  if (abs(min1) < abs(min1_cell)) min1 = min1_cell
!     !  if (abs(min2) < abs(min2_cell)) min2 = min2_cell
!     !  if (abs(min3) < abs(min3_cell)) min3 = min3_cell
!     !  if (abs(max1) < abs(max1_cell)) max1 = max1_cell
!     !  if (abs(max2) < abs(max2_cell)) max2 = max2_cell
!     !  if (abs(max3) < abs(max3_cell)) max3 = max3_cell
!!    !  If the cell is smaller, we redifine new cell to take into acount all atoms
!     !  call destroy_supercell(supercell)
!     !  call init_supercell(natom_uc,(/(max1-min1+1),0,0,  0,(max2-min2+1),0,  0,0,(max3-min3+1)/),&
!&    !                      eff_pot%crystal%rprimd,eff_pot%crystal%typat,&
!&    !                      eff_pot%crystal%xcart,eff_pot%crystal%znucl,supercell)
!
!!    !  Store the information of the supercell of the reference structure into effective potential
!     !  call effective_potential_setSupercell(eff_pot,comm,supercell=supercell)
!     !else
!     !  min1 = min1_cell ; min2 = min2_cell ; min3 = min3_cell
!     !  max1 = max1_cell ; max2 = max2_cell ; max3 = max3_cell
!     !end if
!   !end if

!  Print the new boundary
   write(msg,'(5a,2I3,a,2I3,a,2I3,4a)') ch10,' Bound for ifc (LR):',&
&    ch10,ch10, " x=[",min1_cell,max1_cell,"], y=[",min2_cell,max2_cell,"] and z=[",min3_cell,max3_cell,"]",ch10,ch10,&
&    " Computation of new dipole-dipole interaction."
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

!  Count the new number of ifc
   do i1=min1_cell,max1_cell
     do i2=min2_cell,max2_cell
       do i3=min3_cell,max3_cell
         irpt = irpt +1
         if(i1==0.and.i2==0.and.i3==0) irpt_ref = irpt
       end do
     end do
   end do

   ifc_tmp%nrpt = irpt
   ABI_ALLOCATE(ifc_tmp%cell,(3,ifc_tmp%nrpt))
   ifc_tmp%cell(:,:) = 0

!  Set MPI here and not at the begining because the number of cell is adjust just before
!  Here we store in my_irpt a list with the number of each cell to be treat by this CPU
!  Determine the number of cell for each CPU
   ABI_ALLOCATE(bufsize,(nproc))
   ABI_ALLOCATE(bufdisp,(nproc))

   nrpt_alone = mod(ifc_tmp%nrpt,nproc)
   my_nrpt = int(real(ifc_tmp%nrpt,sp)/nproc)
   if(my_rank >= (nproc-nrpt_alone)) then
     my_nrpt = my_nrpt  + 1
   end if

!  Initialisation of ifc temporary
   ABI_ALLOCATE(buff_ewald,(2,3,natom_uc,3,natom_uc,my_nrpt))
   ABI_ALLOCATE(ifc_tmp%short_atmfrc,(3,natom_uc,3,natom_uc,ifc_tmp%nrpt))
   ABI_ALLOCATE(ifc_tmp%ewald_atmfrc,(3,natom_uc,3,natom_uc,ifc_tmp%nrpt))
   ABI_ALLOCATE(ifc_tmp%atmfrc,(3,natom_uc,3,natom_uc,ifc_tmp%nrpt))
   ABI_ALLOCATE(my_irpt,(my_nrpt))
   ABI_ALLOCATE(my_index_rpt,(3,my_nrpt))

   my_irpt = 0
   my_index_rpt(:,:) = 0
   ifc_tmp%atmfrc(:,:,:,:,:) = zero
   ifc_tmp%short_atmfrc(:,:,:,:,:) = zero
   ifc_tmp%ewald_atmfrc(:,:,:,:,:) = zero
   buff_ewald(:,:,:,:,:,:) = zero

!  Allocation of array
   do irpt = 1,my_nrpt
     if(my_rank >= (nproc-nrpt_alone))then
       my_irpt(irpt)=(int(real(ifc_tmp%nrpt,sp)/nproc))*(my_rank)+&
&                       (my_rank - (nproc-nrpt_alone)) + irpt
     else
       my_irpt(irpt)=(my_nrpt)*(my_rank) + irpt
     end if
   end do

   irpt = 0
   ii = 0
   do i1=min1_cell,max1_cell
     do i2=min2_cell,max2_cell
       do i3=min3_cell,max3_cell
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
   ABI_ALLOCATE(xred_tmp,(3,2*natom_uc))
   ABI_ALLOCATE(xred,(3,supercell%natom))
   ABI_ALLOCATE(zeff_tmp,(3,3,2*natom_uc))
   ABI_ALLOCATE(qdrp_cart,(3,3,3,2*natom_uc))
   ABI_ALLOCATE(dyew,(2,3,2*natom_uc,3,2*natom_uc))
   ABI_ALLOCATE(dyewq0,(2,3,natom_uc,3,natom_uc))

   dyew            = zero
   dyewq0          = zero
   xred(:,:)       = zero
   xred_tmp(:,:)   = zero
   zeff_tmp(:,:,:) = zero
   qdrp_cart       = zero
   sumg0           = 0
   acell           = one

   call matr3inv(supercell%rprimd,gprimd)
   call xcart2xred(supercell%natom,supercell%rprimd,&
&                  supercell%xcart,xred)
   call metric(gmet,gprimd,-1,rmet,supercell%rprimd,ucvol)

!  Fill the atom position of the first cell (reference cell)
   first_coordinate  = ((irpt_ref-1)*natom_uc) + 1
   second_coordinate = first_coordinate + natom_uc-1
   xred_tmp(:,1:natom_uc) = xred(:,first_coordinate:second_coordinate)
!  Fill fake zeff array for ewald9
   zeff_tmp(:,:,1:natom_uc) = eff_pot%harmonics_terms%zeff
   zeff_tmp(:,:,natom_uc+1:2*natom_uc) = eff_pot%harmonics_terms%zeff

   do irpt=1,my_nrpt
     i1=my_index_rpt(1,irpt); i2=my_index_rpt(2,irpt); i3=my_index_rpt(3,irpt)
!    Compute new dipole-dipole interaction
     dyew = zero
     if (i1==0.and.i2==0.and.i3==0) then
       call ewald9(acell,eff_pot%harmonics_terms%epsilon_inf,dyewq0,&
&                  gmet,gprimd,natom_uc,real((/0,0,0/),dp),rmet,&
&                  supercell%rprimd,sumg0,ucvol,xred_tmp(:,1:natom_uc),&
&                  eff_pot%harmonics_terms%zeff,qdrp_cart)
       buff_ewald(:,:,:,:,:,irpt) = dyewq0
     else
       first_coordinate  = ((my_irpt(irpt)-1)*natom_uc) + 1
       second_coordinate = first_coordinate + natom_uc  - 1
       xred_tmp(:,natom_uc+1:2*natom_uc)=&
&              xred(:,first_coordinate:second_coordinate)
       call ewald9(acell,eff_pot%harmonics_terms%epsilon_inf,dyew,gmet,gprimd,&
&                  int(2*natom_uc),real((/0,0,0/),dp),&
&                  rmet,supercell%rprimd,&
&                  sumg0,ucvol,xred_tmp,zeff_tmp,qdrp_cart)
       buff_ewald(:,:,:,:,:,irpt) = &
&           dyew(:,:,1:natom_uc,:,natom_uc+1:2*natom_uc)
     end if
   end do

!  DEALLOCATION OF ARRAYS
   ABI_DEALLOCATE(my_index_rpt)
   ABI_DEALLOCATE(my_irpt)
   ABI_DEALLOCATE(xred_tmp)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(zeff_tmp)
   ABI_DEALLOCATE(qdrp_cart)
   ABI_DEALLOCATE(dyew)
   ABI_DEALLOCATE(dyewq0)

!  Set the bufsize for mpi allgather
   do ii = 1,nproc
     bufsize(ii) = int(real(ifc_tmp%nrpt,sp)/nproc)*3*natom_uc*3*natom_uc
     if(ii > (nproc-nrpt_alone)) then
       bufsize(ii) = bufsize(ii) + 3*natom_uc*3*natom_uc
     end if
   end do

   bufdisp(1) = 0
   do ii = 2,nproc
     bufdisp(ii) = bufdisp(ii-1) + bufsize(ii-1)
   end do

   size_tmp = 3*natom_uc*3*natom_uc*my_nrpt
   call xmpi_allgatherv(buff_ewald(1,:,:,:,:,:),size_tmp,ifc_tmp%ewald_atmfrc,bufsize,bufdisp, comm, ierr)

   ABI_DEALLOCATE(bufsize)
   ABI_DEALLOCATE(bufdisp)
   ABI_DEALLOCATE(buff_ewald)

!  Fill the short range part (calculated previously) only master
!   if(iam_master)then
!     do irpt=1,ifc_tmp%nrpt
!       do irpt2=1,eff_pot%harmonics_terms%ifcs%nrpt
!         if(eff_pot%harmonics_terms%ifcs%cell(1,irpt2)==ifc_tmp%cell(1,irpt).and.&
!&           eff_pot%harmonics_terms%ifcs%cell(2,irpt2)==ifc_tmp%cell(2,irpt).and.&
!&           eff_pot%harmonics_terms%ifcs%cell(3,irpt2)==ifc_tmp%cell(3,irpt).and.&
!&           any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt2)) > tol20)) then
!           ifc_tmp%short_atmfrc(:,:,:,:,irpt) = &
!&                               eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt2)
!         end if
!       end do
!     end do
!   end if


   !call xmpi_bcast(ifc_tmp%short_atmfrc, master, comm, ierr)
   ! Maybe useless
   call xmpi_bcast(eff_pot%harmonics_terms%ifcs%short_atmfrc, master, comm, ierr)

!  Compute total ifc
   !ifc_tmp%atmfrc = ifc_tmp%short_atmfrc + ifc_tmp%ewald_atmfrc
   ! Set the full cell according to the largest box
   full_nrpt = max(ifc_tmp%nrpt,eff_pot%harmonics_terms%ifcs%nrpt)
   ABI_ALLOCATE(full_cell,(3,full_nrpt))
   ABI_CALLOC(full_cell_atmfrc,(3,natom_uc,3,natom_uc,full_nrpt)) ! Allocate and set to 0
   ABI_CALLOC(full_cell_short_atmfrc,(3,natom_uc,3,natom_uc,full_nrpt)) ! Allocate and set to 0
   ABI_CALLOC(full_cell_ewald_atmfrc,(3,natom_uc,3,natom_uc,full_nrpt)) ! Allocate and set to 0
   if ( full_nrpt == ifc_tmp%nrpt ) then
     full_cell = ifc_tmp%cell
   else
     full_cell = eff_pot%harmonics_terms%ifcs%cell
   end if

   ! Copy LR into total_atmfrc
   do irpt=1,ifc_tmp%nrpt ! LR IRPT
     do irpt2=1, full_nrpt
       if(  ifc_tmp%cell(1,irpt)==full_cell(1,irpt2).and.&
&           ifc_tmp%cell(2,irpt)==full_cell(2,irpt2).and.&
&           ifc_tmp%cell(3,irpt)==full_cell(3,irpt2) ) then
         full_cell_atmfrc(:,:,:,:,irpt2) = ifc_tmp%ewald_atmfrc(:,:,:,:,irpt)
         full_cell_ewald_atmfrc(:,:,:,:,irpt2) = ifc_tmp%ewald_atmfrc(:,:,:,:,irpt)
       end if
     end do
   end do
   ! Copy SR into total_atmfrc
   do irpt=1,eff_pot%harmonics_terms%ifcs%nrpt ! SR IRPT
     do irpt2=1, full_nrpt
     if(  eff_pot%harmonics_terms%ifcs%cell(1,irpt)==full_cell(1,irpt2).and.&
&         eff_pot%harmonics_terms%ifcs%cell(2,irpt)==full_cell(2,irpt2).and.&
&         eff_pot%harmonics_terms%ifcs%cell(3,irpt)==full_cell(3,irpt2) ) then
         full_cell_atmfrc(:,:,:,:,irpt2) = full_cell_atmfrc(:,:,:,:,irpt2) + eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt)
         full_cell_short_atmfrc(:,:,:,:,irpt2) = eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt)
       end if
     end do
   end do

!  Count the rpt inferior to the tolerance
   irpt2 = 0
   do irpt=1,full_nrpt
     if(any(abs(full_cell_atmfrc(:,:,:,:,irpt)) > tol20))then
       irpt2 = irpt2 + 1
     end if
   end do

!  Copy ifc into effective potential
!  !!!Warning eff_pot%harmonics_terms%ifcs only contains atmfrc,short_atmfrc,ewald_atmfrc,nrpt
!    and cell!!  rcan,ifc%rpt,wghatm and other quantities
!    are not needed for effective potential!!!
!  Free ifc before copy
   call eff_pot%harmonics_terms%ifcs%free()

!  Fill the effective potential with new atmfr
    eff_pot%harmonics_terms%ifcs%nrpt = irpt2
    ABI_ALLOCATE(eff_pot%harmonics_terms%ifcs%atmfrc,(3,natom_uc,3,natom_uc,irpt2))
    ABI_ALLOCATE(eff_pot%harmonics_terms%ifcs%short_atmfrc,(3,natom_uc,3,natom_uc,irpt2))
    ABI_ALLOCATE(eff_pot%harmonics_terms%ifcs%ewald_atmfrc,(3,natom_uc,3,natom_uc,irpt2))
    ABI_ALLOCATE(eff_pot%harmonics_terms%ifcs%cell,(3,irpt2))

    irpt2 = 0
    do irpt = 1,full_nrpt
     if(any(abs(full_cell_atmfrc(:,:,:,:,irpt)) > tol20))then
       irpt2 = irpt2 + 1
       eff_pot%harmonics_terms%ifcs%atmfrc(:,:,:,:,irpt2) = full_cell_atmfrc(:,:,:,:,irpt)
       eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt2) = full_cell_short_atmfrc(:,:,:,:,irpt)
       eff_pot%harmonics_terms%ifcs%ewald_atmfrc(:,:,:,:,irpt2) = full_cell_ewald_atmfrc(:,:,:,:,irpt)
       eff_pot%harmonics_terms%ifcs%cell(:,irpt2) = full_cell(:,irpt)
     end if
   end do

!  Free temporary ifc
   call ifc_tmp%free()
   !ABI_DEALLOCATE(wghatm)
   !Deallocate temporary arrays
   ABI_DEALLOCATE(full_cell)
   ABI_DEALLOCATE(full_cell_atmfrc) ! Allocate and set to 0
   ABI_DEALLOCATE(full_cell_short_atmfrc) ! Allocate and set to 0
   ABI_DEALLOCATE(full_cell_ewald_atmfrc) ! Allocate and set to 0

 end if

 if(asr >= 0) then
! Impose sum rule
   call harmonics_terms_applySumRule(asr,eff_pot%harmonics_terms%ifcs,natom_uc)
 end if

 write(msg, '(a,(80a),a)' ) ch10,&
&   ('=',ii=1,80)
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

! Free suppercell
 call destroy_supercell(supercell)

end subroutine effective_potential_generateDipDip
!!***

!****f* m_effective_potential/effective_potential_setCoeffs
!!
!! NAME
!! effective_potential_setCoeffs
!!
!! FUNCTION
!! Set the polynomial coefficients of the effective_potential
!!
!! INPUTS
!! coeffs(ncoeff)<type(polynomial_coeff)> = array of polynomial_coeff datatype
!! ncoeff = number of coefficient
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype for effective potential
!!
!! PARENTS
!!      m_effective_potential,m_effective_potential_file,m_fit_polynomial_coeff
!!      mover_effpot
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setCoeffs(coeffs,eff_pot,ncoeff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncoeff
!array
 type(effective_potential_type),intent(inout) :: eff_pot
 type(polynomial_coeff_type),intent(in) :: coeffs(ncoeff)
!Local variables-------------------------------
!scalar
 integer :: ii,jj
 logical :: has_straincoupling=.FALSE.
 character(len=500) :: msg
!array
! *************************************************************************

 if(ncoeff /= size(coeffs))then
   write(msg, '(a)' )&
&       ' ncoeff has not the same size than coeffs array, '
   MSG_BUG(msg)
 end if

! Check if the strain coupling is present
 do ii=1,ncoeff
   do jj=1,coeffs(ii)%nterm
     if (coeffs(ii)%terms(jj)%nstrain > 0) then
       has_straincoupling = .TRUE.
     end if
   end do
 end do

! Set to false the strain coupling from the finite differences
 if(has_straincoupling)then
   if(eff_pot%has_anharmonicsTerms) then
     write(msg, '(8a)' )ch10,&
&        ' --- !WARNING',ch10,&
&        '     There is strain coupling with the fitted coefficients,',ch10,&
&        '     The previous contribution will be set to zero',ch10,&
&        ' ---'
     call wrtout(std_out,msg,"COLL")
   end if
   eff_pot%has_anharmonicsTerms = .FALSE.
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
!! elastics(6,6,6) = 3d order of elastics constant
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype for effective potential
!!
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setElastic3rd(eff_pot,elastics)

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

  if(any(abs(eff_pot%anharmonics_terms%elastic3rd)> tol15)) then
    eff_pot%has_anharmonicsTerms = .TRUE.
  end if

end subroutine effective_potential_setElastic3rd
!!***

!****f* m_effective_potential/effective_potential_setElastic4th
!!
!! NAME
!! effective_potential_setElastic4th
!!
!! FUNCTION
!! Set the 4th order derivative of with respect to 4 strain
!!
!! INPUTS
!! elastics(6,6,6,6) = 4th order of elastics constant
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype for effective potential
!!
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setElastic4th(eff_pot,elastics)

!Arguments ------------------------------------
!scalars
!array
  real(dp),intent(in) :: elastics(6,6,6,6)
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
!array
! *************************************************************************
  call anharmonics_terms_setElastic4th(eff_pot%anharmonics_terms,elastics)

  if(any(abs(eff_pot%anharmonics_terms%elastic4th)> tol15)) then
    eff_pot%has_anharmonicsTerms = .TRUE.
  end if


end subroutine effective_potential_setElastic4th
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
!! strain_phonon(6)<type(ifc_type) = strain-phonon coupling
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype for effective potential
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setStrainPhononCoupling(eff_pot,natom,phonon_strain)

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(ifc_type),intent(in) :: phonon_strain(6)
  type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
  integer :: ii
!array
! *************************************************************************

  call anharmonics_terms_setStrainPhononCoupling(eff_pot%anharmonics_terms,natom,phonon_strain)
  do ii=1,6
    if(any(abs(eff_pot%anharmonics_terms%phonon_strain(ii)%atmfrc)> tol15)) then
      eff_pot%has_anharmonicsTerms = .True.
    end if
  end do

end subroutine effective_potential_setStrainPhononCoupling
!!***

!****f* m_effective_potential/effective_potential_setElasticDispCoupling
!!
!! NAME
!! effective_potential_setElasticDispCoupling
!!
!! FUNCTION
!! Set the elastic constant displacement coupling of the effective_potential
!!
!! INPUTS
!! natom  = number of atoms
!! elastic_displacement = (6,6,3,natom) array with elastic constant displacement coupling
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype for effective potential
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setElasticDispCoupling(eff_pot,natom,elastic_displacement)

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
  if(any(abs(eff_pot%anharmonics_terms%elastic_displacement)> tol15)) then
    eff_pot%has_anharmonicsTerms = .True.
  end if

end subroutine effective_potential_setElasticDispCoupling
!!***

!!****f* m_effective_potential/effective_potential_setConfinement
!!
!! NAME
!! effective_potential_setConfinement
!!
!! FUNCTION
!! Set the confinement in the effective_potential datatype
!!
!! INPUTS
!! cutoff_disp(6) = Cutoff array for the strain
!! cutoff_strain(ndisp) = Cutoff array for the atomic displacement
!! factor_disp = Factor to appy to the polynomial term of the confinement (displacement)
!! factor_strain = Factor to appy to the polynomial term of the confinement (strain)
!! ndisp = Number of displacement (atoms) for the cut off
!! power_disp = Power of the polynome related to the displacement
!! power_strain = Power of the polynome related to the strain
!! need_confinement = optional,Logical related to the necessity of the confinement
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = datatype for effective potential
!!
!! PARENTS
!!      m_effective_potential,multibinit
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setConfinement(cutoff_disp,cutoff_strain,eff_pot,factor_disp,&
&                                             factor_strain,ndisp,power_disp,power_strain,&
&                                             need_confinement)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: power_disp,power_strain,ndisp
 real(dp),intent(in) :: factor_disp,factor_strain
 logical,optional,intent(in)  :: need_confinement
!arrays
 real(dp),intent(in) :: cutoff_disp(ndisp),cutoff_strain(6)
 type(effective_potential_type),intent(inout) :: eff_pot
!Local variables-------------------------------
!scalar
 logical  :: need_confinement_tmp = .FALSE.
!arrays
 character(len=500) :: msg

! *************************************************************************

!Checks
 if (ndisp <= 0) then
   write(msg,'(a,a)')' ndisp can not be inferior or equal to zero'
   MSG_ERROR(msg)
 end if

!First free the type
 call  polynomial_conf_free(eff_pot%confinement)

 if (present(need_confinement)) need_confinement_tmp = need_confinement

 call  polynomial_conf_init(cutoff_disp,cutoff_strain,factor_disp,factor_strain,ndisp,&
&                           eff_pot%confinement,power_disp,power_strain,&
&                           need_confinement=need_confinement_tmp)


end subroutine effective_potential_setConfinement
!!***

!!****f* m_effective_potential/effective_potential_setSupercell
!!
!! NAME
!! effective_potential_setSupercell
!!
!! FUNCTION
!! Set the supercell type in the effective_potential type
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!! comm = MPI communicator
!! ncell(3) = optional, size of the supercell
!! supercell = optional, supercell type to set to eff_pot
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!!
!! PARENTS
!!      m_effective_potential,m_effective_potential_file,mover_effpot
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_setSupercell(eff_pot,comm,ncell,supercell)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
!arrays
 type(effective_potential_type),intent(inout) :: eff_pot
 integer,optional,intent(in) :: ncell(3)
 type(supercell_type),optional,intent(in) :: supercell
!Local variables-------------------------------
!scalar
!arrays
 character(len=500) :: msg

! *************************************************************************

!Checks
 if (.not.present(supercell).and..not.present(ncell)) then
   write(msg,'(a,a)')' You should at least set ncell of supercell type'
   MSG_ERROR(msg)
 end if

 call destroy_supercell(eff_pot%supercell)

 if(present(supercell))then
   call copy_supercell(supercell,eff_pot%supercell)
 else
   call init_supercell(eff_pot%crystal%natom, (/ncell(1),0,0,  0,ncell(2),0,  0,0,ncell(3)/), &
&                      eff_pot%crystal%rprimd,eff_pot%crystal%typat,eff_pot%crystal%xcart,&
&                      eff_pot%crystal%znucl, eff_pot%supercell)
 end if

!Initialisation of new mpi over supercell
 call effective_potential_initmpi(eff_pot,comm)

end subroutine effective_potential_setSupercell
!!***

!****f* m_effective_potential/effective_potential_print
!!
!! NAME
!! effective_potential_print
!!
!! FUNCTION
!! Print the effective_potential in the ouput
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!! option  =  0 no output
!! option  =  1 only useful information
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_print(eff_pot,option,filename)

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
  character(len=500) :: msg
!array
! *************************************************************************

  if(option >= 1) then
    if(present(filename)) then
      write(msg, '(a,a,a,a,a,a)' )ch10,' The file ',trim(filename),&
&     ' contains this effective potential for ',trim(eff_pot%name),':'
    else
      write(msg, '(a,a,a,a)' )ch10,' This effective potential contains ',&
&     trim(eff_pot%name),':'
    end if

    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out,msg,'COLL')

!**********************************************************************
! Write basics values
!**********************************************************************

    write(msg,'(a,F20.10,2a,I3,2a,I4,2a,I4,2a,I3,2a)') &
&     '  - Reference energy:  ',eff_pot%energy ,ch10,&
&     '  - Number of types of atoms:  ',eff_pot%crystal%ntypat ,ch10,&
&     '  - Number of atoms:  ',eff_pot%crystal%natom ,ch10,&
&     '  - Number of cells:  ',eff_pot%harmonics_terms%ifcs%nrpt ,ch10,&
&     '  - Number of qpoints:  ',eff_pot%harmonics_terms%nqpt ,ch10,&
&     '  - Primitive vectors (unit:Bohr):  '
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
    do ii = 1,3
      write(msg,'(3F12.6)') eff_pot%crystal%rprimd(1,ii),&
&                               eff_pot%crystal%rprimd(2,ii),&
&                               eff_pot%crystal%rprimd(3,ii)
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
    end do
    write(msg,'(2a,3F12.6)') '  - acell (unit:Bohr):',ch10,one,one,one

    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
    write(msg,'(a)') '  - Dielectric tensor:  '
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
    do ii=1,3
      write(msg,'(3F12.6)')eff_pot%harmonics_terms%epsilon_inf(1,ii),&
&                              eff_pot%harmonics_terms%epsilon_inf(2,ii),&
&                              eff_pot%harmonics_terms%epsilon_inf(3,ii)
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
    end do
      write(msg,'(a)') '  - Elastic tensor (unit:10^2GPa):  '
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
      fact = HaBohr3_GPa / eff_pot%crystal%ucvol
    do ii=1,6
      write(msg,'(6F12.6)')&
&      eff_pot%harmonics_terms%elastic_constants(1,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(2,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(3,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(4,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(5,ii)*fact/100,&
&      eff_pot%harmonics_terms%elastic_constants(6,ii)*fact/100
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
    end do
    do ia=1,eff_pot%crystal%natom
      write(msg,'(a,I4,2a,F10.4,2a,F10.4,2a,3F12.6,2a)')'  - Atoms',ia,ch10,&
&             "    - atomic number:",eff_pot%crystal%znucl(eff_pot%crystal%typat(ia)),ch10,&
&             "    - atomic mass:",eff_pot%crystal%amu(eff_pot%crystal%typat(ia)),ch10,&
&             "    - cartesian position:",eff_pot%crystal%xcart(:,ia),ch10,&
&             "    - Effective charges:"
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
      do ii = 1,3
        write(msg,'(a,3(F12.6))') "  ",eff_pot%harmonics_terms%zeff(:,ii,ia)
        call wrtout(ab_out,msg,'COLL')
        call wrtout(std_out,msg,'COLL')
      end do
    end do
  end if

end subroutine effective_potential_print
!!***

!****f* m_effective_potential/effective_potential_printSupercell
!!
!! NAME
!! effective_potential_printSupercell
!!
!! FUNCTION
!! Print the supercell of the effective_potential
!! or if present the supercell as input
!! WARNING: need to be consistent with eff_pot
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!! supercell<type(supercell_type)> = optional, supercell type to define
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_printSupercell(eff_pot,supercell)

!Arguments ------------------------------------
!scalars
!array
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(supercell_type),optional,target,intent(in) :: supercell
!Local variables-------------------------------
!scalar
 integer :: iatom,ii
 character(len=500) :: msg
!array
 real(dp), allocatable :: xred(:,:)
 type(supercell_type),pointer :: supercell_tmp

! *************************************************************************

 if(present(supercell)) then
   supercell_tmp => supercell
 else
   supercell_tmp => eff_pot%supercell
 end if

 if(supercell_tmp%natom/= eff_pot%supercell%natom) then
   write(msg, '(3a)' )&
&  ' There is not the same numbers of atoms in the two supercell',ch10,&
&   'Action: modify the code'
   MSG_BUG(msg)
 end if

 ABI_ALLOCATE(xred,(3,supercell_tmp%natom))

!**********************************************************************
! Write basics values
!**********************************************************************

 write (msg, '(4a,I8,a)') ' Structure parameters of the supercell :',ch10,ch10,&
                       '  natom ', supercell_tmp%natom,ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write (msg, '(a)') '  znucl '
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')
 write(msg,*) ''
 do iatom = 1, size(eff_pot%crystal%znucl)
   write (msg, '(a,I5)') trim(msg),int(eff_pot%crystal%znucl(iatom))
   if (mod(iatom,6) == 0) then
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     write(msg,*) ''
   end if
 end do
 write (msg, '(2a)') trim(msg),ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')
 write (msg, '(a,I7,3a)') '  ntypat', size(eff_pot%crystal%znucl),ch10,ch10, '  typat '
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,*) ''
 do iatom = 1, supercell_tmp%natom
   write (msg, '(a,I5)') trim(msg),&
&         supercell_tmp%typat(supercell_tmp%atom_indexing(iatom))
   if (mod(iatom,12) == 0)then
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     write(msg,*) ''
   end if
 end do
 write (msg, '(2a)') trim(msg),ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')
 write (msg, '(3a)') '  acell 1.0 1.0 1.0',ch10,ch10
 write (msg, '(2a)') trim(msg),'  rprim'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 do ii = 1,3
   write(msg,'(3E23.14,3E23.14,3E23.14)') supercell_tmp%rprimd(1,ii),&
&                                         supercell_tmp%rprimd(2,ii),&
&                                         supercell_tmp%rprimd(3,ii)
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end do

  write (msg, '(2a)') ch10,'  xcart'
  call wrtout(ab_out,msg,'COLL')
  call wrtout(std_out,msg,'COLL')
  do iatom = 1, supercell_tmp%natom
    write (msg, '(3E23.14)') supercell_tmp%xcart(1,iatom),&
&                                supercell_tmp%xcart(2,iatom),&
&                                supercell_tmp%xcart(3,iatom)
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
  end do
  call xcart2xred(supercell_tmp%natom,supercell_tmp%rprimd,&
 &                supercell_tmp%xcart,xred)
  write (msg, '(2a)') ch10,'  xred'
  call wrtout(ab_out,msg,'COLL')
  call wrtout(std_out,msg,'COLL')
  do iatom = 1, supercell_tmp%natom
    write (msg, '(3E23.14)') xred(1,iatom),&
&                                xred(2,iatom),&
&                                xred(3,iatom)
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
  end do

  write (msg, '(a)') ''
  call wrtout(ab_out,msg,'COLL')
  call wrtout(std_out,msg,'COLL')

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
!! INPUTS
!! filename = the name of output file
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!! option   = 0 Do nothing
!!          = 1 Generate the XML file with:
!!                - The system definition and the model (Harmonic + Anharmonic)
!!          = 2 Generate two XML files with:
!!                - The system definition and the model (Harmonic)
!!                - The model (Anharmonic)
!!          = 3 Generate one XML files with:
!!                - The system definition and the model (Harmonic)
!!          = 4 Generate one XML files with:
!!                - The model (Anharmonic)
!!
!! OUTPUT
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_writeXML(eff_pot,option,filename,prt_dipdip)

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: option
  character(len=fnlen),optional,intent(in) :: filename
  logical,optional,intent(in) :: prt_dipdip
!arrays
  type(effective_potential_type), intent(in) :: eff_pot

!Local variables-------------------------------
!scalar
 integer :: ii,ia,ib,jj
 integer :: iqpt,irpt,mu,nu
 integer :: unit_xml
 character(len=500) :: msg
 character(len=fnlen) :: namefile
 character(len=10) :: natom
 logical :: new_file = .FALSE.,need_prtdipdip = .FALSE.
!arrays
 real(dp) :: strain(9,6)

! *************************************************************************

 strain(:,1) = (/1,0,0,0,0,0,0,0,0/)
 strain(:,2) = (/0,0,0,0,1,0,0,0,0/)
 strain(:,3) = (/0,0,0,0,0,0,0,0,1/)
 strain(:,4) = half*(/0,0,0,0,0,1,0,1,0/)
 strain(:,5) = half*(/0,0,1,0,0,0,1,0,0/)
 strain(:,6) = half*(/0,1,0,1,0,0,0,0,0/)

 unit_xml = get_unit()
 need_prtdipdip = .TRUE.
 if(present(prt_dipdip)) need_prtdipdip = prt_dipdip
!Print only the reference system in xml format
 if (option ==  1 .or. option == 2 .or. option ==3) then

!  convert natom in character
   write (natom,'(I9)') eff_pot%crystal%natom

!  Compute the name of the XML file
   if(present(filename)) then
     namefile=replace(trim(filename),".out","")
     select case(option)
     case(1)
       namefile=trim(namefile)//"_model.xml"
     case(2)
       namefile=trim(namefile)//"_sys.xml"
     case(3)
       namefile=trim(namefile)//"_sys.xml"
     end select
   else
     namefile='system.xml'
   end if

   call isfile(namefile,'new')

   if (open_file(namefile,msg,unit=unit_xml,form="formatted",&
&      status="new",action="write") /= 0) then
     MSG_ERROR(msg)
   end if

   write(msg,'(a,a,a)')ch10,&
 &       ' Generation of the xml file for the model in ',trim(namefile)

   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

!  Write header
   WRITE(unit_xml,'("<?xml version=""1.0"" ?>")')
   WRITE(unit_xml,'("<System_definition>")')

   WRITE(unit_xml,'("  <energy>")')
   WRITE(unit_xml,'(E23.14)') (eff_pot%energy)
   WRITE(unit_xml,'("  </energy>")')

   WRITE(unit_xml,'("  <unit_cell units=""bohrradius"">")')
   do mu=1,3
     WRITE(unit_xml,'(3(E23.14))') (eff_pot%crystal%rprimd(mu,nu),nu=1,3)
   end do
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
     if(any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt))>tol20)) then
       WRITE(unit_xml,'("  <local_force_constant units=""hartree/bohrradius**2"">")')
       WRITE(unit_xml,'("    <data>")')
       do ia=1,eff_pot%crystal%natom
         do mu=1,3
           do ib=1,eff_pot%crystal%natom
             do  nu=1,3
               WRITE(unit_xml,'(e22.14)', advance="no")&
&                         (eff_pot%harmonics_terms%ifcs%short_atmfrc(mu,ia,nu,ib,irpt))
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
     if(need_prtdipdip)then
       if(all(abs(eff_pot%harmonics_terms%ifcs%atmfrc(:,:,:,:,irpt))<tol20)) then
         if(any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt))>tol20)) then
           write(msg, '(a,a,a,a)' )&
&         ' There is no total range but short range in your effective potential',ch10,&
&         'Action: contact abinit group'
           MSG_BUG(msg)
         end if
       else
         WRITE(unit_xml,'("  <total_force_constant units=""hartree/bohrradius**2"">")')
         WRITE(unit_xml,'("    <data>")')
         do ia=1,eff_pot%crystal%natom
           do mu=1,3
             do ib=1,eff_pot%crystal%natom
               do nu=1,3
                 WRITE(unit_xml,'(e22.14)', advance="no")&
&                   (eff_pot%harmonics_terms%ifcs%atmfrc(mu,ia,nu,ib,irpt))
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
&             (eff_pot%harmonics_terms%strain_coupling(ii,mu,ia))
       end do
       WRITE(unit_xml,'(a)')''
     end do
     WRITE(unit_xml,'("    </correction_force>")')
     if (eff_pot%has_anharmonicsTerms)then
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
&                   (eff_pot%anharmonics_terms%elastic_displacement(ii,nu,mu,ia))
             end do
           end do
           WRITE(unit_xml,'(a)')''
         end do
       end if
       WRITE(unit_xml,'("    </correction_strain_force>")')
       if (eff_pot%anharmonics_terms%has_strain_coupling) then
         do irpt=1,eff_pot%anharmonics_terms%phonon_strain(ii)%nrpt
           WRITE(unit_xml,'("    <correction_force_constant units=""hartree/bohrradius**2"">")')
           WRITE(unit_xml,'("      <data>")')
           do ia=1,eff_pot%crystal%natom
             do mu=1,3
               do ib=1,eff_pot%crystal%natom
                 do  nu=1,3
                   WRITE(unit_xml,'(e22.14)', advance="no")&
&                      (eff_pot%anharmonics_terms%phonon_strain(ii)%atmfrc(mu,ia,nu,ib,irpt))
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
       end if!End if has_straincouplitn
     end if!end Hasstrain_coupling
     WRITE(unit_xml,'("    </strain_coupling>")')
     jj = jj + 1
   end do!end mu

   if(option /=1)  WRITE(unit_xml,'("</System_definition>")')
!  Close file
   CLOSE(unit_xml)

 end if!end option

!Print the coefficients into XML file
 if (option==1 .or. option == 2 .or. option==4) then
!   Compute the name of the XML file
   if(present(filename)) then
     namefile=replace(trim(filename),".out","")
     select case(option)
     case(1)
       new_file = .FALSE.
       namefile=trim(namefile)//"_model.xml"
     case(2)
       new_file = .TRUE.
       namefile=trim(namefile)//"_coeffs.xml"
     case(4)
       new_file = .TRUE.
       namefile=trim(namefile)//"_coeffs.xml"
     end select
   else
     namefile='coeffs.xml'
   end if

   if(eff_pot%anharmonics_terms%ncoeff > 0) then
     call polynomial_coeff_writeXML(eff_pot%anharmonics_terms%coefficients,&
&                                   eff_pot%anharmonics_terms%ncoeff,namefile,unit=unit_xml,&
&                                   newfile=new_file)
   end if
 end if!end option

 if(option==1)then
!  add the end of the file in the case option==1
   open(unit=unit_xml,file=namefile,position="append")
   WRITE(unit_xml,'("</System_definition>")')
!  Close file
   CLOSE(unit_xml)
 end if

end subroutine effective_potential_writeXML
!!***

!!****f* m_effective_potential/effective_potential_writeAbiInput
!! NAME
!! effective_potential_writeAbiInput
!!
!! FUNCTION
!! This routine print the effective potential into input of abinit
!! We can also apply a strain to the structure
!!
!! INPUTS
!! eff_pot<type(effective_potential_type)> = effective_potential datatype
!! filename = the name of input file
!! strain<(strain_type)> = optional,strain datatype if need to apply strain into rprim
!!
!! OUTPUT
!!
!! PARENTS
!!      compute_anharmonics
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_writeAbiInput(eff_pot,filename,strain)

!Arguments ------------------------------------
!scalars
  type(strain_type),optional,intent(in) :: strain
  character(len=fnlen),optional,intent(in) :: filename
!arrays
  type(effective_potential_type), intent(in) :: eff_pot
!Local variables-------------------------------
!scalar
 integer :: unit = 20
 character(len=500) :: msg
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

 if (open_file(namefile,msg,unit=unit,form="formatted",status="new",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

  write(msg,'(a,a,a,a)')ch10,&
 &   ' Generation of the input file in ',namefile,ch10
  call wrtout(ab_out,msg,'COLL')
  call wrtout(std_out,msg,'COLL')

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


!****f* m_effective_potential/effective_potential_evaluate
!!
!! NAME
!! effective_potential_evaluate
!!
!! FUNCTION
!! evaluate the energy, forces and stresses
!! of a structure with the effective potential.
!! For the structure you can only provide xred and rprimd,
!! this routine will compute the displacement and the strain.
!! Otherwise, you can provide the strain and the displacement
!!
!! INPUTS
!! eff_pot <type(effective_potential_type)> = effective potential datatype
!! natom = number of atoms in the supercell
!! rprimd(3,3)   = dimensional primitive translations for real space
!!                 of the perturbed structure
!! xred(3,natom) = optional, reduced coordinates of the atoms of the structure
!! strain(6) = optional, strain apply wrt to the reference in eff_pot
!! displacement(3,natom) = optional, atomic displacement wrt to the reference in eff_pot
!! du_delta(6,3,natom) = optional, variation of displacement wrt to the strain
!! compute_anharmonic  = optional, TRUE if the anhamornic part is needed
!! verbose = optional, FALSE for no verbose
!!
!! OUTPUTS
!! energy =  energy of the structure (Ha)
!! fcart(3,natom) =  forces in cartesian coordinates (Ha/Bohr)
!! fred(3,natom)  =  forces in reduced coordinates
!! strten(6) = stress tensor (Ha/Bohr^3)
!!
!! PARENTS
!!      m_effective_potential,m_fit_data,m_fit_polynomial_coeff,mover
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                       displacement,du_delta,strain,xred,&
&                                       compute_anharmonic,verbose)

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom
!array
  type(effective_potential_type),intent(in) :: eff_pot
  real(dp),intent(out) :: energy
  real(dp),intent(out) :: fcart(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  real(dp),intent(out) :: strten(6)
  real(dp),intent(in)  :: rprimd(3,3)
  real(dp),intent(in),optional  :: xred(3,natom)
  real(dp),intent(in),optional :: strain(6)
  real(dp),intent(in),optional :: displacement(3,natom)
  real(dp),intent(in),optional :: du_delta(6,3,natom)
  logical,intent(in),optional :: verbose,compute_anharmonic
!Local variables-------------------------------
!scalar
  integer :: alpha,ii,ia,mu,ncell,comm
  real(dp):: energy_part,ucvol
  logical :: has_strain = .FALSE.,need_verbose
  logical :: iam_master,need_anharmonic
  integer, parameter:: master = 0
!array
  type(strain_type) :: strain_t
  integer :: sc_size(3)
  real(dp) :: disp_tmp(3,natom)
  real(dp) :: du_delta_tmp(6,3,natom)
  real(dp) :: fcart_part(3,natom)
  real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
  real(dp) :: strain_tmp(6),strten_part(6)
  real(dp) :: energy_coeff_part(eff_pot%anharmonics_terms%ncoeff)
  real(dp),allocatable :: xcart(:,:)
  character(len=500) :: msg

! *************************************************************************

  !MPI variables
  comm = eff_pot%mpi_ifc%comm
  iam_master = (eff_pot%mpi_ifc%my_rank == master)

! Set variables
  ncell          = eff_pot%supercell%ncells
  do ii = 1, 3
    sc_size(ii) = eff_pot%supercell%rlatt(ii,ii)
  end do

  need_verbose = .TRUE.
  if(present(verbose)) then
    need_verbose = verbose
  end if

  need_anharmonic = .TRUE.
  if(present(compute_anharmonic))then
    need_anharmonic = compute_anharmonic
  end if

! Check some variables
  if (natom /= eff_pot%supercell%natom) then
    write(msg,'(a,I7,a,I7,a)')' The number of atoms is not correct :',natom,&
&   ' in argument istead of ',eff_pot%supercell%natom, ' in supercell'
    MSG_ERROR(msg)
  end if

  if (present(displacement))then
    if(size(displacement(1,:)) /= eff_pot%supercell%natom) then
      write(msg,'(a,I7,a,I7,a)')' The number of atoms is not correct :',size(displacement(1,:)),&
&      ' in displacement array instead of ',eff_pot%supercell%natom, ' in supercell'
      MSG_ERROR(msg)
    end if
  end if
  if (present(du_delta))then
    if(size(du_delta,3) /= eff_pot%supercell%natom) then
      write(msg,'(a,I7,a,I7,a)')' The number of atoms is not correct :',size(du_delta,3),&
&      ' in du_delta array instead of ',eff_pot%supercell%natom, ' in supercell'
      MSG_ERROR(msg)
    end if
  end if
  do ii=1,3
    if(eff_pot%supercell%rlatt(ii,ii)<0.or.eff_pot%supercell%rlatt(ii,ii)>50)then
      write(msg, '(a,i0,a,i2,a,a,a,i0,a)' )&
&     'eff_pot%supercell%rlatt(',ii,') is ',int(eff_pot%supercell%rlatt(ii,ii)),&
&     ', which is lower than 0 of superior than 10.',ch10,'Action: correct ncell(',ii,').'
      MSG_ERROR(msg)
    end if
  end do

  if(need_verbose)then
    write(msg, '(a,a,a)' ) ch10,' Calculation of the energy with effective potential'
    call wrtout(ab_out,msg,'COLL')
  end if

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!--------------------------------------------
! 1 - Set the perturbations and intialisation
!--------------------------------------------

! Get strain
  strain_tmp(:) = zero
  if (present(strain)) then
    strain_tmp(:) = strain(:)
    call strain_get(strain_t,mat_delta=reshape((/&
&                                      strain_tmp(1),strain_tmp(6)/two,strain_tmp(5)/two,&
&                                      strain_tmp(6)/two,strain_tmp(2),strain_tmp(4)/two,&
&                                      strain_tmp(5)/two,strain_tmp(4)/two,strain_tmp(3)/),(/3,3/)))
    has_strain = .TRUE.
  else
!   Compute the strain
    call strain_get(strain_t,rprim=eff_pot%supercell%rprimd,&
&                   rprim_def=rprimd,symmetrized=.FALSE.)
    if(need_verbose)then
      write(msg,'(80a)')('-',mu=1,80)
      call wrtout(std_out,msg,'COLL')
      call strain_print(strain_t)
    end if
    if (strain_t%name /= "reference")  then
      has_strain = .TRUE.
      do ii=1,3
        strain_tmp(ii) = strain_t%strain(ii,ii)
      end do
      strain_tmp(4) = (strain_t%strain(2,3) + strain_t%strain(3,2))
      strain_tmp(5) = (strain_t%strain(3,1) + strain_t%strain(1,3))
      strain_tmp(6) = (strain_t%strain(2,1) + strain_t%strain(1,2))
    else
      strain_tmp(:) = zero
    end if
  end if
! Get displacement and the variation of the displacmeent wr to strain
  ABI_ALLOCATE(xcart,(3,natom))
  disp_tmp(:,:) = zero
  du_delta_tmp(:,:,:) = zero
  if((.not.present(displacement).or..not.present(du_delta)).and.present(xred))then
!   Compute the displacement
    call xred2xcart(natom, rprimd, xcart, xred)
    call effective_potential_getDisp(disp_tmp,du_delta_tmp,natom,rprimd,&
&                                    eff_pot%supercell%rprimd,comm,xcart_hist=xcart,&
&                                    xcart_ref=eff_pot%supercell%xcart,&
&                                    compute_displacement = .not.present(displacement),&
&                                    compute_duDelta =  .not.present(du_delta))
  end if

! or set ftrom the arguments
  if(present(displacement)) disp_tmp(:,:)  = displacement(:,:)
  if(present(du_delta))du_delta_tmp(:,:,:) = du_delta(:,:,:)
!Set to zero the outputs
  energy         = zero
  fcart(:,:)     = zero
  strten(:)      = zero

  if(need_verbose)then
    write(msg, '(80a,2a)' ) ('-',mu=1,80),&
&       ch10,' Components of total energy (in Hartree) :'
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
  end if

!------------------------------------
! 2 - Transfert the reference values
!------------------------------------

! Set the value of the energy
  energy = eff_pot%energy * ncell

  if(need_verbose)then
    write(msg, '(a,a,1ES24.16,a)' ) ch10,' Energy of the reference structure         :',&
&                                          energy,' Hartree'
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
  end if

!Set the value of the initial strees (should be zero if the system is relaxed)
!According to the original definition of the model, the reference should always be
! a critical point of the PES so the first order derivative should not be take into account...
!  do ii=1,6
!   energy = energy + ncell * eff_pot%strten(ii)*strain_tmp(ii)
! end do
! strten(:) = ncell * eff_pot%strten(:)
!Set the value of the initial forces (should be zero if the system is relaxed)
! ii = 1
! This next part is not working, need to apply strain to the original fcart
! do ia=1,eff_pot%supercell%natom
!   fcart(:,ia) = eff_pot%fcart(:,ii)
!   do mu=1,3
!     energy = energy + eff_pot%fcart(mu,ii)*disp_tmp(mu,ii)
!   end do
!   ii = ii + 1
!   if(ii > eff_pot%crystal%natom) ii = 1
! end do

!------------------------------------
! 3 - Computation of the IFC part :
!------------------------------------

  energy_part    = zero
  fcart_part(:,:)= zero

  call harmonics_terms_evaluateIFC(eff_pot%harmonics_terms%ifcs%atmfrc(:,:,:,:,:),disp_tmp,&
&                                  energy_part,fcart_part,eff_pot%supercell%natom,&
&                                  eff_pot%crystal%natom,eff_pot%mpi_ifc%my_ncell,&
&                                  eff_pot%mpi_ifc%my_nrpt,eff_pot%mpi_ifc%my_atmrpt_index,&
&                                  eff_pot%mpi_ifc%my_index_cells,sc_size,eff_pot%mpi_ifc%my_rpt,&
&                                  eff_pot%mpi_ifc%comm)

  if(need_verbose)then
    write(msg, '(a,1ES24.16,a)' ) ' Energy of the ifc part                    :',&
&                                     energy_part,' Hartree'
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')

    if(abs(energy_part) < tol16 .and. eff_pot%anharmonics_terms%ncoeff == 0)then
      write(msg, '(8a)' )ch10,&
&        ' --- !WARNING!',ch10,&
&        '        The harmonic part is negative, the simulation will diverge',ch10,&
&        '        if the anharmonic part is not used',ch10,&
&        ' ---'
      call wrtout(std_out,msg,"COLL")
    end if
  end if

  energy = energy + energy_part
  fcart(:,:)= fcart(:,:) + fcart_part(:,:)


!----------------------------------------------------
! 4 - Computation of the elastic part of the energy :
!----------------------------------------------------

  energy_part    = zero
  fcart_part(:,:)= zero
  strten_part(:) = zero

  call harmonics_terms_evaluateElastic(eff_pot%harmonics_terms%elastic_constants,disp_tmp,energy_part,&
&                                fcart_part,eff_pot%supercell%natom,eff_pot%crystal%natom,&
&                                ncell,eff_pot%harmonics_terms%strain_coupling,strten_part,strain_tmp)


  if(has_strain.and.need_verbose)then
    write(msg, '(a,1ES24.16,a)' ) ' Energy of the elastic part                :',&
&                                       energy_part,' Hartree'
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
  end if

  energy = energy + energy_part
  fcart(:,:) = fcart(:,:)  + fcart_part(:,:)
  strten(:) = strten(:) + strten_part(:)

!------------------------------------
! 5 - Treat 3rd order strain-coupling:
!------------------------------------
  if (need_anharmonic.and.eff_pot%has_anharmonicsTerms) then
    energy_part    = zero
    strten_part(:) = zero
    fcart_part(:,:)= zero

!   1-Part due to the anharmonic strain
    if(eff_pot%anharmonics_terms%has_elastic3rd.or.&
&      eff_pot%anharmonics_terms%has_elastic4th.or.&
&      eff_pot%anharmonics_terms%has_elastic_displ)then
      call anharmonics_terms_evaluateElastic(disp_tmp,energy_part,fcart_part,&
&                                    eff_pot%supercell%natom,&
&                                   eff_pot%crystal%natom,ncell,strten_part,strain_tmp,&
&                                   elastic3rd=eff_pot%anharmonics_terms%elastic3rd,&
&                                   elastic4th=eff_pot%anharmonics_terms%elastic4th,&
&                                   elastic_displacement=eff_pot%anharmonics_terms%elastic_displacement)

      if(need_verbose)then
        write(msg, '(a,1ES24.16,a)' ) ' Energy of the anharmonic elastic part     :',&
&                                          energy_part,' Hartree'
        call wrtout(ab_out,msg,'COLL')
        call wrtout(std_out,msg,'COLL')
      end if
      energy = energy + energy_part
      strten(:) = strten(:) + strten_part(:)
    end if

!   2-Part due to the strain-phonon coupling
    if (eff_pot%anharmonics_terms%has_strain_coupling) then
      energy_part    = zero
      strten_part(:) = zero
      fcart_part(:,:)= zero

      call anharmonics_terms_evaluateIFCStrainCoupling(eff_pot%anharmonics_terms%phonon_strain,&
&                                                      disp_tmp,energy_part,fcart_part,&
&                                                      eff_pot%supercell%natom,&
&                                                      eff_pot%crystal%natom,&
&                                                      sc_size,strain_tmp,&
&                                                      strten_part,eff_pot%mpi_ifc%my_cells,&
&                                                      eff_pot%mpi_ifc%my_ncell,&
&                                                      eff_pot%mpi_ifc%my_index_cells,&
&                                                      eff_pot%mpi_ifc%comm)

      if(need_verbose)then
        write(msg, '(a,1ES24.16,a)' ) ' Energy of the 3rd (strain-phonon coupling):',&
&                                         energy_part,' Hartree'
        call wrtout(ab_out,msg,'COLL')
        call wrtout(std_out,msg,'COLL')
      end if
      energy = energy + energy_part
      fcart  = fcart  + fcart_part
      strten = strten + strten_part
    end if

  end if

!----------------------------------
! 6 - Treat polynomial coefficient:
!----------------------------------
  if(need_anharmonic.and.eff_pot%anharmonics_terms%ncoeff > zero)then
    energy_part = zero
    fcart_part(:,:)  = zero
    strten_part(:) = zero
    call polynomial_coeff_evaluate(eff_pot%anharmonics_terms%coefficients,disp_tmp,&
&                                  energy_part,energy_coeff_part,fcart_part,eff_pot%supercell%natom,&
&                                  eff_pot%crystal%natom,eff_pot%anharmonics_terms%ncoeff,&
&                                  sc_size,strain_tmp,strten_part,eff_pot%mpi_coeff%my_ncell,&
&                                  eff_pot%mpi_coeff%my_index_cells,eff_pot%mpi_coeff%comm)

    if(need_verbose)then
      write(msg, '(a,1ES24.16,a)' ) ' Energy of the fitted coefficient          :',&
&                                       energy_part,' Hartree'
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
    end if
    energy = energy + energy_part
    fcart(:,:) = fcart(:,:) + fcart_part(:,:)
    strten(:) = strten(:) + strten_part(:)
  end if

!---------------------------------
! 7 - Compute confinement
!---------------------------------
  if(eff_pot%confinement%need_confinement) then

    energy_part = zero

    call polynomial_conf_evaluate(disp_tmp,eff_pot%confinement%cutoff_disp,energy_part,&
&                                 eff_pot%confinement%factor_disp,&
&                                 eff_pot%confinement%factor_strain,fcart_part,strain_tmp,&
&                                 eff_pot%confinement%cutoff_strain,strten_part,&
&                                 eff_pot%confinement%power_disp,eff_pot%confinement%power_strain,&
&                                 eff_pot%mpi_coeff%my_cells,&
&                                 eff_pot%supercell%natom,eff_pot%crystal%natom,&
&                                 eff_pot%mpi_coeff%my_ncell,eff_pot%mpi_coeff%my_index_cells,&
&                                 eff_pot%mpi_coeff%comm)

      energy = energy + energy_part

    if(abs(energy_part) > tol10 .and. need_verbose )then
      write(msg, '(a,1ES24.16,a)' ) ' Energy of the confinement part            :',&
&                                       energy_part,' Hartree'
      call wrtout(ab_out,msg,'COLL')
      call wrtout(std_out,msg,'COLL')
    end if
  end if

!-----------------------------------
! 8 - Add variation of the atomic
!     displacement due to the strain
!-----------------------------------
  if (has_strain) then
    strten_part(:) = zero
    do ia = 1,eff_pot%supercell%natom
      do mu = 1,3
        do alpha=1,6
          strten_part(alpha) = strten_part(alpha) + fcart(mu,ia) * du_delta_tmp(alpha,mu,ia)
        end do
      end do
    end do
    strten(:) = strten(:) + strten_part(:)
  end if

!---------------------------------
! 8 - Apply factors
!---------------------------------

! divide stess tensor by ucvol
  strten = strten / ucvol

! multiply forces by -1
  fcart = -1 * fcart

! Redistribute the residuale of the forces
  call effective_potential_distributeResidualForces(eff_pot,fcart,eff_pot%supercell%natom)

  call fcart2fred(fcart,fred,rprimd,natom)

!------------------------------------
! 9 - Final Print:
!------------------------------------

  if(need_verbose)then
    write(msg, '(2a,es21.14)' ) ch10,&
&     '    >>>>>>>>> Etotal= ',energy
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')

    write(msg,'(2a,1p,e15.7,a)') ch10,' Unit cell volume ucvol=',ucvol+tol10,' bohr^3'
    call wrtout(std_out,  msg,'COLL')

    write(msg, '(a,80a,3a)' ) ch10,('-',mu=1,80),ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')

    write(msg, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')
    write(msg, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')
    write(msg, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')

! Also output the pressure (minus one third the trace of the stress
! tensor.
    write(msg, '(a,a,es12.4,a)' ) ch10,&
&   '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&   -(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa]'

    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')

    write(msg, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(1 1)=',strten(1)*HaBohr3_GPa,&
&   '  sigma(3 2)=',strten(4)*HaBohr3_GPa
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')
    write(msg, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(2 2)=',strten(2)*HaBohr3_GPa,&
&   '  sigma(3 1)=',strten(5)*HaBohr3_GPa
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')
    write(msg, '(a,1p,e16.8,a,1p,e16.8,a)' ) &
&   '- sigma(3 3)=',strten(3)*HaBohr3_GPa,&
&   '  sigma(2 1)=',strten(6)*HaBohr3_GPa
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,  msg,'COLL')
    write(msg, '(80a,a)' ) ('-',mu=1,80),ch10
    call wrtout(ab_out,msg,'COLL')
    call wrtout(std_out,msg,'COLL')
  end if



  ABI_DEALLOCATE(xcart)

end subroutine effective_potential_evaluate
!!***

!****f* m_effective_potential/effective_potential_getDisp
!!
!! NAME
!! effective_potential_getDisp
!!
!! FUNCTION
!! Compute cartesian coordinates of the displacment
!! between two configurations.
!! You can provide reduce or cartesian coordinates for both perturbed structure
!! and reference structure
!!
!! INPUTS
!! natom = number of atoms
!! rprimd_hist = dimensional primitive translations for real space of the perturbed structure
!! rprimd_ref  = dimensional primitive translations for real space of the reference structure
!! comm = mpi communicator
!! xcart_hist(3,natom) = optional, cartesian coordinates of the atoms in the perturbed structure
!! xred_hist(3,natom)  = optional, reduced coordinates of the atoms in the perturbed structure
!! xred_ref(3,natom)   = optional, reduced coordinates of the atoms in the reference structure
!! xccart_ref(3,natom) = optional, cartesian coordinates of the atoms in the reference structure
!! compute_displacement= optional, flag to compute the displacement array (default is true)
!! compute_duDelta = optional, flag to compute the du_delta array (default is true)
!!                             Be careful, if you specify compute_displacement=.false.,
!!                             you need to provide to this routine the good displacement array!
!!
!!
!! OUTPUT
!! displacement(3,natom) = cartesian atomic displacement between two configurations (Bohr)
!! du_delta(6,3,natom)   = variation of the displacmeent wr to strain
!!
!! PARENTS
!!      m_effective_potential,m_fit_data
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_getDisp(displacement,du_delta,natom,rprimd_hist,rprimd_ref,comm,&
&                                      xcart_hist,xred_hist,xred_ref,xcart_ref,compute_displacement,&
&                                      compute_duDelta)

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: natom,comm
  logical,optional,intent(in) :: compute_displacement,compute_duDelta
!array
  real(dp),intent(in) :: rprimd_ref(3,3),rprimd_hist(3,3)
  real(dp),intent(out) :: displacement(3,natom),du_delta(6,3,natom)
  real(dp),intent(in),optional :: xred_hist(3,natom),xcart_hist(3,natom)
  real(dp),intent(in),optional :: xred_ref(3,natom),xcart_ref(3,natom)
!Local variables-------------------------------
!scalar
  integer :: ii,ia,ib,mu
  integer :: ierr,nproc,my_rank,natom_alone,my_natom
  character(len=500) :: msg
  logical :: has_strain = .FALSE.
  logical :: need_displacement,need_duDelta
!array
  integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
  integer,allocatable :: my_atoms(:)
  type(strain_type) :: strain
  real(dp) :: xcart_hist_tmp(3,natom),xcart_ref_tmp(3,natom)
  real(dp) :: xred_ref_tmp(3,natom)
! *************************************************************************

  if (.not.(present(xred_ref).or.present(xcart_ref))) then
     write(msg, '(3a)' )&
&         'You need at least give xcart_ref or xred_ref '
     MSG_ERROR(msg)
  end if

  if (.not.(present(xred_hist).or.present(xcart_hist))) then
     write(msg, '(3a)' )&
&         'You need at least give xcart_hist or xred_hist '
     MSG_ERROR(msg)
  end if

  need_duDelta = .TRUE.
  need_displacement = .TRUE.
  if(present(compute_duDelta)) need_duDelta = compute_duDelta
  if(present(compute_displacement)) need_displacement = compute_displacement
  if(.not.need_duDelta .and. .not.need_displacement) return

!--------------------------------------------
! 0 - Set the MPI
!--------------------------------------------
  nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
  natom_alone = mod(natom,nproc)
  my_natom = int(aint(real(natom,sp)/(nproc)))
  if(my_rank >= (nproc-natom_alone)) then
    my_natom = my_natom  + 1
  end if
  ABI_ALLOCATE(my_atoms,(my_natom))
  my_atoms = 0
  do ii=1,my_natom
    if(my_rank >= (nproc-natom_alone))then
      my_atoms(ii)=(int(aint(real(natom,sp)/nproc)))*(my_rank)+&
&                              (my_rank - (nproc-natom_alone)) + ii
    else
      my_atoms(ii)=(my_natom)*(my_rank)  + ii
    end if
  end do

!--------------------------------------------
! 1 - Get the strain for this step
!--------------------------------------------

  call strain_get(strain,rprim=rprimd_ref,rprim_def=rprimd_hist)
  if (strain%name /= "reference")  then
    has_strain = .TRUE.
  end if

! fill the history position
  if(present(xcart_hist)) then
    xcart_hist_tmp(:,:) = xcart_hist(:,:)
  else
    call xred2xcart(natom, rprimd_hist, xcart_hist_tmp, xred_hist)
  end if

! Fill the reference position and change the cartesian coordinates
! if the rprimd is different
  if(has_strain) then
    if(present(xcart_ref)) then
      call xcart2xred(natom, rprimd_ref,  xcart_ref,     xred_ref_tmp)
      call xred2xcart(natom, rprimd_hist, xcart_ref_tmp, xred_ref_tmp)
    else
      call xred2xcart(natom, rprimd_hist, xcart_ref_tmp, xred_ref)
    end if
  else
    if(present(xcart_ref)) then
      xcart_ref_tmp(:,:) = xcart_ref(:,:)
    else
      call xred2xcart(natom, rprimd_ref, xcart_ref_tmp, xred_ref)
    end if
  end if

! Compute displacement
  if(need_displacement)then
    displacement(:,:) = zero
    do ii = 1, natom
      displacement(:,ii) = xcart_hist_tmp(:,ii) - xcart_ref_tmp(:,ii)
    end do
  end if

! Get also the variation of the displacmeent wr to strain
  if(has_strain.and.need_duDelta) then
    du_delta(:,:,:)   = zero
!   Compute displacmeent wr to strain
!   within the metric tensor formulation
!   \frac{\partial u_{k}}{\partial \epsilon_{\alpha\beta}}=
!    \frac{1}{2} (\delta_{\alpha,k} u_{\beta}+ \delta_{\beta,k} u_{\alpha})
    do ia=1,my_natom
      ib = my_atoms(ia)
      do ii=1,6
        do mu=1,3
          if(alpha(ii)==mu)then
            du_delta(ii,mu,ib) = du_delta(ii,mu,ib) + half * displacement(beta(ii),ib)
          end if
          if(beta(ii)==mu)then
            du_delta(ii,mu,ib) = du_delta(ii,mu,ib) + half * displacement(alpha(ii),ib)
          end if
        end do
      end do
    end do
    call xmpi_sum(du_delta , comm, ierr)
  end if

  ABI_DEALLOCATE(my_atoms)

end subroutine effective_potential_getDisp
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
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_distributeResidualForces(eff_pot,fcart,natom)

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
    sum_mass = sum_mass + eff_pot%crystal%amu(eff_pot%supercell%typat(ia))
  end do

  do ia=1,natom
    mass_ia = eff_pot%crystal%amu(eff_pot%supercell%typat(ia))
    fcart(:,ia) = fcart(:,ia) - (mass_ia/sum_mass) * sum_f(:)
  end do


end subroutine effective_potential_distributeResidualForces
!!***

!AM_EXPERIMENTAL SECTION
!!****f* m_effective_potential/equal
!! NAME
!!  equal
!!
!! FUNCTION
!! compare two effective potential
!!
!! INPUTS
!! e1<type(effective_potential_type)> = effective_potential datatype
!! e2<type(effective_potential_type)> = effective_potential datatype
!!
!! OUTPUT
!!
!! SOURCE

pure function effective_potential_compare(e1,e2) result (res)

!Arguments ------------------------------------
  type(effective_potential_type), intent(in) :: e1,e2
  logical :: res
! *************************************************************************
  res = .false.
  if(e1%crystal%natom==e2%crystal%natom.and.&
&     e1%harmonics_terms%ifcs%nrpt==e2%harmonics_terms%ifcs%nrpt.and.&
&     e1%crystal%ntypat==e2%crystal%ntypat.and.&
&     e1%harmonics_terms%nqpt==e2%harmonics_terms%nqpt.and.&
&     abs(e1%energy-e2%energy)<tol16.and.&
&     abs(e1%crystal%ucvol-e2%crystal%ucvol)<tol16) then
    res = .true.
  end if

end function effective_potential_compare
!!***

! !****f* m_effective_potential/effective_potential_effpot2ddb
! !!
! !! NAME
! !! effective_potential_effpot2ddb
! !!
! !! FUNCTION
! !! Convert eff_pot into ddb datatype
! !!
! !! INPUTS
! !! eff_pot = effective potential datatype
! !!
! !! OUTPUT
! !! ddb   = ddb with all information
! !!
! !! PARENTS
! !!      m_effective_potential
! !!
! !! CHILDREN
! !!      ab_define_var,isfile,wrtout
! !!
! !! SOURCE

! subroutine effective_potential_effpot2ddb(ddb,crystal,eff_pot,ncell,nph1l,option,qph1l)


! !Arguments ------------------------------------
! !scalars
!   integer,intent(in) :: nph1l,option
! !array
!   integer,intent(in) :: ncell(3)
!   real(dp),intent(in):: qph1l(3,nph1l)
!   type(effective_potential_type),intent(inout) :: eff_pot
!   type(ddb_type),intent(out) :: ddb
!   type(crystal_t),intent(out) :: crystal
! !Local variables-------------------------------
! !scalar
!   integer :: ii,jj,msym
!   real(dp):: ucvol

! ! type(anaddb_dataset_type) :: inp
! !array
!   real(dp) :: gmet(3,3),rmet(3,3)
!   real(dp) :: gprimd(3,3),rprimd(3,3)
!   real(dp),allocatable :: xred(:,:)
! !  character :: title(eff_pot%crystal%ntypat)
!   integer,allocatable :: symrel(:,:,:),symafm(:)
!   real(dp),allocatable :: tnons(:,:)

! ! *************************************************************************

!   ! Number of 2dte blocks in present object
! !  integer,allocatable :: flg(:,:)
!   ! flg(msize,nblok)
!   ! flag to indicate presence of a given block
! !  integer,allocatable :: typ(:)
!   ! typ(nblok)
!   ! type of each block - ddk, dde, phonon etc...
! !  real(dp),allocatable :: amu(:)
!   ! amu(ntypat)
!   ! mass of the atoms (atomic mass unit)
! !  real(dp),allocatable :: nrm(:,:)
!   ! nrm(3,nblok)
!   ! norm of the q-points for each block - can be 0 to indicate a direction of approach to gamma
! !  real(dp),allocatable :: qpt(:,:)
!   ! qpt(9,nblok)
!   ! q-point vector in reciprocal space (reduced lattice coordinates) for each block
! !  real(dp),allocatable :: val(:,:,:)
!   ! val(2,msize,nblok)
!   ! values of the second energy derivatives in each block

! ! Useless value
!    ddb%nblok = -1

!   !option = 1 just print ddb for 1 1 1 cell
!   if(option==1) then
! !   Compute different matrices in real and reciprocal space, also
! !   checks whether ucvol is positive.
!     call metric(gmet,gprimd,-1,rmet,eff_pot%crystal%rprimd,ucvol)

! !   Convert to rprimd
!     do ii=1,3
!       do jj=1,3
!         rprimd(ii,jj)=eff_pot%crystal%rprimd(ii,jj)
!       end do
!     end do

! !   Obtain reciprocal space primitive transl g from inverse trans of r
! !   (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
!     call matr3inv(rprimd,gprimd)

! !   transfert basic values
!     ddb%natom  = eff_pot%crystal%natom
!     ddb%mpert  = ddb%natom+6
!     ddb%msize  = 3*ddb%mpert*3*ddb%mpert;
!     ddb%ntypat = eff_pot%crystal%ntypat
!     ddb%occopt = 3 ! default value
!     ddb%prtvol = 0 ! default value
!     ddb%rprim  = rprimd ! dimensioless real space primitive vectors
!     ddb%gprim  = gprimd ! dimensioless reciprocal space primitive vectors
!     ddb%acell  = one
!     msym   = 1
! !  Setup crystal type
!     ABI_ALLOCATE(xred,(3,ddb%natom))
! !    call xcar2xred(ddb%natom,eff_pot%crystal%rprimd,eff_pot%crystal%xcart,xred)
! !Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
!     ABI_ALLOCATE(symrel,(3,3,msym))
!     ABI_ALLOCATE(symafm,(msym))
!     ABI_ALLOCATE(tnons,(3,msym))

! !    call crystal_init(ddb%amu,crystal,1,ddb%natom,size(eff_pot%crystal%znucl),eff_pot%crystal%ntypat,1,&
! !&       eff_pot%crystal%rprimd,eff_pot%crystal%typat,xred,eff_pot%crystal%znucl,&
! !&       eff_pot%crystal%znucl,0,.FALSE.,.FALSE.,title)!,&
! !&       symrel=symrel,tnons=tnons,symafm=symafm)
! !    call crystal_print(crystal)
! !    stop
! !TEST_AM
!     ABI_DEALLOCATE(symrel)
!     ABI_DEALLOCATE(symafm)
!     ABI_DEALLOCATE(tnons)

!     ABI_DEALLOCATE(xred)

!    else  if (option==2) then
! !   Compute different matrices in real and reciprocal space, also
! !   checks whether ucvol is positive.
!     call metric(gmet,gprimd,-1,rmet,eff_pot%supercell%rprimd,ucvol)

! !   Convert to rprim (dimensionless)
!     do ii=1,3
!       do jj=1,3
!         rprimd(ii,jj)=eff_pot%supercell%rprimd(ii,jj)
!       end do
!     end do

! !   Obtain reciprocal space primitive transl g from inverse trans of r
! !   (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
!     call matr3inv(rprimd,gprimd)

! !   transfert basic values
!     ddb%natom  = eff_pot%supercell%natom
!     ddb%ntypat = eff_pot%crystal%ntypat
!     ddb%mpert  = ddb%natom+6
!     ddb%msize  = 3*ddb%mpert*3*ddb%mpert;
!     ddb%occopt = 3 ! default value
!     ddb%prtvol = 0 ! default value
!     ddb%rprim  = rprimd ! dimensioless real space primitive vectors
!     ddb%gprim  = gprimd ! dimensioless reciprocal space primitive vectors
!     ddb%acell  = one

!    end if
! !TEST_AM
!     !write(std_out,*)"natom ",ddb%natom
!     !write(std_out,*)"ntypat",ddb%ntypat
!     !write(std_out,*)"mpert",ddb%mpert
!     !write(std_out,*)"msize",ddb%msize
!     !write(std_out,*)"occopt",ddb%occopt
!     !write(std_out,*)"prtvol",ddb%prtvol
!     !write(std_out,*)"rprim",ddb%rprim
!     !write(std_out,*)"gprim",ddb%gprim
!     !write(std_out,*)"acell",ddb%acell
! !TEST_AM

!  end subroutine effective_potential_effpot2ddb
! !!***


! !****f* m_effective_potential/effective_potential_printPDOS
! !!
! !! NAME
! !! effective_potential_printPDOS
! !!
! !! FUNCTION
! !! Apply the acoustic sum rule on the effective potential
! !!
! !! INPUTS
! !! eff_pot = effective potential datatype
! !! option  = 0 (default) do nothing
! !!         = 1 print PHFRQ for specific qgrid (need nph1l and qph1l)
! !!         = 2 print PHFRQ for supercell (q=gamma) (need nph1l and qph1l and ncell)
! !! OUTPUT
! !! eff_pot
! !!
! !! PARENTS
! !!
! !! CHILDREN
! !!      ab_define_var,isfile,wrtout
! !!
! !! SOURCE

! subroutine effective_potential_printPDOS(eff_pot,filename,ncell,nph1l,option,qph1l)

! !Arguments ------------------------------------
! !scalars
!   integer,intent(in) :: nph1l,option
! !array
!   integer,intent(in) :: ncell(3)
!   real(dp),intent(in):: qph1l(3,nph1l)
!   type(effective_potential_type),intent(inout) :: eff_pot
!   character(len=fnlen),intent(in) :: filename
! !Local variables-------------------------------
! !scalar
! ! integer :: lenstr
! ! character(len=strlen) :: string
! !array
!  type(crystal_t) :: Crystal
! ! type(anaddb_dataset_type) :: inp
!  type(ddb_type) :: ddb
! ! type(asrq0_t) :: asrq0

! ! *************************************************************************

!  if (option > 0) then

! !  First: transfer into ddb datatype:
!    call effective_potential_effpot2ddb(ddb,Crystal,eff_pot,ncell,nph1l,option,qph1l)

! !  Setup fake anaddb_dataset
! !   string = ''
! !   lenstr = 0
! !   call invars9(inp,lenstr,ddb%natom,string)
! !  fill it with multibinit_dataset values
!  !   inp%prt_ifc = 1
!  !   inp%ifcflag = 1
!  !   inp%qph1l   = qph1l
!  !   inp%nph1l   = nph1l

!  !   ! In case the interatomic forces are not calculated, the
!  !   ! ASR-correction (asrq0%d2asr) has to be determined here from the Dynamical matrix at Gamma.
!  !   if (inp%ifcflag == 0) then
!  !     asrq0 = ddb_get_asrq0(ddb, inp%asr, inp%rfmeth, crystal%xcart)
!  !   end if

!  !  !MG: Note that I'm passing xmpi_comm_self here.
!  !  call mkphbs(eff_pot%harmonics_terms%ifcs,Crystal,inp,ddb,asrq0,filename,xmpi_comm_self)

!  !  call asrq0_free(asrq0)

!   end if

!  end subroutine effective_potential_printPDOS
! !!***

!****f* m_effective_potential/effective_potential_computeGradient
!!
!! NAME
!! effective_potential_computeGradient
!!
!! FUNCTION
!! Compute finate differences on forces to compute dynmical matrix
!! at gamma for supercell
!!
!! INPUTS
!! eff_pot = effective potential datatype
!!
!! OUTPUT
!! dynmat   = ddb with all information
!!
!! PARENTS
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_computeGradient(delta,fcart_out,eff_pot,natom,ncell,option,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option,comm
 real(dp),intent(in) :: delta
!array
 integer,intent(in) :: ncell(3)
 type(effective_potential_type),intent(inout) :: eff_pot
 real(dp),intent(out) :: fcart_out(3,natom)
!Local variables-------------------------------
!scalar
 character(len=500) :: msg
 integer :: ia,ib,ii,mu,nu,npt
 real(dp):: delt,energy
!array
 real(dp):: strten(6)
 real(dp),allocatable :: disp(:,:),diff(:)
 real(dp),allocatable :: fred(:,:),fcart(:,:),xred(:,:)

! *************************************************************************

!Do Some check
 do ii=1,3
   if(abs(eff_pot%supercell%qphon(ii)-ncell(ii)) > tol16)then
     call effective_potential_setSupercell(eff_pot,comm,ncell)
   end if
 end do

  write(msg,'(a,(80a),3a)') ch10,('-',ii=1,80),ch10,' Generation of the dynmical matrix by ',&
&  'finite differences'
!  call wrtout(ab_out,msg,'COLL')
!  call wrtout(std_out,msg,'COLL')

  npt = 3 !default
  select case (option)
  case (1)
!    write(msg,'(2a)') ch10,' Finite differences on 1 points  '
    npt = 2
  case (2)
!    write(msg,'(2a)') ch10,' Finite differences on 3 points  '
    npt = 3
  case (3)
!    write(msg,'(2a)') ch10,' Finite differences on 5 points  '
    npt = 5
  end select

!  call wrtout(ab_out,msg,'COLL')
!  call wrtout(std_out,msg,'COLL')

! Allocation of forces arrays

 ABI_ALLOCATE(disp,(3,natom))
 ABI_ALLOCATE(diff,(npt))
 ABI_ALLOCATE(fred,(3,natom))
 ABI_ALLOCATE(fcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))

 fcart_out = zero

 call xcart2xred(eff_pot%supercell%natom,eff_pot%supercell%rprimd,&
&                eff_pot%supercell%xcart,xred)

 do ia=1,eff_pot%supercell%natom
   do mu=1,3
     diff = zero
     do ii=1,npt
       delt = (-(npt/2+1)+ii) * delta
       disp = zero
       disp(mu,ia) = delt * eff_pot%supercell%rprimd(mu,mu)
       call effective_potential_evaluate(eff_pot,energy,fcart,fred,&
&                                        strten,natom,eff_pot%supercell%rprimd,&
&                                        displacement=disp,&
&                                        compute_anharmonic=.FALSE.,verbose=.false.)

       !       diff(ii,:,:) = fred(:,:)
       diff(ii) = energy
     end do

     select case (option)
     case (1)
       fcart_out(mu,ia) = (diff(1)-diff(2)) / (delta)
     case (2)
       fcart_out(mu,ia) = (diff(3)-diff(1)) / (2*delta)
     case (3)
       fcart_out(mu,ia) = (-diff(5)+8*diff(4)-8*diff(2)+diff(1)) / (12*delta)
     end select
   end do
 end do

!TEST_AM
!Write the phonon  into ddb format wavevector
!write(999, '(a,3es16.8,f6.1)' )' qpt',real((/0,0,0/),dp),1.0
!Write the matrix elements
 do ib=1,eff_pot%supercell%natom
   do nu=1,3
     do ia=1,eff_pot%supercell%natom
       do mu=1,3
!          write(999,'(4i4,2d22.14)')nu,ib,mu,ia,&
! &             dynmat(1,mu,ia,nu,ib),dynmat(2,mu,ia,nu,ib)
       end do
     end do
   end do
 end do
!TEST_AM

! Deallocation of arrays
 ABI_DEALLOCATE(disp)
 ABI_DEALLOCATE(diff)
 ABI_DEALLOCATE(fred)
 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(xred)


 end subroutine effective_potential_computeGradient
!!***

!!****f* m_effective_potential/effective_potential_checkDEV
!!
!! NAME
!! effective_potential_checkDEV
!!
!! FUNCTION
!! Routine for develloper Check by finite differences the equations in
!! effective_potential_evaluate need to provide HIST file, so you need to
!! activate the fit_process or bound_process to activate the reading of the HIST
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = The history of the MD
!! natom = number of atom
!! ntime = number of time in the hist
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

 subroutine effective_potential_checkDEV(eff_pot,hist,natom,ntime)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,ntime
!arrays
 type(effective_potential_type),intent(in) :: eff_pot
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalar
 integer :: ii,jj,ia,mu,npt,istep
! integer :: ifirst
 real(dp):: energy,delt,delta,ucvol
 !arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),mat_def(3,3),identity(3,3)
 real(dp):: fcart(3,natom),fred(3,natom),strten(6),rprimd(3,3)
 real(dp):: rprimd_def(3,3),rprimd_ref(3,3),deltalist(5)
 real(dp):: disp(3,natom),disp_red(3,natom),strain(6),du_delta(6,3,natom),diff(5)
 real(dp),allocatable :: xred(:,:)
 integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
 character(len=500) :: msg

! *************************************************************************

 !Do some checks
 if(ntime /= hist%mxhist)then
   write(msg,'(a)')'ntime is not correct'
   MSG_BUG(msg)
 end if

 if(natom /= size(hist%xred,2)) then
   write(msg,'(a)')'natom is not correct'
   MSG_BUG(msg)
 end if


 ABI_ALLOCATE(xred,(3,natom))
 xred = zero

!option 1 => set the reference for the test
! call xcart2xred(eff_pot%supercell%natom,eff_pot%supercell%rprimd,&
!&                eff_pot%supercell%xcart,xred)
! rprimd =  eff_pot%supercell%rprimd

!option 2 => set a specific step for the test
 istep = 4
 xred = hist%xred(:,:,istep)
 rprimd =  hist%rprimd(:,:,istep)

 rprimd_ref =  eff_pot%supercell%rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 npt=5
 delta = 0.001
 deltalist = (/-2*delta,-delta,real(0.0,dp),delta,2*delta/)
 strain = zero

   do ia=1,natom
     do mu=1,3
       write(std_out,*) "atm: ",ia," dir: ",mu
       do ii=1,npt
         delt = deltalist(ii)

!        Get the initial displacement
         call effective_potential_getDisp(disp,du_delta,natom,rprimd,&
&                                         eff_pot%supercell%rprimd,1,xred_hist=xred,&
&                                         xcart_ref=eff_pot%supercell%xcart,&
&                                         compute_displacement = .true.,compute_duDelta = .true.)

!        Add the delta
         call xcart2xred(natom, rprimd, disp, disp_red)
         disp_red(mu,ia) = disp_red(mu,ia) + delt
         call xred2xcart(natom, rprimd, disp, disp_red)

         call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                          xred=xred,du_delta=du_delta,&
&                                          displacement=disp,compute_anharmonic=.true.,verbose=.false.)
         diff(ii) = energy

       end do

!  Get the initial displacement
   call effective_potential_getDisp(disp,du_delta,natom,rprimd,&
&                                   eff_pot%supercell%rprimd,1,xred_hist=xred,&
&                                   xcart_ref=eff_pot%supercell%xcart,&
&                                   compute_displacement = .true.,compute_duDelta = .true.)

   call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                    xred=xred,du_delta=du_delta,&
&                                    displacement=disp,compute_anharmonic=.true.,verbose=.false.)

   write(std_out,*) "Analyti:",fred(mu,ia)
   write(std_out,*) "FD     :",(-diff(5)+8*diff(4)-8*diff(2)+diff(1)) / (12*delta)
   write(std_out,*) "Diff(%):",abs(100*(fred(mu,ia)-((-diff(5)+8*diff(4)-8*diff(2)+diff(1))&
&             / (12*delta) )) / ((-diff(5)+8*diff(4)-8*diff(2)+diff(1)) / (12*delta) ))

 end do
end do


!  Fill the identity matrix
identity = zero
forall(ii=1:3)identity(ii,ii)=1

 npt=5
 delta = 0.0005
 deltalist = (/-2*delta,-delta,real(0.0,dp),delta,2*delta/)

 do jj=1,6
   write(std_out,*) "strain ",jj
   do ii=1,npt
     strain = zero
     delt = deltalist(ii)
     mat_def = zero
     strain(jj) = strain(jj) + delt

     mat_def(alpha(jj),beta(jj)) = mat_def(alpha(jj),beta(jj)) + half * strain(jj)
     mat_def(beta(jj),alpha(jj)) = mat_def(beta(jj),alpha(jj)) + half * strain(jj)

     rprimd_def =  matmul(identity(:,:)+mat_def(:,:),rprimd)

! The two options should give the same result
! Option 1 => compute the disps and provide them to evaluate
!      call effective_potential_getDisp(disp,du_delta,natom,rprimd_def,&
! &                                     rprimd_ref,1,xred_hist=xred,&
! &                                     xcart_ref=eff_pot%supercell%xcart,&
! &                                     compute_displacement = .true.,compute_duDelta = .true.)

!      call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd_def,&
! &                                      xred=xred,du_delta=du_delta,&
! &                                      displacement=disp,strain=strain,&
! &                                      compute_anharmonic=.true.,verbose=.false.)

!   Option 2 => compute the disps within evaluate
    call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd_def,&
&                                     xred=xred,compute_anharmonic=.true.,verbose=.false.)



     diff(ii) = energy

   end do

!  The two options should give the same result
!  Option 1 => compute the disps and provide them to evaluate
!    call effective_potential_getDisp(disp,du_delta,natom,rprimd,&
! &                                   rprimd_ref,1,xred_hist=xred,&
! &                                   xcart_ref=eff_pot%supercell%xcart,&
! &                                   compute_displacement = .true.,compute_duDelta = .true.)

!    call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
! &                                    xred=xred,du_delta=du_delta,&
! &                                    displacement=disp,&
! &                                    compute_anharmonic=.true.,verbose=.false.)

!  Option 2 => compute the disps within evaluate
   call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                    xred=xred,compute_anharmonic=.true.,verbose=.false.)

 write(std_out,*) "Analyti:",strten(jj)
 write(std_out,*) "FD     :",(-diff(5)+8*diff(4)-8*diff(2)+diff(1)) / (12*delta) / ucvol
 write(std_out,*) "Diff(%):",abs(100*(strten(jj)-((-diff(5)+8*diff(4)-8*diff(2)+diff(1))&
&   / (12*delta) / ucvol)) / ((-diff(5)+8*diff(4)-8*diff(2)+diff(1)) / (12*delta) / ucvol))

end do

end subroutine effective_potential_checkDEV
!!***

!!****f* m_effective_potential/effective_potential_writeNETCDF
!! NAME
!! effective_potential_writeNETCDF
!!
!! FUNCTION
!! This routine print the effective potential into netcdf format
!! Several options are available
!!
!! INPUTS
!! filename = the name of output file
!! eff_pot  = datatype contains the effective potential
!! option   = option for the format of the xml file
!!            1 print the xml for a system
!!
!! OUTPUT
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      ab_define_var,isfile,wrtout
!!
!! SOURCE

subroutine effective_potential_writeNETCDF(eff_pot,option,filename)

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
 character(len=500) :: msg
 character(len=fnlen) :: namefile
!arrays
 real(dp) :: strain(9,6)

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

   write(msg,'(a,a,a)')ch10,&
&   ' Generation of the xml file for the reference structure in ',namefile

   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

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

!AM_END_EXPERIMENTAL SECTION

end module m_effective_potential
!!***
