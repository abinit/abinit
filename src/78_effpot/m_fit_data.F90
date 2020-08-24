!!****m* ABINIT/m_fit_data
!!
!! NAME
!! m_fit_data
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2010-2020 ABINIT group (AM)
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

module m_fit_data

 use defs_basis
 use m_errors
 use m_abicore

 use m_geometry,     only : metric
 
 implicit none
!!***

!!****t* m_fit_data/training_set_type
!! NAME
!! training_set_type
!!
!! FUNCTION
!! datatype for with all the information for the fit process
!! This structure contains all the informations of the training set:
!!   - ntime
!!   - displacement
!!   - du_delta
!!   - strain
!!   - sqomega
!!   - ucvol
!! 
!! SOURCE

 type, public :: training_set_type

   integer :: ntime
!    Number of time in the training set
   
   integer :: natom
!    Number of atoms in the training set
   
   real(dp), allocatable :: displacement(:,:,:)
!    displacement(3,natom,ntime)
!    displacement array, difference of the position in cartisian coordinates
!    between training set and reference

   real(dp), allocatable :: du_delta(:,:,:,:)
!    du_delta(6,3,natom,ntime)
!    du_delta array, variation of displacement wrt to the strain

   real(dp), allocatable :: strain(:,:)
!    strain(6,ntime)
!    strain array, strain in the training set

   real(dp), allocatable :: sqomega(:)
!    sqomega(ntime)
!    sqomega(itime) = (((ucvol(itime)**(-2.))* ((natom_sc)**(0.5)))**(-1.0/3.0))**2

   real(dp), allocatable :: ucvol(:)
!    ucvol(ntime)
!    ucvol array, volume of the cell in the training set

 end type training_set_type
 
!routine for training_set
 public :: training_set_init
 public :: training_set_free 
!!***

!----------------------------------------------------------------------
 
!!****t* m_fit_data/fit_data_type
!! NAME
!! fit_data_type
!!
!! FUNCTION
!! 
!! SOURCE

 type, public :: fit_data_type

   integer :: ntime
!   Number of time in the training set
   
   integer :: natom
!   Number of atoms in the training set

   real(dp),allocatable :: energy_diff(:)
!   energy(ntime)
!   Array with the diffence between energy from training set and energy from initial model.
!   The model constains only harmonic part

   real(dp),allocatable :: fcart_diff(:,:,:)
!   fcart_diff(3,natom,ntime)
!   Array with the diffence between cartesian forces from training set and forces from initial model.
!   The model constains only harmonic part

   real(dp),allocatable :: strten_diff(:,:)
!   strten_diff(6,ntime)
!   Array with the diffence between strain from training set and strain from initial model.
!   The model constains only harmonic part

   type(training_set_type) :: training_set
!    datatype with the informations of the training set
   
 end type fit_data_type

!routine for fit_data
 public :: fit_data_compute
 public :: fit_data_init
 public :: fit_data_free 
!!***

CONTAINS  !===========================================================================================

!!****f* m_fit_data/fit_data_init
!!
!! NAME
!! fit_data_init
!!
!! FUNCTION
!! Initialize fit_data datatype
!!
!! INPUTS
!! energy_diff(3,natom,ntime) = Difference of energy between DFT calculation and 
!!                             fixed part of the model (more often harmonic part)
!! fcart_diff(3,natom,ntime) = Difference of cartesian forces between DFT calculation and 
!!                             fixed part of the model (more often harmonic part)
!! natom = Number of atoms
!! ntime = Number of time (number of snapshot, number of md step...)
!! strten_diff(6,natom) = Difference of stress tensor between DFT calculation and 
!!                        fixed part of the model (more often harmonic part)
!! sqomega(ntime) =  Sheppard and al Factors \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]]
!! ucvol(ntime) = Volume of the system for each time
!! ts<training_set_type> = datatype with the information about the training set
!!
!! OUTPUT
!! fit_data<fit_data_type> = fit_data datatype to be initialized
!!
!! PARENTS
!!      m_fit_data
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_data_init(fit_data,energy_diff,fcart_diff,natom,ntime,strten_diff,ts)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntime
!arrays
 real(dp),intent(in) :: energy_diff(ntime),fcart_diff(3,natom,ntime)
 real(dp),intent(in) :: strten_diff(6,ntime)
 type(training_set_type),intent(in) :: ts
 type(fit_data_type),intent(inout) :: fit_data
!Local variables-------------------------------
!scalar
!arrays
 character(len=500) :: message
! *************************************************************************

 if(natom /= ts%natom)then
   write(message, '(a)')&
&      ' The number of atoms does not correspond to the training set'
   MSG_BUG(message)
 end if

 if(ntime /= ts%ntime)then
   write(message, '(a)')&
&      ' The number of time does not correspond to the training set'
   MSG_BUG(message)
 end if
 
!Free the output 
 call fit_data_free(fit_data)

!Set integer values 
 fit_data%ntime = ntime
 fit_data%natom = natom

!allocate arrays
 ABI_ALLOCATE(fit_data%fcart_diff,(3,natom,ntime))
 fit_data%fcart_diff(:,:,:) = fcart_diff(:,:,:)
 
 ABI_ALLOCATE(fit_data%strten_diff,(6,ntime))
 fit_data%strten_diff(:,:) = strten_diff(:,:)

 ABI_ALLOCATE(fit_data%energy_diff,(ntime))
 fit_data%energy_diff(:) = energy_diff
 
 call training_set_init(fit_data%training_set,ts%displacement,ts%du_delta,&
&                       natom,ntime,ts%strain,ts%sqomega,ts%ucvol)
 
end subroutine fit_data_init
!!***

!!****f* m_fit_data/fit_data_free
!!
!! NAME
!! fit_data_free
!!
!! FUNCTION
!! Free the fit_data datatype
!!
!! INPUTS
!! fit_data<fit_data_type> = fit_data to be free
!! OUTPUT
!!
!! PARENTS
!!      m_fit_data,m_fit_polynomial_coeff,m_opt_effpot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_data_free(fit_data)

 implicit none
  
!Arguments ------------------------------------
!scalars
!arrays
 type(fit_data_type),intent(inout) :: fit_data
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************

! Reset integer values
  fit_data%ntime = 0
  fit_data%natom = 0

! Deallocate arrays  
  if(allocated(fit_data%energy_diff)) then
    ABI_DEALLOCATE(fit_data%energy_diff)
  end if
  if(allocated(fit_data%fcart_diff)) then
    ABI_DEALLOCATE(fit_data%fcart_diff)
  end if
  if(allocated(fit_data%strten_diff)) then
    ABI_DEALLOCATE(fit_data%strten_diff)
  end if
  call training_set_free(fit_data%training_set)
   
end subroutine fit_data_free
!!***

!!****f* m_fit_data/fit_data_compute
!!
!! NAME
!! fit_data_compute
!!
!! FUNCTION
!! Conpute the strain of each configuration.
!! Compute the displacmeent of each configuration.
!! Compute the variation of the displacement due to strain of each configuration.
!! Compute fixed forces and stresse and get the standard deviation.
!! Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]]
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = The history of the MD (or snapshot of DFT
!! comm = MPI communicator
!! verbose  = optional, flag for the verbose mode
!!
!! OUTPUT
!! fit_data<fit_data_type> = fit_data is now filled
!!
!! PARENTS
!!      m_fit_polynomial_coeff,m_opt_effpot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_data_compute(fit_data,eff_pot,hist,comm,verbose)

 use m_strain,only : strain_type,strain_get
 use m_effective_potential,only : effective_potential_type,effective_potential_evaluate
 use m_effective_potential,only : effective_potential_getDisp
 use m_abihist, only : abihist
 use m_strain,only : strain_type,strain_get
 implicit none
  
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 logical,optional,intent(in) :: verbose
 !arrays
 type(fit_data_type),intent(inout) :: fit_data
 type(effective_potential_type),intent(in) :: eff_pot
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalar
 integer :: ii,itime,natom,ntime
 real(dp):: energy
 logical :: need_verbose
!arrays
 character(len=500) :: message
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: energy_diff(:)
 real(dp),allocatable :: du_delta(:,:,:,:),displacement(:,:,:),strain(:,:)
 real(dp),allocatable :: fcart_diff(:,:,:),fred_fixed(:,:,:),fcart_fixed(:,:,:)
 real(dp),allocatable :: strten_diff(:,:),strten_fixed(:,:),sqomega(:),ucvol(:)
 type(strain_type) :: strain_t
 type(training_set_type) :: ts
! *************************************************************************

!Initialisation of optional arguments
 need_verbose = .TRUE.
 if(present(verbose)) need_verbose = verbose

!Checks
 natom = eff_pot%supercell%natom
 ntime = hist%mxhist
 if(natom /= size(hist%xred,2))then
   write(message, '(5a)' )&
&      'The number of atoms in the hist file does not correspond to the supercell.',ch10,&
&      'You should call the routine fit_polynomial_coeff_mapHistToRef before',ch10,&
&      'Action: Contact abinit group'
   MSG_BUG(message)
 end if

!Allocation of temporary arrays
 ABI_ALLOCATE(displacement,(3,natom,ntime))
 ABI_ALLOCATE(du_delta,(6,3,natom,ntime))
 ABI_ALLOCATE(energy_diff,(ntime))
 ABI_ALLOCATE(fcart_fixed,(3,natom,ntime))
 ABI_ALLOCATE(fcart_diff,(3,natom,ntime))
 ABI_ALLOCATE(fred_fixed,(3,natom,ntime))
 ABI_ALLOCATE(strain,(6,ntime))
 ABI_ALLOCATE(strten_fixed,(6,ntime))
 ABI_ALLOCATE(strten_diff,(6,ntime))
 ABI_ALLOCATE(sqomega,(ntime))
 ABI_ALLOCATE(ucvol,(ntime))

 displacement = zero 
 du_delta = zero 
 strain = zero; 
 fcart_fixed  = zero 
 fred_fixed  = zero
 strain = zero 
 strten_fixed = zero 
 strten_diff = zero
 sqomega = zero 
 ucvol = zero

 do itime=1,ntime
!  Get strain
   call strain_get(strain_t,rprim=eff_pot%supercell%rprimd,&
&                  rprim_def=hist%rprimd(:,:,itime),symmetrized=.FALSE.)
   if (strain_t%name /= "reference")  then
     do ii=1,3
       strain(ii,itime) = strain_t%strain(ii,ii)
     end do
     strain(4,itime) = (strain_t%strain(2,3) + strain_t%strain(3,2))
     strain(5,itime) = (strain_t%strain(3,1) + strain_t%strain(1,3))
     strain(6,itime) = (strain_t%strain(2,1) + strain_t%strain(1,2))
   else
     strain(:,itime) = zero
   end if

!  Get displacement and du_delta
   call effective_potential_getDisp(displacement(:,:,itime),du_delta(:,:,:,itime),natom,&
&                                   hist%rprimd(:,:,itime),eff_pot%supercell%rprimd,comm,&
&                                   xred_hist=hist%xred(:,:,itime),xcart_ref=eff_pot%supercell%xcart,&
&                                   compute_displacement=.TRUE.,compute_duDelta=.TRUE.)

!  Get forces and stresses from harmonic part (fixed part)     
   call effective_potential_evaluate(eff_pot,energy,fcart_fixed(:,:,itime),fred_fixed(:,:,itime),&
&                                    strten_fixed(:,itime),natom,hist%rprimd(:,:,itime),&
&                                    displacement=displacement(:,:,itime),&
&                                    du_delta=du_delta(:,:,:,itime),strain=strain(:,itime),&
&                                    compute_anharmonic=.true.,verbose=.FALSE.)
   
!  Compute \Omega^{2} and ucvol for each time
   call metric(gmet,gprimd,-1,rmet,hist%rprimd(:,:,itime),ucvol(itime))
!  Formula: sqomega(itime) = (((ucvol(itime)**(-2.))* ((natom)**(0.5)))**(-1.0/3.0))**2
!   Compact form:
   sqomega(itime) = ((ucvol(itime)**(4.0/3.0)) / ((natom)**(1/3.0)))

!  Compute the difference between History and model (fixed part)
   fcart_diff(:,:,itime) =  hist%fcart(:,:,itime) - fcart_fixed(:,:,itime)
   energy_diff(itime)    =  hist%etot(itime) - energy
   strten_diff(:,itime)  =  hist%strten(:,itime) - strten_fixed(:,itime)
 end do ! End Loop itime 
   
!Set the training set
 call training_set_init(ts,displacement,du_delta,natom,ntime,strain,sqomega,ucvol)
!Set the fit_data
 call fit_data_init(fit_data,energy_diff,fcart_diff,natom,ntime,strten_diff,ts)

!Free space
 call training_set_free(ts)
 ABI_DEALLOCATE(displacement)
 ABI_DEALLOCATE(du_delta)
 ABI_DEALLOCATE(energy_diff)
 ABI_DEALLOCATE(fcart_fixed)
 ABI_DEALLOCATE(fcart_diff)
 ABI_DEALLOCATE(fred_fixed)
 ABI_DEALLOCATE(strain)
 ABI_DEALLOCATE(strten_fixed)
 ABI_DEALLOCATE(strten_diff)
 ABI_DEALLOCATE(sqomega)
 ABI_DEALLOCATE(ucvol)

end subroutine fit_data_compute
!!***




!!****f* m_fit_data/training_set_init
!!
!! NAME
!! training_set_init
!!
!! FUNCTION
!! Initialize training_set datatype
!!
!! INPUTS
!! du_delta(6,3,natom,ntime)  = Variation to displacements wrt to the strain (Bohr)
!! displacement(3,natom,ntime)= Atomic displacement wrt to the reference (Bohr)
!! natom = number of atoms
!! ntime = number of time step
!! strain(6,ntime) = Strain
!! sqomega =  Sheppard and al Factors \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]]
!! ucvol(ntime) = Volume of the supercell for each time (Bohr^3)
!!
!! OUTPUT
!! ts<training_set_type> = training set to be initialized
!!
!! PARENTS
!!      m_fit_data
!!
!! CHILDREN
!!
!! SOURCE

subroutine training_set_init(ts,displacement,du_delta,natom,ntime,strain,sqomega,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntime
!arrays
 real(dp),intent(in) :: displacement(3,natom,ntime),du_delta(6,3,natom,ntime)
 real(dp),intent(in) :: strain(6,ntime),sqomega(ntime),ucvol(ntime)
 type(training_set_type),intent(inout) :: ts  
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************

!Free the output 
 call training_set_free(ts)

!Set integer values 
 ts%ntime = ntime
 ts%natom = natom

!allocate arrays
 ABI_ALLOCATE(ts%displacement,(3,natom,ntime))
 ts%displacement(:,:,:) = displacement(:,:,:)
 
 ABI_ALLOCATE(ts%du_delta,(6,3,natom,ntime))
 ts%du_delta(:,:,:,:) = du_delta(:,:,:,:)
 
 ABI_ALLOCATE(ts%strain,(6,ntime))
 ts%strain(:,:) = strain(:,:)
 
 ABI_ALLOCATE(ts%sqomega,(ntime))
 ts%sqomega(:) = sqomega(:)
 
 ABI_ALLOCATE(ts%ucvol,(ntime))
 ts%ucvol(:) = ucvol(:)
 
end subroutine training_set_init
!!***

!!****f* m_fit_data/training_set_free
!!
!! NAME
!! training_set_free
!!
!! FUNCTION
!! Free the training_set datatype
!!
!! INPUTS
!! training_set<training_set_type> = training_set to be free
!! OUTPUT
!!
!! PARENTS
!!      m_fit_data
!!
!! CHILDREN
!!
!! SOURCE

subroutine training_set_free(ts)

 implicit none
  
!Arguments ------------------------------------
!scalars
!arrays
 type(training_set_type),intent(inout) :: ts  
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************

! Reset integer values
  ts%ntime = 0
  ts%natom = 0

! Deallocate arrays  
  if(allocated(ts%displacement)) then
    ABI_DEALLOCATE(ts%displacement)
  end if
  if(allocated(ts%du_delta)) then
    ABI_DEALLOCATE(ts%du_delta)
  end if
  if(allocated(ts%strain)) then
    ABI_DEALLOCATE(ts%strain)
  end if
  if(allocated(ts%sqomega))then
    ABI_DEALLOCATE(ts%sqomega)
  end if
  if(allocated(ts%ucvol)) then
    ABI_DEALLOCATE(ts%ucvol)
  end if
   
end subroutine training_set_free
!!***

end module m_fit_data
!!***
