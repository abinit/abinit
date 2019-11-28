!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_harmonics_terms
!!
!! NAME
!! m_harmonics_term
!!
!! FUNCTION
!! Module with datatype and tools for the harmonics terms
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

module m_harmonics_terms

 use defs_basis
 use m_errors
 use m_abicore
 use m_supercell,only: getPBCIndexes_supercell
 use m_xmpi,only : xmpi_sum
 use m_ifc

 implicit none

 public :: harmonics_terms_init
 public :: harmonics_terms_free
 public :: harmonics_terms_applySumRule
 public :: harmonics_terms_evaluateIFC
 public :: harmonics_terms_evaluateElastic
 public :: harmonics_terms_setEffectiveCharges
 public :: harmonics_terms_setDynmat
 public :: harmonics_terms_setInternalStrain
!!***

!!****t* m_harmonics_terms/harmonics_terms_type
!! NAME
!! harmonics_terms_type
!!
!! FUNCTION
!! datatype for harmonic part of effective potential.
!!
!! SOURCE

 type, public :: harmonics_terms_type

   integer :: nqpt
!   Number of qpoints

   real(dp) :: epsilon_inf(3,3)
!   epsilon_inf(3,3)
!   Dielectric tensor

   real(dp) :: elastic_constants(6,6)
!   elastic_constant(6,6)
!   Elastic tensor Hartree

   real(dp), allocatable :: strain_coupling(:,:,:)
!   strain_coupling(6,3,natom)
!   internal strain tensor

   real(dp), allocatable :: zeff(:,:,:)
!   zeff(3,3,natom) Effective charges

   type(ifc_type) :: ifcs
!   type with ifcs constants (short + ewald)
!   also contains the number of cell and the indexes

   real(dp), allocatable :: qpoints(:,:)
!   qph1l(3,nqpt)
!   List of qpoints wavevectors

   real(dp), allocatable :: dynmat(:,:,:,:,:,:)
!   dynmat(2,3,natom,3,natom,nqpt)
!   dynamical matrix for each q points

   real(dp), allocatable :: phfrq(:,:)
!   phfrq(3*natom,nqpt)
!   array with all phonons frequencies for each q points in Hartree/cm

 end type harmonics_terms_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_harmonics_terms/harmonics_terms_init
!!
!! NAME
!! harmonics_terms_init
!!
!! FUNCTION
!! Initialize harmonics_terms datatype
!!
!! INPUTS
!! ifc<type(ifc_type)> = interatomic forces constants
!! natom = number of atoms in primitive cell
!! nrpt = number rpt (cell) in the ifc
!! dynmat(2,3,natom,3,natom,3,nqpt) = optional, dynamical matricies for each q-point
!! epsilon_inf(3,3) = optional, dielectric tensor
!! elastic_constant(6,6) = optional, elastic constant
!! strain_coupling(6,3,natom) = optional, internal strain coupling parameters
!! nqpt = optional, number of q-points
!! phfrq(3*natom,nqpt) = optional,phonons frequencies for each q points in Hartree/cm
!! qpoints(3,nqpt) = list of qpoints wavevectors
!! zeff(3,3,natom) = optional,effective charges
!!
!! OUTPUT
!! harmonics_terms<type(harmonics_terms_type)> = harmonics_terms datatype to be initialized
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_init(harmonics_terms,ifcs,natom,nrpt,&
&                               dynmat,epsilon_inf,elastic_constants,strain_coupling,&
&                               nqpt,phfrq,qpoints,zeff)

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
 real(dp),optional,intent(in) :: strain_coupling(6,3,natom),zeff(3,3,natom)
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
 ABI_ALLOCATE(harmonics_terms%ifcs%atmfrc,(3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%atmfrc(:,:,:,:,:) = ifcs%atmfrc(:,:,:,:,:)

!Allocation of ewald part of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%ewald_atmfrc,(3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%ewald_atmfrc(:,:,:,:,:) = ifcs%ewald_atmfrc(:,:,:,:,:)

!Allocation of short range part of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%short_atmfrc,(3,natom,3,natom,nrpt))
 harmonics_terms%ifcs%short_atmfrc(:,:,:,:,:) = ifcs%short_atmfrc(:,:,:,:,:)

!Allocation of cell of ifc
 ABI_ALLOCATE(harmonics_terms%ifcs%cell,(3,nrpt))
 harmonics_terms%ifcs%cell(:,:) = ifcs%cell(:,:)

!Allocation of the dynamical matrix
 harmonics_terms%nqpt = 0
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
 ABI_ALLOCATE(harmonics_terms%strain_coupling,(6,3,natom))
 harmonics_terms%strain_coupling = zero
 if (present(strain_coupling)) then
   call harmonics_terms_setInternalStrain(harmonics_terms,natom,strain_coupling)
 end if

end subroutine harmonics_terms_init
!!***


!****f* m_harmonics_terms/harmonics_terms_free
!!
!! NAME
!! harmonics_terms_free
!!
!! FUNCTION
!! deallocate all dynamic memory for this harmonic datatype
!!
!! INPUTS
!! harmonics_terms<type(harmonics_terms_type)> = harmonics_terms datatype to be free
!!
!! OUTPUT
!! harmonics_terms<type(harmonics_terms_type)> = harmonics_terms datatype to be free
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_free(harmonics_terms)

  implicit none

!Arguments ------------------------------------
!scalars
!array
  type(harmonics_terms_type), intent(inout) :: harmonics_terms
!Local variables-------------------------------
!scalars
!array

! *************************************************************************

  harmonics_terms%nqpt = 0
  harmonics_terms%elastic_constants = zero
  harmonics_terms%epsilon_inf       = zero

  if(allocated(harmonics_terms%zeff))then
    harmonics_terms%zeff=zero
    ABI_DEALLOCATE(harmonics_terms%zeff)
  end if

  if(allocated(harmonics_terms%strain_coupling)) then
    harmonics_terms%strain_coupling=zero
    ABI_DEALLOCATE(harmonics_terms%strain_coupling)
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

  call harmonics_terms%ifcs%free()

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
!! strain_coupling(6,3,natom) = internal strain coupling parameters
!!
!! OUTPUT
!! harmonics_terms<type(harmonics_terms_type)> = harmonics_terms datatype
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_setInternalStrain(harmonics_terms,natom,strain_coupling)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: natom
!array
  real(dp),intent(in) :: strain_coupling(:,:,:)
  type(harmonics_terms_type), intent(inout) :: harmonics_terms
!Local variables-------------------------------
!scalars
!array
  character(len=500) :: msg

! *************************************************************************

! 0-Checks inputs
  if(natom /= size(strain_coupling,3)) then
    write(msg, '(a)' )&
&        ' natom has not the same size strain_coupling array. '
    MSG_BUG(msg)
  end if

! 1-deallocate old array
  if(allocated(harmonics_terms%strain_coupling))then
    ABI_DEALLOCATE(harmonics_terms%strain_coupling)
  end if

! 2-allocate and copy the new array
  ABI_ALLOCATE(harmonics_terms%strain_coupling,(6,3,natom))
  harmonics_terms%strain_coupling(:,:,:) = strain_coupling(:,:,:)

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
!! harmonics_terms<type(harmonics_terms_type)> = harmonics_terms datatype
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_setEffectiveCharges(harmonics_terms,natom,zeff)

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
!! Set the dynamical matricies to the harmonics_terms
!!
!! INPUTS
!! natom = number of atoms
!! nqpt  = number of qpoints
!! dynmat(2,3,natom,3,natom,nqpt) = dynamical matrix in cartesian coordinates
!! phfrq(3*natom,nqpt) = frequency in hartree
!! qpoints(3,nqpt) = list of qpoints
!!
!! OUTPUT
!! harmonics_terms<type(harmonics_terms_type)> = harmonics_terms datatype
!!
!! PARENTS
!!      m_effective_potential,m_harmonics_terms
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_setDynmat(dynmat,harmonics_terms,natom,nqpt,phfrq,qpoints)

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
!!****f* m_harmonics_terms/harmonics_terms_evaluateIFC
!! NAME
!!  harmonics_terms_evaluateIFC
!!
!! FUNCTION
!!  This fonction compute the contribution of the ifc harmonic part of
!!  the energy and forces.
!!
!! INPUTS
!!  atmfrc(3,natom_uc,3,natom_uc,nrpt) = atomic force constants
!!  disp(3,natom_sc) = atomics displacement between configuration and the reference
!!  ncell = total number of cell to treat
!!  nrpt  = total number of rpt to treat
!!  natom_sc = number of atoms in the supercell
!!  natom_uc = number of atoms in the unit cell
!!  nrpt  = number of rpt
!!  atmrpt_index(nrpt,cell) = For each cell in the supercell and each rpt,
!!                            give the index of the first atoms in the rpt cell
!!  rpt(nrpt) = index of rpt in  atmfrc (6th dimension)
!!  index_cells(3,ncell) = indexes of the cells into  supercell (-1 -1 -1 ,...,1 1 1)
!!  comm=MPI communicator
!!
!! OUTPUT
!!   energy = contribution of the ifc to the energy
!!   fcart(3,natom) = contribution of the ifc to the forces
!!
!! PARENT
!!   effective_potential_evaluate
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_evaluateIFC(atmfrc,disp,energy,fcart,natom_sc,natom_uc,&
&                                      ncell,nrpt,atmrpt_index,index_cells,sc_size,rpt,comm)

 implicit none

!Arguments -------------------------------
! scalars
  real(dp),intent(out) :: energy
  integer,intent(in) :: natom_uc,natom_sc,ncell,nrpt
  integer,intent(in) :: comm
! array
  integer,intent(in) :: sc_size(3),atmrpt_index(nrpt,ncell)
  integer,intent(in) ::  index_cells(4,ncell),rpt(nrpt)
  real(dp),intent(in) :: atmfrc(3,natom_uc,3,natom_uc,nrpt)
  real(dp),intent(in) :: disp(3,natom_sc)
  real(dp),intent(out) :: fcart(3,natom_sc)

!Local variables-------------------------------
! scalar
  integer :: i1,i2,i3,ia,ib,icell,ierr,irpt,irpt_tmp,ii,jj,kk,ll
  integer :: mu,nu
  real(dp):: disp1,disp2,ifc,tmp_etot1,tmp_etot2
!Variables for separation of short and dipdip ifc contribution 
 !real(dp):: short_ifc,ewald_ifc
 !real(dp):: tmp_ewald1,tmp_ewald2,tmp_short1,tmp_short2
  ! array
  character(500) :: msg

! *************************************************************************

  if (any(sc_size <= 0)) then
    write(msg,'(a,a)')' sc_size can not be inferior or equal to zero'
    MSG_ERROR(msg)
  end if

! Initialisation of variables
  energy   = zero
  fcart(:,:) = zero

  do icell = 1,ncell
    i1 = index_cells(1,icell)
    i2 = index_cells(2,icell)
    i3 = index_cells(3,icell)
!   index of the first atom in the current cell
    ii = index_cells(4,icell)
    do irpt_tmp = 1,nrpt
      irpt = rpt(irpt_tmp)
!     index of the first atom in the irpt cell
      jj = atmrpt_index(irpt_tmp,icell)
!     Loop over the atom in the cell
      do ib = 1, natom_uc
        ll = jj + ib
        do nu=1,3
          disp2 = disp(nu,ll)
          do ia = 1, natom_uc
            kk = ii + ia
            do mu=1,3
              disp1 = disp(mu,kk)
              ifc = atmfrc(mu,ia,nu,ib,irpt)
              
!              if(abs(ifc) > tol10)then
                tmp_etot1  = disp2 * ifc
!               accumule energy
                tmp_etot2  = disp1*tmp_etot1
                energy =  energy + tmp_etot2
!               accumule forces
                fcart(mu,kk) = fcart(mu,kk) + tmp_etot1
!              end if
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

end subroutine harmonics_terms_evaluateIFC
!!***

!!****f* m_harmonics_terms/harmonics_terms_evaluateElastic
!! NAME
!!  harmonics_terms_evaluateElastic
!!
!! FUNCTION
!! Compute the energy, forces and stresses related to the application of strain
!!
!! INPUTS
!!  elastic_constants(6,6) = elastic constants in Hartree
!!  disp(3,natom_sc) = atomics displacement between configuration and the reference
!!  natom = number of atoms in the supercell
!!  natom_uc = number of atoms in the unit cell
!!  ncell = total number of cell
!!  strain_coupling(6,3,natom) = internal strain coupling parameters
!!  strain(6) = strain between configuration and the reference
!!
!! OUTPUT
!!   energy = contribution to the energy
!!   fcart(3,natom) = contribution to the forces
!!   strten(6) = contribution to the stress tensor
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE
!!
subroutine harmonics_terms_evaluateElastic(elastic_constants,disp,energy,fcart,natom,natom_uc,ncell,&
&                                          strain_coupling,strten,strain)

 real(dp),intent(out):: energy
 integer, intent(in) :: natom,natom_uc,ncell
! array
 real(dp),intent(in) :: elastic_constants(6,6),strain_coupling(6,3,natom)
 real(dp),intent(out):: strten(6)
 real(dp),intent(out):: fcart(3,natom)
 real(dp),intent(in) :: disp(3,natom)
 real(dp),intent(in) :: strain(6)

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
     cij = ncell*elastic_constants(alpha,beta)
     energy = energy + half*cij*strain(alpha)*strain(beta)
     strten(alpha) = strten(alpha) + cij*strain(beta)
   end do
 end do

!2-Part due to the internal strain coupling parameters
 ii = 1
 do ia = 1,natom
   do mu = 1,3
     do alpha=1,6
       cij = strain_coupling(alpha,mu,ii)
!      Accumulte for this atom
       energy = energy + half*cij*strain(alpha)*disp(mu,ia)
       fcart(mu,ia)  = fcart(mu,ia)  + half*cij*strain(alpha)
       strten(alpha) = strten(alpha) + half*cij*disp(mu,ia)
     end do
   end do
   ii = ii +1
!  Reset to 1 if the number of atoms is superior than in the initial cell
   if(ii==natom_uc+1) ii = 1
 end do

end subroutine  harmonics_terms_evaluateElastic
!!***

!****f* m_harmonics_terms/harmonics_terms_applySumRule
!!
!! NAME
!! harmonics_terms_applySumRule
!!
!! FUNCTION
!! Apply the acoustic sum rule on the inter-atomic force constants
!!
!! INPUTS
!! ifc<type(ifc_type)> = interatomic forces constants
!! asr   = acoustic sum rule option (see anaddb help)
!! natom = number of atoms
!! option = optional if |no present asr is done on total ifc
!!                       |present and 1 asr is done on short part
!!                       |present and 2 asr is done on ewald part
!!
!! OUTPUT
!! ifc<type(ifc_type)> = interatomic forces constants
!!
!! PARENTS
!!      compute_anharmonics,m_effective_potential
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine harmonics_terms_applySumRule(asr,ifc,natom,option)

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
 character(500) :: msg
!array
 real(dp),pointer :: atmfrc(:,:,:,:,:)
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
   write(msg,'(a,a)')' Unable to find the cell of reference in IFC'
   MSG_ERROR(msg)
 end if

 atmfrc => ifc%atmfrc
 if (present(option)) then
   if (option == 1) then
     nullify(atmfrc)
     atmfrc => ifc%short_atmfrc
     write(msg,'(3a)') ch10," Impose acoustic sum rule on short range"
   else if (option == 2) then
     nullify(atmfrc)
     atmfrc => ifc%ewald_atmfrc
     write(msg,'(3a)') ch10," Impose acoustic sum rule on long range"
   end if
 else
   write(msg,'(3a)') ch10," Impose acoustic sum rule on total ifc"
 end if
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

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
             sum=sum+atmfrc(mu,ia,nu,ib,irpt)
           end do
         else if(asr==2)then
           do irpt=1, ifc%nrpt
              sum=sum+&
&                 (atmfrc(mu,ia,nu,ib,irpt)+&
&                  atmfrc(nu,ia,mu,ib,irpt))/2
            end do
          end if
        end do

!      Correct the self-interaction in order to fulfill the ASR
        atmfrc(mu,ia,nu,ia,irpt_ref)=&
&       atmfrc(mu,ia,nu,ia,irpt_ref)-sum
        if(asr==2)then
          atmfrc(nu,ia,mu,ia,irpt_ref)=&
&         atmfrc(mu,ia,nu,ia,irpt_ref)
        end if
      end do
    end do
  end do

 if (present(option)) then
   if (option == 1) then
     ifc%short_atmfrc = atmfrc(:,:,:,:,:)
   else if (option == 2) then
     ifc%ewald_atmfrc = atmfrc(:,:,:,:,:)
   end if
 else
   ifc%atmfrc(:,:,:,:,:) = atmfrc(:,:,:,:,:)
 end if

 end subroutine harmonics_terms_applySumRule
!!***

end module m_harmonics_terms
!!***
