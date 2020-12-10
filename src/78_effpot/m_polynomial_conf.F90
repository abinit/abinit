!!****f* ABINIT/m_polynomial_conf
!!
!! NAME
!! m_polynomial_conf
!!
!! FUNCTION
!! Module for using a confinement potential
!! Container type is defined, and destruction
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

module m_polynomial_conf

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi,only   : xmpi_sum

 implicit none

 public :: polynomial_conf_init
 public :: polynomial_conf_free
!AM_2017Need to debug this routine
 public :: polynomial_conf_evaluate
!!***

!!****t* m_polynomial_conf/polynomial_conf_type
!! NAME
!! polynomial_conf_type
!!
!! FUNCTION
!! datatype for specific confinement potential
!!
!! SOURCE

 type, public :: polynomial_conf_type

   integer :: ndisp = 0
!  Number of displacement (atoms) for the cut off

   integer :: power_disp = 0
!  Power of the polynome related to the displacement

   integer :: power_strain = 0
!  Power of the polynome related to the strain

   real(dp):: factor_disp = 0
!  Factor to appy to the polynomial term of the confinement (displacement)

   real(dp):: factor_strain = 0
!  Factor to appy to the polynomial term of the confinement (strain)

   real(dp):: cutoff_strain(6)
!  Cutoff array for the strain

   real(dp),allocatable :: cutoff_disp(:)
!  Cutoff array for the atomic displacement

   logical :: need_confinement =.FALSE.
!  Logical related to the necessity of the confinement

 end type polynomial_conf_type
!!***


CONTAINS  !===========================================================================================


!!****f* m_polynomial_conf/polynomial_conf_init
!!
!! NAME
!! polynomial_conf_init
!!
!! FUNCTION
!! Initialize polynomial_conf_init
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
!! polynomial_conf <type(polynomial_conf)> = datatype with the information for the confinement
!!                                           polynomial
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE
!!

subroutine polynomial_conf_init(cutoff_disp,cutoff_strain,factor_disp,factor_strain,ndisp,&
&                               polynomial_conf,power_disp,power_strain,need_confinement)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ndisp,power_disp,power_strain
 real(dp),intent(in) :: factor_disp,factor_strain
 logical,optional,intent(in)  :: need_confinement
!arrays
 real(dp),intent(in) :: cutoff_disp(ndisp),cutoff_strain(6)
 type(polynomial_conf_type),intent(inout) :: polynomial_conf
!Local variables-------------------------------
!scalar
!arrays
 character(len=500) :: msg

! *************************************************************************

!Checks
 if (ndisp <= 0) then
   write(msg,'(a,a)')' ndisp can not be inferior or equal to zero'
   ABI_ERROR(msg)
 end if

!First free the type
 call  polynomial_conf_free(polynomial_conf)

 polynomial_conf%power_disp    = power_disp
 polynomial_conf%power_strain  = power_strain
 polynomial_conf%factor_disp   = factor_disp
 polynomial_conf%factor_strain = factor_strain
 polynomial_conf%need_confinement = .FALSE.

 polynomial_conf%ndisp   = ndisp
 ABI_ALLOCATE(polynomial_conf%cutoff_disp,(polynomial_conf%ndisp))
 polynomial_conf%cutoff_disp(:) = cutoff_disp(:)

 polynomial_conf%cutoff_strain = cutoff_strain(:)
 if (present(need_confinement)) polynomial_conf%need_confinement = need_confinement

end subroutine polynomial_conf_init
!!***


!!****f* m_polynomial_conf/polynomial_conf_free
!!
!! NAME
!! polynomial_conf_free
!!
!! FUNCTION
!! Free polynomial_conf
!!
!! INPUTS
!! polynomial_conf <type(polynomial_conf)> = polynomial_conf datatype to be free
!!
!! OUTPUT
!! polynomial_conf <type(polynomial_conf)> = polynomial_conf datatype to be free
!!
!! PARENTS
!!      m_effective_potential,m_polynomial_conf
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine polynomial_conf_free(polynomial_conf)

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(polynomial_conf_type), intent(inout) :: polynomial_conf
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************

 if(allocated(polynomial_conf%cutoff_disp))then
   ABI_DEALLOCATE(polynomial_conf%cutoff_disp)
 end if

 polynomial_conf%power_disp    = 0
 polynomial_conf%power_strain  = 0
 polynomial_conf%factor_disp   = zero
 polynomial_conf%factor_strain = zero
 polynomial_conf%cutoff_strain = zero
 polynomial_conf%need_confinement = .FALSE.

end subroutine polynomial_conf_free
!!***

!!****f* m_polynomial_conf/polynomial_conf_evaluate
!! NAME
!!  polynomial_conf_evaluate
!!
!! FUNCTION
!!  This fonction evaluate the energy, (soon forces and stresses) with  the confinement potential
!!
!! INPUTS
!!   disp(3,natom_sc) = atomic displacments of a specific patern wrt to reference structure
!!   disp_ref(natom_uc) = Cutoff array for the atomic displacement
!!   factor_disp = Factor to appy to the polynomial term of the confinement (displacement)
!!   factor_strain = Factor to appy to the polynomial term of the confinement (strain)
!!   strain(6) =   strain of a specific structure wrt to reference
!!   strain_ref(6) = Cutoff array for the strain
!!   power_disp = Power of the polynome related to the displacement
!!   power_strain = Power of the polynome related to the strain
!!   natom_sc = number of atoms in the supercell
!!   natom_uc = number of atoms in the unit cell
!!   ncell   = total number of cell to treat by this cpu
!!   cells(ncell) = index of  the cells into the supercell (1,2,3,4,5)
!!   index_cells(3,ncell) = indexes  of the cells into supercell (-1 -1 -1 ,...,1 1 1)
!!   comm=MPI communicator
!!
!! OUTPUT
!!   energy = contribution to the ifc to the energy
!!   fcart(3,natom) = contribution to the ifc to the forces
!!   strten(6) = contribution to the stress tensor
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine polynomial_conf_evaluate(disp,disp_ref,energy,factor_disp,factor_strain,fcart,&
&                                   strain,strain_ref,strten,power_disp,power_strain,cells,&
&                                   natom_sc,natom_uc,ncell,index_cells,comm)

 implicit none

!Arguments -------------------------------
! scalars
  real(dp),intent(out) :: energy
  integer,intent(in) :: natom_uc,natom_sc,ncell
  integer,intent(in) :: power_disp,power_strain
  integer,intent(in) :: comm
  real(dp),intent(in) :: factor_disp,factor_strain
! array
  integer,intent(in) ::  cells(ncell),index_cells(ncell,3)
  real(dp),intent(in) :: disp(3,natom_sc),strain(6)
  real(dp),intent(out) :: fcart(3,natom_sc),strten(6)
  real(dp),intent(in) :: disp_ref(natom_uc),strain_ref(6)
!Local variables-------------------------------
! scalar
  integer :: ia,icell,ierr,ii,kk
  integer :: mu
  real(dp):: diff,diff_tmp
! array

! *************************************************************************

! Initialisation of variables
  energy   = zero
  fcart(:,:) = zero
  strten(:) = zero

  write(std_out,*) factor_strain,index_cells,power_strain,strain,strain_ref
  do icell = 1,ncell
    ii = (cells(icell)-1)*natom_uc
    do ia = 1, natom_uc
      kk = ii + ia
      diff_tmp = zero
      do mu=1,3
        diff_tmp = diff_tmp + disp(mu,kk)**2
      end do
      diff_tmp = diff_tmp**0.5
!     Compute diff between ref and curent displacement
      diff = diff_tmp - disp_ref(ia)
!     Accumule energy
      energy =  energy + (sign(half, diff)+half)*(factor_disp*((diff_tmp/disp_ref(ia))**power_disp))
    end do
  end do

! MPI_SUM
  call xmpi_sum(energy, comm, ierr)

end subroutine polynomial_conf_evaluate
!!***

end module m_polynomial_conf
!!***
