!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_polynomial_conf
!!
!! NAME
!! m_polynomial_conf
!!
!! FUNCTION
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

module m_polynomial_conf

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 public :: polynomial_conf_init
 public :: polynomial_conf_free
!!***

!!****t* m_polynomial_conf/polynomial_conf_type
!! NAME
!! polynomial_conf_type
!!
!! FUNCTION
!! structure for specific confinement potential
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
!  Factor to appy to the polynomial term of the cofinement (displacement)

   real(dp):: factor_strain = 0
!  Factor to appy to the polynomial term of the cofinement (strain)

   real(dp):: cutoff_strain(6)
!  Cutoff array for the strain

   real(dp),allocatable :: cutoff_disp(:)

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
!!
!! OUTPUT
!!
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_conf_init(cutoff_disp,cutoff_strain,factor_disp,factor_strain,ndisp,&
&                               polynomial_conf,power_disp,power_strain,need_confinement)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_conf_init'
!End of the abilint section

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
   MSG_ERROR(msg)
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
!!
!! OUTPUT
!! polynomial_conf = polynomial_conf structure to be free
!!
!! PARENTS
!!      m_effective_potential_file,m_polynomial_coeff,m_polynomial_conf
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_conf_free(polynomial_conf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_conf_free'
!End of the abilint section

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

 polynomial_conf%power_disp    = zero
 polynomial_conf%power_strain  = zero
 polynomial_conf%factor_disp   = zero
 polynomial_conf%factor_strain = zero
 polynomial_conf%cutoff_strain = zero
 polynomial_conf%need_confinement = .FALSE.

end subroutine polynomial_conf_free
!!***

end module m_polynomial_conf
!!***
