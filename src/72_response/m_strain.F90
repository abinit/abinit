!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_strain
!!
!! NAME
!! m_strain
!!
!! FUNCTION
!! Module for get the strain 
!! Container type is defined
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

module m_strain

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi

 implicit none

 private :: strain_def2strain
 private :: strain_strain2def
 private :: strain_print
 public  :: strain_get
 public  :: strain_init
 public  :: strain_apply
!!***
 
!!****t* defs_abitypes/strain_type
!! NAME
!! strain_type
!!
!! FUNCTION
!! structure for a effective potential constructed.
!!
!! SOURCE

 type, public :: strain_type
   character(len=fnlen) :: name
!   name of the strain (iso,uniaxial,shear...)

   real(dp) :: delta
!   Value of the strain

   integer :: direction
!   Direction of the strain (-1 if isostatic)

   real(dp) :: strain(3,3)
!   Matrix representing the strain 

 end type strain_type
!!***

CONTAINS  !===========================================================================================

!****f* m_strain/strain_init
!!
!! NAME
!! strain_init
!!
!! FUNCTION
!! routine to initialize strain structure
!!
!!
!! INPUTS
!!  name = name of the perturbation
!!  direction = direction of the perturbation
!!  delta = delta to apply in the strain (in percent)
!!
!! OUTPUT
!!  strain = structure with all information of strain
!!
!! PARENTS
!!   strain_phonon_coupling
!!
!! CHILDREN
!!
!! SOURCE
 
 subroutine strain_init(strain,delta,direction,name)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
   character(len=fnlen),optional,intent(in) :: name
   real(dp),optional,intent(in) :: delta
   integer,optional,intent(in) :: direction 
!array
   type(strain_type),intent(out) :: strain
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************
  if (present(name)) then
    strain%name = name
  else
    strain%name = ''
  end if

  if (present(delta)) then
    strain%delta = delta
  else
    strain%delta = zero
  end if

  if (present(direction)) then
    strain%direction = direction
  else
    strain%direction = zero
  end if

  call strain_strain2def(strain%strain,strain)

end subroutine strain_init
!!***

!****f* m_strain/strain_free
!!
!! NAME
!! strain_free
!!
!! FUNCTION
!! routine to free strain structure
!!
!!
!! INPUTS
!!  strain = structure with all information of strain
!!
!! OUTPUT
!!
!! PARENTS
!!   strain_phonon_coupling
!!
!! CHILDREN
!!
!! SOURCE
 
 subroutine strain_free(strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 type(strain_type),intent(inout) :: strain

!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************

    strain%name = ''
    strain%delta = zero
    strain%direction = zero
    strain%strain = zero

end subroutine strain_free
!!***


!****f* m_strain/strain_get
!!
!! NAME
!! strain_get
!!
!! FUNCTION
!! Get the strain for structure, compare to reference
!! structure and fill strain type
!!
!!
!! INPUTS
!!
!! OUTPUT
!!  strain = structure with all information of strain
!!
!! PARENTS
!!   strain_phonon_coupling
!!
!! CHILDREN
!!
!! SOURCE
 
 subroutine strain_get(rprim,rprim_def,strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_get'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 real(dp),intent(in) :: rprim(3,3),rprim_def(3,3)
 type(strain_type),intent(inout) :: strain
!Local variables-------------------------------
!scalar
  integer :: i,j
!arrays
  real(dp) :: mat_delta(3,3),rprim_inv(3,3)
  real(dp) :: identity(3,3)
! *************************************************************************
  write(std_out,'(a,a)') ' strain_get: try to get the strain'

  mat_delta = zero
! Fill the identity matrix
  identity = zero
  forall(i=1:3)identity(i,i)=1

  call matr3inv(rprim,rprim_inv)
  mat_delta =  matmul(rprim_inv,rprim_def)-identity

  identity = zero
  do i=1,3
    do j=1,3
      if (abs(mat_delta(i,j))>tol10) then 
        identity(i,j) = ANINT(mat_delta(i,j)*1000)/1000
      end if
    end do
  end do
  
  mat_delta = identity
  call strain_def2strain(mat_delta,strain)
  write(std_out,'(a,a)') ' end strain_get'

end subroutine strain_get
!!***

!****f* m_strain/strain_apply
!!
!! NAME
!! strain_get
!!
!! FUNCTION
!! Get the strain for structure, compare to reference
!! structure and fill strain type
!!
!!
!! INPUTS
!!
!! OUTPUT
!!  strain = structure with all information of strain
!!
!! PARENTS
!!   strain_phonon_coupling
!!
!! CHILDREN
!!
!! SOURCE
 
 subroutine strain_apply(rprim,rprim_def,strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_apply'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 real(dp),intent(in)  :: rprim(3,3)
 real(dp),intent(out) :: rprim_def(3,3)
 type(strain_type),intent(in) :: strain
!Local variables-------------------------------
!scalar
  integer :: i
!arrays
  real(dp) :: identity(3,3)
! *************************************************************************
  write(std_out,'(a,a)') ' strain_apply: try to get the strain'

! Fill the identity matrix
  identity = zero
  forall(i=1:3)identity(i,i)=1

  rprim_def(:,:) = matmul(rprim(:,:),strain%strain(:,:))

  write(std_out,'(a,a)') ' end strain_get'

end subroutine strain_apply
!!***

!****f* m_strain/strain_def2strain
!!
!! NAME
!! strain_matdef2strain
!!
!! FUNCTION
!! transfer deformation matrix in structure strain
!!
!! INPUTS
!! rprim = contains 
!!
!! OUTPUT
!!
!!
!! PARENTS
!!   epigene
!!
!! CHILDREN
!!
!! SOURCE
 
 subroutine strain_def2strain(mat_strain,strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_def2strain'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 real(dp),intent(in) :: mat_strain(3,3)
 type(strain_type),intent(inout) :: strain 
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************
 strain%name = ""
 strain%delta = zero
 strain%direction = zero
 strain%strain = mat_strain
 
 if (all(abs(mat_strain)<tol10)) then 
   strain%name = "reference"
   strain%delta = zero
   strain%direction = zero
   strain%strain = zero
 else
   if(abs(mat_strain(1,1))>tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))<tol10.and.&
&   abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))>tol10.and.abs(mat_strain(2,3))<tol10.and.&
&   abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))>tol10) then
     if(mat_strain(1,1)==mat_strain(2,2).and.mat_strain(1,1)==mat_strain(3,3)) then
       strain%name = "isostatic"
       strain%delta = mat_strain(1,1)
       strain%direction = -1
       strain%strain = mat_strain
     end if
   end if
   if(abs(mat_strain(1,1))>tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))<tol10.and.&
&    abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))<tol10) then
     strain%name = "uniaxial"
     strain%delta = mat_strain(1,1)
     strain%direction = 1
     strain%strain = mat_strain
   end if
   if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))<tol10.and.&
&    abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))>tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))<tol10) then
     strain%name = "uniaxial"
     strain%delta = mat_strain(2,2)
     strain%direction = 2
     strain%strain = mat_strain
   end if
   if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))<tol10.and.&
&    abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))>tol10) then
     strain%name = "uniaxial"
     strain%delta = mat_strain(3,3)
     strain%direction = 3
     strain%strain = mat_strain
    end if
    if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))<tol10.and.&
&    abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))>tol10.and.&
&    abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))>tol10.and.abs(mat_strain(3,3))<tol10) then
      if (mat_strain(3,2)==mat_strain(3,2)) then 
        strain%name = "shear"
        strain%delta = mat_strain(3,2) * 2
        strain%direction = 4
        strain%strain = mat_strain
      end if
    end if
    if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))>tol10.and.&
&    abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))>tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))<tol10) then
      if (mat_strain(3,1)==mat_strain(1,3)) then 
        strain%name = "shear"
        strain%delta = mat_strain(3,1) * 2
        strain%direction = 5
        strain%strain = mat_strain
      end if
    end if
    if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))>tol10.and.abs(mat_strain(1,3))<tol10.and.&
&    abs(mat_strain(2,1))>tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))<tol10) then
      if (mat_strain(1,2)==mat_strain(2,1)) then 
        strain%name = "shear"
        strain%delta = mat_strain(2,1) * 2
        strain%direction = 6
        strain%strain = mat_strain
      end if
    end if
  end if
  
  call  strain_print(strain)
  
end subroutine strain_def2strain
!!***

!****f* m_strain/strain_strain2def
!!
!! NAME
!! strain_matdef2strain
!!
!! FUNCTION
!! transfer deformation matrix in structure strain
!!
!! INPUTS
!! rprim = contains 
!!
!! OUTPUT
!!
!!
!! PARENTS
!!   epigene
!!
!! CHILDREN
!!
!! SOURCE
 
subroutine strain_strain2def(mat_strain,strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_strain2def'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 real(dp),intent(out) :: mat_strain(3,3)
 type(strain_type),intent(in) :: strain 
!Local variables-------------------------------
!scalar
 integer :: i
!arrays
! *************************************************************************

  mat_strain(:,:) = zero
  forall(i=1:3)mat_strain(i,i)=1  
 
  if(strain%name == "uniaxial") then
    select case(strain%direction)
      case(1)
        mat_strain(1,1) = mat_strain(1,1) + strain%delta
      case(2)
        mat_strain(2,2) = mat_strain(2,2) + strain%delta
      case(3)
        mat_strain(3,3) = mat_strain(3,3) + strain%delta
    end select
  else
    if (strain%name == "shear") then
    select case(strain%direction)
      case(4)
        mat_strain(3,2) =  strain%delta / 2
        mat_strain(2,3) =  strain%delta / 2
      case(5)
        mat_strain(3,1) =  strain%delta / 2
        mat_strain(1,3) =  strain%delta / 2
      case(6)
        mat_strain(2,1) =  strain%delta / 2
        mat_strain(1,2) =  strain%delta / 2
      end select
    end if
  end if
end subroutine strain_strain2def
!!***

!****f* m_strain/strain_print
!!
!! NAME
!! strain_print
!!
!! FUNCTION
!! print the structure strain
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
 
 subroutine strain_print(strain)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!array
 type(strain_type),intent(in) :: strain 
!Local variables-------------------------------
!scalar
 character(len=500) :: message
!arrays
! *************************************************************************

 if(strain%name == "reference") then 
   write(message,'(a,a,a)') ch10,' no strain found:',&
&  ' This structure is equivalent to the reference structure'
   call wrtout(std_out,message,'COLL')
 else
   if(strain%name /= "") then 
     write(message,'(a,a,a,a,a,a,I2,a,(ES10.2),a)') ch10,' strain found:',ch10,&
&      ' The strain is ',trim(strain%name),' type in the direction ',&
&      strain%direction,' with delta of ',strain%delta, ':'
     call wrtout(std_out,message,'COLL')
     write(std_out,'(3(1x,es12.2))') strain%strain
     write(ab_out,'(3(1x,es12.2))') strain%strain

   else

     write(message,'(a,a,a)') ch10,' unable to find strain:'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     write(std_out,'(3(1x,es12.2))') strain%strain
     write(ab_out,'(3(1x,es12.2))') strain%strain

   end if
 end if
end subroutine strain_print
!!***

end module m_strain

