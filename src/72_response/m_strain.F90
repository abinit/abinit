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

module m_strain

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use m_symtk,         only : matr3inv

 implicit none

 private :: strain_def2strain
 private :: strain_strain2def
 public  :: strain_print
 public  :: strain_get
 public  :: strain_init
 public  :: strain_apply
!!***

!!****t* m_strain/strain_type
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
!!      compute_anharmonics,m_effective_potential
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine strain_init(strain,delta,direction,name)

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
   strain%direction = 0
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
!!      compute_anharmonics
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine strain_free(strain)

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
 strain%direction = 0
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
!!  symmetrized = (optional) symmetrize the output
!!
!! OUTPUT
!!  strain = structure with all information of strain
!!
!! PARENTS
!!      compute_anharmonics,m_effective_potential,m_fit_data
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine strain_get(strain,rprim,rprim_def,mat_delta,symmetrized)

!Arguments ------------------------------------
!scalars
!array
 type(strain_type),intent(inout) :: strain
 real(dp),optional,intent(in) :: rprim(3,3),rprim_def(3,3), mat_delta(3,3)
 logical,optional,intent(in) :: symmetrized
!Local variables-------------------------------
!scalar
 integer :: i,j
 logical :: symmetrized_in
 character(len=500) :: message
!arrays
 real(dp) :: mat_delta_tmp(3,3),rprim_inv(3,3)
 real(dp) :: identity(3,3)
! *************************************************************************

!check inputs
 symmetrized_in = .FALSE.
 if(present(symmetrized)) then
   symmetrized_in = symmetrized
 end if

 if((present(rprim_def).and..not.present(rprim)).or.&
&   (present(rprim).and..not.present(rprim_def))) then
    write(message, '(a)' )&
&     ' strain_get: should give rprim_def and rprim as input of the routines'
    MSG_BUG(message)
  end if

 if(present(rprim_def).and.present(rprim))then
   mat_delta_tmp = zero
!  Fill the identity matrix
   identity = zero
   forall(i=1:3)identity(i,i)=1

   call matr3inv(rprim,rprim_inv)
   mat_delta_tmp =  matmul(rprim_def,transpose(rprim_inv))-identity
   identity = zero
   do i=1,3
     do j=1,3
       if (abs(mat_delta_tmp(i,j))>tol10) then
         identity(i,j) = mat_delta_tmp(i,j)
       end if
     end do
   end do

   mat_delta_tmp = identity

 else if (present(mat_delta)) then
   mat_delta_tmp = mat_delta

 else
   write(message, '(a)' )&
&     ' strain_get: should give rprim_def or mat_delta as input of the routines'
   MSG_BUG(message)
 end if

 if(symmetrized_in)then
   mat_delta_tmp(2,3) = (mat_delta_tmp(2,3) + mat_delta_tmp(3,2)) / 2
   mat_delta_tmp(3,1) = (mat_delta_tmp(3,1) + mat_delta_tmp(1,3)) / 2
   mat_delta_tmp(2,1) = (mat_delta_tmp(2,1) + mat_delta_tmp(1,2)) / 2

   mat_delta_tmp(3,2) = mat_delta_tmp(2,3)
   mat_delta_tmp(1,3) = mat_delta_tmp(3,1)
   mat_delta_tmp(1,2) = mat_delta_tmp(2,1)

 end if

 call strain_def2strain(mat_delta_tmp,strain)

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
!!   anharmonic_terms_compute,
!!
!! CHILDREN
!!
!! SOURCE

subroutine strain_apply(rprim,rprim_def,strain)

!Arguments ------------------------------------
!scalars
!array
 real(dp),intent(in)  :: rprim(3,3)
 real(dp),intent(out) :: rprim_def(3,3)
 type(strain_type),intent(in) :: strain
!Local variables-------------------------------
!scalar
 !integer :: i
!arrays
! *************************************************************************

 rprim_def(:,:) = zero
! Fill the identity matrix
 rprim_def(:,:) = matmul(strain%strain(:,:),transpose(rprim(:,:)))

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
!!   multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine strain_def2strain(mat_strain,strain)

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
 strain%direction = 0
 strain%strain = mat_strain

 if (all(abs(mat_strain)<tol10)) then
   strain%name = "reference"
   strain%delta = zero
   strain%direction = 0
   strain%strain = zero
 else
   if(abs(mat_strain(1,1))>tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))<tol10.and.&
&   abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))>tol10.and.abs(mat_strain(2,3))<tol10.and.&
&   abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))>tol10) then
     if((mat_strain(1,1)-mat_strain(2,2))< tol10.and.&
&       (mat_strain(1,1)-mat_strain(3,3))< tol10) then
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
      if (abs(mat_strain(3,2)-mat_strain(3,2))<tol10) then
        strain%name = "shear"
        strain%delta = mat_strain(3,2) * 2
        strain%direction = 4
        strain%strain = mat_strain
      end if
    end if
    if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))<tol10.and.abs(mat_strain(1,3))>tol10.and.&
&    abs(mat_strain(2,1))<tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))>tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))<tol10) then
      if (abs(mat_strain(3,1)-mat_strain(1,3))<tol10) then
        strain%name = "shear"
        strain%delta = mat_strain(3,1) * 2
        strain%direction = 5
        strain%strain = mat_strain
      end if
    end if
    if(abs(mat_strain(1,1))<tol10.and.abs(mat_strain(1,2))>tol10.and.abs(mat_strain(1,3))<tol10.and.&
&    abs(mat_strain(2,1))>tol10.and.abs(mat_strain(2,2))<tol10.and.abs(mat_strain(2,3))<tol10.and.&
&    abs(mat_strain(3,1))<tol10.and.abs(mat_strain(3,2))<tol10.and.abs(mat_strain(3,3))<tol10) then
      if (abs(mat_strain(1,2)-mat_strain(2,1))<tol10) then
        strain%name = "shear"
        strain%delta = mat_strain(2,1) * 2
        strain%direction = 6
        strain%strain = mat_strain
      end if
    end if
  end if

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
!!   multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine strain_strain2def(mat_strain,strain)

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

 select case(strain%direction)
 case(1)
   mat_strain(1,1) = mat_strain(1,1) + strain%delta
 case(2)
   mat_strain(2,2) = mat_strain(2,2) + strain%delta
 case(3)
   mat_strain(3,3) = mat_strain(3,3) + strain%delta
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
!!      m_effective_potential
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine strain_print(strain)

!Arguments ------------------------------------
!scalars
!array
 type(strain_type),intent(in) :: strain
!Local variables-------------------------------
!scalar
 integer :: ii
 character(len=500) :: message
!arrays
! *************************************************************************

 if(strain%name == "reference") then
   write(message,'(4a)') ch10,' no strain found:',&
&  ' This structure is equivalent to the reference structure',ch10
   call wrtout(std_out,message,'COLL')
 else
   if(strain%name /= "") then
     write(message,'(3a,I2,a,(ES10.2),a)') &
&      ' The strain is ',trim(strain%name),' type in the direction ',&
&      strain%direction,' with delta of ',strain%delta, ':'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     do ii = 1,3
       write(message,'(3es17.8)') strain%strain(ii,1),strain%strain(ii,2),strain%strain(ii,3)
       call wrtout(std_out,message,'COLL')
     end do
   else
     write(message,'(a)') ' Strain does not correspond to standard strain:'
     call wrtout(std_out,message,'COLL')
     do ii = 1,3
       write(message,'(3es17.8)') strain%strain(ii,1),strain%strain(ii,2),strain%strain(ii,3)
       call wrtout(std_out,message,'COLL')
     end do
   end if
 end if
end subroutine strain_print
!!***

end module m_strain

