!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_polynomial_coeff
!!
!! NAME
!! m_polynomial_coeff
!!
!! FUNCTION
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

module m_polynomial_coeff

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_polynomial_term
 use m_xmpi

 implicit none

 public :: polynomial_coeff_broacast
 public :: polynomial_coeff_free
 public :: polynomial_coeff_init
 public :: polynomial_coeff_setCoefficient
 public :: polynomial_coeff_writeXML
!!***

!!****t* m_polynomial_coeff/polynomial_coeff_type
!! NAME
!! polynomial_coeff_type
!!
!! FUNCTION
!! structure for a polynomial coefficient
!! contains the value of the coefficient and a 
!! list of terms (displacement) relating to the coefficient
!!
!! SOURCE

 type, public :: polynomial_coeff_type

   character(200) :: name
!     Name of the polynomial_coeff (Sr_y-O1_y)^3) for example

   integer :: nterm
!     Number of terms (short range interaction) for this polynomial_coeff

   real(dp) :: coefficient
!     coefficient = value of the coefficient of this term
!     \frac{\partial E^{k}}{\partial \tau^{k}}

   type(polynomial_term_type),dimension(:),allocatable :: terms
!     terms(nterm)
!     contains all the displacements for this coefficient      

 end type polynomial_coeff_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_polynomial_coeff/polynomial_coeff_init
!!
!! NAME
!! polynomial_coeff_init
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, qpoint chosen, and atoms
!!
!! INPUTS
!!  name     = Name of the polynomial_coeff (Sr_y-O1_y)^3) for example
!!  nterm   = Number of terms (short range interaction) for this polynomial_coeff
!!  coefficient  = Value of the coefficient of this term
!!  termstype(nterm) = Polynomial_term_type contains all the displacements for this coefficient 
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be initialized
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

subroutine polynomial_coeff_init(coefficient,name,nterm,polynomial_coeff,terms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nterm
 real(dp),intent(in) :: coefficient
!arrays
 character(200),intent(in) :: name
 type(polynomial_term_type),intent(in) :: terms(nterm)
 type(polynomial_coeff_type), intent(out) :: polynomial_coeff
!Local variables-------------------------------
!scalar
 integer :: ii
!arrays

! *************************************************************************
 
!First free before initilisation
 call polynomial_coeff_free(polynomial_coeff)

!Initilisation
 polynomial_coeff%name = name
 polynomial_coeff%nterm = nterm
 polynomial_coeff%coefficient = coefficient 
 ABI_DATATYPE_ALLOCATE(polynomial_coeff%terms,(polynomial_coeff%nterm))
 do ii = 1,polynomial_coeff%nterm
   call polynomial_term_init(terms(ii)%atindx,terms(ii)%cell,terms(ii)%direction,terms(ii)%ndisp,&
&                            polynomial_coeff%terms(ii),terms(ii)%power,terms(ii)%weight) 
 end do

end subroutine polynomial_coeff_init
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_free
!!
!! NAME
!! polynomial_coeff_free
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!
!!
!! CHILDREN
!!   
!!
!! SOURCE

subroutine polynomial_coeff_free(polynomial_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(polynomial_coeff_type), intent(inout) :: polynomial_coeff
!Local variables-------------------------------
!scalar
 integer :: ii
!arrays

! *************************************************************************

 if(allocated(polynomial_coeff%terms))then
   do ii = 1,polynomial_coeff%nterm
     call polynomial_term_free(polynomial_coeff%terms(ii))
   end do
   ABI_DATATYPE_DEALLOCATE(polynomial_coeff%terms)
 end if
 polynomial_coeff%name = ""
 polynomial_coeff%nterm = zero
 polynomial_coeff%coefficient = zero

end subroutine polynomial_coeff_free
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_setCoefficient
!!
!! NAME
!! polynomial_coeff_setCoefficient
!!
!! FUNCTION
!! set the coefficient for this  polynomial_coeff type
!!
!! INPUTS
!! coefficient = coefficient of this coefficient 
!! 
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!
!!
!! CHILDREN
!!   
!!
!! SOURCE

subroutine polynomial_coeff_setCoefficient(coefficient,polynomial_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_setCoefficient'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: coefficient
!arrays
 type(polynomial_coeff_type), intent(inout) :: polynomial_coeff
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************

 polynomial_coeff%coefficient = coefficient

end subroutine polynomial_coeff_setCoefficient
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_broacast
!! NAME
!! polynomial_coeff_broacast
!!
!! FUNCTION
!!  MPI broadcast all types for the polynomial_coefficent structure
!!
!! INPUTS
!!   master=Rank of Master
!!   comm=MPI communicator
!!
!! SIDE EFFECTS
!!   coefficent<type(polynomial_coefficent_type)>= Input if node is master, 
!!                              other nodes returns with a completely initialized instance.
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_coeff_broacast(coefficients, master, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_broacast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(polynomial_coeff_type),intent(inout) :: coefficients
 integer, intent(in) :: master,comm

!Local variables-------------------------------
!scalars
 integer :: ierr,ii
!arrays

! *************************************************************************

 if (xmpi_comm_size(comm) == 1) return

 DBG_ENTER("COLL")

 ! Transmit variables
  call xmpi_bcast(coefficients%name, master, comm, ierr)
  call xmpi_bcast(coefficients%nterm, master, comm, ierr)
  call xmpi_bcast(coefficients%coefficient, master, comm, ierr)
 
 !Allocate arrays on the other nodes.
  if (xmpi_comm_rank(comm) /= master) then
    ABI_DATATYPE_ALLOCATE(coefficients%terms,(coefficients%nterm))
    do ii=1,1,coefficients%nterm
      call polynomial_term_free(coefficients%terms(ii))
    end do
  end if
! Set the number of term on each node (needed for allocations of array)
  do ii = 1,coefficients%nterm
    call xmpi_bcast(coefficients%terms(ii)%ndisp, master, comm, ierr)
  end do

! Allocate arrays on the other nodes
  if (xmpi_comm_rank(comm) /= master) then
    do ii = 1,coefficients%nterm
      ABI_ALLOCATE(coefficients%terms(ii)%atindx,(2,coefficients%terms(ii)%ndisp))
      coefficients%terms(ii)%atindx = zero
      ABI_ALLOCATE(coefficients%terms(ii)%direction,(coefficients%terms(ii)%ndisp))
      ABI_ALLOCATE(coefficients%terms(ii)%cell,(3,2,coefficients%terms(ii)%ndisp))      
      ABI_ALLOCATE(coefficients%terms(ii)%power,(coefficients%terms(ii)%ndisp))
    end do
  end if

! Transfert value
  do ii = 1,coefficients%nterm
      call xmpi_bcast(coefficients%terms(ii)%weight, master, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%atindx, master, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%direction, master, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%cell, master, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%power, master, comm, ierr)
  end do

 DBG_EXIT("COLL")

end subroutine polynomial_coeff_broacast
!!***

!!****f*m_polynomial_coeff/polynomial_coeff_writeXML
!! NAME
!! polynomial_coeff_writeXML
!!
!! FUNCTION
!! This routine print the coefficents into xml format
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
!! coeffs   = array with all the coefficiens
!! ncoeff   = number of coeffs to print
!!
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_coeff_writeXML(coeffs,ncoeff,filename)

 use defs_basis
 use m_errors
 use m_io_tools, only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_writeXML'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: ncoeff
!arrays
  type(polynomial_coeff_type), intent(in) :: coeffs(ncoeff)
  character(len=fnlen),optional,intent(in) :: filename
!Local variables-------------------------------
!scalar
 integer :: icoeff,idisp,iterm
 integer :: unit_xml=22
 character(len=500) :: message
 character(len=fnlen) :: namefile
 character(len=1)  :: direction
!arrays

! *************************************************************************

!Check the inputs 
  if (size(coeffs) /= ncoeff) then
    write(message,'(a,a)')' The number of coeffs does not correspond to ncoeff' 
    MSG_ERROR(message)
  end if

!Print the coefficients into XML file
  if(ncoeff>0)then
!   convert natom in character
    if(present(filename)) then
      namefile=trim(filename)
    else
      namefile='coefficients.xml'
    end if

    call isfile(namefile,'new')

    if (open_file(namefile,message,unit=unit_xml,form="formatted",&
&       status="new",action="write") /= 0) then
      MSG_ERROR(message)
    end if

    write(message,'(a,a,a)')ch10,&
&       ' Generation of the xml file for the polynomial fitting in ',namefile

    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')

!   Write header
    WRITE(unit_xml,'("<?xml version=""1.0"" ?>")')
    WRITE(unit_xml,'("<Heff_definition>")')

!   Close header
    do icoeff = 1, ncoeff
      WRITE(unit_xml,'("  <coefficient number=""",i3,""" text=""",a,""">")') &
        icoeff,trim(coeffs(icoeff)%name)
      do iterm = 1,coeffs(icoeff)%nterm
        WRITE(unit_xml,'("    <term weight=""",F9.6,""">")') &
          coeffs(icoeff)%terms(iterm)%weight
        do idisp=1,coeffs(icoeff)%terms(iterm)%ndisp           
          if(coeffs(icoeff)%terms(iterm)%direction(idisp) < 0) then
!           Strain case
            WRITE(unit_xml,'("      <strain power=""",i2,""" voigt=""",i2,"""/>")')&
&                   coeffs(icoeff)%terms(iterm)%power(idisp),&
&                   -1 * coeffs(icoeff)%terms(iterm)%direction(idisp)
          else
!           Atomic displacement case
            select case(coeffs(icoeff)%terms(iterm)%direction(idisp))
            case(1)
              direction ="x"
            case(2)
              direction ="y"
            case(3)
              direction ="z"
            end select
            WRITE(unit_xml,'("      <displacement_diff atom_a=""",i2,""" atom_b=""",i2,&
&                            """ direction=""",a,""" power=""",i2,""">")')&
                    coeffs(icoeff)%terms(iterm)%atindx(1,idisp)-1,&
&                   coeffs(icoeff)%terms(iterm)%atindx(2,idisp)-1,direction,&
&                   coeffs(icoeff)%terms(iterm)%power(idisp)
            WRITE(unit_xml,'("        <cell_a>")',advance='no')
            WRITE(unit_xml,'(3(I4))',advance='no')&
&           coeffs(icoeff)%terms(iterm)%cell(:,1,idisp)
            WRITE(unit_xml,'("</cell_a>")')
            WRITE(unit_xml,'("        <cell_b>")',advance='no')
            WRITE(unit_xml,'(3(I4))',advance='no')&
&                coeffs(icoeff)%terms(iterm)%cell(:,2,idisp)
            WRITE(unit_xml,'("</cell_b>")')
            WRITE(unit_xml,'("      </displacement_diff>")')
          end if
        end do
        WRITE(unit_xml,'("    </term>")') 
      end do
      WRITE(unit_xml,'("  </coefficient>")') 
    end do
    WRITE(unit_xml,'("</Heff_definition>")')
!   Close file
    CLOSE(unit_xml)
  end if

end subroutine polynomial_coeff_writeXML
!!***

end module m_polynomial_coeff
!!***
