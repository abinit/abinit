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
 public :: polynomial_coeff_dot
 public :: polynomial_coeff_free
 public :: polynomial_coeff_init
 public :: polynomial_coeff_getName
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

   character(len=100) :: name = ""
!     Name of the polynomial_coeff (Sr_y-O1_y)^3) for example

   integer :: nterm = zero
!     Number of terms (short range interaction) for this polynomial_coeff

   real(dp) :: coefficient = zero
!     coefficient = value of the coefficient of this term
!     \frac{\partial E^{k}}{\partial \tau^{k}}

   type(polynomial_term_type),dimension(:),allocatable :: terms
!     terms(nterm)
!     contains all the displacements for this coefficient      

 end type polynomial_coeff_type
!!***

 interface operator (==)
   module procedure coeffs_compare
 end interface operator (==)

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
!!      m_anharmonics_terms,m_effective_potential_file
!!
!! CHILDREN
!!      isfile,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_init(coefficient,nterm,polynomial_coeff,terms,name)


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
 character(len=100),optional,intent(in) :: name
 type(polynomial_term_type),intent(in) :: terms(nterm)
 type(polynomial_coeff_type), intent(out) :: polynomial_coeff
!Local variables-------------------------------
!scalar
 integer :: iterm1,iterm2
 integer :: ii,nterm_tmp
 real(dp):: coefficient_tmp
!arrays
 real(dp) :: weights(nterm)
 character(len=100) :: name_tmp
! *************************************************************************
 
!First free before initilisation
 call polynomial_coeff_free(polynomial_coeff)

!Check if the list of term is available or contains identical terms
!in this case, remove all the not needed terms
 nterm_tmp = 0
 weights(:) = one
 do iterm1=1,nterm
   if(weights(iterm1)==0)cycle
   weights(iterm1) = terms(iterm1)%weight
   do iterm2=iterm1+1,nterm
     if(weights(iterm2)==0)cycle     
!    if the terms are identical we check the weight
     if(terms(iterm1)==terms(iterm2))then
       weights(iterm1) = weights(iterm1) + terms(iterm2)%weight
       weights(iterm2) = 0
     end if
   end do
   if(weights(iterm1)/=0) weights(iterm1)= anint(weights(iterm1)/weights(iterm1))
 end do

!Count the number of terms
 nterm_tmp = 0
 do iterm1=1,nterm
   if(weights(iterm1) /= 0) nterm_tmp = nterm_tmp + 1
 end do
 
 if (nterm_tmp ==0)then
   coefficient_tmp = 0.0
 else
   coefficient_tmp = coefficient
 end if
 
 if(present(name))then
   name_tmp = name
 else
   name_tmp = ""
 end if

!Initilisation
 polynomial_coeff%name = name_tmp
 polynomial_coeff%nterm = nterm_tmp
 polynomial_coeff%coefficient = coefficient_tmp
 ABI_DATATYPE_ALLOCATE(polynomial_coeff%terms,(polynomial_coeff%nterm))
 iterm1 = 0
 do ii = 1,nterm
   if(weights(ii)/= 0)then
     iterm1 = iterm1 + 1
     call polynomial_term_init(terms(ii)%atindx,terms(ii)%cell,terms(ii)%direction,terms(ii)%ndisp,&
&                              polynomial_coeff%terms(iterm1),terms(ii)%power,terms(ii)%weight) 
   end if
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
!!      m_anharmonics_terms,m_effective_potential_file,m_polynomial_coeff
!!
!! CHILDREN
!!      isfile,wrtout
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
!!      m_effective_potential_file
!!
!! CHILDREN
!!      isfile,wrtout
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


!!****f* m_polynomial_coeff/polynomial_coeff_getName
!!
!! NAME
!! polynomial_coeff_getName
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
!!      m_effective_potential_file
!!
!! CHILDREN
!!      isfile,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_getName(name,atm1,atm2,dir,power,polynomial_coeff,cell_atm1,cell_atm2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_getName'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 character(len=1),optional,intent(in) :: dir,power
 character(len=5),optional,intent(in) :: atm1,atm2
 character(len=100),optional,intent(out):: name
 type(polynomial_coeff_type),optional, intent(in) :: polynomial_coeff
 integer,optional,intent(in) :: cell_atm1(3),cell_atm2(3)
!Local variables-------------------------------
!scalar
!arrays
 character(len=100):: atm1_tmp,atm2_tmp

! *************************************************************************
 name=""
 if (present(polynomial_coeff)) then
   name = polynomial_coeff%name
 else
   if (present(dir).and.present(power).and.present(atm1).and.present(atm2)) then
     if(present(cell_atm1))then
       if (any(cell_atm1(:) /= zero) )then
         write(atm1_tmp,'(4a,I0,a,I0,a,I0,a)')  trim(atm1),"_",dir,"[",cell_atm1(1)," ",&
&                                      cell_atm1(2)," ",cell_atm1(3),"]"
       else
         atm1_tmp = trim(atm1)//"_"//dir
       end if
     else
       atm1_tmp = trim(atm1)//"_"//dir
     end if
     if(present(cell_atm2))then
       if(any(cell_atm2(:) /= zero))then
         write(atm2_tmp,'(4a,I0,a,I0,a,I0,a)')  trim(atm2),"_",dir,"[",cell_atm2(1)," ",&
&                                      cell_atm2(2)," ",cell_atm2(3),"]"
       else
         atm2_tmp = trim(atm2)//"_"//dir
       end if
     else
       atm2_tmp = trim(atm2)//"_"//dir
     end if
     name="("//trim(atm1_tmp)//"-"//trim(atm2_tmp)//")^"//power
   end if
 end if


end subroutine polynomial_coeff_getName
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
!!      m_effective_potential_file
!!
!! CHILDREN
!!      isfile,wrtout
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
!!      m_effective_potential
!!
!! CHILDREN
!!      isfile,wrtout
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
      WRITE(unit_xml,'("  <coefficient number=""",I0,""" text=""",a,""">")') &
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
            WRITE(unit_xml,'("      <displacement_diff atom_a=""",I0,""" atom_b=""",I0,&
&                            """ direction=""",a,""" power=""",I0,""">")')&
                    coeffs(icoeff)%terms(iterm)%atindx(1,idisp)-1,&
&                   coeffs(icoeff)%terms(iterm)%atindx(2,idisp)-1,direction,&
&                   coeffs(icoeff)%terms(iterm)%power(idisp)
            WRITE(unit_xml,'("        <cell_a>")',advance='no')
            WRITE(unit_xml,'(3(I0,a,I0,a,I0))',advance='no')&
&           coeffs(icoeff)%terms(iterm)%cell(1,1,idisp)," ",&
&           coeffs(icoeff)%terms(iterm)%cell(2,1,idisp)," ",&
&           coeffs(icoeff)%terms(iterm)%cell(3,1,idisp)
            WRITE(unit_xml,'("</cell_a>")')
            WRITE(unit_xml,'("        <cell_b>")',advance='no')
            WRITE(unit_xml,'(3(I0,a,I0,a,I0))',advance='no')&
&           coeffs(icoeff)%terms(iterm)%cell(1,2,idisp)," ",&
&           coeffs(icoeff)%terms(iterm)%cell(2,2,idisp)," ",&
&           coeffs(icoeff)%terms(iterm)%cell(3,2,idisp)
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

!!****f* m_polynomial_coeff/polynomial_coeff_dot
!! NAME
!!  polynomial_coeff_dot
!!
!! FUNCTION
!!  return the multiplication of two coefficients
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polynomial_coeff_dot(coeff1_in,coeff2_in,coeffs_out,natom,power,symbols)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_dot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalar
 integer, intent(in) :: natom,power
!arrays
 type(polynomial_coeff_type), intent(in) :: coeff1_in,coeff2_in
 type(polynomial_coeff_type), intent(out) :: coeffs_out(power)
 character(len=5), intent(in) :: symbols(natom)
!local variables------------------------------
!scalar
 integer :: idisp1,ipower,iterm1,iterm2
 integer :: jpower,kpower,mu,ndisp,nterm_max
 real(dp):: weight
 logical :: res
!arrays
 integer,allocatable :: atindx(:,:),cell(:,:,:),dir_int(:),powers(:)
 character(len=1) :: powerchar
 character(len=1) :: mutodir(3) = (/"x","y","z"/)
 character(len=100):: name,text
 type(polynomial_term_type),dimension(:),allocatable :: terms

! *************************************************************************
  res = .false.


! Get the number of term for the new coefficient
! and allocate the new array terms

! case 1 : the two coefficients are identical
  if(coeff1_in==coeff2_in)then

    do ipower=1,power
      do jpower=1,ipower
        
        kpower = power - ipower + 1

        nterm_max = coeff1_in%nterm
        ABI_DATATYPE_ALLOCATE(terms,(nterm_max))

        do mu=1,3
          ndisp = coeff1_in%terms(1)%ndisp

!       Allocation of the new array for this term
        ABI_ALLOCATE(atindx,(2,ndisp))
        ABI_ALLOCATE(cell,(3,2,ndisp))
        ABI_ALLOCATE(dir_int,(ndisp))
        ABI_ALLOCATE(powers,(ndisp))
    
        do iterm1=1,coeff1_in%nterm
          do iterm2=iterm1,coeff1_in%nterm
!       1-copy the displacement from the first coefficient
            do idisp1=1,coeff1_in%terms(iterm1)%ndisp
              atindx(:,idisp1) = coeff1_in%terms(iterm1)%atindx(:,idisp1)
              cell(:,:,idisp1) = coeff1_in%terms(iterm1)%cell(:,:,idisp1)
              dir_int( idisp1) = coeff1_in%terms(iterm1)%direction(idisp1)
              powers( idisp1)  = coeff1_in%terms(iterm1)%power(idisp1) +&
&                                coeff2_in%terms(iterm1)%power(idisp1)
              weight           = coeff1_in%terms(iterm1)%weight
            end do!end loop disp
            call polynomial_term_init(atindx,cell,dir_int,ndisp,terms(iterm1),powers,weight)
          end do
        end do
      end do

!   Deallocation of the  array
        ABI_DEALLOCATE(atindx)
        ABI_DEALLOCATE(cell)
        ABI_DEALLOCATE(dir_int)
        ABI_DEALLOCATE(powers)

        name=''
        do idisp1=1,terms(1)%ndisp
          write(powerchar,'(I0)') terms(1)%power(idisp1)
          call polynomial_coeff_getName(text,&
&                                      atm1=symbols(terms(1)%atindx(1,idisp1)),&
&                                      atm2=symbols(terms(1)%atindx(2,idisp1)),&
&                                      dir=mutodir(terms(1)%direction(idisp1)),&
&                                      power=trim(powerchar))
          
          name = trim(name)//trim(text)
        end do
        call polynomial_coeff_init(one,nterm_max,coeffs_out(1),terms,name=name)
    
        do iterm1=1,nterm_max
          call polynomial_term_free(terms(iterm1))
        end do
    
        ABI_DEALLOCATE(terms)
        
      end do
    end do
  end if

end subroutine polynomial_coeff_dot
!!***


!!****f* m_polynomial_coeff/coeffs_compare
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

pure function coeffs_compare(c1,c2) result (res)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'coeffs_compare'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  type(polynomial_coeff_type), intent(in) :: c1,c2
  logical :: res
!local
  integer :: iterm1,iterm2
! *************************************************************************
  res = .false.

  do iterm1=1,c1%nterm
    do iterm2=1,c2%nterm
      if(c1%terms(iterm1)==c2%terms(iterm2)) then
        res = .true.
      end if
    end do
  end do

end function coeffs_compare
!!***


end module m_polynomial_coeff
!!***
