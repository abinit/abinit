!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_polynomial_coeff
!!
!! NAME
!! m_polynomial_coeff
!!
!! FUNCTION
!! Module with the datatype polynomial coefficients
!!
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

module m_polynomial_coeff

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_polynomial_term
 use m_phonon_supercell, only: getPBCIndexes_supercell
 use m_xmpi

 implicit none

 public :: polynomial_coeff_broadcast
 public :: polynomial_coeff_evaluate
 public :: polynomial_coeff_free
 public :: polynomial_coeff_init
 public :: polynomial_coeff_getName
 public :: polynomial_coeff_MPIrecv
 public :: polynomial_coeff_MPIsend
 public :: polynomial_coeff_setName
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
!! list of terms (displacements and/or strain) relating to the coefficient by symmetry
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
!     polynomial_term(nterm)<type(polynomial_term)> 
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
!! Initialize a polynomial_coeff datatype
!!
!! INPUTS
!!  name     = Name of the polynomial_coeff (Sr_y-O1_y)^3) for example
!!  nterm    = Number of terms (short range interaction) for this polynomial_coeff
!!  coefficient  = Value of the coefficient of this term
!!  terms(nterm)<type(polynomial_term)> = array of polynomial_term_type 
!!  check   = TRUE if this list of terms has to be check. We remove the symetric of equivalent terms
!!                  for example:  ((Sr_y-O1_y)^1 and -1*(Sr_y-O1_y)^1 => zero
!!            FALSE, defaut, do nothing
!! OUTPUT
!!   polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype to be initialized
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential_file,m_fit_polynomial_coeff
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_getname,polynomial_coeff_init,polynomial_term_free
!!      polynomial_term_init
!!
!! SOURCE

subroutine polynomial_coeff_init(coefficient,nterm,polynomial_coeff,terms,name,check)


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
 logical,optional,intent(in) :: check
!arrays
 character(len=100),optional,intent(in) :: name
 type(polynomial_term_type),intent(in) :: terms(nterm)
 type(polynomial_coeff_type), intent(out) :: polynomial_coeff
!Local variables-------------------------------
!scalar
 integer :: iterm1,iterm2
 integer :: ii,nterm_tmp
 real(dp):: coefficient_tmp
 logical :: check_in = .false.
!arrays
 real(dp) :: weights(nterm)
 character(len=100) :: name_tmp
! *************************************************************************
 
!First free before initilisation
 call polynomial_coeff_free(polynomial_coeff)

 if(present(check)) check_in = check
 if(check_in)then
!  Check if the list of term is available or contains identical terms
!  in this case, remove all the not needed terms
   nterm_tmp = 0
   weights(:) = one
   do iterm1=1,nterm
     if(weights(iterm1)==0)cycle
     weights(iterm1) = terms(iterm1)%weight
     do iterm2=iterm1+1,nterm
       if(weights(iterm2)==0)cycle     
!      if the terms are identical we check the weight
       if(terms(iterm1)==terms(iterm2))then
         weights(iterm1) = weights(iterm1) + terms(iterm2)%weight
         weights(iterm2) = 0
       end if
     end do
     if(weights(iterm1)/=0) weights(iterm1)= anint(weights(iterm1)/weights(iterm1))
   end do

!  Count the number of terms
   nterm_tmp = 0
   do iterm1=1,nterm
     if(weights(iterm1) /= 0) nterm_tmp = nterm_tmp + 1
   end do
 
   if (nterm_tmp ==0)then
     coefficient_tmp = 0.0
   else
     coefficient_tmp = coefficient
   end if

 else
   nterm_tmp = nterm
   coefficient_tmp = coefficient
   weights(:) = terms(:)%weight
 end if!end Check

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
!! Free polynomial_coeff datatype
!!
!! INPUTS
!! polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype 
!!
!! OUTPUT
!! polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype
!!
!! PARENTS
!!      m_anharmonics_terms,m_effective_potential_file,m_fit_polynomial_coeff
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_getname,polynomial_coeff_init,polynomial_term_free
!!      polynomial_term_init
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
!! set the coefficient for of polynomial_coeff
!!
!! INPUTS
!! coefficient = coefficient of this coefficient 
!! 
!! OUTPUT
!! polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      polynomial_coeff_getname,polynomial_coeff_init,polynomial_term_free
!!      polynomial_term_init
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

!!****f* m_polynomial_coeff/polynomial_coeff_setName
!!
!! NAME
!! polynomial_coeff_setName
!!
!! FUNCTION
!! set the name of a  polynomial_coeff type
!!
!! INPUTS
!! name = name of the coeff
!! 
!! OUTPUT
!! polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_coeff_setName(name,polynomial_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_setName'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 character(len=200),intent(in) :: name
 type(polynomial_coeff_type), intent(inout) :: polynomial_coeff
!Local variables-------------------------------
!scalar
!arrays
! *************************************************************************

 polynomial_coeff%name = name

end subroutine polynomial_coeff_setName
!!***


!!****f* m_polynomial_coeff/polynomial_coeff_getName
!!
!! NAME
!! polynomial_coeff_getName
!!
!! FUNCTION
!! get the name of a polynomial coefficient
!!
!! INPUTS
!! natom = number of atoms
!! polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype
!! symbols(natom)  =  array with the atomic symbol:["Sr","Ru","O1","O2","O3"]
!! recompute = (optional) flag to set if the name has to be recomputed
!! iterm = (optional) number of the term used for the name
!! 
!! OUTPUT
!! name = name xof the coefficients
!!
!! PARENTS
!!      m_fit_polynomial_coeff,m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_getname,polynomial_coeff_init,polynomial_term_free
!!      polynomial_term_init
!!
!! SOURCE

subroutine polynomial_coeff_getName(name,natom,polynomial_coeff,symbols,recompute,iterm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_getName'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 integer,optional,intent(in) :: iterm
!arrays
 character(len=5),intent(in) :: symbols(:)
 character(len=100),intent(out):: name
 type(polynomial_coeff_type),optional, intent(in) :: polynomial_coeff
 logical,optional,intent(in) :: recompute
!Local variables-------------------------------
!scalar
 integer :: ii,idisp,iterm_in
 logical :: need_recompute = .FALSE.
!arrays
 integer :: cell_atm1(3),cell_atm2(3)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=1) :: dir
 character(len=2) :: power,powerchar
 character(len=5) :: atm1,atm2
 character(len=100):: atm1_tmp,atm2_tmp
 character(len=200):: text
 character(len=500):: message
! *************************************************************************

!Reset output
 name=""
 iterm_in = 1

!Set the optional arguments 
 if(present(recompute)) need_recompute = recompute
 if(present(iterm)) then
   iterm_in = iterm
 else
   if(need_recompute)then
     iterm_in = -1
     do ii=1,polynomial_coeff%nterm
!      Find the index of the ref 
       if(iterm_in==-1) then !Need to find the reference term
         do idisp=1,polynomial_coeff%terms(ii)%ndisp
           if(polynomial_coeff%terms(ii)%direction(idisp) > zero) then
             iterm_in = ii
             if(any(polynomial_coeff%terms(ii)%cell(:,1,idisp) /= zero).or.&
&               any(polynomial_coeff%terms(ii)%cell(:,2,idisp) /= zero)) then
               iterm_in = -1
               exit
             end if
           end if
         end do!end do disp
       else
         exit             
       end if
     end do!end do term
!    If not find, we set to the first element
     if(iterm_in==-1) iterm_in = 1
   else
     iterm_in = 1
   end if
 end if

!Do check
 if(iterm_in > polynomial_coeff%nterm.or.iterm_in < zero) then
   write(message, '(5a)')&
&      ' The number of the requested term for the generation of',ch10,&
&      'the name of the coefficient is not possible.',ch10,&
&      'Action: Contact Abinit group.'
   MSG_BUG(message)
 end if

 if(polynomial_coeff%name /= "".and..not.need_recompute)then
   name = polynomial_coeff%name
 else
!  Nedd to recompute
   do idisp=1,polynomial_coeff%terms(iterm_in)%ndisp
     text = ""
     !Fill variables for this displacement
     write(powerchar,'(I0)') polynomial_coeff%terms(iterm_in)%power(idisp)
     power=trim(powerchar)
     if(polynomial_coeff%terms(iterm_in)%direction(idisp)>zero) then
       atm1=symbols(polynomial_coeff%terms(iterm_in)%atindx(1,idisp))
       atm2=symbols(polynomial_coeff%terms(iterm_in)%atindx(2,idisp))
       dir=mutodir(polynomial_coeff%terms(iterm_in)%direction(idisp))
       cell_atm1=polynomial_coeff%terms(iterm_in)%cell(:,1,idisp)
       cell_atm2=polynomial_coeff%terms(iterm_in)%cell(:,2,idisp)
!      Construct ATM1
       if (any(cell_atm1(:) /= zero) )then
         write(atm1_tmp,'(4a,I0,a,I0,a,I0,a)')  trim(atm1),"_",dir,"[",cell_atm1(1)," ",&
&                                               cell_atm1(2)," ",cell_atm1(3),"]"
       else
         atm1_tmp = trim(atm1)//"_"//dir
       end if
!      Construct ATM2
       if(any(cell_atm2(:) /= zero))then
         write(atm2_tmp,'(4a,I0,a,I0,a,I0,a)')  trim(atm2),"_",dir,"[",cell_atm2(1)," ",&
 &                                              cell_atm2(2)," ",cell_atm2(3),"]"
       else
         atm2_tmp = trim(atm2)//"_"//dir
       end if

       text="("//trim(atm1_tmp)//"-"//trim(atm2_tmp)//")^"//power
     else
       !Strain case
      dir=mutodir(3+abs(polynomial_coeff%terms(iterm_in)%direction(idisp)))
      text="("//"eta_"//trim(dir)//")^"//power
     end if     
     name = trim(name)//trim(text)
   end do
 end if

end subroutine polynomial_coeff_getName
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_broadcast
!! NAME
!! polynomial_coeff_broadcast
!!
!! FUNCTION
!!  MPI broadcast  polynomial_coefficent datatype
!!
!! INPUTS
!!  source = rank of source
!!  comm = MPI communicator
!!
!! SIDE EFFECTS
!!  coefficients<type(polynomial_coefficent_type)>= Input if node is source, 
!!                              other nodes returns with a completely initialized instance.
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!      polynomial_coeff_getname,polynomial_coeff_init,polynomial_term_free
!!      polynomial_term_init
!!
!! SOURCE

subroutine polynomial_coeff_broadcast(coefficients, source, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_broadcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(polynomial_coeff_type),intent(inout) :: coefficients
 integer, intent(in) :: source,comm

!Local variables-------------------------------
!scalars
 integer :: ierr,ii
!arrays

! *************************************************************************

 if (xmpi_comm_size(comm) == 1) return

! Free the output
 if (xmpi_comm_rank(comm) /= source) then
   call polynomial_coeff_free(coefficients)
 end if

 ! Transmit variables
  call xmpi_bcast(coefficients%name, source, comm, ierr)
  call xmpi_bcast(coefficients%nterm, source, comm, ierr)
  call xmpi_bcast(coefficients%coefficient, source, comm, ierr)
 
 !Allocate arrays on the other nodes.
  if (xmpi_comm_rank(comm) /= source) then
    ABI_DATATYPE_ALLOCATE(coefficients%terms,(coefficients%nterm))
    do ii=1,1,coefficients%nterm
      call polynomial_term_free(coefficients%terms(ii))
    end do
  end if
! Set the number of term on each node (needed for allocations of array)
  do ii = 1,coefficients%nterm
    call xmpi_bcast(coefficients%terms(ii)%ndisp, source, comm, ierr)
  end do

! Allocate arrays on the other nodes
  if (xmpi_comm_rank(comm) /= source) then
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
      call xmpi_bcast(coefficients%terms(ii)%weight, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%atindx, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%direction, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%cell, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%power, source, comm, ierr)
  end do


end subroutine polynomial_coeff_broadcast
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_MPIsend
!! NAME
!! polynomial_coeff_MPIsend
!!
!! FUNCTION
!!  MPI send the polynomial_coefficent datatype
!!
!! INPUTS
!!   tag = tag of the message to send
!!   dest= rank of Dest
!!   comm= MPI communicator
!!
!! SIDE EFFECTS
!!   polynomial_coeff<type(polynomial_coeff)> = polynomial_coeff datatype
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!!
!! SOURCE

subroutine polynomial_coeff_MPIsend(coefficients, tag, dest, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_MPIsend'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(polynomial_coeff_type),intent(inout) :: coefficients
 integer, intent(in) :: dest,comm,tag

!Local variables-------------------------------
!scalars
 integer :: ierr,ii
 integer :: my_rank
!arrays

! *************************************************************************

 if (xmpi_comm_size(comm) == 1) return

  my_rank = xmpi_comm_rank(comm)
! Transmit variables
  call xmpi_send(coefficients%name, dest, 9*tag+0, comm, ierr)
  call xmpi_send(coefficients%nterm, dest, 9*tag+1, comm, ierr)
  call xmpi_send(coefficients%coefficient, dest, 9*tag+2, comm, ierr)
 
! Set the number of term on each node (needed for allocations of array)
  do ii = 1,coefficients%nterm
    call xmpi_send(coefficients%terms(ii)%ndisp, dest, 9*tag+3, comm, ierr)
  end do

! Transfert value
  do ii = 1,coefficients%nterm
      call xmpi_send(coefficients%terms(ii)%weight, dest, 9*tag+4, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%atindx, dest, 9*tag+5, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%direction, dest, 9*tag+6, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%cell, dest, 9*tag+7, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%power, dest, 9*tag+8, comm, ierr)
  end do

end subroutine polynomial_coeff_MPIsend
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_MPIrecv
!! NAME
!! polynomial_coeff_MPIrecv
!!
!! FUNCTION
!!  MPI receive the polynomial_coefficent datatype
!!
!! INPUTS
!!   tag = tag of the message to receive
!!   source = rank of Source
!!   comm = MPI communicator
!!
!! SIDE EFFECTS
!!   coefficients<type(polynomial_coefficent_type)>=  polynomial_coeff datatype
!!
!! PARENTS
!!      m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine polynomial_coeff_MPIrecv(coefficients, tag, source, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_MPIrecv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(polynomial_coeff_type),intent(inout) :: coefficients
 integer, intent(in) :: source,comm,tag

!Local variables-------------------------------
!scalars
 integer :: ierr,ii
!arrays

! *************************************************************************

 if (xmpi_comm_size(comm) == 1) return


! Free the output
  call polynomial_coeff_free(coefficients)

 ! Transmit variables
  call xmpi_recv(coefficients%name, source, 9*tag+0, comm, ierr)
  call xmpi_recv(coefficients%nterm, source, 9*tag+1, comm, ierr)
  call xmpi_recv(coefficients%coefficient, source, 9*tag+2, comm, ierr)
 
 !Allocate arrays on the other nodes.
  ABI_DATATYPE_ALLOCATE(coefficients%terms,(coefficients%nterm))
  do ii=1,1,coefficients%nterm
    call polynomial_term_free(coefficients%terms(ii))
  end do

! Set the number of term on each node (needed for allocations of array)
  do ii = 1,coefficients%nterm
    call xmpi_recv(coefficients%terms(ii)%ndisp, source, 9*tag+3, comm, ierr)
  end do

! Allocate arrays on the other nodes
  do ii = 1,coefficients%nterm
    ABI_ALLOCATE(coefficients%terms(ii)%atindx,(2,coefficients%terms(ii)%ndisp))
    coefficients%terms(ii)%atindx = zero
    ABI_ALLOCATE(coefficients%terms(ii)%direction,(coefficients%terms(ii)%ndisp))
    ABI_ALLOCATE(coefficients%terms(ii)%cell,(3,2,coefficients%terms(ii)%ndisp))      
    ABI_ALLOCATE(coefficients%terms(ii)%power,(coefficients%terms(ii)%ndisp))
  end do

! Transfert value
  do ii = 1,coefficients%nterm
    call xmpi_recv(coefficients%terms(ii)%weight, source, 9*tag+4, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%atindx, source, 9*tag+5, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%direction, source, 9*tag+6, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%cell, source, 9*tag+7, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%power, source, 9*tag+8, comm, ierr)
  end do

end subroutine polynomial_coeff_MPIrecv
!!***

!!****f*m_polynomial_coeff/polynomial_coeff_writeXML
!! NAME
!! polynomial_coeff_writeXML
!!
!! FUNCTION
!! This routine print the coefficents into XML format
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! coeffs(ncoeffs)<type(polynomial_coeff)> = array of polynomial_coeff datatype
!! ncoeff = number of coeffs to print
!! filename = optional,the name of output file
!!                     default is coefficients.xml
!! unit = optional,unit of the output file
!! newfile = optional, TRUE the coefficients are print in new XML (print the headers)
!!                     FALSE (requieres unit) will not print the headers
!!
!! OUTPUT
!!
!! PARENTS
!!      m_effective_potential,m_fit_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_getname,polynomial_coeff_init,polynomial_term_free
!!      polynomial_term_init
!!
!! SOURCE

subroutine polynomial_coeff_writeXML(coeffs,ncoeff,filename,unit,newfile)

 use defs_basis
 use m_errors
 use m_io_tools, only : open_file,get_unit

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
  integer,optional,intent(in) :: unit
  logical,optional,intent(in) :: newfile
!arrays
  type(polynomial_coeff_type), intent(in) :: coeffs(ncoeff)
  character(len=fnlen),optional,intent(in) :: filename
!Local variables-------------------------------
!scalar
 integer :: icoeff,idisp,iterm
 integer :: unit_xml
 logical :: need_header = .TRUE.
 character(len=500) :: message
 character(len=fnlen) :: namefile
 character(len=1)  :: direction

!arrays

! *************************************************************************

!fill the default
 unit_xml = get_unit()

!Check the inputs 
 if(present(filename))then
   namefile=trim(filename)
 else
   namefile='coefficients.xml'
 end if

 if(present(newfile))then
   if (newfile) then
     unit_xml = get_unit()
     need_header = .TRUE.
     call isfile(namefile,'new')
   else
     if(.not.present(unit))then
       write(message,'(a,a)')' You  need to specified the unit' 
       MSG_ERROR(message)
     else
       need_header = .FALSE.
       unit_xml = unit
     end if
   end if
 end if
 if (size(coeffs) /= ncoeff) then
   write(message,'(a,a)')' The number of coeffs does not correspond to ncoeff' 
   MSG_ERROR(message)
 end if

!Print the coefficients into XML file
 if(ncoeff>0)then
   if(need_header)then
!    open new file
     if (open_file(namefile,message,unit=unit_xml,form="formatted",&
&         status="new",action="write") /= 0) then
       MSG_ERROR(message)
     end if
   else
!     just open the file to append the coefficient
     open(unit=unit_xml,file=namefile,position="append")
   end if

!  Write header
   if (need_header)then
     write(message,'(a,a,a)')ch10,&
&         ' Generation of the xml file for the fitted polynomial in ',namefile
   
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     WRITE(unit_xml,'("<?xml version=""1.0"" ?>")')
   end if
   WRITE(unit_xml,'("<Heff_definition>")')
!   Close header
    do icoeff = 1, ncoeff
      WRITE(unit_xml,'("  <coefficient number=""",I0,""" value=""",E19.10,""" text=""",a,""">")') &
         icoeff,coeffs(icoeff)%coefficient,trim(coeffs(icoeff)%name)
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
            WRITE(unit_xml,'(a,I0,a,I0,3a,I0,a)') "      <displacement_diff atom_a=""",&
&                           coeffs(icoeff)%terms(iterm)%atindx(1,idisp)-1,""" atom_b=""",&
&                           coeffs(icoeff)%terms(iterm)%atindx(2,idisp)-1,""" direction=""",&
&                           direction,""" power=""",coeffs(icoeff)%terms(iterm)%power(idisp),&
&                           """>"
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

!!****f* m_polynomial_coeff/polynomial_coeff_evaluate
!! NAME
!!  polynomial_coeff_evaluate
!!
!! FUNCTION
!!  Compute the energy related to the coefficients from
!!  fitted polynome
!!
!! INPUTS
!!  coefficients(ncoeff)<type(polynomial_coeff_type)> = list of coefficients
!!  disp(3,natom_sc) = atomics displacement between configuration and the reference
!!  natom_sc = number of atoms in the supercell
!!  natom_uc = number of atoms in the unit cell
!!  ncoeff   = number of coefficients
!!  sc_size(3) = size of the supercell (2 2 2 for example)
!!  strain(6) = strain between configuration and the reference 
!!  cells(ncell) = number of the cells into the supercell (1,2,3,4,5)
!!  ncell   = total number of cell to treat by this cpu
!!  index_cells(3,ncell) = indexes of the cells into  supercell (-1 -1 -1 ,...,1 1 1)
!!  comm=MPI communicator
!!
!! OUTPUT
!!  energy = contribution to the energy
!!  fcart(3,natom) = contribution  to the forces
!!  strten(6) = contribution to the stress tensor
!!
!! PARENTS
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      asrq0_free,effective_potential_effpot2ddb,invars9,mkphbs
!!
!! SOURCE
!!
subroutine polynomial_coeff_evaluate(coefficients,disp,energy,fcart,natom_sc,natom_uc,ncoeff,sc_size,&
&                                    strain,strten,cells,ncell,index_cells,comm)

!Arguments ------------------------------------
! scalar

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_evaluate'
!End of the abilint section

  real(dp),intent(out):: energy
  integer, intent(in) :: ncell,ncoeff,natom_sc,natom_uc
  integer, intent(in) :: comm
! array
  real(dp),intent(out):: strten(6)
  real(dp),intent(in) :: strain(6)
  real(dp),intent(out):: fcart(3,natom_sc)
  real(dp),intent(in) :: disp(3,natom_sc)
  integer,intent(in) ::   cells(ncell),index_cells(ncell,3)
  integer,intent(in) :: sc_size(3)
  type(polynomial_coeff_type),intent(in) :: coefficients(ncoeff)
 !Local variables-------------------------------
! scalar
  integer :: i1,i2,i3,ia1,ib1,ia2,ib2,idir1,idir2,ierr,ii
  integer :: icoeff,iterm,idisp1,idisp2,icell,power,weight
  real(dp):: coeff,disp1,disp2,tmp1,tmp2,tmp3
! array
  integer :: cell_atoma1(3),cell_atoma2(3)
  integer :: cell_atomb1(3),cell_atomb2(3)
  character(len=500) :: msg
! *************************************************************************

! Check
  if (any(sc_size <= 0)) then
    write(msg,'(a,a)')' No supercell found for getEnergy'
    MSG_ERROR(msg)
  end if

! Initialisation of variables
  energy     = zero
  fcart(:,:) = zero
  strten(:)  = zero

  do icell = 1,ncell
    ii = (cells(icell)-1)*natom_uc
    i1=index_cells(icell,1); i2=index_cells(icell,2); i3=index_cells(icell,3)
!   Loop over coefficients
    do icoeff=1,ncoeff
!     Set the value of the coefficient
      coeff = coefficients(icoeff)%coefficient
!     Loop over terms of this coefficient
      do iterm=1,coefficients(icoeff)%nterm
!       Set the weight of this term
        weight =coefficients(icoeff)%terms(iterm)%weight
        tmp1 = one
!       Loop over displacement and strain
        do idisp1=1,coefficients(icoeff)%terms(iterm)%ndisp

!         Set to one the acculation of forces and strain
          tmp2 = one
          tmp3 = one

!         Set the power of the displacement:
          power = coefficients(icoeff)%terms(iterm)%power(idisp1)

!         Get the direction of the displacement or strain
          idir1 = coefficients(icoeff)%terms(iterm)%direction(idisp1)

!         Strain case idir => -6, -5, -4, -3, -2 or -1
          if (idir1 < zero)then

            if(abs(strain(abs(idir1))) > tol10)then
!             Accumulate energy fo each displacement (\sum ((A_x-O_x)^Y(A_y-O_c)^Z))
              tmp1 = tmp1 * (strain(abs(idir1)))**power           
              if(power > 1) then
!             Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                tmp3 = tmp3 *  power*(strain(abs(idir1)))**(power-1)
              end if
            else
              tmp1 = zero
              if(power > 1) then
                tmp3 = zero
              end if
            end if
          else
!           Displacement case idir = 1, 2  or 3
!           indexes of the cell of the atom a
            cell_atoma1 = coefficients(icoeff)%terms(iterm)%cell(:,1,idisp1)
            if(cell_atoma1(1)/=0.or.cell_atoma1(2)/=0.or.cell_atoma1(3)/=0) then
!             if the cell is not 0 0 0 we apply PBC:
              cell_atoma1(1) =  (i1-1) + cell_atoma1(1)
              cell_atoma1(2) =  (i2-1) + cell_atoma1(2)
              cell_atoma1(3) =  (i3-1) + cell_atoma1(3)
              call getPBCIndexes_supercell(cell_atoma1(1:3),sc_size(1:3))
!             index of the first atom (position in the supercell if the cell is not 0 0 0)
              ia1 = cell_atoma1(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                   cell_atoma1(2)*sc_size(3)*natom_uc+&
&                   cell_atoma1(3)*natom_uc+&
&                   coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
            else
!             index of the first atom (position in the supercell if the cell is 0 0 0)
              ia1 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
            end if

!           indexes of the cell of the atom b  (with PBC) same as ia1
            cell_atomb1 = coefficients(icoeff)%terms(iterm)%cell(:,2,idisp1)
            if(cell_atomb1(1)/=0.or.cell_atomb1(2)/=0.or.cell_atomb1(3)/=0) then
              cell_atomb1(1) =  (i1-1) + cell_atomb1(1)
              cell_atomb1(2) =  (i2-1) + cell_atomb1(2)
              cell_atomb1(3) =  (i3-1) + cell_atomb1(3)
              call getPBCIndexes_supercell(cell_atomb1(1:3),sc_size(1:3))

!            index of the second atom in the (position in the supercell  if the cell is not 0 0 0) 
              ib1 = cell_atomb1(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                   cell_atomb1(2)*sc_size(3)*natom_uc+&
&                   cell_atomb1(3)*natom_uc+&
&                   coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
            else
!             index of the first atom (position in the supercell if the cell is 0 0 0)
              ib1 = ii + coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
            end if

!           Get the displacement for the both atoms
            disp1 = disp(idir1,ia1)
            disp2 = disp(idir1,ib1)

            if(abs(disp1) > tol10 .or. abs(disp2)> tol10)then
!           Accumulate energy fo each displacement (\sum ((A_x-O_x)^Y(A_y-O_c)^Z))
              tmp1 = tmp1 * (disp1-disp2)**power
              if(power > 1) then
!               Accumulate forces for each displacement (\sum (Y(A_x-O_x)^Y-1(A_y-O_c)^Z+...))
                tmp2 = tmp2 * power*(disp1-disp2)**(power-1)
              end if
            else
              tmp1 = zero
              if(power > 1) then
                tmp2 = zero
              end if
            end if
          end if

          do idisp2=1,coefficients(icoeff)%terms(iterm)%ndisp

            if(idisp2 /= idisp1) then
              idir2 = coefficients(icoeff)%terms(iterm)%direction(idisp2)
              if (idir2 < zero)then
!               Strain case
!               Set the power of the strain:
                power = coefficients(icoeff)%terms(iterm)%power(idisp2)
!               Accumulate energy forces
                tmp2 = tmp2 * (strain(abs(idir2)))**power
!               Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                tmp3 = tmp3 * (strain(abs(idir2)))**power
              else
                cell_atoma2=coefficients(icoeff)%terms(iterm)%cell(:,1,idisp2)
                if(cell_atoma2(1)/=0.or.cell_atoma2(2)/=0.or.cell_atoma2(3)/=0) then
                  cell_atoma2(1) =  (i1-1) + cell_atoma2(1)
                  cell_atoma2(2) =  (i2-1) + cell_atoma2(2)
                  cell_atoma2(3) =  (i3-1) + cell_atoma2(3)
                  call getPBCIndexes_supercell(cell_atoma2(1:3),sc_size(1:3))
!                 index of the first atom (position in the supercell and direction)
!                 if the cell of the atom a is not 0 0 0 (may happen)
                  ia2 = cell_atoma2(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                       cell_atoma2(2)*sc_size(3)*natom_uc+&
&                       cell_atoma2(3)*natom_uc+&
&                       coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                else
!                 index of the first atom (position in the supercell and direction)
                  ia2 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                end if

                cell_atomb2= coefficients(icoeff)%terms(iterm)%cell(:,2,idisp2)

                if(cell_atomb2(1)/=0.or.cell_atomb2(2)/=0.or.cell_atomb2(3)/=0) then
!                 indexes of the cell2 (with PBC)
                  cell_atomb2(1) =  (i1-1) + cell_atomb2(1)
                  cell_atomb2(2) =  (i2-1) + cell_atomb2(2)
                  cell_atomb2(3) =  (i3-1) + cell_atomb2(3)
                  call getPBCIndexes_supercell(cell_atomb2(1:3),sc_size(1:3))

!                 index of the second atom in the (position in the supercell) 
                  ib2 = cell_atomb2(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                       cell_atomb2(2)*sc_size(3)*natom_uc+&
&                       cell_atomb2(3)*natom_uc+&
&                       coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
                else
                  ib2 = ii + coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
                end if

                disp1 = disp(idir2,ia2)
                disp2 = disp(idir2,ib2)

!               Set the power of the displacement:
                power = coefficients(icoeff)%terms(iterm)%power(idisp2)

                tmp2 = tmp2 * (disp1-disp2)**power
                tmp3 = tmp3 * (disp1-disp2)**power

              end if
            end if
          end do

          if(idir1<zero)then
!           Accumule stress tensor
            strten(abs(idir1)) = strten(abs(idir1)) + coeff * weight * tmp3
          else
!           Accumule  forces
            fcart(idir1,ia1) =  fcart(idir1,ia1)  + coeff * weight * tmp2
            fcart(idir1,ib1) =  fcart(idir1,ib1)  - coeff * weight * tmp2
          end if
        end do
        
!       accumule energy
        energy = energy +  coeff * weight * tmp1

      end do
    end do
  end do

! MPI_SUM
  call xmpi_sum(energy, comm, ierr)
  call xmpi_sum(fcart , comm, ierr)
  call xmpi_sum(strten , comm, ierr)

end subroutine polynomial_coeff_evaluate
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
