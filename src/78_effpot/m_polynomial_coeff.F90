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
 use m_sort,only : sort_dp
 use m_crystal,only : crystal_t,symbols_crystal
 use m_supercell, only: getPBCIndexes_supercell,distance_supercell
 use m_xmpi

 implicit none

 public :: polynomial_coeff_broadcast
 public :: polynomial_coeff_evaluate
 public :: polynomial_coeff_free
 public :: polynomial_coeff_init
 public :: polynomial_coeff_getList 
 public :: polynomial_coeff_getName
 public :: polynomial_coeff_getNorder
 public :: polynomial_coeff_getOrder1
 public :: polynomial_coeff_MPIrecv
 public :: polynomial_coeff_MPIsend
 public :: polynomial_coeff_setName
 public :: polynomial_coeff_setCoefficient
 public :: polynomial_coeff_writeXML
 private :: computeNorder
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
!!      m_anharmonics_terms,m_effective_potential_file,m_polynomial_coeff
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
!!      m_anharmonics_terms,m_effective_potential_file,m_polynomial_coeff
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
!!      m_polynomial_coeff,m_polynomial_coeff
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
    do ii=1,coefficients%nterm
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
  do ii=1,coefficients%nterm
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
!!      m_effective_potential,m_polynomial_coeff
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
              cell_atoma1(1) =  i1 + cell_atoma1(1)
              cell_atoma1(2) =  i2 + cell_atoma1(2)
              cell_atoma1(3) =  i3 + cell_atoma1(3)
              call getPBCIndexes_supercell(cell_atoma1(1:3),sc_size(1:3))
!             index of the first atom (position in the supercell if the cell is not 0 0 0)
              ia1 = (cell_atoma1(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                   (cell_atoma1(2)-1)*sc_size(3)*natom_uc+&
&                   (cell_atoma1(3)-1)*natom_uc+&
&                   coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
            else
!             index of the first atom (position in the supercell if the cell is 0 0 0)
              ia1 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
            end if

!           indexes of the cell of the atom b  (with PBC) same as ia1
            cell_atomb1 = coefficients(icoeff)%terms(iterm)%cell(:,2,idisp1)
            if(cell_atomb1(1)/=0.or.cell_atomb1(2)/=0.or.cell_atomb1(3)/=0) then
              cell_atomb1(1) =  i1 + cell_atomb1(1)
              cell_atomb1(2) =  i2 + cell_atomb1(2)
              cell_atomb1(3) =  i3 + cell_atomb1(3)
              call getPBCIndexes_supercell(cell_atomb1(1:3),sc_size(1:3))

!            index of the second atom in the (position in the supercell  if the cell is not 0 0 0) 
              ib1 = (cell_atomb1(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                   (cell_atomb1(2)-1)*sc_size(3)*natom_uc+&
&                   (cell_atomb1(3)-1)*natom_uc+&
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
                  cell_atoma2(1) =  i1 + cell_atoma2(1)
                  cell_atoma2(2) =  i2 + cell_atoma2(2)
                  cell_atoma2(3) =  i3 + cell_atoma2(3)
                  call getPBCIndexes_supercell(cell_atoma2(1:3),sc_size(1:3))
!                 index of the first atom (position in the supercell and direction)
!                 if the cell of the atom a is not 0 0 0 (may happen)
                  ia2 = (cell_atoma2(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                       (cell_atoma2(2)-1)*sc_size(3)*natom_uc+&
&                       (cell_atoma2(3)-1)*natom_uc+&
&                       coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                else
!                 index of the first atom (position in the supercell and direction)
                  ia2 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                end if

                cell_atomb2= coefficients(icoeff)%terms(iterm)%cell(:,2,idisp2)

                if(cell_atomb2(1)/=0.or.cell_atomb2(2)/=0.or.cell_atomb2(3)/=0) then
!                 indexes of the cell2 (with PBC)
                  cell_atomb2(1) =  i1 + cell_atomb2(1)
                  cell_atomb2(2) =  i2 + cell_atomb2(2)
                  cell_atomb2(3) =  i3 + cell_atomb2(3)
                  call getPBCIndexes_supercell(cell_atomb2(1:3),sc_size(1:3))

!                 index of the second atom in the (position in the supercell) 
                  ib2 = (cell_atomb2(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                       (cell_atomb2(2)-1)*sc_size(3)*natom_uc+&
&                       (cell_atomb2(3)-1)*natom_uc+&
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

!!****f* m_polynomial_coeff/polynomial_coeff_getList
!!
!! NAME
!! polynomial_coeff_getList
!!
!! FUNCTION
!! Get the list of all  the possible coefficients for the polynome
!!
!! INPUTS
!! cell(3,nrpt) = indexes of the cells into the supercell (-1 -1 -1, 0 0 0 ...)
!! cutoff = cut-off for the inter atomic forces constants
!! dist(natom,natom,nrpt) = distance between atoms atm1 is in the cell 0 0 0
!!                                                 atm2 is in the nrpt cell (see cell(3,nrpt))
!! crystal<type(crystal_t)> = datatype with all the information for the crystal
!! natom = number of atoms in the unit cell
!! nrpt  = number of cell in the supercell
!!
!! OUTPUT
!! list_symcoeff(6,ncoeff_sym,nsym) = array with the list of the coefficients,
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom 
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!! list_symstr(nstr_sym,nsym) = array with the list of the strain  and the symmetrics
!! nstr_sym = number of coefficient for the strain
!! ncoeff_sym = number of coefficient for the IFC
!!
!!
!! PARENTS
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      destroy_supercell,init_supercell,xred2xcart
!!
!! SOURCE

subroutine polynomial_coeff_getList(cell,crystal,cutoff,dist,list_symcoeff,list_symstr,&
&                                       natom,nstr_sym,ncoeff_sym,nrpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_getList'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 integer,intent(out) :: ncoeff_sym,nstr_sym
 real(dp),intent(in):: cutoff
!arrays
 integer,intent(in) :: cell(3,nrpt)
 real(dp),intent(in):: dist(natom,natom,nrpt)
 type(crystal_t), intent(in) :: crystal
 integer,allocatable,intent(out) :: list_symcoeff(:,:,:),list_symstr(:,:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff2,icoeff_tot,icoeff_tmp,idisy1,idisy2,ii
 integer :: ipesy1,ipesy2,isym,irpt,irpt3,irpt_ref,irpt_sym
 integer :: jj,jsym,mu
 integer :: ncoeff,ncoeff2,ncoeff_max,nu
 integer :: nsym,shift_atm1(3)
 integer :: shift_atm2(3)
 real(dp):: tolsym8
 logical :: found
!arrays
 integer :: sym(3,3)
 integer :: transl(3)
 integer,allocatable :: list(:),list_symcoeff_tmp(:,:,:),list_symcoeff_tmp2(:,:,:)
 integer,allocatable :: list_symstr_tmp(:,:),indsym(:,:,:) ,symrec(:,:,:)
 real(dp),allocatable :: blkval(:,:,:,:,:),tnons(:,:)
 real(dp),allocatable :: wkdist(:),xcart(:,:),xred(:,:)
 real(dp) :: difmin(3)
 real(dp) :: rprimd(3,3)
 real(dp) :: tratom(3)
 character(len=500) :: message

! *************************************************************************


!Initialisation of variables
 nsym   = crystal%nsym
 rprimd = crystal%rprimd
 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 xcart(:,:) = crystal%xcart(:,:)
 xred(:,:)  = crystal%xred(:,:)
 ncoeff_max = nrpt*natom*natom*3*3

!Found the ref cell
 irpt_ref = one 
 do irpt=1,nrpt
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
!     exit
   end if
 end do

!Obtain a list of rotated atom labels:
 ABI_ALLOCATE(indsym,(4,nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 ABI_ALLOCATE(tnons,(3,nsym))
 symrec = crystal%symrec
 tnons  = crystal%tnons

 tolsym8=tol8
 call symatm(indsym,natom,nsym,symrec,tnons,&
&            tolsym8,crystal%typat,crystal%xred)
 ABI_ALLOCATE(blkval,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(list,(natom*nrpt))
 ABI_ALLOCATE(list_symcoeff_tmp,(5,ncoeff_max,nsym))
 ABI_ALLOCATE(wkdist,(natom*nrpt))

!1-Fill strain list
 ABI_ALLOCATE(list_symstr_tmp,(6,nsym))
 list_symstr_tmp = one
 do ia=1,6
   if(list_symstr_tmp(ia,1)==0)cycle
!  Transform the voigt notation
   if(ia<=3)then
     mu=ia;nu=ia
   else
     select case(ia)
     case(4)
       mu=2;nu=3
     case(5)
       mu=1;nu=3
     case(6)
       mu=1;nu=2
     end select
   end if
   do isym=1,nsym
!  Get the symmetry matrix 
     sym(:,:) = crystal%symrel(:,:,isym)
     do idisy1=1,3
       do idisy2=1,3
         if((sym(mu,idisy1)/=0.and.sym(nu,idisy2)/=0).or.&
&           (sym(mu,idisy2)/=0.and.sym(nu,idisy1)/=0))then
!          Transform to the voig notation
           if(idisy1==idisy2)then
             list_symstr_tmp(ia,isym) = idisy1
           else
             if(idisy1==1.or.idisy2==1)then
               if(idisy1==2.or.idisy2==2)then
                 list_symstr_tmp(ia,isym) = 6
               end if
               if(idisy1==3.or.idisy2==3)then
                 list_symstr_tmp(ia,isym) = 5
               end if
             else
               list_symstr_tmp(ia,isym) = 4
             end if
           end if
         end if
       end do
     end do
!    Remove the symetric
!TEST_AM
!      if(list_symstr_tmp(ia,isym) > ia) then
!        list_symstr_tmp(list_symstr_tmp(ia,isym),:) = zero
!      end if
!TEST_AM
   end do
  end do

!Count the number of strain and transfert into the final array
  nstr_sym = zero
  do ia=1,6
    if(list_symstr_tmp(ia,1)/=zero) nstr_sym = nstr_sym + 1
  end do

 if(allocated(list_symstr))then
   ABI_DEALLOCATE(list_symstr)
 end if
 ABI_ALLOCATE(list_symstr,(nstr_sym,nsym))

 icoeff_tmp = 1
 do ia=1,6
   if(list_symstr_tmp(ia,1)/=zero) then
     list_symstr(icoeff_tmp,:) = list_symstr_tmp(ia,:)
     icoeff_tmp = icoeff_tmp + 1
   end if
 end do
!END STRAIN

!Set to one blkval, all the coeff have to be compute
 blkval = one 
 icoeff = one
 icoeff_tot = one
 list_symcoeff_tmp = zero

!2-Fill atom list
!Big loop over generic atom 
 do ia=1,natom
   wkdist(:)=reshape(dist(ia,:,:),(/natom*nrpt/))
   do ii=1,natom*nrpt
     list(ii)=ii
   end do
   call sort_dp(natom*nrpt,wkdist,list,tol8)
   do ii=1,natom*nrpt
!    Get the irpt and ib
     irpt=(list(ii)-1)/natom+1     
     ib=list(ii)-natom*(irpt-1)
     if(dist(ia,ib,irpt) > cutoff ) then
!      If this distance is superior to the cutoff, we don't compute
       blkval(:,ia,:,ib,irpt)= zero
       if(irpt==irpt_ref)blkval(:,ib,:,ia,irpt)= zero 
!      Stop the loop
       exit
     end if
     do mu=1,3
       do nu=1,3
!      Check if : - The coefficient is not yet compute
!                 - The directions are the same
!                 - The atoms are not equivalent
         if (mu/=nu) then 
           blkval(mu,ia,nu,ib,irpt)=zero
           blkval(nu,ia,mu,ib,irpt)=zero
           cycle
         end if
!        Pass if the atoms are identical and in the ref cell
         if(irpt==irpt_ref.and.ia==ib) then
           blkval(mu,ia,nu,ib,irpt)=zero
           blkval(nu,ib,mu,ia,irpt)=zero
           cycle
         end if
         
         if(blkval(mu,ia,nu,ib,irpt)==1)then
!          Loop over symmetries 
           do isym=1,nsym
!            Get the symmetry matrix 
             sym(:,:) = crystal%symrel(:,:,isym)
!            Get the corresponding atom and shift with the symetries 
!            For atom 1
             ipesy1 = indsym(4,isym,ia)
             shift_atm1 = indsym(1:3,isym,ia)
!            And atom 2
             do jj=1,3 ! Apply transformation to original coordinates.
               tratom(jj) = dble(sym(1,jj))*(xred(1,ib)+cell(1,irpt)-tnons(1,isym))&
&                          +dble(sym(2,jj))*(xred(2,ib)+cell(2,irpt)-tnons(2,isym))&
&                          +dble(sym(3,jj))*(xred(3,ib)+cell(3,irpt)-tnons(3,isym))
             end do

!            Find symmetrically equivalent atom
             call symchk(difmin,ipesy2,natom,tratom,transl,crystal%typat(ib),&
&                        crystal%typat,xred(:,:))

!            Put information into array indsym: translations and label
             shift_atm2(:)= transl(:) - shift_atm1(:)

             found = .false.
             do irpt3=1,nrpt
               if(cell(1,irpt3)==shift_atm2(1).and.&
&                 cell(2,irpt3)==shift_atm2(2).and.&
&                 cell(3,irpt3)==shift_atm2(3))then
                 found = .true.
                 irpt_sym = irpt3
               end if
             end do
             
!            Now that a symmetric perturbation has been obtained,
!            including the expression of the symmetry matrix, see
!            if the symmetric perturbations are available
             do idisy1=1,3
               do idisy2=1,3
                 if (idisy1/=idisy2) then
!                  Remove this term (is not computed)
!                  Also remove opposite term... (Srx-Tix) = (Ti-Srx)
                   blkval(idisy1,ipesy1,idisy2,ipesy2,irpt_sym) = 0
                   blkval(idisy2,ipesy1,idisy1,ipesy2,irpt_sym) = 0
                   cycle
                 else
                   if(sym(mu,idisy1)/=0.and.sym(nu,idisy2)/=0)then
                     if(.not.found.or.(irpt_sym==irpt_ref.and.ipesy1==ipesy2)) then
!                      Remove this term (is not computed) Sr-Sr or not include in the cell
!                      Also remove oposite term... (Srx-Tix) = (Ti-Srx)
                       blkval(idisy1,ipesy1,idisy2,ipesy2,irpt_sym) = 0
                       blkval(idisy2,ipesy2,idisy1,ipesy1,irpt_sym) = 0
                       cycle
                     else
!                      Fill the list with the coeff and symetric (need all symetrics)
                       list_symcoeff_tmp(1:4,icoeff,isym)=(/idisy1,ipesy1,ipesy2,irpt_sym/)
                       list_symcoeff_tmp(5,icoeff,isym)= sym(mu,idisy1)
                     end if
                   end if
                 end if
               end do
             end do
           end do ! end loop sym
           icoeff = icoeff + 1 
         end if
!        This coeff is now computed 
         blkval(mu,ia,nu,ib,irpt)= zero
       end do ! end loop nu
     end do ! end loop mu
   end do ! end loop ii
 end do ! end loop ia

!Reset the output
 if(allocated(list_symcoeff))then
   ABI_DEALLOCATE(list_symcoeff)
 end if

!Transfert the final array with all the coefficients
!With this array, we can access to all the terms presents
!ncoeff1 + symetrics
!first dimension is 4 (mu,ia,ib,irpt)
!irpt is the index of the cell of the atom ib in the cell array
!example cell(:,irpt=12) can be (-1 0 -2). The cell of ia is 
!always 0 0 0
!Transfert the final array for the list of irreductible coeff and symetries
!With this array, we can access to the irretuctible coefficients (ncoeff1) and  
!all the symetrics of these coefficients (nsym)
!first dimension is 5 (mu,ia,ib,irpt,icoeff)
!icoeff is the position of this coefficients in the list_fullcoeff array

!1/ step remove the zero coeff in this array
 ncoeff = zero
 do icoeff = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff,1)==zero)))then
     ncoeff = ncoeff + 1
   end if
 end do

 ABI_ALLOCATE(list_symcoeff_tmp2,(6,ncoeff,nsym))
 list_symcoeff_tmp2 = zero
 icoeff = zero
 do icoeff_tmp = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff_tmp,1)==zero)))then
     icoeff = icoeff + 1
     list_symcoeff_tmp2(1:5,icoeff,:) = list_symcoeff_tmp(1:5,icoeff_tmp,:)
   end if
 end do


!2/ set the dimension six of list_symcoeff_tmp2(6,icoeffs,1)
!   and check is a symetric coeff is not coresspondig to an other
!   one, in this case we set this coeff to 0
! ncoeff2 = zero
 do icoeff = 1,ncoeff
!  found the index of each coeff in list_fullcoeff
   do isym = 1,nsym
     icoeff2 = getCoeffFromList(list_symcoeff_tmp2(:,:,1),&
&                               list_symcoeff_tmp2(2,icoeff,isym),&
&                               list_symcoeff_tmp2(3,icoeff,isym),&
&                               list_symcoeff_tmp2(4,icoeff,isym),&
&                               list_symcoeff_tmp2(1,icoeff,isym),&
&                               real(list_symcoeff_tmp2(5,icoeff,isym),dp),ncoeff)
     list_symcoeff_tmp2(6,icoeff,isym) = icoeff2
   end do
 end do

!2.5/do checks
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     if(list_symcoeff_tmp2(6,icoeff,isym)==0)then
       write(message, '(a,i0,a,I0,4a)' )&
&           'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&           'have no equivalent',ch10,&
&           'Action: Contact abinit group'
       MSG_BUG(message)
     else
       if(icoeff /= list_symcoeff_tmp2(6,icoeff,isym))then
         if(list_symcoeff_tmp2(1,icoeff,isym)/=&
&           list_symcoeff_tmp2(1,list_symcoeff_tmp2(6,icoeff,isym),1))then
           write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff_tmp2(6,icoeff,1),ch10,&
&          'because the direction is different:',ch10,&
&          'Action: Contact abinit group'
           MSG_BUG(message)
         end if
         if(list_symcoeff_tmp2(4,icoeff,isym)/=&
&           list_symcoeff_tmp2(4,list_symcoeff_tmp2(6,icoeff,isym),1))then
           write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff_tmp2(6,icoeff,1),ch10,&
&          'because the cell is different',ch10,&
&          'Action: Contact abinit group'
           MSG_BUG(message)
         end if
         if((list_symcoeff_tmp2(2,icoeff,isym)/=&
&            list_symcoeff_tmp2(2,list_symcoeff_tmp2(6,icoeff,isym),1).and.&
&            list_symcoeff_tmp2(3,icoeff,isym)/=&
&            list_symcoeff_tmp2(3,list_symcoeff_tmp2(6,icoeff,isym),1)).and.&
&           (list_symcoeff_tmp2(2,icoeff,isym)/=&
&            list_symcoeff_tmp2(3,list_symcoeff_tmp2(6,icoeff,isym),1).and.&
&            list_symcoeff_tmp2(3,icoeff,isym)/=&
&            list_symcoeff_tmp2(2,list_symcoeff_tmp2(6,icoeff,isym),1)))then
           write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff_tmp2(6,icoeff,1),ch10,&
&          'because the atoms different',ch10,&
&          'Action: Contact abinit group'
           MSG_BUG(message)
         end if
       end if
     end if
   end do
 end do

!3/ Remove useless terms like opposites
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     icoeff2 = list_symcoeff_tmp2(6,icoeff,isym)
     if (icoeff2> icoeff)then
       list_symcoeff_tmp2(:,icoeff2,1) = zero
     end if
     do jsym=1,nsym
       icoeff2 = getCoeffFromList(list_symcoeff_tmp2(:,:,jsym),&
&                                 list_symcoeff_tmp2(3,icoeff,isym),&
&                                 list_symcoeff_tmp2(2,icoeff,isym),&
&                                 list_symcoeff_tmp2(4,icoeff,isym),&
&                                 list_symcoeff_tmp2(1,icoeff,isym),&
&                                 real(list_symcoeff_tmp2(5,icoeff,isym),dp),ncoeff)
       if (icoeff2> icoeff)then
         list_symcoeff_tmp2(:,icoeff2,1) = zero
       end if
     end do
   end do
 end do

!4/ Recount the number of coeff after step 3
 ncoeff2 = zero
 do icoeff = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp2(:,icoeff,1)==zero)))then
     ncoeff2 = ncoeff2 + 1
   end if
 end do

 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_ALLOCATE(list_symcoeff_tmp,(6,ncoeff2,nsym))

 list_symcoeff_tmp = zero
 icoeff = zero

 do icoeff_tmp = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp2(:,icoeff_tmp,1)==zero)))then
     icoeff = icoeff + 1
     list_symcoeff_tmp(:,icoeff,:) = list_symcoeff_tmp2(:,icoeff_tmp,:)
   end if
 end do

!5/ Final transfert
 ABI_ALLOCATE(list_symcoeff,(6,ncoeff2,nsym))
 list_symcoeff = zero
 icoeff = zero
 do icoeff = 1,ncoeff2
   list_symcoeff(1:6,icoeff,:) = list_symcoeff_tmp(1:6,icoeff,:)
   do isym=1,nsym
   end do
 end do

!6/ reset the dimension six of list_symcoeff_tmp2(6,icoeffs,1)
!   and check is a symetric coeff is not coresspondig to an other
!   one, in this case we set this coeff to 0
 do icoeff = 1,ncoeff2
!  found the index of each coeff in list_fullcoeff
   do isym = 1,nsym
     icoeff2 = getCoeffFromList(list_symcoeff(:,:,1),&
&                               list_symcoeff(2,icoeff,isym),&
&                               list_symcoeff(3,icoeff,isym),&
&                               list_symcoeff(4,icoeff,isym),&
&                               list_symcoeff(1,icoeff,isym),&
&                               real(list_symcoeff(5,icoeff,isym),dp),ncoeff)
     list_symcoeff(6,icoeff,isym) = icoeff2
   end do
 end do

!Set the max number of coeff inside list_symcoeff
 ncoeff_sym = ncoeff2

!Deallocation
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(list)
 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_DEALLOCATE(list_symcoeff_tmp2)
 ABI_DEALLOCATE(list_symstr_tmp)
 ABI_DEALLOCATE(indsym) 
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(tnons)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred )
 ABI_DEALLOCATE(wkdist)

end subroutine polynomial_coeff_getList
!!***


!!****f* m_polynomial_coeff/polynomial_coeff_getNorder
!!
!! NAME
!! polynomial_coeff_getNorder
!!
!! FUNCTION
!! Compute and store into the datatype coefficients, all the possible
!! coefficients for given orders
!!
!! INPUTS
!! cutoff = cut-off for the inter atomic forces constants
!! crystal<type(crystal_t)> = datatype with all the information for the crystal
!! powers(2) = array with the minimal and maximal power to be computed
!! option = 0 compute all terms
!!          1 still in development
!! comm = MPI communicator
!! anharmstr = logical, optional : TRUE, the anharmonic strain are computed (\eta)^power ...
!!                                   FALSE, (default) the anharmonic strain are not computed
!!
!! OUTPUT
!! polynomial_coeff<(type(polynomial_coeff_type)>(ncoeff) = array of datatype with the polynomial_coeff
!! ncoeff = number of coefficients
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      destroy_supercell,init_supercell,xred2xcart
!!
!! SOURCE

subroutine polynomial_coeff_getNorder(coefficients,crystal,cutoff,ncoeff,powers,option,comm,&
&                                     anharmstr,spcoupling)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_getNorder'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,comm
 integer,intent(out):: ncoeff
 real(dp),intent(in):: cutoff
 logical,optional,intent(in) :: anharmstr,spcoupling
!arrays
 integer,intent(in) :: powers(2)
 type(crystal_t), intent(inout) :: crystal
 type(polynomial_coeff_type),allocatable,intent(inout) :: coefficients(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff2,ii,irpt,irpt_ref
 integer :: lim1,lim2,lim3
 integer :: natom,ncoeff_max,ncoeff_sym,ncoeff_tot,nrpt,nsym,nstr_sym
 integer :: r1,r2,r3
 logical :: need_anharmstr,need_spcoupling
 
!arrays
 integer :: ncell(3)
 integer,allocatable :: cell(:,:),compatibleCoeffs(:,:)
 integer,allocatable :: list_symcoeff(:,:,:),list_symstr(:,:),list_coeff(:)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: dist(:,:,:),rpt(:,:)
 real(dp),allocatable :: xcart(:,:),xred(:,:)
 character(len=5),allocatable :: symbols(:)
 character(len=500) :: message
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs_tmp
 
! *************************************************************************

!Free the output
 if(allocated(coefficients))then
   do ii =1,size(coefficients)
     call polynomial_coeff_free(coefficients(ii))
   end do
   ABI_DEALLOCATE(coefficients)
 end if

!Check
 if(option > crystal%ntypat)then
   write(message, '(3a)' )&
&       'Option can not be superior to ntypat ',ch10,&
&       'Action: contact abinit group'
   MSG_ERROR(message)
 end if

!Initialisation of variables
 need_anharmstr = .TRUE.
 if(present(anharmstr)) need_anharmstr = anharmstr
 need_spcoupling = .TRUE.
 if(present(spcoupling)) need_spcoupling = spcoupling

 natom  = crystal%natom
 nsym   = crystal%nsym
 rprimd = crystal%rprimd

 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 xcart(:,:) = crystal%xcart(:,:)
 xred(:,:)  = crystal%xred(:,:)

!Set the size of the interaction
 ncell = (/anint(cutoff/rprimd(1,1))+1,&
&          anint(cutoff/rprimd(2,2))+1,&
&          anint(cutoff/rprimd(3,3))+1/)

 lim1=((ncell(1)/2))
 lim2=((ncell(2)/2))
 lim3=((ncell(3)/2))
 if(mod(ncell(1),2)/=0) lim1=lim1+1
 if(mod(ncell(2),2)/=0) lim2=lim2+1
 if(mod(ncell(3),2)/=0) lim3=lim3+1
 nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)

!compute new ncell
 ncell(1) = 2*lim1+1
 ncell(2) = 2*lim2+1
 ncell(3) = 2*lim3+1

 !Build the rpt point
 ABI_ALLOCATE(rpt,(3,nrpt))
 ABI_ALLOCATE(cell,(3,nrpt))
 
!WARNING:
!Put the reference cell into the first element
!the code will first deal with the atoms of the first cell
 irpt = one
 irpt_ref = one 
 rpt(:,1) = zero
 cell(:,irpt)=zero
!Fill other rpt:
 do r1=lim1,-lim1,-1
   do r2=lim2,-lim2,-1
     do r3=lim3,-lim3,-1
       if(r1==0.and.r2==0.and.r3==0) then
         cycle
       end if
       irpt=irpt+1
       rpt(1,irpt)=r1*rprimd(1,1)+r2*rprimd(1,2)+r3*rprimd(1,3)
       rpt(2,irpt)=r1*rprimd(2,1)+r2*rprimd(2,2)+r3*rprimd(2,3)
       rpt(3,irpt)=r1*rprimd(3,1)+r2*rprimd(3,2)+r3*rprimd(3,3)
       cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
     end do
   end do
 end do

 ABI_ALLOCATE(symbols,(natom))
 call symbols_crystal(crystal%natom,crystal%ntypat,crystal%npsp,&
&                     symbols,crystal%typat,crystal%znucl)

!Compute the distances between atoms
!Now dist(ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.
 ABI_ALLOCATE(dist,(natom,natom,nrpt))
 dist = zero
 do ia=1,natom
   do ib=1,natom
     do irpt=1,nrpt
       dist(ia,ib,irpt) = ((xcart(1,ib)-xcart(1,ia)+rpt(1,irpt))**2+&
&                          (xcart(2,ib)-xcart(2,ia)+rpt(2,irpt))**2+&
&                          (xcart(3,ib)-xcart(3,ia)+rpt(3,irpt))**2)**0.5
     end do
   end do
 end do

 call polynomial_coeff_getList(cell,crystal,cutoff,dist,list_symcoeff,list_symstr,&
&                                  natom,nstr_sym,ncoeff_sym,nrpt)
 ncoeff_tot = ncoeff_sym+nstr_sym

!Check the distanceance bewteen coefficients and store integer:
! 0: the mix between these coefficient is not possible 
! 1: the mix between these coefficient is possible 
 ABI_ALLOCATE(compatibleCoeffs,(ncoeff_tot,ncoeff_tot))
 compatibleCoeffs(:,:) = one

 do icoeff=1,ncoeff_tot
   do icoeff2=1,ncoeff_tot     
!    Select case:
!    if both icoeff are displacement => check the distance
!    if both icoeff are strain => check the flag
!    Otherwise cycle (we keep the term)
     if(icoeff>ncoeff_sym.and.icoeff2<=ncoeff_sym)cycle
     if(icoeff2<=ncoeff_sym.and.icoeff2>ncoeff_sym)cycle
     if((icoeff>ncoeff_sym.or.icoeff2>ncoeff_sym).and.&
&       .not.need_anharmstr.and..not.need_spcoupling) then
       compatibleCoeffs(icoeff,icoeff2) = zero
     end if
     if(icoeff2<=ncoeff_sym.and.icoeff2<=ncoeff_sym)then
       if(distance_supercell(xcart(:,list_symcoeff(2,icoeff,1)),&
&                  xcart(:,list_symcoeff(2,icoeff2,1)),rprimd,&
&                  cell(:,1),cell(:,1))>=cutoff.or.&
&         distance_supercell(xcart(:,list_symcoeff(2,icoeff,1)),&
&                  xcart(:,list_symcoeff(3,icoeff2,1)),rprimd,&
&                  cell(:,1),cell(:,list_symcoeff(4,icoeff2,1)))>=cutoff.or.&
&         distance_supercell(xcart(:,list_symcoeff(3,icoeff,1)),&
&                  xcart(:,list_symcoeff(2,icoeff2,1)),rprimd,&
&                  cell(:,list_symcoeff(4,icoeff,1)),cell(:,1))>=cutoff.or.&
&         distance_supercell(xcart(:,list_symcoeff(3,icoeff,1)),&
&                  xcart(:,list_symcoeff(3,icoeff2,1)),rprimd,&
&                  cell(:,list_symcoeff(4,icoeff,1)),&
&                  cell(:,list_symcoeff(4,icoeff2,1)))>=cutoff)then
         compatibleCoeffs(icoeff,icoeff2) = zero
       end if
     end if
   end do
 end do

!first call to this routine in order to count the number of maximum coefficients
 ABI_ALLOCATE(list_coeff,(0))
 ABI_ALLOCATE(coeffs_tmp,(0))
 icoeff  = 1
 icoeff2 = 0

 call computeNorder(cell,coeffs_tmp,compatibleCoeffs,list_symcoeff,list_symstr,list_coeff,&
&                   icoeff,icoeff2,natom,ncoeff_sym,nstr_sym,icoeff,nrpt,nsym,1,powers(1),powers(2),&
&                   symbols,nbody=option,compute=.false.,&
&                   anharmstr=need_anharmstr,spcoupling=need_spcoupling)
 ABI_DEALLOCATE(coeffs_tmp)

!Set to the maximum of possible coefficients
 ncoeff_max =  icoeff2

!Second call to this routine in order to compute the coefficients
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 icoeff  = 1
 icoeff2 = 0
 call computeNorder(cell,coeffs_tmp,compatibleCoeffs,list_symcoeff,list_symstr,list_coeff,&
&               icoeff,icoeff2,natom,ncoeff_sym,nstr_sym,ncoeff_max,nrpt,nsym,1,powers(1),powers(2),&
&               symbols,nbody=option,compute=.true.,anharmstr=need_anharmstr,spcoupling=need_spcoupling)

 ABI_DEALLOCATE(list_coeff)

!Final tranfert
!1- Count the total number of coefficient
 ncoeff = zero
 do icoeff=1,ncoeff_max
   if (coeffs_tmp(icoeff)%coefficient /= zero) then
     ncoeff = ncoeff + 1
   end if
 end do

!2- Transfer
 ABI_ALLOCATE(coefficients,(ncoeff))
 icoeff2 = zero
 do icoeff=1,ncoeff_max
   if (coeffs_tmp(icoeff)%coefficient /= zero) then
     icoeff2 = icoeff2 + 1
     call polynomial_coeff_init(one,coeffs_tmp(icoeff)%nterm,coefficients(icoeff2),&
&                               coeffs_tmp(icoeff)%terms,&
&                               name=coeffs_tmp(icoeff)%name)
   end if
 end do

!Free them all
 do icoeff=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff))
 end do
 if(allocated(coeffs_tmp)) then
   ABI_DEALLOCATE(coeffs_tmp)
 end if

 ABI_DEALLOCATE(cell)
 ABI_DEALLOCATE(dist)
 ABI_DEALLOCATE(compatibleCoeffs)
 ABI_DEALLOCATE(list_symcoeff)
 ABI_DEALLOCATE(list_symstr)
 ABI_DEALLOCATE(rpt)
 ABI_DEALLOCATE(symbols)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred)

end subroutine polynomial_coeff_getNorder
!!***

!!****f* m_polynomial_coeff/computeNorder
!! NAME
!! computeNorder
!!
!! FUNCTION
!! Recursive routine to compute the order N of a all the possible coefficient
!! from the list list_symcoeff and list_symstr.
!!
!! INPUTS
!! cell(3,nrpt) = indexes of the cells into the supercell (-1 -1 -1, 0 0 0 ...)
!! compatibleCoeffs(ncoeff+nstr,ncoeff+nstr) = array with the list of compatible coefficients 0 or 1
!! list_symcoeff(6,ncoeff_sym,nsym) = array with the list of the coefficients,
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom 
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!! list_symstr(nstr_sym,nsym) = array with the list of the strain  and the symmetrics
!! index_coeff_in(power-1) = list of previous coefficients computed (start with 0)
!! icoeff = current indexes of the cofficients (start we 1)
!! icoeff_tot = current number of coefficients already computed (start we 0)
!! natom = number of atoms in the unit cell
!! nstr = number of coefficient for related to the atomic displacment into list_symcoeff
!! nstr = number of coefficient for related to the strain into list_symstr
!! ncoeff_out = number of maximum coefficients
!! nrpt = number of cell 
!! nsym = number of symmetries in the system
!! power = initial power to be computed (can be < power_min, this routine will skip the firts power)
!! power_min = minimal power to be computed
!! power_max = maximum power to be computed
!! symbols(natom) = array with the symbols of each atoms (Sr,O,Ti,...) 
!! nbody = optional, number of body for the coefficients, for example:
!!                   0 => all the terms
!!                   1 => only (Sr_x-T_y)^power and (Sr_x-T_y)^power\eta^power  ...
!! compute = logical, optional: TRUE if we store the coefficients
!!                              FALSE just to count the number of coefficient
!! anharmstr = logical, optional : TRUE, the anharmonic strain are computed
!!                                   FALSE, (default) the anharmonic strain are not computed
!!
!! OUTPUT
!! icoeff = current indexes of the cofficients (start we 1)
!! icoeff_tot = current number of coefficients already computed (start we 0)
!! polynomial_coeff<(type(polynomial_coeff_type)>(ncoeff_out) = array of datatype with 
!!                                                              the polynomial_coeff
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

recursive subroutine computeNorder(cell,coeffs_out,compatibleCoeffs,list_coeff,list_str,&
&                                  index_coeff_in,icoeff,icoeff_tot,natom,ncoeff,nstr,ncoeff_out,&
&                                  nrpt,nsym,power,power_min,power_max,symbols,nbody,&
&                                  compute,anharmstr,spcoupling)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'computeNorder'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalar 
 integer,intent(in) :: natom,ncoeff,power,power_min,power_max,ncoeff_out,nsym,nrpt,nstr
 integer,intent(inout) :: icoeff,icoeff_tot
 logical,optional,intent(in) :: compute,anharmstr,spcoupling
 integer,optional,intent(in) :: nbody
!arrays
 integer,intent(in) :: cell(3,nrpt),compatibleCoeffs(ncoeff+nstr,ncoeff+nstr)
 integer,intent(in) :: list_coeff(6,ncoeff,nsym),list_str(nstr,nsym)
 integer,intent(in) :: index_coeff_in(power-1)
 type(polynomial_coeff_type),intent(inout) :: coeffs_out(ncoeff_out)
 character(len=5),intent(in) :: symbols(natom)
!Local variables ---------------------------------------
!scalar
 integer :: ia,ib,ii,icoeff1,icoeff_tmp,icoeff_str
 integer :: irpt,isym,idisp,iterm,mu,nbody_in,ncoeff_max,ndisp,pa,pb
 integer :: nterm_max
 real(dp):: coefficient,weight
 logical :: need_compute,compatible,possible,need_anharmstr,need_spcoupling
!arrays
 integer,allocatable :: index_coeff(:)
 integer,allocatable :: atindx(:,:)
 integer,allocatable :: cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=100):: name
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)
! *************************************************************************

!Set the inputs
 need_compute = .TRUE.
 need_anharmstr = .TRUE.
 need_spcoupling = .TRUE.
 nbody_in = 0 !all kind of terms
 if(present(compute)) need_compute = compute
 if(present(nbody)) nbody_in = nbody
 if(present(anharmstr)) need_anharmstr = anharmstr
 if(present(spcoupling)) need_spcoupling = spcoupling
 if(power <= power_max)then   
   
!  Initialisation of variables
   nterm_max  = nsym
   ncoeff_max = (ncoeff+nstr)
   ndisp = power
   icoeff_tmp = zero
   ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
   ABI_ALLOCATE(terms,(nterm_max))
   ABI_ALLOCATE(atindx,(2,ndisp))
   ABI_ALLOCATE(cells,(3,2,ndisp))
   ABI_ALLOCATE(dir_int,(ndisp))
   ABI_ALLOCATE(powers,(ndisp))
   ABI_ALLOCATE(index_coeff,(power))

   index_coeff(1:power-1) = index_coeff_in(:)

   do icoeff1=icoeff,ncoeff+nstr
!    If the distance between the 2 coefficients is superior than the cut-off,
!    we cycle
!    If the power is one, we need to set icoeff to icoeff1
     if(power==1) icoeff = icoeff1

     if(compatibleCoeffs(icoeff,icoeff1)==0) cycle

!    Reset the flag compatible and possible
     compatible = .TRUE.
     possible   = .TRUE.

     index_coeff(power) = icoeff1
     iterm = zero
     coefficient = one

     if(power >= power_min) then
!      Loop over symetries
       do isym=1,nsym
!        Treat this coeff
         weight = 1
         do idisp=1,ndisp
!          Get index of this displacement term
           if(index_coeff(idisp)<=ncoeff)then
             mu   = list_coeff(1,index_coeff(idisp),isym)
             ia   = list_coeff(2,index_coeff(idisp),isym)
             ib   = list_coeff(3,index_coeff(idisp),isym)
             irpt = list_coeff(4,index_coeff(idisp),isym)
             weight = weight*list_coeff(5,index_coeff(idisp),isym)
!            Fill First term arrays 
             atindx(1,idisp) = ia; atindx(2,idisp) = ib;
             dir_int(idisp) = mu
             powers(idisp)   = 1
             cells(:,1,idisp) = (/0,0,0/)
             cells(:,2,idisp) = cell(:,irpt)
           else
             icoeff_str = index_coeff(idisp)-ncoeff
             atindx(1,idisp) = 0; atindx(2,idisp) = 0;
             dir_int(idisp) = -1 * list_str(icoeff_str,isym)
             powers(idisp)   = 1
             cells(:,1,idisp) = (/0,0,0/)
             cells(:,2,idisp) = (/0,0,0/)
           end if
         end do
         
         iterm = iterm + 1
         call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,&
&                                  weight,check=.true.)
       end do!end do sym
   
       if(iterm > 0)then
!        Do some checks
!        -------------
!        1-Check if the coefficient is full anharmonic strain and if we need to compute it
         if(all(terms(1)%direction(:) < zero))then
           compatible = (need_anharmstr .or. need_spcoupling)
           possible = need_anharmstr
         end if
!        1-Check if the coefficient is strain-coupling and if we need to compute it         
         if(any(terms(1)%direction(:) < zero).and.any(terms(1)%direction(:) > zero))then
           possible   = need_spcoupling
           compatible = need_spcoupling
         end if
         
!        ------------
!        2-Check if this terms is compatible with nbody
         if(nbody_in > zero)then
           pa = one ; pb = one
           ia = zero ; ib = zero
!          Count the number of terms and the power           
           do ii=1,terms(1)%ndisp
             if(terms(1)%direction(ii) < zero) then
               pb = pb*terms(1)%power(ii)
               ib = ib + 1
             else
               pa = pa*terms(1)%power(ii)
               ia = ia + 1
             end if
           end do
           if(ia <= nbody_in)then
             if(ia==nbody_in.and.mod(pa,2)==zero)then
               if(ib==zero)then
                 compatible = .FALSE.
                 possible   = .TRUE.
               else if (ib==nbody_in.and.mod(pb,2)==zero) then
                 compatible = .FALSE.
                 possible   = .TRUE.               
               else
                possible = .FALSE.
                compatible = .FALSE.                 
               end if
             else
                possible = .FALSE.
                compatible = .FALSE.
             end if
           else
             compatible = .FALSE.
             possible = .FALSE.
           end if
         end if

         if(possible)then
!          increase coefficients and set it
           icoeff_tmp = icoeff_tmp + 1
           icoeff_tot = icoeff_tot + 1
           call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),&
&                                     terms(1:iterm),check=.true.)
         end if
       end if

!      Deallocate the terms
       do iterm=1,nterm_max
         call polynomial_term_free(terms(iterm))
       end do

     end if!end if power < power_min

     if(compatible)then
       call computeNorder(cell,coeffs_out,compatibleCoeffs,list_coeff,list_str,index_coeff,&
&                         icoeff1,icoeff_tot,natom,ncoeff,nstr,ncoeff_out,nrpt,nsym,power+1,&
&                         power_min,power_max,symbols,nbody=nbody_in,compute=need_compute,&
&                         anharmstr=need_anharmstr,spcoupling=need_spcoupling)
     end if
   end do

   ABI_DEALLOCATE(terms)
   ABI_DEALLOCATE(atindx)
   ABI_DEALLOCATE(cells)
   ABI_DEALLOCATE(dir_int)
   ABI_DEALLOCATE(index_coeff)
   ABI_DEALLOCATE(powers)

!  Transfer in the final array
   icoeff1 = zero
   do icoeff_tmp=1,ncoeff_max
     if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
!      Increase icoeff and fill the coeffs_out array
       icoeff_tot = icoeff_tot + 1
       if(need_compute)then
         name = ''
!        Get the name of this coefficient
         call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols,recompute=.TRUE.)
         call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                                     coeffs_out(icoeff_tot),coeffs_tmp(icoeff_tmp)%terms,&
&                                     name=name)
       end if
     end if
   end do
!  Deallocation
   do icoeff1=1,ncoeff_max
     call polynomial_coeff_free(coeffs_tmp(icoeff1))
   end do
   ABI_DEALLOCATE(coeffs_tmp)
 end if

end subroutine computeNorder
!!***

!!****f* m_polynomial_coeff/getCoeffFromList
!!
!! NAME
!! getCoeffFromList
!!
!! FUNCTION
!! get the index of a coefficient into the list_coeff
!!
!! INPUTS
!! list_symcoeff(6,ncoeff_sym,nsym) = array with the list of the coefficients,
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom 
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!! ia = index of the atom 1
!! ib = index of the atom 1
!! irpt = indexes of the cell of the second atom 
!! mu = direction of the IFC 
!! weight =  weight of the term (-1 or 1) 
!! ncoeff = number of total coefficients in the list
!!
!! OUTPUT
!! coeff = index of the coefficient
!!
!! PARENTS
!!      m_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

function getCoeffFromList(list_coeff,ia,ib,irpt,mu,weight,ncoeff) result(coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getCoeffFromList'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalar
 integer,intent(in) :: ia,ib,irpt,mu,ncoeff
 real(dp),intent(in):: weight
 integer :: coeff
!arrays
 integer,intent(in) :: list_coeff(6,ncoeff)
!Local variables-------------------------------
!scalar
 integer :: icoeff
!arrays

! *************************************************************************
 coeff = zero
 do icoeff = 1,ncoeff
   if(mu==list_coeff(1,icoeff).and.&
&     ia==list_coeff(2,icoeff).and.&
&     ib==list_coeff(3,icoeff).and.&
&     irpt==list_coeff(4,icoeff)) then
     coeff = icoeff
     exit
   end if
 end do

end function getCoeffFromList
!!***

!!****f* m_polynomial_coeff/polynomial_coeff_getOrder1
!!
!! NAME
!! polynomial_coeff_getOrder1
!!
!! FUNCTION
!! Compute the first order polynomial coefficients from the list
!!
!! INPUTS
!! cell(3,nrpt) = indexes of the cells into the supercell (-1 -1 -1, 0 0 0 ...)
!! cutoff_in = cut-off for the inter atomic forces constants
!! list_symcoeff(6,ncoeff_sym,nsym) = array with the list of the coefficients,
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom 
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!! list_symstr(nstr_sym,nsym) = array with the list of the strain  and the symmetrics
!! natom = number of atoms in the unit cell
!! nrpt = number of cell
!! nsym = number of symmetries in the system
!! rprimd(3,3) = primitive lattice vectors
!! symbols(natom) = array with the symbols of each atoms (Sr,O,Ti,...) 
!! xcart(3,natom) = cartesian coordinates of the atoms in the unit cell
!! comm = MPI communicator
!!
!! OUTPUT
!! polynomial_coeff<(type(polynomial_coeff_type)>(ncoeff_out) = array of datatype with 
!!                                                              the polynomial_coeff
!! ncoeff_out = number of coefficients
!!
!! PARENTS
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      destroy_supercell,init_supercell,xred2xcart
!!
!! SOURCE

subroutine polynomial_coeff_getOrder1(cell,coeffs_out,cutoff_in,list_symcoeff,list_symstr,&
&                                         natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                         rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_coeff_getOrder1'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff,nsym,nrpt
 integer,intent(out) :: ncoeff_out
 real(dp),intent(in) :: cutoff_in
!arrays
 integer,intent(in) :: cell(3,nrpt)
 integer,intent(in) :: list_symcoeff(6,ncoeff,nsym),list_symstr(6,nsym)
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 character(len=5),intent(in) :: symbols(natom)
 type(polynomial_coeff_type),allocatable,intent(inout) :: coeffs_out(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff_tmp,irpt,irpt_ref
 integer :: isym,iterm,mu,ncoeff_max,ndisp,nterm_max
 real(dp):: coefficient,weight
!arrays
 integer,allocatable :: atindx(:,:),cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=1) :: dir_char(3)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=100):: name
 character(len=500) :: message
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)
!TEST_AM
! character(len=fnlen) :: filename
!TEST_AM
! *************************************************************************

!Initialisation of variables
 nterm_max  = nsym
 ncoeff_max = ncoeff
 ndisp = 1
 
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))


 icoeff_tmp = zero 
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))

!Found the ref cell
 irpt_ref = one 
 do irpt=1,nrpt
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
     exit
   end if
 end do

 write(message,'(3a)') " Irreductible coefficient and associated atom 1, atom 2 and direction:",ch10,&
&                     " for the 1st order"
 call wrtout(std_out,message,'COLL') 

 do icoeff=1,ncoeff
!  Reset counter
   iterm = zero
   coefficient = one
   do isym=1,nsym
!    Get index of this displacement term
     mu   = list_symcoeff(1,icoeff,isym)
     ia   = list_symcoeff(2,icoeff,isym)
     ib   = list_symcoeff(3,icoeff,isym)
     irpt = list_symcoeff(4,icoeff,isym)
     weight = list_symcoeff(5,icoeff,isym)
!    And fill arrays for the initialisation
     atindx(1,1) = ia; atindx(2,1) = ib; dir_char(1) = mutodir(mu);
     dir_int(1)  = mu
     ndisp  = 1
     powers(:)   = one
     cells(:,1,1) = (/0,0,0/)
     cells(:,2,1) = cell(:,irpt)
     iterm = iterm + 1
     call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,weight,check=.true.)
   end do!end do sym

   if(iterm > 0)then
!  increase coefficients and set it
     icoeff_tmp = icoeff_tmp + 1
     call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm),check=.true.)
   end if

!  Deallocate the terms
   do iterm=1,nterm_max
     call polynomial_term_free(terms(iterm))
   end do
 end do!end do coeff_sym

 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cells)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(powers)

!Count the number of terms
 ncoeff_out = zero
 do icoeff_tmp=1,ncoeff_max
   if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
     ncoeff_out = ncoeff_out + 1
   end if
 end do

!Transfer in the final array
 ABI_ALLOCATE(coeffs_out,(ncoeff_out))
 icoeff = zero
 do icoeff_tmp=1,ncoeff_max
   if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols,recompute=.TRUE.)
!    Increase icoeff and fill the coeffs_out array
     icoeff = icoeff + 1
      call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
 &                               coeffs_out(icoeff),coeffs_tmp(icoeff_tmp)%terms,name=name)

     write(message,'(2a)')' ',trim(name)
     call wrtout(std_out,message,'COLL') 
     
     do iterm = 1,coeffs_tmp(icoeff_tmp)%nterm
       write(message,'(a,I0,a,I0,2a)') '    Atom ',coeffs_tmp(icoeff_tmp)%terms(iterm)%atindx(1,1),&
&                       ' and atom ',coeffs_tmp(icoeff_tmp)%terms(iterm)%atindx(2,1),&
&                       ' in the direction ',mutodir(coeffs_tmp(icoeff_tmp)%terms(iterm)%direction(1))
       if(any(coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(:,2,1)/=zero))then
         write(message,'(2a,I0,a,I0,a,I0,a)') trim(message),' in the cell ',&
&                                       coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(1,2,1),' ',&
&                                       coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(2,2,1),' ',&
&                                       coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(3,2,1),'.'
       end if
       call wrtout(std_out,message,'COLL') 
     end do
   end if
 end do

!TEST_AM
! filename = "terms_1st_order.xml"
! call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
!TEST_AM

 write(message,'(a,1x,I0,a)') ch10,&
&       ncoeff_out,' fitted coefficients for the 1st order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine polynomial_coeff_getOrder1
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
