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

module m_polynomial_coeff

 use defs_basis
 use m_errors
 use m_abicore
 use m_polynomial_term
 use m_xmpi
#ifdef HAVE_MPI2
 use mpi
#endif

 use m_sort,      only : sort_dp
 use m_io_tools,  only : open_file, get_unit
 use m_symtk,     only : symchk, symatm
 use m_crystal,   only : crystal_t,symbols_crystal
 use m_supercell, only : getPBCIndexes_supercell,distance_supercell,findBound_supercell
 use m_geometry,  only : xcart2xred,metric
 use m_dtfil,     only : isfile

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
 public :: polynomial_coeff_getEvenAnhaStrain
 private :: computeNorder
 private :: computeCombinationFromList
 private :: computeSymmetricCombinations 
 private :: getCoeffFromList
 private :: generateTermsFromList
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

   character(len=200) :: name = ""
!     Name of the polynomial_coeff (Sr_y-O1_y)^3) for example

   integer :: nterm = 0
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

 interface operator (+)
   module procedure coeffs_list_conc
 end interface operator (+)

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
!!      m_polynomial_coeff,mover_effpot
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_init(coefficient,nterm,polynomial_coeff,terms,name,check)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nterm
 real(dp),intent(in) :: coefficient
 logical,optional,intent(in) :: check
!arrays
 character(len=200),optional,intent(in) :: name
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
 character(len=200) :: name_tmp
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
     if(abs(weights(iterm1)) < tol16)cycle
     weights(iterm1) = terms(iterm1)%weight
     do iterm2=iterm1+1,nterm
       if(abs(weights(iterm2)) < tol16)cycle
!      if the terms are identical we check the weight
       if(terms(iterm1)==terms(iterm2))then
         weights(iterm1) = weights(iterm1) + terms(iterm2)%weight
         weights(iterm2) = 0
       end if
     end do
     if(abs(weights(iterm1)) > tol16) then
       weights(iterm1)= anint(weights(iterm1)/weights(iterm1))
     end if
   end do

!  Count the number of terms
   nterm_tmp = 0
   do iterm1=1,nterm
     if(abs(weights(iterm1)) > tol16)then
       nterm_tmp = nterm_tmp + 1
     end if
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
   if(abs(weights(ii)) > tol16)then
     iterm1 = iterm1 + 1
     call polynomial_term_init(terms(ii)%atindx,terms(ii)%cell,terms(ii)%direction,terms(ii)%ndisp,&
&                              terms(ii)%nstrain,polynomial_coeff%terms(iterm1),terms(ii)%power_disp,&
&                              terms(ii)%power_strain,terms(ii)%strain,terms(ii)%weight)
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
!!      m_polynomial_coeff,mover_effpot
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_free(polynomial_coeff)

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
 polynomial_coeff%nterm = 0
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
!!      m_effective_potential_file,mover_effpot
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_setCoefficient(coefficient,polynomial_coeff)

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
!!      m_effective_potential_file,m_fit_polynomial_coeff,m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_setName(name,polynomial_coeff)

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
!!      m_effective_potential_file,m_fit_polynomial_coeff,m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_getName(name,polynomial_coeff,symbols,recompute,iterm)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: iterm
!arrays
 character(len=5),intent(in) :: symbols(:)
 character(len=200),intent(out):: name
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
 character(len=2) :: power_disp,power_dispchar
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
           if(polynomial_coeff%terms(ii)%direction(idisp) > 0) then
             iterm_in = ii
             if(any(polynomial_coeff%terms(ii)%cell(:,1,idisp) /= 0).or.&
&               any(polynomial_coeff%terms(ii)%cell(:,2,idisp) /= 0)) then
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
 if(iterm_in > polynomial_coeff%nterm.or.iterm_in < 0) then
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
     write(power_dispchar,'(I0)') polynomial_coeff%terms(iterm_in)%power_disp(idisp)
     power_disp=trim(power_dispchar)

     atm1=symbols(polynomial_coeff%terms(iterm_in)%atindx(1,idisp))
     atm2=symbols(polynomial_coeff%terms(iterm_in)%atindx(2,idisp))
     dir=mutodir(polynomial_coeff%terms(iterm_in)%direction(idisp))
     cell_atm1=polynomial_coeff%terms(iterm_in)%cell(:,1,idisp)
     cell_atm2=polynomial_coeff%terms(iterm_in)%cell(:,2,idisp)
!    Construct ATM1
     if (any(cell_atm1(:) /= 0) )then
       write(atm1_tmp,'(4a,I0,a,I0,a,I0,a)')  trim(atm1),"_",dir,"[",cell_atm1(1)," ",&
&                                               cell_atm1(2)," ",cell_atm1(3),"]"
     else
       atm1_tmp = trim(atm1)//"_"//dir
     end if
!      Construct ATM2
     if(any(cell_atm2(:) /= 0))then
       write(atm2_tmp,'(4a,I0,a,I0,a,I0,a)')  trim(atm2),"_",dir,"[",cell_atm2(1)," ",&
 &                                              cell_atm2(2)," ",cell_atm2(3),"]"
     else
       atm2_tmp = trim(atm2)//"_"//dir
     end if

     text="("//trim(atm1_tmp)//"-"//trim(atm2_tmp)//")^"//power_disp
     name = trim(name)//trim(text)
   end do
   !Strain case
   do idisp=1,polynomial_coeff%terms(iterm_in)%nstrain
     write(power_dispchar,'(I0)') polynomial_coeff%terms(iterm_in)%power_strain(idisp)
     power_disp=trim(power_dispchar)
     dir=mutodir(3+polynomial_coeff%terms(iterm_in)%strain(idisp))
     text="("//"eta_"//trim(dir)//")^"//power_disp
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
!!      m_effective_potential_file,m_fit_polynomial_coeff,m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_broadcast(coefficients, source, comm)

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
    call xmpi_bcast(coefficients%terms(ii)%nstrain, source, comm, ierr)
  end do

! Allocate arrays on the other nodes
  if (xmpi_comm_rank(comm) /= source) then
    do ii = 1,coefficients%nterm
      ABI_ALLOCATE(coefficients%terms(ii)%atindx,(2,coefficients%terms(ii)%ndisp))
      coefficients%terms(ii)%atindx = 0
      ABI_ALLOCATE(coefficients%terms(ii)%direction,(coefficients%terms(ii)%ndisp))
      ABI_ALLOCATE(coefficients%terms(ii)%cell,(3,2,coefficients%terms(ii)%ndisp))
      ABI_ALLOCATE(coefficients%terms(ii)%power_disp,(coefficients%terms(ii)%ndisp))
      ABI_ALLOCATE(coefficients%terms(ii)%power_strain,(coefficients%terms(ii)%nstrain))
      ABI_ALLOCATE(coefficients%terms(ii)%strain,(coefficients%terms(ii)%nstrain))
    end do
  end if

! Transfert value
  do ii = 1,coefficients%nterm
      call xmpi_bcast(coefficients%terms(ii)%weight, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%atindx, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%direction, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%cell, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%power_disp, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%power_strain, source, comm, ierr)
      call xmpi_bcast(coefficients%terms(ii)%strain, source, comm, ierr)
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
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_MPIsend(coefficients, tag, dest, comm)

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
    call xmpi_send(coefficients%terms(ii)%nstrain, dest, 9*tag+4, comm, ierr)
  end do

! Transfert value
  do ii = 1,coefficients%nterm
      call xmpi_send(coefficients%terms(ii)%weight, dest, 9*tag+5, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%atindx, dest, 9*tag+6, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%direction, dest, 9*tag+7, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%cell, dest, 9*tag+8, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%power_disp, dest, 9*tag+9, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%power_strain, dest, 9*tag+10, comm, ierr)
      call xmpi_send(coefficients%terms(ii)%strain, dest, 9*tag+11, comm, ierr)
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
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_MPIrecv(coefficients, tag, source, comm)

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
    call xmpi_recv(coefficients%terms(ii)%nstrain, source, 9*tag+4, comm, ierr)
  end do

! Allocate arrays on the other nodes
  do ii = 1,coefficients%nterm
    ABI_ALLOCATE(coefficients%terms(ii)%atindx,(2,coefficients%terms(ii)%ndisp))
    coefficients%terms(ii)%atindx = 0
    ABI_ALLOCATE(coefficients%terms(ii)%direction,(coefficients%terms(ii)%ndisp))
    ABI_ALLOCATE(coefficients%terms(ii)%cell,(3,2,coefficients%terms(ii)%ndisp))
    ABI_ALLOCATE(coefficients%terms(ii)%power_disp,(coefficients%terms(ii)%ndisp))
    ABI_ALLOCATE(coefficients%terms(ii)%power_strain,(coefficients%terms(ii)%nstrain))
    ABI_ALLOCATE(coefficients%terms(ii)%strain,(coefficients%terms(ii)%nstrain))
  end do

! Transfert value
  do ii = 1,coefficients%nterm
    call xmpi_recv(coefficients%terms(ii)%weight, source, 9*tag+5, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%atindx, source, 9*tag+6, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%direction, source, 9*tag+7, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%cell, source, 9*tag+8, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%power_disp, source, 9*tag+9, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%power_strain, source, 9*tag+10, comm, ierr)
    call xmpi_recv(coefficients%terms(ii)%strain, source, 9*tag+11, comm, ierr)
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
!! Copyright (C) 2000-2019 ABINIT group (AM)
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
!! replace = optional, TRUE replace filename if filename exists
!!                     FALSE, default not replace if filename exists
!!
!! OUTPUT
!!
!! PARENTS
!!      m_effective_potential,mover_effpot
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_writeXML(coeffs,ncoeff,filename,unit,newfile,replace)

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in) :: ncoeff
  integer,optional,intent(in) :: unit
  logical,optional,intent(in) :: newfile,replace
!arrays
  type(polynomial_coeff_type), intent(in) :: coeffs(ncoeff)
  character(len=fnlen),optional,intent(in) :: filename
!Local variables-------------------------------
!scalar
 integer :: icoeff,idisp,iterm
 integer :: unit_xml
 logical :: need_header = .TRUE.,need_to_replace = .FALSE.
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

 if(present(replace))then
   need_to_replace = replace
 end if

 if(present(newfile))then
   if (newfile) then
     unit_xml = get_unit()
     need_header = .TRUE.

     if(.not. need_to_replace) call isfile(namefile,'new')
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
!         Atomic displacement case
          select case(coeffs(icoeff)%terms(iterm)%direction(idisp))
          case(1)
            direction ="x"
          case(2)
            direction ="y"
          case(3)
            direction ="z"
          end select
          WRITE(unit_xml,'(a,I0,a,I0,3a,I0,a)') "      <displacement_diff atom_a=""",&
&                         coeffs(icoeff)%terms(iterm)%atindx(1,idisp)-1,""" atom_b=""",&
&                         coeffs(icoeff)%terms(iterm)%atindx(2,idisp)-1,""" direction=""",&
&                         direction,""" power=""",coeffs(icoeff)%terms(iterm)%power_disp(idisp),&
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
        end do
        do idisp=1,coeffs(icoeff)%terms(iterm)%nstrain
!         Strain case
          WRITE(unit_xml,'("      <strain power=""",i2,""" voigt=""",i2,"""/>")')&
&               coeffs(icoeff)%terms(iterm)%power_strain(idisp),&
&               coeffs(icoeff)%terms(iterm)%strain(idisp)
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
!!  energy_coeff(ncoeff) = energy contribution of each anharmonic term 
!!  fcart(3,natom) = contribution  to the forces
!!  strten(6) = contribution to the stress tensor
!!
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE
!!
subroutine polynomial_coeff_evaluate(coefficients,disp,energy,energy_coeff,fcart,natom_sc,natom_uc,ncoeff,sc_size,&
&                                    strain,strten,ncell,index_cells,comm,filename)

!Arguments ------------------------------------
! scalar
  real(dp),intent(out):: energy
  integer, intent(in) :: ncell,ncoeff,natom_sc,natom_uc
  integer, intent(in) :: comm
  character(len=fnlen),optional,intent(in) :: filename 
! array
  real(dp),intent(out):: strten(6)
  real(dp),intent(in) :: strain(6)
  real(dp),intent(out):: fcart(3,natom_sc)
  real(dp),intent(in) :: disp(3,natom_sc)
  real(dp),optional,intent(out):: energy_coeff(ncoeff)
  integer,intent(in) :: index_cells(4,ncell)
  integer,intent(in) :: sc_size(3)
  type(polynomial_coeff_type),intent(in) :: coefficients(ncoeff)
 !Local variables-------------------------------
! scalar
  integer :: i1,i2,i3,ia1,ib1,ia2,ib2,idir1,idir2,ierr,ii
  integer :: icoeff,iterm,idisp1,idisp2,idisp1_strain,idisp2_strain,icell,ndisp
  integer :: nstrain,ndisp_tot,power_disp,power_strain,unit_out
  real(dp):: coeff,disp1,disp2,tmp1,tmp2,tmp3,weight
  logical :: file_opened 
! array
  integer :: cell_atoma1(3),cell_atoma2(3)
  integer :: cell_atomb1(3),cell_atomb2(3)
  character(len=500) :: msg
  character(len=fnlen) :: name_file
! *************************************************************************

! Check
  if (any(sc_size <= 0)) then
    write(msg,'(a,a)')' No supercell found for getEnergy'
    MSG_ERROR(msg)
  end if

  if(present(filename)) name_file = filename

! Initialisation of variables
  energy     = zero
  fcart(:,:) = zero
  strten(:)  = zero
  energy_coeff(:) = zero
  do icell = 1,ncell
    ii = index_cells(4,icell);
    i1=index_cells(1,icell); i2=index_cells(2,icell); i3=index_cells(3,icell)
    ia1 = 0 ; ib1 = 0
!   Loop over coefficients
    do icoeff=1,ncoeff
!     Set the value of the coefficient
      coeff = coefficients(icoeff)%coefficient
!     Loop over terms of this coefficient
      do iterm=1,coefficients(icoeff)%nterm
!       Set the weight of this term
        weight =coefficients(icoeff)%terms(iterm)%weight
        tmp1 = one
        ndisp = coefficients(icoeff)%terms(iterm)%ndisp
        nstrain = coefficients(icoeff)%terms(iterm)%nstrain
        ndisp_tot = ndisp + nstrain
!       Loop over displacement and strain
        do idisp1=1,ndisp_tot
!         Set to one the acculation of forces and strain
          tmp2 = one
          tmp3 = one

!         Strain case idisp > ndisp
          if (idisp1 > ndisp)then
!           Set the power_strain of the strain:
            idisp1_strain = idisp1 - ndisp
            power_strain = coefficients(icoeff)%terms(iterm)%power_strain(idisp1_strain)
!           Get the direction of the displacement or strain
            idir1 = coefficients(icoeff)%terms(iterm)%strain(idisp1_strain)
            if(abs(strain(idir1)) > tol10)then
!             Accumulate energy fo each displacement (\sum ((A_x-O_x)^Y(A_y-O_c)^Z))
              tmp1 = tmp1 * (strain(idir1))**power_strain
              if(power_strain > 1) then
!               Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                tmp3 = tmp3 *  power_strain*(strain(idir1))**(power_strain-1)
              end if
            else
              tmp1 = zero
              if(power_strain > 1) then
                tmp3 = zero
              end if
            end if
          else
!           Set the power_disp of the displacement:
            power_disp = coefficients(icoeff)%terms(iterm)%power_disp(idisp1)
!           Get the direction of the displacement or strain
            idir1 = coefficients(icoeff)%terms(iterm)%direction(idisp1)
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

!             index of the second atom in the (position in the supercell  if the cell is not 0 0 0)
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
              tmp1 = tmp1 * (disp1-disp2)**power_disp
              if(power_disp > 1) then
!               Accumulate forces for each displacement (\sum (Y(A_x-O_x)^Y-1(A_y-O_c)^Z+...))
                tmp2 = tmp2 * power_disp*(disp1-disp2)**(power_disp-1)
              end if
            else
              tmp1 = zero
              if(power_disp > 1) then
                tmp2 = zero
              end if
            end if
          end if

          do idisp2=1,ndisp_tot

            if(idisp2 /= idisp1) then
              if (idisp2 > ndisp)then
                idisp2_strain = idisp2 - ndisp
                idir2 = coefficients(icoeff)%terms(iterm)%strain(idisp2_strain)
!               Strain case
!               Set the power_strain of the strain:
                power_strain = coefficients(icoeff)%terms(iterm)%power_strain(idisp2_strain)
!               Accumulate energy forces
                tmp2 = tmp2 * (strain(idir2))**power_strain
!               Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                tmp3 = tmp3 * (strain(idir2))**power_strain
              else
                idir2 = coefficients(icoeff)%terms(iterm)%direction(idisp2)
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

!               Set the power_disp of the displacement:
                power_disp = coefficients(icoeff)%terms(iterm)%power_disp(idisp2)
                tmp2 = tmp2 * (disp1-disp2)**power_disp
                tmp3 = tmp3 * (disp1-disp2)**power_disp

              end if
            end if
          end do

          if(idisp1 > ndisp)then
!           Accumule stress tensor
            strten(idir1) = strten(idir1) + coeff * weight * tmp3
          else
!           Accumule  forces
            fcart(idir1,ia1) =  fcart(idir1,ia1)  + coeff * weight * tmp2
            fcart(idir1,ib1) =  fcart(idir1,ib1)  - coeff * weight * tmp2
          end if
        end do

        energy_coeff(icoeff) = energy_coeff(icoeff) + coeff * weight * tmp1
!       accumule energy
        energy = energy +  coeff * weight * tmp1

      end do
    end do
  end do


! MPI_SUM
  call xmpi_sum(energy, comm, ierr)
  call xmpi_sum(fcart , comm, ierr)
  call xmpi_sum(strten , comm, ierr)

!Write to anharmonic_energy_terms.out ORIGINAL  
  INQUIRE(FILE=name_file,OPENED=file_opened,number=unit_out)
  if(file_opened .eqv. .TRUE.)then
    do icoeff=1,ncoeff
      call xmpi_sum(energy_coeff(icoeff), comm, ierr)
     ! Marcus write energy contributions of anharmonic terms to file 
      if(icoeff <ncoeff)then      
        write(unit_out,'(A,1ES24.16)',advance='no')  '    ',energy_coeff(icoeff)
      else if(icoeff==ncoeff)then  
        write(unit_out,'(A,1ES24.16)',advance='yes') '    ',energy_coeff(icoeff)
      end if     
    enddo 
  end if

  
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
!! dist(3,natom,natom,nrpt) = distance between atoms atm1 is in the cell 0 0 0
!!                                                   atm2 is in the nrpt cell (see cell(3,nrpt))
!!                            for each component x,y and z
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
!! range_ifc(3) = maximum cut-off for the inter atomic forces constants in each direction
!! sc_size(3) = optional,size of the supercell used for the fit.
!!               For example if you want to fit 2x2x2 cell the interation
!!               Sr-Ti and Sr-Ti[2 0 0] will be identical for the fit process
!!               If check_pbc is true we remove these kind of terms
!!
!! PARENTS
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_getList(cell,crystal,dist,list_symcoeff,list_symstr,&
&                                   natom,nstr_sym,ncoeff_sym,nrpt,range_ifc,sc_size)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 integer,intent(out) :: ncoeff_sym,nstr_sym
!arrays
 integer,intent(in) :: cell(3,nrpt)
 real(dp),intent(in):: dist(3,natom,natom,nrpt)
 type(crystal_t), intent(in) :: crystal
 integer,allocatable,intent(out) :: list_symcoeff(:,:,:),list_symstr(:,:,:)
 integer,optional,intent(in) :: sc_size(3)
 real(dp),intent(in):: range_ifc(3)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff2,icoeff_tot,icoeff_tmp,idisy1,idisy2,ii
 integer :: ipesy1,ipesy2,isym,irpt,irpt3,irpt_ref,irpt_sym
 integer :: jj,jsym,mu,doubles
 integer :: ncoeff,ncoeff2,ncoeff3,ncoeff_max,nu
 integer :: nsym,shift_atm1(3)
 integer :: shift_atm2(3)
 real(dp):: dist_orig,dist_sym,tolsym8
 logical :: found,check_pbc,possible
!arrays
 integer :: isym_rec(3,3),isym_rel(3,3),sc_size_in(3)
 integer :: transl(3),min_range(3),max_range(3)
 integer,allocatable :: blkval(:,:,:,:,:),list(:),list_symcoeff_tmp(:,:,:),list_symcoeff_tmp2(:,:,:)
 integer,allocatable :: list_symstr_tmp(:,:,:),indsym(:,:,:) ,symrec(:,:,:),symrel(:,:,:),list_symcoeff_tmp3(:,:,:)
 integer,allocatable :: index_irred(:) 
 real(dp),allocatable :: tnons(:,:)
 real(dp),allocatable :: wkdist(:),xcart(:,:),xred(:,:),distance(:,:,:)
 real(dp) :: difmin(3)
 real(dp) :: rprimd(3,3)
 real(dp) :: tratom(3)
 character(len=500) :: message



!Initialisation of variables
 irpt_sym = 0
 nsym   = crystal%nsym
 rprimd = crystal%rprimd
 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 xcart(:,:) = crystal%xcart(:,:)
 xred(:,:)  = crystal%xred(:,:)
 ncoeff_max = nrpt*natom*natom*3*3

!Found the ref cell
 irpt_ref = 1
 do irpt=1,nrpt
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
!     exit
   end if
 end do

 !Set the size of the interaction
 check_pbc = .FALSE.
 sc_size_in = 0
 min_range = 0; max_range = 0
 if(present(sc_size))then
   sc_size_in = sc_size
   do mu=1,3
     call findBound_supercell(min_range(mu),max_range(mu),sc_size_in(mu))
   end do
 end if

!Obtain a list of rotated atom labels:
 ABI_ALLOCATE(indsym,(4,nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 ABI_ALLOCATE(symrel,(3,3,nsym))
 ABI_ALLOCATE(tnons,(3,nsym))
 symrec = crystal%symrec
 symrel = crystal%symrel
 tnons  = crystal%tnons

 tolsym8=tol14
 call symatm(indsym,natom,nsym,symrec,tnons,&
&            tolsym8,crystal%typat,crystal%xred)
 ABI_ALLOCATE(blkval,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(list,(natom*nrpt))
 ABI_ALLOCATE(list_symcoeff_tmp,(5,ncoeff_max,nsym))
 ABI_ALLOCATE(wkdist,(natom*nrpt))

!1-Fill strain list
 ABI_ALLOCATE(list_symstr_tmp,(6,nsym,2))
 list_symstr_tmp = 1
 do ia=1,6
   if(list_symstr_tmp(ia,1,1)==0)cycle
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
!    Get the symmetry matrix
     isym_rel(:,:) = crystal%symrel(:,:,isym)
     do idisy1=1,3
       do idisy2=1,3
         if((isym_rel(mu,idisy1)/=0.and.isym_rel(nu,idisy2)/=0)) then
!          Transform to the voig notation
           if(idisy1==idisy2)then
             list_symstr_tmp(ia,isym,1) = idisy1
             list_symstr_tmp(ia,isym,2) = isym_rel(mu,idisy1)
           else
             if(idisy1==1.or.idisy2==1)then
               if(idisy1==2.or.idisy2==2)then
                 list_symstr_tmp(ia,isym,1) = 6
               end if
               if(idisy1==3.or.idisy2==3)then
                 list_symstr_tmp(ia,isym,1) = 5
               end if
             else
               list_symstr_tmp(ia,isym,1) = 4
             end if
           end if
           list_symstr_tmp(ia,isym,2) = isym_rel(mu,idisy1) * isym_rel(nu,idisy2)
         end if
       end do
     end do
!    Remove the symetric
!     if(list_symstr_tmp(ia,isym,1) > ia) then
!       list_symstr_tmp(list_symstr_tmp(ia,isym,1),:,1) = 0
!     end if
   end do
 end do

!Count the number of strain and transfert into the final array
  nstr_sym = 0
  do ia=1,6
    if(list_symstr_tmp(ia,1,1)/=0) nstr_sym = nstr_sym + 1
  end do

 if(allocated(list_symstr))then
   ABI_DEALLOCATE(list_symstr)
 end if
 ABI_ALLOCATE(list_symstr,(nstr_sym,nsym,2))

 icoeff_tmp = 1
 do ia=1,6
   if(list_symstr_tmp(ia,1,1)/=0) then
     list_symstr(icoeff_tmp,:,:) = list_symstr_tmp(ia,:,:)
     icoeff_tmp = icoeff_tmp + 1
   end if
 end do
!END STRAIN

!Compute the distance between each atoms. Indeed the dist array contains the difference of
!cartesian coordinate for each direction
 ABI_ALLOCATE(distance,(natom,natom,nrpt))
 do ia=1,natom
   do ib=1,natom
     do irpt=1,nrpt
       distance(ia,ib,irpt) = ((dist(1,ia,ib,irpt))**2+(dist(2,ia,ib,irpt))**2+&
&                              (dist(3,ia,ib,irpt))**2)**0.5
     end do
   end do
 end do


!Set to one blkval, all the coeff have to be compute
 blkval = 1
 icoeff = 1
 icoeff_tot = 1
 list_symcoeff_tmp = 0

!2-Fill atom list
!Big loop over generic atom
 do ia=1,natom
   wkdist(:)=reshape(distance(ia,:,:),(/natom*nrpt/))
   do ii=1,natom*nrpt
     list(ii)=ii
   end do
   call sort_dp(natom*nrpt,wkdist,list,tol8)
   do ii=1,natom*nrpt
!    Get the irpt and ib
     irpt=(list(ii)-1)/natom+1
     ib=list(ii)-natom*(irpt-1)
     possible = .true.
!Old way with the cut off
!     if(((dist(1,ia,ib,irpt)**2+dist(2,ia,ib,irpt)**2+dist(3,ia,ib,irpt)**2)**0.5) > 9)then
!       possible = .false.
!     end if
     do jj=1,3
!        if(abs(dist(jj,ia,ib,irpt)) - range_ifc(jj)  > tol10.or.&
! &          abs(abs(dist(jj,ia,ib,irpt)) - range_ifc(jj))  < tol10)then
        if(abs(dist(jj,ia,ib,irpt)) - range_ifc(jj)  > tol10)then
       possible = .false.
       end if
     end do

!    If this distance is superior to the cutoff, we don't compute that term
     if(.not.possible)then
       blkval(:,ia,:,ib,irpt)= 0
       if(irpt==irpt_ref)blkval(:,ib,:,ia,irpt)= 0
!        Stop the loop
       cycle
     end if

!    If this coefficient is not possible, we cycle...
     if (all(blkval(:,ia,:,ib,irpt)==0)) cycle

!    Save the distance between the two atoms for futur checks
     dist_orig = (dist(1,ia,ib,irpt)**2+dist(2,ia,ib,irpt)**2+dist(3,ia,ib,irpt)**2)**0.5

     do mu=1,3
       do nu=1,3
!      Check if : - The coefficient is not yet compute
!                 - The directions are the same
!                 - The atoms are not equivalent
         if (mu/=nu) then
           blkval(mu,ia,nu,ib,irpt)=0
           blkval(nu,ia,mu,ib,irpt)=0
           cycle
         end if
!        Pass if the atoms are identical and in the ref cell
         if(irpt==irpt_ref.and.ia==ib) then
           blkval(mu,ia,nu,ib,irpt)=0
           blkval(nu,ib,mu,ia,irpt)=0
           cycle
         end if
         if(blkval(mu,ia,nu,ib,irpt)==1)then
!          Loop over symmetries
           do isym=1,nsym
!            Get the symmetry matrix for this sym
             isym_rec(:,:)  = crystal%symrec(:,:,isym)
             isym_rel(:,:) = crystal%symrel(:,:,isym)
!            Get the corresponding atom and shift with the symetries
!            For atom 1
             ipesy1 = indsym(4,isym,ia)
             shift_atm1 = indsym(1:3,isym,ia)
!            And atom 2
             do jj=1,3 ! Apply transformation to original coordinates.
              tratom(jj) = dble(isym_rec(1,jj))*(xred(1,ib)+cell(1,irpt)-tnons(1,isym))&
&                         +dble(isym_rec(2,jj))*(xred(2,ib)+cell(2,irpt)-tnons(2,isym))&
&                         +dble(isym_rec(3,jj))*(xred(3,ib)+cell(3,irpt)-tnons(3,isym))

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

!            Check the distance
             dist_sym = (dist(1,ipesy1,ipesy2,irpt_sym)**2+&
&                        dist(2,ipesy1,ipesy2,irpt_sym)**2+&
&                        dist(3,ipesy1,ipesy2,irpt_sym)**2)**0.5
             if(abs(dist_orig - dist_sym) > tol10)then
               write(message, '(a,i0,2a,I0,a,es15.8,2a,es15.8,2a)' )&
&                'The distance between the atoms for the coefficient number ',icoeff,ch10,&
&                'with the symmetry ',isym,' is ',dist_sym,ch10,'but the original distance is',&
&                   dist_orig,ch10,&
&                'Action: Contact abinit group'
               MSG_BUG(message)
             end if
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
                   if(isym_rel(mu,idisy1)/=0.and.isym_rel(nu,idisy2)/=0)then
                     if(.not.found.or.(irpt_sym==irpt_ref.and.ipesy1==ipesy2)) then
!                      Remove this term (is not computed) Sr-Sr or not include in the cell
!                      Also remove oposite term... (Srx-Tix) = (Ti-Srx)
                       blkval(idisy1,ipesy1,idisy2,ipesy2,irpt_sym) = 0
                       blkval(idisy2,ipesy2,idisy1,ipesy1,irpt_sym) = 0
                       cycle
                     else
!                      Fill the list with the coeff and symmetric (need all symmetrics)
                       list_symcoeff_tmp(1:4,icoeff,isym)=(/idisy1,ipesy1,ipesy2,irpt_sym/)
!                      Check the sign
                       if(isym_rel(mu,idisy1)/=isym_rel(nu,idisy2))then
                         write(message, '(a,i0,a,I0,4a)' )&
&                        'The sign of coefficient number ',icoeff,' with the symmetry ',isym,ch10,&
&                        'can not be found... Something is going wrong',ch10,&
&                        'Action: Contact abinit group'
                         MSG_BUG(message)
                       end if
                       list_symcoeff_tmp(5,icoeff,isym)= isym_rel(nu,idisy2)
                     end if
                   end if
                 end if
               end do
             end do
           end do ! end loop sym
           icoeff = icoeff + 1
         end if
!        This coeff is now computed
         blkval(mu,ia,nu,ib,irpt)= 0
       end do ! end loop nu
     end do ! end loop mu
   end do ! end loop ii
 end do ! end loop ia

!Reset the output
 if(allocated(list_symcoeff))then
   ABI_DEALLOCATE(list_symcoeff)
 end if

 ABI_DEALLOCATE(distance)

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
 ncoeff = 0
 do icoeff = 1,ncoeff_max
     if(.not.(all(list_symcoeff_tmp(:,icoeff,1)==0)))then
     ncoeff = ncoeff + 1
   end if
 end do

 ABI_ALLOCATE(list_symcoeff_tmp2,(6,ncoeff,nsym))
 list_symcoeff_tmp2 = 0
 icoeff = 0
 do icoeff_tmp = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff_tmp,1)==0)))then
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
&                               ncoeff)
     list_symcoeff_tmp2(6,icoeff,isym) = icoeff2
   end do
 end do

!2.5/do checks
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     if(list_symcoeff_tmp2(6,icoeff,isym)==0)then
       write(message, '(a,i0,a,I0,4a)' )&
&           'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&           'has no equivalent',ch10,&
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


!Check if the atom 2 is not (with the PBC) is in the same cell than
!the atom 1. For example if you want to fit 2x2x2 cell the interation
!Sr-Ti and Sr-Ti[2 0 0] will be identical for the fit process...
!If check_pbc is true we remove these kind of terms
 do icoeff = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp2(:,icoeff,1)==0)))then
     do mu=1,3
       if(min_range(mu) > cell(mu,list_symcoeff_tmp2(4,icoeff,1)) .or. &
&       cell(mu,list_symcoeff_tmp2(4,icoeff,1)) > max_range(mu))then
         list_symcoeff_tmp2(:,icoeff,:)=0
         exit
       end if
     end do
   end if
 end do

!3/ Remove useless terms like opposites
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     !icoeff2 = list_symcoeff_tmp2(6,icoeff,isym)
     !if (icoeff2> icoeff)then
     !  list_symcoeff_tmp2(:,icoeff2,1) = 0
     !end if
     do jsym=1,nsym
       icoeff2 = getCoeffFromList(list_symcoeff_tmp2(:,:,jsym),&
&                                 list_symcoeff_tmp2(3,icoeff,isym),&
&                                 list_symcoeff_tmp2(2,icoeff,isym),&
&                                 list_symcoeff_tmp2(4,icoeff,isym),&
&                                 list_symcoeff_tmp2(1,icoeff,isym),&
&                                 ncoeff)
       if (icoeff2> icoeff)then
         list_symcoeff_tmp2(:,icoeff2,1) = 0
       end if
     end do
   end do 
   !MS only keep terms with ia == 1 
   !if (list_symcoeff_tmp2(2,icoeff,1) /= 1)then 
   !    list_symcoeff_tmp2(:,icoeff,1) = 0 
   !endif
 end do

!4/ Recount the number of coeff after step 3
 ncoeff2 = 0
 doubles  = 0
 isym = 0
write(std_out,*) "DEBUG ncoeff2 after 2.5.3: ", ncoeff2
 do icoeff = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp2(:,icoeff,1)==0)))then
     ncoeff2 = ncoeff2 + 1
     !do isym = 1,nsym 
     !  icoeff2 = list_symcoeff_tmp2(6,icoeff,isym)
     !  if (icoeff2> icoeff)then
     !    doubles = doubles +1
     !  write(std_out,*) "DEBUG: doubles: ", doubles
     !  end if
     !end do 
   end if
 end do

write(std_out,*) "DEBUG ncoeff2 after 2.5.3: ", ncoeff2

 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_ALLOCATE(list_symcoeff_tmp,(6,ncoeff,nsym))
 list_symcoeff_tmp = list_symcoeff_tmp2

!4.1 Count irreducible terms
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     icoeff2 = list_symcoeff_tmp2(6,icoeff,isym)
     if (icoeff2> icoeff)then
       list_symcoeff_tmp(:,icoeff2,1) = 0
     end if
   end do 
 end do 

ncoeff3 = 0
 do icoeff = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp(:,icoeff,1)==0)))then
     ncoeff3 = ncoeff3 + 1
   end if
 end do
write(std_out,*) "DEBUG ncoeff3 after 2.5.3: ", ncoeff3

!4.2 Put irreducible terms in front of list_symcoeff_tmp3 
!    Store index of irreducible terms
 ABI_ALLOCATE(list_symcoeff_tmp3,(6,ncoeff2,nsym))
 ABI_ALLOCATE(index_irred,(ncoeff3))
icoeff_tmp = 0 
do icoeff = 1,ncoeff    
   if(.not.(all(list_symcoeff_tmp(:,icoeff,1)==0)))then
     icoeff_tmp = icoeff_tmp+1
     index_irred(icoeff_tmp) = icoeff 
     list_symcoeff_tmp3(:,icoeff_tmp,:) = list_symcoeff_tmp(:,icoeff,:) 
   endif 
enddo

!4.3 Put symmetric equivalents behind in list_symcoeff_tmp3
!    Attention icoeff_tmps keeps it's value of loop before 
!    TODO for check should be equal to ncoeff2
do icoeff = 1,ncoeff     
   if(.not.(all(list_symcoeff_tmp2(:,icoeff,1)==0))&
&     .and..not. any(index_irred == icoeff))then
     icoeff_tmp = icoeff_tmp+1
     list_symcoeff_tmp3(:,icoeff_tmp,:) = list_symcoeff_tmp2(:,icoeff,:) 
   endif 
enddo   

!4.4 A little copy round
ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_ALLOCATE(list_symcoeff_tmp,(6,ncoeff2,nsym))
list_symcoeff_tmp = list_symcoeff_tmp3

!5/ Final transfert
 ABI_ALLOCATE(list_symcoeff,(6,ncoeff2,nsym))
 list_symcoeff = 0
 icoeff = 0
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
&                               ncoeff2)
     list_symcoeff(6,icoeff,isym) = icoeff2
!     list_symcoeff(6,icoeff2,isym) = icoeff
   end do
 end do

!Set the max number of coeff inside list_symcoeff
 ncoeff_sym = ncoeff3

!Deallocation
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(list)
 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_DEALLOCATE(list_symcoeff_tmp2)
 ABI_DEALLOCATE(list_symcoeff_tmp3) 
 ABI_DEALLOCATE(list_symstr_tmp)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(symrel)
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
!! power_disps(2) = array with the minimal and maximal power_disp to be computed
!! max_power_strain = maximum order of the strain of the strain phonon coupling
!! option = 0 compute all terms
!!          1 still in development
!! sc_size(3) = size of the supercell used for the fit.
!!               For example if you want to fit 2x2x2 cell the interation
!!               Sr-Ti and Sr-Ti[2 0 0] will be identical for the fit process
!!               If check_pbc is true we remove these kind of terms
!! comm = MPI communicator
!! anharmstr = logical, optional : TRUE, the anharmonic strain is computed (\eta)^power_disp ...
!!                                   FALSE, (default) the anharmonic strain are not computed
!! spcoupling= logical, optional : TRUE(default) the anharmonic strain-phonon coupling is computed
!! only_odd_power = logical, optional : if TRUE return only odd power
!! only_even_power= logical, optional : if TRUe return only even power
!! distributed = logical, optional : True, the coefficients will be distributed on the CPU
!! verbose  = optional, flag for the verbose mode
!!
!! OUTPUT
!! polynomial_coeff<(type(polynomial_coeff_type)>(ncoeff) = array of datatype with the polynomial_coeff
!! ncoeff = number of coefficients for this CPU if distributed == true, all otherwise
!! ncoeff_tot = total number of coefficient over the CPU
!!
!! PARENTS
!!      m_fit_polynomial_coeff,mover_effpot
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_getNorder(coefficients,crystal,cutoff,ncoeff,ncoeff_tot,power_disps,&
&                                     max_power_strain,option,sc_size,comm,anharmstr,spcoupling,&
&                                     distributed,only_odd_power,only_even_power,verbose)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: max_power_strain,option,comm
 integer,intent(out):: ncoeff,ncoeff_tot
 real(dp),intent(in):: cutoff
 logical,optional,intent(in) :: anharmstr,spcoupling,distributed,verbose
 logical,optional,intent(in) :: only_odd_power,only_even_power
!arrays
 integer,intent(in) :: power_disps(2),sc_size(3)
 type(crystal_t), intent(inout) :: crystal
 type(polynomial_coeff_type),allocatable,intent(inout) :: coefficients(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff2,icoeff3,ierr,ii,ij,irpt,irpt_ref,iterm,i0
 integer :: lim1,lim2,lim3,tmp_ind1,tmp_ind2,i,j,k,cnt,idisp,isym,istrain
 integer :: master,my_rank,my_ncoeff,my_newncoeff,natom,ncombination,ncoeff_max,ncoeff_sym
 integer :: ncoeff_symsym
 integer :: ncoeff_alone,ndisp_max,nproc,nrpt,nsym,nterm,nstr_sym,r1,r2,r3,my_size
 integer :: my_icoeff,rank_to_send,rank_to_receive,rank_to_send_save
 real(dp):: norm
 logical :: iam_master,need_anharmstr,need_spcoupling,need_distributed,need_verbose
 logical :: need_only_odd_power,need_only_even_power
!arrays
 integer :: ncell(3)
 integer,allocatable :: buffsize(:),buffdispl(:)
 integer,allocatable :: cell(:,:),compatibleCoeffs(:,:),list_combination_cmp_tmp(:)
 integer,allocatable :: list_symcoeff(:,:,:),list_symstr(:,:,:),list_coeff(:),list_combination(:,:)
 integer,allocatable :: list_combination_tmp(:,:),list_combination_tmp2(:,:)
 integer,allocatable  :: my_coefflist(:),my_coeffindexes(:),my_newcoeffindexes(:)
 real(dp) :: rprimd(3,3),range_ifc(3)
 real(dp),allocatable :: dist(:,:,:,:),rpt(:,:)
 real(dp),allocatable :: xcart(:,:),xred(:,:)
 character(len=5),allocatable :: symbols(:)
 character(len=200):: name
 character(len=500) :: message
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs_tmp
 type(polynomial_term_type),dimension(:),allocatable :: terms
! character(len=fnlen) :: filename
! *************************************************************************

 !MPI variables
 master = 0
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Free the output
 if(allocated(coefficients))then
   do ii =1,size(coefficients)
     call polynomial_coeff_free(coefficients(ii))
   end do
   ABI_DEALLOCATE(coefficients)
 end if

!Check
 if(option > power_disps(2))then
   write(message, '(3a)' )&
&       'Option can not be superior to the maximum order ',ch10,&
&       'Action: contact abinit group'
   MSG_ERROR(message)
 end if

!Initialisation of variables
 need_anharmstr = .TRUE.
 if(present(anharmstr)) need_anharmstr = anharmstr
 need_spcoupling = .TRUE.
 if(present(spcoupling)) need_spcoupling = spcoupling
 need_distributed = .FALSE.
 if(present(distributed)) need_distributed  = distributed
 need_verbose = .TRUE.
 if(present(verbose)) need_verbose = verbose
 need_only_odd_power = .FALSE.
 if(present(only_odd_power)) need_only_odd_power = only_odd_power
 need_only_even_power = .FALSE.
 if(present(only_even_power)) need_only_even_power = only_even_power

 if(need_only_odd_power.and.need_only_even_power)then
      write(message, '(3a)' )&
&       'need_only_odd_power and need_only_even_power are both true',ch10,&
&       'Action: contact abinit group'
   MSG_ERROR(message)
 end if

 natom  = crystal%natom
 nsym   = crystal%nsym
 rprimd = crystal%rprimd

 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 xcart(:,:) = crystal%xcart(:,:)
 xred(:,:)  = crystal%xred(:,:)

!Compute the max range of the ifc with respect to the trainning set
 range_ifc(:) = zero
 do ii=1,3
   norm = sqrt(rprimd(ii,1)**2+ rprimd(ii,2)**2+rprimd(ii,3)**2)
   range_ifc(ii) = range_ifc(ii) + norm * sc_size(ii) / 2.0
 end do

!compute new ncell
 ncell = sc_size
 lim1=((ncell(1)/2)) + 1
 lim2=((ncell(2)/2)) + 1
 lim3=((ncell(3)/2)) + 1
 if(mod(ncell(1),2)/=0) lim1=lim1+1
 if(mod(ncell(2),2)/=0) lim2=lim2+1
 if(mod(ncell(3),2)/=0) lim3=lim3+1
 nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)

 ncell(1) = 2*lim1+1
 ncell(2) = 2*lim2+1
 ncell(3) = 2*lim3+1

 !Build the rpt point
 ABI_ALLOCATE(rpt,(3,nrpt))
 ABI_ALLOCATE(cell,(3,nrpt))

!WARNING:
!Put the reference cell into the first element
!the code will first deal with the atoms of the first cell
 irpt = 1
 irpt_ref = 1
 rpt(:,1) = zero
 cell(:,irpt)=0
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
!Now dist(3,ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.
 ABI_ALLOCATE(dist,(3,natom,natom,nrpt))
 dist = zero
 do ia=1,natom
   do ib=1,natom
     do irpt=1,nrpt
       dist(1,ia,ib,irpt) = xcart(1,ib)-xcart(1,ia)+rpt(1,irpt)
       dist(2,ia,ib,irpt) = xcart(2,ib)-xcart(2,ia)+rpt(2,irpt)
       dist(3,ia,ib,irpt) = xcart(3,ib)-xcart(3,ia)+rpt(3,irpt)
     end do
   end do
 end do


 if(need_verbose)then
   write(message,'(1a)')' Generation of the list of all the possible coefficients'
   call wrtout(std_out,message,'COLL')
 end if
 call polynomial_coeff_getList(cell,crystal,dist,list_symcoeff,list_symstr,&
&                              natom,nstr_sym,ncoeff_sym,nrpt,range_ifc,sc_size=sc_size)

if(iam_master)then 
i0 = 0 
write(std_out,*) "DEBUG shape list_symcoeff(:,:,:):", shape(list_symcoeff) 
write(std_out,*) "DEBUG shape size(list_symcoeff,2):", size(list_symcoeff(1,:,1))
ncoeff_symsym = size(list_symcoeff(1,:,1))
write(std_out,*) "DEBUG list_symcoeff(:,11,:) for coeff 1:" 
do ii=1,nsym 
  write(std_out,*) "******************"
  write(std_out,*) "nsym: ", ii 
  write(std_out,*) "------------------" 
  write(std_out,*) "direction list_symcoeff(1,11,", ii,"): ", list_symcoeff(1,11,ii)  
  write(std_out,*) "atom1     list_symcoeff(2,11,", ii,"): ", list_symcoeff(2,11,ii)  
  write(std_out,*) "atom2     list_symcoeff(3,11,", ii,"): ", list_symcoeff(3,11,ii)  
  write(std_out,*) "irpt      list_symcoeff(4,11,", ii,"): ", list_symcoeff(4,11,ii) !,"cell: ", cell(:,list_symcoeff(4,1,ii))
  write(std_out,*) "weight    list_symcoeff(5,11,", ii,"): ", list_symcoeff(5,11,ii)  
  write(std_out,*) "index sym list_symcoeff(6,11,", ii,"): ", list_symcoeff(6,11,ii) 
  ia = list_symcoeff(6,11,ii)  
  if(ia /= 0)then 
  write(std_out,*) "------------------" 
  write(std_out,*) "Symmetric Term from list_symcoeff:" 
  write(std_out,*) "direction list_symcoeff(1,", ia,",1): ", list_symcoeff(1,ia,1)  
  write(std_out,*) "atom1     list_symcoeff(2,", ia,",1): ", list_symcoeff(2,ia,1)  
  write(std_out,*) "atom2     list_symcoeff(3,", ia,",1): ", list_symcoeff(3,ia,1)  
  write(std_out,*) "irpt      list_symcoeff(4,", ia,",1): ", list_symcoeff(4,ia,1) !,"cell: ", cell(:,list_symcoeff(4,1,ii))
  write(std_out,*) "weight    list_symcoeff(5,", ia,",1): ", list_symcoeff(5,ia,1)  
  write(std_out,*) "index sym list_symcoeff(6,", ia,",",ii,"): ", list_symcoeff(6,ia,ii)  
  write(std_out,*) "index sym list_symcoeff(6,", ia,",1): ", list_symcoeff(6,ia,1)  
  else 
    i0 = i0 + 1
  endif 
end do 

write(std_out,*) "******************"
!write(std_out,*) "i0 is: ", i0 
! 
!write(std_out,*) "DEBUG WHAT ?! " 
!do ii=1,ncoeff_sym
!   write(std_out,*) "******************"
!   write(std_out,*) "icoeff: ", ii 
!   write(std_out,*) "------------------" 
!   write(std_out,*) "direction list_symcoeff(1,", ii,",1): ", list_symcoeff(1,ii,1)  
!   write(std_out,*) "atom1     list_symcoeff(2,", ii,",1): ", list_symcoeff(2,ii,1)  
!   write(std_out,*) "atom2     list_symcoeff(3,", ii,",1): ", list_symcoeff(3,ii,1)  
!   write(std_out,*) "irpt      list_symcoeff(4,", ii,",1): ", list_symcoeff(4,ii,1) !,"cell: ", cell(:,list_symcoeff(4,1,ii))
!   write(std_out,*) "weight    list_symcoeff(5,", ii,",1): ", list_symcoeff(5,ii,1)  
!   write(std_out,*) "index sym list_symcoeff(6,", ii,",1): ", list_symcoeff(6,ii,1) 
!   write(std_out,*) "------------------" 
!enddo
endif

 ABI_DEALLOCATE(dist)
 ABI_DEALLOCATE(rpt)

!Compute the total number of coefficient
 ncoeff_tot = ncoeff_sym+nstr_sym

 if(iam_master)then
!Check the distanceance bewteen coefficients and store integer:
! 0: the mix between these coefficient is not possible
! 1: the mix between these coefficient is possible
   ABI_ALLOCATE(compatibleCoeffs,(ncoeff_tot,ncoeff_tot))
   compatibleCoeffs(:,:) = 1

   if(need_verbose)then
     write(message,'(1a)')' Check the compatible coefficients with respect to the cutoff'
     call wrtout(std_out,message,'COLL')
   end if

write(std_out,*) "DEBUG: What is ncoeff_tot: ", ncoeff_tot
   do icoeff=1,ncoeff_tot
     do icoeff2=1,ncoeff_tot
!      Select case:
!      if both icoeff are displacement => check the distance
!      if both icoeff are strain => check the flag
!      Otherwise cycle (we keep the term)
       if(icoeff>ncoeff_sym.and.icoeff2<=ncoeff_sym)cycle
       if(icoeff2<=ncoeff_sym.and.icoeff2>ncoeff_sym)cycle
       if((icoeff>ncoeff_sym.or.icoeff2>ncoeff_sym).and.&
&       .not.need_anharmstr.and..not.need_spcoupling) then
         compatibleCoeffs(icoeff,icoeff2) = 0
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
!TEST_AM
           compatibleCoeffs(icoeff,icoeff2) = 0
           compatibleCoeffs(icoeff2,icoeff) = 0
!TEST_AM
         end if
       end if
     end do !end  icoeff 
   end do !icoeff2
ii = 0
   do icoeff=1,ncoeff_tot
     do icoeff2=1,ncoeff_tot
         if(compatibleCoeffs(icoeff,icoeff2) .ne. 0)then 
		ii = ii + 1
!               write(std_out,*) "ii: ", ii
         end if  
     end do 
   end do 

write(std_out,*) "DEBUG: Number of Compatible Coeffs: ", ii

!  Compute all the combination of coefficient up to the given order  (get the number)
   if(need_verbose)then
     write(message,'(1a)')' Compute the number of possible combinations'
     call wrtout(std_out,message,'COLL')
   end if

   ABI_ALLOCATE(list_coeff,(0))
   ABI_ALLOCATE(list_combination,(0,0))
   icoeff  = 1
   icoeff2 = 0
   call computeCombinationFromList(cell,compatibleCoeffs,list_symcoeff,list_symstr,&
&                   list_coeff,list_combination,icoeff,max_power_strain,icoeff2,natom,ncoeff_sym,&
&                   ncoeff_symsym,nstr_sym,icoeff,nrpt,nsym,1,power_disps(1),power_disps(2),symbols,nbody=option,&
&                   compute=.false.,anharmstr=need_anharmstr,spcoupling=need_spcoupling,&
&                   only_odd_power=need_only_odd_power,only_even_power=need_only_even_power)
   ncombination  = icoeff2
   ABI_DEALLOCATE(list_coeff)
   ABI_DEALLOCATE(list_combination)

!  Compute all the combination of coefficient up to the given order
   if(need_verbose)then
     write(message,'(1a)')' Compute the combinations'
     call wrtout(std_out,message,'COLL')
   end if
 
   write(std_out,*) "DEBUG combinations: ", ncombination
   ABI_ALLOCATE(list_coeff,(0))
   ABI_ALLOCATE(list_combination_tmp,(power_disps(2),ncombination))
   write(std_out,*) "DEBUG shape list_combination_tmp: ", shape(list_combination_tmp)
   icoeff  = 1
   icoeff2 = 0
   list_combination_tmp = 0
   write(std_out,*) "DEBUG: list_coeff before first call to CCL:", list_coeff
   write(std_out,*) "DEBUG: ncoeff:", ncoeff_sym
   call computeCombinationFromList(cell,compatibleCoeffs,list_symcoeff,list_symstr,&
&                   list_coeff,list_combination_tmp,icoeff,max_power_strain,icoeff2,natom,&
&                   ncoeff_sym,ncoeff_symsym,nstr_sym,ncombination,nrpt,nsym,1,power_disps(1),power_disps(2),&
&                   symbols,nbody=option,compute=.true.,&
&                   anharmstr=need_anharmstr,spcoupling=need_spcoupling,&
&                   only_odd_power=need_only_odd_power,only_even_power=need_only_even_power)
   ABI_DEALLOCATE(list_coeff)
   ABI_DEALLOCATE(compatibleCoeffs)
!   write(std_out,*) "DEBUG: list_combination after call to CCL:", list_combination_tmp
 else
   ABI_ALLOCATE(list_combination_tmp,(1,1))
 end if

write(std_out,*) "DEBUG: reduce zero combinations"
 !Reduce zero combinations
 i0 = 0 
 do ii=1,ncombination
   if(any(list_combination_tmp(:,ii) /= 0))then 
      i0 = i0 + 1
   endif  
 enddo

 ABI_ALLOCATE(list_combination_tmp2,(power_disps(2),i0))
 i0 = 0 
 do ii=1,ncombination
   if(any(list_combination_tmp(:,ii) /= 0))then 
      i0 = i0+1
      list_combination_tmp2(:,i0) = list_combination_tmp(:,ii)
   endif  
 enddo
 ncombination = i0
 ABI_DEALLOCATE(list_combination_tmp)
 ABI_ALLOCATE(list_combination_tmp,(power_disps(2),i0))
 list_combination_tmp = list_combination_tmp2 
 ABI_DEALLOCATE(list_combination_tmp2)
write(std_out,*) "DEBUG: finish reduce zero combinations"
!write(std_out,*) "DEBUG: shape list_combination after reduction:", shape(list_combination_tmp)
 
i0=0
write(std_out,*) "DEBUG: start reduce doubles"
write(std_out,*) "DEBUG: start sorting combinations"
!Reduce double combinations 
!1. sort all arrays
do i=1,ncombination 
!   write(std_out,*) "DEBUG list_combination_tmp before sort:", list_combination_tmp(:,i)
   j=2 
   k=2
   do while(j <= size(list_combination_tmp(:,i))) 
      k = j
      cnt = 1
      do while(k >= 2 .and. cnt == 1) 
         if(list_combination_tmp(k-1,i) > list_combination_tmp(k,i) .and. list_combination_tmp(k,i) > 0)then 
              tmp_ind1 = list_combination_tmp(k-1,i) 
              tmp_ind2 = list_combination_tmp(k,i)
              list_combination_tmp(k,i) = tmp_ind1
              list_combination_tmp(k-1,i) = tmp_ind2
              k=k-1
         else 
           cnt = cnt +1
         end if
      end do
      j = j+1
   end do 
!   write(std_out,*) "DEBUG list_combination_tmp after sort:", list_combination_tmp(:,i)
enddo
write(std_out,*) "DEBUG: finish sorting combinations"
write(std_out,*) "DEBUG: start counting doubles"
!2.Figure out equal ones 
ABI_ALLOCATE(list_combination_cmp_tmp,(power_disps(2)))
i0 = 0
do i=1,ncombination 
   if(any(list_combination_tmp(:,i) /= 0))then 
     do j=i+1,ncombination
        !If term j is equivalent to term i delet it
        if(all(list_combination_tmp(:,i) == list_combination_tmp(:,j)))then
           list_combination_tmp(:,j) = 0
           i0 = i0 + 1
        else !else loop over symmetries to fina symmetry operation that projects term j on i 
          isym = 2
          do while(isym <= nsym)
             list_combination_cmp_tmp(1)=list_combination_tmp(1,j)
             !Get equivalent term indexes for symmetry isym
             do idisp=2,power_disps(2)
                if(list_combination_tmp(idisp,j) /= 0 .and. list_combination_tmp(idisp,j) <= ncoeff_symsym )then
                   list_combination_cmp_tmp(idisp)=list_symcoeff(6,list_combination_tmp(idisp,j),isym)
                else if(list_combination_tmp(idisp,j) > ncoeff_symsym)then
                   istrain = list_combination_tmp(idisp,j) - ncoeff_symsym  
                   list_combination_cmp_tmp(idisp)=list_symstr(istrain,isym,1)
                else 
                   list_combination_cmp_tmp(idisp) = 0 
                endif 
             enddo
!   write(std_out,*) "DEBUG list_combination_tmp before sort:", list_combination_tmp(:,i)
             !Sort the new symmetric indexes for comparision
             ij=2 
             k=2
             do while(ij <= size(list_combination_cmp_tmp(:))) 
                k = ij
                cnt = 1
                do while(k >= 2 .and. cnt == 1) 
                   if(list_combination_cmp_tmp(k-1) > list_combination_cmp_tmp(k) .and. list_combination_cmp_tmp(k) > 0)then 
                        tmp_ind1 = list_combination_cmp_tmp(k-1) 
                        tmp_ind2 = list_combination_cmp_tmp(k)
                        list_combination_cmp_tmp(k) = tmp_ind1
                        list_combination_cmp_tmp(k-1) = tmp_ind2
                        k=k-1
                   else 
                     cnt = cnt +1
                   end if
                end do ! whlie k>=2
                ij = ij+1
             end do ! while ij< size(list..)
             !Compare. If equivalent delete term j
             if(all(list_combination_tmp(:,i) == list_combination_cmp_tmp(:)))then
                list_combination_tmp(:,j) = 0
                i0 = i0 + 1
                isym = nsym +1 
             else 
                isym = isym + 1
             endif
          enddo !isym 
        endif ! all(list_combinaton...)
     enddo ! j=i+1,ncombination
   end if ! any(list_combination_tmp /= 0 ) 
end do !i=1,ncombination 
ABI_DEALLOCATE(list_combination_cmp_tmp)
write(std_out,*) "DEBUG: finish counting doubles"
write(std_out,*) "DEBUG: Transfer irreducible ones"
 ABI_ALLOCATE(list_combination_tmp2,(power_disps(2),ncombination-i0))
 i0 = 0 
 do ii=1,ncombination
   if(any(list_combination_tmp(:,ii) /= 0))then 
      i0 = i0+1
      list_combination_tmp2(:,i0) = list_combination_tmp(:,ii)
   endif  
 enddo
 write(std_out,*) "DEBUG: finish reduce doubles"
 ncombination = i0
 ABI_DEALLOCATE(list_combination_tmp)
 ABI_ALLOCATE(list_combination_tmp,(power_disps(2),i0))
 list_combination_tmp = list_combination_tmp2 
 ABI_DEALLOCATE(list_combination_tmp2)
 


 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred)

!MPI
 if(need_verbose .and. nproc > 1)then
   write(message,'(1a)')' Distribute the combinations over the CPU'
   call wrtout(std_out,message,'COLL')
 end if

 call xmpi_bcast(ncombination, master, comm, ierr)

 ncoeff_alone = mod(ncombination,nproc)
 my_ncoeff = int(aint(real(ncombination,sp)/(nproc)))

 if(my_rank >= (nproc-ncoeff_alone)) then
   my_ncoeff = my_ncoeff  + 1
 end if

!Set the buffsize for mpi scatterv
 ABI_ALLOCATE(buffsize,(nproc))
 ABI_ALLOCATE(buffdispl,(nproc))
 do ii = 1,nproc
   buffsize(ii) = int(aint(real(ncombination,sp)/(nproc))*power_disps(2))
   if(ii > (nproc-ncoeff_alone)) then
     buffsize(ii) = buffsize(ii) + power_disps(2)
   end if
 end do

 buffdispl(1) = 0
 do ii = 2,nproc
   buffdispl(ii) = buffdispl(ii-1) + buffsize(ii-1)
 end do

 ABI_ALLOCATE(list_combination,(power_disps(2),my_ncoeff))
 list_combination = 0

 my_size = my_ncoeff*power_disps(2)
 call xmpi_scatterv(list_combination_tmp,buffsize,buffdispl,list_combination,my_size,master,&
&                   comm,ierr)

 ABI_DEALLOCATE(buffdispl)
 ABI_DEALLOCATE(buffsize)
 ABI_DEALLOCATE(list_combination_tmp)

write(std_out,*) "DEBUG shape list_symcoeff(:,:,:) before termsfromlist:", shape(list_symcoeff) 
 !Compute the coefficients from the list of combination
 if(need_verbose .and. nproc > 1)then
   write(message,'(1a)')' Compute the coefficients'
   call wrtout(std_out,message,'COLL')
 end if
 ABI_ALLOCATE(coeffs_tmp,(my_ncoeff))
 nterm      = nsym
 ndisp_max  = power_disps(2)
 ncoeff_max = my_ncoeff
 ABI_ALLOCATE(terms,(nterm))
 do ii=1,my_ncoeff
   write(std_out,*) "DEBUG list_combination(:,",ii,"): ", list_combination(:,ii)
   call generateTermsFromList(cell,list_combination(:,ii),list_symcoeff,list_symstr,ncoeff_symsym,&
&                             ndisp_max,nrpt,nstr_sym,nsym,nterm,terms)
   call polynomial_coeff_init(one,nterm,coeffs_tmp(ii),terms(1:nterm),check=.true.)
   !  Free the terms array
   do iterm=1,nterm
     call polynomial_term_free(terms(iterm))
   end do
 end do
 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(cell)


 ABI_DEALLOCATE(list_combination)
 ABI_DEALLOCATE(list_symcoeff)
 ABI_DEALLOCATE(list_symstr)

!Final tranfert
!1- Count the total number of coefficient
 ncoeff = 0
 do icoeff=1,ncoeff_max
   if (abs(coeffs_tmp(icoeff)%coefficient) >tol16) then
     ncoeff = ncoeff + 1
   end if
 end do

!Get the total number of coefficients
!ncoeff_max is the number of total coefficients before the symetries check
!ncoeff_tot is the number of total coefficients after the symetries check
 ncoeff_tot = ncoeff!set the output
 call xmpi_sum(ncoeff_tot,comm,ierr)
 call xmpi_sum(ncoeff_max,comm,ierr)


!Need to redistribute the coefficients over the CPU
!Get the list with the number of coeff on each CPU
!In order to be abble to compute the my_coeffindexes array which is for example:
! if CPU0 has 200  Coeff and CPU1 has 203 Coeff then
! for CPU0:my_coeffindexes=>1-200 and for CPU1:my_coeffindexes=>201-403
 if(need_verbose .and. nproc > 1)then
   write(message,'(1a)')' Redistribute the coefficients over the CPU'
   call wrtout(std_out,message,'COLL')
 end if

 ABI_ALLOCATE(buffdispl,(nproc))
 buffdispl = 0
 buffdispl(my_rank+1) = my_ncoeff
 call xmpi_sum(buffdispl,comm,ierr)
 ABI_ALLOCATE(my_coeffindexes,(my_ncoeff))
 ABI_ALLOCATE(my_coefflist,(my_ncoeff))
 my_coeffindexes = 0
 my_coefflist = 0
 do icoeff=1,my_ncoeff
   my_coefflist(icoeff) = icoeff
   if(my_rank==0) then
     my_coeffindexes(icoeff) = icoeff
   else
     my_coeffindexes(icoeff) = sum(buffdispl(1:my_rank)) + icoeff
   end if
 end do
 ABI_DEALLOCATE(buffdispl)

!Compute the new number of coefficient per CPU
 if(need_distributed) then
   ncoeff_alone = mod(ncoeff_tot,nproc)
   my_newncoeff = int(aint(real(ncoeff_tot,sp)/(nproc)))
   if(my_rank >= (nproc-ncoeff_alone)) then
     my_newncoeff = my_newncoeff  + 1
   end if
 else
   my_newncoeff = ncoeff_tot
 end if

 ncoeff = my_newncoeff ! Set the output

!2:compute the number of coefficients and the list of the corresponding
!  coefficients for each CPU.
 ABI_ALLOCATE(my_newcoeffindexes,(my_newncoeff))
 if(need_distributed) then
   do icoeff=1,my_newncoeff
     if(my_rank >= (nproc-ncoeff_alone))then
       my_newcoeffindexes(icoeff)=int(aint(real(ncoeff_tot,sp)/(nproc)))*(my_rank)+&
&                              (my_rank - (nproc-ncoeff_alone)) + icoeff
     else
       my_newcoeffindexes(icoeff)=(my_newncoeff)*(my_rank)  + icoeff
     end if
   end do
 else
   do icoeff=1,my_newncoeff
     my_newcoeffindexes(icoeff) = icoeff
   end do
 end if

!2- Transfer
 if(.not.need_distributed)then
   if(.not.allocated(coefficients))then
     ABI_DATATYPE_ALLOCATE(coefficients,(my_newncoeff))
   end if
 end if
 icoeff  = 0! icoeff is the current index in the total list of coefficients
 icoeff2 = 0! icoeff2 is the current index in the output coefficients array on each CPU
 icoeff3 = 0! icoeff3 is the current index in total new list of coefficients
 rank_to_send_save = 0

 do icoeff=1,ncoeff_max
!  Need to send the rank with the chosen coefficient
   rank_to_send = 0
   my_icoeff = 0
   do ii=1,my_ncoeff
     if (my_coeffindexes(ii)==icoeff) then
       my_icoeff = ii
       if (abs(coeffs_tmp(my_icoeff)%coefficient) > tol16)then
         rank_to_send = my_rank
       else
         rank_to_send = -1
!        Free the coefficient
         call polynomial_coeff_free(coeffs_tmp(ii))
       end if
       exit
     end if
   end do
   call xmpi_sum(rank_to_send, comm, ierr)
!  This coefficient is not compute
   if (rank_to_send == -1) cycle

!  increase icoeff3
   icoeff3 = icoeff3 + 1

!  Find the receiver CPU
   rank_to_receive = 0
   do ii=1,my_newncoeff
     if (my_newcoeffindexes(ii)==icoeff3) then
       rank_to_receive = my_rank
     end if
   end do
   call xmpi_sum(rank_to_receive, comm, ierr)

   if(need_distributed.and.rank_to_send /= rank_to_send_save) then
     if(my_rank == rank_to_send_save)then
       ABI_DATATYPE_DEALLOCATE(coeffs_tmp)!Free memory if the current CPU has already distribute
                                          !all its own coefficients
     end if
     rank_to_send_save = rank_to_send
   end if

   if(need_distributed.and.my_rank == rank_to_receive)then
     if(.not.allocated(coefficients))then
       ABI_DATATYPE_ALLOCATE(coefficients,(my_newncoeff))
     end if
   end if


   if (need_distributed)then
     if(my_rank==rank_to_send)then
       if(any(my_newcoeffindexes(:)==icoeff3))then
         icoeff2 = icoeff2 + 1
!        Get the name of this coefficient
         call polynomial_coeff_getName(name,coeffs_tmp(my_icoeff),symbols,recompute=.TRUE.)
         call polynomial_coeff_init(one,coeffs_tmp(my_icoeff)%nterm,coefficients(icoeff2),&
&                                  coeffs_tmp(my_icoeff)%terms,name=name,check=.false.)
       else
         call polynomial_coeff_MPIsend(coeffs_tmp(my_icoeff), icoeff, rank_to_receive, comm)
       end if
!      Free the coefficient
       call polynomial_coeff_free(coeffs_tmp(my_icoeff))
     else
       if(any(my_newcoeffindexes(:)==icoeff3))then
         icoeff2 = icoeff2 + 1
         call polynomial_coeff_MPIrecv(coefficients(icoeff2), icoeff, rank_to_send, comm)
         call polynomial_coeff_getName(name,coefficients(icoeff2),symbols,recompute=.TRUE.)
         call polynomial_coeff_SetName(name,coefficients(icoeff2))
       end if
     end if
   else
     icoeff2 = icoeff2 + 1
!    Get the name of this coefficient
     if(my_rank==rank_to_send)then
       call polynomial_coeff_getName(name,coeffs_tmp(my_icoeff),symbols,recompute=.TRUE.)
       call polynomial_coeff_init(one,coeffs_tmp(my_icoeff)%nterm,coefficients(icoeff2),&
&                                 coeffs_tmp(my_icoeff)%terms,name=name,check=.false.)
!      Free the coefficient
       call polynomial_coeff_free(coeffs_tmp(my_icoeff))
     end if
     call polynomial_coeff_broadcast(coefficients(icoeff2),rank_to_send, comm)
   end if
 end do

!  if(iam_master)then
!    filename = "terms_set.xml"
!    call polynomial_coeff_writeXML(coefficients,my_newncoeff,filename=filename)
!  end if
  
 if(need_verbose)then
   write(message,'(1x,I0,2a)') ncoeff_tot,' coefficients generated ',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
 

!Final deallocation
 ABI_DEALLOCATE(symbols)
 ABI_DEALLOCATE(my_coeffindexes)
 ABI_DEALLOCATE(my_newcoeffindexes)
 ABI_DEALLOCATE(my_coefflist)
 if(allocated(coeffs_tmp))then
   ABI_DATATYPE_DEALLOCATE(coeffs_tmp)
 end if
 
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
!! index_coeff_in(power_disp-1) = list of previous coefficients computed (start with 0)
!! icoeff = current indexes of the cofficients (start we 1)
!! icoeff_tot = current number of coefficients already computed (start we 0)
!! natom = number of atoms in the unit cell
!! nstr = number of coefficient for related to the atomic displacment into list_symcoeff
!! nstr = number of coefficient for related to the strain into list_symstr
!! ncoeff_out = number of maximum coefficients
!! nrpt = number of cell
!! nsym = number of symmetries in the system
!! power_disp = initial power_disp to be computed (can be < power_disp_min,
!!              this routine will skip the firts power_disp)
!! power_disp_min = minimal power_disp to be computed
!! power_disp_max = maximum power_disp to be computed
!! symbols(natom) = array with the symbols of each atoms (Sr,O,Ti,...)
!! nbody = optional, number of body for the coefficients, for example:
!!                   0 => all the terms
!!                   1 => only (Sr_x-T_y)^power_disp and (Sr_x-T_y)^power_disp\eta^power_disp  ...
!! compute = logical, optional: TRUE if we store the coefficients
!!                              FALSE just to count the number of coefficient
!! anharmstr = logical, optional : TRUE, the anharmonic strain are computed
!!                                   FALSE, (default) the anharmonic strain are not computed
!! distributed = logical, optional : True, the coefficients will be distributed on the CPU
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
&                                  nrpt,nsym,power_disp,power_disp_min,power_disp_max,symbols,nbody,&
&                                  compute,anharmstr,spcoupling,distributed)

 implicit none

!Arguments ---------------------------------------------
!scalar
 integer,intent(in) :: natom,ncoeff,power_disp,power_disp_min,power_disp_max,ncoeff_out,nsym,nrpt,nstr
 integer,intent(inout) :: icoeff,icoeff_tot
 logical,optional,intent(in) :: compute,anharmstr,spcoupling,distributed
 integer,optional,intent(in) :: nbody
!arrays
 integer,intent(in) :: cell(3,nrpt),compatibleCoeffs(ncoeff+nstr,ncoeff+nstr)
 integer,intent(in) :: list_coeff(6,ncoeff,nsym),list_str(nstr,nsym,2)
 integer,intent(in) :: index_coeff_in(power_disp-1)
 type(polynomial_coeff_type),intent(inout) :: coeffs_out(ncoeff_out)
 character(len=5),intent(in) :: symbols(natom)
!Local variables ---------------------------------------
!scalar
 integer :: ia,ib,ii,icoeff1,icoeff_tmp
 integer :: iterm,nbody_in,ncoeff_max,pa,pb
 integer :: ndisp_max,nterm_max
 real(dp):: coefficient
 logical :: need_compute,compatible,possible,need_anharmstr,need_spcoupling,need_distributed
!arrays
 integer,allocatable :: index_coeff(:)
 character(len=200):: name
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)
! *************************************************************************

!Set the inputs
 need_compute = .TRUE.
 need_anharmstr = .TRUE.
 need_spcoupling = .TRUE.
 need_distributed = .FALSE.
 nbody_in = 0 !all kind of terms
 if(present(compute)) need_compute = compute
 if(present(nbody)) nbody_in = nbody
 if(present(anharmstr)) need_anharmstr = anharmstr
 if(present(spcoupling)) need_spcoupling = spcoupling
 if(present(distributed)) need_distributed  = distributed
 if(power_disp <= power_disp_max)then

!  Initialisation of variables
   nterm_max  = nsym
   ncoeff_max = (ncoeff+nstr)
   ndisp_max = power_disp
   icoeff_tmp = 0
   ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
   ABI_ALLOCATE(terms,(nterm_max))
   ABI_ALLOCATE(index_coeff,(power_disp))

   index_coeff(1:power_disp-1) = index_coeff_in(:)

   do icoeff1=icoeff,ncoeff+nstr
!    If the distance between the 2 coefficients is superior than the cut-off,
!    we cycle
!    If the power_disp is one, we need to set icoeff to icoeff1
     if(power_disp==1) icoeff = icoeff1

     if(compatibleCoeffs(icoeff,icoeff1)==0) cycle

!    Reset the flag compatible and possible
     compatible = .TRUE.
     possible   = .TRUE.

     index_coeff(power_disp) = icoeff1
     iterm = 0
     coefficient = one

     if(power_disp >= power_disp_min) then
       call generateTermsFromList(cell,index_coeff,list_coeff,list_str,ncoeff,&
&                                 ndisp_max,nrpt,nstr,nsym,iterm,terms)

       if(iterm > 0)then
!        Do some checks
!        -------------
!        1-Check if the coefficient is full anharmonic strain and if we need to compute it
         if(terms(1)%ndisp == 0)then
           compatible = (need_anharmstr .or. need_spcoupling)
           possible = need_anharmstr
         end if
!        1-Check if the coefficient is strain-coupling and if we need to compute it
         if(terms(1)%nstrain > 0.and.terms(1)%ndisp > 0)then
           possible   = need_spcoupling
           compatible = need_spcoupling
         end if

!        ------------
!        2-Check if this terms is compatible with nbody
         if(nbody_in > 0)then
           pa = 1 ; pb = 1
           ia = 0 ; ib = 0
!          Count the number of terms and the power_disp
           do ii=1,terms(1)%ndisp
             if(terms(1)%nstrain > 0) then
               pb = pb*terms(1)%power_disp(ii)
               ib = ib + 1
             else
               pa = pa*terms(1)%power_disp(ii)
               ia = ia + 1
             end if
           end do
           if(ia <= nbody_in)then
             if(ia==nbody_in.and.abs(mod(pa,2)) < tol16)then
               if(ib==0)then
                 compatible = .FALSE.
                 possible   = .TRUE.
               else if (ib==nbody_in.and.abs(mod(pb,2)) < tol16) then
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

     end if!end if power_disp < power_disp_min

     if(compatible)then
       call computeNorder(cell,coeffs_out,compatibleCoeffs,list_coeff,list_str,index_coeff,&
&                         icoeff1,icoeff_tot,natom,ncoeff,nstr,ncoeff_out,nrpt,nsym,power_disp+1,&
&                         power_disp_min,power_disp_max,symbols,nbody=nbody_in,compute=need_compute,&
&                         anharmstr=need_anharmstr,spcoupling=need_spcoupling)
     end if
   end do

   ABI_DEALLOCATE(terms)
   ABI_DEALLOCATE(index_coeff)

!  Transfer in the final array
   icoeff1 = 0
   do icoeff_tmp=1,ncoeff_max
     if (abs(coeffs_tmp(icoeff_tmp)%coefficient) > tol16)then
!      Increase icoeff and fill the coeffs_out array
       icoeff_tot = icoeff_tot + 1
       if(need_compute)then
         name = ''
!        Get the name of this coefficient
         call polynomial_coeff_getName(name,coeffs_tmp(icoeff_tmp),symbols,recompute=.TRUE.)
         call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                                   coeffs_out(icoeff_tot),coeffs_tmp(icoeff_tmp)%terms,&
&                                   name=name)
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


!!****f* m_polynomial_coeff/computeCombinationFromList
!! NAME
!! computeCombinationFromList
!!
!! FUNCTION
!! Recursive routine to compute the order N of a all the possible coefficient
!! from the list list_symcoeff and list_symstr.
!!
!! INPUTS
!! cell(3,nrpt) = indexes of the cells into the supercell (-1 -1 -1, 0 0 0 ...)
!! compatibleCoeffs(ncoeff+nstr,ncoeff+nstr) = array with the list of compatible coefficients 0 or 1
!! list_coeff(6,ncoeff_sym,nsym) = array with the list of the coefficients,
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!! list_str(nstr_sym,nsym) = array with the list of the strain  and the symmetrics
!! index_coeff_in(power_disp-1) = list of previous coefficients computed (start with 0)
!! icoeff = current indexes of the combination (start we 1)
!! max_power_strain = maximum order of the strain of the strain phonon coupling
!! nmodel_tot = current number of combination already computed (start we 0)
!! natom = number of atoms in the unit cell
!! ncoeff = number of coefficient for related to the atomic displacment into list_symcoeff
!! nstr = number of coefficient for related to the strain into list_symstr
!! nmodel = number of maximum models
!! nrpt = number of cell
!! nsym = number of symmetries in the system
!!     For example, the sum of all the term like (Sr_y-O_y)^odd, are 0 by symetrie in cubic system.
!!     Here, we build a list with: 0 this term is not allowed for odd
!!                                 1 this term is allowed for odd
!! power_disp = initial power_disp to be computed (can be < power_disp_min,
!!              this routine will skip the firts power_disp)
!! power_disp_min = minimal power_disp to be computed
!! power_disp_max = maximum power_disp to be computed
!! symbols(natom) = array with the symbols of each atoms (Sr,O,Ti,...)
!! nbody = optional, number of body for the coefficients, for example:
!!                   0 => all the terms
!!                   1 => only (Sr_x-T_y)^power_disp and (Sr_x-T_y)^power_disp\eta^power_disp  ...
!! compute = logical, optional: TRUE if we store the coefficients
!!                              FALSE just to count the number of coefficient
!! anharmstr = logical, optional : TRUE, the anharmonic strain are computed
!!                                   FALSE, (default) the anharmonic strain are not computed
!! distributed = logical, optional : True, the coefficients will be distributed on the CPU
!! only_odd_power = logical, optional : if TRUE return only odd power
!! only_even_power= logical, optional : if TRUe return only even power
!!
!! OUTPUT
!! icoeff = current indexes of the cofficients (start we 1)
!! nmodel_tot = current number of coefficients already computed (start we 0)
!! list_combination = list of the possible combination of coefficients
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

recursive subroutine computeCombinationFromList(cell,compatibleCoeffs,list_coeff,list_str,&
&                                  index_coeff_in,list_combination,icoeff,max_power_strain,nmodel_tot,&
&                                  natom,ncoeff,ncoeff_sym,nstr,nmodel,nrpt,nsym,power_disp,power_disp_min,&
&                                  power_disp_max,symbols,nbody,only_odd_power,only_even_power,&
&                                  compute,anharmstr,spcoupling)

 implicit none

!Arguments ---------------------------------------------
!scalar
 integer,intent(in) :: natom,ncoeff,ncoeff_sym,power_disp,power_disp_min,power_disp_max
 integer,intent(in) :: max_power_strain,nmodel,nsym,nrpt,nstr
 integer,intent(inout) :: icoeff,nmodel_tot
 logical,optional,intent(in) :: compute,anharmstr,spcoupling
 integer,optional,intent(in) :: nbody
 logical,optional,intent(in) :: only_odd_power,only_even_power
!arrays
 integer,intent(in) :: cell(3,nrpt),compatibleCoeffs(ncoeff+nstr,ncoeff+nstr)
 integer,intent(in) :: list_coeff(6,ncoeff_sym,nsym),list_str(nstr,nsym,2)
 integer,intent(in) :: index_coeff_in(power_disp-1)
 integer,intent(out) :: list_combination(power_disp_max,nmodel)
 character(len=5),intent(in) :: symbols(natom)
!Local variables ---------------------------------------
!scalar
 integer :: icoeff1,icoeff2,nbody_in,ii,jj,nbody_count,nmodel_tot_test
 integer :: isym_in_test,idisp_in_test,ndisp_test,ndisp_out,nstrain
 logical :: need_compute,compatible,possible,need_anharmstr,need_spcoupling
 logical :: need_only_odd_power,need_only_even_power,compute_test
!arrays
 integer :: powers(power_disp)
 integer,allocatable :: index_coeff(:),dummylist(:)
! *************************************************************************

!Set the inputs
 need_compute = .TRUE.
 need_anharmstr = .TRUE.
 need_spcoupling = .TRUE.
 need_only_odd_power = .FALSE.
 need_only_even_power = .FALSE.
 nbody_in = 0 !all kind of terms
 if(present(compute)) need_compute = compute
 if(present(nbody)) nbody_in = nbody
 if(present(anharmstr)) need_anharmstr = anharmstr
 if(present(spcoupling)) need_spcoupling = spcoupling
 if(present(only_odd_power)) need_only_odd_power = only_odd_power
 if(present(only_even_power)) need_only_even_power = only_even_power

 if(power_disp <= power_disp_max)then

!  Initialisation of variables
   ABI_ALLOCATE(index_coeff,(power_disp))
   index_coeff(1:power_disp-1) = index_coeff_in(:)
!   write(std_out,*) "DEBUG index_coff in recursive CCL: ", index_coeff
!  Loop over ncoeff+nstr
   do icoeff1=icoeff,ncoeff+nstr

!    Reset the flag compatible and possible
     compatible = .TRUE.
     possible   = .TRUE.

!    If the power_disp is one, we need to set icoeff to icoeff1
     if(power_disp==1) then
       icoeff = icoeff1
       if(compatibleCoeffs(icoeff,icoeff1)==0)then
         compatible = .FALSE.
       end if
     end if
!    If the distance between the 2 coefficients is superior than the cut-off, we cycle.
    do icoeff2=1,power_disp-1
!      write(std_out,*) "icoeff1: ", icoeff1 
!      write(std_out,*) "icoeff2: ", icoeff2, "index_icoeff2: ", index_coeff(icoeff2)
      if(compatibleCoeffs(index_coeff(icoeff2),icoeff1)==0)then
        compatible = .FALSE.
      end if
    end do

     if (.not.compatible) cycle !The distance is not compatible

!    Set the index of the new coeff in the list
     index_coeff(power_disp) = icoeff1
!    Do some checks
!    -------------
!    1-Check if the coefficient is full anharmonic strain and if we need to compute it
     if(all(index_coeff > ncoeff))then
       compatible = (need_anharmstr .or. need_spcoupling)
       possible = need_anharmstr
     end if
!    2-Check if the coefficient is strain-coupling and if we need to compute it
     if(any(index_coeff < ncoeff) .and. any(index_coeff > ncoeff))then
       possible   = need_spcoupling
       compatible = need_spcoupling
       if(count(index_coeff > ncoeff) > max_power_strain)then
         possible = .false.
         compatible = .false.
       end if
     end if
!    3-Count number of Strain and number of displacements for compute symmetric terms
     nstrain = 0 
     ndisp_out = 0 
     do ii=1,power_disp
        if(index_coeff(ii) > 0 .and. index_coeff(ii) <= ncoeff)then 
           ndisp_out = ndisp_out + 1
        else 
           nstrain = nstrain +1 
           index_coeff(ii) = index_coeff(ii) - ncoeff + ncoeff_sym
        end if
     end do 

     if(power_disp >= power_disp_min) then

!      count the number of body
       powers(:) = 1
       do ii=1,power_disp
         do jj=ii+1,power_disp
           if (powers(jj) == 0) cycle
           if(index_coeff(ii)==index_coeff(jj))then
             powers(ii) = powers(ii) + 1
             powers(jj) = 0
           end if
         end do
       end do
       nbody_count = count(powers /= 0) 
!       write(std_out,*) "powers: ", powers
!       write(std_out,*) "nbody_count: ", nbody_count
!      check the only_odd and only_even flags
       if(any(mod(powers(1:power_disp),2) /=0) .and. need_only_even_power) then
         possible = .false.
       end if
       if(any(mod(powers(1:power_disp),2) ==0) .and. need_only_odd_power)then
         possible = .false.
       end if
!      Check the nbody flag
       if(nbody_in /= 0)then
         if(power_disp-count(powers==0) > nbody_in) then
           possible = .false.
           compatible = .false.
         end if
       end if

       if(possible) then
!        increase coefficients and set it
!         nmodel_tot = nmodel_tot + 1
!         if(need_compute)then
!           list_combination(1:power_disp,nmodel_tot) = index_coeff
!         end if
         !nmodel_tot_test = 0 
!         !Start from second symmetry in Symmetric Combinations
!         isym_in_test = 2 
!         idisp_in_test = power_disp 
!         ndisp_test = power_disp
!         index_coeff_tmp = index_coeff
         !Count anharmonic strain terms
         if(ndisp_out == 0 .and. nstrain > 0)then 
            nmodel_tot = nmodel_tot + 1 
            if(need_compute)list_combination(1:power_disp,nmodel_tot) = index_coeff
         else !Else construct symmetric terms of atomic displacement (pure disp or disp/strain)
            ABI_ALLOCATE(dummylist,(0))
            write(std_out,*) "DEBUG index_coeff: ", index_coeff
            write(std_out,*) "DEBUG ndisp_out: ", ndisp_out
            write(std_out,*) "DEBUG nstrain: ", nstrain
            compute_test = need_compute
            call computeSymmetricCombinations(nmodel_tot,list_combination,list_coeff,1,1,ndisp_out,nsym,&
                                              dummylist,index_coeff,power_disp_max,nmodel,ncoeff_sym,&
&                                             nstrain,compute_test) 
            write(std_out,*) "DEBUG nmodel_tot: ", nmodel_tot
            !nmodel_tot = nmodel_tot + nmodel_tot_test 
            ABI_DEALLOCATE(dummylist)
         end if !ndisp_out == 0 .and.n nstrain >0 
       end if!possible
     end if!end if power_disp < power_disp_min

     !Change back to irreducible terms ncoeff_limit
     do ii=1,power_disp
        if(index_coeff(ii) > ncoeff_sym)then 
           index_coeff(ii) = index_coeff(ii) + ncoeff - ncoeff_sym
        end if
     end do 

!    If the model is still compatbile with the input flags, we continue.
     if(compatible)then
       call computeCombinationFromList(cell,compatibleCoeffs,list_coeff,list_str,&
&                                     index_coeff,list_combination,icoeff1,max_power_strain,&
&                                     nmodel_tot,natom,ncoeff,ncoeff_sym,nstr,nmodel,nrpt,nsym,power_disp+1,&
&                                     power_disp_min,power_disp_max,symbols,nbody=nbody_in,&
&                                     compute=need_compute,anharmstr=need_anharmstr,&
&                                     spcoupling=need_spcoupling,only_odd_power=need_only_odd_power,&
&                                     only_even_power=need_only_even_power)
     end if
   end do
   ABI_DEALLOCATE(index_coeff)
 end if

end subroutine computeCombinationFromList
!!***

!!****f* m_polynomial_coeff/computeSymmetricCombinations
!! NAME
!! computeSymmetricCombinations
!!
!! FUNCTION
!! Compute all symmetric possible combinations for a given displacement power and combinations of terms 
!!
!! INPUTS
!! list_symcoeff: List irreducible coefficients and their symmetric counterparts  
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!!
!! nsym = number of symmetries 
!! isym = number of current symmetry to iterate from (recursevily changed)
!! idisp = number of current displacement on which we iterate (recursevily changed)  
!! ndisp = current power of displacements (order of the term) 
!! compute = if true computes the terms and stores them if falls counts the possibilites
!! 
!!
!! OUTPUT
!! isym = current index of the symmetry
!! idisp = current index of the displacement
!! list_combination = list of the possible combination of coefficients
!! nmodel_tot = total combinations of terms possible 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

recursive subroutine computeSymmetricCombinations(ncombi,list_combination,list_symcoeff,isym_in,&
&                                               idisp_in,ndisp,nsym,index_isym_in,index_coeff_in,&
&                                               ndisp_max,ncombinations,ncoeff,nstrain,compute)
 
 implicit none

!Arguments ------------------------------------
integer,intent(inout) :: ncombi
integer,intent(in)    :: ndisp,nsym,ndisp_max,ncombinations,ncoeff,nstrain
integer,intent(in)    :: isym_in,idisp_in
logical,intent(in)    :: compute 
!scalar
!arrays
integer,intent(inout) :: list_combination(ndisp_max,ncombinations)
integer,intent(in)    :: list_symcoeff(6,ncoeff,nsym),index_coeff_in(ndisp+nstrain)
integer,intent(in)    :: index_isym_in(idisp_in-1)
!Local variables-------------------------------
!scalar
integer :: isym,idisp
logical :: need_compute 
!arrays
integer :: index_coeff_tmp(ndisp)
integer,allocatable :: index_isym(:)
!Source
! *************************************************************************

!Only start the function if start-symmetry is smaller than maximum symmetry
!and start displacement is smaller than maximum displacement
!otherwise pass through
if(isym_in <= nsym .and. idisp_in <= ndisp)then 
  need_compute = compute  
  !Allocate index_isym 
  !Copy index of symmetries done before. Symmetry of first term is always one 
  ABI_ALLOCATE(index_isym,(idisp_in))
  index_isym(1) = 1
  if(idisp_in > 2)then  
    index_isym(2:idisp_in) = index_isym_in 
  endif

  do isym = isym_in,nsym 
     ! Put the current symmetry to the current displacement 
     index_isym(idisp_in) = isym
     if(idisp_in == ndisp)then  
       if(ndisp == 1 .and. isym ==1 .or. ndisp > 1)then !If term is just one body just store one time       
         ncombi = ncombi + 1  
         if(need_compute)then 
           !loop over displacements in term 
           do idisp=1,ndisp
              index_coeff_tmp(idisp) = list_symcoeff(6,index_coeff_in(idisp),index_isym(idisp))
           end do !idisp=1,ndisp
           if(any(index_coeff_tmp == 0))then ! If symmetry doesn't point to another term write zeros to filter after
              list_combination(:,ncombi) = 0 
           else
             list_combination(:ndisp,ncombi) = index_coeff_tmp 
             !write(std_out,*) "DEBUG, index_coeff_in,: ", index_coeff_in,"ndisp: ", ndisp
             !write(std_out,*) "DEBUG, index_coeff_tmp,: ", index_coeff_in,"ndisp: ", ndisp            
             if(nstrain /= 0)then !If SP coupling copy strain index
                list_combination(ndisp+1:ndisp+nstrain,ncombi) = index_coeff_in(ndisp+1:ndisp+nstrain)
             end if
           end if
         end if ! need compute
       end if !ndisp == 1 .and isym == 1
     end if !(idisp_in == ndisp)
     
     call computeSymmetricCombinations(ncombi,list_combination,list_symcoeff,isym,idisp_in+1,ndisp,& 
&                                      nsym,index_isym(2:),index_coeff_in,ndisp_max,ncombinations,&
&                                      ncoeff,nstrain,compute)
     
  enddo !isym=isym_in,nsym
  ABI_DEALLOCATE(index_isym)
endif!(isym <= nsym)


end subroutine computeSymmetricCombinations

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

function getCoeffFromList(list_coeff,ia,ib,irpt,mu,ncoeff) result(coeff)

 implicit none

!Arguments ------------------------------------
!scalar
 integer,intent(in) :: ia,ib,irpt,mu,ncoeff
 integer :: coeff
!arrays
 integer,intent(in) :: list_coeff(6,ncoeff)
!Local variables-------------------------------
!scalar
 integer :: icoeff
!arrays

! *************************************************************************
 coeff = 0
 !write(std_out,*) "DEBUG shape inside FromList: ", shape(list_coeff)
 do icoeff = 1,ncoeff
   if(mu==list_coeff(1,icoeff).and.&
&     ia==list_coeff(2,icoeff).and.&
&     ib==list_coeff(3,icoeff).and.&
&     irpt==list_coeff(4,icoeff))then!.and.&
!&     abs(weight-list_coeff(5,icoeff)) < tol16) then
     coeff = icoeff
     exit
   end if
 end do

end function getCoeffFromList
!!***


!!****f* m_polynomial_coeff/generateTermsFromList
!!
!! NAME
!! generateTermsFromList
!!
!! FUNCTION
!! Compute for a given list of index the correspondig set of terms
!!
!! INPUTS
!! cell(3,nrpt) = indexes of the cells into the supercell (-1 -1 -1, 0 0 0 ...)
!! index_coeff_in(ndisp) = list of coefficients to be computed
!! list_symcoeff(6,ncoeff,nsym) = array with the list of the coefficients,
!!                                    for each coefficients (ncoeff_sym), we store the symmetrics(nsym)
!!                                    the 6th first dimensions are :
!!                                       1 = direction of the IFC
!!                                       2 = index of the atom number 1 (1=>natom)
!!                                       3 = index of the atom number 2 (1=>natom)
!!                                       4 = indexes of the cell of the second atom
!!                                           (the atom number 1 is always in the cell 0 0 0)
!!                                       5 = weight of the term (-1 or 1)
!!                                       6 = indexes of the symmetric
!! list_symstr(nstr,nsym) = array with the list of the strain  and the symmetrics
!! ncoeff = number of maximum coefficients in the list_symcoeff
!! ndisp = number of maximum diplacement (phonon + strain)
!! nrpt = number of cell
!! nsym = number of symmetries in the system
!!
!! OUTPUT
!! terms<(type(polynomial_term_type)>(nterm)  = list of terms
!! nterm = number of ouput terms
!!
!! PARENTS
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine generateTermsFromList(cell,index_coeff,list_coeff,list_str,ncoeff,ndisp_max,&
&                                nrpt,nstr,nsym,nterm,terms)

 implicit none

!Arguments ------------------------------------
!scalar
 integer,intent(in) :: ndisp_max,ncoeff,nrpt,nstr,nsym
 integer,intent(out):: nterm
!arrays
 integer,intent(in) :: index_coeff(ndisp_max)
 integer,intent(in) :: cell(3,nrpt),list_coeff(6,ncoeff,nsym)
 integer,intent(in) :: list_str(nstr,nsym,2)
 type(polynomial_term_type),intent(out) :: terms(nsym)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff_str,idisp,irpt
 integer :: isym,ndisp,nstrain,mu
 real(dp):: weight
!arrays
 integer :: atindx(2,ndisp_max),cells(3,2,ndisp_max),dir_int(ndisp_max),strain(ndisp_max)
 integer :: power_disps(ndisp_max),power_strain(ndisp_max)

! *************************************************************************
 nterm = 0
!Loop over symetries
 do isym=1,nsym
!Treat this coeff
   weight = 1
   ndisp = 0
   nstrain = 0
   do idisp=1,ndisp_max
!    Get index of this displacement term
!    Check if the index is not zero
     if(index_coeff(idisp)==0)cycle
     if(index_coeff(idisp)<=ncoeff)then
       ndisp = ndisp + 1
       mu   = list_coeff(1,index_coeff(idisp),isym)
       ia   = list_coeff(2,index_coeff(idisp),isym)
       ib   = list_coeff(3,index_coeff(idisp),isym)
       irpt = list_coeff(4,index_coeff(idisp),isym)
       weight = weight*list_coeff(5,index_coeff(idisp),isym)
!      Fill First term arrays
       atindx(1,idisp) = ia; atindx(2,idisp) = ib;
       dir_int(idisp) = mu
       power_disps(idisp)   = 1
       cells(:,1,idisp) = (/0,0,0/)
       cells(:,2,idisp) = cell(:,irpt)
     else
       nstrain = nstrain + 1
       icoeff_str = index_coeff(idisp)-ncoeff
       strain(nstrain) = list_str(icoeff_str,isym,1)
       power_strain(nstrain)  = 1
       weight = weight*list_str(icoeff_str,isym,2)
     end if
   end do
   nterm = nterm + 1
   call polynomial_term_init(atindx,cells,dir_int,ndisp,nstrain,terms(nterm),power_disps,&
&                            power_strain,strain,weight,check=.true.)
 end do!end do sym


end subroutine generateTermsFromList
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
!! natom = number of atoms in the unit cell
!! nrpt = number of cell
!! nsym = number of symmetries in the system
!! symbols(natom) = array with the symbols of each atoms (Sr,O,Ti,...)
!! comm = MPI communicator
!!
!! OUTPUT
!! polynomial_coeff<(type(polynomial_coeff_type)>(ncoeff_out) = array of datatype with
!!                                                              the polynomial_coeff
!! ncoeff_out = number of coefficients
!!
!! PARENTS
!!
!! CHILDREN
!!      polynomial_coeff_free,polynomial_coeff_getname,polynomial_coeff_init
!!      polynomial_term_free,polynomial_term_init,wrtout
!!
!! SOURCE

subroutine polynomial_coeff_getOrder1(cell,coeffs_out,list_symcoeff,&
&                                     natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                     symbols)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff,nsym,nrpt
 integer,intent(out) :: ncoeff_out
!arrays
 integer,intent(in) :: cell(3,nrpt)
 integer,intent(in) :: list_symcoeff(6,ncoeff,nsym)
 character(len=5),intent(in) :: symbols(natom)
 type(polynomial_coeff_type),allocatable,intent(inout) :: coeffs_out(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff_tmp,irpt,irpt_ref
 integer :: isym,iterm,mu,ncoeff_max,ndisp,nstrain,nterm_max
 real(dp):: coefficient,weight
!arrays
 integer,allocatable :: atindx(:,:),cells(:,:,:),dir_int(:)
 integer,allocatable :: power_disps(:),power_strain(:),strain(:)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=200):: name
 character(len=500) :: message
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)
!TEST_AM
 character(len=fnlen) :: filename
!TEST_AM
! *************************************************************************

!Initialisation of variables
 nterm_max  = nsym
 ncoeff_max = ncoeff
 ndisp = 1
 nstrain = 0
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))


 icoeff_tmp = 0
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(power_disps,(ndisp))
 ABI_ALLOCATE(power_strain,(nstrain))
 ABI_ALLOCATE(strain,(nstrain))

!Found the ref cell
 irpt_ref = 1
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
   iterm = 0
   coefficient = one
   do isym=1,nsym
     ndisp   = 1
     nstrain = 0
     mu   = list_symcoeff(1,icoeff,isym)
     ia   = list_symcoeff(2,icoeff,isym)
     ib   = list_symcoeff(3,icoeff,isym)
     irpt = list_symcoeff(4,icoeff,isym)
     weight = list_symcoeff(5,icoeff,isym)
!    Fill First term arrays
     atindx(1,1) = ia; atindx(2,1) = ib;
     dir_int(1) = mu
     power_disps(1)   = 1
     cells(:,1,1) = (/0,0,0/)
     cells(:,2,1) = cell(:,irpt)
     iterm = iterm + 1
     call polynomial_term_init(atindx,cells,dir_int,ndisp,nstrain,terms(iterm),&
&                              power_disps,power_strain,strain,weight,check=.true.)
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
 ABI_DEALLOCATE(power_disps)
 ABI_DEALLOCATE(power_strain)
 ABI_DEALLOCATE(strain)

!Count the number of terms
 ncoeff_out = 0
 do icoeff_tmp=1,ncoeff_max
   if (abs(coeffs_tmp(icoeff_tmp)%coefficient) > tol16)then
     ncoeff_out = ncoeff_out + 1
   end if
 end do

!Transfer in the final array
 ABI_ALLOCATE(coeffs_out,(ncoeff_out))
 icoeff = 0
 do icoeff_tmp=1,ncoeff_max
   if (abs(coeffs_tmp(icoeff_tmp)%coefficient) > tol16)then
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,coeffs_tmp(icoeff_tmp),symbols,recompute=.TRUE.)
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
       if(any(coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(:,2,1)/=0))then
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
 filename = "terms_1st_order.xml"
 call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
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

!!****f* m_polynomial_coeff/polynomial_coeff_getEvenAnhaStrain
!!
!! NAME
!! polynomial_coeff_getEvenAnhaStrain
!!
!! FUNCTION
!! Get even anharmonic strain terms in defined range of order
!!
!! INPUTS
!!
!!
!! OUTPUT
!! polynomial_coeff<(type(polynomial_coeff_type)>(ncoeff_out) = array of datatype with
!!                                                              the polynomial_coeff
!! ncoeff_out = number of coefficients
!!
!! PARENTS
!!
!! CHILDREN
!! polynomial_coeff_getNorder 
!!
!! SOURCE
subroutine polynomial_coeff_getEvenAnhaStrain(strain_terms,crystal,irred_ncoeff,power_strain,comm)

 implicit none

!Arguments ------------------------------------
type(polynomial_coeff_type),allocatable,intent(inout) :: strain_terms(:)
type(crystal_t), intent(inout) :: crystal
integer,intent(out) :: irred_ncoeff 
integer,intent(in) :: power_strain(2)
integer,intent(in) :: comm
!scalars
!arrays
!Local variables-------------------------------
real(dp) :: cutoff 
integer :: ncoeff
integer :: power_strph
integer :: option 
integer :: icoeff1,icoeff2,start
integer :: ncoeff_out 
character(len=fnlen) :: fname
logical,allocatable :: same(:)
type(polynomial_coeff_type),allocatable :: strain_terms_tmp(:)
!scalar
!arrays
integer :: sc_size(3) 
! *************************************************************************

!Initialize Variables for call to polynomial_coeff_getNorder 
cutoff = zero
power_strph = zero
option = 0 
sc_size = (/1,1,1/)


call polynomial_coeff_getNorder(strain_terms_tmp,crystal,cutoff,ncoeff,ncoeff_out,power_strain,& 
&                               power_strph,option,sc_size,comm,anharmstr=.true.,spcoupling=.false.,&
&                               only_odd_power=.false.,only_even_power=.true.,verbose=.false.) 


!TODO Probably put in one routine
!Get irreducible strain
ABI_ALLOCATE(same,(ncoeff_out))
same = .false.
do icoeff1 =1,ncoeff_out
        start = icoeff1 + 1
        !write(*,*) "name coeff_",icoeff1,": ", strain_terms_tmp(icoeff1)%name
        do icoeff2=start,ncoeff_out
                if(.not.same(icoeff2)) same(icoeff2) = coeffs_compare(strain_terms_tmp(icoeff1),strain_terms_tmp(icoeff2))
        enddo 
        !write(*,*) "same(",icoeff1,"): ", same(icoeff1)
enddo 

irred_ncoeff = ncoeff_out - count(same)
ABI_DATATYPE_ALLOCATE(strain_terms,(irred_ncoeff))

!Transfer irreducible strain to output array
icoeff2=0
icoeff1=0
do icoeff1=1,ncoeff_out
        if(.not.same(icoeff1))then 
                icoeff2=icoeff2 + 1
                strain_terms(icoeff2) = strain_terms_tmp(icoeff1)
        endif 
enddo

!TEST MS write coefficients to xml to check 
!fname="EvenAnhaStrain.xml"
!call polynomial_coeff_writeXML(strain_terms,irred_ncoeff,filename=fname,newfile=.true.)
!write(*,*) "Strain terms to xml done"


!Deallocateion
ABI_DATATYPE_DEALLOCATE(strain_terms_tmp)
ABI_DEALLOCATE(same)

end subroutine polynomial_coeff_getEvenAnhaStrain 
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
 implicit none

!Arguments ------------------------------------
  type(polynomial_coeff_type), intent(in) :: c1,c2
  logical :: res
!local
!variable
  integer :: iterm1,iterm2
!array
  integer :: blkval(2,c1%nterm)
! *************************************************************************
  res = .false.
  blkval = 0
  do iterm1=1,c1%nterm
    if(blkval(1,iterm1)==1)cycle!already found
    do iterm2=1,c2%nterm
      if(blkval(2,iterm2)==1)cycle!already found
      if(c1%terms(iterm1)==c2%terms(iterm2)) then
        blkval(1,iterm1) = 1
        blkval(2,iterm2) = 1
      end if
    end do
  end do
  if(.not.any(blkval(:,:)==0))res = .true.

end function coeffs_compare
!!***


!!****f* m_polynomial_coeff/coeffs_list_conc
!! NAME
!! coeff_list_conc
!!
!! FUNCTION
!! 
!! Concatenate list1 and list2 of type polynomial_coeff and store it in list_out
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function coeffs_list_conc(coeff_list1,coeff_list2) result (coeff_list_out)
!Arguments ------------------------------------
 implicit none

!Arguments ------------------------------------
  type(polynomial_coeff_type), intent(in) :: coeff_list1(:),coeff_list2(:)
  type(polynomial_coeff_type),allocatable :: coeff_list_out(:)
!local
!variable
  integer :: ncoeff1,ncoeff2,ncoeff_out 
!array
! *************************************************************************

!Free output 
 if(allocated(coeff_list_out))then
   ABI_DATATYPE_DEALLOCATE(coeff_list_out) 
 endif

!Get sizes of coeff_list1/2
 ncoeff1 = size(coeff_list1) 
 ncoeff2 = size(coeff_list2)
 ncoeff_out = ncoeff1 + ncoeff2

 !Allocate output 
 ABI_DATATYPE_ALLOCATE(coeff_list_out,(ncoeff_out)) 

 !Copy input list into output lists 
 coeff_list_out(1:ncoeff1) = coeff_list1 
 coeff_list_out(ncoeff1+1:) = coeff_list2   
 
end function coeffs_list_conc
!!***

end module m_polynomial_coeff
!!***
