!!****m* ABINIT/m_paw_energies
!! NAME
!!  m_paw_energies
!!
!! FUNCTION
!!  This module contains the definition of the paw_energies_type structured datatype,
!!  as well as related functions and methods.
!!  paw_energies_type variables define several contributions to PAW on-site ENERGIES
!!
!! COPYRIGHT
!! Copyright (C) 2013-2025 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_energies

 USE_DEFS
 USE_MSG_HANDLING

 implicit none

 private

!public parameter
 integer, public, parameter :: n_paw_energies=6

!!***

!----------------------------------------------------------------------

!!****t* m_paw_denpot/paw_energies_type
!! NAME
!! paw_energies_type
!!
!! FUNCTION
!! This structured datatype contains all parts of the PAW contribution to energy
!!
!! SOURCE

 type, public :: paw_energies_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  real(dp) :: epaw
   ! total on-site PAW energy (direct scheme)

  real(dp) :: epaw_dc
   ! total on-site PAW energy (double counting scheme)

  real(dp) :: epaw_core
   ! core contribution to PAW energy (direct scheme)

  real(dp) :: epaw_core_dc
   ! core contribution to PAW energy (double counting scheme)

  real(dp) :: epaw_xc
   ! exchange-correlation on-site contribution to PAW energy

  real(dp) :: entropy_paw
   ! on-site PAW contribution to total entropy

 end type paw_energies_type

!public procedures
 public :: paw_energies_setzero  ! Set all energies in a paw_energies datastructure to zero
 public :: paw_energies_copy     ! Copy a paw_energies_type object into another 
 public :: paw_energies_to_array ! Transfer a paw_energies datastructure into/from a single array 
 public :: paw_energies_print    ! Printout of the object
!!***

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_energies/paw_energies_setzero
!! NAME
!!  paw_energies_setzero
!!
!! FUNCTION
!!  Set all energy contributions to zero in a paw_energies structure
!!
!! SIDE EFFECTS
!!  Paw_energies<type(paw_energies_type)>=content set to zero
!!
!! SOURCE

subroutine paw_energies_setzero(Paw_energies)

!Arguments ------------------------------------
!arrays
 type(Paw_energies_type),intent(inout) :: Paw_energies

!Local variables-------------------------------

! *************************************************************************

 !@Paw_energies_type

 ! === Reset all energies ===

 Paw_energies%epaw         = zero
 Paw_energies%epaw_dc      = zero
 Paw_energies%epaw_core    = zero
 Paw_energies%epaw_core_dc = zero
 Paw_energies%epaw_xc      = zero
 Paw_energies%entropy_paw  = zero

end subroutine paw_energies_setzero
!!***

!----------------------------------------------------------------------

!!****f* m_paw_energies/paw_energies_copy
!!
!! NAME
!! paw_energies_copy
!!
!! FUNCTION
!! Copy a paw_energies_type object into another
!!
!! INPUTS
!!   paw_energies_in <type(paw_energies_type)>=input values (to copy)
!!
!! OUTPUT
!!   paw_energies_out <type(paw_energies_type)>=output values
!!
!! SOURCE

 subroutine paw_energies_copy(paw_energies_in, paw_energies_out)

!Arguments ------------------------------------
!scalars
 type(paw_energies_type),intent(in)  :: paw_energies_in
 type(paw_energies_type),intent(out) :: paw_energies_out

!*************************************************************************

!@paw_energies_type

 paw_energies_out%epaw         = paw_energies_in%epaw
 paw_energies_out%epaw_dc      = paw_energies_in%epaw_dc
 paw_energies_out%epaw_core    = paw_energies_in%epaw_core
 paw_energies_out%epaw_core_dc = paw_energies_in%epaw_core_dc
 paw_energies_out%epaw_xc      = paw_energies_in%epaw_xc
 paw_energies_out%entropy_paw  = paw_energies_in%entropy_paw

end subroutine paw_energies_copy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_energies/paw_energies_to_array
!!
!! NAME
!! paw_energies_to_array
!!
!! FUNCTION
!! Transfer a paw_energies datastructure into a single array or
!! transfer an array into a paw_energies datastructure
!!
!! INPUTS
!!   option= 1: copy paw_energies datastructure into an array
!!   option=-1: copy an array into a paw_energies datastructure
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   paw_energies <type(paw_energies_type)>=energies stored in a datastructure
!!   paw_energies_array=energies stored in a single array
!!
!! SOURCE

 subroutine paw_energies_to_array(paw_energies,paw_energies_array,option)

!Arguments ------------------------------------
!scalars
 type(paw_energies_type),intent(inout)  :: paw_energies
 integer,intent(in) :: option
!arrays
 real(dp),intent(inout) :: paw_energies_array(:)

!Local variables-------------------------------
!scalars
 character(len=100) :: msg

!*************************************************************************

!@paw_energies_type

 if (n_paw_energies<6) then
   msg='error on number of paw_energies!'
   LIBPAW_BUG(msg)
 end if
 if (size(paw_energies_array)<n_paw_energies) then
   msg='error on paw_energies_array size!'
   LIBPAW_BUG(msg)
 end if
 
 if (option==1) then
   paw_energies_array(1)=paw_energies%epaw
   paw_energies_array(2)=paw_energies%epaw_dc
   paw_energies_array(3)=paw_energies%epaw_core
   paw_energies_array(4)=paw_energies%epaw_core_dc
   paw_energies_array(5)=paw_energies%epaw_xc
   paw_energies_array(6)=paw_energies%entropy_paw
 end if

 if (option==-1) then
   paw_energies%epaw         = paw_energies_array(1)
   paw_energies%epaw_dc      = paw_energies_array(2)
   paw_energies%epaw_core    = paw_energies_array(3)
   paw_energies%epaw_core_dc = paw_energies_array(4)
   paw_energies%epaw_xc      = paw_energies_array(5)
   paw_energies%entropy_paw  = paw_energies_array(6)
 end if

end subroutine paw_energies_to_array
!!***

!----------------------------------------------------------------------

!!****f* m_paw_energies/paw_energies_print
!! NAME
!! paw_energies_print
!!
!! FUNCTION
!!  Print out the content of a paw_energies datastructure
!!
!! INPUTS
!!  Paw_energies<paw_energies_type> = PAW energy contributions
!!
!! OUTPUT
!!  Only writing
!!
!! SOURCE

subroutine paw_energies_print(Paw_energies,unit,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
 type(Paw_energies_type), intent(in) :: Paw_energies

!Local variables-------------------------------
!scalars
 integer :: my_unt
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt   =ab_out ; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(msg,'(6a)')&
&  ' ============================== ',ch10,&
&  ' ==== Info on PAW ENERGIES ==== ',ch10,&
&  ' ============================== ',ch10
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a)')'                                 '
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a)')'  ****************************** '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,i4)')'  Total on-site PAW energy (direct) ................ ',Paw_energies%epaw
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i4)')'  Total on-site PAW energy (dble-counting) ......... ',Paw_energies%epaw_dc
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i4)')'  Core contribution to PAW energy (direct) ......... ',Paw_energies%epaw_core
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i4)')'  Core contribution to PAW energy (dble-counting) .. ',Paw_energies%epaw_core_dc
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i4)')'  XC contribution to PAW energy .................... ',Paw_energies%epaw_xc
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i4)')'  Contribution to PAW entropy ...................... ',Paw_energies%entropy_paw
 call wrtout(my_unt,msg,my_mode)

end subroutine paw_energies_print
!!***

!----------------------------------------------------------------------

END MODULE m_paw_energies
!!***
