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
   ! Exchange-correlation on-site contribution to PAW energy

  real(dp) :: spaw
   ! on-site PAW contribution to total entropy

 end type paw_energies_type

!public procedures
 public :: paw_energies_setzero ! Set all energies in a paw_energies datastructure to zero
 public :: paw_energies_print   ! Printout of the object
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
 Paw_energies%spaw         = zero

end subroutine paw_energies_setzero
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
 write(msg,'(a,i4)')'  Contribution to PAW entropy ...................... ',Paw_energies%spaw
 call wrtout(my_unt,msg,my_mode)

end subroutine paw_energies_print
!!***

!----------------------------------------------------------------------

END MODULE m_paw_energies
!!***
