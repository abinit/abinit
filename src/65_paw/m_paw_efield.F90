!!****m* ABINIT/m_paw_efield
!! NAME
!!  m_paw_efield
!!
!! FUNCTION
!!  This module contains routines related to the treatment of electric field in the PAW approach.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (FJ, PH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_efield

 use defs_basis
 use m_abicore
 use m_errors
 use m_time, only : timab
 use m_xmpi, only : xmpi_sum

 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type

 implicit none

 private

!public procedures.
 public :: pawpolev ! Compute the PAW on-site term for polarization

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/pawpolev
!! NAME
!! pawpolev
!!
!! FUNCTION
!! Compute the PAW term for polarization, named expected value term
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (FJ, PH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  ntypat = number of atom types
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  pelev(3)= electronic polarisation. expectation value term (PAW only)
!!
!! PARENTS
!!      berryphase_new
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine pawpolev(my_natom,natom,ntypat,pawrhoij,pawtab,pelev,&
&                   comm_atom) ! optional argument (parallelism)

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat
 integer,optional,intent(in) :: comm_atom
!arrays
 real(dp),intent(out) :: pelev(3)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)


!Local variables ---------------------------------------
!scalars
 integer :: iatom,idir,ierr,irhoij,ispden,itypat,jrhoij,klmn
 logical :: paral_atom
 real(dp) :: c1,ro_dlt
!arrays
 integer,dimension(3) :: idirindx = (/4,2,3/)
 real(dp) :: tsec(2)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(560,1,tsec)

 if (my_natom>0) then
   ABI_CHECK(pawrhoij(1)%qphase==1,'pawpolev not supposed to be called with qphase/=1!')
 end if

!Check for parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))

!note that when vector r is expanded in real spherical harmonics, the factor
!sqrt(four_pi/three) appears, as in the following
!x = sqrt(four_pi/three)*r*S_{1,1}
!y = sqrt(four_pi/three)*r*S_{1,-1}
!z = sqrt(four_pi/three)*r*S_{1,0}
!
!the moments pawtab()%qijl do not include such normalization factors
!see pawinit.F90 for their definition and computation

 c1=sqrt(four_pi/three)

 pelev=zero
 do idir=1,3
   do iatom=1,my_natom
     itypat=pawrhoij(iatom)%itypat
     do ispden=1,pawrhoij(iatom)%nspden
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         pelev(idir)=pelev(idir)+ro_dlt*c1*pawtab(itypat)%qijl(idirindx(idir),klmn)
         jrhoij=jrhoij+pawrhoij(iatom)%cplex_rhoij
       end do
     end do
   end do
 end do

 if (paral_atom) then
   call xmpi_sum(pelev,comm_atom,ierr)
 end if

 call timab(560,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawpolev
!!***

!----------------------------------------------------------------------

END MODULE m_paw_efield
!!***
