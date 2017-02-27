!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_fc
!! NAME
!! calc_fc
!!
!! FUNCTION
!! calculation and output of Fermi-contact term at each atomic site
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (JWZ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nspden=number of spin density components
!!  ntypat=number of atom types
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!  usepaw=1 if PAW is activated
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,make_fc_paw,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine calc_fc(my_natom,natom,nspden,ntypat,pawrad,pawrhoij,pawtab,typat,usepaw,&
 &                  mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_fc'
 use interfaces_14_hidewrite
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars 
 integer,intent(in) :: my_natom,natom,nspden,ntypat,usepaw
 integer,optional,intent(in) :: comm_atom
!arrays 
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,my_comm_atom
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: message
!arrays 
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: fc(:,:)

!***********************************************************************


!Compatibility tests
 if (usepaw /= 1) then
   message = ' usepaw /= 1 but Fermi-contact calculation requires PAW '
   MSG_ERROR(message)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Initialization
 ABI_ALLOCATE(fc,(nspden,natom))

!Computation
 if (paral_atom) then
   call make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab)
 end if

!Printing
 write(message,'(a,a,a)' ) ch10,' Fermi-contact Term Calculation ',ch10
 call wrtout(ab_out,message,'COLL')

 do iatom = 1, natom
   if (nspden == 2) then
     write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC total = ',&
&     fc(1,iatom)+fc(2,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC up - down = ',&
&     fc(1,iatom)-fc(2,iatom)
     call wrtout(ab_out,message,'COLL')
   else
     write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC = ',&
&     fc(1,iatom)
     call wrtout(ab_out,message,'COLL')
   end if
 end do

 write(message,'(3a)')ch10,ch10,ch10
 call wrtout(ab_out,message,'COLL')

!Memory deallocation
 ABI_DEALLOCATE(fc)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine calc_fc
!!***
