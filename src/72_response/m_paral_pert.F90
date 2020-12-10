!!****m* ABINIT/m_paral_pert
!! NAME
!!  m_paral_pert
!!
!! FUNCTION
!!  This module provides routines to manage the parallelisation/distribution
!!  over perturbations
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MT,FJ,MD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paral_pert

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use defs_abitypes, only : MPI_type
 use m_time,      only : timab
 use m_copy,      only : deep_copy
 use m_paw_an,    only : paw_an_type, paw_an_free, paw_an_redistribute
 use m_paw_ij,    only : paw_ij_type, paw_ij_free, paw_ij_redistribute
 use m_pawfgrtab, only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_redistribute
 use m_pawrhoij,  only : pawrhoij_type, pawrhoij_free, pawrhoij_redistribute
 use m_paral_atom,only : get_atm_proc
 use m_mpinfo,    only : initmpi_atom

 implicit none

 private

!public procedures.
 public :: set_pert_comm
 public :: unset_pert_comm
 public :: set_pert_paw
 public :: unset_pert_paw

!private procedures.
 private :: get_exchatom_list
 private :: get_exchatom_list1
!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_paral_pert/set_pert_comm
!! NAME
!! set_pert_comm
!!
!! FUNCTION
!! Set the MPI communicators over the perturbed cell
!!
!! INPUTS
!!  nppert=number of processes for the parallelization over perturbations
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      get_atm_proc,xmpi_bcast,xmpi_comm_translate_ranks
!!
!! SOURCE

subroutine set_pert_comm(mpi_enreg,nppert)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nppert
 type(MPI_type), intent(inout) :: mpi_enreg

! *************************************************************************

 if (mpi_enreg%paral_pert==0) return

!New communicators for the perturbed cell
 mpi_enreg%comm_cell =mpi_enreg%comm_cell_pert
 mpi_enreg%me_cell   =xmpi_comm_rank(mpi_enreg%comm_cell)
 mpi_enreg%nproc_cell=mpi_enreg%nproc/nppert

!Adjust other communicators
 mpi_enreg%comm_kpt    =mpi_enreg%comm_cell
 mpi_enreg%me_kpt      =mpi_enreg%me_cell
 mpi_enreg%comm_kptband=mpi_enreg%comm_cell
 mpi_enreg%comm_wvl    =mpi_enreg%comm_cell

end subroutine set_pert_comm
!!***

!----------------------------------------------------------------------

!!****f* m_paral_pert/unset_pert_comm
!! NAME
!! unset_pert_comm
!!
!! FUNCTION
!! Unset the MPI communicators over the perturbed cell; restore the global communicators.
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      get_atm_proc,xmpi_bcast,xmpi_comm_translate_ranks
!!
!! SOURCE

subroutine unset_pert_comm(mpi_enreg)

!Arguments ---------------------------------------------
!scalars
 type(MPI_type), intent(inout) :: mpi_enreg

! *************************************************************************

 if (mpi_enreg%paral_pert==0) return

!New communicators for the cell
 mpi_enreg%comm_cell =mpi_enreg%comm_world
 mpi_enreg%nproc_cell=mpi_enreg%nproc
 mpi_enreg%me_cell   =mpi_enreg%me

!Adjust other communicators
 mpi_enreg%comm_kpt    =mpi_enreg%comm_cell
 mpi_enreg%me_kpt      =mpi_enreg%me_cell
 mpi_enreg%comm_kptband=mpi_enreg%comm_cell
 mpi_enreg%comm_wvl    =mpi_enreg%comm_cell

end  subroutine unset_pert_comm
!!***

!----------------------------------------------------------------------

!!****f* m_paral_pert/set_pert_paw
!! NAME
!! set_pert_paw
!!
!! FUNCTION
!! Set PAW data for the parallelization over perturbations:
!!  Set MPI communicator over atomic sites
!!  Redistribute PAW on-site data
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!!  paw_an_out(my_natom)<type(paw_an_type)>= --optional--
!!              redistributed PAW arrays given on angular mesh
!!              if not present, paw_an is redistributed "in-place"
!!  paw__ij_out(my_natom)<type(paw_ij_type)>= --optional--
!!              redistributed PAW arrays given on (i,j) channels
!!              if not present, paw_ij is redistributed "in-place"
!!  pawfgrtab_out(my_natom)<type(pawfgrtab_type)>= --optional--
!!              redistributed PAW atomic data given on fine grid
!!              if not present, pawfgrtab is redistributed "in-place"
!!  pawrhoij_out(my_natom)<type(pawrhoij_type)>= --optional--
!!              redistributed PAW rhoij occupancies
!!              if not present, pawrhoij is redistributed "in-place"
!!  old_atmtab=save the indexes of the atoms treated by current proc
!!  old_comm_atom=save the identifier of the MPI communicator
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  paw_an(my_natom)<type(paw_an_type)>=PAW arrays given on angular mesh
!!  paw_ij(my_natom)<type(paw_ij_type)>=PAW arrays given on (i,j) channels
!!  pawfgrtab(my_natom)<type(pawfgrtab_type)>=PAW atomic data given on fine grid
!!  pawrhoij(my_natom)<type(pawrhoij_type)>=PAW rhoij occupancies
!!
!! PARENTS
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      get_atm_proc,xmpi_bcast,xmpi_comm_translate_ranks
!!
!! SOURCE

subroutine set_pert_paw(dtset,mpi_enreg,my_natom,old_atmtab,old_comm_atom,&
&                       paw_an,paw_ij,pawfgrtab,pawrhoij,&
&                       paw_an_out,paw_ij_out,pawfgrtab_out,pawrhoij_out)

!Arguments ---------------------------------------------
!scalars
 integer,intent(inout) :: my_natom
 integer,intent(inout) :: old_comm_atom
 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer,pointer,intent(out) :: old_atmtab(:)
 type(paw_ij_type),allocatable,target,intent(inout) :: paw_ij(:)
 type(paw_ij_type),optional,pointer,intent(inout) :: paw_ij_out(:)
 type(paw_an_type),allocatable,target,intent(inout) :: paw_an(:)
 type(paw_an_type),optional,pointer,intent(inout) :: paw_an_out(:)
 type(pawfgrtab_type),allocatable,target,intent(inout) :: pawfgrtab(:)
 type(pawfgrtab_type),optional,pointer,intent(inout) :: pawfgrtab_out(:)
 type(pawrhoij_type),allocatable,target,intent(inout) :: pawrhoij(:)
 type(pawrhoij_type),optional,pointer,intent(inout) :: pawrhoij_out(:)

!Local variables ---------------------------------------
!scalars
!Type of algo used for communications: 1-brute force, 2-asynchronous
 integer,parameter :: algo_option=2
 logical :: paral_atom
!arrays
 integer,allocatable :: SendAtomProc(:), SendAtomList(:)
 integer,allocatable :: RecvAtomProc(:), RecvAtomList(:)
 real(dp) :: tsec(2)

! *************************************************************************

 call timab(593,1,tsec)

 paral_atom=(dtset%natom/=my_natom)

!Default value when parallelization over atomic sites is not activated
 if ((mpi_enreg%paral_pert==0).or.(.not.paral_atom).or.(dtset%usepaw==0)) then
   old_comm_atom=mpi_enreg%comm_atom
   call deep_copy(mpi_enreg%my_atmtab,old_atmtab)
   if (present(pawrhoij_out))  pawrhoij_out=>pawrhoij
   if (present(paw_ij_out))    paw_ij_out=>paw_ij
   if (present(paw_an_out))    paw_an_out=>paw_an
   if (present(pawfgrtab_out)) pawfgrtab_out=>pawfgrtab
   return
 end if

!Redefine communicator for atoms and store old communicator

 old_comm_atom=mpi_enreg%comm_atom
 call deep_copy(mpi_enreg%my_atmtab,old_atmtab)
 call initmpi_atom(dtset,mpi_enreg)
 my_natom=mpi_enreg%my_natom

!If "asynchronous algo", determine lists of atoms to be exchanged
!  between me and other processes
 if (algo_option==2) then
   call timab(594,1,tsec)
   call get_exchatom_list(old_comm_atom,mpi_enreg%comm_atom,old_atmtab,&
&       mpi_enreg%my_atmtab,dtset%natom,SendAtomProc,&
&       SendAtomList,RecvAtomProc,RecvAtomList)
   call timab(594,2,tsec)
 end if

!Redistribute PAW on-site data

!pawrhoij datastructure
 call timab(595,1,tsec)
 if (present(pawrhoij_out)) then
   nullify(pawrhoij_out)
   if (algo_option==1) then
     call pawrhoij_redistribute(pawrhoij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,pawrhoij_out=pawrhoij_out)
   else if (algo_option==2) then
     call pawrhoij_redistribute(pawrhoij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,pawrhoij_out=pawrhoij_out, &
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 else
   if (algo_option==1) then
     call pawrhoij_redistribute(pawrhoij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call pawrhoij_redistribute(pawrhoij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(595,2,tsec)

!paw_ij datastructure
 call timab(596,1,tsec)
 if (present(paw_ij_out)) then
   nullify(paw_ij_out)
   if (algo_option==1) then
     call paw_ij_redistribute(paw_ij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,paw_ij_out=paw_ij_out)
   else if (algo_option==2) then
     call paw_ij_redistribute(paw_ij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,paw_ij_out=paw_ij_out,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 else
   if (algo_option==1) then
     call paw_ij_redistribute(paw_ij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call paw_ij_redistribute(paw_ij,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(596,2,tsec)

!paw_an datastructure
 call timab(597,1,tsec)
 if (present(paw_an_out)) then
   nullify(paw_an_out)
   if (algo_option==1) then
     call paw_an_redistribute(paw_an,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,paw_an_out=paw_an_out)
   else if (algo_option==2) then
     call paw_an_redistribute(paw_an,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,paw_an_out=paw_an_out,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 else
   if (algo_option==1) then
     call paw_an_redistribute(paw_an,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call paw_an_redistribute(paw_an,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(597,2,tsec)

!pawfgrtab datastructure
 call timab(598,1,tsec)
 if (present(pawfgrtab_out)) then
   nullify(pawfgrtab_out)
   if (algo_option==1) then
     call pawfgrtab_redistribute(pawfgrtab,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,pawfgrtab_out=pawfgrtab_out)
   else if (algo_option==2) then
     call pawfgrtab_redistribute(pawfgrtab,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,pawfgrtab_out=pawfgrtab_out,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 else
   if (algo_option==1) then
     call pawfgrtab_redistribute(pawfgrtab,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call pawfgrtab_redistribute(pawfgrtab,old_comm_atom,mpi_enreg%comm_atom,&
&         mpi_atmtab_in=old_atmtab,mpi_atmtab_out=mpi_enreg%my_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(598,2,tsec)

 if (algo_option==2) then
   ABI_FREE(SendAtomProc)
   ABI_FREE(SendAtomList)
   ABI_FREE(RecvAtomProc)
   ABI_FREE(RecvAtomList)
 end if

 call timab(593,2,tsec)

end  subroutine set_pert_paw
!!***

!----------------------------------------------------------------------

!!****f* m_paral_pert/unset_pert_paw
!! NAME
!! unset_pert_paw
!!
!! FUNCTION
!! Unset PAW data after the parallelization over perturbations:
!!  Restore MPI communicator over atomic sites
!!  Restore PAW on-site data
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  old_atmtab=index of atoms to restore
!!  old_comm_atom=MPI communicator to restore
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!
!! OUTPUT
!!  paw_an_out(my_natom)<type(paw_an_type)>= --optional--
!!              redistributed PAW arrays given on angular mesh
!!              if not present, paw_an is restored "in-place"
!!  paw__ij_out(my_natom)<type(paw_ij_type)>= --optional--
!!              redistributed PAW arrays given on (i,j) channels
!!              if not present, paw_ij is restored "in-place"
!!  pawfgrtab_out(my_natom)<type(pawfgrtab_type)>= --optional--
!!              redistributed PAW atomic data given on fine grid
!!              if not present, pawfgrtab is restored "in-place"
!!  pawrhoij_out(my_natom)<type(pawrhoij_type)>= --optional--
!!              redistributed PAW rhoij occupancies
!!              if not present, pawrhoij is restored "in-place"
!! old_atmtab=save the indexes of atoms treated by current proc
!! old_comm_atom=save the identifier of the MPI communicator
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  paw_an(my_natom)<type(paw_an_type)>=PAW arrays given on angular mesh
!!  paw_ij(my_natom)<type(paw_ij_type)>=PAW arrays given on (i,j) channels
!!  pawfgrtab(my_natom)<type(pawfgrtab_type)>=PAW atomic data given on fine grid
!!  pawrhoij(my_natom)<type(pawrhoij_type)>=PAW rhoij occupancies
!!
!! PARENTS
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      get_atm_proc,xmpi_bcast,xmpi_comm_translate_ranks
!!
!! SOURCE

subroutine unset_pert_paw(dtset,mpi_enreg,my_natom,old_atmtab,old_comm_atom,&
&                       paw_an,paw_ij,pawfgrtab,pawrhoij,&
&                       paw_an_out,paw_ij_out,pawfgrtab_out,pawrhoij_out)

!Arguments ---------------------------------------------
!scalars
 integer,intent(inout) :: my_natom
 integer,intent(inout) :: old_comm_atom
 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer,pointer,intent(inout) :: old_atmtab(:)
 type(paw_ij_type),allocatable,target,intent(inout) :: paw_ij(:)
 type(paw_ij_type),optional,pointer,intent(inout) :: paw_ij_out(:)
 type(paw_an_type),allocatable,target,intent(inout) :: paw_an(:)
 type(paw_an_type),optional,pointer,intent(inout) :: paw_an_out(:)
 type(pawfgrtab_type),allocatable,target,intent(inout) :: pawfgrtab(:)
 type(pawfgrtab_type),optional,pointer,intent(inout) :: pawfgrtab_out(:)
 type(pawrhoij_type),allocatable,target,intent(inout) :: pawrhoij(:)
 type(pawrhoij_type),optional,pointer,intent(inout) :: pawrhoij_out(:)

!Local variables ---------------------------------------
!scalars
!Type of algo used for communications: 1-brute force, 2-asynchronous
 integer,parameter :: algo_option=2
integer :: my_natom_old
 logical :: exchange,paral_atom
!arrays
 integer,allocatable :: SendAtomProc(:), SendAtomList(:)
 integer,allocatable :: RecvAtomProc(:), RecvAtomList(:)
 real(dp) :: tsec(2)

! *************************************************************************

!Nothing to do when parallelization over atomic sites is not activated
 if ((mpi_enreg%paral_pert==0).or.(dtset%usepaw==0)) return
 if (.not.associated(old_atmtab)) return
 paral_atom=(dtset%natom/=size(old_atmtab))
 if (.not.paral_atom) return

!If "asynchronous algo", determine lists of atoms to be exchanged
!  between me and other processes
 exchange=.true.
 if (present(pawrhoij_out).and.present(paw_ij_out).and. &
&    present(paw_an_out)  .and.present(pawfgrtab_out)) exchange=.false.
 if (algo_option==2.and.exchange) then
   call timab(594,1,tsec)
   my_natom_old=size(old_atmtab)
   call get_exchatom_list1(mpi_enreg%comm_atom,old_comm_atom,mpi_enreg%my_atmtab,&
&       old_atmtab,dtset%natom,SendAtomProc,&
&       SendAtomList,RecvAtomProc,RecvAtomList)
   call timab(594,2,tsec)
 end if

!Redistribute PAW on-site data (if in-place storage)
!or destroy them (out-of-place storage)

!pawrhoij datastructure
 call timab(595,1,tsec)
 if (present(pawrhoij_out)) then
   call pawrhoij_free(pawrhoij_out)
   ABI_FREE(pawrhoij_out)
 else
   if (algo_option==1) then
     call pawrhoij_redistribute(pawrhoij,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
       call pawrhoij_redistribute(pawrhoij,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(595,2,tsec)

!paw_ij datastructure
 call timab(596,1,tsec)
 if (present(paw_ij_out)) then
   call paw_ij_free(paw_ij_out)
   ABI_FREE(paw_ij_out)
 else
   if (algo_option==1) then
     call paw_ij_redistribute(paw_ij,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call paw_ij_redistribute(paw_ij,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(596,2,tsec)

!paw_an datastructure
 call timab(597,1,tsec)
 if (present(paw_an_out)) then
   call paw_an_free(paw_an_out)
   ABI_FREE(paw_an_out)
 else
   if (algo_option==1) then
     call paw_an_redistribute(paw_an,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call paw_an_redistribute(paw_an,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(597,2,tsec)

!pawfgrtab datastructure
 call timab(598,1,tsec)
 if (present(pawfgrtab_out)) then
   call pawfgrtab_free(pawfgrtab_out)
   ABI_FREE(pawfgrtab_out)
 else
   if (algo_option==1) then
     call pawfgrtab_redistribute(pawfgrtab,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom)
   else if (algo_option==2) then
     call pawfgrtab_redistribute(pawfgrtab,mpi_enreg%comm_atom,old_comm_atom,&
&         mpi_atmtab_in=mpi_enreg%my_atmtab,mpi_atmtab_out=old_atmtab,&
&         natom=dtset%natom,&
&         SendAtomproc=SendAtomProc,SendAtomList=SendAtomList,&
&         RecvAtomProc=RecvAtomProc,RecvAtomList=RecvAtomList)
   end if
 end if
 call timab(598,2,tsec)

!Release some memory
 if (algo_option==2.and.exchange) then
   ABI_FREE(SendAtomProc)
   ABI_FREE(SendAtomList)
   ABI_FREE(RecvAtomProc)
   ABI_FREE(RecvAtomList)
 end if

!Restore communicator for atoms
 ABI_FREE(old_atmtab)
 old_comm_atom=mpi_enreg%comm_atom
 call deep_copy(mpi_enreg%my_atmtab,old_atmtab)
 call initmpi_atom(dtset,mpi_enreg)
 my_natom=mpi_enreg%my_natom

end  subroutine unset_pert_paw
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_exchatom_list
!! NAME
!! get_exchatom_list
!!
!! FUNCTION
!! This routine determine the list of atoms to be exchanged by current process and other processes.
!! The initial communicator (mpicomm_in) should be a common ancestor of (mpicomm_in and mpicomm_out);
!! It is splitted in "nbgroup" final communicators (mpicomm_out).
!!
!! INPUTS
!!  mpicomm_in=input MPI communicator (over atoms) among which atoms are initially distributed
!!  mpicomm_out=MPI communicator (over atoms) among which we want to distribute atoms
!!  my_atmtab_in= Index of atoms initially treated by current proc (communicator mpicomm_in)
!!  my_atmtab_out= Index of atoms finally treated by current proc (communicator mpicomm_out)
!!  natom = total number of atoms
!!  nbgroup = # of mpicomm_out communicators included in mpicomm_in
!!
!! OUTPUT
!! For current proc :
!! RecvAtomProc(:)= rank of processor from which I expect atom (in mpicomm_in)
!! RecvAtomList(:)= indexes of atoms to be received by me
!!   RecvAtomList(irecv) are the atoms I expect from RecvAtomProc(irecv)
!! SendAtomProc(:)= ranks of process destination of atom (in mpicomm_in)
!! SendAtomList(:)= indexes of atoms to be sent by me
!!   SendAtomList(isend) are the atoms sent to SendAtomProc(isend)
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!      get_atm_proc,xmpi_bcast,xmpi_comm_translate_ranks
!!
!! SOURCE


 subroutine get_exchatom_list(mpicomm_in,mpicomm_out,my_atmtab_in,my_atmtab_out,natom, &
&                             SendAtomProc,SendAtomList,RecvAtomProc,RecvAtomList)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpicomm_in,mpicomm_out,natom
!arrays
 integer,intent(in) :: my_atmtab_in(:),my_atmtab_out(:)
 integer,allocatable,intent(out) :: RecvAtomProc(:),RecvAtomList(:)
 integer,allocatable,intent(out) :: SendAtomList(:),SendAtomProc(:)

!Local variables ---------------------------------------
!scalars
 integer :: igroup,i1,ierr,ii,me_in,me_in_out,me_out,me_out0_in
 integer :: my_natom_in,my_natom_out,nbgroup,nbsend,nproc_in,nproc_max,nproc_out
!arrays
 integer :: buf_int(3),rank0(1),ranks0_out_in(1)
 integer, allocatable :: buf_int_all(:),group(:),master_commout(:),procs(:),ranks(:)
 integer, allocatable  :: ranks1(:,:),sizecomm(:)

! *************************************************************************

 nproc_in=xmpi_comm_size(mpicomm_in)
 me_in=xmpi_comm_rank(mpicomm_in)
 my_natom_in=size(my_atmtab_in)
 nproc_out=xmpi_comm_size(mpicomm_out)
 me_out=xmpi_comm_rank(mpicomm_out)
 my_natom_out=size(my_atmtab_out)

 rank0(1)=0
 call xmpi_comm_translate_ranks(mpicomm_out,1,rank0(1),mpicomm_in,ranks0_out_in(1))
 me_out0_in=ranks0_out_in(1)

 ABI_MALLOC(ranks,(1:nproc_in))
 ABI_MALLOC(sizecomm,(1:nproc_in))
 ABI_MALLOC(master_commout,(1:nproc_in))
 ABI_MALLOC(group,(0:nproc_in-1))
 buf_int(1)=me_out; buf_int(2)=me_out0_in; buf_int(3)=nproc_out;

 ABI_MALLOC(buf_int_all,(3*nproc_in))
 call xmpi_allgather(buf_int,3,buf_int_all,mpicomm_in,ierr)
 nbgroup=0;
 nproc_max=0
 ranks(:)=-1;group(:)=-1;sizecomm(:)=-1;master_commout(:)=-1

 do ii=1,nproc_in
   ranks(ii)=buf_int_all(3*ii-2)    !me_out
   master_commout(ii)=buf_int_all(3*ii-1) !rank of me_out=0 of mpicomm_out expressed in mpicomm_in

   if (ranks(ii)==0) then
     nbgroup=nbgroup+1
     sizecomm(nbgroup)=buf_int_all(3*ii) !nproc_out
     group(master_commout(ii))=nbgroup
     if (sizecomm(nbgroup)>nproc_max) nproc_max=sizecomm(nbgroup)
   end if

 enddo
 ABI_FREE(buf_int_all)


 ABI_MALLOC(ranks1,(nbgroup,0:nproc_max-1))
 ranks1(:,:)=-1
 do ii=1,nproc_in
   me_in_out=ranks(ii)
   igroup=group(master_commout(ii))
   ranks1(igroup,me_in_out)=ii-1 !numbering of procs
 enddo

!Send
 nbsend=0
 if (my_natom_in>0) then
   ABI_MALLOC(SendAtomList,(my_natom_in*nbgroup))
   ABI_MALLOC(SendAtomProc,(my_natom_in*nbgroup))
   ABI_MALLOC(procs,(my_natom_in))
   do igroup=1,nbgroup
     call get_atm_proc(my_atmtab_in,natom,sizecomm(igroup),procs)
     do i1=1,my_natom_in
       nbsend=nbsend+1
       SendAtomProc(nbsend)=ranks1(igroup,procs(i1))
       SendAtomList(nbsend)=my_atmtab_in(i1)
     end do
   end do
   ABI_FREE(procs)
 else
  ABI_MALLOC(SendAtomList,(0))
  ABI_MALLOC(SendAtomProc,(0))
end if

!recv
 if (my_natom_out>0) then !no return before because of xmpi_allgather and of the sending operation
   ABI_MALLOC(RecvAtomProc,(my_natom_out))
   ABI_MALLOC(RecvAtomList,(my_natom_out))
   RecvAtomList(:)=my_atmtab_out(:)
!the atoms are put in increasing order,see get_my_atmtab so the procs are sorted by growing process
   call get_atm_proc(RecvAtomList,natom,nproc_in,RecvAtomProc)
 else
   ABI_MALLOC(RecvAtomList,(0))
   ABI_MALLOC(RecvAtomProc,(0))
 end if

 ABI_FREE(master_commout)
 ABI_FREE(ranks1)
 ABI_FREE(group)
 ABI_FREE(ranks)
 ABI_FREE(sizecomm)

end subroutine get_exchatom_list
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_exchatom_list1
!! NAME
!! get_exchatom_list1
!!
!! FUNCTION
!! This routine determine the list of atoms to be exchanged by current process and other processes.
!! this function redistribute the atoms of one mpicomm_in among mpicomm_out
!!
!! INPUTS
!!  mpicomm_in=input MPI communicator (over atoms) among which atoms are initially distributed
!!  mpicomm_out=MPI communicator (over atoms) among which we want to distribute atoms
!!  my_atmtab_in= Index of atoms initially treated by current proc (communicator mpicomm_in)
!!  my_atmtab_out= Index of atoms finally treated by current proc (communicator mpicomm_out)
!!  natom = total number of atoms
!!
!! OUTPUT
!! For current proc:
!! RecvAtomProc(:)= rank of processor from which I expect atom (in mpicomm_out)
!! RecvAtomList(:)= indexes of atoms to be received by me
!!   RecvAtomList(irecv) are the atoms I expect from RecvAtomProc(irecv)
!! SendAtomProc(:)= ranks of process destination of atom (in mpicomm_out)
!! SendAtomList(:)= indexes of atoms to be sent by me
!!   SendAtomList(isend) are the atoms sent to SendAtomProc(isend)
!!
!! NOTES
!! Previously mpicomm_out has be split in several mpicomm_in.
!! In our purpose, we only need to redistribute the atoms of one mpicomm_in among mpicomm_out
!! because all structures of atoms we have to exchange have the same value in each mpicomm_in.
!! The mpicomm_in we choose is the one in which the processor 0 of mpicomm_out belong
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!      get_atm_proc,xmpi_bcast,xmpi_comm_translate_ranks
!!
!! SOURCE

subroutine get_exchatom_list1(mpicomm_in,mpicomm_out,my_atmtab_in,my_atmtab_out,natom, &
&                             SendAtomProc,SendAtomList,RecvAtomProc,RecvAtomList)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpicomm_in,mpicomm_out,natom
!arrays
 integer,intent(in) :: my_atmtab_in(:),my_atmtab_out(:)
 integer,allocatable,intent(out) :: RecvAtomProc(:),RecvAtomList(:)
 integer,allocatable,intent(out) :: SendAtomList(:),SendAtomProc(:)

!Local variables ---------------------------------------
!scalars
 integer :: i1,ier,me_out,my_natom_in,my_natom_out,nproc_in,nproc_out
 logical :: sender
!arrays
 integer,allocatable :: ranks_in(:),ranks_in_out(:)

! *************************************************************************

 me_out=xmpi_comm_rank(mpicomm_out)
 my_natom_in=size(my_atmtab_in)
 nproc_out=xmpi_comm_size(mpicomm_out)
 my_natom_out=size(my_atmtab_out)
 nproc_in=xmpi_comm_size(mpicomm_in)

 call xmpi_bcast(nproc_in,0,mpicomm_out,ier)
 ABI_MALLOC(ranks_in_out,(0:nproc_in-1))

!All atoms are distributed among each mpicomm_in
!redistribute the atoms of one mpicomm_in among mpicomm_out

!Look for the communicator mpicomm_in from which me_out=0 belong
!Get ranks of all processors of mpicomm_in expressed in mpicomm_out
 if (me_out==0) then
   ABI_MALLOC(ranks_in,(0:nproc_in-1))
   ranks_in=(/ (i1,i1=0,nproc_in-1 )/)
   call xmpi_comm_translate_ranks(mpicomm_in,nproc_in,ranks_in,mpicomm_out,ranks_in_out)
   ABI_FREE(ranks_in)
 end if
 call xmpi_bcast(ranks_in_out,0,mpicomm_out,ier)

!Check if me_out is one of the sending proc.
!(ie belongs to the mpicomm_in from which me_out belong)
 sender = .false.
 do i1=0,nproc_in-1
   if (me_out==ranks_in_out(i1)) then
     sender = .true.
     exit
   end if
 end do

!Send
 if (my_natom_in>0.and.sender) then
   ABI_MALLOC(SendAtomList,(my_natom_in))
   ABI_MALLOC(SendAtomProc,(my_natom_in))
   SendAtomList(:)=my_atmtab_in(:)
!  The atoms are put in increasing order,see get_my_atmtab
!  so the procs are sorted by growing process
   call get_atm_proc(SendAtomList,natom,nproc_out,SendAtomProc)
 else
   ABI_MALLOC(SendAtomList,(0))
   ABI_MALLOC(SendAtomProc,(0))
 end if

!Recv
 if (my_natom_out>0) then
   ABI_MALLOC(RecvAtomProc,(my_natom_out))
   ABI_MALLOC(RecvAtomList,(my_natom_out))
   RecvAtomList(:)=my_atmtab_out(:)
!  The atoms are put in increasing order,see get_my_atmtab
!  so the procs are sorted by growing process
   call get_atm_proc(RecvAtomList,natom,nproc_in,RecvAtomProc)
   RecvAtomProc(:)=ranks_in_out(RecvAtomProc(:))
 else
   ABI_MALLOC(RecvAtomList,(0))
   ABI_MALLOC(RecvAtomProc,(0))
 end if

 ABI_FREE(ranks_in_out)

end subroutine get_exchatom_list1
!!***

!----------------------------------------------------------------------

END MODULE m_paral_pert
!!***
