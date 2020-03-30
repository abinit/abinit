!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paral_atom
!! NAME
!!  m_paral_atom
!!
!! FUNCTION
!!  This module provides routines and methods used to manage the parallelisation/distribution
!!  of PAW data over atomic sites
!!
!! COPYRIGHT
!! Copyright (C) 2012-2020 ABINIT group (MT, MD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paral_atom

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 implicit none

 private

!public procedures.
 public :: get_my_natom
 public :: get_my_atmtab
 public :: free_my_atmtab
 public :: get_proc_atmtab
 public :: get_atm_proc
!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_my_natom
!! NAME
!! get_my_natom
!!
!! FUNCTION
!! Given the total number of atoms, return the number of atoms treated by current process
!!
!! INPUTS
!!  comm_atom=communicator over atoms
!!  natom=total number of atoms
!!
!! OUTPUT
!!  my_natom=number of atoms treated by current process
!!
!! PARENTS
!!      initmpi_atom,m_paw_an,m_paw_ij,m_pawfgrtab,m_pawrhoij
!!
!! CHILDREN
!!
!! SOURCE


subroutine get_my_natom(comm_atom,my_natom,natom)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm_atom,natom
 integer,intent(out) :: my_natom
!arrays

!Local variables ---------------------------------------
!scalars
 integer ::  me,nproc
!arrays

! *************************************************************************

 my_natom=natom
 if (comm_atom/=xmpi_comm_self.and.comm_atom/=xmpi_comm_null)  then
   nproc=xmpi_comm_size(comm_atom)
   me=xmpi_comm_rank(comm_atom)
   my_natom=natom/nproc
   if (me<=(mod(natom,nproc)-1)) my_natom=natom/nproc + 1
 endif

end subroutine get_my_natom
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_my_atmtab
!! NAME
!! get_my_atmtab
!!
!! FUNCTION
!! Given the total number of atoms and a MPI communicator return a table
!! containing the indexes of atoms treated by current processor.
!!
!! INPUTS
!!  comm_atom=communicator over atoms
!!  my_natom_ref= --optional-- a reference value for the local number of atoms
!!                (just for checking purposes)
!!  natom=total number of atoms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  my_atmtab(:)=indexes of atoms treated by current process
!!               nothing is done if my_atmtab(:) is already associated to a target
!!  my_atmtab_allocated=true if my_atmtab is allocated
!!  paral_atom=flag controlling parallelism over atoms
!!
!! PARENTS
!!      calc_efg,calc_fc,denfgr,dfpt_accrho,dfpt_ewald,elt_ewald,eltxccore
!!      initmpi_atom,initrhoij,m_hamiltonian,m_paw_an,m_paw_ij,m_paw_pwaves_lmn
!!      m_pawdij,m_pawfgrtab,m_pawrhoij,make_efg_onsite,make_fc_paw,newfermie1
!!      nhatgrid,outscfcv,paw_mknewh0,pawaccrhoij,pawdenpot,pawdfptenergy
!!      pawgrnl,pawmkaewf,pawmknhat,pawmknhat_psipsi,pawnhatfr,pawprt,pawsushat
!!      pawuj_red,setnoccmmp,setrhoijpbe0
!!
!! CHILDREN
!!
!! SOURCE


subroutine get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
&                        my_natom_ref) ! optional argument

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm_atom,natom
 integer,intent(in),optional :: my_natom_ref
 logical,intent(inout) :: my_atmtab_allocated,paral_atom
!arrays
 integer,pointer :: my_atmtab(:)
!Local variables ---------------------------------------
!scalars
 integer :: iatom,me,my_natom,natom_bef,nmod,nproc
 character(len=500) :: msg
!arrays

! *************************************************************************

 my_atmtab_allocated=.false.
 if (.not.paral_atom) return

 if (comm_atom==xmpi_comm_self.or.comm_atom==xmpi_comm_null) paral_atom=.false.
#ifdef DEV_MJV
print *, 'get_my_atmtab, paral_atom ', paral_atom
#endif
 if (paral_atom)  then
   nproc=xmpi_comm_size(comm_atom)
   paral_atom=(nproc>1)
   if (paral_atom) then
     if (.not.associated(my_atmtab)) then
!      Get local number of atoms
       me=xmpi_comm_rank(comm_atom)
       my_natom=natom/nproc
       if (me<=(mod(natom,nproc)-1)) my_natom=natom/nproc + 1
!      Get table of atoms
       if (my_natom>0) then
         LIBPAW_POINTER_ALLOCATE(my_atmtab,(my_natom))
         my_atmtab_allocated=.true.
         if (my_natom==natom) then
           my_atmtab(1:my_natom)=(/(iatom,iatom=1,natom)/)
         else
!          The atoms are distributed contigously by egal part
!          (the rest is distributed on all the procs)
           nmod=mod(natom,nproc)
           if (me<=(nmod-1)) then
             natom_bef=me*(natom/nproc)+me
           else
             natom_bef=me*(natom/nproc)+nmod
           endif
           do iatom=1,my_natom
             my_atmtab(iatom)=iatom+natom_bef
           enddo
         end if
#ifdef DEV_MJV
print *, 'get_my_atmtab, my_atmtab, natom_bef, nmod, nproc, me ', &
                         my_atmtab, natom_bef, nmod, nproc, me
#endif
       end if
     else
       my_natom=size(my_atmtab)
     end if
     if (present(my_natom_ref).and.(my_natom>0)) then
       if (my_natom_ref/=size(my_atmtab)) then
         msg='my_atmtab should have a size equal to my_natom !'
         MSG_BUG(msg)
       end if
     end if
   end if
 endif
#ifdef DEV_MJV
print *, 'get_my_atmtab, my_natom ', my_natom
#endif

end subroutine get_my_atmtab
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/free_my_atmtab
!! NAME
!! free_my_atmtab
!!
!! FUNCTION
!! Cleanly deallocate a table of atom indexes (my_atmtab)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  my_atmtab_allocated=true if my_atmtab is allocated
!!  my_atmtab(:)=indexes of atoms treated by current process
!!               nothing is done if my_atmtab(:) is already associated to a target
!!
!! PARENTS
!!      calc_efg,calc_fc,denfgr,dfpt_ewald,elt_ewald,eltxccore,initrhoij
!!      m_hamiltonian,m_paw_an,m_paw_ij,m_paw_pwaves_lmn,m_pawdij,m_pawfgrtab
!!      m_pawrhoij,make_efg_onsite,make_fc_paw,newfermie1,nhatgrid,outscfcv
!!      paw_mknewh0,pawaccrhoij,pawdenpot,pawdfptenergy,pawgrnl,pawmkaewf
!!      pawmknhat,pawmknhat_psipsi,pawnhatfr,pawprt,pawsushat,pawuj_red
!!      setnoccmmp,setrhoijpbe0
!!
!! CHILDREN
!!
!! SOURCE

subroutine free_my_atmtab(my_atmtab,my_atmtab_allocated)

!Arguments ---------------------------------------------
!scalars
 logical,intent(inout) :: my_atmtab_allocated
!arrays
 integer,pointer :: my_atmtab(:)
!Local variables ---------------------------------------
!scalars
!arrays

! *************************************************************************

 if (my_atmtab_allocated) then
   LIBPAW_POINTER_DEALLOCATE(my_atmtab)
   nullify(my_atmtab)
   my_atmtab_allocated=.false.
 end if

end subroutine free_my_atmtab
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_proc_atmtab
!! NAME
!! get_proc_atmtab
!!
!! FUNCTION
!!  Given the total number of atoms and the size of a communicator,
!!  return a table containing the indexes of atoms treated by a processor (with index iproc)
!!
!! INPUTS
!!  comm_atom_size= size of communicator (over atoms)
!!  iproc= rank of the processor
!!  natom= total number of atoms
!!
!! OUTPUT
!!  atmtab(natom_out)= indexes of atoms treated by process iproc
!!  natom_out= number of atoms treated by process iproc
!!
!! NOTES
!! In case of modification of the distribution of atom over proc,
!!   get_atmtab must be modified accordingly
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine get_proc_atmtab(iproc,atmtab,natom_out,natom,comm_atom_size)

!Arguments ---------------------------------------------
!scalars
 integer, intent(in) :: comm_atom_size,iproc,natom
!arrays
 integer,intent(out) :: natom_out
 integer,allocatable, intent(out):: atmtab(:)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,natom_bef,nmod,nproc
!arrays

! *************************************************************************

 nproc=comm_atom_size

 natom_out=natom/nproc ; if (iproc<=(mod(natom,nproc)-1)) natom_out=natom/nproc+1

! Get table of atoms
 if (natom_out>0) then
   LIBPAW_ALLOCATE(atmtab,(natom_out))
!  The atoms are distributed contigously by egal part
!  The rest is distributed on all the procs
!  (see get_my_atmtab)
   nmod=mod(natom,nproc)
   if (iproc<=(nmod-1)) then
     natom_bef=iproc*(natom/nproc)+iproc
   else
     natom_bef=iproc*(natom/nproc)+nmod
   end if
   do iatom=1,natom_out
     atmtab(iatom)=iatom+natom_bef
   end do

 else
   natom_out=0
   LIBPAW_ALLOCATE(atmtab,(0))
 end if

end subroutine get_proc_atmtab
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_atm_proc
!! NAME
!! get_atm_proc
!!
!! FUNCTION
!!  Given a list of atoms and a MPI communicator size, return a table
!!  containing the corresponding processor indexes.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2020 ABINIT group (MD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atom_list(:)= index of atoms
!!  nproc=size of communicator over atoms
!!  natom=total number of atoms
!!
!! OUTPUT
!! proc_list(:) = index of procs
!!
!! NOTES
!!  The atoms are distributed contigously by egal part; the rest is distributed
!!  on all the procs (see get_my_atmtab).
!!  In case of modification of the distribution of atom over proc, this routine
!!  must be modified accordingly.
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!
!! SOURCE

 subroutine get_atm_proc(atom_list,natom,nproc,proc_list)

!Arguments ---------------------------------------------
!scalars
 integer, intent(in) :: natom,nproc
!arrays
 integer, intent(in) :: atom_list(:)
 integer, intent(out) :: proc_list(:)

!Local variables ---------------------------------------
!scalars
 integer :: nb_atom,dn,dn1,iatom,natomlim,iatm,jproclim,nmod
!arrays

! *************************************************************************

 nmod=mod(natom,nproc);nb_atom=size(atom_list)

 if (nmod==0) then
   dn=natom/nproc
   do iatm =1, nb_atom
     iatom=atom_list(iatm)
     proc_list(iatm)=(iatom-1)/dn
   end do
 else
   dn=natom/nproc
   dn1=natom/nproc + 1
!  Under (jproclim+1), 1 atome by proc is added
!  The rest nmod is distributed among jproclim+1 first procs
   jproclim=nmod -1
   natomlim=dn1*(jproclim+1)
   do iatm=1,nb_atom
     iatom=atom_list(iatm)
     if (iatom<=natomlim) then
       proc_list(iatm)=(iatom -1 )/dn1
     else
       proc_list(iatm)=jproclim + 1 + (iatom - 1 -(natomlim))/dn
     end if
   enddo
 end if

end subroutine get_atm_proc
!!***

!----------------------------------------------------------------------

END MODULE m_paral_atom
!!***
